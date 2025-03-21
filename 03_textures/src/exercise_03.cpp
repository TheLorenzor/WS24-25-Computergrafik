#include <cglib/rt/intersection.h>
#include <cglib/rt/object.h>
#include <cglib/rt/ray.h>
#include <cglib/rt/raytracing_context.h>
#include <cglib/rt/render_data.h>
#include <cglib/rt/renderer.h>
#include <cglib/rt/texture.h>
#include <cglib/rt/texture_mapping.h>
#include <cglib/core/thread_local_data.h>

#include <cglib/core/image.h>
#include <cglib/core/glmstream.h>
#include <cglib/core/assert.h>

#include <algorithm>


// -----------------------------------------------------------------------------

/*
 * Evaluates a texture for the given uv-coordinate without filtering.
 *
 * This method transformes the uv-coordinate into a st-coordinate and
 * rounds down to integer pixel values.
 *
 * The parameter level in [0, mip_levels.size()-1] is the miplevel of
 * the texture that is to be evaluated.
 */
glm::vec4 ImageTexture::
evaluate_nearest(int level, glm::vec2 const& uv) const
{
	cg_assert(level >= 0 && level < static_cast<int>(mip_levels.size()));
	cg_assert(mip_levels[level]);

	// TODO: compute the st-coordinates for the given uv-coordinates and mipmap level
	int width = mip_levels[level]->getWidth();
    int height = mip_levels[level]->getHeight();
    float s = width * uv.x;
    float t = height * uv.y;

	// get the value of pixel (s, t) of miplevel level
	return get_texel(level, (int)s, (int)t);
}

// -----------------------------------------------------------------------------

/*
 * Implement clamping here.
 *
 * The input "val" may be arbitrary, including negative and very large positive values.
 * The method shall always return a value in [0, size).
 * Out-of-bounds inputs must be clamped to the nearest boundary.
 */
int ImageTexture::
wrap_clamp(int val, int size)
{
	cg_assert(size > 0);
    if (val < 0) {
        return 0;
    }
    if (val > size-1) {
      return size-1;
    }
	return val;
}

/*
 * Implement repeating here.
 *
 * The input "val" may be arbitrary, including negative and very large positive values.
 * The method shall always return a value in [0, size).
 * Out-of-bounds inputs must be mapped back into [0, size) so that 
 * the texture repeats infinitely.
 */
int ImageTexture::
wrap_repeat(int val, int size)
{
	cg_assert(size > 0);
    int returnVal = val%size;

	if (returnVal < 0) {
		returnVal = returnVal+size;
	}
	return returnVal;
}


// -----------------------------------------------------------------------------


/*
 * Implement bilinear filtering here.
 *
 * Use mip_levels[level] as the image to filter.
 * The input uv coordinates may be arbitrary, including values outside [0, 1).
 *
 * Callers of this method must guarantee that level is valid, and
 * that mip_levels[level] is properly initialized. Compare the assertions.
 *
 * The color of a texel is to be interpreted as the color at the texel center.
 */
glm::vec4 ImageTexture::
evaluate_bilinear(int level, glm::vec2 const& uv) const
{
	cg_assert(level >= 0 && level < static_cast<int>(mip_levels.size()));
	cg_assert(mip_levels[level]);

   	int width = mip_levels[level]->getWidth();
    int height = mip_levels[level]->getHeight();

    float s = width * uv.x;
    float t = height * uv.y;
	int sInt = std::floor(s);
	int tInt = std::floor(t);

	glm::vec4 t3 =  get_texel(level, sInt, tInt);
	glm::vec4 t4 =  get_texel(level, sInt+1, tInt);

	glm::vec4 t1 = get_texel(level, sInt, tInt+1);
	glm::vec4 t2 = get_texel(level, sInt+1, tInt+1);

	float a = std::fabs(s-sInt);
	float b = std::fabs(tInt+1-t);

	return t1*(1-a)*(1-b) + t2*a*(1-b) + t3*(1-a)*b + t4 * a*b;

}

// -----------------------------------------------------------------------------

/*
 * This method creates a mipmap hierarchy for
 * the texture.
 *
 * This is done by iteratively reducing the
 * dimenison of a mipmap level and averaging
 * pixel values until the size of a mipmap
 * level is [1, 1].
 *
 * The original data of the texture is stored
 * in mip_levels[0].
 *
 * You can allocale memory for a new mipmap level 
 * with dimensions (size_x, size_y) using
 *		mip_levels.emplace_back(new Image(size_x, size_y));
 */
void ImageTexture::
create_mipmap()
{
	/* this are the dimensions of the original texture/image */
	int size_x = mip_levels[0]->getWidth();
	int size_y = mip_levels[0]->getHeight();

	cg_assert("must be power of two" && !(size_x & (size_x - 1)));
	cg_assert("must be power of two" && !(size_y & (size_y - 1)));

    int i =0;
    float times = 0.0f;
    if (size_x > size_y) {
		times = std::log2(size_y);
    } else {
    	times = std::log2(size_x);
    }

    for (int y = 0; y < int(times); y++) {
    	size_x =size_x / 2;
    	size_y =size_y / 2;
        i++;
    	mip_levels.emplace_back(new Image(size_x, size_y));
    	for (int j =0; j < size_x; j++) {
    		for (int k = 0; k < size_y; k++) {
				glm::vec4 pixel_avg = (
                	mip_levels[i-1]->getPixel((j*2),(k*2))+
                	mip_levels[i-1]->getPixel((j*2)+1,(k*2))+
                	mip_levels[i-1]->getPixel((j*2),(k*2)+1)+
                	mip_levels[i-1]->getPixel((j*2)+1,(k*2)+1));
            	pixel_avg /= 4;
				mip_levels[i]->setPixel(j, k,  pixel_avg);
        	}
		}
    }
}

/*
 * Compute the dimensions of the pixel footprint's AABB in uv-space.
 *
 * First intersect the four rays through the pixel corners with
 * the tangent plane at the given intersection.
 *
 * Then the given code computes uv-coordinates for these
 * intersection points.
 *
 * Finally use the uv-coordinates and compute the AABB in
 * uv-space. 
 *
 * Return width (du) and height (dv) of the AABB.
 *
 */
glm::vec2 Object::
compute_uv_aabb_size(const Ray rays[4], Intersection const& isect)
{
	glm::vec3 intersection_positions[4] = {
		isect.position, isect.position, isect.position, isect.position
	};

	for (int i = 0; i < 4; ++i) {
        Intersection temp_isect;
        if (geo->intersect(rays[i], &temp_isect)) {
            intersection_positions[i] = temp_isect.position;
        }
	}

	// compute uv coordinates from intersection positions
	glm::vec2 intersection_uvs[4];
	get_intersection_uvs(intersection_positions, isect, intersection_uvs);


    glm::vec2 uv_min = intersection_uvs[0];
    glm::vec2 uv_max = intersection_uvs[0];
    for (int i = 1; i < 4; ++i) {
        uv_min = glm::min(uv_min, intersection_uvs[i]);
        uv_max = glm::max(uv_max, intersection_uvs[i]);
    }

	return uv_max - uv_min;
}

/*
 * Implement trilinear filtering at a given uv-position.
 *
 * Transform the AABB dimensions dudv in st-space and
 * take the maximal coordinate as the 1D footprint size T.
 *
 * Determine mipmap levels i and i+1 such that
 *		texel_size(i) <= T <= texel_size(i+1)
 *
 *	Hint: use std::log2(T) for that.
 *
 *	Perform bilinear filtering in both mipmap levels and
 *	linearly interpolate the resulting values.
 *
 */
glm::vec4 ImageTexture::
evaluate_trilinear(glm::vec2 const& uv, glm::vec2 const& dudv) const
{
  	float width_s = dudv.x * this->mip_levels[0]->getWidth();
    float height_t = dudv.y * this->mip_levels[0]->getHeight();
    float T = std::log2(std::max(width_s, height_t));
	if ( T > mip_levels.size() - 1) {
          return evaluate_bilinear(mip_levels.size() - 1, uv);
	} else if (T < 0) {
    	return evaluate_bilinear(0, uv);
    } else {
    	int t1= int(T);
        int t2 = int(T)+1;
        glm::vec4 vec1 = evaluate_bilinear(t1, uv);
        glm::vec4 vec2 = evaluate_bilinear(t2, uv);
        T = T-int(T);
        return vec1*(1-T)+T*vec2;
    }
}

// -----------------------------------------------------------------------------

/*
 * Transform the given direction d using the matrix transform.
 *
 * The output direction must be normalized, even if d is not.
 */
glm::vec3 transform_direction(glm::mat4 const& transform, glm::vec3 const& d)
{
	glm::vec4 v4 = glm::vec4(d,0.f);
    v4 =transform * v4;
    glm::vec3 newVec = glm::vec3(v4[0],v4[1],v4[2]);
	return glm::normalize(newVec);
}

/*
 * Transform the given position p using the matrix transform.
 */
glm::vec3 transform_position(glm::mat4 const& transform, glm::vec3 const& p)
{
	glm::vec4 v4 = glm::vec4(p,1);
	v4 =transform * v4;
	glm::vec3 newVec = glm::vec3(v4[0],v4[1],v4[2]);
	return newVec;
}

/*
 * Intersect with the ray, but do so in object space.
 *
 * First, transform ray into object space. Use the methods you have
 * implemented for this.
 * Then, intersect the object with the transformed ray.
 * Finally, make sure you transform the intersection back into world space.
 *
 * isect is guaranteed to be a valid pointer.
 * The method shall return true if an intersection was found and false otherwise.
 *
 * isect must be filled properly if an intersection was found.
 */
bool Object::
intersect(Ray const& ray, Intersection* isect) const
{
	cg_assert(isect);

	if (RaytracingContext::get_active()->params.transform_objects) {

		// TODO: transform ray, intersect object, transform intersection
    	Ray newRay = transform_ray(ray,this->transform_world_to_object);
		Intersection isectLoc;
        if (geo->intersect(newRay, &isectLoc)) {
        	*isect = transform_intersection(isectLoc,this->transform_object_to_world,this->transform_object_to_world_normal);
			isect->t = glm::length(isect->position - ray.origin);
        	return true;
        }
        return false;
	}
	return geo->intersect(ray, isect);
}

/*
 * Transform a direction from tangent space to object space.
 *
 * Tangent space is a right-handed coordinate system where
 * the tangent is your thumb, the normal is the index finger, and
 * the bitangent is the middle finger.
 *
 * normal, tangent, and bitangent are given in object space.
 * Build a matrix that rotates d from tangent space into object space.
 * Then, transform d with this matrix to obtain the result.
 * 
 * You may assume that normal, tangent, and bitangent are normalized
 * to length 1.
 *
 * The output vector must be normalized to length 1, even if d is not.
 */
glm::vec3 transform_direction_to_object_space(
	glm::vec3 const& d, 
	glm::vec3 const& normal, 
	glm::vec3 const& tangent, 
	glm::vec3 const& bitangent)
{
	cg_assert(std::fabs(glm::length(normal)    - 1.0f) < 1e-4f);
	cg_assert(std::fabs(glm::length(tangent)   - 1.0f) < 1e-4f);
	cg_assert(std::fabs(glm::length(bitangent) - 1.0f) < 1e-4f);

	glm::mat3 rotation_Matrix = glm::mat3(tangent,normal,bitangent);
    glm::vec3 newDir = rotation_Matrix * d;
    newDir = glm::normalize(newDir);

	return newDir;
}

// -----------------------------------------------------------------------------

// Optional bonus assignment! The fourier scene uses evaluate_bilinear(..)
// for displaying amplitude, phase and result images so you need a working
// solution for the bilinear filtering!
//
// Spectrum and reconstruction are 2D memory blocks of the same size [sx, sy]
// in row major layout, i.e one row of size sx after the other. 
//
// Spectrum contains the complex fourier-coefficients \hat{x}_{kl}.
//
// Reconstruct the original grayscale image using fourier transformation!
void DiscreteFourier2D::reconstruct(
	int M, int N,
	std::complex<float> const* spectrum,
	std::complex<float>	* reconstruction)
{
}
// CG_REVISION 8c58412a25ac2367c053bfa840ee81b320bdd315
