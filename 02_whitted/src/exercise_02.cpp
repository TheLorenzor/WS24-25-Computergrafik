#include <cglib/rt/renderer.h>
#include <cglib/rt/intersection_tests.h>
#include <cglib/rt/raytracing_context.h>
#include <cglib/rt/intersection.h>
#include <cglib/rt/ray.h>
#include <cglib/rt/scene.h>
#include <cglib/rt/light.h>
#include <cglib/rt/material.h>
#include <cglib/rt/render_data.h>
#include <math.h>

/*
 * The sphere is defined by its center and radius.
 *
 * Return true, if (and only if) the ray intersects the sphere.
 * In this case, also fill the parameter t with the distance such that
 *    ray_origin + t * ray_direction
 * is the intersection point.
 */
bool intersect_sphere(
    glm::vec3 const& ray_origin,    // starting point of the ray
    glm::vec3 const& ray_direction, // direction of the ray
    glm::vec3 const& center,        // position of the sphere
    float radius,                   // radius of the sphere
    float* t)                       // output parameter which contains distance to the hit point
{
    cg_assert(t);
	cg_assert(std::fabs(glm::length(ray_direction) - 1.f) < EPSILON);

	float a = glm::dot(ray_direction, ray_direction);
	float b = glm::dot(2.f * ray_direction,ray_origin - center);
	float c = glm::dot(ray_origin - center, ray_origin - center)-radius*radius;

	float discr = b*b - 4*a*c;
	if (discr < 0) {
		return false;
	}
	float t_2 = (-b - sqrt(discr))/(2*a);
	float t_1 = (-b + sqrt(discr))/(2*a);
 	if (discr>0 ) {
		if (t_2 <0 || t_1 < 0) {
			return false;
		}
	}
	if (t_1<t_2)
	{
		*t = t_1;
	} else
	{
		*t = t_2;
	}
	return true;

}

/*
 * emission characteristic of a spotlight
 */
glm::vec3 SpotLight::getEmission(
		glm::vec3 const& omega // world space direction
		) const
{
	cg_assert(std::fabs(glm::length(omega) - 1.f) < EPSILON);
	float cos = glm::dot(omega,this->direction);


	return glm::vec3(this->power*(this->falloff+2)*glm::pow(std::max(0.f,cos),this->falloff));
}

glm::vec3 evaluate_phong(
	RenderData &data,			// class containing raytracing information
	MaterialSample const& mat,	// the material at position
	glm::vec3 const& P,			// world space position
	glm::vec3 const& N,			// normal at the position (already normalized)
	glm::vec3 const& V)			// view vector (already normalized)
{
	cg_assert(std::fabs(glm::length(N) - 1.f) < EPSILON);
	cg_assert(std::fabs(glm::length(V) - 1.f) < EPSILON);

	glm::vec3 contribution(0.f);

	// iterate over lights and sum up their contribution
	for (auto& light_uptr : data.context.get_active_scene()->lights) 
	{
		const Light *light = light_uptr.get();
		glm::vec3 L = glm::normalize(light->getPosition()-P);
		glm::vec3 Rl =2* glm::dot(L,N) * N - L;

		float visibility = 1.f;
		if (data.context.params.shadows) {
			if( !visible(data,light->getPosition(),P)) {
                 visibility = 0.f;
            }
		}

        float cos_theta = glm::dot(N,L);
		float oOfx = 0.f;
		if (cos_theta > 0.f) {
			oOfx = 1.f;
		}
		float euklDist = glm::pow(glm::length(P - light->getPosition()),2);
		float prev_term = light->getEmission(-L).x * visibility * oOfx/ euklDist;


		glm::vec3 diffuse(0.f);
		if (data.context.params.diffuse) {
			diffuse = mat.k_d * std::max(0.f,cos_theta);
		}

		glm::vec3 specular(0.f);
		if (data.context.params.specular) {
			float cos_varPhi =  glm::dot(Rl,V);
			specular = mat.k_s * std::pow(std::max(0.f,cos_varPhi),mat.n);
		}
		glm::vec3 ambient = data.context.params.ambient ? mat.k_a : glm::vec3(0.0f);

		contribution += diffuse*prev_term;
		contribution += specular*prev_term;

		contribution += ambient * (light->getPower()/euklDist);
	}

	return contribution;
}

glm::vec3 evaluate_reflection(
	RenderData &data,			// class containing raytracing information
	int depth,					// the current recursion depth
	glm::vec3 const& P,			// world space position
	glm::vec3 const& N,			// normal at the position (already normalized)
	glm::vec3 const& V)			// view vector (already normalized)
{
	glm::vec3 R =reflect(V,N);

	const Ray ray = Ray(P + R* data.context.params.ray_epsilon,R);

	return trace_recursive(data,ray,depth+1);
}

glm::vec3 evaluate_transmission(
	RenderData &data,			// class containing raytracing information
	int depth,					// the current recursion depth
	glm::vec3 const& P,			// world space position
	glm::vec3 const& N,			// normal at the position (already normalized)
	glm::vec3 const& V,			// view vector (already normalized)
	float eta)					// the relative refraction index
{
	glm::vec3 contribution(0.f);

	glm::vec3 newDir = glm::vec3(0,0,0);
	if (refract(V,N,eta,&newDir)) {
		const Ray ray = Ray(P + newDir * data.context.params.ray_epsilon,newDir);
		return trace_recursive(data,ray,depth+1);
	}
	return contribution;
}

glm::vec3 handle_transmissive_material_single_ior(
	RenderData &data,			// class containing raytracing information
	int depth,					// the current recursion depth
	glm::vec3 const& P,			// world space position
	glm::vec3 const& N,			// normal at the position (already normalized)
	glm::vec3 const& V,			// view vector (already normalized)
	float eta)					// the relative refraction index
{
	if (data.context.params.fresnel) {

		float F = fresnel(V,N,eta);

		return evaluate_transmission(data, depth, P, N, V, eta)*(1-F) + evaluate_reflection(data, depth, P, N, V)*F;

	} else {
		// just regular transmission
		return evaluate_transmission(data, depth, P, N, V, eta);
	}
}

glm::vec3 handle_transmissive_material(
	RenderData &data,					// class containing raytracing information
	int depth,							// the current recursion depth
	glm::vec3 const& P,					// world space position
	glm::vec3 const& N,					// normal at the position (already normalized)
	glm::vec3 const& V,					// view vector (already normalized)
	glm::vec3 const& eta_of_channel)	// relative refraction index of red, green and blue color channel
{
	if (data.context.params.dispersion && !(eta_of_channel[0] == eta_of_channel[1] && eta_of_channel[0] == eta_of_channel[2])) {

		glm::vec3 contribution(0.f);

		for (int i =0 ; i< 3; i++)
		{
			contribution +=1.f/3.f * handle_transmissive_material_single_ior(data,depth,P,N,V,eta_of_channel[i]);
		}
		return contribution;
	} else {
		// dont handle transmission, take average refraction index instead.
		const float eta = 1.f/3.f*(eta_of_channel[0]+eta_of_channel[1]+eta_of_channel[2]);
		return handle_transmissive_material_single_ior(data, depth, P, N, V, eta);
	}
	return glm::vec3(0.f);
}
// CG_REVISION 8c58412a25ac2367c053bfa840ee81b320bdd315
