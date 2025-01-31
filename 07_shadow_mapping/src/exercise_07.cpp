#include <cglib/gl/renderer.h>

#include <glm/gtx/transform.hpp>
#include <glm/gtc/matrix_inverse.hpp>

glm::mat4
compute_view_projection_light(
		const glm::vec3 &light_position,  // position of the light source in world space
		const glm::vec3 &light_direction, // light view direction
		float near_clip_plane,            // distance to the near clip plane
		float far_clip_plane,             // distance to the far clip plane
		float field_of_view,              // field of view in radians
		const glm::vec3 &up_vector)       // vector pointing upwards in world space
{


	glm::mat4 view_matrix = glm::lookAt(light_position, light_position + light_direction, up_vector);
	glm::mat4 projection_matrix = glm::perspective(field_of_view, 1.0f, near_clip_plane, far_clip_plane);

	return projection_matrix * view_matrix;
}

void
initialize_alpha_blending()
{
	// TODO: Set up the blend equation and blend function for 
	// alpha blending.
	glBlendEquation(GL_FUNC_ADD);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
}

void
initialize_additive_blending()
{
	// TODO: Set up the blend equation and blend function for 
	// additive blending.
	glBlendFunc(GL_ONE_MINUS_SRC_ALPHA, GL_ONE);
	glBlendEquation(GL_FUNC_ADD);

}

void
initialize_premultiplied_blending()
{
	// TODO: Set up the blend equation and blend function for 
	// blending with premultiplied alpha.
	glBlendFunc(GL_ONE, GL_ONE_MINUS_SRC_ALPHA);
	glBlendEquation(GL_FUNC_ADD);
}

// CG_REVISION d4ab32bd208749f2d2b1439e25d16e642b039298
