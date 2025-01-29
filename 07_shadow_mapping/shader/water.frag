#version 330

uniform sampler2D ReflectionMap;

uniform mat4 RVP_to_tex;                 // transforms points and directions from world space to homogeneous reflection map coordinates
uniform vec3 cam_world_pos;              // camera position in world coordinates
uniform vec3 surface_transmission_color; // the amount of light transmitted through the water surface
uniform vec3 water_transmission_coeff;   // for OPTIONAL task: absorption coefficients to compute the amount of light transmitted through the water volume. Transmission = exp(-coeff * length)
uniform float water_height;              // world-space height of the water plane
uniform float reflection_perturbation;
uniform int enable_fresnel;

in vec3 world_position;             // world-space position of the current water surface fragment
in vec3 world_normal_interpolated;  // unnormalized world-space normal of the current water surface fragment

layout(location = 0, index = 0) out vec4 add_color; // additive output color: rgb are added to the framebuffer (after multiplication by mul_color), alpha is ignored
layout(location = 0, index = 1) out vec4 mul_color; // multiplicative output color: the original framebuffer color is multiplied with rgb, alpha is ignored

vec3 restrict_normal_to_planar_reflection_info(vec3 actual_normal, vec3 cam_dir)
{
	// cam dir projected onto ground plane
	vec3 flat_cam_dir = vec3(cam_dir.x, 0.0, cam_dir.z);

	// fix up reflection normal to never reflect below the water plane (clamp to never point away from camera)
	vec3 reflection_normal = actual_normal;
	reflection_normal = reflection_normal - flat_cam_dir 
		              * max(dot(reflection_normal, flat_cam_dir), 0.0) 
					  / dot(flat_cam_dir, flat_cam_dir);
	reflection_normal = normalize(reflection_normal);
	return reflection_normal;
}

vec3 compute_attenuated_reflection_dir(vec3 world_normal, vec3 cam_pos, vec3 cam_dir)
{
	// restrict the reflection normal to reflect to valid parts of the reflection map
	vec3 reflection_normal = restrict_normal_to_planar_reflection_info(world_normal, cam_dir);

	// compute perturbed & unperturbed reflection directions and blend between them
	// (reflection map does not contain all necessary information for
	// arbitrary diretions, therefore reduce perturbation with distance)
	vec3 perturbed_reflect_dir = reflect(cam_dir, reflection_normal);
	vec3 planar_reflect_dir = reflect(cam_dir, vec3(0.0, 1.0, 0.0));
	float attenuated_reflection_perturbation = 
		mix(reflection_perturbation, 0.05 * reflection_perturbation, 
			min(length(world_position - cam_world_pos) / 500.0, 1.0));
	vec3 reflect_dir = normalize(
		mix(planar_reflect_dir, perturbed_reflect_dir, 
			attenuated_reflection_perturbation) );
	return reflect_dir;
}

void main (void)
{
	vec3 world_normal = normalize(world_normal_interpolated);
	vec3 cam_dir = normalize(world_position - cam_world_pos);
	vec3 reflection_lookup_dir = 
		compute_attenuated_reflection_dir(world_normal, cam_world_pos, cam_dir);

	// TODO: transform reflection_lookup_dir to homogeneous reflection map space,
	// perform dehomogenization,
	// and look up the reflected light in the reflection map
	vec3 reflected_light = vec3(0.5);

	float reflectivity = 0.5;
	if(enable_fresnel > 0) {
		// TODO: compute reflectivity using Schlick's approximation
		// reflectivity = ...
	}
	
	// TODO: compute additive light and multiplicative blend factor
	// using the computed reflectivity, the looked-up reflected light
	// and the water surface transmission color
	// NOTE: in this OpenGL blending setup, the w component of the
	// colors is irrelevant. However, you can output your fresnel
	// reflectivity coefficient in add_color.w, in which case we
	// may be able to give you more points when some part of the
	// calculation is wrong.
	add_color = vec4(0.0); // do sth. sensible here
	mul_color = vec4(smoothstep(0.97, 1.1, world_normal.y)); // do sth. sensible here
	
	// OPTIONAL: handle camera under water
	if (cam_world_pos.y < water_height)
	{
		// OPTIONAL: this will not be graded, but you may insert code to handle
		// cases where the camera is submerged below the water surface here.
		// For this, you need to compute the distance from the camera
		// to the water surface, compute the corresponding transmission
		// and then apply it to both the additive and multiplicative
		// colors to account for absorption of all light travelling from
		// the surface down to the camera.
	}
}
