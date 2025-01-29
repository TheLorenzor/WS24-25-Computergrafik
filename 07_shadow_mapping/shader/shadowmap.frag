#version 330

uniform sampler2D shadow_map; // contains the depth texture rendered from the
                              // lights point of view

uniform vec3 world_light_position;       // position of the light source in world space
uniform vec3 world_light_view_direction; // view direction of the spot light in world space

uniform vec3 K_d;                  // diffuse color of the object

uniform float spotlight_cos_inner; // precomputed constants for the spotlight falloff
uniform float spotlight_cos_outer;

uniform float shadow_bias; // shadow bias to prevent shadow acne

uniform vec3 light_emission; // power of the light source

in vec3 world_normal_interpolated; // interpolated normal in world space
in vec4 pos_shadowmap_space;       // position of the fragment in shadow map space
in vec3 world_position;            // position of the fragment in world space

out vec4 frag_color;

float
spot_light(
		vec3 L,
		vec3 spot_light_view_dir)
{
	float a = dot(L, -spot_light_view_dir);
	float x = (a - spotlight_cos_outer)
	        / (spotlight_cos_inner - spotlight_cos_outer);
	return clamp(x, 0.0, 1.0);
}

vec3 
shade(vec3 world_normal, 
	  vec3 world_position, 
	  vec3 world_light_position)
{
	vec3 light_vec  = world_light_position - world_position;
	vec3 diffuse = K_d * max(dot(world_normal, normalize(light_vec)), 0.0);
	vec3 light = light_emission / dot(light_vec, light_vec)
		* spot_light(normalize(light_vec), world_light_view_direction);
	return diffuse * light;
}

void main (void)
{
	// TODO: project position in light space

	// TODO: determine visibility for the fragment
	float visibility = 1.0;

	vec3 color = vec3(0);
	if (visibility > 0.0) {
		color = shade(normalize(world_normal_interpolated),
			  world_position,	
			  world_light_position);
	}
	frag_color = vec4(color, 1.0);
} 
