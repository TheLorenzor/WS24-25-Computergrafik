#version 330

uniform samplerCube EnvironmentTextureDiffuse;  // environment map, prefiltered for lookup of diffuse lighting
uniform samplerCube EnvironmentTextureGlossy;   // environment map, prefiltered for lookup of specular lighting

uniform vec3 cam_world_pos;    // the camera position in world space
uniform vec3 diffuse_color;    // k_d
uniform vec3 specular_color;   // k_s

in vec3 world_position;              // position of the vertex in world space
in vec3 world_normal_interpolated;   // normal of the vertex in world space

out vec4 frag_color;

void main (void)
{
	// TODO: compute lighting using the prefiltered environment maps
	vec3 v_0 = normalize(world_position-cam_world_pos);
	vec3 r   = normalize(reflect(v_0,world_normal_interpolated));
	vec4 l_d = texture(EnvironmentTextureDiffuse,world_normal_interpolated);
	vec4 l_s = texture(EnvironmentTextureGlossy,r);
	frag_color = vec4(l_d.x*diffuse_color.x,l_d.y*diffuse_color.y,l_d.z*diffuse_color.z,l_d.w);
	frag_color = frag_color + vec4(l_s.x*specular_color.x,l_s.y*specular_color.y,l_s.z*specular_color.z,l_s.w);
}
