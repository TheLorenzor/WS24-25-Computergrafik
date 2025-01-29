#version 330

uniform sampler2D tex; // the texture containing the diffuse material parameter k_d

#define MAX_LIGHT_COUNT 16
uniform int light_count;                        // the number of lights (always smaller than MAX_LIGHT_COUNT
uniform vec3 light_world_pos[MAX_LIGHT_COUNT];  // the position of the lights in world space
uniform vec3 light_intensity[MAX_LIGHT_COUNT];  // the intensity of the lights

uniform vec3 k_s; // the specular material parameter
uniform float n;  // phong exponent of the material

uniform vec3 cam_world_pos; // the camera position in world space

in vec3 world_position;            // the (interpolated) world space position corresponding to the fragment
in vec3 world_normal_interpolated; // the (interpolated) world space normal
in vec2 tex_coord;                 // the (interpolated) uv coordinate

out vec4 frag_color; // the resulting color value (will be written into the framebuffer)

void
main()
{
	// TODO: read k_d from texture

	vec4 k_d = texture(tex,tex_coord);
	// TODO: iterate over all lights and accumulate illumination
	// according the the phong illumination model on the exercise sheet
	vec3 N = normalize(world_normal_interpolated);

	for (int i = 0; i < light_count; ++i)
	{
		// Vector calculation
		vec3 L = normalize(light_world_pos[i]-world_position);
		vec3 V = normalize(cam_world_pos-world_position);
		vec3 R = reflect(V,N);
		// I_l calculation
		vec3 x = light_world_pos[i]-world_position;
		float i_x = light_intensity[i].x / pow(x.x,2);
		float i_y = light_intensity[i].y / pow(x.y,2);
		float i_z = light_intensity[i].z / pow(x.z,2);

		vec4 t = (k_d*max(0,dot(L,N)) + vec4(k_s,1.0) * pow(max(0,dot(R,L)),n));

		frag_color = frag_color + vec4(t.x*i_x,t.y*i_y,t.z*i_z,1.0);
	}

}
