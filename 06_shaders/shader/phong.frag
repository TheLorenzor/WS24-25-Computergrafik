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
	frag_color = vec4(0.0,0.0,0.0,1.0);
	// TODO: iterate over all lights and accumulate illumination
	// according the the phong illumination model on the exercise sheet
	for (int i = 0; i < light_count; ++i)
	{
		vec3 L = normalize(light_world_pos[i]-world_position);
		vec3 N = normalize(world_normal_interpolated);
		//kd calcualtion
		vec4 k_ds = k_d * max(0,dot(L,N));


		vec3 m_V =world_position-cam_world_pos;
		vec3 R = reflect(m_V,world_normal_interpolated);
		vec3 kss = k_s * pow(max(0,dot(R,L)),2);
		vec3 x = light_world_pos[i]-world_position;
		float i_x = light_intensity[i].x / pow(x.x,2);
		float i_y = light_intensity[i].y / pow(x.y,2);
		float i_z = light_intensity[i].z / pow(x.z,2);
		vec3 I = vec3(i_x,i_y,i_z);
		frag_color =frag_color+ (k_ds+vec4(kss,0.0))*vec4(I,0.0);
	}

}
