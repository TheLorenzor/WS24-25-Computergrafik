#version 330

uniform vec3 cam_world_pos;
uniform vec3 sun_world_dir;
uniform vec3 sun_color;

in vec3 world_position;
in vec3 world_normal_interpolated;
in vec2 texcoord;

out vec4 frag_color;

void main (void)
{
	vec3 world_normal = normalize(world_normal_interpolated);
	vec3 from_cam = normalize(world_position - cam_world_pos);
	vec3 reflection_dir = reflect(from_cam, world_normal);

	// phong ...
	float nDotL = dot(-sun_world_dir, world_normal);
	float powRDotL = 0.0;
	if (nDotL > 0.0)
		powRDotL = pow(max(dot(-sun_world_dir, reflection_dir), 0.0), 48 );

	vec3 color = max(0.25 + 0.75 * nDotL, 0.0) * vec3(1.0, 0.65, 0.5) + vec3(6.0 * powRDotL);
	frag_color = vec4(color * sun_color, 1.0);
} 
