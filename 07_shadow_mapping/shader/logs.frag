#version 330

in vec3 world_position;
in vec3 world_normal_interpolated;

out vec4 frag_color;

void main (void)
{
	vec3 light_color = vec3(0.8, 0.8, 1.0) * 0.1;
	vec3 light_dir = vec3(0.0, -1.0, 0.0);

	vec3 world_normal = normalize(world_normal_interpolated);

	float nDotL = dot(-light_dir, world_normal);

	vec3 color = max(0.1 + 0.9 * nDotL, 0.0) * vec3(82, 13, 2) / 255.0 * mix(.1, 1., smoothstep(.0, 6.0, world_position.y));
	frag_color = vec4(color * light_color, 1.0);
} 

