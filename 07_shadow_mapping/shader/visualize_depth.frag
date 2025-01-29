#version 330

//uniform sampler2DShadow shadow_map;
uniform sampler2D shadow_map;

in vec2 tex_coord;

out vec4 frag_color;

uniform float near;
uniform float far;

void main (void)
{
	float depth = texture(shadow_map, tex_coord).r;

	frag_color = vec4((2.0 * near) / (far + near - depth * (far - near)));
} 
