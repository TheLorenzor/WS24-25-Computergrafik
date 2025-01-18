#version 330

in vec3 color;     // interpolated color passed from vertex shader
in vec2 tex_coord; // interpolated uv coordinates passed from vertex shader

out vec4 frag_color; // the final color of the fragment

uniform sampler2D tex; // a simple 2D texture

void main (void)
{
	frag_color = vec4(texture(tex, tex_coord).rgb, 1.0);
}
