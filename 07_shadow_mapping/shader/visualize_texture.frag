#version 330

uniform sampler2D Texture;
uniform bool monochrome;

in vec2 tex_coord;
out vec4 frag_color;

void main (void)
{
	vec3 color = texture(Texture, tex_coord).rgb;
	frag_color = vec4(monochrome ? color.rrr : color, 1.);
} 
