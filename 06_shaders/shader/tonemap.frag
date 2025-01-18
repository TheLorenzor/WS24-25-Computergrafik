#version 330

uniform sampler2D HDRTexture;

in vec2 texcoord;
out vec4 frag_color;

void main (void)
{
	vec4 hdrColor = texture(HDRTexture, texcoord);
	frag_color.rgb = hdrColor.rgb / (dot(hdrColor.rgb, vec3(0.33333)) + 1.0);
	frag_color.a = hdrColor.a;
} 
