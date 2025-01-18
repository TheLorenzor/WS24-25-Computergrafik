#version 330

uniform samplerCube EnvironmentTexture;

in vec3 viewdir_interpolated;
out vec4 frag_color;

void main (void)
{
	vec3 light = texture(EnvironmentTexture, viewdir_interpolated).rgb;
	frag_color = vec4(light, 1.);
} 
