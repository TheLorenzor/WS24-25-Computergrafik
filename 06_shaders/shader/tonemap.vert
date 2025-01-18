#version 330

layout(location = 0) in vec3 POSITION;

out vec2 texcoord;

void main(void)
{
	gl_Position = vec4(POSITION.xy, 0.0, 1.0);
	texcoord = POSITION.xy * 0.5 + vec2(0.5);
}
