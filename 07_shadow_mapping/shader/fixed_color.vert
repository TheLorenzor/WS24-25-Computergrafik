#version 330

layout(location = 0) in vec3 POSITION;

uniform mat4 MVP;
void main(void)
{
	gl_Position = MVP * vec4(POSITION, 1.0);
}
