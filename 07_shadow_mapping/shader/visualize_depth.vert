#version 330

layout(location = 0) in vec3 POSITION;

out vec2 tex_coord;

uniform mat4 MVP;

void main(void)
{
	gl_Position = MVP * vec4(POSITION, 1.0);
	tex_coord = POSITION.xy * 0.5 + 0.5;
}
