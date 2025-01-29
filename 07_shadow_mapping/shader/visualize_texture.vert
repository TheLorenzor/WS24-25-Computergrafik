#version 330

layout(location = 0) in vec3 POSITION;

uniform mat4 VPI;

out vec2 tex_coord;

void main(void)
{
	gl_Position = vec4(POSITION.xy, 0.0, 1.0);
	tex_coord = POSITION.xy * 0.5 + vec2(0.5);
}
