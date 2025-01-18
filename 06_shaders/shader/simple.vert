#version 330

layout(location = 0) in vec3 POSITION; // position of the vertex in object space
layout(location = 2) in vec2 TEXCOORD; // uv coordinates of the vertex

uniform mat4 MVP; // model-view-projection matrix

out vec3 color;
out vec2 tex_coord;

void main(void)
{
	tex_coord = TEXCOORD;
	// TODO: correctly set gl_Position and color
	gl_Position = vec4(1.0);
	color = vec3(1.0);
}
