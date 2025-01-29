#version 330

layout(location = 0) in vec3 POSITION;
layout(location = 1) in vec3 NORMAL;
layout(location = 2) in vec2 TEXCOORD;

uniform mat4 MVP;
uniform mat4 model_matrix;
uniform mat3 normal_matrix;

out vec3 world_position;
out vec3 world_normal_interpolated;
out vec2 texcoord;

void main(void)
{
	gl_Position = MVP * vec4(POSITION, 1.0);
	world_position = vec3(model_matrix * vec4(POSITION, 1.0));
	world_normal_interpolated = normal_matrix * NORMAL;
	texcoord = TEXCOORD;
}
