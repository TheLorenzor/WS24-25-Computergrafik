#version 330

layout(location = 0) in vec3 POSITION;
layout(location = 1) in vec3 NORMAL;
layout(location = 2) in vec2 TEXCOORD;

uniform mat4 MVP;
uniform mat4 light_matrix;
uniform mat4 model_matrix;
uniform mat3 normal_matrix;

out vec3 color;
out vec3 world_normal_interpolated;
out vec4 pos_shadowmap_space;
out vec3 light_vec;
out vec3 world_position;

void main(void)
{
	gl_Position = MVP * vec4(POSITION, 1.0);
	vec4 pos_ws = model_matrix * vec4(POSITION, 1.0);
	pos_shadowmap_space = light_matrix * pos_ws;
	world_position = pos_ws.xyz;
	world_normal_interpolated = normal_matrix * NORMAL;
}
