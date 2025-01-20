#version 330

layout(location = 0) in vec3 POSITION; // position of the vertex in object space
layout(location = 1) in vec3 NORMAL;   // normal of the vertex in object space
layout(location = 2) in vec2 TEXCOORD; // uv coordinates of the vertex

uniform mat4 MVP;              // the model-view-projection matrix
uniform mat4 model_matrix;     // the world transformation of the model
uniform mat3 normal_matrix;    // the world transformation for normals

out vec3 world_position;              // position of the vertex in world space
out vec3 world_normal_interpolated;   // normal of the vertex in world space
out vec2 tex_coord;                   // uv coordinates of the vertex

void
main()
{
	// TODO: correctly set gl_Position, world_position, world_normal_interpolated, tex_coord
	gl_Position = MVP * vec4(POSITION, 1.0);
	vec4 local_pos = vec4(POSITION,1);
	local_pos = model_matrix * local_pos;
	world_normal_interpolated = vec3(1.0);

	world_normal_interpolated = normal_matrix * NORMAL;

	world_position = vec3(local_pos.x/local_pos.w,local_pos.y/local_pos.w,local_pos.z / local_pos.w);
	tex_coord = TEXCOORD;
}
