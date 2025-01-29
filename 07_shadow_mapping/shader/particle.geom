#version 330

uniform mat4 MVP; // model view projection matrix that transforms points from world space to clip space

uniform vec3 cam_world_pos; // position of the camera in world coordinates

layout(points) in; // particle geometry shader retrieves single points as inputs
layout(triangle_strip, max_vertices = 4) out; // particle geometry shader outputs 4 triangle strip vertices to form a camera-aligned quad

in vec3 world_position[1]; // position of the particle in world coordinates
in float world_radius[1]; // radius of the particle in world coordinates

out vec3 world_pos; // position of each particle vertex in world coordinates
out vec2 tex_coord; // texture coordinate of each particle vertex

void
main()
{
	float randomAngle = fract(world_radius[0]) * 6.28;

	// TOOD: compute camera-aligned particle coordinate frame
	
	// TOOD: compute 4 triangle strip vertices to form a camera-aligned particle quad
	world_pos = vec3(0.0); // fill with sensible data for each vertex
	gl_Position = vec4(0.0); // fill with sensible data for each vertex
	tex_coord = vec2(0.0); // fill with sensible data for each vertex
}
