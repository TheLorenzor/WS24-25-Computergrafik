#version 330

layout(location = 0) in vec4 POSITION_RADIUS; // xyz = world-space particle position, w = world-space particle radius

out vec3 world_position; // world-space particle position
out float world_radius; // world-space particle radius

void
main()
{
	// simply pass vertex attributes through to geometry shader
	world_position = POSITION_RADIUS.xyz;
	world_radius = POSITION_RADIUS.w;
}
