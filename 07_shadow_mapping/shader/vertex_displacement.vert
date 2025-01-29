#version 330

layout(location = 0) in vec2 POSITION; // position of the current vertex in the 2-dimensional patch grid

uniform ivec2 instanced_patch_mask_and_shift; // used to compute the global vertex coordinate from the patch instance id
uniform float instanced_patch_size;           // local size of one patch instance, used to compute the global vertex coordinate from the patch instance id

uniform mat4 MVP;             // terrain space to clip space
uniform mat4 model_matrix;    // terrain space to world space
uniform mat3 normal_matrix;   // terrain space to world space

uniform sampler2D HeightMap;           // height map with unscaled height values, stored in the red component
uniform vec2 one_over_height_map_size; // 1 / textureSize(HeightMap)
uniform float height_scaling;          // height scaling, guaranteed to be identical in terrain space and world space

out vec3 world_position;             // world-space position
out vec3 world_normal_interpolated;  // world-space normal
out vec2 grid_coord;                 // passed-through height field coordinate for fragment shader, do not modify

void main(void)
{
	// height field coordinate in grid units, where one height map pixel corresponds to one square unit
	grid_coord = POSITION + instanced_patch_size * vec2(gl_InstanceID & instanced_patch_mask_and_shift.x, gl_InstanceID >> instanced_patch_mask_and_shift.y);

	// position in local terrain space
	vec3 local_position = vec3(grid_coord.x, 0.0, grid_coord.y);

	// TODO: compute vertically displaced position in local terrain space

	// normal in local terrain space
	vec3 local_normal = vec3(0.0, 1.0, 0.0);

	// TODO: compute the actual normal of the displaced surface in local terrain space

	// TODO: transform from terrain space to required spaces
	gl_Position = vec4(100.0); // do sth. sensible
	world_position = vec3(100.0); // do sth. sensible
	world_normal_interpolated = vec3(100.0); // do sth. sensible here
}
