#version 330

uniform sampler2D GrassTex;  // grass color texture
uniform sampler2D RockTex;   // rock color texture
uniform sampler2D SandTex;   // sand color texture
uniform sampler2D HeightMap; // height map

uniform vec3 cam_world_pos;            // camera position in world coordinates
uniform vec3 sun_world_dir;            // sun light direction in world coordinates
uniform vec3 sun_color;                // sun color
uniform vec2 one_over_height_map_size; // 1 / textureSize(HeightMap)
uniform float height_scaling;          // world-space terrain height scaling
uniform float water_height;            // water height in world coordinates
uniform float beach_margin;            // world-space height margin above the water inside which there is only sand
uniform float min_rock_threshold;      // minimum threshold where rock starts to fade in
uniform float rock_blend_margin;       // margin inside which rock is gradually blended in
uniform vec3 water_transmission_coeff; // volumetric transmission coefficients for rgb light components. Transmission = exp(-coeff * length)
uniform bool under_water;              // for OPTIONAL task: specifies if camera under water (may not be deduced from cam position alone!)

uniform bool use_textures;             // enables texturing. Already handled by the existing code

in vec3 world_position;            // world-space position of the current terrain fragment
in vec3 world_normal_interpolated; // unnormalized world-space terrain normal
in vec2 grid_coord;                // height field coordinate of the current terrain fragment, in grid units

out vec4 frag_color; // output color for the current fragment

// samples the given color texture with two scaling factors to add larger-scale variation to the highly repetitive textures
vec3 color_varied(sampler2D tex, vec2 coord)
{
	// sample upscaled texture at a lower resolution to get some color variation
	vec3 colorVariation = texture(tex, 0.07 * coord, 2.0).rgb;
	// somewhat normalize with respect to brightness
	colorVariation *= 3.0 / sqrt(dot(colorVariation, vec3(0.333)) + 0.0001);
	// multiply by detail texture
	return colorVariation * texture(tex, coord).rgb;
}

void main (void)
{
	vec3 n = normalize(world_normal_interpolated);

	// sample color textures
	vec3 grass_color = use_textures 
		?  vec3(0.9, 0.65, 0.9) * color_varied(GrassTex, 0.2 * world_position.xz).rgb 
		: vec3(0.65f);
	vec3 sand_color = use_textures 
		? 0.35 * color_varied(SandTex, 0.25 * world_position.xz).rgb 
		: vec3(0.65f);
	vec3 rock_color = use_textures 
		? color_varied(RockTex, 0.05 * world_position.xz).rgb 
		: vec3(0.65f);


	float sinval = sqrt(1-pow(dot(n,vec3(0,1,0)),2));
	float value = smoothstep(min_rock_threshold,min_rock_threshold+rock_blend_margin,sinval);
	vec3 color = rock_color;
	color = rock_color * value+ grass_color* (1-value);
	float blend_alpha =world_position.y-(water_height+beach_margin);
	if (blend_alpha <= 0.0 ) {
		color = sand_color;
	}
	if (blend_alpha <= 1.0 && blend_alpha >0.0) {
		color = sand_color*(1-blend_alpha)+ blend_alpha*color;
	}

	// fake ambient occlusion
	vec2 heightmap_tex_coord = (grid_coord + vec2(.5)) 
		                     * one_over_height_map_size;
	float height = height_scaling 
		         * textureLod(HeightMap, heightmap_tex_coord, 0).x;
	float filtered_height1 = height_scaling 
		                   * textureLod(HeightMap, heightmap_tex_coord, 4).x;
	float filtered_height2 = height_scaling 
		                   * textureLod(HeightMap, heightmap_tex_coord, 8).x;
	// difference between actual height & mip-filtered heights gives
	// some sense of occlusion by surrounding geometry
	float ambient_occlusion = 1.5 / (1.0 + 0.1 * max(filtered_height1 - height, 0.0));
	ambient_occlusion *= 1.0 / (1.0 + 0.015 * max(filtered_height2 - height, 0.0));

	// basic lambert shading w/ ambient occlusion
	color *= ambient_occlusion * sun_color * max(dot(n, -sun_world_dir), 0.02);

	// TODO: compute light absorption by the volume below the water surface
	// for the underwater segments of both the camera and the sun light rays


	float height_diff =  water_height-world_position.y;
	float height_cam_diff = normalize(cam_world_pos-world_position).y- world_position.y;
	float height_light_diff = normalize(sun_world_dir-world_position).y;
	float underwater_camray_length = abs(height_diff/height_cam_diff); // compute sth. sensible
	float underwater_sun_lightray_length = abs(height_diff/height_light_diff); // compute sth. sensible
	if (height_diff>0.0) {
		color = exp(-water_transmission_coeff * (underwater_camray_length+underwater_sun_lightray_length))*color;
	}

	// OPTIONAL: handle camera under water
	if (under_water) {
		// OPTIONAL: this will not be graded, but you may insert code to handle
		// cases where the camera is submerged below the water surface here.
		// In this case, you need to fix computation of the length of the
		// underwater camera ray segment.
	}

	frag_color = vec4(color, 1.0);
}
