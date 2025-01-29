#version 330

uniform vec3 sun_world_dir;
uniform vec3 sun_color;
uniform bool night_time;
uniform uint seed;

in vec4 world_pt_near_cam_hom;
in vec4 world_pt_on_nearplane_hom;
//in vec3 viewdir_interpolated;
out vec4 frag_color;

vec4 noiseFromPixelPosition(uvec2 position, uvec2 resolution, uint seed)
{
	uint idx = position.x + position.y * resolution.x;
	uint area = position.x * position.y;

	uvec4 idcs = seed + idx * uvec4(41, 29, 53, 43);
	uvec4 areas = seed + area  * uvec4(23, 59, 47, 37);

	return vec4((idcs ^ areas) % 661u) / 330.5f - 1.0f;
}

void main (void)
{
	// dehomogenization needs to be performed _here_, no longer linear w/ oblique near plane clipping matrix (-> reflection matrix)
	vec3 viewdir_interpolated = world_pt_on_nearplane_hom.xyz / world_pt_on_nearplane_hom.w - world_pt_near_cam_hom.xyz / world_pt_near_cam_hom.w;
	vec3 viewdir = normalize(viewdir_interpolated);

	// 'plausible' horizon blend factors
	float horizonBlend = min(1.0 - viewdir.y, 1.0);
	horizonBlend *= horizonBlend;
	horizonBlend = 1.0 - horizonBlend;
	
	vec3 light;
	
	// blend between horizon and zenith
	if (night_time)
		light = mix( vec3(221, 164, 121) * 0.15, vec3(1, 9, 71), horizonBlend ) / 4000.0;
	else
		light = mix( vec3(121, 164, 221), vec3(1, 9, 71), horizonBlend ) / 255.0;
	
	// trace ray against sun / moon
	{
		// that's why they call it physically-based: real-world scaling w/ artistic exaggeration :-)
		vec3 center = -sun_world_dir * 380000.0; // moon distance / 1000
		float radius = night_time ? 5.0 * 1800.0 : 4.0 * 3600.0; // real-world moon & sun radii scaled to our artificial distance, w/ artistic scaling

		// intersect view ray w/ sphere
		float c = dot(center, center) - radius * radius;
		float b = dot(viewdir, -center);
		float d = b * b - c;

	    if (d >= 0.0)
		{
			float d = -b - sqrt(d);
			vec3 p = viewdir * d;
			vec3 n = normalize(p - center);

			if (night_time)
				// lambertian additive moon
				light += vec3( max(0.15 + 0.85 * dot(n, normalize(vec3(2.0, -1.5, -1.6))), 0.0015) ) * sun_color / dot(vec3(.9), sun_color);
			else
				// glowing additive sun
				light += vec3( pow(max(dot(n, -viewdir), 0.0), 8) ) * sun_color;
		}
    }

	// dithering
	if (night_time)
		light += noiseFromPixelPosition(uvec2(gl_FragCoord.xy), uvec2(2048), seed).xyz * 0.0004;

	frag_color = vec4(light, 1.);
} 
