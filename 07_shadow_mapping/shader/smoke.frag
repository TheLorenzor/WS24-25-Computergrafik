#version 330

flat in float heat;
flat in float age;
flat in mat2 R_tex;
out vec4 frag_color;
uniform sampler2D fire_texture;
uniform int render_mode;

void main (void)
{
	// Rotate texture coordinates to make sprites rotate.
	vec2 tc = 1.42 * (
		vec2(0.5, 0.5) + R_tex * (gl_PointCoord - vec2(0.5, 0.5))
	);

	// Compute alpha values for smoke and fire from the texture.
	vec4 color = texture(fire_texture, tc);
	float smoke_alpha = color.a * 0.05 * pow(sin(pow(age, 0.3) * 3.124), 2.0);
	float fire_alpha  = color.a * 0.14 * (1.0-age) * smoothstep(0.1, 0.5, age);

	// Compute color values from the texture.
	float luminance  = dot(vec3(0.2, 0.7, 0.1), color.rgb);
	vec3 smoke_color = 0.2 * vec3(luminance);
	vec3 fire_color  = 4.0 * luminance * color.rgb;

	switch (render_mode)
	{
		case 0: // Just smoke.
			frag_color = vec4(smoke_color, smoke_alpha);
			break;

		case 1: // Just fire.
			frag_color = vec4(fire_color * fire_alpha, fire_alpha);
			break;

		case 2: // Both.
		default:
			frag_color = mix(
				vec4(smoke_color * smoke_alpha, smoke_alpha),
				vec4(fire_color  * fire_alpha, 0.0),
				heat
			);
			break;

	}
} 

