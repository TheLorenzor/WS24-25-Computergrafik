#version 330

uniform sampler2D Texture;

in vec3 world_pos;
in vec2 tex_coord;
out vec4 frag_color;

void
main()
{
	vec4 color = texture(Texture, tex_coord);
	// premultiplied alpha, add
	color.rgb *= color.a;
	// reduce intensity
	color.rgb *= 0.1;
	// cooling (very hacky particle coloring from here on)
	float cooling = 1.0 - clamp(world_pos.y / 700.0, 0.0, 1.0);
	color.a *= cooling * cooling;
	color.rgb = mix(color.rgb, vec3(0.05 * color.a), cooling);
	color.rgb *= 1.0 - cooling;
	frag_color = color;
}
