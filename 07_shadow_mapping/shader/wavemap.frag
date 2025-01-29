#version 330

uniform sampler2D WaveTexture;
uniform float time;

in vec2 signed_normalized_coord;
out vec4 wavemap;

void main(void)
{
	vec2 ncoord = gl_FragCoord.xy * 0.0078125;
	vec2 ndispl = 0.021 * vec2(0.1, 0.9) * time;
	float height = -abs(texture(WaveTexture, ncoord + ndispl).x - texture(WaveTexture, ncoord.yx - ndispl).x );
	wavemap = vec4(height);
} 
