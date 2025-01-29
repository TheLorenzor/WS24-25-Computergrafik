#version 330

uniform sampler2D RandomTexture;

uniform float lacunarity;
uniform float weight_factor;
uniform mat2 rotation_factor;
uniform float valley_smoothness;
uniform vec2 resolution_divisor;

in vec2 signed_normalized_coord;
out vec4 heightmap;

float unsignedToAbs(float us)
{
	return abs(2.0 * us - 1.0);
}

float noise(vec2 coord)
{
	return textureLod(RandomTexture, coord, 0).x;
}

void main(void)
{
	const mat2 rm = mat2(0.80,  0.60, -0.60,  0.80);

	vec2 planePos = gl_FragCoord.xy;
	vec2 scaledPos = planePos / 2048.0 / 16.0 * resolution_divisor;

	float weight = 1.0;

//	scaledPos *= rotation_factor;
	float envelope_height = weight * noise(scaledPos);
	scaledPos *= rotation_factor;
	envelope_height += weight * noise(scaledPos);
	envelope_height *= 0.5;
	scaledPos *= rotation_factor * lacunarity; weight *= weight_factor;
	
	envelope_height += weight * noise(scaledPos); scaledPos *= rotation_factor * lacunarity; weight *= weight_factor;
	envelope_height = 0.5 * pow(envelope_height, 3);

	float detail_height = weight * noise(scaledPos); scaledPos *= rotation_factor * lacunarity; weight *= weight_factor;
	detail_height += weight * noise(scaledPos); scaledPos *= rotation_factor * lacunarity; weight *= weight_factor;
	detail_height += weight * noise(scaledPos); scaledPos *= rotation_factor * lacunarity; weight *= weight_factor;
	detail_height += weight * noise(scaledPos); scaledPos *= rotation_factor * lacunarity; weight *= weight_factor;
	detail_height += weight * noise(scaledPos); scaledPos *= rotation_factor * lacunarity; weight *= weight_factor;
	detail_height += weight * noise(scaledPos); scaledPos *= rotation_factor * lacunarity; weight *= weight_factor;

	float height = envelope_height + mix(1.0, envelope_height, valley_smoothness) * detail_height;
//	height /= mix( 1.0, 10.0, pow(dot(signed_normalized_coord, signed_normalized_coord), 16) );

	heightmap = vec4(height);
} 
