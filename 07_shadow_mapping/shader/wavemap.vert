#version 330

layout(location = 0) in vec3 POSITION;

out vec2 signed_normalized_coord;

void main(void)
{
	gl_Position = vec4(POSITION.xy, 0.0, 1.0);
	signed_normalized_coord = POSITION.xy;
}
