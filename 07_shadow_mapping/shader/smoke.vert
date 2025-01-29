#version 330

layout(location = 0) in vec3 WORLD_POSITION;
layout(location = 1) in vec4 HEAT_AGE_RAND_SIZE;  
flat out float heat;
flat out float age;
flat out mat2 R_tex;

uniform mat4 MVP;
uniform float particle_px_size;

void main(void)
{
	heat = HEAT_AGE_RAND_SIZE.x;
	age = HEAT_AGE_RAND_SIZE.y;
	float rand = HEAT_AGE_RAND_SIZE.z;

	float theta = 2.0 * (rand-0.5) * 0.4 * (age + rand) * 2.0 * 3.124;
	float st = sin(theta);
	float ct = cos(theta);
	R_tex = mat2(ct, st, -st, ct);

	gl_Position   = MVP * vec4(WORLD_POSITION, 1.0);
	gl_PointSize  = particle_px_size * HEAT_AGE_RAND_SIZE.w / gl_Position.w;
}

