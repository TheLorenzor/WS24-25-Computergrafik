#version 330

layout(location = 0) in vec3 POSITION;

uniform mat4 VPI;

out vec4 world_pt_near_cam_hom;
out vec4 world_pt_on_nearplane_hom;

void main(void)
{
	gl_Position = vec4(POSITION.xy, 0.999999, 1.0);
	world_pt_near_cam_hom = VPI * vec4(POSITION.xy, -10.0, 1.0);
	world_pt_on_nearplane_hom = VPI * vec4(POSITION.xy, -1.0, 1.0);
}
