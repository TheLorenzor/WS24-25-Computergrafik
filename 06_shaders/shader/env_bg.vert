#version 330

layout(location = 0) in vec3 POSITION;

uniform mat4 VPI;

out vec3 viewdir_interpolated;

void main(void)
{
	gl_Position = vec4(POSITION.xy, 0.0, 1.0);
	vec4 world_pt_near_cam_hom = VPI * vec4(POSITION.xy, -1.1, 1.0);
	vec4 world_pt_on_nearplane_hom = VPI * vec4(POSITION.xy, -1.0, 1.0);
	viewdir_interpolated = world_pt_on_nearplane_hom.xyz / world_pt_on_nearplane_hom.w - world_pt_near_cam_hom.xyz / world_pt_near_cam_hom.w;
}
