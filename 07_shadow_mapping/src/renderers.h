#include <cglib/gl/renderer.h>

void initialize_alpha_blending();
void initialize_additive_blending();
void initialize_premultiplied_blending();

class ShadowmapRenderer : public Renderer
{
protected:
	std::shared_ptr<GLObjModel> donuts;
	std::shared_ptr<GLObjModel> fullscreen_quad;
	std::shared_ptr<FBO> fbo_shadow_map;
	glm::vec3 light_position = glm::vec3(17.1996, 18.5689, -3.53712);
	glm::vec3 light_direction = glm::vec3(-0.741693f, -0.66015f, 0.11872f);
	unsigned int buffer_cube_vert, buffer_cube_idx, vao_cube;
	float light_near = 10.0f;
	float light_far  = 65.0f;
	float light_fov  = 45.0f;
	float bias       = 0.003f;

	const static MapNameProgramDefinition program_definitions;
public:
	ShadowmapRenderer() :
		Renderer(program_definitions)
	{}

	virtual int  _initialize() override;
	virtual const char * get_name() override { return "Shadow Mapping"; }
	virtual void _draw() override;
	virtual void _destroy() override;
};

class ProceduralLandscapeRenderer : public Renderer
{
private:
	std::shared_ptr<GLObjModel> fullscreen_tri;
	std::shared_ptr<GLObjModel> shuttle;
	glm::mat4 shuttle_transform;
	std::shared_ptr<FBO> height_map;
	std::shared_ptr<FBO> wave_map;
	std::shared_ptr<FBO> reflection_map;

	int plane_patch_width, plane_idx_count;
	unsigned int buffer_plane_vert, buffer_plane_idx, vao_plane;

	int landscape_quality_reduction = 3;

	glm::vec2 resolution_divisor;
	glm::ivec2 patches_per_axis;

	std::vector<glm::vec4> fire_particles;
	unsigned int buffer_fire_vert, vao_fire;

	unsigned grassTex, rockTex, sandTex, randTex, waterTex, fireTex;

public:
	glm::vec3 derived_sun_direction;		// calculated from night_time
	glm::vec3 derived_sun_color;			// calculated from night_time
	bool night_time					= false;

	float water_height				= 30.f;
	float wave_scaling				= 1.4f;
	glm::vec3 water_color			= glm::vec3(0.6f, 0.8f, 0.8f);
	float water_absorption			= 0.2f;
	float water_surface_absorption	= 1.5f;
	glm::vec3 derived_water_absorption;		// calculated from water_color, water_surface_absorption
	float reflection_perturbation	= 1.0f;
	glm::vec3 derived_water_surface_color;	// calculated from water_color, water_absorption
	glm::vec3 derived_water_scattering;		// set in _update_parameters

	float height_scaling			= 100.0f;
	float beach_margin				= 2.0f;
	float min_rock_threshold		= 0.5f;
	float rock_blend_margin			= 0.25f;

	float lacunarity				= 2.01f;
	float weight_factor				= 0.5f;
	float terrain_angle_factor		= 121.5f; // degrees
	glm::mat2 derived_rotation_factor;		// calculated from angle_factor
	float valley_smoothness			= 0.65f;
	bool update_heightmap			= true;

	bool visualize_heightmap = false;
	bool visualize_reflection = false;

	bool no_water_surface = false;
	bool no_textures = false;

	bool wireframe = false;

	bool enable_fresnel = true;

	float time = 0.0f;

	const static MapNameProgramDefinition program_definitions;
public:
	ProceduralLandscapeRenderer() :
		Renderer(program_definitions)
	{}

	virtual int  _initialize() override;
	virtual const char * get_name() override { return "Procedural Landscape"; }
	virtual void _draw() override;
	virtual void _update_objects(double time_step) override;
	virtual void _resize(int width, int height) override;
	virtual void _destroy() override;

	void update_parameters();
	void drawLandscape(glm::mat4 const& view, glm::mat4 const& projection, bool under_water);
	void drawParticles(glm::mat4 const& view, glm::mat4 const& projection);
	void drawShuttle(glm::mat4 const& view, glm::mat4 const& projection);
};

class FirePlaceRenderer : public Renderer
{
private:
	struct Particle
	{
		float t;
		float max_t;
		float angle;
		float rot_speed;
		float rand;
		glm::vec3 position;
		glm::vec3 pos_offset;
		float age;
		float scale;
		float heat;
	};

private:
	std::shared_ptr<GLObjModel> fullscreen_tri;
	std::shared_ptr<GLObjModel> logs;
	float particle_px_size = 128;

	std::vector<Particle>  particles[3];
	std::vector<glm::vec3> particle_pos[3];
	unsigned               particle_pos_buffer[3];
	std::vector<glm::vec4> particle_data[3];
	unsigned               particle_data_buffer[3];
	unsigned               particle_vao[3];
	unsigned               fire_tex;

	int num_smoke_particles = 256;
	int num_fire_particles = 1024;
	const static MapNameProgramDefinition program_definitions;
public:
	FirePlaceRenderer() :
		Renderer(program_definitions)
	{}

	virtual int _initialize() override;
	void init_particles();
	void destroy_particles();
	virtual const char * get_name() override { return "FirePlace"; }
	virtual void _draw() override;
	virtual void _update_objects(double time_step) override;
	virtual void _resize(int width, int height) override;
	virtual void _destroy() override;

	void drawLogs(glm::mat4 const& view, glm::mat4 const& projection);
	void drawParticles(glm::mat4 const& view, glm::mat4 const& projection);

};

glm::mat4
compute_view_projection_light(
		const glm::vec3 &light_position,  // position of the light source in world space
		const glm::vec3 &light_direction, // light view direction
		float near_clip_plane,            // distance to the near clip plane
		float far_clip_plane,             // distance to the far clip plane
		float field_of_view,              // field of view in radians
		const glm::vec3 &up_vector);      // vector pointing upwards in world space 

glm::mat4
compute_light_matrix(
	const glm::mat4 &view_projection_light, // View Projection Matrix for the light
	float shadow_bias);                     // shadow bias to prevent shadow acne

// CG_REVISION d4ab32bd208749f2d2b1439e25d16e642b039298
