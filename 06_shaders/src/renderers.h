#pragma once

#include <cglib/gl/renderer.h>

class CubeRenderer : public Renderer
{
public:
	CubeRenderer()
		: Renderer(program_definitions)
	{
	}

	const char * get_name() override { return "Cube"; }
protected:
	void _draw() override;
	void _destroy() override;
	void _activate() override;
	void _deactivate() override;
	int  _initialize() override;

	unsigned int buffer_cube_vert, buffer_cube_idx, vao_cube;

	static const MapNameProgramDefinition program_definitions;
};

class MonkeyRenderer : public Renderer
{
public:
	MonkeyRenderer()
		: Renderer(program_definitions)
	{
	}
	void _draw() override;
	void _activate() override;
	void _deactivate() override;
	int  _initialize() override;
	void _destroy() override;

	const char * get_name() override { return "Monkeys"; }
protected:
	std::shared_ptr<GLObjModel> monkey, fullscreen_tri;
	unsigned envTex, envTexGlossy, envTexDiffuse;
	const static MapNameProgramDefinition program_definitions;
};

class HDRRenderer
{
public:
	HDRRenderer() :
		shader_manager(program_definitions)
	{}

	void resize(int width, int height);
	void destroy();
	void tonemap();
	void begin_hdr_rendering();
	void end_hdr_rendering();
	int  initialize();

	ShaderManager shader_manager;
	std::shared_ptr<GLObjModel> fullscreen_tri;
	std::shared_ptr<FBO> fbo;
	const static MapNameProgramDefinition program_definitions;
};

class SponzaRenderer : public Renderer
{
public:
	SponzaRenderer() :
		Renderer(program_definitions)
	{}

	const char * get_name() override { return "Sponza"; }

	int  _initialize() override;
	void _draw() override;
	void _activate() override;
	void _deactivate() override;
	void _destroy() override;
	void _update_objects(double time_step) override;
	void _resize(int width, int height) override;

	std::shared_ptr<GLObjModel> sponza;
	const static MapNameProgramDefinition program_definitions;

	float sponza_k_s = 0.5f;

	HDRRenderer hdrRenderer;

	enum { LIGHT_COUNT = 9 };
	glm::vec3 light_pos[LIGHT_COUNT];
	glm::vec3 light_intensity[LIGHT_COUNT];
};

class SphereFlakeRenderer : public Renderer
{
public:
	bool enable_animation = false;
private:
	std::shared_ptr<SceneGraphNode> sceneGraphRoot;
	std::shared_ptr<GLObjModel> sphereModel;
	float animationSpeedFactor = 0.1f;
	int sphere_flake_recursion_depth = 0;
	float sphere_flake_size_factor = 0.5f;

	std::vector<TransformedModel> transformedRenderingObjects;

	const static MapNameProgramDefinition program_definitions;
public:
	SphereFlakeRenderer() :
		Renderer(program_definitions)
	{}

	const char * get_name() override { return "Sphere Flake"; }

	virtual void _activate() override;
	virtual void _deactivate() override;
	virtual void _destroy() override;
	virtual int  _initialize() override;
	virtual void _draw() override;
	virtual void _update_objects(double time_step) override;
};
// CG_REVISION d4ab32bd208749f2d2b1439e25d16e642b039298
