#include <cglib/core/assert.h>
#include <cglib/gl/fbo.h>

#include <cglib/gl/renderer.h>
#include <cglib/gl/glmodel.h>
#include <cglib/rt/aabb.h>
#include <cglib/rt/texture.h>
#include <cglib/core/image.h>
#include <cglib/gl/prefilter_envmap.h>
#include <cglib/gl/scene_graph.h>
#include <cglib/gl/util.h>

#include <glm/gtx/transform.hpp>
#include <glm/gtc/matrix_inverse.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <cglib/core/glmstream.h>

#include <random>

#include <cglib/imgui/imgui.h>

#include "renderers.h"

const MapNameProgramDefinition CubeRenderer::program_definitions = {
	{	"simple",
		{ { "shader/simple.vert" }
		, { "shader/simple.frag" }
		, { "gl_Position", "color" }
		}
	},
};

const MapNameProgramDefinition MonkeyRenderer::program_definitions = {
	{	"env_light",
		{ { "shader/env_light.vert" }
		, { "shader/env_light.frag" }
		, { "world_position", "world_normal_interpolated", "gl_Position" }
		}
	},
	{	"env_bg",
		{ { "shader/env_bg.vert" }
		, { "shader/env_bg.frag" }
		, { "viewdir_interpolated", "gl_Position" }
		}
	},
};

const MapNameProgramDefinition HDRRenderer::program_definitions = {
	{	"tonemap",
		{ { "shader/tonemap.vert" }
		, { "shader/tonemap.frag" }
		}
	},
};

const MapNameProgramDefinition SponzaRenderer::program_definitions = {
	{	"simple_texture",
		{ { "shader/simple_texture.vert" }
		, { "shader/simple_texture.frag" }
		, { "normal", "gl_Position" }
		}
	},
	{	"phong",
		{ { "shader/phong.vert" }
		, { "shader/phong.frag" }
		, { "world_position", "world_normal_interpolated", "tex_coord", "gl_Position" }
		}
	},
};

const MapNameProgramDefinition SphereFlakeRenderer::program_definitions = {
	{	"simple",
		{ { "shader/simple.vert" }
		, { "shader/simple_tex.frag" }
		, { }
		}
	},
};

int CubeRenderer::
_initialize()
{
	create_cube(&vao_cube, &buffer_cube_vert, &buffer_cube_idx);

	glm::vec3 cam_pos(3.0f, 3.f, 3.0f);
	glm::vec3 cam_center(0.f, 0.f, 0.0f);

	camera = std::make_shared<LookAroundCamera>(cam_pos, cam_center, 0.0f);

	return 0;
}

void CubeRenderer::
_destroy()
{
	glDeleteVertexArrays(1, &vao_cube);
	glDeleteBuffers(1, &buffer_cube_vert);
	glDeleteBuffers(1, &buffer_cube_idx);
}

void CubeRenderer::
_draw()
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	// bind shader and set uniforms
	shader_manager["simple"]
		.bind()
		.uniform("MVP", projection * camera->get_view_matrix(Camera::Mono))
		;

	glBindVertexArray(vao_cube);
	glDrawElements(GL_TRIANGLES, 36, GL_UNSIGNED_INT, 0);
	glBindVertexArray(0);

}

void CubeRenderer::
_activate()
{
	glEnable(GL_DEPTH_TEST);
	glDisable(GL_CULL_FACE);
}

void CubeRenderer::
_deactivate()
{
	glDisable(GL_DEPTH_TEST);
}

int MonkeyRenderer::
_initialize()
{
	glm::vec3 cam_pos(0.0f, 0.1f, -3.0f);
	glm::vec3 cam_center(0.f, 0.f, 0.5f);

	camera = std::make_shared<LookAroundCamera>(cam_pos, cam_center, 0.0f);
	monkey = std::make_shared<GLObjModel>("assets/suzanne.obj");
	fullscreen_tri = std::make_shared<GLObjModel>("assets/fstri.obj");

	{
		Image img;
		img.load("assets/appartment.jpg", 2.2f);

		ImageTexture tex(img, TextureFilterMode::TRILINEAR, TextureWrapMode::CLAMP);
		{
			tex.create_mipmap();
			glGenTextures(1, &envTex);

			glBindTexture(GL_TEXTURE_CUBE_MAP, envTex);

			resample_to_cubemap(tex);

			glGenerateMipmap(GL_TEXTURE_CUBE_MAP);
			glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
			glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
		}

		{
			auto reflf = prefilter_environment_specular(*tex.get_mip_levels()[5], 40.0f);
			auto& img_ = *reflf;
			img_.save("prefiltered_envmap_specular.png", 2.2f);
			glGenTextures(1, &envTexGlossy);

			glBindTexture(GL_TEXTURE_CUBE_MAP, envTexGlossy);

			ImageTexture reflftex(*reflf, TextureFilterMode::BILINEAR, TextureWrapMode::CLAMP);
			resample_to_cubemap(reflftex);

			glGenerateMipmap(GL_TEXTURE_CUBE_MAP);
			glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
			glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
		}

		{
			auto difff = prefilter_environment_diffuse(*tex.get_mip_levels()[6]);
			auto& img_ = *difff;
			img_.save("prefiltered_envmap_diffuse.png", 2.2f);
			glGenTextures(1, &envTexDiffuse);

			glBindTexture(GL_TEXTURE_CUBE_MAP, envTexDiffuse);

			ImageTexture diffftex(*difff, TextureFilterMode::BILINEAR, TextureWrapMode::CLAMP);
			resample_to_cubemap(diffftex);

			glGenerateMipmap(GL_TEXTURE_CUBE_MAP);
			glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
			glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
		}
	}

	return 0;
}

void MonkeyRenderer::
_draw()
{
	glEnable(GL_FRAMEBUFFER_SRGB);
	glEnable(GL_TEXTURE_CUBE_MAP_SEAMLESS);

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	// bind environment textures
	glActiveTexture(GL_TEXTURE1);
	glBindTexture(GL_TEXTURE_CUBE_MAP, envTex);
	glActiveTexture(GL_TEXTURE2);
	glBindTexture(GL_TEXTURE_CUBE_MAP, envTexGlossy);
	glActiveTexture(GL_TEXTURE3);
	glBindTexture(GL_TEXTURE_CUBE_MAP, envTexDiffuse);
	glActiveTexture(GL_TEXTURE0);

	auto VP = projection * camera->get_view_matrix(Camera::Mono);
	auto VPI = glm::inverse(VP);

	glDisable(GL_DEPTH_TEST);
	glDisable(GL_CULL_FACE);

	shader_manager["env_bg"]
		.bind()
		.uniform("VPI", VPI)
		.uniform("EnvironmentTexture", 1)
		;

	fullscreen_tri->draw();

	glEnable(GL_DEPTH_TEST);
	glEnable(GL_CULL_FACE);

	float diff_norm = 1.0f / float(M_PI);
	float spec_norm = (40.0f + 2.0f) / (2.0f * float(M_PI));

	// bind shader and set uniforms
	auto&& shader_transform = [&](Shader& shader, glm::mat4 const& model_matrix, glm::vec3 diffuse, glm::vec3 specular)
	{
		shader
			.uniform("MVP", VP * model_matrix)
			.uniform("model_matrix", model_matrix)
			.uniform("normal_matrix", glm::inverseTranspose(glm::mat3(model_matrix)))
			.uniform("diffuse_color", diff_norm * diffuse)
			.uniform("specular_color", spec_norm * specular)
			;
	};

	auto &env_light = shader_manager["env_light"];
	env_light
		.bind()
		.uniform("cam_world_pos", camera->get_position(Camera::Mono))
		.uniform("EnvironmentTextureGlossy", 2)
		.uniform("EnvironmentTextureDiffuse", 3)
		;

	// glossy
	shader_transform(env_light, glm::translate(glm::mat4(1.f), glm::vec3(3.0f, 0.0f, 0.0f)) * glm::rotate(glm::mat4(1.f),float(M_PI), glm::vec3(0.0f, 1.0f, 0.0f)), glm::vec3(0.0f), glm::vec3(0.7f, 0.6f, 0.5f));
	monkey->draw();

	// fresnel
	shader_transform(env_light, glm::translate(glm::mat4(1.f),glm::vec3(0.0f, 0.0f, 0.0f)) * glm::rotate(glm::mat4(1.f),float(M_PI), glm::vec3(0.0f, 1.0f, 0.0f)), glm::vec3(1.0f, 0.2f, 0.1f), glm::vec3(0.7f, 0.6f, 0.5f));
	monkey->draw();

	// diffuse
	shader_transform(env_light, glm::translate(glm::mat4(1.f),glm::vec3(-3.0f, 0.0f, 0.0f)) * glm::rotate(glm::mat4(1.f),float(M_PI), glm::vec3(0.0f, 1.0f, 0.0f)), glm::vec3(1.0f, 0.2f, 0.1f), glm::vec3(0.0f));
	monkey->draw();

	glDisable(GL_FRAMEBUFFER_SRGB);
}

void MonkeyRenderer::
_activate()
{
	glEnable(GL_DEPTH_TEST);
}

void MonkeyRenderer::
_deactivate()
{
	glDisable(GL_DEPTH_TEST);
}

void MonkeyRenderer::
_destroy()
{
	glDeleteTextures(1, &envTex);
	glDeleteTextures(1, &envTexGlossy);
	glDeleteTextures(1, &envTexDiffuse);

	monkey = nullptr;
	fullscreen_tri = nullptr;
}

int HDRRenderer::
initialize()
{
	fullscreen_tri = std::make_shared<GLObjModel>("assets/fstri.obj");
	return 0;
}

void HDRRenderer::
destroy()
{
	fullscreen_tri = nullptr;
}

void HDRRenderer::
resize(int width, int height)
{
	FBOInit fbodesc(width, height);
	fbodesc.attach_color(0, GL_RGBA16F);
	fbodesc.attach_depth(GL_DEPTH_COMPONENT32);
	fbo = std::make_shared<FBO>(fbodesc);
}

void HDRRenderer::
tonemap()
{
	glDisable(GL_DEPTH_TEST);
	glDisable(GL_CULL_FACE);

	auto was_bound = fbo->is_bound;
	if (was_bound) fbo->unbind();

	glActiveTexture(GL_TEXTURE1);
	glBindTexture(GL_TEXTURE_2D, fbo->tex_color[0]);
	shader_manager["tonemap"].bind()
		.uniform("HDRTexture", 1);
	glActiveTexture(GL_TEXTURE0);
	fullscreen_tri->draw();

	glBindTexture(GL_TEXTURE_2D, 0);
	if (was_bound) fbo->bind();

	glEnable(GL_DEPTH_TEST);
	glEnable(GL_CULL_FACE);
}

void HDRRenderer::
begin_hdr_rendering()
{
	fbo->bind();
}

void HDRRenderer::
end_hdr_rendering()
{
	fbo->unbind();
}

void SponzaRenderer::
_activate()
{
	glEnable(GL_DEPTH_TEST);
	glEnable(GL_CULL_FACE);
}

void SponzaRenderer::
_deactivate()
{
	glDisable(GL_DEPTH_TEST);
	glDisable(GL_CULL_FACE);
}

int SponzaRenderer::
_initialize()
{
	glm::vec3 cam_pos(-3.50705, 7.71856, 0.784408);
	glm::vec3 cam_view_dir(0.931884, -0.297391, -0.207728);

	camera = std::make_shared<RTFreeFlightCamera>(cam_pos, cam_view_dir, 0.0f);
	sponza = std::make_shared<GLObjModel>("assets/crytek-sponza/sponza_triag.obj");

	hdrRenderer.initialize();

	return 0;
}

void SponzaRenderer::
_destroy()
{
	hdrRenderer.destroy();
	sponza = nullptr;
}

void SponzaRenderer::
_resize(int width, int height)
{
	hdrRenderer.resize(width, height);
}

void SponzaRenderer::
_draw()
{
	enable_animation_prev = enable_animation;
	if (interactive)
	{
		ImGui::Text("Sponza");
		ImGui::Checkbox("Animate", &enable_animation);
		if (!enable_animation)
		{
			if (enable_animation_prev)
				slider_time = 28.5f;
			ImGui::DragFloat("Time", &slider_time, 0.1f, 0.f, 100.f);
		}
		else
		{
			slider_time = current_time;
		}
	}
	else if (!enable_animation) // only in non-interactive mode!
	{
		current_time = 20.5f;
	}

	glEnable(GL_DEPTH_TEST);
	glEnable(GL_CULL_FACE);
	glEnable(GL_FRAMEBUFFER_SRGB);

	hdrRenderer.begin_hdr_rendering();

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	auto VP = projection * camera->get_view_matrix(Camera::Mono);
	auto model_matrix = glm::scale(glm::mat4(1.f),glm::vec3(0.01f));
	auto normal_matrix = glm::inverseTranspose(glm::mat3(model_matrix));
	auto MVP = VP * model_matrix;

	// bind shader and set uniforms
	shader_manager["phong"]
		.bind()
		.uniform("MVP", MVP)
		.uniform("model_matrix", model_matrix)
		.uniform("normal_matrix", normal_matrix)
		.uniform("cam_world_pos", camera->get_position(Camera::Mono))
		.uniform("light_count", int(LIGHT_COUNT))
		.uniform("light_world_pos[0]", light_pos, LIGHT_COUNT)
		.uniform("light_intensity[0]", light_intensity, LIGHT_COUNT)
		.uniform("k_s", glm::vec3(0.4f))
		.uniform("n", 80.0f)
		;

	glActiveTexture(GL_TEXTURE0);
	sponza->draw();


	hdrRenderer.end_hdr_rendering();
	hdrRenderer.tonemap();

	glDisable(GL_FRAMEBUFFER_SRGB);
}

void SponzaRenderer::
_update_objects(double time_step)
{
	if (interactive)
	{
		if (!enable_animation)
		{
			current_time = slider_time;
		}
		enable_animation_prev = enable_animation;
	}
	else
	{
		current_time = 20.5;
	}

	for (int i = 0; i < LIGHT_COUNT; ++i)
	{
		float fi = (float) i;
		light_pos[i] = glm::vec3( 8.0f * sin(0.05f * (1.0f + fi) * current_time + 5.0f * fi), 2.6f + 2.4f * sin(0.02f * (1.0f + fi) * current_time + 3.0f * fi), 2.1f * sin(0.07f * (1.0f + fi) * current_time + 7.f * fi) );
		light_intensity[i] = 1.5f * (glm::vec3(0.5f) + 0.5f * glm::vec3( sin(fi * 2.f), sin(fi * 4.f), sin(fi * 3.f) ));
	}

	int j = 0;
	for (int i = 0; i < LIGHT_COUNT / 3; ++i, ++j) {
		float fi = (float) i;
		float angle = sin( 0.025f * (1.0f + fi) * current_time ) * float(M_PI) + 5.0f * fi;
		light_pos[j] = glm::vec3( 10.6f * sin(angle), 1.8f + 1.6f * sin(0.2f * (1.0f + fi) * current_time + 3.0f * fi), 4.9f * cos(angle) );
	}
	for (int i = 0; i < LIGHT_COUNT / 3; ++i, ++j) {
		float fi = (float) i;
		float angle = sin( 0.02f * (5.0f + fi) * current_time) * float(M_PI) + 4.0f * fi ;
		light_pos[j] = glm::vec3( 11.4f * sin(angle), 6.5f + 1.7f * sin(0.3f * (1.0f + fi) * current_time + 3.0f * fi), 5.5f * cos(angle) );
	}
}

void SphereFlakeRenderer::
_activate()
{
	glEnable(GL_DEPTH_TEST);
	glEnable(GL_CULL_FACE);
}

void SphereFlakeRenderer::
_deactivate()
{
	glDisable(GL_DEPTH_TEST);
	glDisable(GL_CULL_FACE);
}

int SphereFlakeRenderer::
_initialize()
{
	glm::vec3 cam_pos(-3.82484f, 2.31422f, -1.91606f);
	glm::vec3 cam_view_dir(0.792147f, -0.478072f, 0.379408f);

	camera = std::make_shared<RTFreeFlightCamera>(cam_pos, cam_view_dir, 0.0f);

	sphereModel = std::make_shared<GLObjModel>("assets/sphere.obj");

	sceneGraphRoot = buildSphereFlakeSceneGraph(sphereModel, 1.0f, sphere_flake_size_factor, sphere_flake_recursion_depth);

	return 0;
}

void SphereFlakeRenderer::
_destroy()
{
	sceneGraphRoot = nullptr;
	sphereModel = nullptr;
}

void SphereFlakeRenderer::
_draw()
{
	if (interactive)
	{
		ImGui::Text("Sphere Flake");
		ImGui::Checkbox("Animate", &enable_animation);
		if(ImGui::SliderInt("recursion depth", &sphere_flake_recursion_depth, 0, 5)) {
			sceneGraphRoot = buildSphereFlakeSceneGraph(sphereModel, 1.0f, sphere_flake_size_factor, sphere_flake_recursion_depth);
		}
	}


	glEnable(GL_DEPTH_TEST);
	glEnable(GL_CULL_FACE);
	glEnable(GL_FRAMEBUFFER_SRGB);
	GLERROR;

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	cg_assert(camera);

	transformedRenderingObjects.clear();
	sceneGraphRoot->collectTransformedModels(transformedRenderingObjects);

	glActiveTexture(GL_TEXTURE0);

	for (const auto& obj : transformedRenderingObjects) {
		shader_manager["simple"]
			.bind()
			.uniform("MVP",
				projection
				* camera->get_view_matrix(Camera::Mono)
				* obj.object_to_world)
			.uniform("tex", 0)
			;
		obj.model->draw();
	}

	glDisable(GL_FRAMEBUFFER_SRGB);
}

void SphereFlakeRenderer::
_update_objects(double time_step)
{
	if (enable_animation)
		animateSphereFlake(*sceneGraphRoot, time_step * animationSpeedFactor);
}


// CG_REVISION d4ab32bd208749f2d2b1439e25d16e642b039298
