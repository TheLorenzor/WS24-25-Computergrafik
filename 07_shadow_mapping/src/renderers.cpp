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
#include <cglib/imgui/imgui_orient.h>

#include "renderers.h"

const MapNameProgramDefinition ShadowmapRenderer::program_definitions = {
	{	"shadowmap",
		{ { "shader/shadowmap.vert" }
		, { "shader/shadowmap.frag" }
		, { }
		}
	},

	{	"visualize_depth_map",
		{ { "shader/visualize_depth.vert" }
		, { "shader/visualize_depth.frag" }
		, { }
		}
	},

	{	"fixed_color",
		{ { "shader/fixed_color.vert" }
		, { "shader/fixed_color.frag" }
		, { }
		}
	},
};

const MapNameProgramDefinition FirePlaceRenderer::program_definitions = {
	{	"simple",
		{ { "shader/simple.vert" }
		, { "shader/simple.frag" }
		, { }
		}
	},
	{	"logs",
		{ { "shader/logs.vert" }
		, { "shader/logs.frag" }
		, { }
		}
	},
	{	"smoke",
		{ { "shader/smoke.vert" }
		, { "shader/smoke.frag" }
		, { }
		}
	},
};

const MapNameProgramDefinition ProceduralLandscapeRenderer::program_definitions = {
	{	"heightmap",
		{ { "shader/heightmap.vert" }
		, { "shader/heightmap.frag" }
		, { }
		}
	},

	{	"wavemap",
		{ { "shader/wavemap.vert" }
		, { "shader/wavemap.frag" }
		, { }
		}
	},

	{	"terrain",
		{ { "shader/vertex_displacement.vert" }
		, { "shader/terrain.frag" }
		, { }
		}
	},

	{	"sky_bg",
		{ { "shader/sky_bg.vert" }
		, { "shader/sky_bg.frag" }
		, { }
		}
	},

	{	"water",
		{ { "shader/vertex_displacement.vert" }
		, { "shader/water.frag" }
		, { }
		}
	},

	{	"fire",
		{ { "shader/particle.vert" }
		, { "shader/fire.frag" }
		, { }
		, { "shader/particle.geom" }
		}
	},

	{	"shuttle",
		{ { "shader/shuttle.vert" }
		, { "shader/shuttle.frag" }
		, { }
		}
	},

	{	"visualize_texture",
		{ { "shader/visualize_texture.vert" }
		, { "shader/visualize_texture.frag" }
		}
	},
};

int ShadowmapRenderer::
_initialize()
{
	glm::vec3 cam_pos(22.2307, 25.586, 18.1663);
	glm::vec3 cam_view_dir(-0.630628, -0.625181, -0.459846);

	camera = std::make_shared<RTFreeFlightCamera>(cam_pos, cam_view_dir, 0.0f);
	donuts = std::make_shared<GLObjModel>("assets/donuts.obj");
	fullscreen_quad = std::make_shared<GLObjModel>("assets/fsquad.obj");
	fbo_shadow_map = std::make_shared<FBO>(
			FBOInit(1024, 1024).attach_depth(GL_DEPTH_COMPONENT32));

	glBindTexture(GL_TEXTURE_2D, fbo_shadow_map->tex_depth);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_BORDER);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_BORDER);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glBindTexture(GL_TEXTURE_2D, 0);

	create_cube_lines(&vao_cube, &buffer_cube_vert, &buffer_cube_idx);

	return 0;
}

void ShadowmapRenderer::
_destroy()
{
	glDeleteVertexArrays(1, &vao_cube);
	glDeleteBuffers(1, &buffer_cube_idx);
	glDeleteBuffers(1, &buffer_cube_vert);
	donuts = nullptr;
	fullscreen_quad = nullptr;
	fbo_shadow_map = nullptr;
}

void ShadowmapRenderer::
_draw()
{
	ImGui::DragFloat("Shadow Bias", &bias, 0.001f);
	ImGui::DragFloat("Light near clip plane", &light_near, 0.1f);
	ImGui::DragFloat("Light far clip plane", &light_far, 0.1f);
	ImGui::DragFloat("Light Fov", &light_fov, 1.f);
	ImGui::DirectionGizmo("Light Direction", reinterpret_cast<ImVec3 &>(light_direction));

	glEnable(GL_DEPTH_TEST);
	// compute the view_projection matrix for the light
	// does not contain model matrix since it will be used
	// later on to project world space positions into the
	// shadow map
	auto view_projection_light = compute_view_projection_light(
			light_position,
			light_direction,
			light_near,
			light_far,
			glm::radians(light_fov),
			glm::vec3(0.0f, 1.0f, 0.0f));

	glm::mat4 model_matrix(1.f); // identity matrix for our scene, we only have one model

	{	// render scene into shadow map
		fbo_shadow_map->bind();   // draw into the shadow map
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT); // clear the shadow map
		shader_manager["fixed_color"]
			.bind()
			.uniform("MVP", view_projection_light * model_matrix)
			;
		donuts->draw();           // draw the scene
		fbo_shadow_map->unbind(); // stop drawing into the shadow map
	}

	glEnable(GL_FRAMEBUFFER_SRGB);

	{	// draw the scene with shadows
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

		// bind the depth texture for shadow lookup to texture unit 1
		glActiveTexture(GL_TEXTURE1);
		glBindTexture(GL_TEXTURE_2D, fbo_shadow_map->tex_depth);
		// make texture unit 0 active again
		glActiveTexture(GL_TEXTURE0);
		shader_manager["shadowmap"]
			.bind()
			.uniform("MVP",
					projection
					* camera->get_view_matrix(Camera::Mono))
			.uniform("model_matrix", model_matrix)
			.uniform("normal_matrix", glm::mat3(glm::transpose(glm::inverse(model_matrix))))
			.uniform("light_matrix", view_projection_light)
			.uniform("world_light_position", light_position)
			.uniform("world_light_view_direction", light_direction)
			.uniform("shadow_map", 1)
			.uniform("spotlight_cos_outer", std::cos(glm::radians(0.5f * light_fov)))
			.uniform("spotlight_cos_inner", std::cos(glm::radians(0.5f * 2.0f)))
			.uniform("shadow_bias", bias)
			.uniform("light_emission", glm::vec3(450.0f))
			.uniform("K_d", glm::vec3(0.8f))
			;
		donuts->draw(); // draw the scene
	}

	{	// draw the light frustum
		shader_manager["fixed_color"]
			.bind()
			.uniform("MVP",
					projection
					* camera->get_view_matrix(Camera::Mono)
					* glm::inverse(view_projection_light)
					)
			.uniform("color", glm::vec4(1.0))
			;
		glBindVertexArray(vao_cube);
		glDrawElements(GL_LINES, 24, GL_UNSIGNED_INT, 0);
		glBindVertexArray(0);
	}

	glDisable(GL_FRAMEBUFFER_SRGB);

	{	// visualize the depth image
		shader_manager["visualize_depth_map"].bind()
			.uniform("shadow_map", 1)
			.uniform("near", light_near)
			.uniform("far",  light_far)
			.uniform("MVP",
					projection
					* camera->get_view_matrix(Camera::Mono)
					* glm::inverse(view_projection_light)
					* glm::translate(glm::vec3(0, 0, -1))
					)
			;
		glActiveTexture(GL_TEXTURE1);
		glBindTexture(GL_TEXTURE_2D, fbo_shadow_map->tex_depth);
		glActiveTexture(GL_TEXTURE0);
		fullscreen_quad->draw();
	}
}

int ProceduralLandscapeRenderer::
_initialize()
{
	glm::vec3 cam_pos(432.104f, 92.8843f, 969.301f);
	glm::vec3 cam_view_dir(-0.326702f, -0.0303592f, -0.94464f);

	camera = std::make_shared<RTFreeFlightCamera>(cam_pos, cam_view_dir, 0.0f);
	camera->m_speed *= 40.0f;

	fullscreen_tri = std::make_shared<GLObjModel>("assets/fstri.obj");
	shuttle = std::make_shared<GLObjModel>("assets/tintin_rocket_tri.obj");

	grassTex = loadGLTexture("assets/grass.png", 2.2f);
	rockTex  = loadGLTexture("assets/rock.tga",  2.2f);
	sandTex  = loadGLTexture("assets/sand.jpg",  2.2f);
	fireTex  = loadGLTexture("assets/fire.png",  2.2f);
	waterTex = loadGLTexture("assets/water.jpg", 1.0f);

	height_map = std::make_shared<FBO>(
			FBOInit(2048, 2048).attach_color(0, GL_R16F).use_color_mip());
	glActiveTexture(GL_TEXTURE0);
	glBindTexture(GL_TEXTURE_2D, height_map->tex_color[0]);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
	glBindTexture(GL_TEXTURE_2D, 0);

	wave_map = std::make_shared<FBO>(
			FBOInit(512, 512).attach_color(0, GL_R16F));
	glActiveTexture(GL_TEXTURE0);
	glBindTexture(GL_TEXTURE_2D, wave_map->tex_color[0]);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
	glBindTexture(GL_TEXTURE_2D, 0);

	reflection_map = std::make_shared<FBO>(
			FBOInit(1024, 1024).attach_color(0, GL_RGBA16F).attach_depth(GL_DEPTH_COMPONENT32));

	// Random image
	{
		int randomNumbersPerAxis = 256;
		std::vector<float> randImg(randomNumbersPerAxis * randomNumbersPerAxis);
		std::mt19937 mt;
		std::uniform_real_distribution<float> uniform_float;
		for (int y = 0; y < randomNumbersPerAxis; ++y) {
			for (int x = 0; x < randomNumbersPerAxis; ++x) {
				randImg[x + y * randomNumbersPerAxis] = uniform_float(mt);
			}
		}

		glGenTextures(1, &randTex);
		glBindTexture(GL_TEXTURE_2D, randTex);
		glTexImage2D(GL_TEXTURE_2D, 0, GL_R16, randomNumbersPerAxis, randomNumbersPerAxis, 0, GL_RED, GL_FLOAT, randImg.data());
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
		glBindTexture(GL_TEXTURE_2D, 0);
	}

	// Plane
	glGenVertexArrays(1, &vao_plane);
	glBindVertexArray(vao_plane);

	int plane_patch_vertices_per_axis = 256;
	{
		glGenBuffers(1, &buffer_plane_vert);
		glBindBuffer(GL_ARRAY_BUFFER, buffer_plane_vert);

		std::vector<glm::vec2> vertices(plane_patch_vertices_per_axis * plane_patch_vertices_per_axis);

		for (int i = 0; i < plane_patch_vertices_per_axis; ++i) {
			for (int j = 0; j < plane_patch_vertices_per_axis; ++j) {
				vertices[i + j * plane_patch_vertices_per_axis] = glm::vec2(i, j);
			}
		}

		glBufferData(GL_ARRAY_BUFFER, sizeof(*vertices.data()) * vertices.size(), vertices.data(), GL_STATIC_DRAW);
		glEnableVertexAttribArray(0);
		glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 0, 0);
	}

	{
		glGenBuffers(1, &buffer_plane_idx);
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, buffer_plane_idx);

		plane_patch_width = plane_patch_vertices_per_axis - 1;
		plane_idx_count = 6 * plane_patch_width * plane_patch_width;
		std::vector<glm::ivec3> indices(plane_idx_count / 3);

		for (int i = 0; i < plane_patch_width; ++i) {
			for (int j = 0; j < plane_patch_width; ++j) {
				int baseTri = 2 * (i + j * plane_patch_width);
				int baseVtx = i + j * plane_patch_vertices_per_axis;
				indices[baseTri] = glm::ivec3(baseVtx, baseVtx + plane_patch_vertices_per_axis, baseVtx + 1);
				indices[baseTri + 1] = glm::ivec3(baseVtx + 1, baseVtx + plane_patch_vertices_per_axis, baseVtx + plane_patch_vertices_per_axis + 1);
			}
		}

		glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(*indices.data()) * indices.size(), indices.data(), GL_STATIC_DRAW);
	}

	glBindVertexArray(0);

	// Smoke
	{
		std::mt19937 mt;
		std::uniform_real_distribution<float> uniform_float;

		auto startPos = glm::vec3(-449.921f, 145.108f, -701.308f);

		auto currentPos = startPos;
		auto currentDir = glm::vec3(0.0f, 1.0f, 0.0f);
		auto height = 800.0f;
		auto offset = 2.5f;
		auto steps = 5000;

		fire_particles.reserve(steps);

		for (int i = 0; i < steps; ++i)
		{
			// approx time since launch
			float age = 1.0f - float(i) / float(steps);
			age *= age;

			auto particlePos = currentPos;
			particlePos += age * currentPos.y * uniform_float(mt) * glm::vec3(2.0f * uniform_float(mt) - 1.0f, -1.0f, 2.0f * uniform_float(mt) - 1.0f);
			float particleRadius = (15.0f + 60.0f * uniform_float(mt)) * 0.55f;
			fire_particles.push_back(glm::vec4(particlePos, particleRadius));

			currentPos += currentDir * ((height - startPos.y) / steps);
			currentDir.x += offset / steps;
		}

		{
			glm::vec3 up = glm::normalize(currentDir);
			glm::vec3 right = glm::normalize(glm::cross(up, glm::vec3(0.0f, 0.0f, 1.0f)));
			glm::vec3 look = glm::cross(right, up);
			shuttle_transform = glm::translate(currentPos) * glm::scale(glm::vec3(5.0f)) * glm::mat4(glm::mat3(right, up, look));
		}

		glGenVertexArrays(1, &vao_fire);
		glBindVertexArray(vao_fire);

		glGenBuffers(1, &buffer_fire_vert);
		glBindBuffer(GL_ARRAY_BUFFER, buffer_fire_vert);
		glBufferData(GL_ARRAY_BUFFER, sizeof(*fire_particles.data()) * fire_particles.size(), fire_particles.data(), GL_DYNAMIC_DRAW);
		glEnableVertexAttribArray(0);
		glVertexAttribPointer(0, 4, GL_FLOAT, GL_FALSE, 0, 0);

		glBindVertexArray(0);
	}

	update_parameters();

	return 0;
}

void ProceduralLandscapeRenderer::
_destroy()
{
	glDeleteTextures(1, &randTex);
	glDeleteTextures(1, &grassTex);
	glDeleteTextures(1, &rockTex);
	glDeleteTextures(1, &sandTex);
	glDeleteTextures(1, &fireTex);
	glDeleteTextures(1, &waterTex);
	glDeleteVertexArrays(1, &vao_plane);
	glDeleteBuffers(1, &buffer_plane_idx);
	glDeleteBuffers(1, &buffer_plane_vert);
	glDeleteVertexArrays(1, &vao_fire);
	glDeleteBuffers(1, &buffer_fire_vert);

	fullscreen_tri = nullptr;
	shuttle = nullptr;
	height_map = nullptr;
	wave_map = nullptr;
	reflection_map = nullptr;
}

static inline float
sgn(float a)
{
    if (a > 0.0f) return 1.0f;
    if (a < 0.0f) return -1.0f;
    return 0.0F;
}

// http://www.terathon.com/code/oblique.html
static glm::mat4
add_clip_plane_to_proj_matrix(const glm::mat4& in_proj, const glm::vec4& clipPlane)
{
	glm::mat4 proj = in_proj;
	glm::vec4 q;

    // Calculate the clip-space corner point opposite the clipping plane
    // as (sgn(clipPlane.x), sgn(clipPlane.y), 1, 1) and
    // transform it into camera space by multiplying it
    // by the inverse of the projection matrix

    q.x = (sgn(clipPlane.x) + proj[2][0]) / proj[0][0];
    q.y = (sgn(clipPlane.y) + proj[2][1]) / proj[1][1];
    q.z = -1.0f;
    q.w = (1.0f + proj[2][2]) / proj[3][2];

    // Calculate the scaled plane vector
    auto c = clipPlane * (2.0f / glm::dot(clipPlane, q));

    // Replace the third row of the projection matrix
    proj[0][2] = c.x;
    proj[1][2] = c.y;
    proj[2][2] = c.z + 1.0f;
    proj[3][2] = c.w;

	return proj;
}

static int
highest_bit_pos(int x)
{
	int c = 0;
	while ((x = x >> 1)) ++c;
	return c;
}

void ProceduralLandscapeRenderer::
_draw()
{
	bool cfg_changed = false;

	if (ImGui::CollapsingHeader("Rendering Settings"))
	{
		cfg_changed |= ImGui::Checkbox("Heightmap Visualization", &visualize_heightmap);
		cfg_changed |= ImGui::Checkbox("Reflection Visualization", &visualize_reflection);
		cfg_changed |= ImGui::Checkbox("No Water Surface", &no_water_surface);
		cfg_changed |= ImGui::Checkbox("No Textures", &no_textures);
		cfg_changed |= ImGui::SliderInt("Reduce Quality", &landscape_quality_reduction, 0, 3);
		cfg_changed |= ImGui::Checkbox("Wireframe", &wireframe);
	}

	if (ImGui::CollapsingHeader("Scene Settings"))
	{
		cfg_changed |= ImGui::Checkbox("Night Time", &night_time);

		if (ImGui::CollapsingHeader("Water"))
		{
			cfg_changed |= ImGui::Checkbox("Enable Fresnel", &enable_fresnel);
			cfg_changed |= ImGui::DragFloat("Water Height", &water_height, 1.f);
			cfg_changed |= ImGui::DragFloat("Wave Scaling", &wave_scaling, 0.1f);
			cfg_changed |= ImGui::DragFloat("Reflection Perturbation", &reflection_perturbation, 0.01f, 0.f, 1.f);
			cfg_changed |= ImGui::ColorEdit3("Water Color", &water_color[0]);
			cfg_changed |= ImGui::DragFloat("Water Absorption", &water_absorption, 0.01f, 0.f);
			cfg_changed |= ImGui::DragFloat("Water Surface Absorption", &water_surface_absorption, 0.01f, 0.f);
		}

		if (ImGui::CollapsingHeader("Terrain"))
		{
			cfg_changed |= ImGui::DragFloat("Height Scaling", &height_scaling, 0.1f, 0.f);
			cfg_changed |= ImGui::DragFloat("Beach Margin", &beach_margin, 0.01f);
			cfg_changed |= ImGui::DragFloat("Min Rock Threshold", &min_rock_threshold, 0.01f, 0.f);
			cfg_changed |= ImGui::DragFloat("Rock Blend Margin", &rock_blend_margin, 0.01f, 0.f, 1.f);
			cfg_changed |= ImGui::DragFloat("Lacunarity", &lacunarity, 0.001f, 1.f);
			cfg_changed |= ImGui::DragFloat("Weight Factor", &weight_factor, 0.001f, 0.f);
			cfg_changed |= ImGui::DragFloat("Angle Factor", &terrain_angle_factor, 0.1f);
			cfg_changed |= ImGui::DragFloat("Valley Smoothness", &valley_smoothness, 0.001f, 0.f, 1.f);
		}
	}

	if (cfg_changed)
		update_parameters();


	// render height map
	glDisable(GL_DEPTH_TEST);
	glDisable(GL_CULL_FACE);
	glDisable(GL_FRAMEBUFFER_SRGB);

	// height field
	if (update_heightmap)
	{
		height_map->bind();

		glClear(GL_COLOR_BUFFER_BIT);

		glActiveTexture(GL_TEXTURE1);
		glBindTexture(GL_TEXTURE_2D, randTex);
		shader_manager["heightmap"].bind()
			.uniform("lacunarity", lacunarity)
			.uniform("weight_factor", weight_factor)
			.uniform("rotation_factor", derived_rotation_factor)
			.uniform("valley_smoothness", valley_smoothness)
			.uniform("resolution_divisor", resolution_divisor)
			.uniform("RandomTexture", 1);
		glActiveTexture(GL_TEXTURE0);
		fullscreen_tri->draw();

		height_map->unbind();

		// generate mipmaps
		glActiveTexture(GL_TEXTURE0);
		glBindTexture(GL_TEXTURE_2D, height_map->tex_color[0]);
		glGenerateMipmap(GL_TEXTURE_2D);
		glBindTexture(GL_TEXTURE_2D, 0);

		update_heightmap = false;
	}

	// wave field
	if (!no_water_surface)
	{
		wave_map->bind();

		glClear(GL_COLOR_BUFFER_BIT);

		glActiveTexture(GL_TEXTURE1);
		glBindTexture(GL_TEXTURE_2D, waterTex);
		shader_manager["wavemap"].bind()
			.uniform("time", time)
			.uniform("WaveTexture", 1);
		glActiveTexture(GL_TEXTURE0);
		fullscreen_tri->draw();

		wave_map->unbind();
	}

	glEnable(GL_FRAMEBUFFER_SRGB);

	auto view = camera->get_view_matrix(Camera::Mono);
	auto camPos = camera->get_position(Camera::Mono);
	bool under_water = camPos.y < water_height;

	// reflect cam below water plane
	auto worldReflect = glm::translate(glm::vec3(0.0f, 2.0f * water_height, 0.0f))
		* glm::mat4(glm::mat3(glm::vec3(1.0f, 0.0f, 0.0f), glm::vec3(0.0f, -1.0f, 0.0f), glm::vec3(0.0f, 0.0f, 1.0f)));
	auto reflected_view = view * worldReflect;
	auto reflection_plane = glm::vec4(0.0f, 1.0f, 0.0f, -water_height) * glm::affineInverse(reflected_view);
	auto reflected_proj = add_clip_plane_to_proj_matrix(projection, under_water ?  -reflection_plane : reflection_plane);

	if (!no_water_surface)
	{
		// draw reflected landscape
		reflection_map->bind();
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

		// reflected landscape has reversed winding order!
		glFrontFace(GL_CW);
		drawLandscape(reflected_view, reflected_proj, under_water);
		drawShuttle(reflected_view, reflected_proj);
		drawParticles(reflected_view, reflected_proj);
		glFrontFace(GL_CCW);

		reflection_map->unbind();
	}

	// draw primary landscape
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	drawLandscape(view, projection, under_water);
	drawShuttle(view, projection);

	// render water
	glEnable(GL_DEPTH_TEST);

	if (!no_water_surface)
	{
		glActiveTexture(GL_TEXTURE0);
		glBindTexture(GL_TEXTURE_2D, wave_map->tex_color[0]);
		glActiveTexture(GL_TEXTURE1);
		glBindTexture(GL_TEXTURE_2D, reflection_map->tex_color[0]);
		glActiveTexture(GL_TEXTURE0);

		glDisable(GL_CULL_FACE);
		glEnable(GL_BLEND);
		glBlendEquation(GL_FUNC_ADD);
		glBlendFunc(GL_ONE, GL_SRC1_COLOR);

		auto waterExtent = resolution_divisor * glm::vec2(patches_per_axis * plane_patch_width);
		auto model_matrix = glm::translate(glm::vec3(-0.5f * waterExtent.x, water_height, -0.5f * waterExtent.y)) * glm::scale(glm::vec3(resolution_divisor.x, 1.0f, resolution_divisor.y));
		auto normal_matrix = glm::inverseTranspose(glm::mat3(model_matrix));
		auto VP = projection * view;
		auto MVP = VP * model_matrix;

		shader_manager["water"]
			.bind()
			.uniform("instanced_patch_mask_and_shift", glm::ivec2(patches_per_axis.x - 1, highest_bit_pos(patches_per_axis.y)))
			.uniform("instanced_patch_size", float(plane_patch_width))
			.uniform("height_scaling", wave_scaling)
			.uniform("one_over_height_map_size", 1.0f / glm::vec2(wave_map->width, wave_map->height))
			.uniform("MVP", MVP)
			.uniform("RVP_to_tex", glm::translate(glm::vec3(.5f)) * glm::scale(glm::vec3(.5f)) * reflected_proj * reflected_view)
			.uniform("water_height", water_height)
			.uniform("model_matrix", model_matrix)
			.uniform("normal_matrix", normal_matrix)
			.uniform("cam_world_pos", camPos)
			.uniform("surface_transmission_color", derived_water_surface_color)
			.uniform("water_transmission_coeff", derived_water_absorption)
			.uniform("reflection_perturbation", reflection_perturbation)
			.uniform("HeightMap", 0)
			.uniform("ReflectionMap", 1)
			.uniform("enable_fresnel", int(enable_fresnel))
			;

		glBindVertexArray(vao_plane);
		glDrawElementsInstanced(GL_TRIANGLES, plane_idx_count, GL_UNSIGNED_INT, 0, patches_per_axis.x * patches_per_axis.y);
		glBindVertexArray(0);

		glDisable(GL_BLEND);
		glEnable(GL_CULL_FACE);
	}

	drawParticles(view, projection);

	glDisable(GL_DEPTH_TEST);
	glDisable(GL_CULL_FACE);

	if (visualize_heightmap)
	{
		glDisable(GL_FRAMEBUFFER_SRGB);

		glActiveTexture(GL_TEXTURE1);
		glBindTexture(GL_TEXTURE_2D, height_map->tex_color[0]);
		shader_manager["visualize_texture"].bind()
			.uniform("monochrome", 1)
			.uniform("Texture", 1);
		glActiveTexture(GL_TEXTURE0);
		fullscreen_tri->draw();
	}
	else if (visualize_reflection)
	{
		glEnable(GL_FRAMEBUFFER_SRGB);

		glActiveTexture(GL_TEXTURE1);
		glBindTexture(GL_TEXTURE_2D, reflection_map->tex_color[0]);
		shader_manager["visualize_texture"].bind()
			.uniform("monochrome", 0)
			.uniform("Texture", 1);
		glActiveTexture(GL_TEXTURE0);
		fullscreen_tri->draw();
	}

	glDisable(GL_FRAMEBUFFER_SRGB);
}

void ProceduralLandscapeRenderer::
drawLandscape(glm::mat4 const& view, glm::mat4 const& projection_, bool under_water)
{
	auto terrainExtent = resolution_divisor * glm::vec2(patches_per_axis * plane_patch_width);

	auto VP = projection_ * view;
	auto VPI = glm::inverse(VP);
	auto model_matrix = glm::translate(-0.5f * glm::vec3(terrainExtent.x, 0.0f, terrainExtent.y)) * glm::scale(glm::vec3(resolution_divisor.x, 1.0f, resolution_divisor.y));
	auto normal_matrix = glm::inverseTranspose(glm::mat3(model_matrix));
	auto MVP = VP * model_matrix;

	// render terrain
	glEnable(GL_DEPTH_TEST);
	
	if (wireframe)
	{
		glDisable(GL_CULL_FACE);
		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	}
	else
		glEnable(GL_CULL_FACE);

	glActiveTexture(GL_TEXTURE0);
	glBindTexture(GL_TEXTURE_2D, grassTex);
	glActiveTexture(GL_TEXTURE1);
	glBindTexture(GL_TEXTURE_2D, rockTex);
	glActiveTexture(GL_TEXTURE2);
	glBindTexture(GL_TEXTURE_2D, sandTex);
	glActiveTexture(GL_TEXTURE3);
	glBindTexture(GL_TEXTURE_2D, height_map->tex_color[0]);
	glActiveTexture(GL_TEXTURE0);

	auto camPos = glm::vec3(glm::affineInverse(view)[3]);

	shader_manager["terrain"]
		.bind()
		.uniform("instanced_patch_mask_and_shift", glm::ivec2(patches_per_axis.x - 1, highest_bit_pos(patches_per_axis.y)))
		.uniform("instanced_patch_size", float(plane_patch_width))
		.uniform("MVP", MVP)
		.uniform("model_matrix", model_matrix)
		.uniform("normal_matrix", normal_matrix)
		.uniform("cam_world_pos", camPos)
		.uniform("water_height", water_height)
		.uniform("under_water", int(under_water && camPos.y < water_height)) // need double check: only under water if submerged & primary cam
		.uniform("height_scaling", height_scaling)
		.uniform("one_over_height_map_size", 1.0f / glm::vec2(height_map->width, height_map->height))
		.uniform("beach_margin", beach_margin)
		.uniform("min_rock_threshold", min_rock_threshold)
		.uniform("rock_blend_margin", rock_blend_margin)
		.uniform("sun_world_dir", derived_sun_direction)
		.uniform("sun_color", derived_sun_color)
		.uniform("water_transmission_coeff", derived_water_absorption)
		.uniform("use_textures", int(!no_textures))
		.uniform("GrassTex", 0)
		.uniform("RockTex", 1)
		.uniform("SandTex", 2)
		.uniform("HeightMap", 3)
		;

	glBindVertexArray(vao_plane);
	glDrawElementsInstanced(GL_TRIANGLES, plane_idx_count, GL_UNSIGNED_INT, 0, patches_per_axis.x * patches_per_axis.y);
	glBindVertexArray(0);
	
	if (wireframe)
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

	// render sky
	glDepthMask(GL_FALSE);
	glDisable(GL_CULL_FACE);

	shader_manager["sky_bg"]
		.bind()
		.uniform("VPI", VPI)
		.uniform("sun_world_dir", derived_sun_direction)
		.uniform("sun_color", derived_sun_color)
		.uniform("night_time", int(night_time))
		.uniform("seed", unsigned(rand()));
	fullscreen_tri->draw();

	glEnable(GL_CULL_FACE);
	glDepthMask(GL_TRUE);
}

void ProceduralLandscapeRenderer::
drawParticles(glm::mat4 const& view, glm::mat4 const& projection_)
{
	glEnable(GL_CULL_FACE);
	glEnable(GL_DEPTH_TEST);
	glDepthMask(GL_FALSE);
	glEnable(GL_BLEND);
	glBlendEquation(GL_FUNC_SUBTRACT);
	glBlendFunc(GL_ZERO, GL_ZERO);
	initialize_premultiplied_blending();

	glActiveTexture(GL_TEXTURE0);
	glBindTexture(GL_TEXTURE_2D, fireTex);
	glActiveTexture(GL_TEXTURE0);

	auto VI = glm::affineInverse(view);
	auto cam_pos = glm::vec3(VI[3]);
	auto cam_dir = -glm::vec3(VI[2]);

	shader_manager["fire"]
		.bind()
		.uniform("MVP", projection_ * view)
		.uniform("cam_world_pos", glm::vec3(VI[3]))
		.uniform("Texture", 0)
		;

	std::sort(fire_particles.begin(), fire_particles.end(), [=](const glm::vec4& a, const glm::vec4& b) {
		return glm::dot(glm::vec3(a) - cam_pos, cam_dir) > glm::dot(glm::vec3(b) - cam_pos, cam_dir);
	});

	{
		glBindBuffer(GL_ARRAY_BUFFER, buffer_fire_vert);
		auto dataSize = sizeof(*fire_particles.data()) * fire_particles.size();
		glBufferSubData(GL_ARRAY_BUFFER, 0, dataSize, fire_particles.data());
	}

	auto particleCount = int(fire_particles.size());
	glBindVertexArray(vao_fire);
	glDrawArrays(GL_POINTS, 0, particleCount);
	glBindVertexArray(0);

	glDisable(GL_BLEND);
	glDepthMask(GL_TRUE);
}

void ProceduralLandscapeRenderer::
drawShuttle(glm::mat4 const& view, glm::mat4 const& projection_)
{
	// render terrain
	glEnable(GL_DEPTH_TEST);
	glEnable(GL_CULL_FACE);
//	glDisable(GL_CULL_FACE); // model broken

	shader_manager["shuttle"]
		.bind()
		.uniform("MVP", projection_ * view * shuttle_transform)
		.uniform("model_matrix", shuttle_transform)
		.uniform("normal_matrix", glm::inverseTranspose(glm::mat3(shuttle_transform)))
		.uniform("cam_world_pos", glm::vec3(glm::affineInverse(view)[3]))
		.uniform("sun_world_dir", derived_sun_direction)
		.uniform("sun_color", derived_sun_color)
		;

	glActiveTexture(GL_TEXTURE0);
	shuttle->draw();
}

void ProceduralLandscapeRenderer::
update_parameters()
{
	int int_resolution_shift = glm::max(landscape_quality_reduction, 0);
	auto full_quality_patches_per_axis = glm::ivec2(height_map->width, height_map->height) / plane_patch_width;
	auto lower_quality_patches_per_axis = glm::max(glm::ivec2(1), full_quality_patches_per_axis >> int_resolution_shift);
	patches_per_axis = lower_quality_patches_per_axis;
	resolution_divisor = glm::vec2(full_quality_patches_per_axis) / glm::vec2(lower_quality_patches_per_axis);

	derived_sun_direction = night_time ? -glm::vec3(0.122453f, 0.225675f, -0.966476f) : glm::normalize(glm::vec3(-0.5f, -0.8f, 0.9f));
	derived_sun_color = night_time ? 1.5f * glm::vec3(0.07f, 0.07f, 0.04f) : glm::vec3(1.0f);

	derived_water_surface_color = glm::pow(water_color, glm::vec3(water_surface_absorption));
	derived_water_absorption = glm::log(glm::pow(water_color, glm::vec3(-water_absorption)));
	derived_water_scattering = glm::vec3(0.1f); // todo?

	derived_rotation_factor = glm::mat2(glm::rotate(glm::radians(terrain_angle_factor), glm::vec3(0.0f, 0.0f, 1.0f)));
	update_heightmap = true;
}

void ProceduralLandscapeRenderer::
_update_objects(double time_step)
{
	time += time_step;
}

void ProceduralLandscapeRenderer::
_resize(int width, int height)
{
	// override default projection matrix
	projection = glm::perspective(
			glm::radians(fov),
			float(width) / float(height),
			1.0f, 2500.0f);
}

/*
 * Fire place renderer
 */
int FirePlaceRenderer::_initialize()
{
	glm::vec3 cam_pos(0.f, 20.f, 30.f);
	glm::vec3 cam_view_dir(glm::normalize(glm::vec3(0.f, -0.25f, -1.f)));
	camera = std::make_shared<RTFreeFlightCamera>(cam_pos, cam_view_dir, 0.0f);
	//camera->m_speed *= 40.0f;

	fullscreen_tri = std::make_shared<GLObjModel>("assets/fstri.obj");
	logs = std::make_shared<GLObjModel>("assets/logs.obj");

	init_particles();

	return 0;
}

void FirePlaceRenderer::init_particles()
{
	std::mt19937 mt;
	std::uniform_real_distribution<float> uniform_float;

	fire_tex = loadGLTexture("assets/fire.png", 2.2f);
	glBindTexture(GL_TEXTURE_2D, fire_tex);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_BORDER);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_BORDER);
	glBindTexture(GL_TEXTURE_2D, 0);

	// SMOKE
	particles[0].resize(num_smoke_particles);
	particle_pos[0].resize(num_smoke_particles);
	particle_data[0].resize(num_smoke_particles);

	for (auto it = particles[0].begin(); it != particles[0].end(); ++it)
	{
		it->max_t = 1.f + 0.25f * std::sin(uniform_float(mt) * 2.f * float(M_PI));
		it->max_t /= 1.3f;
		it->t = uniform_float(mt) * it->max_t;
		it->angle = uniform_float(mt) * 2.f * float(M_PI);
		it->rot_speed = std::sin(uniform_float(mt) * 2.f * float(M_PI));
		it->rand = uniform_float(mt);
		it->pos_offset = glm::vec3(0.f, 3.f, 0.f);
	}

	// FIRE
	particles[1].resize(num_fire_particles);
	particle_pos[1].resize(num_fire_particles);
	particle_data[1].resize(num_fire_particles);

	for (auto it = particles[1].begin(); it != particles[1].end(); ++it)
	{
		it->max_t = 3.f + 2.f * std::sin(uniform_float(mt) * 2.f * float(M_PI));
		it->max_t /= 8.0f;
		it->t = uniform_float(mt) * it->max_t;
		it->angle = uniform_float(mt) * 2.f * float(M_PI);
		it->rot_speed = 2.0f * std::sin(uniform_float(mt) * 2.f * float(M_PI));
		it->rand = uniform_float(mt);
		it->pos_offset = glm::vec3(2.f * uniform_float(mt), 3.f, 2.f * uniform_float(mt));
	}

	// SMOKE + FIRE
	particle_pos[2].resize(num_smoke_particles+num_fire_particles);
	particle_data[2].resize(num_smoke_particles+num_fire_particles);

	glGenVertexArrays(3, particle_vao);
	glGenBuffers(3, particle_pos_buffer);
	glGenBuffers(3, particle_data_buffer);
	for (int i = 0; i < 3; ++i)
	{
		glBindVertexArray(particle_vao[i]);

		// pos.
		glBindBuffer(GL_ARRAY_BUFFER, particle_pos_buffer[i]);
		auto const& pos = particle_pos[i];
		glBufferData(GL_ARRAY_BUFFER, sizeof(*pos.data()) * pos.size(), 
			pos.data(), GL_DYNAMIC_DRAW);
		glEnableVertexAttribArray(0);
		glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, 0);

		// data.
		glBindBuffer(GL_ARRAY_BUFFER, particle_data_buffer[i]);
		auto const& data = particle_data[i];
		glBufferData(GL_ARRAY_BUFFER, sizeof(*data.data()) * data.size(), 
			data.data(), GL_DYNAMIC_DRAW);
		glEnableVertexAttribArray(1);
		glVertexAttribPointer(1, 4, GL_FLOAT, GL_FALSE, 0, 0);
	}
	glBindVertexArray(0);

}

void FirePlaceRenderer::
_update_objects(double time_step)
{
	// Smoke first, then fire.
	float const speed[]   = { 0.25f, .3f };
	float const delta_t[] = { (float)(speed[0] * time_step), (float)(speed[1] * time_step) };
	float const max_t     = 1.f;

	// Update particle positions.
	for (int m = 0; m < 2; ++m)
	{
		for (std::size_t i = 0; i < particles[m].size(); ++i)
		{
			Particle& p = particles[m][i];
			p.t += delta_t[m];

			// Die / spawn a new particle.
			if (p.t > max_t * p.max_t)
			{
				p.t = std::max(0.f, p.t - max_t * p.max_t);
			}

			p.angle += time_step;
			p.age = p.t / (max_t * p.max_t);
			float const base_scale = 2.f;

			switch (m)
			{
				// Smoke
				case 0:
					p.scale = base_scale * 1.5f * (1.f + p.age * 2.f);
					p.heat = 0.f;
					p.position = glm::vec3(
						5.f * std::sin(p.t*3.f*p.rand) / std::max(0.00001f, p.t*5.f),
						17.f*p.t,
						0.f
					);
					break;

				// Fire
				case 1:
					p.scale = base_scale * (1.f-std::pow(p.age, 1.5f));
					p.heat = 1.f;
					p.position = glm::vec3(
						p.rand * 2.f * std::sin(7.0f * p.t*5.f*p.rand) / std::max(1.0f, p.t*10.f),
						9.0f * (2.f+p.rand) * p.t - 0.8f,
						0.f
					);
					break;
			}

			p.position += p.pos_offset;
			float const angle = p.angle * p.rot_speed;
			glm::mat4 const R = glm::rotate(angle, glm::vec3(0.f, 1.f, 0.f));
			p.position = glm::vec3(R * glm::vec4(p.position, 1.f));
		}
	}

	// Make combined array.
	particles[2].resize(particles[0].size() + particles[1].size());
	for (std::size_t i = 0; i < particles[0].size(); ++i)
	{
		particles[2][i] = particles[0][i];
	}
	for (std::size_t i = 0; i < particles[1].size(); ++i)
	{
		particles[2][i+particles[0].size()] = particles[1][i];
	}


	auto const cam_pos = camera->get_position(Camera::Mono);
	auto const cam_dir = glm::normalize(camera->get_direction());
	for (int m = 0; m < 3; ++m)
	{
		std::sort(particles[m].begin(), particles[m].end(), 
			[=](Particle const& p1, Particle const& p2) 
		{
			return glm::dot(glm::vec3(p1.position) - cam_pos, cam_dir) 
				> glm::dot(glm::vec3(p2.position) - cam_pos, cam_dir);
		});
	}
	
	// Update buffers.
	for (int m = 0; m < 3; ++m)
	{
		for (std::size_t i = 0; i < particles[m].size(); ++i)
		{
			Particle& p = particles[m][i];
			particle_pos[m][i] = p.position;
			particle_data[m][i] = glm::vec4(p.heat, p.age, p.rand, p.scale);
		}
	}


	// Upload buffers.
	for (int i = 0; i < 3; ++i)
	{
		{
			auto const& src = particle_pos[i];
			glBindBuffer(GL_ARRAY_BUFFER, particle_pos_buffer[i]);
			auto dataSize = sizeof(*src.data()) * src.size();
			glBufferSubData(GL_ARRAY_BUFFER, 0, dataSize, src.data());
		}

		{
			auto const& src = particle_data[i];
			glBindBuffer(GL_ARRAY_BUFFER, particle_data_buffer[i]);
			auto dataSize = sizeof(*src.data()) * src.size();
			glBufferSubData(GL_ARRAY_BUFFER, 0, dataSize, src.data());
		}
	}
}

void FirePlaceRenderer::
_resize(int width, int height)
{
	particle_px_size = std::max(float(width), float(2 * height));
}

void FirePlaceRenderer::
_destroy()
{
	fullscreen_tri = nullptr;
	logs = nullptr;

	destroy_particles();
}

void FirePlaceRenderer::destroy_particles()
{
	particles[0].clear();
	particles[1].clear();
	particle_pos[0].clear();
	particle_pos[1].clear();
	particle_pos[2].clear();
	particle_data[0].clear();
	particle_data[1].clear();
	particle_data[2].clear();
	glDeleteVertexArrays(3, particle_vao);
	glDeleteBuffers(3, particle_pos_buffer);
	glDeleteBuffers(3, particle_data_buffer);
	glDeleteTextures(1, &fire_tex);
}

void FirePlaceRenderer::
_draw()
{
	bool cfg_changed = false;
	cfg_changed |= ImGui::SliderInt("Fire Particles", &num_fire_particles, 0, 2048);
	cfg_changed |= ImGui::SliderInt("Smoke Particles", &num_smoke_particles, 0, 1024);
	if (cfg_changed)
	{
		destroy_particles();
		init_particles();
	}

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glEnable(GL_FRAMEBUFFER_SRGB);
	auto view = camera->get_view_matrix(Camera::Mono);

	glEnable(GL_DEPTH_TEST);
	glEnable(GL_CULL_FACE);

	drawLogs(view, projection);

	glDisable(GL_DEPTH_TEST);
	glDisable(GL_CULL_FACE);

	glEnable(GL_DEPTH_TEST);

	drawParticles(view, projection);

	glDisable(GL_DEPTH_TEST);
	glDisable(GL_FRAMEBUFFER_SRGB);
}


void FirePlaceRenderer::
drawLogs(
	glm::mat4 const& view, 
	glm::mat4 const& projection_)
{
	// render terrain
	glEnable(GL_DEPTH_TEST);
	glEnable(GL_CULL_FACE);

	glm::mat4 logs_transform = glm::translate(glm::vec3(-15.f, 0.f, 0.f))
				   * glm::rotate(1.5f, glm::vec3(0.f, 1.f, 0.f));
	shader_manager["logs"]
		.bind()
		.uniform("MVP", projection_ * view * logs_transform)
		.uniform("model_matrix", logs_transform)
		.uniform("normal_matrix", glm::inverseTranspose(glm::mat3(logs_transform)))
		;

	logs->draw();

	logs_transform = glm::translate(glm::vec3(0.f, 0.f, 0.f));
	shader_manager["logs"]
		.bind()
		.uniform("MVP", projection_ * view * logs_transform)
		.uniform("model_matrix", logs_transform)
		.uniform("normal_matrix", glm::inverseTranspose(glm::mat3(logs_transform)))
		;

	logs->draw();

	logs_transform = glm::translate(glm::vec3(15.f, 0.f, 0.f))
				   * glm::rotate(0.75f, glm::vec3(0.f, 1.f, 0.f));
	shader_manager["logs"]
		.bind()
		.uniform("MVP", projection_ * view * logs_transform)
		.uniform("model_matrix", logs_transform)
		.uniform("normal_matrix", glm::inverseTranspose(glm::mat3(logs_transform)))
		;

	logs->draw();
}

void FirePlaceRenderer::
drawParticles(
	glm::mat4 const& view,
	glm::mat4 const& projection_)
{
	float const offset[] = {-15.f, 0.f, 15.f};

	glActiveTexture(GL_TEXTURE0);
	glBindTexture(GL_TEXTURE_2D, fire_tex);

	for (int i = 0; i < 3; ++i)
	{
		auto M = glm::translate(glm::vec3(offset[i], 0.f, 0.f));
		auto MVP = projection_ * view * M;

		glDisable(GL_CULL_FACE);
		glDepthMask(GL_FALSE);
		glEnable(GL_DEPTH_TEST);
		glEnable(GL_PROGRAM_POINT_SIZE);

		glEnable(GL_BLEND);
		glBlendEquation(GL_FUNC_SUBTRACT);
		glBlendFunc(GL_ZERO, GL_ZERO);

		switch(i) {
		case 0:  initialize_alpha_blending();         break;
		case 1:  initialize_additive_blending();      break;
		default: initialize_premultiplied_blending(); break;
		}

		shader_manager["smoke"]
			.bind()
			.uniform("MVP", MVP)
			.uniform("fire_texture", 0)
			.uniform("particle_px_size", particle_px_size)
			.uniform("render_mode", i)
			;

		glBindVertexArray(particle_vao[i]);
		glDrawArrays(GL_POINTS, 0, static_cast<int>(particle_pos[i].size()));
		glBindVertexArray(0);
		glDisable(GL_BLEND);
		glDisable(GL_PROGRAM_POINT_SIZE);
		glDisable(GL_DEPTH_TEST);
		glDepthMask(GL_TRUE);
	}

	glBindTexture(GL_TEXTURE_2D, 0);
}
// CG_REVISION d4ab32bd208749f2d2b1439e25d16e642b039298
