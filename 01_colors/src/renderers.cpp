#include <cglib/colors/exercise.h>
#include <cglib/colors/physics.h>
#include <cglib/colors/convert.h>
#include <cglib/colors/cmf.h>
#include <glm/gtx/transform.hpp>
#include <glm/gtc/matrix_inverse.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <cglib/colors/font.h>
#include <cglib/imgui/imgui.h>
#include "renderers.h"

RGBCubeRenderer::
RGBCubeRenderer()
{
	cube_vertices = std::vector < glm::vec3 > {
		glm::vec3(0, 0, 0), glm::vec3(0, 0, 1), glm::vec3(0, 1, 1),
        glm::vec3(0, 0, 0), glm::vec3(0, 1, 1), glm::vec3(0, 1, 0),
        glm::vec3(0, 0, 0), glm::vec3(1, 0, 1), glm::vec3(0, 0, 1),
        glm::vec3(0, 0, 0), glm::vec3(1, 0, 0), glm::vec3(1, 0, 1),
        glm::vec3(0, 0, 0), glm::vec3(0, 1, 0), glm::vec3(1, 1, 0),
        glm::vec3(0, 0, 0), glm::vec3(1, 1, 0), glm::vec3(1, 0, 0),
        glm::vec3(1, 1, 1), glm::vec3(1, 1, 0), glm::vec3(0, 1, 0),
        glm::vec3(1, 1, 1), glm::vec3(0, 1, 0), glm::vec3(0, 1, 1),
        glm::vec3(1, 1, 1), glm::vec3(0, 1, 1), glm::vec3(0, 0, 1),
        glm::vec3(1, 1, 1), glm::vec3(0, 0, 1), glm::vec3(1, 0, 1),
        glm::vec3(1, 1, 1), glm::vec3(1, 0, 1), glm::vec3(1, 0, 0),
        glm::vec3(1, 1, 1), glm::vec3(1, 0, 0), glm::vec3(1, 1, 0)
    };

	for (auto& a : cube_vertices) a += glm::vec3(-1.0, -1.0, 0);

	cube_colors = std::vector < glm::vec3 > {
		glm::vec3(0, 0, 0), glm::vec3(0, 0, 1), glm::vec3(0, 1, 1),
        glm::vec3(0, 0, 0), glm::vec3(0, 1, 1), glm::vec3(0, 1, 0),
        glm::vec3(0, 0, 0), glm::vec3(1, 0, 1), glm::vec3(0, 0, 1),
        glm::vec3(0, 0, 0), glm::vec3(1, 0, 0), glm::vec3(1, 0, 1),
        glm::vec3(0, 0, 0), glm::vec3(0, 1, 0), glm::vec3(1, 1, 0),
        glm::vec3(0, 0, 0), glm::vec3(1, 1, 0), glm::vec3(1, 0, 0),
        glm::vec3(1, 1, 1), glm::vec3(1, 1, 0), glm::vec3(0, 1, 0),
        glm::vec3(1, 1, 1), glm::vec3(0, 1, 0), glm::vec3(0, 1, 1),
        glm::vec3(1, 1, 1), glm::vec3(0, 1, 1), glm::vec3(0, 0, 1),
        glm::vec3(1, 1, 1), glm::vec3(0, 0, 1), glm::vec3(1, 0, 1),
        glm::vec3(1, 1, 1), glm::vec3(1, 0, 1), glm::vec3(1, 0, 0),
        glm::vec3(1, 1, 1), glm::vec3(1, 0, 0), glm::vec3(1, 1, 0)
    };

	cube_wire_colors = std::vector <glm::vec3>(std::vector<glm::vec3>(cube_vertices.size(), glm::vec3(0.f)));

	glm::vec3 d(1.f, 1.f, -1.f);
	glm::vec3 cam_center(-1.f, -1.f, 0.0f);
	glm::vec3 cam_pos = cam_center + 1.0f * d;
	camera = std::make_shared<LookAroundCamera>(cam_pos, cam_center, 0.0f);
}

void RGBCubeRenderer::
_draw()
{
	ImGui::Checkbox("Show Wireframe", &showWireframe);
	// Set up the projection matrix. We will never change the camera after this
	// (rotation is done by rotating the world), so we can do this here.
	glMatrixMode(GL_PROJECTION);
	glLoadMatrixf((const float *)&projection);

	// From now on, matrix operations shall affect the model-view matrix.
	glMatrixMode(GL_MODELVIEW);

	// We generally want to use the z-buffer.
	glEnable(GL_DEPTH_TEST);

	glClearColor(0.01f, 0.01f, 0.01f, 1.f);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	// Enable Gamma Correction
	glEnable(GL_FRAMEBUFFER_SRGB);

	// Store current transformation to avoid side effects.
	glPushMatrix();

	//enforce_aspect_ratio(params.screen_width, params.screen_height);

	// Rotate cube according to user input on the GUI.
	glMultMatrixf(glm::value_ptr(camera->get_view_matrix(Camera::Mono)));
	
	// center the cube
	glTranslatef(-0.5, -0.5, -0.5);

	// Draw the cube as a simple list of triangles. We use the vertex coordinates
	// as colors here to show a visualization of the RGB color space.
	draw_triangles(cube_vertices, cube_colors);

	// Draw a black wireframe on top.
	if (showWireframe)
	{
		wireframe_overlay([&]() {
			draw_triangles(cube_vertices, cube_wire_colors);
		});
	}

	// Restore transformation.
	glPopMatrix();

	// Disable Gamma Correction
	glDisable(GL_FRAMEBUFFER_SRGB);
}

CGFontRenderer::CGFontRenderer()
{
	build_font_triangles(font_c_vertices, font_d_vertices);

	font_c_vertices_colors.resize(font_c_vertices.size());
	font_d_vertices_colors.resize(font_d_vertices.size());

    for (size_t i = 0; i < font_c_vertices.size(); ++i) {
		font_c_wire_colors.push_back(glm::vec3(1));
    }

    for (size_t i = 0; i < font_d_vertices.size(); ++i) {
		font_d_wire_colors.push_back(glm::vec3(1));
    }

	glm::vec3 d(1.f, 1.f, 1.f);
	glm::vec3 cam_center(0.0f, 0.0f, 0.0f);
	glm::vec3 cam_pos = cam_center + 0.05f * d;
	camera = std::make_shared<LookAroundCamera>(cam_pos, cam_center, 0.0f);
}

void CGFontRenderer::
_draw()
{
	ImGui::Checkbox("Show Wireframe", &showWireframe);

	// Set up the projection matrix. We will never change the camera after this
	// (rotation is done by rotating the world), so we can do this here.
	glMatrixMode(GL_PROJECTION);
	glLoadMatrixf((const float *)&projection);

	// From now on, matrix operations shall affect the model-view matrix.
	glMatrixMode(GL_MODELVIEW);

	// We generally want to use the z-buffer.
	glEnable(GL_DEPTH_TEST);

	glClearColor(0.01f, 0.01f, 0.01f, 1.f);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	// Disable Gamma Correction
	glEnable(GL_FRAMEBUFFER_SRGB);

	// Store current transformation to avoid side effects.
	glPushMatrix();

	//enforce_aspect_ratio(params.screen_width, params.screen_height);

	// Move a bit back so that we can see the whole cube.
	glTranslatef(0, 0, -1);

	// Rotate cube according to user input on the GUI.
	glMultMatrixf(glm::value_ptr(camera->get_view_matrix(Camera::Mono)));

	update_vertex_light_color(font_c_vertices, font_c_vertices_colors, camera->get_position(Camera::Mono), glm::vec3(1.f));
	update_vertex_light_color(font_d_vertices, font_d_vertices_colors, camera->get_position(Camera::Mono), glm::vec3(0.4f, 0.4f, 1.f));

	// Draw the cube as a simple list of triangles. We use the vertex coordinates
	// as colors here to show a visualization of the RGB color space.
	draw_triangles(font_c_vertices, font_c_vertices_colors);
	draw_triangles(font_d_vertices, font_d_vertices_colors);

	// Draw a black wireframe on top.
	if (showWireframe)
	{
		wireframe_overlay([&]() {
			draw_triangles(font_c_vertices, font_c_wire_colors);
			draw_triangles(font_d_vertices, font_d_wire_colors);
		});
	}

	// Restore transformation.
	glPopMatrix();

	// Enable Gamma Correction
	glDisable(GL_FRAMEBUFFER_SRGB);
}

void CGFontRenderer::
update_vertex_light_color(
		const std::vector<glm::vec3>& vertices,
		std::vector<glm::vec3>& colors,
		const glm::vec3& light_position,
	 	const glm::vec3& color) const
{
	cg_assert(vertices.size() == colors.size());
	cg_assert(vertices.size() % 3 == 0);

	glm::vec3 norm_light_position = glm::normalize(light_position);

	for (size_t i = 0; i < vertices.size(); i += 3) {
		glm::vec3 a = vertices[i + 0];
		glm::vec3 b = vertices[i + 1];
		glm::vec3 c = vertices[i + 2];

		glm::vec3 normal = glm::normalize(glm::cross(c - a, b - a));

		float diffuse_a = glm::dot(normal, a - norm_light_position);
		float diffuse_b = glm::dot(normal, b - norm_light_position);
		float diffuse_c = glm::dot(normal, c - norm_light_position);

		colors[i + 0] = color * diffuse_a;
		colors[i + 1] = color * diffuse_b;
		colors[i + 2] = color * diffuse_c;
	}
}

ColorMatchingRenderer::
ColorMatchingRenderer()
{
	update_spectrum(tesselation, temperature);

	glm::vec3 cam_pos(3.0f, 3.f, 3.0f);
	glm::vec3 cam_center(0.f, 0.f, 0.0f);

	camera = std::make_shared<LookAroundCamera>(cam_pos, cam_center, 0.0f);
}

/*
 * Render the blackbody spectrum mode.
 */
void ColorMatchingRenderer::
_draw()
{
	ImGui::Checkbox("Show Wireframe", &showWireframe);
	if(ImGui::SliderFloat("Temperature", &temperature, 400.0f, 8000.0f)
	|  ImGui::SliderInt("Tessellation", &tesselation, 1, 128))
	{
		update_spectrum(tesselation, temperature);
	}
	// Set up the projection matrix. We will never change the camera after this
	// (rotation is done by rotating the world), so we can do this here.
	glMatrixMode(GL_PROJECTION);
	glLoadMatrixf((const float *)&projection);

	// From now on, matrix operations shall affect the model-view matrix.
	glMatrixMode(GL_MODELVIEW);

	// We generally want to use the z-buffer.
	glEnable(GL_DEPTH_TEST);

	glClearColor(0.01f, 0.01f, 0.01f, 1.f);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	
	// Enable Gamma Correction
	glEnable(GL_FRAMEBUFFER_SRGB);

	glMatrixMode(GL_PROJECTION);
	glLoadMatrixf((const float *)&projection);
	// From now on, matrix operations shall affect the model-view matrix.
	glMatrixMode(GL_MODELVIEW);

	glPushMatrix();
	glTranslatef(0.f, 0.f, -1.f); // Makes sure 2d things at z=0 are shown, and zooms out a little.

	// Spectrum.
	{
		glPushMatrix();

		// Transform so that the spectrum can be drawn in [0, 1]^2.
		glTranslatef(-0.5f, -0.15f, 0.f);
		glScalef(1.f, 0.6f, 1.f);

		draw_triangle_strip(spectrum_vertices, spectrum_colors);

		glPopMatrix();
	}

	// Gradient.
	{
		glPushMatrix();

		// Transform so that the temperature gradient can be drawn in [0, 1]^2.
		glTranslatef(-0.5f, -0.5f, 0.f);
		glScalef(1.f, 0.2f, 1.f);

		// Current temperature marker.
		wireframe_overlay([&](){
			float const rel_temp = (temperature - BLACKBODY_TMIN)
								 / (BLACKBODY_TMAX-BLACKBODY_TMIN);
			float oldLineWidth;
			glGetFloatv(GL_LINE_WIDTH, &oldLineWidth);
			glLineWidth(10.f);
			glBegin(GL_LINES);
				glColor3f(1.f, 1.f, 1.f);
				glVertex2f(rel_temp, 0.f - 0.1f);
				glVertex2f(rel_temp, 1.f + 0.1f);
			glEnd();
			glLineWidth(oldLineWidth);
		});

		// The actual gradient.
		draw_triangle_strip(blackbody_vertices, blackbody_colors);

		// Gradient wireframe.
		if (showWireframe)
		{
			wireframe_overlay([&](){
				draw_triangle_strip(blackbody_vertices, blackbody_wire_colors);
			});
		}

		glPopMatrix();
	}
	
	glPopMatrix();

	// Disable Gamma Correction
	glDisable(GL_FRAMEBUFFER_SRGB);
}
/*
 * Below, we tesselate the spectrum visualization and the blackbody gradient, and update their
 * colors.
 */
void ColorMatchingRenderer::update_spectrum(
	std::uint32_t resolution,
	float T)
{
	// Generate the spectrum vertices and color if necessary.
	std::uint32_t const spectrum_res = cmf::x.size()-1;
	if (spectrum_vertices.size() != (spectrum_res*2) + 2)
	{
		generate_strip(spectrum_res, &spectrum_vertices);
		spectrum_colors.resize(spectrum_vertices.size());
		for (std::size_t i = 0; i < spectrum_vertices.size(); ++i)
		{
			std::size_t const idx = i / 2;
			if (idx >= cmf::x.size())
			{
				continue;
			}
			spectrum_colors[i] = convert::xyz_to_rgb(glm::vec3(cmf::x[idx], cmf::y[idx], cmf::z[idx]));
		}
	}

	// Update the spectrum shape.
	{
		std::vector<float> spectrum;
		blackbody_spectrum(&spectrum, T);
		// Normalize spectrum to maximum 1 for color mapping.
		float max_val = 0;
		for (std::size_t i = 0; i < spectrum.size(); ++i)
		{
			max_val = std::max<float>(max_val, spectrum[i]);
		}

		for (std::size_t i = 0; i < spectrum_vertices.size(); ++i)
		{
			if (i % 2 == 1 && i/2 < spectrum.size())
			{
				spectrum_vertices[i - 1][1] = spectrum[i/2] / max_val;
			}
		}
	}

	// Finally, regenerate the gradient strip.
	if (blackbody_vertices.size() != (resolution*2) + 2)
	{
		generate_strip(resolution, &blackbody_vertices);
		blackbody_colors.resize(blackbody_vertices.size());

		// Update the gradient colors.
		std::vector<float> spectrum;
		for (std::uint32_t i = 0; i < blackbody_vertices.size(); ++i)
		{
			auto& p = blackbody_vertices.at(i);

			float const temp = BLACKBODY_TMIN + p[0] * (BLACKBODY_TMAX-BLACKBODY_TMIN);
			blackbody_spectrum(&spectrum, temp);

			glm::vec3 const rgb = spectrum_to_rgb(spectrum);

			// Color map so we can see all of the blackbody_colors.
			glm::vec3 xyY = convert::xyz_to_xyy(convert::rgb_to_xyz(rgb));
			xyY[2] = colormap_luminance(xyY[2]);
			blackbody_colors.at(i) = convert::xyz_to_rgb(convert::xyy_to_xyz(xyY));
		}

		// We show the wireframe with inverted colors for better visibility.
		blackbody_wire_colors.resize(blackbody_colors.size());
		for (std::size_t i = 0; i < blackbody_colors.size(); ++i)
		{
			blackbody_wire_colors[i] = glm::vec3(1.f, 1.f, 1.f) - blackbody_colors[i];
		}
	}
}

GravityFieldRenderer::GravityFieldRenderer()
{
	update_gravity_field(128, mass0Center, mass0, mass1Center, mass1);

	glm::vec3 d = normalize(glm::vec3(1.f, 1.5f, 1.f));
	glm::vec3 cam_center(0.5f, -0.5f, 0.5f);
	glm::vec3 cam_pos = cam_center + 1.5f * d;

	camera = std::make_shared<LookAroundCamera>(cam_pos, cam_center, 0.0f);
}

void GravityFieldRenderer::
_draw()
{
	ImGui::Checkbox("Show Wireframe", &showWireframe);
	if(ImGui::SliderFloat("Mass 1", &mass0, 0.0f, 1.0f)
	|  ImGui::SliderFloat("Mass 2", &mass1, 0.0f, 1.0f)
	|  ImGui::SliderInt("Tessellation", &tesselation, 1, 128))
	{
		update_gravity_field(tesselation, mass0Center, mass0, mass1Center, mass1);
	}

	// Set up the projection matrix. We will never change the camera after this
	// (rotation is done by rotating the world), so we can do this here.
	glMatrixMode(GL_PROJECTION);
	glLoadMatrixf((const float *)&projection);

	// From now on, matrix operations shall affect the model-view matrix.
	glMatrixMode(GL_MODELVIEW);
	
	// Enable Gamma Correction
	glEnable(GL_FRAMEBUFFER_SRGB);

	// We generally want to use the z-buffer.
	glEnable(GL_DEPTH_TEST);

	glClearColor(0.01f, 0.01f, 0.01f, 1.f);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	// Store current transformation to avoid side effects.
	glPushMatrix();

	// Rotate field according to user input on the GUI.
	glMultMatrixf(glm::value_ptr(camera->get_view_matrix(Camera::Mono)));

	// Draw the gravity field as an indexed grid of triangles.
	draw_indexed_triangles(gravity_field_vertices, gravity_field_colors, gravity_field_indices);

	// Draw a black wireframe overlay.
	if (showWireframe)
	{
		wireframe_overlay([&]() {
			draw_indexed_triangles(gravity_field_vertices, gravity_wire_colors, gravity_field_indices);
		});
	}

	// Restore matrices.
	glPopMatrix();
	
	// Disable Gamma Correction
	glDisable(GL_FRAMEBUFFER_SRGB);
}

/*
 * The following two functions take care of tesselating the height field grid and
 * updating its height.
 */
void GravityFieldRenderer::update_gravity_field(
	std::uint32_t resolution,
	glm::vec3 const& p1,
	float m1,
	glm::vec3 const& p2,
	float m2)
{
	cg_assert(resolution >= 1);

	// Only regenerate the grid if necessary.
	std::uint32_t const num_vertices = (resolution + 1) * (resolution + 1);
	if (num_vertices != gravity_field_vertices.size())
	{
		generate_grid(resolution, &gravity_field_vertices, &gravity_field_indices);

        // bring grid from xy-plane to xz-plane
        for (auto &v : gravity_field_vertices) {
            v = glm::vec3(glm::rotate(90.f, glm::vec3(1.f, 0.f, 0.f)) * glm::vec4(v, 1.f));
        }

		gravity_wire_colors = std::vector<glm::vec3>(gravity_field_vertices.size(), glm::vec3(0.f));
		gravity_field_colors.resize(gravity_field_vertices.size());
	}

	// If there still is no grid, exit now. This will be the case when nothing has
	// been implemented.
	if (num_vertices != gravity_field_vertices.size())
	{
		return;
	}

	// Update heightfield and colors.

	// Since the potential is proportional to 1/r, it is unbounded (infinite at the point masses).
	// This leads to rendering artifacts, and so we clamp it to minimum -2.
	float const min_potential = -2;

	for (std::uint32_t y = 0; y < resolution+1; ++y)
	{
		for (std::uint32_t x = 0; x < resolution+1; ++x)
		{
			auto const index = y * (resolution+1) + x;
			auto& vertex     = gravity_field_vertices.at(index);
			auto& color      = gravity_field_colors.at(index);

			glm::vec2 const delta1(vertex[0]-p1[0], vertex[2]-p1[2]);
			float dist1 = glm::length(delta1);
			glm::vec2 const delta2(vertex[0]-p2[0], vertex[2]-p2[2]);
			float dist2 = glm::length(delta2);

			// Sum up potentials of both point masses.
			vertex.y = std::max(min_potential, gravitational_potential(m1, dist1))
				     + std::max(min_potential, gravitational_potential(m2, dist2));

			// Below is just color mapping. We color p1 blue and p2 red, and
			// generate a gradient of hues between them.
			float hue;
			float const hue1   = 0.6f;
			float const hue2   = 1.0f;
			float const weight = std::min(dist1, dist2) / std::max(dist1, dist2);
			if (dist1 > dist2) { hue = (1.f-weight) * hue1 + weight * 0.5f * (hue1+hue2); }
			else               { hue = (1.f-weight) * hue2 + weight * 0.5f * (hue1+hue2); }

			// Note how the brightness increases as we approach one of the point masses.
			color = convert::hsv_to_rgb(glm::vec3(hue, 1, vertex.y / min_potential));
		}
	}
}
// CG_REVISION 8c58412a25ac2367c053bfa840ee81b320bdd315
