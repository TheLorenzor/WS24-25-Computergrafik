#include <cglib/core/assert.h>
#include <cglib/core/camera.h>
#include <cglib/gl/shader.h>

#include <glm/glm.hpp>
#include <cglib/gl/renderer.h>

/*
 * Call the drawing functor, but first set up the OpenGL context so that
 * all geometry is drawn as a wireframe overlay.
 */
template <typename DrawFunc>
void wireframe_overlay(DrawFunc&& draw)
{
	// Store a copy of the current OpenGL state for all attributes that
	// we are going to modify.
	glPushAttrib(GL_DEPTH_BUFFER_BIT | GL_POLYGON_BIT);

	// Disable depth testing and writing (this is an overlay), and draw as lines.
	glDisable(GL_DEPTH_TEST);
	glDepthMask(GL_FALSE);
	glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

	// Slightly thicker lines for
	//glLineWidth(2.f);

	// Draw using the provided functor.
	draw();

	// Restore OpenGL state.
	glPopAttrib();
}

class RGBCubeRenderer : public Renderer
{
public:
	RGBCubeRenderer();

	virtual void _draw() override;
	const char * get_name() override { return "RGB Cube"; }

	std::vector<glm::vec3> cube_vertices;
	std::vector<glm::vec3> cube_colors;
	std::vector<glm::vec3> cube_wire_colors;

	bool showWireframe = false;
};

class CGFontRenderer : public Renderer
{
public:
	CGFontRenderer();

	const char * get_name() override { return "CG Font"; }

	virtual void _draw() override;

	std::vector<glm::vec3> font_c_vertices;
	std::vector<glm::vec3> font_c_vertices_colors;
	std::vector<glm::vec3> font_c_wire_colors;

	std::vector<glm::vec3> font_d_vertices;
	std::vector<glm::vec3> font_d_vertices_colors;
	std::vector<glm::vec3> font_d_wire_colors;

	bool showWireframe = false;

private:
	void update_vertex_light_color(const std::vector<glm::vec3>& vertices, std::vector<glm::vec3>& colors, const glm::vec3& light_position, const glm::vec3& color) const;
};

class GravityFieldRenderer : public Renderer
{
public:
	GravityFieldRenderer();

	const char * get_name() override { return "Gravity Field"; }

	virtual void _draw() override;

	// Gravitational potential height field.
	std::vector<glm::vec3>  gravity_field_vertices;
	std::vector<glm::vec3>  gravity_field_colors;
	std::vector<glm::vec3>  gravity_wire_colors;
	std::vector<glm::uvec3> gravity_field_indices;

	bool showWireframe = false;
	glm::vec3 mass0Center = glm::vec3(0.25f, 0.f, 0.25f);
	glm::vec3 mass1Center = glm::vec3(0.75f, 0.f, 0.75f);
	float mass0 = 0.050f;
	float mass1 = 0.020f;
	int tesselation = 16;

	void update_gravity_field(
		// Grid resolution.
		std::uint32_t resolution,
		// Point mass 1.
		glm::vec3 const& p1,
		float m1,
		// Point mass 2.
		glm::vec3 const& p2,
		float m2);
};

class ColorMatchingRenderer : public Renderer
{
public:
	ColorMatchingRenderer();

	const char * get_name() override { return "Color Matching"; }

	virtual void _draw() override;

	// Triangle strip used for the spectrum visualization.
	std::vector<glm::vec3> spectrum_vertices;
	std::vector<glm::vec3> spectrum_colors;
	std::vector<glm::vec3> blackbody_vertices;
	std::vector<glm::vec3> blackbody_colors;
	std::vector<glm::vec3> blackbody_wire_colors;

	bool showWireframe = false;
	float temperature = 6500.f;
	int tesselation = 16;

	void update_spectrum(
		// Gradient strip resolution.
		std::uint32_t resolution,
		// Current temperature.
		float T);
};
// CG_REVISION 8c58412a25ac2367c053bfa840ee81b320bdd315
