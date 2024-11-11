#include <cglib/colors/exercise.h>
#include <cglib/colors/convert.h>
#include <cglib/colors/cmf.h>

#include <cglib/core/glheaders.h>
#include <cglib/core/glmstream.h>

#include <cglib/core/assert.h>
#include <iostream>

/*
 * Draw the given vertices directly as GL_TRIANGLES.
 * For each vertex, also set the corresponding color.
 */
void draw_triangles(
	std::vector<glm::vec3> const& vertices,
	std::vector<glm::vec3> const& colors)
{
	cg_assert(vertices.size() == colors.size());
	cg_assert(vertices.size() % 3 == 0);

	for (int i =0 ; i < static_cast<int>(vertices.size())/3; ++i)
	{
		glBegin(GL_TRIANGLES);
		glColor3fv(&colors[3*i][0]);
		glVertex3fv(&vertices[3*i][0]);

		glColor3fv(&colors[3*i+1][0]);
		glVertex3fv(&vertices[3*i+1][0]);

		glColor3fv(&colors[3*i+2][0]);
		glVertex3fv(&vertices[3*i+2][0]);

		glEnd();
	}

}

/*
 * Generate a regular grid of resolution N*N (2*N*N triangles) in the xy-plane (z=0).
 * Store the grid in vertex-index form.
 *
 * The vertices of the triangles should be in counter clock-wise order.
 *
 * The grid must fill exactly the square [0, 1]x[0, 1], and must
 * be generated as an Indexed Face Set (Shared Vertex representation).
 * 
 * The first vertex should be at position (0,0,0) and the last
 * vertex at position (1,1,0)
 *
 * An example for N = 3:
 *
 *   ^
 *   |  ----------
 *   |  |\ |\ |\ |
 *   |  | \| \| \|
 *   |  ----------
 *   |  |\ |\ |\ |
 * y |  | \| \| \|
 *   |  ----------
 *   |  |\ |\ |\ |
 *   |  | \| \| \|
 *   |  ----------
 *   |
 *   |-------------->
 *          x
 *
 */
void generate_grid(
	std::uint32_t N,
	std::vector<glm::vec3>* vertices,
	std::vector<glm::uvec3>* indices)
{
	cg_assert(N >= 1);
	cg_assert(vertices);
	cg_assert(indices);

	vertices->clear();
	indices->clear();
	// Code here

	float subPoint = 1.0f/ static_cast<float>(N);
	for (uint32_t i =0; i <= N; ++i) {
		for (uint32_t j =0; j <= N; ++j) {
			glm::vec3 pointInSpace = glm::vec3(j*subPoint,i*subPoint,0.0);
			vertices->push_back(pointInSpace);
		}
	}
	for (uint32_t i =0; i < N; ++i) {
		for (uint32_t j =0; j < N; ++j) {
			uint32_t n = N+1;
			glm::uvec3 upsideTriangle = glm::uvec3(i*n+j,i*n+j+1,(i+1)*n+j);
			glm::uvec3 downsideTriangle = glm::uvec3(i*n+j+1,(i+1)*n+j+1,(i+1)*n+j);
			indices->push_back(upsideTriangle);
			indices->push_back(downsideTriangle);
		}
	}
}

/*
 * Draw the given vertices as indexed GL_TRIANGLES using glDrawElements.
 * For each vertex, also set the corresponding color.
 *
 * Don't forget to enable the correct client states. After drawing
 * the triangles, you need to disable the client states again.
 */
void draw_indexed_triangles(
	std::vector<glm::vec3>  const& vertices,
	std::vector<glm::vec3>  const& colors,
	std::vector<glm::uvec3> const& indices)
{
	cg_assert(vertices.size() == colors.size());

	glEnableClientState(GL_VERTEX_ARRAY);
	glEnableClientState(GL_COLOR_ARRAY);

	glVertexPointer(3, GL_FLOAT, 0, vertices.data());
	glColorPointer(3, GL_FLOAT, 0, colors.data());
	glDrawElements(GL_TRIANGLES, indices.size()*3, GL_UNSIGNED_INT, indices.data());

	glDisableClientState(GL_VERTEX_ARRAY);
	glDisableClientState(GL_COLOR_ARRAY);
}

/*
 * Generate a triangle strip with N segments (2*N triangles)
 * in the xy plane (z=0).
 *
 * The vertices of the triangles should be in counter clock-wise order.
 *
 * The triangle strip must fill exactly the square [0, 1]x[0, 1].
 *
 * The first vertex should be at position (0,1,0) and the last
 * vertex at position (1,0,0)
 *
 * An example for N = 3:
 *
 *   ^  (0,1,0)
 *   |  ----------
 *   |  | /| /| /|
 * y |  |/ |/ |/ |
 *   |  ---------- (1,0,0)
 *   |
 *   |-------------->
 *           x
 *
 */
void generate_strip(
	std::uint32_t N,
	std::vector<glm::vec3>* vertices)
{
	cg_assert(N >= 1);
	cg_assert(vertices);

	vertices->clear();
	//N=1;
	float subPointDistance = 1.0f / static_cast<float>(N);
	vertices->push_back(glm::vec3(0,1,0));
	vertices->push_back(glm::vec3(0,0,0));
	vertices->push_back(glm::vec3(subPointDistance,1,0));
	for (uint32_t i =0; i < N-1; ++i) {
		vertices->push_back(glm::vec3(static_cast<float>(i+1)*subPointDistance,0,0));
		vertices->push_back(glm::vec3(static_cast<float>(i+2)*subPointDistance,1,0));
	}
	vertices->push_back(glm::vec3(1,0,0));
}

/*
 * Draw the given vertices as a triangle strip.
 * For each vertex, also set the corresponding color.
 */
void draw_triangle_strip(
	std::vector<glm::vec3> const& vertices,
	std::vector<glm::vec3> const& colors)
{
	cg_assert(vertices.size() == colors.size());

	glBegin(GL_TRIANGLE_STRIP);
	for (int i =0; i < static_cast<int>(vertices.size()); ++i)
	{
		glVertex3f(vertices[i].x,vertices[i].y,vertices[i].z);
		glColor3f(colors[i].x,colors[i].y,colors[i].z);
	}
	glEnd();
}

/*
 * Integrate the given piecewise linear function 
 * using trapezoidal integration.
 *
 * The function is given at points
 *     x[0], ..., x[N]
 * and its corresponding values are
 *     y[0], ..., y[N]
 */
float integrate_trapezoidal_student(
	std::vector<float> const& x,
	std::vector<float> const& y)
{
	cg_assert(x.size() == y.size());
	cg_assert(x.size() > 1);

	float sumOfAll = 0.0;

	for (uint32_t i =0; i+1 < x.size(); ++i)
	{
		sumOfAll += y[i] * y[i+1] /2.0*(x[i+1] - x[i]);
	}
	return sumOfAll;
}

/*
 * Convert the given spectrum to RGB using your
 * implementation of integrate_trapezoidal(...)
 *
 * The color matching functions and the wavelengths
 * for which they are given can be found in
 *     cglib/colors/cmf.h
 * and
 *     cglib/src/colors/cmf.cpp
 *
 * The wavelengths corresponding to the spectral values 
 * given in spectrum are defined in cmf::wavelengths
 */
glm::vec3 spectrum_to_rgb(std::vector<float> const& spectrum)
{
	cg_assert(spectrum.size() == cmf::wavelengths.size());
	float x = integrate_trapezoidal_student(cmf::x,spectrum);
	float y = integrate_trapezoidal_student(cmf::y,spectrum);
	float z = integrate_trapezoidal_student(cmf::z,spectrum);

	return convert::xyy_to_rgb(glm::vec3(x,y,z));
}
// CG_REVISION 8c58412a25ac2367c053bfa840ee81b320bdd315

