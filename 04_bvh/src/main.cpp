#include <cglib/rt/bvh.h>
#include <cglib/rt/host_render.h>
#include <cglib/rt/intersection_tests.h>
#include <cglib/rt/ray.h>
#include <cglib/rt/raytracing_context.h>
#include <cglib/rt/render_data.h>
#include <cglib/rt/renderer.h>
#include <cglib/rt/sampling_patterns.h>
#include <cglib/rt/scene.h>
#include <cglib/rt/raytracing_parameters.h>

#include <cglib/core/glmstream.h>
#include <cglib/core/image.h>
#include <cglib/core/parameters.h>
#include <cglib/core/thread_local_data.h>
#include <cglib/core/thread_local_data.h>
#include <cglib/core/timer.h>

#include <glm/gtc/matrix_transform.hpp>

#include <cglib/core/assert.h>
#include <cmath>
#include <sstream>
#include <chrono>
#include <iostream>
#include <stack>
#include <complex>

using std::cout;
using std::cerr;
using std::endl;


/*
 * This is the main rendering kernel.
 *
 * It is supposed to compute the RGB color of the given pixel (x,y).
 *
 * RenderData contains data relevant for the computation of the color
 * for one pixel. Thread-local data is referenced by this struct, aswell. The
 * pointer tld is guaranteed to be valid (not nullptr).
 */
glm::vec3 render_pixel(int x, int y, RaytracingContext const& context, RenderData &data)
{
	cg_assert(data.tld);

	if(dynamic_cast<GaussScene *>(context.get_active_scene())) {
		float u = float(x)/context.params.image_width;
		float v = float(y)/context.params.image_height;
		const auto &scene = context.get_active_scene();
		const char *tex_names[RaytracingParameters::GAUSS_MODE_COUNT] = { "input", "filtered_naive", "filtered_separable" };
		const auto &texture = scene->textures.find(tex_names[context.params.gauss_mode])->second;
		auto &mip = texture->get_mip_levels()[0];
		float a = float(mip->getWidth()) / float(mip->getHeight());
		return texture->evaluate_bilinear(0, a > 1.0
				? glm::vec2(u, v * a - (a - 1.0f) * 0.5f)
				: glm::vec2(u / a - (1.0f - a) * 0.5f, v));
	}

	std::vector<glm::vec2> samples;
	int spp = data.context.params.spp;

	if(spp > 1) {
		int grid_size = int(sqrtf(static_cast<float>(spp)));
		if(data.context.params.stratified)
			generate_stratified_samples(&samples, grid_size, grid_size, data.tld);
		else
			generate_random_samples(&samples, grid_size, grid_size, data.tld);
		glm::vec3 accum(0.0f);

		for(size_t i = 0; i < samples.size(); i++) {
			float fx = float(x) + samples[i].x;
			float fy = float(y) + samples[i].y;

			data.x = fx;
			data.y = fy;

			Ray ray = createPrimaryRay(data, fx, fy);
			accum += trace_recursive(data, ray, 0/*depth*/);
		}

		return accum / float(samples.size());
	}
	else {
		float fx = float(x) + 0.5f;
		float fy = float(y) + 0.5f;

		data.x = fx;
		data.y = fy;

		Ray ray = createPrimaryRay(data, fx, fy);
		return trace_recursive(data, ray, 0/*depth*/);
	}
}

static const std::string image_prefix = "assignment_images/";

static void render_triangles(std::string const& output_name, int num_triangles)
{
    RaytracingContext context;
	context.params.interactive = 0;
	context.params.image_width  = 512;
	context.params.image_height = 512;
	context.params.num_triangles = num_triangles;
	context.params.output_file_name = output_name;

	context.params.output_file_name = image_prefix + output_name;
	context.add_scene(std::make_shared<TriangleScene>(context.params));
	HostRender::run(context, render_pixel);
}

static void render_monkey(std::string const& output_name)
{
    RaytracingContext context;
	context.params.interactive = 0;
	context.params.image_width  = 512;
	context.params.image_height = 512;
	context.params.render_mode = RaytracingParameters::NORMAL;

	context.params.output_file_name = image_prefix + output_name;
	context.add_scene(std::make_shared<MonkeyScene>(context.params));
	HostRender::run(context, render_pixel);
}

static void render_sponza(std::string const& output_name)
{
    RaytracingContext context;
	context.params.interactive = 0;
	context.params.image_width  = 128;
	context.params.image_height = 128;
	context.params.render_mode = RaytracingParameters::NORMAL;

	context.params.output_file_name = image_prefix + output_name;
	context.add_scene(std::make_shared<SponzaScene>(context.params));
	HostRender::run(context, render_pixel);
}

void gauss_filter()
{
	// load the image
	Image img;
	img.load("assets/lion.png", 1.f);

	int radius = 9;
	int kernel_size = radius * 2 + 1;
	float sigma = radius * 0.5;

	std::cout << "Filter naive..." << std::endl;
	Image filtered_naive(img.getWidth(), img.getHeight());
	img.filter_gaussian(&filtered_naive, sigma, kernel_size);

	std::cout << "Filter separable..." << std::endl;
	Image filtered_seperable(img.getWidth(), img.getHeight());
	img.filter_gaussian_separable(&filtered_seperable, sigma, kernel_size);

	std::cout << "Save results..." << std::endl;
	filtered_naive.save(image_prefix+"gauss_filtered_naive.png", 1.f);
	filtered_seperable.save(image_prefix+"gauss_filtered_seperable.png", 1.f);
}

void create_images()
{
	render_triangles("triangle.png", 1);
	render_triangles("triangles.png", 10);
	render_monkey("monkey.png");
	render_sponza("sponza.png");
}

int
main(int argc, char const**argv)
{
    RaytracingContext context;
    if (!context.params.parse_command_line(argc, argv)) {
        std::cerr << "invalid command line argument" << std::endl;
        return -1;
    }
	
	if(context.params.create_images) {
		create_images();
		return 0;
	}
	if(context.params.gauss) {
		gauss_filter();
		return 0;
	}

	context.add_scene(std::make_shared<TriangleScene>(context.params));
	context.add_scene(std::make_shared<MonkeyScene>(context.params));
	context.add_scene(std::make_shared<SponzaScene>(context.params));
	context.add_scene(std::make_shared<GaussScene>(context.params));

	return HostRender::run(context, render_pixel);
}

// CG_REVISION 8c58412a25ac2367c053bfa840ee81b320bdd315
