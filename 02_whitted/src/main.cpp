#include <cglib/rt/intersection_tests.h>
#include <cglib/rt/ray.h>
#include <cglib/rt/raytracing_context.h>
#include <cglib/rt/render_data.h>
#include <cglib/rt/renderer.h>
#include <cglib/rt/scene.h>
#include <cglib/rt/host_render.h>

#include <cglib/core/glmstream.h>
#include <cglib/core/heatmap.h>
#include <cglib/core/image.h>
#include <cglib/core/parameters.h>
#include <cglib/core/thread_local_data.h>
#include <cglib/core/thread_local_data.h>

#include <glm/gtc/matrix_transform.hpp>

#include <cglib/core/assert.h>
#include <cmath>
#include <sstream>
#include <chrono>

static const std::string image_prefix = "assignment_images/";

enum {
	RES_WIDTH = 512,
	RES_HEIGHT = 512
};

enum {
	SPHERES      = (1 << 0),
	SHADOWS      = (1 << 1),
	AMBIENT      = (1 << 2),
	DIFFUSE      = (1 << 3),
	SPECULAR     = (1 << 4),
	REFLECTION   = (1 << 5),
	TRANSMISSION = (1 << 6),
	DISPERSION   = (1 << 7),
	FRESNEL      = (1 << 8),
	PHONG        = (1 << 9),
	SPOTLIGHT    = (1 << 10),
};

enum {
	SCENE_BOX,
	SCENE_SPHERES
};

std::string get_scene_prefix(int scene)
{
	switch(scene) {
		case SCENE_BOX: return "box";
		case SCENE_SPHERES: return "spheres";
		default: cg_assert(!"invalid scene");
	}
	return "";
}

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
	float fx = x + 0.5f;
	float fy = y + 0.5f;

	data.x = fx;
	data.y = fy;

	Ray ray = createPrimaryRay(data, fx, fy);
	return trace_recursive(data, ray, 0/*depth*/);
}

static void
render_image(const std::string &output_name, int params, int scene, int max_depth = 3)
{
    RaytracingContext context;
	context.params.max_depth = max_depth;
	context.params.draw_spheres = !!(params & SPHERES);
	context.params.phong =        !!(params & PHONG);
	context.params.spot_light =   !!(params & SPOTLIGHT);
	context.params.shadows =      !!(params & SHADOWS);
	context.params.ambient =      !!(params & AMBIENT);
	context.params.diffuse =      !!(params & DIFFUSE);
	context.params.specular =     !!(params & SPECULAR);
	context.params.reflection =   !!(params & REFLECTION);
	context.params.transmission = !!(params & TRANSMISSION);
	context.params.dispersion =   !!(params & DISPERSION);
	context.params.fresnel =      !!(params & FRESNEL);
	context.params.interactive = 0;
	context.params.image_width  = RES_WIDTH;
	context.params.image_height = RES_HEIGHT;
	context.params.output_file_name = output_name;

	switch(scene) {
	case SCENE_BOX:
		context.add_scene(std::make_shared<CornellBox>(context.params));
		break;
	case SCENE_SPHERES:
		context.add_scene(std::make_shared<SpherePortrait>(context.params));
		break;
	default:
		cg_assert(!"invalid scene");
	}

	HostRender::run(context, render_pixel);
}

static void
create_images()
{
	render_image(image_prefix+get_scene_prefix(SCENE_BOX)+"_after_a.png", SPHERES | AMBIENT, SCENE_BOX, 4);
	render_image(image_prefix+get_scene_prefix(SCENE_BOX)+"_after_b_no_shadow.png", PHONG | AMBIENT | DIFFUSE | SPECULAR | SPHERES, SCENE_BOX, 4);
	render_image(image_prefix+get_scene_prefix(SCENE_BOX)+"_after_b.png", SHADOWS | PHONG | AMBIENT | DIFFUSE | SPECULAR | SPHERES, SCENE_BOX, 4);
	render_image(image_prefix+get_scene_prefix(SCENE_BOX)+"_after_d.png", REFLECTION | SHADOWS | PHONG | AMBIENT | DIFFUSE | SPECULAR | SPHERES, SCENE_BOX, 4);
	render_image(image_prefix+get_scene_prefix(SCENE_BOX)+"_after_e.png", TRANSMISSION | REFLECTION | SHADOWS | PHONG | AMBIENT | DIFFUSE | SPECULAR | SPHERES, SCENE_BOX, 4);
	render_image(image_prefix+get_scene_prefix(SCENE_BOX)+"_after_f.png", FRESNEL | TRANSMISSION | REFLECTION | SHADOWS | PHONG | AMBIENT | DIFFUSE | SPECULAR | SPHERES, SCENE_BOX, 4);
	render_image(image_prefix+get_scene_prefix(SCENE_BOX)+"_after_g.png", DISPERSION | FRESNEL | TRANSMISSION | REFLECTION | SHADOWS | PHONG | AMBIENT | DIFFUSE | SPECULAR | SPHERES, SCENE_BOX, 4);

	render_image(image_prefix+get_scene_prefix(SCENE_SPHERES)+"_after_a.png", SPHERES | AMBIENT, SCENE_SPHERES, 4);
	render_image(image_prefix+get_scene_prefix(SCENE_SPHERES)+"_after_b.png", SHADOWS | PHONG | AMBIENT | DIFFUSE | SPECULAR | SPHERES, SCENE_SPHERES, 4);
	render_image(image_prefix+get_scene_prefix(SCENE_SPHERES)+"_after_c.png", SPOTLIGHT | SHADOWS | PHONG | AMBIENT | DIFFUSE | SPECULAR | SPHERES, SCENE_SPHERES, 4);
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

    context.add_scene(std::make_shared<CornellBox>(context.params));
    context.add_scene(std::make_shared<SphereReflection>(context.params));
    context.add_scene(std::make_shared<SpherePortrait>(context.params));

	return HostRender::run(context, render_pixel);
}

// CG_REVISION 8c58412a25ac2367c053bfa840ee81b320bdd315
