#include <cglib/core/glheaders.h>
#include <cglib/core/assert.h>

#include <cglib/gl/device_rendering_context.h>
#include <cglib/gl/device_render.h>
#include <cglib/gl/renderer.h>

#include <sstream>
#include <fstream>
#include <string>
#include <iostream>
# include <type_traits>

#include "renderers.h"

enum Renderers { CUBE, MONKEY, SPHERE_FLAKE, SPONZA, NUM_RENDERERS };


enum {
	FLAGS_DEFAULT = 0,
	IMG_COMPARE_NO_SPONZA_K_S = (1 << 0)
};

std::string assignment_images_path(std::string const& filename)
{
//	ensure_directory(dir);
	std::ostringstream os;
	os << "assignment_images";
	os << "/";
	os << filename;
	return os.str();
}

void create_image(int renderer, int flags)
{
	DeviceRenderingContext context;
	context.params.interactive = false;
	//context.params.enable_animation = false;
	std::string img_name;
	switch(renderer) {
	case CUBE: {
			auto r = std::make_shared<CubeRenderer>();
			context.renderers.push_back(r);
			img_name = "cube.png";
			break;
		}
	case MONKEY: {
			auto r = std::make_shared<MonkeyRenderer>();
			context.renderers.push_back(r);
			img_name = "monkey.png";
			break;
		}
	case SPHERE_FLAKE: {
			auto r = std::make_shared<SphereFlakeRenderer>();
			r->enable_animation = false;
			r->set_interactive(false);//this disables drawing imgui
			context.renderers.push_back(r);
			img_name = "sphere_flake.png";
			break;
		}
	case SPONZA: {
			auto r = std::make_shared<SponzaRenderer>();
			r->enable_animation = false;
			context.renderers.push_back(r);
			img_name = "sponza.png";
			r->set_interactive(false);
			context.renderers.push_back(r);
			break;
		}
	}

	context.params.output_file_name = assignment_images_path(img_name);

    DeviceRender::run(context, GUI::DEFAULT_FLAGS);
}

void create_images()
{
	std::cout << "Render Monkey" << std::endl;
	create_image(MONKEY, FLAGS_DEFAULT);
	std::cout << "Render Cube" << std::endl;
	create_image(CUBE, FLAGS_DEFAULT);
	std::cout << "Render Sphere Flake" << std::endl;
	create_image(SPHERE_FLAKE, FLAGS_DEFAULT);
	std::cout << "Render Sponza" << std::endl;
	create_image(SPONZA, FLAGS_DEFAULT);
}

int
main(int argc, char const**argv)
{
	DeviceRenderingContext *context = new DeviceRenderingContext;
    if (!context->params.parse_command_line(argc, argv)) {
        std::cerr << "invalid command line argument" << std::endl;
        return -1;
    }

    if (context->params.create_images)
    {
	    delete context;
	    create_images();
	    return 0;
    }

	context->renderers.push_back(std::make_shared<CubeRenderer>());
	context->renderers.push_back(std::make_shared<MonkeyRenderer>());
	context->renderers.push_back(std::make_shared<SphereFlakeRenderer>());
	context->renderers.push_back(std::make_shared<SponzaRenderer>());

    return DeviceRender::run(*context, GUI::DEFAULT_FLAGS);
}

// CG_REVISION d4ab32bd208749f2d2b1439e25d16e642b039298
