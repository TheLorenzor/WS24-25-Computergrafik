#include <cglib/core/glheaders.h>
#include <cglib/core/assert.h>

#include <cglib/gl/device_rendering_context.h>
#include <cglib/gl/device_render.h>
#include <cglib/gl/renderer.h>
#include "renderers.h"

int
main(int argc, char const**argv)
{
	DeviceRenderingContext context;
    if (!context.params.parse_command_line(argc, argv)) {
        std::cerr << "invalid command line argument" << std::endl;
        return -1;
    }

	context.renderers.push_back(std::make_shared<ShadowmapRenderer>());
	context.renderers.push_back(std::make_shared<FirePlaceRenderer>());
	context.renderers.push_back(std::make_shared<ProceduralLandscapeRenderer>());

    return DeviceRender::run(context, GUI::DEFAULT_FLAGS);
}
// CG_REVISION d4ab32bd208749f2d2b1439e25d16e642b039298
