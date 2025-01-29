cmake_minimum_required (VERSION 2.8) 
project (cg) 

set(EXERCISE_NUMBER 7)

set(EXERCISE_SOURCE_SHADERS
	shader/fire.frag
	shader/fixed_color.frag
	shader/fixed_color.vert
	shader/heightmap.frag
	shader/heightmap.vert
	shader/logs.frag
	shader/logs.vert
	shader/particle.geom
	shader/particle.vert
	shader/shadowmap.frag
	shader/shadowmap.vert
	shader/shuttle.frag
	shader/shuttle.vert
	shader/sky_bg.frag
	shader/sky_bg.vert
	shader/smoke.frag
	shader/smoke.vert
	shader/terrain.frag
	shader/vertex_displacement.vert
	shader/visualize_depth.frag
	shader/visualize_depth.vert
	shader/visualize_texture.frag
	shader/visualize_texture.vert
	shader/water.frag
	shader/wavemap.frag
	shader/wavemap.vert
)

set(EXERCISE_SOURCE
	src/exercise_07.cpp
	src/main.cpp
	src/renderers.cpp
	src/renderers.h
)

if (NOT DEFINED CGLIB_DIR)
	set(CGLIB_DIR ${CMAKE_SOURCE_DIR}/../cglib)
endif()
get_filename_component(CGLIB_DIR ${CGLIB_DIR} ABSOLUTE)

include(${CGLIB_DIR}/CMakeListsCommon.cmake)

