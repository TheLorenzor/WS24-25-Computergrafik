cmake_minimum_required (VERSION 2.8) 
project (cg) 

set(EXERCISE_NUMBER 6)

set(EXERCISE_SOURCE_SHADERS
		shader/simple.vert
		shader/simple.frag
		shader/env_light.frag
		shader/env_light.vert
		shader/phong.frag
		shader/phong.vert
)

set(EXERCISE_SOURCE
	src/main.cpp
	src/exercise_06.cpp
	src/renderers.cpp
)

if (NOT DEFINED CGLIB_DIR)
	set(CGLIB_DIR ${CMAKE_SOURCE_DIR}/../cglib)
endif()
get_filename_component(CGLIB_DIR ${CGLIB_DIR} ABSOLUTE)

include(${CGLIB_DIR}/CMakeListsCommon.cmake)
