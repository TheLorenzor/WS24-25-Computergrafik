cmake_minimum_required (VERSION 2.8) 
project (cg) 

set(EXERCISE_NUMBER 2)

set(EXERCISE_SOURCE_SHADERS
)

set(EXERCISE_SOURCE
	src/main.cpp
	src/exercise_02.cpp
)

if (NOT DEFINED CGLIB_DIR)
	set(CGLIB_DIR ${CMAKE_SOURCE_DIR}/../cglib)
endif()
get_filename_component(CGLIB_DIR ${CGLIB_DIR} ABSOLUTE)

include(${CGLIB_DIR}/CMakeListsCommon.cmake)
