#******************************************************************************
# Copyright (C) Siarhei Uzunbajakau, 2023.
#
# This program is free software. You can use, modify, and redistribute it under
# the terms of the GNU Lesser General Public License as published by the Free
# Software Foundation, either version 3 or (at your option) any later version.
# This program is distributed without any warranty.
#
# Refer to COPYING.LESSER for more details.
#*****************************************************************************/

set(PRJ_NAME "mms-vt-ii")

cmake_minimum_required(VERSION 2.8.12)

include_directories(${DEAL_II_DIR} /usr/include/suitesparse ./include
../shared/include)

set(SOURCE_FILES
	"include/solver.hpp"
	"include/exact_solution.hpp"
	"include/settings.hpp"
	"../shared/include/static_vector_solver_i.hpp"
	"../shared/include/static_vector_solver_ii.hpp"
	"../shared/include/static_vector_input.hpp"
	"../shared/include/project_Hgrad_to_Hdiv.hpp"
	"../shared/include/project_Hcurl_to_L2.hpp"
	"../shared/include/project_T_to_J.hpp"
	"../shared/include/project_A_to_B.hpp"
	"../shared/include/misc.hpp"
	"../shared/include/constants.hpp"
	"src/main.cpp"
	"src/solver.cpp"
	"src/exact_solution.cpp"
	"src/static_scalar_input.cpp"
	"src/static_vector_input.cpp"
	"../shared/src/misc.cpp"
	"../shared/src/constants.cpp"
	)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

find_package(deal.II 9.3.2 QUIET
	HINTS ${deal.II_DIR} ${DEAL_II_DIR} ../ ../../ $ENV{DEAL_II_DIR}
)
if(NOT ${deal.II_FOUND})
message(FATAL_ERROR "\n"
	"*** Could not locate a (sufficiently recent) version of deal.II. ***\n\n"
	"You may want to either pass a flag -DDEAL_II_DIR=/path/to/deal.II to cmake\n"
	"or set an environment variable \"DEAL_II_DIR\" that contains this path."
)
endif()
DEAL_II_INITIALIZE_CACHED_VARIABLES()

project(${PRJ_NAME})

set(TARGET_CIRCLE "${PRJ_NAME}-circle")
set(TARGET_SQUARE "${PRJ_NAME}-square")

foreach(TARGET IN LISTS TARGET_CIRCLE TARGET_SQUARE)

add_executable(${TARGET} ${SOURCE_FILES})
DEAL_II_SETUP_TARGET(${TARGET})

set_target_properties(${TARGET}
	PROPERTIES RUNTIME_OUTPUT_DIRECTORY
	"${PROJECT_SOURCE_DIR}/bin/$<CONFIG>")

message(STATUS "TARGET=${TARGET}")

endforeach(TARGET)

target_compile_options(${TARGET_CIRCLE} PRIVATE -DDOMAIN__=0)
target_compile_options(${TARGET_SQUARE} PRIVATE -DDOMAIN__=1)

