cmake_minimum_required(VERSION 3.11)

project(YAMD)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
set(CMAKE_CXX_STANDARD 17)

include(tools)

add_external_package(Eigen3 VERSION 3.4.0 CONFIG)

set(YAMD_HDRS
        # Add header files here
        ../HeaderFiles/verlet.h
        ../HeaderFiles/types.h
        ../HeaderFiles/atoms.h
        ../HeaderFiles/lj_direct_summation.h
        )

set(YAMD_SRCS
        # Add source files here
        ../SourceFiles/verlet.cpp
        ../SourceFiles/lj_direct_summation.cpp
        )

set(YAMD_executable_SRCS
        main.cpp
        )

add_subdirectory(tests)

# Simulation code
add_executable(YAMD ${YAMD_SRCS} ${YAMD_HDRS} ${YAMD_executable_SRCS})
target_link_libraries(YAMD PUBLIC Eigen3::Eigen)
