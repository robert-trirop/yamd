include(FetchContent)
FetchContent_Declare(
        googletest
        URL https://github.com/google/googletest/archive/eaf9a3fd77869cf95befb87455a2e2a2e85044ff.zip
)
FetchContent_MakeAvailable(googletest)
# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

set(YAMD_tests_SRCS
        test_verlet.cpp
        test_lj_direct_summation.cpp
        test_berendsen_thermostat.cpp
        test_neighbors.cpp)

# Tests
include_directories(YAMD_tests ${CMAKE_SOURCE_DIR})
add_executable(YAMD_tests ${YAMD_SRCS} ${YAMD_HDRS} ${YAMD_tests_SRCS})
target_link_libraries(YAMD_tests PUBLIC Eigen3::Eigen gtest gtest_main)

include(GoogleTest)
gtest_discover_tests(YAMD_tests)

add_test(NAME YAMD_tests COMMAND YAMD_tests)
