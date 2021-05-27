include(FetchContent)
FetchContent_Declare(
        googletest
        URL https://github.com/google/googletest/archive/eaf9a3fd77869cf95befb87455a2e2a2e85044ff.zip
)
FetchContent_MakeAvailable(googletest)
# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

set(MYPROJECT_tests_SRCS
        test_verlet.cpp)

# Tests
include_directories(yamd_tests ${CMAKE_SOURCE_DIR})
add_executable(MYPROJECT_tests ${MYPROJECT_SRCS} ${MYPROJECT_HDRS} ${MYPROJECT_tests_SRCS})
target_link_libraries(MYPROJECT_tests PUBLIC Eigen3::Eigen gtest gtest_main)

include(GoogleTest)
gtest_discover_tests(MYPROJECT_tests)

add_test(NAME MYPROJECT_tests COMMAND MYPROJECT_tests)
