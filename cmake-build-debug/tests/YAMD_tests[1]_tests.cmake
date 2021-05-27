add_test( VerletTest.BasicAssertions /home/trirop/yamd/cmake-build-debug/tests/YAMD_tests [==[--gtest_filter=VerletTest.BasicAssertions]==] --gtest_also_run_disabled_tests)
set_tests_properties( VerletTest.BasicAssertions PROPERTIES WORKING_DIRECTORY /home/trirop/yamd/cmake-build-debug/tests)
set( YAMD_tests_TESTS VerletTest.BasicAssertions)
