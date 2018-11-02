cmake_minimum_required(VERSION 3.12)
project(hermite)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_COMPILER /usr/bin/clang++)

set(CS_INCLUDE $ENV{HOME}/library/cs)
set(GTEST_INCLUDE $ENV{HOME}/opt/googletest-1.8.0/googletest/include)
set(GTEST_LINK $ENV{HOME}/opt/googletest-1.8.0/build/googlemock/gtest)

add_subdirectory(demo)
add_subdirectory(test)
enable_testing()
add_test(NAME CPP_TEST_SOLVER4th
        COMMAND gtest_solver4th)