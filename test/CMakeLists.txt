find_package(GTest REQUIRED)
set(CMAKE_CXX_FLAGS "-Wall -O0 -Wextra")
link_directories(${GTEST_LINK})
add_executable(gtest_solver4th gtest_solver4th.cc)
target_link_libraries(gtest_solver4th
  PRIVATE
    GTest::GTest
    GTest::Main
)
