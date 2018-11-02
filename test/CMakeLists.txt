set(CMAKE_CXX_FLAGS "-Wall -O0 -Wextra")

include_directories(${GTEST_INCLUDE} ${CS_INCLUDE})
link_directories(${GTEST_LINK})

add_executable(gtest_solver4th gtest_solver4th.cc)
target_link_libraries(gtest_solver4th pthread gtest)