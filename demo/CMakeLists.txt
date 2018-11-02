include_directories(${CS_INCLUDE})
set(CMAKE_CXX_FLAGS "-g -Wall -O3 -Wextra -pthread")

add_executable(twobody twobody.cc)
target_link_libraries(twobody pthread)

configure_file(./analyze_twobody.ipynb ../../build/demo/)
