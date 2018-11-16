cmake_minimum_required(VERSION 3.12)
project(stat906A3)

set(CMAKE_CXX_STANDARD 14)
include_directories(./include)
find_package(PythonLibs 2.7)
add_executable(stat906A3 main.cpp include/q2.h include/q3.h include/global.h)

target_include_directories(stat906A3 PRIVATE /usr/include/python2.7)
target_link_libraries(stat906A3 python2.7)
