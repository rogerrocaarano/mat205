cmake_minimum_required(VERSION 3.24)
project(mat205)

set(CMAKE_CXX_STANDARD 23)

add_executable(mat205 main.cpp lib/matrix.h lib/matrix.cpp)
