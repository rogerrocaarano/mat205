cmake_minimum_required(VERSION 3.24)
project(mat205)

set(CMAKE_CXX_STANDARD 23)

add_executable(
        mat_205_ecuacionesLineales
        ecuaciones-lineales.cpp lib/matrix.h lib/matrix.cpp
)
add_executable(
        mat_205_interpolacion
        interpolacion-main.cpp
        lib/interpolacion.h lib/interpolacion.cpp
        lib/polinomio.cpp lib/polinomio.h
        lib/matriz.h lib/matriz.cpp
)
add_executable(mat_205_integracionTrapecios integracionTrapecios.cpp)
add_executable(mat_205_integracionSimpson integracionSimpson.cpp)