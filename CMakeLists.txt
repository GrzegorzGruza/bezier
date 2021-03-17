cmake_minimum_required(VERSION 3.17)
project(jnp_task7)

set(CMAKE_CXX_STANDARD 14)

include_directories(.)

add_executable(bezier_not_my
        bezier_.h
        bezier_example_.cc)

add_executable(bezier
        bezier.h
        bezier_example.cc)
