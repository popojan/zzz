cmake_minimum_required(VERSION 3.24)
project(zzz C)

set(CMAKE_C_STANDARD 11)

set(CMAKE_C_FLAGS "-O3 -funroll-loops")

find_library(ARB arb)

include_directories(/usr/include/arb/)

link_libraries(arb flint gmp)

add_executable(zzz main.c)
