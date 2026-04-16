cmake_minimum_required(VERSION 3.5)
project(zzz C)

set(CMAKE_C_STANDARD 11)

set(CMAKE_C_FLAGS "-O3 -funroll-loops")

set(CMAKE_LIB)

option(STATIC_LINKING "static linking dependencies" OFF)

IF(STATIC_LINKING)
    include_directories(${PROJECT_SOURCE_DIR}/deps/)

    link_directories(
        ${PROJECT_SOURCE_DIR}/deps/gmp-6.2.1/.libs
        ${PROJECT_SOURCE_DIR}/deps/flint-2.9.0
        ${PROJECT_SOURCE_DIR}/deps/mpfr-4.2.0/src/.libs
        ${PROJECT_SOURCE_DIR}/deps/arb
    )

    link_libraries(
        m
        libarb.a
        libflint.a
        libmpfr.a
        libgmp.a
    )
ELSE()
    find_library(GMP gmp)
    find_library(FLINT flint)
    find_library(MPFR mpfr)
    link_libraries(gmp flint mpfr)
ENDIF()

add_library(ghy STATIC ghy.c)

add_executable(zzz main.c)
add_executable(zproxy zproxy.c)
add_executable(zghy zghy.c)
add_executable(zhad zhad.c)
target_link_libraries(zghy ghy m)
target_link_libraries(zhad ghy m)
