cmake_minimum_required(VERSION 3.5)
project(zzz C)

set(CMAKE_C_STANDARD 11)

set(CMAKE_C_FLAGS "-O3 -funroll-loops")

option(STATIC_LINKING "static linking dependencies" OFF)

# All zeta-zero sources live under src/ (psi.c / pellreg*.c in the tree root are
# unrelated and not built here).
set(ZZZ_SRC ${PROJECT_SOURCE_DIR}/src)
include_directories(${ZZZ_SRC})

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

add_library(ghy STATIC ${ZZZ_SRC}/ghy.c ${ZZZ_SRC}/weil.c ${ZZZ_SRC}/loop.c)

add_executable(zzz     ${ZZZ_SRC}/main.c)
add_executable(zproxy  ${ZZZ_SRC}/zproxy.c)
add_executable(zghy    ${ZZZ_SRC}/zghy.c)
add_executable(zhad    ${ZZZ_SRC}/zhad.c)
add_executable(zhybrid ${ZZZ_SRC}/zhybrid.c)
target_link_libraries(zzz     ghy m)
target_link_libraries(zghy    ghy m)
target_link_libraries(zhad    ghy m)
target_link_libraries(zhybrid ghy m)

# ---- experimental: Vulkan-accelerated batch method-B (opt-in) --------------
option(BUILD_VULKAN "build the experimental Vulkan batch method-B tool (zvk)" OFF)
if(BUILD_VULKAN)
    add_subdirectory(src/vulkan)
endif()
