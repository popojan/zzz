cmake_minimum_required(VERSION 3.10)
project(zzz C)

set(CMAKE_C_STANDARD 11)
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Release (default), Debug, RelWithDebInfo" FORCE)
endif()

option(STATIC_LINKING "link FLINT/MPFR/GMP as static archives from ZZZ_DEPS_PREFIX (see deps/README.md)" OFF)
set(ZZZ_DEPS_PREFIX "${PROJECT_SOURCE_DIR}/deps/local" CACHE PATH
    "install prefix holding include/{flint,gmp.h,mpfr.h} and lib/lib{flint,mpfr,gmp}.a")

# All zeta-zero sources live under src/ (psi.c / pellreg*.c in the tree root are
# unrelated and not built here).
set(ZZZ_SRC ${PROJECT_SOURCE_DIR}/src)
include_directories(${ZZZ_SRC})

# FLINT >= 3 (arb is part of it); GMP and MPFR underneath.
if(STATIC_LINKING)
    include_directories(${ZZZ_DEPS_PREFIX}/include)
    set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
    set(_zzz_libdirs ${ZZZ_DEPS_PREFIX}/lib ${ZZZ_DEPS_PREFIX}/lib64)
    find_library(FLINT flint PATHS ${_zzz_libdirs} NO_DEFAULT_PATH)
    find_library(MPFR  mpfr  PATHS ${_zzz_libdirs} NO_DEFAULT_PATH)
    find_library(GMP   gmp   PATHS ${_zzz_libdirs} NO_DEFAULT_PATH)
else()
    find_library(FLINT flint)
    find_library(MPFR  mpfr)
    find_library(GMP   gmp)
endif()
foreach(_lib FLINT MPFR GMP)
    if(NOT ${_lib})
        message(FATAL_ERROR "${_lib} not found (STATIC_LINKING=${STATIC_LINKING}, ZZZ_DEPS_PREFIX=${ZZZ_DEPS_PREFIX}); see deps/README.md")
    endif()
endforeach()
set(ZZZ_LIBS ${FLINT} ${MPFR} ${GMP} m)   # dependency order matters for static archives

add_library(ghy STATIC ${ZZZ_SRC}/ghy.c ${ZZZ_SRC}/weil.c ${ZZZ_SRC}/loop.c)

add_executable(zzz     ${ZZZ_SRC}/main.c)
add_executable(zproxy  ${ZZZ_SRC}/zproxy.c)
add_executable(zghy    ${ZZZ_SRC}/zghy.c)
add_executable(zhad    ${ZZZ_SRC}/zhad.c)
add_executable(zhybrid ${ZZZ_SRC}/zhybrid.c)
target_link_libraries(ghy PUBLIC ${ZZZ_LIBS})
target_link_libraries(zproxy ${ZZZ_LIBS})
foreach(_t zzz zghy zhad zhybrid)
    target_link_libraries(${_t} ghy)
endforeach()

# -Wno-missing-field-initializers: the argp option table relies on trailing zeros
foreach(_t ghy zzz zproxy zghy zhad zhybrid)
    target_compile_options(${_t} PRIVATE -funroll-loops -Wall -Wextra -Wno-missing-field-initializers)
endforeach()

# ---- smoke tests: `ctest --test-dir build` (~12 s) -------------------------
enable_testing()
foreach(_case zero1 ghy boot weil grids loop)
    add_test(NAME smoke.${_case}
             COMMAND sh ${PROJECT_SOURCE_DIR}/tests/smoke.sh ${_case}
                     $<TARGET_FILE_DIR:zzz> ${PROJECT_SOURCE_DIR}/doc/ghy/first-20-zeros.txt)
endforeach()

# ---- experimental: Vulkan-accelerated batch method-B (opt-in) --------------
option(BUILD_VULKAN "build the experimental Vulkan batch method-B tool (zvk)" OFF)
if(BUILD_VULKAN)
    add_subdirectory(src/vulkan)
endif()
