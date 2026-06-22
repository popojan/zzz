# Experimental Vulkan batch method-B tool (zvk).  Enabled only with
# -DBUILD_VULKAN=ON from the top-level build.  Prereqs (Linux):
#   * Vulkan loader + headers   (vulkan-tools, libvulkan-dev)
#   * glslangValidator          (glslang-tools)
#   * kompute                   (fetched below) -- needs network at configure time
#   * FLINT/GMP/MPFR            (already required by zzz)

enable_language(CXX)   # top-level project is C-only; zvk.cpp needs C++ (scoped here)

find_package(Vulkan REQUIRED)
find_program(GLSLANG NAMES glslangValidator glslang REQUIRED)  # name varies by version

include(FetchContent)
FetchContent_Declare(
    kompute
    GIT_REPOSITORY https://github.com/KomputeProject/kompute.git
    GIT_TAG        v0.9.0
)
# keep kompute's chatter off our stdout (ordinates only); set before MakeAvailable
set(KOMPUTE_OPT_LOG_LEVEL "Off" CACHE STRING "" FORCE)
FetchContent_MakeAvailable(kompute)

# compile the compute shader to SPIR-V next to the binary
set(SPV ${CMAKE_CURRENT_BINARY_DIR}/zeromb.comp.spv)
add_custom_command(
    OUTPUT  ${SPV}
    COMMAND ${GLSLANG} -V ${CMAKE_CURRENT_SOURCE_DIR}/zeromb.comp -o ${SPV}
    DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/zeromb.comp
    COMMENT "glslang: zeromb.comp -> SPIR-V")
add_custom_target(zvk_shader DEPENDS ${SPV})

add_executable(zvk zvk.cpp)
add_dependencies(zvk zvk_shader)
set_target_properties(zvk PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED ON)
target_compile_definitions(zvk PRIVATE ZVK_SPV="${SPV}")
target_link_libraries(zvk PRIVATE kompute Vulkan::Vulkan flint gmp mpfr m)
