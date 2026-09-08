cmake_minimum_required(VERSION 4.0)

project(ShellHomogenization CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_BUILD_TYPE "Release")


# gmsh https://gitlab.onelab.info/gmsh/ddm/-/blob/master/CMakeLists.txt?ref_type=heads
find_library(GMSH_LIB gmsh)
find_path(GMSH_INC gmsh.h)
if(NOT GMSH_INC)
  message(FATAL_ERROR "Could not find Gmsh header 'gmsh.h' or library")
endif()
include_directories(${GMSH_INC})

# Eigen
find_package (Eigen3 5 REQUIRED NO_MODULE)

# libigl https://github.com/libigl/libigl-example-project/tree/main
add_subdirectory(ext/libigl)
igl_include(glfw)

# tactile
add_library(tactile STATIC ext/tactile/src/tiling.cpp)
target_include_directories(tactile PUBLIC ext/tactile/src)

# clipper2
find_package(Clipper2 REQUIRED)

# prebuilt libtorch uses its own openmp, which is in conflict with openmp in gmsh
# works if installed via homebrew
find_package(Torch REQUIRED)



add_library(homogenization)

# to include header files
include_directories(${CMAKE_SOURCE_DIR}/include)

file(GLOB HEADERS "include/*.h" "include/*/*.h")
file(GLOB SOURCES "src/*.cpp" "src/*/*.cpp")
# file(GLOB HEADERS 
#   "include/Element/Tet10.h" 
#   "include/Element/Shell.h" 
#   "include/Mesh/Meshlib.h" 
#   "include/Mesh/Isohedral.h"
#   "include/Material/Elasticity.h" 
#   "include/Material/Neural.h" 
#   "include/Newton.h" 
#   "include/Homogenization.h" 
#   "include/CylindricalProj.h" 
#   "include/Compression.h"
# )
# file(GLOB SOURCES 
#   "src/Element/Tet10.cpp" 
#   "src/Element/Shell.cpp" 
#   "src/Mesh/Meshlib.cpp" 
#   "src/Mesh/Isohedral.cpp" 
#   "src/Material/Elasticity.cpp" 
#   "src/Material/Neural.cpp" 
#   "src/Newton.cpp" 
#   "src/Homogenization.cpp" 
#   "src/CylindricalProj.cpp" 
#   "src/Compression.cpp"
# )
target_sources(homogenization PRIVATE ${HEADERS} ${SOURCES})
target_link_libraries(homogenization
  PUBLIC
  tactile
  Clipper2::Clipper2
  ${GMSH_LIB}
  Eigen3::Eigen
  igl::glfw
  ${TORCH_LIBRARIES}
)


add_executable(run_native simulation/run_native.cpp)
target_link_libraries(run_native PRIVATE homogenization)

add_executable(run_homogenized simulation/run_homogenized.cpp)
target_link_libraries(run_homogenized PRIVATE homogenization)
