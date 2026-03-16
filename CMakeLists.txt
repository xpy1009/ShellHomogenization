cmake_minimum_required(VERSION 3.18)

project(zh0818 CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_BUILD_TYPE "Release")

# to include header files
include_directories(${CMAKE_SOURCE_DIR}/include)

file(GLOB HEADERS "include/*.h" "include/*/*.h")
file(GLOB SOURCES "src/*.cpp" "src/*/*.cpp")
add_executable(zh0818 ${HEADERS} ${SOURCES})

# gmsh 2025.08.14
# use homebrew to install gmsh
# set up cmake as in https://gitlab.onelab.info/gmsh/ddm/-/blob/master/CMakeLists.txt?ref_type=heads
find_library(GMSH_LIB gmsh)
if(NOT GMSH_LIB)
  message(FATAL_ERROR "Could not find Gmsh library")
endif()
find_path(GMSH_INC gmsh.h)
if(NOT GMSH_INC)
  message(FATAL_ERROR "Could not find Gmsh header 'gmsh.h'")
endif()
include_directories(${GMSH_INC})
target_link_libraries(${PROJECT_NAME} ${GMSH_LIB})

# Eigen 
# libigl automatically download Eigen 3.4, while Eigen 5 is already available
find_package (Eigen3 REQUIRED NO_MODULE)
target_link_libraries(${PROJECT_NAME} Eigen3::Eigen)

# libigl https://github.com/libigl/libigl-example-project/tree/main
# include(libigl)
add_subdirectory(ext/libigl)
igl_include(glfw)
target_link_libraries(${PROJECT_NAME} 
    igl::glfw
)

# polyscope (igl::imgui causes build fail for polyscope)
# include(FetchContent)
# FetchContent_Declare(
#    polyscope
#    GIT_REPOSITORY https://github.com/nmwsharp/polyscope.git
#    GIT_TAG        master
# )
# FetchContent_MakeAvailable(polyscope)
add_subdirectory(ext/polyscope)
target_link_libraries(${PROJECT_NAME} polyscope)

# tbb 
# https://www.intel.com/content/www/us/en/docs/onetbb/get-started-guide/2021-11/integrate-onetbb.html
find_package(TBB REQUIRED)
target_link_libraries(${PROJECT_NAME} TBB::tbb)

# tactile
add_library(tactile STATIC ext/tactile/src/tiling.cpp)
target_include_directories(tactile PUBLIC ext/tactile/src)
target_link_libraries(${PROJECT_NAME} tactile)
target_link_libraries(tactile glm)

# clipper2
# install via homebrew
# https://vcpkg.link/ports/clipper2
find_package(Clipper2 REQUIRED)
target_link_libraries(${PROJECT_NAME} Clipper2::Clipper2)
