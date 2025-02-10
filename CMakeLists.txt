# cmake version to use
cmake_minimum_required(VERSION 3.12)

# project name
project(mdrom)

# some c++ flags
set(CMAKE_CXX_STANDARD 17)  # Need c++ 17 for inline variables
set(CMAKE_CXX_STANDARD_REQUIRED True)

# optimize for speed
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -pthread -fopenmp")
set(CMAKE_CXX_FLAGS_RELEASE "-DNDEBUG")  # Disable debugging in release builds

# where to look for packages and libraries
SET(yaml-cpp_DIR /home/sergio/repos/yaml-cpp_install/lib/cmake/yaml-cpp/)

# include headers
include_directories(inc)

# find the relevant packages
find_package(Eigen3)
find_package(yaml-cpp REQUIRED)

# add the exectuable files
add_executable(a src/main.cpp src/dynamics.cpp src/control.cpp)
add_executable(s src/main_3D.cpp src/dynamics_3D.cpp src/control_3D.cpp)

# link your executables
target_link_libraries(a Eigen3::Eigen yaml-cpp::yaml-cpp)
target_link_libraries(s Eigen3::Eigen yaml-cpp::yaml-cpp)

# also add the mdrom as a library for the whole project
add_library(control_3D SHARED src/control_3D.cpp)
target_link_libraries(control_3D Eigen3::Eigen yaml-cpp::yaml-cpp)
