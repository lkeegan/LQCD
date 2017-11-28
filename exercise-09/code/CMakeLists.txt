cmake_minimum_required(VERSION 3.1)
project (ex09)

set(CMAKE_VERBOSE_MAKEFILE ON)
set(CMAKE_CXX_STANDARD 11)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -pedantic -Wextra")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --coverage")

find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})

# lattice library
set(LATTICE_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/lattice)
add_library(Lattice INTERFACE)
target_include_directories(Lattice INTERFACE ${LATTICE_INCLUDE_DIR})

# Make ising executable
add_executable(MC_SU2 mc_SU2.cpp)
target_link_libraries(MC_SU2 Lattice)
