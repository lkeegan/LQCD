cmake_minimum_required(VERSION 3.1)
project (ex05)

set(CMAKE_VERBOSE_MAKEFILE ON)
set(CMAKE_CXX_STANDARD 11)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -pedantic -Wextra")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --coverage")

# lattice library
set(LATTICE_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/lattice)
add_library(Lattice INTERFACE)
target_include_directories(Lattice INTERFACE ${LATTICE_INCLUDE_DIR})

# Make ising executable
add_executable(MC_Ising mc_ising.cpp)
target_link_libraries(MC_Ising Lattice)