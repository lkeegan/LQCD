cmake_minimum_required(VERSION 3.1)
project (ex02)

set(CMAKE_VERBOSE_MAKEFILE ON)
set(CMAKE_CXX_STANDARD 11)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -pedantic -Wextra")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --coverage")

# Catch unit testing library
set(CATCH_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/catch)
add_library(Catch INTERFACE)
target_include_directories(Catch INTERFACE ${CATCH_INCLUDE_DIR})

# lattice library
set(LATTICE_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/lattice)
add_library(Lattice INTERFACE)
target_include_directories(Lattice INTERFACE ${LATTICE_INCLUDE_DIR})

# Make unit test executable
add_executable(tests test/main.cpp)
target_link_libraries(tests Catch Lattice)

# Make ising executable
add_executable(MC_Oscillator mc_anharmonic_oscillator.cpp)
target_link_libraries(MC_Oscillator Lattice)

#enable_testing()
add_test(NAME Addition COMMAND tests)

#Make running the unit tests part of the build, so failed test = failed build
add_custom_command(TARGET tests POST_BUILD COMMAND tests -d yes)
