cmake_minimum_required(VERSION 3.5)
project(cartpole_examples LANGUAGES CXX)
add_definitions(-std=c++14)

find_package(casadi REQUIRED)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

add_library(cartpole_solver src/DirectTranscriptSolver.cpp src/DirectCollocationSolver.cpp) 
target_link_libraries(cartpole_solver casadi)

add_executable(mytest src/test.cpp)
target_link_libraries(mytest cartpole_solver)

add_executable(simpletest src/simple_test.cpp)
target_link_libraries(simpletest cartpole_solver)