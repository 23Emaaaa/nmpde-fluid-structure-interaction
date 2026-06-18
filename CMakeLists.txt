cmake_minimum_required(VERSION 3.12.0)
project(PDE_project LANGUAGES CXX C)

include(common/cmake-common.cmake)

add_executable(exerciseFSI src/exerciseFSI.cpp src/FSIProblem.cpp)
deal_ii_setup_target(exerciseFSI)
