cmake_minimum_required(VERSION 3.12.0)
project(project_title LANGUAGES CXX C)

# FIX: Comment this out until you actually add the file to the repo
# include(../common/cmake-common.cmake)

add_executable(lab-01 src/lab-01.cpp src/Poisson1D.cpp)

# FIX: Comment this out too (it depends on the file above)
# deal_ii_setup_target(lab-01)
