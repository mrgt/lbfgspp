project(lbfgspp)
cmake_minimum_required(VERSION 2.8)
set(CMAKE_BUILD_TYPE Release)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
add_library(lbfgspp SHARED src/lbfgs.cpp)
