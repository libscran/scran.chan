cmake_minimum_required(VERSION 3.14)

project(static_scran
    VERSION 1.0.0
    DESCRIPTION "Static version of the scran C++ library"
    LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)

add_library(static_scran 
    src/PerCellQCMetrics.cpp)

target_include_directories(static_scran PUBLIC include/)

add_subdirectory(extern)

target_link_libraries(static_scran scran)
