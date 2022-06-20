cmake_minimum_required(VERSION 3.14)

project(dependencies
    VERSION 1.0.0
    DESCRIPTION "C++ dependencies for the scran.chan shared library"
    LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)

include(FetchContent)

FetchContent_Declare(
  scran
  GIT_REPOSITORY https://github.com/LTLA/libscran
  GIT_TAG 899efd34c08ef76545ad5a9b6b5eec5860281788
)
FetchContent_MakeAvailable(scran)

FetchContent_Declare(
  mnncorrect 
  GIT_REPOSITORY https://github.com/LTLA/CppMnnCorrect
  GIT_TAG ba8c2f5987e1a0eb94b6e0859040301b13e14acf
)
FetchContent_MakeAvailable(mnncorrect)

FetchContent_Declare(
  qdtsne
  GIT_REPOSITORY https://github.com/LTLA/qdtsne
  GIT_TAG 2a6fd4d78180625389dd75637939d74929b5668d
)
FetchContent_MakeAvailable(qdtsne)

FetchContent_Declare(
  umappp 
  GIT_REPOSITORY https://github.com/LTLA/umappp
  GIT_TAG e47f81b55c6adda777d22f361ec5def89248e5d7
)
FetchContent_MakeAvailable(umappp)

add_library(dummy INTERFACE)

target_link_libraries(dummy INTERFACE scran mnncorrect qdtsne umappp igraph)

set_property(TARGET igraph PROPERTY POSITION_INDEPENDENT_CODE ON)
