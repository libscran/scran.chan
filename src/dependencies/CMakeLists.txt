cmake_minimum_required(VERSION 3.14)

project(dependencies
    VERSION 1.0.0
    DESCRIPTION "C++ dependencies for the scran.chan shared library"
    LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(IGRAPH_WARNINGS_AS_ERRORS OFF CACHE BOOL "" FORCE)

include(FetchContent)

FetchContent_Declare(
  scran
  GIT_REPOSITORY https://github.com/LTLA/libscran
  GIT_TAG f0dafa6c79bfb96ee1804f81a4f51d7109fb072f
)
FetchContent_MakeAvailable(scran)

FetchContent_Declare(
  mnncorrect 
  GIT_REPOSITORY https://github.com/LTLA/CppMnnCorrect
  GIT_TAG 7c00002b315d6a7f9d6bead1fce307d539e911ce
)
FetchContent_MakeAvailable(mnncorrect)

FetchContent_Declare(
  qdtsne
  GIT_REPOSITORY https://github.com/LTLA/qdtsne
  GIT_TAG 9309c2e24180218363a5d898572051aa74a54377
)
FetchContent_MakeAvailable(qdtsne)

FetchContent_Declare(
  umappp 
  GIT_REPOSITORY https://github.com/LTLA/umappp
  GIT_TAG ebcf025444fc709c76485cbd1d943d8e13b3507b 
)
FetchContent_MakeAvailable(umappp)

add_library(dummy INTERFACE)

target_link_libraries(dummy INTERFACE scran mnncorrect qdtsne umappp igraph)

set_property(TARGET igraph PROPERTY POSITION_INDEPENDENT_CODE ON)
