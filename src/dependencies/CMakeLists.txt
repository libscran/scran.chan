cmake_minimum_required(VERSION 3.14)

project(dependencies
    VERSION 1.0.0
    DESCRIPTION "C++ dependencies for the scran.chan shared library"
    LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(IGRAPH_WARNINGS_AS_ERRORS OFF CACHE BOOL "" FORCE)

# Avoid missing symbol errors when CMake gets tricked into using the system libraries.
set(IGRAPH_USE_INTERNAL_LAPACK ON CACHE BOOL "" FORCE) 
set(IGRAPH_USE_INTERNAL_ARPACK ON CACHE BOOL "" FORCE) 
set(CMAKE_POSITION_INDEPENDENT_CODE ON) # required to avoid ARPACK linking errors for shared libraries.

include(FetchContent)

FetchContent_Declare(
  tatami
  GIT_REPOSITORY https://github.com/LTLA/tatami
  GIT_TAG 573f1fb772db0580ff80f0318006af1555c0b24e
)
FetchContent_MakeAvailable(tatami)

FetchContent_Declare(
  scran
  GIT_REPOSITORY https://github.com/LTLA/libscran
  GIT_TAG dfbd83397d7ed2110c0ee046ac4b1749a5c76ac4
)
FetchContent_MakeAvailable(scran)

FetchContent_Declare(
  mnncorrect 
  GIT_REPOSITORY https://github.com/LTLA/CppMnnCorrect
  GIT_TAG 5ba5c790f01b3a676420892151791786dfb0a8d6
)
FetchContent_MakeAvailable(mnncorrect)

FetchContent_Declare(
  qdtsne
  GIT_REPOSITORY https://github.com/LTLA/qdtsne
  GIT_TAG 3113fca1d5d88307ab1dbafa6ca0c6e42cfd51ed
)
FetchContent_MakeAvailable(qdtsne)

FetchContent_Declare(
  umappp 
  GIT_REPOSITORY https://github.com/LTLA/umappp
  GIT_TAG e086d6811b300878902d02ce164a7721902bfef3
)
FetchContent_MakeAvailable(umappp)

add_library(dummy INTERFACE)

target_link_libraries(dummy INTERFACE scran mnncorrect qdtsne umappp igraph)

set_property(TARGET igraph PROPERTY POSITION_INDEPENDENT_CODE ON)
