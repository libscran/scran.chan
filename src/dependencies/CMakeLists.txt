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
  GIT_TAG 142468ce7cdc17bbaa390b595d63276483499e64
)
FetchContent_MakeAvailable(scran)

FetchContent_Declare(
  mnncorrect 
  GIT_REPOSITORY https://github.com/LTLA/CppMnnCorrect
  GIT_TAG 3da4e85887beaf10d7197b0949c5ac477bf0470b
)
FetchContent_MakeAvailable(mnncorrect)

FetchContent_Declare(
  qdtsne
  GIT_REPOSITORY https://github.com/LTLA/qdtsne
  GIT_TAG 00de30a3b2d45edbab06f34a39771f4333fee1a7 
)
FetchContent_MakeAvailable(qdtsne)

FetchContent_Declare(
  umappp 
  GIT_REPOSITORY https://github.com/LTLA/umappp
  GIT_TAG 46e169b8df80725abca76e8007840ec17c54c819
)
FetchContent_MakeAvailable(umappp)

add_library(dummy INTERFACE)

target_link_libraries(dummy INTERFACE scran mnncorrect qdtsne umappp igraph)

set_property(TARGET igraph PROPERTY POSITION_INDEPENDENT_CODE ON)
