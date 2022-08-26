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
  GIT_TAG ae74e0345303a2d2c6e70d599d72c0e02d346fb6
)
FetchContent_MakeAvailable(scran)

FetchContent_Declare(
  mnncorrect 
  GIT_REPOSITORY https://github.com/LTLA/CppMnnCorrect
  GIT_TAG 8ae1bc5686809a9287d0a9030ffb9b99d8baf1ff
)
FetchContent_MakeAvailable(mnncorrect)

FetchContent_Declare(
  qdtsne
  GIT_REPOSITORY https://github.com/LTLA/qdtsne
  GIT_TAG 72c81a00eda97cae2011abe06acce4525ecdd82a
)
FetchContent_MakeAvailable(qdtsne)

FetchContent_Declare(
  umappp 
  GIT_REPOSITORY https://github.com/LTLA/umappp
  GIT_TAG b65785b631c73e34008253f6dc2c04a19c8c6795
)
FetchContent_MakeAvailable(umappp)

add_library(dummy INTERFACE)

target_link_libraries(dummy INTERFACE scran mnncorrect qdtsne umappp igraph)

set_property(TARGET igraph PROPERTY POSITION_INDEPENDENT_CODE ON)
