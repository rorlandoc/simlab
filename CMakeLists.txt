cmake_minimum_required(VERSION 3.20)
project(simlab VERSION 0.0.1 LANGUAGES CXX)

include(cmake/CPM.cmake)
CPMAddPackage(
  NAME eigen
  VERSION 3.4.0
  URL "https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.gz"
  DOWNLOAD_ONLY True
)
if(eigen_ADDED)
  add_library(eigen INTERFACE IMPORTED)
  target_include_directories(eigen INTERFACE ${eigen_SOURCE_DIR})
endif()

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(simlab 
  ${CMAKE_CURRENT_SOURCE_DIR}/src/simlab/material.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/simlab/integration.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/simlab/shape_function.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/simlab/element.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/simlab/mesh.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/simlab/solver.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/main.cpp)
target_link_libraries(simlab PRIVATE eigen)
target_include_directories(simlab PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
