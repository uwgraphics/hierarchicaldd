cmake_minimum_required(VERSION 3.0)
set (CMAKE_CXX_STANDARD 11)
project (Hierarchical_Domain_Decomposition_Demo)
set(BUILD_SHARED_LIBS ON) 

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin-demos)

set(EIGEN3_INCLUDE_DIR "${CMAKE_SOURCE_DIR}/external_libraries/eigen-3.3.3" CACHE PATH "Path to Eigen3 include directory")
set(MEMKIND_LIBRARY "" CACHE FILEPATH "Path to libmemkind install." )
set(MEMKIND_INCLUDE_DIR "" CACHE PATH "Path to libmemkind header files." )

find_package( CUDA REQUIRED)

add_definitions(-DCOMPILE_WITHOUT_RLE_SUPPORT)
add_definitions(-DCOMPILE_WITHOUT_ZLIB_SUPPORT)
add_definitions(-DCOMPILE_WITHOUT_DYADIC_SUPPORT)
add_definitions(-DCOMPILE_WITHOUT_DOUBLE_SUPPORT)
add_definitions(-DUSE_PTHREADS)
add_definitions(-DCOMPILE_ID_TYPES_AS_INT)

add_subdirectory(external_libraries/Physbam)
add_subdirectory(tools/opengl_3d)
add_subdirectory(demos)

if(NOT "${CMAKE_CXX_COMPILER_ID}" MATCHES "Intel")
  message(FATAL_ERROR "Error: This demo requires the Intel icc compiler.")
endif()


