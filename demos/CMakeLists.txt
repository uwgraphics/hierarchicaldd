include_directories (${CMAKE_SOURCE_DIR}/external_libraries/Physbam/Public_Library/)
include_directories (${CMAKE_SOURCE_DIR}/external_libraries/Physbam/Common_Libraries/)

find_package(OpenMP REQUIRED)

if(OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

add_subdirectory(SPGrid_Create_Simple_Domain)

set( ENABLE_CPU_DEMO ON CACHE bool "Build the demo that only uses the CPU" )
set( ENABLE_TWO_LEVEL_CPU_DEMO ON CACHE bool "Build the two-level demo that only uses the CPU" )
set( ENABLE_CUDA_DEMO OFF CACHE bool "Build the demo that uses the GPU" )
set( ENABLE_KNL_DEMO OFF CACHE bool "Build the demo that uses Intel's Knights Landing Architecture" )

if( ENABLE_CPU_DEMO )
  add_subdirectory(Domain_Decomposition_CPU_Only)
endif()

if( ENABLE_TWO_LEVEL_CPU_DEMO )
  add_subdirectory(Domain_Decomposition_CPU_Only_Two_Levels)
endif()

if( ENABLE_CUDA_DEMO )
  add_subdirectory(Domain_Decomposition_CUDA_Matrix_Free)
endif()

if( ENABLE_KNL_DEMO )
  add_subdirectory(Domain_Decomposition_KNL)
endif()
