add_library( IS_DD_CPU_ONLY
  Interface_Solver/EIGEN_INTERFACE.cpp
  Interface_Solver/GRID_HIERARCHY_INSTANTIATION_HELPER.cpp
  Interface_Solver/GRID_HIERARCHY_ITERATOR_INSTANTIATION_HELPER.cpp
  Interface_Solver/PROLONGATION_MATRIX_HELPER.cpp
  Interface_Solver/SIGMA_MULTIGRID_SOLVER.cpp
  Interface_Solver/SIGMA_SMOOTHER.cpp
  )
target_link_libraries( IS_DD_CPU_ONLY
  PUBLIC
  pthread
  )

target_include_directories(IS_DD_CPU_ONLY
  PUBLIC
  ${EIGEN3_INCLUDE_DIR}
)

target_compile_options( IS_DD_CPU_ONLY
  PUBLIC
  "-mkl"
  "-pthread"
)


add_executable( Test3D_Domain_Decomposition_CPU_Only 
  Test_Bench_3D/main.cpp
  SPGrid_Domain_Decomposition_Solver/CG_VECTOR.cpp
  SPGrid_Domain_Decomposition_Solver/CG_SYSTEM.cpp
  SPGrid_Domain_Decomposition_Solver/CONJUGATE_GRADIENT_DUMPING.cpp
)

target_link_libraries( Test3D_Domain_Decomposition_CPU_Only 
  IS_DD_CPU_ONLY
  SPGrid
  SPGrid_Fluids
  Threading_Tools
)

target_compile_options(Test3D_Domain_Decomposition_CPU_Only 
  PUBLIC
  "-pthread"
)

set_target_properties(Test3D_Domain_Decomposition_CPU_Only 
  PROPERTIES
  LINK_FLAGS
  "-mkl"
)

add_executable( Test2D_Domain_Decomposition_CPU_Only 
  Test_Bench_2D/main.cpp
  SPGrid_Domain_Decomposition_Solver/CG_VECTOR.cpp
  SPGrid_Domain_Decomposition_Solver/CG_SYSTEM.cpp
  SPGrid_Domain_Decomposition_Solver/CONJUGATE_GRADIENT_DUMPING.cpp
)

target_link_libraries( Test2D_Domain_Decomposition_CPU_Only 
  IS_DD_CPU_ONLY
  SPGrid
  SPGrid_Fluids
  Threading_Tools
)

target_compile_options(Test2D_Domain_Decomposition_CPU_Only 
  PUBLIC
  "-pthread"
)

set_target_properties(Test2D_Domain_Decomposition_CPU_Only 
  PROPERTIES
  LINK_FLAGS
  "-mkl"
)
