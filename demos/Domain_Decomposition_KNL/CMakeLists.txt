add_library( IS_DD_KNL
  Interface_Solver/EIGEN_INTERFACE.cpp
  Interface_Solver/GRID_HIERARCHY_INSTANTIATION_HELPER.cpp
  Interface_Solver/GRID_HIERARCHY_ITERATOR_INSTANTIATION_HELPER.cpp
  Interface_Solver/PROLONGATION_MATRIX_HELPER.cpp
  Interface_Solver/SIGMA_MULTIGRID_SOLVER.cpp
  Interface_Solver/SIGMA_SMOOTHER.cpp
  )
target_link_libraries( IS_DD_KNL
  PUBLIC
  pthread
  )

target_include_directories(IS_DD_KNL
  PUBLIC
  ${EIGEN3_INCLUDE_DIR}
)

target_compile_options( IS_DD_KNL
  PUBLIC
  "-mkl"
  "-pthread"
)


add_executable( Test3D_Domain_Decomposition_KNL 
  Test_Bench_3D/main.cpp
  SPGrid_Domain_Decomposition_Solver/CG_VECTOR.cpp
  SPGrid_Domain_Decomposition_Solver/CG_SYSTEM.cpp
  KNL_Libraries/KNL_OVERSEER.cpp
  KNL_Libraries/Kernels/Clear_Helper_Phi.cpp
  KNL_Libraries/Kernels/Correction_Helper_Phi.cpp
  KNL_Libraries/Kernels/Laplace_Helper_Phi.cpp
  KNL_Libraries/Kernels/Minus_Laplace_Helper_Phi.cpp
  KNL_Libraries/Kernels/Residual_Helper_Phi.cpp
  KNL_Libraries/Kernels/Restriction_Helper_Phi.cpp
  KNL_Libraries/Kernels/Prolongation_Helper_Phi.cpp
)

target_link_libraries( Test3D_Domain_Decomposition_KNL 
  IS_DD_KNL
  SPGrid
  SPGrid_Fluids
  Threading_Tools
  ${MEMKIND_LIBRARY}
)

target_include_directories( Test3D_Domain_Decomposition_KNL 
  PUBLIC
  ${MEMKIND_INCLUDE_DIR}
)

target_compile_options(Test3D_Domain_Decomposition_KNL 
  PUBLIC
  "-pthread"
)

set_target_properties(Test3D_Domain_Decomposition_KNL 
  PROPERTIES
  LINK_FLAGS
  "-mkl"
)
