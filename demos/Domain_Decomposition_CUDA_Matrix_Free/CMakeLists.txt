add_library( IS_DD_CUDA
  Interface_Solver/EIGEN_INTERFACE.cpp
  Interface_Solver/GRID_HIERARCHY_INSTANTIATION_HELPER.cpp
  Interface_Solver/GRID_HIERARCHY_ITERATOR_INSTANTIATION_HELPER.cpp
  Interface_Solver/PROLONGATION_MATRIX_HELPER.cpp
  Interface_Solver/SIGMA_MULTIGRID_SOLVER.cpp
  Interface_Solver/SIGMA_SMOOTHER.cpp
  )
target_link_libraries( IS_DD_CUDA
  PUBLIC
  pthread
  )

target_include_directories(IS_DD_CUDA
  PUBLIC
  ${EIGEN3_INCLUDE_DIR}
)

target_compile_options( IS_DD_CUDA
  PUBLIC
  "-mkl"
  "-pthread"
)

set(
    CUDA_NVCC_FLAGS
    ${CUDA_NVCC_FLAGS};
    -arch=compute_35 -code=sm_35 -Xcompiler --ptxas-options=-allow-expensive-optimizations=true,-O3
    )

cuda_add_library( NODES_CUDA
  Nodes/Overseer.cpp
  Nodes/Communicator.cpp
  Nodes/Client_Node_CUDA.cpp
)

target_compile_options( NODES_CUDA
  PUBLIC
  "-D_GLIBCXX_USE_NANOSLEEP"
  "-D__NVCC__"
  "-DCUDA_ENABLED"
)


cuda_add_library( DD_CUDA_KERNELS
  CUDA_Kernels/Clear_Helper_CUDA.cu
  CUDA_Kernels/Correction_Helper_CUDA.cu
  CUDA_Kernels/Interface_Block_Helper_CUDA.cu
  CUDA_Kernels/Laplace_Helper_CUDA.cu
  CUDA_Kernels/Minus_Laplace_Helper_CUDA.cu
  CUDA_Kernels/Norm_Helper_CUDA.cu
  CUDA_Kernels/Prolongation_Helper_CUDA.cu
  CUDA_Kernels/Residual_Helper_CUDA.cu
  CUDA_Kernels/Restriction_Helper_CUDA.cu
  CUDA_Kernels/Smoother_Helper_CUDA.cu
  CUDA_Kernels/Smoother_In_Cache_Helper_CUDA.cu
)


cuda_add_executable( Test3D_Domain_Decomposition_CUDA 
  Test_Bench_3D/main.cpp
  SPGrid_Domain_Decomposition_Solver/CG_SYSTEM.cpp
  SPGrid_Domain_Decomposition_Solver/CG_VECTOR.cpp
  )


target_link_libraries( Test3D_Domain_Decomposition_CUDA
  IS_DD_CUDA
  DD_CUDA_KERNELS
  NODES_CUDA
  SPGrid
  SPGrid_Fluids
  Threading_Tools
)

target_compile_options( Test3D_Domain_Decomposition_CUDA
  PUBLIC
  "-mkl"
  "-pthread"
)

set_target_properties(Test3D_Domain_Decomposition_CUDA 
  PROPERTIES
  LINK_FLAGS
  "-mkl"
)
