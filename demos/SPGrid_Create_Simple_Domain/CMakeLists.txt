add_executable(Create_2D_Domain
  Create_2D_Domain/main.cpp
)

target_link_libraries( Create_2D_Domain
  PUBLIC
  SPGrid
  SPGrid_Fluids
  Threading_Tools
)


add_executable(Create_3D_Domain
  Create_3D_Domain/main.cpp
)

target_link_libraries( Create_3D_Domain
  PUBLIC
  SPGrid
  SPGrid_Fluids
  Threading_Tools
)
