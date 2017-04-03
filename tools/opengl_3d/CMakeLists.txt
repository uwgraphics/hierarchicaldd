add_executable(opengl_3d_viewer
main.cpp)

find_package(OpenGL REQUIRED)
#find_package(GLEW REQUIRED)
find_package(GLUT)
find_package(X11 REQUIRED)

target_include_directories (opengl_3d_viewer PUBLIC ${CMAKE_SOURCE_DIR}/external_libraries/Physbam/Public_Library/)
target_include_directories (opengl_3d_viewer PUBLIC ${CMAKE_SOURCE_DIR}/external_libraries/Physbam/Common_Libraries/)

target_link_libraries( opengl_3d_viewer
  Visualization
  ${OPENGL_LIBRARIES}
#  GLEW
  ${X11_X11_LIB}
  ${GLUT_LIBRARIES}
)


