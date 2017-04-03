//#####################################################################
// Copyright (c) 2014, Mridul Aanjaneya.
// Distributed under the FreeBSD license (see license.txt)
//#####################################################################
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT_PARTICLES_3D_DEFINITIONS.h>
#include <SPGrid_Fluids/Particles/VORTEX_PARTICLES.h>

template class OPENGL_COMPONENT_PARTICLES_3D<float,VORTEX_PARTICLES<VECTOR<float,3> >,float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class OPENGL_COMPONENT_PARTICLES_3D<double,VORTEX_PARTICLES<VECTOR<double,3> >,double>;
#endif
