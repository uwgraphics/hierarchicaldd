//#####################################################################
// Copyright (c) 2014, Raj Setaluri, Mridul Aanjneya
// Distributed under the FreeBSD license (see license.txt)
//#####################################################################
#include <SPGrid_Fluids/Visualization/VISUALIZE_TRIANGULATED_SURFACE.h>
#include <SPGrid_Fluids/Simulation/FLUIDS_SIMULATION_DATA.h>
using namespace PhysBAM;
//#####################################################################
// Visualize_Triangulated_Surface
//#####################################################################
template<class T_STRUCT,class T,int d> void VISUALIZE_TRIANGULATED_SURFACE<T_STRUCT,T,d>::
Visualize_Triangulated_Surface(STREAM_TYPE stream_type,T_HIERARCHY& hierarchy,const std::string filename,const int axis,const T h,const bool flat,const T scale)
{PHYSBAM_FATAL_ERROR("Not valid for 3D! Only for 2D!");}
//#####################################################################
template class VISUALIZE_TRIANGULATED_SURFACE<FLUIDS_SIMULATION_DATA<float>,float,3>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class VISUALIZE_TRIANGULATED_SURFACE<FLUIDS_SIMULATION_DATA<double>,double,3>;
#endif
