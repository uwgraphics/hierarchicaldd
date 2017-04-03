//#####################################################################
// Copyright (c) 2014, Raj Setaluri, Mridul Aanjneya
// Distributed under the FreeBSD license (see license.txt)
//#####################################################################
#include <SPGrid_Fluids/Visualization/VISUALIZE_HEIGHTFIELD.h>
#include <SPGrid_Fluids/Simulation/FLUIDS_SIMULATION_DATA.h>

using namespace PhysBAM;
//#####################################################################
// Visualize_Heightfield
//#####################################################################
template<class T_STRUCT,class T,int d> void VISUALIZE_HEIGHTFIELD<T_STRUCT,T,d>::
Visualize_Heightfield(STREAM_TYPE stream_type,T_HIERARCHY& hierarchy,T T_STRUCT::* data_channel,const std::string output_directory,const int frame,const T scale)
{PHYSBAM_FATAL_ERROR();}
//#####################################################################
template class VISUALIZE_HEIGHTFIELD<FLUIDS_SIMULATION_DATA<float>,float,3>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class VISUALIZE_HEIGHTFIELD<FLUIDS_SIMULATION_DATA<double>,double,3>;
#endif
