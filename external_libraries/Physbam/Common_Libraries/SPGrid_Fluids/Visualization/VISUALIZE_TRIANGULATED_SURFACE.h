//#####################################################################
// Copyright (c) 2014, Raj Setaluri, Mridul Aanjneya
// Distributed under the FreeBSD license (see license.txt)
//#####################################################################
#ifndef __VISUALIZE_TRIANGULATED_SURFACE_h__
#define __VISUALIZE_TRIANGULATED_SURFACE_h__

#include <SPGrid_Fluids/Grids/GRID_HIERARCHY.h>
#include <PhysBAM_Tools/Data_Structures/TRIPLE.h>

namespace PhysBAM{
using namespace SPGrid;
//#####################################################################
// Class VISUALIZE_TRIANGULATED_SURFACE
//#####################################################################
template<class T_STRUCT,class T,int d>
class VISUALIZE_TRIANGULATED_SURFACE
{ 
    typedef GRID_HIERARCHY<T_STRUCT,T,d> T_HIERARCHY;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<unsigned>::type Flag_array_type;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<const unsigned>::type Const_flag_array_type;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<T>::type Data_array_type;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<const T>::type Const_data_array_type;
    typedef typename Flag_array_type::MASK Flag_array_mask;

    typedef VECTOR<T,d> TV;
    typedef VECTOR<int,d> TV_INT;
    typedef TRIPLE<int,int,TV_INT> CELL_INFO;
//#############################################################################
public:
    static void Visualize_Triangulated_Surface(STREAM_TYPE stream_type,T_HIERARCHY& hierarchy,const std::string output_directory,const int axis,const T h,const bool flat,const T scale=(T)1.);
//#############################################################################
};
}
#endif
