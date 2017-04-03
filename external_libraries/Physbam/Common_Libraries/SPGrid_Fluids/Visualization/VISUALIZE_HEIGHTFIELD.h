//#####################################################################
// Copyright (c) 2014, Raj Setaluri, Mridul Aanjneya
// Distributed under the FreeBSD license (see license.txt)
//#####################################################################
#ifndef __VISUALIZE_HEIGHTFIELD_h__
#define __VISUALIZE_HEIGHTFIELD_h__

#include <SPGrid_Fluids/Grids/GRID_HIERARCHY.h>

namespace PhysBAM{
using namespace SPGrid;
//#####################################################################
// Class VISUALIZE_HEIGHTFIELD
//#####################################################################
template<class T_STRUCT,class T,int d>
class VISUALIZE_HEIGHTFIELD
{ 
    typedef GRID_HIERARCHY<T_STRUCT,T,d> T_HIERARCHY;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<unsigned>::type Flag_array_type;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<const unsigned>::type Const_flag_array_type;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<T>::type Data_array_type;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<const T>::type Const_data_array_type;
    
    typedef VECTOR<T,d> TV;
    typedef VECTOR<int,d> TV_INT;
// #############################################################################
public:
    static void Visualize_Heightfield(STREAM_TYPE stream_type,T_HIERARCHY& hierarchy,T T_STRUCT::* data_channel,const std::string output_directory,const int frame,const T scale=(T)1.);    
// #############################################################################
};
}
#endif
