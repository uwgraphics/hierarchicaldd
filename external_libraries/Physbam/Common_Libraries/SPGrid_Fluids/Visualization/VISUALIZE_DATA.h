//#####################################################################
// Copyright (c) 2012, Eftychios Sifakis, Sean Bauer
// Distributed under the FreeBSD license (see license.txt)
//#####################################################################
#ifndef __VISUALIZE_DATA_h__
#define __VISUALIZE_DATA_h__

#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <SPGrid_Fluids/Grids/GRID_HIERARCHY.h>

namespace PhysBAM{
using namespace SPGrid;
//#####################################################################
// Class VISUALIZE_DATA
//#####################################################################
template<class T_STRUCT,class T,int d>
class VISUALIZE_DATA
{
    typedef GRID_HIERARCHY<T_STRUCT,T,d> T_HIERARCHY;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<unsigned>::type Flag_array_type;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<const unsigned>::type Const_flag_array_type;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<T>::type Data_array_type;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<const T>::type Const_data_array_type;

    typedef GRID<VECTOR<T,d> > T_GRID;
    typedef VECTOR<int,d> T_INDEX;

// #############################################################################
public:
    static void Visualize_Data(T_HIERARCHY& hierarchy,const int frame,const std::string directory="output_data");
    static void Write_Output(STREAM_TYPE stream_type,T_HIERARCHY& hierarchy,const int frame,const std::string directory,const int level);
// #############################################################################
};
}
#endif
