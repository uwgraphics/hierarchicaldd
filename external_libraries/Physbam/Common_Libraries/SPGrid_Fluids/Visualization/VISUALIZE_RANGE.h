//#####################################################################
// Copyright (c) 2012, Eftychios Sifakis, Sean Bauer
// Distributed under the FreeBSD license (see license.txt)
//#####################################################################
#ifndef __VISUALIZE_RANGE_h__
#define __VISUALIZE_RANGE_h__

#include <SPGrid_Fluids/Grids/GRID_HIERARCHY.h>
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/FACE_ARRAYS.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/GRID_ARRAYS_POLICY_UNIFORM.h>

namespace PhysBAM{
using namespace SPGrid;
//#####################################################################
// Class VISUALIZE_RANGE
//#####################################################################
template<class T_STRUCT,class T,int d>
class VISUALIZE_RANGE
{ 
    typedef float RW;
    typedef GRID_HIERARCHY<T_STRUCT,T,d> T_HIERARCHY;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<unsigned>::type Flag_array_type;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<const unsigned>::type Const_flag_array_type;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<T>::type Data_array_type;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<const T>::type Const_data_array_type;
    
    typedef VECTOR<int,d> T_INDEX;
    typedef VECTOR<T,d> TV;
    typedef GRID<TV> T_GRID;
    typedef ARRAY<T,T_INDEX> T_ARRAYS_SCALAR;
    typedef ARRAY<T,FACE_INDEX<d> > T_FACE_ARRAYS;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::FACE_ARRAYS::template REBIND<bool>::TYPE T_FACE_ARRAYS_BOOL;
    typedef ARRAY<bool,T_INDEX> T_ARRAYS_BOOL;
    /* typedef typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_SCALAR T_ARRAYS_SCALAR; */
    /* typedef typename GRID_ARRAYS_POLICY<T_GRID>::FACE_ARRAYS T_FACE_ARRAYS_SCALAR; */
    /* typedef typename T_FACE_ARRAYS_SCALAR::template REBIND<bool>::TYPE T_FACE_ARRAYS_BOOL; */
    /* typedef typename T_ARRAYS_SCALAR::template REBIND<bool>::TYPE T_ARRAYS_BOOL; */

// #############################################################################
public:
    static void Visualize_Range(const RANGE<TV>& range,const ARRAY<int>& frames,const std::string input_dir="output",const std::string output_dir="output_range");
    static void Visualize_Range_Helper(T_HIERARCHY& hierarchy,const RANGE<TV>& working_range,const ARRAY<RANGE<T_INDEX> >& working_cell_domains,const int output_frame,const std::string output_dir);
// #############################################################################
};
}
#endif
