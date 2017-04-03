//#####################################################################
// Copyright (c) 2014, Raj Setaluri, Mridul Aanjneya
// Distributed under the FreeBSD license (see license.txt)
//#####################################################################
#ifndef __SMOOTHER_h__
#define __SMOOTHER_h__

#include <SPGrid_Fluids/Grids/GRID_HIERARCHY.h>
#include <SPGrid_Fluids/Grids/GRID_TOPOLOGY_HELPER.h>

namespace PhysBAM{
using namespace SPGrid;
//#####################################################################
// Class SMOOTHER
//#####################################################################
template<class T_STRUCT,class T,int d>
class SMOOTHER
{ 
    typedef GRID_HIERARCHY<T_STRUCT,T,d> T_HIERARCHY;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<unsigned>::type Flag_array_type;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<const unsigned>::type Const_flag_array_type;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<T>::type Data_array_type;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<const T>::type Const_data_array_type;
    enum{faces_per_cell=GRID_TOPOLOGY_HELPER<typename Flag_array_type::MASK>::faces_per_cell};
    
    typedef VECTOR<T,d> TV;
    typedef VECTOR<int,d> TV_INT;
// #############################################################################
public:
    static void Compute_Residual(T_HIERARCHY& hierarchy,T T_STRUCT::* u_channel,T T_STRUCT::* b_channel,T T_STRUCT::* r_channel);
    static void Multiply_With_System_Matrix(T_HIERARCHY& hierarchy,T T_STRUCT::* u_channel,T T_STRUCT::* Lu_channel);
    static void Single_Level_Interior_Jacobi_Smoother(T_HIERARCHY& hierarchy,T T_STRUCT::* x_channel,T T_STRUCT::* b_channel,T T_STRUCT::* diagonal_channel,
        T T_STRUCT::* temp_channel,const int iterations,const T omega,const bool diagonal_is_inverted=false,const int* const smoothing_level=0);
    static void Single_Level_Boundary_Jacobi_Smoother(T_HIERARCHY& hierarchy,T T_STRUCT::* x_channel,T T_STRUCT::* b_channel,T T_STRUCT::* diagonal_channel,
        T T_STRUCT::* temp_channel,const int iterations,const T omega,const unsigned boundary_mask,const bool diagonal_is_inverted=false,const int* const smoothing_level=0);
    static void Single_Level_Jacobi_Smoother(T_HIERARCHY& hierarchy,T T_STRUCT::* x_channel,T T_STRUCT::* b_channel,T T_STRUCT::* delta_channel,T T_STRUCT::* dinv_channel,const int level,const int iterations,const T omega,const unsigned long mask,const std::pair<const unsigned long*,unsigned>& blocks);
// #############################################################################
};
}
#endif
