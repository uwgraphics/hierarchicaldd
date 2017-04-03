//#####################################################################
// Copyright (c) 2014, Raj Setaluri, Mridul Aanjneya
// Distributed under the FreeBSD license (see license.txt)
//#####################################################################
#ifndef __MULTIGRID_SOLVER_h__
#define __MULTIGRID_SOLVER_h__

#include <SPGrid_Fluids/Grids/GRID_HIERARCHY.h>

namespace PhysBAM{
using namespace SPGrid;
//#####################################################################
// Class MULTIGRID_SOLVER
//#####################################################################
template<class T_STRUCT,class T,int d>
class MULTIGRID_SOLVER
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
    static void V_Cycle(ARRAY<T_HIERARCHY*>& multigrid_hierarchy,T T_STRUCT::* x_channel,T T_STRUCT::* b_channel,T T_STRUCT::* residual_channel,
        T T_STRUCT::* cg_q_channel,T T_STRUCT::* cg_r_channel,T T_STRUCT::* cg_s_channel,T T_STRUCT::* cg_k_channel,T T_STRUCT::* cg_z_channel,
        VECTOR<T T_STRUCT::*,d> L_channels,T T_STRUCT::* ic_diag_channel,const ARRAY<int>& substitution_partitions,
        T T_STRUCT::* mg_diag_channel,T T_STRUCT::* temp_channel,const int interior_smoother_iterations,const int boundary_smoother_iterations,
        const T interior_omega,const T boundary_omega,const bool diagonal_is_inverted=false);
// #############################################################################
};
}
#endif
