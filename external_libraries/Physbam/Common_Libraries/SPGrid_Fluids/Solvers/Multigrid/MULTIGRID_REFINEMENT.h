//#####################################################################
// Copyright (c) 2014, Raj Setaluri, Mridul Aanjneya
// Distributed under the FreeBSD license (see license.txt)
//#####################################################################
#ifndef __MULTIGRID_REFINEMENT_h__
#define __MULTIGRID_REFINEMENT_h__

#include <SPGrid_Fluids/Grids/GRID_HIERARCHY.h>
#include <SPGrid_Fluids/Grids/GRID_TOPOLOGY_HELPER.h>

namespace PhysBAM{
using namespace SPGrid;
//#####################################################################
// Class MULTIGRID_REFINEMENT
//#####################################################################
template<class T_STRUCT,class T,int d>
class MULTIGRID_REFINEMENT
{ 
    typedef GRID_HIERARCHY<T_STRUCT,T,d> T_HIERARCHY;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<unsigned>::type Flag_array_type;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<const unsigned>::type Const_flag_array_type;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<T>::type Data_array_type;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<const T>::type Const_data_array_type;
    typedef VECTOR<T,d> TV;
    typedef VECTOR<int,d> TV_INT;
    enum{nodes_per_cell=GRID_TOPOLOGY_HELPER<typename Flag_array_type::MASK>::nodes_per_cell};
    enum{parents_per_cell=nodes_per_cell};
    enum{restriction_stencil_denominator=(d==2)?16:64};
// #############################################################################
public:
    static void Restrict(T_HIERARCHY& fine_hierarchy,T_HIERARCHY& coarse_hierarchy,T T_STRUCT::* fine_data_channel,T T_STRUCT::* coarse_data_channel,VECTOR<int,2> finest_active_level);
    static void Prolongate(T_HIERARCHY& fine_hierarchy,T_HIERARCHY& coarse_hierarchy,T T_STRUCT::* fine_data_channel,T T_STRUCT::* coarse_data_channel,VECTOR<int,2> finest_active_level);
// #############################################################################
};
}
#endif
