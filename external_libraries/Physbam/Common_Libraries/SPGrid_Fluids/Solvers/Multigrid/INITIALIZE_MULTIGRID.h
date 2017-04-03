//#####################################################################
// Copyright (c) 2014, Raj Setaluri, Mridul Aanjneya
// Distributed under the FreeBSD license (see license.txt)
//#####################################################################
#ifndef __INITIALIZE_MULTIGRID_h__
#define __INITIALIZE_MULTIGRID_h__

#include <SPGrid_Fluids/Grids/GRID_HIERARCHY.h>

namespace PhysBAM{
using namespace SPGrid;
//#####################################################################
// Class INITIALIZE_MULTIGRID
//#####################################################################
template<class T_STRUCT,class T,int d>
class INITIALIZE_MULTIGRID
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
    static void Initialize_Multigrid_Flags(T_HIERARCHY& hierarchy,unsigned T_STRUCT::* flags_channel,unsigned T_STRUCT::* multigrid_flags_channel);
    static void Initialize_Multigrid_Flags(T_HIERARCHY& base_hierarchy,ARRAY<T_HIERARCHY*>& multigrid_hierarchy,const int mg_levels);
    static void Initialize_Boundary_Flags(const ARRAY<T_HIERARCHY*>& multigrid_hierarchy,const int boundary_radius,const unsigned mask);
// #############################################################################
};
}
#endif
