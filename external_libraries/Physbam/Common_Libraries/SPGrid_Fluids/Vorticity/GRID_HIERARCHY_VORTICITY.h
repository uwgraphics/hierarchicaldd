//#####################################################################
// Copyright 2014, Mridul Aanjaneya, Raj Setaluri.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class GRID_HIERARCHY_VORTICITY
//#####################################################################
#ifndef __GRID_HIERARCHY_VORTICITY_h__
#define __GRID_HIERARCHY_VORTICITY_h__

#include <SPGrid_Fluids/Grids/GRID_HIERARCHY.h>
#include <SPGrid_Fluids/Grids/GRID_TOPOLOGY_HELPER.h>

namespace PhysBAM{

template<class TV> class RIGID_GEOMETRY_COLLECTION;

template<class T_STRUCT,class T,int d>
class GRID_HIERARCHY_VORTICITY
{
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<unsigned>::type Flag_array_type;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<const unsigned>::type Const_flag_array_type;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<T>::type Data_array_type;
    typedef typename Flag_array_type::MASK Flag_array_mask;
    typedef GRID_HIERARCHY<T_STRUCT,T,d> T_HIERARCHY;
    typedef VECTOR<T,d> TV;
    typedef VECTOR<int,d> TV_INT;

    enum {
        nodes_per_cell=GRID_TOPOLOGY_HELPER<typename Data_array_type::MASK>::nodes_per_cell,
        nodes_per_face=GRID_TOPOLOGY_HELPER<typename Data_array_type::MASK>::nodes_per_face
    };

//#####################################################################
public:
    static void Compute_Vorticity(T_HIERARCHY& hierarchy,const VECTOR<T T_STRUCT::*,d> node_velocities,unsigned T_STRUCT::* const flags_field,const typename VECTOR<T T_STRUCT::*,d>::SPIN vorticity);
//#####################################################################
};
}
#endif
