//#####################################################################
// Copyright 2012-2013, Sean Bauer, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class GRID_HIERARCHY_INTERPOLATION
//#####################################################################
#ifndef __GRID_HIERARCHY_INTERPOLATION_h__
#define __GRID_HIERARCHY_INTERPOLATION_h__

#include <SPGrid_Fluids/Grids/GRID_HIERARCHY.h>
#include <SPGrid_Fluids/Grids/GRID_TOPOLOGY_HELPER.h>

namespace PhysBAM{

template<class T_STRUCT,class T,int d>
class GRID_HIERARCHY_INTERPOLATION
{
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<unsigned>::type Flag_array_type;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<T>::type Data_array_type;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<const unsigned>::type Const_flag_array_type;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<const T>::type Const_data_array_type;
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
    static T Cell_Interpolation_Helper(const T_HIERARCHY& hierarchy,unsigned long nodes_of_cell_offsets[nodes_per_cell],const int level,
        const unsigned long offset,const TV weights,T T_STRUCT::* density_field,T T_STRUCT::* node_density_field);

    static T Face_Interpolation_Helper(const T_HIERARCHY& hierarchy,unsigned long nodes_of_face_offsets[nodes_per_face],const int axis,const int level,const unsigned long offset,
        const VECTOR<T,d-1> weights,unsigned T_STRUCT::* flags_field,VECTOR<T T_STRUCT::*,d> face_velocities,
        VECTOR<T T_STRUCT::*,d> node_velocities); 
//#####################################################################
};
}
#endif
