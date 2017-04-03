//#####################################################################
// Copyright 2013, Raj Setaluri, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Face_Advection_Helper
//#####################################################################
#ifndef __Face_Advection_Helper_h__
#define __Face_Advection_Helper_h__

#include <PhysBAM_Geometry/Solids_Geometry/RIGID_GEOMETRY_COLLECTION.h>
#include <SPGrid/Core/SPGrid_Allocator.h>
#include <SPGrid_Fluids/Grids/GRID_HIERARCHY.h>
#include <SPGrid_Fluids/Grids/GRID_TOPOLOGY_HELPER.h>

namespace PhysBAM{

template<class T_STRUCT,class T_DATA,class T_FLAGS,int d>
class Face_Advection_Helper
{
    typedef T_DATA T;
    typedef GRID_HIERARCHY<T_STRUCT,T,d> T_HIERARCHY;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<unsigned>::type Flag_array_type;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<const T_FLAGS>::type Const_flag_array_type;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<T>::type Data_array_type;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<const T>::type Const_data_array_type;
    typedef typename Flag_array_type::MASK Flag_array_mask;
    typedef VECTOR<int,d> TV_INT;
    typedef VECTOR<T,d> TV;
    enum{nodes_per_face=GRID_TOPOLOGY_HELPER<typename Data_array_type::MASK>::nodes_per_face};
    enum{elements_per_block=Flag_array_mask::elements_per_block};
    enum{data_size=Data_array_type::MASK::field_size};

    const T_HIERARCHY& hierarchy;
    const VECTOR<T T_STRUCT::*,d>& face_velocities;
    const VECTOR<T T_STRUCT::*,d>& node_velocities;
    T T_STRUCT::* result_field;
    T_FLAGS T_STRUCT::* flags_field;
    const VECTOR<void*,d>& node_velocity_ptrs;
    unsigned long (&nodes_of_face_offsets)[nodes_per_face];
    const VECTOR<unsigned long,d>& other_face_offsets;
    const TV intra_cell_dX;
    const T dt;
    const unsigned face_advect_mask;
    const int level;
    const int axis;
    const RIGID_GEOMETRY_COLLECTION<TV>* const rigid_geometry_collection;

public:

    Face_Advection_Helper(const T_HIERARCHY& hierarchy_input,const VECTOR<T T_STRUCT::*,d>& face_velocities_input,const VECTOR<T T_STRUCT::*,d>& node_velocities_input,
        T T_STRUCT::* result_field_input,T_FLAGS T_STRUCT::* flags_field_input,const VECTOR<void*,d>& node_velocity_ptrs_input,
        unsigned long (&nodes_of_face_offsets_input)[nodes_per_face],const VECTOR<unsigned long,d>& other_face_offsets_input,const TV intra_cell_dX_input,const T dt_input,
        const unsigned face_advect_mask_input,const int level_input,const int axis_input,const RIGID_GEOMETRY_COLLECTION<TV>* const rigid_geometry_collection_input)
        :hierarchy(hierarchy_input),face_velocities(face_velocities_input),node_velocities(node_velocities_input),result_field(result_field_input),
         flags_field(flags_field_input),node_velocity_ptrs(node_velocity_ptrs_input),nodes_of_face_offsets(nodes_of_face_offsets_input),
         other_face_offsets(other_face_offsets_input),intra_cell_dX(intra_cell_dX_input),dt(dt_input),face_advect_mask(face_advect_mask_input),
         level(level_input),axis(axis_input),rigid_geometry_collection(rigid_geometry_collection_input)
    {}

    Face_Advection_Helper(SPGrid_Allocator<T_STRUCT,d>& allocator,const std::pair<const unsigned long*,unsigned>& blocks,
        const T_HIERARCHY& hierarchy_input,const VECTOR<T T_STRUCT::*,d>& face_velocities_input,const VECTOR<T T_STRUCT::*,d>& node_velocities_input,
        T T_STRUCT::* result_field_input,T_FLAGS T_STRUCT::* flags_field_input,const VECTOR<void*,d>& node_velocity_ptrs_input,
        unsigned long (&nodes_of_face_offsets_input)[nodes_per_face],const VECTOR<unsigned long,d>& other_face_offsets_input,const TV intra_cell_dX_input,const T dt_input,
        const unsigned face_advect_mask_input,const int level_input,const int axis_input,const RIGID_GEOMETRY_COLLECTION<TV>* const rigid_geometry_collection_input)
        :hierarchy(hierarchy_input),face_velocities(face_velocities_input),node_velocities(node_velocities_input),result_field(result_field_input),
         flags_field(flags_field_input),node_velocity_ptrs(node_velocity_ptrs_input),nodes_of_face_offsets(nodes_of_face_offsets_input),
         other_face_offsets(other_face_offsets_input),intra_cell_dX(intra_cell_dX_input),dt(dt_input),face_advect_mask(face_advect_mask_input),
         level(level_input),axis(axis_input),rigid_geometry_collection(rigid_geometry_collection_input)
    {Run(allocator,blocks);}

    void Run(SPGrid_Allocator<T_STRUCT,d>& allocator,const std::pair<const unsigned long*,unsigned>& blocks) const;
//#####################################################################
};
}
#endif
