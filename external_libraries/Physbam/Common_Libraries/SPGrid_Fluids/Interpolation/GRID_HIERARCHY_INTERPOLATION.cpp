//#####################################################################
// Copyright 2012-2013, Sean Bauer, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class GRID_HIERARCHY_INTERPOLATION
//#####################################################################
#include <SPGrid_Fluids/Interpolation/GRID_HIERARCHY_INTERPOLATION.h>

#include <SPGrid_Fluids/Simulation/FLUIDS_SIMULATION_DATA.h>
#include <PhysBAM_Tools/Interpolation/LINEAR_INTERPOLATION.h>
#include <SPGrid_Fluids/Grids/GRID_HIERARCHY_LOOKUP.h>

using namespace PhysBAM;

//#####################################################################
// Function Cell_Interpolation_Helper
//#####################################################################
template<class T_STRUCT, class T,int d> inline T GRID_HIERARCHY_INTERPOLATION<T_STRUCT,T,d>::
Cell_Interpolation_Helper(const T_HIERARCHY& hierarchy,unsigned long nodes_of_cell_offsets[nodes_per_cell],const int level,
    const unsigned long offset,const TV weights,T T_STRUCT::* density_field,T T_STRUCT::* node_density_field)
{
    T density_node_array[nodes_per_cell];    
    T interpolated_cell_value;
    
    for(int node=0;node<nodes_per_cell;node++){
        T node_val=hierarchy.Array(level,node_density_field)(Flag_array_mask::Packed_Add(offset,nodes_of_cell_offsets[node]));
        interpolated_cell_value+=node_val;
        density_node_array[node]=node_val;}
    interpolated_cell_value/=(T)nodes_per_cell;

    // weights are reversed to account for "reverse" ordering of nodes_of_cell_offsets array
    T multilinear_from_nodes=LINEAR_INTERPOLATION<T,T>::Linear(density_node_array,weights.Reversed());
    T phi=(T)2*min(weights.Min(),((T)1.-weights).Min());
    T actual_cell_value=hierarchy.Array(level,density_field)(offset);

    return multilinear_from_nodes + phi*(actual_cell_value-interpolated_cell_value);
}
//#####################################################################
// Function Face_Interpolation_Helper
//#####################################################################
template<class T_STRUCT, class T,int d> inline T GRID_HIERARCHY_INTERPOLATION<T_STRUCT,T,d>::
Face_Interpolation_Helper(const T_HIERARCHY& hierarchy,unsigned long nodes_of_face_offsets[nodes_per_face],const int axis,const int level,const unsigned long offset,
    const VECTOR<T,d-1> weights,unsigned T_STRUCT::* flags_field,VECTOR<T T_STRUCT::*,d> face_velocities,
    VECTOR<T T_STRUCT::*,d> node_velocities)
{    
    unsigned face_valid_mask=GRID_TOPOLOGY_HELPER<typename Data_array_type::MASK>::Face_Valid_Mask(axis);

    int face_level=level;
    unsigned long face_offset=offset;
    VECTOR<T,d-1> face_weights=weights;
    const bool face_found=GRID_HIERARCHY_LOOKUP<T_STRUCT,T,d>::Face_Lookup(hierarchy,face_offset,face_level,face_weights,face_valid_mask,axis);
	//PHYSBAM_ASSERT(face_found); // debug
    
    T velocity_node_array[nodes_per_face];
    T interpolated_face_value=T();

    for(int node=0;node<nodes_per_face;node++){
        const T node_val=hierarchy.Array(face_level,node_velocities(axis))(Flag_array_mask::Packed_Add(face_offset,nodes_of_face_offsets[node]));
        interpolated_face_value+=node_val;
        velocity_node_array[node]=node_val;}
    interpolated_face_value/=(T)nodes_per_face;
        
    // weights are reversed to account for "reverse" ordering of nodes_of_face_offsets array
    const T multilinear_from_nodes=LINEAR_INTERPOLATION<T,T>::Linear(velocity_node_array,face_weights.Reversed());
    const T phi=(T)2*min(face_weights.Min(),((T)1.-face_weights).Min());

    //PHYSBAM_ASSERT(hierarchy.Set(face_level).Is_Set(face_offset,face_valid_mask)); // debug
    const T actual_face_value=hierarchy.Array(face_level,face_velocities(axis))(face_offset);

    return multilinear_from_nodes + phi*(actual_face_value-interpolated_face_value);
}
//#####################################################################
template class GRID_HIERARCHY_INTERPOLATION<FLUIDS_SIMULATION_DATA<float>,float,2>;
template class GRID_HIERARCHY_INTERPOLATION<FLUIDS_SIMULATION_DATA<float>,float,3>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class GRID_HIERARCHY_INTERPOLATION<FLUIDS_SIMULATION_DATA<double>,double,2>;
template class GRID_HIERARCHY_INTERPOLATION<FLUIDS_SIMULATION_DATA<double>,double,3>;
#endif

