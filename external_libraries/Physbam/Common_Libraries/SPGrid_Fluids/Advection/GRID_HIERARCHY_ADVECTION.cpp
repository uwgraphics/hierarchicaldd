//#####################################################################
// Copyright 2012-2013, Sean Bauer, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class GRID_HIERARCHY_ADVECTION
//#####################################################################
#include <SPGrid_Fluids/Advection/GRID_HIERARCHY_ADVECTION.h>
#include <PhysBAM_Geometry/Implicit_Objects/IMPLICIT_OBJECT.h>
#include <SPGrid_Fluids/Simulation/FLUIDS_SIMULATION_DATA.h>
#include <SPGrid/Tools/SPGrid_Clear.h>
#include <SPGrid/Tools/SPGrid_Copy.h>
#include <SPGrid/Tools/SPGrid_Threading_Helper.h>
#include <SPGrid_Fluids/Advection/Density_Advection_Helper.h>
#include <SPGrid_Fluids/Advection/Face_Advection_Helper.h>

using namespace PhysBAM;

namespace PhysBAM{
extern int PhysBAM_number_of_threads;
}

//#####################################################################
// Function Advect_Face_Velocities
//#####################################################################
template<class T_STRUCT, class T,int d> void GRID_HIERARCHY_ADVECTION<T_STRUCT,T,d>::
Advect_Face_Velocities(T_HIERARCHY& hierarchy,const VECTOR<T T_STRUCT::*,d> face_velocities,
    const VECTOR<T T_STRUCT::*,d> node_velocities,unsigned T_STRUCT::* const flags_field,T T_STRUCT::* const temp_field,
    const T dt,const RIGID_GEOMETRY_COLLECTION<TV>* const rigid_geometry_collection)
{
    // get d oppsite face offsets
    VECTOR<unsigned long,d> other_face_offsets;
    for(int v=1;v<=d;v++) other_face_offsets(v)=Flag_array_mask::Linear_Offset(std_array<int,d>(TV_INT::Axis_Vector(v)));

    // advect velocities for each axis
    for(int axis=1;axis<=d;axis++){
        unsigned face_advect_mask=GRID_TOPOLOGY_HELPER<typename Data_array_type::MASK>::Face_Active_Mask(axis);
        unsigned long nodes_of_face_offsets[nodes_per_face];
        GRID_TOPOLOGY_HELPER<typename Data_array_type::MASK>::Nodes_Of_Face_Offsets(nodes_of_face_offsets,axis);
        TV intra_cell_dX=(T).5*(TV::All_Ones_Vector()-TV::Axis_Vector(axis));
        
        // clear temp field
        for(int level=1;level<=hierarchy.Levels();level++)
            if(PhysBAM_number_of_threads)
                SPGrid_Computations::Threading_Helper<T_STRUCT,d>(hierarchy.Allocator(level),hierarchy.Blocks(level)).Run_Parallel(
                    SPGrid_Computations::Clear<T_STRUCT,T,d>(temp_field),PhysBAM_number_of_threads);
            else
                SPGrid_Computations::Clear<T_STRUCT,T,d>(hierarchy.Allocator(level),hierarchy.Blocks(level),temp_field);        

        // do advection for each level
        for(int level=1;level<=hierarchy.Levels();level++){
            VECTOR<void*,d> node_velocity_ptrs;
            for(int v=1;v<=d;v++) node_velocity_ptrs(v)=hierarchy.Array(level,node_velocities(v)).Get_Data_Ptr();
            if(PhysBAM_number_of_threads)
                SPGrid_Computations::Threading_Helper<T_STRUCT,d>(hierarchy.Allocator(level),hierarchy.Blocks(level)).Run_Parallel(
                    Face_Advection_Helper<T_STRUCT,T,unsigned,d>(hierarchy,face_velocities,node_velocities,temp_field,flags_field,node_velocity_ptrs,
                        nodes_of_face_offsets,other_face_offsets,intra_cell_dX,dt,face_advect_mask,level,axis,rigid_geometry_collection),PhysBAM_number_of_threads);
            else
                Face_Advection_Helper<T_STRUCT,T,unsigned,d>(hierarchy.Allocator(level),hierarchy.Blocks(level),
                    hierarchy,face_velocities,node_velocities,temp_field,flags_field,node_velocity_ptrs,
                    nodes_of_face_offsets,other_face_offsets,intra_cell_dX,dt,face_advect_mask,level,axis,rigid_geometry_collection);}

        // for each level, copy temp-->face_velocities(axis)
        for(int level=1;level<=hierarchy.Levels();level++)
            if(PhysBAM_number_of_threads)
                SPGrid_Computations::Threading_Helper<T_STRUCT,d>(hierarchy.Allocator(level),hierarchy.Blocks(level)).Run_Parallel(
                    SPGrid_Computations::Masked_Copy<T_STRUCT,T,unsigned,d>(temp_field,face_velocities(axis),&T_STRUCT::flags,face_advect_mask),PhysBAM_number_of_threads);
            else
                SPGrid_Computations::Masked_Copy<T_STRUCT,T,unsigned,d>(hierarchy.Allocator(level),hierarchy.Blocks(level),temp_field,face_velocities(axis),&T_STRUCT::flags,face_advect_mask);
    }
    
}
//#####################################################################
// Function Advect_Densities
//#####################################################################
template<class T_STRUCT, class T,int d> void GRID_HIERARCHY_ADVECTION<T_STRUCT,T,d>::
Advect_Densities(T_HIERARCHY& hierarchy,const VECTOR<T T_STRUCT::*,d> face_velocities,
    const VECTOR<T T_STRUCT::*,d> node_velocities,T T_STRUCT::* density_field,T T_STRUCT::* node_density_field,
    unsigned T_STRUCT::* const flags_field,T T_STRUCT::* const temp_field,
    const T dt,const RIGID_GEOMETRY_COLLECTION<TV>* const rigid_geometry_collection)
{
    // where to advect density?
    unsigned density_advect_mask=SPGrid_Cell_Type_Interior;

    // offset to center of cell
    const TV intra_cell_dX=(T).5*TV::All_Ones_Vector();

    // compute offsets for nodes of cell
    unsigned long nodes_of_cell_offsets[nodes_per_cell];
    GRID_TOPOLOGY_HELPER<typename Data_array_type::MASK>::Nodes_Of_Cell_Offsets(nodes_of_cell_offsets);

    // clear temp field
    for(int level=1;level<=hierarchy.Levels();level++)
        if(PhysBAM_number_of_threads)
            SPGrid_Computations::Threading_Helper<T_STRUCT,d>(hierarchy.Allocator(level),hierarchy.Blocks(level)).Run_Parallel(
                SPGrid_Computations::Clear<T_STRUCT,T,d>(temp_field),PhysBAM_number_of_threads);
        else
            SPGrid_Computations::Clear<T_STRUCT,T,d>(hierarchy.Allocator(level),hierarchy.Blocks(level),temp_field);        

    // do advection for each level
    for(int level=1;level<=hierarchy.Levels();level++){
        // node velocity pointers
        VECTOR<void*,d> node_velocity_ptrs;
        for(int v=1;v<=d;v++) node_velocity_ptrs(v)=hierarchy.Array(level,node_velocities(v)).Get_Data_Ptr();
        if(PhysBAM_number_of_threads)
            SPGrid_Computations::Threading_Helper<T_STRUCT,d>(hierarchy.Allocator(level),hierarchy.Blocks(level)).Run_Parallel(
                Density_Advection_Helper<T_STRUCT,T,unsigned,d>(hierarchy,density_field,node_density_field,temp_field,flags_field,
                    node_velocity_ptrs,nodes_of_cell_offsets,intra_cell_dX,dt,density_advect_mask,level,rigid_geometry_collection),PhysBAM_number_of_threads);
        else
            Density_Advection_Helper<T_STRUCT,T,unsigned,d>(hierarchy.Allocator(level),hierarchy.Blocks(level),
                hierarchy,density_field,node_density_field,temp_field,flags_field,node_velocity_ptrs,
                nodes_of_cell_offsets,intra_cell_dX,dt,density_advect_mask,level,rigid_geometry_collection);}

    // for each level, copy from temp to density
	for(int level=1;level<=hierarchy.Levels();level++)
        if(PhysBAM_number_of_threads)
            SPGrid_Computations::Threading_Helper<T_STRUCT,d>(hierarchy.Allocator(level),hierarchy.Blocks(level)).Run_Parallel(
                SPGrid_Computations::Masked_Copy<T_STRUCT,T,unsigned,d>(temp_field,density_field,&T_STRUCT::flags,density_advect_mask),PhysBAM_number_of_threads);
        else
            SPGrid_Computations::Masked_Copy<T_STRUCT,T,unsigned,d>(hierarchy.Allocator(level),hierarchy.Blocks(level),temp_field,density_field,&T_STRUCT::flags,density_advect_mask);
}
//#####################################################################
template class GRID_HIERARCHY_ADVECTION<FLUIDS_SIMULATION_DATA<float>,float,2>;
template class GRID_HIERARCHY_ADVECTION<FLUIDS_SIMULATION_DATA<float>,float,3>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class GRID_HIERARCHY_ADVECTION<FLUIDS_SIMULATION_DATA<double>,double,2>;
template class GRID_HIERARCHY_ADVECTION<FLUIDS_SIMULATION_DATA<double>,double,3>;
#endif
