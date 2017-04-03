//#####################################################################
// Copyright 2012-2013, Sean Bauer, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class GRID_HIERARCHY_AVERAGING
//#####################################################################
#include <SPGrid_Fluids/Interpolation/GRID_HIERARCHY_AVERAGING.h>
#include <SPGrid/Tools/SPGrid_Block_Iterator.h>
#include <SPGrid/Tools/SPGrid_Threading_Helper.h>
#include <SPGrid/Tools/SPGrid_Clear.h>
#include <SPGrid/Tools/SPGrid_Masked_Average_Offset_Grid.h>
#include <SPGrid/Tools/SPGrid_Downsample_Accumulate_Shared.h>
#include <SPGrid/Tools/SPGrid_Masked_Normalize.h>
#include <SPGrid/Tools/SPGrid_Upsample_Inject_Shared.h>
#include <SPGrid_Fluids/Grids/GRID_TOPOLOGY_HELPER.h>
#include <SPGrid_Fluids/Simulation/FLUIDS_SIMULATION_DATA.h>
#include <SPGrid_Fluids/Projection/GRID_HIERARCHY_PROJECTION.h>
#include <SPGrid_Fluids/Interpolation/CONSTRAIN_T_JUNCTION_NODES.h>

using namespace PhysBAM;
namespace PhysBAM{
extern int PhysBAM_number_of_threads;
}
//#####################################################################
// Function Average_Face_Velocities_To_Nodes
//#####################################################################
template<class T_STRUCT, class T,int d> void GRID_HIERARCHY_AVERAGING<T_STRUCT,T,d>::
Average_Face_Velocities_To_Nodes(T_HIERARCHY& hierarchy,const VECTOR<T T_STRUCT::*,d> face_velocities,
    const VECTOR<T T_STRUCT::*,d> node_velocities,unsigned T_STRUCT::* const flags_field,T T_STRUCT::* const weight_field)
{
    static const int nodes_per_face=GRID_TOPOLOGY_HELPER<typename Data_array_type::MASK>::nodes_per_face;

    for(int axis=1;axis<=d;axis++){
        
        // Compute offsets for nodes of each face
        int nodes_of_face_offsets[nodes_per_face];
        GRID_TOPOLOGY_HELPER<typename Data_array_type::MASK>::Nodes_Of_Face_Shadow_Grid_Offsets(nodes_of_face_offsets,axis);
        unsigned face_valid_mask=GRID_TOPOLOGY_HELPER<typename Data_array_type::MASK>::Face_Valid_Mask(axis);

        // Clear nodal values and weights, perform in-level averaging

        for(int level=1;level<=hierarchy.Levels();level++){

            const T weight=(T)(1u<<(1*(hierarchy.Levels()-level)));

            if(PhysBAM_number_of_threads)
                SPGrid_Computations::Threading_Helper<T_STRUCT,d>(hierarchy.Allocator(level),hierarchy.Blocks(level)).Run_Parallel(
                    SPGrid_Computations::Clear<T_STRUCT,T,d,2>(node_velocities(axis),weight_field),PhysBAM_number_of_threads);
            else
                SPGrid_Computations::Clear<T_STRUCT,T,d,2>(hierarchy.Allocator(level),hierarchy.Blocks(level),
                    node_velocities(axis),weight_field);

            if(PhysBAM_number_of_threads){
                SPGrid_Computations::Threading_Helper<T_STRUCT,d>(hierarchy.Allocator(level),hierarchy.Blocks(level)).Run_Parallel(
                    SPGrid_Computations::Masked_Average_Offset_Grid<nodes_per_face,T_STRUCT,T,unsigned,d>(face_velocities(axis),node_velocities(axis),weight_field,flags_field,
                        nodes_of_face_offsets,face_valid_mask,weight),hierarchy.Red_Partition(level));
                SPGrid_Computations::Threading_Helper<T_STRUCT,d>(hierarchy.Allocator(level),hierarchy.Blocks(level)).Run_Parallel(
                    SPGrid_Computations::Masked_Average_Offset_Grid<nodes_per_face,T_STRUCT,T,unsigned,d>(face_velocities(axis),node_velocities(axis),weight_field,flags_field,
                        nodes_of_face_offsets,face_valid_mask,weight),hierarchy.Black_Partition(level));}
            else
                SPGrid_Computations::Masked_Average_Offset_Grid<nodes_per_face,T_STRUCT,T,unsigned,d>(hierarchy.Allocator(level),hierarchy.Blocks(level),
                    face_velocities(axis),node_velocities(axis),weight_field,flags_field,nodes_of_face_offsets,face_valid_mask,weight);
        }

        // Accumulate weights, normalize and update across levels

        for(int level=1;level<hierarchy.Levels();level++)
            SPGrid_Computations::Downsample_Accumulate_Shared(hierarchy.Allocator(level),hierarchy.Allocator(level+1),
                hierarchy.Blocks(level),node_velocities(axis),weight_field,flags_field,(unsigned)SPGrid_Node_Coarse_Shared);

        for(int level=1;level<=hierarchy.Levels();level++)
            if(PhysBAM_number_of_threads)
                SPGrid_Computations::Threading_Helper<T_STRUCT,d>(hierarchy.Allocator(level),hierarchy.Blocks(level)).Run_Parallel(
                    SPGrid_Computations::Masked_Normalize<T_STRUCT,T,unsigned,d>(node_velocities(axis),weight_field,flags_field,(unsigned)SPGrid_Node_Active),PhysBAM_number_of_threads);
            else
                SPGrid_Computations::Masked_Normalize<T_STRUCT,T,unsigned,d>(hierarchy.Allocator(level),hierarchy.Blocks(level),
                    node_velocities(axis),weight_field,flags_field,(unsigned)SPGrid_Node_Active);

        for(int level=hierarchy.Levels()-1;level>=1;level--)
            SPGrid_Computations::Upsample_Inject_Shared(hierarchy.Allocator(level),hierarchy.Allocator(level+1),
                hierarchy.Blocks(level),node_velocities(axis),flags_field,(unsigned)SPGrid_Node_Coarse_Shared);
    }
}
//#####################################################################
// Function Average_Cell_Density_To_Nodes
//#####################################################################
template<class T_STRUCT, class T,int d> void GRID_HIERARCHY_AVERAGING<T_STRUCT,T,d>::
Average_Cell_Density_To_Nodes(T_HIERARCHY& hierarchy,T T_STRUCT::* const cell_density,
    T T_STRUCT::* const node_density,unsigned T_STRUCT::* const flags_field,T T_STRUCT::* const weight_field)
{
    static const int nodes_per_cell=GRID_TOPOLOGY_HELPER<typename Data_array_type::MASK>::nodes_per_cell;
    int nodes_of_cell_offsets[nodes_per_cell];
    GRID_TOPOLOGY_HELPER<typename Data_array_type::MASK>::Nodes_Of_Cell_Shadow_Grid_Offsets(nodes_of_cell_offsets);

    // Clear nodal values and weights, perform in-level averaging

    for(int level=1;level<=hierarchy.Levels();level++){

        const T weight=(T)(1u<<(1*(hierarchy.Levels()-level))); // CHECK THIS WEIGHT

        if(PhysBAM_number_of_threads)
            SPGrid_Computations::Threading_Helper<T_STRUCT,d>(hierarchy.Allocator(level),hierarchy.Blocks(level)).Run_Parallel(
                SPGrid_Computations::Clear<T_STRUCT,T,d,2>(node_density,weight_field),PhysBAM_number_of_threads);
        else
            SPGrid_Computations::Clear<T_STRUCT,T,d,2>(hierarchy.Allocator(level),hierarchy.Blocks(level),
                node_density,weight_field);

        if(PhysBAM_number_of_threads){
            SPGrid_Computations::Threading_Helper<T_STRUCT,d>(hierarchy.Allocator(level),hierarchy.Blocks(level)).Run_Parallel(
                SPGrid_Computations::Masked_Average_Offset_Grid<nodes_per_cell,T_STRUCT,T,unsigned,d>(cell_density,node_density,weight_field,flags_field,
                    nodes_of_cell_offsets,(unsigned)SPGrid_Cell_Type_Interior,weight),hierarchy.Red_Partition(level));
            SPGrid_Computations::Threading_Helper<T_STRUCT,d>(hierarchy.Allocator(level),hierarchy.Blocks(level)).Run_Parallel(
                SPGrid_Computations::Masked_Average_Offset_Grid<nodes_per_cell,T_STRUCT,T,unsigned,d>(cell_density,node_density,weight_field,flags_field,
                    nodes_of_cell_offsets,(unsigned)SPGrid_Cell_Type_Interior,weight),hierarchy.Black_Partition(level));}
        else
            SPGrid_Computations::Masked_Average_Offset_Grid<nodes_per_cell,T_STRUCT,T,unsigned,d>(hierarchy.Allocator(level),hierarchy.Blocks(level),
                cell_density,node_density,weight_field,flags_field,nodes_of_cell_offsets,(unsigned)SPGrid_Cell_Type_Interior,weight);
    }

    // Accumulate weights, normalize and update across levels
    
    for(int level=1;level<hierarchy.Levels();level++)
        SPGrid_Computations::Downsample_Accumulate_Shared(hierarchy.Allocator(level),hierarchy.Allocator(level+1),
            hierarchy.Blocks(level),node_density,weight_field,flags_field,(unsigned)SPGrid_Node_Coarse_Shared);

    for(int level=1;level<=hierarchy.Levels();level++)
        if(PhysBAM_number_of_threads)
            SPGrid_Computations::Threading_Helper<T_STRUCT,d>(hierarchy.Allocator(level),hierarchy.Blocks(level)).Run_Parallel(
                SPGrid_Computations::Masked_Normalize<T_STRUCT,T,unsigned,d>(node_density,weight_field,flags_field,(unsigned)SPGrid_Node_Active),PhysBAM_number_of_threads);
        else
            SPGrid_Computations::Masked_Normalize<T_STRUCT,T,unsigned,d>(hierarchy.Allocator(level),hierarchy.Blocks(level),
                node_density,weight_field,flags_field,(unsigned)SPGrid_Node_Active);

    for(int level=hierarchy.Levels()-1;level>=1;level--)
        SPGrid_Computations::Upsample_Inject_Shared(hierarchy.Allocator(level),hierarchy.Allocator(level+1),
            hierarchy.Blocks(level),node_density,flags_field,(unsigned)SPGrid_Node_Coarse_Shared);

    // Fix T-Junctions
    for(int level=1;level<=hierarchy.Levels();level++) // Don't need to do technically, but being safe, TODO: thread
        SPGrid_Computations::Masked_Clear<T_STRUCT,T,unsigned,d>(hierarchy.Allocator(level),hierarchy.Blocks(level),node_density,flags_field,(unsigned)SPGrid_Node_T_Junction);
    for(int level=hierarchy.Levels()-1;level>=1;level--)
        CONSTRAIN_T_JUNCTION_NODES<T_STRUCT,T,unsigned,d>::CONSTRAIN_T_JUNCTION_NODES(hierarchy,node_density,flags_field,level).Run(hierarchy.Allocator(level),hierarchy.Blocks(level));
}
//#####################################################################
// Function Average_Cell_Density_To_Faces
//#####################################################################
template<class T_STRUCT, class T,int d> void GRID_HIERARCHY_AVERAGING<T_STRUCT,T,d>::
Average_Cell_Density_To_Vertical_Faces(T_HIERARCHY& hierarchy,T T_STRUCT::* const cell_density,
    T T_STRUCT::* const vertical_face_density,unsigned T_STRUCT::* const flags_field,const int vertical_axis)
{
    // NOTE: Assuming that all valid faces which lie between different resolutions are marked as scaled !!

    static const T scale_uniform=(T).5;
    static const T scale_nonuniform=(T)(2./3.);
    
    const unsigned vertical_face_valid_mask=GRID_TOPOLOGY_HELPER<typename Data_array_type::MASK>::Face_Valid_Mask(vertical_axis);
    const unsigned vertical_face_scaled_mask=GRID_TOPOLOGY_HELPER<typename Data_array_type::MASK>::Face_Minus_Scaled_Mask(vertical_axis);
    const unsigned long negative_vertical_axis_offset=GRID_TOPOLOGY_HELPER<typename Data_array_type::MASK>::Negative_Axis_Vector_Offset(vertical_axis);

    GRID_HIERARCHY_PROJECTION<T_STRUCT,T,d>::Propagate_Ghost_Values(hierarchy,flags_field,cell_density); 
        
    for(int level=1;level<=hierarchy.Levels();level++){
        Const_data_array_type d=hierarchy.Allocator(level).Get_Const_Array(cell_density);
        Data_array_type vertical_face_d=hierarchy.Allocator(level).Get_Array(vertical_face_density);
        Const_flag_array_type flags=hierarchy.Allocator(level).Get_Const_Array(flags_field);
        for(SPGrid_Block_Iterator<typename Data_array_type::MASK> iterator(hierarchy.Blocks(level));iterator.Valid();iterator.Next()){
            const unsigned flag=iterator.Data(flags);
            if(flag & vertical_face_valid_mask){
                T result;

#if 0
                // debug
                // PHYSBAM_ASSERT(hierarchy.Set(level).Is_Set(Flag_array_type::MASK::Packed_Add(negative_vertical_axis_offset,iterator.Offset()),(SPGrid_Cell_Type_Interior|SPGrid_Cell_Type_Ghost|SPGrid_Cell_Type_Dirichlet)));
                // if(!(hierarchy.Set(level).Is_Set(iterator.Offset(),(SPGrid_Cell_Type_Interior|SPGrid_Cell_Type_Ghost|SPGrid_Cell_Type_Dirichlet))))
                //     {LOG::cout<<"level="<<level<<", index="<<iterator.Index()<<std::endl; exit(0);}                
                PHYSBAM_ASSERT((hierarchy.Set(level).Is_Set(Flag_array_type::MASK::Packed_Add(negative_vertical_axis_offset,iterator.Offset()),(SPGrid_Cell_Type_Interior|SPGrid_Cell_Type_Ghost|SPGrid_Cell_Type_Dirichlet)))
                               ||
                               (hierarchy.Set(level).Is_Set(iterator.Offset(),(SPGrid_Cell_Type_Interior|SPGrid_Cell_Type_Ghost|SPGrid_Cell_Type_Dirichlet))));
                // debug
#endif
                
                T cell_value=iterator.Data(d);
                T other_cell_value=iterator.Data(d,negative_vertical_axis_offset);
                if(flag & vertical_face_scaled_mask){
                    if(flag & SPGrid_Cell_Type_Ghost){
                        result = scale_nonuniform*other_cell_value + ((T)1.-scale_nonuniform)*cell_value;}
                    else{
#if 0
                        // debug
                        PHYSBAM_ASSERT(hierarchy.Set(level).Is_Set(Flag_array_type::MASK::Packed_Add(negative_vertical_axis_offset,iterator.Offset()),SPGrid_Cell_Type_Ghost));
                        // debug
#endif                        
                        result = scale_nonuniform*cell_value + ((T)1.-scale_nonuniform)*other_cell_value;}
                } else{
                    result = scale_uniform*(cell_value + other_cell_value);
                }
                iterator.Data(vertical_face_d) = result;}
        }
    }    
}
//#####################################################################
// Function Average_Node_Density_To_Faces
//#####################################################################
template<class T_STRUCT, class T,int d> void GRID_HIERARCHY_AVERAGING<T_STRUCT,T,d>::
Average_Node_Density_To_Vertical_Faces(T_HIERARCHY& hierarchy,T T_STRUCT::* const node_density,T T_STRUCT::* const vertical_face_density,unsigned T_STRUCT::* const flags_field,const int vertical_axis,unsigned mask)
{
    const unsigned long neighbor_offset=GRID_TOPOLOGY_HELPER<typename Flag_array_type::MASK>::Axis_Vector_Offset(vertical_axis);

    for(int level=1;level<=hierarchy.Levels();level++){
        Const_data_array_type node_d=hierarchy.Allocator(level).Get_Const_Array(node_density);
        Data_array_type vertical_face_d=hierarchy.Allocator(level).Get_Array(vertical_face_density);
        Const_flag_array_type flags=hierarchy.Allocator(level).Get_Const_Array(flags_field);
        for(SPGrid_Block_Iterator<typename Data_array_type::MASK> iterator(hierarchy.Blocks(level));iterator.Valid();iterator.Next())
            if(iterator.Data(flags)&mask)
                iterator.Data(vertical_face_d)=(T).5*(iterator.Data(node_d)+iterator.Data(node_d,neighbor_offset));}
}
//#####################################################################
template class GRID_HIERARCHY_AVERAGING<FLUIDS_SIMULATION_DATA<float>,float,2>;
template class GRID_HIERARCHY_AVERAGING<FLUIDS_SIMULATION_DATA<float>,float,3>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class GRID_HIERARCHY_AVERAGING<FLUIDS_SIMULATION_DATA<double>,double,2>;
template class GRID_HIERARCHY_AVERAGING<FLUIDS_SIMULATION_DATA<double>,double,3>;
#endif

