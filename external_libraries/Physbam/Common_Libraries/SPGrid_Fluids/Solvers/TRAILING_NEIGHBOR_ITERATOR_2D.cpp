//#####################################################################
// Copyright 2013, Sean Bauer, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class TRAILING_NEIGHBOR_ITERATOR
//#####################################################################
#include <SPGrid_Fluids/Solvers/TRAILING_NEIGHBOR_ITERATOR.h>
#include <SPGrid_Fluids/Simulation/FLUIDS_SIMULATION_DATA.h>
#include <SPGrid_Fluids/Simulation/FLUIDS_SIMULATION_FLAGS.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T_STRUCT,class T> TRAILING_NEIGHBOR_ITERATOR<T_STRUCT,T,2>::
TRAILING_NEIGHBOR_ITERATOR(Hierarchy_type& hierarchy_in,CELL_ID cid_in,unsigned T_STRUCT::* flags_field_in,T T_STRUCT::* data_field_in,T T_STRUCT::* diag_field_in,T T_STRUCT::* Ux_field_in,T T_STRUCT::* Uy_field_in) 
    :hierarchy(hierarchy_in)
{
    cid = cid_in;
    flags_field = flags_field_in;
    data_field = data_field_in;
    diag_field = diag_field_in;
    Ux_field = Ux_field_in;
    Uy_field = Uy_field_in;

    // Initialize list of neighbors
    Enumerate_Neighbors();

    index = 1; // index within neighbor array
}
//#####################################################################
// Function Enumerate_Neighbors
//#####################################################################
template<class T_STRUCT,class T> void TRAILING_NEIGHBOR_ITERATOR<T_STRUCT,T,2>::
Enumerate_Neighbors()
{
    unsigned flags = hierarchy.Array(cid.level, flags_field)(cid.offset);

    // Find same or coarser neighbors first
    if(flags & SPGrid_Face_Minus_X_Active)
    {
        int neighbor_level = cid.level;
        unsigned long neighbor_offset = Flag_array_mask::template Packed_OffsetXdim<-1>(cid.offset);
        unsigned long coef_offset = cid.offset;
        bool use_neighbor = (neighbor_offset>cid.offset) || (hierarchy.Array(cid.level, flags_field)(neighbor_offset)&SPGrid_Cell_Type_Ghost);
        if(use_neighbor)
        {
            T* data_ptr = &hierarchy.Array(cid.level, data_field)(neighbor_offset);
            T* coef_ptr = &hierarchy.Array(cid.level, Ux_field)(coef_offset);
            neighbors.Append(NEIGHBOR(neighbor_level, neighbor_offset, data_ptr, coef_ptr));
        }
    }

    if(flags & SPGrid_Face_Plus_X_Active)
    {
        int neighbor_level = cid.level;
        unsigned long neighbor_offset = Flag_array_mask::template Packed_OffsetXdim<1>(cid.offset);
        unsigned long coef_offset = neighbor_offset;
        
        bool use_neighbor = (neighbor_offset>cid.offset) || (hierarchy.Array(cid.level, flags_field)(neighbor_offset)&SPGrid_Cell_Type_Ghost);
        if(use_neighbor)
        {
            T* data_ptr = &hierarchy.Array(cid.level, data_field)(neighbor_offset);
            T* coef_ptr = &hierarchy.Array(cid.level, Ux_field)(coef_offset);
            neighbors.Append(NEIGHBOR(neighbor_level, neighbor_offset, data_ptr, coef_ptr));
        }
    }
    
    if(flags & SPGrid_Face_Minus_Y_Active)
    {
        int neighbor_level = cid.level;
        unsigned long neighbor_offset = Flag_array_mask::template Packed_OffsetYdim<-1>(cid.offset);
        unsigned long coef_offset = cid.offset;
        bool use_neighbor = (neighbor_offset>cid.offset) || (hierarchy.Array(cid.level, flags_field)(neighbor_offset)&SPGrid_Cell_Type_Ghost);
        if(use_neighbor)
        {
            T* data_ptr = &hierarchy.Array(cid.level, data_field)(neighbor_offset);
            T* coef_ptr = &hierarchy.Array(cid.level, Uy_field)(coef_offset);
            neighbors.Append(NEIGHBOR(neighbor_level, neighbor_offset, data_ptr, coef_ptr));
        }
    }

    if(flags & SPGrid_Face_Plus_Y_Active)
    {
        int neighbor_level = cid.level;
        unsigned long neighbor_offset = Flag_array_mask::template Packed_OffsetYdim<1>(cid.offset);
        unsigned long coef_offset = neighbor_offset;
        bool use_neighbor = (neighbor_offset>cid.offset) || (hierarchy.Array(cid.level, flags_field)(neighbor_offset)&SPGrid_Cell_Type_Ghost);
        if(use_neighbor)
        {
            T* data_ptr = &hierarchy.Array(cid.level, data_field)(neighbor_offset);
            T* coef_ptr = &hierarchy.Array(cid.level, Uy_field)(coef_offset);
            neighbors.Append(NEIGHBOR(neighbor_level, neighbor_offset, data_ptr, coef_ptr));
        }
    }
}
//#####################################################################
// Function Coarse_Cell_Finder
//#####################################################################
template<class T_STRUCT,class T> void TRAILING_NEIGHBOR_ITERATOR<T_STRUCT,T,2>::
Coarse_Cell_Finder(int& neighbor_level,unsigned long& neighbor_offset)
{
    // Also check for ghost
    while(neighbor_level<=hierarchy.Levels())
    {
        unsigned flags = hierarchy.Array(neighbor_level, flags_field)(neighbor_offset);
        if(flags & SPGrid_Cell_Type_Interior)
            return;

        if(flags & SPGrid_Cell_Type_Ghost)
        {
            neighbor_offset = Flag_array_mask::DownsampleOffset(neighbor_offset);
            neighbor_level++;
            continue;
        }
        PHYSBAM_FATAL_ERROR("Inconsistent edge flag and cell topology");
    }
}
//#####################################################################
template class TRAILING_NEIGHBOR_ITERATOR<FLUIDS_SIMULATION_DATA<float>,float,2>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class TRAILING_NEIGHBOR_ITERATOR<FLUIDS_SIMULATION_DATA<double>,double,2>;
#endif
