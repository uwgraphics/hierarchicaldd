//#####################################################################
// Copyright 2013, Sean Bauer, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class TRAILING_NGN_ITERATOR
//#####################################################################
#include <SPGrid_Fluids/Solvers/TRAILING_NGN_ITERATOR.h>
#include <SPGrid_Fluids/Simulation/FLUIDS_SIMULATION_DATA.h>
#include <SPGrid_Fluids/Simulation/FLUIDS_SIMULATION_FLAGS.h>
#include <algorithm>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T_STRUCT, class T> TRAILING_NGN_ITERATOR<T_STRUCT,T,3>::
TRAILING_NGN_ITERATOR(const bool discard_dirichlet_neighbors_in,Hierarchy_type& hierarchy_in,CELL_ID cid_in,unsigned T_STRUCT::* flags_field_in,T T_STRUCT::* data_field_in,T T_STRUCT::* diag_field_in,T T_STRUCT::* Ux_field_in,T T_STRUCT::* Uy_field_in,T T_STRUCT::* Uz_field_in)
    :discard_dirichlet_neighbors(discard_dirichlet_neighbors_in),hierarchy(hierarchy_in)
{
    cid = cid_in;
    flags_field = flags_field_in;
    data_field = data_field_in;
    diag_field = diag_field_in;
    Ux_field = Ux_field_in;
    Uy_field = Uy_field_in;
    Uz_field = Uz_field_in;

    // Initialize list of neighbors
    Enumerate_Neighbors();

    // Sort list by level, then offset
    Sort_Neighbors();

    index = 1; // index within neighbor array
}
//#####################################################################
// Function Enumerate_Neighbors
//#####################################################################
template<class T_STRUCT, class T> void TRAILING_NGN_ITERATOR<T_STRUCT,T,3>::
Enumerate_Neighbors()
{
    // Add my own cell
    T* data_ptr = &hierarchy.Array(cid.level, data_field)(cid.offset);
    T* coef_ptr = &hierarchy.Array(cid.level, diag_field)(cid.offset);
    if(discard_dirichlet_neighbors && (hierarchy.Array(cid.level,flags_field)(cid.offset) & SPGrid_Cell_Type_Dirichlet));
    else
        neighbors.Append(NEIGHBOR(cid, data_ptr, coef_ptr));

    unsigned flags = hierarchy.Array(cid.level, flags_field)(cid.offset);

    // Find same or coarser neighbors first
    if(flags & SPGrid_Face_Minus_X_Active)
    {
        int neighbor_level = cid.level;
        unsigned long neighbor_offset = Flag_array_mask::template Packed_OffsetXdim<-1>(cid.offset);
        unsigned long coef_offset = cid.offset;
        Coarse_Cell_Finder(neighbor_level, neighbor_offset);

        T* data_ptr = &hierarchy.Array(cid.level, data_field)(neighbor_offset);
        T* coef_ptr;
        CELL_ID nid(neighbor_level,neighbor_offset);
        if(cid < nid)
        {
            coef_ptr = &hierarchy.Array(cid.level, Ux_field)(coef_offset);
            if(discard_dirichlet_neighbors && (hierarchy.Array(neighbor_level,flags_field)(neighbor_offset) & SPGrid_Cell_Type_Dirichlet));
            else
                neighbors.Append(NEIGHBOR(neighbor_level, neighbor_offset, data_ptr, coef_ptr));
        }
    }
    if(flags & SPGrid_Face_Plus_X_Active)
    {
        int neighbor_level = cid.level;
        unsigned long neighbor_offset = Flag_array_mask::template Packed_OffsetXdim<1>(cid.offset);
        unsigned long coef_offset = neighbor_offset;
        Coarse_Cell_Finder(neighbor_level, neighbor_offset);

        T* data_ptr = &hierarchy.Array(cid.level, data_field)(neighbor_offset);
        T* coef_ptr;
        CELL_ID nid(neighbor_level,neighbor_offset);
        if(cid < nid)
        {
            coef_ptr = &hierarchy.Array(cid.level, Ux_field)(coef_offset);
            if(discard_dirichlet_neighbors && (hierarchy.Array(neighbor_level,flags_field)(neighbor_offset) & SPGrid_Cell_Type_Dirichlet));
            else
                neighbors.Append(NEIGHBOR(neighbor_level, neighbor_offset, data_ptr, coef_ptr));
        }
    }
    
    if(flags & SPGrid_Face_Minus_Y_Active)
    {
        int neighbor_level = cid.level;
        unsigned long neighbor_offset = Flag_array_mask::template Packed_OffsetYdim<-1>(cid.offset);
        unsigned long coef_offset = cid.offset;
        Coarse_Cell_Finder(neighbor_level, neighbor_offset);

        T* data_ptr = &hierarchy.Array(cid.level, data_field)(neighbor_offset);
        T* coef_ptr;
        CELL_ID nid(neighbor_level,neighbor_offset);
        if(cid < nid)
        {
            coef_ptr = &hierarchy.Array(cid.level, Uy_field)(coef_offset);
            if(discard_dirichlet_neighbors && (hierarchy.Array(neighbor_level,flags_field)(neighbor_offset) & SPGrid_Cell_Type_Dirichlet));
            else
                neighbors.Append(NEIGHBOR(neighbor_level, neighbor_offset, data_ptr, coef_ptr));
        }
    }
    if(flags & SPGrid_Face_Plus_Y_Active)
    {
        int neighbor_level = cid.level;
        unsigned long neighbor_offset = Flag_array_mask::template Packed_OffsetYdim<1>(cid.offset);
        unsigned long coef_offset = neighbor_offset;
        Coarse_Cell_Finder(neighbor_level, neighbor_offset);

        T* data_ptr = &hierarchy.Array(cid.level, data_field)(neighbor_offset);
        T* coef_ptr;
        CELL_ID nid(neighbor_level,neighbor_offset);
        if(cid < nid)
        {
            coef_ptr = &hierarchy.Array(cid.level, Uy_field)(coef_offset);
            if(discard_dirichlet_neighbors && (hierarchy.Array(neighbor_level,flags_field)(neighbor_offset) & SPGrid_Cell_Type_Dirichlet));
            else
                neighbors.Append(NEIGHBOR(neighbor_level, neighbor_offset, data_ptr, coef_ptr));
        }
    }
    if(flags & SPGrid_Face_Minus_Z_Active)
    {
        int neighbor_level = cid.level;
        unsigned long neighbor_offset = Flag_array_mask::template Packed_OffsetZdim<-1>(cid.offset);
        unsigned long coef_offset = cid.offset;
        Coarse_Cell_Finder(neighbor_level, neighbor_offset);

        T* data_ptr = &hierarchy.Array(cid.level, data_field)(neighbor_offset);
        T* coef_ptr;
        CELL_ID nid(neighbor_level,neighbor_offset);
        if(cid < nid)
        {
            coef_ptr = &hierarchy.Array(cid.level, Uz_field)(coef_offset);
            if(discard_dirichlet_neighbors && (hierarchy.Array(neighbor_level,flags_field)(neighbor_offset) & SPGrid_Cell_Type_Dirichlet));
            else
                neighbors.Append(NEIGHBOR(neighbor_level, neighbor_offset, data_ptr, coef_ptr));
        }
    }
    if(flags & SPGrid_Face_Plus_Z_Active)
    {
        int neighbor_level = cid.level;
        unsigned long neighbor_offset = Flag_array_mask::template Packed_OffsetZdim<1>(cid.offset);
        unsigned long coef_offset = neighbor_offset;
        Coarse_Cell_Finder(neighbor_level, neighbor_offset);

        T* data_ptr = &hierarchy.Array(cid.level, data_field)(neighbor_offset);
        T* coef_ptr;
        CELL_ID nid(neighbor_level,neighbor_offset);
        if(cid < nid)
        {
            coef_ptr = &hierarchy.Array(cid.level, Uz_field)(coef_offset);
            if(discard_dirichlet_neighbors && (hierarchy.Array(neighbor_level,flags_field)(neighbor_offset) & SPGrid_Cell_Type_Dirichlet));
            else
                neighbors.Append(NEIGHBOR(neighbor_level, neighbor_offset, data_ptr, coef_ptr));
        }
    }
}
//#####################################################################
// Function Sort_Neighbors
//#####################################################################
template<class T_STRUCT, class T> void TRAILING_NGN_ITERATOR<T_STRUCT,T,3>::
Sort_Neighbors()
{
    // Implement sort by level, then by offset within level
    std::sort(neighbors.begin(), neighbors.end());
}
//#####################################################################
// Function Coarse_Cell_Finder
//#####################################################################
template<class T_STRUCT, class T> void TRAILING_NGN_ITERATOR<T_STRUCT,T,3>::
Coarse_Cell_Finder(int& neighbor_level, unsigned long& neighbor_offset)
{
    // Also check for ghost
    while(neighbor_level<=hierarchy.Levels())
    {
        unsigned flags = hierarchy.Array(neighbor_level, flags_field)(neighbor_offset);
        if(flags & (SPGrid_Cell_Type_Interior|SPGrid_Cell_Type_Dirichlet))
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
template class TRAILING_NGN_ITERATOR<FLUIDS_SIMULATION_DATA<float>,float,3>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class TRAILING_NGN_ITERATOR<FLUIDS_SIMULATION_DATA<double>,double,3>;
#endif
