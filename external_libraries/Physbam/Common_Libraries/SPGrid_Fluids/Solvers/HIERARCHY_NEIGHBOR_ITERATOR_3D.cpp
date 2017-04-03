//#####################################################################
// Copyright 2013, Sean Bauer, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class HIERARCHY_NEIGHBOR_ITERATOR
//#####################################################################
#include <SPGrid_Fluids/Solvers/HIERARCHY_NEIGHBOR_ITERATOR.h>
#include <SPGrid_Fluids/Simulation/FLUIDS_SIMULATION_DATA.h>
#include <SPGrid_Fluids/Simulation/FLUIDS_SIMULATION_FLAGS.h>
#include <algorithm>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T_STRUCT, class T> HIERARCHY_NEIGHBOR_ITERATOR<T_STRUCT,T,3>::
HIERARCHY_NEIGHBOR_ITERATOR(const bool discard_dirichlet_neighbors_in,Hierarchy_type& hierarchy_in,CELL_ID cid_in,unsigned T_STRUCT::* flags_field_in,T T_STRUCT::* data_field_in,T T_STRUCT::* diag_field_in,T T_STRUCT::* Lx_field_in,T T_STRUCT::* Ly_field_in,T T_STRUCT::* Lz_field_in,T T_STRUCT::* Ux_field_in,T T_STRUCT::* Uy_field_in,T T_STRUCT::* Uz_field_in)
    :discard_dirichlet_neighbors(discard_dirichlet_neighbors_in),hierarchy(hierarchy_in)
{
    cid = cid_in;
    flags_field = flags_field_in;
    data_field = data_field_in;
    diag_field = diag_field_in;
    Lx_field = Lx_field_in;
    Ly_field = Ly_field_in;
    Lz_field = Lz_field_in;
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
// Constructor
//#####################################################################
template<class T_STRUCT, class T> HIERARCHY_NEIGHBOR_ITERATOR<T_STRUCT,T,3>::
HIERARCHY_NEIGHBOR_ITERATOR(const bool discard_dirichlet_neighbors_in,Hierarchy_type& hierarchy_in,CELL_ID cid_in,unsigned T_STRUCT::* flags_field_in,T T_STRUCT::* data_field_in,T T_STRUCT::* diag_field_in,VECTOR<T T_STRUCT::*,3> L_fields_in,VECTOR<T T_STRUCT::*,3> U_fields_in) 
    :discard_dirichlet_neighbors(discard_dirichlet_neighbors_in),hierarchy(hierarchy_in)
{
    cid = cid_in;
    flags_field = flags_field_in;
    data_field = data_field_in;
    diag_field = diag_field_in;
    Lx_field = L_fields_in(1);
    Ly_field = L_fields_in(2);
    Lz_field = L_fields_in(3);
    Ux_field = U_fields_in(1);
    Uy_field = U_fields_in(2);
    Uz_field = U_fields_in(3);

    // Initialize list of neighbors
    Enumerate_Neighbors();

    // Sort list by level, then offset
    Sort_Neighbors();

    index = 1; // index within neighbor array
}
//#####################################################################
// Function Enumerate_Neighbors
//#####################################################################
template<class T_STRUCT, class T> void HIERARCHY_NEIGHBOR_ITERATOR<T_STRUCT,T,3>::
Enumerate_Neighbors()
{
    std::stack<CHILD_CELL> cell_stack;
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
            coef_ptr = &hierarchy.Array(cid.level, Ux_field)(coef_offset);
        else
            coef_ptr = &hierarchy.Array(cid.level, Lx_field)(coef_offset);

        if(discard_dirichlet_neighbors && (hierarchy.Array(neighbor_level,flags_field)(neighbor_offset) & SPGrid_Cell_Type_Dirichlet));
        else
            neighbors.Append(NEIGHBOR(neighbor_level, neighbor_offset, data_ptr, coef_ptr));
    }
    else
        cell_stack.push(CHILD_CELL(CELL_ID(cid.level, cid.offset), SPGrid_Face_Minus_X_Active));


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
            coef_ptr = &hierarchy.Array(cid.level, Ux_field)(coef_offset);
        else
            coef_ptr = &hierarchy.Array(cid.level, Lx_field)(coef_offset);
        
        if(discard_dirichlet_neighbors && (hierarchy.Array(neighbor_level,flags_field)(neighbor_offset) & SPGrid_Cell_Type_Dirichlet));
        else
            neighbors.Append(NEIGHBOR(neighbor_level, neighbor_offset, data_ptr, coef_ptr));
    }
    else
        cell_stack.push(CHILD_CELL(CELL_ID(cid.level, cid.offset), SPGrid_Face_Plus_X_Active));
    
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
            coef_ptr = &hierarchy.Array(cid.level, Uy_field)(coef_offset);
        else
            coef_ptr = &hierarchy.Array(cid.level, Ly_field)(coef_offset);
         
        if(discard_dirichlet_neighbors && (hierarchy.Array(neighbor_level,flags_field)(neighbor_offset) & SPGrid_Cell_Type_Dirichlet));
        else
            neighbors.Append(NEIGHBOR(neighbor_level, neighbor_offset, data_ptr, coef_ptr));
    }
    else
        cell_stack.push(CHILD_CELL(CELL_ID(cid.level, cid.offset), SPGrid_Face_Minus_Y_Active));


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
            coef_ptr = &hierarchy.Array(cid.level, Uy_field)(coef_offset);
        else
            coef_ptr = &hierarchy.Array(cid.level, Ly_field)(coef_offset);
         
        if(discard_dirichlet_neighbors && (hierarchy.Array(neighbor_level,flags_field)(neighbor_offset) & SPGrid_Cell_Type_Dirichlet));
        else
            neighbors.Append(NEIGHBOR(neighbor_level, neighbor_offset, data_ptr, coef_ptr));
    }
    else
        cell_stack.push(CHILD_CELL(CELL_ID(cid.level, cid.offset), SPGrid_Face_Plus_Y_Active));
    
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
            coef_ptr = &hierarchy.Array(cid.level, Uz_field)(coef_offset);
        else
            coef_ptr = &hierarchy.Array(cid.level, Lz_field)(coef_offset);
         
        if(discard_dirichlet_neighbors && (hierarchy.Array(neighbor_level,flags_field)(neighbor_offset) & SPGrid_Cell_Type_Dirichlet));
        else
            neighbors.Append(NEIGHBOR(neighbor_level, neighbor_offset, data_ptr, coef_ptr));
    }
    else
        cell_stack.push(CHILD_CELL(CELL_ID(cid.level, cid.offset), SPGrid_Face_Minus_Z_Active));


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
            coef_ptr = &hierarchy.Array(cid.level, Uz_field)(coef_offset);
        else
            coef_ptr = &hierarchy.Array(cid.level, Lz_field)(coef_offset);
         
        if(discard_dirichlet_neighbors && (hierarchy.Array(neighbor_level,flags_field)(neighbor_offset) & SPGrid_Cell_Type_Dirichlet));
        else
            neighbors.Append(NEIGHBOR(neighbor_level, neighbor_offset, data_ptr, coef_ptr));
    }
    else
        cell_stack.push(CHILD_CELL(CELL_ID(cid.level, cid.offset), SPGrid_Face_Plus_Z_Active));

    // Add my own cell
    T* data_ptr = &hierarchy.Array(cid.level, data_field)(cid.offset);
    T* coef_ptr = &hierarchy.Array(cid.level, diag_field)(cid.offset);
     
    if(discard_dirichlet_neighbors && (hierarchy.Array(cid.level,flags_field)(cid.offset) & SPGrid_Cell_Type_Dirichlet));
    else
        neighbors.Append(NEIGHBOR(cid, data_ptr, coef_ptr));

    // Now process cell_stack & generate finer neighbors
    while(!cell_stack.empty())
    {
        //std::cout<<"grabbing element off cell stack"<<std::endl;
        CHILD_CELL child = cell_stack.top();
        cell_stack.pop();
        
        CELL_ID cid = child.x;
        unsigned mask = child.y;

        if(hierarchy.Set(cid.level).Is_Set(cid.offset,mask))
        {
            // Found active face, add neighbor to neighbors list
            unsigned long offset;
            T* data_ptr;
            T* coef_ptr;

            switch(mask)
            {
                case SPGrid_Face_Minus_X_Active:
                    offset = Flag_array_mask::template Packed_Offset<-1, 0, 0>(cid.offset);
                    data_ptr = &hierarchy.Array(cid.level, data_field)(offset);
                    coef_ptr = &hierarchy.Array(cid.level, Lx_field)(cid.offset);
                    break;
                case SPGrid_Face_Plus_X_Active:
                    offset = Flag_array_mask::template Packed_Offset< 1, 0, 0>(cid.offset);
                    data_ptr = &hierarchy.Array(cid.level, data_field)(offset);
                    coef_ptr = &hierarchy.Array(cid.level, Lx_field)(offset);
                    break;
                case SPGrid_Face_Minus_Y_Active:
                    offset = Flag_array_mask::template Packed_Offset< 0,-1, 0>(cid.offset);
                    data_ptr = &hierarchy.Array(cid.level, data_field)(offset);
                    coef_ptr = &hierarchy.Array(cid.level, Ly_field)(cid.offset);
                    break;
                case SPGrid_Face_Plus_Y_Active:
                    offset = Flag_array_mask::template Packed_Offset< 0, 1, 0>(cid.offset);
                    data_ptr = &hierarchy.Array(cid.level, data_field)(offset);
                    coef_ptr = &hierarchy.Array(cid.level, Ly_field)(offset);
                    break;
                case SPGrid_Face_Minus_Z_Active:
                    offset = Flag_array_mask::template Packed_Offset< 0, 0,-1>(cid.offset);
                    data_ptr = &hierarchy.Array(cid.level, data_field)(offset);
                    coef_ptr = &hierarchy.Array(cid.level, Lz_field)(cid.offset);
                    break;
                case SPGrid_Face_Plus_Z_Active:
                    offset = Flag_array_mask::template Packed_Offset< 0, 0, 1>(cid.offset);
                    data_ptr = &hierarchy.Array(cid.level, data_field)(offset);
                    coef_ptr = &hierarchy.Array(cid.level, Lz_field)(offset);
                    break;
                default:
                    std::cout<<"Reached default case of cell-splitter"<<std::endl;
                    exit(1);
            }
             
            if(discard_dirichlet_neighbors && (hierarchy.Array(cid.level,flags_field)(offset) & SPGrid_Cell_Type_Dirichlet));
            else
                neighbors.Append(NEIGHBOR(cid.level, offset, data_ptr, coef_ptr));
        }
        else
        {
            //std::cout<<"Splitting ghost cell"<<std::endl;
            int level = cid.level-1;
            if(level==0) continue;
            unsigned long base_offset = Flag_array_mask::UpsampleOffset(cid.offset);
            unsigned long offsetA, offsetB, offsetC, offsetD;

            // Split up cell, add two children to stack
            switch(mask)
            {
                case SPGrid_Face_Minus_X_Active:
                    offsetA = base_offset;
                    offsetB = Flag_array_mask::template Packed_Offset<0,1,0>(base_offset);
                    offsetC = Flag_array_mask::template Packed_Offset<0,0,1>(base_offset);
                    offsetD = Flag_array_mask::template Packed_Offset<0,1,1>(base_offset);
                    break;
                case SPGrid_Face_Plus_X_Active:
                    offsetA = Flag_array_mask::template Packed_Offset<1,0,0>(base_offset);
                    offsetB = Flag_array_mask::template Packed_Offset<1,1,0>(base_offset);
                    offsetC = Flag_array_mask::template Packed_Offset<1,0,1>(base_offset);
                    offsetD = Flag_array_mask::template Packed_Offset<1,1,1>(base_offset);
                    break;
                case SPGrid_Face_Minus_Y_Active:
                    offsetA = base_offset;
                    offsetB = Flag_array_mask::template Packed_Offset<0,0,1>(base_offset);
                    offsetC = Flag_array_mask::template Packed_Offset<1,0,0>(base_offset);
                    offsetD = Flag_array_mask::template Packed_Offset<1,0,1>(base_offset);
                    break;
                case SPGrid_Face_Plus_Y_Active:
                    offsetA = Flag_array_mask::template Packed_Offset<0,1,0>(base_offset);
                    offsetB = Flag_array_mask::template Packed_Offset<1,1,0>(base_offset);
                    offsetC = Flag_array_mask::template Packed_Offset<0,1,1>(base_offset);
                    offsetD = Flag_array_mask::template Packed_Offset<1,1,1>(base_offset);
                    break;
                case SPGrid_Face_Minus_Z_Active:
                    offsetA = base_offset;
                    offsetB = Flag_array_mask::template Packed_Offset<0,1,0>(base_offset);
                    offsetC = Flag_array_mask::template Packed_Offset<1,0,0>(base_offset);
                    offsetD = Flag_array_mask::template Packed_Offset<1,1,0>(base_offset);
                    break;
                case SPGrid_Face_Plus_Z_Active:
                    offsetA = Flag_array_mask::template Packed_Offset<0,0,1>(base_offset);
                    offsetB = Flag_array_mask::template Packed_Offset<0,1,1>(base_offset);
                    offsetC = Flag_array_mask::template Packed_Offset<1,0,1>(base_offset);
                    offsetD = Flag_array_mask::template Packed_Offset<1,1,1>(base_offset);
                    break;
                default:
                    std::cout<<"Reached default case of cell-splitter"<<std::endl;
                    exit(1);
            }

            if(hierarchy.Set(level).Is_Set(offsetA,SPGrid_Cell_Type_Ghost)) cell_stack.push(CHILD_CELL(CELL_ID(level,offsetA),mask));
            if(hierarchy.Set(level).Is_Set(offsetB,SPGrid_Cell_Type_Ghost)) cell_stack.push(CHILD_CELL(CELL_ID(level,offsetB),mask));
            if(hierarchy.Set(level).Is_Set(offsetC,SPGrid_Cell_Type_Ghost)) cell_stack.push(CHILD_CELL(CELL_ID(level,offsetC),mask));
            if(hierarchy.Set(level).Is_Set(offsetD,SPGrid_Cell_Type_Ghost)) cell_stack.push(CHILD_CELL(CELL_ID(level,offsetD),mask));
        }
    }
}
//#####################################################################
// Function Sort_Neighbors
//#####################################################################
template<class T_STRUCT, class T> void HIERARCHY_NEIGHBOR_ITERATOR<T_STRUCT,T,3>::
Sort_Neighbors()
{
    // Implement sort by level, then by offset within level
    std::sort(neighbors.begin(), neighbors.end());
}
//#####################################################################
// Function Coarse_Cell_Finder
//#####################################################################
template<class T_STRUCT, class T> void HIERARCHY_NEIGHBOR_ITERATOR<T_STRUCT,T,3>::
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
template class HIERARCHY_NEIGHBOR_ITERATOR<FLUIDS_SIMULATION_DATA<float>,float,3>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class HIERARCHY_NEIGHBOR_ITERATOR<FLUIDS_SIMULATION_DATA<double>,double,3>;
#endif
