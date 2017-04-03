//#####################################################################
// Copyright 2014, Raj Setaluri.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class CONSTRAIN_T_JUNCTION_NODES
//#####################################################################
#ifndef __CONSTRAIN_T_JUNCTION_NODES__
#define __CONSTRAIN_T_JUNCTION_NODES__

#include <SPGrid/Core/SPGrid_Allocator.h>
#include <SPGrid_Fluids/Grids/GRID_HIERARCHY.h>

using namespace SPGrid;

namespace PhysBAM{

template<class T_MASK,int d> struct CONSTRAIN_T_JUNCTION_NODES_HELPER;

template<class T_MASK>
struct CONSTRAIN_T_JUNCTION_NODES_HELPER<T_MASK,2>
{
    enum{d=2};
    enum{odd_bits=T_MASK::template LinearOffset<1,1>::value,
        Parity_00=T_MASK::template LinearOffset<0,0>::value,
        Parity_01=T_MASK::template LinearOffset<0,1>::value,
        Parity_10=T_MASK::template LinearOffset<1,0>::value,
        Parity_11=T_MASK::template LinearOffset<1,1>::value};
    enum{
        plus_x  = T_MASK::template LinearOffset< 1, 0>::value,
        plus_y  = T_MASK::template LinearOffset< 0, 1>::value,
        minus_x = T_MASK::template LinearOffset<-1, 0>::value,
        minus_y = T_MASK::template LinearOffset< 0,-1>::value};
    static inline void Find_Neighboring_Nodes(const unsigned long offset,unsigned long* neighbor_offsets,int& m)
    {
        const unsigned long parity=(offset&odd_bits);
        const unsigned long coarse_offset=T_MASK::DownsampleOffset(offset);
        if(parity==Parity_11){
            m=1;
            neighbor_offsets[0]=coarse_offset;}
        else if(parity==Parity_10){
            m=2;
            neighbor_offsets[0]=coarse_offset;
            neighbor_offsets[1]=T_MASK::Packed_Add(coarse_offset,plus_y);}
        else if(parity==Parity_01){
            m=2;
            neighbor_offsets[0]=coarse_offset;
            neighbor_offsets[1]=T_MASK::Packed_Add(coarse_offset,plus_x);}
        else if(parity==Parity_00) PHYSBAM_FATAL_ERROR("Not valid arrangement for T-Junction node!");
        else PHYSBAM_FATAL_ERROR("Shouldn't reach this point.");
    }
};

template<class T_MASK>
struct CONSTRAIN_T_JUNCTION_NODES_HELPER<T_MASK,3>
{
    enum{d=3};
    enum{odd_bits=T_MASK::template LinearOffset<1,1,1>::value,
         Parity_000=T_MASK::template LinearOffset<0,0,0>::value,
         Parity_001=T_MASK::template LinearOffset<0,0,1>::value,
         Parity_010=T_MASK::template LinearOffset<0,1,0>::value,
         Parity_011=T_MASK::template LinearOffset<0,1,1>::value,
         Parity_100=T_MASK::template LinearOffset<1,0,0>::value,
         Parity_101=T_MASK::template LinearOffset<1,0,1>::value,
         Parity_110=T_MASK::template LinearOffset<1,1,0>::value,
         Parity_111=T_MASK::template LinearOffset<1,1,1>::value};
    enum{
        plus_x        = T_MASK::template LinearOffset< 1, 0, 0>::value,
        plus_y        = T_MASK::template LinearOffset< 0, 1, 0>::value,
        plus_z        = T_MASK::template LinearOffset< 0, 0, 1>::value,
        plus_x_and_y  = T_MASK::template LinearOffset< 1, 1, 0>::value,
        plus_x_and_z  = T_MASK::template LinearOffset< 1, 0, 1>::value,
        plus_y_and_z  = T_MASK::template LinearOffset< 0, 1, 1>::value,
        minus_x       = T_MASK::template LinearOffset<-1, 0, 0>::value,
        minus_y       = T_MASK::template LinearOffset< 0,-1, 0>::value,
        minus_z       = T_MASK::template LinearOffset< 0, 0,-1>::value,
        minus_x_and_y = T_MASK::template LinearOffset<-1,-1, 0>::value,
        minus_x_and_z = T_MASK::template LinearOffset<-1, 0,-1>::value,
        minus_y_and_z = T_MASK::template LinearOffset< 0,-1,-1>::value};
    static inline void Find_Neighboring_Nodes(const unsigned long offset,unsigned long* neighbor_offsets,int& m)
    {
        const unsigned long parity=(offset&odd_bits);
        const unsigned long coarse_offset=T_MASK::DownsampleOffset(offset);
        if(parity==Parity_111){
            m=1;
            neighbor_offsets[0]=coarse_offset;}
        else if(parity==Parity_110){
            m=2;
            neighbor_offsets[0]=coarse_offset;
            neighbor_offsets[1]=T_MASK::Packed_Add(coarse_offset,plus_z);}
        else if(parity==Parity_101){ 
            m=2;
            neighbor_offsets[0]=coarse_offset;
            neighbor_offsets[1]=T_MASK::Packed_Add(coarse_offset,plus_y);}
        else if(parity==Parity_011){
            m=2;
            neighbor_offsets[0]=coarse_offset;
            neighbor_offsets[1]=T_MASK::Packed_Add(coarse_offset,plus_x);}
        else if(parity==Parity_100){
            m=4;
            neighbor_offsets[0]=coarse_offset;
            neighbor_offsets[1]=T_MASK::Packed_Add(coarse_offset,plus_y);
            neighbor_offsets[1]=T_MASK::Packed_Add(coarse_offset,plus_z);
            neighbor_offsets[1]=T_MASK::Packed_Add(coarse_offset,plus_y_and_z);}
        else if(parity==Parity_010){
            m=4;
            neighbor_offsets[0]=coarse_offset;
            neighbor_offsets[1]=T_MASK::Packed_Add(coarse_offset,plus_x);
            neighbor_offsets[1]=T_MASK::Packed_Add(coarse_offset,plus_z);
            neighbor_offsets[1]=T_MASK::Packed_Add(coarse_offset,plus_x_and_z);}
        else if(parity==Parity_001){
            m=4;
            neighbor_offsets[0]=coarse_offset;
            neighbor_offsets[1]=T_MASK::Packed_Add(coarse_offset,plus_x);
            neighbor_offsets[1]=T_MASK::Packed_Add(coarse_offset,plus_y);
            neighbor_offsets[1]=T_MASK::Packed_Add(coarse_offset,plus_x_and_y);}
        else if(parity==Parity_000) PHYSBAM_FATAL_ERROR("Not valid arrangement for T-Junction node!");
        else PHYSBAM_FATAL_ERROR("Shouldn't reach this point.");
    }
};

template<class T_STRUCT,class T_DATA,class T_FLAGS,int d>
class CONSTRAIN_T_JUNCTION_NODES
{
    typedef T_DATA T;
    typedef GRID_HIERARCHY<T_STRUCT,T,d> T_HIERARCHY;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<unsigned>::type Flag_array_type;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<const T_FLAGS>::type Const_flag_array_type;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<T>::type Data_array_type;
    typedef typename Flag_array_type::MASK Flag_array_mask;

    static const int max_neighbors=(d==2)?2:4;

    T_HIERARCHY& hierarchy;
    T T_STRUCT::* node_field;
    T_FLAGS T_STRUCT::* flags_field;
    const int level;
    unsigned long* neighbor_offsets;

public:
    CONSTRAIN_T_JUNCTION_NODES(T_HIERARCHY& hierarchy_input,T T_STRUCT::* node_field_input,T_FLAGS T_STRUCT::* flags_field_input,const int level_input)
        :hierarchy(hierarchy_input),node_field(node_field_input),flags_field(flags_field_input),level(level_input),neighbor_offsets(new unsigned long[max_neighbors])
    {}

    ~CONSTRAIN_T_JUNCTION_NODES()
    {delete neighbor_offsets;}

    CONSTRAIN_T_JUNCTION_NODES(SPGrid_Allocator<T_STRUCT,d>& allocator,const std::pair<const unsigned long*,unsigned>& blocks,
        T_HIERARCHY& hierarchy_input,T T_STRUCT::* node_field_input,T_FLAGS T_STRUCT::* flags_field_input,const int level_input)
        :hierarchy(hierarchy_input),node_field(node_field_input),flags_field(flags_field_input),level(level_input),neighbor_offsets(new unsigned long[max_neighbors])
    {Run(allocator,blocks);}

    void Run(SPGrid_Allocator<T_STRUCT,d>& allocator,const std::pair<const unsigned long*,unsigned>& blocks) const
    {Data_array_type node_array=allocator.Get_Array(node_field);
    Const_flag_array_type flags=allocator.Get_Const_Array(flags_field);
    for(SPGrid_Block_Iterator<typename Data_array_type::MASK> iterator(blocks);iterator.Valid();iterator.Next())
        if(iterator.Data(flags)&SPGrid_Node_T_Junction){
            int m=0;
            CONSTRAIN_T_JUNCTION_NODES_HELPER<Flag_array_mask,d>::Find_Neighboring_Nodes(iterator.Offset(),neighbor_offsets,m);
            double value=(double)0;
            for(int node=0;node<m;node++) value+=(double)hierarchy.Array(level+1,node_field)(neighbor_offsets[node]);
            value/=(double)m;
            iterator.Data(node_array)=value;}}
//#####################################################################
};
}
#endif
