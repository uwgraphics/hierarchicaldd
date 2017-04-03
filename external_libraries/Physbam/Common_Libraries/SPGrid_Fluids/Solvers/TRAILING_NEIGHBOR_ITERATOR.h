//#####################################################################
// Copyright 2013, Sean Bauer, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class TRAILING_NEIGHBOR_ITERATOR
//#####################################################################
#ifndef __TRAILING_NEIGHBOR_ITERATOR__
#define __TRAILING_NEIGHBOR_ITERATOR__

#include <stack>
#include <SPGrid_Fluids/Grids/GRID_HIERARCHY.h>
#include "Parity_Helper.h"
#include "NEIGHBOR_STRUCT.h"
#include <PhysBAM_Tools/Data_Structures/PAIR.h>

namespace PhysBAM{

template<class T_STRUCT, class T, int d> class TRAILING_NEIGHBOR_ITERATOR;

typedef PAIR<CELL_ID,unsigned> CHILD_CELL;

template<class T_STRUCT, class T>
class TRAILING_NEIGHBOR_ITERATOR<T_STRUCT,T,2>
{
    enum{d=2};
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<unsigned>::type Flag_array_type;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<T>::type Data_array_type;
    typedef typename Flag_array_type::MASK Flag_array_mask;
    typedef typename Data_array_type::MASK Data_array_mask;
    typedef GRID_HIERARCHY<T_STRUCT,T,d> Hierarchy_type;
    typedef Parity_Helper<Flag_array_mask,d> Parity;
    typedef NEIGHBOR_STRUCT<T> NEIGHBOR;

private:
    CELL_ID cid;
    ARRAY<NEIGHBOR> neighbors;
    int index;
    Hierarchy_type& hierarchy;
    unsigned T_STRUCT::* flags_field;
    T T_STRUCT::* data_field;
    T T_STRUCT::* diag_field;
    T T_STRUCT::* Ux_field;
    T T_STRUCT::* Uy_field;

public:
    TRAILING_NEIGHBOR_ITERATOR(Hierarchy_type& hierarchy_in,CELL_ID cid_in,unsigned T_STRUCT::* flags_field_in,T T_STRUCT::* data_field_in,T T_STRUCT::* diag_field_in,T T_STRUCT::* Ux_field_in,T T_STRUCT::* Uy_field_in); 
    
    // Wrong dim constructor
    TRAILING_NEIGHBOR_ITERATOR(Hierarchy_type& h,CELL_ID c,unsigned T_STRUCT::* f,T T_STRUCT::* d1,T T_STRUCT::* d2,T T_STRUCT::* d3,T T_STRUCT::* d4,T T_STRUCT::* d5)
        :hierarchy(h)
    {PHYSBAM_FATAL_ERROR("Dimension mismatch in TRAILING_NEIGHBOR_ITERATOR");}
    
    void Enumerate_Neighbors();

    void Coarse_Cell_Finder(int& neighbor_level, unsigned long& neighbor_offset);

    inline NEIGHBOR& Neighbor() {return neighbors(index);}

    inline T& Coefficient() {return *neighbors(index).coefficient;}
    inline T& Data()        {return *neighbors(index).data;}
    inline void Next()      {index++;}
    inline bool Valid()     {return index<=neighbors.Size();}
    inline int Size()       {return neighbors.Size();}
};

template<class T_STRUCT, class T>
class TRAILING_NEIGHBOR_ITERATOR<T_STRUCT,T,3>
{
    enum{d=3};
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<unsigned>::type Flag_array_type;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<T>::type Data_array_type;
    typedef typename Flag_array_type::MASK Flag_array_mask;
    typedef typename Data_array_type::MASK Data_array_mask;
    typedef GRID_HIERARCHY<T_STRUCT,T,d> Hierarchy_type;
    typedef Parity_Helper<Flag_array_mask,d> Parity;
    typedef NEIGHBOR_STRUCT<T> NEIGHBOR;

private:
    CELL_ID cid;
    ARRAY<NEIGHBOR> neighbors;
    int index;
    Hierarchy_type& hierarchy;
    unsigned T_STRUCT::* flags_field;
    T T_STRUCT::* data_field;
    T T_STRUCT::* diag_field;
    T T_STRUCT::* Ux_field;
    T T_STRUCT::* Uy_field;
    T T_STRUCT::* Uz_field;

public:
    TRAILING_NEIGHBOR_ITERATOR(Hierarchy_type& hierarchy_in,CELL_ID cid_in,unsigned T_STRUCT::* flags_field_in,T T_STRUCT::* data_field_in,T T_STRUCT::* diag_field_in,T T_STRUCT::* Ux_field_in,T T_STRUCT::* Uy_field_in,T T_STRUCT::* Uz_field_in);
    
    // Wrong dim constructor
    TRAILING_NEIGHBOR_ITERATOR(Hierarchy_type& h,CELL_ID c,unsigned T_STRUCT::* f,T T_STRUCT::* d1,T T_STRUCT::* d2,T T_STRUCT::* d3,T T_STRUCT::* d4)
        :hierarchy(h)
    {PHYSBAM_FATAL_ERROR("Dimension mismatch in TRAILING_NEIGHBOR_ITERATOR");}
    
    void Enumerate_Neighbors();

    void Coarse_Cell_Finder(int& neighbor_level,unsigned long& neighbor_offset);

    inline NEIGHBOR& Neighbor() {return neighbors(index);}

    inline T& Coefficient() {return *neighbors(index).coefficient;}
    inline T& Data()        {return *neighbors(index).data;}
    inline void Next()      {index++;}
    inline bool Valid()     {return index<=neighbors.Size();}
    inline int Size()       {return neighbors.Size();}
};
}
#endif
