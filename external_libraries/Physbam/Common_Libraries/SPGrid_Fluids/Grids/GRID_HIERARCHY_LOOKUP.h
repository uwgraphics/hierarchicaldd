//#####################################################################
// Copyright 2013, Raj Setaluri, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __GRID_HIERARCHY_LOOKUP_h__
#define __GRID_HIERARCHY_LOOKUP_h__

#include <SPGrid/Core/SPGrid_Mask.h>
#include <SPGrid_Fluids/Grids/GRID_HIERARCHY.h>

namespace PhysBAM{
using namespace SPGrid;
template<class T_STRUCT,class T,int d>
//#####################################################################
// Class GRID_HIERARCHY_LOOKUP
//#####################################################################
class GRID_HIERARCHY_LOOKUP
{
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<unsigned>::type::MASK Flag_array_mask;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<unsigned>::type Flag_array_type;
    typedef GRID_HIERARCHY<T_STRUCT,T,d> Hierarchy_type;
    typedef VECTOR<T,d> TV;
    typedef VECTOR<T,d-1> TV2;
    typedef VECTOR<int,d> TV_INT;

    template<class T_MASK> struct Cell_Parity_Helper {};        
    template<int log2_struct,int log2_field>
    struct Cell_Parity_Helper<SPGrid_Mask<log2_struct,log2_field,2> >
    {
        enum{
            x_mask=SPGrid_Mask<log2_struct,log2_field,2>::template LinearOffset<1,0>::value,
            y_mask=SPGrid_Mask<log2_struct,log2_field,2>::template LinearOffset<0,1>::value,
            z_mask=0UL
        };
    };    
    template<int log2_struct,int log2_field>
    struct Cell_Parity_Helper<SPGrid_Mask<log2_struct,log2_field,3> >
    {
        enum{
            x_mask=SPGrid_Mask<log2_struct,log2_field,3>::template LinearOffset<1,0,0>::value,
            y_mask=SPGrid_Mask<log2_struct,log2_field,3>::template LinearOffset<0,1,0>::value,
            z_mask=SPGrid_Mask<log2_struct,log2_field,3>::template LinearOffset<0,0,1>::value
        };
    };

// #############################################################################
public:
    static bool Cell_Lookup(const Hierarchy_type& hierarchy,const TV& X,unsigned long &flags_offset,int &level);
    static bool Cell_Lookup(const Hierarchy_type& hierarchy,const TV& X,unsigned long &flags_offset,int &level,TV& weights);
    static bool Cell_Lookup(const Hierarchy_type& hierarchy,const TV& X,unsigned long &flags_offset,TV_INT &index,int &level,TV& weights);
    static bool Cell_Lookup(const Hierarchy_type& hierarchy,unsigned long& flags_offset,int& level,TV& weights);
    static bool Cell_Lookup_Coarser(const Hierarchy_type& hierarchy,unsigned long& flags_offset,int& level,TV& weights);
    static bool Cell_Lookup_Finer(const Hierarchy_type& hierarchy,unsigned long& flags_offset,int& level,TV& weights);
    static bool Face_Lookup(const Hierarchy_type& hierarchy,unsigned long& flags_offset,int& level,TV2& weights,const unsigned mask,const int axis);
    static bool Face_Lookup_Finer(const Hierarchy_type& hierarchy,unsigned long& flags_offset,int& level,TV2& weights,const unsigned mask,const int axis);
// #############################################################################
};
}
#endif
