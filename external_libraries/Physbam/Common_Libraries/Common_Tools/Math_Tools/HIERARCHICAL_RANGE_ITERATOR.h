//#####################################################################
// Copyright 2012, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class HIERARCHICAL_RANGE_ITERATOR
//#####################################################################
#ifndef __HIERARCHICAL_RANGE_ITERATOR__
#define __HIERARCHICAL_RANGE_ITERATOR__

#include <PhysBAM_Tools/Data_Structures/STACK.h>
#include <PhysBAM_Tools/Math_Tools/RANGE.h>
#include <Common_Tools/Log/ADAPTIVE_PROGRESS_INDICATOR.h>

namespace PhysBAM{

template<int d,class RANGE_FUNCTOR>
class HIERARCHICAL_RANGE_ITERATOR
{
    typedef VECTOR<int,d> T_INDEX;
    typedef RANGE<T_INDEX> T_RANGE;

    STACK<T_RANGE> stack;
    ADAPTIVE_PROGRESS_INDICATOR indicator;
    RANGE_FUNCTOR& functor;

public:

    HIERARCHICAL_RANGE_ITERATOR(const T_RANGE range,RANGE_FUNCTOR& functor_input)
        :indicator(VECTOR<unsigned long,d>(range.Edge_Lengths()).Product(),0x100000),functor(functor_input)
    {stack.Push(range);}

    bool Valid() const
    {return !stack.Empty();}

    void Next()
    {
        T_RANGE range=stack.Pop();
        bool recurse=functor.Consume(range);
        if(!recurse) {indicator.Progress(VECTOR<unsigned long,d>(range.Edge_Lengths()).Product());return;}
        T_INDEX edge_lengths=range.Edge_Lengths();
        // The order of splitting here matters based on size padding
        int split_axis=1;for(;edge_lengths(split_axis)<edge_lengths(d);split_axis++);
        T_RANGE lower_range(range),upper_range(range);
        int split=(range.max_corner(split_axis)+range.min_corner(split_axis))/2;
        lower_range.max_corner(split_axis)=split;upper_range.min_corner(split_axis)=split;
        stack.Push(lower_range);stack.Push(upper_range);
    }

//#####################################################################
};
}
#endif
