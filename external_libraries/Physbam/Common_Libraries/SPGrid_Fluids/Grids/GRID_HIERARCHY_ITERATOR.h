//#####################################################################
// Copyright 2012, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class GRID_HIERARCHY_ITERATOR
//#####################################################################
#ifndef __GRID_HIERARCHY_ITERATOR__
#define __GRID_HIERARCHY_ITERATOR__

#include <PhysBAM_Tools/Data_Structures/PAIR.h>
#include <PhysBAM_Tools/Data_Structures/STACK.h>
#include <PhysBAM_Tools/Math_Tools/RANGE.h>
#include <Common_Tools/Log/ADAPTIVE_PROGRESS_INDICATOR.h>
#include <math.h>

namespace PhysBAM{

template<int d,class CELL_FUNCTOR>
class GRID_HIERARCHY_ITERATOR
{
    typedef VECTOR<int,d> T_INDEX;
    typedef RANGE<T_INDEX> T_RANGE;
    typedef PAIR<int,T_INDEX> T_CELL; // (level,cell_index)

    STACK<T_CELL> stack;
    ADAPTIVE_PROGRESS_INDICATOR indicator;
    CELL_FUNCTOR& functor;

public:
    GRID_HIERARCHY_ITERATOR(const T_RANGE& range,const int level,CELL_FUNCTOR& functor_input);

    bool Valid() const
    {return !stack.Empty();}

    void Next();
//#####################################################################
};
}
#endif
