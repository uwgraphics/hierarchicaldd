//#####################################################################
// Copyright 2014, Raj Setaluri.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class HIERARCHICAL_ADAPTIVE_FINE_RANGE_RASTERIZER
//#####################################################################
#ifndef __HIERARCHICAL_ADAPTIVE_FINE_RANGE_RASTERIZER__
#define __HIERARCHICAL_ADAPTIVE_FINE_RANGE_RASTERIZER__

#include <SPGrid_Fluids/Grids/Rasterizers/HIERARCHICAL_RASTERIZER.h>
#include <SPGrid_Fluids/Grids/GRID_HIERARCHY.h>
#include <SPGrid_Fluids/Simulation/FLUIDS_SIMULATION_FLAGS.h>

using namespace SPGrid;

namespace PhysBAM{

template<class T_STRUCT,class T,int d>
class HIERARCHICAL_ADAPTIVE_FINE_RANGE_RASTERIZER:public HIERARCHICAL_RASTERIZER<T_STRUCT,T,d>
{
    typedef HIERARCHICAL_RASTERIZER<T_STRUCT,T,d> BASE;
    typedef VECTOR<T,d> TV;
    typedef VECTOR<int,d> T_INDEX;
    typedef RANGE<TV> T_SCALAR_RANGE;
    typedef PAIR<int,T_INDEX> T_CELL;

    const ARRAY<T_SCALAR_RANGE> ranges;

public:
    using BASE::hierarchy;

    HIERARCHICAL_ADAPTIVE_FINE_RANGE_RASTERIZER(GRID_HIERARCHY<T_STRUCT,T,d>& hierarchy_input,const ARRAY<T_SCALAR_RANGE>& ranges_input)
        :BASE(hierarchy_input),ranges(ranges_input)
    {}

    bool Consume(const T_CELL& cell)
    {
        const int level=cell.x;
        const T_INDEX& index=cell.y;

        PHYSBAM_ASSERT(hierarchy.Grid(level).Domain_Indices().Lazy_Inside(index));
        PHYSBAM_ASSERT(level>=1);

        T_SCALAR_RANGE scalar_range=hierarchy.Grid(level).Cell_Domain(index);
        const TV X=scalar_range.Center();

        if(level==1){hierarchy.Activate_Cell(level,index,SPGrid_Cell_Type_Interior);return false;} // if at bottom then activate

        for(int i=1;i<=ranges.m;++i) if(ranges(i).Lazy_Inside(X)) return true; // if inside a range then recurse
        
        hierarchy.Activate_Cell(level,index,SPGrid_Cell_Type_Interior); // if not at bottom, and not inside a range then activate
        return false;
    }
//#####################################################################
};
}
#endif
