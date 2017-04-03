//#####################################################################
// Copyright 2014, Raj Setaluri.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class HIERARCHICAL_ADAPTIVE_RANDOM_RASTERIZER
//#####################################################################
#ifndef __HIERARCHICAL_ADAPTIVE_RANDOM_RASTERIZER__
#define __HIERARCHICAL_ADAPTIVE_RANDOM_RASTERIZER__

#include <SPGrid_Fluids/Grids/Rasterizers/HIERARCHICAL_RASTERIZER.h>
#include <SPGrid_Fluids/Grids/GRID_HIERARCHY.h>
#include <SPGrid_Fluids/Simulation/FLUIDS_SIMULATION_FLAGS.h>
#include <PhysBAM_Tools/Random_Numbers/RANDOM_NUMBERS.h>

using namespace SPGrid;

namespace PhysBAM{

template<class T_STRUCT,class T,int d>
class HIERARCHICAL_ADAPTIVE_RANDOM_RASTERIZER:public HIERARCHICAL_RASTERIZER<T_STRUCT,T,d>
{
    typedef HIERARCHICAL_RASTERIZER<T_STRUCT,T,d> BASE;
    typedef VECTOR<T,d> TV;
    typedef VECTOR<int,d> T_INDEX;
    typedef RANGE<TV> T_SCALAR_RANGE;
    typedef PAIR<int,T_INDEX> T_CELL;

    RANDOM_NUMBERS<T> random;

public:
    using BASE::hierarchy;

    HIERARCHICAL_ADAPTIVE_RANDOM_RASTERIZER(GRID_HIERARCHY<T_STRUCT,T,d>& hierarchy_input)
        :BASE(hierarchy_input)
    {
        random.Set_Seed(1);
    }

    bool Consume(const T_CELL& cell)
    {
        const int level=cell.x;
        const T_INDEX& index=cell.y;

        PHYSBAM_ASSERT(hierarchy.Grid(level).Domain_Indices().Lazy_Inside(index));
        PHYSBAM_ASSERT(level>=1);

        if(level==1){hierarchy.Activate_Cell(level,index,SPGrid_Cell_Type_Interior);return false;}
        if(random.Get_Uniform_Integer(0,1)) return true; // randomly choose

        hierarchy.Activate_Cell(level,index,SPGrid_Cell_Type_Interior);
        return false;
    }
//#####################################################################
};
}
#endif
