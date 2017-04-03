//#####################################################################
// Copyright 2016, Mridul Aanjaneya.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class HIERARCHICAL_UNIFORM_INVERSE_GEOMETRY_RASTERIZER
//#####################################################################
#ifndef __HIERARCHICAL_UNIFORM_INVERSE_GEOMETRY_RASTERIZER__
#define __HIERARCHICAL_UNIFORM_INVERSE_GEOMETRY_RASTERIZER__

#include <SPGrid_Fluids/Grids/Rasterizers/HIERARCHICAL_RASTERIZER.h>
#include <Common_Tools/Math_Tools/RANGE_ITERATOR.h>
#include <SPGrid_Fluids/Grids/GRID_HIERARCHY.h>
#include <SPGrid_Fluids/Simulation/FLUIDS_SIMULATION_FLAGS.h>

using namespace SPGrid;

namespace PhysBAM{
template<class T_STRUCT,class T,int d>
class HIERARCHICAL_UNIFORM_INVERSE_GEOMETRY_RASTERIZER:public HIERARCHICAL_RASTERIZER<T_STRUCT,T,d>
{
    typedef HIERARCHICAL_RASTERIZER<T_STRUCT,T,d> BASE;
    typedef VECTOR<T,d> TV;typedef VECTOR<int,d> T_INDEX;
    typedef RANGE<TV> T_SCALAR_RANGE;typedef PAIR<int,T_INDEX> T_CELL;

    const RIGID_GEOMETRY<TV>* const rigid_geometry;

public:
    using BASE::hierarchy;

    HIERARCHICAL_UNIFORM_INVERSE_GEOMETRY_RASTERIZER(GRID_HIERARCHY<T_STRUCT,T,d>& hierarchy_input,const RIGID_GEOMETRY<TV>* const rigid_geometry_input)
        :BASE(hierarchy_input),rigid_geometry(rigid_geometry_input)
    {}

    bool Consume(const T_CELL& cell)
    {
        const int level=cell.x;
        const T_INDEX& index=cell.y;
        
        TV center=hierarchy.Grid(level).Center(index);

        const int cell_size=1<<(level-1);
        RANGE<T_INDEX> fine_indices((index-1)*cell_size+1,index*cell_size);
        PHYSBAM_ASSERT(hierarchy.Grid(1).Domain_Indices().Lazy_Inside(fine_indices.min_corner));
        PHYSBAM_ASSERT(hierarchy.Grid(1).Domain_Indices().Lazy_Inside(fine_indices.max_corner));

        if(level==1 && rigid_geometry->Implicit_Geometry_Extended_Value(center)>(T)0.)
            hierarchy.Activate_Cell(level,index,SPGrid_Cell_Type_Interior);
        return (level>1);
    }
//#####################################################################
};
}
#endif
