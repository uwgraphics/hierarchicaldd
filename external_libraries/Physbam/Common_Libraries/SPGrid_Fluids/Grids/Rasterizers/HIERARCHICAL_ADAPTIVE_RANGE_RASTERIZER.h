//#####################################################################
// Copyright 2012, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class HIERARCHICAL_ADAPTIVE_RANGE_RASTERIZER
//#####################################################################
#ifndef __HIERARCHICAL_ADAPTIVE_RANGE_RASTERIZER__
#define __HIERARCHICAL_ADAPTIVE_RANGE_RASTERIZER__

#include <SPGrid_Fluids/Grids/Rasterizers/HIERARCHICAL_RASTERIZER.h>
#include <SPGrid_Fluids/Grids/GRID_HIERARCHY.h>
#include <SPGrid_Fluids/Simulation/FLUIDS_SIMULATION_FLAGS.h>

using namespace SPGrid;

namespace PhysBAM{

template<class T_STRUCT,class T,int d>
class HIERARCHICAL_ADAPTIVE_RANGE_RASTERIZER:public HIERARCHICAL_RASTERIZER<T_STRUCT,T,d>
{
    typedef HIERARCHICAL_RASTERIZER<T_STRUCT,T,d> BASE;
    typedef VECTOR<T,d> TV;
    typedef VECTOR<int,d> T_INDEX;
    typedef RANGE<TV> T_SCALAR_RANGE;
    typedef PAIR<int,T_INDEX> T_CELL;

    const ARRAY<T_SCALAR_RANGE> ranges;
    const TV center;
    const T radius;
    const T radius_squared;
    const RIGID_GEOMETRY<TV>* const sphere_rigid_geometry;

public:
    using BASE::hierarchy;using BASE::Inside_Object;using BASE::Maximum_Distance_Squared;

    HIERARCHICAL_ADAPTIVE_RANGE_RASTERIZER(GRID_HIERARCHY<T_STRUCT,T,d>& hierarchy_input,const ARRAY<T_SCALAR_RANGE>& ranges_input,const TV& center_input,const T radius_input,const RIGID_GEOMETRY<TV>* const sphere_rigid_geometry_input=0)
        :BASE(hierarchy_input),ranges(ranges_input),center(center_input),radius(radius_input),radius_squared(radius*radius),sphere_rigid_geometry(sphere_rigid_geometry_input)
    {}

    T Minimum_Distance(const T_SCALAR_RANGE& scalar_range,const ARRAY<T_SCALAR_RANGE>& ranges)
    {
        const TV scalar_range_center=scalar_range.Center();
        T min_distance=std::numeric_limits<T>::max();
        for(int i=1;i<=ranges.m;++i){
            T current_distance=(ranges(i).Clamp(scalar_range_center)-scalar_range.Clamp(ranges(i).Center())).Max_Abs();
            min_distance=min(min_distance,current_distance);}
        return min_distance;
    }

    bool Consume(const T_CELL& cell)
    {
        const int level=cell.x;
        const T_INDEX& index=cell.y;

        PHYSBAM_ASSERT(hierarchy.Grid(level).Domain_Indices().Lazy_Inside(index));
        PHYSBAM_ASSERT(level>=1);

        T_SCALAR_RANGE scalar_range=hierarchy.Grid(level).Cell_Domain(index);
        const TV X=scalar_range.Center();

        if(level==1){hierarchy.Activate_Cell(level,index,SPGrid_Cell_Type_Interior);return false;}

        for(int i=1;i<=ranges.m;++i) if(ranges(i).Lazy_Inside(X)) return true;
        int shell_width=8;
        const T minimum_distance=Minimum_Distance(scalar_range,ranges);
        T limit=shell_width*scalar_range.Edge_Lengths()(1);
        if(minimum_distance>limit){ 
            hierarchy.Activate_Cell(level,index,SPGrid_Cell_Type_Interior);
            return false;}

        return true;
    }
//#####################################################################
};
}
#endif
