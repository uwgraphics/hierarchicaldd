//#####################################################################
// Copyright 2012, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class HIERARCHICAL_ADAPTIVE_SPHERE_RASTERIZER
//#####################################################################
#ifndef __HIERARCHICAL_ADAPTIVE_SPHERE_RASTERIZER__
#define __HIERARCHICAL_ADAPTIVE_SPHERE_RASTERIZER__

#include <SPGrid_Fluids/Grids/Rasterizers/HIERARCHICAL_RASTERIZER.h>
#include <Common_Tools/Math_Tools/RANGE_ITERATOR.h>
#include <SPGrid_Fluids/Grids/GRID_HIERARCHY.h>
#include <SPGrid_Fluids/Simulation/FLUIDS_SIMULATION_FLAGS.h>

using namespace SPGrid;

namespace PhysBAM{

template<class T_STRUCT,class T,int d>
class HIERARCHICAL_ADAPTIVE_SPHERE_RASTERIZER:public HIERARCHICAL_RASTERIZER<T_STRUCT,T,d>
{
    typedef HIERARCHICAL_RASTERIZER<T_STRUCT,T,d> BASE;
    typedef VECTOR<T,d> TV;
    typedef VECTOR<int,d> T_INDEX;
    typedef RANGE<TV> T_SCALAR_RANGE;
    typedef PAIR<int,T_INDEX> T_CELL;

    const TV center;
    const T radius;
    const T radius_squared;
    const int shell_width;
    const TV prolate_factor;
    const RIGID_GEOMETRY<TV>* const sphere_rigid_geometry;
    bool initialize_inside;
    int max_level;
    const int grade;

public:
    using BASE::hierarchy;using BASE::Inside_Object;using BASE::Maximum_Distance_Squared;

    HIERARCHICAL_ADAPTIVE_SPHERE_RASTERIZER(GRID_HIERARCHY<T_STRUCT,T,d>& hierarchy_input,const TV& center_input,
        const T radius_input,const int shell_width_input,const TV& prolate_factor_input,
        const RIGID_GEOMETRY<TV>* const sphere_rigid_geometry_input=0,const bool initialize_inside_input=false,const int max_level_input=-1,const int grade_input=1)
        :BASE(hierarchy_input),center(center_input),radius(radius_input),radius_squared(radius*radius),
         shell_width(shell_width_input),prolate_factor(prolate_factor_input),sphere_rigid_geometry(sphere_rigid_geometry_input),initialize_inside(initialize_inside_input),
         max_level(max_level_input),grade(grade_input)
    {
        if(max_level==-1)
            max_level = hierarchy.Levels();
    }

    T Minimum_Distance_Squared(const T_SCALAR_RANGE& scalar_range)
    {
        return (scalar_range.Clamp(center)-center).Magnitude_Squared();
    }

    bool Consume(const T_CELL& cell)
    {
        const int level=cell.x;
        const T_INDEX& index=cell.y;

        PHYSBAM_ASSERT(hierarchy.Grid(level).Domain_Indices().Lazy_Inside(index));
        PHYSBAM_ASSERT(level>=1);

        T_SCALAR_RANGE scalar_range=hierarchy.Grid(level).Cell_Domain(index);

        T minimum_distance_squared=Minimum_Distance_Squared(scalar_range);

        // Simple base cases
        if(sphere_rigid_geometry){
            if(Inside_Object(scalar_range,*sphere_rigid_geometry)) return false;} // Block is inside the sphere, no need to recurse
        else{
            if(Maximum_Distance_Squared(scalar_range,center)<radius_squared && !initialize_inside) return false;} // Block is inside the sphere, no need to recurse                
        if(level==1) { hierarchy.Activate_Cell(level,index,SPGrid_Cell_Type_Interior); return false; }
        if(minimum_distance_squared<radius_squared) return true; // Block spans sphere edge, recurse
        if( grade>2 && (level%grade - 1) != 0) return true; // If not along graded levels then recurse

        T limit = radius + shell_width*scalar_range.Edge_Lengths()(1);
        T limit_squared = limit*limit;

        TV dX=scalar_range.Clamp(center)-center;
        for(int v=1;v<=d;v++) if(dX(v)>0) dX(v)/=prolate_factor(v);
        minimum_distance_squared=dX.Magnitude_Squared();

        if(minimum_distance_squared > limit_squared) 
        { 
            if(level<=max_level) hierarchy.Activate_Cell(level,index,SPGrid_Cell_Type_Interior);
            return false; 
        }
        return true;
    }
//#####################################################################
};
}
#endif
