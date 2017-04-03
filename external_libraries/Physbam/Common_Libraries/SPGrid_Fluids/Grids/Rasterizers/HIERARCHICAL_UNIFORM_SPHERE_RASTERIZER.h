//#####################################################################
// Copyright 2012, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class HIERARCHICAL_UNIFORM_SPHERE_RASTERIZER
//#####################################################################
#ifndef __HIERARCHICAL_UNIFORM_SPHERE_RASTERIZER__
#define __HIERARCHICAL_UNIFORM_SPHERE_RASTERIZER__

#include <SPGrid_Fluids/Grids/Rasterizers/HIERARCHICAL_RASTERIZER.h>
#include <Common_Tools/Math_Tools/RANGE_ITERATOR.h>
#include <SPGrid_Fluids/Grids/GRID_HIERARCHY.h>
#include <SPGrid_Fluids/Simulation/FLUIDS_SIMULATION_FLAGS.h>

using namespace SPGrid;

namespace PhysBAM{

template<class T_STRUCT,class T,int d>
class HIERARCHICAL_UNIFORM_SPHERE_RASTERIZER:public HIERARCHICAL_RASTERIZER<T_STRUCT,T,d>
{
    typedef HIERARCHICAL_RASTERIZER<T_STRUCT,T,d> BASE;
    typedef VECTOR<T,d> TV;
    typedef VECTOR<int,d> T_INDEX;
    typedef RANGE<TV> T_SCALAR_RANGE;
    typedef PAIR<int,T_INDEX> T_CELL;

    const TV center;
    const T radius;
    const T radius_squared;
    const RIGID_GEOMETRY<TV>* const sphere_rigid_geometry;

public:
    using BASE::hierarchy;using BASE::Inside_Object;using BASE::Maximum_Distance_Squared;

    HIERARCHICAL_UNIFORM_SPHERE_RASTERIZER(GRID_HIERARCHY<T_STRUCT,T,d>& hierarchy_input,const TV& center_input,const T radius_input,const RIGID_GEOMETRY<TV>* const sphere_rigid_geometry_input=0)
        :BASE(hierarchy_input),center(center_input),radius(radius_input),radius_squared(radius*radius),sphere_rigid_geometry(sphere_rigid_geometry_input)
    {}

    bool Consume(const T_CELL& cell)
    {
        const int level=cell.x;
        const T_INDEX& index=cell.y;
        
        T_SCALAR_RANGE scalar_range=hierarchy.Grid(level).Cell_Domain(index);       

        const int cell_size=1<<(level-1);
        RANGE<T_INDEX> fine_indices((index-1)*cell_size+1,index*cell_size);
        PHYSBAM_ASSERT(hierarchy.Grid(1).Domain_Indices().Lazy_Inside(fine_indices.min_corner));
        PHYSBAM_ASSERT(hierarchy.Grid(1).Domain_Indices().Lazy_Inside(fine_indices.max_corner));

        if(sphere_rigid_geometry){
            if(level==1 && !(Inside_Object(scalar_range,*sphere_rigid_geometry)))
                hierarchy.Activate_Cell(level,index,SPGrid_Cell_Type_Interior);}
        else{
            if(level==1 && !(Maximum_Distance_Squared(scalar_range,center)<radius_squared)){
                hierarchy.Activate_Cell(level,index,SPGrid_Cell_Type_Interior);}}
        return (level>1);
    }
//#####################################################################
};
}
#endif
