//#####################################################################
// Copyright 2014, Mridul Aanjaneya, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class HIERARCHICAL_ADAPTIVE_FLASK_RASTERIZER
//#####################################################################
#ifndef __HIERARCHICAL_ADAPTIVE_FLASK_RASTERIZER__
#define __HIERARCHICAL_ADAPTIVE_FLASK_RASTERIZER__

#include <PhysBAM_Geometry/Basic_Geometry/SPHERE.h>
#include <PhysBAM_Geometry/Basic_Geometry/CYLINDER.h>
#include <PhysBAM_Geometry/Solids_Geometry/RIGID_GEOMETRY_COLLECTION.h>
#include <SPGrid_Fluids/Grids/Rasterizers/HIERARCHICAL_RASTERIZER.h>
#include <Common_Tools/Math_Tools/RANGE_ITERATOR.h>
#include <SPGrid_Fluids/Grids/GRID_HIERARCHY.h>
#include <SPGrid_Fluids/Simulation/FLUIDS_SIMULATION_FLAGS.h>
using namespace SPGrid;

namespace PhysBAM{
template<class T_STRUCT,class T,int d>
class HIERARCHICAL_ADAPTIVE_FLASK_RASTERIZER:public HIERARCHICAL_RASTERIZER<T_STRUCT,T,d>
{
    typedef HIERARCHICAL_RASTERIZER<T_STRUCT,T,d> BASE;
    typedef VECTOR<T,d> TV;
    typedef VECTOR<int,d> T_INDEX;
    typedef RANGE<TV> T_SCALAR_RANGE;
    typedef PAIR<int,T_INDEX> T_CELL;

    const int shell_width;
    const TV prolate_factor;
    const ARRAY<SPHERE<TV> >& spheres;
    const RIGID_GEOMETRY_COLLECTION<TV>& rigid_geometry_collection;
    int max_level;
    const int grade;

public:
    using BASE::hierarchy;using BASE::Intersect_Object;using BASE::Maximum_Distance_Squared;

    HIERARCHICAL_ADAPTIVE_FLASK_RASTERIZER(GRID_HIERARCHY<T_STRUCT,T,d>& hierarchy_input,
        const int shell_width_input,const TV& prolate_factor_input,const ARRAY<SPHERE<TV> >& spheres_input,const RIGID_GEOMETRY_COLLECTION<TV>& rigid_geometry_collection_input,
        const int max_level_input=-1,const int grade_input=1)
        :BASE(hierarchy_input),shell_width(shell_width_input),prolate_factor(prolate_factor_input),spheres(spheres_input),rigid_geometry_collection(rigid_geometry_collection_input),max_level(max_level_input),grade(grade_input)
    {
        if(max_level==-1) max_level=hierarchy.Levels();
    }

    T Distance_From_Surface(const T_SCALAR_RANGE& scalar_range,const int sphere_index)
    {
        const TV& center=spheres(sphere_index).center;
        return spheres(sphere_index).radius-(scalar_range.Clamp(center)-center).Magnitude();
    }

    bool Consume(const T_CELL& cell)
    {
        const int level=cell.x;
        const T_INDEX& index=cell.y;

        PHYSBAM_ASSERT(hierarchy.Grid(level).Domain_Indices().Lazy_Inside(index));
        PHYSBAM_ASSERT(level>=1);

        T_SCALAR_RANGE scalar_range=hierarchy.Grid(level).Cell_Domain(index);

        int sphere_index=0,cylinder_index=0;
        for(int i=1;i<=spheres.Size();++i)
            if(Intersect_Object(scalar_range,rigid_geometry_collection.Rigid_Geometry(i))){sphere_index=i;break;}
        for(int i=spheres.Size()+1;i<=rigid_geometry_collection.particles.array_collection->Size();++i)
            if(Intersect_Object(scalar_range,rigid_geometry_collection.Rigid_Geometry(i))){cylinder_index=i;break;}
        if(sphere_index==0 && cylinder_index==0) return false;

        if(level==1){hierarchy.Activate_Cell(level,index,SPGrid_Cell_Type_Interior);return false;}
        if(cylinder_index) return true;
        if( grade>2 && (level%grade - 1) != 0) return true; // If not along graded levels then recurse

        T limit=shell_width*scalar_range.Edge_Lengths()(1);
        T distance_from_surface=Distance_From_Surface(scalar_range,sphere_index);

        if(distance_from_surface>limit)
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
