//#####################################################################
// Copyright 2012, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class HIERARCHICAL_RASTERIZER
//#####################################################################
#ifndef __HIERARCHICAL_RASTERIZER__
#define __HIERARCHICAL_RASTERIZER__

#include <SPGrid_Fluids/Grids/GRID_HIERARCHY.h>
#include <Common_Tools/Math_Tools/RANGE_ITERATOR.h>
#include <PhysBAM_Geometry/Solids_Geometry/RIGID_GEOMETRY.h>

using namespace SPGrid;

namespace PhysBAM{

template<class T_STRUCT,class T,int d>
class HIERARCHICAL_RASTERIZER
{

    typedef VECTOR<T,d> TV;
    typedef VECTOR<int,d> T_INDEX;
    typedef RANGE<TV> T_SCALAR_RANGE;
    typedef PAIR<int,T_INDEX> T_CELL;

public:
    GRID_HIERARCHY<T_STRUCT,T,d>& hierarchy;

    HIERARCHICAL_RASTERIZER(GRID_HIERARCHY<T_STRUCT,T,d>& hierarchy_input)
        :hierarchy(hierarchy_input)
    {}
    
    bool Inside_Object(const T_SCALAR_RANGE& scalar_range,const RIGID_GEOMETRY<TV>& rigid_geometry)
    {
        for(RANGE_ITERATOR<d> iterator(RANGE<T_INDEX>::Unit_Box());iterator.Valid();iterator.Next()){
            TV corner;
            for(int v=1;v<=d;v++) corner(v)=iterator.Index()(v)?scalar_range.max_corner(v):scalar_range.min_corner(v);
            if(!(rigid_geometry.Implicit_Geometry_Lazy_Inside(corner))) return false;}
        return true;
    }

    bool Intersect_Object(const T_SCALAR_RANGE& scalar_range,const RIGID_GEOMETRY<TV>& rigid_geometry)
    {
        for(RANGE_ITERATOR<d> iterator(RANGE<T_INDEX>::Unit_Box());iterator.Valid();iterator.Next()){
            TV corner;
            for(int v=1;v<=d;v++) corner(v)=iterator.Index()(v)?scalar_range.max_corner(v):scalar_range.min_corner(v);
            if(rigid_geometry.Implicit_Geometry_Lazy_Inside(corner)) return true;}
        return false;
    }

    T Maximum_Distance_Squared(const T_SCALAR_RANGE& scalar_range,const TV& center)
    {
        T distance_squared=0.f;
        for(RANGE_ITERATOR<d> iterator(RANGE<T_INDEX>::Unit_Box());iterator.Valid();iterator.Next()){
            TV corner;
            for(int v=1;v<=d;v++) corner(v)=iterator.Index()(v)?scalar_range.max_corner(v):scalar_range.min_corner(v);
            distance_squared=std::max(distance_squared,(corner-center).Magnitude_Squared());}
        return distance_squared;
    }

    virtual bool Consume(const T_CELL& cell)=0;
//#####################################################################
};
}
#endif
