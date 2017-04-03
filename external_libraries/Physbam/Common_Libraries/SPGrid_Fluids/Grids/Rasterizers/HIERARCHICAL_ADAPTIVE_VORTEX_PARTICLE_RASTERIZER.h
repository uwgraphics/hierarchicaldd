//#####################################################################
// Copyright 2014, Mridul Aanjaneya.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class HIERARCHICAL_ADAPTIVE_VORTEX_PARTICLE_RASTERIZER
//#####################################################################
#ifndef __HIERARCHICAL_ADAPTIVE_VORTEX_PARTICLE_RASTERIZER__
#define __HIERARCHICAL_ADAPTIVE_VORTEX_PARTICLE_RASTERIZER__

#include <PhysBAM_Geometry/Spatial_Acceleration/PARTICLE_HIERARCHY.h>
#include <SPGrid_Fluids/Grids/Rasterizers/HIERARCHICAL_RASTERIZER.h>
#include <SPGrid_Fluids/Grids/GRID_HIERARCHY.h>
#include <SPGrid_Fluids/Particles/VORTEX_PARTICLES.h>
using namespace SPGrid;

namespace PhysBAM{
template<class T_STRUCT,class T,int d>
class HIERARCHICAL_ADAPTIVE_VORTEX_PARTICLE_RASTERIZER:public HIERARCHICAL_RASTERIZER<T_STRUCT,T,d>
{
    typedef HIERARCHICAL_RASTERIZER<T_STRUCT,T,d> BASE;
    typedef VECTOR<T,d> TV;
    typedef VECTOR<int,d> T_INDEX;
    typedef RANGE<TV> T_SCALAR_RANGE;
    typedef PAIR<int,T_INDEX> T_CELL;

    const VORTEX_PARTICLES<TV>& vortex_particles;
    PARTICLE_HIERARCHY<TV,ARRAY_VIEW<TV> > particle_hierarchy;
    T particle_radius;
    const RIGID_GEOMETRY<TV>* const rigid_geometry;
    const int rigid_geometry_shell_width;
    const int particle_shell_width;

public:
    using BASE::hierarchy;using BASE::Inside_Object;

    HIERARCHICAL_ADAPTIVE_VORTEX_PARTICLE_RASTERIZER(GRID_HIERARCHY<T_STRUCT,T,d>& hierarchy_input,const VORTEX_PARTICLES<TV>& vortex_particles_input,const RIGID_GEOMETRY<TV>* const rigid_geometry_input=0,const int rigid_geometry_shell_width_input=4,const int particle_shell_width_input=5);

    bool Consume(const T_CELL& cell);
};
}
#endif
