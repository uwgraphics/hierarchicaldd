//#####################################################################
// Copyright 2014, Mridul Aanjaneya.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class HIERARCHICAL_ADAPTIVE_VORTEX_PARTICLE_RASTERIZER
//#####################################################################
#include <SPGrid_Fluids/Grids/Rasterizers/HIERARCHICAL_ADAPTIVE_VORTEX_PARTICLE_RASTERIZER.h>
#include <SPGrid_Fluids/Simulation/FLUIDS_SIMULATION_FLAGS.h>
#include <SPGrid_Fluids/Simulation/FLUIDS_SIMULATION_DATA.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T_STRUCT,class T,int d> HIERARCHICAL_ADAPTIVE_VORTEX_PARTICLE_RASTERIZER<T_STRUCT,T,d>::
HIERARCHICAL_ADAPTIVE_VORTEX_PARTICLE_RASTERIZER(GRID_HIERARCHY<T_STRUCT,T,d>& hierarchy_input,const VORTEX_PARTICLES<VECTOR<T,d> >& vortex_particles_input,const RIGID_GEOMETRY<TV>* const rigid_geometry_input,const int rigid_geometry_shell_width_input,const int particle_shell_width_input)
    :BASE(hierarchy_input),vortex_particles(vortex_particles_input),particle_hierarchy(vortex_particles.X,true,1),rigid_geometry(rigid_geometry_input),rigid_geometry_shell_width(rigid_geometry_shell_width_input),particle_shell_width(particle_shell_width_input)
{particle_radius=particle_shell_width*hierarchy.Grid(1).dX(1);}
//#####################################################################
// 
//#####################################################################
template<class T_STRUCT,class T,int d> bool HIERARCHICAL_ADAPTIVE_VORTEX_PARTICLE_RASTERIZER<T_STRUCT,T,d>::
Consume(const PAIR<int,VECTOR<int,d> >& cell)
{
    const int level=cell.x;const T_INDEX& index=cell.y;

    PHYSBAM_ASSERT(hierarchy.Grid(level).Domain_Indices().Lazy_Inside(index));
    PHYSBAM_ASSERT(level>=1);

    T_SCALAR_RANGE scalar_range=hierarchy.Grid(level).Cell_Domain(index);
    T particle_limit_squared=sqr(particle_radius+particle_shell_width*scalar_range.Edge_Lengths()(1));
    T rigid_geometry_limit_squared=sqr(rigid_geometry_shell_width*scalar_range.Edge_Lengths()(1));
    const TV X=scalar_range.Center();

    if(rigid_geometry && Inside_Object(scalar_range,*rigid_geometry)) return false; // Block is inside the sphere, no need to recurse
    if(level==1){hierarchy.Activate_Cell(level,index,SPGrid_Cell_Type_Interior);return false;}

    T distance_to_rigid_geometry_squared=rigid_geometry?sqr(rigid_geometry->Implicit_Geometry_Extended_Value(X)):std::numeric_limits<T>::max();
    T particle_scale=(T)0.;
    for(int level=1;level<=hierarchy.Levels();++level) particle_scale+=particle_shell_width*hierarchy.Grid(level).dX(1);

    ARRAY<int> intersection_list;
    particle_hierarchy.Intersection_List(scalar_range,intersection_list,particle_scale);
    T minimum_distance_squared=std::numeric_limits<T>::max();
    for(int i=1;i<=intersection_list.m;i++){int p=intersection_list(i);
        const TV& Xp=vortex_particles.X(p);const TV X_closest=scalar_range.Clamp(Xp);
        T distance_squared=(Xp-X_closest).Magnitude_Squared();
        if(distance_squared<sqr(particle_radius)) return true;
        minimum_distance_squared=min(minimum_distance_squared,distance_squared);}

    if(minimum_distance_squared>particle_limit_squared && distance_to_rigid_geometry_squared>rigid_geometry_limit_squared){
        hierarchy.Activate_Cell(level,index,SPGrid_Cell_Type_Interior);
        return false;}

    return true;
}
//#####################################################################
template class HIERARCHICAL_ADAPTIVE_VORTEX_PARTICLE_RASTERIZER<FLUIDS_SIMULATION_DATA<float>,float,2>;
template class HIERARCHICAL_ADAPTIVE_VORTEX_PARTICLE_RASTERIZER<FLUIDS_SIMULATION_DATA<float>,float,3>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class HIERARCHICAL_ADAPTIVE_VORTEX_PARTICLE_RASTERIZER<FLUIDS_SIMULATION_DATA<double>,double,2>;
template class HIERARCHICAL_ADAPTIVE_VORTEX_PARTICLE_RASTERIZER<FLUIDS_SIMULATION_DATA<double>,double,3>;
#endif
