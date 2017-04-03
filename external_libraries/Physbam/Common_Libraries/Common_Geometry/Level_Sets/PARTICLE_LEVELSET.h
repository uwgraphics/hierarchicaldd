//#####################################################################
// Copyright 2016, Mridul Aanjaneya.
// This file is part of PhysBAM whose distribution is governed by the license 
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class PARTICLE_LEVELSET
//#####################################################################
#ifndef __PARTICLE_LEVELSET__
#define __PARTICLE_LEVELSET__

#include <Common_Geometry/Particles/PARTICLE_LEVELSET_PARTICLES.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/FAST_LEVELSET.h>
#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/FACE_ARRAYS.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Grids_Uniform_Boundaries/BOUNDARY_POLICY_UNIFORM.h>
#include <PhysBAM_Tools/Grids_Uniform_Interpolation/FACE_LOOKUP_UNIFORM.h>
#include <PhysBAM_Tools/Grids_Uniform_Interpolation/LINEAR_INTERPOLATION_UNIFORM.h>
#include <PhysBAM_Tools/Random_Numbers/RANDOM_NUMBERS.h>

namespace PhysBAM{
template<class TV>
class PARTICLE_LEVELSET
{
    typedef typename TV::SCALAR T;typedef VECTOR<int,TV::dimension> T_INDEX;
    typedef GRID<TV> T_GRID;typedef typename T_GRID::CELL_ITERATOR CELL_ITERATOR;
    typedef ARRAY<T,T_INDEX> T_SCALAR_FIELD;typedef ARRAY<T,FACE_INDEX<TV::dimension> > T_FACE_SCALAR_FIELD;
    typedef ARRAY<PARTICLE_LEVELSET_PARTICLES<TV>*,T_INDEX> T_ARRAYS_PARTICLE_LEVELSET_PARTICLES;
    typedef typename BOUNDARY_POLICY<T_GRID>::BOUNDARY_SCALAR T_BOUNDARY_SCALAR;
    typedef FACE_LOOKUP_UNIFORM<T_GRID> T_FACE_LOOKUP;
    enum {d=TV::dimension};
  public:
    int number_of_particles_per_cell,number_of_ghost_cells;
    T half_band_width;
    T minimum_particle_radius,maximum_particle_radius;
    T outside_particle_distance_multiplier; // how far outside a particle needs to be before it counts as outside
    int maximum_iterations_for_attraction;
    RANDOM_NUMBERS<T> random;

    FAST_LEVELSET<T_GRID> levelset;
    T_FACE_SCALAR_FIELD& face_velocities;
    T_BOUNDARY_SCALAR& boundary;
    T_ARRAYS_PARTICLE_LEVELSET_PARTICLES positive_particles,negative_particles;
    LINEAR_INTERPOLATION_UNIFORM<T_GRID,TV> interpolation;

    PARTICLE_LEVELSET(GRID<TV>& grid,T_SCALAR_FIELD& phi,T_FACE_SCALAR_FIELD& face_velocities_input,T_BOUNDARY_SCALAR& boundary_input,const int number_of_ghost_cells_input=3)
        :number_of_ghost_cells(number_of_ghost_cells_input),levelset(grid,phi,number_of_ghost_cells),face_velocities(face_velocities_input),boundary(boundary_input)
    {
        random.Set_Seed(0);
        Set_Outside_Particle_Distance();
	Set_Maximum_Iterations_For_Attraction();
        Set_Number_Of_Particles_Per_Cell(16);
    }

    virtual ~PARTICLE_LEVELSET();

    void Set_Number_Of_Particles_Per_Cell(const int number)
    {number_of_particles_per_cell=number;}

    void Set_Minimum_Particle_Radius(const T radius)
    {minimum_particle_radius=radius;}

    void Set_Maximum_Particle_Radius(const T radius)
    {maximum_particle_radius=radius;}

    void Set_Band_Width(const int number_of_cells=6)
    {half_band_width=number_of_cells*levelset.grid.dX.Max()*(T).5;}

    void Set_Outside_Particle_Distance(const T fraction_of_the_radius=1)
    {outside_particle_distance_multiplier=fraction_of_the_radius;}

    void Set_Maximum_Iterations_For_Attraction(const int maximum_iterations=15)
    {maximum_iterations_for_attraction=maximum_iterations;}

//#####################################################################
    void Initialize_Particle_Levelset_Grid_Values();
    void Add_Particles_To_Both_Sides_Of_Interface_Band();
    void Attract_Particles_To_Interface_And_Adjust_Radii();
    void Euler_Step_Particles(const T dt,const T time);
    void Delete_Particles_Outside_Grid();
    void Modify_Levelset_Using_Escaped_Particles();
    void Adjust_Particle_Radii();
    void Reseed_Particles();
  protected:
    void Update_Cell_Containers(T_ARRAYS_PARTICLE_LEVELSET_PARTICLES& particles);
    void Add_Particles_To_Both_Sides_Of_Interface_Band(T_ARRAYS_PARTICLE_LEVELSET_PARTICLES& particles);
    void Attract_Particles_To_Interface_And_Adjust_Radii(T_ARRAYS_PARTICLE_LEVELSET_PARTICLES& particles,const int sign);
    void Attract_Individual_Particle_To_Interface_And_Adjust_Radius(PARTICLE_LEVELSET_PARTICLES<TV>& particles,const T& phi_min,const T& phi_max,const int k);
    void Euler_Step_Particles(T_ARRAYS_PARTICLE_LEVELSET_PARTICLES& particles,T_FACE_SCALAR_FIELD& face_velocities_ghost,const T dt);
    void Delete_Particles_Outside_Grid(T_ARRAYS_PARTICLE_LEVELSET_PARTICLES& particles);
    void Modify_Levelset_Using_Escaped_Particles(T_SCALAR_FIELD& phi,T_ARRAYS_PARTICLE_LEVELSET_PARTICLES& particles,const int sign);
    void Adjust_Particle_Radii(T_ARRAYS_PARTICLE_LEVELSET_PARTICLES& particles,const int sign);
//#####################################################################
};
}
#endif
