//#####################################################################
// Copyright 2016, Mridul Aanjaneya.
// This file is part of PhysBAM whose distribution is governed by the license 
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Common_Geometry/Level_Sets/PARTICLE_LEVELSET.h>
#include <Common_Tools/Math_Tools/RANGE_ITERATOR.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Math_Tools/exchange.h>
using namespace PhysBAM;
//#####################################################################
// Destructor
//#####################################################################
template<class TV> PARTICLE_LEVELSET<TV>::
~PARTICLE_LEVELSET()
{
    positive_particles.Delete_Pointers_And_Clean_Memory();
    negative_particles.Delete_Pointers_And_Clean_Memory();
}
//#####################################################################
// Initialize_Particle_Levelset_Grid_Values
//#####################################################################
template<class TV> void PARTICLE_LEVELSET<TV>::
Initialize_Particle_Levelset_Grid_Values()
{
    positive_particles.Resize(levelset.grid.Cell_Indices());
    negative_particles.Resize(levelset.grid.Cell_Indices());
    Set_Minimum_Particle_Radius((T).1*levelset.grid.dX.Max());
    Set_Maximum_Particle_Radius((T).5*levelset.grid.dX.Max());
    Set_Band_Width();
}
//#####################################################################
// Add_Particles_To_Both_Sides_Of_Interface_Band
//#####################################################################
template<class TV> void PARTICLE_LEVELSET<TV>::
Add_Particles_To_Both_Sides_Of_Interface_Band()
{
    Add_Particles_To_Both_Sides_Of_Interface_Band(positive_particles);
    Add_Particles_To_Both_Sides_Of_Interface_Band(negative_particles);
}
//#####################################################################
// Add_Particles_To_Both_Sides_Of_Interface_Band
//#####################################################################
template<class TV> void PARTICLE_LEVELSET<TV>::
Add_Particles_To_Both_Sides_Of_Interface_Band(T_ARRAYS_PARTICLE_LEVELSET_PARTICLES& particles)
{
    for(CELL_ITERATOR iterator(levelset.grid);iterator.Valid();iterator.Next()){
        const T_INDEX& index=iterator.Cell_Index();ARRAY<TV> node_locations(T_GRID::number_of_nodes_per_cell);bool inside=false;
        levelset.grid.Node_Locations_In_Cell_From_Minimum_Corner_Node(index,node_locations);
        for(int i=1;i<=node_locations.Size();++i) if(fabs(levelset.Extended_Phi(node_locations(i)))<=half_band_width){inside=true;break;}

        if(inside){TV min_corner=levelset.grid.Node(index),max_corner=min_corner+levelset.grid.dX;
            if(!particles(index)) particles(index)=new PARTICLE_LEVELSET_PARTICLES<TV>();
            for(int k=1;k<=number_of_particles_per_cell;++k){TV location;
                for(int v=1;v<=d;++v) location(v)=random.Get_Uniform_Number(min_corner(v),max_corner(v));
                particles(index)->array_collection->Add_Elements(1);
                particles(index)->X(particles(index)->X.Size())=location;
                particles(index)->radius(particles(index)->X.Size())=maximum_particle_radius;}}}
}
//#####################################################################
// Attract_Particles_To_Interface_And_Adjust_Radii
//#####################################################################
template<class TV> void PARTICLE_LEVELSET<TV>::
Attract_Particles_To_Interface_And_Adjust_Radii()
{
    levelset.Compute_Normals();
    Attract_Particles_To_Interface_And_Adjust_Radii(negative_particles,-1);
    Attract_Particles_To_Interface_And_Adjust_Radii(positive_particles,1);
    Update_Cell_Containers(negative_particles);
    Update_Cell_Containers(positive_particles);
}
//#####################################################################
// Attract_Particles_To_Interface_And_Adjust_Radii
//#####################################################################
template<class TV> void PARTICLE_LEVELSET<TV>::
Attract_Particles_To_Interface_And_Adjust_Radii(T_ARRAYS_PARTICLE_LEVELSET_PARTICLES& particles,const int sign)
{
    T phi_min=sign*minimum_particle_radius,phi_max=sign*half_band_width;
    if(phi_min>phi_max) exchange(phi_min,phi_max);
    for(CELL_ITERATOR iterator(levelset.grid);iterator.Valid();iterator.Next()){const T_INDEX& index=iterator.Cell_Index();
        if(particles(index)) for(int k=particles(index)->array_collection->Size();k>=1;--k)
            Attract_Individual_Particle_To_Interface_And_Adjust_Radius(*(particles(index)),phi_min,phi_max,k);}
}
//#####################################################################
// Attract_Individual_Particle_To_Interface_And_Adjust_Radius
//#####################################################################
template<class TV> void PARTICLE_LEVELSET<TV>::
Attract_Individual_Particle_To_Interface_And_Adjust_Radius(PARTICLE_LEVELSET_PARTICLES<TV>& particles,const T& phi_min,const T& phi_max,const int k)
{
    T phi=levelset.Extended_Phi(particles.X(k));
    bool inside=(phi>=phi_min && phi<=phi_max);
    if(!inside){T phi_goal=random.Get_Uniform_Number(phi_min,phi_max);TV X_new;int iteration=0;
        while(!inside && iteration<=maximum_iterations_for_attraction){iteration++;
            TV normal=levelset.Extended_Normal(particles.X(k)).Normalized();
            T distance=phi_goal-phi;T dt=(T)1.;bool inside_domain=false;
            while(!inside_domain){X_new=particles.X(k)+normal*dt*distance;
                if(levelset.grid.domain.Lazy_Outside(X_new)){inside_domain=false;dt*=(T).5;}
                else inside_domain=true;}
            phi=levelset.Extended_Phi(X_new);
            inside=(phi>=phi_min && phi<=phi_max);
            if(!inside){dt*=(T).5;
                X_new=particles.X(k)+normal*dt*distance;
                phi=levelset.Extended_Phi(X_new);
                inside=(phi>=phi_min && phi<=phi_max);}
            particles.X(k)=X_new;}}

    if(!inside) particles.array_collection->Delete_Element(k);
    else{   // adjust radii for particles that are inside
        T new_radius=fabs(phi);
        if(new_radius<minimum_particle_radius) particles.radius(k)=minimum_particle_radius;
        else if(new_radius>maximum_particle_radius) particles.radius(k)=maximum_particle_radius;
        else particles.radius(k)=new_radius;}
}
//#####################################################################
// Update_Cell_Containers
//#####################################################################
template<class TV> void PARTICLE_LEVELSET<TV>::
Update_Cell_Containers(T_ARRAYS_PARTICLE_LEVELSET_PARTICLES& particles)
{
    for(CELL_ITERATOR iterator(levelset.grid);iterator.Valid();iterator.Next()){const T_INDEX& index=iterator.Cell_Index();
        if(particles(index) && particles(index)->array_collection->Size()==0){delete particles(index);particles(index)=0;}
        else if(particles(index)) for(int i=particles(index)->array_collection->Size();i>=1;--i){const TV& location=particles(index)->X(i);
            T_INDEX new_cell=levelset.grid.Clamp_To_Cell(location);
            if(new_cell!=index){
                if(!particles(new_cell)) particles(new_cell)=new PARTICLE_LEVELSET_PARTICLES<TV>();
                particles(new_cell)->array_collection->Add_Elements(1);
                const int size=particles(new_cell)->array_collection->Size();
                particles(new_cell)->X(size)=location;
                particles(new_cell)->radius(size)=particles(index)->radius(i);
                if(particles(index)->array_collection->Size()>1) particles(index)->array_collection->Delete_Element(i);
                else{delete particles(index);particles(index)=0;break;}}}}
}
//#####################################################################
// Euler_Step_Particles
//#####################################################################
template<class TV> void PARTICLE_LEVELSET<TV>::
Euler_Step_Particles(const T dt,const T time)
{
    T_FACE_SCALAR_FIELD face_velocities_ghost(levelset.grid,number_of_ghost_cells);
    boundary.Fill_Ghost_Cells_Face(levelset.grid,face_velocities,face_velocities_ghost,time+dt,number_of_ghost_cells);
    Euler_Step_Particles(negative_particles,face_velocities_ghost,dt);
    Euler_Step_Particles(positive_particles,face_velocities_ghost,dt);
}
//#####################################################################
// Euler_Step_Particles
//#####################################################################
template<class TV> void PARTICLE_LEVELSET<TV>::
Euler_Step_Particles(T_ARRAYS_PARTICLE_LEVELSET_PARTICLES& particles,T_FACE_SCALAR_FIELD& face_velocities,const T dt)
{
    for(CELL_ITERATOR iterator(levelset.grid);iterator.Valid();iterator.Next()){const T_INDEX& index=iterator.Cell_Index();
        if(particles(index)) for(int k=1;k<=particles(index)->array_collection->Size();++k){
            TV fluid_velocity=interpolation.Clamped_To_Array_Face(levelset.grid,T_FACE_LOOKUP(face_velocities),particles(index)->X(k));
            particles(index)->X(k)+=dt*fluid_velocity;}}
}
//#####################################################################
// Delete_Particles_Outside_Grid
//#####################################################################
template<class TV> void PARTICLE_LEVELSET<TV>::
Delete_Particles_Outside_Grid()
{
    Delete_Particles_Outside_Grid(negative_particles);
    Delete_Particles_Outside_Grid(positive_particles);
    Update_Cell_Containers(negative_particles);
    Update_Cell_Containers(positive_particles);
}
//#####################################################################
// Delete_Particles_Outside_Grid
//#####################################################################
template<class TV> void PARTICLE_LEVELSET<TV>::
Delete_Particles_Outside_Grid(T_ARRAYS_PARTICLE_LEVELSET_PARTICLES& particles)
{
    for(CELL_ITERATOR iterator(levelset.grid);iterator.Valid();iterator.Next()){const T_INDEX& index=iterator.Cell_Index();
        if(particles(index)){for(int k=particles(index)->array_collection->Size();k>=1;--k){
            T_INDEX new_cell=levelset.grid.Clamp_To_Cell(particles(index)->X(k));
            if(!levelset.grid.Inside_Domain(new_cell)) particles(index)->array_collection->Delete_Element(k);}
            if(particles(index)->array_collection->Size()==0){delete particles(index);particles(index)=0;}}}
}
//#####################################################################
// Modify_Levelset_Using_Escaped_Particles
//#####################################################################
template<class TV> void PARTICLE_LEVELSET<TV>::
Modify_Levelset_Using_Escaped_Particles()
{
    T_SCALAR_FIELD phi_minus(levelset.phi),phi_plus(levelset.phi);
    Modify_Levelset_Using_Escaped_Particles(phi_minus,negative_particles,-1);
    Modify_Levelset_Using_Escaped_Particles(phi_plus,positive_particles,1);

    for(CELL_ITERATOR iterator(levelset.grid);iterator.Valid();iterator.Next()){const T_INDEX& index=iterator.Cell_Index();
        levelset.phi(index)=minmag(phi_minus(index),phi_plus(index));}
}
//#####################################################################
// Modify_Levelset_Using_Escaped_Particles
//#####################################################################
template<class TV> void PARTICLE_LEVELSET<TV>::
Modify_Levelset_Using_Escaped_Particles(T_SCALAR_FIELD& phi,T_ARRAYS_PARTICLE_LEVELSET_PARTICLES& particles,const int sign)
{
    for(CELL_ITERATOR iterator(levelset.grid);iterator.Valid();iterator.Next()){const T_INDEX& index=iterator.Cell_Index();
        if(particles(index)) for(int k=particles(index)->array_collection->Size();k>=1;--k){T phi_particle=levelset.Extended_Phi(particles(index)->X(k));
            if(-sign*phi_particle>outside_particle_distance_multiplier*particles(index)->radius(k)){T_INDEX base_index=index;
                for(int axis=1;axis<=d;++axis) if(particles(index)->X(k)(axis)<levelset.grid.Center(index)(axis)){
                    int factor=(base_index(axis)==1-number_of_ghost_cells)?0:1;base_index(axis)-=factor;}
                for(RANGE_ITERATOR<d> local_iterator(RANGE<T_INDEX>::Unit_Box());local_iterator.Valid();local_iterator.Next()){
                    T_INDEX current_index=base_index+local_iterator.Index();
                    T value=-sign*((levelset.grid.Center(current_index)-particles(index)->X(k)).Magnitude()-particles(index)->radius(k));
                    if(sign==-1) phi(current_index)=min(value,phi(current_index));
                    else phi(current_index)=max(value,phi(current_index));}}}}
}
//#####################################################################
// Adjust_Particle_Radii
//#####################################################################
template<class TV> void PARTICLE_LEVELSET<TV>::
Adjust_Particle_Radii()
{
    Adjust_Particle_Radii(negative_particles,-1);
    Adjust_Particle_Radii(positive_particles,1);
}
//#####################################################################
// Adjust_Particle_Radii
//#####################################################################
template<class TV> void PARTICLE_LEVELSET<TV>::
Adjust_Particle_Radii(T_ARRAYS_PARTICLE_LEVELSET_PARTICLES& particles,const int sign)
{
    for(CELL_ITERATOR iterator(levelset.grid);iterator.Valid();iterator.Next()){const T_INDEX& index=iterator.Cell_Index();
        if(particles(index)) for(int k=particles(index)->array_collection->Size();k>=1;--k){
            // the new radius will be negative if the particle is on the wrong side of the interface
            T new_radius=sign*levelset.Extended_Phi(particles(index)->X(k));
            if(new_radius<minimum_particle_radius) particles(index)->radius(k)=minimum_particle_radius;
            else if(new_radius>maximum_particle_radius) particles(index)->radius(k)=maximum_particle_radius;
            else particles(index)->radius(k)=new_radius;}}
}
//#####################################################################
// Reseed_Particles
//#####################################################################
template<class TV> void PARTICLE_LEVELSET<TV>::
Reseed_Particles()
{
    for(CELL_ITERATOR iterator(levelset.grid);iterator.Valid();iterator.Next()){const T_INDEX& index=iterator.Cell_Index();
        if(positive_particles(index)){delete positive_particles(index);positive_particles(index)=0;}
        if(negative_particles(index)){delete negative_particles(index);negative_particles(index)=0;}}
    Add_Particles_To_Both_Sides_Of_Interface_Band();
    Attract_Particles_To_Interface_And_Adjust_Radii();
}
//#####################################################################
template class PARTICLE_LEVELSET<VECTOR<float,2> >;
template class PARTICLE_LEVELSET<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class PARTICLE_LEVELSET<VECTOR<double,2> >;
template class PARTICLE_LEVELSET<VECTOR<double,3> >;
#endif
