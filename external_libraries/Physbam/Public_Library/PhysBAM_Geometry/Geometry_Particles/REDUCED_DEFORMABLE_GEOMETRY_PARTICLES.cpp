//#####################################################################
// Copyright 2013, Rahul Sheth.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Tools/Vectors/VECTOR_EXPRESSION.h>
#include <PhysBAM_Geometry/Geometry_Particles/REDUCED_DEFORMABLE_GEOMETRY_PARTICLES.h>
#include <PhysBAM_Geometry/Solids_Geometry/REDUCED_DEFORMABLE_REST_GEOMETRY.h>
#include <PhysBAM_Geometry/Solids_Geometry/REDUCED_DEFORMABLE_BASIS.h>
namespace PhysBAM{
//#####################################################################
// Function Resize
//#####################################################################
template<class TV> void REDUCED_DEFORMABLE_GEOMETRY_PARTICLES<TV>::
Resize(const int new_size)
{
    for(int p=new_size+1;p<=array_collection->Size();p++){
        if(rest_geometry(p)) Remove_Geometry(p);
        if(reduced_basis(p)) Remove_Reduced_Basis(p);
        if(full_displacements(p)) delete full_displacements(p);
        if(full_velocities(p)) delete full_velocities(p);}
    array_collection->Resize(new_size);
}
//#####################################################################
// Function Remove_Geometry
//#####################################################################
template<class TV> void REDUCED_DEFORMABLE_GEOMETRY_PARTICLES<TV>::
Remove_Geometry(const int p)
{
    delete rest_geometry(p);
}
//#####################################################################
// Function Remove_Reduced_Basis
//#####################################################################
template<class TV> inline void REDUCED_DEFORMABLE_GEOMETRY_PARTICLES<TV>::
Remove_Reduced_Basis(const int p)
{
    delete reduced_basis(p);
}
//#####################################################################
// Function Clean_Memory
//#####################################################################
template<class TV> void REDUCED_DEFORMABLE_GEOMETRY_PARTICLES<TV>::
Clean_Memory()
{
    for(int p=1;p<=array_collection->Size();p++){
        if(rest_geometry(p)) Remove_Geometry(p);
        if(reduced_basis(p)) Remove_Reduced_Basis(p);}
    array_collection->Clean_Memory();
}
//#####################################################################
// Function Delete_All_Particles
//#####################################################################
template<class TV> void REDUCED_DEFORMABLE_GEOMETRY_PARTICLES<TV>::
Delete_All_Particles()
{
    for(int p=1;p<=array_collection->Size();p++){
        if(rest_geometry(p)) Remove_Geometry(p);
        if(reduced_basis(p)) Remove_Reduced_Basis(p);}
    array_collection->Delete_All_Elements();
}
//#####################################################################
// Function Set_State
//#####################################################################
template<class TV> void REDUCED_DEFORMABLE_GEOMETRY_PARTICLES<TV>::
Set_State(const REDUCED_DEFORMABLE_GEOMETRY_STATE<TV>& state_in,const int index)
{
    state_in.Get_Translation(this->X(index));
    state_in.Get_Translational_Velocity(this->V(index));
    state_in.Get_Rotation(this->rotation(index));
    state_in.Get_Angular_Velocity(this->angular_velocity(index));
    state_in.Get_Reduced_Displacements(this->reduced_displacements(index));
    state_in.Get_Reduced_Velocities(this->reduced_velocities(index));
}
//#####################################################################
// Function Update_Reduced_Basis
//#####################################################################
template<class TV> void REDUCED_DEFORMABLE_GEOMETRY_PARTICLES<TV>::
Update_Reduced_Basis(const MATRIX_MXN<T>& basis_in,const int index)
{
    if(reduced_basis(index)) delete reduced_basis(index);
    reduced_basis(index)=new REDUCED_DEFORMABLE_BASIS<TV>(basis_in);
    this->reduced_displacements(index).Resize(basis_in.Columns());
    this->reduced_velocities(index).Resize(basis_in.Columns());
}
//#####################################################################
// Function Object_Space_Position_At_Particle
//#####################################################################
template<class TV> TV REDUCED_DEFORMABLE_GEOMETRY_PARTICLES<TV>::
Object_Space_Position_At_Particle(const int body_index,const int particle_index)
{
    return Object_Space_Position_At_Particle_Helper(body_index,particle_index,(T)1.0);
}
//#####################################################################
// Function Object_Space_Velocity_At_Particle
//#####################################################################
template<class TV> TV REDUCED_DEFORMABLE_GEOMETRY_PARTICLES<TV>::
Object_Space_Velocity_At_Particle(const int body_index,const int particle_index)
{
    return Object_Space_Velocity_At_Particle_Helper(body_index,particle_index,(T)1.0);
}
//#####################################################################
// Function World_Space_Position_At_Particle
//#####################################################################
template<class TV> TV REDUCED_DEFORMABLE_GEOMETRY_PARTICLES<TV>::
World_Space_Position_At_Particle(const int body_index,const int particle_index)
{
    return World_Space_Position_At_Particle_Helper(body_index,particle_index,(T)1.0);
}
//#####################################################################
// Function World_Space_Velocity_At_Particle
//#####################################################################
template<class TV> TV REDUCED_DEFORMABLE_GEOMETRY_PARTICLES<TV>::
World_Space_Velocity_At_Particle(const int body_index,const int particle_index)
{
    return World_Space_Velocity_At_Particle_Helper(body_index,particle_index,(T)1.0);
}
//#####################################################################
// Function World_Space_Moment_Arm_At_Particle
//#####################################################################
template<class TV> TV REDUCED_DEFORMABLE_GEOMETRY_PARTICLES<TV>::
World_Space_Moment_Arm_At_Particle(const int body_index,const int particle_index)
{
    return World_Space_Moment_Arm_At_Particle_Helper(body_index,particle_index,(T)1.0);
}
//#####################################################################
// Function Object_Space_Position_At_Particle_Helper
//#####################################################################
template<class TV> TV REDUCED_DEFORMABLE_GEOMETRY_PARTICLES<TV>::
Object_Space_Position_At_Particle_Helper(const int body_index,const int particle_index,const T mass)
{
    int offset=rest_geometry(body_index)->rest_particle_indices(1)-1;
    if(full_displacements(body_index)==NULL){full_displacements(body_index)=new LAZY_MATRIX_VECTOR_PRODUCT<MATRIX_MXN<T>,VECTOR_ND<T> >(reduced_basis(body_index)->operator*(reduced_displacements(body_index)));}
    TV position;
    for(int i=1;i<=TV::dimension;i++) {
        position(i)=(scaling_needed?scaling_factors(body_index)*sqrt((T)reduced_basis(body_index)->reduced_basis.Rows())/sqrt(mass):(T)1.0)*full_displacements(body_index)->operator()((particle_index-1)*TV::dimension+i); }
    position += rest_geometry(body_index)->rest_geometry_particles->X(particle_index+offset);
    return position;
}
//#####################################################################
// Function Object_Space_Velocity_At_Particle_Helper
//#####################################################################
template<class TV> TV REDUCED_DEFORMABLE_GEOMETRY_PARTICLES<TV>::
Object_Space_Velocity_At_Particle_Helper(const int body_index,const int particle_index,const T mass)
{
    if(full_velocities(body_index)==NULL){full_velocities(body_index)=new LAZY_MATRIX_VECTOR_PRODUCT<MATRIX_MXN<T>,VECTOR_ND<T> >(reduced_basis(body_index)->operator*(reduced_velocities(body_index)));}
    TV velocity;
    for(int i=1;i<=TV::dimension;i++) {
        velocity(i)=(scaling_needed?scaling_factors(body_index)*sqrt((T)reduced_basis(body_index)->reduced_basis.Rows())/sqrt(mass):(T)1.0)*full_velocities(body_index)->operator()((particle_index-1)*TV::dimension+i); }
    return velocity;
}
//#####################################################################
// Function World_Space_Position_At_Particle_Helper
//#####################################################################
template<class TV> TV REDUCED_DEFORMABLE_GEOMETRY_PARTICLES<TV>::
World_Space_Position_At_Particle_Helper(const int body_index,const int particle_index,const T mass)
{
    return FRAME<TV>(X(body_index),rotation(body_index))*Object_Space_Position_At_Particle_Helper(body_index,particle_index,mass);
}
//#####################################################################
// Function World_Space_Velocity_At_Particle_Helper
//#####################################################################
template<class TV> TV REDUCED_DEFORMABLE_GEOMETRY_PARTICLES<TV>::
World_Space_Velocity_At_Particle_Helper(const int body_index,const int particle_index,const T mass)
{
    return rotation(body_index).Rotate(Object_Space_Velocity_At_Particle_Helper(body_index,particle_index,mass)) + V(body_index) \
        + TV::Cross_Product(angular_velocity(body_index),World_Space_Moment_Arm_At_Particle_Helper(body_index,particle_index,mass));
}
//#####################################################################
// Function World_Space_Moment_Arm_At_Particle_Helper
//#####################################################################
template<class TV> TV REDUCED_DEFORMABLE_GEOMETRY_PARTICLES<TV>::
World_Space_Moment_Arm_At_Particle_Helper(const int body_index,const int particle_index,const T mass)
{
    return FRAME<TV>(TV(),rotation(body_index))*Object_Space_Position_At_Particle_Helper(body_index,particle_index,mass);
}
//#####################################################################
template class REDUCED_DEFORMABLE_GEOMETRY_PARTICLES<VECTOR<float,1> >;
template class REDUCED_DEFORMABLE_GEOMETRY_PARTICLES<VECTOR<float,2> >;
template class REDUCED_DEFORMABLE_GEOMETRY_PARTICLES<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class REDUCED_DEFORMABLE_GEOMETRY_PARTICLES<VECTOR<double,1> >;
template class REDUCED_DEFORMABLE_GEOMETRY_PARTICLES<VECTOR<double,2> >;
template class REDUCED_DEFORMABLE_GEOMETRY_PARTICLES<VECTOR<double,3> >;
#endif
}
