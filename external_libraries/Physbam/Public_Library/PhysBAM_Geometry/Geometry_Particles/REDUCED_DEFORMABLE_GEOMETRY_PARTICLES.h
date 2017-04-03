//#####################################################################
// Copyright 2013, Rahul Sheth.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class REDUCED_DEFORMABLE_GEOMETRY_PARTICLES
//#####################################################################
#ifndef __REDUCED_DEFORMABLE_GEOMETRY_PARTICLES__
#define __REDUCED_DEFORMABLE_GEOMETRY_PARTICLES__

#include <PhysBAM_Tools/Matrices/ROTATION.h>
#include <PhysBAM_Geometry/Geometry_Particles/GEOMETRY_PARTICLES.h>
#include <PhysBAM_Geometry/Solids_Geometry/REDUCED_DEFORMABLE_BASIS.h>
#include <PhysBAM_Geometry/Solids_Geometry/REDUCED_DEFORMABLE_GEOMETRY_STATE.h>

namespace PhysBAM{

template<class TV> class REDUCED_DEFORMABLE_REST_GEOMETRY;

template<class TV>
class REDUCED_DEFORMABLE_GEOMETRY_PARTICLES : public CLONEABLE<REDUCED_DEFORMABLE_GEOMETRY_PARTICLES<TV>,GEOMETRY_PARTICLES<TV> >
{
    typedef typename TV::SCALAR T;
    typedef typename TV::SPIN T_SPIN;
    typedef CLONEABLE<REDUCED_DEFORMABLE_GEOMETRY_PARTICLES<TV>,GEOMETRY_PARTICLES<TV> > BASE;
public:
    using BASE::array_collection;
    //t and t_dot
    using BASE::X;using BASE::V;
    //R and R_dot
    ARRAY_VIEW<ROTATION<TV> > rotation;
    ARRAY_VIEW<T_SPIN> angular_velocity;
    ARRAY_VIEW<T_SPIN> angular_acceleration;
    //S, q, and q_dot
    ARRAY_VIEW<REDUCED_DEFORMABLE_BASIS<TV>*> reduced_basis;
    ARRAY_VIEW<VECTOR_ND<T> > reduced_displacements;
    ARRAY_VIEW<VECTOR_ND<T> > reduced_velocities;
    //Rest Geometry - x_o
    ARRAY_VIEW<REDUCED_DEFORMABLE_REST_GEOMETRY<TV>*> rest_geometry;
    //Body-specific basis scaling factor
    ARRAY_VIEW<T> scaling_factors;
    //Lazy full displacements and velocities
    ARRAY_VIEW<const LAZY_MATRIX_VECTOR_PRODUCT<MATRIX_MXN<T>,VECTOR_ND<T> >*> full_displacements;
    ARRAY_VIEW<const LAZY_MATRIX_VECTOR_PRODUCT<MATRIX_MXN<T>,VECTOR_ND<T> >*> full_velocities;
    bool scaling_needed;

    REDUCED_DEFORMABLE_GEOMETRY_PARTICLES(ARRAY_COLLECTION* array_collection_in,bool initialize=true)
        :rotation(0,0),angular_velocity(0,0),angular_acceleration(0,0),reduced_basis(0,0),reduced_displacements(0,0),reduced_velocities(0,0),rest_geometry(0,0),scaling_factors(0,0),full_displacements(0,0),full_velocities(0,0),scaling_needed(false)
    {if(initialize){if(array_collection){delete array_collection;}array_collection=array_collection_in;Initialize_Array_Collection();}}

    REDUCED_DEFORMABLE_GEOMETRY_PARTICLES(bool initialize=true)
        :rotation(0,0),angular_velocity(0,0),angular_acceleration(0,0),reduced_basis(0,0),reduced_displacements(0,0),reduced_velocities(0,0),rest_geometry(0,0),scaling_factors(0,0),full_displacements(0,0),full_velocities(0,0),scaling_needed(false)
    {if(initialize){Initialize_Array_Collection();}}

    void Initialize_Array_Collection()
    {BASE::Initialize_Array_Collection(); //X
     this->Store_Velocity(); //always store V
     array_collection->Add_Array(ATTRIBUTE_ID_ROTATION,&rotation);
     array_collection->Add_Array(ATTRIBUTE_ID_ANGULAR_VELOCITY,&angular_velocity);
     array_collection->Add_Array(ATTRIBUTE_ID_ANGULAR_ACCELERATION,&angular_acceleration);
     array_collection->Add_Array(ATTRIBUTE_ID_REDUCED_BASIS,&reduced_basis);
     array_collection->Add_Array(ATTRIBUTE_ID_REDUCED_DISPLACEMENTS,&reduced_displacements);
     array_collection->Add_Array(ATTRIBUTE_ID_REDUCED_VELOCITIES,&reduced_velocities);
     array_collection->Add_Array(ATTRIBUTE_ID_REDUCED_REST_GEOMETRY,&rest_geometry);
     array_collection->Add_Array(ATTRIBUTE_ID_REDUCED_SCALING_FACTOR,&scaling_factors);
     array_collection->Add_Array(ATTRIBUTE_ID_REDUCED_LAZY_DISPLACEMENTS,&full_displacements);
     array_collection->Add_Array(ATTRIBUTE_ID_REDUCED_LAZY_VELOCITIES,&full_velocities);}
    
    ~REDUCED_DEFORMABLE_GEOMETRY_PARTICLES()
    {Clean_Memory();}
//#####################################################################
    void Resize(const int new_size);
    virtual void Remove_Geometry(const int p);
    inline void Remove_Reduced_Basis(const int p);
    void Clean_Memory();
    void Delete_All_Particles();
    void Set_State(const REDUCED_DEFORMABLE_GEOMETRY_STATE<TV>& state_in,const int index);
    void Update_Reduced_Basis(const MATRIX_MXN<T>& basis_in,const int index);
    virtual TV Object_Space_Position_At_Particle(const int body_index,const int particle_index);
    virtual TV Object_Space_Velocity_At_Particle(const int body_index,const int particle_index);
    virtual TV World_Space_Position_At_Particle(const int body_index,const int particle_index);
    virtual TV World_Space_Velocity_At_Particle(const int body_index,const int particle_index);
    virtual TV World_Space_Moment_Arm_At_Particle(const int body_index,const int particle_index);
    TV Object_Space_Position_At_Particle_Helper(const int body_index,const int particle_index,const T mass);
    TV Object_Space_Velocity_At_Particle_Helper(const int body_index,const int particle_index,const T mass);
    TV World_Space_Position_At_Particle_Helper(const int body_index,const int particle_index,const T mass);
    TV World_Space_Velocity_At_Particle_Helper(const int body_index,const int particle_index,const T mass);
    TV World_Space_Moment_Arm_At_Particle_Helper(const int body_index,const int particle_index,const T mass);
//#####################################################################
};
}
#endif
