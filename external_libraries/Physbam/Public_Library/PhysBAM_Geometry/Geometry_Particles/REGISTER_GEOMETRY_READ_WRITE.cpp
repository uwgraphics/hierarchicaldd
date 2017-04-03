//#####################################################################
// Copyright 2009-2013, Geoffrey Irving, Michael Lentine, Rahul Sheth.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class REGISTER_GEOMETRY_READ_WRITE
//#####################################################################
#include <PhysBAM_Tools/Arrays/ARRAY_COLLECTION.h>
#include <PhysBAM_Tools/Matrices/FRAME.h>
#include <PhysBAM_Tools/Point_Clouds/POINT_CLOUD_FORWARD.h>
#include <PhysBAM_Tools/Read_Write/Arrays/READ_WRITE_ARRAY_COLLECTION.h>
#include <PhysBAM_Tools/Vectors/TWIST.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Geometry/Geometry_Particles/GEOMETRY_PARTICLES_FORWARD.h>
#include <PhysBAM_Geometry/Solids_Geometry/RIGID_GEOMETRY.h>
#include <PhysBAM_Geometry/Solids_Geometry/REDUCED_DEFORMABLE_BASIS.h>
#include <PhysBAM_Geometry/Solids_Geometry/REDUCED_DEFORMABLE_REST_GEOMETRY.h>
namespace PhysBAM{
bool Register_Implicit_Object_Transformed(); // TODO(jontg): Move these out of here
bool Register_Multibody_Levelset_Implicit_Object();
bool Register_Analytic_Implicit_Object();
void Register_Attribute_Name(const ATTRIBUTE_ID id,const char* name);
void Register_Free_Particles();

void Initialize_Geometry_Particle()
{
    static bool done=false;if(done) return;done=true;
    Register_Attribute_Name(ATTRIBUTE_ID_ROTATION,"rotation");
    Register_Attribute_Name(ATTRIBUTE_ID_RIGID_GEOMETRY,"rigid_geometry");
    Register_Attribute_Name(ATTRIBUTE_ID_X,"X");
    Register_Attribute_Name(ATTRIBUTE_ID_V,"V");
    Register_Attribute_Name(ATTRIBUTE_ID_ANGULAR_VELOCITY,"angular_velocity");
    Register_Attribute_Name(ATTRIBUTE_ID_STRUCTURE_IDS,"structure_ids");
    Register_Attribute_Name(ATTRIBUTE_ID_ID,"id");
    Register_Attribute_Name(ATTRIBUTE_ID_COLOR,"color");
    Register_Attribute_Name(ATTRIBUTE_ID_RADIUS,"radius");
    Register_Attribute_Name(ATTRIBUTE_ID_THICKNESS,"thickness");
    Register_Attribute_Name(ATTRIBUTE_ID_OLD_ID,"old_id");
    Register_Attribute_Name(ATTRIBUTE_ID_GUID,"guid");
    Register_Attribute_Name(ATTRIBUTE_ID_TID,"tid");
    Register_Attribute_Name(ATTRIBUTE_ID_LOCAL_INDEX,"local_index");
    Register_Attribute_Name(ATTRIBUTE_ID_THREAD_BOUNDARY,"thread_boundary");
    Register_Attribute_Name(ATTRIBUTE_ID_THREAD_BOUNDARY_VECTOR,"thread boundary vector");
    Register_Attribute_Name(ATTRIBUTE_ID_INCOMPRESSIBILITY,"incompressibility");
    Register_Attribute_Name(ATTRIBUTE_ID_DISTANCE_TO_RIM,"distance to rim");
    Register_Attribute_Name(ATTRIBUTE_ID_DISTANCE_TO_VOLUME,"distance to volume");
    Register_Attribute_Name(ATTRIBUTE_ID_BUBBLE,"bubble");
    Register_Attribute_Name(ATTRIBUTE_ID_SMALL_REGION_ID,"small region");
    Register_Attribute_Name(ATTRIBUTE_ID_SURFACE_TENSION_FORCE,"surface tension force");
    Register_Attribute_Name(ATTRIBUTE_ID_PRESSURE_FORCE,"pressure force");
    Register_Attribute_Name(ATTRIBUTE_ID_BODY_FORCE,"body force");
    Register_Attribute_Name(ATTRIBUTE_ID_VISCOSITY_FORCE,"viscosity force");
    Register_Attribute_Name(ATTRIBUTE_ID_AREA_WEIGHTED_NORMAL,"area weighted normal");
    Register_Attribute_Name(ATTRIBUTE_ID_DIMENSION,"dimension");
    Register_Attribute_Name(ATTRIBUTE_ID_DIMENSIONAL_MASS,"dimensional mass");
    Register_Attribute_Name(ATTRIBUTE_ID_DIMENSIONAL_SIZE,"dimensional size");
    Register_Attribute_Name(ATTRIBUTE_ID_DISPLAY_SIZE,"display_size");
    Register_Attribute_Name(ATTRIBUTE_ID_BRADIUS,"bradius");
    Register_Attribute_Name(ATTRIBUTE_ID_BMASS,"bmass");
    Register_Attribute_Name(ATTRIBUTE_ID_PRESSURE,"pressure");
    Register_Attribute_Name(ATTRIBUTE_ID_SIGNED_DISTANCE,"signed distance");
    Register_Attribute_Name(ATTRIBUTE_ID_RADIAL_VELOCITY,"radial velocity");
    Register_Attribute_Name(ATTRIBUTE_ID_WEIGHTS,"weights");
    Register_Attribute_Name(ATTRIBUTE_ID_REDUCED_BASIS,"reduced basis");
    Register_Attribute_Name(ATTRIBUTE_ID_REDUCED_DISPLACEMENTS,"reduced displacements");
    Register_Attribute_Name(ATTRIBUTE_ID_REDUCED_VELOCITIES,"reduced velocities");
    Register_Attribute_Name(ATTRIBUTE_ID_REDUCED_REST_GEOMETRY,"reduced rest geometry");
    Register_Attribute_Name(ATTRIBUTE_ID_RGB_COLOR,"rgb color");

    #define READ_WRITE_VECTOR_HELPER(T,RW,d) \
        Read_Write<ARRAY_COLLECTION,RW>::Register_Read_Write<VECTOR<T,d> >(); \
        Read_Write<ARRAY_COLLECTION,RW>::Register_Read_Write<ARRAY<PAIR<VECTOR<int,d>,T> > >(); \
        Read_Write<ARRAY_COLLECTION,RW>::Register_Read_Write<ROTATION<VECTOR<T,d> > >(); \
        Read_Write<ARRAY_COLLECTION,RW>::Register_Read_Write<TWIST<VECTOR<T,d> > >(); \
        Read_Write<ARRAY_COLLECTION,RW>::Register_Read_Write<RIGID_GEOMETRY<VECTOR<T,d> >*>(); \
        Read_Write<ARRAY_COLLECTION,RW>::Register_Read_Write<REDUCED_DEFORMABLE_REST_GEOMETRY<VECTOR<T,d> >*>(); \
        Read_Write<ARRAY_COLLECTION,RW>::Register_Read_Write<REDUCED_DEFORMABLE_BASIS<VECTOR<T,d> >*>();

    #define READ_WRITE_SCALAR_HELPER(T,RW) \
        Read_Write<ARRAY_COLLECTION,RW>::Register_Read_Write<T>(); \
        Read_Write<ARRAY_COLLECTION,RW>::Register_Read_Write<VECTOR<T,0> >(); \
        Read_Write<ARRAY_COLLECTION,RW>::Register_Read_Write<VECTOR_ND<T> >(); \
        READ_WRITE_VECTOR_HELPER(T,RW,1);READ_WRITE_VECTOR_HELPER(T,RW,2);READ_WRITE_VECTOR_HELPER(T,RW,3);

    #define READ_WRITE_HELPER(RW) \
        Read_Write<ARRAY_COLLECTION,RW>::Register_Read_Write<VECTOR<int,1> >(); \
        Read_Write<ARRAY_COLLECTION,RW>::Register_Read_Write<VECTOR<int,2> >(); \
        Read_Write<ARRAY_COLLECTION,RW>::Register_Read_Write<VECTOR<int,3> >(); \
        Read_Write<ARRAY_COLLECTION,RW>::Register_Read_Write<int>(); \
        Read_Write<ARRAY_COLLECTION,RW>::Register_Read_Write<VECTOR<bool,1> >(); \
        Read_Write<ARRAY_COLLECTION,RW>::Register_Read_Write<VECTOR<bool,2> >(); \
        Read_Write<ARRAY_COLLECTION,RW>::Register_Read_Write<VECTOR<bool,3> >(); \
        Read_Write<ARRAY_COLLECTION,RW>::Register_Read_Write<ARRAY<int> >();

    READ_WRITE_SCALAR_HELPER(float,float);
    READ_WRITE_SCALAR_HELPER(float,double);
    READ_WRITE_HELPER(float);
    #ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
    READ_WRITE_SCALAR_HELPER(double,float);
    READ_WRITE_SCALAR_HELPER(double,double);
    READ_WRITE_HELPER(double);
    #endif

    Register_Implicit_Object_Transformed();
    Register_Multibody_Levelset_Implicit_Object();
    Register_Analytic_Implicit_Object();
    Register_Free_Particles();
}
}
