//#####################################################################
// Copyright (c) 2014, Mridul Aanjaneya.
// Distributed under the FreeBSD license (see license.txt)
//#####################################################################
// Class READ_WRITE_PARTICLES
//#####################################################################
#include <PhysBAM_Tools/Arrays/ARRAY_COLLECTION.h>
#include <PhysBAM_Tools/Read_Write/Arrays/READ_WRITE_ARRAY_COLLECTION.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Geometry/Geometry_Particles/GEOMETRY_PARTICLES_FORWARD.h>
#include <PhysBAM_Geometry/Geometry_Particles/REGISTER_GEOMETRY_READ_WRITE.h>
#include <SPGrid_Fluids/Particles/PARTICLES_FORWARD.h>
#include <SPGrid_Fluids/Read_Write/Particles/READ_WRITE_PARTICLES.h>
namespace PhysBAM{
    
void Initialize_Particles()
{
    static bool done=false;if(done) return;done=true;
    Register_Attribute_Name(ATTRIBUTE_ID_VORTICITY,"vorticity");

    #define READ_WRITE_HELPER(RW) \
        Read_Write<ARRAY_COLLECTION,RW>::Register_Read_Write<unsigned short>();

    READ_WRITE_HELPER(float);
    #ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
    READ_WRITE_HELPER(double);
    #endif

    Initialize_Geometry_Particle();
}
}
