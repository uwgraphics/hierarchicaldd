//#####################################################################
// Copyright 2016, Mridul Aanjaneya.
// This file is part of PhysBAM whose distribution is governed by the license 
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Common_Geometry/Particles/PARTICLE_LEVELSET_PARTICLES.h>
#include <Common_Geometry/Particles/PARTICLE_LEVELSET_PARTICLES_FORWARD.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> PARTICLE_LEVELSET_PARTICLES<TV>::
PARTICLE_LEVELSET_PARTICLES()
    :radius(0,0)
{
    array_collection->Add_Array(ATTRIBUTE_ID_RADIUS,&radius);
}
//#####################################################################
// Resize
//#####################################################################
template<class TV> void PARTICLE_LEVELSET_PARTICLES<TV>::
Resize(const int new_size)
{
    array_collection->Resize(new_size);
}
//#####################################################################
template class PARTICLE_LEVELSET_PARTICLES<VECTOR<float,2> >;
template class PARTICLE_LEVELSET_PARTICLES<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class PARTICLE_LEVELSET_PARTICLES<VECTOR<double,2> >;
template class PARTICLE_LEVELSET_PARTICLES<VECTOR<double,3> >;
#endif
