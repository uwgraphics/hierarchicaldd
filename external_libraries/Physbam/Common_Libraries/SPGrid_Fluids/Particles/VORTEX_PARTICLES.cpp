//#####################################################################
// Copyright (c) 2004-2014, Mridul Aanjaneya, Ron Fedkiw, Geoffrey Irving, Michael Lentine, Andrew Selle.
// Distributed under the FreeBSD license (see license.txt)
//#####################################################################
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <SPGrid_Fluids/Particles/PARTICLES_FORWARD.h>
#include <SPGrid_Fluids/Particles/VORTEX_PARTICLES.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> VORTEX_PARTICLES<TV>::
VORTEX_PARTICLES()
    :vorticity(0,0),radius(0,0)
{
    array_collection->Add_Array(ATTRIBUTE_ID_VORTICITY,&vorticity);
    array_collection->Add_Array(ATTRIBUTE_ID_RADIUS,&radius);
}
//#####################################################################
template class VORTEX_PARTICLES<VECTOR<float,1> >;
template class VORTEX_PARTICLES<VECTOR<float,2> >;
template class VORTEX_PARTICLES<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class VORTEX_PARTICLES<VECTOR<double,1> >;
template class VORTEX_PARTICLES<VECTOR<double,2> >;
template class VORTEX_PARTICLES<VECTOR<double,3> >;
#endif
