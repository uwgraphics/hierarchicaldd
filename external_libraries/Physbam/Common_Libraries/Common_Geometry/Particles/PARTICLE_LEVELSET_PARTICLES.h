//#####################################################################
// Copyright 2016, Mridul Aanjaneya.
// This file is part of PhysBAM whose distribution is governed by the license 
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class PARTICLE_LEVELSET_PARTICLES
//#####################################################################
#ifndef __PARTICLE_LEVELSET_PARTICLES__
#define __PARTICLE_LEVELSET_PARTICLES__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Clone/CLONEABLE.h>
#include <PhysBAM_Tools/Point_Clouds/POINT_CLOUD.h>

namespace PhysBAM{
template<class TV>
class PARTICLE_LEVELSET_PARTICLES: public CLONEABLE<PARTICLE_LEVELSET_PARTICLES<TV>,POINT_CLOUD<TV> >
{
    typedef typename TV::SCALAR T;
    typedef CLONEABLE<PARTICLE_LEVELSET_PARTICLES<TV>,POINT_CLOUD<TV> > BASE;
  public:
    using BASE::array_collection;using BASE::X;

    ARRAY_VIEW<T> radius;
    
    PARTICLE_LEVELSET_PARTICLES(); 
    ~PARTICLE_LEVELSET_PARTICLES() {}

//#####################################################################
    void Resize(const int new_size);
//#####################################################################
};
}
#endif
