//#####################################################################
// Copyright (c) 2004-2014, Mridul Aanjaneya, Ron Fedkiw, Geoffrey Irving, Michael Lentine, Andrew Selle.
// Distributed under the FreeBSD license (see license.txt)
//#####################################################################
// Class VORTEX_PARTICLES
//#####################################################################
#ifndef __VORTEX_PARTICLES__
#define __VORTEX_PARTICLES__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Clone/CLONEABLE.h>
#include <PhysBAM_Tools/Point_Clouds/POINT_CLOUD.h>
namespace PhysBAM{

template<class TV>
class VORTEX_PARTICLES:public CLONEABLE<VORTEX_PARTICLES<TV>,POINT_CLOUD<TV> >
{
    typedef typename TV::SCALAR T;
    typedef CLONEABLE<VORTEX_PARTICLES<TV>,POINT_CLOUD<TV> > BASE;
public:
    using BASE::array_collection;

    ARRAY_VIEW<typename TV::SPIN> vorticity;
    ARRAY_VIEW<T> radius;

    VORTEX_PARTICLES();
    ~VORTEX_PARTICLES() {}
};
}
#endif
