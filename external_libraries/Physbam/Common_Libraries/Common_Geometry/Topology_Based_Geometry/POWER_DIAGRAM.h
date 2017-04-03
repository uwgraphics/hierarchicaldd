//#####################################################################
// Copyright 2015, Mridul Aanjaneya.
// This file is part of PhysBAM whose distribution is governed by the license 
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class POWER_DIAGRAM
//#####################################################################
#ifndef __POWER_DIAGRAM__
#define __POWER_DIAGRAM__

#include <Common_Geometry/Topology_Based_Geometry/VORONOI_DIAGRAM.h>
#include <PhysBAM_Tools/Arrays/ARRAY.h>

namespace PhysBAM{
template<class TV>
class POWER_DIAGRAM: public VORONOI_DIAGRAM<TV>
{
  public:
    typedef VORONOI_DIAGRAM<TV> BASE;
    typedef typename TV::SCALAR T;
    using BASE::center;using BASE::face_vertices;using BASE::neighbors;

    ARRAY<TV> centroid;

    POWER_DIAGRAM() {}

    POWER_DIAGRAM(const int number_of_sites)
        :BASE(number_of_sites)
    {
        centroid.Resize(number_of_sites);
    }

    ~POWER_DIAGRAM() {}

    void Clean_Memory()
    {
        centroid.Clean_Memory();
    }
};
}
#endif
