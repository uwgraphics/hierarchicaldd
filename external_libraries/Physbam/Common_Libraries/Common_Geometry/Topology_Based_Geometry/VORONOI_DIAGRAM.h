//#####################################################################
// Copyright 2015, Mridul Aanjaneya.
// This file is part of PhysBAM whose distribution is governed by the license 
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class VORONOI_DIAGRAM
//#####################################################################
#ifndef __VORONOI_DIAGRAM__
#define __VORONOI_DIAGRAM__

#include <PhysBAM_Tools/Arrays/ARRAY.h>

namespace PhysBAM{
template<class TV>
class VORONOI_DIAGRAM
{
  public:
    typedef typename TV::SCALAR T;

    ARRAY<TV> center;
    ARRAY<ARRAY<ARRAY<TV> > > face_vertices;
    ARRAY<ARRAY<TV> > face_barycenters;
    ARRAY<ARRAY<T> > face_areas;
    ARRAY<ARRAY<int> > neighbors;
    ARRAY<T> volume;

    VORONOI_DIAGRAM() {}

    VORONOI_DIAGRAM(const int number_of_sites)
    {
        center.Resize(number_of_sites);
        face_vertices.Resize(number_of_sites);
        neighbors.Resize(number_of_sites);
        volume.Resize(number_of_sites);
    }

    ~VORONOI_DIAGRAM() {}

    void Clean_Memory()
    {
        center.Clean_Memory();
        face_vertices.Clean_Memory();
        neighbors.Clean_Memory();
        volume.Clean_Memory();
    }

//#####################################################################
    void Compute_Face_Barycenters();
//#####################################################################
};
}
#endif
