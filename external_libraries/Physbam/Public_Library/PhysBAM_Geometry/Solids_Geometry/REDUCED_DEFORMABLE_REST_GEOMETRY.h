//#####################################################################
// Copyright 2013, Rahul Sheth.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class REDUCED_DEFORMABLE_REST_GEOMETRY
//#####################################################################
#ifndef __REDUCED_DEFORMABLE_REST_GEOMETRY__
#define __REDUCED_DEFORMABLE_REST_GEOMETRY__

#include <PhysBAM_Tools/Utilities/Find_Type.h>
#include <PhysBAM_Tools/Utilities/TYPE_UTILITIES.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Geometry/Implicit_Objects/IMPLICIT_OBJECT.h>
#include <PhysBAM_Geometry/Solids_Geometry/REDUCED_DEFORMABLE_GEOMETRY_COLLECTION.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/STRUCTURE.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TOPOLOGY_BASED_SIMPLEX_POLICY.h>

namespace PhysBAM{

template<class TV>
class REDUCED_DEFORMABLE_REST_GEOMETRY : public NONCOPYABLE
{
    typedef typename TV::SCALAR T;
    typedef typename TV::SPIN T_SPIN;
    typedef typename TOPOLOGY_BASED_SIMPLEX_POLICY<TV,TV::m-1>::OBJECT T_SIMPLICIAL_OBJECT;

  public:
    REDUCED_DEFORMABLE_GEOMETRY_COLLECTION<TV>* const reduced_deformable_geometry_collection;
    GEOMETRY_PARTICLES<TV>* rest_geometry_particles; //shared particles
    T_SIMPLICIAL_OBJECT* simplicial_object;
    IMPLICIT_OBJECT<TV>* implicit_object;
    int parent_particle_index; //parent particle index
    ARRAY<STRUCTURE<TV>*> structures;
    ARRAY<int> rest_particle_indices; //indices of the particles that this body owns
    bool own_structures;

    REDUCED_DEFORMABLE_REST_GEOMETRY(REDUCED_DEFORMABLE_GEOMETRY_COLLECTION<TV>* const reduced_deformable_geometry_collection_in,int parent_particle_index_in=0);
  protected:
    REDUCED_DEFORMABLE_REST_GEOMETRY(REDUCED_DEFORMABLE_REST_GEOMETRY<TV>* rest_geometry_in);
  public:
    virtual ~REDUCED_DEFORMABLE_REST_GEOMETRY();

    template<class T_STRUCTURE>
    T_STRUCTURE Find_Structure(const int index=1) const
    {return Find_Type<T_STRUCTURE>(structures,index);}
//#####################################################################
    template<int d> typename ENABLE_IF<d==1>::TYPE Add_Structure(STRUCTURE<VECTOR<T,d> >& structure);
    template<int d> typename ENABLE_IF<d==2>::TYPE Add_Structure(STRUCTURE<VECTOR<T,d> >& structure);
    template<int d> typename ENABLE_IF<d==3>::TYPE Add_Structure(STRUCTURE<VECTOR<T,d> >& structure);
    void Remove_Structure(STRUCTURE<TV>* structure);
    void Update_Rest_Particle_Indices();
//#####################################################################
};

}
#endif
