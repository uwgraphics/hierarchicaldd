//#####################################################################
// Copyright 2013, Linhai Qiu
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class FAST_MARCHING_METHOD_SURFACE_MESH  
//#####################################################################
#ifndef __FAST_MARCHING_METHOD_SURFACE_MESH__
#define __FAST_MARCHING_METHOD_SURFACE_MESH__

#include <PhysBAM_Tools/Arrays/ARRAYS_FORWARD.h>
#include <PhysBAM_Geometry/Solids_Geometry/RIGID_GEOMETRY.h>
#include <PhysBAM_Geometry/Level_Sets/FAST_MARCHING.h>

namespace PhysBAM{

template<class TV>
class FAST_MARCHING_METHOD_SURFACE_MESH:public NONCOPYABLE
{
    typedef typename TV::SCALAR T;
    typedef typename TV::template REBIND<int>::TYPE TV_INT;
    typedef typename TOPOLOGY_BASED_SIMPLEX_POLICY<TV,TV::m-1>::OBJECT T_SIMPLICIAL_OBJECT;
    const T_SIMPLICIAL_OBJECT& simplicial_object;
    const ARRAY_VIEW<TV>& normal;
public:
    FAST_MARCHING_METHOD_SURFACE_MESH(const T_SIMPLICIAL_OBJECT& simplicial_object_input,const ARRAY_VIEW<TV>& normal_input);
    ~FAST_MARCHING_METHOD_SURFACE_MESH();

//#####################################################################
    void Fast_Marching_Method(ARRAY_VIEW<T>& phi);
private:
    void Update_Or_Add_Neighbor(ARRAY_VIEW<T>& phi,ARRAY<bool>& done,ARRAY<int>& close_k,ARRAY<int>& heap,int& heap_length,const int& neighbor);
    void Initialize_Interface(ARRAY_VIEW<T>& phi,ARRAY<bool>& done,ARRAY<int>& close_k,ARRAY<int>& heap,int& heap_length);
    void Update_Close_Point(ARRAY_VIEW<T>& phi,ARRAY<bool>& done,const int& index);
    void Add_To_Initial(ARRAY<bool>& done,ARRAY<int>& close_k,const int& index);
//#####################################################################
};
}
#endif
