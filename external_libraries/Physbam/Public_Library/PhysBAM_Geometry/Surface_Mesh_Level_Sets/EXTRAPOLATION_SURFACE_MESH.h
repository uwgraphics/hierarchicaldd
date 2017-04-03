//#####################################################################
// Copyright 2013, Yue Yu.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class EXTRAPOLATION_SURFACE_MESH  
//##################################################################### 
//
// Extrapolates the values of u in the normal direction from phi <= 0 to phi > 0,
// overwriting u where phi > 0, or possibly some real cells when using the Isobaric Fix.
//
//#####################################################################
#ifndef __EXTRAPOLATION_SURFACE_MESH__
#define __EXTRAPOLATION_SURFACE_MESH__

#include <PhysBAM_Tools/Arrays/ARRAYS_FORWARD.h>
#include <PhysBAM_Geometry/Topology/SEGMENT_MESH.h>
#include <PhysBAM_Geometry/Topology/TRIANGLE_MESH.h>
#include <PhysBAM_Geometry/Solids_Geometry/RIGID_GEOMETRY.h>
namespace PhysBAM{

template<class TV,class T2>
class EXTRAPOLATION_SURFACE_MESH
{
    typedef typename TV::SCALAR T;
    typedef typename TV::template REBIND<int>::TYPE TV_INT;
    typedef typename IF<TV::m==2,SEGMENT_MESH,TRIANGLE_MESH>::TYPE SURFACE_MESH;
    const SURFACE_MESH& mesh;
public:
    const ARRAY_VIEW<T>& phi;
    const ARRAY_VIEW<TV>& X;
    const ARRAY_VIEW<TV>& normal;
    ARRAY_VIEW<T2>& u;
    T max_dX;
    ARRAY<T> local_average_dX;
    T band_width; // band for extrapolation near the interface
public:

    EXTRAPOLATION_SURFACE_MESH(const SURFACE_MESH& mesh_input,const ARRAY_VIEW<T>& phi_input,const ARRAY_VIEW<TV>& X_input,const ARRAY_VIEW<TV>& normal_input,ARRAY_VIEW<T2>& u_input);
    ~EXTRAPOLATION_SURFACE_MESH();
    void Set_Band_Width(const T number_of_cells=(T)3)
    {band_width=number_of_cells*max_dX;}

//#####################################################################
    void Extrapolate(ARRAY<ARRAY<PAIR<int,T> > >* weights=NULL);
private:
    void Initialize(ARRAY<bool>& done,ARRAY<int>& close_k,ARRAY<int>& heap,int& heap_length);
    void Update_Close_Point(ARRAY<bool>& done,int index,ARRAY<int>& closest_node_indices);
//#####################################################################
}; 
}
#endif
