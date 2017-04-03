//#####################################################################
// Copyright 2013, Wenlong Lu, Yue Yu.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ADVECTION_SEMI_LAGRANGIAN_SURFACE_MESH
//#####################################################################
#ifndef __ADVECTION_SEMI_LAGRANGIAN_SURFACE_MESH__
#define __ADVECTION_SEMI_LAGRANGIAN_SURFACE_MESH__

#include <PhysBAM_Geometry/Topology/SEGMENT_MESH.h>
#include <PhysBAM_Geometry/Topology/TRIANGLE_MESH.h>
#include <PhysBAM_Geometry/Solids_Geometry/RIGID_GEOMETRY.h>
#include <PhysBAM_Tools/Arrays/ARRAY_VIEW.h>
#include <PhysBAM_Tools/Data_Structures/TRIPLE.h>
#include <PhysBAM_Tools/Random_Numbers/RANDOM_NUMBERS.h>

namespace PhysBAM{

template<class TV,class T2>
class ADVECTION_SEMI_LAGRANGIAN_SURFACE_MESH
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::dimension> TV_INT;
    typedef typename IF<TV::m==3,TRIANGLE_MESH,SEGMENT_MESH>::TYPE SURFACE_MESH_TYPE;
public:
    ARRAY_VIEW<TV>* shadow_X;
    bool rotate_before_barycentric_interpolation;

    ADVECTION_SEMI_LAGRANGIAN_SURFACE_MESH(bool rotate=true):shadow_X(0),rotate_before_barycentric_interpolation(rotate){}


    template<int d> typename ENABLE_IF<d==2,void>::TYPE
    Update_Advection_Equation_Node(const SURFACE_MESH_TYPE& mesh,
            ARRAY_VIEW<T2>& Z,const ARRAY_VIEW<T2>& Z_ghost,
            const ARRAY_VIEW<VECTOR<T,d> >& X,const ARRAY_VIEW<VECTOR<T,d> >& V,const ARRAY_VIEW<VECTOR<T,d> >& normal,
            const RIGID_GEOMETRY<TV>& rigid_geometry,const FRAME<TV>& last_step_frame,
            const T dt,const T time,ARRAY<ARRAY<TRIPLE<int,T,T2> > >* weights_to_cell=0,
            bool forward=false);
    template<int d> typename ENABLE_IF<d==3,void>::TYPE
    Update_Advection_Equation_Node(const SURFACE_MESH_TYPE& mesh,
            ARRAY_VIEW<T2>& Z,const ARRAY_VIEW<T2>& Z_ghost,
            const ARRAY_VIEW<VECTOR<T,d> >& X,const ARRAY_VIEW<VECTOR<T,d> >& V,const ARRAY_VIEW<VECTOR<T,d> >& normal,
            const RIGID_GEOMETRY<TV>& rigid_geometry,const FRAME<TV>& last_step_frame,
            const T dt,const T time,ARRAY<ARRAY<TRIPLE<int,T,T2> > >* weights_to_cell=0,
            bool forward=false);
    static bool Intersection_Segments(int plane,const TV& p1,const TV& p2,const TV& p3,const TV& p4,T& s,T& t,bool& too_short);
    static bool Intersection_Segments(const TV& normal,const TV& p1,const TV& p2,const TV& p3,const TV& p4,T& s,T& t,bool& too_short);
private:
    static T tolerance;
    RANDOM_NUMBERS<T> rng;
};

template<class TV,class T2>
typename TV::SCALAR ADVECTION_SEMI_LAGRANGIAN_SURFACE_MESH<TV,T2>::tolerance=(T)1e-7;
}
#endif
