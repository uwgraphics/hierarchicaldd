//#####################################################################
// Copyright 2013, Yue Yu.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// namespace TESSELLATION
//##################################################################### 
#ifndef __CONICAL_FRUSTUM_TESSELLATION__
#define __CONICAL_FRUSTUM_TESSELLATION__
 
namespace PhysBAM{
template<class TV> class CONICAL_FRUSTUM;
template<class T> class TRIANGULATED_SURFACE;
template<class T> class SEGMENTED_CURVE_2D;
template<class T,int d> class VECTOR;

namespace TESSELLATION{
//#####################################################################
    template<class T> TRIANGULATED_SURFACE<T>* Generate_Triangles(const CONICAL_FRUSTUM<VECTOR<T,3> >& cylinder,const int radius_height=4,const int resolution_radius=16);
    template<class T> SEGMENTED_CURVE_2D<T>* Generate_Triangles(const CONICAL_FRUSTUM<VECTOR<T,2> >& cylinder,const int radius_height=-1,const int resolution_radius=-1);
//#####################################################################
}
}
#endif
