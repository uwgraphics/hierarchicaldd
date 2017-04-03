//#####################################################################
// Copyright 2009-2013, Jon Gretarsson, Yue Yu.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// namespace TESSELLATION
//##################################################################### 
#ifndef __RANGE_TESSELLATION__
#define __RANGE_TESSELLATION__
 
namespace PhysBAM{
template<class TV> class RANGE;
template<class T> class TRIANGULATED_SURFACE;
template<class T> class SEGMENTED_CURVE_2D;
template<class T,int d> class VECTOR;

namespace TESSELLATION{
//#####################################################################
    template<class T> TRIANGULATED_SURFACE<T>* Generate_Triangles(const RANGE<VECTOR<T,3> >& box);
    template<class T> TRIANGULATED_SURFACE<T>* Generate_Triangles(const RANGE<VECTOR<T,3> >& box,const T suggested_dx);
    template<class T> SEGMENTED_CURVE_2D<T>* Generate_Triangles(const RANGE<VECTOR<T,2> >& box);
//#####################################################################
}
}
#endif
