//#####################################################################
// Copyright 2013, Yue Yu.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// namespace TESSELLATION
//##################################################################### 
#include <PhysBAM_Tools/Matrices/FRAME.h>
#include <PhysBAM_Tools/Matrices/MATRIX_3X3.h>
#include <PhysBAM_Geometry/Basic_Geometry/CONICAL_FRUSTUM.h>
#include <PhysBAM_Geometry/Tessellation/CONICAL_FRUSTUM_TESSELLATION.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/SEGMENTED_CURVE_2D.h>

namespace PhysBAM{
namespace TESSELLATION{
//#####################################################################
// Function Generate_Triangles
//#####################################################################
template<class T> TRIANGULATED_SURFACE<T>* Generate_Triangles(const CONICAL_FRUSTUM<VECTOR<T,3> >& cylinder,const int resolution_height,const int resolution_radius)
{
    typedef VECTOR<T,3> TV;
    TRIANGULATED_SURFACE<T>* surface=TRIANGULATED_SURFACE<T>::Create();
    GEOMETRY_PARTICLES<TV>& particles=surface->particles;
    surface->Initialize_Conical_Frustum_Mesh_And_Particles(resolution_height,resolution_radius,cylinder.height,cylinder.radius1,cylinder.radius2,true);
    FRAME<TV> frame(cylinder.plane1.x1,ROTATION<TV>::From_Rotated_Vector(TV(1,0,0),cylinder.plane2.x1-cylinder.plane1.x1));
    particles.X=frame*particles.X;
    return surface;
}
template<class T> SEGMENTED_CURVE_2D<T>* Generate_Triangles(const CONICAL_FRUSTUM<VECTOR<T,2> >& cylinder,const int,const int)
{
    typedef VECTOR<T,2> TV;
    SEGMENTED_CURVE_2D<T>* surface=SEGMENTED_CURVE_2D<T>::Create();
    GEOMETRY_PARTICLES<TV>& particles=surface->particles;
    surface->Initialize_Conical_Frustum_Mesh_And_Particles(cylinder.height,cylinder.radius1,cylinder.radius2,true);
    FRAME<TV> frame(cylinder.plane1.x1,ROTATION<TV>::From_Rotated_Vector(TV(1,0),cylinder.plane2.x1-cylinder.plane1.x1));
    particles.X=frame*particles.X;
    return surface;
}
//#####################################################################
template TRIANGULATED_SURFACE<float>* Generate_Triangles(const CONICAL_FRUSTUM<VECTOR<float,3> >&,const int,const int);
template SEGMENTED_CURVE_2D<float>* Generate_Triangles(const CONICAL_FRUSTUM<VECTOR<float,2> >&,const int,const int);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template TRIANGULATED_SURFACE<double>* Generate_Triangles(const CONICAL_FRUSTUM<VECTOR<double,3> >&,const int,const int);
template SEGMENTED_CURVE_2D<double>* Generate_Triangles(const CONICAL_FRUSTUM<VECTOR<double,2> >&,const int,const int);
#endif
}
}
