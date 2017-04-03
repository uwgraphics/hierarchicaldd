//#####################################################################
// Copyright 2015, Mridul Aanjaneya.
// This file is part of PhysBAM whose distribution is governed by the license 
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Common_Geometry/Topology_Based_Geometry/VORONOI_DIAGRAM.h>
#include <PhysBAM_Tools/Math_Tools/constants.h>
#include <PhysBAM_Tools/Vectors/VECTOR_2D.h>
#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>
using namespace PhysBAM;
//#####################################################################
// Compute_Face_Barycenters
//#####################################################################
template<class T>
void Face_Barycenter_Helper(const ARRAY<VECTOR<T,2> >& face,VECTOR<T,2>& centroid,T& area)
{
    PHYSBAM_ASSERT(face.Size()==2);
    centroid=(face(1)+face(2))*(T).5;
    area=(face(1)-face(2)).Magnitude();
}
template<class T>
void Face_Barycenter_Helper(const ARRAY<VECTOR<T,3> >& face,VECTOR<T,3>& centroid,T& area)
{
    typedef VECTOR<T,3> TV;
    TV base_vertex=face(1);
    for(int k=1;k<=face.Size();++k){int next=(k==face.Size())?1:k+1;
        T current_area=TV::Cross_Product(face(k)-base_vertex,face(next)-base_vertex).Magnitude()*(T).5;
        centroid+=current_area*(base_vertex+face(k)+face(next))*one_third;
        area+=current_area;}
    centroid/=area;
}
template<class TV> void VORONOI_DIAGRAM<TV>::
Compute_Face_Barycenters()
{
    face_areas.Clean_Memory();
    face_barycenters.Clean_Memory();
    for(int i=1;i<=face_vertices.Size();++i){
        ARRAY<TV> barycenters;ARRAY<T> areas;
        for(int j=1;j<=face_vertices(i).Size();++j){T area=(T)0.;TV centroid;
            Face_Barycenter_Helper(face_vertices(i)(j),centroid,area);
            barycenters.Append(centroid);
            areas.Append(area);}
        face_barycenters.Append(barycenters);
        face_areas.Append(areas);}
}
//#####################################################################
template class VORONOI_DIAGRAM<VECTOR<float,2> >;
template class VORONOI_DIAGRAM<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class VORONOI_DIAGRAM<VECTOR<double,2> >;
template class VORONOI_DIAGRAM<VECTOR<double,3> >;
#endif
