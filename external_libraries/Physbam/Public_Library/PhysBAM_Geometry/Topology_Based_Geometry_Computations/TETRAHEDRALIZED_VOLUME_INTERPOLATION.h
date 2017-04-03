//#####################################################################
// Copyright 2013, Rahul Sheth.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class TETRAHEDRALIZED_VOLUME_INTERPOLATION
//#####################################################################
#ifndef __TETRAHEDRALIZED_VOLUME_INTERPOLATION__
#define __TETRAHEDRALIZED_VOLUME_INTERPOLATION__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>
#include <PhysBAM_Geometry/Topology/TRIANGLE_MESH.h>
namespace PhysBAM {
template<class T> class TETRAHEDRALIZED_VOLUME;
template<class T> class TRIANGULATED_SURFACE;

namespace TOPOLOGY_BASED_GEOMETRY_COMPUTATIONS {
//#####################################################################
// Function Compute_Mean_Value_Coordinates_For_Position
//#####################################################################
template<class T>
ARRAY<T> Compute_Mean_Value_Coordinates_For_Position(const VECTOR<T,3>& v,const int num_particles,
    const ARRAY_VIEW<VECTOR<T,3> >& tri_node_positions,const ARRAY<int>& tri_node_indices,TRIANGLE_MESH tri_mesh_in)
{
    const T epsilon=(T)1e-10;
    const T cage_face_epsilon=(T)1e-5;
    ARRAY<T> coords(num_particles);

    ARRAY<T> d(num_particles);
    ARRAY<VECTOR<T,3> > u(num_particles);
    for(int j=1;j<=tri_node_indices.Size();j++){
        const VECTOR<T,3>& pj=tri_node_positions(tri_node_indices(j));
        d(tri_node_indices(j))=(pj-v).Magnitude();
        if(d(tri_node_indices(j))<epsilon){coords(tri_node_indices(j))=(T)1.0;return coords;}
        u(tri_node_indices(j))=(pj-v)/d(tri_node_indices(j));}
    for(int i=1;i<=tri_mesh_in.elements.Size();i++){
        int x1,x2,x3;
        tri_mesh_in.elements(i).Get(x1,x2,x3);
        T sign_local=-sign_nonzero(VECTOR<T,3>::Dot_Product(VECTOR<T,3>::Cross_Product(u(x1),u(x2)),u(x3)));
        T theta1=(T)2.0*asin(min((u(x2)-u(x3)).Magnitude()/(T)2.0,(T)1.0));
        T theta2=(T)2.0*asin(min((u(x3)-u(x1)).Magnitude()/(T)2.0,(T)1.0));
        T theta3=(T)2.0*asin(min((u(x1)-u(x2)).Magnitude()/(T)2.0,(T)1.0));
        T h = (theta1+theta2+theta3)/(T)2.0;
        if(fabs(pi-h)<cage_face_epsilon) {
            ARRAYS_COMPUTATIONS::Fill(coords,(T)0);
            coords(x1)=sin(theta1)*d(x3)*d(x2);
            coords(x2)=sin(theta2)*d(x1)*d(x3);
            coords(x3)=sin(theta3)*d(x2)*d(x1);
            T bary_total=coords(x1)+coords(x2)+coords(x3);
            coords(x1)/=bary_total;coords(x2)/=bary_total;coords(x3)/=bary_total;
            return coords;}
        T c1=((T)2.0*sin(h)*sin(h-theta1))/(sin(theta2)*sin(theta3))-(T)1.0;
        T c2=((T)2.0*sin(h)*sin(h-theta2))/(sin(theta3)*sin(theta1))-(T)1.0;
        T c3=((T)2.0*sin(h)*sin(h-theta3))/(sin(theta1)*sin(theta2))-(T)1.0;
        c1=clamp(c1,(T)-1,(T)1);c2=clamp(c2,(T)-1,(T)1);c3=clamp(c3,(T)-1,(T)1);
        T s1=sign_local*sqrt((T)1.0-c1*c1);
        T s2=sign_local*sqrt((T)1.0-c2*c2);
        T s3=sign_local*sqrt((T)1.0-c3*c3);
        if(fabs(s1)<epsilon||fabs(s2)<epsilon||fabs(s3)<epsilon) continue;
        T w1=(theta1-c2*theta3-c3*theta2)/(d(x1)*sin(theta2)*s3);
        T w2=(theta2-c3*theta1-c1*theta3)/(d(x2)*sin(theta3)*s1);
        T w3=(theta3-c1*theta2-c2*theta1)/(d(x3)*sin(theta1)*s2);
        coords(x1)+=w1;coords(x2)+=w2;coords(x3)+=w3;}
    T total_w=(T)0;
    for(int i=1;i<=num_particles;i++){total_w+=coords(i);}
    for(int i=1;i<=num_particles;i++){coords(i)/=total_w;}
    return coords;
}
//#####################################################################
// Function Compute_Mean_Value_Coordinates
//#####################################################################
template<class T>
void Compute_Mean_Value_Coordinates(TETRAHEDRALIZED_VOLUME<T>& tv,ARRAY<ARRAY<T> >& mv_coords)
{
    ARRAY<int> tv_nodes(tv.mesh.number_nodes);
    if(tv.triangulated_surface==NULL){tv.Initialize_Triangulated_Surface();tv.mesh.Initialize_Boundary_Nodes();}
    mv_coords.Exact_Resize(tv_nodes.Size());
    ARRAY<int>& boundary_nodes=*tv.mesh.boundary_nodes;
    for(int i=1;i<=tv_nodes.Size();i++){
        mv_coords(i)=Compute_Mean_Value_Coordinates_For_Position(tv.particles.X(i),tv.particles.array_collection->Size(),
            tv.particles.X,boundary_nodes,*tv.mesh.boundary_mesh);}
}
} }
#endif
