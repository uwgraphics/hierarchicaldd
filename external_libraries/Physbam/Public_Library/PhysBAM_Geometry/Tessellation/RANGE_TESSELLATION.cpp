//#####################################################################
// Copyright 2009-2013, Jon Gretarsson, Yue Yu.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// namespace TESSELLATION
//##################################################################### 
#include <PhysBAM_Tools/Math_Tools/RANGE.h>
#include <PhysBAM_Tools/Arrays_Computations/SORT.h>
#include <PhysBAM_Geometry/Tessellation/RANGE_TESSELLATION.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/SEGMENTED_CURVE_2D.h>

namespace PhysBAM{
namespace TESSELLATION{
//#####################################################################
// Function Generate_Triangles
//#####################################################################
template<class T> TRIANGULATED_SURFACE<T>* Generate_Triangles(const RANGE<VECTOR<T,3> >& box)
{
    typedef VECTOR<T,3> TV;
    TRIANGULATED_SURFACE<T>* surface=TRIANGULATED_SURFACE<T>::Create();
    GEOMETRY_PARTICLES<TV>& particles=surface->particles;particles.array_collection->Add_Elements(8);
    particles.X(1)=TV(box.min_corner.x,box.min_corner.y,box.max_corner.z);particles.X(2)=TV(box.max_corner.x,box.min_corner.y,box.max_corner.z);
    particles.X(3)=TV(box.max_corner.x,box.max_corner.y,box.max_corner.z);particles.X(4)=TV(box.min_corner.x,box.max_corner.y,box.max_corner.z);
    particles.X(5)=TV(box.min_corner.x,box.min_corner.y,box.min_corner.z);particles.X(6)=TV(box.max_corner.x,box.min_corner.y,box.min_corner.z);
    particles.X(7)=TV(box.max_corner.x,box.max_corner.y,box.min_corner.z);particles.X(8)=TV(box.min_corner.x,box.max_corner.y,box.min_corner.z);
    TRIANGLE_MESH& mesh=surface->mesh;mesh.number_nodes=8;mesh.elements.Preallocate(12);
    mesh.elements.Append(VECTOR<int,3>(1,4,8));mesh.elements.Append(VECTOR<int,3>(8,5,1));
    mesh.elements.Append(VECTOR<int,3>(2,6,7));mesh.elements.Append(VECTOR<int,3>(7,3,2));
    mesh.elements.Append(VECTOR<int,3>(5,2,1));mesh.elements.Append(VECTOR<int,3>(2,5,6));
    mesh.elements.Append(VECTOR<int,3>(7,8,3));mesh.elements.Append(VECTOR<int,3>(4,3,8));
    mesh.elements.Append(VECTOR<int,3>(5,8,7));mesh.elements.Append(VECTOR<int,3>(7,6,5));
    mesh.elements.Append(VECTOR<int,3>(1,2,3));mesh.elements.Append(VECTOR<int,3>(1,3,4));
    return surface;
}
template<class T> TRIANGULATED_SURFACE<T>* Generate_Triangles(const RANGE<VECTOR<T,3> >& box,const T suggested_dx)
{
    typedef VECTOR<T,3> TV;
    typedef VECTOR<int,3> TV_INT;
    TV box_suggested_counts=box.Edge_Lengths()/suggested_dx;
    for(int d=1;d<=3;++d)if(box_suggested_counts(d)-(T)1<=(T)-1e-5){std::stringstream ss;ss<<"suggested dx must be smaller than all of "<<box.Edge_Lengths();PHYSBAM_FATAL_ERROR(ss.str());}
    TV_INT box_counts;VECTOR<bool,3> even;even.Fill(false);
    for(int d=1;d<=3;++d){
        if(box_suggested_counts(d)-int(box_suggested_counts(d))<(T)1e-5)box_counts(d)=int(box_suggested_counts(d));
        else if(box_suggested_counts(d)-(int(box_suggested_counts(d))+1)>(T)-1e-5)box_counts(d)=int(box_suggested_counts(d))+1;
        if(box_counts(d))even(d)=(box_counts(d)%2==0);}
    int even_count=0;for(int d=1;d<=3;++d)if(even(d))++even_count;
    while(even_count<2){
        T min_distortion=2;int min_distortion_dimension=0;int closest_even_count=0;
        for(int d=1;d<=3;++d)if(!box_counts(d)){
            closest_even_count=int(box_suggested_counts(d))+int(box_suggested_counts(d))%2;
            T distortion=abs(closest_even_count-box_suggested_counts(d))/box_suggested_counts(d);
            if(distortion<min_distortion){min_distortion=distortion;min_distortion_dimension=d;}}
        if(min_distortion_dimension){box_counts(min_distortion_dimension)=closest_even_count;even(min_distortion_dimension)=true;++even_count;continue;}
        min_distortion=2;min_distortion_dimension=0;
        for(int d=1;d<=3;++d)if(!even(d)){
            T distortion=(T)1/(T)box_counts(d);
            if(distortion<min_distortion){min_distortion=distortion;min_distortion_dimension=d;}}
        ++box_counts(min_distortion_dimension);even(min_distortion_dimension)=true;++even_count;}
    for(int d=1;d<=3;++d)if(!box_counts(d))box_counts(d)=int((T).5+box_suggested_counts(d));
    int odd_dimension=0;
    ARRAY<PAIR<int,int> > dimension_counts;for(int d=1;d<=3;++d){if(even(d))dimension_counts.Append(Tuple(d,box_counts(d)));else odd_dimension=d;}
    Sort(dimension_counts,Field_Comparison(&PAIR<int,int>::y));
    if(even_count==3)odd_dimension=dimension_counts(3).x;
    int medium_dimension=dimension_counts(2).x;int short_dimension=dimension_counts(1).x;
    TV_INT counts(box_counts(odd_dimension),box_counts(medium_dimension),box_counts(short_dimension));// sorted box_counts
    TV_INT sorted_d(odd_dimension,medium_dimension,short_dimension);
    int   odd_medium_nodes=(counts(2)-1)*counts(1)-(counts(2)/2-1);
    int    odd_short_nodes=(counts(3)-1)*counts(1)-(counts(3)/2-1);
    int medium_short_nodes=(counts(3)-1)*counts(2)-(counts(3)/2-1);
    int total_number_of_nodes=4*(box_counts.Sum()-1)+2*(odd_medium_nodes+odd_short_nodes+medium_short_nodes);
    int total_number_of_elements=2*(counts(2)*(2*counts(1)+1)+counts(3)*(2*counts(1)+1)+counts(3)*(2*counts(2)+1));

    TRIANGULATED_SURFACE<T>* surface=TRIANGULATED_SURFACE<T>::Create();
    GEOMETRY_PARTICLES<TV>& particles=surface->particles;particles.array_collection->Add_Elements(total_number_of_nodes);
    TRIANGLE_MESH& mesh=surface->mesh;mesh.number_nodes=total_number_of_nodes;mesh.elements.Preallocate(total_number_of_elements);

    HASHTABLE<TV_INT,int> nodes_to_particle;
    TV unit_lengths=(T).5*box.Edge_Lengths()/(TV)box_counts;
    // add nodes on edges and corners
    for(int d=1;d<=3;++d)
        for(int iside=0;iside<=1;++iside)
            for(int jside=0;jside<=1;++jside)
                for(int index=(d==1?0:1);index<=box_counts(d)-(d==1?0:1);++index){
                    TV_INT node_index=TV_INT::Axis_Vector(d)*index+TV_INT::Axis_Vector(d%3+1)*iside*box_counts(d%3+1)+TV_INT::Axis_Vector((d+1)%3+1)*jside*box_counts((d+1)%3+1);
                    node_index*=2;
                    particles.X(nodes_to_particle.Size()+1)=box.min_corner+unit_lengths*(TV)node_index;
                    nodes_to_particle.Set(node_index,nodes_to_particle.Size()+1);}
    // add nodes on faces
    for(int sd=1;sd<=3;++sd){
        int xd=1,yd=2;if(sd==1){xd=2;yd=3;}else if(sd==2){xd=1;yd=3;}
        for(int side=0;side<=1;++side)
            for(int j=2;j<counts(yd)*2;j+=2)
                for(int i=2-(j/2)%2;i<counts(xd)*2;i+=2){
                    TV_INT node_index=TV_INT::Axis_Vector(sorted_d(sd))*side*counts(sd)*2+TV_INT::Axis_Vector(sorted_d(xd))*i+TV_INT::Axis_Vector(sorted_d(yd))*j;
                    particles.X(nodes_to_particle.Size()+1)=box.min_corner+unit_lengths*(TV)node_index;
                    nodes_to_particle.Insert(node_index,nodes_to_particle.Size()+1);}}
    // add triangles that face outward
    for(int sd=1;sd<=3;++sd){
        int xd=1,yd=2;if(sd==1){xd=2;yd=3;}else if(sd==2){xd=1;yd=3;}
        for(int side=0;side<=1;++side){
            int orientation=sorted_d(yd)-sorted_d(xd);if(orientation==-2)orientation=1;else if(orientation==2)orientation=-1;
            orientation*=2*side-1;
            TV_INT side_offset=TV_INT::Axis_Vector(sorted_d(sd))*side*counts(sd)*2;
            for(int j=0;j<counts(yd)*2;j+=2)
                for(int i=0;i<=counts(xd)*2;++i){
                    int ij_up=(1-i%2*2)*(1-(j/2)%2*2);
                    int i_offset=orientation*ij_up;
                    TV_INT node_index_2=side_offset+TV_INT::Axis_Vector(sorted_d(xd))*i+TV_INT::Axis_Vector(sorted_d(yd))*(j+1-ij_up);
                    TV_INT node_index_1=side_offset+TV_INT::Axis_Vector(sorted_d(xd))*min(counts(xd)*2,max(0,i-i_offset))+TV_INT::Axis_Vector(sorted_d(yd))*(j+1+ij_up);
                    TV_INT node_index_3=side_offset+TV_INT::Axis_Vector(sorted_d(xd))*min(counts(xd)*2,max(0,i+i_offset))+TV_INT::Axis_Vector(sorted_d(yd))*(j+1+ij_up);
                    mesh.elements.Append(TV_INT(nodes_to_particle.Get(node_index_1),nodes_to_particle.Get(node_index_2),nodes_to_particle.Get(node_index_3)));}}}
    return surface;
}
template<class T> SEGMENTED_CURVE_2D<T>* Generate_Triangles(const RANGE<VECTOR<T,2> >& box)
{
    typedef VECTOR<T,2> TV;
    SEGMENTED_CURVE_2D<T>* surface=new SEGMENTED_CURVE_2D<T>(*new SEGMENT_MESH,*new GEOMETRY_PARTICLES<TV>);
    //surface->need_destroy_mesh=surface->need_destroy_particles=true;
    GEOMETRY_PARTICLES<TV>& particles=surface->particles;particles.array_collection->Add_Elements(4);
    particles.X(1)=TV(box.min_corner.x,box.min_corner.y);particles.X(2)=TV(box.max_corner.x,box.min_corner.y);
    particles.X(3)=TV(box.max_corner.x,box.max_corner.y);particles.X(4)=TV(box.min_corner.x,box.max_corner.y);
    SEGMENT_MESH& mesh=surface->mesh;mesh.number_nodes=4;mesh.elements.Preallocate(4);
    mesh.elements.Append(VECTOR<int,2>(1,2));mesh.elements.Append(VECTOR<int,2>(2,3));
    mesh.elements.Append(VECTOR<int,2>(3,4));mesh.elements.Append(VECTOR<int,2>(4,1));
    return surface;
}
//#####################################################################
}
template TRIANGULATED_SURFACE<float>* TESSELLATION::Generate_Triangles(const RANGE<VECTOR<float,3> >&);
template TRIANGULATED_SURFACE<float>* TESSELLATION::Generate_Triangles(const RANGE<VECTOR<float,3> >&,const float);
template SEGMENTED_CURVE_2D<float>* TESSELLATION::Generate_Triangles(const RANGE<VECTOR<float,2> >&);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template TRIANGULATED_SURFACE<double>* TESSELLATION::Generate_Triangles(const RANGE<VECTOR<double,3> >&);
template TRIANGULATED_SURFACE<double>* TESSELLATION::Generate_Triangles(const RANGE<VECTOR<double,3> >&,const double);
template SEGMENTED_CURVE_2D<double>* TESSELLATION::Generate_Triangles(const RANGE<VECTOR<double,2> >&);
#endif
}
