//####################################################################
// Copyright 2002-2005, Ronald Fedkiw, Eran Guendelman, Geoffrey Irving, Frank Losasso, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class EXTRAPOLATION_SURFACE_MESH  
//##################################################################### 
#include <PhysBAM_Geometry/Surface_Mesh_Level_Sets/EXTRAPOLATION_SURFACE_MESH.h>
#include <PhysBAM_Geometry/Level_Sets/FAST_MARCHING.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV,class T2> EXTRAPOLATION_SURFACE_MESH<TV,T2>::
EXTRAPOLATION_SURFACE_MESH(const SURFACE_MESH& mesh_input,const ARRAY_VIEW<T>& phi_input,const ARRAY_VIEW<TV>& X_input,const ARRAY_VIEW<TV>& normal_input,ARRAY_VIEW<T2>& u_input)
    :mesh(mesh_input),phi(phi_input),X(X_input),normal(normal_input),u(u_input),max_dX((T)0),local_average_dX(X_input.m)
{
    for(int index=1;index<=X.m;index++){
        const ARRAY<int>& neighbor_nodes=mesh.neighbor_nodes->operator()(index);
            for(int n=1;n<=neighbor_nodes.m;n++){
                if(index>neighbor_nodes(n))continue; //no duplicated checks
                T edge=(X(index)-X(neighbor_nodes(n))).Magnitude();
                if(edge>max_dX)max_dX=edge;}}
    band_width=(T)1*max_dX;
    for(int index=1;index<=X.m;index++){
        const ARRAY<int>& neighbor_nodes=mesh.neighbor_nodes->operator()(index);
        T sum(0);for(int n=1;n<=neighbor_nodes.m;n++)sum+=(X(index)-X(neighbor_nodes(n))).Magnitude();
        local_average_dX(index)=sum/T(neighbor_nodes.m);}
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV,class T2> EXTRAPOLATION_SURFACE_MESH<TV,T2>::
~EXTRAPOLATION_SURFACE_MESH()
{} 
template<class T> void Add_To_Heap(const ARRAY_VIEW<T>& phi,ARRAY<int>& heap,int& heap_length,ARRAY<int>& close_k,const int index)
{heap_length++;heap(heap_length)=index;close_k(index)=heap_length;FAST_MARCHING<T>::Up_Heap(phi,close_k,heap,heap_length);}
//#####################################################################
// Function Extrapolate
//#####################################################################
template<class TV,class T2> void EXTRAPOLATION_SURFACE_MESH<TV,T2>::
Extrapolate(ARRAY<ARRAY<PAIR<int,T> > >* weights)
{
    int heap_length=0;
    ARRAY<int> close_k(phi.m);close_k.Fill(0);
    ARRAY<int> heap(phi.m,false);
    ARRAY<bool> done(phi.m);done.Fill(false);
    Initialize(done,close_k,heap,heap_length);
    if(weights){weights->Remove_All();weights->Resize(phi.m);}
    while(heap_length && phi(heap(1)) <= band_width){
        int index=heap(1);
        done(index)=true;close_k(index)=0;//add to done, remove from close
        FAST_MARCHING<T>::Down_Heap(phi,close_k,heap,heap_length);heap_length--;
        ARRAY<int> closest_node_indices;
        Update_Close_Point(done,index,closest_node_indices);
        if(weights){
            if(closest_node_indices.m!=1){std::stringstream ss;ss<<"using weights, index="<<index<<" closest_node_indices="<<closest_node_indices;PHYSBAM_FATAL_ERROR(ss.str());}
            ARRAY<ARRAY<PAIR<int,T> > >& w=*weights;
            w(index).Append(Tuple(closest_node_indices(1),(T)1.0)); //TODO handle >1 nodes after mpi port

            const int w_index_size_pre_insert=w(index).m;
            ARRAY<bool> node_to_delete(w_index_size_pre_insert);node_to_delete.Fill(false);
            for(int j=1;j<=w_index_size_pre_insert;++j){const int& jndex=w(index)(j).x;
                if(jndex==index)PHYSBAM_FATAL_ERROR("jndex==index, inf loop expr");
                if(w(jndex).m){
                    node_to_delete(j)=true;
                    for(int k=1;k<=w(jndex).m;++k){const int& kndex=w(jndex)(k).x;
                        w(index).Append(Tuple(kndex,w(index)(j).y*w(jndex)(k).y));}}}
            // delete the marked compound jndex from w(index) and merge duplicated kndex's in w(index)
            for(int j=w_index_size_pre_insert;j>=1;--j)if(node_to_delete(j))w(index).Remove_Index(j);
            for(int k=1;k<=w(index).m;++k){const int& kndex=w(index)(k).x;
                for(int h=w(index).m;h>k;--h){const int& hndex=w(index)(h).x;
                    if(hndex==kndex){
                        w(index)(k).y+=w(index)(h).y;
                        w(index).Remove_Index(h);}}}}
        const ARRAY<int>& neighbor_nodes=mesh.neighbor_nodes->operator()(index);
        for(int n=1;n<=neighbor_nodes.m;n++)
            if(!done(neighbor_nodes(n))&&!close_k(neighbor_nodes(n)))
                Add_To_Heap(phi,heap,heap_length,close_k,neighbor_nodes(n));}

    //TODO boundary->Apply_Boundary_Condition(node_grid,u,time);
}
//#####################################################################
// Function Initialize
//#####################################################################
// pass heap_length by reference
template<class TV,class T2> void EXTRAPOLATION_SURFACE_MESH<TV,T2>::
Initialize(ARRAY<bool>& done,ARRAY<int>& close_k,ARRAY<int>& heap,int& heap_length)
{
    for(int index=1;index<=phi.m;index++){if(phi(index)<=0)done(index)=true;}
    for(int index=1;index<=phi.m;index++) if(done(index)){
        const ARRAY<int>& neighbor_nodes=mesh.neighbor_nodes->operator()(index);
        for(int n=1;n<=neighbor_nodes.m;n++)
            if(!done(neighbor_nodes(n))&&!close_k(neighbor_nodes(n)))
                Add_To_Heap(phi,heap,heap_length,close_k,neighbor_nodes(n));}
}
//#####################################################################
// Function Update_Close_Point
//##################################################################### 
// needs done=0 around the outside of the domain
template<class T,int d> void 
Rotate_Helper(T& z,const ARRAY_VIEW<VECTOR<T,d> >& X,const ARRAY_VIEW<VECTOR<T,d> >& normal,int index,int closest_node_index)
{
}
template<class T> void 
Rotate_Helper(VECTOR<T,2>& z,const ARRAY_VIEW<VECTOR<T,2> >& X,const ARRAY_VIEW<VECTOR<T,2> >& normal,int index,int closest_node_index)
{
    typedef VECTOR<T,2> TV;
    const TV destination_Z=normal(index);
    const TV destination_X(-destination_Z.y,destination_Z.x);
    const TV source_Z=normal(closest_node_index);
    const TV source_X(-source_Z.y,source_Z.x);
    MATRIX<T,2> M_source,M_destination;
    M_source.Column(1)=source_Z;M_source.Column(2)=source_X;
    M_destination.Column(1)=destination_Z;M_destination.Column(2)=destination_X;
    ROTATION<TV> R_source=ROTATION<TV>(M_source);
    ROTATION<TV> R_destination=ROTATION<TV>(M_destination);
    z=(R_destination*R_source.Inverse()).Rotate(z);
}
template<class T> void 
Rotate_Helper(VECTOR<T,3>& z,const ARRAY_VIEW<VECTOR<T,3> >& X,const ARRAY_VIEW<VECTOR<T,3> >& normal,int index,int closest_node_index)
{
    typedef VECTOR<T,3> TV;

    const TV destination_Z=normal(index);
    TV destination_X=X(index)-X(closest_node_index);
    const TV destination_Y=TV::Cross_Product(destination_Z,destination_X).Normalized();
    destination_X=TV::Cross_Product(destination_Y,destination_Z);

    const TV source_Z=normal(closest_node_index);
    TV source_X=X(index)-X(closest_node_index);
    const TV source_Y=TV::Cross_Product(source_Z,source_X).Normalized();
    source_X=TV::Cross_Product(source_Y,source_Z);

    MATRIX<T,3> M_source,M_destination;
    M_source.Column(1)=source_X;M_source.Column(2)=source_Y;M_source.Column(3)=source_Z;
    M_destination.Column(1)=destination_X;M_destination.Column(2)=destination_Y;M_destination.Column(3)=destination_Z;
    ROTATION<TV> R_source=ROTATION<TV>(M_source);
    ROTATION<TV> R_destination=ROTATION<TV>(M_destination);
    z=(R_destination*R_source.Inverse()).Rotate(z);
}

template<class TV,class T2> void EXTRAPOLATION_SURFACE_MESH<TV,T2>::
Update_Close_Point(ARRAY<bool>& done,int index,ARRAY<int>& closest_node_indices)
{
    typedef typename TV::SCALAR T;
    T tolerance=(T).01*local_average_dX(index);
    T smallest_neighbor_phi(FLT_MAX);closest_node_indices.Remove_All();
    const ARRAY<int>& neighbor_nodes=mesh.neighbor_nodes->operator()(index);
    for(int n=1;n<=neighbor_nodes.m;n++)if(done(neighbor_nodes(n))){
        if(phi(neighbor_nodes(n))-smallest_neighbor_phi<-tolerance){
            smallest_neighbor_phi=phi(neighbor_nodes(n));
            closest_node_indices.Remove_All();
            closest_node_indices.Append(neighbor_nodes(n));}
        else if(phi(neighbor_nodes(n))-smallest_neighbor_phi<tolerance){
            closest_node_indices.Append(neighbor_nodes(n));}}
    T2 sum((T)0);
    for(int i=1;i<=closest_node_indices.m;++i){
        T2 neighbor_z=u(closest_node_indices(i));
        Rotate_Helper(neighbor_z,X,normal,index,closest_node_indices(i));
        sum+=neighbor_z;}
    u(index)=sum/T(closest_node_indices.m);
    //LOG::cout<<"extrapolate: u("<<index<<")=u("<<closest_node_indices<<")"<<std::endl;
}
//#####################################################################
template class EXTRAPOLATION_SURFACE_MESH<VECTOR<float,2>,float>;
template class EXTRAPOLATION_SURFACE_MESH<VECTOR<float,3>,float>;
template class EXTRAPOLATION_SURFACE_MESH<VECTOR<float,2>,VECTOR<float,2> >;
template class EXTRAPOLATION_SURFACE_MESH<VECTOR<float,3>,VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class EXTRAPOLATION_SURFACE_MESH<VECTOR<double,2>,double>;
template class EXTRAPOLATION_SURFACE_MESH<VECTOR<double,3>,double>;
template class EXTRAPOLATION_SURFACE_MESH<VECTOR<double,2>,VECTOR<double,2> >;
template class EXTRAPOLATION_SURFACE_MESH<VECTOR<double,3>,VECTOR<double,3> >;
#endif
