//#####################################################################
// Copyright 2013, Linhai Qiu
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Geometry/Surface_Mesh_Level_Sets/FAST_MARCHING_METHOD_SURFACE_MESH.h>
#include <PhysBAM_Geometry/Surface_Mesh_Advection/ADVECTION_SEMI_LAGRANGIAN_SURFACE_MESH.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/POINT_SIMPLICES_1D.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/SEGMENTED_CURVE.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> FAST_MARCHING_METHOD_SURFACE_MESH<TV>::
FAST_MARCHING_METHOD_SURFACE_MESH(const T_SIMPLICIAL_OBJECT& simplicial_object_input,const ARRAY_VIEW<TV>& normal_input)
:simplicial_object(simplicial_object_input),normal(normal_input)
{}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> FAST_MARCHING_METHOD_SURFACE_MESH<TV>::
~FAST_MARCHING_METHOD_SURFACE_MESH()
{}
//#####################################################################
// Fast_Marching_Method
//#####################################################################
template<class TV> void FAST_MARCHING_METHOD_SURFACE_MESH<TV>::
Fast_Marching_Method(ARRAY_VIEW<T>& phi)
{
    int heap_length=0;
    ARRAY<int> close_k(phi.m);close_k.Fill(0);
    ARRAY<int> heap(phi.m,false);
    ARRAY<bool> done(phi.m);done.Fill(false);
    Initialize_Interface(phi,done,close_k,heap,heap_length);
    while(heap_length!=0){
        int index=heap(1);
        done(index)=true;close_k(index)=0;//add to done, remove from close
        FAST_MARCHING<T>::Down_Heap(phi,close_k,heap,heap_length);heap_length--;
        ARRAY<int>& neighbor_nodes=simplicial_object.mesh.neighbor_nodes->operator()(index);
        for(int n=1;n<=neighbor_nodes.m;n++){
            if(!done(neighbor_nodes(n)))
                Update_Or_Add_Neighbor(phi,done,close_k,heap,heap_length,neighbor_nodes(n));}}
}
//#####################################################################
// Update_Or_Add_Neighbor
//#####################################################################
template<class TV> void FAST_MARCHING_METHOD_SURFACE_MESH<TV>::
Update_Or_Add_Neighbor(ARRAY_VIEW<T>& phi,ARRAY<bool>& done,ARRAY<int>& close_k,ARRAY<int>& heap,int& heap_length,const int& neighbor)
{
    if(close_k(neighbor)){
        Update_Close_Point(phi,done,neighbor);
        FAST_MARCHING<T>::Up_Heap(phi,close_k,heap,close_k(neighbor));}
    else{
        heap_length++;heap(heap_length)=neighbor;close_k(neighbor)=heap_length;
        Update_Close_Point(phi,done,neighbor);
        FAST_MARCHING<T>::Up_Heap(phi,close_k,heap,heap_length);}
}
//#####################################################################
// Initialize_Interface
//#####################################################################
template<class T> void Initialize_Interface_Helper(T& closest_distance,const ARRAY<VECTOR<T,1> >& value){PHYSBAM_NOT_IMPLEMENTED();}
template<class T> void Initialize_Interface_Helper(T& closest_distance,const ARRAY<VECTOR<T,2> >& value)
{
    if(value.m==1) closest_distance=min(closest_distance,value(1).Magnitude());
}
template<class T> void Initialize_Interface_Helper(T& closest_distance,const ARRAY<VECTOR<T,3> >& value)
{
    typedef VECTOR<T,3> TV;
    if(value.m==1) closest_distance=min(closest_distance,value(1).Magnitude());
    else if(value.m==2){
        TV interface_edge=value(2)-value(1);
        TV parallelogram_area=TV::Cross_Product(value(1),value(2));
        T normal_distance=parallelogram_area.Magnitude()/interface_edge.Magnitude();
        TV normal_path=TV::Cross_Product(interface_edge,parallelogram_area).Normalized()*normal_distance;
        T t=TV::Dot_Product(normal_path-value(1),interface_edge)/interface_edge.Magnitude_Squared();
        if(t<0) closest_distance=min(closest_distance,value(1).Magnitude());
        else if(t>1) closest_distance=min(closest_distance,value(2).Magnitude());
        else{closest_distance=min(closest_distance,normal_distance);}
    }
}
template<class TV> void FAST_MARCHING_METHOD_SURFACE_MESH<TV>::
Initialize_Interface(ARRAY_VIEW<T>& phi,ARRAY<bool>& done,ARRAY<int>& close_k,ARRAY<int>& heap,int& heap_length)
{
    const T tolerance=(T)1e-7;
    const ARRAY<TV_INT>& elements=simplicial_object.mesh.elements;
    for(int element_index=1;element_index<=elements.m;element_index++){
        TV_INT element=elements(element_index);
        bool some_in=false,some_out=false;
        for(int n=1;n<=TV_INT::m;n++){
            if(phi(element(n))>-tolerance) some_out=true;
            if(phi(element(n))<tolerance) some_in=true;}
        if(some_in&&some_out){
            for(int n=1;n<=TV_INT::m;n++)Add_To_Initial(done,close_k,element(n));}}

    ARRAY<T> phi_new(phi.m);
    for(int index=1;index<=phi.m;index++)if(done(index)){
        ARRAY<int>& neighbor_elements=simplicial_object.mesh.incident_elements->operator()(index);
        ARRAY_VIEW<TV>& X=simplicial_object.particles.X;
        T current_phi=phi(index);
        if(abs(current_phi)<tolerance){phi_new(index)=current_phi;continue;}
        T closest_distance(FLT_MAX);
        bool must_find_interface_neighbor=false;
        for(int i=1;i<=neighbor_elements.m;i++){
            int element_index=neighbor_elements(i);
            TV_INT element=simplicial_object.mesh.elements(element_index);
            int n=1;for(;n<=TV::m;n++)if(element(n)==index)break;
            ARRAY<TV> value;
            for(int j=0;j<=TV::m-2;++j){int neighbor_index=element((n+j)%TV::m+1);
                T neighbor_phi=phi(neighbor_index);
                if(done(neighbor_index) && ((current_phi>0 && neighbor_phi<tolerance) || (current_phi<0 && neighbor_phi>-tolerance))){
                    value.Append(LEVELSET_UTILITIES<T>::Theta(current_phi,neighbor_phi)*(X(neighbor_index)-X(index)));}}
            Initialize_Interface_Helper(closest_distance,value);
            if(value.m)must_find_interface_neighbor=true;
        }
        if(!must_find_interface_neighbor)PHYSBAM_FATAL_ERROR("must_find_interface_neighbor");
        phi_new(index)=closest_distance*LEVELSET_UTILITIES<T>::Sign(current_phi);
    }
    for(int index=1;index<=phi.m;index++)if(done(index)){phi(index)=phi_new(index);}

    // initialize close points
    for(int index=1;index<=phi.m;index++) if(close_k(index)){
        Update_Close_Point(phi,done,index);
        heap_length++;heap(heap_length)=index;
        FAST_MARCHING<T>::Up_Heap(phi,close_k,heap,heap_length);}
}
//#####################################################################
// Update_Close_Point
//#####################################################################
template<class T> void Update_Close_Point_Helper(const VECTOR<int,1>& element,const ARRAY_VIEW<VECTOR<T,1> >& X,const ARRAY_VIEW<VECTOR<T,1> >& normal,const POINT_SIMPLEX_MESH& mesh,ARRAY_VIEW<T>& phi,const ARRAY<bool>& done,const int& index){PHYSBAM_NOT_IMPLEMENTED();}
template<class T> void Update_Close_Point_Helper(const VECTOR<int,2>& element,const ARRAY_VIEW<VECTOR<T,2> >& X,const ARRAY_VIEW<VECTOR<T,2> >& normal,const SEGMENT_MESH& mesh,ARRAY_VIEW<T>& phi,const ARRAY<bool>& done,const int& index){
    int n=1;for(;n<=2;n++)if(element(n)==index)break;
    int a=element(3-n);
    if(!done(a))return;
    T length_a=(X(a)-X(element(n))).Magnitude();
    T phi_from_a=phi(a)>0?phi(a)+length_a:phi(a)-length_a;
    phi(index)=minmag(phi(index),phi_from_a);
}
template<class T> void Update_Close_Point_Helper_Known_Indices(const ARRAY_VIEW<VECTOR<T,3> >& X,ARRAY_VIEW<T>& phi,const ARRAY<bool>& done,const int& index,int a,int b)
{
    typedef VECTOR<T,3> TV;
    int c=index;
    if(!done(a)&&!done(b))return;
    if(abs(phi(a))>abs(phi(b))){int temp=a;a=b;b=temp;}

    TV pos_a=X(a),pos_b=X(b),pos_c=X(c);
    T length_a=(pos_b-pos_c).Magnitude();
    T length_b=(pos_a-pos_c).Magnitude();
    T phi_from_a=phi(a)>0?phi(a)+length_b:phi(a)-length_b,
      phi_from_b=phi(b)>0?phi(b)+length_a:phi(b)-length_a;
    if(!done(a)){phi(index)=minmag(phi(index),phi_from_b);return;}
    if(!done(b)){phi(index)=minmag(phi(index),phi_from_a);return;}
    assert(LEVELSET_UTILITIES<T>::Sign(phi(a))==LEVELSET_UTILITIES<T>::Sign(phi(b)));
    T length_c=(pos_a-pos_b).Magnitude();
    T cos_theta=(length_a*length_a+length_b*length_b-length_c*length_c)/((T)2*length_a*length_b);
    T u=abs(phi(b))-abs(phi(a));
    T q_a=length_c*length_c;
    T q_b=(T)2*length_b*u*(length_a*cos_theta-length_b);
    T q_c=length_b*length_b*(u*u-length_a*length_a*(1-cos_theta*cos_theta));
    T b24ac=q_b*q_b-(T)4*q_a*q_c;
    if(b24ac<0){
        if(b24ac>(T)-1e-07)b24ac=(T)0;
        else{phi(index)=minmag(phi(index),phi_from_a,phi_from_b);return;}}
    T t=(-q_b+sqrt(b24ac))/((T)2*q_a);
    T condition=length_b*(t-u)/(t*length_a);
    if(u<t&&condition>cos_theta&&(cos_theta<1e-5||condition<(T)1/cos_theta))
        phi(index)=minmag(phi(index),phi(a)>0?phi(a)+t:phi(a)-t);
    else
        phi(index)=minmag(phi(index),phi_from_a,phi_from_b);
}
template<class T> void Update_Close_Point_Helper(const VECTOR<int,3>& element,const ARRAY_VIEW<VECTOR<T,3> >& X,const ARRAY_VIEW<VECTOR<T,3> >& normal,const TRIANGLE_MESH& mesh,ARRAY_VIEW<T>& phi,const ARRAY<bool>& done,const int& index)
{
    typedef VECTOR<T,3> TV;
    int n=1;for(;n<=3;++n)if(element(n)==index)break;
    int a=element(n%3+1),b=element((n+1)%3+1),c=index;
    TV edge_ca=X(a)-X(c),edge_cb=X(b)-X(c);
    if(TV::Dot_Product(edge_ca,edge_cb)>=(T)0) Update_Close_Point_Helper_Known_Indices(X,phi,done,index,a,b);
    else{ // point c is obtuse
        T length_a=(X(b)-X(c)).Magnitude(),length_b=(X(a)-X(c)).Magnitude(),length_c=(X(a)-X(b)).Magnitude();
        T cos_theta_a=(length_b*length_b+length_c*length_c-length_a*length_a)/((T)2*length_b*length_c);
        T cos_theta_b=(length_a*length_a+length_c*length_c-length_b*length_b)/((T)2*length_a*length_c);
        TV edge_ab=X(b)-X(a);
        ARRAY<TV> pos(2),ray(2);
        pos(1)=X(a)+edge_ab.Normalized()*(edge_ca.Magnitude()/cos_theta_a);
        pos(2)=X(b)-edge_ab.Normalized()*(edge_cb.Magnitude()/cos_theta_b);
        ray(1)=(pos(1)-X(c)).Normalized();ray(2)=(pos(2)-X(c)).Normalized();
        const T tolerance=(T)1e-7;
        int virtual_node_index=0;
        bool debug=false;
        if(debug)LOG::cout<<"starting obtuse edge: "<<a<<", "<<b<<std::endl;
        for(VECTOR<int,2> current_edge(a,b);!virtual_node_index;){
            if(debug)LOG::cout<<"\tfollowing obtuse edge: "<<current_edge<<std::endl;
            TV current_edge_ray=(X(current_edge(2))-X(current_edge(1))).Normalized();
            int opposite_node_index=0;
            ARRAY<int>& incident_elements=mesh.incident_elements->operator()(current_edge(2));
            for(int incident_index=1;incident_index<=incident_elements.Size();incident_index++){
                const VECTOR<int,3>& incident_element=mesh.elements(incident_elements(incident_index));
                int i=1;for(;i<=3;i++)if(incident_element(i)==current_edge(2))break;
                if(incident_element(i%3+1)==current_edge(1)){opposite_node_index=incident_element((i+1)%3+1);break;}}
            TV new_edge_12=X(opposite_node_index)-X(current_edge(1));
            TV new_triangle_normal=TV::Cross_Product(new_edge_12,current_edge_ray).Normalized();
            ARRAY<T> orientations(2);ARRAY<TV> folded_ray(2);
            ARRAY<ARRAY<TV> > edge_to_nodes(3);for(int t=1;t<=3;++t){edge_to_nodes(t).Resize(2);}
            for(int side=1;side<=2;++side){
                TV tangent_part=TV::Dot_Product(current_edge_ray,ray(side))*current_edge_ray;
                TV folded_normal=(new_edge_12-TV::Dot_Product(current_edge_ray,new_edge_12)*current_edge_ray).Normalized();
                folded_ray(side)=tangent_part+(ray(side)-tangent_part).Magnitude()*folded_normal;
                edge_to_nodes(1)(side)=X(current_edge(1))-pos(side);
                edge_to_nodes(2)(side)=X(opposite_node_index)-pos(side);
                edge_to_nodes(3)(side)=X(current_edge(2))-pos(side);
                orientations(side)=TV::Dot_Product(TV::Cross_Product(folded_ray(side),edge_to_nodes(2)(side).Normalized()),new_triangle_normal);}
            int check_interesection_index=0;VECTOR<int,2> next_edge;
            if(orientations(1)<tolerance){
                if(orientations(2)>-tolerance) virtual_node_index=opposite_node_index;
                else{check_interesection_index=2;
                    next_edge=VECTOR<int,2>(opposite_node_index,current_edge(2));}}
            else{
                if(orientations(2)>-tolerance){check_interesection_index=1;
                    next_edge=VECTOR<int,2>(current_edge(1),opposite_node_index);}
                else PHYSBAM_FATAL_ERROR("obtuse search paths crossed: on two edges");}
            if(check_interesection_index){VECTOR<T,2> t;
                for(int side=1;side<=2;++side){
                    T s=0;bool too_short=false;
                    ADVECTION_SEMI_LAGRANGIAN_SURFACE_MESH<TV,T>::Intersection_Segments(new_triangle_normal,TV(),folded_ray(side),
                        edge_to_nodes(check_interesection_index)(side),
                        edge_to_nodes(check_interesection_index+1)(side),s,t(side),too_short);
                    TV edge_to_new_pos=((T)1-t(side))*edge_to_nodes(check_interesection_index)(side)+t(side)*edge_to_nodes(check_interesection_index+1)(side);
                    pos(side)+=edge_to_new_pos;
                    ray(side)=edge_to_new_pos.Normalized();}
                if(t(1)<t(2))PHYSBAM_FATAL_ERROR("obtuse search paths crossed: on the same edge");}
            current_edge=next_edge;
        }
        if(debug)LOG::cout<<"\tfound node: "<<virtual_node_index<<" index:"<<index<<", a:"<<a<<", b:"<<b<<", virt:"<<virtual_node_index<<std::endl;
        if(debug)LOG::cout<<"\t(cont)"<<virtual_node_index<<" before:"<<phi(index)<<","<<phi(a)<<","<<phi(b)<<","<<phi(virtual_node_index)<<std::endl;
        Update_Close_Point_Helper_Known_Indices(X,phi,done,index,a,virtual_node_index);
        T result_from_a=phi(index);
        Update_Close_Point_Helper_Known_Indices(X,phi,done,index,b,virtual_node_index);
        phi(index)=minmag(result_from_a,phi(index));
        if(debug)LOG::cout<<"\t(cont)"<<virtual_node_index<<" after:"<<phi(index)<<","<<phi(a)<<","<<phi(b)<<","<<phi(virtual_node_index)<<std::endl;
        }
}
template<class TV> void FAST_MARCHING_METHOD_SURFACE_MESH<TV>::
Update_Close_Point(ARRAY_VIEW<T>& phi,ARRAY<bool>& done,const int& index)
{
    phi(index)=(T)FLT_MAX;
    ARRAY<int>& neighbor_elements=simplicial_object.mesh.incident_elements->operator()(index);
    ARRAY_VIEW<TV>& X=simplicial_object.particles.X;
    for(int i=1;i<=neighbor_elements.m;i++){
        int element_index=neighbor_elements(i);
        TV_INT element=simplicial_object.mesh.elements(element_index);
        Update_Close_Point_Helper(element,X,normal,simplicial_object.mesh,phi,done,index);
    }
}
//#####################################################################
// Add_To_Initial
//#####################################################################
template<class TV> void FAST_MARCHING_METHOD_SURFACE_MESH<TV>::
Add_To_Initial(ARRAY<bool>& done,ARRAY<int>& close_k,const int& index)
{
    done(index)=true;close_k(index)=0; // add to done, remove from close
    ARRAY<int>& neighbor_nodes=simplicial_object.mesh.neighbor_nodes->operator()(index);
    for(int n=1;n<=neighbor_nodes.m;n++){
        if(!done(neighbor_nodes(n))) close_k(neighbor_nodes(n))=1;}
}
//#####################################################################
template class FAST_MARCHING_METHOD_SURFACE_MESH<VECTOR<float,1> >;
template class FAST_MARCHING_METHOD_SURFACE_MESH<VECTOR<float,2> >;
template class FAST_MARCHING_METHOD_SURFACE_MESH<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class FAST_MARCHING_METHOD_SURFACE_MESH<VECTOR<double,1> >;
template class FAST_MARCHING_METHOD_SURFACE_MESH<VECTOR<double,2> >;
template class FAST_MARCHING_METHOD_SURFACE_MESH<VECTOR<double,3> >;
#endif
