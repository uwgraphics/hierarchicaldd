//#####################################################################
// Copyright 2013, Wenlong Lu, Yue Yu.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Geometry/Surface_Mesh_Advection/ADVECTION_SEMI_LAGRANGIAN_SURFACE_MESH.h>
#include <PhysBAM_Tools/Matrices/MATRIX_2X2.h>
#include <PhysBAM_Tools/Matrices/MATRIX_3X3.h>
#include <PhysBAM_Geometry/Basic_Geometry/SEGMENT_2D.h>
#include <PhysBAM_Geometry/Basic_Geometry/TRIANGLE_3D.h>
#include <PhysBAM_Tools/Log/LOG.h>
namespace PhysBAM{

template<class TV,class T2> bool ADVECTION_SEMI_LAGRANGIAN_SURFACE_MESH<TV,T2>::
Intersection_Segments(const TV& normal,const TV& p1,const TV& p2,const TV& p3,const TV& p4,T& s,T& t,bool& too_short)
{
    if(TV::m!=3) PHYSBAM_NOT_IMPLEMENTED();

    T max_normal=(T)0;int intersection_plane=0;
    for(int n=1;n<=3;n++)
        if(abs(normal(n))>max_normal){max_normal=abs(normal(n));intersection_plane=n;}
    return Intersection_Segments(intersection_plane,p1,p2,p3,p4,s,t,too_short);
}

template<class TV,class T2> bool ADVECTION_SEMI_LAGRANGIAN_SURFACE_MESH<TV,T2>::
Intersection_Segments(int plane,const TV& p1, const TV& p2,const TV& p3, const TV& p4,T& s, T& t,bool& too_short)
{
    if(TV::m!=3) PHYSBAM_NOT_IMPLEMENTED();

    if(TV::Cross_Product((p2-p1).Normalized(),(p4-p3).Normalized()).Magnitude()<tolerance){
        s=t=std::numeric_limits<T>::quiet_NaN();return true;}

    T x1,x2,x3,x4,y1,y2,y3,y4; switch(plane){
        case 1:
            x1=p1(2);y1=p1(3);
            x2=p2(2);y2=p2(3);
            x3=p3(2);y3=p3(3);
            x4=p4(2);y4=p4(3);
            break;
        case 2:
            x1=p1(3);y1=p1(1);
            x2=p2(3);y2=p2(1);
            x3=p3(3);y3=p3(1);
            x4=p4(3);y4=p4(1);
            break;
        case 3:
        default:
            x1=p1(1);y1=p1(2);
            x2=p2(1);y2=p2(2);
            x3=p3(1);y3=p3(2);
            x4=p4(1);y4=p4(2);
            break;
    }
    T d_x12=x1*y2-y1*x2;
    T d_x34=x3*y4-y3*x4;
    T x1_x2=x1-x2;
    T x3_x4=x3-x4;
    T y1_y2=y1-y2;
    T y3_y4=y3-y4;
    T divide=x1_x2*y3_y4-y1_y2*x3_x4;

    T x=(d_x12*x3_x4-x1_x2*d_x34)/divide;
    T y=(d_x12*y3_y4-y1_y2*d_x34)/divide;
    
    if(abs(x2-x1)>abs(y2-y1)){
        if(abs(x2-x1)>tolerance) s=(x-x1)/(x2-x1);
        else{s=std::numeric_limits<T>::quiet_NaN();too_short=true;}}
    else{
        if(abs(y2-y1)>tolerance) s=(y-y1)/(y2-y1);
        else{s=std::numeric_limits<T>::quiet_NaN();too_short=true;}}
    if(abs(x4-x3)>abs(y4-y3)){
        if(abs(x4-x3)>tolerance) t=(x-x3)/(x4-x3);
        else{t=std::numeric_limits<T>::quiet_NaN();too_short=true;}}
    else{
        if(abs(y4-y3)>tolerance) t=(y-y3)/(y4-y3);
        else{t=std::numeric_limits<T>::quiet_NaN();too_short=true;}}
    return false;
}

template<class T,int d> void 
Rotate_Helper(T& z,const ROTATION<VECTOR<T,d> >& rotation)
{
}
template<class T,int d> void 
Rotate_Helper(VECTOR<T,d>& z,const ROTATION<VECTOR<T,d> >& rotation)
{
    z=rotation.Rotate(z);
}

template<class TV,class T2> template<int d> typename ENABLE_IF<d==2,void>::TYPE ADVECTION_SEMI_LAGRANGIAN_SURFACE_MESH<TV,T2>::
Update_Advection_Equation_Node(const SURFACE_MESH_TYPE& mesh,
            ARRAY_VIEW<T2>& Z,const ARRAY_VIEW<T2>& Z_ghost,
            const ARRAY_VIEW<VECTOR<T,d> >& X,const ARRAY_VIEW<VECTOR<T,d> >& V,const ARRAY_VIEW<VECTOR<T,d> >& normal,
            const RIGID_GEOMETRY<TV>& rigid_geometry,const FRAME<TV>& last_step_frame,
            const T dt,const T time,ARRAY<ARRAY<TRIPLE<int,T,T2> > >* weights_to_cell,
            bool forward)
{
    ROTATION<TV> current_rigid_rotation=rigid_geometry.Frame().Rotation();
    ROTATION<TV> last_rigid_rotation=last_step_frame.Rotation();
    for(int index=1;index<=Z.Size();index++){
        bool debug=false;
        if(debug)LOG::cout<<"index="<<index<<std::endl;
        TV mesh_velocity=rigid_geometry.Object_Space_Vector(rigid_geometry.Pointwise_Object_Velocity(rigid_geometry.Frame()*X(index)));
        if(forward)mesh_velocity*=(T)-1;
        TV current_movement=(rigid_geometry.Object_Space_Vector(V(index))-mesh_velocity)*(-dt);
        current_movement=current_movement-TV::Dot_Product(normal(index),current_movement)*normal(index);
        VECTOR<int,2> current_edge(index,index); 
        T current_t=1;
        const TV destination_Z=normal(index);
        const TV destination_X(-destination_Z.y,destination_Z.x);
        if(debug)LOG::cout<<"Initial Condition: "<<X(index)<<" "<<current_movement<<" "<<current_movement.Magnitude()<<" "<<current_edge<<std::endl;
        while(current_movement.Magnitude()>tolerance){
            TV projected_movement=current_movement-TV::Dot_Product(normal(current_edge(2)),current_movement)*normal(current_edge(2));
            if(projected_movement.Magnitude()<tolerance){
                if(debug)LOG::cout<<"projected_movement vanish="<<projected_movement<<std::endl;
                current_t=1;
                break;}
            if(debug)LOG::cout<<"Stepping phase starts"<<std::endl;
            ARRAY<int> incident_elements=mesh.incident_elements->operator()(current_edge(2));
            bool must_find_an_edge=false;
            for(int incident_index=1;incident_index<=incident_elements.Size();incident_index++){
                VECTOR<int,2> incident_element=mesh.elements(incident_elements(incident_index));
                int n=1;for(;n<=2;++n)if(incident_element(n)==current_edge(2))break;
                VECTOR<T,1> last_orientation=TV::Cross_Product(normal(current_edge(2)),current_movement);
                VECTOR<T,1> current_orientation=TV::Cross_Product(normal(current_edge(2)),X(incident_element(3-n))-X(current_edge(2)));
                if(debug){LOG::cout<<"normal(current_edge(2))="<<normal(current_edge(2))<<",current_movement="<<current_movement<<",edge="<<X(incident_element(3-n))-X(current_edge(2))<<std::endl;
                    LOG::cout<<"last_orientation="<<last_orientation<<" current_orientation="<<current_orientation<<std::endl;}
                if(last_orientation(1)*current_orientation(1)>0){
                    current_edge(1)=current_edge(2);current_edge(2)=incident_element(3-n);
                    must_find_an_edge=true;
                    break;
                }
            }
            if(!must_find_an_edge)PHYSBAM_FATAL_ERROR();
            TV current_edge_vector=X(current_edge(2))-X(current_edge(1));
            current_t=current_movement.Magnitude()/current_edge_vector.Magnitude();
            if(current_t>1){
                TV tangent(-normal(current_edge(2)).y,normal(current_edge(2)).x);
                if(TV::Dot_Product(tangent,current_edge_vector)<0)tangent=-tangent;
                current_movement=tangent*(current_movement.Magnitude()-current_edge_vector.Magnitude()); // TODO debug this
                if(debug)LOG::cout<<"Move Rotate: "<<current_edge(1)<<" "<<current_movement<<" "<<current_movement.Magnitude()<<" "<<current_edge<<std::endl;
            }
            else{
                current_movement=TV();
                if(debug)LOG::cout<<"Move Done: "<<X(current_edge(1))<<" "<<current_movement<<" "<<current_movement.Magnitude()<<" "<<current_edge<<std::endl;
                break;
            }
            if(debug)LOG::cout<<"Stepping phase ends"<<std::endl;
        }
        // barycentric weight, interpolation, done
        if(shadow_X) (*shadow_X)(index)=(1-current_t)*X(current_edge(1))+current_t*X(current_edge(2));
        if(debug)LOG::cout<<"index="<<index<<" t="<<current_t<<" "<<Z_ghost(current_edge(1))<<" "<<Z_ghost(current_edge(2))<<std::endl;
        assert(current_t>-tolerance && current_t-(T)1<tolerance);
        T2 Z_ghost_1=Z_ghost(current_edge(1)),Z_ghost_2=Z_ghost(current_edge(2));
        const TV source_Z=((1-current_t)*normal(current_edge(1))+current_t*normal(current_edge(2))).Normalized();
        const TV source_X(-source_Z.y,source_Z.x);
        MATRIX<T,2> M_source,M_destination;
        M_source.Column(1)=source_Z;M_source.Column(2)=source_X;
        M_destination.Column(1)=destination_Z;M_destination.Column(2)=destination_X;
        ROTATION<TV> R_source=ROTATION<TV>(M_source);
        ROTATION<TV> R_destination=ROTATION<TV>(M_destination);
        MATRIX<T,2> M_node_1,M_node_2;
        M_node_1.Column(1)=normal(current_edge(1));M_node_1.Column(2)=TV(-normal(current_edge(1)).y,normal(current_edge(1)).x);
        M_node_2.Column(1)=normal(current_edge(2));M_node_2.Column(2)=TV(-normal(current_edge(2)).y,normal(current_edge(2)).x);
        ROTATION<TV> R_node_1=ROTATION<TV>(M_node_1);
        ROTATION<TV> R_node_2=ROTATION<TV>(M_node_2);
        if(rotate_before_barycentric_interpolation){
            Rotate_Helper(Z_ghost_1,R_source*R_node_1.Inverse()*last_rigid_rotation.Inverse());
            Rotate_Helper(Z_ghost_2,R_source*R_node_2.Inverse()*last_rigid_rotation.Inverse());}
        T weight1=(T)1-current_t,weight2=current_t;
        if(current_t<(T)0){weight2=(T)0;weight1=(T)1;}
        if(current_t>(T)1){weight2=(T)1;weight1=(T)0;}
        if(weights_to_cell){
            if(!forward){
                T2 Z_ghost_1_rotated_to_destination=Z_ghost_1;
                T2 Z_ghost_2_rotated_to_destination=Z_ghost_2;
                Rotate_Helper(Z_ghost_1_rotated_to_destination,current_rigid_rotation*R_destination*R_source.Inverse());
                Rotate_Helper(Z_ghost_2_rotated_to_destination,current_rigid_rotation*R_destination*R_source.Inverse());
                weights_to_cell->operator()(current_edge(1)).Append(TRIPLE<int,T,T2>(index,weight1,Z_ghost_1_rotated_to_destination));
                weights_to_cell->operator()(current_edge(2)).Append(TRIPLE<int,T,T2>(index,weight2,Z_ghost_2_rotated_to_destination));}
            else{
                T2 Z_ghost_1_rotated_to_destination=Z_ghost(index);
                T2 Z_ghost_2_rotated_to_destination=Z_ghost(index);
                Rotate_Helper(Z_ghost_1_rotated_to_destination,current_rigid_rotation*R_node_1*R_destination.Inverse()*last_rigid_rotation.Inverse());
                Rotate_Helper(Z_ghost_2_rotated_to_destination,current_rigid_rotation*R_node_2*R_destination.Inverse()*last_rigid_rotation.Inverse());
                weights_to_cell->operator()(current_edge(1)).Append(TRIPLE<int,T,T2>(index,weight1,Z_ghost_1_rotated_to_destination));
                weights_to_cell->operator()(current_edge(2)).Append(TRIPLE<int,T,T2>(index,weight2,Z_ghost_2_rotated_to_destination));}}
        Z(index)=weight1*Z_ghost_1+weight2*Z_ghost_2;
        Rotate_Helper(Z(index),current_rigid_rotation*R_destination*R_source.Inverse());
        
        if(debug)LOG::cout<<"index done="<<index<<std::endl;
    }
}

template<class T>
bool Find_Backtrace_Face(bool force_to_find_an_element,const TRIANGLE_MESH& mesh,const ARRAY_VIEW<VECTOR<T,3> >& X,int close_vertex,const VECTOR<T,3>& current_pos,VECTOR<T,3>& normal_of_close_vertex,T tolerance,RANDOM_NUMBERS<T>& rng,
                        VECTOR<T,3>& current_movement,bool& projected_movement_vanished,VECTOR<int,2>& current_edge,T& current_t,VECTOR<int,3>& current_element,int& current_element_vertex_index)
{
    typedef VECTOR<T,3> TV;
    bool debug=false;//(close_vertex==3759||close_vertex==4636);
    ARRAY<int> incident_elements=mesh.incident_elements->operator()(close_vertex);
    TV projected_movement=current_movement-TV::Dot_Product(normal_of_close_vertex,current_movement)*normal_of_close_vertex;
    if(projected_movement.Magnitude()<tolerance){
        if(debug)LOG::cout<<"projected_movement vanish on vertex="<<projected_movement<<std::endl;
        current_movement=TV();
        projected_movement_vanished=true;
        return true;}
    TV projected_movement_unit=projected_movement;
    projected_movement_unit.Normalize();
    if(debug)LOG::cout<<"just projected (close vertex case) projected_movement="<<projected_movement<<std::endl;
    //LOG::cout<<"test1="<<TV::Dot_Product(normal_of_close_vertex,projected_movement)<<std::endl;
    for(int incident_index=1;incident_index<=incident_elements.Size();incident_index++){
        VECTOR<int,3> incident_element=mesh.elements(incident_elements(incident_index));
        VECTOR<int,3> test_triangle;test_triangle(1)=close_vertex; 
        int current_element_index_n=-1;
        for(int n=1;n<=3;n++){
            if(incident_element(n)==close_vertex){
                test_triangle(2)=incident_element(n%3+1);
                test_triangle(3)=incident_element((n+1)%3+1);
                current_element_index_n=n%3+1;
                break;}}
        VECTOR<TV,3> test_triangle_X;
        for(int n=1;n<=3;n++)test_triangle_X(n)=X(test_triangle(n));
        //LOG::cout<<incident_element<<test_triangle<<test_triangle_X<<std::endl;
        TV edge02=test_triangle_X(2)-current_pos;
        TV edge03=test_triangle_X(3)-current_pos;
        if(TV::Dot_Product(TV::Cross_Product(edge02,edge03).Normalized(),normal_of_close_vertex)<(T)0)continue;
        TV projected_edge02=edge02-TV::Dot_Product(normal_of_close_vertex,edge02)*normal_of_close_vertex;
        TV projected_edge03=edge03-TV::Dot_Product(normal_of_close_vertex,edge03)*normal_of_close_vertex;
        if(debug)LOG::cout<<"cross02="<<TV::Dot_Product(TV::Cross_Product(projected_edge02,projected_movement_unit),normal_of_close_vertex)
                         <<" cross03="<<TV::Dot_Product(TV::Cross_Product(projected_movement_unit,projected_edge03),normal_of_close_vertex)<<std::endl;
        //LOG::cout<<"test2="<<TV::Dot_Product(TV::Cross_Product(projected_edge02,projected_movement_unit),normal_of_close_vertex)<<std::endl;
        //LOG::cout<<"test3="<<TV::Dot_Product(TV::Cross_Product(projected_movement_unit,projected_edge03),normal_of_close_vertex)<<std::endl;
        if(TV::Dot_Product(TV::Cross_Product(projected_edge02,projected_movement_unit),normal_of_close_vertex)<(T)0
        || TV::Dot_Product(TV::Cross_Product(projected_movement_unit,projected_edge03),normal_of_close_vertex)<(T)0){
            //LOG::cout<<"not inside incident element "<<incident_element<<" "<<normal_of_close_vertex<<" "<<projected_movement_unit<<std::endl;
            //LOG::cout<<std::endl;
            continue;}
        else{
            current_edge(1)=test_triangle(2);
            current_edge(2)=close_vertex;
            current_t=1; // normal pulled
            //LOG::cout<<"inside incident element "<<incident_element<<std::endl;
            //current_pos=test_triangle_X(1); // pull to the point DONE disable position pull
            current_element_vertex_index=current_element_index_n;
            current_element=incident_element;
            return true;}
    }
    if(force_to_find_an_element){
        TV face_normal_as_new_vertex_normal;
        LOG::cout<<"forced mode, close_vertex="<<close_vertex<<" X(close_vertex)="<<X(close_vertex)<<std::endl;
        T max_movement_projected_on_face(0);bool found_triangle=false;
        for(int attempt=1;attempt<=3;++attempt){
            LOG::cout<<"attempt "<<attempt<<std::endl;int attempt_3_id=1+int(rng.Get_Number()*(incident_elements.Size()));
            for(int incident_index=1;incident_index<=incident_elements.Size();incident_index++){
                if(attempt==3&&incident_index!=attempt_3_id)continue;
                VECTOR<int,3> incident_element=mesh.elements(incident_elements(incident_index));
                VECTOR<int,3> test_triangle;test_triangle(1)=close_vertex; 
                int current_element_index_n=-1;
                for(int n=1;n<=3;n++){
                    if(incident_element(n)==close_vertex){
                        test_triangle(2)=incident_element(n%3+1);
                        test_triangle(3)=incident_element((n+1)%3+1);
                        current_element_index_n=n%3+1;
                        break;}}
                VECTOR<TV,3> test_triangle_X;
                for(int n=1;n<=3;n++)test_triangle_X(n)=X(test_triangle(n));
                //LOG::cout<<incident_element<<test_triangle<<test_triangle_X<<std::endl;
                TV edge02=test_triangle_X(2)-current_pos;
                TV edge03=test_triangle_X(3)-current_pos;
                TV edge23=test_triangle_X(3)-test_triangle_X(2);
                TV face_normal=TV::Cross_Product(edge02,edge03).Normalized();
                TV ray0to23=TV::Cross_Product(edge23,face_normal).Normalized();
                TV plane_normal=TV::Cross_Product(current_movement,normal_of_close_vertex).Normalized();
                if(attempt==1 && TV::Dot_Product(edge02,plane_normal)*TV::Dot_Product(edge03,plane_normal)>(T)0){continue;}
                else if(attempt==3 || TV::Dot_Product(current_movement,ray0to23)>max_movement_projected_on_face){
                    LOG::cout<<"test_triangle="<<test_triangle<<" face_normal="<<face_normal<<" found triangle in forced mode"<<std::endl;
                    //build tet using current_movement,edge02+edge03, find height ray.
                    TV height_ray=(current_movement.Normalized()-(edge02+edge03).Normalized()).Normalized();
                    face_normal_as_new_vertex_normal=TV::Dot_Product(height_ray,face_normal)>(T)0?height_ray:-height_ray;
                    found_triangle=true;
                    max_movement_projected_on_face=TV::Dot_Product(current_movement,ray0to23);
                    current_edge(1)=test_triangle(2);
                    current_edge(2)=close_vertex;
                    current_t=1;
                    current_element_vertex_index=current_element_index_n;
                    current_element=incident_element;}}
            if(found_triangle){normal_of_close_vertex=face_normal_as_new_vertex_normal;return true;}}}
    return false;
}
template<class TV,class T2> template<int d> typename ENABLE_IF<d==3,void>::TYPE ADVECTION_SEMI_LAGRANGIAN_SURFACE_MESH<TV,T2>::
Update_Advection_Equation_Node(const SURFACE_MESH_TYPE& mesh,
            ARRAY_VIEW<T2>& Z,const ARRAY_VIEW<T2>& Z_ghost,
            const ARRAY_VIEW<VECTOR<T,d> >& X,const ARRAY_VIEW<VECTOR<T,d> >& V,const ARRAY_VIEW<VECTOR<T,d> >& normal,
            const RIGID_GEOMETRY<TV>& rigid_geometry,const FRAME<TV>& last_step_frame,
            const T dt,const T time,ARRAY<ARRAY<TRIPLE<int,T,T2> > >* weights_to_cell,
            bool forward)
{
    ROTATION<TV> current_rigid_rotation=rigid_geometry.Frame().Rotation();
    ROTATION<TV> last_rigid_rotation=last_step_frame.Rotation();
    for(int index=1;index<=Z.Size();index++){
        bool debug=false;//(index==3759||index==4636);
        if(debug)LOG::cout<<"index="<<index<<std::endl;
        TV current_pos=X(index);
        TV mesh_velocity=rigid_geometry.Object_Space_Vector(rigid_geometry.Pointwise_Object_Velocity(rigid_geometry.Frame()*X(index)));
        if(forward)mesh_velocity*=(T)-1;
        TV current_movement=(rigid_geometry.Object_Space_Vector(V(index))-mesh_velocity)*(-dt);
        current_movement=current_movement-TV::Dot_Product(normal(index),current_movement)*normal(index);
        VECTOR<int,3> current_element(0,0,0); 
        VECTOR<int,2> current_edge(index,index);
        T current_t=1;

        const TV destination_Z=normal(index);
        TV destination_X;
        if(current_movement.Magnitude()<tolerance){destination_X=TV(1,0,0);if(TV::Cross_Product(destination_Z,destination_X).Magnitude()<tolerance)destination_X=TV(0,1,0);}
        else destination_X=current_movement.Normalized();
        const TV destination_Y=TV::Cross_Product(destination_Z,destination_X).Normalized();
        destination_X=TV::Cross_Product(destination_Y,destination_Z);
        TV current_rotated_direction=destination_X;

        if(debug)LOG::cout<<"Initial Condition: "<<current_pos<<" "<<current_movement<<" "<<current_movement.Magnitude()<<" "<<current_element<<std::endl;
        
        bool projected_movement_vanished=true;
        int close_vertex=index;
        while(current_movement.Magnitude()>tolerance){
            projected_movement_vanished=false;
            bool found_close_vertex=false;
            int current_element_vertex_index=-1; 
            if(debug)LOG::cout<<"Stepping phase starts"<<std::endl;
            // if current_pos is close to one of the triangle vertices, 
            if(close_vertex<0){
                //T close_tolerance(1e-3);
                if(abs(current_t)<tolerance){
                    close_vertex=current_edge(1);
                    if(debug)LOG::cout<<"Search and find close vertex="<<close_vertex<<std::endl;
                }
                else if(abs((T)1-current_t)<tolerance){
                    close_vertex=current_edge(2);
                    if(debug)LOG::cout<<"Search and find close vertex="<<close_vertex<<std::endl;
                }
            }
            // search one ring of that vertex, choose may be a new triangle to start with
            TV vertex_normal;
            if(close_vertex>0){
                found_close_vertex=true;
                vertex_normal=normal(close_vertex);
                if(!Find_Backtrace_Face(true,mesh,X,close_vertex,current_pos,vertex_normal,tolerance,rng,current_movement,projected_movement_vanished,current_edge,current_t,current_element,current_element_vertex_index)){
                        std::stringstream ss;ss<<"Tracing outside the mesh, close_vertex="<<close_vertex<<" X(close_vertex)="<<X(close_vertex);PHYSBAM_FATAL_ERROR(ss.str());
                }
                if(debug)LOG::cout<<"Close Vertex Rotate: "<<current_pos<<" "<<current_movement<<" "<<current_movement.Magnitude()<<" "<<current_element<<std::endl;
                if(projected_movement_vanished)break;}

            close_vertex=-1;
            //LOG::cout<<"current_element="<<current_element<<X(current_element(1))<<X(current_element(2))<<X(current_element(3))<<std::endl;
            TV edge_normal=found_close_vertex?vertex_normal:((1-current_t)*normal(current_edge(1))+current_t*normal(current_edge(2))).Normalized();

            const TV face_normal=TV::Cross_Product(X(current_element(2))-X(current_element(1)),X(current_element(3))-X(current_element(1))).Normalized();
            const TV current_edge_direction=(X(current_edge(2))-X(current_edge(1))).Normalized();
            if(debug)LOG::cout<<"edge_normal="<<edge_normal<<" current_edge_direction="<<current_edge_direction<<" face_normal="<<face_normal<<" dot_product="<<TV::Dot_Product(edge_normal,face_normal)<<std::endl;
            if(!found_close_vertex && (TV::Dot_Product(current_movement,TV::Cross_Product(current_edge_direction,edge_normal).Normalized())<tolerance || TV::Dot_Product(edge_normal,face_normal)<tolerance)){
                VECTOR<int,3> last_element(0,0,0);
                ARRAY<int> incident_elements=mesh.incident_elements->operator()(current_edge(2));
                bool must_found_a_triangle=false;
                for(int incident_index=1;incident_index<=incident_elements.Size();incident_index++){
                    VECTOR<int,3> incident_element=mesh.elements(incident_elements(incident_index));
                    int i;for(i=1;i<=3;i++)
                        if(incident_element(i)==current_edge(1))break;
                    if(current_edge(2)==incident_element(i%3+1)){
                        last_element=incident_element;
                        must_found_a_triangle=true;
                        break;}}
                if(!must_found_a_triangle)PHYSBAM_FATAL_ERROR("no last element!");
                const TV last_face_normal=TV::Cross_Product(X(last_element(2))-X(last_element(1)),X(last_element(3))-X(last_element(1))).Normalized();

                if(debug)LOG::cout<<"use averaged edge normal (before: "<<edge_normal;
                TV tangent_part_edge_normal=current_edge_direction*TV::Dot_Product(edge_normal,current_edge_direction);
                edge_normal=tangent_part_edge_normal+(edge_normal-tangent_part_edge_normal).Magnitude()*((face_normal+last_face_normal).Normalized());
                if(debug)LOG::cout<<" after: "<<edge_normal<<") for current_edge="<<current_edge<<" X(current_edge(1))="<<X(current_edge(1))<<std::endl;
                if(TV::Dot_Product(edge_normal,face_normal)<tolerance)PHYSBAM_FATAL_ERROR("bad edge normal");
            }

            TV projected_movement=current_movement-TV::Dot_Product(edge_normal,current_movement)*edge_normal;
            if(projected_movement.Magnitude()<tolerance){
                if(debug)LOG::cout<<"projected_movement vanish on edge="<<projected_movement<<std::endl;
                current_movement=TV();
                projected_movement_vanished=true;
                break;
            }
            bool must_intersect_with_an_edge=false;
            for(int n=1;n<=3;n++){
                if(found_close_vertex && current_element_vertex_index!=n) continue;
                if(!found_close_vertex && current_element(n%3+1)==current_edge(1) && current_element(n)==current_edge(2)){if(debug)LOG::cout<<"regular edge: skip intersecting dup edge"<<std::endl;continue;}
                // DONE read current_edge here, and project everyone by that for s and t
                if(debug)LOG::cout<<"Intersect segment "<<current_pos<<" "<<current_pos+current_movement<<" "<<X(current_element(n))<<" " <<X(current_element(n%3+1))<<std::endl;
                TV edge12=X(current_element(n))-current_pos;
                TV edge13=X(current_element(n%3+1))-current_pos;
                TV projected_edge12=edge12-TV::Dot_Product(edge_normal,edge12)*edge_normal;
                TV projected_edge13=edge13-TV::Dot_Product(edge_normal,edge13)*edge_normal;
                if(debug)LOG::cout<<"just projected projected_movement="<<projected_movement<<std::endl;
                T s,t;bool too_short=false;
                bool parallel=Intersection_Segments(edge_normal,TV::Constant_Vector((T)0),projected_movement,projected_edge12,projected_edge13,s,t,too_short);
                if(debug)LOG::cout<<"Intersect segment n="<<n<<" s="<<s<<" t="<<t<<std::endl;
                if(!found_close_vertex && !too_short){
                    if(parallel){if(debug){LOG::cout<<"skipped parallel edge"<<std::endl;}continue;}
                    if(t<-tolerance||t>(T)1+tolerance){continue;}
                }
                t=clamp(t,(T)0,(T)1);
                //LOG::cout<<"Intersect segment clamped n="<<n<<" s="<<s<<" t="<<t<<std::endl;
                //if(s<(T)0 && !too_short){std::stringstream ss;ss<<"negative s! current_pos="<<current_pos<<" current_edge="<<current_edge<<" index="<<index<<" X(index)="<<X(index);PHYSBAM_FATAL_ERROR(ss.str());}
                must_intersect_with_an_edge=true;
                if(s>(T)0 && s<=(T)1 && !too_short){
                    // move to next element
                    // update pos movement triangle
                    current_pos=X(current_element(n))*((T)1-t)+X(current_element(n%3+1))*t;
                    if(debug)LOG::cout<<"current pos is "<<current_pos<<", original starting point is "<<X(index)<<std::endl;
                    current_movement=current_movement*((T)1-s);
                    ARRAY<int> incident_elements=mesh.incident_elements->operator()(current_element(n));
                    VECTOR<int,3> last_element=current_element;
                    bool must_found_a_triangle=false;
                    for(int incident_index=1;incident_index<=incident_elements.Size();incident_index++){
                        VECTOR<int,3> incident_element=mesh.elements(incident_elements(incident_index));
                        int i;for(i=1;i<=3;i++)
                            if(incident_element(i)==current_element(n))break;
                        if(current_element(n%3+1)==incident_element((i+1)%3+1)){
                            if(debug)LOG::cout<<"move element from "<<current_element<<" to "<<incident_element<<std::endl;
                            current_element=incident_element;
                            must_found_a_triangle=true;
                            break;
                        }
                    }
                    assert(must_found_a_triangle);
                    current_edge(1)=last_element(n);
                    current_edge(2)=last_element(n%3+1);
                    current_t=t;
                    if(debug)LOG::cout<<"current_edge="<<current_edge<<" current_t="<<current_t<<std::endl;
                    TV rotated_direction=(1-t)*edge12+t*edge13; // TODO in the close vertex case, project to the new current_edge plane...
                    rotated_direction.Normalize();
                    current_rotated_direction=rotated_direction;
                    current_movement=current_movement.Magnitude()*rotated_direction;
                    if(debug)LOG::cout<<"Move Rotate: "<<current_pos<<" "<<current_movement<<" "<<current_movement.Magnitude()<<" "<<current_element<<std::endl;
                    break;}
                else{
                    if(s>(T)0)current_pos=current_pos+current_movement;
                    current_movement=TV();
                    if(debug)LOG::cout<<"Move Done: "<<current_pos<<" "<<current_movement<<" "<<current_movement.Magnitude()<<" "<<current_element<<std::endl;
                    break;
                }}
            if(!must_intersect_with_an_edge){std::stringstream ss;ss<<"did not intersect with an edge! current_pos="<<current_pos<<" current_edge="<<current_edge<<" index="<<index<<" X(index)="<<X(index);PHYSBAM_FATAL_ERROR(ss.str());}
            assert(must_intersect_with_an_edge);
            if(debug)LOG::cout<<"Stepping phase ends"<<std::endl;
        }
        // barycentric weigth, interpolation, done
        TV source_normal;
        TV barycentric;
        ARRAY<T2> Z_ghost_node;ARRAY<int> node_index;
        // shadow, normal, source Z's before rotation
        if(projected_movement_vanished){
            if(shadow_X) (*shadow_X)(index)=current_pos;
            if(debug)LOG::cout<<"index="<<index<<" t="<<current_t<<" "<<Z_ghost(current_edge(1))<<" "<<Z_ghost(current_edge(2))<<std::endl;
            assert(current_t>-tolerance && current_t-(T)1<tolerance);
            TV edge_normal=(1-current_t)*normal(current_edge(1))+current_t*normal(current_edge(2));edge_normal.Normalize();
            source_normal=edge_normal;
            Z_ghost_node.Append(Z_ghost(current_edge(1)));Z_ghost_node.Append(Z_ghost(current_edge(2)));
            node_index.Append(current_edge(1));node_index.Append(current_edge(2));
        }
        else{
            if(shadow_X) (*shadow_X)(index)=current_pos;
            barycentric=TRIANGLE_3D<T>::Barycentric_Coordinates(
                    current_pos,X(current_element(1)),X(current_element(2)),X(current_element(3)));
            if(debug)LOG::cout<<"index="<<index<<" Barycentri Weight="<<barycentric<<" "<<Z_ghost(current_element(1))<<" "<<Z_ghost(current_element(2))<<" "<<Z_ghost(current_element(3))<<std::endl;
            if(debug)LOG::cout<<"current_pos "<<current_pos<<" "<<-X(current_element(1))*barycentric(1)-X(current_element(2))*barycentric(2)-X(current_element(3))*barycentric(3)<<std::endl;
            //assert((current_pos-X(current_element(1))*barycentric(1)-X(current_element(2))*barycentric(2)-X(current_element(3))*barycentric(3)).Magnitude()<tolerance);
            TV face_normal=normal(current_element(1))*barycentric(1)+normal(current_element(2))*barycentric(2)+normal(current_element(3))*barycentric(3);face_normal.Normalize();
            source_normal=face_normal;
            Z_ghost_node.Append(Z_ghost(current_element(1)));Z_ghost_node.Append(Z_ghost(current_element(2)));Z_ghost_node.Append(Z_ghost(current_element(3)));
            node_index.Append(current_element(1));node_index.Append(current_element(2));node_index.Append(current_element(3));
        }
        ARRAY<ROTATION<TV> > R_center(Z_ghost_node.Size());
        ARRAY<ROTATION<TV> > R_node(Z_ghost_node.Size());
        if(rotate_before_barycentric_interpolation) for(int i=1;i<=Z_ghost_node.Size();++i){
            const TV center_Z=source_normal;
            TV center_X=(X(node_index(i))-current_pos).Normalized(); // TODO use in-simplex position instead of current_pos
            const TV center_Y=TV::Cross_Product(center_Z,center_X).Normalized();
            center_X=TV::Cross_Product(center_Y,center_Z);
            const TV node_Z=normal(node_index(i));
            TV node_X=(X(node_index(i))-current_pos).Normalized();
            const TV node_Y=TV::Cross_Product(node_Z,node_X).Normalized();
            node_X=TV::Cross_Product(node_Y,node_Z);
            MATRIX<T,3> M_center,M_node;
            M_center.Column(1)=center_X;M_center.Column(2)=center_Y;M_center.Column(3)=center_Z;
            M_node.Column(1)=node_X;M_node.Column(2)=node_Y;M_node.Column(3)=node_Z;
            R_center(i)=ROTATION<TV>(M_center);
            R_node(i)=ROTATION<TV>(M_node);
            Rotate_Helper(Z_ghost_node(i),R_center(i)*R_node(i).Inverse()*last_rigid_rotation.Inverse());
        }
        const TV source_Z=source_normal;
        TV source_X=current_rotated_direction.Normalized();
        const TV source_Y=TV::Cross_Product(source_Z,source_X).Normalized();
        source_X=TV::Cross_Product(source_Y,source_Z);
        MATRIX<T,3> M_source,M_destination;
        M_source.Column(1)=source_X;M_source.Column(2)=source_Y;M_source.Column(3)=source_Z;
        M_destination.Column(1)=destination_X;M_destination.Column(2)=destination_Y;M_destination.Column(3)=destination_Z;
        ROTATION<TV> R_source=ROTATION<TV>(M_source);
        ROTATION<TV> R_destination=ROTATION<TV>(M_destination);

        if(projected_movement_vanished){
            T weight1=(T)1-current_t,weight2=current_t;
            if(current_t<(T)0){weight2=(T)0;weight1=(T)1;}
            if(current_t>(T)1){weight2=(T)1;weight1=(T)0;}
            Z(index)=weight1*Z_ghost_node(1)+weight2*Z_ghost_node(2);
            if(weights_to_cell){
                if(!forward){
                    T2 z_ghost_node_1_rotated_to_destination=Z_ghost_node(1);
                    T2 z_ghost_node_2_rotated_to_destination=Z_ghost_node(2);
                    Rotate_Helper(z_ghost_node_1_rotated_to_destination,current_rigid_rotation*R_destination*R_source.Inverse());
                    Rotate_Helper(z_ghost_node_2_rotated_to_destination,current_rigid_rotation*R_destination*R_source.Inverse());
                    weights_to_cell->operator()(current_edge(1)).Append(TRIPLE<int,T,T2>(index,weight1,z_ghost_node_1_rotated_to_destination));
                    weights_to_cell->operator()(current_edge(2)).Append(TRIPLE<int,T,T2>(index,weight2,z_ghost_node_2_rotated_to_destination));}
                else{
                    T2 z_ghost_node_1_rotated_to_destination=Z_ghost(index);
                    T2 z_ghost_node_2_rotated_to_destination=Z_ghost(index);
                    Rotate_Helper(z_ghost_node_1_rotated_to_destination,current_rigid_rotation*R_node(1)*R_center(1).Inverse()*R_source*R_destination.Inverse()*last_rigid_rotation.Inverse());
                    Rotate_Helper(z_ghost_node_2_rotated_to_destination,current_rigid_rotation*R_node(2)*R_center(2).Inverse()*R_source*R_destination.Inverse()*last_rigid_rotation.Inverse());
                    weights_to_cell->operator()(current_edge(1)).Append(TRIPLE<int,T,T2>(index,weight1,z_ghost_node_1_rotated_to_destination));
                    weights_to_cell->operator()(current_edge(2)).Append(TRIPLE<int,T,T2>(index,weight2,z_ghost_node_2_rotated_to_destination));}}}
        else{
            T weight1=barycentric(1),weight2=barycentric(2),weight3=barycentric(3);
            weight1=clamp(weight1,(T)0,(T)1);weight2=clamp(weight2,(T)0,(T)1);weight3=clamp(weight3,(T)0,(T)1);
            T sum_weight=weight1+weight2+weight3;
            weight1/=sum_weight;weight2/=sum_weight;weight3/=sum_weight;
            Z(index)=Z_ghost_node(1)*weight1+Z_ghost_node(2)*weight2+Z_ghost_node(3)*weight3;
            if(weights_to_cell){
                if(!forward){
                    T2 z_ghost_node_1_rotated_to_destination=Z_ghost_node(1);
                    T2 z_ghost_node_2_rotated_to_destination=Z_ghost_node(2);
                    T2 z_ghost_node_3_rotated_to_destination=Z_ghost_node(3);
                    Rotate_Helper(z_ghost_node_1_rotated_to_destination,current_rigid_rotation*R_destination*R_source.Inverse());
                    Rotate_Helper(z_ghost_node_2_rotated_to_destination,current_rigid_rotation*R_destination*R_source.Inverse());
                    Rotate_Helper(z_ghost_node_3_rotated_to_destination,current_rigid_rotation*R_destination*R_source.Inverse());
                    weights_to_cell->operator()(current_element(1)).Append(TRIPLE<int,T,T2>(index,weight1,z_ghost_node_1_rotated_to_destination));
                    weights_to_cell->operator()(current_element(2)).Append(TRIPLE<int,T,T2>(index,weight2,z_ghost_node_2_rotated_to_destination));
                    weights_to_cell->operator()(current_element(3)).Append(TRIPLE<int,T,T2>(index,weight3,z_ghost_node_3_rotated_to_destination));}
                else{
                    T2 z_ghost_node_1_rotated_to_destination=Z_ghost(index);
                    T2 z_ghost_node_2_rotated_to_destination=Z_ghost(index);
                    T2 z_ghost_node_3_rotated_to_destination=Z_ghost(index);
                    Rotate_Helper(z_ghost_node_1_rotated_to_destination,current_rigid_rotation*R_node(1)*R_center(1).Inverse()*R_source*R_destination.Inverse()*last_rigid_rotation.Inverse());
                    Rotate_Helper(z_ghost_node_2_rotated_to_destination,current_rigid_rotation*R_node(2)*R_center(2).Inverse()*R_source*R_destination.Inverse()*last_rigid_rotation.Inverse());
                    Rotate_Helper(z_ghost_node_3_rotated_to_destination,current_rigid_rotation*R_node(3)*R_center(3).Inverse()*R_source*R_destination.Inverse()*last_rigid_rotation.Inverse());
                    weights_to_cell->operator()(current_element(1)).Append(TRIPLE<int,T,T2>(index,weight1,z_ghost_node_1_rotated_to_destination));
                    weights_to_cell->operator()(current_element(2)).Append(TRIPLE<int,T,T2>(index,weight2,z_ghost_node_2_rotated_to_destination));
                    weights_to_cell->operator()(current_element(3)).Append(TRIPLE<int,T,T2>(index,weight3,z_ghost_node_3_rotated_to_destination));}
                }
            }

        Rotate_Helper(Z(index),current_rigid_rotation*R_destination*R_source.Inverse());

        if(debug)LOG::cout<<"index done="<<index<<std::endl;
    }
}
#define INSTANTIATION_HELPER(T,d) \
template class ADVECTION_SEMI_LAGRANGIAN_SURFACE_MESH<VECTOR<T,d>,T>; \
template class ADVECTION_SEMI_LAGRANGIAN_SURFACE_MESH<VECTOR<T,d>,VECTOR<T,d> >; \
template  void ADVECTION_SEMI_LAGRANGIAN_SURFACE_MESH<VECTOR<T,d>,T           >::Update_Advection_Equation_Node<d>(const SURFACE_MESH_TYPE&, ARRAY_VIEW<T           >&,const ARRAY_VIEW<T           >&, const ARRAY_VIEW<VECTOR<T,d> >&,const ARRAY_VIEW<VECTOR<T,d> >&,const ARRAY_VIEW<VECTOR<T,d> >&,const RIGID_GEOMETRY<VECTOR<T,d> >&,const FRAME<VECTOR<T,d> >&,const T,const T,ARRAY<ARRAY<TRIPLE<int,T,T> > >*,bool); \
template  void ADVECTION_SEMI_LAGRANGIAN_SURFACE_MESH<VECTOR<T,d>,VECTOR<T,d> >::Update_Advection_Equation_Node<d>(const SURFACE_MESH_TYPE&, ARRAY_VIEW<VECTOR<T,d> >&,const ARRAY_VIEW<VECTOR<T,d> >&, const ARRAY_VIEW<VECTOR<T,d> >&,const ARRAY_VIEW<VECTOR<T,d> >&,const ARRAY_VIEW<VECTOR<T,d> >&,const RIGID_GEOMETRY<VECTOR<T,d> >&,const FRAME<VECTOR<T,d> >&,const T,const T,ARRAY<ARRAY<TRIPLE<int,T,VECTOR<T,d> > > >*,bool);

INSTANTIATION_HELPER(float,2)
INSTANTIATION_HELPER(float,3)
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
INSTANTIATION_HELPER(double,2)
INSTANTIATION_HELPER(double,3)
#endif
}
