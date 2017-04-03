//#####################################################################
// Copyright 2015, Mridul Aanjaneya, Ming Gao.
// This file is part of PhysBAM whose distribution is governed by the license 
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Common_Tools/Data_Structures/QUADTREE.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Geometry/Level_Sets/LEVELSET_UTILITIES.h>
#include <PhysBAM_Geometry/Basic_Geometry/SEGMENT_2D.h>
#include <PhysBAM_Geometry/Basic_Geometry/TRIANGLE_2D.h>
#include <PhysBAM_Geometry/Basic_Geometry/POLYGON.h>
#include <queue>
using namespace PhysBAM;
//#####################################################################
// Compute_Two_Level_Subdivision
//#####################################################################
template<class T> void QUADTREE<T>::
Compute_Two_Level_Subdivision(const RANGE<TV>& domain_input,const int scale,bool uniform_input)
{
    uniform=uniform_input;
    domain=domain_input;
    int resolution=1;TV dx=(T).5*domain.Edge_Lengths();
    QUADTREE_NODE<T> *node=new QUADTREE_NODE<T>(domain.Center(),dx);
    nodes.Append(node);node->index=1;
    
    std::queue<QUADTREE_NODE<T>*> list;list.push(node);
    while(resolution<scale){resolution*=2;TV current_dx=dx/resolution;
        int counter=0;
        // subdivide the quadtree
        while(counter<sqr(resolution/2)){++counter;
            QUADTREE_NODE<T> *top=list.front();
            list.pop();
            // create four new corners
            TV north_west_corner=top->center+TV(-current_dx(1),current_dx(2));
            TV north_east_corner=top->center+TV(current_dx(1),current_dx(2));
            TV south_west_corner=top->center+TV(-current_dx(1),-current_dx(2));
            TV south_east_corner=top->center+TV(current_dx(1),-current_dx(2));
            // create four new quadtree nodes
            QUADTREE_NODE<T> *north_east=new QUADTREE_NODE<T>(north_east_corner,current_dx);north_east->parent_index=top->index;
            QUADTREE_NODE<T> *north_west=new QUADTREE_NODE<T>(north_west_corner,current_dx);north_west->parent_index=top->index;
            QUADTREE_NODE<T> *south_east=new QUADTREE_NODE<T>(south_east_corner,current_dx);south_east->parent_index=top->index;
            QUADTREE_NODE<T> *south_west=new QUADTREE_NODE<T>(south_west_corner,current_dx);south_west->parent_index=top->index;
            // reassign pointers
            top->north_east_index=nodes.Size()+1;
            top->north_west_index=nodes.Size()+2;
            top->south_east_index=nodes.Size()+3;
            top->south_west_index=nodes.Size()+4;
            // append in nodes
            nodes.Append(north_east);north_east->index=nodes.Size();
            nodes.Append(north_west);north_west->index=nodes.Size();
            nodes.Append(south_east);south_east->index=nodes.Size();
            nodes.Append(south_west);south_west->index=nodes.Size();
            // push in queue
            list.push(north_east);
            list.push(north_west);
            list.push(south_east);
            list.push(south_west);}}





    // list only contains cells of given resolution; subdivide the left side


    // specific to uniform boundary cells
    T minx=1.,maxx=0.,miny=1.,maxy=0.;
    for(int i=1;i<=nodes.Size();++i)
        if(nodes(i)->north_east_index==0){
            if(nodes(i)->center.x>maxx) maxx=nodes(i)->center.x;
            if(nodes(i)->center.x<minx) minx=nodes(i)->center.x;
            if(nodes(i)->center.y>maxy) maxy=nodes(i)->center.y;
            if(nodes(i)->center.y<miny) miny=nodes(i)->center.y;
        }
    //LOG::cout << minx << "," << maxx << "," << miny << "," << maxy << std::endl;
    ghost.Exact_Resize(nodes.Size());
    for(int i=1;i<=nodes.Size();++i)
    {
        ghost(i)=false;
        if(nodes(i)->north_east_index==0)
            if(nodes(i)->center.x==maxx || nodes(i)->center.x==minx || nodes(i)->center.y==maxy || nodes(i)->center.y==miny) ghost(i)=true;
    }
    //for(int i=1;i<=nodes.Size();++i) if(ghost(i)) LOG::cout << "ghost " << i << std::endl;
    
    if(!uniform){
        while(!list.empty()){QUADTREE_NODE<T> *top=list.front();
            TV current_dx=top->dx_over_two/2;
            list.pop();
            //if( !ghost(top->index) && ((top->center.x < node->center.x &&  top->center.y < node->center.y) || (top->center.x > node->center.x && top->center.y > node->center.y)))
            //if( (top->center.x < node->center.x&& top->center.y < node->center.y) || (top->center.x > node->center.x && top->center.y > node->center.y))
            if(!ghost(top->index) && std::fabs((top->center-TV(.5,.5)).Magnitude())>0.25111111111)
            {
                // create four new corners
                TV north_west_corner=top->center+TV(-current_dx(1),current_dx(2));
                TV north_east_corner=top->center+TV(current_dx(1),current_dx(2));
                TV south_west_corner=top->center+TV(-current_dx(1),-current_dx(2));
                TV south_east_corner=top->center+TV(current_dx(1),-current_dx(2));
                // create four new quadtree nodes
                QUADTREE_NODE<T> *north_east=new QUADTREE_NODE<T>(north_east_corner,current_dx);north_east->parent_index=top->index;
                QUADTREE_NODE<T> *north_west=new QUADTREE_NODE<T>(north_west_corner,current_dx);north_west->parent_index=top->index;
                QUADTREE_NODE<T> *south_east=new QUADTREE_NODE<T>(south_east_corner,current_dx);south_east->parent_index=top->index;
                QUADTREE_NODE<T> *south_west=new QUADTREE_NODE<T>(south_west_corner,current_dx);south_west->parent_index=top->index;
                // reassign pointers
                top->north_east_index=nodes.Size()+1;
                top->north_west_index=nodes.Size()+2;
                top->south_east_index=nodes.Size()+3;
                top->south_west_index=nodes.Size()+4;
                // append in nodes
                nodes.Append(north_east);north_east->index=nodes.Size();
                nodes.Append(north_west);north_west->index=nodes.Size();
                nodes.Append(south_east);south_east->index=nodes.Size();
                nodes.Append(south_west);south_west->index=nodes.Size();}}
    }

    int ghost_size=ghost.Size();
    for(int i=1;i<=nodes.Size()-ghost_size;++i)
        ghost.Append(false);

}
//#####################################################################
// Compute_Neighbors
//#####################################################################
template<class T> void QUADTREE<T>::
Compute_Neighbors()
{
   for(int i=1;i<=nodes.Size();++i)
    {
        ARRAY<PAIR<int,int> > indices(4);
        indices(1)=East_Neighbor_Leaf(i);
        indices(2)=North_Neighbor_Leaf(i);
        indices(3)=West_Neighbor_Leaf(i);
        indices(4)=South_Neighbor_Leaf(i);
        //indices(5)=East_Neighbor(i)==0?0:North_Neighbor(East_Neighbor(i));
        //indices(6)=West_Neighbor(i)==0?0:North_Neighbor(West_Neighbor(i));
        //indices(7)=East_Neighbor(i)==0?0:South_Neighbor(East_Neighbor(i));
        //indices(8)=West_Neighbor(i)==0?0:South_Neighbor(West_Neighbor(i));
        neighbor_indices.Append(indices);
    }
}
//#####################################################################
// North_Neighbor
//#####################################################################
template<class T> int QUADTREE<T>::
North_Neighbor_Same_Level(const int index) const
{
    if(index==1) return 0;
    int parent_index=nodes(index)->parent_index;
    if(index==nodes(parent_index)->south_west_index) return nodes(parent_index)->north_west_index;
    if(index==nodes(parent_index)->south_east_index) return nodes(parent_index)->north_east_index;
    int u=North_Neighbor_Same_Level(parent_index);
    if(u==0 || nodes(u)->north_east_index==0) return u;
    else if(index==nodes(parent_index)->north_west_index) return nodes(u)->south_west_index;
         else return nodes(u)->south_east_index;
    return 0;
}
template<class T> PAIR<int,int> QUADTREE<T>::
North_Neighbor_Leaf(const int index)
{
    if(nodes(index)->north_east_index!=0) return PAIR<int,int>(-1,-1);
    int neighbor=North_Neighbor_Same_Level(index);
    if(neighbor==0) return PAIR<int,int>(0,0);
    if(nodes(neighbor)->north_east_index==0) return PAIR<int,int>(0,neighbor);
    return PAIR<int,int>(nodes(neighbor)->south_east_index,nodes(neighbor)->south_west_index);
}
//#####################################################################
// South_Neighbor
//#####################################################################
template<class T> int QUADTREE<T>::
South_Neighbor_Same_Level(const int index) const
{
    if(index==1) return 0;
    int parent_index=nodes(index)->parent_index;
    if(index==nodes(parent_index)->north_west_index) return nodes(parent_index)->south_west_index;
    if(index==nodes(parent_index)->north_east_index) return nodes(parent_index)->south_east_index;
    int u=South_Neighbor_Same_Level(parent_index);
    if(u==0 || nodes(u)->south_east_index==0) return u;
    else if(index==nodes(parent_index)->south_west_index) return nodes(u)->north_west_index;
         else return nodes(u)->north_east_index;
    return 0;
}
template<class T> PAIR<int,int> QUADTREE<T>::
South_Neighbor_Leaf(const int index)
{
    if(nodes(index)->south_east_index!=0) return PAIR<int,int>(-1,-1);
    int neighbor=South_Neighbor_Same_Level(index);
    if(neighbor==0) return PAIR<int,int>(0,0);
    if(nodes(neighbor)->south_east_index==0) return PAIR<int,int>(0,neighbor);
    return PAIR<int,int>(nodes(neighbor)->north_west_index,nodes(neighbor)->north_east_index);
}
//#####################################################################
// West_Neighbor
//#####################################################################
template<class T> int QUADTREE<T>::
West_Neighbor_Same_Level(const int index) const
{
    if(index==1) return 0;
    int parent_index=nodes(index)->parent_index;
    if(index==nodes(parent_index)->north_east_index) return nodes(parent_index)->north_west_index;
    if(index==nodes(parent_index)->south_east_index) return nodes(parent_index)->south_west_index;
    int u=West_Neighbor_Same_Level(parent_index);
    if(u==0 || nodes(u)->north_east_index==0) return u;
    else if(index==nodes(parent_index)->north_west_index) return nodes(u)->north_east_index;
         else return nodes(u)->south_east_index;
    return 0;
}
template<class T> PAIR<int,int> QUADTREE<T>::
West_Neighbor_Leaf(const int index)
{
    if(nodes(index)->north_east_index!=0) return PAIR<int,int>(-1,-1);
    int neighbor=West_Neighbor_Same_Level(index);
    if(neighbor==0) return PAIR<int,int>(0,0);
    if(nodes(neighbor)->north_east_index==0) return PAIR<int,int>(0,neighbor);
    return PAIR<int,int>(nodes(neighbor)->north_east_index,nodes(neighbor)->south_east_index);
}
//#####################################################################
// East_Neighbor
//#####################################################################
template<class T> int QUADTREE<T>::
East_Neighbor_Same_Level(const int index) const
{
    if(index==1) return 0;
    int parent_index=nodes(index)->parent_index;
    if(index==nodes(parent_index)->north_west_index) return nodes(parent_index)->north_east_index;
    if(index==nodes(parent_index)->south_west_index) return nodes(parent_index)->south_east_index;
    int u=East_Neighbor_Same_Level(parent_index);
    if(u==0 || nodes(u)->north_west_index==0) return u;
    else if(index==nodes(parent_index)->north_east_index) return nodes(u)->north_west_index;
         else return nodes(u)->south_west_index;
    return 0;
}
template<class T> PAIR<int,int> QUADTREE<T>::
East_Neighbor_Leaf(const int index)
{
    if(nodes(index)->north_east_index!=0) return PAIR<int,int>(-1,-1);
    int neighbor=East_Neighbor_Same_Level(index);
    if(neighbor==0) return PAIR<int,int>(0,0);
    if(nodes(neighbor)->north_east_index==0) return PAIR<int,int>(0,neighbor);
    return PAIR<int,int>(nodes(neighbor)->south_west_index,nodes(neighbor)->north_west_index);
}
//#####################################################################
// Check_Child_Number
//#####################################################################
template<class T> int QUADTREE<T>::
Check_Child_Number(int index) const
{
    int parent_index=nodes(index)->parent_index;
    assert(parent_index>0);
    if(nodes(parent_index)->north_east_index==index) return 1;
    if(nodes(parent_index)->north_west_index==index) return 2;
    if(nodes(parent_index)->south_west_index==index) return 3;
    if(nodes(parent_index)->south_east_index==index) return 4;
    return 0;
}
//#####################################################################
// Corner_Neighbor
//#####################################################################
template<class T> int QUADTREE<T>::
Corner_Neighbor_Leaf_Coarse(int index) const
{
    int parent_index=nodes(index)->parent_index;assert(parent_index>0);
    int index2;
    if(nodes(parent_index)->north_east_index==index){index2=North_Neighbor_Same_Level(East_Neighbor_Same_Level(parent_index));
        if(nodes(index2)->north_west_index==0) return index2;else return nodes(index2)->south_west_index;}
    if(nodes(parent_index)->north_west_index==index){index2=North_Neighbor_Same_Level(West_Neighbor_Same_Level(parent_index));
        if(nodes(index2)->north_west_index==0) return index2;else return nodes(index2)->south_east_index;}
    if(nodes(parent_index)->south_east_index==index){index2=South_Neighbor_Same_Level(East_Neighbor_Same_Level(parent_index));
        if(nodes(index2)->north_west_index==0) return index2;else return nodes(index2)->north_west_index;}
    if(nodes(parent_index)->south_west_index==index){index2=South_Neighbor_Same_Level(West_Neighbor_Same_Level(parent_index));
        if(nodes(index2)->north_west_index==0) return index2;else return nodes(index2)->north_east_index;}
    return -1;
}
//#####################################################################
// Tree_Index_2_Array_Index
//#####################################################################
template<class T> int QUADTREE<T>::
Tree_Index_2_Array_Index(int tree_index)
{
    if(tree_index_2_array_index.Size()==0)
    {
        Array_Index_2_Tree_Index(1);
        for(int i=1;i<=array_index_2_tree_index.Size();++i)
            tree_index_2_array_index.Insert(Array_Index_2_Tree_Index(i),i);
    }
    return tree_index_2_array_index.Get(tree_index);
}
//#####################################################################
// Tree_Index_2_Array_Index
//#####################################################################
template<class T> int QUADTREE<T>::
Array_Index_2_Tree_Index(int array_index)
{
    if(array_index_2_tree_index.Size()==0)
    {
        for(int i=1;i<=nodes.Size();++i)
            if(nodes(i)->north_east_index==0 && !ghost(i))
                array_index_2_tree_index.Append(nodes(i)->index);
    }
    return array_index_2_tree_index(array_index);
}
//#####################################################################
// Multiply
//#####################################################################
template<class T> void QUADTREE<T>::
Multiply(const ARRAY<T>& x,ARRAY<T>& b) 
{
    if(uniform) Multiply_Uniform(x,b);
    else Multiply_Non_Uniform(x,b);
}
template<class T> void QUADTREE<T>::
Multiply_Uniform(const ARRAY<T>& x,ARRAY<T>& b) // Ax=b
{
    //b.Exact_Resize(x.Size());
    for(int j1=1;j1<=x.Size();++j1)
    {
        b(j1)=0;
        T p1=x(j1),p2;
        int index1=Array_Index_2_Tree_Index(j1);
        if(nodes(index1)->levelset>0 && nodes(index1)->Neumann_of_cell==false) continue;
        if(ghost(index1)){LOG::cout<<"sth is wrong!should not be ghost!"<<std::endl;exit(0);}
        T dx=nodes(index1)->dx_over_two(1)*2;
        for(int i=1;i<=4;++i){
            T scale=1.;
            T distance=dx;
            int index2=neighbor_indices(index1)(i).y;// for uniform we know y is the real index
            if(ghost(index2)) p2=0.;
            else{
                int j2=Tree_Index_2_Array_Index(index2);
                if(nodes(index2)->levelset>0 && nodes(index2)->Neumann_of_cell==false){
                    p2=p1;
                }
                else{
                    p2=x(j2);
                    scale=std::fabs(nodes(index1)->dof_theta_of_edge(i));
                }
            }
            T g=(p2-p1)/distance;                
            b(j1)+=g*dx*scale;
        }
    }
}
template<class T> void QUADTREE<T>::
Multiply_Non_Uniform(const ARRAY<T>& x,ARRAY<T>& b) // Ax=b
{
    b.Exact_Resize(x.Size());
    for(int j=1;j<=x.Size();++j)
    {
        b(j)=0;
        T p1=x(j),p2;
        int index=Array_Index_2_Tree_Index(j);
        if(nodes(index)->levelset>0 && nodes(index)->Neumann_of_cell==false) continue;
        if(ghost(index)){LOG::cout<<"sth is wrong!should be no ghost!"<<std::endl;exit(0);}
        int counter=1;
        for(int i=1;i<=4;++i){
            int index1=neighbor_indices(index)(i).x;
            int index2=neighbor_indices(index)(i).y;
            T distance=(nodes(index)->center-nodes(index2)->center).Magnitude();
            int j1,j2;
            T scale=1.;
            if(index1!=0){
                if(ghost(index1)) p2=0.;
                else{
                    j1=Tree_Index_2_Array_Index(index1);
                    if(nodes(index1)->levelset>0 && nodes(index1)->Neumann_of_cell==false){
                        p2=p1;
                    }
                    else
                    {
                        p2=x(j1);
                        scale=std::fabs(nodes(index)->dof_theta_of_edge(counter));
                    }

                }
                T g=(p2-p1)/distance;                
                b(j)+=g*nodes(index)->area_of_face(counter++)*scale;
                scale=1.;
                if(ghost(index2)) p2=0.;
                else{
                    j2=Tree_Index_2_Array_Index(index2);
                    if(nodes(index2)->levelset>0 && nodes(index2)->Neumann_of_cell==false){
                        p2=p1;
                    }
                    else{
                        p2=x(j2);
                        scale=std::fabs(nodes(index)->dof_theta_of_edge(counter));
                    }
                }
                g=(p2-p1)/distance;
                b(j)+=g*nodes(index)->area_of_face(counter++)*scale;
            }
            else{
                if(ghost(index2)) p2=0.;
                else{
                    j2=Tree_Index_2_Array_Index(index2);
                    if(nodes(index2)->levelset>0 && nodes(index2)->Neumann_of_cell==false){
                        p2=p1;
                    }
                    else{
                        p2=x(j2);
                        scale=std::fabs(nodes(index)->dof_theta_of_edge(counter));
                    }
                }
                T g=(p2-p1)/distance;
                b(j)+=g*nodes(index)->area_of_face(counter++)*scale;
            }

        }
        assert(counter==nodes(index)->area_of_face.Size()+1 && index);
    }
}
//#####################################################################
// Project
//#####################################################################
template<class T> void QUADTREE<T>::
Project(ARRAY<T>& x)
{
    if(uniform) Project_Uniform(x);
    else Project_Non_Uniform(x);
}
template<class T> void QUADTREE<T>::
Project_Uniform(ARRAY<T>& x)
{
    for(int j1=1;j1<=x.Size();++j1)
    {
        int index1=Array_Index_2_Tree_Index(j1);
        if(nodes(index1)->levelset>0 && nodes(index1)->Neumann_of_cell==false) 
        {
            x(j1)=0.;
        }
    }
}
template<class T> void QUADTREE<T>::
Project_Non_Uniform(ARRAY<T>& x)
{
    for(int j1=1;j1<=x.Size();++j1)
    {
        int index1=Array_Index_2_Tree_Index(j1);
        if(nodes(index1)->levelset>0 && nodes(index1)->Neumann_of_cell==false) 
            x(j1)=0.;
    }
}
//#####################################################################
// Verify_First_of_Two_Finer_Neighbors
//#####################################################################
template<class T> bool QUADTREE<T>::
Verify_First_of_Two_Finer_Neighbors(int coarser,int finer)
{
    for(int i=1;i<=4;++i)
    {
        int index1,index2;
        index1=neighbor_indices(coarser)(i).x;
        index2=neighbor_indices(coarser)(i).y;
        if(index1==0) continue;
        if(index1==finer) return true;
        if(index2==finer) return false;
    }
    LOG::cout << "Verify_First_of_Two_Finer_Neighbors failed!" << std::endl;
    exit(0);
    return false;
}
//#####################################################################
// Phi_Circle
//#####################################################################
template<class T> T QUADTREE<T>::
Phi_Circle(TV position)
{
    // >0, air
    // <=0, dof
    return 0.25111111111-(position-TV(0.5,0.5)).Magnitude();
}
//#####################################################################
// Compute_Volume_And_Area_For_Pressure
//#####################################################################
template<class T> void QUADTREE<T>::
Compute_Volume_And_Area_For_Pressure()
{

#if 0
    // if neighbor is higher, lose volume; else gain volume
    for(int i=1;i<=nodes.Size();++i)
    {
        QUADTREE_NODE<T>* node=nodes(i);
        if(node->north_east_index==0 && !ghost(i))
        {
            bool flag=false;
            node->volume_pressure=(node->dx_over_two*2).Product();
            for(int j=1;j<=4;++j)
            {
                int index1,index2;
                index1=neighbor_indices(i)(j).x;
                index2=neighbor_indices(i)(j).y;
                if(index1!=0){
                    T dx=node->dx_over_two(1); // edge length of smaller cells
                    node->normal_of_face.Append(-(node->center-nodes(index1)->center).Normalized());
                    node->center_of_face.Append(node->center+node->normal_of_face(node->normal_of_face.Size())*4*dx/sqrt(10));
                    node->normal_of_face.Append(-(node->center-nodes(index2)->center).Normalized());
                    node->center_of_face.Append(node->center+node->normal_of_face(node->normal_of_face.Size())*4*dx/sqrt(10));
                }
                else{
                    if(nodes(index2)->dx_over_two!=node->dx_over_two){
                        T dx=node->dx_over_two(1)*2;  // edge length of smaller cells
                        node->normal_of_face.Append(-(node->center-nodes(index2)->center).Normalized());
                        node->center_of_face.Append(nodes(index2)->center-node->normal_of_face(node->normal_of_face.Size())*4*dx/sqrt(10));
                    }
                    else{
                        node->normal_of_face.Append(-(node->center-nodes(index2)->center).Normalized());
                        node->center_of_face.Append(node->center+node->normal_of_face(node->normal_of_face.Size())*node->dx_over_two(1));
                    }
                }
            }

            for(int j=1;j<=4;++j)
            {
                int index1,index2;
                index1=neighbor_indices(i)(j).x;
                index2=neighbor_indices(i)(j).y;
                if(index1!=0){
                    T dx=node->dx_over_two(1); // edge length of smaller cells
                    T dv=.5*dx*dx/3;
                    node->volume_pressure+=2*dv;
                    T area=dx/3*sqrt(10);
                    node->area_of_face.Append(area);
                    node->area_of_face.Append(area);
                }
                else{
                    if(nodes(index2)->dx_over_two!=node->dx_over_two){
                        T dx=node->dx_over_two(1)*2;  // edge length of smaller cells
                        T dv=.5*dx*dx/3;
                        node->volume_pressure-=dv;
                        T area1=dx/3*sqrt(10);
                        T area2=dx/3*2;
                        node->area_of_face.Append(area1);
                        if(!Verify_First_of_Two_Finer_Neighbors(index2,i)){
                            if(j!=4){
                                node->area_of_face.Append(area2);++j;}
                            else
                                node->area_of_face(1)=area2;
                        }
                        else{
                            if(j>1) 
                            {
                                node->area_of_face(node->area_of_face.Size()-1)=area2;
                                if(j==4) flag=false; // more than 1 coarser neighbors 
                            }
                            else flag=true;
                        }
                    }
                    else{
                        T area=node->dx_over_two(1)*2;
                        node->area_of_face.Append(area);
                    }
                }
            }
            if(flag==true) node->area_of_face(4)=node->dx_over_two(1)*2/3*2;
        }
    }
#else
    // if neighbor is higher, lose volume; else gain volume
    ARRAY<TV> diagonal_vector,axis_vector;
    diagonal_vector.Append(TV(1.,-1.));
    diagonal_vector.Append(TV(1.,1.));
    diagonal_vector.Append(TV(-1.,1.));
    diagonal_vector.Append(TV(-1.,-1.));
    axis_vector.Append(TV(1.,0.));
    axis_vector.Append(TV(0.,1.));
    axis_vector.Append(TV(-1.,0.));
    axis_vector.Append(TV(0.,-1.));
    for(int i=1;i<=nodes.Size();++i)
    {
        QUADTREE_NODE<T>* node=nodes(i);
        if(node->north_east_index==0 && !ghost(i))
        {
            for(int j=1;j<=4;++j)
            {
                int k1=j,k2=j==4?1:j+1;

                int index1,index2;
                index1=neighbor_indices(i)(j).x;
                index2=neighbor_indices(i)(j).y;
                if(index1!=0){
                    T dx=node->dx_over_two(1); // edge length of smaller cells
                    node->normal_of_face.Append(-(node->center-nodes(index1)->center).Normalized());
                    node->center_of_face.Append(node->center+node->normal_of_face(node->normal_of_face.Size())*4*dx/sqrt(10));
                    node->normal_of_face.Append(-(node->center-nodes(index2)->center).Normalized());
                    node->center_of_face.Append(node->center+node->normal_of_face(node->normal_of_face.Size())*4*dx/sqrt(10));
                    // nodes 
                }
                else{
                    if(nodes(index2)->dx_over_two!=node->dx_over_two){
                        T dx=node->dx_over_two(1)*2;  // edge length of smaller cells
                        node->normal_of_face.Append(-(node->center-nodes(index2)->center).Normalized());
                        node->center_of_face.Append(nodes(index2)->center-node->normal_of_face(node->normal_of_face.Size())*4*dx/sqrt(10));

                    }
                    else{
                        T dx=node->dx_over_two(1);  
                        node->normal_of_face.Append(-(node->center-nodes(index2)->center).Normalized());
                        node->center_of_face.Append(node->center+node->normal_of_face(node->normal_of_face.Size())*node->dx_over_two(1));
                    }
                }
            }
            bool flag=false;
            for(int j=1;j<=4;++j)
            {
                int k1=j,k2=j==4?1:j+1;

                int index1,index2;
                index1=neighbor_indices(i)(j).x;
                index2=neighbor_indices(i)(j).y;
                if(index1!=0){
                    T dx=node->dx_over_two(1); // edge length of smaller cells
                    // nodes 
                    node->position_of_node.Append(node->center+diagonal_vector(j)*dx);
                    node->position_of_node.Append(node->center+axis_vector(j)*dx*4/3);
                }
                else{
                    if(nodes(index2)->dx_over_two!=node->dx_over_two){
                        T dx=node->dx_over_two(1)*2;  // edge length of smaller cells

                        if(Verify_First_of_Two_Finer_Neighbors(index2,i)){
                            //if(i==172) LOG::cout<<"got here 172" << std::endl;
                            node->position_of_node.Append(node->center+diagonal_vector(j)*dx*.5-axis_vector(j)*dx/3.);
                        }
                        else{
                            //if(i==173) LOG::cout<<"got here 173" << std::endl;
                            if(j<4){
                                node->position_of_node.Append(node->center+diagonal_vector(j)*dx*.5);
                                node->position_of_node.Append(node->center+diagonal_vector(j+1)*dx*.5-axis_vector(j)*dx/3.);
                                j++;
                                flag=false;
                            }
                            else{
                                node->position_of_node.Append(node->center+diagonal_vector(j)*dx*.5);
                                flag=true;
                            }

                        }
                    }
                    else{
                        T dx=node->dx_over_two(1);  
                        node->position_of_node.Append(node->center+diagonal_vector(j)*dx);
                    }
                }
            }
            if(flag) node->position_of_node(1)=node->center+diagonal_vector(1)*node->dx_over_two(1)-axis_vector(4)*node->dx_over_two(1)*2./3.;

            flag=false;
            for(int j=1;j<=4;++j)
            {
                int index1,index2;
                index1=neighbor_indices(i)(j).x;
                index2=neighbor_indices(i)(j).y;
                if(index1!=0){
                    T dx=node->dx_over_two(1); // edge length of smaller cells
                    T dv=.5*dx*dx/3;
                    //node->volume_pressure+=2*dv;
                    T area=dx/3*sqrt(10);
                    node->area_of_face.Append(area);
                    node->area_of_face.Append(area);
                }
                else{
                    if(nodes(index2)->dx_over_two!=node->dx_over_two){
                        T dx=node->dx_over_two(1)*2;  // edge length of smaller cells
                        T dv=.5*dx*dx/3;
                        //node->volume_pressure-=dv;
                        T area1=dx/3*sqrt(10);
                        T area2=dx/3*2;
                        node->area_of_face.Append(area1);
                        if(!Verify_First_of_Two_Finer_Neighbors(index2,i)){
                            if(j!=4){
                                node->area_of_face.Append(area2);++j;}
                            else
                                node->area_of_face(1)=area2;
                        }
                        else{
                            if(j>1) 
                            {
                                node->area_of_face(node->area_of_face.Size()-1)=area2;
                                if(j==4) flag=false; // more than 1 coarser neighbors 
                            }
                            else flag=true;
                        }
                    }
                    else{
                        T area=node->dx_over_two(1)*2;
                        node->area_of_face.Append(area);
                    }
                }
            }
            if(flag==true) node->area_of_face(4)=node->dx_over_two(1)*2/3*2;

            node->Neumann_of_cell=false;
            ARRAY<TV> point;
            for(int j=1;j<=node->position_of_node.Size();++j) // here j iterates over edges
            {
                node->dof_theta_of_edge.Append(0.);
                int k1=j,k2=j==node->position_of_node.Size()?1:j+1;
                T phi1=Phi_Circle(node->position_of_node(k1)), phi2=Phi_Circle(node->position_of_node(k2));
                // draw a table for This
                // >0 =0 <0
                if(phi1<0 && phi2>0){
                    node->Neumann_of_cell=true;
                    node->dof_theta_of_edge(j)=LEVELSET_UTILITIES<T>::Theta(phi1,phi2);
                    node->position_of_Neumann_node.Append_Unique(node->position_of_node(k1));
                    node->position_of_Neumann_node.Append_Unique(node->position_of_node(k1)+(node->position_of_node(k2)-node->position_of_node(k1))*node->dof_theta_of_edge(j));
                    point.Append(node->position_of_node(k1)+(node->position_of_node(k2)-node->position_of_node(k1))*node->dof_theta_of_edge(j));
                    node->air_node=node->position_of_node(k2);
                }
                else if(phi1>0 && phi2<0)
                {
                    node->Neumann_of_cell=true;
                    node->dof_theta_of_edge(j)=-LEVELSET_UTILITIES<T>::Theta(phi2,phi1);
                    node->position_of_Neumann_node.Append_Unique(node->position_of_node(k2)-(node->position_of_node(k1)-node->position_of_node(k2))*node->dof_theta_of_edge(j));
                    node->position_of_Neumann_node.Append_Unique(node->position_of_node(k2));
                    point.Append(node->position_of_node(k2)-(node->position_of_node(k1)-node->position_of_node(k2))*node->dof_theta_of_edge(j));
                    node->air_node=node->position_of_node(k1);
                }
                else if(phi1<0 && phi2<0)
                {
                    node->dof_theta_of_edge(j)=1.;
                    node->position_of_Neumann_node.Append_Unique(node->position_of_node(k1));
                    node->position_of_Neumann_node.Append_Unique(node->position_of_node(k2));
                }
                else if(phi1>0 && phi2>0)
                {
                    node->dof_theta_of_edge(j)=0.;
                }
            }
                //LOG::cout << i << "," << point.Size() <<std::endl;
                //for(int j=1;j<=node->position_of_node.Size();++j) LOG::cout << "position " <<j<<"," <<node->position_of_node(j) << std::endl;
                //for(int j=1;j<=node->position_of_node.Size();++j) LOG::cout << "theta " <<j<<"," <<node->dof_theta_of_edge(j) << std::endl;
                //for(int j=1;j<=node->position_of_node.Size();++j){
                    //int k1=j,k2=j==node->position_of_node.Size()?1:j+1;
                    //LOG::cout << i<< " diff " <<j<<"," <<node->area_of_face(j)-(node->position_of_node(k1)-node->position_of_node(k2)).Magnitude() << std::endl;
                    //LOG::cout << i<< " Magnitude " <<j<<"," <<(node->position_of_node(k1)-node->position_of_node(k2)).Magnitude() << std::endl;
                    //LOG::cout << i<< " area_of_face " <<j<<"," <<node->area_of_face(j) << std::endl;
                //}
            if(node->Neumann_of_cell)
            {
                assert(point.Size()==2);
                SEGMENT_2D<T> segment(point(1),point(2));
                node->new_point=(point(1)+point(2))*.5;
                node->new_normal=segment.Normal();
                if(TV::Dot_Product((node->air_node-node->new_point),node->new_normal)<0)
                    node->new_normal*=-1.;

                node->new_length=segment.Size();

                POLYGON<TV> polygon(node->position_of_Neumann_node.Size());
                polygon.X=node->position_of_Neumann_node;
                //LOG::cout << i << "," << node->position_of_Neumann_node.Size() << std::endl;
                node->volume_pressure=polygon.Area();
                //LOG::cout << i << "," << node->position_of_Neumann_node.Size() << "," << node->volume_pressure << std::endl;
            }
            else{
                POLYGON<TV> polygon(node->position_of_node.Size());
                polygon.X=node->position_of_node;
                node->volume_pressure=polygon.Area();
            }

                //LOG::cout << "point " << node->new_point << std::endl;
                //LOG::cout << "normal " << node->new_normal << std::endl;
                //LOG::cout << "length " << node->new_length << std::endl;
                //LOG::cout << "volume before " << node->volume_pressure << std::endl;

        }
    }
#endif

#if 0
    // decide each cell to be Neumann boundary or not
    // uniform only now
    // need to update
    for(int i=1;i<=nodes.Size();++i)
    {
        QUADTREE_NODE<T>* node=nodes(i);
        if(node->north_east_index==0 && !ghost(i))
        {
            T x=node->dx_over_two.x;
            T y=node->dx_over_two.y;
            node->position_of_node.Append(node->center+TV(x,-y));
            node->position_of_node.Append(node->center+TV(x,y));
            node->position_of_node.Append(node->center+TV(-x,y));
            node->position_of_node.Append(node->center+TV(-x,-y));
        
            node->Neumann_of_cell=false;
            for(int j=1;j<=4;++j) // here j iterates over edges
            {
                node->dof_theta_of_edge.Append(0.);
                int k1=j,k2=j==4?1:j+1;
                T phi1=Phi_Circle(node->position_of_node(k1)), phi2=Phi_Circle(node->position_of_node(k2));
                // draw a table for This
                // >0 =0 <0
                if(phi1*phi2<0){
                    node->Neumann_of_cell=true;
                    if(phi1<0)
                        node->dof_theta_of_edge(j)=LEVELSET_UTILITIES<T>::Theta(phi1,phi2);
                    else
                        node->dof_theta_of_edge(j)=-LEVELSET_UTILITIES<T>::Theta(phi2,phi1);
                }
                else if((phi1<0 && phi2<=0) || (phi1<=0 && phi2<0))
                    node->dof_theta_of_edge(j)=1.;
            }
        }
    }
    for(int i=1;i<=nodes.Size();++i)
    {
        QUADTREE_NODE<T>* node=nodes(i);
        if(node->north_east_index==0 && !ghost(i))
        {
            if(node->Neumann_of_cell) 
            {
                ARRAY<TV> point;
                for(int j=1;j<=4;++j) // here j iterates over edges
                {
                    int k1=j;int k2=j==4?1:j+1;
                    if(node->dof_theta_of_edge(j)!=0 && node->dof_theta_of_edge(j)!=1)
                    {
                        TV tmp;
                        if(node->dof_theta_of_edge(j)>0)
                            tmp=node->position_of_node(k1)+(node->position_of_node(k2)-node->position_of_node(k1))*node->dof_theta_of_edge(j);
                        else
                            tmp=node->position_of_node(k1)+(node->position_of_node(k2)-node->position_of_node(k1))*(1.+node->dof_theta_of_edge(j));
                        point.Append(tmp);
                    }
                    if(node->dof_theta_of_edge(j)==1)
                    {
                        T phi1=Phi_Circle(node->position_of_node(k1)), phi2=Phi_Circle(node->position_of_node(k2));
                        if(phi1==0){LOG::cout<<"AAAAAAAA"<<std::endl; point.Append(node->position_of_node(k1));}
                        if(phi2==0){LOG::cout<<"AAAAAAAA"<<std::endl; point.Append(node->position_of_node(k2));}
                    }
                }
                assert(point.Size()==2);
                SEGMENT_2D<T> segment(point(1),point(2));
                node->new_point=(point(1)+point(2))*.5;
                node->new_normal=segment.Normal();
                node->new_length=segment.Size();

                int number_of_air_node=0;
                int test_node_air=0;
                int test_node_dof=0;
                for(int j=1;j<=4;++j)
                {
                    if(Phi_Circle(node->position_of_node(j))>0){
                        number_of_air_node++;
                        test_node_air=j;
                    }
                    else test_node_dof=j;
                }
                assert(number_of_air_node!=0);

                if(TV::Dot_Product((node->position_of_node(test_node_air)-node->new_point),node->new_normal)<0)
                    node->new_normal*=-1.;
                
                LOG::cout << i << " --- number_of_air_node " << number_of_air_node << std::endl;
                LOG::cout << "theta1 " << node->dof_theta_of_edge(1) << std::endl;
                LOG::cout << "theta2 " << node->dof_theta_of_edge(2) << std::endl;
                LOG::cout << "theta3 " << node->dof_theta_of_edge(3) << std::endl;
                LOG::cout << "theta4 " << node->dof_theta_of_edge(4) << std::endl;
                LOG::cout << "point " << node->new_point << std::endl;
                LOG::cout << "normal " << node->new_normal << std::endl;
                LOG::cout << "length " << node->new_length << std::endl;
                LOG::cout << "volume before " << node->volume_pressure << std::endl;
                // fix volume
                if(number_of_air_node==1)
                {
                    node->volume_pressure-=TRIANGLE_2D<T>::Area(point(1),point(2),node->position_of_node(test_node_air));
                }
                else if(number_of_air_node==2)
                {
                    T tmp=0.;
                    for(int j=1;j<=4;++j)
                    {
                        if(node->dof_theta_of_edge(j)!=1) tmp+=std::fabs(node->dof_theta_of_edge(j))/2.;
                    }
                    node->volume_pressure*=tmp;
                }
                else if(number_of_air_node==3)
                {
                    node->volume_pressure=TRIANGLE_2D<T>::Area(point(1),point(2),node->position_of_node(test_node_dof));
                }
                else{LOG::cout << "bad happend" << std::endl;exit(0);}
                LOG::cout << "volume after " << node->volume_pressure << std::endl;
                
                // use polygon to compute volume
                ARRAY<TV> polygon_nodes;
                for(int j=1;j<=4;++j)
                {
                    int k1=j,k2=j==4?1:j+1;
                    T phi1=Phi_Circle(node->position_of_node(k1)), phi2=Phi_Circle(node->position_of_node(k2));
                    if(phi1*phi2<0)
                    {
                            TV tmp;
                            if(node->dof_theta_of_edge(j)>0)
                                tmp=node->position_of_node(k1)+(node->position_of_node(k2)-node->position_of_node(k1))*node->dof_theta_of_edge(j);
                            else
                                tmp=node->position_of_node(k1)+(node->position_of_node(k2)-node->position_of_node(k1))*(1.+node->dof_theta_of_edge(j));
                        if(phi1>0){
                            polygon_nodes.Append(tmp);
                            polygon_nodes.Append(node->position_of_node(k2));
                        }else{
                            polygon_nodes.Append(node->position_of_node(k1));
                            polygon_nodes.Append(tmp);
                        }
                    }
                    if(phi1<0 && phi2<0)
                    {
                            polygon_nodes.Append(node->position_of_node(k1));
                            polygon_nodes.Append(node->position_of_node(k2));
                    }
                }
                POLYGON<TV> polygon(polygon_nodes.Size());
                polygon.X=polygon_nodes;
                LOG::cout << "volume after2 " << polygon.Area() << std::endl;
                LOG::cout << "volume diff " << polygon.Area()-node->volume_pressure << std::endl;
            }
        }
    }

                //LOG::cout << "theta1 " << nodes(63)->dof_theta_of_edge(1) << std::endl;
                //LOG::cout << "theta2 " << nodes(63)->dof_theta_of_edge(2) << std::endl;
                //LOG::cout << "theta3 " << nodes(63)->dof_theta_of_edge(3) << std::endl;
                //LOG::cout << "theta4 " << nodes(63)->dof_theta_of_edge(4) << std::endl;
#endif
}
//#####################################################################
template class QUADTREE_NODE<float>;
template class QUADTREE<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class QUADTREE_NODE<double>;
template class QUADTREE<double>;
#endif
