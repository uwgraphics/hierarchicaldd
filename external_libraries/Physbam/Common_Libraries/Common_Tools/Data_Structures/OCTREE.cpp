//#####################################################################
// Copyright 2015, Mridul Aanjaneya, Ming Gao.
// This file is part of PhysBAM whose distribution is governed by the license 
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Common_Tools/Data_Structures/OCTREE.h>
#include <PhysBAM_Tools/Math_Tools/cube.h>
#include <PhysBAM_Tools/Data_Structures/PAIR.h>
#include <PhysBAM_Tools/Math_Tools/exchange.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Geometry/Level_Sets/LEVELSET_UTILITIES.h>
#include <queue>
#include <vector>
//#include <voro++/voro++.hh>
#include <PhysBAM_Geometry/Basic_Geometry/POLYGON.h>
#include <PhysBAM_Geometry/Basic_Geometry/TRIANGLE_3D.h>
#include <PhysBAM_Geometry/Basic_Geometry/PLANE.h>
#include <Common_Geometry/Read_Write/Topology_Based_Geometry/READ_WRITE_VORONOI_DIAGRAM.h>

//#include <Eigen/Core>
//#include <Eigen/Dense>
//#include <Eigen/Sparse>
//#include <Eigen/Eigenvalues>

using namespace PhysBAM;
//using namespace voro;
//#####################################################################
// Destructor
//#####################################################################
/*template<class T> OCTREE<T>::*/
//~OCTREE()
//{
    //for(int i=1;i<=nodes.Size();++i) if(nodes(i)) delete nodes(i);
/*}*/
//#####################################################################
// Compute_Two_Level_Subdivision
//#####################################################################
template<class T> void OCTREE<T>::
Compute_A_Simple_Subdivision(const RANGE<TV>& domain_input,bool fine)
{
    domain=domain_input;
    TV dx=(T).5*domain.Edge_Lengths();
    //LOG::cout << dx << std::endl;exit(0);
    OCTREE_NODE<T> *node=new OCTREE_NODE<T>(domain.Center(),dx);
    nodes.Append(node);
    node->index=1;
    node->parent_index=0;

    int scale=4;
    std::queue<OCTREE_NODE<T>*> list;list.push(node);
    int resolution=1;
    while(resolution<scale/2){resolution*=2;TV current_dx=dx/resolution;
        int counter=0;
        // subdivide the quadtree
        while(counter<cube(resolution/2)){++counter;
            OCTREE_NODE<T> *top=list.front();
            list.pop();
            // create 8 new corners
            ARRAY<TV> new_corners(8);
            new_corners(1)=top->center+TV(current_dx(1),-current_dx(2),current_dx(3));
            new_corners(2)=top->center+TV(current_dx(1),current_dx(2),current_dx(3));
            new_corners(3)=top->center+TV(-current_dx(1),current_dx(2),current_dx(3));
            new_corners(4)=top->center+TV(-current_dx(1),-current_dx(2),current_dx(3));
            new_corners(5)=top->center+TV(current_dx(1),-current_dx(2),-current_dx(3));
            new_corners(6)=top->center+TV(current_dx(1),current_dx(2),-current_dx(3));
            new_corners(7)=top->center+TV(-current_dx(1),current_dx(2),-current_dx(3));
            new_corners(8)=top->center+TV(-current_dx(1),-current_dx(2),-current_dx(3));
            // create 8 new octree nodes
            ARRAY<OCTREE_NODE<T>*> new_nodes(8);
            for(int i=1;i<=8;++i){
                new_nodes(i)=new OCTREE_NODE<T>(new_corners(i),current_dx);
                new_nodes(i)->parent_index=top->index;
                new_nodes(i)->n_of_8=i;
            }
            // reassign pointers
            for(int i=1;i<=8;++i) top->children(i)=nodes.Size()+i;
            // append in nodes
            for(int i=1;i<=8;++i){
                nodes.Append(new_nodes(i));
                new_nodes(i)->index=nodes.Size();}
            // push in queue
            for(int i=1;i<=8;++i) list.push(new_nodes(i));}}

        while(!list.empty()){OCTREE_NODE<T> *top=list.front();
            //LOG::cout << top->index << std::endl;
            TV current_dx=top->dx_over_two/2;
            list.pop();
            bool kkkk;
            if(!fine) kkkk=top->index==4 || top->index==3 || top->index==8 || top->index==7 || top->index==2 || top->index==6 || top->index==9 
                || top->center.y==0.875
                ;
            else kkkk=top->index==4 || top->index==3 || top->index==8 || top->index==7 || top->index==2 || top->index==6 || top->index==9 
                || top->center.y==0.875
                || top->center==TV(0.625,0.625,0.375)
                ;
            if(kkkk){
                // create 8 new corners
                ARRAY<TV> new_corners(8);
                // create 8 new octree nodes
                ARRAY<OCTREE_NODE<T>*> new_nodes(8);
                new_corners(1)=top->center+TV(current_dx(1),-current_dx(2),current_dx(3));
                new_corners(2)=top->center+TV(current_dx(1),current_dx(2),current_dx(3));
                new_corners(3)=top->center+TV(-current_dx(1),current_dx(2),current_dx(3));
                new_corners(4)=top->center+TV(-current_dx(1),-current_dx(2),current_dx(3));
                new_corners(5)=top->center+TV(current_dx(1),-current_dx(2),-current_dx(3));
                new_corners(6)=top->center+TV(current_dx(1),current_dx(2),-current_dx(3));
                new_corners(7)=top->center+TV(-current_dx(1),current_dx(2),-current_dx(3));
                new_corners(8)=top->center+TV(-current_dx(1),-current_dx(2),-current_dx(3));
                for(int i=1;i<=8;++i){
                    new_nodes(i)=new OCTREE_NODE<T>(new_corners(i),current_dx);
                    new_nodes(i)->parent_index=top->index;
                    new_nodes(i)->n_of_8=i;
                }
                // reassign pointers
                for(int i=1;i<=8;++i) top->children(i)=nodes.Size()+i;
                // append in nodes
                for(int i=1;i<=8;++i){
                    nodes.Append(new_nodes(i));
                    new_nodes(i)->index=nodes.Size();}
                // push in queue
                for(int i=1;i<=8;++i) list.push(new_nodes(i));}}
        

        //LOG::cout << "number of nodes " << nodes.Size() << std::endl;

}
template<class T> void OCTREE<T>::
Compute_Two_Level_Subdivision(const RANGE<TV>& domain_input,const int scale,bool uniform_input)
{
    uniform=uniform_input;
    domain=domain_input;
    int resolution=1;TV dx=(T).5*domain.Edge_Lengths();
    OCTREE_NODE<T> *node=new OCTREE_NODE<T>(domain.Center(),dx);
    nodes.Append(node);node->index=1;node->parent_index=0;

    std::queue<OCTREE_NODE<T>*> list;list.push(node);
    while(resolution<scale/2){resolution*=2;TV current_dx=dx/resolution;
        int counter=0;
        // subdivide the quadtree
        while(counter<cube(resolution/2)){++counter;
            OCTREE_NODE<T> *top=list.front();
            list.pop();
            // create 8 new corners
            ARRAY<TV> new_corners(8);
            new_corners(1)=top->center+TV(current_dx(1),-current_dx(2),current_dx(3));
            new_corners(2)=top->center+TV(current_dx(1),current_dx(2),current_dx(3));
            new_corners(3)=top->center+TV(-current_dx(1),current_dx(2),current_dx(3));
            new_corners(4)=top->center+TV(-current_dx(1),-current_dx(2),current_dx(3));
            new_corners(5)=top->center+TV(current_dx(1),-current_dx(2),-current_dx(3));
            new_corners(6)=top->center+TV(current_dx(1),current_dx(2),-current_dx(3));
            new_corners(7)=top->center+TV(-current_dx(1),current_dx(2),-current_dx(3));
            new_corners(8)=top->center+TV(-current_dx(1),-current_dx(2),-current_dx(3));
            // create 8 new octree nodes
            ARRAY<OCTREE_NODE<T>*> new_nodes(8);
            for(int i=1;i<=8;++i){
                new_nodes(i)=new OCTREE_NODE<T>(new_corners(i),current_dx);
                new_nodes(i)->parent_index=top->index;
                new_nodes(i)->n_of_8=i;
            }
            // reassign pointers
            for(int i=1;i<=8;++i) top->children(i)=nodes.Size()+i;
            // append in nodes
            for(int i=1;i<=8;++i){
                nodes.Append(new_nodes(i));
                new_nodes(i)->index=nodes.Size();}
            // push in queue
            for(int i=1;i<=8;++i) list.push(new_nodes(i));}}

    // specific to uniform boundary cells
    T minx=1.,maxx=0.,miny=1.,maxy=0.,minz=1.,maxz=0.;
    for(int i=1;i<=nodes.Size();++i)
        if(nodes(i)->children(1)==0){
            if(nodes(i)->center.x>maxx) maxx=nodes(i)->center.x;
            if(nodes(i)->center.x<minx) minx=nodes(i)->center.x;
            if(nodes(i)->center.y>maxy) maxy=nodes(i)->center.y;
            if(nodes(i)->center.y<miny) miny=nodes(i)->center.y;
            if(nodes(i)->center.z>maxz) maxz=nodes(i)->center.z;
            if(nodes(i)->center.z<minz) minz=nodes(i)->center.z;
        }
    //LOG::cout << minx << "," << maxx << "," << miny << "," << maxy << std::endl;
    ghost.Exact_Resize(nodes.Size());
    for(int i=1;i<=nodes.Size();++i)
    {
        ghost(i)=false;
        if(nodes(i)->children(1)==0)
            if(nodes(i)->center.x==maxx || nodes(i)->center.x==minx || nodes(i)->center.y==maxy || nodes(i)->center.y==miny || nodes(i)->center.z==maxz || nodes(i)->center.z==minz) ghost(i)=true;
    }


    // list only contains cells of given resolution; subdivide the left side
    if(!uniform){
        while(!list.empty()){OCTREE_NODE<T> *top=list.front();
            TV current_dx=top->dx_over_two/2;
            list.pop();
            if(!ghost(top->index) && !(top->center.x < node->center.x && top->center.y < node->center.y)){
            //if(!ghost(top->index) && top->center.x < node->center.x){
            //if(!ghost(top->index) && (
                        //(std::fabs((top->center-TV(.5,.5,.5)).Magnitude())>0.15111111111 && std::fabs((top->center-TV(.5,.5,.5)).Magnitude())<0.25111111111) ||
                        //(std::fabs((top->center-TV(.5,.5,.5)).Magnitude())>0.35111111111)
                        //)){
                // create 8 new corners
                ARRAY<TV> new_corners(8);
                // create 8 new octree nodes
                ARRAY<OCTREE_NODE<T>*> new_nodes(8);
                new_corners(1)=top->center+TV(current_dx(1),-current_dx(2),current_dx(3));
                new_corners(2)=top->center+TV(current_dx(1),current_dx(2),current_dx(3));
                new_corners(3)=top->center+TV(-current_dx(1),current_dx(2),current_dx(3));
                new_corners(4)=top->center+TV(-current_dx(1),-current_dx(2),current_dx(3));
                new_corners(5)=top->center+TV(current_dx(1),-current_dx(2),-current_dx(3));
                new_corners(6)=top->center+TV(current_dx(1),current_dx(2),-current_dx(3));
                new_corners(7)=top->center+TV(-current_dx(1),current_dx(2),-current_dx(3));
                new_corners(8)=top->center+TV(-current_dx(1),-current_dx(2),-current_dx(3));
                for(int i=1;i<=8;++i){
                    new_nodes(i)=new OCTREE_NODE<T>(new_corners(i),current_dx);
                    new_nodes(i)->parent_index=top->index;
                    new_nodes(i)->n_of_8=i;
                }
                // reassign pointers
                for(int i=1;i<=8;++i) top->children(i)=nodes.Size()+i;
                // append in nodes
                for(int i=1;i<=8;++i){
                    nodes.Append(new_nodes(i));
                    new_nodes(i)->index=nodes.Size();}
                }}}

    int ghost_size=ghost.Size();
    for(int i=1;i<=nodes.Size()-ghost_size;++i)
        ghost.Append(false);

    number_of_leaves_including_ghosts=0;
    for(int i=1;i<=nodes.Size();++i)
    {
        if(nodes(i)->children(1)==0){
            number_of_leaves_including_ghosts++;
            nodes(i)->voronoi_index=number_of_leaves_including_ghosts;
            voronoi_index_2_tree_index.Append(nodes(i)->index);
        }
    }

}
//#####################################################################
// Compute_Neighbors
//#####################################################################
template<class T> void OCTREE<T>::
Compute_Neighbors()
{
   for(int i=1;i<=nodes.Size();++i)
    {
        ARRAY<VECTOR<int,4> > indices(6);
        for(int j=1;j<=6;++j) indices(j)=Neighbor_Leaf(i,j);
        neighbor_indices.Append(indices);
    }
}
//#####################################################################
// Corner_Neighbor_Leaf_Fine
//#####################################################################
template<class T> void OCTREE<T>::
Corner_Neighbor_Leaf_Fine(const int index,const int n_of_6,ARRAY<int>& corner_indices,ARRAY<int>& corner_codes) // can be 1 or 2 corner neighbors
{
    corner_indices.Resize(0);
    // -1 coarse
    // -2 fine
    corner_codes.Resize(0);
    int index2=Neighbor_Same_Level(index,n_of_6);

    if(n_of_6==1)
    {
        int u=Neighbor_Same_Level(index2,3);
        if(u==0 || nodes(u)->children(1)==0){
            corner_indices.Append(u);
            corner_codes.Append(-1);
        }
        else{
            corner_indices.Append(nodes(u)->children(4));
            corner_indices.Append(nodes(u)->children(8));
            corner_codes.Append(-2);
        }

        u=Neighbor_Same_Level(index2,4);
        if(u==0 || nodes(u)->children(1)==0){
            corner_indices.Append(u);
            corner_codes.Append(-1);
        }
        else{
            corner_indices.Append(nodes(u)->children(3));
            corner_indices.Append(nodes(u)->children(7));
            corner_codes.Append(-2);
        }

        u=Neighbor_Same_Level(index2,5);
        if(u==0 || nodes(u)->children(1)==0){
            corner_indices.Append(u);
            corner_codes.Append(-1);
        }
        else{
            corner_indices.Append(nodes(u)->children(7));
            corner_indices.Append(nodes(u)->children(8));
            corner_codes.Append(-2);
        }

        u=Neighbor_Same_Level(index2,6);
        if(u==0 || nodes(u)->children(1)==0){
            corner_indices.Append(u);
            corner_codes.Append(-1);
        }
        else{
            corner_indices.Append(nodes(u)->children(3));
            corner_indices.Append(nodes(u)->children(4));
            corner_codes.Append(-2);
        }

        u=Neighbor_Same_Level(Neighbor_Same_Level(index2,3),5);
        if(u==0 || nodes(u)->children(1)==0){
            corner_indices.Append(u);
            corner_codes.Append(-1);
        }
        else{
            corner_indices.Append(nodes(u)->children(8));
            corner_codes.Append(-2);
        }

        u=Neighbor_Same_Level(Neighbor_Same_Level(index2,5),4);
        if(u==0 || nodes(u)->children(1)==0){
            corner_indices.Append(u);
            corner_codes.Append(-1);
        }
        else{
            corner_indices.Append(nodes(u)->children(7));
            corner_codes.Append(-2);
        }

        u=Neighbor_Same_Level(Neighbor_Same_Level(index2,4),6);
        if(u==0 || nodes(u)->children(1)==0){
            corner_indices.Append(u);
            corner_codes.Append(-1);
        }
        else{
            corner_indices.Append(nodes(u)->children(3));
            corner_codes.Append(-2);
        }

        u=Neighbor_Same_Level(Neighbor_Same_Level(index2,6),3);
        if(u==0 || nodes(u)->children(1)==0){
            corner_indices.Append(u);
            corner_codes.Append(-1);
        }
        else{
            corner_indices.Append(nodes(u)->children(4));
            corner_codes.Append(-2);
        }
    }

    if(n_of_6==2)
    {
        int u=Neighbor_Same_Level(index2,3);
        if(u==0 || nodes(u)->children(1)==0){
            corner_indices.Append(u);
            corner_codes.Append(-1);
        }
        else{
            corner_indices.Append(nodes(u)->children(1));
            corner_indices.Append(nodes(u)->children(5));
            corner_codes.Append(-2);
        }

        u=Neighbor_Same_Level(index2,4);
        if(u==0 || nodes(u)->children(1)==0){
            corner_indices.Append(u);
            corner_codes.Append(-1);
        }
        else{
            corner_indices.Append(nodes(u)->children(2));
            corner_indices.Append(nodes(u)->children(6));
            corner_codes.Append(-2);
        }

        u=Neighbor_Same_Level(index2,5);
        if(u==0 || nodes(u)->children(1)==0){
            corner_indices.Append(u);
            corner_codes.Append(-1);
        }
        else{
            corner_indices.Append(nodes(u)->children(5));
            corner_indices.Append(nodes(u)->children(6));
            corner_codes.Append(-2);
        }

        u=Neighbor_Same_Level(index2,6);
        if(u==0 || nodes(u)->children(1)==0){
            corner_indices.Append(u);
            corner_codes.Append(-1);
        }
        else{
            corner_indices.Append(nodes(u)->children(1));
            corner_indices.Append(nodes(u)->children(2));
            corner_codes.Append(-2);
        }

        u=Neighbor_Same_Level(Neighbor_Same_Level(index2,3),5);
        if(u==0 || nodes(u)->children(1)==0){
            corner_indices.Append(u);
            corner_codes.Append(-1);
        }
        else{
            corner_indices.Append(nodes(u)->children(5));
            corner_codes.Append(-2);
        }

        u=Neighbor_Same_Level(Neighbor_Same_Level(index2,5),4);
        if(u==0 || nodes(u)->children(1)==0){
            corner_indices.Append(u);
            corner_codes.Append(-1);
        }
        else{
            corner_indices.Append(nodes(u)->children(6));
            corner_codes.Append(-2);
        }

        u=Neighbor_Same_Level(Neighbor_Same_Level(index2,4),6);
        if(u==0 || nodes(u)->children(1)==0){
            corner_indices.Append(u);
            corner_codes.Append(-1);
        }
        else{
            corner_indices.Append(nodes(u)->children(2));
            corner_codes.Append(-2);
        }

        u=Neighbor_Same_Level(Neighbor_Same_Level(index2,6),3);
        if(u==0 || nodes(u)->children(1)==0){
            corner_indices.Append(u);
            corner_codes.Append(-1);
        }
        else{
            corner_indices.Append(nodes(u)->children(1));
            corner_codes.Append(-2);
        }
    }

    if(n_of_6==3)
    {
        int u=Neighbor_Same_Level(index2,1);
        if(u==0 || nodes(u)->children(1)==0){
            corner_indices.Append(u);
            corner_codes.Append(-1);
        }
        else{
            corner_indices.Append(nodes(u)->children(4));
            corner_indices.Append(nodes(u)->children(8));
            corner_codes.Append(-2);
        }

        u=Neighbor_Same_Level(index2,2);
        if(u==0 || nodes(u)->children(1)==0){
            corner_indices.Append(u);
            corner_codes.Append(-1);
        }
        else{
            corner_indices.Append(nodes(u)->children(1));
            corner_indices.Append(nodes(u)->children(5));
            corner_codes.Append(-2);
        }

        u=Neighbor_Same_Level(index2,5);
        if(u==0 || nodes(u)->children(1)==0){
            corner_indices.Append(u);
            corner_codes.Append(-1);
        }
        else{
            corner_indices.Append(nodes(u)->children(5));
            corner_indices.Append(nodes(u)->children(8));
            corner_codes.Append(-2);
        }

        u=Neighbor_Same_Level(index2,6);
        if(u==0 || nodes(u)->children(1)==0){
            corner_indices.Append(u);
            corner_codes.Append(-1);
        }
        else{
            corner_indices.Append(nodes(u)->children(1));
            corner_indices.Append(nodes(u)->children(4));
            corner_codes.Append(-2);
        }

        u=Neighbor_Same_Level(Neighbor_Same_Level(index2,1),5);
        if(u==0 || nodes(u)->children(1)==0){
            corner_indices.Append(u);
            corner_codes.Append(-1);
        }
        else{
            corner_indices.Append(nodes(u)->children(8));
            corner_codes.Append(-2);
        }

        u=Neighbor_Same_Level(Neighbor_Same_Level(index2,5),2);
        if(u==0 || nodes(u)->children(1)==0){
            corner_indices.Append(u);
            corner_codes.Append(-1);
        }
        else{
            corner_indices.Append(nodes(u)->children(5));
            corner_codes.Append(-2);
        }

        u=Neighbor_Same_Level(Neighbor_Same_Level(index2,2),6);
        if(u==0 || nodes(u)->children(1)==0){
            corner_indices.Append(u);
            corner_codes.Append(-1);
        }
        else{
            corner_indices.Append(nodes(u)->children(1));
            corner_codes.Append(-2);
        }

        u=Neighbor_Same_Level(Neighbor_Same_Level(index2,6),1);
        if(u==0 || nodes(u)->children(1)==0){
            corner_indices.Append(u);
            corner_codes.Append(-1);
        }
        else{
            corner_indices.Append(nodes(u)->children(4));
            corner_codes.Append(-2);
        }
    }


    if(n_of_6==4)
    {
        int u=Neighbor_Same_Level(index2,1);
        if(u==0 || nodes(u)->children(1)==0){
            corner_indices.Append(u);
            corner_codes.Append(-1);
        }
        else{
            corner_indices.Append(nodes(u)->children(3));
            corner_indices.Append(nodes(u)->children(7));
            corner_codes.Append(-2);
        }

        u=Neighbor_Same_Level(index2,2);
        if(u==0 || nodes(u)->children(1)==0){
            corner_indices.Append(u);
            corner_codes.Append(-1);
        }
        else{
            corner_indices.Append(nodes(u)->children(2));
            corner_indices.Append(nodes(u)->children(6));
            corner_codes.Append(-2);
        }

        u=Neighbor_Same_Level(index2,5);
        if(u==0 || nodes(u)->children(1)==0){
            corner_indices.Append(u);
            corner_codes.Append(-1);
        }
        else{
            corner_indices.Append(nodes(u)->children(6));
            corner_indices.Append(nodes(u)->children(7));
            corner_codes.Append(-2);
        }

        u=Neighbor_Same_Level(index2,6);
        if(u==0 || nodes(u)->children(1)==0){
            corner_indices.Append(u);
            corner_codes.Append(-1);
        }
        else{
            corner_indices.Append(nodes(u)->children(2));
            corner_indices.Append(nodes(u)->children(3));
            corner_codes.Append(-2);
        }

        u=Neighbor_Same_Level(Neighbor_Same_Level(index2,1),5);
        if(u==0 || nodes(u)->children(1)==0){
            corner_indices.Append(u);
            corner_codes.Append(-1);
        }
        else{
            corner_indices.Append(nodes(u)->children(7));
            corner_codes.Append(-2);
        }

        u=Neighbor_Same_Level(Neighbor_Same_Level(index2,5),2);
        if(u==0 || nodes(u)->children(1)==0){
            corner_indices.Append(u);
            corner_codes.Append(-1);
        }
        else{
            corner_indices.Append(nodes(u)->children(6));
            corner_codes.Append(-2);
        }

        u=Neighbor_Same_Level(Neighbor_Same_Level(index2,2),6);
        if(u==0 || nodes(u)->children(1)==0){
            corner_indices.Append(u);
            corner_codes.Append(-1);
        }
        else{
            corner_indices.Append(nodes(u)->children(2));
            corner_codes.Append(-2);
        }

        u=Neighbor_Same_Level(Neighbor_Same_Level(index2,6),1);
        if(u==0 || nodes(u)->children(1)==0){
            corner_indices.Append(u);
            corner_codes.Append(-1);
        }
        else{
            corner_indices.Append(nodes(u)->children(3));
            corner_codes.Append(-2);
        }
    }

    if(n_of_6==5)
    {
        int u=Neighbor_Same_Level(index2,1);
        if(u==0 || nodes(u)->children(1)==0){
            corner_indices.Append(u);
            corner_codes.Append(-1);
        }
        else{
            corner_indices.Append(nodes(u)->children(7));
            corner_indices.Append(nodes(u)->children(8));
            corner_codes.Append(-2);
        }

        u=Neighbor_Same_Level(index2,2);
        if(u==0 || nodes(u)->children(1)==0){
            corner_indices.Append(u);
            corner_codes.Append(-1);
        }
        else{
            corner_indices.Append(nodes(u)->children(5));
            corner_indices.Append(nodes(u)->children(6));
            corner_codes.Append(-2);
        }

        u=Neighbor_Same_Level(index2,3);
        if(u==0 || nodes(u)->children(1)==0){
            corner_indices.Append(u);
            corner_codes.Append(-1);
        }
        else{
            corner_indices.Append(nodes(u)->children(5));
            corner_indices.Append(nodes(u)->children(8));
            corner_codes.Append(-2);
        }

        u=Neighbor_Same_Level(index2,4);
        if(u==0 || nodes(u)->children(1)==0){
            corner_indices.Append(u);
            corner_codes.Append(-1);
        }
        else{
            corner_indices.Append(nodes(u)->children(6));
            corner_indices.Append(nodes(u)->children(7));
            corner_codes.Append(-2);
        }

        u=Neighbor_Same_Level(Neighbor_Same_Level(index2,1),3);
        if(u==0 || nodes(u)->children(1)==0){
            corner_indices.Append(u);
            corner_codes.Append(-1);
        }
        else{
            corner_indices.Append(nodes(u)->children(8));
            corner_codes.Append(-2);
        }

        u=Neighbor_Same_Level(Neighbor_Same_Level(index2,3),2);
        if(u==0 || nodes(u)->children(1)==0){
            corner_indices.Append(u);
            corner_codes.Append(-1);
        }
        else{
            corner_indices.Append(nodes(u)->children(5));
            corner_codes.Append(-2);
        }

        u=Neighbor_Same_Level(Neighbor_Same_Level(index2,2),4);
        if(u==0 || nodes(u)->children(1)==0){
            corner_indices.Append(u);
            corner_codes.Append(-1);
        }
        else{
            corner_indices.Append(nodes(u)->children(6));
            corner_codes.Append(-2);
        }

        u=Neighbor_Same_Level(Neighbor_Same_Level(index2,4),1);
        if(u==0 || nodes(u)->children(1)==0){
            corner_indices.Append(u);
            corner_codes.Append(-1);
        }
        else{
            corner_indices.Append(nodes(u)->children(7));
            corner_codes.Append(-2);
        }
    }

    if(n_of_6==6)
    {
        int u=Neighbor_Same_Level(index2,1);
        if(u==0 || nodes(u)->children(1)==0){
            corner_indices.Append(u);
            corner_codes.Append(-1);
        }
        else{
            corner_indices.Append(nodes(u)->children(3));
            corner_indices.Append(nodes(u)->children(4));
            corner_codes.Append(-2);
        }

        u=Neighbor_Same_Level(index2,2);
        if(u==0 || nodes(u)->children(1)==0){
            corner_indices.Append(u);
            corner_codes.Append(-1);
        }
        else{
            corner_indices.Append(nodes(u)->children(1));
            corner_indices.Append(nodes(u)->children(2));
            corner_codes.Append(-2);
        }

        u=Neighbor_Same_Level(index2,3);
        if(u==0 || nodes(u)->children(1)==0){
            corner_indices.Append(u);
            corner_codes.Append(-1);
        }
        else{
            corner_indices.Append(nodes(u)->children(1));
            corner_indices.Append(nodes(u)->children(4));
            corner_codes.Append(-2);
        }

        u=Neighbor_Same_Level(index2,4);
        if(u==0 || nodes(u)->children(1)==0){
            corner_indices.Append(u);
            corner_codes.Append(-1);
        }
        else{
            corner_indices.Append(nodes(u)->children(2));
            corner_indices.Append(nodes(u)->children(3));
            corner_codes.Append(-2);
        }

        u=Neighbor_Same_Level(Neighbor_Same_Level(index2,1),3);
        if(u==0 || nodes(u)->children(1)==0){
            corner_indices.Append(u);
            corner_codes.Append(-1);
        }
        else{
            corner_indices.Append(nodes(u)->children(8));
            corner_codes.Append(-2);
        }

        u=Neighbor_Same_Level(Neighbor_Same_Level(index2,3),2);
        if(u==0 || nodes(u)->children(1)==0){
            corner_indices.Append(u);
            corner_codes.Append(-1);
        }
        else{
            corner_indices.Append(nodes(u)->children(5));
            corner_codes.Append(-2);
        }

        u=Neighbor_Same_Level(Neighbor_Same_Level(index2,2),4);
        if(u==0 || nodes(u)->children(1)==0){
            corner_indices.Append(u);
            corner_codes.Append(-1);
        }
        else{
            corner_indices.Append(nodes(u)->children(6));
            corner_codes.Append(-2);
        }

        u=Neighbor_Same_Level(Neighbor_Same_Level(index2,4),1);
        if(u==0 || nodes(u)->children(1)==0){
            corner_indices.Append(u);
            corner_codes.Append(-1);
        }
        else{
            corner_indices.Append(nodes(u)->children(7));
            corner_codes.Append(-2);
        }
    }











/*


    if(n_of_6==2)
    {
        int u=Neighbor_Same_Level(index2,3);
        if(u==0 || nodes(u)->children(1)==0) corner_indices.Append(u);
        else{
            corner_indices.Append(nodes(u)->children(1));
            corner_indices.Append(nodes(u)->children(5));}
        u=Neighbor_Same_Level(index2,4);
        if(u==0 || nodes(u)->children(1)==0) corner_indices.Append(u);
        else{
            corner_indices.Append(nodes(u)->children(2));
            corner_indices.Append(nodes(u)->children(6));}
        u=Neighbor_Same_Level(index2,5);
        if(u==0 || nodes(u)->children(1)==0) corner_indices.Append(u);
        else{
            corner_indices.Append(nodes(u)->children(5));
            corner_indices.Append(nodes(u)->children(6));}
        u=Neighbor_Same_Level(index2,6);
        if(u==0 || nodes(u)->children(1)==0) corner_indices.Append(u);
        else{
            corner_indices.Append(nodes(u)->children(1));
            corner_indices.Append(nodes(u)->children(2));}
    }
    if(n_of_6==3)
    {
        int u=Neighbor_Same_Level(index2,1);
        if(u==0 || nodes(u)->children(1)==0) corner_indices.Append(u);
        else{
            corner_indices.Append(nodes(u)->children(4));
            corner_indices.Append(nodes(u)->children(8));}
        u=Neighbor_Same_Level(index2,2);
        if(u==0 || nodes(u)->children(1)==0) corner_indices.Append(u);
        else{
            corner_indices.Append(nodes(u)->children(1));
            corner_indices.Append(nodes(u)->children(5));}
        u=Neighbor_Same_Level(index2,5);
        if(u==0 || nodes(u)->children(1)==0) corner_indices.Append(u);
        else{
            corner_indices.Append(nodes(u)->children(5));
            corner_indices.Append(nodes(u)->children(8));}
        u=Neighbor_Same_Level(index2,6);
        if(u==0 || nodes(u)->children(1)==0) corner_indices.Append(u);
        else{
            corner_indices.Append(nodes(u)->children(1));
            corner_indices.Append(nodes(u)->children(4));}
    }
    if(n_of_6==4)
    {
        int u=Neighbor_Same_Level(index2,1);
        if(u==0 || nodes(u)->children(1)==0) corner_indices.Append(u);
        else{
            corner_indices.Append(nodes(u)->children(3));
            corner_indices.Append(nodes(u)->children(7));}
        u=Neighbor_Same_Level(index2,2);
        if(u==0 || nodes(u)->children(1)==0) corner_indices.Append(u);
        else{
            corner_indices.Append(nodes(u)->children(2));
            corner_indices.Append(nodes(u)->children(6));}
        u=Neighbor_Same_Level(index2,5);
        if(u==0 || nodes(u)->children(1)==0) corner_indices.Append(u);
        else{
            corner_indices.Append(nodes(u)->children(6));
            corner_indices.Append(nodes(u)->children(7));}
        u=Neighbor_Same_Level(index2,6);
        if(u==0 || nodes(u)->children(1)==0) corner_indices.Append(u);
        else{
            corner_indices.Append(nodes(u)->children(2));
            corner_indices.Append(nodes(u)->children(3));}
    }
    if(n_of_6==5)
    {
        int u=Neighbor_Same_Level(index2,1);
        if(u==0 || nodes(u)->children(1)==0) corner_indices.Append(u);
        else{
            corner_indices.Append(nodes(u)->children(7));
            corner_indices.Append(nodes(u)->children(8));}
        u=Neighbor_Same_Level(index2,2);
        if(u==0 || nodes(u)->children(1)==0) corner_indices.Append(u);
        else{
            corner_indices.Append(nodes(u)->children(5));
            corner_indices.Append(nodes(u)->children(6));}
        u=Neighbor_Same_Level(index2,3);
        if(u==0 || nodes(u)->children(1)==0) corner_indices.Append(u);
        else{
            corner_indices.Append(nodes(u)->children(5));
            corner_indices.Append(nodes(u)->children(8));}
        u=Neighbor_Same_Level(index2,4);
        if(u==0 || nodes(u)->children(1)==0) corner_indices.Append(u);
        else{
            corner_indices.Append(nodes(u)->children(6));
            corner_indices.Append(nodes(u)->children(7));}
    }
    if(n_of_6==6)
    {
        int u=Neighbor_Same_Level(index2,1);
        if(u==0 || nodes(u)->children(1)==0) corner_indices.Append(u);
        else{
            corner_indices.Append(nodes(u)->children(3));
            corner_indices.Append(nodes(u)->children(4));}
        u=Neighbor_Same_Level(index2,2);
        if(u==0 || nodes(u)->children(1)==0) corner_indices.Append(u);
        else{
            corner_indices.Append(nodes(u)->children(1));
            corner_indices.Append(nodes(u)->children(2));}
        u=Neighbor_Same_Level(index2,3);
        if(u==0 || nodes(u)->children(1)==0) corner_indices.Append(u);
        else{
            corner_indices.Append(nodes(u)->children(1));
            corner_indices.Append(nodes(u)->children(4));}
        u=Neighbor_Same_Level(index2,4);
        if(u==0 || nodes(u)->children(1)==0) corner_indices.Append(u);
        else{
            corner_indices.Append(nodes(u)->children(2));
            corner_indices.Append(nodes(u)->children(3));}
    }
    */
}
//#####################################################################
// All_26_Neighbors
//#####################################################################
template<class T> void OCTREE<T>::
All_26_Neighbors(const int index,ARRAY<int>& all_neighbor_indices,ARRAY<int>& all_neighbor_codes,bool exclude_vertex_neighbor)
{
    // code 1 : 1 level higher resolution (split)
    // code 2 : same resolution 
    // code 3 : 1 level lower resolution 
    // 6 face neighbors 
    // 12 edge neighbors 
    // 8 vertex neighbors 
    
    T dx=nodes(index)->dx_over_two(1);
    int n_of_8=nodes(index)->n_of_8;

    all_neighbor_indices.Exact_Resize(0);all_neighbor_indices.Append(index);
    all_neighbor_codes.Exact_Resize(0);all_neighbor_codes.Append(n_of_8);
    // 6 face neighbors 
    { 
        int directions[]={2,1,4,3,6,5};
        for(int i=1;i<=6;++i)
        {
            int index2=S_FACE_NEIGHBOR(index,directions[i-1]);
            if(index2==0){all_neighbor_codes.Append(0);continue;}
            if(nodes(index2)->children(1)==0)
            {
                all_neighbor_indices.Append_Unique(index2);
                if(nodes(index2)->dx_over_two(1)==dx)
                    all_neighbor_codes.Append(2);
                else 
                    all_neighbor_codes.Append(3);
            }
            else
            {
                all_neighbor_codes.Append(1);
                for(int j=1;j<=8;++j) assert(nodes(nodes(index2)->children(j))->children(1)==0);
                    //all_neighbor_indices.Append_Unique(nodes(index2)->children(j));
                switch(i){
                    case 2:all_neighbor_indices.Append_Unique(nodes(index2)->children(3));all_neighbor_indices.Append_Unique(nodes(index2)->children(4));
                    all_neighbor_indices.Append_Unique(nodes(index2)->children(7));all_neighbor_indices.Append_Unique(nodes(index2)->children(8));break;

                    case 1:all_neighbor_indices.Append_Unique(nodes(index2)->children(1));all_neighbor_indices.Append_Unique(nodes(index2)->children(2));
                    all_neighbor_indices.Append_Unique(nodes(index2)->children(5));all_neighbor_indices.Append_Unique(nodes(index2)->children(6));break;

                    case 4:all_neighbor_indices.Append_Unique(nodes(index2)->children(1));all_neighbor_indices.Append_Unique(nodes(index2)->children(4));
                    all_neighbor_indices.Append_Unique(nodes(index2)->children(5));all_neighbor_indices.Append_Unique(nodes(index2)->children(8));break;

                    case 3:all_neighbor_indices.Append_Unique(nodes(index2)->children(3));all_neighbor_indices.Append_Unique(nodes(index2)->children(2));
                    all_neighbor_indices.Append_Unique(nodes(index2)->children(7));all_neighbor_indices.Append_Unique(nodes(index2)->children(6));break;

                    case 6:all_neighbor_indices.Append_Unique(nodes(index2)->children(5));all_neighbor_indices.Append_Unique(nodes(index2)->children(6));
                    all_neighbor_indices.Append_Unique(nodes(index2)->children(7));all_neighbor_indices.Append_Unique(nodes(index2)->children(8));break;

                    case 5:default:all_neighbor_indices.Append_Unique(nodes(index2)->children(1));all_neighbor_indices.Append_Unique(nodes(index2)->children(2));
                    all_neighbor_indices.Append_Unique(nodes(index2)->children(3));all_neighbor_indices.Append_Unique(nodes(index2)->children(4));break;
                }
            }
        }
    }
    // 12 edge neighbors 
    { 
        int directions[]={24,23,26,25,14,13,16,15,46,45,36,35};
        for(int i=1;i<=12;++i)
        {
            int index2=S_EDGE_NEIGHBOR(index,directions[i-1]);
            if(index2==0){all_neighbor_codes.Append(0);continue;}
            if(nodes(index2)->children(1)==0)
            {
                all_neighbor_indices.Append_Unique(index2);
                if(nodes(index2)->dx_over_two(1)==dx)
                    all_neighbor_codes.Append(2);
                else 
                    all_neighbor_codes.Append(3);
            }
            else
            {
                all_neighbor_codes.Append(1);
                for(int j=1;j<=8;++j)
                {
                    assert(nodes(nodes(index2)->children(j))->children(1)==0);
                    //all_neighbor_indices.Append_Unique(nodes(index2)->children(j));
                }
                switch(i){
                    case 1:all_neighbor_indices.Append_Unique(nodes(index2)->children(2));all_neighbor_indices.Append_Unique(nodes(index2)->children(6));break;
                    case 2:all_neighbor_indices.Append_Unique(nodes(index2)->children(1));all_neighbor_indices.Append_Unique(nodes(index2)->children(5));break;
                    case 3:all_neighbor_indices.Append_Unique(nodes(index2)->children(1));all_neighbor_indices.Append_Unique(nodes(index2)->children(2));break;
                    case 4:all_neighbor_indices.Append_Unique(nodes(index2)->children(5));all_neighbor_indices.Append_Unique(nodes(index2)->children(6));break;

                    case 5:all_neighbor_indices.Append_Unique(nodes(index2)->children(3));all_neighbor_indices.Append_Unique(nodes(index2)->children(7));break;
                    case 6:all_neighbor_indices.Append_Unique(nodes(index2)->children(4));all_neighbor_indices.Append_Unique(nodes(index2)->children(8));break;
                    case 7:all_neighbor_indices.Append_Unique(nodes(index2)->children(3));all_neighbor_indices.Append_Unique(nodes(index2)->children(4));break;
                    case 8:all_neighbor_indices.Append_Unique(nodes(index2)->children(7));all_neighbor_indices.Append_Unique(nodes(index2)->children(8));break;

                    case 9:all_neighbor_indices.Append_Unique(nodes(index2)->children(2));all_neighbor_indices.Append_Unique(nodes(index2)->children(3));break;
                    case 10:all_neighbor_indices.Append_Unique(nodes(index2)->children(7));all_neighbor_indices.Append_Unique(nodes(index2)->children(6));break;
                    case 11:all_neighbor_indices.Append_Unique(nodes(index2)->children(1));all_neighbor_indices.Append_Unique(nodes(index2)->children(4));break;
                    case 12:default:all_neighbor_indices.Append_Unique(nodes(index2)->children(5));all_neighbor_indices.Append_Unique(nodes(index2)->children(8));break;
                }
            }
        }
    }
    // 8 vertex neighbors
    if(exclude_vertex_neighbor==false)
    { 
        int directions[]={246,245,236,235,146,145,136,135};
        for(int i=1;i<=8;++i)
        {
            int index2=S_VERTEX_NEIGHBOR(index,directions[i-1]);
            if(index2==0){all_neighbor_codes.Append(0);continue;}
            if(nodes(index2)->children(1)==0)
            {
                //all_neighbor_indices.Append_Unique(index2);
                if(nodes(index2)->dx_over_two(1)==dx)
                    all_neighbor_codes.Append(2);
                else 
                    all_neighbor_codes.Append(3);
            }
            else
            {
                all_neighbor_codes.Append(1);
                for(int j=1;j<=8;++j)
                {
                    all_neighbor_indices.Append_Unique(nodes(index2)->children(j));
                    assert(nodes(nodes(index2)->children(j))->children(1)==0);
                }
            }
        }
    }
}
//#####################################################################
// Check_Neighbor_Type 
// 1. face_neighbor
// 2. edge_neighbor
//#####################################################################
template<class T> int OCTREE<T>::
Check_Neighbor_Type(const int index,const int index_neighbor,int& direction)
{
    T dx=nodes(index)->dx_over_two(1);
    int n_of_8=nodes(index)->n_of_8;

    // face
    {
        int directions[]={2,1,4,3,6,5};
        for(int i=1;i<=6;++i)
        {
            direction=directions[i-1];
            int index2=S_FACE_NEIGHBOR(index,directions[i-1]);
            if(index2==0)continue;
            if(nodes(index2)->children(1)==0){
                if(index2==index_neighbor) return 1;}
            else{
                switch(i){
                    case 2: if(index_neighbor==nodes(index2)->children(3) || index_neighbor==nodes(index2)->children(4) ||
                                index_neighbor==nodes(index2)->children(7) || index_neighbor==nodes(index2)->children(8)) return 1;
                    case 1: if(index_neighbor==nodes(index2)->children(1) || index_neighbor==nodes(index2)->children(2) ||
                                index_neighbor==nodes(index2)->children(5) || index_neighbor==nodes(index2)->children(6)) return 1;
                    case 4: if(index_neighbor==nodes(index2)->children(1) || index_neighbor==nodes(index2)->children(4) ||
                                index_neighbor==nodes(index2)->children(5) || index_neighbor==nodes(index2)->children(8)) return 1;
                    case 3: if(index_neighbor==nodes(index2)->children(3) || index_neighbor==nodes(index2)->children(2) ||
                                index_neighbor==nodes(index2)->children(7) || index_neighbor==nodes(index2)->children(6)) return 1;
                    case 6: if(index_neighbor==nodes(index2)->children(5) || index_neighbor==nodes(index2)->children(6) ||
                                index_neighbor==nodes(index2)->children(7) || index_neighbor==nodes(index2)->children(8)) return 1;
                    case 5:default: if(index_neighbor==nodes(index2)->children(1) || index_neighbor==nodes(index2)->children(2) ||
                                index_neighbor==nodes(index2)->children(3) || index_neighbor==nodes(index2)->children(4)) return 1;
                }
            }
        }
    }
    // edge
    { 
        int directions[]={24,23,26,25,14,13,16,15,46,45,36,35};
        for(int i=1;i<=12;++i)
        {
            direction=directions[i-1];
            int index2=S_EDGE_NEIGHBOR(index,directions[i-1]);
            if(index2==0)continue;
            if(nodes(index2)->children(1)==0){
                if(index2==index_neighbor) return 2;}
            else
            {
                switch(i){
                    case 1: if(index_neighbor==nodes(index2)->children(2) || index_neighbor==nodes(index2)->children(6)) return 2; 
                    case 2: if(index_neighbor==nodes(index2)->children(1) || index_neighbor==nodes(index2)->children(5)) return 2; 
                    case 3: if(index_neighbor==nodes(index2)->children(1) || index_neighbor==nodes(index2)->children(2)) return 2; 
                    case 4: if(index_neighbor==nodes(index2)->children(5) || index_neighbor==nodes(index2)->children(6)) return 2; 
                    case 5: if(index_neighbor==nodes(index2)->children(3) || index_neighbor==nodes(index2)->children(7)) return 2; 
                    case 6: if(index_neighbor==nodes(index2)->children(4) || index_neighbor==nodes(index2)->children(8)) return 2; 
                    case 7: if(index_neighbor==nodes(index2)->children(3) || index_neighbor==nodes(index2)->children(4)) return 2; 
                    case 8: if(index_neighbor==nodes(index2)->children(7) || index_neighbor==nodes(index2)->children(8)) return 2; 
                    case 9: if(index_neighbor==nodes(index2)->children(2) || index_neighbor==nodes(index2)->children(3)) return 2; 
                    case 10: if(index_neighbor==nodes(index2)->children(7) || index_neighbor==nodes(index2)->children(6)) return 2; 
                    case 11: if(index_neighbor==nodes(index2)->children(1) || index_neighbor==nodes(index2)->children(4)) return 2; 
                    case 12:default: if(index_neighbor==nodes(index2)->children(5) || index_neighbor==nodes(index2)->children(8)) return 2; 
                }
            }
        }
    }
    return 0;
}
//#####################################################################
// Corner_Neighbor_Leaf_Coarse
//#####################################################################
template<class T> int OCTREE<T>::
Corner_Neighbor_Leaf_Coarse(int index,VECTOR<int,2> plane,bool& coarse)
{
    int parent_index=nodes(index)->parent_index;
    int index2;
    coarse=false;
    switch(nodes(index)->n_of_8)
    {
        case 1:
            if(plane==VECTOR<int,2>(1,3)) // xz plane
            {
                index2=Neighbor_Same_Level(Neighbor_Same_Level(parent_index,1),5);
                if(nodes(index2)->children(1)==0) return index2;
                else{coarse=true; return nodes(index2)->children(8);}
            }
            else if(plane==VECTOR<int,2>(1,2)) // xy plane
            {
                index2=Neighbor_Same_Level(Neighbor_Same_Level(parent_index,1),4);
                if(nodes(index2)->children(1)==0) return index2;
                else{coarse=true; return nodes(index2)->children(3);}
            }
            else // yz plane
            {
                index2=Neighbor_Same_Level(Neighbor_Same_Level(parent_index,4),5);
                if(nodes(index2)->children(1)==0) return index2;
                else{coarse=true; return nodes(index2)->children(6);}
            }
            break;
        case 2:
            if(plane==VECTOR<int,2>(1,3)) // xz plane
            {
                index2=Neighbor_Same_Level(Neighbor_Same_Level(parent_index,1),5);
                if(nodes(index2)->children(1)==0) return index2;
                else{coarse=true; return nodes(index2)->children(7);}
            }
            else if(plane==VECTOR<int,2>(1,2)) // xy plane
            {
                index2=Neighbor_Same_Level(Neighbor_Same_Level(parent_index,1),3);
                if(nodes(index2)->children(1)==0) return index2;
                else{coarse=true; return nodes(index2)->children(4);}
            }
            else // yz plane
            {
                index2=Neighbor_Same_Level(Neighbor_Same_Level(parent_index,3),5);
                if(nodes(index2)->children(1)==0) return index2;
                else{coarse=true; return nodes(index2)->children(5);}
            }
            break;
        case 3:
            if(plane==VECTOR<int,2>(1,3)) // xz plane
            {
                index2=Neighbor_Same_Level(Neighbor_Same_Level(parent_index,2),5);
                if(nodes(index2)->children(1)==0) return index2;
                else{coarse=true; return nodes(index2)->children(6);}
            }
            else if(plane==VECTOR<int,2>(1,2)) // xy plane
            {
                index2=Neighbor_Same_Level(Neighbor_Same_Level(parent_index,2),3);
                if(nodes(index2)->children(1)==0) return index2;
                else{coarse=true; return nodes(index2)->children(1);}
            }
            else // yz plane
            {
                index2=Neighbor_Same_Level(Neighbor_Same_Level(parent_index,3),5);
                if(nodes(index2)->children(1)==0) return index2;
                else{coarse=true; return nodes(index2)->children(8);}
            }
            break;
        case 4:
            if(plane==VECTOR<int,2>(1,3)) // xz plane
            {
                index2=Neighbor_Same_Level(Neighbor_Same_Level(parent_index,2),5);
                if(nodes(index2)->children(1)==0) return index2;
                else{coarse=true; return nodes(index2)->children(5);}
            }
            else if(plane==VECTOR<int,2>(1,2)) // xy plane
            {
                index2=Neighbor_Same_Level(Neighbor_Same_Level(parent_index,2),4);
                if(nodes(index2)->children(1)==0) return index2;
                else{coarse=true; return nodes(index2)->children(2);}
            }
            else // yz plane
            {
                index2=Neighbor_Same_Level(Neighbor_Same_Level(parent_index,4),5);
                if(nodes(index2)->children(1)==0) return index2;
                else{coarse=true; return nodes(index2)->children(7);}
            }
            break;
        case 5:
            if(plane==VECTOR<int,2>(1,3)) // xz plane
            {
                index2=Neighbor_Same_Level(Neighbor_Same_Level(parent_index,1),6);
                if(nodes(index2)->children(1)==0) return index2;
                else{coarse=true; return nodes(index2)->children(4);}
            }
            else if(plane==VECTOR<int,2>(1,2)) // xy plane
            {
                index2=Neighbor_Same_Level(Neighbor_Same_Level(parent_index,1),4);
                if(nodes(index2)->children(1)==0) return index2;
                else{coarse=true; return nodes(index2)->children(7);}
            }
            else // yz plane
            {
                index2=Neighbor_Same_Level(Neighbor_Same_Level(parent_index,4),6);
                if(nodes(index2)->children(1)==0) return index2;
                else{coarse=true; return nodes(index2)->children(2);}
            }
            break;
        case 6:
            if(plane==VECTOR<int,2>(1,3)) // xz plane
            {
                index2=Neighbor_Same_Level(Neighbor_Same_Level(parent_index,1),6);
                if(nodes(index2)->children(1)==0) return index2;
                else{coarse=true; return nodes(index2)->children(3);}
            }
            else if(plane==VECTOR<int,2>(1,2)) // xy plane
            {
                index2=Neighbor_Same_Level(Neighbor_Same_Level(parent_index,1),3);
                if(nodes(index2)->children(1)==0) return index2;
                else{coarse=true; return nodes(index2)->children(8);}
            }
            else // yz plane
            {
                index2=Neighbor_Same_Level(Neighbor_Same_Level(parent_index,3),6);
                if(nodes(index2)->children(1)==0) return index2;
                else{coarse=true; return nodes(index2)->children(1);}
            }
            break;
        case 7:
            if(plane==VECTOR<int,2>(1,3)) // xz plane
            {
                index2=Neighbor_Same_Level(Neighbor_Same_Level(parent_index,2),6);
                if(nodes(index2)->children(1)==0) return index2;
                else{coarse=true; return nodes(index2)->children(2);}
            }
            else if(plane==VECTOR<int,2>(1,2)) // xy plane
            {
                index2=Neighbor_Same_Level(Neighbor_Same_Level(parent_index,2),3);
                if(nodes(index2)->children(1)==0) return index2;
                else{coarse=true; return nodes(index2)->children(5);}
            }
            else // yz plane
            {
                index2=Neighbor_Same_Level(Neighbor_Same_Level(parent_index,3),6);
                if(nodes(index2)->children(1)==0) return index2;
                else{coarse=true; return nodes(index2)->children(4);}
            }
            break;
        case 8:default:
            if(plane==VECTOR<int,2>(1,3)) // xz plane
            {
                index2=Neighbor_Same_Level(Neighbor_Same_Level(parent_index,2),6);
                if(nodes(index2)->children(1)==0) return index2;
                else{coarse=true; return nodes(index2)->children(1);}
            }
            else if(plane==VECTOR<int,2>(1,2)) // xy plane
            {
                index2=Neighbor_Same_Level(Neighbor_Same_Level(parent_index,2),4);
                if(nodes(index2)->children(1)==0) return index2;
                else{coarse=true; return nodes(index2)->children(6);}
            }
            else // yz plane
            {
                index2=Neighbor_Same_Level(Neighbor_Same_Level(parent_index,4),6);
                if(nodes(index2)->children(1)==0) return index2;
                else{coarse=true; return nodes(index2)->children(3);}
            }
    }
}
//#####################################################################
// Tree_Index_2_Array_Index
// Array_Index_2_Tree_Index
//#####################################################################
template<class T> int OCTREE<T>::
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
template<class T> int OCTREE<T>::
Array_Index_2_Tree_Index(int array_index)
{
    if(array_index_2_tree_index.Size()==0)
    {
        for(int i=1;i<=nodes.Size();++i)
            if(nodes(i)->children(1)==0 && !ghost(i))
                array_index_2_tree_index.Append(nodes(i)->index);
    }
    return array_index_2_tree_index(array_index);
}
//#####################################################################
// Neighbor_Leaf
//#####################################################################
template<class T> VECTOR<int,4> OCTREE<T>::
Neighbor_Leaf(const int index, const int n_of_6)
{
    if(n_of_6==1)
    {
        if(nodes(index)->children(1)!=0) return VECTOR<int,4>(-1,-1,-1,-1);
        int neighbor=Neighbor_Same_Level(index,n_of_6);
        if(neighbor==0) return VECTOR<int,4>(0,0,0,0);
        if(nodes(neighbor)->children(1)==0) return VECTOR<int,4>(0,0,0,neighbor);
        return VECTOR<int,4>(nodes(neighbor)->children(4),nodes(neighbor)->children(3),nodes(neighbor)->children(7),nodes(neighbor)->children(8));
    }
    if(n_of_6==2)
    {
        if(nodes(index)->children(1)!=0) return VECTOR<int,4>(-1,-1,-1,-1);
        int neighbor=Neighbor_Same_Level(index,n_of_6);
        if(neighbor==0) return VECTOR<int,4>(0,0,0,0);
        if(nodes(neighbor)->children(1)==0) return VECTOR<int,4>(0,0,0,neighbor);
        return VECTOR<int,4>(nodes(neighbor)->children(5),nodes(neighbor)->children(6),nodes(neighbor)->children(2),nodes(neighbor)->children(1));
    }
    if(n_of_6==3)
    {
        if(nodes(index)->children(1)!=0) return VECTOR<int,4>(-1,-1,-1,-1);
        int neighbor=Neighbor_Same_Level(index,n_of_6);
        if(neighbor==0) return VECTOR<int,4>(0,0,0,0);
        if(nodes(neighbor)->children(1)==0) return VECTOR<int,4>(0,0,0,neighbor);
        return VECTOR<int,4>(nodes(neighbor)->children(5),nodes(neighbor)->children(1),nodes(neighbor)->children(4),nodes(neighbor)->children(8));
    }
    if(n_of_6==4)
    {
        if(nodes(index)->children(1)!=0) return VECTOR<int,4>(-1,-1,-1,-1);
        int neighbor=Neighbor_Same_Level(index,n_of_6);
        if(neighbor==0) return VECTOR<int,4>(0,0,0,0);
        if(nodes(neighbor)->children(1)==0) return VECTOR<int,4>(0,0,0,neighbor);
        return VECTOR<int,4>(nodes(neighbor)->children(2),nodes(neighbor)->children(6),nodes(neighbor)->children(7),nodes(neighbor)->children(3));
    }
    if(n_of_6==5)
    {
        if(nodes(index)->children(1)!=0) return VECTOR<int,4>(-1,-1,-1,-1);
        int neighbor=Neighbor_Same_Level(index,n_of_6);
        if(neighbor==0) return VECTOR<int,4>(0,0,0,0);
        if(nodes(neighbor)->children(1)==0) return VECTOR<int,4>(0,0,0,neighbor);
        return VECTOR<int,4>(nodes(neighbor)->children(8),nodes(neighbor)->children(7),nodes(neighbor)->children(6),nodes(neighbor)->children(5));
    }
    if(n_of_6==6)
    {
        if(nodes(index)->children(1)!=0) return VECTOR<int,4>(-1,-1,-1,-1);
        int neighbor=Neighbor_Same_Level(index,n_of_6);
        if(neighbor==0) return VECTOR<int,4>(0,0,0,0);
        if(nodes(neighbor)->children(1)==0) return VECTOR<int,4>(0,0,0,neighbor);
        return VECTOR<int,4>(nodes(neighbor)->children(1),nodes(neighbor)->children(2),nodes(neighbor)->children(3),nodes(neighbor)->children(4));
    }
    return VECTOR<int,4>(0,0,0,0);

}
//#####################################################################
// Neighbor_Same_Level 
//#####################################################################
template<class T> int OCTREE<T>::
Neighbor_Same_Level(const int index,const int n_of_6,int& increment)
{
    if(n_of_6==1) // offset=(1,0,0) or east neighbor
    {
        if(index==1) return 0;
        int parent_index=nodes(index)->parent_index;
        if(nodes(index)->n_of_8==3) return nodes(parent_index)->children(2);
        if(nodes(index)->n_of_8==4) return nodes(parent_index)->children(1);
        if(nodes(index)->n_of_8==7) return nodes(parent_index)->children(6);
        if(nodes(index)->n_of_8==8) return nodes(parent_index)->children(5);
        int u=Neighbor_Same_Level(parent_index,n_of_6);
        if(u==0 || nodes(u)->children(1)==0){increment++;return u;}
        else{
            if(nodes(index)->n_of_8==2) return nodes(u)->children(3);
            if(nodes(index)->n_of_8==1) return nodes(u)->children(4);
            if(nodes(index)->n_of_8==6) return nodes(u)->children(7);
            if(nodes(index)->n_of_8==5) return nodes(u)->children(8);
        }
    }
    if(n_of_6==2) // west
    {
        if(index==1) return 0;
        int parent_index=nodes(index)->parent_index;
        if(nodes(index)->n_of_8==2) return nodes(parent_index)->children(3);
        if(nodes(index)->n_of_8==1) return nodes(parent_index)->children(4);
        if(nodes(index)->n_of_8==6) return nodes(parent_index)->children(7);
        if(nodes(index)->n_of_8==5) return nodes(parent_index)->children(8);
        int u=Neighbor_Same_Level(parent_index,n_of_6);
        if(u==0 || nodes(u)->children(1)==0){increment++;return u;}
        else{
            if(nodes(index)->n_of_8==3) return nodes(u)->children(2);
            if(nodes(index)->n_of_8==4) return nodes(u)->children(1);
            if(nodes(index)->n_of_8==7) return nodes(u)->children(6);
            if(nodes(index)->n_of_8==8) return nodes(u)->children(5);
        }
    }
    if(n_of_6==3) // north 
    {
        if(index==1) return 0;
        int parent_index=nodes(index)->parent_index;
        if(nodes(index)->n_of_8==1) return nodes(parent_index)->children(2);
        if(nodes(index)->n_of_8==4) return nodes(parent_index)->children(3);
        if(nodes(index)->n_of_8==5) return nodes(parent_index)->children(6);
        if(nodes(index)->n_of_8==8) return nodes(parent_index)->children(7);
        int u=Neighbor_Same_Level(parent_index,n_of_6);
        if(u==0 || nodes(u)->children(1)==0){increment++;return u;} 
        else{
            if(nodes(index)->n_of_8==2) return nodes(u)->children(1);
            if(nodes(index)->n_of_8==3) return nodes(u)->children(4);
            if(nodes(index)->n_of_8==6) return nodes(u)->children(5);
            if(nodes(index)->n_of_8==7) return nodes(u)->children(8);
        }
    }
    if(n_of_6==4) // south
    {
        if(index==1) return 0;
        int parent_index=nodes(index)->parent_index;
        if(nodes(index)->n_of_8==2) return nodes(parent_index)->children(1);
        if(nodes(index)->n_of_8==3) return nodes(parent_index)->children(4);
        if(nodes(index)->n_of_8==6) return nodes(parent_index)->children(5);
        if(nodes(index)->n_of_8==7) return nodes(parent_index)->children(8);
        int u=Neighbor_Same_Level(parent_index,n_of_6);
        if(u==0 || nodes(u)->children(1)==0){increment++;return u;}
        else{
            if(nodes(index)->n_of_8==1) return nodes(u)->children(2);
            if(nodes(index)->n_of_8==4) return nodes(u)->children(3);
            if(nodes(index)->n_of_8==5) return nodes(u)->children(6);
            if(nodes(index)->n_of_8==8) return nodes(u)->children(7);
        }
    }
    if(n_of_6==5) // front
    {
        if(index==1) return 0;
        int parent_index=nodes(index)->parent_index;
        if(nodes(index)->n_of_8==5) return nodes(parent_index)->children(1);
        if(nodes(index)->n_of_8==6) return nodes(parent_index)->children(2);
        if(nodes(index)->n_of_8==7) return nodes(parent_index)->children(3);
        if(nodes(index)->n_of_8==8) return nodes(parent_index)->children(4);
        int u=Neighbor_Same_Level(parent_index,n_of_6);
        if(u==0 || nodes(u)->children(1)==0){increment++;return u;}
        else{
            if(nodes(index)->n_of_8==1) return nodes(u)->children(5);
            if(nodes(index)->n_of_8==2) return nodes(u)->children(6);
            if(nodes(index)->n_of_8==3) return nodes(u)->children(7);
            if(nodes(index)->n_of_8==4) return nodes(u)->children(8);
        }
    }
    if(n_of_6==6) // back
    {
        if(index==1) return 0;
        int parent_index=nodes(index)->parent_index;
        if(nodes(index)->n_of_8==1) return nodes(parent_index)->children(5);
        if(nodes(index)->n_of_8==2) return nodes(parent_index)->children(6);
        if(nodes(index)->n_of_8==3) return nodes(parent_index)->children(7);
        if(nodes(index)->n_of_8==4) return nodes(parent_index)->children(8);
        int u=Neighbor_Same_Level(parent_index,n_of_6);
        if(u==0 || nodes(u)->children(1)==0){increment++;return u;}
        else{
            if(nodes(index)->n_of_8==5) return nodes(u)->children(1);
            if(nodes(index)->n_of_8==6) return nodes(u)->children(2);
            if(nodes(index)->n_of_8==7) return nodes(u)->children(3);
            if(nodes(index)->n_of_8==8) return nodes(u)->children(4);
        }
    }
    return 1;
}
template<class T> int OCTREE<T>::
Neighbor_Same_Level(const int index,const int n_of_6)
{
    if(n_of_6==1) // offset=(1,0,0) or east neighbor
    {
        if(index==1) return 0;
        int parent_index=nodes(index)->parent_index;
        if(nodes(index)->n_of_8==3) return nodes(parent_index)->children(2);
        if(nodes(index)->n_of_8==4) return nodes(parent_index)->children(1);
        if(nodes(index)->n_of_8==7) return nodes(parent_index)->children(6);
        if(nodes(index)->n_of_8==8) return nodes(parent_index)->children(5);
        int u=Neighbor_Same_Level(parent_index,n_of_6);
        if(u==0 || nodes(u)->children(1)==0) return u;
        else{
            if(nodes(index)->n_of_8==2) return nodes(u)->children(3);
            if(nodes(index)->n_of_8==1) return nodes(u)->children(4);
            if(nodes(index)->n_of_8==6) return nodes(u)->children(7);
            if(nodes(index)->n_of_8==5) return nodes(u)->children(8);
        }
    }
    if(n_of_6==2) // west
    {
        if(index==1) return 0;
        int parent_index=nodes(index)->parent_index;
        if(nodes(index)->n_of_8==2) return nodes(parent_index)->children(3);
        if(nodes(index)->n_of_8==1) return nodes(parent_index)->children(4);
        if(nodes(index)->n_of_8==6) return nodes(parent_index)->children(7);
        if(nodes(index)->n_of_8==5) return nodes(parent_index)->children(8);
        int u=Neighbor_Same_Level(parent_index,n_of_6);
        if(u==0 || nodes(u)->children(1)==0) return u;
        else{
            if(nodes(index)->n_of_8==3) return nodes(u)->children(2);
            if(nodes(index)->n_of_8==4) return nodes(u)->children(1);
            if(nodes(index)->n_of_8==7) return nodes(u)->children(6);
            if(nodes(index)->n_of_8==8) return nodes(u)->children(5);
        }
    }
    if(n_of_6==3) // north 
    {
        if(index==1) return 0;
        int parent_index=nodes(index)->parent_index;
        if(nodes(index)->n_of_8==1) return nodes(parent_index)->children(2);
        if(nodes(index)->n_of_8==4) return nodes(parent_index)->children(3);
        if(nodes(index)->n_of_8==5) return nodes(parent_index)->children(6);
        if(nodes(index)->n_of_8==8) return nodes(parent_index)->children(7);
        int u=Neighbor_Same_Level(parent_index,n_of_6);
        if(u==0 || nodes(u)->children(1)==0) return u;
        else{
            if(nodes(index)->n_of_8==2) return nodes(u)->children(1);
            if(nodes(index)->n_of_8==3) return nodes(u)->children(4);
            if(nodes(index)->n_of_8==6) return nodes(u)->children(5);
            if(nodes(index)->n_of_8==7) return nodes(u)->children(8);
        }
    }
    if(n_of_6==4) // south
    {
        if(index==1) return 0;
        int parent_index=nodes(index)->parent_index;
        if(nodes(index)->n_of_8==2) return nodes(parent_index)->children(1);
        if(nodes(index)->n_of_8==3) return nodes(parent_index)->children(4);
        if(nodes(index)->n_of_8==6) return nodes(parent_index)->children(5);
        if(nodes(index)->n_of_8==7) return nodes(parent_index)->children(8);
        int u=Neighbor_Same_Level(parent_index,n_of_6);
        if(u==0 || nodes(u)->children(1)==0) return u;
        else{
            if(nodes(index)->n_of_8==1) return nodes(u)->children(2);
            if(nodes(index)->n_of_8==4) return nodes(u)->children(3);
            if(nodes(index)->n_of_8==5) return nodes(u)->children(6);
            if(nodes(index)->n_of_8==8) return nodes(u)->children(7);
        }
    }
    if(n_of_6==5) // front
    {
        if(index==1) return 0;
        int parent_index=nodes(index)->parent_index;
        if(nodes(index)->n_of_8==5) return nodes(parent_index)->children(1);
        if(nodes(index)->n_of_8==6) return nodes(parent_index)->children(2);
        if(nodes(index)->n_of_8==7) return nodes(parent_index)->children(3);
        if(nodes(index)->n_of_8==8) return nodes(parent_index)->children(4);
        int u=Neighbor_Same_Level(parent_index,n_of_6);
        if(u==0 || nodes(u)->children(1)==0) return u;
        else{
            if(nodes(index)->n_of_8==1) return nodes(u)->children(5);
            if(nodes(index)->n_of_8==2) return nodes(u)->children(6);
            if(nodes(index)->n_of_8==3) return nodes(u)->children(7);
            if(nodes(index)->n_of_8==4) return nodes(u)->children(8);
        }
    }
    if(n_of_6==6) // back
    {
        if(index==1) return 0;
        int parent_index=nodes(index)->parent_index;
        if(nodes(index)->n_of_8==1) return nodes(parent_index)->children(5);
        if(nodes(index)->n_of_8==2) return nodes(parent_index)->children(6);
        if(nodes(index)->n_of_8==3) return nodes(parent_index)->children(7);
        if(nodes(index)->n_of_8==4) return nodes(parent_index)->children(8);
        int u=Neighbor_Same_Level(parent_index,n_of_6);
        if(u==0 || nodes(u)->children(1)==0) return u;
        else{
            if(nodes(index)->n_of_8==5) return nodes(u)->children(1);
            if(nodes(index)->n_of_8==6) return nodes(u)->children(2);
            if(nodes(index)->n_of_8==7) return nodes(u)->children(3);
            if(nodes(index)->n_of_8==8) return nodes(u)->children(4);
        }
    }
    return 1;
}
//#####################################################################
// Multiply
//#####################################################################
template<class T> void OCTREE<T>::
Multiply(const ARRAY<T>& x,ARRAY<T>& b) 
{
    if(uniform) Multiply_Uniform(x,b);
    else Multiply_Non_Uniform(x,b);
}
template<class T> void OCTREE<T>::
Multiply_Uniform(const ARRAY<T>& x,ARRAY<T>& b) // Ax=b
{
    //b.Exact_Resize(x.Size());
    for(int j1=1;j1<=x.Size();++j1)
    {
        b(j1)=0;
        T p1=x(j1),p2;
        int index1=Array_Index_2_Tree_Index(j1);
        if(nodes(index1)->levelset>0 && nodes(index1)->Neumann_face_indices.Size()==0) continue;
        if(ghost(index1)){LOG::cout<<"sth is wrong!should not be ghost!"<<std::endl;exit(0);}
        T dx=nodes(index1)->dx_over_two(1)*2;
        for(int i=1;i<=6;++i){
            T distance=dx;
            int index2=neighbor_indices(index1)(i)(4);// for uniform we know y is the real index
            if(ghost(index2)) p2=0.;
            else{
                int j2=Tree_Index_2_Array_Index(index2);
                if(nodes(index2)->levelset>0 && nodes(index2)->Neumann_face_indices.Size()==0){
                    p2=p1;}
                else{
                    p2=x(j2);}
            }
            T g=(p2-p1)/distance;                
            b(j1)+=g*nodes(index1)->face_areas(i);
        }
    }
}
template<class T> void OCTREE<T>::
Multiply_Non_Uniform(const ARRAY<T>& x,ARRAY<T>& b) // Ax=b
{
#if 0
    for(int j=1;j<=x.Size();++j)
    {
        b(j)=0;
        T p1=x(j),p2;
        int index=Array_Index_2_Tree_Index(j);
        if(ghost(index)){LOG::cout<<"sth is wrong!should be no ghost!"<<std::endl;exit(0);}
        for(int i=1;i<=6;++i){
            VECTOR<int,4> indices=neighbor_indices(index)(i);
            int j2;
            if(indices(1)!=0){
                for(int k=1;k<=4;++k)
                {
                    int index2=indices(k);
                    if(ghost(index2)) p2=0.;
                    else{
                        j2=Tree_Index_2_Array_Index(index2);
                        p2=x(j2);
                    }
                    T g=(p2-p1)/(nodes(index)->dx_over_two(1)+nodes(index2)->dx_over_two(1));
                    b(j)+=g*nodes(index)->dx_over_two(1)*nodes(index)->dx_over_two(1);
                }
            }
            else{
                int index2=indices(4);
                if(nodes(index)->dx_over_two(1)!=nodes(index2)->dx_over_two(1))
                {
                    int index11,index22,index33,index44;
                    int ii=0;
                    switch(i)
                    {
                        case 1: ii=2;break;
                        case 2: ii=1;break;
                        case 3: ii=4;break;
                        case 4: ii=3;break;
                        case 5: ii=6;break;
                        case 6:default:ii=5;
                    }
                    index11=neighbor_indices(indices(4))(ii)(1);
                    index22=neighbor_indices(indices(4))(ii)(2);
                    index33=neighbor_indices(indices(4))(ii)(3);
                    index44=neighbor_indices(indices(4))(ii)(4);
                    T p1_tmp=x(Tree_Index_2_Array_Index(index11))/4+x(Tree_Index_2_Array_Index(index22))/4+x(Tree_Index_2_Array_Index(index33))/4+x(Tree_Index_2_Array_Index(index44))/4;
   
                    if(ghost(index2)) p2=0.;
                    else 
                    {
                        j2=Tree_Index_2_Array_Index(index2);
                        p2=x(j2);
                    }
                    T g=(p2-p1_tmp)/(nodes(index)->dx_over_two(1)+nodes(index2)->dx_over_two(1));
                    b(j)+=g*nodes(index)->dx_over_two(1)*2*nodes(index)->dx_over_two(1)*2;
                }
                else
                {
                    if(ghost(index2)) p2=0.;
                    else{
                        j2=Tree_Index_2_Array_Index(index2);
                        p2=x(j2);
                    }
                    T g=(p2-p1)/(nodes(index)->dx_over_two(1)+nodes(index2)->dx_over_two(1));
                    b(j)+=g*nodes(index)->dx_over_two(1)*2*nodes(index)->dx_over_two(1)*2;
                }
            }
        }
    }
   

#else
    for(int j=1;j<=x.Size();++j)
    {
        b(j)=0;
        T p1=x(j),p2;
        int index=Array_Index_2_Tree_Index(j);
        OCTREE_NODE<T>* node=nodes(index);
        if(nodes(index)->levelset>0 && nodes(index)->Neumann_face_indices.Size()==0) continue;
        else if((nodes(index)->center-TV(.5,.5,.5)).Magnitude()<0.3)
        {
            for(int i=1;i<=voronoi_diagram->neighbors(node->voronoi_index).Size();++i){
                int index2=voronoi_index_2_tree_index(voronoi_diagram->neighbors(node->voronoi_index)(i));
                T distance=(nodes(index)->center-nodes(index2)->center).Magnitude();
                if(ghost(index2)) p2=0.;
                else{
                    if(nodes(index2)->levelset>0 && nodes(index2)->Neumann_face_indices.Size()==0)
                        p2=p1;
                    else
                        p2=x(Tree_Index_2_Array_Index(index2));
                }
                T g=(p2-p1)/distance;                
                b(j)+=g*node->face_areas(i);
            }

        }
        else
        {
                for(int i=1;i<=6;++i){
                    VECTOR<int,4> indices=neighbor_indices(index)(i);
                    int j2;
                    if(indices(1)!=0){
                        for(int k=1;k<=4;++k)
                        {
                            int index2=indices(k);
                            if(ghost(index2)) p2=0.;
                            else{
                                j2=Tree_Index_2_Array_Index(index2);
                                p2=x(j2);
                            }
                            T g=(p2-p1)/(nodes(index)->dx_over_two(1)+nodes(index2)->dx_over_two(1));
                            b(j)+=g*nodes(index)->dx_over_two(1)*nodes(index)->dx_over_two(1);
                        }
                    }
                    else{
                        int index2=indices(4);
                        if(nodes(index)->dx_over_two(1)!=nodes(index2)->dx_over_two(1))
                        {
                            int index11,index22,index33,index44;
                            int ii=0;
                            switch(i)
                            {
                                case 1: ii=2;break;
                                case 2: ii=1;break;
                                case 3: ii=4;break;
                                case 4: ii=3;break;
                                case 5: ii=6;break;
                                case 6:default:ii=5;
                            }
                            index11=neighbor_indices(indices(4))(ii)(1);
                            index22=neighbor_indices(indices(4))(ii)(2);
                            index33=neighbor_indices(indices(4))(ii)(3);
                            index44=neighbor_indices(indices(4))(ii)(4);
                            T p1_tmp=x(Tree_Index_2_Array_Index(index11))/4+x(Tree_Index_2_Array_Index(index22))/4+x(Tree_Index_2_Array_Index(index33))/4+x(Tree_Index_2_Array_Index(index44))/4;
   
                            if(ghost(index2)) p2=0.;
                            else 
                            {
                                j2=Tree_Index_2_Array_Index(index2);
                                p2=x(j2);
                            }
                            T g=(p2-p1_tmp)/(nodes(index)->dx_over_two(1)+nodes(index2)->dx_over_two(1));
                            b(j)+=g*nodes(index)->dx_over_two(1)*2*nodes(index)->dx_over_two(1)*2;
                        }
                        else
                        {
                            if(ghost(index2)) p2=0.;
                            else{
                                j2=Tree_Index_2_Array_Index(index2);
                                p2=x(j2);
                            }
                            T g=(p2-p1)/(nodes(index)->dx_over_two(1)+nodes(index2)->dx_over_two(1));
                            b(j)+=g*nodes(index)->dx_over_two(1)*2*nodes(index)->dx_over_two(1)*2;
                        }
                    }
                }
        }
    }
#endif
}
//#####################################################################
// Project
//#####################################################################
template<class T> void OCTREE<T>::
Project(ARRAY<T>& x)
{
    if(uniform) Project_Uniform(x);
    else Project_Non_Uniform(x);
}
template<class T> void OCTREE<T>::
Project_Uniform(ARRAY<T>& x)
{
    for(int j1=1;j1<=x.Size();++j1)
    {
        int index1=Array_Index_2_Tree_Index(j1);
        if(nodes(index1)->levelset>0 && nodes(index1)->Neumann_face_indices.Size()==0) x(j1)=0.;
    }
}
template<class T> void OCTREE<T>::
Project_Non_Uniform(ARRAY<T>& x)
{
    for(int j1=1;j1<=x.Size();++j1)
    {
        int index1=Array_Index_2_Tree_Index(j1);
        if(nodes(index1)->levelset>0 && nodes(index1)->Neumann_face_indices.Size()==0) x(j1)=0.;
    }
}
//#####################################################################
// Update_Volume_And_Area_For_Neumann
//#####################################################################
template<class T> T OCTREE<T>::
Phi_Circle(TV position)
{
    return 0.15111111111-(position-TV(.5,.5,.5)).Magnitude();
    //return -1.;
}
template<class T>
void Face_Barycenter_Helper(const ARRAY<VECTOR<T,3> >& face,VECTOR<T,3>& centroid,T& area)
{
    centroid=VECTOR<T,3>(0.,0.,0.);
    area=0.;
    typedef VECTOR<T,3> TV;
    TV base_vertex=face(1);
    for(int k=1;k<=face.Size();++k){int next=(k==face.Size())?1:k+1;
        T current_area=TV::Cross_Product(face(k)-base_vertex,face(next)-base_vertex).Magnitude()*(T).5;
        centroid+=current_area*(base_vertex+face(k)+face(next))*one_third;
        area+=current_area;}
    centroid/=area;
}
template<class T> void OCTREE<T>::
Update_Volume_And_Area_For_Neumann()
{
    for(int i=1;i<=nodes.Size();++i)
    {
        if(nodes(i)->children(1)!=0) continue;
        if(ghost(i)) continue;
        OCTREE_NODE<T>* node=nodes(i);
        node->new_volume=node->volume;

        ARRAY<TV> intersecting_points;
        ARRAY<VECTOR<TV,2> > edges_visited;
        ARRAY<TV> centroids_for_volume(voronoi_diagram->face_vertices(node->voronoi_index).Size());
        ARRAY<TV> normals_for_volume(voronoi_diagram->face_vertices(node->voronoi_index).Size());
        ARRAY<T> areas_for_volume(voronoi_diagram->face_vertices(node->voronoi_index).Size());
        ARRAY<TV> centroids_for_volume2(voronoi_diagram->face_vertices(node->voronoi_index).Size());
        ARRAY<TV> normals_for_volume2(voronoi_diagram->face_vertices(node->voronoi_index).Size());
        ARRAY<T> areas_for_volume2(voronoi_diagram->face_vertices(node->voronoi_index).Size());
        ARRAY<VECTOR<TV,2> > intersecting_pairs;
        for(int j=1;j<=voronoi_diagram->face_vertices(node->voronoi_index).Size();++j)
        {
            ARRAY<TV> updated_face;
            ARRAY<TV> updated_face2;
            const ARRAY<TV>& vertices=voronoi_diagram->face_vertices(node->voronoi_index)(j);
            ARRAY<T> phis(vertices.Size()); 
            for(int k=1;k<=vertices.Size();++k) phis(k)=Phi_Circle(vertices(k));

            int counter=1;
            VECTOR<TV,2> pair;
            for(int k=1;k<=vertices.Size();++k)
            {
                int e1=k,e2=k==vertices.Size()?1:k+1;
                VECTOR<TV,2> tmp1=VECTOR<TV,2>(vertices(e1),vertices(e2));;
                VECTOR<TV,2> tmp2=VECTOR<TV,2>(vertices(e2),vertices(e1));;
                bool flag=false;
                for(int l=1;l<=edges_visited.Size();++l)
                    if(tmp1==edges_visited(l) || tmp2==edges_visited(l)){flag=true;break;}
               if(flag) edges_visited.Append(tmp1);

                T phi1=phis(e1),phi2=phis(e2);
                if(phi1>0) updated_face.Append_Unique(vertices(e1));
                else updated_face2.Append_Unique(vertices(e1));
                if((phi1<0 && phi2>0) || (phi1>0 && phi2<0))
                {
                    node->Neumann_face_indices.Append_Unique(j);
                    T theta=LEVELSET_UTILITIES<T>::Theta(phi1,phi2);
                    TV tmp=vertices(e1)+(vertices(e2)-vertices(e1))*theta;
                    updated_face.Append(tmp);
                    updated_face2.Append(tmp);
                    pair(counter++)=tmp;
                }
                if(phi2>0) updated_face.Append_Unique(vertices(e2));
                else updated_face2.Append_Unique(vertices(e2));
            }
            if(counter==3) intersecting_pairs.Append(pair);
            if(updated_face.Size()>2)
            {
                Face_Barycenter_Helper(updated_face,centroids_for_volume(j),areas_for_volume(j));
                int index2=voronoi_index_2_tree_index((voronoi_diagram->neighbors(node->voronoi_index))(j));
                normals_for_volume(j)=(nodes(index2)->center-node->center).Normalized();
            }
            if(updated_face2.Size()>2)
            {
                Face_Barycenter_Helper(updated_face2,centroids_for_volume2(j),areas_for_volume2(j));
                int index2=voronoi_index_2_tree_index((voronoi_diagram->neighbors(node->voronoi_index))(j));
                normals_for_volume2(j)=(nodes(index2)->center-node->center).Normalized();
            }
            node->face_areas(j)-=areas_for_volume(j);
            node->face_centers(j)=centroids_for_volume2(j);
        }
        if(intersecting_pairs.Size()==0) continue;
        // make intersecting_points ordered
        intersecting_points.Append(intersecting_pairs(1).x);
        intersecting_points.Append(intersecting_pairs(1).y);
        TV current=intersecting_pairs(1).y;
        intersecting_pairs.Remove_Index_Lazy(1);
        while(intersecting_pairs.Size()>1)
        {
            for(int j=1;j<=intersecting_pairs.Size();++j)
            {
                if(std::fabs(TV::Dot_Product(intersecting_pairs(j).x,current)-current.Magnitude_Squared())<1e-15*current.Magnitude_Squared()){
                    intersecting_points.Append(intersecting_pairs(j).y);
                    current=intersecting_pairs(j).y;
                    intersecting_pairs.Remove_Index_Lazy(j);
                    break;
                }
                if(std::fabs(TV::Dot_Product(intersecting_pairs(j).y,current)-current.Magnitude_Squared())<1e-15*current.Magnitude_Squared()){
                    intersecting_points.Append(intersecting_pairs(j).x);
                    current=intersecting_pairs(j).x;
                    intersecting_pairs.Remove_Index_Lazy(j);
                    break;
                }
            }
        }
        if(intersecting_points.Size()>2)
        {
            T cutted_volume=0;
            for(int j=1;j<=normals_for_volume.Size();++j)
                cutted_volume+=TV::Dot_Product(centroids_for_volume(j),normals_for_volume(j))*areas_for_volume(j);

            ARRAY<TV> tri(3);
            tri(1)=TV();
            for(int j=1;j<=intersecting_points.Size();++j) tri(1)+=intersecting_points(j);
            tri(1)/=intersecting_points.Size();
            for(int j=1;j<=intersecting_points.Size();++j)
            {
                int e1=j,e2=j==intersecting_points.Size()?1:j+1;
                tri(2)=intersecting_points(e1);
                tri(3)=intersecting_points(e2);
                TV center;T area;
                Face_Barycenter_Helper(tri,center,area);
                node->new_face_center.Append(center);
                node->new_face_area.Append(area);
                TV normal=TV::Cross_Product(tri(2)-tri(1),tri(3)-tri(1)).Normalized();
                if(TV::Dot_Product(normal,tri(1)-TV(.5,.5,.5))>0.) normal*=-1.;
                node->new_face_normal.Append(normal);
                cutted_volume+=TV::Dot_Product(center,-normal)*area;
            }
            cutted_volume/=3.;
            node->new_volume=node->volume-cutted_volume;
        }
#if 0
        //if(node->Neumann_face_indices.Size())
        if(false)
        {
            node->new_face_centers=node->face_centers;
            int n=intersecting_points.Size();
            if(bad>0) LOG::cout << "bad " << bad << "/" << intersecting_points.Size() << std::endl;
            if(n<3){LOG::cout << "????" << std::endl;exit(0);}
            Eigen::Matrix<T,Eigen::Dynamic,3> X(n,3);
            Eigen::Matrix<T,1,3> u(1,3);
            for(int jj=0;jj<3;++jj) u(0,jj)=0.;
            for(int ii=0;ii<n;++ii)
            {
                for(int jj=0;jj<3;++jj){
                    X(ii,jj)=intersecting_points(ii+1)(jj+1);
                    u(0,jj)+=X(ii,jj);
                }
            }
            for(int jj=0;jj<3;++jj) u(0,jj)/=n;
            Eigen::Matrix<T,Eigen::Dynamic,3> B(n,3);
            Eigen::Matrix<T,Eigen::Dynamic,1> h(n,1);
            for(int ii=0;ii<n;++ii) h(ii,0)=1.;
            B=X-h*u;
            Eigen::Matrix<T,3,3> C(3,3);
            C=B.transpose()*B/(n-1);
            Eigen::EigenSolver<Eigen::Matrix<T,3,3> > evd_solver(C);
            Eigen::Matrix<T,3,3>  V = evd_solver.eigenvectors().real();
            Eigen::Matrix<T,3,1>  eigs = evd_solver.eigenvalues().real();
            int smallest=0;
            for(int ii=1;ii<3;++ii) if(eigs(smallest,0)>eigs(ii,0)) smallest=ii;
            TV normal=TV(V.col(smallest)(0),V.col(smallest)(1),V.col(smallest)(2)); 
            TV point=TV(u(0,0),u(0,1),u(0,2));
            if(TV::Dot_Product(normal,point-TV(.5,.5,.5))>0.) normal*=-1.;
            node->new_face_normal=normal.Normalized();
            PLANE<T> plane(normal,point);

            // for neumann faces, recalculate 
            intersecting_points.Resize(0);
            edges_visited.Resize(0);
            ARRAY<TV> centroids_for_volume(voronoi_diagram->face_vertices(node->voronoi_index).Size());
            ARRAY<TV> normals_for_volume(voronoi_diagram->face_vertices(node->voronoi_index).Size());
            ARRAY<T> areas_for_volume(voronoi_diagram->face_vertices(node->voronoi_index).Size());
            ARRAY<TV> centroids_for_volume2(voronoi_diagram->face_vertices(node->voronoi_index).Size());
            ARRAY<TV> normals_for_volume2(voronoi_diagram->face_vertices(node->voronoi_index).Size());
            ARRAY<T> areas_for_volume2(voronoi_diagram->face_vertices(node->voronoi_index).Size());
            ARRAY<VECTOR<TV,2> > intersecting_pairs;
            for(int j=1;j<=voronoi_diagram->face_vertices(node->voronoi_index).Size();++j)
            {
                ARRAY<TV> updated_face;
                ARRAY<TV> updated_face2;
                const ARRAY<TV>& vertices=voronoi_diagram->face_vertices(node->voronoi_index)(j);
                ARRAY<T> phis(vertices.Size()); 
                for(int k=1;k<=vertices.Size();++k) phis(k)=Phi_Circle(vertices(k));
                
                int counter=1;
                VECTOR<TV,2> pair;
                for(int k=1;k<=vertices.Size();++k)
                {
                    int e1=k,e2=k==vertices.Size()?1:k+1;
                    VECTOR<TV,2> tmp1=VECTOR<TV,2>(vertices(e1),vertices(e2));;
                    VECTOR<TV,2> tmp2=VECTOR<TV,2>(vertices(e2),vertices(e1));;
                    bool flag=false;
                    for(int l=1;l<=edges_visited.Size();++l)
                        if(tmp1==edges_visited(l) || tmp2==edges_visited(l)){flag=true;break;}
                    if(!flag) edges_visited.Append(tmp1);

                    T phi1=phis(e1),phi2=phis(e2);
                    if(phi1>0) updated_face.Append_Unique(vertices(e1));
                    else updated_face2.Append_Unique(vertices(e1));
                    
                    if((phi1<0 && phi2>0) || (phi1>0 && phi2<0))
                    {
                        T theta;
                        assert(plane.Segment_Plane_Intersection(vertices(e1),vertices(e2),theta));
                        //plane.Segment_Plane_Intersection(vertices(e1),vertices(e2),theta);

                        LOG::cout << "new theta " << theta << std::endl;
                        LOG::cout << "old theta " <<LEVELSET_UTILITIES<T>::Theta(phi1,phi2) << std::endl;
                        TV tmp=vertices(e1)+(vertices(e2)-vertices(e1))*theta;
                        updated_face.Append(tmp);
                        updated_face2.Append(tmp);
                        pair(counter++)=tmp;
                    }
                    if(phi2>0) updated_face.Append_Unique(vertices(e2));
                    else updated_face2.Append_Unique(vertices(e2));
                }
                if(counter==3) intersecting_pairs.Append(pair);
                if(updated_face.Size()>2)
                {
                    Face_Barycenter_Helper(updated_face,centroids_for_volume(j),areas_for_volume(j));
                    node->face_areas(j)-=areas_for_volume(j);
                    int index2=voronoi_index_2_tree_index((voronoi_diagram->neighbors(node->voronoi_index))(j));
                    normals_for_volume(j)=(nodes(index2)->center-node->center).Normalized();
                }
                if(updated_face2.Size()>2)
                {
                    Face_Barycenter_Helper(updated_face2,centroids_for_volume2(j),areas_for_volume2(j));
                    node->new_face_centers(j)=centroids_for_volume2(j);
                    int index2=voronoi_index_2_tree_index((voronoi_diagram->neighbors(node->voronoi_index))(j));
                    normals_for_volume2(j)=(nodes(index2)->center-node->center).Normalized();
                }
            }
            // make intersecting_points ordered
            intersecting_points.Append(intersecting_pairs(1).x);
            intersecting_points.Append(intersecting_pairs(1).y);
            TV current=intersecting_pairs(1).y;
            intersecting_pairs.Remove_Index_Lazy(1);
            while(intersecting_pairs.Size()>1)
            {
                for(int j=1;j<=intersecting_pairs.Size();++j)
                {
                    if(std::fabs(TV::Dot_Product(intersecting_pairs(j).x,current)-current.Magnitude_Squared())<1e-15*current.Magnitude_Squared()){
                        intersecting_points.Append(intersecting_pairs(j).y);
                        current=intersecting_pairs(j).y;
                        intersecting_pairs.Remove_Index_Lazy(j);
                        break;
                    }
                    if(std::fabs(TV::Dot_Product(intersecting_pairs(j).y,current)-current.Magnitude_Squared())<1e-15*current.Magnitude_Squared()){
                        intersecting_points.Append(intersecting_pairs(j).x);
                        current=intersecting_pairs(j).x;
                        intersecting_pairs.Remove_Index_Lazy(j);
                        break;
                    }
                }
            }
            if(intersecting_points.Size()>2)
            {
                Face_Barycenter_Helper(intersecting_points,node->new_face_center,node->new_face_area);
                assert(node->new_face_area>0);
                T cutted_volume=0;
                for(int j=1;j<=normals_for_volume.Size();++j)
                    cutted_volume+=TV::Dot_Product(centroids_for_volume(j),normals_for_volume(j))*areas_for_volume(j);
                cutted_volume+=TV::Dot_Product(node->new_face_center,-node->new_face_normal)*node->new_face_area;
                cutted_volume/=3.;
                node->new_volume=node->volume-cutted_volume;
            }
        }
#endif
    }
    return;
}
//#####################################################################
// Functions from Hanan Samet's paper
//#####################################################################
template<class T> int OCTREE<T>::
S_FATHER(const int index)
{
    return nodes(index)->parent_index;
}
template<class T> int OCTREE<T>::
S_SON(const int index,const int octant)
{
    if(nodes(index)->children(1)!=0)
    {
        switch(octant)
        {
            case 145:return nodes(index)->children(1);break;
            case 135:return nodes(index)->children(2);break;
            case 235:return nodes(index)->children(3);break;
            case 245:return nodes(index)->children(4);break;
            case 146:return nodes(index)->children(5);break;
            case 136:return nodes(index)->children(6);break;
            case 236:return nodes(index)->children(7);break;
            case 246:return nodes(index)->children(8);break;
            default:LOG::cout<<"error S_SON"<<std::endl;
        }
    }
    return index;
}
template<class T> int OCTREE<T>::
S_SONTYPE(const int index)
{
    int sontype=nodes(index)->n_of_8;
    switch(sontype)
    {
        case 1:return 145;break;
        case 2:return 135;break;
        case 3:return 235;break;
        case 4:return 245;break;
        case 5:return 146;break;
        case 6:return 136;break;
        case 7:return 236;break;
        case 8:return 246;break;
        default:LOG::cout<<"error S_SONTYPE"<<std::endl;
    }
    return 0;
}
template<class T> bool OCTREE<T>::
S_ADJ(const int direction,const int octant)
{
    switch(direction)
    {
        case 2:
            switch(octant){
                case 246:return true ;break;case 245:return true ;break;case 236:return true ;break;case 235:return true ;break;
                case 146:return false;break;case 145:return false;break;case 136:return false;break;case 135:return false;break;default:return false;}break;
        case 1:
            switch(octant){
                case 246:return false;break;case 245:return false;break;case 236:return false;break;case 235:return false;break;
                case 146:return true ;break;case 145:return true ;break;case 136:return true ;break;case 135:return true ;break;default:return false;}break;
        case 4:
            switch(octant){
                case 246:return true ;break;case 245:return true ;break;case 236:return false;break;case 235:return false;break;
                case 146:return true ;break;case 145:return true ;break;case 136:return false;break;case 135:return false;break;default:return false;}break;
        case 3:
            switch(octant){
                case 246:return false;break;case 245:return false;break;case 236:return true ;break;case 235:return true ;break;
                case 146:return false;break;case 145:return false;break;case 136:return true ;break;case 135:return true ;break;default:return false;}break;
        case 6:
            switch(octant){
                case 246:return true ;break;case 245:return false;break;case 236:return true ;break;case 235:return false;break;
                case 146:return true ;break;case 145:return false;break;case 136:return true ;break;case 135:return false;break;default:return false;}break;
        case 5:
            switch(octant){
                case 246:return false;break;case 245:return true ;break;case 236:return false;break;case 235:return true ;break;
                case 146:return false;break;case 145:return true ;break;case 136:return false;break;case 135:return true ;break;default:return false;}break;
        case 24:
            switch(octant){
                case 246:return true ;break;case 245:return true ;break;case 236:return false;break;case 235:return false;break;
                case 146:return false;break;case 145:return false;break;case 136:return false;break;case 135:return false;break;default:return false;}break;
        case 23:
            switch(octant){
                case 246:return false;break;case 245:return false;break;case 236:return true ;break;case 235:return true ;break;
                case 146:return false;break;case 145:return false;break;case 136:return false;break;case 135:return false;break;default:return false;}break;
        case 26:
            switch(octant){
                case 246:return true ;break;case 245:return false;break;case 236:return true ;break;case 235:return false;break;
                case 146:return false;break;case 145:return false;break;case 136:return false;break;case 135:return false;break;default:return false;}break;
        case 25:
            switch(octant){
                case 246:return false;break;case 245:return true ;break;case 236:return false;break;case 235:return true ;break;
                case 146:return false;break;case 145:return false;break;case 136:return false;break;case 135:return false;break;default:return false;}break;
        case 14:
            switch(octant){
                case 246:return false;break;case 245:return false;break;case 236:return false;break;case 235:return false;break;
                case 146:return true ;break;case 145:return true ;break;case 136:return false;break;case 135:return false;break;default:return false;}break;
        case 13:
            switch(octant){
                case 246:return false;break;case 245:return false;break;case 236:return false;break;case 235:return false;break;
                case 146:return false;break;case 145:return false;break;case 136:return true ;break;case 135:return true ;break;default:return false;}break;
        case 16:
            switch(octant){
                case 246:return false;break;case 245:return false;break;case 236:return false;break;case 235:return false;break;
                case 146:return true ;break;case 145:return false;break;case 136:return true ;break;case 135:return false;break;default:return false;}break;
        case 15:
            switch(octant){
                case 246:return false;break;case 245:return false;break;case 236:return false;break;case 235:return false;break;
                case 146:return false;break;case 145:return true ;break;case 136:return false;break;case 135:return true ;break;default:return false;}break;
        case 46:
            switch(octant){
                case 246:return true ;break;case 245:return false;break;case 236:return false;break;case 235:return false;break;
                case 146:return true ;break;case 145:return false;break;case 136:return false;break;case 135:return false;break;default:return false;}break;
        case 45:
            switch(octant){
                case 246:return false;break;case 245:return true ;break;case 236:return false;break;case 235:return false;break;
                case 146:return false;break;case 145:return true ;break;case 136:return false;break;case 135:return false;break;default:return false;}break;
        case 36:
            switch(octant){
                case 246:return false;break;case 245:return false;break;case 236:return true ;break;case 235:return false;break;
                case 146:return false;break;case 145:return false;break;case 136:return true ;break;case 135:return false;break;default:return false;}break;
        case 35:
            switch(octant){
                case 246:return false;break;case 245:return false;break;case 236:return false;break;case 235:return true ;break;
                case 146:return false;break;case 145:return false;break;case 136:return false;break;case 135:return true ;break;default:return false;}break;
        case 246:
            switch(octant){
                case 246:return true ;break;case 245:return false;break;case 236:return false;break;case 235:return false;break;
                case 146:return false;break;case 145:return false;break;case 136:return false;break;case 135:return false;break;default:return false;}break;
        case 245:
            switch(octant){
                case 246:return false;break;case 245:return true ;break;case 236:return false;break;case 235:return false;break;
                case 146:return false;break;case 145:return false;break;case 136:return false;break;case 135:return false;break;default:return false;}break;
        case 236:
            switch(octant){
                case 246:return false;break;case 245:return false;break;case 236:return true ;break;case 235:return false;break;
                case 146:return false;break;case 145:return false;break;case 136:return false;break;case 135:return false;break;default:return false;}break;
        case 235:
            switch(octant){
                case 246:return false;break;case 245:return false;break;case 236:return false;break;case 235:return true ;break;
                case 146:return false;break;case 145:return false;break;case 136:return false;break;case 135:return false;break;default:return false;}break;
        case 146:
            switch(octant){
                case 246:return false;break;case 245:return false;break;case 236:return false;break;case 235:return false;break;
                case 146:return true ;break;case 145:return false;break;case 136:return false;break;case 135:return false;break;default:return false;}break;
        case 145:
            switch(octant){
                case 246:return false;break;case 245:return false;break;case 236:return false;break;case 235:return false;break;
                case 146:return false;break;case 145:return true ;break;case 136:return false;break;case 135:return false;break;default:return false;}break;
        case 136:
            switch(octant){
                case 246:return false;break;case 245:return false;break;case 236:return false;break;case 235:return false;break;
                case 146:return false;break;case 145:return false;break;case 136:return true ;break;case 135:return false;break;default:return false;}break;
        case 135:
            switch(octant){
                case 246:return false;break;case 245:return false;break;case 236:return false;break;case 235:return false;break;
                case 146:return false;break;case 145:return false;break;case 136:return false;break;case 135:return true ;break;default:return false;}break;
        default:return false;
    }
    return false;
}
template<class T> int OCTREE<T>::
S_REFLECT(const int direction,const int octant)
{
    switch(direction)
    {
        case 2:
            switch(octant){
                case 246:return 146;break;case 245:return 145;break;case 236:return 136;break;case 235:return 135;break;
                case 146:return 246;break;case 145:return 245;break;case 136:return 236;break;case 135:return 235;break;default:LOG::cout<<"error S_COMMON_FACE"<<std::endl;}break;
        case 1:
            switch(octant){
                case 246:return 146;break;case 245:return 145;break;case 236:return 136;break;case 235:return 135;break;
                case 146:return 246;break;case 145:return 245;break;case 136:return 236;break;case 135:return 235;break;default:LOG::cout<<"error S_COMMON_FACE"<<std::endl;}break;
        case 4:
            switch(octant){
                case 246:return 236;break;case 245:return 235;break;case 236:return 246;break;case 235:return 245;break;
                case 146:return 136;break;case 145:return 135;break;case 136:return 146;break;case 135:return 145;break;default:LOG::cout<<"error S_COMMON_FACE"<<std::endl;}break;
        case 3:
            switch(octant){
                case 246:return 236;break;case 245:return 235;break;case 236:return 246;break;case 235:return 245;break;
                case 146:return 136;break;case 145:return 135;break;case 136:return 146;break;case 135:return 145;break;default:LOG::cout<<"error S_COMMON_FACE"<<std::endl;}break;
        case 6:
            switch(octant){
                case 246:return 245;break;case 245:return 246;break;case 236:return 235;break;case 235:return 236;break;
                case 146:return 145;break;case 145:return 146;break;case 136:return 135;break;case 135:return 136;break;default:LOG::cout<<"error S_COMMON_FACE"<<std::endl;}break;
        case 5:
            switch(octant){
                case 246:return 245;break;case 245:return 246;break;case 236:return 235;break;case 235:return 236;break;
                case 146:return 145;break;case 145:return 146;break;case 136:return 135;break;case 135:return 136;break;default:LOG::cout<<"error S_COMMON_FACE"<<std::endl;}break;
        case 24:
            switch(octant){
                case 246:return 136;break;case 245:return 135;break;case 236:return 146;break;case 235:return 145;break;
                case 146:return 236;break;case 145:return 235;break;case 136:return 246;break;case 135:return 245;break;default:LOG::cout<<"error S_COMMON_FACE"<<std::endl;}break;
        case 23:
            switch(octant){
                case 246:return 136;break;case 245:return 135;break;case 236:return 146;break;case 235:return 145;break;
                case 146:return 236;break;case 145:return 235;break;case 136:return 246;break;case 135:return 245;break;default:LOG::cout<<"error S_COMMON_FACE"<<std::endl;}break;
        case 26:
            switch(octant){
                case 246:return 145;break;case 245:return 146;break;case 236:return 135;break;case 235:return 136;break;
                case 146:return 245;break;case 145:return 246;break;case 136:return 235;break;case 135:return 236;break;default:LOG::cout<<"error S_COMMON_FACE"<<std::endl;}break;
        case 25:
            switch(octant){
                case 246:return 145;break;case 245:return 146;break;case 236:return 135;break;case 235:return 136;break;
                case 146:return 245;break;case 145:return 246;break;case 136:return 235;break;case 135:return 236;break;default:LOG::cout<<"error S_COMMON_FACE"<<std::endl;}break;
        case 14:
            switch(octant){
                case 246:return 136;break;case 245:return 135;break;case 236:return 146;break;case 235:return 145;break;
                case 146:return 236;break;case 145:return 235;break;case 136:return 246;break;case 135:return 245;break;default:LOG::cout<<"error S_COMMON_FACE"<<std::endl;}break;
        case 13:
            switch(octant){
                case 246:return 136;break;case 245:return 135;break;case 236:return 146;break;case 235:return 145;break;
                case 146:return 236;break;case 145:return 235;break;case 136:return 246;break;case 135:return 245;break;default:LOG::cout<<"error S_COMMON_FACE"<<std::endl;}break;
        case 16:
            switch(octant){
                case 246:return 145;break;case 245:return 146;break;case 236:return 135;break;case 235:return 136;break;
                case 146:return 245;break;case 145:return 246;break;case 136:return 235;break;case 135:return 236;break;default:LOG::cout<<"error S_COMMON_FACE"<<std::endl;}break;
        case 15:
            switch(octant){
                case 246:return 145;break;case 245:return 146;break;case 236:return 135;break;case 235:return 136;break;
                case 146:return 245;break;case 145:return 246;break;case 136:return 235;break;case 135:return 236;break;default:LOG::cout<<"error S_COMMON_FACE"<<std::endl;}break;
        case 46:
            switch(octant){
                case 246:return 235;break;case 245:return 236;break;case 236:return 245;break;case 235:return 246;break;
                case 146:return 135;break;case 145:return 136;break;case 136:return 145;break;case 135:return 146;break;default:LOG::cout<<"error S_COMMON_FACE"<<std::endl;}break;
        case 45:
            switch(octant){
                case 246:return 235;break;case 245:return 236;break;case 236:return 245;break;case 235:return 246;break;
                case 146:return 135;break;case 145:return 136;break;case 136:return 145;break;case 135:return 146;break;default:LOG::cout<<"error S_COMMON_FACE"<<std::endl;}break;
        case 36:
            switch(octant){
                case 246:return 235;break;case 245:return 236;break;case 236:return 245;break;case 235:return 246;break;
                case 146:return 135;break;case 145:return 136;break;case 136:return 145;break;case 135:return 146;break;default:LOG::cout<<"error S_COMMON_FACE"<<std::endl;}break;
        case 35:
            switch(octant){
                case 246:return 235;break;case 245:return 236;break;case 236:return 245;break;case 235:return 246;break;
                case 146:return 135;break;case 145:return 136;break;case 136:return 145;break;case 135:return 146;break;default:LOG::cout<<"error S_COMMON_FACE"<<std::endl;}break;
        case 246:
            switch(octant){
                case 246:return 135;break;case 245:return 136;break;case 236:return 145;break;case 235:return 146;break;
                case 146:return 235;break;case 145:return 236;break;case 136:return 245;break;case 135:return 246;break;default:LOG::cout<<"error S_COMMON_FACE"<<std::endl;}break;
        case 245:
            switch(octant){
                case 246:return 135;break;case 245:return 136;break;case 236:return 145;break;case 235:return 146;break;
                case 146:return 235;break;case 145:return 236;break;case 136:return 245;break;case 135:return 246;break;default:LOG::cout<<"error S_COMMON_FACE"<<std::endl;}break;
        case 236:
            switch(octant){
                case 246:return 135;break;case 245:return 136;break;case 236:return 145;break;case 235:return 146;break;
                case 146:return 235;break;case 145:return 236;break;case 136:return 245;break;case 135:return 246;break;default:LOG::cout<<"error S_COMMON_FACE"<<std::endl;}break;
        case 235:
            switch(octant){
                case 246:return 135;break;case 245:return 136;break;case 236:return 145;break;case 235:return 146;break;
                case 146:return 235;break;case 145:return 236;break;case 136:return 245;break;case 135:return 246;break;default:LOG::cout<<"error S_COMMON_FACE"<<std::endl;}break;
        case 146:
            switch(octant){
                case 246:return 135;break;case 245:return 136;break;case 236:return 145;break;case 235:return 146;break;
                case 146:return 235;break;case 145:return 236;break;case 136:return 245;break;case 135:return 246;break;default:LOG::cout<<"error S_COMMON_FACE"<<std::endl;}break;
        case 145:
            switch(octant){
                case 246:return 135;break;case 245:return 136;break;case 236:return 145;break;case 235:return 146;break;
                case 146:return 235;break;case 145:return 236;break;case 136:return 245;break;case 135:return 246;break;default:LOG::cout<<"error S_COMMON_FACE"<<std::endl;}break;
        case 136:
            switch(octant){
                case 246:return 135;break;case 245:return 136;break;case 236:return 145;break;case 235:return 146;break;
                case 146:return 235;break;case 145:return 236;break;case 136:return 245;break;case 135:return 246;break;default:LOG::cout<<"error S_COMMON_FACE"<<std::endl;}break;
        case 135:
            switch(octant){
                case 246:return 135;break;case 245:return 136;break;case 236:return 145;break;case 235:return 146;break;
                case 146:return 235;break;case 145:return 236;break;case 136:return 245;break;case 135:return 246;break;default:LOG::cout<<"error S_COMMON_FACE"<<std::endl;}break;
        default:LOG::cout<<"error S_REFLECT"<<std::endl;
    }
    return 0;
}
template<class T> int OCTREE<T>::
S_COMMON_FACE(const int direction,const int octant)
{
    switch(direction)
    {
        case 24:
            switch(octant){
                case 246:return 0;break;case 245:return 0;break;case 236:return 2;break;case 235:return 2;break;
                case 146:return 4;break;case 145:return 4;break;case 136:return 0;break;case 135:return 0;break;default:return 0;}break;
        case 23:
            switch(octant){
                case 246:return 2;break;case 245:return 2;break;case 236:return 0;break;case 235:return 0;break;
                case 146:return 0;break;case 145:return 0;break;case 136:return 3;break;case 135:return 3;break;default:return 0;}break;
        case 26:
            switch(octant){
                case 246:return 0;break;case 245:return 2;break;case 236:return 0;break;case 235:return 2;break;
                case 146:return 6;break;case 145:return 0;break;case 136:return 6;break;case 135:return 0;break;default:return 0;}break;
        case 25:
            switch(octant){
                case 246:return 2;break;case 245:return 0;break;case 236:return 2;break;case 235:return 0;break;
                case 146:return 0;break;case 145:return 5;break;case 136:return 0;break;case 135:return 5;break;default:return 0;}break;
        case 14:
            switch(octant){
                case 246:return 4;break;case 245:return 4;break;case 236:return 0;break;case 235:return 0;break;
                case 146:return 0;break;case 145:return 0;break;case 136:return 1;break;case 135:return 1;break;default:return 0;}break;
        case 13:
            switch(octant){
                case 246:return 0;break;case 245:return 0;break;case 236:return 3;break;case 235:return 3;break;
                case 146:return 1;break;case 145:return 1;break;case 136:return 0;break;case 135:return 0;break;default:return 0;}break;
        case 16:
            switch(octant){
                case 246:return 6;break;case 245:return 0;break;case 236:return 6;break;case 235:return 0;break;
                case 146:return 0;break;case 145:return 1;break;case 136:return 0;break;case 135:return 1;break;default:return 0;}break;
        case 15:
            switch(octant){
                case 246:return 0;break;case 245:return 5;break;case 236:return 0;break;case 235:return 5;break;
                case 146:return 1;break;case 145:return 0;break;case 136:return 1;break;case 135:return 0;break;default:return 0;}break;
        case 46:
            switch(octant){
                case 246:return 0;break;case 245:return 4;break;case 236:return 6;break;case 235:return 0;break;
                case 146:return 0;break;case 145:return 4;break;case 136:return 6;break;case 135:return 0;break;default:return 0;}break;
        case 45:
            switch(octant){
                case 246:return 4;break;case 245:return 0;break;case 236:return 0;break;case 235:return 5;break;
                case 146:return 4;break;case 145:return 0;break;case 136:return 0;break;case 135:return 5;break;default:return 0;}break;
        case 36:
            switch(octant){
                case 246:return 6;break;case 245:return 0;break;case 236:return 0;break;case 235:return 3;break;
                case 146:return 6;break;case 145:return 0;break;case 136:return 0;break;case 135:return 3;break;default:return 0;}break;
        case 35:
            switch(octant){
                case 246:return 0;break;case 245:return 5;break;case 236:return 3;break;case 235:return 0;break;
                case 146:return 0;break;case 145:return 5;break;case 136:return 3;break;case 135:return 0;break;default:return 0;}break;
        case 246:
            switch(octant){
                case 246:return 0;break;case 245:return 0;break;case 236:return 0;break;case 235:return 2;break;
                case 146:return 0;break;case 145:return 4;break;case 136:return 6;break;case 135:return 0;break;default:return 0;}break;
        case 245:
            switch(octant){
                case 246:return 0;break;case 245:return 0;break;case 236:return 2;break;case 235:return 0;break;
                case 146:return 4;break;case 145:return 0;break;case 136:return 0;break;case 135:return 5;break;default:return 0;}break;
        case 236:
            switch(octant){
                case 246:return 0;break;case 245:return 2;break;case 236:return 0;break;case 235:return 0;break;
                case 146:return 6;break;case 145:return 0;break;case 136:return 0;break;case 135:return 3;break;default:return 0;}break;
        case 235:
            switch(octant){
                case 246:return 2;break;case 245:return 0;break;case 236:return 0;break;case 235:return 0;break;
                case 146:return 0;break;case 145:return 5;break;case 136:return 3;break;case 135:return 0;break;default:return 0;}break;
        case 146:
            switch(octant){
                case 246:return 0;break;case 245:return 4;break;case 236:return 6;break;case 235:return 0;break;
                case 146:return 0;break;case 145:return 0;break;case 136:return 0;break;case 135:return 1;break;default:return 0;}break;
        case 145:
            switch(octant){
                case 246:return 4;break;case 245:return 0;break;case 236:return 0;break;case 235:return 5;break;
                case 146:return 0;break;case 145:return 0;break;case 136:return 1;break;case 135:return 0;break;default:return 0;}break;
        case 136:
            switch(octant){
                case 246:return 6;break;case 245:return 0;break;case 236:return 0;break;case 235:return 3;break;
                case 146:return 0;break;case 145:return 1;break;case 136:return 0;break;case 135:return 0;break;default:return 0;}break;
        case 135:
            switch(octant){
                case 246:return 0;break;case 245:return 5;break;case 236:return 3;break;case 235:return 0;break;
                case 146:return 1;break;case 145:return 0;break;case 136:return 0;break;case 135:return 0;break;default:return 0;}break;
        default:return 0;
    }
    return 0;
}
template<class T> int OCTREE<T>::
S_COMMON_EDGE(const int direction,const int octant)
{
    switch(direction)
    {
        case 246:
            switch(octant){
                case 246:return  0;break;case 245:return 24;break;case 236:return 26;break;case 235:return  0;break;
                case 146:return 46;break;case 145:return  0;break;case 136:return  0;break;case 135:return  0;break;default:return 0;}break;
        case 245:
            switch(octant){
                case 246:return 24;break;case 245:return  0;break;case 236:return  0;break;case 235:return 25;break;
                case 146:return  0;break;case 145:return 45;break;case 136:return  0;break;case 135:return  0;break;default:return 0;}break;
        case 236:
            switch(octant){
                case 246:return 26;break;case 245:return  0;break;case 236:return  0;break;case 235:return 23;break;
                case 146:return  0;break;case 145:return  0;break;case 136:return 36;break;case 135:return  0;break;default:return 0;}break;
        case 235:
            switch(octant){
                case 246:return  0;break;case 245:return 25;break;case 236:return 23;break;case 235:return  0;break;
                case 146:return  0;break;case 145:return  0;break;case 136:return  0;break;case 135:return 35;break;default:return 0;}break;
        case 146:
            switch(octant){
                case 246:return 46;break;case 245:return  0;break;case 236:return  0;break;case 235:return  0;break;
                case 146:return  0;break;case 145:return 14;break;case 136:return 16;break;case 135:return  0;break;default:return 0;}break;
        case 145:
            switch(octant){
                case 246:return  0;break;case 245:return 45;break;case 236:return  0;break;case 235:return  0;break;
                case 146:return 14;break;case 145:return  0;break;case 136:return  0;break;case 135:return 15;break;default:return 0;}break;
        case 136:
            switch(octant){
                case 246:return  0;break;case 245:return  0;break;case 236:return 36;break;case 235:return  0;break;
                case 146:return 16;break;case 145:return  0;break;case 136:return  0;break;case 135:return 13;break;default:return 0;}break;
        case 135:
            switch(octant){
                case 246:return  0;break;case 245:return  0;break;case 236:return  0;break;case 235:return 35;break;
                case 146:return  0;break;case 145:return 15;break;case 136:return 13;break;case 135:return  0;break;default:return 0;}break;
        default:return 0;
    }
    return 0;
}
template<class T> int OCTREE<T>::
S_FACE_NEIGHBOR(const int index,const int face_direction)
{
    assert(face_direction<10);
    int P=index,I=face_direction,Q;
    if(S_FATHER(P) && S_ADJ(I,S_SONTYPE(P))) Q=S_FACE_NEIGHBOR(S_FATHER(P),I);
    else Q=S_FATHER(P);
    if(Q) return S_SON(Q,S_REFLECT(I,S_SONTYPE(P)));
    else return Q;
}
template<class T> int OCTREE<T>::
S_EDGE_NEIGHBOR(const int index,const int edge_direction)
{
    assert(edge_direction>10 && edge_direction<100);
    int P=index,I=edge_direction,Q;
    if(S_FATHER(P)==0) Q=0;
    else if(S_ADJ(I,S_SONTYPE(P))) Q=S_EDGE_NEIGHBOR(S_FATHER(P),I);
    else if(S_COMMON_FACE(I,S_SONTYPE(P))!=0) Q=S_FACE_NEIGHBOR(S_FATHER(P),S_COMMON_FACE(I,S_SONTYPE(P)));
    else Q=S_FATHER(P);
    //LOG::cout << Q << std::endl;
    if(Q) return S_SON(Q,S_REFLECT(I,S_SONTYPE(P)));
    else return Q;
}
template<class T> int OCTREE<T>::
S_VERTEX_NEIGHBOR(const int index,const int vertex_direction)
{
    assert(vertex_direction>100);
    int P=index,I=vertex_direction,Q;
    if(S_FATHER(P)==0) Q=0;
    else if(S_ADJ(I,S_SONTYPE(P))) Q=S_VERTEX_NEIGHBOR(S_FATHER(P),I);
    else if(S_COMMON_EDGE(I,S_SONTYPE(P))!=0) Q=S_EDGE_NEIGHBOR(S_FATHER(P),S_COMMON_EDGE(I,S_SONTYPE(P)));
    else if(S_COMMON_FACE(I,S_SONTYPE(P))!=0) Q=S_FACE_NEIGHBOR(S_FATHER(P),S_COMMON_FACE(I,S_SONTYPE(P)));
    else Q=S_FATHER(P);
    //LOG::cout << Q << std::endl;
    if(Q) return S_SON(Q,S_REFLECT(I,S_SONTYPE(P)));
    else return Q;   
}
//#####################################################################
template class OCTREE_NODE<float>;
template class OCTREE<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class OCTREE_NODE<double>;
template class OCTREE<double>;
#endif
