//#####################################################################
// Copyright 2015, Mridul Aanjaneya, Ming Gao.
// This file is part of PhysBAM whose distribution is governed by the license 
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class QUADTREE
//#####################################################################
#ifndef __QUADTREE__
#define __QUADTREE__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Data_Structures/PAIR.h>
#include <PhysBAM_Tools/Data_Structures/HASHTABLE.h>
#include <PhysBAM_Tools/Math_Tools/RANGE.h>
#include <PhysBAM_Tools/Vectors/VECTOR_2D.h>

namespace PhysBAM{
template<class T>
class QUADTREE_NODE
{
    typedef VECTOR<T,2> TV;
  public:
    TV center;
    TV dx_over_two;
    int north_east_index,north_west_index,south_east_index,south_west_index;
    int parent_index;
    int index;
    T levelset;
    T volume_pressure;
    bool Neumann_of_cell;
    // dof_theta_of_edge
    // 1: full dof
    // 0: full Neumann boundary 
    // 0-1: the partial of dof, from node1
    // -1-0: the partial of dof, from node2
    ARRAY<T> dof_theta_of_edge; 
    TV new_point; 
    TV new_normal; 
    T new_length; 
    ARRAY<TV> position_of_node;
    ARRAY<TV> position_of_Neumann_node;
    TV air_node;
    ARRAY<TV> center_of_face;
    ARRAY<TV> normal_of_face;
    ARRAY<T> area_of_face;

    QUADTREE_NODE(const TV& center_input=TV(),const TV& dx_over_two_input=TV(),int north_east_index_input=0,int north_west_index_input=0,int south_east_index_input=0,int south_west_index_input=0,int parent_index_input=0,int index_input=0)
        :center(center_input),dx_over_two(dx_over_two_input),north_east_index(north_east_index_input),north_west_index(north_west_index_input),south_east_index(south_east_index_input),
        south_west_index(south_west_index_input),parent_index(parent_index_input),index(index_input)
    {}

};

template<class T>
class QUADTREE
{
    typedef VECTOR<T,2> TV;
  public:
    RANGE<TV> domain;
    bool uniform;
    ARRAY<QUADTREE_NODE<T>*> nodes;
    // 1-4 east north west south (north_east north_west south_east south_west)
    // pair (i,j) 
    //      when i=0, j is the only useful index, meaning only 1 neighbor is valid
    //      when i,j both are not zero: counterclockwise
    ARRAY<ARRAY<PAIR<int,int> > > neighbor_indices; 
    ARRAY<bool> ghost;
    HASHTABLE<int,int> tree_index_2_array_index;
    ARRAY<int> array_index_2_tree_index;
    HASHTABLE<int,T> Dirichlet_boundary;

    QUADTREE() {}

    void Clean_Memory()
    {nodes.Clean_Memory();}

//#####################################################################
    void Compute_Two_Level_Subdivision(const RANGE<TV>& domain,const int scale,bool uniform=false);
    void Compute_Neighbors();
    void Compute_Volume_And_Area_For_Pressure();
    int North_Neighbor_Same_Level(const int index) const;
    int South_Neighbor_Same_Level(const int index) const;
    int West_Neighbor_Same_Level(const int index) const;
    int East_Neighbor_Same_Level(const int index) const;
    PAIR<int,int> North_Neighbor_Leaf(const int index);
    PAIR<int,int> South_Neighbor_Leaf(const int index);
    PAIR<int,int> West_Neighbor_Leaf(const int index);
    PAIR<int,int> East_Neighbor_Leaf(const int index);
    int Check_Child_Number(int index) const;
    int Corner_Neighbor_Leaf_Coarse(int index) const;
    int Tree_Index_2_Array_Index(int tree_index);
    int Array_Index_2_Tree_Index(int array_index);
    bool Verify_First_of_Two_Finer_Neighbors(int coarser,int finer);
//#####################################################################
    void Project(ARRAY<T>& x);
    void Project_Uniform(ARRAY<T>& x);
    void Project_Non_Uniform(ARRAY<T>& x);
    void Multiply(const ARRAY<T>& x,ARRAY<T>& b);
    void Multiply_Uniform(const ARRAY<T>& x,ARRAY<T>& b);
    void Multiply_Non_Uniform(const ARRAY<T>& x,ARRAY<T>& b);
    T Phi_Circle(TV position);

};
}
#endif
