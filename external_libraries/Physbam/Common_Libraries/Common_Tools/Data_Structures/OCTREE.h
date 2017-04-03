//#####################################################################
// Copyright 2015, Mridul Aanjaneya, Ming Gao.
// This file is part of PhysBAM whose distribution is governed by the license 
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OCTREE
//#####################################################################
#ifndef __OCTREE__
#define __OCTREE__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Math_Tools/RANGE.h>
#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Tools/Data_Structures/HASHTABLE.h>
#include <Common_Geometry/Topology_Based_Geometry/VORONOI_DIAGRAM.h>

namespace PhysBAM{
template<class T>
class OCTREE_NODE
{
    typedef VECTOR<T,3> TV;
  public:
    TV center,dx_over_two;
    // 8 children
    // front 4 : south_east_index -> counter clock wise 
    // back 4 : south_east_index -> counter clock wise 
    ARRAY<int> children;
    int parent_index;
    int index;
    int n_of_8;
    T levelset;

    // power diagram
    int voronoi_index; // starting from 1
    ARRAY<TV> face_centers;
    //ARRAY<TV> new_face_centers;
    ARRAY<T> face_areas;
    T volume;
    T new_volume;
    ARRAY<int> Neumann_face_indices;
    ARRAY<TV> new_face_normal;
    ARRAY<T> new_face_area;
    ARRAY<TV> new_face_center;

    // fast marching
    ARRAY<VECTOR<int,3> > tet_list;
    ARRAY<int> fm_neighbor_list;
    ARRAY<int> code;


    OCTREE_NODE()
    {children.Resize(8);children.Fill(0);}

    OCTREE_NODE(const TV& center_input,const TV& dx_over_two_input)
        :center(center_input),dx_over_two(dx_over_two_input)
    {children.Resize(8);children.Fill(0);}

    OCTREE_NODE(const TV& center_input,const TV& dx_over_two_input,const ARRAY<int>& children_input)
        :center(center_input),dx_over_two(dx_over_two_input)
    {
        assert(children_input.Size()==8);
        children.Resize(8);
        for(int i=1;i<=8;++i) children(i)=children_input(i);
    }
};

template<class T>
class OCTREE
{
    typedef VECTOR<T,3> TV;
  public:
    RANGE<TV> domain;
    bool uniform;
    ARRAY<OCTREE_NODE<T>*> nodes;
    bool READ_VORONOI_FROM_FILE;
    T distance;

    int number_of_leaves_including_ghosts;

    // 0 0 0 index or 
    // index1 index2 index3 index4
    ARRAY<ARRAY<VECTOR<int,4> > > neighbor_indices; 
    ARRAY<bool> ghost;
    HASHTABLE<int,int> tree_index_2_array_index;
    ARRAY<int> array_index_2_tree_index;
    HASHTABLE<int,T> Dirichlet_boundary;

    VORONOI_DIAGRAM<TV> *voronoi_diagram;
    ARRAY<int> voronoi_index_2_tree_index;
    ARRAY<T> dxs;

    OCTREE() {}
    //~OCTREE();

    void Clean_Memory()
    {nodes.Clean_Memory();}

//#####################################################################
    void Compute_Two_Level_Subdivision(const RANGE<TV>& domain,const int scale,bool uniform);
    void Compute_A_Simple_Subdivision(const RANGE<TV>& domain,bool fine);
    void Compute_Neighbors();
    //virtual void Compute_Volume_And_Area_For_Pressure(){}
    void Update_Volume_And_Area_For_Neumann();
    T Phi_Circle(TV position);
    int Neighbor_Same_Level(const int index,const int n_of_6);
    int Neighbor_Same_Level(const int index,const int n_of_6,int& increment);
    VECTOR<int,4> Neighbor_Leaf(const int index, const int n_of_6);
    int Tree_Index_2_Array_Index(int tree_index);
    int Array_Index_2_Tree_Index(int array_index);
    int Corner_Neighbor_Leaf_Coarse(int index,VECTOR<int,2> plane,bool& coarse);
    void Corner_Neighbor_Leaf_Fine(const int index,const int n_of_6,ARRAY<int>& corner_indices,ARRAY<int>& corner_codes); // can be 1 or 2 corner neighbors
    void All_26_Neighbors(const int index,ARRAY<int>& all_neighbor_indices,ARRAY<int>& all_neighbor_codes,bool exclude_vertex_neighbor=false); // can be 1 or 2 corner neighbors
    int Check_Neighbor_Type(const int index,const int index_neighbor,int& direction);
//#####################################################################
    void Project(ARRAY<T>& x);
    void Project_Uniform(ARRAY<T>& x);
    void Project_Non_Uniform(ARRAY<T>& x);
    void Multiply(const ARRAY<T>& x,ARRAY<T>& b);
    void Multiply_Uniform(const ARRAY<T>& x,ARRAY<T>& b);
    void Multiply_Non_Uniform(const ARRAY<T>& x,ARRAY<T>& b);
//#####################################################################
    int S_FATHER(const int index);
    int S_SON(const int index,const int octant);
    int S_SONTYPE(const int index);
    bool S_ADJ(const int direction,const int octant);
    int S_REFLECT(const int direction,const int octant);
    int S_COMMON_FACE(const int direction,const int octant);
    int S_COMMON_EDGE(const int direction,const int octant);
    int S_FACE_NEIGHBOR(const int index,const int face_direction);
    int S_EDGE_NEIGHBOR(const int index,const int edge_direction);
    int S_VERTEX_NEIGHBOR(const int index,const int vertex_direction);
};
}
#endif
