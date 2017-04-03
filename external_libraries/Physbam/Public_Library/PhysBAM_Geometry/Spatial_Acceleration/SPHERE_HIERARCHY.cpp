//#####################################################################
// Copyright 2013, Matthew Cong.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Arrays/INDIRECT_ARRAY.h>
#include <PhysBAM_Tools/Data_Structures/KD_TREE.h>
#include <PhysBAM_Tools/Math_Tools/RANGE.h>
#include <PhysBAM_Geometry/Basic_Geometry/SPHERE.h>
#include <PhysBAM_Geometry/Spatial_Acceleration/SPHERE_HIERARCHY.h>
#include <PhysBAM_Geometry/Spatial_Acceleration/SPHERE_HIERARCHY_DEFINITION.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> SPHERE_HIERARCHY<TV>::
SPHERE_HIERARCHY()
{}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> SPHERE_HIERARCHY<TV>::
~SPHERE_HIERARCHY()
{}
//#####################################################################
// Function Clean_Memory
//#####################################################################
template<class TV> void SPHERE_HIERARCHY<TV>::
Clean_Memory()
{
    leaves=root=0;
    parents.Clean_Memory();children.Clean_Memory();sphere_hierarchy.Clean_Memory();
}
//#####################################################################
// Function Set_Leaf_Spheres
//#####################################################################
template<class TV> void SPHERE_HIERARCHY<TV>::
Set_Leaf_Spheres(const ARRAY<SPHERE<TV> >& spheres,const bool reinitialize)
{
    if(reinitialize){sphere_hierarchy=spheres;Initialize_Hierarchy_Using_KD_Tree();}
    else sphere_hierarchy.Prefix(leaves)=spheres;
    Update_Nonleaf_Spheres();
}
//#####################################################################
// Function Thicken_Leaf_Spheres
//#####################################################################
template<class TV> void SPHERE_HIERARCHY<TV>::
Thicken_Leaf_Spheres(const T extra_thickness)
{
    for(int k=1;k<=leaves;k++) sphere_hierarchy(k).radius+=extra_thickness;
}
//#####################################################################
// Function Initialize_Hierarchy_Using_KD_Tree
//#####################################################################
template<class TV> void SPHERE_HIERARCHY<TV>::
Initialize_Hierarchy_Using_KD_Tree()
{
    KD_TREE<TV> kd_tree(false);
    ARRAY<TV> centroids(sphere_hierarchy.m);
    for(int l=1;l<=sphere_hierarchy.m;l++)centroids(l)=sphere_hierarchy(l).center;
    kd_tree.Create_Left_Balanced_KD_Tree(centroids);
    leaves=sphere_hierarchy.m;parents.Resize(leaves);children.Remove_All();
    root=Initialize_Hierarchy_Using_KD_Tree_Helper(kd_tree.root_node);assert(root==2*leaves-1);
    sphere_hierarchy.Resize(root);
}
//#####################################################################
// Function Initialize_Hierarchy_Using_KD_Tree_Helper
//#####################################################################
template<class TV> int SPHERE_HIERARCHY<TV>::
Initialize_Hierarchy_Using_KD_Tree_Helper(KD_TREE_NODE<T>* node)
{
    if(!node->left&&!node->right)return node->node_index;
    int left_child=Initialize_Hierarchy_Using_KD_Tree_Helper(node->left);
    if(!node->right)return left_child;
    int right_child=Initialize_Hierarchy_Using_KD_Tree_Helper(node->right);
    children.Append(VECTOR<int,2>(left_child,right_child));parents.Append(0);
    return parents(left_child)=parents(right_child)=children.m+leaves;
}
//#####################################################################
// Function Update_Nonleaf_Spheres
//#####################################################################
template<class TV> void SPHERE_HIERARCHY<TV>::
Update_Nonleaf_Spheres()
{
    for(int k=leaves+1;k<=sphere_hierarchy.m;k++) sphere_hierarchy(k)=SPHERE<TV>::Combine(sphere_hierarchy(children(k-leaves)(1)),sphere_hierarchy(children(k-leaves)(2)));
}
//#####################################################################
// Function Update_Modified_Nonleaf_Spheres
//#####################################################################
template<class TV> void SPHERE_HIERARCHY<TV>::
Update_Modified_Nonleaf_Spheres(ARRAY<bool>& modified)
{
    for(int k=leaves+1;k<=sphere_hierarchy.m;k++){
        const VECTOR<int,2>& child=children(k-leaves);
        modified(k)=modified(child[1]) || modified(child[2]);
        if(modified(k)) sphere_hierarchy(k)=SPHERE<TV>::Combine(sphere_hierarchy(child[1]),sphere_hierarchy(child[2]));}
}
//#####################################################################
#define INSTANTIATION_HELPER(T,d) \
    template class SPHERE_HIERARCHY<VECTOR<T,d> >;

INSTANTIATION_HELPER(float,1)
INSTANTIATION_HELPER(float,2)
INSTANTIATION_HELPER(float,3)
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
INSTANTIATION_HELPER(double,1)
INSTANTIATION_HELPER(double,2)
INSTANTIATION_HELPER(double,3)
#endif

