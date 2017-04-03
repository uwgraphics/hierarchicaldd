//#####################################################################
// Copyright 2013, Matthew Cong.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SPHERE_HIERARCHY
//##################################################################### 
#ifndef __SPHERE_HIERARCHY__
#define __SPHERE_HIERARCHY__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Data_Structures/STACK.h>
#include <PhysBAM_Tools/Utilities/NONCOPYABLE.h>
#include <PhysBAM_Tools/Utilities/PHYSBAM_OVERRIDE.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Geometry/Basic_Geometry/SPHERE.h>
namespace PhysBAM{

template<class T> class KD_TREE_NODE;

struct SPHERE_VISITOR_TRIVIAL
{
    ARRAY<ARRAY<int> >& intersection_list;

    SPHERE_VISITOR_TRIVIAL(ARRAY<ARRAY<int> >& intersection_list)
        :intersection_list(intersection_list)
    {}

    bool Cull_Self(const int self_sphere_index) const
    {return false;}

    bool Cull(const int self_sphere_index,const int other_sphere_index) const
    {return false;}

    void Store(const int self_sphere_index,const int other_sphere_index) const
    {intersection_list(self_sphere_index).Append(other_sphere_index);}
};

template<class T_NESTED_VISITOR>
struct SPHERE_VISITOR_MPI
{
    T_NESTED_VISITOR nested_visitor;
    const ARRAY<char> &processors1,&processors2; // entries are 2*has_ours + has_greater 

    SPHERE_VISITOR_MPI(T_NESTED_VISITOR nested_visitor,const ARRAY<char>& processors1,const ARRAY<char>& processors2)
        :nested_visitor(nested_visitor),processors1(processors1),processors2(processors2)
    {}

    bool Cull_Self(const int sphere) const
    {return processors1(sphere)<2 || nested_visitor.Cull_Self(sphere);} // has_ours

    bool Cull(const int sphere1,const int sphere2) const
    {return (processors1(sphere1)<<processors2(sphere2))<4 || nested_visitor.Cull(sphere1,sphere2);} // has_ours_1 and (has_ours_2 or has_greater_2) or has_ours_2 and has_greater_1

    void Store(const int self_sphere_index,const int other_sphere_index)
    {nested_visitor.Store(self_sphere_index,other_sphere_index);}
};

template<class TV>
class SPHERE_HIERARCHY:public NONCOPYABLE
{
private:
    typedef typename TV::SCALAR T;
public:
    int leaves,root;
    ARRAY<int> parents;
    ARRAY<VECTOR<int,2> > children;
    ARRAY<SPHERE<TV> > sphere_hierarchy;

    SPHERE_HIERARCHY();
    virtual ~SPHERE_HIERARCHY();

    bool Leaf(const int sphere) const
    {return sphere<=leaves;}

//#####################################################################
    virtual void Clean_Memory();
protected:
    int Initialize_Hierarchy_Using_KD_Tree_Helper(KD_TREE_NODE<T>* node);
public:
    virtual void Initialize_Hierarchy_Using_KD_Tree();
    void Set_Leaf_Spheres(const ARRAY<SPHERE<TV> >& spheres,const bool reinitialize=false);
    void Thicken_Leaf_Spheres(const T extra_thickness);
    void Update_Nonleaf_Spheres();
    void Update_Modified_Nonleaf_Spheres(ARRAY<bool>& modified);
private:
    template<class T_VISITOR> void Intersection_List(const SPHERE_HIERARCHY<TV>& other_hierarchy,T_VISITOR& visitor,const int self_sphere,const int other_sphere,const T extra_thickness=(T)0.) const;
    template<class T_VISITOR> void Intersection_List(const FRAME<TV>& frame,const FRAME<TV>& other_frame,const SPHERE_HIERARCHY<TV>& other_hierarchy,T_VISITOR& visitor,const int self_sphere,const int other_sphere,const T extra_thickness=(T)0.) const;
    template<class T_VISITOR> void Swept_Intersection_List(const VECTOR<FRAME<TV>,2>& frame,const VECTOR<FRAME<TV>,2>& other_frame,const SPHERE_HIERARCHY<TV>& other_hierarchy,T_VISITOR& visitor,const int self_sphere,const int other_sphere,const T extra_thickness=(T)0.) const;
public:
    template<class T_VISITOR> void Intersection_List(const SPHERE_HIERARCHY<TV>& other_hierarchy,T_VISITOR& visitor,const T extra_thickness=(T)0.) const;
    template<class T_VISITOR> void Intersection_List(const FRAME<TV>& frame,const FRAME<TV>& other_frame,const SPHERE_HIERARCHY<TV>& other_hierarchy,T_VISITOR& visitor,const T extra_thickness=(T)0.) const;
    template<class T_VISITOR> void Swept_Intersection_List(const VECTOR<FRAME<TV>,2>& frame,const VECTOR<FRAME<TV>,2>& other_frame,const SPHERE_HIERARCHY<TV>& other_hierarchy,T_VISITOR& visitor,const T extra_thickness=(T)0.) const;
    template<class T_VISITOR> void Intersection_List(T_VISITOR& visitor) const;
//#####################################################################
};
}
#endif

