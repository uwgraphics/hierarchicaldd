//#####################################################################
// Copyright 2013, Matthew Cong.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Geometry/Spatial_Acceleration/SPHERE_HIERARCHY.h>
#include <PhysBAM_Tools/Matrices/FRAME.h>
using namespace PhysBAM;
//#####################################################################
// Function Intersection_List
//#####################################################################
template<class TV> template<class T_VISITOR> void SPHERE_HIERARCHY<TV>::
Intersection_List(const FRAME<TV>& frame,const FRAME<TV>& other_frame,const SPHERE_HIERARCHY<TV>& other_hierarchy,T_VISITOR& visitor,const int self_sphere,const int other_sphere,const T extra_thickness) const
{
    STACK<VECTOR<int,2> > stack;stack.Push(VECTOR<int,2>(self_sphere,other_sphere));

    const int self_leaves=leaves,other_leaves=other_hierarchy.leaves;
    ARRAY_VIEW<const SPHERE<TV> > self_sphere_hierarchy(sphere_hierarchy),other_sphere_hierarchy(other_hierarchy.sphere_hierarchy);
    ARRAY_VIEW<const VECTOR<int,2> > self_children(children),other_children(other_hierarchy.children);

    FRAME<TV> transform=frame.Inverse()*other_frame;
    while(!stack.Empty()){
        int self_hierarchy_sphere,other_hierarchy_sphere;stack.Pop().Get(self_hierarchy_sphere,other_hierarchy_sphere);
        if(visitor.Cull(self_hierarchy_sphere,other_hierarchy_sphere)) continue;
        SPHERE<TV> other_sphere(other_sphere_hierarchy(other_hierarchy_sphere).center+transform.X(),other_sphere_hierarchy(other_hierarchy_sphere).radius);
        if(!self_sphere_hierarchy(self_hierarchy_sphere).Intersection(other_sphere,extra_thickness)) continue;
        if(self_hierarchy_sphere<=self_leaves && other_hierarchy_sphere<=other_leaves) visitor.Store(self_hierarchy_sphere,other_hierarchy_sphere);
        else if(self_hierarchy_sphere<=self_leaves){int child_sphere1,child_sphere2;other_children(other_hierarchy_sphere-other_leaves).Get(child_sphere1,child_sphere2);
            stack.Push(VECTOR<int,2>(self_hierarchy_sphere,child_sphere1));stack.Push(VECTOR<int,2>(self_hierarchy_sphere,child_sphere2));}
        else if(other_hierarchy_sphere<=other_leaves){int child_sphere1,child_sphere2;self_children(self_hierarchy_sphere-self_leaves).Get(child_sphere1,child_sphere2);
            stack.Push(VECTOR<int,2>(child_sphere1,other_hierarchy_sphere));stack.Push(VECTOR<int,2>(child_sphere2,other_hierarchy_sphere));}
        else {int self_child_sphere1,self_child_sphere2,other_child_sphere1,other_child_sphere2;
            self_children(self_hierarchy_sphere-self_leaves).Get(self_child_sphere1,self_child_sphere2);other_children(other_hierarchy_sphere-other_leaves).Get(other_child_sphere1,other_child_sphere2);
            stack.Push(VECTOR<int,2>(self_child_sphere1,other_child_sphere1));stack.Push(VECTOR<int,2>(self_child_sphere1,other_child_sphere2));
            stack.Push(VECTOR<int,2>(self_child_sphere2,other_child_sphere1));stack.Push(VECTOR<int,2>(self_child_sphere2,other_child_sphere2));}}
}
//#####################################################################
// Function Intersection_List
//#####################################################################
template<class TV> template<class T_VISITOR> void SPHERE_HIERARCHY<TV>::
Swept_Intersection_List(const VECTOR<FRAME<TV>,2>& frame,const VECTOR<FRAME<TV>,2>& other_frame,const SPHERE_HIERARCHY<TV>& other_hierarchy,T_VISITOR& visitor,const int self_sphere,const int other_sphere,const T extra_thickness) const
{
    STACK<VECTOR<int,2> > stack;stack.Push(VECTOR<int,2>(self_sphere,other_sphere));

    // make some fields local
    const int self_leaves=leaves,other_leaves=other_hierarchy.leaves;
    ARRAY_VIEW<const SPHERE<TV> > self_sphere_hierarchy(sphere_hierarchy),other_sphere_hierarchy(other_hierarchy.sphere_hierarchy);
    ARRAY_VIEW<const VECTOR<int,2> > self_children(children),other_children(other_hierarchy.children);

    while(!stack.Empty()){
        int self_hierarchy_sphere,other_hierarchy_sphere;stack.Pop().Get(self_hierarchy_sphere,other_hierarchy_sphere);
        if(visitor.Cull(self_hierarchy_sphere,other_hierarchy_sphere)) continue;
        SPHERE<TV> other_sphere1=SPHERE<TV>(other_sphere_hierarchy(other_hierarchy_sphere),other_frame(1));
        SPHERE<TV> other_sphere2=SPHERE<TV>(other_sphere_hierarchy(other_hierarchy_sphere),other_frame(2));
        SPHERE<TV> other_sphere=SPHERE<TV>::Combine(other_sphere1,other_sphere2);
        SPHERE<TV> self_sphere1=SPHERE<TV>(self_sphere_hierarchy(self_hierarchy_sphere),frame(1));
        SPHERE<TV> self_sphere2=SPHERE<TV>(self_sphere_hierarchy(self_hierarchy_sphere),frame(2));
        SPHERE<TV> self_sphere=SPHERE<TV>::Combine(self_sphere1,self_sphere2);
        if(!self_sphere.Intersection(other_sphere,extra_thickness)) continue;
        if(self_hierarchy_sphere<=self_leaves && other_hierarchy_sphere<=other_leaves) visitor.Store(self_hierarchy_sphere,other_hierarchy_sphere);
        else if(self_hierarchy_sphere<=self_leaves){int child_sphere1,child_sphere2;other_children(other_hierarchy_sphere-other_leaves).Get(child_sphere1,child_sphere2);
            stack.Push(VECTOR<int,2>(self_hierarchy_sphere,child_sphere1));stack.Push(VECTOR<int,2>(self_hierarchy_sphere,child_sphere2));}
        else if(other_hierarchy_sphere<=other_leaves){int child_sphere1,child_sphere2;self_children(self_hierarchy_sphere-self_leaves).Get(child_sphere1,child_sphere2);
            stack.Push(VECTOR<int,2>(child_sphere1,other_hierarchy_sphere));stack.Push(VECTOR<int,2>(child_sphere2,other_hierarchy_sphere));}
        else {int self_child_sphere1,self_child_sphere2,other_child_sphere1,other_child_sphere2;
            self_children(self_hierarchy_sphere-self_leaves).Get(self_child_sphere1,self_child_sphere2);other_children(other_hierarchy_sphere-other_leaves).Get(other_child_sphere1,other_child_sphere2);
            stack.Push(VECTOR<int,2>(self_child_sphere1,other_child_sphere1));stack.Push(VECTOR<int,2>(self_child_sphere1,other_child_sphere2));
            stack.Push(VECTOR<int,2>(self_child_sphere2,other_child_sphere1));stack.Push(VECTOR<int,2>(self_child_sphere2,other_child_sphere2));}}
}
//#####################################################################
// Function Intersection_List
//#####################################################################
template<class TV> template<class T_VISITOR> void SPHERE_HIERARCHY<TV>::
Intersection_List(const SPHERE_HIERARCHY<TV>& other_hierarchy,T_VISITOR& visitor,const int self_sphere,const int other_sphere,const T extra_thickness) const
{
    STACK<VECTOR<int,2> > stack;stack.Push(VECTOR<int,2>(self_sphere,other_sphere));

    // make some fields local
    const int self_leaves=leaves,other_leaves=other_hierarchy.leaves;
    ARRAY_VIEW<const SPHERE<TV> > self_sphere_hierarchy(sphere_hierarchy),other_sphere_hierarchy(other_hierarchy.sphere_hierarchy);
    ARRAY_VIEW<const VECTOR<int,2> > self_children(children),other_children(other_hierarchy.children);

    while(!stack.Empty()){
        int self_hierarchy_sphere,other_hierarchy_sphere;stack.Pop().Get(self_hierarchy_sphere,other_hierarchy_sphere);
        if(visitor.Cull(self_hierarchy_sphere,other_hierarchy_sphere)) continue;
        if(!self_sphere_hierarchy(self_hierarchy_sphere).Intersection(other_sphere_hierarchy(other_hierarchy_sphere),extra_thickness)) continue;
        if(self_hierarchy_sphere<=self_leaves && other_hierarchy_sphere<=other_leaves) visitor.Store(self_hierarchy_sphere,other_hierarchy_sphere);
        else if(self_hierarchy_sphere<=self_leaves){int child_sphere1,child_sphere2;other_children(other_hierarchy_sphere-other_leaves).Get(child_sphere1,child_sphere2);
            stack.Push(VECTOR<int,2>(self_hierarchy_sphere,child_sphere1));stack.Push(VECTOR<int,2>(self_hierarchy_sphere,child_sphere2));}
        else if(other_hierarchy_sphere<=other_leaves){int child_sphere1,child_sphere2;self_children(self_hierarchy_sphere-self_leaves).Get(child_sphere1,child_sphere2);
            stack.Push(VECTOR<int,2>(child_sphere1,other_hierarchy_sphere));stack.Push(VECTOR<int,2>(child_sphere2,other_hierarchy_sphere));}
        else {int self_child_sphere1,self_child_sphere2,other_child_sphere1,other_child_sphere2;
            self_children(self_hierarchy_sphere-self_leaves).Get(self_child_sphere1,self_child_sphere2);other_children(other_hierarchy_sphere-other_leaves).Get(other_child_sphere1,other_child_sphere2);
            stack.Push(VECTOR<int,2>(self_child_sphere1,other_child_sphere1));stack.Push(VECTOR<int,2>(self_child_sphere1,other_child_sphere2));
            stack.Push(VECTOR<int,2>(self_child_sphere2,other_child_sphere1));stack.Push(VECTOR<int,2>(self_child_sphere2,other_child_sphere2));}}
}
//#####################################################################
// Function Intersection_List
//#####################################################################
template<class TV> template<class T_VISITOR> void SPHERE_HIERARCHY<TV>::
Intersection_List(const FRAME<TV>& frame,const FRAME<TV>& other_frame,const SPHERE_HIERARCHY<TV>& other_hierarchy,T_VISITOR& visitor,const T extra_thickness) const
{
    if(&other_hierarchy!=this || frame!=other_frame)
        return Intersection_List(frame,other_frame,other_hierarchy,visitor,root,other_hierarchy.root,extra_thickness);

    STACK<int> stack;stack.Push(root);
    while(!stack.Empty()){
        int sphere=stack.Pop();
        if(sphere>leaves && !visitor.Cull_Self(sphere)){
            int child_sphere1,child_sphere2;children(sphere-leaves).Get(child_sphere1,child_sphere2);
            stack.Push(child_sphere1);stack.Push(child_sphere2);
            Intersection_List(*this,visitor,child_sphere1,child_sphere2,extra_thickness);}}
}
//#####################################################################
// Function Intersection_List
//#####################################################################
template<class TV> template<class T_VISITOR> void SPHERE_HIERARCHY<TV>::
Swept_Intersection_List(const VECTOR<FRAME<TV>,2>& frame,const VECTOR<FRAME<TV>,2>& other_frame,const SPHERE_HIERARCHY<TV>& other_hierarchy,T_VISITOR& visitor,const T extra_thickness) const
{
    return Swept_Intersection_List(frame,other_frame,other_hierarchy,visitor,root,other_hierarchy.root,extra_thickness);
}
//#####################################################################
// Function Intersection_List
//#####################################################################
template<class TV> template<class T_VISITOR> void SPHERE_HIERARCHY<TV>::
Intersection_List(const SPHERE_HIERARCHY<TV>& other_hierarchy,T_VISITOR& visitor,const T extra_thickness) const
{
    Intersection_List(FRAME<TV>(),FRAME<TV>(),other_hierarchy,visitor,extra_thickness);
}
//#####################################################################
// Function Intersection_List
//#####################################################################
template<class TV> template<class T_VISITOR> void SPHERE_HIERARCHY<TV>::
Intersection_List(T_VISITOR& visitor) const
{
    STACK<int> stack;stack.Push(root);
    while(!stack.Empty()){int sphere=stack.Pop();
        if(visitor.Cull(sphere)) continue;
        if(Leaf(sphere)) visitor.Store(sphere);
        else{
            int child_sphere1,child_sphere2;children(sphere-leaves).Get(child_sphere1,child_sphere2);
            stack.Push(child_sphere1);stack.Push(child_sphere2);}}
}
//#####################################################################

