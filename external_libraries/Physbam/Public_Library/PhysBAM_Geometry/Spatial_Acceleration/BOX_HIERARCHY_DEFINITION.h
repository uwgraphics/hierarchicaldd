//#####################################################################
// Copyright 2004-2011, Zhaosheng Bao, Ron Fedkiw, Jon Gretarsson, Geoffrey Irving, Michael Lentine, Andrew Selle, Joseph Teran, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Geometry/Basic_Geometry/ORIENTED_BOX.h>
#include <PhysBAM_Geometry/Spatial_Acceleration/BOX_HIERARCHY.h>
#include <PhysBAM_Tools/Matrices/FRAME.h>
using namespace PhysBAM;
//#####################################################################
// Function Intersection_List
//#####################################################################
template<class TV> template<class T_VISITOR,class T_THICKNESS> void BOX_HIERARCHY<TV>::
Intersection_List(const FRAME<TV>& frame,const FRAME<TV>& other_frame,const BOX_HIERARCHY<TV>& other_hierarchy,T_VISITOR& visitor,const int self_box,const int other_box,const T_THICKNESS extra_thickness) const
{
    assert((IS_SAME<T_THICKNESS,ZERO>::value || !!extra_thickness));

    // borrow stack ownership to improve aliasing semantics
    STACK<VECTOR<int,2> > stack;stack.Exchange(dual_traversal_stack); // borrow stack ownership to improve aliasing semantics
    assert(stack.Empty());stack.Push(VECTOR<int,2>(self_box,other_box));

    // make some fields local
    const int self_leaves=leaves,other_leaves=other_hierarchy.leaves;
    ARRAY_VIEW<const RANGE<TV> > self_box_hierarchy(box_hierarchy),other_box_hierarchy(other_hierarchy.box_hierarchy);
    ARRAY_VIEW<const VECTOR<int,2> > self_children(children),other_children(other_hierarchy.children);

    FRAME<TV> transform=frame.Inverse()*other_frame;
    while(!stack.Empty()){
        int self_hierarchy_box,other_hierarchy_box;stack.Pop().Get(self_hierarchy_box,other_hierarchy_box);
        if(visitor.Cull(self_hierarchy_box,other_hierarchy_box)) continue;
        RANGE<TV> other_box=ORIENTED_BOX<TV>(other_box_hierarchy(other_hierarchy_box),transform).Axis_Aligned_Bounding_Box();
        if(!self_box_hierarchy(self_hierarchy_box).Intersection(other_box,extra_thickness)) continue;
        if(self_hierarchy_box<=self_leaves && other_hierarchy_box<=other_leaves) visitor.Store(self_hierarchy_box,other_hierarchy_box);
        else if(self_hierarchy_box<=self_leaves){int child_box1,child_box2;other_children(other_hierarchy_box-other_leaves).Get(child_box1,child_box2);
            stack.Push(VECTOR<int,2>(self_hierarchy_box,child_box1));stack.Push(VECTOR<int,2>(self_hierarchy_box,child_box2));}
        else if(other_hierarchy_box<=other_leaves){int child_box1,child_box2;self_children(self_hierarchy_box-self_leaves).Get(child_box1,child_box2);
            stack.Push(VECTOR<int,2>(child_box1,other_hierarchy_box));stack.Push(VECTOR<int,2>(child_box2,other_hierarchy_box));}
        else {int self_child_box1,self_child_box2,other_child_box1,other_child_box2;
            self_children(self_hierarchy_box-self_leaves).Get(self_child_box1,self_child_box2);other_children(other_hierarchy_box-other_leaves).Get(other_child_box1,other_child_box2);
            stack.Push(VECTOR<int,2>(self_child_box1,other_child_box1));stack.Push(VECTOR<int,2>(self_child_box1,other_child_box2));
            stack.Push(VECTOR<int,2>(self_child_box2,other_child_box1));stack.Push(VECTOR<int,2>(self_child_box2,other_child_box2));}}
    stack.Exchange(dual_traversal_stack); // release stack ownership
}
//#####################################################################
// Function Intersection_List
//#####################################################################
template<class TV> template<class T_VISITOR,class T_THICKNESS> void BOX_HIERARCHY<TV>::
Swept_Intersection_List(const VECTOR<FRAME<TV>,2>& frame,const VECTOR<FRAME<TV>,2>& other_frame,const BOX_HIERARCHY<TV>& other_hierarchy,T_VISITOR& visitor,const int self_box,const int other_box,const T_THICKNESS extra_thickness) const
{
    assert((IS_SAME<T_THICKNESS,ZERO>::value || !!extra_thickness));

    // borrow stack ownership to improve aliasing semantics
    STACK<VECTOR<int,2> > stack;stack.Exchange(dual_traversal_stack); // borrow stack ownership to improve aliasing semantics
    assert(stack.Empty());stack.Push(VECTOR<int,2>(self_box,other_box));

    // make some fields local
    const int self_leaves=leaves,other_leaves=other_hierarchy.leaves;
    ARRAY_VIEW<const RANGE<TV> > self_box_hierarchy(box_hierarchy),other_box_hierarchy(other_hierarchy.box_hierarchy);
    ARRAY_VIEW<const VECTOR<int,2> > self_children(children),other_children(other_hierarchy.children);

    while(!stack.Empty()){
        int self_hierarchy_box,other_hierarchy_box;stack.Pop().Get(self_hierarchy_box,other_hierarchy_box);
        if(visitor.Cull(self_hierarchy_box,other_hierarchy_box)) continue;
        RANGE<TV> other_box1=ORIENTED_BOX<TV>(other_box_hierarchy(other_hierarchy_box),other_frame(1)).Axis_Aligned_Bounding_Box();
        RANGE<TV> other_box2=ORIENTED_BOX<TV>(other_box_hierarchy(other_hierarchy_box),other_frame(2)).Axis_Aligned_Bounding_Box();
        RANGE<TV> other_box=RANGE<TV>::Combine(other_box1,other_box2);
        RANGE<TV> self_box1=ORIENTED_BOX<TV>(self_box_hierarchy(self_hierarchy_box),frame(1)).Axis_Aligned_Bounding_Box();
        RANGE<TV> self_box2=ORIENTED_BOX<TV>(self_box_hierarchy(self_hierarchy_box),frame(2)).Axis_Aligned_Bounding_Box();
        RANGE<TV> self_box=RANGE<TV>::Combine(self_box1,self_box2);
        if(!self_box.Intersection(other_box,extra_thickness)) continue;
        if(self_hierarchy_box<=self_leaves && other_hierarchy_box<=other_leaves) visitor.Store(self_hierarchy_box,other_hierarchy_box);
        else if(self_hierarchy_box<=self_leaves){int child_box1,child_box2;other_children(other_hierarchy_box-other_leaves).Get(child_box1,child_box2);
            stack.Push(VECTOR<int,2>(self_hierarchy_box,child_box1));stack.Push(VECTOR<int,2>(self_hierarchy_box,child_box2));}
        else if(other_hierarchy_box<=other_leaves){int child_box1,child_box2;self_children(self_hierarchy_box-self_leaves).Get(child_box1,child_box2);
            stack.Push(VECTOR<int,2>(child_box1,other_hierarchy_box));stack.Push(VECTOR<int,2>(child_box2,other_hierarchy_box));}
        else {int self_child_box1,self_child_box2,other_child_box1,other_child_box2;
            self_children(self_hierarchy_box-self_leaves).Get(self_child_box1,self_child_box2);other_children(other_hierarchy_box-other_leaves).Get(other_child_box1,other_child_box2);
            stack.Push(VECTOR<int,2>(self_child_box1,other_child_box1));stack.Push(VECTOR<int,2>(self_child_box1,other_child_box2));
            stack.Push(VECTOR<int,2>(self_child_box2,other_child_box1));stack.Push(VECTOR<int,2>(self_child_box2,other_child_box2));}}
    stack.Exchange(dual_traversal_stack); // release stack ownership
}

//#####################################################################
// Function Intersection_List
//#####################################################################
template<class TV> template<class T_VISITOR,class T_THICKNESS> void BOX_HIERARCHY<TV>::
Intersection_List(const BOX_HIERARCHY<TV>& other_hierarchy,T_VISITOR& visitor,const int self_box,const int other_box,const T_THICKNESS extra_thickness) const
{
    assert((IS_SAME<T_THICKNESS,ZERO>::value || !!extra_thickness));

    // borrow stack ownership to improve aliasing semantics
    STACK<VECTOR<int,2> > stack;stack.Exchange(dual_traversal_stack); // borrow stack ownership to improve aliasing semantics
    assert(stack.Empty());stack.Push(VECTOR<int,2>(self_box,other_box));

    // make some fields local
    const int self_leaves=leaves,other_leaves=other_hierarchy.leaves;
    ARRAY_VIEW<const RANGE<TV> > self_box_hierarchy(box_hierarchy),other_box_hierarchy(other_hierarchy.box_hierarchy);
    ARRAY_VIEW<const VECTOR<int,2> > self_children(children),other_children(other_hierarchy.children);

    while(!stack.Empty()){
        int self_hierarchy_box,other_hierarchy_box;stack.Pop().Get(self_hierarchy_box,other_hierarchy_box);
        if(visitor.Cull(self_hierarchy_box,other_hierarchy_box)) continue;
        if(!self_box_hierarchy(self_hierarchy_box).Intersection(other_box_hierarchy(other_hierarchy_box),extra_thickness)) continue;
        if(self_hierarchy_box<=self_leaves && other_hierarchy_box<=other_leaves) visitor.Store(self_hierarchy_box,other_hierarchy_box);
        else if(self_hierarchy_box<=self_leaves){int child_box1,child_box2;other_children(other_hierarchy_box-other_leaves).Get(child_box1,child_box2);
            stack.Push(VECTOR<int,2>(self_hierarchy_box,child_box1));stack.Push(VECTOR<int,2>(self_hierarchy_box,child_box2));}
        else if(other_hierarchy_box<=other_leaves){int child_box1,child_box2;self_children(self_hierarchy_box-self_leaves).Get(child_box1,child_box2);
            stack.Push(VECTOR<int,2>(child_box1,other_hierarchy_box));stack.Push(VECTOR<int,2>(child_box2,other_hierarchy_box));}
        else {int self_child_box1,self_child_box2,other_child_box1,other_child_box2;
            self_children(self_hierarchy_box-self_leaves).Get(self_child_box1,self_child_box2);other_children(other_hierarchy_box-other_leaves).Get(other_child_box1,other_child_box2);
            stack.Push(VECTOR<int,2>(self_child_box1,other_child_box1));stack.Push(VECTOR<int,2>(self_child_box1,other_child_box2));
            stack.Push(VECTOR<int,2>(self_child_box2,other_child_box1));stack.Push(VECTOR<int,2>(self_child_box2,other_child_box2));}}
    stack.Exchange(dual_traversal_stack); // release stack ownership
}
//#####################################################################
// Function Intersection_List
//#####################################################################
template<class TV> template<class T_VISITOR,class T_THICKNESS> void BOX_HIERARCHY<TV>::
Intersection_List(const FRAME<TV>& frame,const FRAME<TV>& other_frame,const BOX_HIERARCHY<TV>& other_hierarchy,T_VISITOR& visitor,const T_THICKNESS extra_thickness) const
{
    if(!IS_SAME<T_THICKNESS,ZERO>::value && !extra_thickness){Intersection_List(frame,other_frame,other_hierarchy,visitor,ZERO());return;}

    if(frame != other_frame) return Intersection_List(frame,other_frame,other_hierarchy,visitor,root,other_hierarchy.root,extra_thickness);
    else if(&other_hierarchy != this) return Intersection_List(other_hierarchy,visitor,root,other_hierarchy.root,extra_thickness);

    // borrow stack ownership to improve aliasing semantics
    STACK<int> stack;stack.Exchange(traversal_stack);
    assert(stack.Empty());stack.Push(root);

    // make some fields local
    const int self_leaves=leaves;
    ARRAY_VIEW<const RANGE<TV> > self_box_hierarchy(box_hierarchy);
    ARRAY_VIEW<const VECTOR<int,2> > self_children(children);

    while(!stack.Empty()){
        int box=stack.Pop();
        if(box>self_leaves && !visitor.Cull_Self(box)){
            int child_box1,child_box2;self_children(box-self_leaves).Get(child_box1,child_box2);
            stack.Push(child_box1);stack.Push(child_box2);
            Intersection_List(*this,visitor,child_box1,child_box2,extra_thickness);}}
    stack.Exchange(traversal_stack); // release stack ownership
}
//#####################################################################
// Function Intersection_List
//#####################################################################
template<class TV> template<class T_VISITOR,class T_THICKNESS> void BOX_HIERARCHY<TV>::
Swept_Intersection_List(const VECTOR<FRAME<TV>,2>& frame,const VECTOR<FRAME<TV>,2>& other_frame,const BOX_HIERARCHY<TV>& other_hierarchy,T_VISITOR& visitor,const T_THICKNESS extra_thickness) const
{
    if(!IS_SAME<T_THICKNESS,ZERO>::value && !extra_thickness){Swept_Intersection_List(frame,other_frame,other_hierarchy,visitor,ZERO());return;}
    return Swept_Intersection_List(frame,other_frame,other_hierarchy,visitor,root,other_hierarchy.root,extra_thickness);
}
//#####################################################################
// Function Intersection_List
//#####################################################################
template<class TV> template<class T_VISITOR,class T_THICKNESS> void BOX_HIERARCHY<TV>::
Intersection_List(const BOX_HIERARCHY<TV>& other_hierarchy,T_VISITOR& visitor,const T_THICKNESS extra_thickness) const
{
    Intersection_List(FRAME<TV>(),FRAME<TV>(),other_hierarchy,visitor,extra_thickness);
}
//#####################################################################
// Function Intersection_List
//#####################################################################
template<class TV> template<class T_VISITOR> void BOX_HIERARCHY<TV>::
Intersection_List(T_VISITOR& visitor) const
{
    STACK<int> stack;stack.Exchange(traversal_stack); // borrow stack ownership to improve aliasing semantics
    stack.Remove_All();stack.Push(root);
    while(!stack.Empty()){int box=stack.Pop();
        if(visitor.Cull(box)) continue;
        if(Leaf(box)) visitor.Store(box);
        else{
            int child_box1,child_box2;children(box-leaves).Get(child_box1,child_box2);
            stack.Push(child_box1);stack.Push(child_box2);}}
    stack.Exchange(traversal_stack); // release stack ownership
}
//#####################################################################

