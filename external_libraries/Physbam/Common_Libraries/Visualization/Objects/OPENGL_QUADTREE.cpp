//#####################################################################
// Copyright 2015, Mridul Aanjaneya, Ming Gao.
// This file is part of PhysBAM whose distribution is governed by the license 
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_PREFERENCES.h>
#include <Visualization/Objects/OPENGL_QUADTREE.h>
using namespace PhysBAM;
//#####################################################################
// Display
//#####################################################################
template<class T> void OPENGL_QUADTREE<T>::
Display(const int in_color) const
{
    if(!draw)return;
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    Send_Transform_To_GL_Pipeline();

    glPushAttrib(GL_ENABLE_BIT);
    glDisable(GL_LIGHTING);
    glDisable(GL_DEPTH_TEST);

    GLint mode;
    glGetIntegerv(GL_RENDER_MODE, &mode);

    if(mode == GL_SELECT){
        glPushAttrib(GL_ENABLE_BIT | GL_POINT_BIT);
        glDisable(GL_CULL_FACE);

        // Draw cells for selection
        glPushName(1);
        for(int i=1;i<=tree.nodes.Size();++i){QUADTREE_NODE<T>* node=tree.nodes(i);
            if(node->north_east_index==0){TV min_corner=node->center-node->dx_over_two,max_corner=node->center+node->dx_over_two;
                glPushName(i);
                OpenGL_Begin(GL_QUADS);
                OpenGL_Quad_2D(min_corner,max_corner);
                OpenGL_End();
                glPopName();}}
        glPopName();glPopAttrib();}
    else{OPENGL_COLOR color;
        if(draw_centers){color=OPENGL_COLOR::Yellow();
            color.Send_To_GL_Pipeline();
            OpenGL_Begin(GL_POINTS);
            for(int i=1;i<=tree.nodes.Size();++i) OpenGL_Vertex(tree.nodes(i)->center);
            OpenGL_End();}
    
        color=OPENGL_COLOR::Gray();
        color.Send_To_GL_Pipeline();
        OpenGL_Begin(GL_LINES);

        // domain box
        TV edge_lengths=tree.domain.Edge_Lengths();
        for(int axis=1;axis<=2;++axis){OpenGL_Line(tree.domain.min_corner,tree.domain.min_corner+TV::Axis_Vector(axis)*edge_lengths(axis));
            OpenGL_Line(tree.domain.max_corner,tree.domain.max_corner-TV::Axis_Vector(axis)*edge_lengths(axis));}

        for(int i=1;i<=tree.nodes.Size();++i){
            if(tree.nodes(i)->north_east_index>0){
            for(int axis=1;axis<=2;++axis){
                OpenGL_Line(tree.nodes(i)->center,tree.nodes(i)->center+tree.nodes(i)->dx_over_two(axis)*TV::Axis_Vector(axis));
                OpenGL_Line(tree.nodes(i)->center,tree.nodes(i)->center-tree.nodes(i)->dx_over_two(axis)*TV::Axis_Vector(axis));}}}
        OpenGL_End();
    
        // Highlight current selection
        if(current_selection){
            if(current_selection->type == OPENGL_SELECTION::GRID_CELL_2D){
                OPENGL_SELECTION_QUADTREE_CELL<T>* real_selection=(OPENGL_SELECTION_QUADTREE_CELL<T>*)current_selection;
                int index=real_selection->index;
                TV min_corner=tree.nodes(index)->center-tree.nodes(index)->dx_over_two,max_corner=tree.nodes(index)->center+tree.nodes(index)->dx_over_two;
                OPENGL_SELECTION::Draw_Highlighted_Quad(min_corner,max_corner);}}}

    glPopAttrib();
    glPopMatrix();
}

//#####################################################################
// Update_Domain
//#####################################################################
template<class T> void OPENGL_QUADTREE<T>::
Update_Domain()
{
    domain=RANGE<VECTOR<float,2> >();
    for(int i=1;i<=tree.nodes.Size();++i) for(int axis=1;axis<=2;++axis){
        domain.min_corner(axis)=min(domain.min_corner(axis),tree.nodes(i)->center(axis));
        domain.max_corner(axis)=max(domain.min_corner(axis),tree.nodes(i)->center(axis));}
}
//#####################################################################
// Bounding_Box
//#####################################################################
template<class T> RANGE<VECTOR<float,3> > OPENGL_QUADTREE<T>::
Bounding_Box() const
{
    return World_Space_Box(domain);
}
//#####################################################################
// Function Get_Selection
//#####################################################################
template<class T> OPENGL_SELECTION *OPENGL_QUADTREE<T>::
Get_Selection(GLuint* buffer,int buffer_size)
{
    OPENGL_SELECTION* selection=0;
    if(buffer_size == 2){
        if(buffer[0] == 1) selection=new OPENGL_SELECTION_QUADTREE_CELL<T>(this,buffer[1]);}
    return selection;
}
//#####################################################################
// Function Highlight_Selection
//#####################################################################
template<class T> void OPENGL_QUADTREE<T>::
Highlight_Selection(OPENGL_SELECTION *selection)
{
    delete current_selection;current_selection=0;
    if (selection->type == OPENGL_SELECTION::GRID_CELL_2D){
        OPENGL_SELECTION_QUADTREE_CELL<T>* real_selection=(OPENGL_SELECTION_QUADTREE_CELL<T>*)selection;
        current_selection=new OPENGL_SELECTION_QUADTREE_CELL<T>(this,real_selection->index);}
}
//#####################################################################
// Function Clear_Highlight
//#####################################################################
template<class T> void OPENGL_QUADTREE<T>::
Clear_Highlight()
{
    delete current_selection;current_selection=0;
}
//#####################################################################
// Print_Selection_Info
//#####################################################################
template<class T> void OPENGL_QUADTREE<T>::
Print_Selection_Info(std::ostream& stream,OPENGL_SELECTION* selection) const
{
    if(current_selection && current_selection->type==OPENGL_SELECTION::GRID_CELL_2D){
        int index=((OPENGL_SELECTION_QUADTREE_CELL<T>*)current_selection)->index;
        stream<<"Index: "<<index<<std::endl;}
}
//#####################################################################
// Function Bounding_Box
//#####################################################################
template<class T> RANGE<VECTOR<float,3> > OPENGL_SELECTION_QUADTREE_CELL<T>::
Bounding_Box() const
{
    PHYSBAM_ASSERT(object);
    TV center=((OPENGL_QUADTREE<T>*)object)->tree.nodes(index)->center,dx_over_two=((OPENGL_QUADTREE<T>*)object)->tree.nodes(index)->dx_over_two;
    TV min_corner=center-dx_over_two,max_corner=center+dx_over_two;
    RANGE<VECTOR<T,2> > box(min_corner,max_corner);
    return object->World_Space_Box(RANGE<VECTOR<float,2> >(box));
}
//#####################################################################
template class OPENGL_QUADTREE<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class OPENGL_QUADTREE<double>;
#endif
