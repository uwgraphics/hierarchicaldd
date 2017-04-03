//#####################################################################
// Copyright 2015, Mridul Aanjaneya, Ming Gao.
// This file is part of PhysBAM whose distribution is governed by the license 
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_PREFERENCES.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_UNIFORM_SLICE.h>
#include <Visualization/Objects/OPENGL_OCTREE.h>
using namespace PhysBAM;
//#####################################################################
// Display
//#####################################################################
template<class T> void OPENGL_OCTREE<T>::
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

    OPENGL_UNIFORM_SLICE* slice=(OPENGL_UNIFORM_SLICE*)this->slice;
    OPENGL_COLOR color;color=OPENGL_COLOR::Gray(0.5);
    color.Send_To_GL_Pipeline();

    if(!slice || slice->mode==OPENGL_SLICE::NO_SLICE || slice->mode == OPENGL_SLICE::NODE_SLICE){
        for(int i=1;i<=tree.nodes.Size();++i){
            if(tree.nodes(i)->children(1)==0){ARRAY<TV> nodes(8);const TV dx=(T)2.*tree.nodes(i)->dx_over_two;
                nodes(1)=tree.nodes(i)->center-tree.nodes(i)->dx_over_two;
                nodes(2)=nodes(1)+dx(1)*TV::Axis_Vector(1);
                nodes(3)=nodes(1)+dx(1)*TV::Axis_Vector(1)+dx(2)*TV::Axis_Vector(2);
                nodes(4)=nodes(1)+dx(2)*TV::Axis_Vector(2);
                nodes(5)=nodes(1)+dx(3)*TV::Axis_Vector(3);
                nodes(6)=nodes(2)+dx(3)*TV::Axis_Vector(3);
                nodes(7)=nodes(3)+dx(3)*TV::Axis_Vector(3);
                nodes(8)=nodes(4)+dx(3)*TV::Axis_Vector(3);
                Draw_Cube(nodes);}}}
    else if(slice->mode == OPENGL_SLICE::CELL_SLICE){int axis=slice->axis;
        T_INDEX slice_index=T_INDEX::All_Ones_Vector();slice_index(axis)=slice->index;
        TV location=slice->grid.Center(slice_index);

        if(mode==GL_SELECT){
            glPushAttrib(GL_ENABLE_BIT);
            glDisable(GL_CULL_FACE);
            glPushName(1);

            TV x_vector(1,0,0),y_vector(0,1,0),z_vector(0,0,1);
            TV axis_1,axis_2,axis_3;
            if(axis==1){axis_1=y_vector;axis_2=z_vector;axis_3=x_vector;}
            else if(axis==2){axis_1=z_vector;axis_2=x_vector;axis_3=y_vector;}
            else if(axis==3){axis_1=x_vector;axis_2=y_vector;axis_3=z_vector;}

            for(int i=1;i<=tree.nodes.Size();++i){
                if(tree.nodes(i)->children(1)==0){TV dx=(T)2.*tree.nodes(i)->dx_over_two;
                    if(tree.nodes(i)->center(axis)-tree.nodes(i)->dx_over_two(axis)<location(axis) && tree.nodes(i)->center(axis)+tree.nodes(i)->dx_over_two(axis)>location(axis)){
                        glPushName(i);
                        TV left_corner=tree.nodes(i)->center-tree.nodes(i)->dx_over_two;
                        for(int v=1;v<=3;++v){
                            if(axis_1(v)>(T)0.) axis_1(v)=dx(v);
                            if(axis_2(v)>(T)0.) axis_2(v)=dx(v);
                            if(axis_3(v)>(T)0.) axis_3(v)=dx(v);}
                        OpenGL_Begin(GL_QUADS);
                        OpenGL_Quad(left_corner,axis_1,axis_2);
                        OpenGL_Quad(left_corner+axis_3,axis_1,axis_2);
                        OpenGL_End();
                        glPopName();}}}
            glPopName();glPopAttrib();}
        else{for(int i=1;i<=tree.nodes.Size();++i){
                if(tree.nodes(i)->children(1)==0){ARRAY<TV> nodes(8);TV dx=(T)2.*tree.nodes(i)->dx_over_two;
                    if(tree.nodes(i)->center(axis)-tree.nodes(i)->dx_over_two(axis)<location(axis) && tree.nodes(i)->center(axis)+tree.nodes(i)->dx_over_two(axis)>location(axis)){
                        nodes(1)=tree.nodes(i)->center-tree.nodes(i)->dx_over_two;
                        nodes(2)=nodes(1)+dx(1)*TV::Axis_Vector(1);
                        nodes(3)=nodes(1)+dx(1)*TV::Axis_Vector(1)+dx(2)*TV::Axis_Vector(2);
                        nodes(4)=nodes(1)+dx(2)*TV::Axis_Vector(2);
                        nodes(5)=nodes(1)+dx(3)*TV::Axis_Vector(3);
                        nodes(6)=nodes(2)+dx(3)*TV::Axis_Vector(3);
                        nodes(7)=nodes(3)+dx(3)*TV::Axis_Vector(3);
                        nodes(8)=nodes(4)+dx(3)*TV::Axis_Vector(3);
                        Draw_Cube(nodes);}}}}}

    // Highlight current selection
    if(current_selection){
        if(current_selection->type == OPENGL_SELECTION::GRID_CELL_3D){
            OPENGL_SELECTION_OCTREE_CELL<T>* real_selection=(OPENGL_SELECTION_OCTREE_CELL<T>*)current_selection;
            int i=real_selection->index;
            TV dx=(T)2.*tree.nodes(i)->dx_over_two;
            TV min_corner=tree.nodes(i)->center-(T).5*dx,max_corner=tree.nodes(i)->center+(T).5*dx;
            OPENGL_SELECTION::Draw_Highlighted_Box(min_corner,max_corner);}}

    glPopAttrib();
    glPopMatrix();
}
//#####################################################################
// Draw_Cube
//#####################################################################
template<class T> void OPENGL_OCTREE<T>::
Draw_Cube(const ARRAY<TV>& nodes) const
{
    OpenGL_Begin(GL_LINES);
    OpenGL_Line(nodes(1),nodes(2));
    OpenGL_Line(nodes(2),nodes(3));
    OpenGL_Line(nodes(3),nodes(4));
    OpenGL_Line(nodes(4),nodes(1));
    OpenGL_Line(nodes(5),nodes(6));
    OpenGL_Line(nodes(6),nodes(7));
    OpenGL_Line(nodes(7),nodes(8));
    OpenGL_Line(nodes(8),nodes(5));
    OpenGL_Line(nodes(1),nodes(5));
    OpenGL_Line(nodes(2),nodes(6));
    OpenGL_Line(nodes(3),nodes(7));
    OpenGL_Line(nodes(4),nodes(8));
    OpenGL_End();
}
//#####################################################################
// Update_Domain
//#####################################################################
template<class T> void OPENGL_OCTREE<T>::
Update_Domain()
{
    domain=RANGE<VECTOR<float,3> >();
    for(int i=1;i<=tree.nodes.Size();++i) for(int axis=1;axis<=3;++axis){
        domain.min_corner(axis)=min(domain.min_corner(axis),tree.nodes(i)->center(axis));
        domain.max_corner(axis)=max(domain.min_corner(axis),tree.nodes(i)->center(axis));}
}
//#####################################################################
// Get_Selection
//#####################################################################
template<class T> OPENGL_SELECTION* OPENGL_OCTREE<T>::
Get_Selection(GLuint *buffer, int buffer_size)
{
    OPENGL_SELECTION* selection=0;
    if(buffer_size==2){
        if(buffer[0]==1) selection=new OPENGL_SELECTION_OCTREE_CELL<T>(this,int(buffer[1]));}
    return selection;
}
//#####################################################################
// Highlight_Selection
//#####################################################################
template<class T> void OPENGL_OCTREE<T>::
Highlight_Selection(OPENGL_SELECTION *selection)
{
    delete current_selection;current_selection=0;
    if(selection->type==OPENGL_SELECTION::GRID_CELL_3D){
        OPENGL_SELECTION_OCTREE_CELL<T>* real_selection=(OPENGL_SELECTION_OCTREE_CELL<T>*)selection;
        current_selection=new OPENGL_SELECTION_OCTREE_CELL<T>(this,real_selection->index);}
}
//#####################################################################
// Clear_Highlight
//#####################################################################
template<class T> void OPENGL_OCTREE<T>::
Clear_Highlight()
{
    delete current_selection;current_selection=0;
}
//#####################################################################
// Print_Selection_Info
//#####################################################################
template<class T> void OPENGL_OCTREE<T>::
Print_Selection_Info(std::ostream& stream,OPENGL_SELECTION* selection) const
{
    if(current_selection && current_selection->type==OPENGL_SELECTION::GRID_CELL_3D){
        int i=((OPENGL_SELECTION_OCTREE_CELL<T>*)current_selection)->index;
        stream<<"Selected cell "<<i<<" ("<<tree.nodes(i)->center<<")"<<std::endl;}
}
//#####################################################################
// Bounding_Box
//#####################################################################
template<class T> RANGE<VECTOR<float,3> > OPENGL_OCTREE<T>::
Bounding_Box() const
{
    return World_Space_Box(domain);
}
//#####################################################################
// Bounding_Box
//#####################################################################
template<class T> RANGE<VECTOR<float,3> > OPENGL_SELECTION_OCTREE_CELL<T>::
Bounding_Box() const
{
    PHYSBAM_ASSERT(object);
    const OCTREE<T>& tree=((OPENGL_OCTREE<T>*)object)->tree;
    RANGE<TV> box(tree.nodes(index)->center-tree.nodes(index)->dx_over_two,tree.nodes(index)->center+tree.nodes(index)->dx_over_two);
    return object->World_Space_Box(RANGE<VECTOR<float,3> >(box));
}
//#####################################################################
template class OPENGL_OCTREE<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class OPENGL_OCTREE<double>;
#endif
