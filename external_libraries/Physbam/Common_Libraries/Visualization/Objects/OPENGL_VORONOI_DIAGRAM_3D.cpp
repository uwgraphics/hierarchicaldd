//#####################################################################
// Copyright 2015, Mridul Aanjaneya.
// This file is part of PhysBAM whose distribution is governed by the license 
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_CONSTANT_COLOR_MAP.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_PREFERENCES.h>
#include <Visualization/Objects/OPENGL_VORONOI_DIAGRAM_3D.h>
using namespace PhysBAM;
//#####################################################################
// Display
//#####################################################################
template<class T> void OPENGL_VORONOI_DIAGRAM_3D<T>::
Display(const int in_color) const
{
    if(!draw) return;
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    Send_Transform_To_GL_Pipeline();

    glPushAttrib(GL_ENABLE_BIT);
    glDisable(GL_LIGHTING);
    glDisable(GL_DEPTH_TEST);

    color.Send_To_GL_Pipeline();

    GLint mode;
    glGetIntegerv(GL_RENDER_MODE, &mode);

    if(mode == GL_SELECT){}
    else{OpenGL_Begin(GL_LINES);
        for(int i=1;i<=voronoi_diagram.face_vertices.Size();++i) for(int j=1;j<=voronoi_diagram.face_vertices(i).Size();++j)
            for(int k=1;k<=voronoi_diagram.face_vertices(i)(j).Size();++k){
                const TV& X1=voronoi_diagram.face_vertices(i)(j)(k);
                const TV& X2=(k==voronoi_diagram.face_vertices(i)(j).Size())?voronoi_diagram.face_vertices(i)(j)(1):voronoi_diagram.face_vertices(i)(j)(k+1);
                OpenGL_Line<T,3>(X1,X2);}
        OpenGL_End();}

    glPopAttrib();
    glPopMatrix();
}
//#####################################################################
// Update_Domain
//#####################################################################
template<class T> void OPENGL_VORONOI_DIAGRAM_3D<T>::
Update_Domain()
{
    domain=RANGE<VECTOR<float,3> >();
    for(int i=1;i<=voronoi_diagram.face_vertices.Size();++i) for(int j=1;j<=voronoi_diagram.face_vertices(i).Size();++j)
        for(int k=1;k<=voronoi_diagram.face_vertices(i)(j).Size();++k){const TV& X=voronoi_diagram.face_vertices(i)(j)(k);
            for(int axis=1;axis<=3;++axis){
                domain.min_corner(axis)=min(domain.min_corner(axis),(float)X(axis));
                domain.max_corner(axis)=max(domain.max_corner(axis),(float)X(axis));}}
}
//#####################################################################
// Bounding_Box
//#####################################################################
template<class T> RANGE<VECTOR<float,3> > OPENGL_VORONOI_DIAGRAM_3D<T>::
Bounding_Box() const
{
    return World_Space_Box(domain);
}
//#####################################################################
template class OPENGL_VORONOI_DIAGRAM_3D<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class OPENGL_VORONOI_DIAGRAM_3D<double>;
#endif
