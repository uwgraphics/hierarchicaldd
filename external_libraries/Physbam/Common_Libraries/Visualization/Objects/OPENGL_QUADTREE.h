//#####################################################################
// Copyright 2015, Mridul Aanjaneya, Ming Gao.
// This file is part of PhysBAM whose distribution is governed by the license 
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_QUADTREE
//#####################################################################
#ifndef __OPENGL_QUADTREE__
#define __OPENGL_QUADTREE__

#include <Common_Tools/Data_Structures/QUADTREE.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_COLOR.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_OBJECT.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_SELECTION.h>

namespace PhysBAM{
template<class T>
class OPENGL_QUADTREE: public OPENGL_OBJECT
{
  public:
    typedef VECTOR<T,2> TV;

    QUADTREE<T> tree;
    bool draw,draw_centers;
    int current_frame;
    RANGE<VECTOR<float,2> > domain;
  private:
    OPENGL_SELECTION *current_selection;

  public:
    OPENGL_QUADTREE()
        :draw(true),draw_centers(false),current_frame(-1),current_selection(0)
    {}

//#####################################################################
    virtual RANGE<VECTOR<float,3> > Bounding_Box() const PHYSBAM_OVERRIDE;
    void Display(const int in_color=1) const PHYSBAM_OVERRIDE;
    virtual OPENGL_SELECTION *Get_Selection(GLuint *buffer, int buffer_size);
    void Highlight_Selection(OPENGL_SELECTION *selection) PHYSBAM_OVERRIDE;
    void Clear_Highlight() PHYSBAM_OVERRIDE;
    void Print_Selection_Info(std::ostream& stream,OPENGL_SELECTION* selection) const PHYSBAM_OVERRIDE;
    void Update_Domain();
//#####################################################################
};

template<class T>
class OPENGL_SELECTION_QUADTREE_CELL: public OPENGL_SELECTION
{
    typedef VECTOR<T,2> TV;
public:
    int index;
    OPENGL_SELECTION_QUADTREE_CELL(OPENGL_OBJECT *object,const int index) 
        : OPENGL_SELECTION(OPENGL_SELECTION::GRID_CELL_2D, object),index(index) {}

    RANGE<VECTOR<float,3> > Bounding_Box() const PHYSBAM_OVERRIDE;
};
}
#endif
