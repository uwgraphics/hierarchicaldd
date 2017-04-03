//#####################################################################
// Copyright 2015, Mridul Aanjaneya, Ming Gao.
// This file is part of PhysBAM whose distribution is governed by the license 
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_OCTREE
//#####################################################################
#ifndef __OPENGL_OCTREE__
#define __OPENGL_OCTREE__

#include <Common_Tools/Data_Structures/OCTREE.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_COLOR.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_OBJECT.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_SELECTION.h>

namespace PhysBAM{
template<class T>
class OPENGL_OCTREE: public OPENGL_OBJECT
{
  public:
    typedef VECTOR<T,3> TV;typedef VECTOR<int,3> T_INDEX;

    OCTREE<T> tree;
    bool draw;
    int current_frame;
    RANGE<VECTOR<float,3> > domain;
  private:
    OPENGL_SELECTION *current_selection;

  public:
    OPENGL_OCTREE()
        :draw(true),current_frame(-1),current_selection(0)
    {}

//#####################################################################
    virtual RANGE<VECTOR<float,3> > Bounding_Box() const PHYSBAM_OVERRIDE;
    void Display(const int in_color=1) const PHYSBAM_OVERRIDE;
    void Draw_Cube(const ARRAY<TV>& nodes) const;
    void Update_Domain();
    virtual OPENGL_SELECTION *Get_Selection(GLuint *buffer, int buffer_size);
    void Highlight_Selection(OPENGL_SELECTION *selection) PHYSBAM_OVERRIDE;
    void Clear_Highlight() PHYSBAM_OVERRIDE;
    void Print_Selection_Info(std::ostream& stream,OPENGL_SELECTION* selection) const PHYSBAM_OVERRIDE;
//#####################################################################
};

template<class T>
class OPENGL_SELECTION_OCTREE_CELL: public OPENGL_SELECTION
{
    typedef VECTOR<T,3> TV;
  public:
    int index;

    OPENGL_SELECTION_OCTREE_CELL(OPENGL_OBJECT *object,const int& index_input=0) 
        :OPENGL_SELECTION(OPENGL_SELECTION::GRID_CELL_3D,object),index(index_input)
    {}

    RANGE<VECTOR<float,3> > Bounding_Box() const PHYSBAM_OVERRIDE;
};
}
#endif
