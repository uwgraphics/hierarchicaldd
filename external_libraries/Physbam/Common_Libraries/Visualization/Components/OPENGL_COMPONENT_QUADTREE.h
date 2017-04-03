//#####################################################################
// Copyright 2015, Mridul Aanjaneya, Ming Gao.
// This file is part of PhysBAM whose distribution is governed by the license 
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_COMPONENT_QUADTREE
//#####################################################################
#ifndef __OPENGL_COMPONENT_QUADTREE__
#define __OPENGL_COMPONENT_QUADTREE__

#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT.h>
#include <Visualization/Objects/OPENGL_QUADTREE.h>

namespace PhysBAM{
template<class T,class RW=T>
class OPENGL_COMPONENT_QUADTREE: public OPENGL_COMPONENT
{
    typedef VECTOR<T,2> TV;
  public:
    OPENGL_QUADTREE<T> opengl_quadtree;
  private:
    std::string filename;
    int frame_loaded;
    bool valid;
  public:

    OPENGL_COMPONENT_QUADTREE(const std::string& filename_input,const int frame_input);

    bool Is_Up_To_Date(int frame) const PHYSBAM_OVERRIDE {return valid && frame_loaded==frame;}

    bool Use_Bounding_Box() const PHYSBAM_OVERRIDE {return draw && valid;}

    virtual OPENGL_SELECTION *Get_Selection(GLuint *buffer, int buffer_size)
    {return opengl_quadtree.Get_Selection(buffer,buffer_size);}

    virtual void Highlight_Selection(OPENGL_SELECTION *selection) PHYSBAM_OVERRIDE
    {opengl_quadtree.Highlight_Selection(selection);}

    virtual void Clear_Highlight() PHYSBAM_OVERRIDE
    {opengl_quadtree.Clear_Highlight();}

//#####################################################################
    bool Valid_Frame(int frame_input) const PHYSBAM_OVERRIDE;
    void Set_Frame(int frame_input) PHYSBAM_OVERRIDE;
    void Set_Draw(bool draw_input=true) PHYSBAM_OVERRIDE;
    void Display(const int in_color=1) const PHYSBAM_OVERRIDE;
    virtual RANGE<VECTOR<float,3> > Bounding_Box() const PHYSBAM_OVERRIDE;
    void Toggle_Draw_Centers();
    void Print_Selection_Info(std::ostream& output_stream,OPENGL_SELECTION* selection) const PHYSBAM_OVERRIDE;

    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_QUADTREE, Toggle_Draw_Centers, "Toggle centers");

  private:
    void Reinitialize();
//#####################################################################
};
}
#endif
