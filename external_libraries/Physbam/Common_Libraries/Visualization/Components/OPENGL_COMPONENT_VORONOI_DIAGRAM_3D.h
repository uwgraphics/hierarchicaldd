//#####################################################################
// Copyright 2015, Mridul Aanjaneya.
// This file is part of PhysBAM whose distribution is governed by the license 
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_COMPONENT_VORONOI_DIAGRAM_3D
//#####################################################################
#ifndef __OPENGL_COMPONENT_VORONOI_DIAGRAM_3D__
#define __OPENGL_COMPONENT_VORONOI_DIAGRAM_3D__

#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT.h>
#include <Visualization/Objects/OPENGL_VORONOI_DIAGRAM_3D.h>

namespace PhysBAM{
template<class T,class RW=T>
class OPENGL_COMPONENT_VORONOI_DIAGRAM_3D: public OPENGL_COMPONENT
{
    typedef VECTOR<T,3> TV;
  public:
    OPENGL_VORONOI_DIAGRAM_3D<T> opengl_voronoi_diagram;
  private:
    std::string filename;
    int frame_loaded;
    bool valid;
  public:

    OPENGL_COMPONENT_VORONOI_DIAGRAM_3D(const std::string& filename_input,const int frame_input);

    bool Is_Up_To_Date(int frame) const PHYSBAM_OVERRIDE {return valid && frame_loaded==frame;}

    bool Use_Bounding_Box() const PHYSBAM_OVERRIDE {return draw && valid;}

//#####################################################################
    bool Valid_Frame(int frame_input) const PHYSBAM_OVERRIDE;
    void Set_Frame(int frame_input) PHYSBAM_OVERRIDE;
    void Set_Draw(bool draw_input=true) PHYSBAM_OVERRIDE;
    void Display(const int in_color=1) const PHYSBAM_OVERRIDE;
    virtual RANGE<VECTOR<float,3> > Bounding_Box() const PHYSBAM_OVERRIDE;
  private:
    void Reinitialize();
//#####################################################################
};
}
#endif
