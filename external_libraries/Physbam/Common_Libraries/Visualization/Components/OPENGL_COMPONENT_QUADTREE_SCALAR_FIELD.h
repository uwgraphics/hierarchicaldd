//#####################################################################
// Copyright 2015, Mridul Aanjaneya.
// This file is part of PhysBAM whose distribution is governed by the license 
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_COMPONENT_QUADTREE_SCALAR_FIELD
//##################################################################### 
#ifndef __OPENGL_COMPONENT_QUADTREE_SCALAR_FIELD__
#define __OPENGL_COMPONENT_QUADTREE_SCALAR_FIELD__

#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT.h>
#include <Visualization/Objects/OPENGL_QUADTREE_SCALAR_FIELD.h>
#include <string>

namespace PhysBAM{
template<class T,class T2=T,class RW=T>
class OPENGL_COMPONENT_QUADTREE_SCALAR_FIELD: public OPENGL_COMPONENT
{
    typedef VECTOR<T,2> TV;
public:
    OPENGL_COMPONENT_QUADTREE_SCALAR_FIELD(OPENGL_QUADTREE<T>& opengl_quadtree,const std::string &filename_input,OPENGL_COLOR_MAP<T2>* color_map_input);
    ~OPENGL_COMPONENT_QUADTREE_SCALAR_FIELD() {}

    bool Valid_Frame(int frame_input) const PHYSBAM_OVERRIDE;
    bool Is_Up_To_Date(int frame) const PHYSBAM_OVERRIDE {return valid && frame_loaded==frame;}

    void Set_Frame(int frame_input) PHYSBAM_OVERRIDE;
    void Set_Draw(bool draw_input=true) PHYSBAM_OVERRIDE;

    void Display(const int in_color=1) const PHYSBAM_OVERRIDE;
    bool Use_Bounding_Box() const PHYSBAM_OVERRIDE {return draw&&valid;}
    virtual RANGE<VECTOR<float,3> > Bounding_Box() const PHYSBAM_OVERRIDE;

    void Print_Selection_Info(std::ostream& output_stream,OPENGL_SELECTION* current_selection) const PHYSBAM_OVERRIDE;

    void Toggle_Color_Map();
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_QUADTREE_SCALAR_FIELD, Toggle_Color_Map, "Toggle color map");

private:
    void Reinitialize();

public:
    OPENGL_QUADTREE_SCALAR_FIELD<T,T2> opengl_scalar_field;

private:
    std::string filename;
    int frame_loaded;
    bool valid;
};
}
#endif
