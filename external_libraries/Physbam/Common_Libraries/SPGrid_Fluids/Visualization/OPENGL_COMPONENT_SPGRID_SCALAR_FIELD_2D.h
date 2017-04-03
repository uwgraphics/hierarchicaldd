//#####################################################################
// Copyright (c) 2014, Mridul Aanjaneya.
// Distributed under the FreeBSD license (see license.txt)
//#####################################################################
// Class OPENGL_COMPONENT_SPGRID_SCALAR_FIELD_2D
//##################################################################### 
#ifndef __OPENGL_COMPONENT_SPGRID_SCALAR_FIELD_2D__
#define __OPENGL_COMPONENT_SPGRID_SCALAR_FIELD_2D__

#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT.h>
#include <SPGrid_Fluids/Grids/GRID_HIERARCHY.h>
#include <SPGrid_Fluids/Simulation/FLUIDS_SIMULATION_DATA.h>
#include <SPGrid_Fluids/Visualization/OPENGL_SPGRID_2D.h>
#include <SPGrid_Fluids/Visualization/OPENGL_SPGRID_SCALAR_FIELD_2D.h>
#include <string>

namespace PhysBAM{
template<class T,class T2=T,class RW=T>
class OPENGL_COMPONENT_SPGRID_SCALAR_FIELD_2D : public OPENGL_COMPONENT
{
    typedef VECTOR<T,2> TV;typedef FLUIDS_SIMULATION_DATA<T> T_STRUCT;
public:
    OPENGL_COMPONENT_SPGRID_SCALAR_FIELD_2D(GRID<TV> &fine_mac_grid,OPENGL_SPGRID_2D<T>& opengl_spgrid,const std::string &channel_filename_input,T T_STRUCT::* field_input,OPENGL_COLOR_MAP<T2>* color_map_input,bool& draw_separate_levels,int& level_counter);
    ~OPENGL_COMPONENT_SPGRID_SCALAR_FIELD_2D() {}

    bool Valid_Frame(int frame_input) const PHYSBAM_OVERRIDE;
    bool Is_Up_To_Date(int frame) const PHYSBAM_OVERRIDE {return valid && frame_loaded==frame;}

    void Set_Frame(int frame_input) PHYSBAM_OVERRIDE;
    void Set_Draw(bool draw_input=true) PHYSBAM_OVERRIDE;

    void Display(const int in_color=1) const PHYSBAM_OVERRIDE;
    bool Use_Bounding_Box() const PHYSBAM_OVERRIDE {return draw&&valid;}
    virtual RANGE<VECTOR<float,3> > Bounding_Box() const PHYSBAM_OVERRIDE;

    void Print_Selection_Info(std::ostream& output_stream,OPENGL_SELECTION* current_selection) const PHYSBAM_OVERRIDE;
    void Toggle_Color_Map();
    void Toggle_Draw_Mode();
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_SPGRID_SCALAR_FIELD_2D, Toggle_Color_Map, "Toggle color map");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_SPGRID_SCALAR_FIELD_2D, Toggle_Draw_Mode, "Toggle draw mode");

private:
    void Reinitialize();

public:
    OPENGL_SPGRID_SCALAR_FIELD_2D<T,T2> opengl_scalar_field;

private:
    std::string channel_filename;
    T T_STRUCT::* field;
    int frame_loaded;
    bool valid;
};
}
#endif
