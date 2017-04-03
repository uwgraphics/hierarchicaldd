//#####################################################################
// Copyright (c) 2014, Mridul Aanjaneya.
// Distributed under the FreeBSD license (see license.txt)
//#####################################################################
// Class OPENGL_COMPONENT_SPGRID_2D
//##################################################################### 
#ifndef __OPENGL_COMPONENT_SPGRID_2D__
#define __OPENGL_COMPONENT_SPGRID_2D__

#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT.h>
#include <SPGrid_Fluids/Visualization/OPENGL_SPGRID_2D.h>

namespace PhysBAM{
template<class T,class RW=T>
class OPENGL_COMPONENT_SPGRID_2D : public OPENGL_COMPONENT
{
    typedef FLUIDS_SIMULATION_DATA<T> T_STRUCT;typedef VECTOR<T,2> TV;
public:
    GRID<TV>& fine_grid;
    OPENGL_SPGRID_2D<T> opengl_spgrid;
private:
    std::string filename;
    int frame_loaded;
    bool valid;
public:

    virtual OPENGL_SELECTION *Get_Selection(GLuint *buffer, int buffer_size)
    {return opengl_spgrid.Get_Selection(buffer,buffer_size);}

    virtual void Highlight_Selection(OPENGL_SELECTION *selection) PHYSBAM_OVERRIDE
    {opengl_spgrid.Highlight_Selection(selection);}

    virtual void Clear_Highlight() PHYSBAM_OVERRIDE
    {opengl_spgrid.Clear_Highlight();}

//#####################################################################
    OPENGL_COMPONENT_SPGRID_2D(GRID<TV> &fine_grid_input,const std::string& filename_input,bool& draw_separate_levels,int& level_counter,bool& draw_nodal_flags,int& nodal_flag_counter,const int frame_input);
    bool Valid_Frame(int frame_input) const PHYSBAM_OVERRIDE;
    bool Is_Up_To_Date(int frame) const PHYSBAM_OVERRIDE {return valid && frame_loaded==frame;}
    void Set_Frame(int frame_input) PHYSBAM_OVERRIDE;
    void Set_Draw(bool draw_input=true) PHYSBAM_OVERRIDE;
    void Display(const int in_color=1) const PHYSBAM_OVERRIDE;
    bool Use_Bounding_Box() const PHYSBAM_OVERRIDE {return draw && valid;}
    virtual RANGE<VECTOR<float,3> > Bounding_Box() const PHYSBAM_OVERRIDE;
    void Print_Selection_Info(std::ostream& output_stream,OPENGL_SELECTION* selection) const PHYSBAM_OVERRIDE;
private:
    void Reinitialize();
//#####################################################################
};
}
#endif
