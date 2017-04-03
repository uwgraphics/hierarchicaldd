//#####################################################################
// Copyright (c) 2014, Mridul Aanjaneya.
// Distributed under the FreeBSD license (see license.txt)
//#####################################################################
// Class OPENGL_COMPONENT_SPGRID_MAC_VELOCITY_FIELD_3D
//##################################################################### 
#ifndef __OPENGL_COMPONENT_SPGRID_MAC_VELOCITY_FIELD_3D__
#define __OPENGL_COMPONENT_SPGRID_MAC_VELOCITY_FIELD_3D__

#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT.h>
#include <SPGrid_Fluids/Grids/GRID_HIERARCHY.h>
#include <SPGrid_Fluids/Simulation/FLUIDS_SIMULATION_DATA.h>
#include <SPGrid_Fluids/Rendering/OpenGL/OPENGL_SPGRID_3D.h>
#include <SPGrid_Fluids/Rendering/OpenGL/OPENGL_SPGRID_MAC_VELOCITY_FIELD_3D.h>
#include <string>
namespace PhysBAM{
template<class T,class RW=T>
class OPENGL_COMPONENT_SPGRID_MAC_VELOCITY_FIELD_3D : public OPENGL_COMPONENT
{
    typedef VECTOR<T,3> TV;typedef FLUIDS_SIMULATION_DATA<T> T_STRUCT;
public:
    OPENGL_SPGRID_MAC_VELOCITY_FIELD_3D<T> opengl_mac_velocity_field;
private:
    std::string channel_filename1,channel_filename2,channel_filename3;
    T T_STRUCT::* field1;
    T T_STRUCT::* field2;
    T T_STRUCT::* field3;
    int frame_loaded;
    bool valid;

public:
    OPENGL_COMPONENT_SPGRID_MAC_VELOCITY_FIELD_3D(GRID<TV> &fine_mac_grid,OPENGL_SPGRID_3D<T>& opengl_spgrid,std::string& filename1_input,std::string& filename2_input,std::string& filename3_input,T T_STRUCT::* field1_input,T T_STRUCT::* field2_input,T T_STRUCT::* field3_input,bool& draw_separate_levels,int& level_counter);
    ~OPENGL_COMPONENT_SPGRID_MAC_VELOCITY_FIELD_3D() {}

    bool Valid_Frame(int frame_input) const PHYSBAM_OVERRIDE;
    bool Is_Up_To_Date(int frame) const PHYSBAM_OVERRIDE {return valid && frame_loaded==frame;}
    void Set_Frame(int frame_input) PHYSBAM_OVERRIDE;
    void Set_Draw(bool draw_input=true) PHYSBAM_OVERRIDE;
    void Display(const int in_color=1) const PHYSBAM_OVERRIDE;
    void Print_Selection_Info(std::ostream& stream,OPENGL_SELECTION* selection) const PHYSBAM_OVERRIDE;
    virtual void Set_Slice(OPENGL_SLICE *slice_input) PHYSBAM_OVERRIDE {slice=slice_input;opengl_mac_velocity_field.Set_Slice(slice_input);}
    bool Use_Bounding_Box() const PHYSBAM_OVERRIDE {return draw && valid;}
    virtual RANGE<VECTOR<float,3> > Bounding_Box() const PHYSBAM_OVERRIDE;

    void Increase_Vector_Size();
    void Decrease_Vector_Size();
    void Toggle_Arrowhead();

    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_SPGRID_MAC_VELOCITY_FIELD_3D, Toggle_Draw, "Toggle draw");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_SPGRID_MAC_VELOCITY_FIELD_3D, Increase_Vector_Size, "Increase vector size");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_SPGRID_MAC_VELOCITY_FIELD_3D, Decrease_Vector_Size, "Decrease vector size");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_SPGRID_MAC_VELOCITY_FIELD_3D, Toggle_Arrowhead, "Toggle arrowhead");

private:
    void Reinitialize();
};
}
#endif
