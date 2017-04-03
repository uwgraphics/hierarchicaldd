//#####################################################################
// Copyright (c) 2014, Mridul Aanjaneya.
// Distributed under the FreeBSD license (see license.txt)
//#####################################################################
#include <PhysBAM_Tools/Read_Write/Grids_Uniform/READ_WRITE_GRID.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <SPGrid_Fluids/Rendering/OpenGL/OPENGL_COMPONENT_SPGRID_3D.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T,class RW> OPENGL_COMPONENT_SPGRID_3D<T,RW>::
OPENGL_COMPONENT_SPGRID_3D(GRID<TV> &fine_grid_input,const std::string& filename_input,bool& draw_separate_levels,int& level_counter,bool& draw_nodal_flags,int& nodal_flag_counter,const int frame_input)
    :OPENGL_COMPONENT("SPGrid 3D"),fine_grid(fine_grid_input),opengl_spgrid(fine_grid,draw_separate_levels,level_counter,draw_nodal_flags,nodal_flag_counter,filename_input,OPENGL_COLOR::Gray(.5),frame_input),filename(filename_input),frame_loaded(-1),valid(false)
{
    is_animation = FILE_UTILITIES::Is_Animated(filename);
    Set_Frame(frame_input);
}
//#####################################################################
// Function Valid_Frame
//#####################################################################
template<class T,class RW> bool OPENGL_COMPONENT_SPGRID_3D<T,RW>::
Valid_Frame(int frame_input) const
{
    return FILE_UTILITIES::Frame_File_Exists(filename,frame_input);
}
//#####################################################################
// Function Set_Frame
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_SPGRID_3D<T,RW>::
Set_Frame(int frame_input)
{
    OPENGL_COMPONENT::Set_Frame(frame_input);
    Reinitialize();
}
//#####################################################################
// Function Set_Draw
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_SPGRID_3D<T,RW>::
Set_Draw(bool draw_input)
{
    OPENGL_COMPONENT::Set_Draw(draw_input);
    Reinitialize();
}
//#####################################################################
// Function Display
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_SPGRID_3D<T,RW>::
Display(const int in_color) const
{
    if (valid && draw) opengl_spgrid.Display(in_color);
}
//#####################################################################
// Function Bounding_Box
//#####################################################################
template<class T,class RW> RANGE<VECTOR<float,3> > OPENGL_COMPONENT_SPGRID_3D<T,RW>::
Bounding_Box() const
{
    if (valid && draw) return opengl_spgrid.Bounding_Box();
    else return RANGE<VECTOR<float,3> >::Centered_Box();
}
//#####################################################################
// Function Print_Selection_Info
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_SPGRID_3D<T,RW>::
Print_Selection_Info(std::ostream& output_stream,OPENGL_SELECTION* current_selection) const
{
    if(Is_Up_To_Date(frame)){
        opengl_spgrid.Print_Selection_Info(output_stream,current_selection);}
}
//#####################################################################
// Function Reinitialize
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_SPGRID_3D<T,RW>::
Reinitialize()
{
    if(draw){
        if((is_animation && frame_loaded != frame) || (!is_animation && frame_loaded<0)){
            valid=false;
            std::string directory_filename=FILE_UTILITIES::Get_Frame_Filename(filename,frame);
            std::string flags_filename=directory_filename+"/flags";
            if(FILE_UTILITIES::File_Exists(flags_filename)){
                if(opengl_spgrid.current_frame!=frame){
                    opengl_spgrid.current_frame=frame;
                    opengl_spgrid.Initialize_Hierarchy(frame);}}
            else return;

            frame_loaded=frame;
            valid=true;}}
}
//#####################################################################
template class OPENGL_COMPONENT_SPGRID_3D<float,float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class OPENGL_COMPONENT_SPGRID_3D<double,double>;
#endif
