//#####################################################################
// Copyright (c) 2014, Mridul Aanjaneya.
// Distributed under the FreeBSD license (see license.txt)
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Read_Write/Grids_Uniform_Arrays/READ_WRITE_ARRAYS.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <SPGrid_Fluids/Visualization/OPENGL_COMPONENT_SPGRID_SCALAR_FIELD_2D.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T,class T2,class RW> OPENGL_COMPONENT_SPGRID_SCALAR_FIELD_2D<T,T2,RW>::
OPENGL_COMPONENT_SPGRID_SCALAR_FIELD_2D(GRID<TV> &fine_mac_grid,OPENGL_SPGRID_2D<T>& opengl_spgrid,const std::string &channel_filename_input,T T_STRUCT::* field_input,OPENGL_COLOR_MAP<T2>* color_map_input,bool& draw_separate_levels,int& level_counter)
    : OPENGL_COMPONENT("SPGrid Scalar Field 2D"),opengl_scalar_field(fine_mac_grid,opengl_spgrid,color_map_input,draw_separate_levels,level_counter,field_input),
      channel_filename(channel_filename_input),field(field_input),frame_loaded(-1),valid(false)
{
    is_animation = FILE_UTILITIES::Is_Animated(channel_filename);
}
//#####################################################################
// Function Valid_Frame
//#####################################################################
template<class T,class T2,class RW> bool OPENGL_COMPONENT_SPGRID_SCALAR_FIELD_2D<T,T2,RW>::
Valid_Frame(int frame_input) const
{
    return FILE_UTILITIES::Frame_File_Exists(channel_filename,frame_input);
}
//#####################################################################
// Function Set_Frame
//#####################################################################
template<class T,class T2,class RW> void OPENGL_COMPONENT_SPGRID_SCALAR_FIELD_2D<T,T2,RW>::
Set_Frame(int frame_input)
{
    OPENGL_COMPONENT::Set_Frame(frame_input);
    Reinitialize();
}
//#####################################################################
// Function Set_Draw
//#####################################################################
template<class T,class T2,class RW> void OPENGL_COMPONENT_SPGRID_SCALAR_FIELD_2D<T,T2,RW>::
Set_Draw(bool draw_input)
{
    OPENGL_COMPONENT::Set_Draw(draw_input);
    Reinitialize();
}
//#####################################################################
// Function Display
//#####################################################################
template<class T,class T2,class RW> void OPENGL_COMPONENT_SPGRID_SCALAR_FIELD_2D<T,T2,RW>::
Display(const int in_color) const
{
    if (valid && draw) opengl_scalar_field.Display(in_color);
}
//#####################################################################
// Function Bounding_Box
//#####################################################################
template<class T,class T2,class RW> RANGE<VECTOR<float,3> > OPENGL_COMPONENT_SPGRID_SCALAR_FIELD_2D<T,T2,RW>::
Bounding_Box() const
{
    if (valid && draw) return opengl_scalar_field.Bounding_Box();
    else return RANGE<VECTOR<float,3> >::Centered_Box();
}
//#####################################################################
// Function Print_Selection_Info
//#####################################################################
template<class T,class T2,class RW> void OPENGL_COMPONENT_SPGRID_SCALAR_FIELD_2D<T,T2,RW>::
Print_Selection_Info(std::ostream& output_stream,OPENGL_SELECTION* current_selection) const
{
    if(Is_Up_To_Date(frame)){
        output_stream<<component_name<<": ";
        opengl_scalar_field.Print_Selection_Info(output_stream,current_selection);}
}
//#####################################################################
// Function Reinitialize
//#####################################################################
template<class T,class T2,class RW> void OPENGL_COMPONENT_SPGRID_SCALAR_FIELD_2D<T,T2,RW>::
Reinitialize()
{
    if(draw){
        if((is_animation && frame_loaded != frame) || (!is_animation && frame_loaded<0)){
            valid=false;
            std::string filename=FILE_UTILITIES::Get_Frame_Filename(channel_filename,frame);
            if(FILE_UTILITIES::File_Exists(filename)){
                if(opengl_scalar_field.opengl_spgrid.current_frame!=frame){
                    opengl_scalar_field.opengl_spgrid.current_frame=frame;
                    opengl_scalar_field.opengl_spgrid.Initialize_Hierarchy(frame);}
                opengl_scalar_field.opengl_spgrid.hierarchy->Read_Data_Channel(filename,field);}
            else return;

            frame_loaded=frame;
            valid=true;}}
}
//#####################################################################
// Function Toggle_Color_Map
//#####################################################################
template<class T,class T2,class RW> void OPENGL_COMPONENT_SPGRID_SCALAR_FIELD_2D<T,T2,RW>::
Toggle_Color_Map()
{
    opengl_scalar_field.Toggle_Color_Map();
}
//#####################################################################
// Function Toggle_Draw_Mode
//#####################################################################
template<class T,class T2,class RW> void OPENGL_COMPONENT_SPGRID_SCALAR_FIELD_2D<T,T2,RW>::
Toggle_Draw_Mode()
{
    opengl_scalar_field.Toggle_Draw_Mode();
}
//#####################################################################
template class OPENGL_COMPONENT_SPGRID_SCALAR_FIELD_2D<float,int,float>;
template class OPENGL_COMPONENT_SPGRID_SCALAR_FIELD_2D<float,bool,float>;
template class OPENGL_COMPONENT_SPGRID_SCALAR_FIELD_2D<float,float,float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class OPENGL_COMPONENT_SPGRID_SCALAR_FIELD_2D<double,int,double>;
template class OPENGL_COMPONENT_SPGRID_SCALAR_FIELD_2D<double,bool,double>;
template class OPENGL_COMPONENT_SPGRID_SCALAR_FIELD_2D<double,double,double>;
#endif
