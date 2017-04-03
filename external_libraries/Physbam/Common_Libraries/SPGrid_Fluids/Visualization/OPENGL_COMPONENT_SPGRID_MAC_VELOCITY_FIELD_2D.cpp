//#####################################################################
// Copyright (c) 2014, Mridul Aanjaneya.
// Distributed under the FreeBSD license (see license.txt)
//#####################################################################
#include <PhysBAM_Tools/Math_Tools/RANGE.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <SPGrid_Fluids/Projection/GRID_HIERARCHY_PROJECTION.h>
#include <SPGrid_Fluids/Visualization/OPENGL_COMPONENT_SPGRID_MAC_VELOCITY_FIELD_2D.h>
using namespace PhysBAM;
//#####################################################################
// Function Constructor
//#####################################################################
template<class T,class RW> OPENGL_COMPONENT_SPGRID_MAC_VELOCITY_FIELD_2D<T,RW>::
OPENGL_COMPONENT_SPGRID_MAC_VELOCITY_FIELD_2D(GRID<TV> &fine_mac_grid,OPENGL_SPGRID_2D<T>& opengl_spgrid,std::string& filename1_input,std::string& filename2_input,T T_STRUCT::* field1_input,T T_STRUCT::* field2_input,bool& draw_separate_levels,int& level_counter)
    : OPENGL_COMPONENT("SPGrid MAC Velocity Field 2D"),opengl_mac_velocity_field(fine_mac_grid,opengl_spgrid,field1_input,field2_input,draw_separate_levels,level_counter),channel_filename1(filename1_input),channel_filename2(filename2_input),field1(field1_input),field2(field2_input)
{
    is_animation=FILE_UTILITIES::Is_Animated(channel_filename1);
    frame_loaded=-1;
}
//#####################################################################
// Function Valid_Frame
//#####################################################################
template<class T,class RW> bool OPENGL_COMPONENT_SPGRID_MAC_VELOCITY_FIELD_2D<T,RW>::
Valid_Frame(int frame_input) const
{
    return FILE_UTILITIES::Frame_File_Exists(channel_filename1,frame_input);
}
//#####################################################################
// Function Print_Selection_Info
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_SPGRID_MAC_VELOCITY_FIELD_2D<T,RW>::
Print_Selection_Info(std::ostream& stream,OPENGL_SELECTION* selection) const
{
    if(Is_Up_To_Date(frame)){
        stream<<component_name<<": "<<std::endl;
        opengl_mac_velocity_field.Print_Selection_Info(stream,selection);}
}
//#####################################################################
// Function Set_Frame
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_SPGRID_MAC_VELOCITY_FIELD_2D<T,RW>::
Set_Frame(int frame_input)
{
    OPENGL_COMPONENT::Set_Frame(frame_input);
    Reinitialize();
}
//#####################################################################
// Function Set_Frame
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_SPGRID_MAC_VELOCITY_FIELD_2D<T,RW>::
Set_Draw(bool draw_input)
{
    OPENGL_COMPONENT::Set_Draw(draw_input);
    Reinitialize();
}
//#####################################################################
// Function Display
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_SPGRID_MAC_VELOCITY_FIELD_2D<T,RW>::
Display(const int in_color) const
{
    if(valid && draw) opengl_mac_velocity_field.Display(in_color);
}
//#####################################################################
// Function Set_Frame
//#####################################################################
template<class T,class RW> RANGE<VECTOR<float,3> > OPENGL_COMPONENT_SPGRID_MAC_VELOCITY_FIELD_2D<T,RW>::
Bounding_Box() const
{
    if (valid && draw) return opengl_mac_velocity_field.Bounding_Box();
    else return RANGE<VECTOR<float,3> >::Centered_Box();
}
//#####################################################################
// Function Increase_Vector_Size
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_SPGRID_MAC_VELOCITY_FIELD_2D<T,RW>::
Increase_Vector_Size()
{
    opengl_mac_velocity_field.Scale_Vector_Size(1.1);
}
//#####################################################################
// Function Decrease_Vector_Size
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_SPGRID_MAC_VELOCITY_FIELD_2D<T,RW>::
Decrease_Vector_Size()
{
    opengl_mac_velocity_field.Scale_Vector_Size(1/1.1);
}
//#####################################################################
// Function Toggle_Arrowhead
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_SPGRID_MAC_VELOCITY_FIELD_2D<T,RW>::
Toggle_Arrowhead()
{
    opengl_mac_velocity_field.draw_arrowhead = !opengl_mac_velocity_field.draw_arrowhead;
}
//#####################################################################
// Function Reinitialize
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_SPGRID_MAC_VELOCITY_FIELD_2D<T,RW>::
Reinitialize()
{
    if(draw){
        if((is_animation && frame_loaded!=frame) || (!is_animation && frame_loaded<0)){
            valid=false;
            ARRAY<std::string> channel_filenames;
            channel_filenames.Append(channel_filename1);
            channel_filenames.Append(channel_filename2);
            ARRAY<T T_STRUCT::*> fields;
            fields.Append(field1);
            fields.Append(field2);
            if(FILE_UTILITIES::File_Exists(FILE_UTILITIES::Get_Frame_Filename(channel_filenames(1),frame)) && opengl_mac_velocity_field.opengl_spgrid.current_frame!=frame){
                opengl_mac_velocity_field.opengl_spgrid.current_frame=frame;
                opengl_mac_velocity_field.opengl_spgrid.Initialize_Hierarchy(frame);}
            for(int i=1;i<=channel_filenames.m;++i){
                std::string filename=FILE_UTILITIES::Get_Frame_Filename(channel_filenames(i),frame);
                if(FILE_UTILITIES::File_Exists(filename))
                    opengl_mac_velocity_field.opengl_spgrid.hierarchy->Read_Data_Channel(filename,fields(i));
                else return;}
            opengl_mac_velocity_field.divergence_channel=&T_STRUCT::ch7;
            GRID_HIERARCHY_PROJECTION<T_STRUCT,T,2>::Compute_Divergence(*(opengl_mac_velocity_field.opengl_spgrid.hierarchy),&T_STRUCT::flags,fields(1),fields(2),opengl_mac_velocity_field.divergence_channel);

            frame_loaded=frame;valid=true;}}
}
//#####################################################################
template class OPENGL_COMPONENT_SPGRID_MAC_VELOCITY_FIELD_2D<float,float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class OPENGL_COMPONENT_SPGRID_MAC_VELOCITY_FIELD_2D<double,double>;
#endif
