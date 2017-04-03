//#####################################################################
// Copyright 2015, Mridul Aanjaneya, Ming Gao.
// This file is part of PhysBAM whose distribution is governed by the license 
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Common_Tools/Read_Write/Data_Structures/READ_WRITE_OCTREE.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <Visualization/Components/OPENGL_COMPONENT_OCTREE.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T,class RW> OPENGL_COMPONENT_OCTREE<T,RW>::
OPENGL_COMPONENT_OCTREE(const std::string& filename_input,const int frame_input)
    :OPENGL_COMPONENT("Octree"),filename(filename_input),frame_loaded(-1),valid(false)
{
    is_animation = FILE_UTILITIES::Is_Animated(filename);
    Set_Frame(frame_input);
}
//#####################################################################
// Valid_Frame
//#####################################################################
template<class T,class RW> bool OPENGL_COMPONENT_OCTREE<T,RW>::
Valid_Frame(int frame_input) const
{
    return FILE_UTILITIES::Frame_File_Exists(filename,frame_input);
}
//#####################################################################
// Set_Frame
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_OCTREE<T,RW>::
Set_Frame(int frame_input)
{
    OPENGL_COMPONENT::Set_Frame(frame_input);
    Reinitialize();
}
//#####################################################################
// Set_Draw
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_OCTREE<T,RW>::
Set_Draw(bool draw_input)
{
    OPENGL_COMPONENT::Set_Draw(draw_input);
    Reinitialize();
}
//#####################################################################
// Display
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_OCTREE<T,RW>::
Display(const int in_color) const
{
    if(valid && draw) opengl_octree.Display(in_color);
}
//#####################################################################
// Bounding_Box
//#####################################################################
template<class T,class RW> RANGE<VECTOR<float,3> > OPENGL_COMPONENT_OCTREE<T,RW>::
Bounding_Box() const
{
    if(valid && draw) return opengl_octree.Bounding_Box();
    else return RANGE<VECTOR<float,3> >::Centered_Box();
}
//#####################################################################
// Print_Selection_Info
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_OCTREE<T,RW>::
Print_Selection_Info(std::ostream& output_stream,OPENGL_SELECTION* selection) const
{   
    if(Is_Up_To_Date(frame)){
        opengl_octree.Print_Selection_Info(output_stream,selection);}
}
//#####################################################################
// Reinitialize
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_OCTREE<T,RW>::
Reinitialize()
{
    if(draw){
        if((is_animation && frame_loaded != frame) || (!is_animation && frame_loaded<0)){
            valid=false;
            std::string current_filename=FILE_UTILITIES::Get_Frame_Filename(filename,frame);
            if(FILE_UTILITIES::File_Exists(current_filename)){
                if(opengl_octree.current_frame!=frame){
                    opengl_octree.current_frame=frame;
                    FILE_UTILITIES::Read_From_File<RW>(current_filename,opengl_octree.tree);
                    opengl_octree.Update_Domain();}}
            else return;

            frame_loaded=frame;
            valid=true;}}
}
//#####################################################################
template class OPENGL_COMPONENT_OCTREE<float,float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class OPENGL_COMPONENT_OCTREE<double,double>;
#endif
