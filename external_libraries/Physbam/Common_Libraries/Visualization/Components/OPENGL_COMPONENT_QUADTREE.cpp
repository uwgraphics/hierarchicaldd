//#####################################################################
// Copyright 2015, Mridul Aanjaneya, Ming Gao.
// This file is part of PhysBAM whose distribution is governed by the license 
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Common_Tools/Read_Write/Data_Structures/READ_WRITE_QUADTREE.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <Visualization/Components/OPENGL_COMPONENT_QUADTREE.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T,class RW> OPENGL_COMPONENT_QUADTREE<T,RW>::
OPENGL_COMPONENT_QUADTREE(const std::string& filename_input,const int frame_input)
    :OPENGL_COMPONENT("Quadtree"),filename(filename_input),frame_loaded(-1),valid(false)
{
    is_animation = FILE_UTILITIES::Is_Animated(filename);
    Set_Frame(frame_input);
}
//#####################################################################
// Valid_Frame
//#####################################################################
template<class T,class RW> bool OPENGL_COMPONENT_QUADTREE<T,RW>::
Valid_Frame(int frame_input) const
{
    return FILE_UTILITIES::Frame_File_Exists(filename,frame_input);
}
//#####################################################################
// Set_Frame
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_QUADTREE<T,RW>::
Set_Frame(int frame_input)
{
    OPENGL_COMPONENT::Set_Frame(frame_input);
    Reinitialize();
}
//#####################################################################
// Set_Draw
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_QUADTREE<T,RW>::
Set_Draw(bool draw_input)
{
    OPENGL_COMPONENT::Set_Draw(draw_input);
    Reinitialize();
}
//#####################################################################
// Toggle_Draw_Centers
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_QUADTREE<T,RW>::
Toggle_Draw_Centers()
{
    opengl_quadtree.draw_centers=!opengl_quadtree.draw_centers;
}
//#####################################################################
// Display
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_QUADTREE<T,RW>::
Display(const int in_color) const
{
    if(valid && draw) opengl_quadtree.Display(in_color);
}
//#####################################################################
// Bounding_Box
//#####################################################################
template<class T,class RW> RANGE<VECTOR<float,3> > OPENGL_COMPONENT_QUADTREE<T,RW>::
Bounding_Box() const
{
    if(valid && draw) return opengl_quadtree.Bounding_Box();
    else return RANGE<VECTOR<float,3> >::Centered_Box();
}
//#####################################################################
// Function Print_Selection_Info
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_QUADTREE<T,RW>::
Print_Selection_Info(std::ostream& output_stream,OPENGL_SELECTION* current_selection) const
{
    if(Is_Up_To_Date(frame)){
        opengl_quadtree.Print_Selection_Info(output_stream,current_selection);}
}
//#####################################################################
// Reinitialize
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_QUADTREE<T,RW>::
Reinitialize()
{
    if(draw){
        if((is_animation && frame_loaded != frame) || (!is_animation && frame_loaded<0)){
            valid=false;
            std::string current_filename=FILE_UTILITIES::Get_Frame_Filename(filename,frame);
            if(FILE_UTILITIES::File_Exists(current_filename)){
                if(opengl_quadtree.current_frame!=frame){
                    opengl_quadtree.current_frame=frame;
                    FILE_UTILITIES::Read_From_File<RW>(current_filename,opengl_quadtree.tree);
                    opengl_quadtree.Update_Domain();}}
            else return;

            frame_loaded=frame;
            valid=true;}}
}
//#####################################################################
template class OPENGL_COMPONENT_QUADTREE<float,float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class OPENGL_COMPONENT_QUADTREE<double,double>;
#endif
