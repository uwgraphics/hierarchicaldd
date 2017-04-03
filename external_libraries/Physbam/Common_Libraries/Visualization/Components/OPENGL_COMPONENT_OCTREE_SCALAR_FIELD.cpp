//#####################################################################
// Copyright 2015, Mridul Aanjaneya.
// This file is part of PhysBAM whose distribution is governed by the license 
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Read_Write/Data_Structures/READ_WRITE_HASHTABLE.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <Visualization/Components/OPENGL_COMPONENT_OCTREE_SCALAR_FIELD.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T,class T2,class RW> OPENGL_COMPONENT_OCTREE_SCALAR_FIELD<T,T2,RW>::
OPENGL_COMPONENT_OCTREE_SCALAR_FIELD(OPENGL_OCTREE<T>& opengl_octree,const std::string &filename_input,OPENGL_COLOR_MAP<T2>* color_map_input)
    :OPENGL_COMPONENT("Octree Scalar Field"),opengl_scalar_field(opengl_octree,color_map_input),filename(filename_input),frame_loaded(-1),valid(false)
{
    is_animation = FILE_UTILITIES::Is_Animated(filename);
}
//#####################################################################
// Function Valid_Frame
//#####################################################################
template<class T,class T2,class RW> bool OPENGL_COMPONENT_OCTREE_SCALAR_FIELD<T,T2,RW>::
Valid_Frame(int frame_input) const
{
    return FILE_UTILITIES::Frame_File_Exists(filename,frame_input);
}
//#####################################################################
// Function Set_Frame
//#####################################################################
template<class T,class T2,class RW> void OPENGL_COMPONENT_OCTREE_SCALAR_FIELD<T,T2,RW>::
Set_Frame(int frame_input)
{
    OPENGL_COMPONENT::Set_Frame(frame_input);
    Reinitialize();
}
//#####################################################################
// Function Set_Draw
//#####################################################################
template<class T,class T2,class RW> void OPENGL_COMPONENT_OCTREE_SCALAR_FIELD<T,T2,RW>::
Set_Draw(bool draw_input)
{
    OPENGL_COMPONENT::Set_Draw(draw_input);
    Reinitialize();
}
//#####################################################################
// Function Display
//#####################################################################
template<class T,class T2,class RW> void OPENGL_COMPONENT_OCTREE_SCALAR_FIELD<T,T2,RW>::
Display(const int in_color) const
{
    if (valid && draw) opengl_scalar_field.Display(in_color);
}
//#####################################################################
// Function Bounding_Box
//#####################################################################
template<class T,class T2,class RW> RANGE<VECTOR<float,3> > OPENGL_COMPONENT_OCTREE_SCALAR_FIELD<T,T2,RW>::
Bounding_Box() const
{
    if (valid && draw) return opengl_scalar_field.Bounding_Box();
    else return RANGE<VECTOR<float,3> >::Centered_Box();
}
//#####################################################################
// Function Print_Selection_Info
//#####################################################################
template<class T,class T2,class RW> void OPENGL_COMPONENT_OCTREE_SCALAR_FIELD<T,T2,RW>::
Print_Selection_Info(std::ostream& output_stream,OPENGL_SELECTION* current_selection) const
{
    if(Is_Up_To_Date(frame)){
        output_stream<<component_name<<": ";
        opengl_scalar_field.Print_Selection_Info(output_stream,current_selection);}
}
//#####################################################################
// Function Reinitialize
//#####################################################################
template<class T,class T2,class RW> void OPENGL_COMPONENT_OCTREE_SCALAR_FIELD<T,T2,RW>::
Reinitialize()
{
    if(draw){
        if((is_animation && frame_loaded != frame) || (!is_animation && frame_loaded<0)){
            valid=false;
            std::string data_filename=FILE_UTILITIES::Get_Frame_Filename(filename,frame);
            if(FILE_UTILITIES::File_Exists(data_filename)){
                if(opengl_scalar_field.opengl_octree.current_frame!=frame){
                    opengl_scalar_field.opengl_octree.current_frame=frame;}
                FILE_UTILITIES::Read_From_File<RW>(data_filename,opengl_scalar_field.data);}
            else return;

            frame_loaded=frame;
            valid=true;}}
}
//#####################################################################
// Function Toggle_Color_Map
//#####################################################################
template<class T,class T2,class RW> void OPENGL_COMPONENT_OCTREE_SCALAR_FIELD<T,T2,RW>::
Toggle_Color_Map()
{
    opengl_scalar_field.Toggle_Color_Map();
}
//#####################################################################
template class OPENGL_COMPONENT_OCTREE_SCALAR_FIELD<float,int,float>;
template class OPENGL_COMPONENT_OCTREE_SCALAR_FIELD<float,bool,float>;
template class OPENGL_COMPONENT_OCTREE_SCALAR_FIELD<float,float,float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class OPENGL_COMPONENT_OCTREE_SCALAR_FIELD<double,int,double>;
template class OPENGL_COMPONENT_OCTREE_SCALAR_FIELD<double,bool,double>;
template class OPENGL_COMPONENT_OCTREE_SCALAR_FIELD<double,double,double>;
#endif
