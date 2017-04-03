//#####################################################################
// Copyright 2015, Mridul Aanjaneya.
// This file is part of PhysBAM whose distribution is governed by the license 
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Common_Geometry/Read_Write/Topology_Based_Geometry/READ_WRITE_VORONOI_DIAGRAM.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <Visualization/Components/OPENGL_COMPONENT_VORONOI_DIAGRAM_3D.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T,class RW> OPENGL_COMPONENT_VORONOI_DIAGRAM_3D<T,RW>::
OPENGL_COMPONENT_VORONOI_DIAGRAM_3D(const std::string& filename_input,const int frame_input)
    :OPENGL_COMPONENT("Voronoi 3D"),filename(filename_input),frame_loaded(-1),valid(false)
{
    is_animation = FILE_UTILITIES::Is_Animated(filename);
    Set_Frame(frame_input);
}
//#####################################################################
// Valid_Frame
//#####################################################################
template<class T,class RW> bool OPENGL_COMPONENT_VORONOI_DIAGRAM_3D<T,RW>::
Valid_Frame(int frame_input) const
{
    return FILE_UTILITIES::Frame_File_Exists(filename,frame_input);
}
//#####################################################################
// Set_Frame
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_VORONOI_DIAGRAM_3D<T,RW>::
Set_Frame(int frame_input)
{
    OPENGL_COMPONENT::Set_Frame(frame_input);
    Reinitialize();
}
//#####################################################################
// Set_Draw
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_VORONOI_DIAGRAM_3D<T,RW>::
Set_Draw(bool draw_input)
{
    OPENGL_COMPONENT::Set_Draw(draw_input);
    Reinitialize();
}
//#####################################################################
// Display
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_VORONOI_DIAGRAM_3D<T,RW>::
Display(const int in_color) const
{
    if(valid && draw) opengl_voronoi_diagram.Display(in_color);
}
//#####################################################################
// Bounding_Box
//#####################################################################
template<class T,class RW> RANGE<VECTOR<float,3> > OPENGL_COMPONENT_VORONOI_DIAGRAM_3D<T,RW>::
Bounding_Box() const
{
    if(valid && draw) return opengl_voronoi_diagram.Bounding_Box();
    else return RANGE<VECTOR<float,3> >::Centered_Box();
}
//#####################################################################
// Reinitialize
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_VORONOI_DIAGRAM_3D<T,RW>::
Reinitialize()
{
    if(draw){
        if((is_animation && frame_loaded != frame) || (!is_animation && frame_loaded<0)){
            valid=false;
            std::string current_filename=FILE_UTILITIES::Get_Frame_Filename(filename,frame);
            if(FILE_UTILITIES::File_Exists(current_filename)){
                if(opengl_voronoi_diagram.current_frame!=frame){
                    opengl_voronoi_diagram.current_frame=frame;
                    FILE_UTILITIES::Read_From_File<RW>(current_filename,opengl_voronoi_diagram.voronoi_diagram);
                    opengl_voronoi_diagram.Update_Domain();}}
            else return;

            frame_loaded=frame;
            valid=true;}}
}
//#####################################################################
template class OPENGL_COMPONENT_VORONOI_DIAGRAM_3D<float,float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class OPENGL_COMPONENT_VORONOI_DIAGRAM_3D<double,double>;
#endif
