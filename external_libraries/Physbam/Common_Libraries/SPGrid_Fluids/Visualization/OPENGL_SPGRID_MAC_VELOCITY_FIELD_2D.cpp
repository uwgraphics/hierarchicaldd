//#####################################################################
// Copyright (c) 2014, Mridul Aanjaneya.
// Distributed under the FreeBSD license (see license.txt)
//#####################################################################
#include <PhysBAM_Tools/Math_Tools/RANGE.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_SHAPES.h>
#include <SPGrid_Fluids/Grids/GRID_TOPOLOGY_HELPER.h>
#include <SPGrid_Fluids/Simulation/FLUIDS_SIMULATION_FLAGS.h>
#include <SPGrid_Fluids/Visualization/OPENGL_SPGRID_2D.h>
#include <SPGrid_Fluids/Visualization/OPENGL_SPGRID_MAC_VELOCITY_FIELD_2D.h>
using namespace PhysBAM;
//#####################################################################
// Function Constructor
//#####################################################################
template<class T> OPENGL_SPGRID_MAC_VELOCITY_FIELD_2D<T>::
OPENGL_SPGRID_MAC_VELOCITY_FIELD_2D(GRID<TV> &fine_mac_grid_input,OPENGL_SPGRID_2D<T>& opengl_spgrid_input,T T_STRUCT::* field1_input,T T_STRUCT::* field2_input,bool& draw_separate_levels_input,int& level_counter_input,const OPENGL_COLOR &color,double size,bool draw_arrowhead)
    : fine_mac_grid(fine_mac_grid_input),opengl_spgrid(opengl_spgrid_input),field1(field1_input),field2(field2_input),draw_separate_levels(draw_separate_levels_input),level_counter(level_counter_input),vector_color(color),size(size),draw_arrowhead(draw_arrowhead),draw(true)
{}
//#####################################################################
// Function Bounding_Box
//#####################################################################
template<class T> RANGE<VECTOR<float,3> > OPENGL_SPGRID_MAC_VELOCITY_FIELD_2D<T>::
Bounding_Box() const
{
    return RANGE<VECTOR<float,3> >(fine_mac_grid.domain.min_corner.x,fine_mac_grid.domain.max_corner.x,fine_mac_grid.domain.min_corner.y,fine_mac_grid.domain.max_corner.y,0,0);
}
//#####################################################################
// Function Scale_Vector_Size
//#####################################################################
template<class T> void OPENGL_SPGRID_MAC_VELOCITY_FIELD_2D<T>::
Scale_Vector_Size(const T scale)
{
    size*=scale;
}
//#####################################################################
// Function Draw_Vector
//#####################################################################
template<class T> void OPENGL_SPGRID_MAC_VELOCITY_FIELD_2D<T>::
Draw_Vector(const TV& vector_field,const TV& vector_location,const double& current_size) const
{
    if(draw_arrowhead) OPENGL_SHAPES::Draw_Arrow(vector_location,vector_location+(T)current_size*vector_field);
    else OpenGL_Line(vector_location,vector_location+(T)current_size*vector_field);
}
//#####################################################################
// Function Display
//#####################################################################
template<class T> void OPENGL_SPGRID_MAC_VELOCITY_FIELD_2D<T>::
Display(const int in_color) const
{
    if(!draw)return;
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    Send_Transform_To_GL_Pipeline();
    
    glPushAttrib(GL_ENABLE_BIT | GL_LIGHTING_BIT | GL_TEXTURE_BIT);
    glDisable(GL_LIGHTING);glDisable(GL_TEXTURE_2D);glDisable(GL_DEPTH_TEST);

    vector_color.Send_To_GL_Pipeline();

    VECTOR<T T_STRUCT::*,2> face_velocity_channels;
    face_velocity_channels(1)=field1;face_velocity_channels(2)=field2;

    OpenGL_Begin(GL_LINES);
    double current_size=size;
    for(int level=opengl_spgrid.hierarchy->Levels();level>=1;--level){if(!draw_separate_levels || (draw_separate_levels && level==level_counter)){
        const GRID<TV> current_grid=opengl_spgrid.hierarchy->Grid(level).Get_MAC_Grid();
        VECTOR<void*,2> face_velocity_data_ptrs;
        for(int v=1;v<=TV::dimension;v++) face_velocity_data_ptrs(v)=opengl_spgrid.hierarchy->Array(level,face_velocity_channels(v)).Get_Data_Ptr();
        Const_flag_array_type flags=opengl_spgrid.hierarchy->Allocator(level).Get_Const_Array(&T_STRUCT::flags);    
        for(SPGrid_Block_Iterator<typename Data_array_type::MASK> iterator(opengl_spgrid.hierarchy->Blocks(level));iterator.Valid();iterator.Next()){
            TV_INT index=iterator.Index().template Cast<TV_INT>();
            for(int axis=1;axis<=TV::dimension;++axis){
                unsigned face_valid_mask=GRID_TOPOLOGY_HELPER<typename Flag_array_type::MASK>::Face_Valid_Mask(axis);
                if(iterator.Data(flags) & face_valid_mask){
                    TV vector_field=iterator.template Data<Data_array_type>(face_velocity_data_ptrs(axis))*TV::Axis_Vector(axis);
                    TV vector_location=current_grid.X(index)-(T).5*current_grid.dX(axis)*TV::Axis_Vector(axis);
                    Draw_Vector(vector_field,vector_location,current_size);}}}}
        current_size*=(T).5;}
    OpenGL_End();

    glPopAttrib();
    glPopMatrix();
}
//#####################################################################
// Print_Selection_Info
//#####################################################################
template<class T> void OPENGL_SPGRID_MAC_VELOCITY_FIELD_2D<T>::
Print_Selection_Info(std::ostream& output_stream,OPENGL_SELECTION* selection) const
{
    if(selection && selection->type==OPENGL_SELECTION::GRID_CELL_2D){
        VECTOR<T T_STRUCT::*,2> face_velocity_channels;
        VECTOR<unsigned long, GRID<TV>::dimension> neighbor_offsets;
        for(int v=1;v<=TV::dimension;v++) neighbor_offsets(v)=GRID_TOPOLOGY_HELPER<typename Flag_array_type::MASK>::Axis_Vector_Offset(v);
        face_velocity_channels(1)=&T_STRUCT::ch0;face_velocity_channels(2)=&T_STRUCT::ch1;
        int level=((OPENGL_SELECTION_SPGRID_CELL_2D<T>*)selection)->level;
        VECTOR<int,2> index=((OPENGL_SELECTION_SPGRID_CELL_2D<T>*)selection)->index;
        const unsigned long offset=Flag_array_type::MASK::Linear_Offset(std_array<int,2>(index));
        const unsigned long other_x_offset=Flag_array_type::MASK::Packed_Add(offset,neighbor_offsets(1));
        const unsigned long other_y_offset=Flag_array_type::MASK::Packed_Add(offset,neighbor_offsets(2));
        if(opengl_spgrid.hierarchy->Set(level).Is_Set(offset,GRID_TOPOLOGY_HELPER<typename Flag_array_type::MASK>::Face_Valid_Mask(1)))
            output_stream<<"    u left = "<<opengl_spgrid.hierarchy->Array(level,face_velocity_channels(1))(offset);
        else output_stream<<"    u left = N/A";
        if(opengl_spgrid.hierarchy->Set(level).Is_Set(other_x_offset,GRID_TOPOLOGY_HELPER<typename Flag_array_type::MASK>::Face_Valid_Mask(1)))
            output_stream<<", right = "<<opengl_spgrid.hierarchy->Array(level,face_velocity_channels(1))(other_x_offset)<<std::endl;
        else output_stream<<", right = N/A"<<std::endl;
        if(opengl_spgrid.hierarchy->Set(level).Is_Set(offset,GRID_TOPOLOGY_HELPER<typename Flag_array_type::MASK>::Face_Valid_Mask(2)))
            output_stream<<"    v bottom = "<<opengl_spgrid.hierarchy->Array(level,face_velocity_channels(2))(offset);
        else output_stream<<"    v bottom = N/A";
        if(opengl_spgrid.hierarchy->Set(level).Is_Set(other_y_offset,GRID_TOPOLOGY_HELPER<typename Flag_array_type::MASK>::Face_Valid_Mask(2)))
            output_stream<<", top = "<<opengl_spgrid.hierarchy->Array(level,face_velocity_channels(2))(other_y_offset)<<std::endl;
        else output_stream<<", top = N/A"<<std::endl;
        output_stream<<"    Divergence: "<<opengl_spgrid.hierarchy->Array(level,divergence_channel)(index[1],index[2])<<std::endl;}
}
//#####################################################################
template class OPENGL_SPGRID_MAC_VELOCITY_FIELD_2D<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class OPENGL_SPGRID_MAC_VELOCITY_FIELD_2D<double>;
#endif
