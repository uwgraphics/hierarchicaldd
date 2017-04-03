//#####################################################################
// Copyright (c) 2014, Mridul Aanjaneya.
// Distributed under the FreeBSD license (see license.txt)
//#####################################################################
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_COLOR_RAMP.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_UNIFORM_SLICE.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_WORLD.h>
#include <SPGrid_Fluids/Simulation/FLUIDS_SIMULATION_FLAGS.h>
#include <SPGrid_Fluids/Rendering/OpenGL/OPENGL_SPGRID_3D.h>
#include <SPGrid_Fluids/Rendering/OpenGL/OPENGL_SPGRID_SCALAR_FIELD_3D.h>
using namespace PhysBAM;
//#####################################################################
// Function Constructor
//#####################################################################
template<class T,class T2> OPENGL_SPGRID_SCALAR_FIELD_3D<T,T2>::
OPENGL_SPGRID_SCALAR_FIELD_3D(GRID<TV> &fine_mac_grid_input,OPENGL_SPGRID_3D<T>& opengl_spgrid_input,OPENGL_COLOR_MAP<T2>* color_map_input,bool& draw_separate_levels_input,int& level_counter_input,T T_STRUCT::* field_input,DRAW_MODE draw_mode_input)
    :fine_mac_grid(fine_mac_grid_input),opengl_spgrid(opengl_spgrid_input),draw_separate_levels(draw_separate_levels_input),level_counter(level_counter_input),field(field_input),draw_mode(draw_mode_input),current_color_map(1),scale_range(false)
{
    PHYSBAM_ASSERT(color_map_input);
    Initialize_Color_Maps(color_map_input);
    Set_Draw_Mode(draw_mode_input);
}
//#####################################################################
// Function Destructor
//#####################################################################
template<class T,class T2> OPENGL_SPGRID_SCALAR_FIELD_3D<T,T2>::
~OPENGL_SPGRID_SCALAR_FIELD_3D()
{
    color_maps.Delete_Pointers_And_Clean_Memory();
}
//#####################################################################
// Function Initialize_Color_Maps
//#####################################################################
template<> void OPENGL_SPGRID_SCALAR_FIELD_3D<float,bool>::
Initialize_Color_Maps(OPENGL_COLOR_MAP<bool>* color_map_input)
{
    color_maps.Append(color_map_input);
}
//#####################################################################
// Function Initialize_Color_Maps
//#####################################################################
template<> void OPENGL_SPGRID_SCALAR_FIELD_3D<double,bool>::
Initialize_Color_Maps(OPENGL_COLOR_MAP<bool>* color_map_input)
{
    color_maps.Append(color_map_input);
}
//#####################################################################
// Function Initialize_Color_Maps
//#####################################################################
template<> void OPENGL_SPGRID_SCALAR_FIELD_3D<float,int>::
Initialize_Color_Maps(OPENGL_COLOR_MAP<int>* color_map_input)
{
    color_maps.Append(color_map_input);
}
//#####################################################################
// Function Initialize_Color_Maps
//#####################################################################
template<> void OPENGL_SPGRID_SCALAR_FIELD_3D<double,int>::
Initialize_Color_Maps(OPENGL_COLOR_MAP<int>* color_map_input)
{
    color_maps.Append(color_map_input);
}
//#####################################################################
// Function Initialize_Color_Maps
//#####################################################################
template<class T,class T2> void OPENGL_SPGRID_SCALAR_FIELD_3D<T,T2>::
Initialize_Color_Maps(OPENGL_COLOR_MAP<T2>* color_map_input)
{
    color_maps.Append(color_map_input);
    color_maps.Append(OPENGL_COLOR_RAMP<T2>::Matlab_Jet(8e4,1e6));
    color_maps.Append(OPENGL_COLOR_RAMP<T2>::Matlab_Hot(8e4,1e6));
}
//#####################################################################
// Function Set_Scale_Range
//#####################################################################
template<> void OPENGL_SPGRID_SCALAR_FIELD_3D<float,bool>::
Set_Scale_Range(const bool range_min,const bool range_max)
{PHYSBAM_FATAL_ERROR();}
//#####################################################################
// Function Set_Scale_Range
//#####################################################################
template<> void OPENGL_SPGRID_SCALAR_FIELD_3D<double,bool>::
Set_Scale_Range(const bool range_min,const bool range_max)
{PHYSBAM_FATAL_ERROR();}
//#####################################################################
// Function Set_Scale_Range
//#####################################################################
template<> void OPENGL_SPGRID_SCALAR_FIELD_3D<float,int>::
Set_Scale_Range(const int range_min,const int range_max)
{PHYSBAM_FATAL_ERROR();}
//#####################################################################
// Function Set_Scale_Range
//#####################################################################
template<> void OPENGL_SPGRID_SCALAR_FIELD_3D<double,int>::
Set_Scale_Range(const int range_min,const int range_max)
{PHYSBAM_FATAL_ERROR();}
//#####################################################################
// Function Set_Scale_Range
//#####################################################################
template<class T,class T2> void OPENGL_SPGRID_SCALAR_FIELD_3D<T,T2>::
Set_Scale_Range(const T2 range_min,const T2 range_max)
{
    scale_range=true;
    scale_range_min=range_min;
    T2 range_length=(range_max-range_min);
    scale_range_dx=range_length>1e-10?(T2)1/range_length:(T2)0;
}
//#####################################################################
// Function Reset_Scale_Range
//#####################################################################
template<class T,class T2> void OPENGL_SPGRID_SCALAR_FIELD_3D<T,T2>::
Reset_Scale_Range()
{
    scale_range=false;
}
//#####################################################################
// Function Pre_Map_Value
//#####################################################################
template<> bool OPENGL_SPGRID_SCALAR_FIELD_3D<float,bool>::
Pre_Map_Value(const bool value) const
{
    return value;
}
//#####################################################################
// Function Pre_Map_Value
//#####################################################################
template<> bool OPENGL_SPGRID_SCALAR_FIELD_3D<double,bool>::
Pre_Map_Value(const bool value) const
{
    return value;
}
//#####################################################################
// Function Pre_Map_Value
//#####################################################################
template<class T,class T2> T2 OPENGL_SPGRID_SCALAR_FIELD_3D<T,T2>::
Pre_Map_Value(const T2 value) const
{
    if(!scale_range) return value;
    else return (value-scale_range_min)*scale_range_dx; 
}
//#####################################################################
// Function Display_Points
//#####################################################################
template<class T,class T2> void OPENGL_SPGRID_SCALAR_FIELD_3D<T,T2>::
Display_Points(const int in_color) const
{
    glPushAttrib(GL_LIGHTING_BIT | GL_TEXTURE_BIT | GL_POINT_BIT);
    glPointSize(5);glDisable(GL_LIGHTING);glDisable(GL_TEXTURE_2D);
    OPENGL_UNIFORM_SLICE* slice=(OPENGL_UNIFORM_SLICE*)this->slice;
    int axis;TV location;
    if(slice){axis=slice->axis;
        location=opengl_spgrid.hierarchy->Grid(1).Get_MAC_Grid().domain.min_corner+opengl_spgrid.hierarchy->Grid(1).Get_MAC_Grid().dX*(TV(slice->index)-TV::All_Ones_Vector()*(T).5);}

    OpenGL_Begin(GL_POINTS);
    for(int level=1;level<=opengl_spgrid.hierarchy->Levels();++level) if(!draw_separate_levels || (draw_separate_levels && level==level_counter)){
        const TV& current_dX=opengl_spgrid.hierarchy->Grid(level).dX;
        Const_data_array_type spgrid_density=opengl_spgrid.hierarchy->Allocator(level).Get_Const_Array(field);
        Const_flag_array_type flags=opengl_spgrid.hierarchy->Allocator(level).Get_Const_Array(&T_STRUCT::flags);    
        for(SPGrid_Block_Iterator<typename Flag_array_type::MASK> iterator(opengl_spgrid.hierarchy->Blocks(level));iterator.Valid();iterator.Next()){
            TV_INT index=iterator.Index().template Cast<TV_INT>();
            TV left_corner(index-TV_INT::All_Ones_Vector());
            unsigned cell_flags=iterator.Data(flags);
            if(cell_flags&(SPGrid_Cell_Type_Interior)){
                if(!slice || ((slice->mode == OPENGL_SLICE::CELL_SLICE && location[axis]>=left_corner[axis]*current_dX[axis] && location[axis]<=((left_corner+TV::All_Ones_Vector())*current_dX)[axis]) || slice->mode != OPENGL_SLICE::CELL_SLICE)){
                    color_maps(current_color_map)->Lookup(iterator.Data(spgrid_density)).Send_To_GL_Pipeline();
                    OpenGL_Vertex(opengl_spgrid.hierarchy->Grid(level).Get_MAC_Grid().X(index));}}}}
    OpenGL_End();
    glPopAttrib();
}
//#####################################################################
// Function Display
//#####################################################################
template<class T,class T2> void OPENGL_SPGRID_SCALAR_FIELD_3D<T,T2>::
Display(const int in_color) const
{ 
    if(draw_mode==DRAW_TEXTURE) Display_Texture();
    else if(draw_mode==DRAW_POINTS) Display_Points(in_color);
}
//#####################################################################
// Function Display_Texture
//#####################################################################
template<class T,class T2> void OPENGL_SPGRID_SCALAR_FIELD_3D<T,T2>::
Display_Texture() const
{
    glPushAttrib(GL_ENABLE_BIT|GL_DEPTH_BUFFER_BIT|GL_CURRENT_BIT);
    glDisable(GL_LIGHTING);
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_BLEND);
    glDepthMask(GL_FALSE);

    VECTOR<float,3> view_forward,view_up,view_right;
    OPENGL_WORLD::Singleton()->Get_View_Frame(view_forward,view_up,view_right);
    int dominant_axis=view_forward.Dominant_Axis();

    OPENGL_UNIFORM_SLICE* slice=(OPENGL_UNIFORM_SLICE*)this->slice;
    int axis;TV location;
    if(slice){axis=slice->axis;
        location=opengl_spgrid.hierarchy->Grid(1).Get_MAC_Grid().domain.min_corner+opengl_spgrid.hierarchy->Grid(1).Get_MAC_Grid().dX*(TV(slice->index)-TV::All_Ones_Vector()*(T).5);}

    OpenGL_Begin(GL_QUADS);
    for(int level=1;level<=opengl_spgrid.hierarchy->Levels();++level) if(!draw_separate_levels || (draw_separate_levels && level==level_counter)){
        const TV& current_dX=opengl_spgrid.hierarchy->Grid(level).dX;
        Const_data_array_type spgrid_density=opengl_spgrid.hierarchy->Allocator(level).Get_Const_Array(field);
        Const_flag_array_type flags=opengl_spgrid.hierarchy->Allocator(level).Get_Const_Array(&T_STRUCT::flags);    
        for(SPGrid_Block_Iterator<typename Flag_array_type::MASK> iterator(opengl_spgrid.hierarchy->Blocks(level));iterator.Valid();iterator.Next()){
            TV_INT index=iterator.Index().template Cast<TV_INT>();
            TV left_corner(index-TV_INT::All_Ones_Vector());
            TV pos=opengl_spgrid.hierarchy->Grid(level).Get_MAC_Grid().X(index);
            unsigned cell_flags=iterator.Data(flags);
            if(cell_flags&(SPGrid_Cell_Type_Interior)){
                if(!slice || ((slice->mode == OPENGL_SLICE::CELL_SLICE && location[axis]>=left_corner[axis]*current_dX[axis] && location[axis]<=((left_corner+TV::All_Ones_Vector())*current_dX)[axis]) || slice->mode != OPENGL_SLICE::CELL_SLICE)){
                    color_maps(current_color_map)->Lookup(iterator.Data(spgrid_density)).Send_To_GL_Pipeline();
                    if(dominant_axis==1){
                        if(view_forward[1]>0){
                            glVertex3f(pos.x,pos.y-0.5*current_dX.y,pos.z-0.5*current_dX.z);
                            glVertex3f(pos.x,pos.y-0.5*current_dX.y,pos.z+0.5*current_dX.z);
                            glVertex3f(pos.x,pos.y+0.5*current_dX.y,pos.z+0.5*current_dX.z);
                            glVertex3f(pos.x,pos.y+0.5*current_dX.y,pos.z-0.5*current_dX.z);}
                        else{
                            glVertex3f(pos.x,pos.y-0.5*current_dX.y,pos.z+0.5*current_dX.z);
                            glVertex3f(pos.x,pos.y-0.5*current_dX.y,pos.z-0.5*current_dX.z);
                            glVertex3f(pos.x,pos.y+0.5*current_dX.y,pos.z-0.5*current_dX.z);
                            glVertex3f(pos.x,pos.y+0.5*current_dX.y,pos.z+0.5*current_dX.z);}}
                    else if(dominant_axis==2){
                        if(view_forward[2]>0){
                            glVertex3f(pos.x-0.5*current_dX.x,pos.y,pos.z-0.5*current_dX.z);
                            glVertex3f(pos.x+0.5*current_dX.x,pos.y,pos.z-0.5*current_dX.z);
                            glVertex3f(pos.x+0.5*current_dX.x,pos.y,pos.z+0.5*current_dX.z);
                            glVertex3f(pos.x-0.5*current_dX.x,pos.y,pos.z+0.5*current_dX.z);}
                        else{
                           glVertex3f(pos.x-0.5*current_dX.x,pos.y,pos.z+0.5*current_dX.z);
                           glVertex3f(pos.x+0.5*current_dX.x,pos.y,pos.z+0.5*current_dX.z);
                           glVertex3f(pos.x+0.5*current_dX.x,pos.y,pos.z-0.5*current_dX.z);
                           glVertex3f(pos.x-0.5*current_dX.x,pos.y,pos.z-0.5*current_dX.z);}}
                    else if(dominant_axis==3){
                        if(view_forward[3]>0){
                            glVertex3f(pos.x-0.5*current_dX.x,pos.y-0.5*current_dX.y,pos.z);
                            glVertex3f(pos.x-0.5*current_dX.x,pos.y+0.5*current_dX.y,pos.z);
                            glVertex3f(pos.x+0.5*current_dX.x,pos.y+0.5*current_dX.y,pos.z);
                            glVertex3f(pos.x+0.5*current_dX.x,pos.y-0.5*current_dX.y,pos.z);}
                        else{
                            glVertex3f(pos.x-0.5*current_dX.x,pos.y+0.5*current_dX.y,pos.z);
                            glVertex3f(pos.x-0.5*current_dX.x,pos.y-0.5*current_dX.y,pos.z);
                            glVertex3f(pos.x+0.5*current_dX.x,pos.y-0.5*current_dX.y,pos.z);
                            glVertex3f(pos.x+0.5*current_dX.x,pos.y+0.5*current_dX.y,pos.z);}}}}}}
    OpenGL_End();
    glPopAttrib();
}
//#####################################################################
// Function Bounding_Box
//#####################################################################
template<class T,class T2> RANGE<VECTOR<float,3> > OPENGL_SPGRID_SCALAR_FIELD_3D<T,T2>::
Bounding_Box() const
{
    // May not be the exact bounds, but close enough...
    return World_Space_Box(RANGE<VECTOR<float,3> >(fine_mac_grid.domain.min_corner.x,fine_mac_grid.domain.max_corner.x,fine_mac_grid.domain.min_corner.y,fine_mac_grid.domain.max_corner.y,0,0));
}
//#####################################################################
// Function Set_Draw_Mode
//#####################################################################
template<class T,class T2> void OPENGL_SPGRID_SCALAR_FIELD_3D<T,T2>::
Set_Draw_Mode(DRAW_MODE draw_mode_input)
{
    draw_mode = draw_mode_input;
}
//#####################################################################
// Function Print_Selection_Info
//#####################################################################
template<class T,class T2> void OPENGL_SPGRID_SCALAR_FIELD_3D<T,T2>::
Print_Selection_Info(std::ostream& output_stream,OPENGL_SELECTION* current_selection) const
{
    if(current_selection && current_selection->type==OPENGL_SELECTION::GRID_CELL_3D){
        int level=((OPENGL_SELECTION_SPGRID_CELL_3D<T>*)current_selection)->level;
        TV_INT index=((OPENGL_SELECTION_SPGRID_CELL_3D<T>*)current_selection)->index;
        output_stream<<opengl_spgrid.hierarchy->Array(level,field)(index[1],index[2],index[3])<<std::endl;}
}
//#####################################################################
// Function Toggle_Draw_Mode
//#####################################################################
template<class T,class T2> void OPENGL_SPGRID_SCALAR_FIELD_3D<T,T2>::
Toggle_Draw_Mode()
{
    DRAW_MODE new_draw_mode=(DRAW_MODE)(((int)draw_mode+1)%2);
    Set_Draw_Mode(new_draw_mode);
}
//#####################################################################
// Function Toggle_Color_Map
//#####################################################################
template<class T,class T2> void OPENGL_SPGRID_SCALAR_FIELD_3D<T,T2>::
Toggle_Color_Map()
{
    current_color_map=current_color_map%color_maps.m+1;
}
//#####################################################################
template class OPENGL_SPGRID_SCALAR_FIELD_3D<float,int>;
template class OPENGL_SPGRID_SCALAR_FIELD_3D<float,bool>;
template class OPENGL_SPGRID_SCALAR_FIELD_3D<float,float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class OPENGL_SPGRID_SCALAR_FIELD_3D<double,int>;
template class OPENGL_SPGRID_SCALAR_FIELD_3D<double,bool>;
template class OPENGL_SPGRID_SCALAR_FIELD_3D<double,double>;
#endif
