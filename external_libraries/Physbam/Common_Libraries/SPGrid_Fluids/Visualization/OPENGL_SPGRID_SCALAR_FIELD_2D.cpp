//#####################################################################
// Copyright (c) 2014, Mridul Aanjaneya.
// Distributed under the FreeBSD license (see license.txt)
//#####################################################################
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_COLOR_RAMP.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_TEXTURED_RECT.h>
#include <SPGrid_Fluids/Simulation/FLUIDS_SIMULATION_FLAGS.h>
#include <SPGrid_Fluids/Visualization/OPENGL_SPGRID_2D.h>
#include <SPGrid_Fluids/Visualization/OPENGL_SPGRID_SCALAR_FIELD_2D.h>
using namespace PhysBAM;
//#####################################################################
// Function Constructor
//#####################################################################
template<class T,class T2> OPENGL_SPGRID_SCALAR_FIELD_2D<T,T2>::
OPENGL_SPGRID_SCALAR_FIELD_2D(GRID<TV> &fine_mac_grid_input,OPENGL_SPGRID_2D<T>& opengl_spgrid_input,OPENGL_COLOR_MAP<T2>* color_map_input,bool& draw_separate_levels_input,int& level_counter_input,T T_STRUCT::* field_input,DRAW_MODE draw_mode_input)
    :fine_mac_grid(fine_mac_grid_input),opengl_spgrid(opengl_spgrid_input),draw_separate_levels(draw_separate_levels_input),level_counter(level_counter_input),field(field_input),draw_mode(draw_mode_input),current_color_map(1),opengl_textured_rect(0),scale_range(false)
{
    PHYSBAM_ASSERT(color_map_input);
    Initialize_Color_Maps(color_map_input);
    Set_Draw_Mode(draw_mode_input);
}
//#####################################################################
// Function Destructor
//#####################################################################
template<class T,class T2> OPENGL_SPGRID_SCALAR_FIELD_2D<T,T2>::
~OPENGL_SPGRID_SCALAR_FIELD_2D()
{
    if(opengl_textured_rect){
        delete opengl_textured_rect->texture;
        delete opengl_textured_rect;}
    color_maps.Delete_Pointers_And_Clean_Memory();
}
//#####################################################################
// Function Initialize_Color_Maps
//#####################################################################
template<> void OPENGL_SPGRID_SCALAR_FIELD_2D<float,bool>::
Initialize_Color_Maps(OPENGL_COLOR_MAP<bool>* color_map_input)
{
    color_maps.Append(color_map_input);
}
//#####################################################################
// Function Initialize_Color_Maps
//#####################################################################
template<> void OPENGL_SPGRID_SCALAR_FIELD_2D<double,bool>::
Initialize_Color_Maps(OPENGL_COLOR_MAP<bool>* color_map_input)
{
    color_maps.Append(color_map_input);
}
//#####################################################################
// Function Initialize_Color_Maps
//#####################################################################
template<> void OPENGL_SPGRID_SCALAR_FIELD_2D<float,int>::
Initialize_Color_Maps(OPENGL_COLOR_MAP<int>* color_map_input)
{
    color_maps.Append(color_map_input);
}
//#####################################################################
// Function Initialize_Color_Maps
//#####################################################################
template<> void OPENGL_SPGRID_SCALAR_FIELD_2D<double,int>::
Initialize_Color_Maps(OPENGL_COLOR_MAP<int>* color_map_input)
{
    color_maps.Append(color_map_input);
}
//#####################################################################
// Function Initialize_Color_Maps
//#####################################################################
template<class T,class T2> void OPENGL_SPGRID_SCALAR_FIELD_2D<T,T2>::
Initialize_Color_Maps(OPENGL_COLOR_MAP<T2>* color_map_input)
{
    color_maps.Append(color_map_input);
    color_maps.Append(OPENGL_COLOR_RAMP<T2>::Matlab_Jet(0,1));
    color_maps.Append(OPENGL_COLOR_RAMP<T2>::Matlab_Hot(0,1));
}
//#####################################################################
// Function Set_Scale_Range
//#####################################################################
template<> void OPENGL_SPGRID_SCALAR_FIELD_2D<float,bool>::
Set_Scale_Range(const bool range_min,const bool range_max)
{PHYSBAM_FATAL_ERROR();}
//#####################################################################
// Function Set_Scale_Range
//#####################################################################
template<> void OPENGL_SPGRID_SCALAR_FIELD_2D<double,bool>::
Set_Scale_Range(const bool range_min,const bool range_max)
{PHYSBAM_FATAL_ERROR();}
//#####################################################################
// Function Set_Scale_Range
//#####################################################################
template<> void OPENGL_SPGRID_SCALAR_FIELD_2D<float,int>::
Set_Scale_Range(const int range_min,const int range_max)
{PHYSBAM_FATAL_ERROR();}
//#####################################################################
// Function Set_Scale_Range
//#####################################################################
template<> void OPENGL_SPGRID_SCALAR_FIELD_2D<double,int>::
Set_Scale_Range(const int range_min,const int range_max)
{PHYSBAM_FATAL_ERROR();}
//#####################################################################
// Function Set_Scale_Range
//#####################################################################
template<class T,class T2> void OPENGL_SPGRID_SCALAR_FIELD_2D<T,T2>::
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
template<class T,class T2> void OPENGL_SPGRID_SCALAR_FIELD_2D<T,T2>::
Reset_Scale_Range()
{
    scale_range=false;
}
//#####################################################################
// Function Pre_Map_Value
//#####################################################################
template<> bool OPENGL_SPGRID_SCALAR_FIELD_2D<float,bool>::
Pre_Map_Value(const bool value) const
{
    return value;
}
//#####################################################################
// Function Pre_Map_Value
//#####################################################################
template<> bool OPENGL_SPGRID_SCALAR_FIELD_2D<double,bool>::
Pre_Map_Value(const bool value) const
{
    return value;
}
//#####################################################################
// Function Pre_Map_Value
//#####################################################################
template<class T,class T2> T2 OPENGL_SPGRID_SCALAR_FIELD_2D<T,T2>::
Pre_Map_Value(const T2 value) const
{
    if(!scale_range) return value;
    else return (value-scale_range_min)*scale_range_dx; 
}
//#####################################################################
// Function Display_Points
//#####################################################################
template<class T,class T2> void OPENGL_SPGRID_SCALAR_FIELD_2D<T,T2>::
Display_Points(const int in_color) const
{
    glPushAttrib(GL_LIGHTING_BIT | GL_TEXTURE_BIT | GL_POINT_BIT);
    glPointSize(5);glDisable(GL_LIGHTING);glDisable(GL_TEXTURE_2D);

    OpenGL_Begin(GL_POINTS);
    for(int level=1;level<=opengl_spgrid.hierarchy->Levels();++level) if(!draw_separate_levels || (draw_separate_levels && level==level_counter)){
        Const_data_array_type spgrid_density=opengl_spgrid.hierarchy->Allocator(level).Get_Const_Array(field);
        Const_flag_array_type flags=opengl_spgrid.hierarchy->Allocator(level).Get_Const_Array(&T_STRUCT::flags);    
        for(SPGrid_Block_Iterator<typename Flag_array_type::MASK> iterator(opengl_spgrid.hierarchy->Blocks(level));iterator.Valid();iterator.Next()){
            TV_INT index=iterator.Index().template Cast<TV_INT>();
            unsigned cell_flags=iterator.Data(flags);
            if(cell_flags&(SPGrid_Cell_Type_Interior)){
                color_maps(current_color_map)->Lookup(iterator.Data(spgrid_density)).Send_To_GL_Pipeline();
                OpenGL_Vertex(opengl_spgrid.hierarchy->Grid(level).Get_MAC_Grid().X(index));}}}
    OpenGL_End();
    glPopAttrib();
}
//#####################################################################
// Function Display
//#####################################################################
template<class T,class T2> void OPENGL_SPGRID_SCALAR_FIELD_2D<T,T2>::
Display(const int in_color) const
{ 
    if(draw_mode==DRAW_TEXTURE){
        PHYSBAM_ASSERT(opengl_textured_rect);
        Update();
        opengl_textured_rect->Display(in_color);}
    else if(draw_mode==DRAW_POINTS) Display_Points(in_color);
}
//#####################################################################
// Function Bounding_Box
//#####################################################################
template<class T,class T2> RANGE<VECTOR<float,3> > OPENGL_SPGRID_SCALAR_FIELD_2D<T,T2>::
Bounding_Box() const
{
    // May not be the exact bounds, but close enough...
    return World_Space_Box(RANGE<VECTOR<float,3> >(fine_mac_grid.domain.min_corner.x,fine_mac_grid.domain.max_corner.x,fine_mac_grid.domain.min_corner.y,fine_mac_grid.domain.max_corner.y,0,0));
}
//#####################################################################
// Function Set_Draw_Mode
//#####################################################################
template<class T,class T2> void OPENGL_SPGRID_SCALAR_FIELD_2D<T,T2>::
Set_Draw_Mode(DRAW_MODE draw_mode_input)
{
    draw_mode = draw_mode_input;

    if(draw_mode!=DRAW_TEXTURE && opengl_textured_rect){
        delete opengl_textured_rect->texture;delete opengl_textured_rect;opengl_textured_rect=0;}

    if(draw_mode==DRAW_TEXTURE){
        if(!opengl_textured_rect) opengl_textured_rect=new OPENGL_TEXTURED_RECT;}

}
//#####################################################################
// Function Update
//#####################################################################
template<class T,class T2> void OPENGL_SPGRID_SCALAR_FIELD_2D<T,T2>::
Update() const
{
    TV_INT start_index(1,1);TV_INT end_index=fine_mac_grid.Domain_Indices().max_corner;
    if(draw_mode==DRAW_TEXTURE) Update_Texture(start_index,end_index);
}
//#####################################################################
// Function Print_Selection_Info
//#####################################################################
template<class T,class T2> void OPENGL_SPGRID_SCALAR_FIELD_2D<T,T2>::
Print_Selection_Info(std::ostream& output_stream,OPENGL_SELECTION* current_selection) const
{
    if(current_selection && current_selection->type==OPENGL_SELECTION::GRID_CELL_2D){
        int level=((OPENGL_SELECTION_SPGRID_CELL_2D<T>*)current_selection)->level;
        VECTOR<int,2> index=((OPENGL_SELECTION_SPGRID_CELL_2D<T>*)current_selection)->index;
        output_stream<<opengl_spgrid.hierarchy->Array(level,field)(index[1],index[2])<<std::endl;}
}
//#####################################################################
// Function Update_Texture
//#####################################################################
template<class T,class T2> void OPENGL_SPGRID_SCALAR_FIELD_2D<T,T2>::
Update_Texture(const VECTOR<int,2>& start_index,const VECTOR<int,2>& end_index) const
{
    PHYSBAM_ASSERT(opengl_textured_rect);

    if ((end_index.x-start_index.x+1) == 0 || (end_index.y-start_index.y+1) == 0)
    {
        delete opengl_textured_rect->texture; opengl_textured_rect->texture = 0; return;
    }

    // Handle values arrays which are not (1,m)(1,n)
    VECTOR<T,2> half_dX=(T)0.5*fine_mac_grid.dX;
    RANGE<VECTOR<T,2> > domain(fine_mac_grid.X(start_index)-half_dX,fine_mac_grid.X(end_index)+half_dX);

    // Set underlying OPENGL_OBJECT's transformation
    opengl_textured_rect->frame->t=VECTOR<float,3>(Convert_2d_To_3d(domain.Center()));

    // Set OPENGL_TEXTURED_RECT's dimensions
    opengl_textured_rect->width=domain.Edge_Lengths().x;
    opengl_textured_rect->height=domain.Edge_Lengths().y;

    int tex_width = end_index.x-start_index.x+1;
    int tex_height = end_index.y-start_index.y+1;

    if (!opengl_textured_rect->texture || opengl_textured_rect->texture->width != tex_width || opengl_textured_rect->texture->height != tex_height)
    {
        delete opengl_textured_rect->texture;
        opengl_textured_rect->texture = new OPENGL_TEXTURE();
        opengl_textured_rect->texture->Initialize(tex_width, tex_height);
    }

    OPENGL_COLOR* bitmap = new OPENGL_COLOR[tex_width*tex_height];
    OPENGL_COLOR_MAP<T2>* color_map=color_maps(current_color_map);
    for(int i=1;i<=tex_width;++i)for(int j=1;j<=tex_height;++j){
        int idx=(j-1)*tex_width+i-1;
        TV_INT index=TV_INT(start_index.x+i-1,start_index.y+j-1);
        unsigned long offset=Flag_array_type::MASK::Linear_Offset(std_array<int,2>(index));
        int level=0;
        for(int l=1;l<=opengl_spgrid.hierarchy->Levels();++l){
            if(opengl_spgrid.hierarchy->Set(l).Is_Set(offset,SPGrid_Cell_Type_Interior)){level=l;break;}
            offset=Flag_array_type::MASK::DownsampleOffset(offset);}

        T2 value=(level>0 && (!draw_separate_levels || (draw_separate_levels && level==level_counter)))?Pre_Map_Value(opengl_spgrid.hierarchy->Array(level,field)(offset)):(T2)0;
        OPENGL_COLOR color_value=color_map->Lookup(value);
        bitmap[idx]=color_value;}

    opengl_textured_rect->texture->Update_Texture(bitmap);
    delete[] bitmap;
}
//#####################################################################
// Function Toggle_Draw_Mode
//#####################################################################
template<class T,class T2> void OPENGL_SPGRID_SCALAR_FIELD_2D<T,T2>::
Toggle_Draw_Mode()
{
    DRAW_MODE new_draw_mode=(DRAW_MODE)(((int)draw_mode+1)%2);
    Set_Draw_Mode(new_draw_mode);
}
//#####################################################################
// Function Toggle_Color_Map
//#####################################################################
template<class T,class T2> void OPENGL_SPGRID_SCALAR_FIELD_2D<T,T2>::
Toggle_Color_Map()
{
    current_color_map=current_color_map%color_maps.m+1;
}
//#####################################################################
template class OPENGL_SPGRID_SCALAR_FIELD_2D<float,int>;
template class OPENGL_SPGRID_SCALAR_FIELD_2D<float,bool>;
template class OPENGL_SPGRID_SCALAR_FIELD_2D<float,float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class OPENGL_SPGRID_SCALAR_FIELD_2D<double,int>;
template class OPENGL_SPGRID_SCALAR_FIELD_2D<double,bool>;
template class OPENGL_SPGRID_SCALAR_FIELD_2D<double,double>;
#endif
