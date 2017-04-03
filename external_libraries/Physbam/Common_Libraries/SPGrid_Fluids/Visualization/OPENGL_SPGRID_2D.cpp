//#####################################################################
// Copyright (c) 2014, Mridul Aanjaneya.
// Distributed under the FreeBSD license (see license.txt)
//#####################################################################
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_CONSTANT_COLOR_MAP.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_PREFERENCES.h>
#include <SPGrid/Tools/SPGrid_Block_Iterator.h>
#include <SPGrid_Fluids/Grids/GRID_TOPOLOGY_HELPER.h>
#include <SPGrid_Fluids/Simulation/FLUIDS_SIMULATION_FLAGS.h>
#include <SPGrid_Fluids/Visualization/OPENGL_SPGRID_2D.h>
using namespace PhysBAM;
//#####################################################################
// Function Display
//#####################################################################
template<class T> void OPENGL_SPGRID_2D<T>::
Display(const int in_color) const
{
    if(!draw)return;
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    Send_Transform_To_GL_Pipeline();

    glPushAttrib(GL_ENABLE_BIT);
    glDisable(GL_LIGHTING);
    glDisable(GL_DEPTH_TEST);

    color.Send_To_GL_Pipeline();

    GLint mode;
    glGetIntegerv(GL_RENDER_MODE, &mode);

    if(mode == GL_SELECT){
        glPushAttrib(GL_ENABLE_BIT | GL_POINT_BIT);
        glDisable(GL_CULL_FACE);

        const int levels=hierarchy->Levels();
        RANGE<TV> domain=hierarchy->Grid(levels).domain;
        const TV epsilon=(T).5*hierarchy->Grid(levels).dX*T();

        // Draw cells for selection
        glPushName(1);
        for(int level=levels;level>=1;--level) if(!draw_separate_levels || (draw_separate_levels && level==level_counter)){
            glPushName(level);
            Const_flag_array_type flags=hierarchy->Allocator(level).Get_Const_Array(&T_STRUCT::flags);
            for(SPGrid_Block_Iterator<typename Flag_array_type::MASK> iterator(hierarchy->Blocks(level));iterator.Valid();iterator.Next()){
                TV_INT index=iterator.Index().template Cast<TV_INT>();
                unsigned cell_flags=iterator.Data(flags);
                if(cell_flags&(SPGrid_Cell_Type_Interior) || cell_flags&(SPGrid_Cell_Type_Ghost)){
                    const TV& current_dX=hierarchy->Grid(level).dX;const RANGE<TV>& current_domain=hierarchy->Grid(level).domain;
                    TV left_corner(index-TV_INT::All_Ones_Vector());
                    glPushName(index[1]);glPushName(index[2]);
                    OpenGL_Begin(GL_QUADS);
                    OpenGL_Quad_2D(current_domain.min_corner+left_corner*current_dX,current_domain.min_corner+(left_corner+TV::All_Ones_Vector())*current_dX);
                    OpenGL_End();
                    glPopName();glPopName();}}
            glPopName();}
        glPopName();glPopAttrib();}
    else{const int levels=hierarchy->Levels();
        RANGE<TV> domain=hierarchy->Grid(levels).domain;
        const TV epsilon=(T).5*hierarchy->Grid(levels).dX*T();

        OpenGL_Begin(GL_LINES);
        for(int level=levels;level>=1;--level) if(!draw_separate_levels || (draw_separate_levels && level==level_counter)){
            Const_flag_array_type flags=hierarchy->Allocator(level).Get_Const_Array(&T_STRUCT::flags);
            for(SPGrid_Block_Iterator<typename Flag_array_type::MASK> iterator(hierarchy->Blocks(level));iterator.Valid();iterator.Next()){
                TV_INT index=iterator.Index().template Cast<TV_INT>();
                unsigned cell_flags=iterator.Data(flags);
                if(cell_flags&(SPGrid_Cell_Type_Interior)){
                    const TV& current_dX=hierarchy->Grid(level).dX;const RANGE<TV>& current_domain=hierarchy->Grid(level).domain;
                    for(int axis=1;axis<=2;++axis){
                        int other_axis=axis%2+1;
                        TV left_corner(index-TV_INT::All_Ones_Vector());
                        OpenGL_Line<T,2>(current_domain.min_corner+left_corner*current_dX,current_domain.min_corner+(left_corner+TV::Axis_Vector(axis))*current_dX);
                        OpenGL_Line<T,2>(current_domain.min_corner+(left_corner+TV::Axis_Vector(other_axis))*current_dX,current_domain.min_corner+(left_corner+TV::Axis_Vector(other_axis)+TV::Axis_Vector(axis))*current_dX);}}}}
        OpenGL_End();

        if(draw_separate_levels){
            glPushAttrib(GL_LINE_BIT);glLineWidth(2*OPENGL_PREFERENCES::line_width);
            OpenGL_Begin(GL_LINES);
            int level=level_counter;
            Const_flag_array_type flags=hierarchy->Allocator(level).Get_Const_Array(&T_STRUCT::flags);
            for(SPGrid_Block_Iterator<typename Flag_array_type::MASK> iterator(hierarchy->Blocks(level));iterator.Valid();iterator.Next()){
                TV_INT index=iterator.Index().template Cast<TV_INT>();
                unsigned cell_flags=iterator.Data(flags);
                if(cell_flags&(SPGrid_Cell_Type_Ghost)){
                    const TV& current_dX=hierarchy->Grid(level).dX;const RANGE<TV>& current_domain=hierarchy->Grid(level).domain;
                    for(int axis=1;axis<=2;++axis){
                        int other_axis=axis%2+1;
                        TV left_corner(index-TV_INT::All_Ones_Vector());
                        OpenGL_Line<T,2>(current_domain.min_corner+left_corner*current_dX,current_domain.min_corner+(left_corner+TV::Axis_Vector(axis))*current_dX);
                        OpenGL_Line<T,2>(current_domain.min_corner+(left_corner+TV::Axis_Vector(other_axis))*current_dX,current_domain.min_corner+(left_corner+TV::Axis_Vector(other_axis)+TV::Axis_Vector(axis))*current_dX);}}}
            OpenGL_End();
            glPopAttrib();}
    
        // Highlight current selection
        if(current_selection){
            if(current_selection->type == OPENGL_SELECTION::GRID_CELL_2D){
                OPENGL_SELECTION_SPGRID_CELL_2D<T>* real_selection=(OPENGL_SELECTION_SPGRID_CELL_2D<T>*)current_selection;
                TV min_corner=hierarchy->Grid(real_selection->level).Node(real_selection->index),max_corner=min_corner+hierarchy->Grid(real_selection->level).dX;
                OPENGL_SELECTION::Draw_Highlighted_Quad(min_corner,max_corner);}}}

    glPopAttrib();
    if(draw_boundary_conditions){
        Display_Psi_D_Boundary_Conditions();
        Display_Psi_N_Boundary_Conditions();}
    if(draw_nodal_flags)
        Display_Nodal_Flags();
    glPopMatrix();
}
//#####################################################################
// Function Display_Psi_D_Boundary_Conditions
//#####################################################################
template<class T> void OPENGL_SPGRID_2D<T>::
Display_Psi_D_Boundary_Conditions() const
{
    glPushAttrib(GL_LIGHTING_BIT | GL_TEXTURE_BIT | GL_POINT_BIT);
    glPointSize(5);glDisable(GL_LIGHTING);glDisable(GL_TEXTURE_2D);

    OPENGL_CONSTANT_COLOR_MAP<bool> color_map(OPENGL_COLOR::Magenta());
    OpenGL_Begin(GL_POINTS);
    for(int level=1;level<=hierarchy->Levels();++level) if(!draw_separate_levels || (draw_separate_levels && level==level_counter)){
        Const_flag_array_type flags=hierarchy->Allocator(level).Get_Const_Array(&T_STRUCT::flags);    
        for(SPGrid_Block_Iterator<typename Flag_array_type::MASK> iterator(hierarchy->Blocks(level));iterator.Valid();iterator.Next()){
            TV_INT index=iterator.Index().template Cast<TV_INT>();
            unsigned cell_flags=iterator.Data(flags);
            if((cell_flags&(SPGrid_Cell_Type_Interior) && !(cell_flags&(SPGrid_Cell_Type_Active))) || (cell_flags&(SPGrid_Cell_Type_Dirichlet))){
                color_map.Lookup(true).Send_To_GL_Pipeline();
                OpenGL_Vertex(hierarchy->Grid(level).Get_MAC_Grid().X(index));}}}
    OpenGL_End();
    glPopAttrib();
}
//#####################################################################
// Function Display_Psi_N_Boundary_Conditions
//#####################################################################
template<class T> void OPENGL_SPGRID_2D<T>::
Display_Psi_N_Boundary_Conditions() const
{
    glPushAttrib(GL_LIGHTING_BIT | GL_TEXTURE_BIT | GL_POINT_BIT);
    glPointSize(5);glDisable(GL_LIGHTING);glDisable(GL_TEXTURE_2D);

    OPENGL_CONSTANT_COLOR_MAP<bool> color_map(OPENGL_COLOR::Cyan());
    OpenGL_Begin(GL_POINTS);
    for(int level=hierarchy->Levels();level>=1;--level) if(!draw_separate_levels || (draw_separate_levels && level==level_counter)){
        const GRID<TV> current_grid=hierarchy->Grid(level).Get_MAC_Grid();
        Const_flag_array_type flags=hierarchy->Allocator(level).Get_Const_Array(&T_STRUCT::flags);    
        for(SPGrid_Block_Iterator<typename Data_array_type::MASK> iterator(hierarchy->Blocks(level));iterator.Valid();iterator.Next()){
            TV_INT index=iterator.Index().template Cast<TV_INT>();
            for(int axis=1;axis<=TV::dimension;++axis){
                unsigned face_valid_mask=GRID_TOPOLOGY_HELPER<typename Flag_array_type::MASK>::Face_Valid_Mask(axis);
                unsigned face_active_mask=GRID_TOPOLOGY_HELPER<typename Flag_array_type::MASK>::Face_Active_Mask(axis);
                if((iterator.Data(flags)&face_valid_mask) && !(iterator.Data(flags)&face_active_mask)){
                    color_map.Lookup(true).Send_To_GL_Pipeline();
                    OpenGL_Vertex(current_grid.X(index)-(T).5*current_grid.dX(axis)*TV::Axis_Vector(axis));}}}}

    OpenGL_End();
    glPopAttrib();
}
//#####################################################################
// Function Display_Nodal_Flags
//#####################################################################
template<class T> void OPENGL_SPGRID_2D<T>::
Display_Nodal_Flags() const
{
    glPushAttrib(GL_LIGHTING_BIT | GL_TEXTURE_BIT | GL_POINT_BIT);
    glPointSize(5);glDisable(GL_LIGHTING);glDisable(GL_TEXTURE_2D);

    unsigned node_mask;
    if(nodal_flag_counter==ACTIVE_NODES) node_mask=SPGrid_Node_Active;
    else if(nodal_flag_counter==COARSE_SHARED_NODES) node_mask=SPGrid_Node_Coarse_Shared;
    else if(nodal_flag_counter==T_JUNCTION_NODES) node_mask=SPGrid_Node_T_Junction;

    OPENGL_CONSTANT_COLOR_MAP<bool> color_map(OPENGL_COLOR::Yellow());
    OpenGL_Begin(GL_POINTS);
    for(int level=hierarchy->Levels();level>=1;--level) if(!draw_separate_levels || (draw_separate_levels && level==level_counter)){
        const GRID<TV> current_grid=hierarchy->Grid(level).Get_MAC_Grid();
        Const_flag_array_type flags=hierarchy->Allocator(level).Get_Const_Array(&T_STRUCT::flags);    
        for(SPGrid_Block_Iterator<typename Data_array_type::MASK> iterator(hierarchy->Blocks(level));iterator.Valid();iterator.Next())
            if(iterator.Data(flags)&node_mask){
                TV_INT index=iterator.Index().template Cast<TV_INT>();
                color_map.Lookup(true).Send_To_GL_Pipeline();
                OpenGL_Vertex(current_grid.Node(index));}}

    OpenGL_End();
    glPopAttrib();
}
//#####################################################################
// Function Initialize_Hierarchy
//#####################################################################
template<class T> void OPENGL_SPGRID_2D<T>::
Initialize_Hierarchy(const int frame)
{
    if(hierarchy) delete hierarchy;
    std::string directory_filename=FILE_UTILITIES::Get_Frame_Filename(filename,frame);
    int levels=0;
    FILE_UTILITIES::template Read_From_Text_File<int>(STRING_UTILITIES::string_sprintf("%s/levels",directory_filename.c_str()),levels);
    hierarchy=new GRID_HIERARCHY<T_STRUCT,T,2>(fine_mac_grid,levels);
    hierarchy->Initialize_Sets();
    hierarchy->Read_Block_Offsets(STRING_UTILITIES::string_sprintf("%s/block_offsets",directory_filename.c_str()));
    hierarchy->Read_Flags_Channel(STRING_UTILITIES::string_sprintf("%s/flags",directory_filename.c_str()));
}
//#####################################################################
// Function Set_Frame
//#####################################################################
template<class T> void OPENGL_SPGRID_2D<T>::
Set_Frame(int frame_input)
{
    frame=frame_input;
    return;
}
//#####################################################################
// Function Bounding_Box
//#####################################################################
template<class T> RANGE<VECTOR<float,3> > OPENGL_SPGRID_2D<T>::
Bounding_Box() const
{
    return World_Space_Box(RANGE<VECTOR<float,2> >(fine_mac_grid.domain));
}
//#####################################################################
// Function Get_Selection
//#####################################################################
template<class T> OPENGL_SELECTION *OPENGL_SPGRID_2D<T>::
Get_Selection(GLuint* buffer,int buffer_size)
{
    OPENGL_SELECTION* selection=0;
    if(buffer_size == 4){
        if(buffer[0] == 1) selection=new OPENGL_SELECTION_SPGRID_CELL_2D<T>(this,buffer[1],VECTOR<int,2>(buffer[2],buffer[3]));}
    return selection;
}
//#####################################################################
// Function Highlight_Selection
//#####################################################################
template<class T> void OPENGL_SPGRID_2D<T>::
Highlight_Selection(OPENGL_SELECTION *selection)
{
    delete current_selection;current_selection=0;
    if (selection->type == OPENGL_SELECTION::GRID_CELL_2D){
        OPENGL_SELECTION_SPGRID_CELL_2D<T>* real_selection=(OPENGL_SELECTION_SPGRID_CELL_2D<T>*)selection;
        current_selection=new OPENGL_SELECTION_SPGRID_CELL_2D<T>(this,real_selection->level,real_selection->index);}
}
//#####################################################################
// Function Clear_Highlight
//#####################################################################
template<class T> void OPENGL_SPGRID_2D<T>::
Clear_Highlight()
{
    delete current_selection;current_selection=0;
}
//#####################################################################
// Function Toggle_Draw_Separate_Levels
//#####################################################################
template<class T> void OPENGL_SPGRID_2D<T>::
Toggle_Draw_Separate_Levels()
{
    if(!draw_separate_levels){draw_separate_levels=true;level_counter=1;}
    else if(level_counter<hierarchy->Levels()) ++level_counter;
    else draw_separate_levels=false;
}
//#####################################################################
// Function Toggle_Draw_Psi_D_Boundary_Conditions
//#####################################################################
template<class T> void OPENGL_SPGRID_2D<T>::
Toggle_Draw_Boundary_Conditions()
{
    draw_boundary_conditions=!draw_boundary_conditions;
}
//#####################################################################
// Function Toggle_Draw_Nodal_Flags
//#####################################################################
template<class T> void OPENGL_SPGRID_2D<T>::
Toggle_Draw_Nodal_Flags()
{
    if(!draw_nodal_flags){draw_nodal_flags=true;nodal_flag_counter=1;}
    else if(nodal_flag_counter<NUMBER_OF_NODE_TYPES) ++nodal_flag_counter;
    else draw_nodal_flags=false;
}
//#####################################################################
// Print_Selection_Info
//#####################################################################
template<class T> void OPENGL_SPGRID_2D<T>::
Print_Selection_Info(std::ostream& stream,OPENGL_SELECTION* selection) const
{
    if(current_selection && current_selection->type==OPENGL_SELECTION::GRID_CELL_2D){
        int level=((OPENGL_SELECTION_SPGRID_CELL_2D<T>*)current_selection)->level;
        VECTOR<int,2> index=((OPENGL_SELECTION_SPGRID_CELL_2D<T>*)current_selection)->index;
        stream<<"Level: "<<level<<std::endl<<"Selected cell: "<<index<<" ("<<hierarchy->Grid(level).Get_MAC_Grid().X(index)<<")"<<std::endl;
        if(hierarchy->Array(level,&T_STRUCT::flags)(index[1],index[2])&SPGrid_Cell_Type_Ghost) stream<<"Ghost"<<std::endl;
        if(hierarchy->Array(level,&T_STRUCT::flags)(index[1],index[2])&SPGrid_Cell_Type_Dirichlet) stream<<"Dirichlet"<<std::endl;
        if(hierarchy->Array(level,&T_STRUCT::flags)(index[1],index[2])&SPGrid_Cell_Type_Active) stream<<"Active"<<std::endl;
        if(hierarchy->Array(level,&T_STRUCT::flags)(index[1],index[2])&SPGrid_Cell_Type_Interior) stream<<"Interior"<<std::endl;}
}
//#####################################################################
// Function Bounding_Box
//#####################################################################
template<class T> RANGE<VECTOR<float,3> > OPENGL_SELECTION_SPGRID_CELL_2D<T>::
Bounding_Box() const
{
    PHYSBAM_ASSERT(object);
    const GRID<TV>& grid=((OPENGL_SPGRID_2D<T>*)object)->hierarchy->Grid(level);
    RANGE<VECTOR<T,2> > box(grid.Node(index),grid.Node(index.x+1,index.y+1));
    return object->World_Space_Box(RANGE<VECTOR<float,2> >(box));
}
//#####################################################################
template class OPENGL_SPGRID_2D<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class OPENGL_SPGRID_2D<double>;
#endif
