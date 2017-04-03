//#####################################################################
// Copyright (c) 2014, Mridul Aanjaneya.
// Distributed under the FreeBSD license (see license.txt)
//#####################################################################
// Class OPENGL_SPGRID_2D
//##################################################################### 
#ifndef __OPENGL_SPGRID_2D__
#define __OPENGL_SPGRID_2D__
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_CALLBACK.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_COLOR.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_OBJECT.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_SELECTION.h>
#include <SPGrid_Fluids/Grids/GRID_HIERARCHY.h>
#include <SPGrid_Fluids/Simulation/FLUIDS_SIMULATION_DATA.h>

namespace PhysBAM{

template<class T>
class OPENGL_SPGRID_2D : public OPENGL_OBJECT
{
public:
    typedef VECTOR<T,2> TV;typedef VECTOR<int,2> TV_INT;
    typedef FLUIDS_SIMULATION_DATA<T> T_STRUCT;typedef GRID<TV> T_GRID;
    typedef typename SPGrid_Allocator<T_STRUCT,2>::template Array<unsigned>::type Flag_array_type;
    typedef typename SPGrid_Allocator<T_STRUCT,2>::template Array<const unsigned>::type Const_flag_array_type;
    typedef typename SPGrid_Allocator<T_STRUCT,2>::template Array<T>::type Data_array_type;
    typedef typename SPGrid_Allocator<T_STRUCT,2>::template Array<const T>::type Const_data_array_type;
    GRID<TV>& fine_mac_grid;
    GRID_HIERARCHY<T_STRUCT,T,2>* hierarchy;
    OPENGL_COLOR color;
    bool draw,draw_boundary_conditions;
    bool& draw_separate_levels;
    int& level_counter;
    bool& draw_nodal_flags;
    int& nodal_flag_counter;
    int current_frame;
private:
    OPENGL_SELECTION *current_selection;
    std::string filename;
    int frame;

    enum{ACTIVE_NODES=1,COARSE_SHARED_NODES,T_JUNCTION_NODES} NODE_TYPE;
    enum{NUMBER_OF_NODE_TYPES=3};

public:
    OPENGL_SPGRID_2D(GRID<TV> &fine_mac_grid_input,bool& draw_separate_levels_input,int& level_counter_input,bool& draw_nodal_flags_input,int& nodal_flag_counter_input,const std::string& filename_input,const OPENGL_COLOR &color_input=OPENGL_COLOR::White(),const int frame_input=0)
        :fine_mac_grid(fine_mac_grid_input),hierarchy(0),color(color_input),draw(true),draw_boundary_conditions(false),draw_separate_levels(draw_separate_levels_input),level_counter(level_counter_input),draw_nodal_flags(draw_nodal_flags_input),nodal_flag_counter(nodal_flag_counter_input),current_frame(-1),current_selection(0),filename(filename_input),frame(frame_input)
    {}

    void Display(const int in_color=1) const PHYSBAM_OVERRIDE;
    void Display_Psi_D_Boundary_Conditions() const;
    void Display_Psi_N_Boundary_Conditions() const;
    void Display_Nodal_Flags() const;
    virtual void Set_Frame(int frame_input);
    virtual RANGE<VECTOR<float,3> > Bounding_Box() const PHYSBAM_OVERRIDE;

    virtual OPENGL_SELECTION *Get_Selection(GLuint *buffer, int buffer_size);
    void Highlight_Selection(OPENGL_SELECTION *selection) PHYSBAM_OVERRIDE;
    void Clear_Highlight() PHYSBAM_OVERRIDE;
    void Initialize_Hierarchy(const int frame);

    void Toggle_Draw_Separate_Levels();
    void Toggle_Draw_Boundary_Conditions();
    void Toggle_Draw_Nodal_Flags();
    void Print_Selection_Info(std::ostream& stream,OPENGL_SELECTION* selection) const PHYSBAM_OVERRIDE;

    DEFINE_CALLBACK_CREATOR(OPENGL_SPGRID_2D, Toggle_Draw_Separate_Levels);
    DEFINE_CALLBACK_CREATOR(OPENGL_SPGRID_2D, Toggle_Draw_Boundary_Conditions);
    DEFINE_CALLBACK_CREATOR(OPENGL_SPGRID_2D, Toggle_Draw_Nodal_Flags);
};

template<class T>
class OPENGL_SELECTION_SPGRID_CELL_2D : public OPENGL_SELECTION
{
private:
    typedef VECTOR<T,2> TV;typedef VECTOR<int,2> TV_INT;
public:
    int level;
    VECTOR<int,2> index;
    OPENGL_SELECTION_SPGRID_CELL_2D(OPENGL_OBJECT *object,const int level=1,const VECTOR<int,2> &index=TV_INT()) 
        : OPENGL_SELECTION(OPENGL_SELECTION::GRID_CELL_2D, object),level(level),index(index) {}

    RANGE<VECTOR<float,3> > Bounding_Box() const PHYSBAM_OVERRIDE;
};
}
#endif
