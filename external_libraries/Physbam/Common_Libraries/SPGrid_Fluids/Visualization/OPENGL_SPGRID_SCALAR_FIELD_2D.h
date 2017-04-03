//#####################################################################
// Copyright (c) 2014, Mridul Aanjaneya.
// Distributed under the FreeBSD license (see license.txt)
//#####################################################################
// Class OPENGL_SPGRID_SCALAR_FIELD_2D
//##################################################################### 
#ifndef __OPENGL_SPGRID_SCALAR_FIELD_2D__
#define __OPENGL_SPGRID_SCALAR_FIELD_2D__

#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_CALLBACK.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_COLOR_MAP.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_OBJECT.h>
#include <SPGrid_Fluids/Grids/GRID_HIERARCHY.h>
#include <SPGrid_Fluids/Simulation/FLUIDS_SIMULATION_DATA.h>
#include <SPGrid_Fluids/Visualization/OPENGL_SPGRID_2D.h>

namespace PhysBAM{
class OPENGL_TEXTURED_RECT;
template<class T,class T2=T>
class OPENGL_SPGRID_SCALAR_FIELD_2D : public OPENGL_OBJECT
{
public:
    typedef VECTOR<T,2> TV;typedef VECTOR<int,2> TV_INT;
    typedef FLUIDS_SIMULATION_DATA<T> T_STRUCT;typedef GRID<TV> T_GRID;
    typedef typename SPGrid_Allocator<T_STRUCT,2>::template Array<unsigned>::type Flag_array_type;
    typedef typename SPGrid_Allocator<T_STRUCT,2>::template Array<const unsigned>::type Const_flag_array_type;
    typedef typename SPGrid_Allocator<T_STRUCT,2>::template Array<T>::type Data_array_type;
    typedef typename SPGrid_Allocator<T_STRUCT,2>::template Array<const T>::type Const_data_array_type;
    GRID<TV>& fine_mac_grid;
    OPENGL_SPGRID_2D<T>& opengl_spgrid;

    enum DRAW_MODE {DRAW_TEXTURE, DRAW_POINTS};
    ARRAY<OPENGL_COLOR_MAP<T2>*> color_maps; // all owned by us
    bool& draw_separate_levels;
    int& level_counter;
    T T_STRUCT::* field;
    DRAW_MODE draw_mode;
    int current_color_map;
private:
    mutable OPENGL_TEXTURED_RECT *opengl_textured_rect;
    bool scale_range;
    T2 scale_range_min,scale_range_dx;

public:
    OPENGL_SPGRID_SCALAR_FIELD_2D(GRID<TV> &fine_mac_grid_input,OPENGL_SPGRID_2D<T>& opengl_spgrid_input,OPENGL_COLOR_MAP<T2>* color_map_input,bool& draw_separate_levels_input,int& level_counter_input,T T_STRUCT::* field_input,DRAW_MODE draw_mode_input=DRAW_TEXTURE);
    ~OPENGL_SPGRID_SCALAR_FIELD_2D();

    void Set_Scale_Range(const T2 range_min,const T2 range_max);
    void Reset_Scale_Range();
    T2 Pre_Map_Value(const T2 value) const;

    void Display(const int in_color=1) const PHYSBAM_OVERRIDE;
    virtual RANGE<VECTOR<float,3> > Bounding_Box() const PHYSBAM_OVERRIDE;

    void Set_Draw_Mode(DRAW_MODE draw_mode);
    virtual void Update() const;  // Call when attributes have changed
    void Display_Points(const int in_color) const;

    void Print_Selection_Info(std::ostream& output_stream,OPENGL_SELECTION* current_selection) const PHYSBAM_OVERRIDE;
    void Toggle_Color_Map();
    void Toggle_Draw_Mode();

private:
    void Initialize_Color_Maps(OPENGL_COLOR_MAP<T2>* color_map_input);
    void Update_Texture(const VECTOR<int,2>& start_index,const VECTOR<int,2>& end_index) const;
};
}
#endif
