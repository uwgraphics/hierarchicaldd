//#####################################################################
// Copyright (c) 2014, Mridul Aanjaneya.
// Distributed under the FreeBSD license (see license.txt)
//#####################################################################
// Class OPENGL_SPGRID_MAC_VELOCITY_FIELD_2D
//##################################################################### 
#ifndef __OPENGL_SPGRID_MAC_VELOCITY_FIELD_2D__
#define __OPENGL_SPGRID_MAC_VELOCITY_FIELD_2D__

#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_COLOR.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_OBJECT.h>
#include <SPGrid_Fluids/Grids/GRID_HIERARCHY.h>
#include <SPGrid_Fluids/Simulation/FLUIDS_SIMULATION_DATA.h>
#include <SPGrid_Fluids/Visualization/OPENGL_SPGRID_2D.h>
namespace PhysBAM{
template<class T>
class OPENGL_SPGRID_MAC_VELOCITY_FIELD_2D:public OPENGL_OBJECT
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
    T T_STRUCT::* field1;
    T T_STRUCT::* field2;
    bool& draw_separate_levels;
    int& level_counter;
    OPENGL_COLOR vector_color;
    double size;
    bool draw_arrowhead,draw;
    T T_STRUCT::* divergence_channel;

    OPENGL_SPGRID_MAC_VELOCITY_FIELD_2D(GRID<TV> &fine_mac_grid_input,OPENGL_SPGRID_2D<T>& opengl_spgrid_input,T T_STRUCT::* field1_input,T T_STRUCT::* field2_input,bool& draw_separate_levels_input,int& level_counter_input,const OPENGL_COLOR &color=OPENGL_COLOR::White(),double size=0.025,bool draw_arrowhead=true);
    ~OPENGL_SPGRID_MAC_VELOCITY_FIELD_2D() {}

    void Display(const int in_color=1) const PHYSBAM_OVERRIDE;
    void Print_Selection_Info(std::ostream& stream,OPENGL_SELECTION* selection) const PHYSBAM_OVERRIDE;
    virtual RANGE<VECTOR<float,3> > Bounding_Box() const PHYSBAM_OVERRIDE;
    void Scale_Vector_Size(const T scale);
    void Draw_Vector(const TV& vector_field,const TV& vector_location,const double& current_size) const;
};
}
#endif
