//#####################################################################
// Copyright 2012-2013, Sean Bauer, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class GRID_HIERARCHY_AVERAGING
//#####################################################################
#ifndef __GRID_HIERARCHY_AVERAGING_h__
#define __GRID_HIERARCHY_AVERAGING_h__

#include <SPGrid_Fluids/Grids/GRID_HIERARCHY.h>

namespace PhysBAM{

template<class T_STRUCT,class T,int d>
class GRID_HIERARCHY_AVERAGING
{
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<unsigned>::type Flag_array_type;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<T>::type Data_array_type;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<const unsigned>::type Const_flag_array_type;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<const T>::type Const_data_array_type;
    typedef typename Flag_array_type::MASK Flag_array_mask;
    typedef GRID_HIERARCHY<T_STRUCT,T,d> T_HIERARCHY;

//#####################################################################
public:
    static void Average_Face_Velocities_To_Nodes(T_HIERARCHY& hierarchy,const VECTOR<T T_STRUCT::*,d> face_velocities,
        const VECTOR<T T_STRUCT::*,d> node_velocities,unsigned T_STRUCT::* const flags_field,T T_STRUCT::* const weight_field);
    static void Average_Cell_Density_To_Nodes(T_HIERARCHY& hierarchy,T T_STRUCT::* const cell_density,
        T T_STRUCT::* const node_density,unsigned T_STRUCT::* const flags_field,T T_STRUCT::* const weight_field);
    static void Average_Cell_Density_To_Vertical_Faces(T_HIERARCHY& hierarchy,T T_STRUCT::* const cell_density,
        T T_STRUCT::* const vertical_face_density,unsigned T_STRUCT::* const flags_field,const int vertical_axis);
    static void Average_Node_Density_To_Vertical_Faces(T_HIERARCHY& hierarchy,T T_STRUCT::* const node_density,
        T T_STRUCT::* const vertical_face_density,unsigned T_STRUCT::* const flags_field,const int vertical_axis,unsigned mask);
//#####################################################################
};
}
#endif
