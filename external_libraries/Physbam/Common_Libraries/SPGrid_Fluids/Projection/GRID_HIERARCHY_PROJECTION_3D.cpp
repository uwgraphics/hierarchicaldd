//#####################################################################
// Copyright 2012-2013, Sean Bauer, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class GRID_HIERARCHY_PROJECTION
//#####################################################################

//#include <math.h>

#include "GRID_HIERARCHY_PROJECTION.h"

using namespace PhysBAM;

namespace PhysBAM{
extern int PhysBAM_number_of_threads;
}

template<class T_STRUCT, class T,int d>
void GRID_HIERARCHY_PROJECTION<T_STRUCT,T,d>::Compute_Divergence(
    Hierarchy_type& hierarchy,
    unsigned T_STRUCT::* flags_field, 
    T T_STRUCT::* u_field, 
    T T_STRUCT::* v_field, 
    T T_STRUCT::* w_field, 
    T T_STRUCT::* d_field)
{

    for(int level=1;level<=hierarchy.Levels();level++)
    {
        Divergence_Helper<T,NextLogTwo<sizeof(T_STRUCT)>::value,d> helper(
            (T*)hierarchy.Array(level, d_field).Get_Data_Ptr(),
            (T*)hierarchy.Array(level, u_field).Get_Data_Ptr(),
            (T*)hierarchy.Array(level, v_field).Get_Data_Ptr(),
            (T*)hierarchy.Array(level, w_field).Get_Data_Ptr(),
            (unsigned*)hierarchy.Array(level, flags_field).Get_Data_Ptr(),
            hierarchy.Blocks(level).first,
            hierarchy.Blocks(level).second,
            Divergence_Scale_Uniform(hierarchy,level)
            );
        if(PhysBAM_number_of_threads)
            helper.Run_Parallel(PhysBAM_number_of_threads);
        else
            helper.Run();
    }
    
    Accumulate_Ghost_Values(hierarchy, flags_field, d_field);
}

template<class T_STRUCT, class T,int d>
void GRID_HIERARCHY_PROJECTION<T_STRUCT,T,d>::Compute_Gradient(
    Hierarchy_type& hierarchy,
    unsigned T_STRUCT::* flags_field, 
    T T_STRUCT::* gu_field, 
    T T_STRUCT::* gv_field, 
    T T_STRUCT::* gw_field, 
    T T_STRUCT::* d_field)
{
    Propagate_Ghost_Values(hierarchy, flags_field, d_field);

    for(int level=1;level<=hierarchy.Levels();level++)
    {
        Gradient_Helper<T,NextLogTwo<sizeof(T_STRUCT)>::value,d> helper(
            (T*)hierarchy.Array(level, gu_field).Get_Data_Ptr(),
            (T*)hierarchy.Array(level, gv_field).Get_Data_Ptr(),
            (T*)hierarchy.Array(level, gw_field).Get_Data_Ptr(),
            (T*)hierarchy.Array(level, d_field).Get_Data_Ptr(),
            (unsigned*)hierarchy.Array(level, flags_field).Get_Data_Ptr(),
            hierarchy.Blocks(level).first,
            hierarchy.Blocks(level).second,
            Gradient_Scale_Uniform(hierarchy,level),
            Gradient_Scale_Nonuniform(hierarchy,level)
            );
        if(PhysBAM_number_of_threads)
            helper.Run_Parallel(PhysBAM_number_of_threads);
        else
            helper.Run();
    }
}

template class GRID_HIERARCHY_PROJECTION<FLUIDS_SIMULATION_DATA<float>,float,3>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class GRID_HIERARCHY_PROJECTION<FLUIDS_SIMULATION_DATA<float>,float,3>;
#endif

