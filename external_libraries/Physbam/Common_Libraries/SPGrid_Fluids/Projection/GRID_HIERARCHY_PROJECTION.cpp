//#####################################################################
// Copyright 2012-2013, Sean Bauer, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class GRID_HIERARCHY_PROJECTION
//#####################################################################
#include <SPGrid_Fluids/Projection/GRID_HIERARCHY_PROJECTION.h>

#include <SPGrid_Fluids/Simulation/FLUIDS_SIMULATION_DATA.h>
#include <SPGrid_Fluids/Projection/Ghost_Value_Propagate.h>
#include <SPGrid_Fluids/Projection/Ghost_Value_Accumulate.h>
#include <SPGrid_Fluids/Projection/SPGRID_LAPLACE.h>
#include <SPGrid_Fluids/Projection/Laplace_Helper.h>
#include <SPGrid_Fluids/Projection/Variable_Beta_Laplace_Helper.h>

//#define SERIAL_LAPLACE

using namespace PhysBAM;

namespace PhysBAM{
extern int PhysBAM_number_of_threads;
}

template<class T_STRUCT, class T,int d>
void GRID_HIERARCHY_PROJECTION<T_STRUCT,T,d>::Propagate_Ghost_Values(
    Hierarchy_type& hierarchy,
    unsigned T_STRUCT::* flags_field, 
    T T_STRUCT::* u_field)
{
    for (int level=hierarchy.Levels()-1;level>=1;level--)
    {
        Ghost_Value_Propagate<T,T_STRUCT,d>
            helper( (unsigned*)hierarchy.Array(level,flags_field).Get_Data_Ptr(),
                    (T*)hierarchy.Array(level, u_field).Get_Data_Ptr(),
                    (T*)hierarchy.Array(level+1, u_field).Get_Data_Ptr(),
                    hierarchy.Blocks(level).first,
                    hierarchy.Blocks(level).second);
        if(PhysBAM_number_of_threads)
            helper.Run_Parallel(PhysBAM_number_of_threads);
        else
            helper.Run();
    }
}

template<class T_STRUCT, class T,int d>
void GRID_HIERARCHY_PROJECTION<T_STRUCT,T,d>::Accumulate_Ghost_Values(
    Hierarchy_type& hierarchy,
    unsigned T_STRUCT::* flags_field, 
    T T_STRUCT::* u_field)
{
    for (int level=1;level<=hierarchy.Levels()-1;level++)
    {
        Ghost_Value_Accumulate<T,T_STRUCT,d>
            helper( (unsigned*)hierarchy.Array(level+1,flags_field).Get_Data_Ptr(),
                    (T*)hierarchy.Array(level, u_field).Get_Data_Ptr(),
                    (T*)hierarchy.Array(level+1, u_field).Get_Data_Ptr(),
                    hierarchy.Blocks(level+1).first,
                    hierarchy.Blocks(level+1).second);
        if(PhysBAM_number_of_threads)
            helper.Run_Parallel(PhysBAM_number_of_threads);
        else
            helper.Run();
    }
}

template<class T_STRUCT, class T,int d>
void GRID_HIERARCHY_PROJECTION<T_STRUCT,T,d>::Compute_Laplacian(
    Hierarchy_type& hierarchy,
    unsigned T_STRUCT::* flags_field, 
    T T_STRUCT::* u_field,
    T T_STRUCT::* result_field)
{
    Propagate_Ghost_Values(hierarchy, flags_field, u_field);

    for(int level=1;level<=hierarchy.Levels();level++){
#ifdef SERIAL_LAPLACE
        SPGRID_LAPLACE<Data_array_type,Flag_array_type>::Compute(
            hierarchy.Blocks(level),
            hierarchy.Array(level, u_field),
            hierarchy.Array(level, result_field),
            hierarchy.Array(level, &T_STRUCT::flags),
            Laplace_Scale_Uniform(hierarchy,level),
            Laplace_Scale_Nonuniform(hierarchy,level)
            );
#else
        Laplace_Helper<T,NextLogTwo<sizeof(T_STRUCT)>::value,d> helper(
            (T*)hierarchy.Array(level, result_field).Get_Data_Ptr(),
            (T*)hierarchy.Array(level, u_field).Get_Data_Ptr(),
            (unsigned*)hierarchy.Array(level, &T_STRUCT::flags).Get_Data_Ptr(),
            hierarchy.Blocks(level).first,
            hierarchy.Blocks(level).second,
            Laplace_Scale_Uniform(hierarchy,level),
            Laplace_Scale_Nonuniform(hierarchy,level)
            );
        if(PhysBAM_number_of_threads)
            helper.Run_Parallel(PhysBAM_number_of_threads);
        else
            helper.Run();
#endif        
    }

    Accumulate_Ghost_Values(hierarchy, flags_field, result_field);
}

template<class T_STRUCT, class T,int d>
void GRID_HIERARCHY_PROJECTION<T_STRUCT,T,d>::Compute_Variable_Beta_Laplacian( // variable beta
    Hierarchy_type& hierarchy,
    unsigned T_STRUCT::* flags_field, 
    T T_STRUCT::* u_field,
    T T_STRUCT::* result_field,
    T T_STRUCT::* variable_beta_field)
{
    Propagate_Ghost_Values(hierarchy, flags_field, u_field);

    for(int level=1;level<=hierarchy.Levels();level++){
#ifdef SERIAL_LAPLACE
        PHYSBAM_NOT_IMPLEMENTED();
#else
        Variable_Beta_Laplace_Helper<T,NextLogTwo<sizeof(T_STRUCT)>::value,d> helper(
            (T*)hierarchy.Array(level, result_field).Get_Data_Ptr(),
            (T*)hierarchy.Array(level, u_field).Get_Data_Ptr(),
            (T*)hierarchy.Array(level, variable_beta_field).Get_Data_Ptr(),
            (unsigned*)hierarchy.Array(level, &T_STRUCT::flags).Get_Data_Ptr(),
            hierarchy.Blocks(level).first,
            hierarchy.Blocks(level).second,
            Laplace_Scale_Uniform(hierarchy,level),
            Laplace_Scale_Nonuniform(hierarchy,level)
            );
        if(PhysBAM_number_of_threads)
            helper.Run_Parallel(PhysBAM_number_of_threads);
        else
            helper.Run();
#endif        
    }

    Accumulate_Ghost_Values(hierarchy, flags_field, result_field);
}

template class GRID_HIERARCHY_PROJECTION<FLUIDS_SIMULATION_DATA<float>,float,2>;
template class GRID_HIERARCHY_PROJECTION<FLUIDS_SIMULATION_DATA<float>,float,3>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class GRID_HIERARCHY_PROJECTION<FLUIDS_SIMULATION_DATA<double>,double,2>;
template class GRID_HIERARCHY_PROJECTION<FLUIDS_SIMULATION_DATA<double>,double,3>;
#endif
