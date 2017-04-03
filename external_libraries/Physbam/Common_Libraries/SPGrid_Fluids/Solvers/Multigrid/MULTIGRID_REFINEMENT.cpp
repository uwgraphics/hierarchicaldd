//#####################################################################
// Copyright (c) 2014, Raj Setaluri, Mridul Aanjneya
// Distributed under the FreeBSD license (see license.txt)
//#####################################################################
#include <SPGrid_Fluids/Solvers/Multigrid/MULTIGRID_REFINEMENT.h>
#include <SPGrid_Fluids/Simulation/FLUIDS_SIMULATION_DATA.h>
#include <SPGrid_Fluids/Simulation/FLUIDS_SIMULATION_FLAGS.h>
#include <SPGrid/Tools/SPGrid_Block_Iterator.h>
#include <SPGrid/Tools/SPGrid_Clear.h>
#include <SPGrid/Tools/SPGrid_Threading_Helper.h>
#include <SPGrid/Tools/SPGrid_Multiple_Allocator_Masked_Plus_Equals_Helper.h>
#include <SPGrid_Fluids/Grids/GRID_HIERARCHY_INITIALIZER.h>
#include <PhysBAM_Tools/Interpolation/LINEAR_INTERPOLATION.h>
#include <SPGrid_Fluids/Projection/Ghost_Value_Propagate.h>
#include <SPGrid_Fluids/Projection/Ghost_Value_Accumulate.h>
#include <SPGrid_Fluids/Solvers/Multigrid/Threaded_Restriction_Stencil_Helper.h>
#include <SPGrid_Fluids/Solvers/Multigrid/Threaded_Prolongation_Stencil_Helper.h>

using namespace PhysBAM;

namespace PhysBAM{
extern int PhysBAM_number_of_threads;
}

//#####################################################################
// Restrict
//#####################################################################
template<class T_STRUCT,class T,int d> void MULTIGRID_REFINEMENT<T_STRUCT,T,d>::
Restrict(T_HIERARCHY& fine_hierarchy,T_HIERARCHY& coarse_hierarchy,T T_STRUCT::* fine_data_channel,T T_STRUCT::* coarse_data_channel,VECTOR<int,2> finest_active_level)
{
    // grabbing levels
    const int fine_finest_active_level=finest_active_level(1);
    const int coarse_finest_active_level=finest_active_level(2);
    const int levels=fine_hierarchy.Levels();
    PHYSBAM_ASSERT(levels==coarse_hierarchy.Levels()); // debug
    PHYSBAM_ASSERT(fine_finest_active_level>=1&&coarse_finest_active_level<=levels); // debug

    // clear coarse data
    for(int level=1;level<=coarse_hierarchy.Levels();level++)
        if(PhysBAM_number_of_threads) SPGrid_Computations::Threading_Helper<T_STRUCT,d>(coarse_hierarchy.Allocator(level),coarse_hierarchy.Blocks(level)).Run_Parallel(
            SPGrid_Computations::Clear<T_STRUCT,T,d>(coarse_data_channel),PhysBAM_number_of_threads);
        else SPGrid_Computations::Clear<T_STRUCT,T,d>(coarse_hierarchy.Allocator(level),coarse_hierarchy.Blocks(level),coarse_data_channel);

    // restrict from fine to coarse for finest level
    {const int fine_level=fine_finest_active_level;const int coarse_level=coarse_finest_active_level;
    Data_array_type coarse_data=coarse_hierarchy.Allocator(coarse_level).Get_Array(coarse_data_channel);
    Const_data_array_type fine_data=fine_hierarchy.Allocator(fine_level).Get_Const_Array(fine_data_channel);
    Const_flag_array_type coarse_flags=coarse_hierarchy.Allocator(coarse_level).Get_Const_Array(&T_STRUCT::flags);
    Const_flag_array_type fine_flags=fine_hierarchy.Allocator(fine_level).Get_Const_Array(&T_STRUCT::flags);
    const T scale=(T)1./(T)restriction_stencil_denominator;
    Threaded_Restriction_Stencil_Helper<T,NextLogTwo<sizeof(T_STRUCT)>::value,d> threaded_restriction_stencil_helper(
        (T*)coarse_data.Get_Data_Ptr(),(T*)fine_data.Get_Data_Ptr(),(unsigned*)coarse_flags.Get_Data_Ptr(),
        coarse_hierarchy.Blocks(coarse_level).first,coarse_hierarchy.Blocks(coarse_level).second,
        scale,(unsigned)(SPGrid_Cell_Type_Active|SPGrid_Cell_Type_Ghost));
    if(PhysBAM_number_of_threads) threaded_restriction_stencil_helper.Run_Parallel(PhysBAM_number_of_threads);
    else threaded_restriction_stencil_helper.Run();}

    // accumulate within coarse
    for(int level=coarse_finest_active_level;level<levels;level++){
        Ghost_Value_Accumulate<T,T_STRUCT,d> helper(
            (unsigned*)coarse_hierarchy.Set(level+1).array.Get_Data_Ptr(),
            (T*)coarse_hierarchy.Array(level,coarse_data_channel).Get_Data_Ptr(),
            (T*)coarse_hierarchy.Array(level+1,coarse_data_channel).Get_Data_Ptr(),
            coarse_hierarchy.Blocks(level+1).first,
            coarse_hierarchy.Blocks(level+1).second);
        if(PhysBAM_number_of_threads) helper.Run_Parallel(PhysBAM_number_of_threads);
        else helper.Run();}

    // copy fine to coarse for all higher levels
    for(int level=levels;level>=coarse_finest_active_level;level--){
        SPGrid_Computations::Multiple_Allocator_Masked_Plus_Equals_Helper<T_STRUCT,T,unsigned,d> helper(
            coarse_hierarchy.Allocator(level),fine_hierarchy.Allocator(level),coarse_hierarchy.Allocator(level),fine_hierarchy.Allocator(level),
            fine_hierarchy.Blocks(level),coarse_data_channel,fine_data_channel,coarse_data_channel,&T_STRUCT::flags,(unsigned)(SPGrid_Cell_Type_Active));
        if(PhysBAM_number_of_threads) helper.Run_Parallel(PhysBAM_number_of_threads);
        else helper.Run();}
}
//#####################################################################
// Prolongate
//#####################################################################
template<class T_STRUCT,class T,int d> void MULTIGRID_REFINEMENT<T_STRUCT,T,d>::
Prolongate(T_HIERARCHY& fine_hierarchy,T_HIERARCHY& coarse_hierarchy,T T_STRUCT::* fine_data_channel,T T_STRUCT::* coarse_data_channel,VECTOR<int,2> finest_active_level)
{
    // auxiliary stuff
    unsigned long nodes_of_cell_offsets[nodes_per_cell];
    GRID_TOPOLOGY_HELPER<typename Flag_array_type::MASK>::Nodes_Of_Cell_Offsets(nodes_of_cell_offsets);
    VECTOR<unsigned long,d> parity_masks;
    for(int v=1;v<=d;v++) parity_masks(v)=GRID_TOPOLOGY_HELPER<typename Flag_array_type::MASK>::Axis_Vector_Offset(v);
    VECTOR<unsigned long,d> negative_axis_vector_offsets;
    for(int v=1;v<=d;v++) negative_axis_vector_offsets(v)=GRID_TOPOLOGY_HELPER<typename Flag_array_type::MASK>::Negative_Axis_Vector_Offset(v);

    // grabbing levels
    const int fine_finest_active_level=finest_active_level(1);
    const int coarse_finest_active_level=finest_active_level(2);
    const int levels=fine_hierarchy.Levels();
    PHYSBAM_ASSERT(levels==coarse_hierarchy.Levels()); // debug
    PHYSBAM_ASSERT(fine_finest_active_level>=1&&coarse_finest_active_level<=levels); // debug

    // clear fine data
    for(int level=1;level<=fine_hierarchy.Levels();level++)
        if(PhysBAM_number_of_threads) SPGrid_Computations::Threading_Helper<T_STRUCT,d>(fine_hierarchy.Allocator(level),fine_hierarchy.Blocks(level)).Run_Parallel(
            SPGrid_Computations::Clear<T_STRUCT,T,d>(fine_data_channel),PhysBAM_number_of_threads);
        else SPGrid_Computations::Clear<T_STRUCT,T,d>(fine_hierarchy.Allocator(level),fine_hierarchy.Blocks(level),fine_data_channel);

    // copy coarse to fine for all higher levels
    for(int level=levels;level>=coarse_finest_active_level;level--){
        SPGrid_Computations::Multiple_Allocator_Masked_Plus_Equals_Helper<T_STRUCT,T,unsigned,d> helper(
            fine_hierarchy.Allocator(level),coarse_hierarchy.Allocator(level),fine_hierarchy.Allocator(level),fine_hierarchy.Allocator(level),
            fine_hierarchy.Blocks(level),fine_data_channel,coarse_data_channel,fine_data_channel,&T_STRUCT::flags,(unsigned)(SPGrid_Cell_Type_Active));
        if(PhysBAM_number_of_threads) helper.Run_Parallel(PhysBAM_number_of_threads);
        else helper.Run();}

    // propagate within coarse
    for(int level=levels;level>coarse_finest_active_level;level--){
        Ghost_Value_Propagate<T,T_STRUCT,d> helper(
            (unsigned*)coarse_hierarchy.Set(level-1).array.Get_Data_Ptr(),
            (T*)coarse_hierarchy.Array(level-1,coarse_data_channel).Get_Data_Ptr(),
            (T*)coarse_hierarchy.Array(level,coarse_data_channel).Get_Data_Ptr(),
            coarse_hierarchy.Blocks(level-1).first,
            coarse_hierarchy.Blocks(level-1).second);
        if(PhysBAM_number_of_threads) helper.Run_Parallel(PhysBAM_number_of_threads);
        else helper.Run();}

    // prolongate from coarse to fine for finest level
    {const int fine_level=fine_finest_active_level;const int coarse_level=coarse_finest_active_level;
    Const_data_array_type coarse_data=coarse_hierarchy.Allocator(coarse_level).Get_Const_Array(coarse_data_channel);
    Data_array_type fine_data=fine_hierarchy.Allocator(fine_level).Get_Array(fine_data_channel);
    Const_flag_array_type fine_flags=fine_hierarchy.Allocator(fine_level).Get_Const_Array(&T_STRUCT::flags);
    Threaded_Prolongation_Stencil_Helper<T,NextLogTwo<sizeof(T_STRUCT)>::value,d> threaded_prolongation_stencil_helper(
        (T*)fine_data.Get_Data_Ptr(),(T*)coarse_data.Get_Data_Ptr(),(unsigned*)fine_flags.Get_Data_Ptr(),
        fine_hierarchy.Blocks(fine_level).first,fine_hierarchy.Blocks(fine_level).second,
        (unsigned)(SPGrid_Cell_Type_Active));
    if(PhysBAM_number_of_threads) threaded_prolongation_stencil_helper.Run_Parallel(PhysBAM_number_of_threads);
    else threaded_prolongation_stencil_helper.Run();}
}
//#####################################################################
template class MULTIGRID_REFINEMENT<FLUIDS_SIMULATION_DATA<float>,float,2>;
template class MULTIGRID_REFINEMENT<FLUIDS_SIMULATION_DATA<float>,float,3>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class MULTIGRID_REFINEMENT<FLUIDS_SIMULATION_DATA<double>,double,2>;
template class MULTIGRID_REFINEMENT<FLUIDS_SIMULATION_DATA<double>,double,3>;
#endif
