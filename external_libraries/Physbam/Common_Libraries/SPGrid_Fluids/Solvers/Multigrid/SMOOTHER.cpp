//#####################################################################
// Copyright (c) 2014, Raj Setaluri, Mridul Aanjneya
// Distributed under the FreeBSD license (see license.txt)
//#####################################################################
#include <SPGrid_Fluids/Solvers/Multigrid/SMOOTHER.h>
#include <SPGrid_Fluids/Solvers/Blocked_Set_Helper.h>
#include <SPGrid_Fluids/Simulation/FLUIDS_SIMULATION_DATA.h>
#include <SPGrid_Fluids/Simulation/FLUIDS_SIMULATION_FLAGS.h>
#include <SPGrid_Fluids/Projection/GRID_HIERARCHY_PROJECTION.h>
#include <SPGrid/Tools/SPGrid_Block_Iterator.h>
#include <SPGrid/Tools/SPGrid_Clear.h>
#include <SPGrid/Tools/SPGrid_Copy.h>
#include <SPGrid/Tools/SPGrid_Scale.h>
#include <SPGrid/Tools/SPGrid_Arithmetic.h>
#include <SPGrid/Tools/SPGrid_Threading_Helper.h>
#include <SPGrid_Fluids/Solvers/Multigrid/Smoother_Helper.h>

using namespace PhysBAM;

namespace PhysBAM{
extern int PhysBAM_number_of_threads;
}

//#####################################################################
// Multiply_With_System_Matrix
//#####################################################################
template<class T_STRUCT,class T,int d> void SMOOTHER<T_STRUCT,T,d>::
Multiply_With_System_Matrix(T_HIERARCHY& hierarchy,T T_STRUCT::* u_channel,T T_STRUCT::* Lu_channel)
{
    GRID_HIERARCHY_PROJECTION<T_STRUCT,T,d>::Compute_Laplacian(hierarchy,&T_STRUCT::flags,u_channel,Lu_channel);
    // clear Lu of non active values
    for(int level=1;level<=hierarchy.Levels();level++){
        Blocked_Set_Helper<T, Data_array_type::MASK::elements_per_block> helper(
            (T*)(hierarchy.Array(level,Lu_channel).Get_Data_Ptr()),
            (T)0,(unsigned*)hierarchy.Array(level,&T_STRUCT::flags).Get_Data_Ptr(),
            hierarchy.Blocks(level).first,hierarchy.Blocks(level).second);
        if(PhysBAM_number_of_threads) helper.Run_Parallel(PhysBAM_number_of_threads);
        else helper.Run();}
}
//#####################################################################
// Compute_Residual
//#####################################################################
template<class T_STRUCT,class T,int d> void SMOOTHER<T_STRUCT,T,d>::
Compute_Residual(T_HIERARCHY& hierarchy,T T_STRUCT::* u_channel,T T_STRUCT::* b_channel,T T_STRUCT::* r_channel)
{
    for(int level=1;level<=hierarchy.Levels();level++)
        if(PhysBAM_number_of_threads) SPGrid_Computations::Threading_Helper<T_STRUCT,d>(hierarchy.Allocator(level),hierarchy.Blocks(level)).Run_Parallel(
            SPGrid_Computations::Clear<T_STRUCT,T,d>(r_channel),PhysBAM_number_of_threads);
        else SPGrid_Computations::Clear<T_STRUCT,T,d>(hierarchy.Allocator(level),hierarchy.Blocks(level),r_channel);
    Multiply_With_System_Matrix(hierarchy,u_channel,r_channel);
    for(int level=1;level<=hierarchy.Levels();level++)
        if(PhysBAM_number_of_threads) SPGrid_Computations::Threading_Helper<T_STRUCT,d>(hierarchy.Allocator(level),hierarchy.Blocks(level)).Run_Parallel(
            SPGrid_Computations::Masked_Subtract<T_STRUCT,T,unsigned,d>(b_channel,r_channel,r_channel,&T_STRUCT::flags,(unsigned)SPGrid_Cell_Type_Active),PhysBAM_number_of_threads);
        else SPGrid_Computations::Masked_Subtract<T_STRUCT,T,unsigned,d>(hierarchy.Allocator(level),hierarchy.Blocks(level),b_channel,r_channel,r_channel,&T_STRUCT::flags,(unsigned)SPGrid_Cell_Type_Active);
}
//#####################################################################
// Single_Level_Jacobi_Smoother
//#####################################################################
template<class T_STRUCT,class T,int d> void SMOOTHER<T_STRUCT,T,d>::
Single_Level_Interior_Jacobi_Smoother(T_HIERARCHY& hierarchy,T T_STRUCT::* x_channel,T T_STRUCT::* b_channel,T T_STRUCT::* diagonal_channel,T T_STRUCT::* temp_channel,const int iterations,const T omega,const bool diagonal_is_inverted,const int* const smoothing_level)
{
    for(int level=1;level<=hierarchy.Levels();level++)
        SPGrid_Computations::Clear<T_STRUCT,T,d>(hierarchy.Allocator(level),hierarchy.Blocks(level),temp_channel);
    for(int i=1;i<=iterations;i++){
        Compute_Residual(hierarchy,x_channel,b_channel,temp_channel);
        if(diagonal_is_inverted){for(int level=1;level<=hierarchy.Levels();level++)
                SPGrid_Computations::Masked_Multiply<T_STRUCT,T,unsigned,d>(hierarchy.Allocator(level),hierarchy.Blocks(level),temp_channel,diagonal_channel,temp_channel,&T_STRUCT::flags,(unsigned)SPGrid_Cell_Type_Active);}
        else{for(int level=1;level<=hierarchy.Levels();level++)
                SPGrid_Computations::Masked_Divide<T_STRUCT,T,unsigned,d>(hierarchy.Allocator(level),hierarchy.Blocks(level),temp_channel,diagonal_channel,temp_channel,&T_STRUCT::flags,(unsigned)SPGrid_Cell_Type_Active);}
        for(int level=1;level<=hierarchy.Levels();level++)
            if(!smoothing_level || (*smoothing_level)==level)
                SPGrid_Computations::Masked_Saxpy<T_STRUCT,T,unsigned,d>(hierarchy.Allocator(level),hierarchy.Blocks(level),omega,temp_channel,x_channel,x_channel,&T_STRUCT::flags,(unsigned)SPGrid_Cell_Type_Active);}
}
//#####################################################################
// Single_Level_Boundary_Jacobi_Smoother
//#####################################################################
template<class T_STRUCT,class T,int d> void SMOOTHER<T_STRUCT,T,d>::
Single_Level_Boundary_Jacobi_Smoother(T_HIERARCHY& hierarchy,T T_STRUCT::* x_channel,T T_STRUCT::* b_channel,T T_STRUCT::* diagonal_channel,T T_STRUCT::* temp_channel,const int iterations,const T omega,const unsigned boundary_mask,const bool diagonal_is_inverted,const int* const smoothing_level)
{
    for(int level=1;level<=hierarchy.Levels();level++)
        SPGrid_Computations::Clear<T_STRUCT,T,d>(hierarchy.Allocator(level),hierarchy.Blocks(level),temp_channel);
    for(int i=1;i<=iterations;i++){
        Compute_Residual(hierarchy,x_channel,b_channel,temp_channel);
        if(diagonal_is_inverted){for(int level=1;level<=hierarchy.Levels();level++)
                SPGrid_Computations::Masked_Multiply<T_STRUCT,T,unsigned,d>(hierarchy.Allocator(level),hierarchy.Blocks(level),temp_channel,diagonal_channel,temp_channel,&T_STRUCT::flags,(unsigned)SPGrid_Cell_Type_Active);}
        else{for(int level=1;level<=hierarchy.Levels();level++)
                SPGrid_Computations::Masked_Divide<T_STRUCT,T,unsigned,d>(hierarchy.Allocator(level),hierarchy.Blocks(level),temp_channel,diagonal_channel,temp_channel,&T_STRUCT::flags,(unsigned)SPGrid_Cell_Type_Active);}
        for(int level=1;level<=hierarchy.Levels();level++)
            if(!smoothing_level || (*smoothing_level)==level)
                SPGrid_Computations::Masked_Saxpy<T_STRUCT,T,unsigned,d>(hierarchy.Allocator(level),hierarchy.Blocks(level),omega,temp_channel,x_channel,x_channel,&T_STRUCT::flags,(unsigned)(boundary_mask));}
}
//#####################################################################
// Single_Level_Boundary_Jacobi_Smoother
//#####################################################################
template<class T_STRUCT,class T,int d> void SMOOTHER<T_STRUCT,T,d>::
Single_Level_Jacobi_Smoother(T_HIERARCHY& hierarchy,T T_STRUCT::* x_channel,T T_STRUCT::* b_channel,T T_STRUCT::* delta_channel,T T_STRUCT::* dinv_channel,const int level,const int iterations,const T omega,const unsigned long mask,const std::pair<const unsigned long*,unsigned>& blocks)
{
    const T laplace_scale_uniform=GRID_HIERARCHY_PROJECTION<T_STRUCT,T,d>::Laplace_Scale_Uniform(hierarchy,level);
    const T laplace_scale_nonuniform=two_thirds*laplace_scale_uniform;
    // debug
    for(int l=level-1;l>=1;l--) PHYSBAM_ASSERT(hierarchy.Blocks(l).second==0);
    // debug
    unsigned long face_neighbor_offsets[faces_per_cell];
    GRID_TOPOLOGY_HELPER<typename Flag_array_type::MASK>::Face_Neighbor_Offsets(face_neighbor_offsets);
    if(PhysBAM_number_of_threads) SPGrid_Computations::Threading_Helper<T_STRUCT,d>(hierarchy.Allocator(level),blocks).Run_Parallel(SPGrid_Computations::Clear<T_STRUCT,T,d>(delta_channel),PhysBAM_number_of_threads);
    else SPGrid_Computations::Clear<T_STRUCT,T,d>(hierarchy.Allocator(level),blocks,delta_channel);
    Const_flag_array_type flags=hierarchy.Allocator(level).Get_Const_Array(&T_STRUCT::flags);
    Const_data_array_type b=hierarchy.Allocator(level).Get_Const_Array(b_channel);
    Const_data_array_type dinv=hierarchy.Allocator(level).Get_Const_Array(dinv_channel);
    Data_array_type x=hierarchy.Allocator(level).Get_Array(x_channel);
    Data_array_type delta=hierarchy.Allocator(level).Get_Array(delta_channel);
    for(int iteration=1;iteration<=iterations;iteration++){
        Smoother_Helper<T,NextLogTwo<sizeof(T_STRUCT)>::value,d> smoother_helper((T*)x.Get_Data_Ptr(),(T*)b.Get_Data_Ptr(),(T*)delta.Get_Data_Ptr(),(T*)dinv.Get_Data_Ptr(),(unsigned*)flags.Get_Data_Ptr(),blocks.first,blocks.second,laplace_scale_uniform,laplace_scale_nonuniform,omega,mask);
        if(PhysBAM_number_of_threads) smoother_helper.Run_Parallel(PhysBAM_number_of_threads);
        else smoother_helper.Run();
        if(PhysBAM_number_of_threads) SPGrid_Computations::Threading_Helper<T_STRUCT,d>(hierarchy.Allocator(level),blocks).Run_Parallel(
            SPGrid_Computations::Masked_Add<T_STRUCT,T,unsigned,d>(x_channel,delta_channel,x_channel,&T_STRUCT::flags,mask),PhysBAM_number_of_threads);
        else SPGrid_Computations::Masked_Add<T_STRUCT,T,unsigned,d>(hierarchy.Allocator(level),blocks,x_channel,delta_channel,x_channel,&T_STRUCT::flags,mask);}
}
//#####################################################################
template class SMOOTHER<FLUIDS_SIMULATION_DATA<float>,float,2>;
template class SMOOTHER<FLUIDS_SIMULATION_DATA<float>,float,3>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class SMOOTHER<FLUIDS_SIMULATION_DATA<double>,double,2>;
template class SMOOTHER<FLUIDS_SIMULATION_DATA<double>,double,3>;
#endif
