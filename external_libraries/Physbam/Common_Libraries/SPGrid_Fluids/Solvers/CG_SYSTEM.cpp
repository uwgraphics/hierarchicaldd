//#####################################################################
// Copyright 2011, Sean Bauer, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class CG_SYSTEM
//#####################################################################
#include "CG_SYSTEM.h"
#include "CG_VECTOR.h"
#include <iomanip>

#include <SPGrid_Fluids/Solvers/Blocked_Set_Helper.h>
#include <SPGrid_Fluids/Solvers/Inner_Product_Helper.h>
#include <SPGrid_Fluids/Solvers/Convergence_Norm_Helper.h>

#include <SPGrid_Fluids/Simulation/FLUIDS_SIMULATION_DATA.h>
#include <SPGrid_Fluids/Projection/GRID_HIERARCHY_PROJECTION.h>
#include <SPGrid_Fluids/Solvers/HIERARCHY_PRECONDITIONER.h>
#include <SPGrid_Fluids/Solvers/Multigrid/MULTIGRID_SOLVER.h>
#include <SPGrid/Tools/SPGrid_Copy.h>
#include <SPGrid/Tools/SPGrid_Clear.h>

//#define TIMING

using namespace PhysBAM;

namespace PhysBAM{
extern int PhysBAM_number_of_threads;
}

//#####################################################################
// Constructor
//#####################################################################
template<class T_STRUCT,class T,int d> CG_SYSTEM<T_STRUCT,T,d>::
CG_SYSTEM(Hierarchy_type& hierarchy_input,ARRAY<Hierarchy_type*>* multigrid_hierarchy_in,T T_STRUCT::* mg_residual_channel_in,
    T T_STRUCT::* mg_diag_channel_in,T T_STRUCT::* mg_temp_channel_in,T T_STRUCT::* mg_x_channel_in,T T_STRUCT::* mg_b_channel_in,T T_STRUCT::* mg_cg_q_channel_in,
    T T_STRUCT::* mg_cg_r_channel_in,T T_STRUCT::* mg_cg_s_channel_in,T T_STRUCT::* mg_cg_k_channel_in,
    T T_STRUCT::* mg_cg_z_channel_in,const int interior_smoother_iterations_in,
    const int boundary_smoother_iterations_in,const T interior_omega_in,const T boundary_omega_in,
    const bool mg_diag_is_inverted_in,VECTOR<T T_STRUCT::*,d> L_channels_in,T T_STRUCT::* diag_channel_in,
    const ARRAY<int>& substitution_partitions_in)
    :BASE(true,false),hierarchy(hierarchy_input),A(T_MATRIX()),multigrid_hierarchy(multigrid_hierarchy_in),mg_residual_channel(mg_residual_channel_in),
    mg_diag_channel(mg_diag_channel_in),mg_temp_channel(mg_temp_channel_in),mg_x_channel(mg_x_channel_in),mg_b_channel(mg_b_channel_in),
    mg_cg_q_channel(mg_cg_q_channel_in),mg_cg_r_channel(mg_cg_r_channel_in),mg_cg_s_channel(mg_cg_s_channel_in),mg_cg_k_channel(mg_cg_k_channel_in),mg_cg_z_channel(mg_cg_z_channel_in),    
    interior_smoother_iterations(interior_smoother_iterations_in),boundary_smoother_iterations(boundary_smoother_iterations_in),
    interior_omega(interior_omega_in),boundary_omega(boundary_omega_in),mg_diag_is_inverted(mg_diag_is_inverted_in),L_channels(L_channels_in),diag_channel(diag_channel_in),
    ic_preconditioner(false),substitution_partitions(substitution_partitions_in),use_variable_beta(false)
{}
template<class T_STRUCT,class T,int d> CG_SYSTEM<T_STRUCT,T,d>::
CG_SYSTEM(Hierarchy_type& hierarchy_input,VECTOR<T T_STRUCT::*,d> L_channels_in,T T_STRUCT::* diag_channel_in,const ARRAY<int>& substitution_partitions_in)
    :BASE(true,false),hierarchy(hierarchy_input),multigrid_hierarchy(0),A(T_MATRIX()),L_channels(L_channels_in),diag_channel(diag_channel_in),
     interior_smoother_iterations(0),boundary_smoother_iterations(0),interior_omega((T)0.),boundary_omega((T)0.),mg_diag_is_inverted(false),
     ic_preconditioner(true),substitution_partitions(substitution_partitions_in),use_variable_beta(false)
{}
//#####################################################################
// Function Multiply
//#####################################################################
template<class T_STRUCT,class T,int d> void CG_SYSTEM<T_STRUCT,T,d>::
Multiply(const VECTOR_BASE& v,VECTOR_BASE& result) const
{
#ifdef TIMING
    LOG::SCOPE scope("CG_SYSTEM::Multiply");
#endif

    T T_STRUCT::* v_field = CG_VECTOR<T_STRUCT,T,d>::Cg_Vector(v).field;
    T T_STRUCT::* result_field = CG_VECTOR<T_STRUCT,T,d>::Cg_Vector(result).field;

    if(use_variable_beta)
        GRID_HIERARCHY_PROJECTION<T_STRUCT,T,d>::Compute_Variable_Beta_Laplacian(
            hierarchy,
            &T_STRUCT::flags,
            v_field,
            result_field,
            variable_beta_channel);
    else
        GRID_HIERARCHY_PROJECTION<T_STRUCT,T,d>::Compute_Laplacian(
            hierarchy,
            &T_STRUCT::flags,
            v_field,
            result_field);
}
//#####################################################################
// Function Inner_Product
//#####################################################################
template<class T_STRUCT,class T,int d> 
double CG_SYSTEM<T_STRUCT,T,d>::
Inner_Product(const VECTOR_BASE& v1,const VECTOR_BASE& v2) const
{
#ifdef TIMING
    LOG::SCOPE scope("CG_SYSTEM::Inner_Product");
#endif

    const Hierarchy_type& v1_hierarchy=CG_VECTOR<T_STRUCT,T,d>::Hierarchy(v1);
    const Hierarchy_type& v2_hierarchy=CG_VECTOR<T_STRUCT,T,d>::Hierarchy(v2);
    T T_STRUCT::* const v1_field = CG_VECTOR<T_STRUCT,T,d>::Cg_Vector(v1).field;
    T T_STRUCT::* const v2_field = CG_VECTOR<T_STRUCT,T,d>::Cg_Vector(v2).field;
    PHYSBAM_ASSERT(&hierarchy==&v1_hierarchy);
    PHYSBAM_ASSERT(&hierarchy==&v2_hierarchy);

    // Take dot-product of hierarchy, use doubles for temporaries
    double sum = 0;

    for(int level=1; level<=hierarchy.Levels();level++)
    {
        double level_sum = 0;
        Const_data_array_type d1 = v1_hierarchy.Allocator(level).Get_Array(v1_field);
        Const_data_array_type d2 = v2_hierarchy.Allocator(level).Get_Array(v2_field);
#ifdef SERIAL_INNER_PRODUCT
        for(SPGrid_Block_Iterator<typename Data_array_type::MASK> iterator(hierarchy.Blocks(level));iterator.Valid();iterator.Next())
            level_sum += iterator.Data(d1)*iterator.Data(d2);
#else
        Inner_Product_Helper<T,Data_array_type::MASK::elements_per_block> helper(
            (T*)d1.Get_Data_Ptr(),
            (T*)d2.Get_Data_Ptr(),
            (unsigned*)hierarchy.Array(level, &T_STRUCT::flags).Get_Data_Ptr(),
            hierarchy.Blocks(level).first,
            hierarchy.Blocks(level).second);
        if(PhysBAM_number_of_threads)
            level_sum = helper.Run_Parallel(PhysBAM_number_of_threads);
        else
            level_sum = helper.Run();
#endif
        sum+=level_sum;
    }

    return sum;
}
//#####################################################################
// Function Convergence_Norm
//#####################################################################
template<class T_STRUCT,class T,int d> T CG_SYSTEM<T_STRUCT,T,d>::
Convergence_Norm(const VECTOR_BASE& v) const
{
#ifdef TIMING
    LOG::SCOPE scope("CG_SYSTEM::Convergence_Norm");
#endif

    // Take maximum value of channel
    T T_STRUCT::* v_field = CG_VECTOR<T_STRUCT,T,d>::Cg_Vector(v).field;
    T maximum_value=0;
    T temp_value=0;

    for(int level=1; level<=hierarchy.Levels();level++)
    {
        Data_array_type d1 = hierarchy.Allocator(level).Get_Array(v_field);
#ifdef SERIAL_CONVERGENCE_NORM
        for(SPGrid_Block_Iterator<typename Data_array_type::MASK> iterator(hierarchy.Blocks(level));iterator.Valid();iterator.Next())
            if(iterator.Data(hierarchy.Set(level).array) & SPGrid_Cell_Type_Interior)
                maximum_value = PhysBAM::maxabs(maximum_value,iterator.Data(d1));
#else
        Convergence_Norm_Helper<T,Data_array_type::MASK::elements_per_block> helper(
            (T*)d1.Get_Data_Ptr(),
            hierarchy.Blocks(level).first,
            hierarchy.Blocks(level).second);
        if(PhysBAM_number_of_threads)
            maximum_value = PhysBAM::maxabs(maximum_value,helper.Run_Parallel(PhysBAM_number_of_threads));
        else
            maximum_value = PhysBAM::maxabs(maximum_value,helper.Run());
#endif
    }

    return maximum_value;
}
//#####################################################################
// Function Project
//#####################################################################
template<class T_STRUCT,class T,int d> void CG_SYSTEM<T_STRUCT,T,d>::
Project(VECTOR_BASE& v) const
{
    if(!ic_preconditioner){
        LOG::SCOPE scope("CG_SYSTEM::Apply_MG_Preconditioner");
    }
#ifdef TIMING
    LOG::SCOPE scope("CG_SYSTEM::Project");
#endif

    // Set all non-Interior nodes to zero.
    T T_STRUCT::* v_field = CG_VECTOR<T_STRUCT,T,d>::Cg_Vector(v).field;
    for(int level=1; level<=hierarchy.Levels();level++)
    {
        Data_array_type d1 = hierarchy.Allocator(level).Get_Array(v_field);
#ifdef SERIAL_PROJECT
        for(SPGrid_Block_Iterator<typename Data_array_type::MASK> iterator(hierarchy.Blocks(level));iterator.Valid();iterator.Next())
            if( !(iterator.Data(hierarchy.Set(level).array) & SPGrid_Cell_Type_Active) )
                iterator.Data(d1) = (T)0;
#else
        Blocked_Set_Helper<T, Data_array_type::MASK::elements_per_block> helper(
            (T*)d1.Get_Data_Ptr(),
            (T)0,
            (unsigned*)hierarchy.Array(level,&T_STRUCT::flags).Get_Data_Ptr(),
            hierarchy.Blocks(level).first,
            hierarchy.Blocks(level).second);
        if(PhysBAM_number_of_threads)
            helper.Run_Parallel(PhysBAM_number_of_threads);
        else
            helper.Run();
#endif
    }
}
//#####################################################################
// Function Set_Boundary_Conditions
//#####################################################################
template<class T_STRUCT,class T,int d> void CG_SYSTEM<T_STRUCT,T,d>::
Set_Boundary_Conditions(VECTOR_BASE& x) const
{
#ifdef TIMING
    LOG::SCOPE scope("CG_SYSTEM::Set_Boundary_Conditions");
#endif

#if 0
    Project(x);
#endif
}
//#####################################################################
// Function Project_Nullspace
//#####################################################################
template<class T_STRUCT,class T,int d> void CG_SYSTEM<T_STRUCT,T,d>::  
Project_Nullspace(VECTOR_BASE& x) const
{
#ifdef TIMING
    LOG::SCOPE scope("CG_SYSTEM::Project_Nullspace");
#endif

    return;

#ifdef DO_PROJECT_NULLSPACE
    Hierarchy_type& hierarchy=CG_VECTOR<T_STRUCT,T,d>::Hierarchy(x);
    T T_STRUCT::* const field = CG_VECTOR<T_STRUCT,T,d>::Cg_Vector(x).field;

    double sum = 0;
    double total_area = 0;

    for(int level=1; level<=hierarchy.Levels();level++)
    {
        // int cell_area = 1<<(level-1);
        int cell_area=1;
        cell_area *= cell_area;
        Const_data_array_type d1 = hierarchy.Allocator(level).Get_Const_Array(field);
#ifdef SERIAL_PROJECT_NULLSPACE
        for(SPGrid_Block_Iterator<typename Data_array_type::MASK> iterator(hierarchy.Blocks(level));iterator.Valid();iterator.Next())
            if(iterator.Data(hierarchy.Set(level).array) & SPGrid_Cell_Type_Interior)
            {
                sum += iterator.Data(d1)*cell_area;
                total_area += cell_area;
            }
#else
        Project_Nullspace_Helper<T,Data_array_type::MASK::elements_per_block> helper(
            (T*)d1.Get_Data_Ptr(),
            (unsigned*)hierarchy.Array(level, &T_STRUCT::flags).Get_Data_Ptr(),
            hierarchy.Blocks(level).first,
            hierarchy.Blocks(level).second);
        unsigned long count=0;
        sum += helper.Run_Parallel(num_threads,count)*cell_area;
        total_area += count*cell_area;        
#endif
    }

    double average_value = sum/total_area;

    for(int level=1; level<=hierarchy.Levels();level++)
    {
        Data_array_type d1 = hierarchy.Allocator(level).Get_Array(field);
#ifdef SERIAL_PROJECT_NULLSPACE
        for(SPGrid_Block_Iterator<typename Data_array_type::MASK> iterator(hierarchy.Blocks(level));iterator.Valid();iterator.Next())
            if(iterator.Data(hierarchy.Set(level).array) & SPGrid_Cell_Type_Interior)
                iterator.Data(d1) = iterator.Data(d1) - (T)average_value;
#else
        Blocked_SubtractConstant_Helper<T, Data_array_type::MASK::elements_per_block> helper(
            (T*)d1.Get_Data_Ptr(),
            (T)average_value,
            (unsigned*)hierarchy.Array(level, &T_STRUCT::flags).Get_Data_Ptr(),
            hierarchy.Blocks(level).first,
            hierarchy.Blocks(level).second);
        helper.Run_Parallel(num_threads);           
#endif
    }
#endif
}
//#####################################################################
// Function Use_Variable_Beta
//#####################################################################
template<class T_STRUCT,class T,int d> void CG_SYSTEM<T_STRUCT,T,d>::  
Use_Variable_Beta(T T_STRUCT::* variable_beta_channel_input)
{
    use_variable_beta=true;
    variable_beta_channel=variable_beta_channel_input;
}
//#####################################################################
// Function Use_Constant_Beta
//#####################################################################
template<class T_STRUCT,class T,int d> void CG_SYSTEM<T_STRUCT,T,d>::  
Use_Constant_Beta()
{
    use_variable_beta=false;
}
//#####################################################################
// Function Apply_Preconditioner
//#####################################################################
template<class T_STRUCT,class T,int d> void CG_SYSTEM<T_STRUCT,T,d>::  
Apply_Preconditioner(const VECTOR_BASE& r, VECTOR_BASE& z) const
{
    if(!ic_preconditioner){
        LOG::SCOPE scope("CG_SYSTEM::Apply_Preconditioner");
    }
#ifdef TIMING
    LOG::SCOPE scope("CG_SYSTEM::Apply_Preconditioner");
#endif

#ifdef IDENTITY_PRECONDITIONER
    for(int level=1;level<=hierarchy.Levels();level++) SPGrid_Computations::Copy<T_STRUCT,T,d>(hierarchy.Allocator(level),hierarchy.Blocks(level),CG_VECTOR<T_STRUCT,T,d>::Cg_Vector(r).field,CG_VECTOR<T_STRUCT,T,d>::Cg_Vector(z).field);
#else
    // Incomplete Cholesky
    if(ic_preconditioner){
        T T_STRUCT::* const r_channel=CG_VECTOR<T_STRUCT,T,d>::Cg_Vector(r).field;
        T T_STRUCT::* const z_channel=CG_VECTOR<T_STRUCT,T,d>::Cg_Vector(z).field;
        for(int level=1;level<=hierarchy.Levels();level++)
            if(PhysBAM_number_of_threads)SPGrid_Computations::Threading_Helper<T_STRUCT,d>(hierarchy.Allocator(level),hierarchy.Blocks(level)).Run_Parallel(SPGrid_Computations::Clear<T_STRUCT,T,d>(z_channel),PhysBAM_number_of_threads);
            else SPGrid_Computations::Clear<T_STRUCT,T,d>(hierarchy.Allocator(level),hierarchy.Blocks(level),z_channel);
        HIERARCHY_PRECONDITIONER<T_STRUCT,T,d>::Solve_Forward_Substitution (hierarchy,r_channel,z_channel,diag_channel,L_channels,substitution_partitions);
        for(int level=1;level<=hierarchy.Levels();level++) SPGrid_Computations::Masked_Clear(hierarchy.Allocator(level),hierarchy.Blocks(level),z_channel,&T_STRUCT::flags,(unsigned)(SPGrid_Cell_Type_Ghost|SPGrid_Cell_Type_Dirichlet));
        HIERARCHY_PRECONDITIONER<T_STRUCT,T,d>::Solve_Backward_Substitution(hierarchy,z_channel,z_channel,diag_channel,L_channels,substitution_partitions);}
    // MG
    else{T T_STRUCT::* const r_channel=CG_VECTOR<T_STRUCT,T,d>::Cg_Vector(r).field;
        T T_STRUCT::* const z_channel=CG_VECTOR<T_STRUCT,T,d>::Cg_Vector(z).field;
        for(int level=1;level<=hierarchy.Levels();level++) // clear z
            if(PhysBAM_number_of_threads) SPGrid_Computations::Threading_Helper<T_STRUCT,d>(hierarchy.Allocator(level),hierarchy.Blocks(level)).Run_Parallel(SPGrid_Computations::Clear<T_STRUCT,T,d>(z_channel),PhysBAM_number_of_threads);
            else SPGrid_Computations::Clear<T_STRUCT,T,d>(hierarchy.Allocator(level),hierarchy.Blocks(level),z_channel);
        for(int i=1;i<=(*multigrid_hierarchy).m;i++)
            for(int level=1;level<=(*multigrid_hierarchy)(i)->Levels();level++) // clear z (in mg)
                if(PhysBAM_number_of_threads) SPGrid_Computations::Threading_Helper<T_STRUCT,d>((*multigrid_hierarchy)(i)->Allocator(level),(*multigrid_hierarchy)(i)->Blocks(level)).Run_Parallel(SPGrid_Computations::Clear<T_STRUCT,T,d>(mg_x_channel),PhysBAM_number_of_threads);
                else SPGrid_Computations::Clear<T_STRUCT,T,d>((*multigrid_hierarchy)(i)->Allocator(level),(*multigrid_hierarchy)(i)->Blocks(level),mg_x_channel);
        for(int level=1;level<=hierarchy.Levels();level++){ // copy over b
            Const_flag_array_type flags=hierarchy.Allocator(level).Get_Const_Array(&T_STRUCT::flags);
            Const_data_array_type b=hierarchy.Allocator(level).Get_Const_Array(r_channel);
            Data_array_type mg_b=(*multigrid_hierarchy)(1)->Allocator(level).Get_Array(mg_b_channel);
            for(SPGrid_Block_Iterator<typename Data_array_type::MASK> iterator(hierarchy.Blocks(level));iterator.Valid();iterator.Next())
                if(iterator.Data(flags)&SPGrid_Cell_Type_Active) iterator.Data(mg_b)=iterator.Data(b);}
		MULTIGRID_SOLVER<T_STRUCT,T,d>::V_Cycle( // run V-Cycle
            (*multigrid_hierarchy),mg_x_channel,mg_b_channel,mg_residual_channel,
            mg_cg_q_channel,mg_cg_r_channel,mg_cg_s_channel,mg_cg_k_channel,mg_cg_z_channel,
            L_channels,diag_channel,substitution_partitions,
            mg_diag_channel,mg_temp_channel,interior_smoother_iterations,boundary_smoother_iterations,
            interior_omega,boundary_omega,mg_diag_is_inverted);
        for(int level=1;level<=hierarchy.Levels();level++){ // copy back x
            Const_flag_array_type flags=hierarchy.Allocator(level).Get_Const_Array(&T_STRUCT::flags);
            Data_array_type x=hierarchy.Allocator(level).Get_Array(z_channel);
            Const_data_array_type mg_x=(*multigrid_hierarchy)(1)->Allocator(level).Get_Const_Array(mg_x_channel);
            for(SPGrid_Block_Iterator<typename Data_array_type::MASK> iterator(hierarchy.Blocks(level));iterator.Valid();iterator.Next())
                if(iterator.Data(flags)&SPGrid_Cell_Type_Active) iterator.Data(x)=iterator.Data(mg_x);}}
#endif    
}
//#####################################################################
template class CG_SYSTEM<FLUIDS_SIMULATION_DATA<float>,float,2>;
template class CG_SYSTEM<FLUIDS_SIMULATION_DATA<float>,float,3>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class CG_SYSTEM<FLUIDS_SIMULATION_DATA<double>,double,2>;
template class CG_SYSTEM<FLUIDS_SIMULATION_DATA<double>,double,3>;
#endif
