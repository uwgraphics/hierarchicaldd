//#####################################################################
// Copyright (c) 2014, Raj Setaluri, Mridul Aanjneya
// Distributed under the FreeBSD license (see license.txt)
//#####################################################################
#include <SPGrid_Fluids/Solvers/Multigrid/MULTIGRID_SOLVER.h>
#include <SPGrid_Fluids/Simulation/FLUIDS_SIMULATION_DATA.h>
#include <SPGrid_Fluids/Simulation/FLUIDS_SIMULATION_FLAGS.h>
#include <SPGrid_Fluids/Solvers/CG_VECTOR.h>
#include <SPGrid_Fluids/Solvers/CG_SYSTEM.h>
#include <PhysBAM_Tools/Krylov_Solvers/CONJUGATE_GRADIENT.h>
#include <SPGrid_Fluids/Solvers/Multigrid/SMOOTHER.h>
#include <SPGrid_Fluids/Solvers/Multigrid/MULTIGRID_REFINEMENT.h>
#include <SPGrid/Tools/SPGrid_Arithmetic.h>
#include <SPGrid/Tools/SPGrid_Clear.h>
#include <SPGrid/Tools/SPGrid_Copy.h>
#include <SPGrid/Tools/SPGrid_Threading_Helper.h>
#include <SPGrid_Fluids/Visualization/VISUALIZE_HEIGHTFIELD.h>
#include <PhysBAM_Tools/Random_Numbers/RANDOM_NUMBERS.h>
#include <SPGrid_Fluids/Solvers/BUILD_MATRIX.h>
#include <SPGrid_Fluids/Projection/GRID_HIERARCHY_PROJECTION.h>
#include <PhysBAM_Tools/Krylov_Solvers/PCG_SPARSE.h>

using namespace PhysBAM;

namespace PhysBAM{
extern int PhysBAM_number_of_threads;
}

//#####################################################################
// V_Cycle
//#####################################################################
template<class T_STRUCT,class T,int d> void MULTIGRID_SOLVER<T_STRUCT,T,d>::
V_Cycle(ARRAY<T_HIERARCHY*>& multigrid_hierarchy,T T_STRUCT::* x_channel,T T_STRUCT::* b_channel,T T_STRUCT::* residual_channel,
    T T_STRUCT::* cg_q_channel,T T_STRUCT::* cg_r_channel,T T_STRUCT::* cg_s_channel,T T_STRUCT::* cg_k_channel,T T_STRUCT::* cg_z_channel,
    VECTOR<T T_STRUCT::*,d> L_channels,T T_STRUCT::* ic_diag_channel,const ARRAY<int>& substitution_partitions,
    T T_STRUCT::* mg_diag_channel,T T_STRUCT::* temp_channel,const int interior_smoother_iterations,const int boundary_smoother_iterations,
    const T interior_omega,const T boundary_omega,const bool diagonal_is_inverted)
{
    const int multigrid_levels=multigrid_hierarchy.m;
    LOG::SCOPE scope("MULTIGRID_SOLVER::Preconditioner");
    PHYSBAM_ASSERT(diagonal_is_inverted);

    // debug
    for(int i=1;i<=multigrid_levels;i++){
        int finest_active_level=1;
        for(;finest_active_level<=multigrid_hierarchy(i)->Levels();finest_active_level++)
            if(multigrid_hierarchy(i)->Blocks(finest_active_level).second>0) break;
        PHYSBAM_ASSERT(i==finest_active_level);}
    // debug

    // downstroke
    for(int i=1;i<multigrid_levels;i++){
        // clear initial guess into smoother
        SMOOTHER<T_STRUCT,T,d>::Single_Level_Jacobi_Smoother(*multigrid_hierarchy(i),x_channel,b_channel,temp_channel,mg_diag_channel,i,boundary_smoother_iterations,boundary_omega,(unsigned)(MG_Boundary),multigrid_hierarchy(i)->Boundary_Blocks(i));
        SMOOTHER<T_STRUCT,T,d>::Single_Level_Jacobi_Smoother(*multigrid_hierarchy(i),x_channel,b_channel,temp_channel,mg_diag_channel,i,interior_smoother_iterations,interior_omega,(unsigned)(SPGrid_Cell_Type_Active),multigrid_hierarchy(i)->Blocks(i));
        SMOOTHER<T_STRUCT,T,d>::Single_Level_Jacobi_Smoother(*multigrid_hierarchy(i),x_channel,b_channel,temp_channel,mg_diag_channel,i,boundary_smoother_iterations,boundary_omega,(unsigned)(MG_Boundary),multigrid_hierarchy(i)->Boundary_Blocks(i));
        SMOOTHER<T_STRUCT,T,d>::Compute_Residual(*multigrid_hierarchy(i),x_channel,b_channel,residual_channel);        
        MULTIGRID_REFINEMENT<T_STRUCT,T,d>::Restrict(*multigrid_hierarchy(i),*multigrid_hierarchy(i+1),residual_channel,b_channel,VECTOR<int,2>(i,i+1));
        for(int level=1;level<=multigrid_hierarchy(i+1)->Levels();level++)
            if(PhysBAM_number_of_threads) SPGrid_Computations::Threading_Helper<T_STRUCT,d>(multigrid_hierarchy(i+1)->Allocator(level),multigrid_hierarchy(i+1)->Blocks(level)).Run_Parallel(
                SPGrid_Computations::Clear<T_STRUCT,T,d>(x_channel),PhysBAM_number_of_threads);
            else SPGrid_Computations::Clear<T_STRUCT,T,d>(multigrid_hierarchy(i+1)->Allocator(level),multigrid_hierarchy(i+1)->Blocks(level),x_channel);}
    //#define CG_AT_BOTTOM
#ifdef CG_AT_BOTTOM
    // solve at bottom
//#define PCG_SPARSE_SOLVE_AT_BOTTOM
#ifdef PCG_SPARSE_SOLVE_AT_BOTTOM
    {LOG::SCOPE scope("MULTIGRID_SOLVER::BOTTOM_SOLVE");
    T_HIERARCHY* const bottom_hierarchy=multigrid_hierarchy(multigrid_levels);
    typedef PAIR<int,unsigned long> CID;
    VECTOR<T T_STRUCT::*,d> U_channels; // U is temporary -- watch out for channels!
    U_channels(1)=x_channel;
    U_channels(2)=cg_r_channel;
    if(d==3) U_channels(3)=cg_z_channel;
    HASHTABLE<CID,int> active_dof_hash;ARRAY<CID> active_dof_array;
    for(int level=1;level<=bottom_hierarchy->Levels();level++){ // clear L, diag at this level
        for(int v=1;v<=d;v++) SPGrid_Computations::Clear<T_STRUCT,T,d,2>(bottom_hierarchy->Allocator(level),bottom_hierarchy->Blocks(level),L_channels(v),U_channels(v));
        SPGrid_Computations::Clear<T_STRUCT,T,d>(bottom_hierarchy->Allocator(level),bottom_hierarchy->Blocks(level),ic_diag_channel);}
    BUILD_MATRIX<T_STRUCT,T,d>::Generate_Conversion_Structures(active_dof_hash,active_dof_array,*bottom_hierarchy,&T_STRUCT::flags);
    BUILD_MATRIX<T_STRUCT,T,d>::Build_Matrix(*bottom_hierarchy,ic_diag_channel,L_channels,U_channels,&T_STRUCT::flags);
    SPARSE_MATRIX_FLAT_NXN<T> matrix;
    VECTOR_ND<T> vecX;
    VECTOR_ND<T> vecB;
    VECTOR_ND<T> vecQ;
    VECTOR_ND<T> vecS;
    VECTOR_ND<T> vecR;
    VECTOR_ND<T> vecK;
    VECTOR_ND<T> vecZ;
    vecX.Resize(active_dof_hash.Size());
    vecB.Resize(active_dof_hash.Size());
    vecQ.Resize(active_dof_hash.Size());
    vecS.Resize(active_dof_hash.Size());
    vecR.Resize(active_dof_hash.Size());
    vecK.Resize(active_dof_hash.Size());
    vecZ.Resize(active_dof_hash.Size());
    BUILD_MATRIX<T_STRUCT,T,d>::SPGrid_Channels_To_Matrix(*bottom_hierarchy,ic_diag_channel,L_channels,U_channels,matrix,active_dof_hash,&T_STRUCT::flags);
    BUILD_MATRIX<T_STRUCT,T,d>::SPGrid_Channel_To_Vector(active_dof_array,*bottom_hierarchy,b_channel,vecB);
    PCG_SPARSE<T> pcg;
    const int pcg_iterations=1000;const int pcg_restart_iterations=pcg_iterations+1;
    //pcg.Show_Residuals(true);pcg.Show_Results(true);
    pcg.Set_Maximum_Iterations(pcg_iterations);
    pcg.cg_restart_iterations=pcg_restart_iterations;
    pcg.Use_Incomplete_Cholesky();
    pcg.Solve(matrix,vecX,vecB,vecQ,vecS,vecR,vecK,vecZ,(T)1e-9);
    for(int level=1;level<=bottom_hierarchy->Levels();level++) SPGrid_Computations::Clear<T_STRUCT,T,d>(bottom_hierarchy->Allocator(level),bottom_hierarchy->Blocks(level),x_channel);
    BUILD_MATRIX<T_STRUCT,T,d>::Vector_To_SPGrid_Channel(active_dof_array,*bottom_hierarchy,x_channel,vecX);}
#else
    {LOG::SCOPE scope("MULTIGRID_SOLVER::BOTTOM_SOLVE");
    ARRAY<int> substitution_partitions_at_bottom(substitution_partitions);substitution_partitions_at_bottom.Resize(multigrid_hierarchy(multigrid_levels)->Levels());
    CG_SYSTEM<T_STRUCT,T,d> cg_system(*multigrid_hierarchy(multigrid_levels),L_channels,ic_diag_channel,substitution_partitions_at_bottom);
    CONJUGATE_GRADIENT<T> cg;
    for(int level=1;level<=multigrid_hierarchy(multigrid_levels)->Levels();level++){
        if(PhysBAM_number_of_threads){
            SPGrid_Computations::Threading_Helper<T_STRUCT,d>(multigrid_hierarchy(multigrid_levels)->Allocator(level),multigrid_hierarchy(multigrid_levels)->Blocks(level)).Run_Parallel(
                SPGrid_Computations::Clear<T_STRUCT,T,d,2>(cg_q_channel,cg_r_channel),PhysBAM_number_of_threads);
            SPGrid_Computations::Threading_Helper<T_STRUCT,d>(multigrid_hierarchy(multigrid_levels)->Allocator(level),multigrid_hierarchy(multigrid_levels)->Blocks(level)).Run_Parallel(
                SPGrid_Computations::Clear<T_STRUCT,T,d,2>(cg_s_channel,cg_k_channel),PhysBAM_number_of_threads);
            SPGrid_Computations::Threading_Helper<T_STRUCT,d>(multigrid_hierarchy(multigrid_levels)->Allocator(level),multigrid_hierarchy(multigrid_levels)->Blocks(level)).Run_Parallel(
                SPGrid_Computations::Clear<T_STRUCT,T,d,2>(cg_z_channel,x_channel),PhysBAM_number_of_threads);}
        else{SPGrid_Computations::Clear<T_STRUCT,T,d,2>(multigrid_hierarchy(multigrid_levels)->Allocator(level),multigrid_hierarchy(multigrid_levels)->Blocks(level),cg_q_channel,cg_r_channel);
            SPGrid_Computations::Clear<T_STRUCT,T,d,2>(multigrid_hierarchy(multigrid_levels)->Allocator(level),multigrid_hierarchy(multigrid_levels)->Blocks(level),cg_s_channel,cg_k_channel);
            SPGrid_Computations::Clear<T_STRUCT,T,d,2>(multigrid_hierarchy(multigrid_levels)->Allocator(level),multigrid_hierarchy(multigrid_levels)->Blocks(level),cg_z_channel,x_channel);}}
    CG_VECTOR<T_STRUCT,T,d> x_V(*multigrid_hierarchy(multigrid_levels),x_channel);
    CG_VECTOR<T_STRUCT,T,d> b_V(*multigrid_hierarchy(multigrid_levels),b_channel);
    CG_VECTOR<T_STRUCT,T,d> q_V(*multigrid_hierarchy(multigrid_levels),cg_q_channel);
    CG_VECTOR<T_STRUCT,T,d> r_V(*multigrid_hierarchy(multigrid_levels),cg_r_channel);
    CG_VECTOR<T_STRUCT,T,d> s_V(*multigrid_hierarchy(multigrid_levels),cg_s_channel);
    CG_VECTOR<T_STRUCT,T,d> k_V(*multigrid_hierarchy(multigrid_levels),cg_k_channel);
    CG_VECTOR<T_STRUCT,T,d> z_V(*multigrid_hierarchy(multigrid_levels),cg_z_channel);
    const int cg_iterations=1000;const int cg_restart_iterations=cg_iterations+1;
    cg.print_residuals=false;cg.print_diagnostics=false;
    cg.restart_iterations=cg_restart_iterations;
    cg.Solve(cg_system,x_V,b_V,q_V,s_V,r_V,k_V,z_V,(T)1e-5,0,cg_iterations);}
#endif
#else
    //PHYSBAM_FATAL_ERROR("Dont use smoother at bottom!");
    SMOOTHER<T_STRUCT,T,d>::Single_Level_Boundary_Jacobi_Smoother(*multigrid_hierarchy(multigrid_levels),x_channel,b_channel,mg_diag_channel,temp_channel,200,boundary_omega,(unsigned)(MG_Boundary),diagonal_is_inverted);
    SMOOTHER<T_STRUCT,T,d>::Single_Level_Interior_Jacobi_Smoother(*multigrid_hierarchy(multigrid_levels),x_channel,b_channel,mg_diag_channel,temp_channel,200,interior_omega,diagonal_is_inverted);
    SMOOTHER<T_STRUCT,T,d>::Single_Level_Boundary_Jacobi_Smoother(*multigrid_hierarchy(multigrid_levels),x_channel,b_channel,mg_diag_channel,temp_channel,200,boundary_omega,(unsigned)(MG_Boundary),diagonal_is_inverted);
#endif
    // upstroke
    for(int i=multigrid_levels-1;i>=1;i--){
        MULTIGRID_REFINEMENT<T_STRUCT,T,d>::Prolongate(*multigrid_hierarchy(i),*multigrid_hierarchy(i+1),cg_q_channel,x_channel,VECTOR<int,2>(i,i+1));
        for(int level=1;level<=multigrid_hierarchy(i)->Levels();level++)
            if(PhysBAM_number_of_threads) SPGrid_Computations::Threading_Helper<T_STRUCT,d>(multigrid_hierarchy(i)->Allocator(level),multigrid_hierarchy(i)->Blocks(level)).Run_Parallel(
                SPGrid_Computations::Masked_Add<T_STRUCT,T,unsigned,d>(x_channel,cg_q_channel,x_channel,&T_STRUCT::flags,(unsigned)SPGrid_Cell_Type_Active),PhysBAM_number_of_threads);
            else SPGrid_Computations::Masked_Add<T_STRUCT,T,unsigned,d>(multigrid_hierarchy(i)->Allocator(level),multigrid_hierarchy(i)->Blocks(level),x_channel,cg_q_channel,x_channel,&T_STRUCT::flags,(unsigned)SPGrid_Cell_Type_Active);
        GRID_HIERARCHY_PROJECTION<T_STRUCT,T,d>::Propagate_Ghost_Values(*multigrid_hierarchy(i),&T_STRUCT::flags,x_channel);
        SMOOTHER<T_STRUCT,T,d>::Single_Level_Jacobi_Smoother(*multigrid_hierarchy(i),x_channel,b_channel,temp_channel,mg_diag_channel,i,boundary_smoother_iterations,boundary_omega,(unsigned)(MG_Boundary),multigrid_hierarchy(i)->Boundary_Blocks(i));
        SMOOTHER<T_STRUCT,T,d>::Single_Level_Jacobi_Smoother(*multigrid_hierarchy(i),x_channel,b_channel,temp_channel,mg_diag_channel,i,interior_smoother_iterations,interior_omega,(unsigned)(SPGrid_Cell_Type_Active),multigrid_hierarchy(i)->Blocks(i));
        SMOOTHER<T_STRUCT,T,d>::Single_Level_Jacobi_Smoother(*multigrid_hierarchy(i),x_channel,b_channel,temp_channel,mg_diag_channel,i,boundary_smoother_iterations,boundary_omega,(unsigned)(MG_Boundary),multigrid_hierarchy(i)->Boundary_Blocks(i));}
}
//#####################################################################
template class MULTIGRID_SOLVER<FLUIDS_SIMULATION_DATA<float>,float,2>;
template class MULTIGRID_SOLVER<FLUIDS_SIMULATION_DATA<float>,float,3>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class MULTIGRID_SOLVER<FLUIDS_SIMULATION_DATA<double>,double,2>;
template class SMOOTHER<FLUIDS_SIMULATION_DATA<double>,double,3>;
#endif
