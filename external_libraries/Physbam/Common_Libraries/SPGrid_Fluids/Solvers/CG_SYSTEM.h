//#####################################################################
// Copyright 2011, Taylor Patterson, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class CG_SYSTEM
//#####################################################################
#ifndef __CG_SYSTEM__
#define __CG_SYSTEM__

#include <PhysBAM_Tools/Krylov_Solvers/KRYLOV_SYSTEM_BASE.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Vectors/VECTOR_2D.h>
#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>
#include <PhysBAM_Tools/Vectors/VECTOR_ND.h>
#include <PhysBAM_Tools/Matrices/SPARSE_MATRIX_FLAT_NXN.h>

#include <SPGrid_Fluids/Solvers/CG_VECTOR.h>

namespace PhysBAM{

template<class T_STRUCT,class T,int d> class NONLINEAR_ELASTICITY;

template<class T_STRUCT,class T,int d>
class CG_SYSTEM:public KRYLOV_SYSTEM_BASE<T>
{
    typedef KRYLOV_SYSTEM_BASE<T> BASE;
    typedef KRYLOV_VECTOR_BASE<T> VECTOR_BASE;
    typedef SPARSE_MATRIX_FLAT_NXN<T> T_MATRIX;
    typedef VECTOR_ND<T> T_VECTOR;
    

    typedef GRID_HIERARCHY<T_STRUCT,T,d> Hierarchy_type;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<T>::type Data_array_type;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<const T>::type Const_data_array_type;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<unsigned>::type Flag_array_type;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<const unsigned>::type Const_flag_array_type;
    Hierarchy_type& hierarchy;
    T_MATRIX A;
    
    VECTOR<T T_STRUCT::*,d> L_channels;
    T T_STRUCT::* diag_channel;

    ARRAY<Hierarchy_type*>* multigrid_hierarchy;
    T T_STRUCT::* mg_residual_channel;
    T T_STRUCT::* mg_diag_channel;
    T T_STRUCT::* mg_temp_channel;
    T T_STRUCT::* mg_x_channel;
    T T_STRUCT::* mg_b_channel;
    T T_STRUCT::* mg_cg_q_channel;
    T T_STRUCT::* mg_cg_r_channel;
    T T_STRUCT::* mg_cg_s_channel;
    T T_STRUCT::* mg_cg_k_channel;
    T T_STRUCT::* mg_cg_z_channel;
    const int interior_smoother_iterations;
    const int boundary_smoother_iterations;
    const T interior_omega;
    const T boundary_omega;
    const bool mg_diag_is_inverted;

    bool ic_preconditioner;

    T T_STRUCT::* variable_beta_channel;
    bool use_variable_beta;

    const ARRAY<int>& substitution_partitions;

//#####################################################################
public:
    CG_SYSTEM(Hierarchy_type& hierarchy_input,VECTOR<T T_STRUCT::*,d> L_channels_in,T T_STRUCT::* diag_channel_in,const ARRAY<int>& substitution_partitions_in); // IC
    CG_SYSTEM(Hierarchy_type& hierarchy_input,ARRAY<Hierarchy_type*>* multigrid_hierarchy_in,T T_STRUCT::* mg_residual_channel_in, // MG
        T T_STRUCT::* mg_diag_channel_in,T T_STRUCT::* mg_temp_channel_in,T T_STRUCT::* mg_x_channel_in,T T_STRUCT::* mg_b_channel_in,T T_STRUCT::* mg_cg_q_channel_in,
        T T_STRUCT::* mg_cg_r_channel_in,T T_STRUCT::* mg_cg_s_channel_in,T T_STRUCT::* mg_cg_k_channel_in,
        T T_STRUCT::* mg_cg_z_channel_in,const int interior_smoother_iterations_in,
        const int boundary_smoother_iterations_in,const T interior_omega_in,const T boundary_omega_in,
        const bool mg_diag_is_inverted_in,VECTOR<T T_STRUCT::*,d> L_channels_in,T T_STRUCT::* diag_channel_in,
        const ARRAY<int>& substitution_partitions_in);
    void Multiply_With_Matrix(const VECTOR_BASE& v,VECTOR_BASE& result) const;
    void Multiply(const VECTOR_BASE& v,VECTOR_BASE& result) const;
    double Inner_Product(const VECTOR_BASE& x,const VECTOR_BASE& y) const;
    T Convergence_Norm(const VECTOR_BASE& x) const;
    void Project(VECTOR_BASE& x) const;
    void Set_Boundary_Conditions(VECTOR_BASE& x) const;
    void Project_Nullspace(VECTOR_BASE& x) const;
    void Use_Variable_Beta(T T_STRUCT::* variable_beta_channel_input);
    void Use_Constant_Beta();
protected:
    void Apply_Preconditioner(const VECTOR_BASE& r, VECTOR_BASE& z) const;
//#####################################################################
};
}
#endif
