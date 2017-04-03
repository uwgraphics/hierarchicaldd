//#####################################################################
// Copyright 2012, Mridul Aanjaneya.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifdef USE_MPI
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Parallel_Computation/PCG_SPARSE_BUBBLE_SOLID_FLUID_MPI.h>
using namespace PhysBAM;
//#####################################################################
// Function Parallel_Solve
//#####################################################################
//template<class T_GRID> void PCG_SPARSE_BUBBLE_SOLID_FLUID_MPI<T_GRID>::
//Parallel_Solve(SPARSE_MATRIX_FLAT_NXN<T>& A,VECTOR_ND<T>& x,VECTOR_ND<T>& b,const T tolerance,const bool recompute_preconditioner)
//{
//    Initialize_Datatypes();
//    int local_n=A.n,interior_n=partition.interior_indices.Size()+1;
//    int global_n=Global_Sum_Bubbles(interior_n);
//    T global_tolerance=Global_Max_Bubbles(tolerance);
//    int desired_iterations=global_n;if(pcg.enforce_compatibility) desired_iterations--;if(pcg.maximum_iterations) desired_iterations=min(desired_iterations,pcg.maximum_iterations);
//
//    VECTOR_ND<T> temp(local_n,false),p(local_n,false),z_interior(interior_n,false);
//
//    // build interior views of x,b,p,z,temp
//    VECTOR_ND<T> x_interior,b_interior,p_interior,temp_interior;
//    x_interior.Set_Subvector_View(x,partition.interior_indices);
//    b_interior.Set_Subvector_View(b,partition.interior_indices);
//    p_interior.Set_Subvector_View(p,partition.interior_indices);
//    temp_interior.Set_Subvector_View(temp,partition.interior_indices);
//
//    // adjust x for the null space
//    if(pcg.enforce_compatibility && pcg.remove_null_space_solution_component) x_interior-=(T)(Global_Sum_Bubbles(x_interior.Sum_Double_Precision())/global_n);
//
//    // find initial residual, r=b-Ax - reusing b for the residual
//    Fill_Ghost_Cells(x);
//    A.Times(x,temp);b_interior-=temp_interior;
//    if(pcg.enforce_compatibility) b_interior-=(T)(Global_Sum_Bubbles(b_interior.Sum_Double_Precision())/global_n);
//    if(Global_Max_Bubbles(b_interior.Max_Abs())<=global_tolerance){
//#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
//        if(pcg.show_results) LOG::filecout("NO ITERATIONS NEEDED\n");
//#endif
//        return;}
//
//    // find an incomplete cholesky preconditioner - actually an LU that saves square roots, and an inverted diagonal to save on divides
//    if(pcg.incomplete_cholesky && (recompute_preconditioner || !A.C)){
//        delete A.C;A.C=A.Create_Submatrix(partition.interior_indices);
//        A.C->In_Place_Incomplete_Cholesky_Factorization(pcg.modified_incomplete_cholesky,pcg.modified_incomplete_cholesky_coefficient,
//            pcg.preconditioner_zero_tolerance,pcg.preconditioner_zero_replacement);}
//
//    double rho=0,rho_old=0;
//    for(int iteration=1;;iteration++){
//        if(pcg.incomplete_cholesky){
//            // solve Mz=r
//            A.C->Solve_Forward_Substitution(b_interior,temp_interior,true); // diagonal should be treated as the identity
//            A.C->Solve_Backward_Substitution(temp_interior,z_interior,false,true);} // diagonal is inverted to save on divides
//        else z_interior=b_interior; // set z=r when there is no preconditioner
//
//        // for Neumann boundary conditions only, make sure z sums to zero
//        if(pcg.enforce_compatibility) z_interior-=(T)(Global_Sum_Bubbles(z_interior.Sum_Double_Precision())/global_n);
//
//        // update search direction
//        rho_old=rho;rho=Global_Sum_Bubbles(VECTOR_ND<T>::Dot_Product_Double_Precision(z_interior,b_interior));
//        T beta=0;if(iteration==1) p_interior=z_interior;else{beta=(T)(rho/rho_old);for(int i=1;i<=interior_n;i++) p_interior(i)=z_interior(i)+beta*p_interior(i);} // when iteration=1, beta=0
//
//        // update solution and residual
//        Fill_Ghost_Cells(p);
//        A.Times(p,temp);
//        T alpha=(T)(rho/Global_Sum_Bubbles(VECTOR_ND<T>::Dot_Product_Double_Precision(p_interior,temp_interior)));
//        for(int i=1;i<=interior_n;i++){x_interior(i)+=alpha*p_interior(i);b_interior(i)-=alpha*temp_interior(i);}
//
//        // remove null space component of b before computing residual norm because we might have converged up to the null space but have some null space component left due to roundoff
//        if(pcg.enforce_compatibility) b_interior-=(T)(Global_Sum_Bubbles(b_interior.Sum_Double_Precision())/global_n);
//
//#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
//        T residual=Global_Max_Bubbles(b_interior.Max_Abs());
//
//        // check for convergence
//        std::stringstream ss;
//        if(pcg.show_residual) ss<<residual<<std::endl;
//        if(residual<=global_tolerance){if(pcg.show_results) ss<<"NUMBER OF ITERATIONS = "<<iteration<<std::endl;LOG::filecout(ss.str());break;}
//        if(iteration==desired_iterations){if(pcg.show_results) ss<<"DID NOT CONVERGE IN "<<iteration<<" ITERATIONS"<<std::endl;LOG::filecout(ss.str());break;}
//#endif
//    }
//  
//    Fill_Ghost_Cells(x);
//}
//#####################################################################
template class PCG_SPARSE_BUBBLE_SOLID_FLUID_MPI<GRID<VECTOR<float,1> > >;
template class PCG_SPARSE_BUBBLE_SOLID_FLUID_MPI<GRID<VECTOR<float,2> > >;
template class PCG_SPARSE_BUBBLE_SOLID_FLUID_MPI<GRID<VECTOR<float,3> > >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class PCG_SPARSE_BUBBLE_SOLID_FLUID_MPI<GRID<VECTOR<double,1> > >;
template class PCG_SPARSE_BUBBLE_SOLID_FLUID_MPI<GRID<VECTOR<double,2> > >;
template class PCG_SPARSE_BUBBLE_SOLID_FLUID_MPI<GRID<VECTOR<double,3> > >;
#endif
#endif
