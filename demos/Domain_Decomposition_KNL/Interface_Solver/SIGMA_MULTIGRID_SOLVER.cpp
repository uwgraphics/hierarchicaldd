//#####################################################################
//Removed For Anonymity: Copyright Authors of Submission pap173s1
// This file is part of PhysBAM whose distribution is governed by the license 
//Removed For Anonymity: Copyright Authors of Submission pap173s1
//#####################################################################
#include "SIGMA_MULTIGRID_SOLVER.h"
using namespace PhysBAM;

//#####################################################################
// Function Symbolic_Initialize
//#####################################################################
template<typename T,int d> void SIGMA_MULTIGRID_SOLVER<T,d>::
Symbolic_Initialize(){
    smoother.Symbolic_Factorization();
    const int nlevels=solver_data.Get_Number_Of_Levels();
    u.resize(nlevels);
    b.resize(nlevels);
    r.resize(nlevels);
    
    for(int level=0;level<nlevels;++level){
        u[level].resize(solver_data.Arr[level].rows());
        b[level].resize(solver_data.Arr[level].rows());
        r[level].resize(solver_data.Arr[level].rows());}
}
//#####################################################################
// Function Numerical_Initialize
//#####################################################################
template<typename T,int d> void SIGMA_MULTIGRID_SOLVER<T,d>::
Numerical_Initialize(){
    smoother.Numerical_Factorization();
}
//#####################################################################
// Function V_Cycle
//#####################################################################
template<typename T,int d> void SIGMA_MULTIGRID_SOLVER<T,d>::
V_Cycle(){
    const int nlevels=solver_data.Get_Number_Of_Levels();
    for(int level=0;level<nlevels-1;++level){
        for(int i=0;i<smoothing_iteration;++i) smoother.Smooth(level,u[level],r[level],b[level]);
        smoother.Compute_Residual(level,u[level],r[level],b[level]);
        u[level+1]=T_EIGEN_VECTOR::Zero(u[level+1].size());
        b[level+1]=solver_data.Prr[level].transpose()*r[level];}
    for(int i=0;i<200;++i) smoother.Smooth(nlevels-1,u[nlevels-1],r[nlevels-1],b[nlevels-1]);
    for(int level=nlevels-2;level>=0;--level){
        u[level]+=solver_data.Prr[level]*u[level+1];        
        for(int i=0;i<smoothing_iteration;++i) smoother.Smooth(level,u[level],r[level],b[level]);}
}
//#####################################################################
template class SIGMA_MULTIGRID_SOLVER<float,2>;
template class SIGMA_MULTIGRID_SOLVER<float,3>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class SIGMA_MULTIGRID_SOLVER<double,2>;
template class SIGMA_MULTIGRID_SOLVER<double,3>;
#endif

