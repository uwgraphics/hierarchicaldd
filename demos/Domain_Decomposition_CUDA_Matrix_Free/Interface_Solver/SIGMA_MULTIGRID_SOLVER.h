//#####################################################################
//Removed For Anonymity: Copyright Authors of Submission pap173s1
// This file is part of PhysBAM whose distribution is governed by the license 
//Removed For Anonymity: Copyright Authors of Submission pap173s1
//#####################################################################
#ifndef __SIGMA_MULTIGRID_SOLVER_H__
#define __SIGMA_MULTIGRID_SOLVER_H__
#include "SIGMA_SMOOTHER.h"

namespace PhysBAM{
template<typename T,int d>
class SIGMA_MULTIGRID_SOLVER{
    typedef VECTOR<int,d> T_INDEX;
    //In Eigen 3, all dense matrices and vectors are indexed with pointer sized int.
    typedef Eigen::SparseMatrix<T,Eigen::RowMajor,int> T_PROLONGATION_MATRIX;
    typedef Eigen::Matrix<T,Eigen::Dynamic,1> T_EIGEN_VECTOR;
public:
    INTERFACE_SOLVER_DATA<T,d>& solver_data;
    SIGMA_SMOOTHER<T,d> smoother;
    std::vector<T_EIGEN_VECTOR>& u;
    std::vector<T_EIGEN_VECTOR>& r;
    std::vector<T_EIGEN_VECTOR>& b;
    const int smoothing_iteration;
    SIGMA_MULTIGRID_SOLVER(INTERFACE_SOLVER_DATA<T,d>& solver_data_in,
                           std::vector<T_EIGEN_VECTOR>& u_in,
                           std::vector<T_EIGEN_VECTOR>& r_in,
                           std::vector<T_EIGEN_VECTOR>& b_in,
                           const int smoothing_iteration_in=3)
        :solver_data(solver_data_in),u(u_in),r(r_in),b(b_in),smoother(solver_data_in),smoothing_iteration(smoothing_iteration_in){};
    void Symbolic_Initialize();
    void Numerical_Initialize();
    void V_Cycle();
};
}
#endif
