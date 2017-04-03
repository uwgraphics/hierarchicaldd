//#####################################################################
//Removed For Anonymity: Copyright Authors of Submission pap173s1
// This file is part of PhysBAM whose distribution is governed by the license 
//Removed For Anonymity: Copyright Authors of Submission pap173s1
//#####################################################################
#ifndef __SIGMA_SMOOTHER_H__
#define __SIGMA_SMOOTHER_H__
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <Common_Tools/Math_Tools/RANGE_ITERATOR.h>
#include "INTERFACE_SOLVER_DATA.h"
#include <Eigen/Sparse>

namespace PhysBAM{

template<typename T,int d>
class SIGMA_SMOOTHER{
    typedef VECTOR<int,d> T_INDEX;
    //In Eigen 3, all dense matrices and vectors are indexed with pointer sized int.
    typedef Eigen::Matrix<T,Eigen::Dynamic,1> T_EIGEN_VECTOR;
public:
    INTERFACE_SOLVER_DATA<T,d>& solver_data;
    bool symbolic_factorized;
    //##################################################################### 
    SIGMA_SMOOTHER(INTERFACE_SOLVER_DATA<T,d>& solver_data_in)
        :solver_data(solver_data_in),symbolic_factorized(false){};
    ~SIGMA_SMOOTHER(){};
    void Symbolic_Factorization();
    void Numerical_Factorization();
    void Solve_Aii(T_EIGEN_VECTOR& result,T_EIGEN_VECTOR& rhs,T_INDEX subdomain_index=T_INDEX(),int level=0);
    //the smooth function replace u in place.
    void Smooth(const int level,T_EIGEN_VECTOR& u,T_EIGEN_VECTOR& u_next,const T_EIGEN_VECTOR& b)const; 
    void Compute_Residual(const int level,const T_EIGEN_VECTOR& u,T_EIGEN_VECTOR& r,const T_EIGEN_VECTOR& b)const;
    void Multiply(const int level,const T_EIGEN_VECTOR& u,T_EIGEN_VECTOR& result)const;
    //##################################################################### 
};
}
#endif
