//#####################################################################
//Removed For Anonymity: Copyright Authors of Submission pap173s1
// This file is part of PhysBAM whose distribution is governed by the license 
//Removed For Anonymity: Copyright Authors of Submission pap173s1
//#####################################################################
// Class Prolongation_Matrix_Helper.h
//#####################################################################
#ifndef __PROLONGATION_MATRIX_HELPER_H__
#define __PROLONGATION_MATRIX_HELPER_H__
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Data_Structures/HASHTABLE.h>
#include <Eigen/Sparse>
#include <vector>
namespace PhysBAM{
template<typename T,int d>
class PROLONGATION_MATRIX_HELPER{ 
    typedef VECTOR<int,d> T_INDEX;
    typedef Eigen::Triplet<T,int> TRIPLET;
public:
    static void Construct_Prolongation_Matrix(Eigen::SparseMatrix<T,Eigen::RowMajor,int>& prolongation_matrix,
                                              const HASHTABLE<T_INDEX,int>&index_map_nd_to_1d_fine,
                                              const HASHTABLE<T_INDEX,int>&index_map_nd_to_1d_coarse);
};
}
#endif
