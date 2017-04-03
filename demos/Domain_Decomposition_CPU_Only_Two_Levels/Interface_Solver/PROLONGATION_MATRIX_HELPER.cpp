//#####################################################################
//Removed For Anonymity: Copyright Authors of Submission pap173s1
// This file is part of PhysBAM whose distribution is governed by the license 
//Removed For Anonymity: Copyright Authors of Submission pap173s1
//#####################################################################
#include <Common_Tools/Math_Tools/RANGE_ITERATOR.h>
#include <PhysBAM_Tools/Data_Structures/HASHTABLE_ITERATOR.h>
#include "PROLONGATION_MATRIX_HELPER.h"
using namespace PhysBAM;
//#####################################################################
// Function Construct_Prolongation_Matrix
//#####################################################################
template<typename T,int d> void
PROLONGATION_MATRIX_HELPER<T,d>::Construct_Prolongation_Matrix(Eigen::SparseMatrix<T,Eigen::RowMajor,int>& prolongation_matrix,
                                                               const HASHTABLE<T_INDEX,int>& index_map_nd_to_1d_fine,
                                                               const HASHTABLE<T_INDEX,int>& index_map_nd_to_1d_coarse){
    int fine_dof = index_map_nd_to_1d_fine.Size();
    int coarse_dof = index_map_nd_to_1d_coarse.Size();
    prolongation_matrix.resize(fine_dof,coarse_dof);
    std::vector<TRIPLET> triplet_list;
    for(HASHTABLE_ITERATOR<T_INDEX,const int> iterator(index_map_nd_to_1d_fine);iterator.Valid();iterator.Next()){
        const T_INDEX& fine_index = iterator.Key();
        T_INDEX coarse_node; 
        for(int v = 1;v <= d;++v) coarse_node(v) = std::floor(fine_index(v) / 2.0);
        T total_weight = 0;
        //iterate through the neighbor of the corresponding coarse node
        for(RANGE_ITERATOR<d> neighbor_iterator(RANGE<T_INDEX>(coarse_node,coarse_node+1));neighbor_iterator.Valid();neighbor_iterator.Next()){
            T_INDEX coarse_index = neighbor_iterator.Index();
            T weight = T(1);
            for(int v = 1;v <= d;++v) weight *= (T(2) - abs(coarse_index(v) * 2 - fine_index(v))) / 2.0;
            total_weight += weight;
            if(weight > 0){
                int coarse_index_1D=index_map_nd_to_1d_coarse.Get_Default(coarse_index,-1);
                if(coarse_index_1D==-1) continue;
                PHYSBAM_ASSERT(iterator.Data() != -1);
                triplet_list.push_back(TRIPLET(iterator.Data(),coarse_index_1D,weight));}}
        PHYSBAM_ASSERT(total_weight == 1.0);}
    prolongation_matrix.setFromTriplets(triplet_list.begin(),triplet_list.end());
    prolongation_matrix.makeCompressed();
}
//#####################################################################
template class PROLONGATION_MATRIX_HELPER<float,2>;
template class PROLONGATION_MATRIX_HELPER<float,3>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class PROLONGATION_MATRIX_HELPER<double,2>;
template class PROLONGATION_MATRIX_HELPER<double,3>;
#endif
