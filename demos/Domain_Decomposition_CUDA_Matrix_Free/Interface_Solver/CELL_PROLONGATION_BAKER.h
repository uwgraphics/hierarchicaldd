//#####################################################################
//Removed For Anonymity: Copyright Authors of Submission pap173s1
// This file is part of PhysBAM whose distribution is governed by the license 
//Removed For Anonymity: Copyright Authors of Submission pap173s1
//#####################################################################
// Class CELL_PROLONGATION_BAKER.h
//#####################################################################
#ifndef __CELL_PROLONGATION_BAKER_H__
#define __CELL_PROLONGATION_BAKER_H__
#include <Common_Tools/Math_Tools/RANGE_ITERATOR.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <Eigen/Core>
#include <vector>
namespace PhysBAM{
template<typename T,int d>
class CELL_PROLONGATION_BAKER{
public:
    typedef VECTOR<int,d> T_INDEX;
    std::vector<Eigen::Matrix<T,1<<d,1<<d> > prolongation_matrices;
    CELL_PROLONGATION_BAKER():prolongation_matrices(1<<d){}
    ARRAY<int,T_INDEX> index_map;
    void Bake(){
        index_map.Resize(RANGE<T_INDEX>::Unit_Box());
        int index_counter=0;
        for(RANGE_ITERATOR<d> iterator(RANGE<T_INDEX>::Unit_Box());iterator.Valid();iterator.Next(),index_counter++){
            index_map(iterator.Index())=index_counter;}
        int matrix_counter=0;
        for(RANGE_ITERATOR<d> matrix_iterator(RANGE<T_INDEX>::Unit_Box());matrix_iterator.Valid();matrix_iterator.Next(),matrix_counter++){
            prolongation_matrices[matrix_counter].setZero();
            const T_INDEX& fine_offset=matrix_iterator.Index();
            int fine_index=0;
            for(RANGE_ITERATOR<d> fine_iterator(RANGE<T_INDEX>::Unit_Box());fine_iterator.Valid();fine_iterator.Next(),fine_index++){
                const T_INDEX fine_node_index=fine_iterator.Index()+fine_offset;
                T_INDEX coarse_node; 
                for(int v=1;v<=d;++v) coarse_node(v)=std::floor(fine_node_index(v)/2.0);
                T total_weight = 0;
                //iterate through the neighbor of the corresponding coarse node
                for(RANGE_ITERATOR<d> neighbor_iterator(RANGE<T_INDEX>(coarse_node,coarse_node+1));neighbor_iterator.Valid();neighbor_iterator.Next()){
                    const T_INDEX& coarse_node_index = neighbor_iterator.Index();
                    T weight=T(1);
                    for(int v=1;v<=d;++v) weight*=(T(2)-abs(coarse_node_index(v)*2-fine_node_index(v)))/2.0;
                    total_weight+=weight;
                    for(int v=1;v<=d;++v) if(coarse_node_index(v)>1) PHYSBAM_ASSERT(weight==0);
                    if(weight==0) continue;
                    prolongation_matrices[matrix_counter](fine_index,index_map(coarse_node_index))=weight;
                }
                PHYSBAM_ASSERT(total_weight==1);
            }
        }
    }
};
}
#endif
