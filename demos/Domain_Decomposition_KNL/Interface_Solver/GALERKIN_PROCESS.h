//#####################################################################
//Removed For Anonymity: Copyright Authors of Submission pap173s1
// Distributed under the FreeBSD license (see license.txt)
//#####################################################################
#ifndef __GALERKIN_PROCESS_H__
#define __GALERKIN_PROCESS_H__

#include <Common_Tools/Math_Tools/RANGE_ITERATOR.h>
#include <SPGrid_Fluids/Grids/GRID_HIERARCHY.h>
#include <Eigen/Core>

namespace PhysBAM{
using namespace SPGrid;
//#####################################################################
// Class GALERKIN_PROCESS
//#####################################################################
template<class T_STRUCT_SOLVER,class T,int d>
class GALERKIN_PROCESS{
    typedef typename SPGrid_Allocator<T_STRUCT_SOLVER,d>::template Array<unsigned>::type Flag_array_type_solver;
    typedef typename SPGrid_Allocator<T_STRUCT_SOLVER,d>::template Array<unsigned>::mask Flag_array_mask_solver;
    typedef SPGrid_Set<Flag_array_type_solver> Set_type_solver;

    typedef std_array<int,d> coord_t;
    typedef VECTOR<int,d> T_INDEX;
// #############################################################################
public:
    static Eigen::Matrix<T,1<<d,1<<d> Galerkin_Coarsen(const Set_type_solver& set,
                                                       const std::vector<Eigen::Matrix<T,1<<d,1<<d>,Eigen::aligned_allocator<std::pair<const int, Eigen::Matrix<T,1<<d,1<<d> > > >& prolongation_matrices,
                                                       const ARRAY<int,T_INDEX>& index_map,const T_INDEX& cell_index,int level){
        static T edge_weight=1.0/(1<<(d-1));
        Eigen::Matrix<T,1<<d,1<<d> matrix;
        matrix.setZero();
        if(level==1){
            std::array<unsigned,1<<d> node_flag;
            for(RANGE_ITERATOR<d> node_iterator(RANGE<VECTOR<int,d> >(cell_index-1,cell_index));node_iterator.Valid();node_iterator.Next()){
                const T_INDEX& node_index=node_iterator.Index();
                node_flag[index_map(node_index+1-cell_index)]=0;
                unsigned long node_offset=Flag_array_mask_solver::Linear_Offset(std_array<int,d>(node_index));
                if(set.Is_Set(node_offset,SPGrid_Solver_Cell_Type_Active|SPGrid_Solver_Cell_Type_Interface))
                    node_flag[index_map(node_index+1-cell_index)]=SPGrid_Solver_Cell_Type_Active;
                else if(set.Is_Set(node_offset,SPGrid_Solver_Cell_Type_Dirichlet)){
                    node_flag[index_map(node_index+1-cell_index)]=SPGrid_Solver_Cell_Type_Dirichlet;}}
            for(int v=1;v<=d;v++)
                for(RANGE_ITERATOR<d-1> face_iterator(RANGE<VECTOR<int,d-1> >(cell_index.Remove_Index(v)-1,cell_index.Remove_Index(v)));
                    face_iterator.Valid();face_iterator.Next()){
                    T_INDEX lower_index(face_iterator.Index().Insert(cell_index(v)-1,v));
                    T_INDEX upper_index(face_iterator.Index().Insert(cell_index(v),v));
                    unsigned long lower_offset=Flag_array_mask_solver::Linear_Offset(std_array<int,d>(lower_index));
                    unsigned long upper_offset=Flag_array_mask_solver::Linear_Offset(std_array<int,d>(upper_index));
                    const int lower_id=index_map(lower_index+1-cell_index);
                    const int upper_id=index_map(upper_index+1-cell_index);
                    if((node_flag[lower_id]&SPGrid_Solver_Cell_Type_Active)&&(node_flag[upper_id]&SPGrid_Solver_Cell_Type_Active)){
                        matrix(lower_id,upper_id)-=edge_weight;
                        matrix(lower_id,lower_id)+=edge_weight;
                        matrix(upper_id,upper_id)+=edge_weight;
                        matrix(upper_id,lower_id)-=edge_weight;}
                    if((node_flag[lower_id]&SPGrid_Solver_Cell_Type_Active)&&(node_flag[upper_id]&SPGrid_Solver_Cell_Type_Dirichlet)){
                        matrix(lower_id,lower_id)+=edge_weight;}
                    if((node_flag[lower_id]&SPGrid_Solver_Cell_Type_Dirichlet)&&(node_flag[upper_id]&SPGrid_Solver_Cell_Type_Active)){
                        matrix(upper_id,upper_id)+=edge_weight;}
                }
        }else{
            T_INDEX children_base_index=2*cell_index-1;
            for(RANGE_ITERATOR<d> children_iterator(RANGE<T_INDEX>::Unit_Box());children_iterator.Valid();children_iterator.Next()){
                const T_INDEX& children_index=children_iterator.Index();
                matrix+=prolongation_matrices[index_map(children_index)].transpose()*Galerkin_Coarsen(set,prolongation_matrices,index_map,children_base_index+children_index,level-1)*prolongation_matrices[index_map(children_index)];
            }      
        }
        return matrix;
    }
// #############################################################################
};
}
#endif

