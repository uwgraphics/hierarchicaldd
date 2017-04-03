//#####################################################################
//Removed For Anonymity: Copyright Authors of Submission pap173s1
// This file is part of PhysBAM whose distribution is governed by the license 
//Removed For Anonymity: Copyright Authors of Submission pap173s1
//#####################################################################
#include <Common_Tools/Grids_Uniform_PDE_Linear/STENCIL_ITERATOR.h>
#include <Common_Tools/Math_Tools/RANGE_ITERATOR.h>
#include <PhysBAM_Tools/Data_Structures/HASHTABLE_ITERATOR.h>
#include <SPGrid_Fluids/Simulation/FLUIDS_SIMULATION_DATA.h>
#include <SPGrid_Fluids/Solvers/Domain_Decomposition/SPGRID_MULTIGRID_FLAGS.h>
#include "../SPGrid_Domain_Decomposition_Solver/SPGRID_DOMAIN_DECOMPOSITION_DATA.h"
#include "EIGEN_INTERFACE.h"
using namespace PhysBAM;
/////////////////////////////////////////////////////
// Function Copy_To_Eigen_Array
/////////////////////////////////////////////////////
template<typename T_STRUCT,typename T,int d> void 
EIGEN_INTERFACE<T_STRUCT,T,d>::Copy_To_Eigen_Array(T_EIGEN_VECTOR& array_out,
                                                   const T T_STRUCT::* u_field,const T_FLAG T_STRUCT::* flags_field, 
                                                   /*const*/ SPG_Allocator& allocator,/*const*/ SPG_Set_Type& set,
                                                   const std::vector<T_INDEX>& index_map_1d_to_nd){
    SPG_Const_Data_Array_Type u=allocator.Get_Const_Array(u_field);
    SPG_Const_Flags_Array_Type flags=allocator.Get_Const_Array(flags_field);
    
    #pragma omp parallel for
    for(int i=0;i<index_map_1d_to_nd.size();++i){
        //std::cout<<omp_get_num_threads()<<std::endl;
        unsigned long offset=T_MASK::Linear_Offset(std_array<int,d>(index_map_1d_to_nd[i]));
        if(set.IsPageActive(offset)&&(flags(offset)&SPGrid_Solver_Cell_Type_Interface))
            array_out(i)=u(offset);}
}
/////////////////////////////////////////////////////
// Function Copy_From_Eigen_Array
/////////////////////////////////////////////////////
template<typename T_STRUCT,typename T,int d> void 
EIGEN_INTERFACE<T_STRUCT,T,d>::Copy_From_Eigen_Array(const T_EIGEN_VECTOR& array_in,
                                                     T T_STRUCT::* u_field,const T_FLAG T_STRUCT::* flags_field, 
                                                     SPG_Allocator& allocator,/*const*/ SPG_Set_Type& set,
                                                     const std::vector<T_INDEX>& index_map_1d_to_nd){
    SPG_Data_Array_Type u=allocator.Get_Array(u_field);
    SPG_Const_Flags_Array_Type flags=allocator.Get_Const_Array(flags_field);
    #pragma omp parallel for
    for(int i=0;i<index_map_1d_to_nd.size();++i){
        unsigned long offset=T_MASK::Linear_Offset(std_array<int,d>(index_map_1d_to_nd[i]));
        if(set.IsPageActive(offset)&&(flags(offset)&SPGrid_Solver_Cell_Type_Interface)){
            u(offset)=array_in(i);}}
}
//#####################################################################
template class EIGEN_INTERFACE<SPGRID_DOMAIN_DECOMPOSITION_DATA<float>,float,2>;
template class EIGEN_INTERFACE<SPGRID_DOMAIN_DECOMPOSITION_DATA<float>,float,3>;

