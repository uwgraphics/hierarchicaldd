//#####################################################################
//Removed For Anonymity: Copyright Authors of Submission pap173s1
// This file is part of PhysBAM whose distribution is governed by the license 
//Removed For Anonymity: Copyright Authors of Submission pap173s1
//#####################################################################
// Class EIGEN_INTERFACE_STENCIL & EIGEN_INTERFACE_SPGRID
//#####################################################################
#ifndef __EIGEN_INTERFACE_H__
#define __EIGEN_INTERFACE_H__
#include <Common_Tools/Math_Tools/RANGE_ITERATOR.h>
#include <PhysBAM_Tools/Data_Structures/HASHTABLE.h>
#include <PhysBAM_Tools/Data_Structures/HASHTABLE_ITERATOR.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Matrices/SPARSE_MATRIX_NXN.h>
#include <SPGrid/Core/SPGrid_Allocator.h>
#include <SPGrid/Core/SPGrid_Set.h>
#include <Eigen/Sparse>
#include <vector>

using namespace PhysBAM;
using namespace SPGrid;

template<typename T_STRUCT,typename T,int d>
class EIGEN_INTERFACE{
    typedef VECTOR<int,d> T_INDEX;
    typedef unsigned T_FLAG;
    typedef SPGrid_Allocator<T_STRUCT,d> SPG_Allocator;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<const T_FLAG>::type SPG_Const_Flags_Array_Type;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<T_FLAG>::type SPG_Flags_Array_Type;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<T>::type SPG_Data_Array_Type;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<const T>::type SPG_Const_Data_Array_Type;
    typedef SPGrid_Set<SPG_Flags_Array_Type> SPG_Set_Type;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<const T_FLAG>::mask T_MASK;
    typedef Eigen::Matrix<T,Eigen::Dynamic,1> T_EIGEN_VECTOR;
public:
    static void Copy_To_Eigen_Array(T_EIGEN_VECTOR& array_out,const T T_STRUCT::* u_field,const T_FLAG T_STRUCT::* flags_field, 
                                    /*const*/ SPG_Allocator& allocator,/*const*/ SPG_Set_Type& set,
                                    const std::vector<T_INDEX>& index_map_1d_to_nd);
    static void Copy_From_Eigen_Array(const T_EIGEN_VECTOR& array_in,T T_STRUCT::* u_field,const T_FLAG T_STRUCT::* flags_field, 
                                      SPG_Allocator& allocator,/*const*/ SPG_Set_Type& set,
                                      const std::vector<T_INDEX>& index_map_1d_to_nd);
};

#endif

