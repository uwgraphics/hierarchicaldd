#ifndef __SPGRID_DOMAIN_DECOMPOSITION_MPI_H__
#define __SPGRID_DOMAIN_DECOMPOSITION_MPI_H__
#include "SPGrid_Gemini.h"

using namespace PhysBAM;
using namespace SPGrid;

template<typename T,typename T_STRUCT,int d,typename INDEX>
class SPGrid_Domain_Decomposition_MPI{
    typedef SPGrid_Allocator<T_STRUCT,d> SPG_Allocator;
    typedef typename SPG_Allocator::template Array<T>::type SPG_Data_Array_Type;
    typedef typename SPG_Allocator::template Array<const T>::type SPG_Const_Data_Array_Type;
    typedef typename SPG_Allocator::template Array<unsigned>::type SPG_Flags_Array_Type;
    typedef typename SPG_Allocator::template Array<const unsigned>::type SPG_Const_Flags_Array_Type;
    typedef typename SPG_Allocator::template Array<T>::mask T_MASK;
    typedef SPGrid_Set<SPG_Flags_Array_Type> SPG_Set_Type;
    typedef VECTOR<int,d> T_INDEX;
    typedef VECTOR<T,d> TV;
    typedef Eigen::Matrix<T,Eigen::Dynamic,1> T_EIGEN_VECTOR;
public:
    SPGrid_Gemini& spgrid_gemini;
    
}
#endif
