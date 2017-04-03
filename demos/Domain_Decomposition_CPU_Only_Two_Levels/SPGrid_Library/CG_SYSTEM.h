//#####################################################################
//Removed For Anonymity: Copyright Authors of Submission pap173s1
//Removed For Anonymity: Copyright Authors of Submission pap173s1
//#####################################################################
// Class CG_SYSTEM
//#####################################################################
#ifndef __CG_SYSTEM__
#define __CG_SYSTEM__

#include <PhysBAM_Tools/Krylov_Solvers/KRYLOV_SYSTEM_BASE.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Vectors/VECTOR_2D.h>
#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>
#include <PhysBAM_Tools/Vectors/VECTOR_ND.h>
#include <PhysBAM_Tools/Matrices/SPARSE_MATRIX_FLAT_NXN.h>
#include "CG_VECTOR.h"
#include "SPGrid_DD_Wrapper.h"

using namespace SPGrid;

namespace PhysBAM{

template<class T_STRUCT,class T,int d,typename INDEX>
class CG_SYSTEM:public KRYLOV_SYSTEM_BASE<T>
{
    typedef KRYLOV_SYSTEM_BASE<T> BASE;
    typedef KRYLOV_VECTOR_BASE<T> VECTOR_BASE;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<T>::type Data_array_type;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<const T>::type Const_data_array_type;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<const unsigned>::type Const_flag_array_type;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<unsigned>::type Flag_array_type;
    typedef SPGrid_Allocator<T_STRUCT,d> SPG_Allocator;
    typedef SPGrid_Set<Flag_array_type> SPG_Set_Type;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<const unsigned>::mask T_MASK;
    typedef VECTOR<int,d> T_INDEX;

    const int number_of_threads;
    SPGrid_DD_Wrapper<T,T_STRUCT,d,INDEX>& dd_wrapper;
    T T_STRUCT::* r_tmp_field;
    T T_STRUCT::* x_tmp_field;

//#####################################################################
public:
    CG_SYSTEM(SPGrid_DD_Wrapper<T,T_STRUCT,d,INDEX>& dd_wrapper_input,T T_STRUCT::* x_tmp_field,T T_STRUCT::* r_tmp_field, int number_of_threads = 0);
    void Multiply(const VECTOR_BASE& v,VECTOR_BASE& result) const;
    double Inner_Product(const VECTOR_BASE& x,const VECTOR_BASE& y) const;
    T Convergence_Norm(const VECTOR_BASE& x) const;
    void Project(VECTOR_BASE& x) const;
    void Set_Boundary_Conditions(VECTOR_BASE& x) const;
    void Project_Nullspace(VECTOR_BASE& x) const;
protected:
    void Apply_Preconditioner(const VECTOR_BASE& r, VECTOR_BASE& z) const;
//#####################################################################
};
}
#endif
