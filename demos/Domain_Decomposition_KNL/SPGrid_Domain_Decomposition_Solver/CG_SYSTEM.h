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
#include "SPGrid_Domain_Decomposition_Preconditioner.h"
#include "CG_VECTOR.h"

using namespace SPGrid;

namespace PhysBAM{

template<class T_STRUCT,class T,int d,typename T_offset_ptr>
class CG_SYSTEM:public KRYLOV_SYSTEM_BASE<T>
{
    typedef KRYLOV_SYSTEM_BASE<T> BASE;
    typedef KRYLOV_VECTOR_BASE<T> VECTOR_BASE;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<T>::type Data_Array_Type;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<const T>::type Const_Data_Array_Type;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<const unsigned>::type Const_Flag_Array_Type;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<unsigned>::type Flag_Array_Type;
    typedef SPGrid_Allocator<T_STRUCT,d> SPG_Allocator;
    typedef SPGrid_Set<Flag_Array_Type> SPG_Set_Type;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<const unsigned>::mask T_MASK;
    typedef VECTOR<int,d> T_INDEX;

    const int number_of_threads;
    unsigned T_STRUCT::* flags_field;
    T T_STRUCT::* preconditioner_tmp_field;
    const SPGrid_Allocator<T_STRUCT,d>& allocator_flags;
    SPGrid_Domain_Decomposition_Preconditioner<T,T_STRUCT,d,T_offset_ptr>& dd_preconditioner;
//#####################################################################
public:
    CG_SYSTEM(const SPGrid_Allocator<T_STRUCT,d>& allocator_flags,unsigned T_STRUCT::* flags_field,
              SPGrid_Domain_Decomposition_Preconditioner<T,T_STRUCT,d,T_offset_ptr>& dd_preconditioner,T T_STRUCT::* preconditioner_tmp_field,int number_of_threads=0);
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
