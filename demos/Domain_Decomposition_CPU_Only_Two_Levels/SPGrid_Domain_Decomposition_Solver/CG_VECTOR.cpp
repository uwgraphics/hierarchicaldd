//#####################################################################
//Removed For Anonymity: Copyright Authors of Submission pap173s1
//Removed For Anonymity: Copyright Authors of Submission pap173s1
//#####################################################################
// Class CG_VECTOR
//#####################################################################
#include "CG_VECTOR.h"
#include <Common_Tools/Math_Tools/RANGE_ITERATOR.h>
#include <SPGrid/Tools/SPGrid_Block_Iterator.h>
#include <SPGrid/Tools/SPGrid_Threading_Helper.h>
#include <Threading_Tools/PTHREAD_QUEUE.h>
#include "SPGrid_Gemini_Arithmetic_Helpers.h"
#include "SPGRID_DOMAIN_DECOMPOSITION_DATA.h"

using namespace PhysBAM;
using namespace SPGrid;
//#####################################################################
// operator+=
//#####################################################################
template<class T_STRUCT,class T,int d> KRYLOV_VECTOR_BASE<T>& CG_VECTOR<T_STRUCT,T,d>::
operator+=(const BASE& bv)
{
    if(number_of_threads)
        SPGrid_Computations::Threading_Helper<T_STRUCT,d>(Allocator(*this),Set(*this).Get_Blocks()).Run_Parallel(SPGrid_Computations::SPGrid_Plus_Equal<T_STRUCT,T,d>(Allocator(bv),Cg_Vector(bv).field,field),number_of_threads);
    else
        SPGrid_Computations::SPGrid_Plus_Equal<T_STRUCT,T,d>(Allocator(*this),Set(*this).Get_Blocks(),Allocator(bv),Cg_Vector(bv).field,field);
    return *this;
}
//#####################################################################
// operator-=
//#####################################################################
template<class T_STRUCT,class T,int d> KRYLOV_VECTOR_BASE<T>& CG_VECTOR<T_STRUCT,T,d>::
operator-=(const BASE& bv)
{
    if(number_of_threads)
        SPGrid_Computations::Threading_Helper<T_STRUCT,d>(Allocator(*this),Set(*this).Get_Blocks()).Run_Parallel(SPGrid_Computations::SPGrid_Minus_Equal<T_STRUCT,T,d>(Allocator(bv),Cg_Vector(bv).field,field),number_of_threads);
    else
        SPGrid_Computations::SPGrid_Minus_Equal<T_STRUCT,T,d>(Allocator(*this),Set(*this).Get_Blocks(),Allocator(bv),Cg_Vector(bv).field,field);
    return *this;
}
//#####################################################################
// operator*=
//#####################################################################
template<class T_STRUCT,class T,int d> KRYLOV_VECTOR_BASE<T>& CG_VECTOR<T_STRUCT,T,d>::
operator*=(const T a)
{
    if(number_of_threads)
        SPGrid_Computations::Threading_Helper<T_STRUCT,d>(Allocator(*this),Set(*this).Get_Blocks()).Run_Parallel(SPGrid_Computations::SPGrid_Scale<T_STRUCT,T,d>(field,a),number_of_threads);
    else
        SPGrid_Computations::SPGrid_Scale<T_STRUCT,T,d>(Allocator(*this),Set(*this).Get_Blocks(),field,a);
    return *this;
}
//#####################################################################
// Function Copy
//#####################################################################
template<class T_STRUCT,class T,int d> void CG_VECTOR<T_STRUCT,T,d>::
Copy(const T c,const BASE& bv)
{
    if(number_of_threads)
        SPGrid_Computations::Threading_Helper<T_STRUCT,d>(Allocator(*this),Set(*this).Get_Blocks()).Run_Parallel(SPGrid_Computations::SPGrid_Copy<T_STRUCT,T,d>(Allocator(bv),Cg_Vector(bv).field,field,c),number_of_threads);
    else
        SPGrid_Computations::SPGrid_Copy<T_STRUCT,T,d>(Allocator(*this),Set(*this).Get_Blocks(),Allocator(bv),Cg_Vector(bv).field,field,c);
}
//#####################################################################
// Function Copy
//#####################################################################
template<class T_STRUCT,class T,int d> void CG_VECTOR<T_STRUCT,T,d>::
Copy(const T c1,const BASE& bv1,const BASE& bv2)
{
    if(number_of_threads)
        SPGrid_Computations::Threading_Helper<T_STRUCT,d>(Allocator(*this),Set(*this).Get_Blocks()).Run_Parallel(SPGrid_Computations::SPGrid_Copy2<T_STRUCT,T,d>(Allocator(bv1),Allocator(bv2),Cg_Vector(bv1).field,Cg_Vector(bv2).field,field,c1),number_of_threads);
    else
        SPGrid_Computations::SPGrid_Copy2<T_STRUCT,T,d>(Allocator(*this),Set(*this).Get_Blocks(),Allocator(bv1),Allocator(bv2),Cg_Vector(bv1).field,Cg_Vector(bv2).field,field,c1);
}
//#####################################################################
// Function Raw_Size
//#####################################################################
template<class T_STRUCT,class T,int d> int CG_VECTOR<T_STRUCT,T,d>::
Raw_Size() const
{
    PHYSBAM_NOT_IMPLEMENTED();
}
//#####################################################################
// Function Raw_Get
//#####################################################################
template<class T_STRUCT,class T,int d> T& CG_VECTOR<T_STRUCT,T,d>::
Raw_Get(int i)
{
    PHYSBAM_NOT_IMPLEMENTED();
}
//#####################################################################
template class CG_VECTOR<SPGRID_DOMAIN_DECOMPOSITION_DATA<float>,float,2>;
template class CG_VECTOR<SPGRID_DOMAIN_DECOMPOSITION_DATA<float>,float,3>;

