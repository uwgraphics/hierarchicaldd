//#####################################################################
//Removed For Anonymity: Copyright Authors of Submission pap173s1
//Removed For Anonymity: Copyright Authors of Submission pap173s1
//#####################################################################
// Class CG_SYSTEM
//#####################################################################
#include "CG_SYSTEM.h"
#include "CG_VECTOR.h"
#include <SPGrid/Tools/SPGrid_Threading_Helper.h>
#include <SPGrid/Tools/SPGrid_Block_Iterator.h>
#include <Common_Tools/Math_Tools/RANGE_ITERATOR.h>
#include "SPGrid_Gemini_Arithmetic_Helpers.h"
#include "SPGRID_DOMAIN_DECOMPOSITION_DATA.h"

using namespace SPGrid;
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T_STRUCT,class T,int d,typename T_offset_ptr> CG_SYSTEM<T_STRUCT,T,d,T_offset_ptr>::
CG_SYSTEM(const SPGrid_Allocator<T_STRUCT,d>& allocator_flags_input,unsigned T_STRUCT::* flags_field_input,SPGrid_Domain_Decomposition_Preconditioner<T,T_STRUCT,d,T_offset_ptr>& dd_preconditioner_input,T T_STRUCT::* preconditioner_tmp_field_input,int number_of_threads_input)
    :allocator_flags(allocator_flags_input),flags_field(flags_field_input),dd_preconditioner(dd_preconditioner_input),preconditioner_tmp_field(preconditioner_tmp_field_input),number_of_threads(number_of_threads_input),BASE(false,false)
{}

//#####################################################################
// Function Multiply
//#####################################################################
template<class T_STRUCT,class T,int d,typename T_offset_ptr> void CG_SYSTEM<T_STRUCT,T,d,T_offset_ptr>::
Multiply(const VECTOR_BASE& v,VECTOR_BASE& result) const
{
    LOG::SCOPE scope("Multiply");
    T T_STRUCT::* const v_field = CG_VECTOR<T_STRUCT,T,d>::Cg_Vector(v).field;
    T T_STRUCT::* const result_field = CG_VECTOR<T_STRUCT,T,d>::Cg_Vector(result).field;

    const SPG_Allocator& v_allocator=CG_VECTOR<T_STRUCT,T,d>::Allocator(v);
    SPG_Allocator& result_allocator=CG_VECTOR<T_STRUCT,T,d>::Allocator(result);
    const SPG_Set_Type& v_set=CG_VECTOR<T_STRUCT,T,d>::Set(v);
    const SPG_Set_Type& result_set=CG_VECTOR<T_STRUCT,T,d>::Set(result);
    PHYSBAM_ASSERT(&v_set==&result_set);
    if(number_of_threads)
        SPGrid_Computations::Threading_Helper<T_STRUCT,d>(result_allocator,v_set.Get_Blocks()).Run_Parallel(SPGrid_Computations::SPGrid_Laplace<T_STRUCT,T,unsigned,d>(v_allocator,allocator_flags,v_field,result_field,flags_field,SPGrid_Solver_Cell_Type_Active|SPGrid_Solver_Cell_Type_Interface),number_of_threads);
    else
        SPGrid_Computations::SPGrid_Laplace<T_STRUCT,T,unsigned,d>(result_allocator,v_set.Get_Blocks(),v_allocator,allocator_flags,v_field,result_field,flags_field,SPGrid_Solver_Cell_Type_Active|SPGrid_Solver_Cell_Type_Interface);
}
//#####################################################################
// Function Inner_Product
//#####################################################################
template<class T_STRUCT,class T,int d,typename T_offset_ptr> 
double CG_SYSTEM<T_STRUCT,T,d,T_offset_ptr>::
Inner_Product(const VECTOR_BASE& v1,const VECTOR_BASE& v2) const
{
    LOG::SCOPE scope("Inner_Product");
    T T_STRUCT::* const v1_field = CG_VECTOR<T_STRUCT,T,d>::Cg_Vector(v1).field;
    T T_STRUCT::* const v2_field = CG_VECTOR<T_STRUCT,T,d>::Cg_Vector(v2).field;
    const SPG_Allocator& v1_allocator=CG_VECTOR<T_STRUCT,T,d>::Allocator(v1);
    const SPG_Allocator& v2_allocator=CG_VECTOR<T_STRUCT,T,d>::Allocator(v2);
    const SPG_Set_Type& v1_set=CG_VECTOR<T_STRUCT,T,d>::Set(v1);
    const SPG_Set_Type& v2_set=CG_VECTOR<T_STRUCT,T,d>::Set(v2);
    PHYSBAM_ASSERT(&v1_set == &v2_set);

    // Take dot-product of hierarchy, use doubles for temporaries
    double sum = 0;
    Const_Data_Array_Type d1 = v1_allocator.Get_Const_Array(v1_field);
    Const_Data_Array_Type d2 = v2_allocator.Get_Const_Array(v2_field);
    for(SPGrid_Block_Iterator<T_MASK> iterator(v1_set.Get_Blocks());iterator.Valid();iterator.Next_Block()){
        unsigned long data_offset=iterator.Offset();
        T_INDEX base_index=iterator.Index().template Cast<T_INDEX>();
        for(RANGE_ITERATOR<d> iterator(RANGE<T_INDEX>(base_index,base_index+v1_allocator.Block_Size().template Cast<T_INDEX>()-1));
            iterator.Valid();
            iterator.Next(),data_offset+=sizeof(T)){
            sum+=d1(data_offset)*d2(data_offset);}}
    return sum;
}
//#####################################################################
// Function Convergence_Norm
//#####################################################################
template<class T_STRUCT,class T,int d,typename T_offset_ptr> T CG_SYSTEM<T_STRUCT,T,d,T_offset_ptr>::
Convergence_Norm(const VECTOR_BASE& v) const
{
    LOG::SCOPE scope("Convergence_Norm");
    // Take maximum value of channel
    T T_STRUCT::* const v_field = CG_VECTOR<T_STRUCT,T,d>::Cg_Vector(v).field;
    const SPG_Allocator& v_allocator=CG_VECTOR<T_STRUCT,T,d>::Allocator(v);
    const SPG_Set_Type& v_set=CG_VECTOR<T_STRUCT,T,d>::Set(v);
    // Take dot-product of hierarchy, use doubles for temporaries
    T max = 0;
    Const_Data_Array_Type data = v_allocator.Get_Const_Array(v_field);
    for(SPGrid_Block_Iterator<T_MASK> iterator(v_set.Get_Blocks());iterator.Valid();iterator.Next_Block()){
        unsigned long data_offset=iterator.Offset();
        T_INDEX base_index=iterator.Index().template Cast<T_INDEX>();
        for(RANGE_ITERATOR<d> iterator(RANGE<T_INDEX>(base_index,base_index+v_allocator.Block_Size().template Cast<T_INDEX>()-1));
            iterator.Valid();
            iterator.Next(),data_offset+=sizeof(T)){
            max=(fabs(data(data_offset))>max)?fabs(data(data_offset)):max;}}
    return max;
}
//#####################################################################
// Function Project
//#####################################################################
template<class T_STRUCT,class T,int d,typename T_offset_ptr> void CG_SYSTEM<T_STRUCT,T,d,T_offset_ptr>::
Project(VECTOR_BASE& v) const
{
    LOG::SCOPE scope("Project");
    // Set all non-Interior nodes to zero.
    T T_STRUCT::* const v_field = CG_VECTOR<T_STRUCT,T,d>::Cg_Vector(v).field;
    SPG_Allocator& v_allocator=CG_VECTOR<T_STRUCT,T,d>::Allocator(v);
    const SPG_Set_Type& v_set=CG_VECTOR<T_STRUCT,T,d>::Set(v);
    Data_Array_Type data = v_allocator.Get_Array(v_field);
    Const_Flag_Array_Type flags = allocator_flags.Get_Const_Array(flags_field);
    for(SPGrid_Block_Iterator<T_MASK> iterator(v_set.Get_Blocks());iterator.Valid();iterator.Next_Block()){
        unsigned long data_offset=iterator.Offset();
        unsigned long flag_offset=iterator.Offset();
        T_INDEX base_index=iterator.Index().template Cast<T_INDEX>();
        for(RANGE_ITERATOR<d> iterator(RANGE<T_INDEX>(base_index,base_index+v_allocator.Block_Size().template Cast<T_INDEX>()-1));
            iterator.Valid();
            iterator.Next(),flag_offset+=sizeof(unsigned),data_offset+=sizeof(T)){
            if(!(flags(flag_offset)&(SPGrid_Solver_Cell_Type_Active|SPGrid_Solver_Cell_Type_Interface)))
                data(data_offset)=0.f;}}
}
//#####################################################################
// Function Set_Boundary_Conditions
//#####################################################################
template<class T_STRUCT,class T,int d,typename T_offset_ptr> void CG_SYSTEM<T_STRUCT,T,d,T_offset_ptr>::
Set_Boundary_Conditions(VECTOR_BASE& x) const
{
}
//#####################################################################
// Function Project_Nullspace
//#####################################################################
template<class T_STRUCT,class T,int d,typename T_offset_ptr> void CG_SYSTEM<T_STRUCT,T,d,T_offset_ptr>::
Project_Nullspace(VECTOR_BASE& x) const
{
}
//#####################################################################
// Function Apply_Preconditioner
//#####################################################################
template<class T_STRUCT,class T,int d,typename T_offset_ptr> void CG_SYSTEM<T_STRUCT,T,d,T_offset_ptr>::
Apply_Preconditioner(const VECTOR_BASE& r, VECTOR_BASE& z) const
{
    T T_STRUCT::* const r_field = CG_VECTOR<T_STRUCT,T,d>::Cg_Vector(r).field;
    T T_STRUCT::* const z_field = CG_VECTOR<T_STRUCT,T,d>::Cg_Vector(z).field;
    dd_preconditioner.Apply_Preconditioner(z_field,r_field,preconditioner_tmp_field);
    //z=r;
}

//#####################################################################
//template class CG_SYSTEM<SPGRID_DOMAIN_DECOMPOSITION_DATA<float>,float,2,unsigned int>;
template class CG_SYSTEM<SPGRID_DOMAIN_DECOMPOSITION_DATA<float>,float,3,unsigned int>;
template class CG_SYSTEM<SPGRID_DOMAIN_DECOMPOSITION_DATA<float>,float,3,unsigned long>;
