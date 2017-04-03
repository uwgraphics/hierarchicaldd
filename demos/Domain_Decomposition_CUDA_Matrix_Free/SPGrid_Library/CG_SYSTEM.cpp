//#####################################################################
//Removed For Anonymity: Copyright Authors of Submission pap173s1
//Removed For Anonymity: Copyright Authors of Submission pap173s1
//#####################################################################
// Class CG_SYSTEM
//#####################################################################
#include "SPGrid_Laplace_Legacy.h"
#include "CG_SYSTEM.h"
#include "CG_VECTOR.h"
#include <iomanip>
#include <SPGrid/Tools/SPGrid_Threading_Helper.h>
#include <SPGrid/Tools/SPGrid_Block_Iterator.h>
#include <Common_Tools/Math_Tools/RANGE_ITERATOR.h>
#include <SPGrid_Fluids/Solvers/Domain_Decomposition/SPGRID_MULTIGRID_FLAGS.h>
#include "SPGrid_V_Cycle_Helper.h"

using namespace SPGrid;
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T_STRUCT,class T,int d,typename INDEX> CG_SYSTEM<T_STRUCT,T,d,INDEX>::
CG_SYSTEM(SPGrid_DD_Wrapper<T,T_STRUCT,d,INDEX>& dd_wrapper_input,T T_STRUCT::* x_tmp_field_input,T T_STRUCT::* r_tmp_field_input,int number_of_threads_input)
    :dd_wrapper(dd_wrapper_input),x_tmp_field(x_tmp_field_input),r_tmp_field(r_tmp_field_input),BASE(false,false),number_of_threads(number_of_threads_input){}


//#####################################################################
// Function Multiply
//#####################################################################
template<class T_STRUCT,class T,int d,typename INDEX> void CG_SYSTEM<T_STRUCT,T,d,INDEX>::
Multiply(const VECTOR_BASE& v,VECTOR_BASE& result) const
{
    T T_STRUCT::* const v_field = CG_VECTOR<T_STRUCT,T,d>::Cg_Vector(v).field;
    T T_STRUCT::* const result_field = CG_VECTOR<T_STRUCT,T,d>::Cg_Vector(result).field;

    //dd_wrapper.Multiply(v_field,result_field);
    //return;

    const SPG_Allocator& v_allocator=CG_VECTOR<T_STRUCT,T,d>::Allocator(v);
    SPG_Allocator& result_allocator=CG_VECTOR<T_STRUCT,T,d>::Allocator(result);
    const SPG_Set_Type& v_set=CG_VECTOR<T_STRUCT,T,d>::Set(v);
    const SPG_Set_Type& result_set=CG_VECTOR<T_STRUCT,T,d>::Set(result);
    PHYSBAM_ASSERT(&v_allocator == &result_allocator);
    PHYSBAM_ASSERT(&v_set == &result_set);
    if(number_of_threads)
        SPGrid_Computations::Threading_Helper<T_STRUCT,d>(result_allocator,v_set.Get_Blocks()).Run_Parallel(SPGrid_Computations::Laplace<T_STRUCT,T,unsigned,d>(v_field,result_field,&T_STRUCT::flags),number_of_threads);
    else
        SPGrid_Computations::Laplace<T_STRUCT,T,unsigned,d>(result_allocator,v_set.Get_Blocks(),v_field,result_field,&T_STRUCT::flags);
}
//#####################################################################
// Function Inner_Product
//#####################################################################
template<class T_STRUCT,class T,int d,typename INDEX> 
double CG_SYSTEM<T_STRUCT,T,d,INDEX>::
Inner_Product(const VECTOR_BASE& v1,const VECTOR_BASE& v2) const
{
    T T_STRUCT::* const v1_field = CG_VECTOR<T_STRUCT,T,d>::Cg_Vector(v1).field;
    T T_STRUCT::* const v2_field = CG_VECTOR<T_STRUCT,T,d>::Cg_Vector(v2).field;
    const SPG_Allocator& v1_allocator=CG_VECTOR<T_STRUCT,T,d>::Allocator(v1);
    const SPG_Allocator& v2_allocator=CG_VECTOR<T_STRUCT,T,d>::Allocator(v2);
    const SPG_Set_Type& v1_set=CG_VECTOR<T_STRUCT,T,d>::Set(v1);
    const SPG_Set_Type& v2_set=CG_VECTOR<T_STRUCT,T,d>::Set(v2);
    PHYSBAM_ASSERT(&v1_allocator == &v2_allocator);
    PHYSBAM_ASSERT(&v1_set == &v2_set);

    // Take dot-product of hierarchy, use doubles for temporaries
    double sum = 0;
    Const_data_array_type d1 = v1_allocator.Get_Const_Array(v1_field);
    Const_data_array_type d2 = v1_allocator.Get_Const_Array(v2_field);
    for(SPGrid_Block_Iterator<T_MASK> iterator(v1_set.Get_Blocks());iterator.Valid();iterator.Next_Block()){
        unsigned long data_offset=iterator.Offset();
        T_INDEX base_index=iterator.Index().template Cast<T_INDEX>();
        for(RANGE_ITERATOR<d> iterator(RANGE<T_INDEX>(base_index,base_index+v1_allocator.Block_Size().template Cast<T_INDEX>()-1));
            iterator.Valid();
            iterator.Next(),data_offset+=sizeof(T)){
            sum += d1(data_offset) * d2(data_offset); 
        }
    }
    return sum;
}
//#####################################################################
// Function Convergence_Norm
//#####################################################################
template<class T_STRUCT,class T,int d,typename INDEX> T CG_SYSTEM<T_STRUCT,T,d,INDEX>::
Convergence_Norm(const VECTOR_BASE& v) const
{
    // Take maximum value of channel
    T T_STRUCT::* const v_field = CG_VECTOR<T_STRUCT,T,d>::Cg_Vector(v).field;
    const SPG_Allocator& v_allocator=CG_VECTOR<T_STRUCT,T,d>::Allocator(v);
    const SPG_Set_Type& v_set=CG_VECTOR<T_STRUCT,T,d>::Set(v);
    // Take dot-product of hierarchy, use doubles for temporaries
    T max = 0;
    Const_data_array_type data = v_allocator.Get_Const_Array(v_field);
    for(SPGrid_Block_Iterator<T_MASK> iterator(v_set.Get_Blocks());iterator.Valid();iterator.Next_Block()){
        unsigned long data_offset=iterator.Offset();
        T_INDEX base_index=iterator.Index().template Cast<T_INDEX>();
        for(RANGE_ITERATOR<d> iterator(RANGE<T_INDEX>(base_index,base_index+v_allocator.Block_Size().template Cast<T_INDEX>()-1));
            iterator.Valid();
            iterator.Next(),data_offset+=sizeof(T)){
            max = (fabs(data(data_offset)) > max) ? fabs(data(data_offset)) : max; 
        }
    }
    return max;
}
//#####################################################################
// Function Project
//#####################################################################
template<class T_STRUCT,class T,int d,typename INDEX> void CG_SYSTEM<T_STRUCT,T,d,INDEX>::
Project(VECTOR_BASE& v) const
{
    // Set all non-Interior nodes to zero.
    T T_STRUCT::* const v_field = CG_VECTOR<T_STRUCT,T,d>::Cg_Vector(v).field;
    SPG_Allocator& v_allocator=CG_VECTOR<T_STRUCT,T,d>::Allocator(v);
    const SPG_Set_Type& v_set=CG_VECTOR<T_STRUCT,T,d>::Set(v);
    Data_array_type data = v_allocator.Get_Array(v_field);
    Const_flag_array_type flags = v_allocator.Get_Const_Array(&T_STRUCT::flags);
    for(SPGrid_Block_Iterator<T_MASK> iterator(v_set.Get_Blocks());iterator.Valid();iterator.Next_Block()){
        unsigned long data_offset=iterator.Offset();
        unsigned long flag_offset=iterator.Offset();
        T_INDEX base_index=iterator.Index().template Cast<T_INDEX>();
        for(RANGE_ITERATOR<d> iterator(RANGE<T_INDEX>(base_index,base_index+v_allocator.Block_Size().template Cast<T_INDEX>()-1));
            iterator.Valid();
            iterator.Next(),flag_offset+=sizeof(unsigned),data_offset+=sizeof(T)){
            if(!(flags(flag_offset) & (SPGrid_Solver_Cell_Type_Active | SPGrid_Solver_Cell_Type_Interface)))
                data(data_offset) = 0.f;
        }
    }
}
//#####################################################################
// Function Set_Boundary_Conditions
//#####################################################################
template<class T_STRUCT,class T,int d,typename INDEX> void CG_SYSTEM<T_STRUCT,T,d,INDEX>::
Set_Boundary_Conditions(VECTOR_BASE& x) const
{
}
//#####################################################################
// Function Project_Nullspace
//#####################################################################
template<class T_STRUCT,class T,int d,typename INDEX> void CG_SYSTEM<T_STRUCT,T,d,INDEX>::
Project_Nullspace(VECTOR_BASE& x) const
{
}
//#####################################################################
// Function Apply_Preconditioner
//#####################################################################
template<class T_STRUCT,class T,int d,typename INDEX> void CG_SYSTEM<T_STRUCT,T,d,INDEX>::
Apply_Preconditioner(const VECTOR_BASE& r, VECTOR_BASE& z) const
{
    //z=r;
    dd_wrapper.Apply_DD(CG_VECTOR<T_STRUCT,T,d>::Cg_Vector(z).field,CG_VECTOR<T_STRUCT,T,d>::Cg_Vector(r).field,number_of_threads);
        
        //SPGrid_Computations::Clear<T_STRUCT,T,d>(dd_wrapper.allocator,dd_wrapper.set.Get_Blocks(),CG_VECTOR<T_STRUCT,T,d>::Cg_Vector(z).field);
        //SPGrid_Computations::Clear<T_STRUCT,T,d>(dd_wrapper.allocator,dd_wrapper.set.Get_Blocks(),x_tmp_field);
    //const int smoothing_itrs = 0;
    // for(int i = 0;i < smoothing_itrs;++i)
    //     SPGrid_Computations::V_Cycle_Helper<T_STRUCT,T,d>::global_smoothing(CG_VECTOR<T_STRUCT,T,d>::Cg_Vector(z).field,r_tmp_field,CG_VECTOR<T_STRUCT,T,d>::Cg_Vector(r).field,dd_wrapper.flags_field,dd_wrapper.allocator,dd_wrapper.set,number_of_threads);
    

    //CG_VECTOR<T_STRUCT,T,d> b_tmp(dd_wrapper.allocator,dd_wrapper.set,dd_wrapper.b_v_field);
    //b_tmp.Copy(T(1),r);
    //SPGrid_Computations::V_Cycle_Helper<T_STRUCT,T,d>::v_cycle(CG_VECTOR<T_STRUCT,T,d>::Cg_Vector(z).field,dd_wrapper.r_v_field,dd_wrapper.b_v_field,dd_wrapper.flags_field,dd_wrapper.v_cycle_topology,number_of_threads);

    /*SPGrid_Computations::V_Cycle_Helper<T_STRUCT,T,d>::compute_residual_global(CG_VECTOR<T_STRUCT,T,d>::Cg_Vector(z).field,r_tmp_field,CG_VECTOR<T_STRUCT,T,d>::Cg_Vector(r).field,dd_wrapper.flags_field,dd_wrapper.allocator,dd_wrapper.set,number_of_threads);
    dd_wrapper.Apply_DD(x_tmp_field,r_tmp_field,number_of_threads);
    CG_VECTOR<T_STRUCT,T,d> x_tmp(dd_wrapper.allocator,dd_wrapper.set,x_tmp_field);
    z += x_tmp;*/

    //for(int i = 0;i < smoothing_itrs;++i)
    //SPGrid_Computations::V_Cycle_Helper<T_STRUCT,T,d>::global_smoothing(CG_VECTOR<T_STRUCT,T,d>::Cg_Vector(z).field,r_tmp_field,CG_VECTOR<T_STRUCT,T,d>::Cg_Vector(r).field,dd_wrapper.flags_field,dd_wrapper.allocator,dd_wrapper.set,number_of_threads);
}

//#####################################################################
template class CG_SYSTEM<FLUIDS_SIMULATION_DATA<float>,float,2,int>;
template class CG_SYSTEM<FLUIDS_SIMULATION_DATA<float>,float,3,int>;
