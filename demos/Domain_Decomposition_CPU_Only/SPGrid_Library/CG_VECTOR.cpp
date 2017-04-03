//#####################################################################
//Removed For Anonymity: Copyright Authors of Submission pap173s1
//Removed For Anonymity: Copyright Authors of Submission pap173s1
//#####################################################################
// Class CG_VECTOR
//#####################################################################
#include "CG_VECTOR.h"
#include <SPGrid/Tools/SPGrid_Block_Iterator.h>
#include <Common_Tools/Math_Tools/RANGE_ITERATOR.h>

using namespace PhysBAM;
using namespace SPGrid;
//#####################################################################
// operator+=
//#####################################################################
template<class T_STRUCT, class T,int d> KRYLOV_VECTOR_BASE<T>& CG_VECTOR<T_STRUCT,T,d>::
operator+=(const BASE& bv)
{
    T T_STRUCT::* const bv_field = Cg_Vector(bv).field;
    PHYSBAM_ASSERT(&allocator==&Allocator(bv));
    PHYSBAM_ASSERT(&set==&Set(bv));
    
    Data_array_type d1 = allocator.Get_Array(field);
    Const_data_array_type d2 = Allocator(bv).Get_Array(bv_field);
    Const_flag_array_type flags = Allocator(bv).Get_Const_Array(&T_STRUCT::flags);
    for(SPGrid_Block_Iterator<T_MASK> iterator(set.Get_Blocks());iterator.Valid();iterator.Next_Block()){
        unsigned long flag_offset=iterator.Offset();
        unsigned long data_offset=iterator.Offset();
        T_INDEX base_index=iterator.Index().template Cast<T_INDEX>();
        for(RANGE_ITERATOR<d> iterator(RANGE<T_INDEX>(base_index,base_index+allocator.Block_Size().template Cast<T_INDEX>()-1));
            iterator.Valid();
            iterator.Next(),flag_offset+=sizeof(unsigned),data_offset+=sizeof(T)){
            if(flags(flag_offset) & (SPGrid_Solver_Cell_Type_Active | SPGrid_Solver_Cell_Type_Interface))
                d1(data_offset) += d2(data_offset);
        }
    }
}
//#####################################################################
// operator-=
//#####################################################################
template<class T_STRUCT, class T,int d> KRYLOV_VECTOR_BASE<T>& CG_VECTOR<T_STRUCT,T,d>::
operator-=(const BASE& bv)
{
    T T_STRUCT::* const bv_field = Cg_Vector(bv).field;
    PHYSBAM_ASSERT(&allocator==&Allocator(bv));
    PHYSBAM_ASSERT(&set==&Set(bv));

    Data_array_type d1 = allocator.Get_Array(field);
    Const_data_array_type d2 = Allocator(bv).Get_Array(bv_field);
    Const_flag_array_type flags = Allocator(bv).Get_Const_Array(&T_STRUCT::flags);
    for(SPGrid_Block_Iterator<T_MASK> iterator(set.Get_Blocks());iterator.Valid();iterator.Next_Block()){
        unsigned long flag_offset=iterator.Offset();
        unsigned long data_offset=iterator.Offset();
        T_INDEX base_index=iterator.Index().template Cast<T_INDEX>();
        for(RANGE_ITERATOR<d> iterator(RANGE<T_INDEX>(base_index,base_index+allocator.Block_Size().template Cast<T_INDEX>()-1));
            iterator.Valid();
            iterator.Next(),flag_offset+=sizeof(unsigned),data_offset+=sizeof(T)){
            if(flags(flag_offset) & (SPGrid_Solver_Cell_Type_Active | SPGrid_Solver_Cell_Type_Interface))
                d1(data_offset) -= d2(data_offset);
        }
    }
}
//#####################################################################
// operator*=
//#####################################################################
template<class T_STRUCT, class T,int d> KRYLOV_VECTOR_BASE<T>& CG_VECTOR<T_STRUCT,T,d>::
operator*=(const T a)
{
    Data_array_type d1 = allocator.Get_Array(field);
    Const_flag_array_type flags = allocator.Get_Const_Array(&T_STRUCT::flags);
    for(SPGrid_Block_Iterator<T_MASK> iterator(set.Get_Blocks());iterator.Valid();iterator.Next_Block()){
        unsigned long flag_offset=iterator.Offset();
        unsigned long data_offset=iterator.Offset();
        T_INDEX base_index=iterator.Index().template Cast<T_INDEX>();
        for(RANGE_ITERATOR<d> iterator(RANGE<T_INDEX>(base_index,base_index+allocator.Block_Size().template Cast<T_INDEX>()-1));
            iterator.Valid();
            iterator.Next(),flag_offset+=sizeof(unsigned),data_offset+=sizeof(T)){
            if(flags(flag_offset) & (SPGrid_Solver_Cell_Type_Active | SPGrid_Solver_Cell_Type_Interface))
                d1(data_offset) *= a;
        }
    }
    return *this;
}
//#####################################################################
// Function Copy
//#####################################################################
template<class T_STRUCT, class T,int d> void CG_VECTOR<T_STRUCT,T,d>::
Copy(const T c,const BASE& bv)
{
    T T_STRUCT::* const bv_field = Cg_Vector(bv).field;
    PHYSBAM_ASSERT(&allocator==&Allocator(bv));
    PHYSBAM_ASSERT(&set==&Set(bv));

    Data_array_type d1 = allocator.Get_Array(field);
    Const_data_array_type d2 = Allocator(bv).Get_Array(bv_field);
    Const_flag_array_type flags = Allocator(bv).Get_Const_Array(&T_STRUCT::flags);
    for(SPGrid_Block_Iterator<T_MASK> iterator(set.Get_Blocks());iterator.Valid();iterator.Next_Block()){
        unsigned long flag_offset=iterator.Offset();
        unsigned long data_offset=iterator.Offset();
        T_INDEX base_index=iterator.Index().template Cast<T_INDEX>();
        for(RANGE_ITERATOR<d> iterator(RANGE<T_INDEX>(base_index,base_index+allocator.Block_Size().template Cast<T_INDEX>()-1));
            iterator.Valid();
            iterator.Next(),flag_offset+=sizeof(unsigned),data_offset+=sizeof(T)){
            if(flags(flag_offset) & (SPGrid_Solver_Cell_Type_Active | SPGrid_Solver_Cell_Type_Interface))
                d1(data_offset) = c * d2(data_offset);
        }
    }
}
//#####################################################################
// Function Copy
//#####################################################################
template<class T_STRUCT, class T,int d> void CG_VECTOR<T_STRUCT,T,d>::
Copy(const T c1,const BASE& bv1,const BASE& bv2)
{
    T T_STRUCT::* const bv1_field = Cg_Vector(bv1).field;
    T T_STRUCT::* const bv2_field = Cg_Vector(bv2).field;
    PHYSBAM_ASSERT(&allocator==&Allocator(bv1));
    PHYSBAM_ASSERT(&set==&Set(bv1));
    PHYSBAM_ASSERT(&allocator==&Allocator(bv2));
    PHYSBAM_ASSERT(&set==&Set(bv2));

    Data_array_type d1 = allocator.Get_Array(field);
    Const_data_array_type d2 = Allocator(bv1).Get_Array(bv1_field);
    Const_data_array_type d3 = Allocator(bv2).Get_Array(bv2_field);
    Const_flag_array_type flags = allocator.Get_Const_Array(&T_STRUCT::flags);
    for(SPGrid_Block_Iterator<T_MASK> iterator(set.Get_Blocks());iterator.Valid();iterator.Next_Block()){
        unsigned long flag_offset=iterator.Offset();
        unsigned long data_offset=iterator.Offset();
        T_INDEX base_index=iterator.Index().template Cast<T_INDEX>();
        for(RANGE_ITERATOR<d> iterator(RANGE<T_INDEX>(base_index,base_index+allocator.Block_Size().template Cast<T_INDEX>()-1));
            iterator.Valid();
            iterator.Next(),flag_offset+=sizeof(unsigned),data_offset+=sizeof(T)){
            if(flags(flag_offset) & (SPGrid_Solver_Cell_Type_Active | SPGrid_Solver_Cell_Type_Interface))
                d1(data_offset) = c1 * d2(data_offset) + d3(data_offset); 
        }
    }
}
//#####################################################################
// Function Raw_Size
//#####################################################################
template<class T_STRUCT, class T,int d> int CG_VECTOR<T_STRUCT,T,d>::
Raw_Size() const
{
    PHYSBAM_NOT_IMPLEMENTED();
}
//#####################################################################
// Function Raw_Get
//#####################################################################
template<class T_STRUCT, class T,int d> T& CG_VECTOR<T_STRUCT,T,d>::
Raw_Get(int i)
{
    PHYSBAM_NOT_IMPLEMENTED();
}
//#####################################################################
template class CG_VECTOR<FLUIDS_SIMULATION_DATA<float>,float,2>;
template class CG_VECTOR<FLUIDS_SIMULATION_DATA<float>,float,3>;

