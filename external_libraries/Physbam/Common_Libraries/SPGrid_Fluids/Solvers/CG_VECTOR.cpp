//#####################################################################
// Copyright 2011, Sean Bauer, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class CG_VECTOR
//#####################################################################
#include "CG_VECTOR.h"

#include <SPGrid_Fluids/Solvers/Blocked_Add_Helper.h>
#include <SPGrid_Fluids/Solvers/Blocked_Subtract_Helper.h>
#include <SPGrid_Fluids/Solvers/Blocked_Copy_Helper.h>
#include <SPGrid_Fluids/Simulation/FLUIDS_SIMULATION_DATA.h>

//#define TIMING

//#define SERIAL_ADD
//#define SERIAL_SUBTRACT
#define SERIAL_MULTIPLY
//#define SERIAL_COPY_C_V1
//#define SERIAL_COPY_C_V1_V2

using namespace PhysBAM;
using namespace SPGrid;

namespace PhysBAM{
extern int PhysBAM_number_of_threads;
}

//#####################################################################
// operator+=
//#####################################################################
template<class T_STRUCT, class T,int d> KRYLOV_VECTOR_BASE<T>& CG_VECTOR<T_STRUCT,T,d>::
operator+=(const BASE& bv)
{
#ifdef TIMING
    {
    LOG::SCOPE scope("CG_VECTOR::+=");
#endif

    const Hierarchy_type& bv_hierarchy=CG_VECTOR::Hierarchy(bv);
    T T_STRUCT::* const bv_field = Cg_Vector(bv).field;
    PHYSBAM_ASSERT(&hierarchy==&bv_hierarchy);

    for(int level=1; level<=hierarchy.Levels();level++)
    {
        Data_array_type d1 = hierarchy.Allocator(level).Get_Array(field);
        Const_data_array_type d2 = bv_hierarchy.Allocator(level).Get_Array(bv_field);
        Const_flag_array_type flags = hierarchy.Allocator(level).Get_Const_Array(&T_STRUCT::flags);

#ifdef SERIAL_ADD
        for(SPGrid_Block_Iterator<typename Data_array_type::MASK> iterator(hierarchy.Blocks(level));iterator.Valid();iterator.Next())
            if(iterator.Data(flags) & SPGrid_Cell_Type_Interior)
                iterator.Data(d1) += iterator.Data(d2);
#else
        Blocked_Add_Helper<T,Data_array_type::MASK::elements_per_block> helper(
            (T*)d1.Get_Data_Ptr(),
            (T*)d2.Get_Data_Ptr(),
            (unsigned*)flags.Get_Data_Ptr(),
            hierarchy.Blocks(level).first,
            hierarchy.Blocks(level).second);
        if(PhysBAM_number_of_threads)
            helper.Run_Parallel(PhysBAM_number_of_threads);
        else
            helper.Run();
#endif 
    }

#ifdef TIMING
    }
#endif

    return *this;
}
//#####################################################################
// operator-=
//#####################################################################
template<class T_STRUCT, class T,int d> KRYLOV_VECTOR_BASE<T>& CG_VECTOR<T_STRUCT,T,d>::
operator-=(const BASE& bv)
{
#ifdef TIMING
    {
    LOG::SCOPE scope("CG_VECTOR::-=");
#endif

    const Hierarchy_type& bv_hierarchy=CG_VECTOR::Hierarchy(bv);
    T T_STRUCT::* const bv_field = Cg_Vector(bv).field;
    PHYSBAM_ASSERT(&hierarchy==&bv_hierarchy);

    for(int level=1; level<=hierarchy.Levels();level++)
    {
        Data_array_type d1 = hierarchy.Allocator(level).Get_Array(field);
        Const_data_array_type d2 = bv_hierarchy.Allocator(level).Get_Array(bv_field);
        Const_flag_array_type flags = hierarchy.Allocator(level).Get_Const_Array(&T_STRUCT::flags);
         
#ifdef SERIAL_SUBTRACT
        for(SPGrid_Block_Iterator<typename Data_array_type::MASK> iterator(hierarchy.Blocks(level));iterator.Valid();iterator.Next())
            if(iterator.Data(flags) & SPGrid_Cell_Type_Interior)
                iterator.Data(d1) -= iterator.Data(d2);
#else
        Blocked_Subtract_Helper<T,Data_array_type::MASK::elements_per_block> helper(
            (T*)d1.Get_Data_Ptr(),
            (T*)d2.Get_Data_Ptr(),
            (unsigned*)flags.Get_Data_Ptr(),
            hierarchy.Blocks(level).first,
            hierarchy.Blocks(level).second);
        if(PhysBAM_number_of_threads)
            helper.Run_Parallel(PhysBAM_number_of_threads);
        else
            helper.Run();
#endif 
    }

#ifdef TIMING
    }
#endif

    return *this;
}
//#####################################################################
// operator*=
//#####################################################################
template<class T_STRUCT, class T,int d> KRYLOV_VECTOR_BASE<T>& CG_VECTOR<T_STRUCT,T,d>::
operator*=(const T a)
{
#ifdef TIMING
    {
    LOG::SCOPE scope("CG_VECTOR::*=");
#endif

    for(int level=1; level<=hierarchy.Levels();level++)
    {
        Data_array_type d1 = hierarchy.Allocator(level).Get_Array(field);
        Const_flag_array_type flags = hierarchy.Allocator(level).Get_Const_Array(&T_STRUCT::flags); 

#ifdef SERIAL_MULTIPLY
        for(SPGrid_Block_Iterator<typename Data_array_type::MASK> iterator(hierarchy.Blocks(level));iterator.Valid();iterator.Next())
            if(iterator.Data(flags) & SPGrid_Cell_Type_Interior)
                iterator.Data(d1)*=a;
#else
        Blocked_Multiply_Helper<T,Data_array_type::MASK::elements_per_block> helper(
            (T*)d1.Get_Data_Ptr(),
            a,
            (unsigned*)flags.Get_Data_Ptr(),
            hierarchy.Blocks(level).first,
            hierarchy.Blocks(level).second);
        helper.Run_Parallel(num_threads);
#endif 
    }

#ifdef TIMING
    }
#endif

    return *this;
}
//#####################################################################
// Function Copy
//#####################################################################
template<class T_STRUCT, class T,int d> void CG_VECTOR<T_STRUCT,T,d>::
Copy(const T c,const BASE& bv)
{
#ifdef TIMING
    {
    LOG::SCOPE scope("CG_VECTOR::Copy(c,bv)");
#endif

    const Hierarchy_type& bv_hierarchy=CG_VECTOR::Hierarchy(bv);
    T T_STRUCT::* const bv_field = Cg_Vector(bv).field;
    PHYSBAM_ASSERT(&hierarchy==&bv_hierarchy);

    for(int level=1; level<=hierarchy.Levels();level++)
    {
        Data_array_type d1 = hierarchy.Allocator(level).Get_Array(field);
        Const_data_array_type d2 = bv_hierarchy.Allocator(level).Get_Array(bv_field);
        Const_flag_array_type flags = hierarchy.Allocator(level).Get_Const_Array(&T_STRUCT::flags); 

#ifdef SERIAL_COPY_C_V1
        for(SPGrid_Block_Iterator<typename Data_array_type::MASK> iterator(hierarchy.Blocks(level));iterator.Valid();iterator.Next())
            if(iterator.Data(flags) & SPGrid_Cell_Type_Interior)
                iterator.Data(d1) = c * iterator.Data(d2);
#else
        Blocked_Copy2_Helper<T,Data_array_type::MASK::elements_per_block> helper(
            (T*)d1.Get_Data_Ptr(),
            (T*)d2.Get_Data_Ptr(),
            c,
            (unsigned*)flags.Get_Data_Ptr(),
            hierarchy.Blocks(level).first,
            hierarchy.Blocks(level).second);
        if(PhysBAM_number_of_threads)
            helper.Run_Parallel(PhysBAM_number_of_threads);
        else
            helper.Run();
#endif 
    }

#ifdef TIMING
    }
#endif
}
//#####################################################################
// Function Copy
//#####################################################################
template<class T_STRUCT, class T,int d> void CG_VECTOR<T_STRUCT,T,d>::
Copy(const T c1,const BASE& bv1,const BASE& bv2)
{
#ifdef TIMING
    {
    LOG::SCOPE scope("CG_VECTOR::Copy(c1,bv1,bv2)");
#endif

    const Hierarchy_type& bv1_hierarchy=CG_VECTOR::Hierarchy(bv1);
    const Hierarchy_type& bv2_hierarchy=CG_VECTOR::Hierarchy(bv2);
    T T_STRUCT::* const bv1_field = Cg_Vector(bv1).field;
    T T_STRUCT::* const bv2_field = Cg_Vector(bv2).field;
    PHYSBAM_ASSERT(&hierarchy==&bv1_hierarchy);
    PHYSBAM_ASSERT(&hierarchy==&bv2_hierarchy);

    for(int level=1; level<=hierarchy.Levels();level++)
    {
        Data_array_type d1 = hierarchy.Allocator(level).Get_Array(field);
        Const_data_array_type d2 = bv1_hierarchy.Allocator(level).Get_Array(bv1_field);
        Const_data_array_type d3 = bv2_hierarchy.Allocator(level).Get_Array(bv2_field);
        Const_flag_array_type flags = hierarchy.Allocator(level).Get_Const_Array(&T_STRUCT::flags); 

#ifdef SERIAL_COPY_C_V1_V2
        for(SPGrid_Block_Iterator<typename Data_array_type::MASK> iterator(hierarchy.Blocks(level));iterator.Valid();iterator.Next())
            if(iterator.Data(flags) & SPGrid_Cell_Type_Interior)
                iterator.Data(d1) = c1 * iterator.Data(d2) + iterator.Data(d3);
#else
        Blocked_Copy_Helper<T,Data_array_type::MASK::elements_per_block> helper(
            (T*)d1.Get_Data_Ptr(),
            (T*)d2.Get_Data_Ptr(),
            c1,
            (T*)d3.Get_Data_Ptr(),
            (unsigned*)flags.Get_Data_Ptr(),
            hierarchy.Blocks(level).first,
            hierarchy.Blocks(level).second);
        if(PhysBAM_number_of_threads)
            helper.Run_Parallel(PhysBAM_number_of_threads);
        else
            helper.Run();
#endif 
    }

#ifdef TIMING
    }
#endif
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
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class CG_VECTOR<FLUIDS_SIMULATION_DATA<double>,double,2>;
template class CG_VECTOR<FLUIDS_SIMULATION_DATA<double>,double,3>;
#endif

