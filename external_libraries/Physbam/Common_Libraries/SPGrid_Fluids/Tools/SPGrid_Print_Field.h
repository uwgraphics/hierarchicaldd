//#####################################################################
// Scaleright 2013, Raj Setaluri, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Subroutine SPGrid_Print_Field::Print_Field
//#####################################################################
#ifndef __SPGrid_Print_Field_h__
#define __SPGrid_Print_Field_h__

#include <SPGrid_Fluids/Grids/GRID_HIERARCHY.h>
#include <SPGrid/Tools/SPGrid_Block_Iterator.h>

namespace SPGrid_Print_Field{

using namespace SPGrid;

template<class T_STRUCT,class T,int d> void
Print_Field(GRID_HIERARCHY<T_STRUCT,T,d>& hierarchy,T T_STRUCT::* data_field,std::ostream& output)
{
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<const T>::type Const_data_array_type;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<T>::type Data_array_type;
    
    for(int level=1;level<=hierarchy.Levels();level++){
        Const_data_array_type data=hierarchy.Allocator(level).Get_Const_Array(data_field);
        for(SPGrid_Block_Iterator<typename Data_array_type::MASK> iterator(hierarchy.Blocks(level));iterator.Valid();iterator.Next())
            output<<"level="<<level<<", index="<<iterator.Index()<<", value="<<iterator.Data(data)<<std::endl;}
}

template<class T_STRUCT,class T,int d> void
Masked_Print_Field(GRID_HIERARCHY<T_STRUCT,T,d>& hierarchy,T T_STRUCT::* data_field,unsigned T_STRUCT::* flags_field,std::ostream& output,const unsigned mask)
{
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<const T>::type Const_data_array_type;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<const unsigned>::type Const_flag_array_type;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<T>::type Data_array_type;

    for(int level=1;level<=hierarchy.Levels();level++){
        Const_data_array_type data=hierarchy.Allocator(level).Get_Const_Array(data_field);
        Const_flag_array_type flags=hierarchy.Allocator(level).Get_Const_Array(flags_field);
        for(SPGrid_Block_Iterator<typename Data_array_type::MASK> iterator(hierarchy.Blocks(level));iterator.Valid();iterator.Next())
            if(iterator.Data(flags)&mask)
                output<<"level="<<level<<", index="<<iterator.Index()<<", value="<<iterator.Data(data)<<std::endl;}
}
//#####################################################################
}
#endif
