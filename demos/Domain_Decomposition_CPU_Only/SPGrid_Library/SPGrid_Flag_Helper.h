//#####################################################################
//Removed For Anonymity: Copyright Authors of Submission pap173s1
//Removed For Anonymity: Copyright Authors of Submission pap173s1
//#####################################################################
// Subroutine SPGrid_Computations::Flaging
//#####################################################################
#ifndef __SPGrid_Flaging_h__
#define __SPGrid_Flaging_h__

#include <Common_Tools/Math_Tools/RANGE_ITERATOR.h>
#include <SPGrid/Core/SPGrid_Allocator.h>
#include "SPGRID_MULTIGRID_FLAGS.h"

namespace SPGrid_Computations{

using namespace SPGrid;
using namespace PhysBAM;

template<class T_STRUCT,int d>
    class Flaging;

template<class T_STRUCT>
    class Flaging<T_STRUCT,3>
{
    enum{d=3};
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<unsigned>::type Flag_Array_Type;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<unsigned>::mask T_MASK;
    typedef VECTOR<int,d> T_INDEX;

    unsigned T_STRUCT::* flags_field;
    SPGrid_Set<Flag_Array_Type>& set;
public:
     Flaging(unsigned T_STRUCT::* flags_field_input,SPGrid_Set<Flag_Array_Type>& set_input)
         :flags_field(flags_field_input),set(set_input)
    {}
    
     Flaging(SPGrid_Allocator<T_STRUCT,d>& allocator,const std::pair<const unsigned long*,unsigned>& blocks,unsigned T_STRUCT::* flags_field_input,SPGrid_Set<Flag_Array_Type>& set_input)
         :flags_field(flags_field_input),set(set_input)
    {Run(allocator,blocks);}
    
    void Run(SPGrid_Allocator<T_STRUCT,d>& allocator,const std::pair<const unsigned long*,unsigned>& blocks) const
    {
        Flag_Array_Type flags=allocator.Get_Array(flags_field);    
        enum{
            block_xsize = 1u << T_MASK::block_xbits,
            block_ysize = 1u << T_MASK::block_ybits,
            block_zsize = 1u << T_MASK::block_zbits  
        };
        typedef unsigned (&block_flag)[block_xsize][block_ysize][block_zsize];
        for(SPGrid_Block_Iterator<T_MASK> iterator(blocks);iterator.Valid();iterator.Next_Block()){
            unsigned long offset=iterator.Offset();
            T_INDEX base_index=iterator.Index().template Cast<T_INDEX>();
            for(RANGE_ITERATOR<d> iterator(RANGE<T_INDEX>(base_index,base_index+allocator.Block_Size().template Cast<T_INDEX>()-1));
                iterator.Valid();
                iterator.Next(),offset+=sizeof(unsigned)){
                unsigned& flag = flags(offset);
                if(flag & (SPGrid_Solver_Cell_Type_Active | SPGrid_Solver_Cell_Type_Interface)){
                    if(set.Is_Set(T_MASK::Packed_OffsetXdim<-1>(offset),(SPGrid_Solver_Cell_Type_Active | SPGrid_Solver_Cell_Type_Interface | SPGrid_Solver_Cell_Type_Dirichlet)))
                        flag |= SPGrid_Solver_Face_Minus_X_Active;
                    if(set.Is_Set(T_MASK::Packed_OffsetXdim<+1>(offset),(SPGrid_Solver_Cell_Type_Active | SPGrid_Solver_Cell_Type_Interface | SPGrid_Solver_Cell_Type_Dirichlet)))
                        flag |= SPGrid_Solver_Face_Plus_X_Active;
                    if(set.Is_Set(T_MASK::Packed_OffsetYdim<-1>(offset),(SPGrid_Solver_Cell_Type_Active | SPGrid_Solver_Cell_Type_Interface | SPGrid_Solver_Cell_Type_Dirichlet)))
                        flag |= SPGrid_Solver_Face_Minus_Y_Active;
                    if(set.Is_Set(T_MASK::Packed_OffsetYdim<+1>(offset),(SPGrid_Solver_Cell_Type_Active | SPGrid_Solver_Cell_Type_Interface | SPGrid_Solver_Cell_Type_Dirichlet)))
                        flag |= SPGrid_Solver_Face_Plus_Y_Active;
                    if(set.Is_Set(T_MASK::Packed_OffsetZdim<-1>(offset),(SPGrid_Solver_Cell_Type_Active | SPGrid_Solver_Cell_Type_Interface | SPGrid_Solver_Cell_Type_Dirichlet)))
                        flag |= SPGrid_Solver_Face_Minus_Z_Active;
                    if(set.Is_Set(T_MASK::Packed_OffsetZdim<+1>(offset),(SPGrid_Solver_Cell_Type_Active | SPGrid_Solver_Cell_Type_Interface | SPGrid_Solver_Cell_Type_Dirichlet)))
                        flag |= SPGrid_Solver_Face_Plus_Z_Active;
                }
            }
        }
    } 
};
///////////////////////////////////////////////////////////////////////////////////////
// 2D
///////////////////////////////////////////////////////////////////////////////////////
template<class T_STRUCT>
    class Flaging<T_STRUCT,2>
{
    enum{d=2};
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<unsigned>::type Flag_Array_Type;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<unsigned>::mask T_MASK;
    typedef VECTOR<int,d> T_INDEX;

    unsigned T_STRUCT::* flags_field;
    SPGrid_Set<Flag_Array_Type>& set;
public:
     Flaging(unsigned T_STRUCT::* flags_field_input,SPGrid_Set<Flag_Array_Type>& set_input)
         :flags_field(flags_field_input),set(set_input)
    {}
    
     Flaging(SPGrid_Allocator<T_STRUCT,d>& allocator,const std::pair<const unsigned long*,unsigned>& blocks,unsigned T_STRUCT::* flags_field_input,SPGrid_Set<Flag_Array_Type>& set_input)
         :flags_field(flags_field_input),set(set_input)
    {Run(allocator,blocks);}
    
    void Run(SPGrid_Allocator<T_STRUCT,d>& allocator,const std::pair<const unsigned long*,unsigned>& blocks) const
    {
        Flag_Array_Type flags=allocator.Get_Array(flags_field);    
        enum{
            block_xsize = 1u << T_MASK::block_xbits,
            block_ysize = 1u << T_MASK::block_ybits
        };
        typedef unsigned (&block_flag)[block_xsize][block_ysize];
        for(SPGrid_Block_Iterator<T_MASK> iterator(blocks);iterator.Valid();iterator.Next_Block()){
            unsigned long offset=iterator.Offset();
            T_INDEX base_index=iterator.Index().template Cast<T_INDEX>();
            for(RANGE_ITERATOR<d> iterator(RANGE<T_INDEX>(base_index,base_index+allocator.Block_Size().template Cast<T_INDEX>()-1));
                iterator.Valid();
                iterator.Next(),offset+=sizeof(unsigned)){
                unsigned& flag = flags(offset);
                if(flag & (SPGrid_Solver_Cell_Type_Active | SPGrid_Solver_Cell_Type_Interface)){
                    if(set.Is_Set(T_MASK::Packed_OffsetXdim<-1>(offset),(SPGrid_Solver_Cell_Type_Active | SPGrid_Solver_Cell_Type_Interface | SPGrid_Solver_Cell_Type_Dirichlet)))
                        flag |= SPGrid_Solver_Face_Minus_X_Active;
                    if(set.Is_Set(T_MASK::Packed_OffsetXdim<+1>(offset),(SPGrid_Solver_Cell_Type_Active | SPGrid_Solver_Cell_Type_Interface | SPGrid_Solver_Cell_Type_Dirichlet)))
                        flag |= SPGrid_Solver_Face_Plus_X_Active;
                    if(set.Is_Set(T_MASK::Packed_OffsetYdim<-1>(offset),(SPGrid_Solver_Cell_Type_Active | SPGrid_Solver_Cell_Type_Interface | SPGrid_Solver_Cell_Type_Dirichlet)))
                        flag |= SPGrid_Solver_Face_Minus_Y_Active;
                    if(set.Is_Set(T_MASK::Packed_OffsetYdim<+1>(offset),(SPGrid_Solver_Cell_Type_Active | SPGrid_Solver_Cell_Type_Interface | SPGrid_Solver_Cell_Type_Dirichlet)))
                        flag |= SPGrid_Solver_Face_Plus_Y_Active;
                }
            }
        }
    } 
};
//#####################################################################
}
#endif
