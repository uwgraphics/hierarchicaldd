//#####################################################################
//Removed For Anonymity: Copyright Authors of Submission pap173s1
//Removed For Anonymity: Copyright Authors of Submission pap173s1
//#####################################################################
// Subroutine SPGrid_Computations::Linearized_Data_Copy_Hashtable
//#####################################################################
#ifndef __SPGrid_Linearized_Data_Copy_Hashtable_Helper_h__
#define __SPGrid_Linearized_Data_Copy_Hashtable_Helper_h__
#include <Common_Tools/Math_Tools/RANGE_ITERATOR.h>
#include <SPGrid/Core/SPGrid_Allocator.h>
#include <SPGrid_Fluids/Solvers/Domain_Decomposition/SPGRID_MULTIGRID_FLAGS.h>
#include <PhysBAM_Tools/Data_Structures/HASHTABLE.h>
#include <PhysBAM_Tools/Data_Structures/HASHTABLE_ITERATOR.h>
namespace SPGrid_Computations{
using namespace SPGrid;
using namespace PhysBAM;
template<class T_STRUCT,class T,int d,class T_offset_ptr>
class Copy_To_Linearized_Data_With_Interface_Clear
{
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<unsigned>::type Flag_Array_Type;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<const unsigned>::type Const_Flag_Array_Type;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<const T>::type Const_Data_Array_Type;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<unsigned>::mask T_MASK;
    typedef VECTOR<int,d> T_INDEX;
    enum{page_size=4096u};
    unsigned T_STRUCT::* flags_field;
    T T_STRUCT::* data_field;
    char* linearized_data;
    const PhysBAM::HASHTABLE<unsigned long,T_offset_ptr>& offsets_map;

public:
    Copy_To_Linearized_Data_With_Interface_Clear(unsigned T_STRUCT::* flags_field_input,T T_STRUCT::* data_field_input,char* linearized_data_input,const PhysBAM::HASHTABLE<unsigned long,T_offset_ptr>& offsets_map_input)
        :flags_field(flags_field_input),data_field(data_field_input),linearized_data(linearized_data_input),offsets_map(offsets_map_input)
    {}
    
    Copy_To_Linearized_Data_With_Interface_Clear(SPGrid_Allocator<T_STRUCT,d>& allocator,const std::pair<const unsigned long*,unsigned>& blocks,unsigned T_STRUCT::* flags_field_input,T T_STRUCT::* data_field_input,char* linearized_data_input,const PhysBAM::HASHTABLE<unsigned long,T_offset_ptr>& offsets_map_input)
        :flags_field(flags_field_input),data_field(data_field_input),linearized_data(linearized_data_input),offsets_map(offsets_map_input)
    {Run(allocator,blocks);}
    
    void Run(SPGrid_Allocator<T_STRUCT,d>& allocator,const std::pair<const unsigned long*,unsigned>& blocks) const
    {
        Const_Flag_Array_Type flags=allocator.Get_Const_Array(flags_field);
        Const_Data_Array_Type data=allocator.Get_Const_Array(data_field);
        size_t offset_flags=OffsetOfMember(flags_field)<<T_MASK::block_bits;
        size_t offset_data=OffsetOfMember(data_field)<<T_MASK::block_bits;
        for(SPGrid_Block_Iterator<T_MASK> iterator(blocks);iterator.Valid();iterator.Next_Block()){
            unsigned long offset=iterator.Offset();
            unsigned long page_offset=offset;
            T_offset_ptr linearized_page_offset;
            PHYSBAM_ASSERT(offsets_map.Get(page_offset,linearized_page_offset));
            //if(!offsets_map.Get(page_offset,linearized_page_offset)) continue;
            unsigned* linearized_flags_ptr=reinterpret_cast<unsigned*>((unsigned long)(linearized_data)+linearized_page_offset+offset_flags);
            T* linearized_data_ptr=reinterpret_cast<T*>((unsigned long)(linearized_data)+linearized_page_offset+offset_data);
            T_INDEX base_index=iterator.Index().template Cast<T_INDEX>();
            for(RANGE_ITERATOR<d> block_iterator(RANGE<T_INDEX>(base_index,base_index+allocator.Block_Size().template Cast<T_INDEX>()-1));
                block_iterator.Valid();
                block_iterator.Next(),offset+=sizeof(unsigned),++linearized_flags_ptr,++linearized_data_ptr){
                const unsigned& flag = flags(offset);
                const T& d = data(offset);
                if((*linearized_flags_ptr)&SPGrid_Solver_Cell_Type_Interface){            
                    PHYSBAM_ASSERT(flag&SPGrid_Solver_Cell_Type_Interface);
                    (*linearized_data_ptr)=T();}
                else if((*linearized_flags_ptr)&SPGrid_Solver_Cell_Type_Active){
                    PHYSBAM_ASSERT(flag&SPGrid_Solver_Cell_Type_Active);
                    (*linearized_data_ptr)=d;}
            }
        }
    } 
};

template<class T_STRUCT,class T,int d,class T_offset_ptr>
class Copy_From_Linearized_Data
{
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<const unsigned>::type Const_Flag_Array_Type;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<T>::type Data_Array_Type;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<unsigned>::mask T_MASK;
    typedef VECTOR<int,d> T_INDEX;
    enum{page_size=4096u};
    unsigned T_STRUCT::* flags_field;
    T T_STRUCT::* data_field;
    char* linearized_data;
    const PhysBAM::HASHTABLE<unsigned long,T_offset_ptr>& offsets_map;

public:
    Copy_From_Linearized_Data(unsigned T_STRUCT::* flags_field_input,T T_STRUCT::* data_field_input,char* linearized_data_input,const PhysBAM::HASHTABLE<unsigned long,T_offset_ptr>& offsets_map_input)
        :flags_field(flags_field_input),data_field(data_field_input),linearized_data(linearized_data_input),offsets_map(offsets_map_input)
    {}
    
    Copy_From_Linearized_Data(SPGrid_Allocator<T_STRUCT,d>& allocator,const std::pair<const unsigned long*,unsigned>& blocks,unsigned T_STRUCT::* flags_field_input,T T_STRUCT::* data_field_input,char* linearized_data_input,const PhysBAM::HASHTABLE<unsigned long,T_offset_ptr>& offsets_map_input)
        :flags_field(flags_field_input),data_field(data_field_input),linearized_data(linearized_data_input),offsets_map(offsets_map_input)
    {Run(allocator,blocks);}
    
    void Run(SPGrid_Allocator<T_STRUCT,d>& allocator,const std::pair<const unsigned long*,unsigned>& blocks) const
    {
        Const_Flag_Array_Type flags=allocator.Get_Const_Array(flags_field);
        Data_Array_Type data=allocator.Get_Array(data_field);
        size_t offset_flags=OffsetOfMember(flags_field)<<T_MASK::block_bits;
        size_t offset_data=OffsetOfMember(data_field)<<T_MASK::block_bits;
        for(SPGrid_Block_Iterator<T_MASK> iterator(blocks);iterator.Valid();iterator.Next_Block()){
            unsigned long offset=iterator.Offset();
            unsigned long page_offset=offset;
            T_offset_ptr linearized_page_offset;
            PHYSBAM_ASSERT(offsets_map.Get(page_offset,linearized_page_offset));
            //if(!offsets_map.Get(page_offset,linearized_page_offset)) continue;
            const unsigned* linearized_flags_ptr=reinterpret_cast<unsigned*>((unsigned long)(linearized_data)+linearized_page_offset+offset_flags);
            const T* linearized_data_ptr=reinterpret_cast<T*>((unsigned long)(linearized_data)+linearized_page_offset+offset_data);
            T_INDEX base_index=iterator.Index().template Cast<T_INDEX>();
            for(RANGE_ITERATOR<d> block_iterator(RANGE<T_INDEX>(base_index,base_index+allocator.Block_Size().template Cast<T_INDEX>()-1));
                block_iterator.Valid();
                block_iterator.Next(),offset+=sizeof(unsigned),++linearized_flags_ptr,++linearized_data_ptr){
                const unsigned& flag = flags(offset);
                T& d = data(offset);
                if((*linearized_flags_ptr)&SPGrid_Solver_Cell_Type_Active){
                    PHYSBAM_ASSERT(flag&SPGrid_Solver_Cell_Type_Active);
                    d=(*linearized_data_ptr);}
            }
        }
    } 
};

template<class T_STRUCT,class T,int d,class T_offset_ptr>
class Accumulatively_Interface_Substract_From_Linearized_Data
{
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<unsigned>::type Flag_Array_Type;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<const unsigned>::type Const_Flag_Array_Type;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<T>::type Data_Array_Type;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<unsigned>::mask T_MASK;
    typedef VECTOR<int,d> T_INDEX;
    enum{page_size=4096u};
    unsigned T_STRUCT::* flags_field;
    T T_STRUCT::* data_field;
    const char* linearized_data;
    const PhysBAM::HASHTABLE<unsigned long,T_offset_ptr>& offsets_map;

public:
    Accumulatively_Interface_Substract_From_Linearized_Data(unsigned T_STRUCT::* flags_field_input,T T_STRUCT::* data_field_input,const char* linearized_data_input,const PhysBAM::HASHTABLE<unsigned long,T_offset_ptr>& offsets_map_input)
        :flags_field(flags_field_input),data_field(data_field_input),linearized_data(linearized_data_input),offsets_map(offsets_map_input)
    {}
    
    Accumulatively_Interface_Substract_From_Linearized_Data(SPGrid_Allocator<T_STRUCT,d>& allocator,const std::pair<const unsigned long*,unsigned>& blocks,unsigned T_STRUCT::* flags_field_input,T T_STRUCT::* data_field_input,const char* linearized_data_input,const PhysBAM::HASHTABLE<unsigned long,T_offset_ptr>& offsets_map_input)
        :flags_field(flags_field_input),data_field(data_field_input),linearized_data(linearized_data_input),offsets_map(offsets_map_input)
    {Run(allocator,blocks);}
    
    void Run(SPGrid_Allocator<T_STRUCT,d>& allocator,const std::pair<const unsigned long*,unsigned>& blocks) const
    {
        //T norm = 0;
        Const_Flag_Array_Type flags=allocator.Get_Const_Array(flags_field);
        Data_Array_Type data=allocator.Get_Array(data_field);
        size_t offset_flags=OffsetOfMember(flags_field)<<T_MASK::block_bits;
        size_t offset_data=OffsetOfMember(data_field)<<T_MASK::block_bits;
        for(SPGrid_Block_Iterator<T_MASK> iterator(blocks);iterator.Valid();iterator.Next_Block()){
            unsigned long offset=iterator.Offset();
            unsigned long page_offset=offset;
            T_offset_ptr linearized_page_offset;
            PHYSBAM_ASSERT(offsets_map.Get(page_offset,linearized_page_offset));
            //if(!offsets_map.Get(page_offset,linearized_page_offset)) continue;
            const unsigned* linearized_flags_ptr=reinterpret_cast<const unsigned*>(reinterpret_cast<unsigned long>(linearized_data)+linearized_page_offset+offset_flags);
            const T* linearized_data_ptr=reinterpret_cast<const T*>(reinterpret_cast<unsigned long>(linearized_data)+linearized_page_offset+offset_data);
            T_INDEX base_index=iterator.Index().template Cast<T_INDEX>();
            for(RANGE_ITERATOR<d> iterator(RANGE<T_INDEX>(base_index,base_index+allocator.Block_Size().template Cast<T_INDEX>()-1));
                iterator.Valid();
                iterator.Next(),offset+=sizeof(unsigned),++linearized_flags_ptr,++linearized_data_ptr){
                const unsigned& flag = flags(offset);
                T& d = data(offset);
                if((*linearized_flags_ptr)&SPGrid_Solver_Cell_Type_Interface){
                    //PHYSBAM_ASSERT(flag==(*linearized_flags_ptr));
                    PHYSBAM_ASSERT(flag&SPGrid_Solver_Cell_Type_Interface);
                    d-=(*linearized_data_ptr);
                }
            }
        }
    }
};

template<class T_STRUCT,class T,int d,class T_offset_ptr>
class Copy_Interface_To_Linearized_Data
{
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<unsigned>::type Flag_Array_Type;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<const unsigned>::type Const_Flag_Array_Type;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<const T>::type Const_Data_Array_Type;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<unsigned>::mask T_MASK;
    typedef VECTOR<int,d> T_INDEX;
    enum{page_size=4096u};
    unsigned T_STRUCT::* flags_field;
    T T_STRUCT::* data_field;
    char* linearized_data;
    const PhysBAM::HASHTABLE<unsigned long,T_offset_ptr>& offsets_map;

public:
    Copy_Interface_To_Linearized_Data(unsigned T_STRUCT::* flags_field_input,T T_STRUCT::* data_field_input,char* linearized_data_input,const PhysBAM::HASHTABLE<unsigned long,T_offset_ptr>& offsets_map_input)
        :flags_field(flags_field_input),data_field(data_field_input),linearized_data(linearized_data_input),offsets_map(offsets_map_input)
    {}
    
    Copy_Interface_To_Linearized_Data(SPGrid_Allocator<T_STRUCT,d>& allocator,const std::pair<const unsigned long*,unsigned>& blocks,unsigned T_STRUCT::* flags_field_input,T T_STRUCT::* data_field_input,char* linearized_data_input,const PhysBAM::HASHTABLE<unsigned long,T_offset_ptr>& offsets_map_input)
        :flags_field(flags_field_input),data_field(data_field_input),linearized_data(linearized_data_input),offsets_map(offsets_map_input)
    {Run(allocator,blocks);}
    
    void Run(SPGrid_Allocator<T_STRUCT,d>& allocator,const std::pair<const unsigned long*,unsigned>& blocks) const
    {
        Const_Flag_Array_Type flags=allocator.Get_Const_Array(flags_field);
        Const_Data_Array_Type data=allocator.Get_Const_Array(data_field);
        size_t offset_flags=OffsetOfMember(flags_field)<<T_MASK::block_bits;
        size_t offset_data=OffsetOfMember(data_field)<<T_MASK::block_bits;
        for(SPGrid_Block_Iterator<T_MASK> iterator(blocks);iterator.Valid();iterator.Next_Block()){
            unsigned long offset=iterator.Offset();
            unsigned long page_offset=offset;
            T_offset_ptr linearized_page_offset;
            PHYSBAM_ASSERT(offsets_map.Get(page_offset,linearized_page_offset));
            //if(!offsets_map.Get(page_offset,linearized_page_offset)) continue;
            unsigned* linearized_flags_ptr=reinterpret_cast<unsigned*>((unsigned long)(linearized_data)+linearized_page_offset+offset_flags);
            T* linearized_data_ptr=reinterpret_cast<T*>((unsigned long)(linearized_data)+linearized_page_offset+offset_data);
            T_INDEX base_index=iterator.Index().template Cast<T_INDEX>();
            for(RANGE_ITERATOR<d> block_iterator(RANGE<T_INDEX>(base_index,base_index+allocator.Block_Size().template Cast<T_INDEX>()-1));
                block_iterator.Valid();
                block_iterator.Next(),offset+=sizeof(unsigned),++linearized_flags_ptr,++linearized_data_ptr){
                const unsigned& flag = flags(offset);
                const T& d = data(offset);
                if((*linearized_flags_ptr)&SPGrid_Solver_Cell_Type_Interface){            
                    PHYSBAM_ASSERT(flag&SPGrid_Solver_Cell_Type_Interface);
                    (*linearized_data_ptr)=d;}
            }
        }
    } 
};
    ////////////////////////////////////////////////DEBUG/////////////////////////////////////////////
template<class T_STRUCT,class T,int d,class T_offset_ptr>
class Set_One_From_Linearized_Data
{
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<const unsigned>::type Const_Flag_Array_Type;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<T>::type Data_Array_Type;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<unsigned>::mask T_MASK;
    typedef VECTOR<int,d> T_INDEX;
    enum{page_size=4096u};
    unsigned T_STRUCT::* flags_field;
    T T_STRUCT::* data_field;
    char* linearized_data;
    const PhysBAM::HASHTABLE<unsigned long,T_offset_ptr>& offsets_map;

public:
    Set_One_From_Linearized_Data(unsigned T_STRUCT::* flags_field_input,T T_STRUCT::* data_field_input,char* linearized_data_input,const PhysBAM::HASHTABLE<unsigned long,T_offset_ptr>& offsets_map_input)
        :flags_field(flags_field_input),data_field(data_field_input),linearized_data(linearized_data_input),offsets_map(offsets_map_input)
    {}
    
    Set_One_From_Linearized_Data(SPGrid_Allocator<T_STRUCT,d>& allocator,const std::pair<const unsigned long*,unsigned>& blocks,unsigned T_STRUCT::* flags_field_input,T T_STRUCT::* data_field_input,char* linearized_data_input,const PhysBAM::HASHTABLE<unsigned long,T_offset_ptr>& offsets_map_input)
        :flags_field(flags_field_input),data_field(data_field_input),linearized_data(linearized_data_input),offsets_map(offsets_map_input)
    {Run(allocator,blocks);}
    
    void Run(SPGrid_Allocator<T_STRUCT,d>& allocator,const std::pair<const unsigned long*,unsigned>& blocks) const
    {
        Const_Flag_Array_Type flags=allocator.Get_Const_Array(flags_field);
        Data_Array_Type data=allocator.Get_Array(data_field);
        size_t offset_flags=OffsetOfMember(flags_field)<<T_MASK::block_bits;
        size_t offset_data=OffsetOfMember(data_field)<<T_MASK::block_bits;
        for(SPGrid_Block_Iterator<T_MASK> iterator(blocks);iterator.Valid();iterator.Next_Block()){
            unsigned long offset=iterator.Offset();
            unsigned long page_offset=offset;
            T_offset_ptr linearized_page_offset;
            PHYSBAM_ASSERT(offsets_map.Get(page_offset,linearized_page_offset));
            //if(!offsets_map.Get(page_offset,linearized_page_offset)) continue;
            const unsigned* linearized_flags_ptr=reinterpret_cast<unsigned*>((unsigned long)(linearized_data)+linearized_page_offset+offset_flags);
            const T* linearized_data_ptr=reinterpret_cast<T*>((unsigned long)(linearized_data)+linearized_page_offset+offset_data);
            T_INDEX base_index=iterator.Index().template Cast<T_INDEX>();
            for(RANGE_ITERATOR<d> block_iterator(RANGE<T_INDEX>(base_index,base_index+allocator.Block_Size().template Cast<T_INDEX>()-1));
                block_iterator.Valid();
                block_iterator.Next(),offset+=sizeof(unsigned),++linearized_flags_ptr,++linearized_data_ptr){
                const unsigned& flag = flags(offset);
                T& d = data(offset);
                if((*linearized_flags_ptr)&SPGrid_Solver_Cell_Type_Active){
                    PHYSBAM_ASSERT(flag&SPGrid_Solver_Cell_Type_Active);
                    d=T(1);}
            }
        }
    } 
};

}
#endif
