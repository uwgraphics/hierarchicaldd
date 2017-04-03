//#####################################################################
//Removed For Anonymity: Copyright Authors of Submission pap173s1
//Removed For Anonymity: Copyright Authors of Submission pap173s1
//#####################################################################
// Subroutine SPGrid_Computations::Linearized_Data_Copy
//#####################################################################
#ifndef __SPGrid_Linearized_Data_Copy_Helper_h__
#define __SPGrid_Linearized_Data_Copy_Helper_h__

#include <Common_Tools/Math_Tools/RANGE_ITERATOR.h>
#include <SPGrid/Core/SPGrid_Allocator.h>
#include <SPGrid_Fluids/Solvers/Domain_Decomposition/SPGRID_MULTIGRID_FLAGS.h>

namespace SPGrid_Computations{

using namespace SPGrid;
using namespace PhysBAM;

template<class T_STRUCT,class T,int d,class T_offset_ptr>
class Linearized_Data_Copy_To
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
    unsigned flags_to_check;
    const PhysBAM::HASHTABLE<unsigned long,T_offset_ptr>& offsets_map;

public:
    Linearized_Data_Copy_To(unsigned T_STRUCT::* flags_field_input,T T_STRUCT::* data_field_input,unsigned flags_to_check_input,char* linearized_data_input,const PhysBAM::HASHTABLE<unsigned long,T_offset_ptr>& offsets_map_input)
        :flags_field(flags_field_input),data_field(data_field_input),flags_to_check(flags_to_check_input),linearized_data(linearized_data_input),offsets_map(offsets_map_input)
    {}
    
    Linearized_Data_Copy_To(SPGrid_Allocator<T_STRUCT,d>& allocator,const std::pair<const unsigned long*,unsigned>& blocks,unsigned T_STRUCT::* flags_field_input,T T_STRUCT::* data_field_input,unsigned flags_to_check_input,char* linearized_data_input,const PhysBAM::HASHTABLE<unsigned long,T_offset_ptr>& offsets_map_input)
        :flags_field(flags_field_input),data_field(data_field_input),flags_to_check(flags_to_check_input),linearized_data(linearized_data_input),offsets_map(offsets_map_input)
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
            if(!offsets_map.Get(page_offset,linearized_page_offset)) continue;
            unsigned* linearized_flags_ptr=reinterpret_cast<unsigned*>((unsigned long)(linearized_data)+linearized_page_offset+offset_flags);
            T* linearized_data_ptr=reinterpret_cast<T*>((unsigned long)(linearized_data)+linearized_page_offset+offset_data);
            T_INDEX base_index=iterator.Index().template Cast<T_INDEX>();
            for(RANGE_ITERATOR<d> iterator(RANGE<T_INDEX>(base_index,base_index+allocator.Block_Size().template Cast<T_INDEX>()-1));
                iterator.Valid();
                iterator.Next(),offset+=sizeof(unsigned),++linearized_flags_ptr,++linearized_data_ptr){
                const unsigned& flag = flags(offset);
                const T& d = data(offset);
                if(((flag&SPGrid_Solver_PartitionID_Mask)==(flags_to_check&SPGrid_Solver_PartitionID_Mask))||//check subdomain id
                   ((flag&(~SPGrid_Solver_PartitionID_Mask))&(flags_to_check&(~SPGrid_Solver_PartitionID_Mask)))){//check other flags
                    //copy both flag and data
                    (*linearized_flags_ptr)=flag;
                    (*linearized_data_ptr)=d;
                }
            }
        }
    } 
};


template<class T_STRUCT,class T,int d,class T_offset_ptr>
class Linearized_Data_Copy_To_All
{
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<unsigned>::type Flag_Array_Type;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<const unsigned>::type Const_Flag_Array_Type;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<const T>::type Const_Data_Array_Type;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<unsigned>::mask T_MASK;
    typedef VECTOR<int,d> T_INDEX;
    enum{page_size=4096u};
    char* linearized_data;
    unsigned flags_to_check;
    const PhysBAM::HASHTABLE<unsigned long,T_offset_ptr>& offsets_map;

public:
    Linearized_Data_Copy_To_All(char* linearized_data_input,const PhysBAM::HASHTABLE<unsigned long,T_offset_ptr>& offsets_map_input)
        :linearized_data(linearized_data_input),offsets_map(offsets_map_input)
    {}
    
    Linearized_Data_Copy_To_All(SPGrid_Allocator<T_STRUCT,d>& allocator,const std::pair<const unsigned long*,unsigned>& blocks,char* linearized_data_input,const PhysBAM::HASHTABLE<unsigned long,T_offset_ptr>& offsets_map_input)
        :linearized_data(linearized_data_input),offsets_map(offsets_map_input)
    {Run(allocator,blocks);}
    
    void Run(SPGrid_Allocator<T_STRUCT,d>& allocator,const std::pair<const unsigned long*,unsigned>& blocks) const
    {
        //Assume that flags are the first field
        Const_Flag_Array_Type data=allocator.Get_Const_Array(&T_STRUCT::flags);
        for(SPGrid_Block_Iterator<T_MASK> iterator(blocks);iterator.Valid();iterator.Next_Block()){
            unsigned long offset=iterator.Offset();
            unsigned long page_offset=offset;
            T_offset_ptr linearized_page_offset;
            if(!offsets_map.Get(page_offset,linearized_page_offset)) {/*std::cout<<"missing offset encountered during copy!"<<std::endl;*/continue;}
            T* linearized_data_ptr=reinterpret_cast<T*>(reinterpret_cast<unsigned long>(linearized_data)+linearized_page_offset);
            memcpy(linearized_data_ptr,(void*)&data(offset),page_size);
        }
    } 
};

template<class T_STRUCT,class T,int d,class T_offset_ptr>
class Linearized_Data_Copy_From
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
    unsigned flags_to_check;
    const PhysBAM::HASHTABLE<unsigned long,T_offset_ptr>& offsets_map;

public:
    Linearized_Data_Copy_From(unsigned T_STRUCT::* flags_field_input,T T_STRUCT::* data_field_input,unsigned flags_to_check_input,const char* linearized_data_input,const PhysBAM::HASHTABLE<unsigned long,T_offset_ptr>& offsets_map_input)
        :flags_field(flags_field_input),data_field(data_field_input),flags_to_check(flags_to_check_input),linearized_data(linearized_data_input),offsets_map(offsets_map_input)
    {}
    
    Linearized_Data_Copy_From(SPGrid_Allocator<T_STRUCT,d>& allocator,const std::pair<const unsigned long*,unsigned>& blocks,unsigned T_STRUCT::* flags_field_input,T T_STRUCT::* data_field_input,unsigned flags_to_check_input,const char* linearized_data_input,const PhysBAM::HASHTABLE<unsigned long,T_offset_ptr>& offsets_map_input)
        :flags_field(flags_field_input),data_field(data_field_input),flags_to_check(flags_to_check_input),linearized_data(linearized_data_input),offsets_map(offsets_map_input)
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
            if(!offsets_map.Get(page_offset,linearized_page_offset)) continue;
            const unsigned* linearized_flags_ptr=reinterpret_cast<const unsigned*>(reinterpret_cast<unsigned long>(linearized_data)+linearized_page_offset+offset_flags);
            const T* linearized_data_ptr=reinterpret_cast<const T*>(reinterpret_cast<unsigned long>(linearized_data)+linearized_page_offset+offset_data);
            T_INDEX base_index=iterator.Index().template Cast<T_INDEX>();
            for(RANGE_ITERATOR<d> iterator(RANGE<T_INDEX>(base_index,base_index+allocator.Block_Size().template Cast<T_INDEX>()-1));
                iterator.Valid();
                iterator.Next(),offset+=sizeof(unsigned),++linearized_flags_ptr,++linearized_data_ptr){
                const unsigned& flag = flags(offset);
                T& d = data(offset);
                if(((flag&SPGrid_Solver_PartitionID_Mask)>0&&((flag&SPGrid_Solver_PartitionID_Mask)==(flags_to_check&SPGrid_Solver_PartitionID_Mask)))||//check subdomain id
                   ((flag&(~SPGrid_Solver_PartitionID_Mask))&(flags_to_check&(~SPGrid_Solver_PartitionID_Mask)))){//check other flags
                    //copy only data
                    PHYSBAM_ASSERT((*linearized_flags_ptr)==flag);
                    d=(*linearized_data_ptr);
                    //if(fabs(d) > norm) norm = fabs(d);
                }
            }
        }
        //std::cout << "Norm of the copied result: " << norm << std::endl;
    } 
};

template<class T_STRUCT,class T,int d,class T_offset_ptr>
class Linearized_Data_Copy_From_All
{
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<unsigned>::type Flag_Array_Type;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<const unsigned>::type Const_Flag_Array_Type;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<T>::type Data_Array_Type;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<unsigned>::mask T_MASK;
    typedef VECTOR<int,d> T_INDEX;
    enum{page_size=4096u};
    const char* linearized_data;
    const PhysBAM::HASHTABLE<unsigned long,T_offset_ptr>& offsets_map;

public:
    Linearized_Data_Copy_From_All(const char* linearized_data_input,const PhysBAM::HASHTABLE<unsigned long,T_offset_ptr>& offsets_map_input)
        :linearized_data(linearized_data_input),offsets_map(offsets_map_input)
    {}
    
    Linearized_Data_Copy_From_All(SPGrid_Allocator<T_STRUCT,d>& allocator,const std::pair<const unsigned long*,unsigned>& blocks,const char* linearized_data_input,const PhysBAM::HASHTABLE<unsigned long,T_offset_ptr>& offsets_map_input)
        :linearized_data(linearized_data_input),offsets_map(offsets_map_input)
    {Run(allocator,blocks);}
    
    void Run(SPGrid_Allocator<T_STRUCT,d>& allocator,const std::pair<const unsigned long*,unsigned>& blocks) const
    {
        //Assume that flags are the first field
        Flag_Array_Type data=allocator.Get_Array(&T_STRUCT::flags);
        for(SPGrid_Block_Iterator<T_MASK> iterator(blocks);iterator.Valid();iterator.Next_Block()){
            unsigned long offset=iterator.Offset();
            unsigned long page_offset=offset;
            T_offset_ptr linearized_page_offset;
            if(!offsets_map.Get(page_offset,linearized_page_offset)) {/*std::cout<<"missing offset encountered during copy!"<<std::endl;*/continue;}
            T* linearized_data_ptr=reinterpret_cast<T*>((unsigned long)linearized_data+linearized_page_offset);
            /*int entry=0;
            for(RANGE_ITERATOR<d> iterator(RANGE<T_INDEX>(T_INDEX(),allocator.Block_Size().template Cast<T_INDEX>()-1));
                iterator.Valid();
                iterator.Next(),++entry){
                PHYSBAM_ASSERT(((unsigned*)&data(offset))[entry]==((unsigned*)linearized_data_ptr)[entry]);
                    }*/
            memcpy((void*)&data(offset),linearized_data_ptr,page_size);
        }
    }
};

template<class T_STRUCT,class T,int d,class T_offset_ptr>
class Linearized_Data_Copy_Accumulative_Substract_From
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
    unsigned flags_to_check;
    const PhysBAM::HASHTABLE<unsigned long,T_offset_ptr>& offsets_map;

public:
    Linearized_Data_Copy_Accumulative_Substract_From(unsigned T_STRUCT::* flags_field_input,T T_STRUCT::* data_field_input,unsigned flags_to_check_input,const char* linearized_data_input,const PhysBAM::HASHTABLE<unsigned long,T_offset_ptr>& offsets_map_input)
        :flags_field(flags_field_input),data_field(data_field_input),flags_to_check(flags_to_check_input),linearized_data(linearized_data_input),offsets_map(offsets_map_input)
    {}
    
    Linearized_Data_Copy_Accumulative_Substract_From(SPGrid_Allocator<T_STRUCT,d>& allocator,const std::pair<const unsigned long*,unsigned>& blocks,unsigned T_STRUCT::* flags_field_input,T T_STRUCT::* data_field_input,unsigned flags_to_check_input,const char* linearized_data_input,const PhysBAM::HASHTABLE<unsigned long,T_offset_ptr>& offsets_map_input)
        :flags_field(flags_field_input),data_field(data_field_input),flags_to_check(flags_to_check_input),linearized_data(linearized_data_input),offsets_map(offsets_map_input)
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
            if(!offsets_map.Get(page_offset,linearized_page_offset)) continue;
            const unsigned* linearized_flags_ptr=reinterpret_cast<const unsigned*>(reinterpret_cast<unsigned long>(linearized_data)+linearized_page_offset+offset_flags);
            const T* linearized_data_ptr=reinterpret_cast<const T*>(reinterpret_cast<unsigned long>(linearized_data)+linearized_page_offset+offset_data);
            T_INDEX base_index=iterator.Index().template Cast<T_INDEX>();
            for(RANGE_ITERATOR<d> iterator(RANGE<T_INDEX>(base_index,base_index+allocator.Block_Size().template Cast<T_INDEX>()-1));
                iterator.Valid();
                iterator.Next(),offset+=sizeof(unsigned),++linearized_flags_ptr,++linearized_data_ptr){
                const unsigned& flag = flags(offset);
                T& d = data(offset);
                if(((flag&SPGrid_Solver_PartitionID_Mask)>0&&((flag&SPGrid_Solver_PartitionID_Mask)==(flags_to_check&SPGrid_Solver_PartitionID_Mask)))||//check subdomain id
                   ((flag&(~SPGrid_Solver_PartitionID_Mask))&(flags_to_check&(~SPGrid_Solver_PartitionID_Mask)))){//check other flags
                    //copy only data
                    //if((*linearized_flags_ptr)!=flag) LOG::cout<<flag<<" : "<<*linearized_flags_ptr<<std::endl;
                    PHYSBAM_ASSERT((*linearized_flags_ptr)==flag);
                    d-=(*linearized_data_ptr);
                    //if(fabs(d) > norm) norm = fabs(d);
                }
            }
        }
        //LOG::cout<<"Norm after accumulative substract: "<<norm<<std::endl;
    }
};
//###################################################################
// (DEBUG)
//###################################################################
template<class T_STRUCT,class T,int d,class T_offset_ptr>
class Linearized_Data_Compare
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
    unsigned flags_to_check;
    const PhysBAM::HASHTABLE<unsigned long,T_offset_ptr>& offsets_map;

public:
    Linearized_Data_Compare(unsigned T_STRUCT::* flags_field_input,T T_STRUCT::* data_field_input,unsigned flags_to_check_input,const char* linearized_data_input,const PhysBAM::HASHTABLE<unsigned long,T_offset_ptr>& offsets_map_input)
        :flags_field(flags_field_input),data_field(data_field_input),flags_to_check(flags_to_check_input),linearized_data(linearized_data_input),offsets_map(offsets_map_input)
    {}
    
    Linearized_Data_Compare(SPGrid_Allocator<T_STRUCT,d>& allocator,const std::pair<const unsigned long*,unsigned>& blocks,unsigned T_STRUCT::* flags_field_input,T T_STRUCT::* data_field_input,unsigned flags_to_check_input,const char* linearized_data_input,const PhysBAM::HASHTABLE<unsigned long,T_offset_ptr>& offsets_map_input)
        :flags_field(flags_field_input),data_field(data_field_input),flags_to_check(flags_to_check_input),linearized_data(linearized_data_input),offsets_map(offsets_map_input)
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
            if(!offsets_map.Get(page_offset,linearized_page_offset)) continue;
            const unsigned* linearized_flags_ptr=reinterpret_cast<const unsigned*>(reinterpret_cast<unsigned long>(linearized_data)+linearized_page_offset+offset_flags);
            const T* linearized_data_ptr=reinterpret_cast<const T*>(reinterpret_cast<unsigned long>(linearized_data)+linearized_page_offset+offset_data);
            T_INDEX base_index=iterator.Index().template Cast<T_INDEX>();
            for(RANGE_ITERATOR<d> iterator(RANGE<T_INDEX>(base_index,base_index+allocator.Block_Size().template Cast<T_INDEX>()-1));
                iterator.Valid();
                iterator.Next(),offset+=sizeof(unsigned),++linearized_flags_ptr,++linearized_data_ptr){
                const unsigned& flag = flags(offset);
                T& d = data(offset);
                if(((flag&SPGrid_Solver_PartitionID_Mask)>0&&((flag&SPGrid_Solver_PartitionID_Mask)==(flags_to_check&SPGrid_Solver_PartitionID_Mask)))||//check subdomain id
                   ((flag&(~SPGrid_Solver_PartitionID_Mask))&(flags_to_check&(~SPGrid_Solver_PartitionID_Mask)))){//check other flags
                    //copy only data
                    PHYSBAM_ASSERT((*linearized_flags_ptr)==flag);
                    if(d!=(*linearized_data_ptr))
                        std::cerr << "difference found at " << iterator.Index() << " : " << (*linearized_data_ptr) << " vs " << d << std::endl;;
                    //if(fabs(d) > norm) norm = fabs(d);
                }
            }
        }
        //std::cout << "Norm of the copied result: " << norm << std::endl;
    } 
};

}
#endif
