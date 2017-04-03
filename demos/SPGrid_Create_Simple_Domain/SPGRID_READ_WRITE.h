//#####################################################################
//Removed For Anonymity: Copyright Authors of Submission pap173s1
// Distributed under the FreeBSD license (see license.txt)
//#####################################################################
#ifndef __SPGRID_READ_WRITE_H__
#define __SPGRID_READ_WRITE_H__
#include <SPGrid/Core/SPGrid_Allocator.h>
#include <SPGrid/Core/SPGrid_Set.h>
#include <SPGrid/Tools/SPGrid_Block_Iterator.h>
#include <PhysBAM_Tools/Read_Write/Arrays/READ_WRITE_ARRAY_VIEW.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Tools/Read_Write/Vectors/READ_WRITE_VECTOR.h>

namespace PhysBAM{
using namespace SPGrid;
//#####################################################################
// Class SPGRID_READ_WRITE
//#####################################################################
template<class T_STRUCT,class T,int d>
class SPGRID_READ_WRITE
{
    typedef float RW;
    typedef SPGrid_Allocator<T_STRUCT,d> Allocator_type;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<unsigned>::type Flag_array_type;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<const unsigned>::type Const_flag_array_type;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<T>::type Data_array_type;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<const T>::type Const_data_array_type;
    typedef SPGrid_Set<Flag_array_type> Set_type;
    typedef VECTOR<int,d> T_INDEX;
public:
    static void Read_SPGrid(const std::string& base_dir,Allocator_type*& allocator_data,Allocator_type*& allocator_set,Set_type*& set,T T_STRUCT::* field){
        RW rw=RW(); STREAM_TYPE stream_type(rw);
        T_INDEX size;
        FILE_UTILITIES::Read_From_File(stream_type,base_dir+"/size",size);
        LOG::cout<<"Size: "<<size<<std::endl;
        allocator_set=new Allocator_type(std_array<int,d>(size));
        allocator_data=new Allocator_type(std_array<int,d>(size));
        set=new Set_type(allocator_set->Get_Array(&T_STRUCT::flags));
        Read_Block_Offsets(base_dir+"/block_offset",*set);
        Read_Data_Channel(base_dir+"/data",*allocator_data,*set,field);
        Read_Flags_Channel(base_dir+"/flags",*allocator_set,*set);
    }
    static void Write_SPGrid(const std::string& base_dir,Allocator_type& allocator,Set_type& set,const T T_STRUCT::* field){
        FILE_UTILITIES::Create_Directory(base_dir);
        RW rw=RW(); STREAM_TYPE stream_type(rw);
        T_INDEX size=allocator.Padded_Size().template Cast<T_INDEX>();
        FILE_UTILITIES::Write_To_File(stream_type,base_dir+"/size",size);
        Write_Block_Offsets(base_dir+"/block_offset",set);
        Write_Data_Channel(base_dir+"/data",allocator,set,field);
        Write_Flags_Channel(base_dir+"/flags",allocator,set);
    }
    static void Write_Data_Channel(const std::string& filename,Allocator_type& allocator,Set_type& set,const T T_STRUCT::* field){
        std::ostream* output=FILE_UTILITIES::Safe_Open_Output(filename);
        static const unsigned int elements_per_block = Data_array_type::MASK::elements_per_block;
        void* data_ptr=allocator.Get_Array(field).Get_Data_Ptr();
        for(SPGrid_Block_Iterator<typename Const_data_array_type::MASK> iterator(set.Get_Blocks());iterator.Valid();iterator.Next_Block()){
            PhysBAM::ARRAY_VIEW<const T> block_array(elements_per_block,&iterator.template Data<Const_data_array_type>(data_ptr));
            Write_Binary<T>(*output,block_array);}
        delete output;
    }
    static void Read_Data_Channel(const std::string& filename,Allocator_type& allocator,Set_type& set,T T_STRUCT::* field){
        std::istream* input=FILE_UTILITIES::Safe_Open_Input(filename);
        static const unsigned int elements_per_block = Data_array_type::MASK::elements_per_block;
        void* data_ptr=allocator.Get_Array(field).Get_Data_Ptr();
        for(SPGrid_Block_Iterator<typename Data_array_type::MASK> iterator(set.Get_Blocks());iterator.Valid();iterator.Next_Block()){
            PhysBAM::ARRAY_VIEW<T> block_array(elements_per_block,&iterator.template Data<Data_array_type>(data_ptr));
            Read_Binary<T>(*input,block_array);}
        delete input;
    }
    static void Write_Flags_Channel(const std::string& filename,Allocator_type& allocator,Set_type& set){
        std::ostream* output=FILE_UTILITIES::Safe_Open_Output(filename);
        static const unsigned int elements_per_block = Flag_array_type::MASK::elements_per_block;
        void* data_ptr=allocator.Get_Array(&T_STRUCT::flags).Get_Data_Ptr();
        for(SPGrid_Block_Iterator<typename Const_data_array_type::MASK> iterator(set.Get_Blocks());iterator.Valid();iterator.Next_Block()){
            PhysBAM::ARRAY_VIEW<const unsigned> block_array(elements_per_block,&iterator.template Data<Const_flag_array_type>(data_ptr));
            Write_Binary<T>(*output,block_array);}
        delete output;
    }
    static void Read_Flags_Channel(const std::string& filename,Allocator_type& allocator,Set_type& set){
        std::istream* input=FILE_UTILITIES::Safe_Open_Input(filename);
        static const unsigned int elements_per_block = Data_array_type::MASK::elements_per_block;
        void* data_ptr=allocator.Get_Array(&T_STRUCT::flags).Get_Data_Ptr();
        for(SPGrid_Block_Iterator<typename Data_array_type::MASK> iterator(set.Get_Blocks());iterator.Valid();iterator.Next_Block()){
            PhysBAM::ARRAY_VIEW<unsigned> block_array(elements_per_block,&iterator.template Data<Flag_array_type>(data_ptr));
            Read_Binary<T>(*input,block_array);}
        delete input;
    }
    static void Write_Block_Offsets(const std::string& filename,Set_type& set){
        std::ostream* output=FILE_UTILITIES::Safe_Open_Output(filename);
        Write_Binary<T>(*output,set.Get_Blocks().second);
        PhysBAM::ARRAY_VIEW<const unsigned long> blocks(set.Get_Blocks().second,set.Get_Blocks().first);
        Write_Binary<T>(*output,blocks);
        delete output;
    }
    static void Read_Block_Offsets(const std::string& filename,Set_type& set){
        std::istream* input=FILE_UTILITIES::Safe_Open_Input(filename);
        int size;
        Read_Binary<T>(*input,size);
        std::vector<unsigned long>& block_offsets=set.block_offsets;
        block_offsets.resize(size);
        PhysBAM::ARRAY_VIEW<unsigned long> blocks(set.Get_Blocks().second,&block_offsets[0]);
        Read_Binary<T>(*input,blocks);
        set.Clear_Bitmap();
        set.FillBitmapWithBlockOffsets(block_offsets);
        delete input;
    }
};
}
#endif
