//#####################################################################
// Copyright (c) 2012, Eftychios Sifakis, Sean Bauer
// Distributed under the FreeBSD license (see license.txt)
//#####################################################################
#include <PhysBAM_Tools/Read_Write/Arrays/READ_WRITE_ARRAY_VIEW.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <SPGrid_Fluids/Grids/GRID_HIERARCHY.h>
#include <SPGrid_Fluids/Simulation/FLUIDS_SIMULATION_DATA.h>
#include <SPGrid_Fluids/Simulation/FLUIDS_SIMULATION_FLAGS.h>
#include <SPGrid_Fluids/Tools/ARRAY_BLOCK_HASH.h>
#include <SPGrid/Tools/SPGrid_Partitioning_Helper.h>

using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T_STRUCT,class T,int d> GRID_HIERARCHY<T_STRUCT,T,d>::
GRID_HIERARCHY(const T_GRID& base_grid_input,const int levels_input)
    :base_grid(base_grid_input),levels(levels_input)
{
    Initialize_Grids();
    Initialize_Allocators();
}
//#####################################################################
// Destructor
//#####################################################################
template<class T_STRUCT,class T,int d> GRID_HIERARCHY<T_STRUCT,T,d>::
~GRID_HIERARCHY()
{
    sets.Delete_Pointers_And_Clean_Memory();
    allocators.Delete_Pointers_And_Clean_Memory();
}
//#####################################################################
// Function Initialize_Grids
//#####################################################################
template<class T_STRUCT,class T,int d> void GRID_HIERARCHY<T_STRUCT,T,d>::
Initialize_Grids()
{
    grids.Resize(levels);
    grids(1)=base_grid;
    for(int level=2;level<=levels;level++){
        for(int v=1;v<=d;v++) PHYSBAM_ASSERT(grids(level-1).Numbers_Of_Cells()(v)%2==0);
        grids(level)=T_GRID(grids(level-1).Numbers_Of_Cells()/2+1,base_grid.Domain());}
}
//#####################################################################
// Function Initialize_Allocators
//#####################################################################
template<class T_STRUCT,class T,int d> void GRID_HIERARCHY<T_STRUCT,T,d>::
Initialize_Allocators()
{
    allocators.Resize(levels);
    for(int level=1;level<=levels;level++){
        coord_t size(grids(level).Numbers_Of_Cells()+2);
        allocators(level)=new Allocator_type(size);}
}
//#####################################################################
// Function Initialize_Sets
//#####################################################################
template<class T_STRUCT,class T,int d> void GRID_HIERARCHY<T_STRUCT,T,d>::
Initialize_Sets()
{
    sets.Resize(levels);
    for (int level=1;level<=levels;level++)
        sets(level)=new SPGrid_Set<Flag_array_type>(allocators(level)->Get_Array(&T_STRUCT::flags));
}
//#####################################################################
// Initialize_Red_Black_Partition
//#####################################################################
template<class T_STRUCT,class T,int d> void GRID_HIERARCHY<T_STRUCT,T,d>::
Initialize_Red_Black_Partition(const int number_of_partitions)
{
    // copy over blocks
    copy_of_blocks.Resize(Levels());
    for(int level=1;level<=Levels();level++){
        copy_of_blocks(level).resize(Blocks(level).second);
        for(int i=0;i<copy_of_blocks(level).size();i++)
            copy_of_blocks(level)[i]=(Blocks(level).first)[i];}

    // set up read black partition
    red_blocks.Resize(Levels());
    black_blocks.Resize(Levels());
    for(int level=1;level<=Levels();level++){
        SPGrid_Computations::Partitioning_Helper<T_STRUCT,d> partitioning_helper(Allocator(level),Set(level),copy_of_blocks(level));
        red_blocks(level)=std::vector<std::pair<const unsigned long*,unsigned> >();
        black_blocks(level)=std::vector<std::pair<const unsigned long*,unsigned> >();
        partitioning_helper.Generate_Red_Black_Partition(number_of_partitions,red_blocks(level),black_blocks(level));}

    // print information
    for(int level=1;level<=Levels();level++){
        LOG::cout<<"Level : "<<level<<std::endl;
        LOG::cout<<"  Red Blocks : "<<std::endl;
        for(int i=0;i<red_blocks(level).size();i++)
        LOG::cout<<"      # Blocks : "<<red_blocks(level)[i].second<<std::endl;
        LOG::cout<<"  Black Blocks : "<<std::endl;
        for(int i=0;i<black_blocks(level).size();i++)
        LOG::cout<<"      # Blocks : "<<black_blocks(level)[i].second<<std::endl;}
}
//#####################################################################
// Function Initialize_Boundary_Blocks
//#####################################################################
template<class T_STRUCT,class T,int d> void GRID_HIERARCHY<T_STRUCT,T,d>::
Initialize_Boundary_Blocks()
{
    typedef typename Const_flag_array_type::MASK T_MASK;
    boundary_blocks.Resize(levels);
    const unsigned boundary_mask(MG_Boundary);
    for(int level=1;level<=levels;level++) boundary_blocks(level).resize(0);
    for(int level=1;level<=levels;level++){Const_flag_array_type flags=Allocator(level).Get_Const_Array(&T_STRUCT::flags);
        for(SPGrid_Block_Iterator<T_MASK> iterator(Blocks(level));iterator.Valid();iterator.Next_Block()){
            const unsigned long block_offset=iterator.Offset();
            const unsigned* const flags_ptr(&iterator.Data(flags));
            for(int i=0;i<T_MASK::elements_per_block;i++)
                if(flags_ptr[i]&boundary_mask){boundary_blocks(level).push_back(block_offset);break;}}}
}
//#####################################################################
// Function Update_Block_Offsets
//#####################################################################
template<class T_STRUCT,class T,int d> void GRID_HIERARCHY<T_STRUCT,T,d>::
Update_Block_Offsets()
{
    for(int level=1;level<=levels;level++)
        sets(level)->Refresh_Block_Offsets();
}
//#####################################################################
// Function Clear_Bitmaps
//#####################################################################
template<class T_STRUCT,class T,int d> void GRID_HIERARCHY<T_STRUCT,T,d>::
Clear_Bitmaps()
{
    for (int level=1;level<=levels;level++)
        sets(level)->Clear_Bitmap();
}
//#####################################################################
// Function Print_Grids
//#####################################################################
template<class T_STRUCT,class T,int d> void GRID_HIERARCHY<T_STRUCT,T,d>::
Print_Grids(std::ostream& output)
{
    for(int level=1;level<=levels;level++)
        output<<"Level #"<<level<<" : Geometric domain "<<grids(level).Domain()
              <<" with cell indices "<<grids(level).Cell_Indices()<<std::endl;
}
//#####################################################################
// Function Hash
//#####################################################################
template<class T_STRUCT,class T,int d> int GRID_HIERARCHY<T_STRUCT,T,d>::
Hash(T T_STRUCT::* field)
{
    ARRAY<int> hashes;
    if(levels%2==0) hashes.Append(HASH().value);
    
    for(int level=1;level<=levels;level++){
        PHYSBAM_ASSERT(!sets(level) || !sets(level)->dirty);
        hashes.Append(ARRAY_BLOCK_HASH(Array(level,field),Blocks(level)));
        if(hashes.m==3){
            hashes(1)=PhysBAM::HASH(hashes(1),hashes(2),hashes(3)).value;
            hashes.Resize(1);}}
    
    PHYSBAM_ASSERT(hashes.m==1);
    return hashes(1);
}
//#####################################################################
// Function Write_Data_Channel
//#####################################################################
template<class T_STRUCT,class T,int d> void GRID_HIERARCHY<T_STRUCT,T,d>::
Write_Data_Channel(const std::string& filename,const T T_STRUCT::* field) const
{
    std::ostream* output=FILE_UTILITIES::Safe_Open_Output(filename);
    static const unsigned int elements_per_block = Data_array_type::MASK::elements_per_block;
    for(int level=1;level<=levels;level++){
        void* data_ptr=Array(level,field).Get_Data_Ptr();
        for(SPGrid_Block_Iterator<typename Const_data_array_type::MASK> iterator(Blocks(level));iterator.Valid();iterator.Next_Block()){
            PhysBAM::ARRAY_VIEW<const T> block_array(elements_per_block,&iterator.template Data<Const_data_array_type>(data_ptr));
            Write_Binary<T>(*output,block_array);}}
    delete output;
}
//#####################################################################
// Function Read_Data_Channel
//#####################################################################
template<class T_STRUCT,class T,int d> void GRID_HIERARCHY<T_STRUCT,T,d>::
Read_Data_Channel(const std::string& filename,T T_STRUCT::* field)
{
    std::istream* input=FILE_UTILITIES::Safe_Open_Input(filename);
    static const unsigned int elements_per_block = Data_array_type::MASK::elements_per_block;
    for(int level=1;level<=levels;level++){
        void* data_ptr=Array(level,field).Get_Data_Ptr();
        for(SPGrid_Block_Iterator<typename Data_array_type::MASK> iterator(Blocks(level));iterator.Valid();iterator.Next_Block()){
            PhysBAM::ARRAY_VIEW<T> block_array(elements_per_block,&iterator.template Data<Data_array_type>(data_ptr));
            Read_Binary<T>(*input,block_array);}}
    delete input;
}
//#####################################################################
// Function Write_Flags_Channel
//#####################################################################
template<class T_STRUCT,class T,int d> void GRID_HIERARCHY<T_STRUCT,T,d>::
Write_Flags_Channel(const std::string& filename) const
{
    std::ostream* output=FILE_UTILITIES::Safe_Open_Output(filename);
    static const unsigned int elements_per_block = Data_array_type::MASK::elements_per_block;
    for(int level=1;level<=levels;level++){
        void* data_ptr=Array(level,&T_STRUCT::flags).Get_Data_Ptr();
        for(SPGrid_Block_Iterator<typename Const_data_array_type::MASK> iterator(Blocks(level));iterator.Valid();iterator.Next_Block()){
            PhysBAM::ARRAY_VIEW<const unsigned> block_array(elements_per_block,&iterator.template Data<Const_flag_array_type>(data_ptr));
            Write_Binary<T>(*output,block_array);}}
    delete output;
}
//#####################################################################
// Function Read_Flags_Channel
//#####################################################################
template<class T_STRUCT,class T,int d> void GRID_HIERARCHY<T_STRUCT,T,d>::
Read_Flags_Channel(const std::string& filename)
{
    std::istream* input=FILE_UTILITIES::Safe_Open_Input(filename);
    static const unsigned int elements_per_block = Data_array_type::MASK::elements_per_block;
    for(int level=1;level<=levels;level++){
        void* data_ptr=Array(level,&T_STRUCT::flags).Get_Data_Ptr();
        for(SPGrid_Block_Iterator<typename Data_array_type::MASK> iterator(Blocks(level));iterator.Valid();iterator.Next_Block()){
            PhysBAM::ARRAY_VIEW<unsigned> block_array(elements_per_block,&iterator.template Data<Flag_array_type>(data_ptr));
            Read_Binary<T>(*input,block_array);}}
    delete input;
}
//#####################################################################
// Function Write_Block_Offsets
//#####################################################################
template<class T_STRUCT,class T,int d> void GRID_HIERARCHY<T_STRUCT,T,d>::
Write_Block_Offsets(const std::string& filename) const
{
    std::ostream* output=FILE_UTILITIES::Safe_Open_Output(filename);
    for(int level=1;level<=levels;level++){
        Write_Binary<T>(*output,Blocks(level).second);
        PhysBAM::ARRAY_VIEW<const unsigned long> blocks(Blocks(level).second,Blocks(level).first);
        Write_Binary<T>(*output,blocks);}
    delete output;
}
//#####################################################################
// Function Read_Block_Offsets
//#####################################################################
template<class T_STRUCT,class T,int d> void GRID_HIERARCHY<T_STRUCT,T,d>::
Read_Block_Offsets(const std::string& filename) const
{
    std::istream* input=FILE_UTILITIES::Safe_Open_Input(filename);
    for(int level=1;level<=levels;level++){
        int size;
        Read_Binary<T>(*input,size);
        std::vector<unsigned long>& block_offsets=sets(level)->block_offsets;
        block_offsets.resize(size);
        PhysBAM::ARRAY_VIEW<unsigned long> blocks(Blocks(level).second,&block_offsets[0]);
        Read_Binary<T>(*input,blocks);
        sets(level)->Clear_Bitmap();
        sets(level)->FillBitmapWithBlockOffsets(block_offsets);}
    delete input;
}
//#####################################################################
template class GRID_HIERARCHY<FLUIDS_SIMULATION_DATA<float>,float,2>;
template class GRID_HIERARCHY<FLUIDS_SIMULATION_DATA<float>,float,3>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class GRID_HIERARCHY<FLUIDS_SIMULATION_DATA<double>,double,2>;
template class GRID_HIERARCHY<FLUIDS_SIMULATION_DATA<double>,double,3>;
#endif
