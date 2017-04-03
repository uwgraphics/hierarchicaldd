//#####################################################################
//Removed For Anonymity: Copyright Authors of Submission pap173s1
//Removed For Anonymity: Copyright Authors of Submission pap173s1
//#####################################################################
#ifndef __SPGRID_CONVERTOR_H__
#define __SPGRID_CONVERTOR_H__

#include <SPGrid/Core/SPGrid_Allocator.h>
#include <SPGrid/Tools/SPGrid_Block_Iterator.h>
#include <SPGrid_Fluids/Grids/GRID_TOPOLOGY_HELPER.h>
#include <SPGrid_Fluids/Simulation/FLUIDS_SIMULATION_FLAGS.h>
#include <SPGrid_Fluids/Solvers/Domain_Decomposition/SPGRID_MULTIGRID_FLAGS.h>
#include "SPGrid_Gemini.h"

namespace SPGrid_Computations{

using namespace PhysBAM;
using namespace SPGrid;

inline unsigned long Covert_Block_Offset_16_to_4(unsigned long offset_16){return (offset_16>>2)&(0xfffffffffffff000u);}
inline unsigned long Covert_Block_Offset_4_to_16(unsigned long offset_4){return (offset_4&(0xfffffffffffff000u))<<2;}
    
template<class T_STRUCT_16,class T_STRUCT_4,class T_FLAGS,int d>
class SPGrid_Creator_16_to_4{
    //This create TWO 4 channel SPGrid and mark the block offsets from a 16 channel SPGrid
    typedef SPGrid_Allocator<T_STRUCT_16,d> Allocator_Type_16;
    typedef typename Allocator_Type_16::Array<unsigned>::type SPG_Flags_Array_Type_16;
    typedef typename Allocator_Type_16::Array<unsigned>::mask T_MASK_16;
    typedef SPGrid_Set<SPG_Flags_Array_Type_16> SPG_Set_Type_16;

    typedef SPGrid_Allocator<T_STRUCT_4,d> Allocator_Type_4;
    typedef typename Allocator_Type_4::Array<unsigned>::mask T_MASK_4;
    typedef SPGrid_Gemini<T_STRUCT_4,d> Gemini;

public:
    static void Create(Gemini& gemini,Allocator_Type_16& allocator_16,SPG_Set_Type_16& set_16){
        static_assert(d==3,"Creation from 16 channel of 4 channel SPGrid. Only 3D is supported");
        if(gemini.allocator_first!=0) delete gemini.allocator_first;
        gemini.allocator_first=new typename Gemini::SPG_Allocator(allocator_16.Padded_Size());
        
        if(gemini.allocator_second!=0) delete gemini.allocator_second;
        gemini.allocator_second=new typename Gemini::SPG_Allocator(allocator_16.Padded_Size());
        
        typename Gemini::SPG_Gemini_Flags_Array_Type gemini_flags=gemini.allocator_first->Get_Array(&T_STRUCT_4::flags);
        if(gemini.set!=0) delete gemini.set;
        gemini.set=new typename Gemini::SPG_Gemini_Set_Type(gemini_flags);

        for(SPGrid_Block_Iterator<T_MASK_16> iterator(set_16.Get_Blocks());iterator.Valid();iterator.Next_Block()){
            unsigned long offset_16=iterator.Offset();
            unsigned long offset_4=Covert_Block_Offset_16_to_4(offset_16);
            //To convert 16 channel SPGrid to a 4 channel SPGrid. Each 4 channel block is always consist of 4 continues 16 channel blocks.
            gemini.set->MarkPageActive(offset_4);}
        gemini.set->Refresh_Block_Offsets();
    }
};

template<class T_STRUCT_16,class T_STRUCT_4,class T_DATA,int d> class SPGrid_Convertor_16_To_4;
template<class T_STRUCT_16,class T_STRUCT_4,class T_DATA>
class SPGrid_Convertor_16_To_4<T_STRUCT_16,T_STRUCT_4,T_DATA,3>
{
    enum{d=3};
    //This is intend to convert between a 16 channel SPGrid and a 4 channel SPGrid. Do not use it for any other purposes!
    typedef SPGrid_Allocator<T_STRUCT_16,d> Allocator_Type_16;
    typedef typename SPGrid_Allocator<T_STRUCT_16,d>::template Array<const T_DATA>::type Const_Data_Array_Type_16;
    typedef typename SPGrid_Allocator<T_STRUCT_16,d>::Array<unsigned>::type Flags_Array_Type_16;
    typedef SPGrid_Set<Flags_Array_Type_16> Set_Type_16;

    typedef SPGrid_Allocator<T_STRUCT_4,d> Allocator_Type_4;
    typedef typename SPGrid_Allocator<T_STRUCT_4,d>::template Array<T_DATA>::type Data_Array_Type_4;
    typedef typename SPGrid_Allocator<T_STRUCT_4,d>::template Array<T_DATA>::mask T_MASK_4;

    enum{block_4_xsize=1u<<T_MASK_4::block_xbits,
         block_4_ysize=1u<<T_MASK_4::block_ybits,
         block_4_zsize=1u<<T_MASK_4::block_ybits};

    typedef VECTOR<int,d> T_INDEX;

    Allocator_Type_16& allocator_16;
    Set_Type_16& set_16;
    T_DATA T_STRUCT_16::* data_field_16;
    T_DATA T_STRUCT_4::* data_field_4;

public:
    SPGrid_Convertor_16_To_4(Allocator_Type_16& allocator_16_in,Set_Type_16& set_16_in,
                             T_DATA T_STRUCT_16::* data_field_16_input,T_DATA T_STRUCT_4::* data_field_4_input)
        :allocator_16(allocator_16_in),set_16(set_16_in),data_field_16(data_field_16_input),data_field_4(data_field_4_input){}

    SPGrid_Convertor_16_To_4(Allocator_Type_4& allocator_4,const std::pair<const unsigned long*,unsigned>&blocks_4,
                             Allocator_Type_16& allocator_16_in,Set_Type_16& set_16_in,
                             T_DATA T_STRUCT_16::* data_field_16_input,T_DATA T_STRUCT_4::* data_field_4_input)
        :allocator_16(allocator_16_in),set_16(set_16_in),data_field_16(data_field_16_input),data_field_4(data_field_4_input)
    {Run(allocator_4,blocks_4);}
    
    void Run(Allocator_Type_4& allocator_4,const std::pair<const unsigned long*,unsigned>& blocks_4)const{
        Data_Array_Type_4 data_4=allocator_4.Get_Array(data_field_4);
        Const_Data_Array_Type_16 data_16=allocator_16.Get_Const_Array(data_field_16);
        unsigned long offset_16[2][2];
        std_array<unsigned,d> block_size_4=allocator_4.Block_Size();
        std_array<unsigned,d> block_size_16=allocator_16.Block_Size();
        for(SPGrid_Block_Iterator<T_MASK_4> block_iterator(blocks_4);block_iterator.Valid();block_iterator.Next_Block()){
            unsigned long offset_4=block_iterator.Offset();
            offset_16[0][0]=Covert_Block_Offset_4_to_16(offset_4);
            offset_16[0][1]=offset_16[0][0]|(1<<12);
            offset_16[1][0]=offset_16[0][0]|(1<<13);
            offset_16[1][1]=offset_16[1][0]|(1<<12);
            typedef T_DATA (&BLOCK_4)[block_4_xsize][block_4_ysize][block_4_zsize];            
            for(int i=0;i<=1;++i) for(int j=0;j<=1;++j) if(set_16.IsPageActive(offset_16[i][j])){
                unsigned long offset_16_element=offset_16[i][j];
                for(int index_i=0;index_i<block_size_16(0);++index_i)
                for(int index_j=0;index_j<block_size_16(1);++index_j)
                for(int index_k=0;index_k<block_size_16(2);++index_k,offset_16_element+=sizeof(T_DATA)){
                    reinterpret_cast<BLOCK_4>(data_4(offset_4))[index_i][index_j+i*block_size_16(1)][index_k+j*block_size_16(2)]=data_16(offset_16_element);}}}        
    }
};

template<class T_STRUCT_16,class T_STRUCT_4,class T_DATA,int d> class SPGrid_Convertor_4_To_16;
template<class T_STRUCT_16,class T_STRUCT_4,class T_DATA>
class SPGrid_Convertor_4_To_16<T_STRUCT_16,T_STRUCT_4,T_DATA,3>
{
    enum{d=3};
    //This is intend to convert between a 16 channel SPGrid and a 4 channel SPGrid. Do not use it for any other purposes!
    typedef SPGrid_Allocator<T_STRUCT_16,d> Allocator_Type_16;
    typedef typename SPGrid_Allocator<T_STRUCT_16,d>::template Array<T_DATA>::type Data_Array_Type_16;
    typedef typename SPGrid_Allocator<T_STRUCT_16,d>::Array<unsigned>::type Flags_Array_Type_16;
    typedef SPGrid_Set<Flags_Array_Type_16> Set_Type_16;

    typedef SPGrid_Allocator<T_STRUCT_4,d> Allocator_Type_4;
    typedef typename SPGrid_Allocator<T_STRUCT_4,d>::template Array<const T_DATA>::type Const_Data_Array_Type_4;
    typedef typename SPGrid_Allocator<T_STRUCT_4,d>::template Array<T_DATA>::mask T_MASK_4;

    enum{block_4_xsize=1u<<T_MASK_4::block_xbits,
         block_4_ysize=1u<<T_MASK_4::block_ybits,
         block_4_zsize=1u<<T_MASK_4::block_ybits};

    typedef VECTOR<int,d> T_INDEX;

    Allocator_Type_16& allocator_16;
    Set_Type_16& set_16;
    T_DATA T_STRUCT_16::* data_field_16;
    T_DATA T_STRUCT_4::* data_field_4;

public:
    SPGrid_Convertor_4_To_16(Allocator_Type_16& allocator_16_in,Set_Type_16& set_16_in,
                             T_DATA T_STRUCT_16::* data_field_16_input,T_DATA T_STRUCT_4::* data_field_4_input)
        :allocator_16(allocator_16_in),set_16(set_16_in),data_field_16(data_field_16_input),data_field_4(data_field_4_input){}
    
    SPGrid_Convertor_4_To_16(Allocator_Type_4& allocator_4,const std::pair<const unsigned long*,unsigned>&blocks_4,
                             Allocator_Type_16& allocator_16_in,Set_Type_16& set_16_in,
                             T_DATA T_STRUCT_16::* data_field_16_input,T_DATA T_STRUCT_4::* data_field_4_input)
        :allocator_16(allocator_16_in),set_16(set_16_in),data_field_16(data_field_16_input),data_field_4(data_field_4_input)
    {Run(allocator_4,blocks_4);}
    
    void Run(Allocator_Type_4& allocator_4,const std::pair<const unsigned long*,unsigned>&blocks_4)const{
        static_assert(d==3,"Copy from 16 channel of 4 channel SPGrid. Only 3D is supported");
        Const_Data_Array_Type_4 data_4=allocator_4.Get_Const_Array(data_field_4);
        Data_Array_Type_16 data_16=allocator_16.Get_Array(data_field_16);
        unsigned long offset_16[2][2];
        std_array<unsigned,d> block_size_4=allocator_4.Block_Size();
        std_array<unsigned,d> block_size_16=allocator_16.Block_Size();
        for(SPGrid_Block_Iterator<T_MASK_4> block_iterator(blocks_4);block_iterator.Valid();block_iterator.Next_Block()){
            unsigned long offset_4=block_iterator.Offset();
            offset_16[0][0]=Covert_Block_Offset_4_to_16(offset_4);
            offset_16[0][1]=offset_16[0][0]|(1<<12);
            offset_16[1][0]=offset_16[0][0]|(1<<13);
            offset_16[1][1]=offset_16[1][0]|(1<<12);
            typedef const T_DATA (&BLOCK_4)[block_4_xsize][block_4_ysize][block_4_zsize];            
            for(int i=0;i<=1;++i) for(int j=0;j<=1;++j) if(set_16.IsPageActive(offset_16[i][j])){
                unsigned long offset_16_element=offset_16[i][j];
                for(int index_i=0;index_i<block_size_16(0);++index_i)
                for(int index_j=0;index_j<block_size_16(1);++index_j)
                for(int index_k=0;index_k<block_size_16(2);++index_k,offset_16_element+=sizeof(T_DATA)){
                    data_16(offset_16_element)=reinterpret_cast<BLOCK_4>(data_4(offset_4))[index_i][index_j+i*block_size_16(1)][index_k+j*block_size_16(2)];}}}
    }
};

template<class T_STRUCT_16,class T_STRUCT_4,int d> class SPGrid_Flag_Convertor_16_To_4;
template<class T_STRUCT_16,class T_STRUCT_4>
class SPGrid_Flag_Convertor_16_To_4<T_STRUCT_16,T_STRUCT_4,3>
{
    enum{d=3};
    //This is intend to convert between a 16 channel SPGrid and a 4 channel SPGrid. Do not use it for any other purposes!
    typedef SPGrid_Allocator<T_STRUCT_16,d> Allocator_Type_16;
    typedef typename SPGrid_Allocator<T_STRUCT_16,d>::Array<unsigned>::type Flags_Array_Type_16;
    typedef typename SPGrid_Allocator<T_STRUCT_16,d>::Array<const unsigned>::type Const_Flags_Array_Type_16;
    typedef SPGrid_Set<Flags_Array_Type_16> Set_Type_16;

    typedef SPGrid_Allocator<T_STRUCT_4,d> Allocator_Type_4;
    typedef typename SPGrid_Allocator<T_STRUCT_4,d>::Array<unsigned>::type Flags_Array_Type_4;
    typedef typename SPGrid_Allocator<T_STRUCT_4,d>::template Array<unsigned>::mask T_MASK_4;

    enum{block_4_xsize=1u<<T_MASK_4::block_xbits,
         block_4_ysize=1u<<T_MASK_4::block_ybits,
         block_4_zsize=1u<<T_MASK_4::block_ybits};

    typedef VECTOR<int,d> T_INDEX;

    Allocator_Type_16& allocator_16;
    Set_Type_16& set_16;
    unsigned T_STRUCT_16::* flag_field_16;
    unsigned T_STRUCT_4::* flag_field_4;

public:
    SPGrid_Flag_Convertor_16_To_4(Allocator_Type_16& allocator_16_in,Set_Type_16& set_16_in,
                                  unsigned T_STRUCT_16::* flag_field_16_input,unsigned T_STRUCT_4::* flag_field_4_input)
        :allocator_16(allocator_16_in),set_16(set_16_in),flag_field_16(flag_field_16_input),flag_field_4(flag_field_4_input){}

    SPGrid_Flag_Convertor_16_To_4(Allocator_Type_4& allocator_4,const std::pair<const unsigned long*,unsigned>&blocks_4,
                                  Allocator_Type_16& allocator_16_in,Set_Type_16& set_16_in,
                                  unsigned T_STRUCT_16::* flag_field_16_input,unsigned T_STRUCT_4::* flag_field_4_input)
        :allocator_16(allocator_16_in),set_16(set_16_in),flag_field_16(flag_field_16_input),flag_field_4(flag_field_4_input)
    {Run(allocator_4,blocks_4);}
    
    void Run(Allocator_Type_4& allocator_4,const std::pair<const unsigned long*,unsigned>& blocks_4)const{
        Flags_Array_Type_4 data_4=allocator_4.Get_Array(flag_field_4);
        Const_Flags_Array_Type_16 data_16=allocator_16.Get_Const_Array(flag_field_16);
        unsigned long offset_16[2][2];
        std_array<unsigned,d> block_size_4=allocator_4.Block_Size();
        std_array<unsigned,d> block_size_16=allocator_16.Block_Size();
        for(SPGrid_Block_Iterator<T_MASK_4> block_iterator(blocks_4);block_iterator.Valid();block_iterator.Next_Block()){
            unsigned long offset_4=block_iterator.Offset();
            offset_16[0][0]=Covert_Block_Offset_4_to_16(offset_4);
            offset_16[0][1]=offset_16[0][0]|(1<<12);
            offset_16[1][0]=offset_16[0][0]|(1<<13);
            offset_16[1][1]=offset_16[1][0]|(1<<12);
            typedef unsigned (&BLOCK_4)[block_4_xsize][block_4_ysize][block_4_zsize];            
            for(int i=0;i<=1;++i) for(int j=0;j<=1;++j) if(set_16.IsPageActive(offset_16[i][j])){
                unsigned long offset_16_element=offset_16[i][j];
                for(int index_i=0;index_i<block_size_16(0);++index_i)
                for(int index_j=0;index_j<block_size_16(1);++index_j)
                for(int index_k=0;index_k<block_size_16(2);++index_k,offset_16_element+=sizeof(unsigned)){
                    unsigned& flag_4=reinterpret_cast<BLOCK_4>(data_4(offset_4))[index_i][index_j+i*block_size_16(1)][index_k+j*block_size_16(2)];
                    flag_4=0x0u;
                    if(data_16(offset_16_element)&SPGrid_Cell_Type_Dirichlet)
                        flag_4=SPGrid_Solver_Cell_Type_Dirichlet;
                    else if(data_16(offset_16_element)&SPGrid_Cell_Type_Active)
                        flag_4=SPGrid_Solver_Cell_Type_Active;}}}        
    }
};
}
#endif
