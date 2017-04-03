//#####################################################################
//Removed For Anonymity: Copyright Authors of Submission pap173s1
//Removed For Anonymity: Copyright Authors of Submission pap173s1
//#####################################################################
// Subroutine SPGrid_Computations::Restriction
//#####################################################################
#ifndef __SPGrid_Restriction_h__
#define __SPGrid_Restriction_h__

#include <Common_Tools/Math_Tools/RANGE_ITERATOR.h>
#include <SPGrid/Core/SPGrid_Allocator.h>
#include <SPGrid_Fluids/Solvers/Domain_Decomposition/SPGRID_MULTIGRID_FLAGS.h>

namespace SPGrid_Computations{

using namespace SPGrid;
using namespace PhysBAM;

template<class T_STRUCT,class T_DATA,class T_FLAGS,int d>class Restriction;
template<class T_STRUCT,class T_DATA,class T_FLAGS>
class Restriction<T_STRUCT,T_DATA,T_FLAGS,3>
{
    enum{d=3};
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<const T_DATA>::type Const_data_array_type;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<T_DATA>::type Data_array_type;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<const T_FLAGS>::type Const_flag_array_type;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<T_FLAGS>::mask Flag_array_mask;
    typedef VECTOR<int,d> T_INDEX;

    T_DATA T_STRUCT::* r_field;
    T_DATA T_STRUCT::* b_field;
    T_FLAGS T_STRUCT::* flags_field;
    const SPGrid_Allocator<T_STRUCT,d>& fine_allocator;

public:
    Restriction(const SPGrid_Allocator<T_STRUCT,d>& fine_allocator_input,
                T_DATA T_STRUCT::* b_field_input,T_DATA T_STRUCT::* r_field_input,T_FLAGS T_STRUCT::* flags_field_input)
        :fine_allocator(fine_allocator_input),b_field(b_field_input),r_field(r_field_input),flags_field(flags_field_input)
    {}
    
     Restriction(SPGrid_Allocator<T_STRUCT,d>& coarse_allocator,const std::pair<const unsigned long*,unsigned>& coarse_blocks,
                 const SPGrid_Allocator<T_STRUCT,d>& fine_allocator_input,
                 T_DATA T_STRUCT::* b_field_input,T_DATA T_STRUCT::* r_field_input,T_FLAGS T_STRUCT::* flags_field_input)
         :fine_allocator(fine_allocator_input),b_field(b_field_input),r_field(r_field_input),flags_field(flags_field_input)
    {Run(coarse_allocator,coarse_blocks);}
    
    void Run(SPGrid_Allocator<T_STRUCT,d>& coarse_allocator,const std::pair<const unsigned long*,unsigned>& coarse_blocks) const
    {
        float fine_stencil_scales[4][4][4];
        float stencil_1D[4]={1./8.,3./8.,3./8.,1./8.}; 
        for(int i=0;i<=3;++i){
            for(int j=0;j<=3;++j){
                for(int k=0;k<=3;++k){
                    fine_stencil_scales[i][j][k]=stencil_1D[i]*stencil_1D[j]*stencil_1D[k]; 
                }
            }
        }

        Data_array_type b=coarse_allocator.Get_Array(b_field);
        Const_data_array_type r=fine_allocator.Get_Const_Array(r_field);
        Const_flag_array_type flags=coarse_allocator.Get_Const_Array(flags_field);    
        //Const_flag_array_type fine_flags=fine_allocator.Get_Const_Array(flags_field);    
        enum{
            block_xsize = 1u << Const_flag_array_type::MASK::block_xbits,
            block_ysize = 1u << Const_flag_array_type::MASK::block_ybits,
            block_zsize = 1u << Const_flag_array_type::MASK::block_zbits  
        };
        typedef const unsigned (&const_block_flag)[block_xsize][block_ysize][block_zsize];    
        typedef unsigned (&block_flag)[block_xsize][block_ysize][block_zsize];    
        typedef const T_DATA (&const_block_data)[block_xsize][block_ysize][block_zsize];    
        typedef T_DATA (&block_data)[block_xsize][block_ysize][block_zsize];    

        for(SPGrid_Block_Iterator<Flag_array_mask> iterator(coarse_blocks);iterator.Valid();iterator.Next_Block()){
            unsigned long coarse_offset=iterator.Offset();
            T_INDEX base_index=iterator.Index().template Cast<T_INDEX>();
            std_array<int,d> coarse_coor = Const_flag_array_type::MASK::LinearToCoord(coarse_offset);
            for(RANGE_ITERATOR<d> iterator(RANGE<T_INDEX>(base_index,base_index+coarse_allocator.Block_Size().template Cast<T_INDEX>()-1));
                iterator.Valid();
                iterator.Next(),coarse_offset+=sizeof(unsigned)){
                unsigned coarse_flag = flags(coarse_offset);
                if(coarse_flag & SPGrid_Solver_Cell_Type_Active){
                    T_INDEX coarse_index=iterator.Index();
                    std_array<int,d> fine_coor = std_array<int,d>(coarse_index*2-1);
                    const unsigned long fine_offset=Const_flag_array_type::MASK::Linear_Offset(fine_coor);
                    PHYSBAM_ASSERT(fine_offset==Const_flag_array_type::MASK::UpsampleOffset(coarse_offset));
                    
                    unsigned long fine_offset_array[4][4][4];
                    fine_offset_array[0][0][0] = Const_flag_array_type::MASK::Packed_Offset<-1,-1,-1>(fine_offset);
                    fine_offset_array[0][0][1] = Const_flag_array_type::MASK::Packed_Offset<-1,-1,0>(fine_offset);
                    fine_offset_array[0][0][2] = Const_flag_array_type::MASK::Packed_Offset<-1,-1,+1>(fine_offset);
                    fine_offset_array[0][0][3] = Const_flag_array_type::MASK::Packed_Offset<-1,-1,+2>(fine_offset);

                    fine_offset_array[0][1][0] = Const_flag_array_type::MASK::Packed_Offset<-1,0,-1>(fine_offset);
                    fine_offset_array[0][1][1] = Const_flag_array_type::MASK::Packed_Offset<-1,0,0>(fine_offset);
                    fine_offset_array[0][1][2] = Const_flag_array_type::MASK::Packed_Offset<-1,0,+1>(fine_offset);
                    fine_offset_array[0][1][3] = Const_flag_array_type::MASK::Packed_Offset<-1,0,+2>(fine_offset);

                    fine_offset_array[0][2][0] = Const_flag_array_type::MASK::Packed_Offset<-1,1,-1>(fine_offset);
                    fine_offset_array[0][2][1] = Const_flag_array_type::MASK::Packed_Offset<-1,1,0>(fine_offset);
                    fine_offset_array[0][2][2] = Const_flag_array_type::MASK::Packed_Offset<-1,1,+1>(fine_offset);
                    fine_offset_array[0][2][3] = Const_flag_array_type::MASK::Packed_Offset<-1,1,+2>(fine_offset);

                    fine_offset_array[0][3][0] = Const_flag_array_type::MASK::Packed_Offset<-1,2,-1>(fine_offset);
                    fine_offset_array[0][3][1] = Const_flag_array_type::MASK::Packed_Offset<-1,2,0>(fine_offset);
                    fine_offset_array[0][3][2] = Const_flag_array_type::MASK::Packed_Offset<-1,2,+1>(fine_offset);
                    fine_offset_array[0][3][3] = Const_flag_array_type::MASK::Packed_Offset<-1,2,+2>(fine_offset);

                    fine_offset_array[1][0][0] = Const_flag_array_type::MASK::Packed_Offset<0,-1,-1>(fine_offset);
                    fine_offset_array[1][0][1] = Const_flag_array_type::MASK::Packed_Offset<0,-1,0>(fine_offset);
                    fine_offset_array[1][0][2] = Const_flag_array_type::MASK::Packed_Offset<0,-1,+1>(fine_offset);
                    fine_offset_array[1][0][3] = Const_flag_array_type::MASK::Packed_Offset<0,-1,+2>(fine_offset);

                    fine_offset_array[1][1][0] = Const_flag_array_type::MASK::Packed_Offset<0,0,-1>(fine_offset);
                    fine_offset_array[1][1][1] = Const_flag_array_type::MASK::Packed_Offset<0,0,0>(fine_offset);
                    fine_offset_array[1][1][2] = Const_flag_array_type::MASK::Packed_Offset<0,0,+1>(fine_offset);
                    fine_offset_array[1][1][3] = Const_flag_array_type::MASK::Packed_Offset<0,0,+2>(fine_offset);

                    fine_offset_array[1][2][0] = Const_flag_array_type::MASK::Packed_Offset<0,1,-1>(fine_offset);
                    fine_offset_array[1][2][1] = Const_flag_array_type::MASK::Packed_Offset<0,1,0>(fine_offset);
                    fine_offset_array[1][2][2] = Const_flag_array_type::MASK::Packed_Offset<0,1,+1>(fine_offset);
                    fine_offset_array[1][2][3] = Const_flag_array_type::MASK::Packed_Offset<0,1,+2>(fine_offset);

                    fine_offset_array[1][3][0] = Const_flag_array_type::MASK::Packed_Offset<0,2,-1>(fine_offset);
                    fine_offset_array[1][3][1] = Const_flag_array_type::MASK::Packed_Offset<0,2,0>(fine_offset);
                    fine_offset_array[1][3][2] = Const_flag_array_type::MASK::Packed_Offset<0,2,+1>(fine_offset);
                    fine_offset_array[1][3][3] = Const_flag_array_type::MASK::Packed_Offset<0,2,+2>(fine_offset);

                    fine_offset_array[2][0][0] = Const_flag_array_type::MASK::Packed_Offset<1,-1,-1>(fine_offset);
                    fine_offset_array[2][0][1] = Const_flag_array_type::MASK::Packed_Offset<1,-1,0>(fine_offset);
                    fine_offset_array[2][0][2] = Const_flag_array_type::MASK::Packed_Offset<1,-1,+1>(fine_offset);
                    fine_offset_array[2][0][3] = Const_flag_array_type::MASK::Packed_Offset<1,-1,+2>(fine_offset);

                    fine_offset_array[2][1][0] = Const_flag_array_type::MASK::Packed_Offset<1,0,-1>(fine_offset);
                    fine_offset_array[2][1][1] = Const_flag_array_type::MASK::Packed_Offset<1,0,0>(fine_offset);
                    fine_offset_array[2][1][2] = Const_flag_array_type::MASK::Packed_Offset<1,0,+1>(fine_offset);
                    fine_offset_array[2][1][3] = Const_flag_array_type::MASK::Packed_Offset<1,0,+2>(fine_offset);

                    fine_offset_array[2][2][0] = Const_flag_array_type::MASK::Packed_Offset<1,1,-1>(fine_offset);
                    fine_offset_array[2][2][1] = Const_flag_array_type::MASK::Packed_Offset<1,1,0>(fine_offset);
                    fine_offset_array[2][2][2] = Const_flag_array_type::MASK::Packed_Offset<1,1,+1>(fine_offset);
                    fine_offset_array[2][2][3] = Const_flag_array_type::MASK::Packed_Offset<1,1,+2>(fine_offset);

                    fine_offset_array[2][3][0] = Const_flag_array_type::MASK::Packed_Offset<1,2,-1>(fine_offset);
                    fine_offset_array[2][3][1] = Const_flag_array_type::MASK::Packed_Offset<1,2,0>(fine_offset);
                    fine_offset_array[2][3][2] = Const_flag_array_type::MASK::Packed_Offset<1,2,+1>(fine_offset);
                    fine_offset_array[2][3][3] = Const_flag_array_type::MASK::Packed_Offset<1,2,+2>(fine_offset);

                    fine_offset_array[3][0][0] = Const_flag_array_type::MASK::Packed_Offset<2,-1,-1>(fine_offset);
                    fine_offset_array[3][0][1] = Const_flag_array_type::MASK::Packed_Offset<2,-1,0>(fine_offset);
                    fine_offset_array[3][0][2] = Const_flag_array_type::MASK::Packed_Offset<2,-1,+1>(fine_offset);
                    fine_offset_array[3][0][3] = Const_flag_array_type::MASK::Packed_Offset<2,-1,+2>(fine_offset);

                    fine_offset_array[3][1][0] = Const_flag_array_type::MASK::Packed_Offset<2,0,-1>(fine_offset);
                    fine_offset_array[3][1][1] = Const_flag_array_type::MASK::Packed_Offset<2,0,0>(fine_offset);
                    fine_offset_array[3][1][2] = Const_flag_array_type::MASK::Packed_Offset<2,0,+1>(fine_offset);
                    fine_offset_array[3][1][3] = Const_flag_array_type::MASK::Packed_Offset<2,0,+2>(fine_offset);

                    fine_offset_array[3][2][0] = Const_flag_array_type::MASK::Packed_Offset<2,1,-1>(fine_offset);
                    fine_offset_array[3][2][1] = Const_flag_array_type::MASK::Packed_Offset<2,1,0>(fine_offset);
                    fine_offset_array[3][2][2] = Const_flag_array_type::MASK::Packed_Offset<2,1,+1>(fine_offset);
                    fine_offset_array[3][2][3] = Const_flag_array_type::MASK::Packed_Offset<2,1,+2>(fine_offset);

                    fine_offset_array[3][3][0] = Const_flag_array_type::MASK::Packed_Offset<2,2,-1>(fine_offset);
                    fine_offset_array[3][3][1] = Const_flag_array_type::MASK::Packed_Offset<2,2,0>(fine_offset);
                    fine_offset_array[3][3][2] = Const_flag_array_type::MASK::Packed_Offset<2,2,+1>(fine_offset);
                    fine_offset_array[3][3][3] = Const_flag_array_type::MASK::Packed_Offset<2,2,+2>(fine_offset);

                    b(coarse_offset) = 0.f;      
                    float test = 0.f;                  
                    for(int i2=0;i2<4;++i2)
                    for(int j2=0;j2<4;++j2)
                    for(int k2=0;k2<4;++k2){
                        // if(!(fine_flags(fine_offset_array[i2][j2][k2])&SPGrid_Solver_Cell_Type_Active) &&
                        //    !(fine_flags(fine_offset_array[i2][j2][k2])&SPGrid_Solver_Cell_Type_Interface)){
                        //     if(r(fine_offset_array[i2][j2][k2]) != T_DATA(0))
                        //         std::cout<<"errrrr!!!!!"<<std::endl;                            
                        //  }
                        test += fine_stencil_scales[i2][j2][k2];
                        b(coarse_offset) += fine_stencil_scales[i2][j2][k2]*r(fine_offset_array[i2][j2][k2]);                        
                    }
                    //std::cout<<"Test!"<<test<<std::endl;
                    PHYSBAM_ASSERT(test == 1.f);
                }
            }
        }
    } 
};

/////////////////////////////////////////////////////////////////////////
// 2D
/////////////////////////////////////////////////////////////////////////
template<class T_STRUCT,class T_DATA,class T_FLAGS>
class Restriction<T_STRUCT,T_DATA,T_FLAGS,2>
{
    enum{d=2};
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<const T_DATA>::type Const_data_array_type;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<T_DATA>::type Data_array_type;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<const T_FLAGS>::type Const_flag_array_type;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<T_FLAGS>::mask Flag_array_mask;
    typedef VECTOR<int,d> T_INDEX;

    T_DATA T_STRUCT::* r_field;
    T_DATA T_STRUCT::* b_field;
    T_FLAGS T_STRUCT::* flags_field;
    const SPGrid_Allocator<T_STRUCT,d>& fine_allocator;

public:
    Restriction(
                const SPGrid_Allocator<T_STRUCT,d>& fine_allocator_input,
                T_DATA T_STRUCT::* b_field_input,T_DATA T_STRUCT::* r_field_input,T_FLAGS T_STRUCT::* flags_field_input)
        :fine_allocator(fine_allocator_input),b_field(b_field_input),r_field(r_field_input),flags_field(flags_field_input)
    {}
    
     Restriction(SPGrid_Allocator<T_STRUCT,d>& coarse_allocator,const std::pair<const unsigned long*,unsigned>& coarse_blocks,
                 const SPGrid_Allocator<T_STRUCT,d>& fine_allocator_input,
                 T_DATA T_STRUCT::* b_field_input,T_DATA T_STRUCT::* r_field_input,T_FLAGS T_STRUCT::* flags_field_input)
         :fine_allocator(fine_allocator_input),b_field(b_field_input),r_field(r_field_input),flags_field(flags_field_input)
    {Run(coarse_allocator,coarse_blocks);}
    
    void Run(SPGrid_Allocator<T_STRUCT,d>& coarse_allocator,const std::pair<const unsigned long*,unsigned>& coarse_blocks) const
    {
        float fine_stencil_scales[4][4];
        float stencil_1D[4]={1./8.,3./8.,3./8.,1./8.}; 
        for(int i=0;i<=3;++i)
            for(int j=0;j<=3;++j)
                fine_stencil_scales[i][j]=stencil_1D[i]*stencil_1D[j]; 
        
        Data_array_type b=coarse_allocator.Get_Array(b_field);
        Const_data_array_type r=fine_allocator.Get_Const_Array(r_field);
        Const_flag_array_type flags=coarse_allocator.Get_Const_Array(flags_field);    
        enum{
            block_xsize = 1u << Const_flag_array_type::MASK::block_xbits,
            block_ysize = 1u << Const_flag_array_type::MASK::block_ybits
        };
        typedef const unsigned (&const_block_flag)[block_xsize][block_ysize];    
        typedef unsigned (&block_flag)[block_xsize][block_ysize];    
        typedef const T_DATA (&const_block_data)[block_xsize][block_ysize];    
        typedef T_DATA (&block_data)[block_xsize][block_ysize];    

        for(SPGrid_Block_Iterator<Flag_array_mask> iterator(coarse_blocks);iterator.Valid();iterator.Next_Block()){
            unsigned long coarse_offset=iterator.Offset();
            T_INDEX base_index=iterator.Index().template Cast<T_INDEX>();
            std_array<int,d> coarse_coor = Const_flag_array_type::MASK::LinearToCoord(coarse_offset);
            for(RANGE_ITERATOR<d> iterator(RANGE<T_INDEX>(base_index,base_index+coarse_allocator.Block_Size().template Cast<T_INDEX>()-1));
                iterator.Valid();
                iterator.Next(),coarse_offset+=sizeof(unsigned)){
                unsigned coarse_flag = flags(coarse_offset);
                if(coarse_flag & SPGrid_Solver_Cell_Type_Active){
                    T_INDEX coarse_index=iterator.Index();
                    std_array<int,d> fine_coor = std_array<int,d>(coarse_index*2-1);
                    const unsigned long fine_offset=Const_flag_array_type::MASK::Linear_Offset(fine_coor);
                    PHYSBAM_ASSERT(fine_offset==Const_flag_array_type::MASK::UpsampleOffset(coarse_offset));
                    
                    unsigned long fine_offset_array[4][4];
                    fine_offset_array[0][0] = Const_flag_array_type::MASK::Packed_Offset<-1,-1>(fine_offset);
                    fine_offset_array[0][1] = Const_flag_array_type::MASK::Packed_Offset<-1,0>(fine_offset);
                    fine_offset_array[0][2] = Const_flag_array_type::MASK::Packed_Offset<-1,+1>(fine_offset);
                    fine_offset_array[0][3] = Const_flag_array_type::MASK::Packed_Offset<-1,+2>(fine_offset);

                    fine_offset_array[1][0] = Const_flag_array_type::MASK::Packed_Offset<0,-1>(fine_offset);
                    fine_offset_array[1][1] = Const_flag_array_type::MASK::Packed_Offset<0,0>(fine_offset);
                    fine_offset_array[1][2] = Const_flag_array_type::MASK::Packed_Offset<0,+1>(fine_offset);
                    fine_offset_array[1][3] = Const_flag_array_type::MASK::Packed_Offset<0,+2>(fine_offset);

                    fine_offset_array[2][0] = Const_flag_array_type::MASK::Packed_Offset<1,-1>(fine_offset);
                    fine_offset_array[2][1] = Const_flag_array_type::MASK::Packed_Offset<1,0>(fine_offset);
                    fine_offset_array[2][2] = Const_flag_array_type::MASK::Packed_Offset<1,+1>(fine_offset);
                    fine_offset_array[2][3] = Const_flag_array_type::MASK::Packed_Offset<1,+2>(fine_offset);

                    fine_offset_array[3][0] = Const_flag_array_type::MASK::Packed_Offset<2,-1>(fine_offset);
                    fine_offset_array[3][1] = Const_flag_array_type::MASK::Packed_Offset<2,0>(fine_offset);
                    fine_offset_array[3][2] = Const_flag_array_type::MASK::Packed_Offset<2,+1>(fine_offset);
                    fine_offset_array[3][3] = Const_flag_array_type::MASK::Packed_Offset<2,+2>(fine_offset);

                    b(coarse_offset) = 0.f;                        
                    for(int i2=0;i2<4;++i2)
                    for(int j2=0;j2<4;++j2){
                        b(coarse_offset) += fine_stencil_scales[i2][j2]*r(fine_offset_array[i2][j2]); 
                    }
                }
            }
        }
    } 
};
//#####################################################################
}
#endif
