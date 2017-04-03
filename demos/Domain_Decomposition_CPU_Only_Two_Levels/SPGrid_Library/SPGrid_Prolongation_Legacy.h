//#####################################################################
//Removed For Anonymity: Copyright Authors of Submission pap173s1
//Removed For Anonymity: Copyright Authors of Submission pap173s1
//#####################################################################
// Subroutine SPGrid_Computations::Prolongation
//#####################################################################
#ifndef __SPGrid_Prolongation_h__
#define __SPGrid_Prolongation_h__

#include <SPGrid/Core/SPGrid_Allocator.h>
#include <SPGrid_Fluids/Solvers/Domain_Decomposition/SPGRID_MULTIGRID_FLAGS.h>

namespace SPGrid_Computations{

using namespace SPGrid;
using namespace PhysBAM;

template<class T_STRUCT,class T_DATA,class T_FLAGS,int d> class Prolongation;
template<class T_STRUCT,class T_DATA,class T_FLAGS>
class Prolongation<T_STRUCT,T_DATA,T_FLAGS,3>
{
    enum{d = 3};
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<const T_DATA>::type Const_data_array_type;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<T_DATA>::type Data_array_type;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<const T_FLAGS>::type Const_flag_array_type;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<T_FLAGS>::mask Flag_array_mask;
    typedef VECTOR<int,d> T_INDEX;

    T_DATA T_STRUCT::* u_field;
    T_FLAGS T_STRUCT::* flags_field;
    const SPGrid_Allocator<T_STRUCT,d>& coarse_allocator;

public:
    Prolongation(
                const SPGrid_Allocator<T_STRUCT,d>& coarse_allocator_input,
                T_DATA T_STRUCT::* u_field_input,T_FLAGS T_STRUCT::* flags_field_input)
        :coarse_allocator(coarse_allocator_input),u_field(u_field_input),flags_field(flags_field_input)
    {}
    
     Prolongation(SPGrid_Allocator<T_STRUCT,d>& fine_allocator,const std::pair<const unsigned long*,unsigned>& fine_blocks,
                  const SPGrid_Allocator<T_STRUCT,d>& coarse_allocator_input,
                  T_DATA T_STRUCT::* u_field_input,T_FLAGS T_STRUCT::* flags_field_input)
         :coarse_allocator(coarse_allocator_input),u_field(u_field_input),flags_field(flags_field_input)
    {Run(fine_allocator,fine_blocks);}
    
    void Run(SPGrid_Allocator<T_STRUCT,d>& fine_allocator,const std::pair<const unsigned long*,unsigned>& fine_blocks) const
    {
        float fine_stencil_scales[2][2][2];//first three indicates the fine cell's position in the stencil. the last three indicates the coarse cell's position.
        float stencil_1D[2]={3./4.f,1./4.f};
        for(int i=0;i<=1;++i)
            for(int j=0;j<=1;++j)
                for(int k=0;k<=1;++k)
                    fine_stencil_scales[i][j][k]=stencil_1D[i]*stencil_1D[j]*stencil_1D[k] * 4.0f; 

        Const_data_array_type u_coarse=coarse_allocator.Get_Const_Array(u_field);
        Data_array_type u_fine=fine_allocator.Get_Array(u_field);
        Const_flag_array_type flags=fine_allocator.Get_Const_Array(flags_field);    
        
        for(SPGrid_Block_Iterator<Flag_array_mask> fine_iterator(fine_blocks);fine_iterator.Valid();fine_iterator.Next_Block()){
            unsigned long fine_offset=fine_iterator.Offset();
            std_array<int,d> fine_coor = Const_flag_array_type::MASK::LinearToCoord(fine_offset);
            T_INDEX base_index=fine_iterator.Index().template Cast<T_INDEX>();
            for(RANGE_ITERATOR<d> iterator(RANGE<T_INDEX>(base_index,base_index+fine_allocator.Block_Size().template Cast<T_INDEX>()-1));
                iterator.Valid();
                iterator.Next(),fine_offset+=sizeof(unsigned)){
                unsigned fine_flag = flags(fine_offset);
                if(fine_flag & SPGrid_Solver_Cell_Type_Active){
                    T_INDEX fine_index=iterator.Index();
                    std_array<int,d> coarse_coor = std_array<int,d>((fine_index+1)/2);;
                    const unsigned long coarse_offset=Const_flag_array_type::MASK::DownsampleOffset(fine_offset);
                    //const unsigned long tmp = Flag_array_mask::Linear_Offset(coarse_coor);
                    //PHYSBAM_ASSERT(tmp == coarse_offset);
                    std_array<int,d> fine_coor_ref = std_array<int,d>(coarse_coor.template Cast<T_INDEX>()*2-1);
                    const int parity_x=1-(fine_index(1)-fine_coor_ref(0));
                    const int parity_y=1-(fine_index(2)-fine_coor_ref(1));
                    const int parity_z=1-(fine_index(3)-fine_coor_ref(2));

                    unsigned long coarse_offset_array[3][3][3];
                    coarse_offset_array[0][0][0] = Const_flag_array_type::MASK::Packed_Offset<-1,-1,-1>(coarse_offset);
                    coarse_offset_array[0][0][1] = Const_flag_array_type::MASK::Packed_Offset<-1,-1,0>(coarse_offset);
                    coarse_offset_array[0][0][2] = Const_flag_array_type::MASK::Packed_Offset<-1,-1,1>(coarse_offset);

                    coarse_offset_array[0][1][0] = Const_flag_array_type::MASK::Packed_Offset<-1,0,-1>(coarse_offset);
                    coarse_offset_array[0][1][1] = Const_flag_array_type::MASK::Packed_Offset<-1,0,0>(coarse_offset);
                    coarse_offset_array[0][1][2] = Const_flag_array_type::MASK::Packed_Offset<-1,0,1>(coarse_offset);

                    coarse_offset_array[0][2][0] = Const_flag_array_type::MASK::Packed_Offset<-1,1,-1>(coarse_offset);
                    coarse_offset_array[0][2][1] = Const_flag_array_type::MASK::Packed_Offset<-1,1,0>(coarse_offset);
                    coarse_offset_array[0][2][2] = Const_flag_array_type::MASK::Packed_Offset<-1,1,1>(coarse_offset);

                    coarse_offset_array[1][0][0] = Const_flag_array_type::MASK::Packed_Offset<0,-1,-1>(coarse_offset);
                    coarse_offset_array[1][0][1] = Const_flag_array_type::MASK::Packed_Offset<0,-1,0>(coarse_offset);
                    coarse_offset_array[1][0][2] = Const_flag_array_type::MASK::Packed_Offset<0,-1,1>(coarse_offset);

                    coarse_offset_array[1][1][0] = Const_flag_array_type::MASK::Packed_Offset<0,0,-1>(coarse_offset);
                    coarse_offset_array[1][1][1] = Const_flag_array_type::MASK::Packed_Offset<0,0,0>(coarse_offset);
                    coarse_offset_array[1][1][2] = Const_flag_array_type::MASK::Packed_Offset<0,0,1>(coarse_offset);

                    coarse_offset_array[1][2][0] = Const_flag_array_type::MASK::Packed_Offset<0,1,-1>(coarse_offset);
                    coarse_offset_array[1][2][1] = Const_flag_array_type::MASK::Packed_Offset<0,1,0>(coarse_offset);
                    coarse_offset_array[1][2][2] = Const_flag_array_type::MASK::Packed_Offset<0,1,1>(coarse_offset);

                    coarse_offset_array[2][0][0] = Const_flag_array_type::MASK::Packed_Offset<1,-1,-1>(coarse_offset);
                    coarse_offset_array[2][0][1] = Const_flag_array_type::MASK::Packed_Offset<1,-1,0>(coarse_offset);
                    coarse_offset_array[2][0][2] = Const_flag_array_type::MASK::Packed_Offset<1,-1,1>(coarse_offset);

                    coarse_offset_array[2][1][0] = Const_flag_array_type::MASK::Packed_Offset<1,0,-1>(coarse_offset);
                    coarse_offset_array[2][1][1] = Const_flag_array_type::MASK::Packed_Offset<1,0,0>(coarse_offset);
                    coarse_offset_array[2][1][2] = Const_flag_array_type::MASK::Packed_Offset<1,0,1>(coarse_offset);

                    coarse_offset_array[2][2][0] = Const_flag_array_type::MASK::Packed_Offset<1,1,-1>(coarse_offset);
                    coarse_offset_array[2][2][1] = Const_flag_array_type::MASK::Packed_Offset<1,1,0>(coarse_offset);
                    coarse_offset_array[2][2][2] = Const_flag_array_type::MASK::Packed_Offset<1,1,1>(coarse_offset);

                    float test = 0.f;
                    for(int i=0;i<2;++i)
                    for(int j=0;j<2;++j)
                    for(int k=0;k<2;++k){
                        test += fine_stencil_scales[i^parity_x][j^parity_y][k^parity_z];
                        u_fine(fine_offset)+=fine_stencil_scales[i^parity_x][j^parity_y][k^parity_z]*
                                                   u_coarse(coarse_offset_array[1+i-parity_x][1+j-parity_y][1+k-parity_z]);
                    }
                    //std::cout<<"Test!"<<test<<std::endl;
                    PHYSBAM_ASSERT(test == 4.0f);
                }
            }
        }
    }
};

///////////////////////////////////////////////////////////////////////////////////
// 2D
///////////////////////////////////////////////////////////////////////////////////
template<class T_STRUCT,class T_DATA,class T_FLAGS>
class Prolongation<T_STRUCT,T_DATA,T_FLAGS,2>
{
    enum{d = 2};
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<const T_DATA>::type Const_data_array_type;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<T_DATA>::type Data_array_type;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<const T_FLAGS>::type Const_flag_array_type;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<T_FLAGS>::mask Flag_array_mask;
    typedef VECTOR<int,d> T_INDEX;

    T_DATA T_STRUCT::* u_field;
    T_FLAGS T_STRUCT::* flags_field;
    const SPGrid_Allocator<T_STRUCT,d>& coarse_allocator;

public:
    Prolongation(
                const SPGrid_Allocator<T_STRUCT,d>& coarse_allocator_input,
                T_DATA T_STRUCT::* u_field_input,T_FLAGS T_STRUCT::* flags_field_input)
        :coarse_allocator(coarse_allocator_input),u_field(u_field_input),flags_field(flags_field_input)
    {}
    
     Prolongation(SPGrid_Allocator<T_STRUCT,d>& fine_allocator,const std::pair<const unsigned long*,unsigned>& fine_blocks,
                  const SPGrid_Allocator<T_STRUCT,d>& coarse_allocator_input,
                  T_DATA T_STRUCT::* u_field_input,T_FLAGS T_STRUCT::* flags_field_input)
         :coarse_allocator(coarse_allocator_input),u_field(u_field_input),flags_field(flags_field_input)
    {Run(fine_allocator,fine_blocks);}
    
    void Run(SPGrid_Allocator<T_STRUCT,d>& fine_allocator,const std::pair<const unsigned long*,unsigned>& fine_blocks) const
    {
        float fine_stencil_scales[2][2];//first three indicates the fine cell's position in the stencil. the last three indicates the coarse cell's position.
        float stencil_1D[2]={3./4.f,1./4.f};
        for(int i=0;i<=1;++i)
            for(int j=0;j<=1;++j)
                fine_stencil_scales[i][j]=stencil_1D[i]*stencil_1D[j] * 4.f; 

        Const_data_array_type u_coarse=coarse_allocator.Get_Const_Array(u_field);
        Data_array_type u_fine=fine_allocator.Get_Array(u_field);
        Const_flag_array_type flags=fine_allocator.Get_Const_Array(flags_field);    
        
        for(SPGrid_Block_Iterator<Flag_array_mask> fine_iterator(fine_blocks);fine_iterator.Valid();fine_iterator.Next_Block()){
            unsigned long fine_offset=fine_iterator.Offset();
            std_array<int,d> fine_coor = Const_flag_array_type::MASK::LinearToCoord(fine_offset);
            T_INDEX base_index=fine_iterator.Index().template Cast<T_INDEX>();
            for(RANGE_ITERATOR<d> iterator(RANGE<T_INDEX>(base_index,base_index+fine_allocator.Block_Size().template Cast<T_INDEX>()-1));
                iterator.Valid();
                iterator.Next(),fine_offset+=sizeof(unsigned)){
                unsigned fine_flag = flags(fine_offset);
                if(fine_flag & SPGrid_Solver_Cell_Type_Active){
                    T_INDEX fine_index=iterator.Index();
                    std_array<int,d> coarse_coor = std_array<int,d>((fine_index+1)/2);;
                    const unsigned long coarse_offset=Const_flag_array_type::MASK::DownsampleOffset(fine_offset);
                    //const unsigned long tmp = Flag_array_mask::Linear_Offset(coarse_coor);
                    //PHYSBAM_ASSERT(tmp == coarse_offset);
                    std_array<int,d> fine_coor_ref = std_array<int,d>(coarse_coor.template Cast<T_INDEX>()*2-1);
                    const int parity_x=1-(fine_index(1)-fine_coor_ref(0));
                    const int parity_y=1-(fine_index(2)-fine_coor_ref(1));

                    unsigned long coarse_offset_array[3][3];
                    coarse_offset_array[0][0] = Const_flag_array_type::MASK::Packed_Offset<-1,-1>(coarse_offset);
                    coarse_offset_array[0][1] = Const_flag_array_type::MASK::Packed_Offset<-1,0>(coarse_offset);
                    coarse_offset_array[0][2] = Const_flag_array_type::MASK::Packed_Offset<-1,1>(coarse_offset);

                    coarse_offset_array[1][0] = Const_flag_array_type::MASK::Packed_Offset<0,-1>(coarse_offset);
                    coarse_offset_array[1][1] = Const_flag_array_type::MASK::Packed_Offset<0,0>(coarse_offset);
                    coarse_offset_array[1][2] = Const_flag_array_type::MASK::Packed_Offset<0,1>(coarse_offset);

                    coarse_offset_array[2][0] = Const_flag_array_type::MASK::Packed_Offset<1,-1>(coarse_offset);
                    coarse_offset_array[2][1] = Const_flag_array_type::MASK::Packed_Offset<1,0>(coarse_offset);
                    coarse_offset_array[2][2] = Const_flag_array_type::MASK::Packed_Offset<1,1>(coarse_offset);

                    //float test = 0.f;
                    for(int i=0;i<2;++i)
                    for(int j=0;j<2;++j){
                        u_fine(fine_offset)+=fine_stencil_scales[i^parity_x][j^parity_y]*
                            u_coarse(coarse_offset_array[1+i-parity_x][1+j-parity_y]);
                    }
                    //PHYSBAM_ASSERT(test == 4.f);
                }
            }
        }
    }
};
//#####################################################################
}
#endif
