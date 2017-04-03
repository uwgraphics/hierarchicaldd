//#####################################################################
//Removed For Anonymity: Copyright Authors of Submission pap173s1
// Distributed under the FreeBSD license (see license.txt)
//#####################################################################
#include "Restriction_Helper_Phi.h"
#include <iostream>
#include <Common_Tools/Math_Tools/RANGE_ITERATOR.h>

using namespace SPGrid;
using namespace PhysBAM;
//#####################################################################
// Function Run_Index_Range
//#####################################################################
// T_MASK corresponds to the mask for the data (not the mask channel)
namespace{
template<int block_ysize,int block_zsize>
struct GatherHelper
{
    enum {no_straddle=0, z_straddle=1, y_straddle=2, yz_straddle=3, x_straddle=4, xz_straddle=5, xy_straddle=6, xyz_straddle=7};
    static_assert(block_zsize==4 || block_zsize==8,"Unsupported size");
    static const int offsets[4][16];
    static const uint16_t masks[4][4];
};
const int GatherHelper<8,8>::offsets[4][16] __attribute__((aligned(64))) = {
    {  0,  1,  0xdeadbeef,  2,  0xdeadbeef,  3,  0xdeadbeef,  4,  8,  9,  0xdeadbeef, 10, 0xdeadbeef, 11, 0xdeadbeef, 12},
    { 0xdeadbeef, 0xdeadbeef, 0xdeadbeef, 0xdeadbeef, 0xdeadbeef, 0xdeadbeef, 0xdeadbeef, -4, 0xdeadbeef, 0xdeadbeef, 0xdeadbeef, 0xdeadbeef, 0xdeadbeef, 0xdeadbeef, 0xdeadbeef, 4 },
    { 0xdeadbeef, 0xdeadbeef, 0xdeadbeef, 0xdeadbeef, 0xdeadbeef, 0xdeadbeef, 0xdeadbeef, 0xdeadbeef, -56, -55, 0xdeadbeef, -54, 0xdeadbeef, -53, 0xdeadbeef, -52 },
    { 0xdeadbeef, 0xdeadbeef, 0xdeadbeef, 0xdeadbeef, 0xdeadbeef, 0xdeadbeef, 0xdeadbeef, 0xdeadbeef, 0xdeadbeef, 0xdeadbeef, 0xdeadbeef, 0xdeadbeef, 0xdeadbeef, 0xdeadbeef, 0xdeadbeef, -60 }
};
const uint16_t GatherHelper<8,8>::masks[4][4] = {
    { 0xabab, 0x0000, 0x0000, 0x0000 },
    { 0x2b2b, 0x8080, 0x0000, 0x0000 },
    { 0x00ab, 0x0000, 0xab00, 0x0000 },
    { 0x002b, 0x0080, 0x2b00, 0x8000 }
};

template<int block_zsize> void Bilinear_Transpose(__m512&,__m512);
template<>
void Bilinear_Transpose<8>(__m512& Vcoarse,__m512 Vfine)
{
    __m512 Vone_quarter=_mm512_set1_ps(.25f/2.f);
    __m512 Vthree_quarters=_mm512_set1_ps(.75f/2.f);
    __m512 Vtmp1=_mm512_mul_ps(Vfine,Vthree_quarters);
    Vtmp1=_mm512_fmadd_ps(_mm512_swizzle_ps(Vfine,_MM_SWIZ_REG_CDAB),Vone_quarter,Vtmp1);
    __m512 Vtmp2=_mm512_permute4f128_ps(Vtmp1,_MM_PERM_BADC);
    Vtmp1=_mm512_mul_ps(Vtmp1,Vthree_quarters);
    Vtmp1=_mm512_fmadd_ps(Vtmp2,Vone_quarter,Vtmp1);

    Vcoarse=Vtmp1;

    __m512 Vtmp=_mm512_permute4f128_ps(Vcoarse,_MM_PERM_CDAB);
    __mmask16 fourth_and_twelfth_elements=_mm512_int2mask(0x0808);
    Vcoarse=_mm512_mask_add_ps(Vcoarse,fourth_and_twelfth_elements,Vcoarse,_mm512_swizzle_ps(Vtmp,_MM_SWIZ_REG_AAAA));
    __mmask16 second_element_in_lane=_mm512_int2mask(0x2222);
    Vcoarse=_mm512_mask_add_ps(Vcoarse,second_element_in_lane,Vcoarse,_mm512_swizzle_ps(Vcoarse,_MM_SWIZ_REG_CCCC));
}
}
void _mm512_mask_i32accumulate_ps (void* base_addr, __mmask16 k, __m512i vindex, __m512 a, int scale_input)
{
    enum {scale=4};
    if(scale!=scale_input) abort();
    __m512 Vtmp=_mm512_mask_i32gather_ps(_mm512_setzero_ps(),k,vindex,base_addr,scale);
    //scale by a factor of 8 to balance the Prolongtion/Restriction factor.
    Vtmp=_mm512_mask_add_ps(Vtmp,k,Vtmp,a);
    _mm512_mask_i32scatter_ps(base_addr,k,vindex,Vtmp,scale);
}

template<int block_xsize,int block_ysize,int block_zsize,int fine_parity,int coarse_parity>
void RestrictHelper(
                    float* coarse_block_input,
                    int x,
                    int y,
                    __m512 fine_low_v,
                    __m512 fine_high_v
                    ){
    enum{
        SIMD_width=16,
        block_ystep=SIMD_width/block_zsize
    };
    typedef float (&block_type)[block_xsize][block_ysize][block_zsize];
    block_type coarse_block =reinterpret_cast<block_type>(coarse_block_input[0]);

    __m512i center_offsets_v;
    __mmask16 center_mask_v;
    __mmask16 center_special_mask_v;
    if(coarse_parity==GatherHelper<block_ysize,block_zsize>::no_straddle || coarse_parity==GatherHelper<block_ysize,block_zsize>::x_straddle){
        center_offsets_v=_mm512_load_epi32(GatherHelper<block_ysize,block_zsize>::offsets[GatherHelper<block_ysize,block_zsize>::no_straddle]);
        center_mask_v=_mm512_int2mask(GatherHelper<block_ysize,block_zsize>::masks[fine_parity&0x3][GatherHelper<block_ysize,block_zsize>::no_straddle]);
        center_special_mask_v=_mm512_int2mask(GatherHelper<block_ysize,block_zsize>::masks[fine_parity&GatherHelper<block_ysize,block_zsize>::z_straddle][GatherHelper<block_ysize,block_zsize>::no_straddle]);}

    __m512i z_plus_offsets_v;
    __mmask16 z_plus_mask_v;
    __mmask16 z_special_mask_v;
    if(coarse_parity==GatherHelper<block_ysize,block_zsize>::z_straddle || coarse_parity==GatherHelper<block_ysize,block_zsize>::xz_straddle){
        z_plus_offsets_v=_mm512_load_epi32(GatherHelper<block_ysize,block_zsize>::offsets[GatherHelper<block_ysize,block_zsize>::z_straddle]);
        z_plus_mask_v=_mm512_int2mask(GatherHelper<block_ysize,block_zsize>::masks[fine_parity&0x3][GatherHelper<block_ysize,block_zsize>::z_straddle]);
        z_special_mask_v=_mm512_int2mask(GatherHelper<block_ysize,block_zsize>::masks[fine_parity&GatherHelper<block_ysize,block_zsize>::z_straddle][GatherHelper<block_ysize,block_zsize>::z_straddle]);}

    __m512i y_plus_offsets_v;
    __mmask16 y_plus_mask_v;
    if(coarse_parity==GatherHelper<block_ysize,block_zsize>::y_straddle || coarse_parity==GatherHelper<block_ysize,block_zsize>::xy_straddle){
        y_plus_offsets_v=_mm512_load_epi32(GatherHelper<block_ysize,block_zsize>::offsets[GatherHelper<block_ysize,block_zsize>::y_straddle]);
        y_plus_mask_v=_mm512_int2mask(GatherHelper<block_ysize,block_zsize>::masks[fine_parity&0x3][GatherHelper<block_ysize,block_zsize>::y_straddle]);}

    __m512i yz_plus_offsets_v;
    __mmask16 yz_plus_mask_v;
    if(coarse_parity==GatherHelper<block_ysize,block_zsize>::yz_straddle || coarse_parity==GatherHelper<block_ysize,block_zsize>::xyz_straddle){
        yz_plus_offsets_v=_mm512_load_epi32(GatherHelper<block_ysize,block_zsize>::offsets[GatherHelper<block_ysize,block_zsize>::yz_straddle]);
        yz_plus_mask_v=_mm512_int2mask(GatherHelper<block_ysize,block_zsize>::masks[fine_parity&0x3][GatherHelper<block_ysize,block_zsize>::yz_straddle]);}

    __m512 one_quarter_v=_mm512_set1_ps(.25f/2.f);
    __m512 three_quarters_v=_mm512_set1_ps(.75f/2.f);

    __m512 coarse_x_low_v = _mm512_mul_ps(three_quarters_v,fine_low_v);
    coarse_x_low_v = _mm512_fmadd_ps(one_quarter_v,fine_high_v,coarse_x_low_v);
    __m512 coarse_x_high_v = _mm512_mul_ps(three_quarters_v,fine_high_v);
    coarse_x_high_v = _mm512_fmadd_ps(one_quarter_v,fine_low_v,coarse_x_high_v);

    __m512 coarse_parents_low_v=_mm512_setzero_ps();
    __m512 coarse_parents_high_v=_mm512_setzero_ps();

    Bilinear_Transpose<block_zsize>(coarse_parents_low_v,coarse_x_low_v);
    Bilinear_Transpose<block_zsize>(coarse_parents_high_v,coarse_x_high_v);

    int cx = ((fine_parity&GatherHelper<block_ysize,block_zsize>::x_straddle)?(x+block_xsize)/2:x/2);
    int cy = ((fine_parity&GatherHelper<block_ysize,block_zsize>::y_straddle)?(y+block_ysize)/2:y/2);
    int cz = ((fine_parity&GatherHelper<block_ysize,block_zsize>::z_straddle)?(block_zsize)/2:0);

    if(y + block_ystep >= block_ysize){
        if(coarse_parity==GatherHelper<block_ysize,block_zsize>::no_straddle)
            _mm512_mask_i32accumulate_ps(&coarse_block[cx][cy][cz],center_mask_v,center_offsets_v,coarse_parents_low_v,4);
        if(coarse_parity==GatherHelper<block_ysize,block_zsize>::z_straddle)
            _mm512_mask_i32accumulate_ps(&coarse_block[cx][cy][cz],z_plus_mask_v,z_plus_offsets_v,coarse_parents_low_v,4);
        if(coarse_parity==GatherHelper<block_ysize,block_zsize>::y_straddle)
            _mm512_mask_i32accumulate_ps(&coarse_block[cx][cy][cz],y_plus_mask_v,y_plus_offsets_v,coarse_parents_low_v,4);
        if(coarse_parity==GatherHelper<block_ysize,block_zsize>::yz_straddle)
            _mm512_mask_i32accumulate_ps(&coarse_block[cx][cy][cz],yz_plus_mask_v,yz_plus_offsets_v,coarse_parents_low_v,4);
    }else{
        //the whole thing is still inside the block regarding y direction
        if(coarse_parity==GatherHelper<block_ysize,block_zsize>::no_straddle)
            _mm512_mask_i32accumulate_ps(&coarse_block[cx][cy][cz],center_special_mask_v,center_offsets_v,coarse_parents_low_v,4);
        if(coarse_parity==GatherHelper<block_ysize,block_zsize>::z_straddle)
            _mm512_mask_i32accumulate_ps(&coarse_block[cx][cy][cz],z_special_mask_v,z_plus_offsets_v,coarse_parents_low_v,4);
    }
    if(cx+1 < block_xsize){
        if(y + block_ystep >= block_ysize){
            if(coarse_parity==GatherHelper<block_ysize,block_zsize>::no_straddle)
                _mm512_mask_i32accumulate_ps(&coarse_block[cx+1][cy][cz],center_mask_v,center_offsets_v,coarse_parents_high_v,4);
            if(coarse_parity==GatherHelper<block_ysize,block_zsize>::z_straddle)
                _mm512_mask_i32accumulate_ps(&coarse_block[cx+1][cy][cz],z_plus_mask_v,z_plus_offsets_v,coarse_parents_high_v,4);
            if(coarse_parity==GatherHelper<block_ysize,block_zsize>::y_straddle)
                _mm512_mask_i32accumulate_ps(&coarse_block[cx+1][cy][cz],y_plus_mask_v,y_plus_offsets_v,coarse_parents_high_v,4);
            if(coarse_parity==GatherHelper<block_ysize,block_zsize>::yz_straddle)
                _mm512_mask_i32accumulate_ps(&coarse_block[cx+1][cy][cz],yz_plus_mask_v,yz_plus_offsets_v,coarse_parents_high_v,4);            
        }else{
            //the whole thing is still inside the block regarding y direction
            if(coarse_parity==GatherHelper<block_ysize,block_zsize>::no_straddle)
                _mm512_mask_i32accumulate_ps(&coarse_block[cx+1][cy][cz],center_special_mask_v,center_offsets_v,coarse_parents_high_v,4);
            if(coarse_parity==GatherHelper<block_ysize,block_zsize>::z_straddle)
                _mm512_mask_i32accumulate_ps(&coarse_block[cx+1][cy][cz],z_special_mask_v,z_plus_offsets_v,coarse_parents_high_v,4);
        }
    }else{
        if(y + block_ystep >= block_ysize){
            if(coarse_parity==GatherHelper<block_ysize,block_zsize>::x_straddle)
                _mm512_mask_i32accumulate_ps(&coarse_block[0][cy][cz],center_mask_v,center_offsets_v,coarse_parents_high_v,4);
            if(coarse_parity==GatherHelper<block_ysize,block_zsize>::xz_straddle)
                _mm512_mask_i32accumulate_ps(&coarse_block[0][cy][cz],z_plus_mask_v,z_plus_offsets_v,coarse_parents_high_v,4);
            if(coarse_parity==GatherHelper<block_ysize,block_zsize>::xy_straddle)
                _mm512_mask_i32accumulate_ps(&coarse_block[0][cy][cz],y_plus_mask_v,y_plus_offsets_v,coarse_parents_high_v,4);
            if(coarse_parity==GatherHelper<block_ysize,block_zsize>::xyz_straddle)
                _mm512_mask_i32accumulate_ps(&coarse_block[0][cy][cz],yz_plus_mask_v,yz_plus_offsets_v,coarse_parents_high_v,4);              
        }else{
            //the whole thing is still inside the block regarding y direction
            if(coarse_parity==GatherHelper<block_ysize,block_zsize>::x_straddle)
                _mm512_mask_i32accumulate_ps(&coarse_block[0][cy][cz],center_special_mask_v,center_offsets_v,coarse_parents_high_v,4);
            if(coarse_parity==GatherHelper<block_ysize,block_zsize>::xz_straddle)
                _mm512_mask_i32accumulate_ps(&coarse_block[0][cy][cz],z_special_mask_v,z_plus_offsets_v,coarse_parents_high_v,4);
        }
    }
}
template<int block_xsize,int block_ysize,int block_zsize,int fine_parity,int coarse_parity>
void ResctrictBlockPair(float* coarse_block_input,
                          const float* fine_block_input){
    enum{
        SIMD_width=16,
        block_ystep=SIMD_width/block_zsize
    };
   typedef const float (&block_type)[block_xsize][block_ysize][block_zsize];
    block_type fine_block=reinterpret_cast<block_type>(fine_block_input[0]);
    for(int x=0;x<block_xsize;x+=2)
    for(int y=0;y<block_ysize;y+=block_ystep){
        __m512 fine_result_low_v = _mm512_load_ps(&fine_block[x][y][0]);
        __m512 fine_result_high_v = _mm512_load_ps(&fine_block[x+1][y][0]);
       RestrictHelper<block_xsize,block_ysize,block_zsize,fine_parity,coarse_parity>(coarse_block_input,x,y,fine_result_low_v,fine_result_high_v);
    }
}
template<class T,int block_xsize,int block_ysize,int block_zsize,class T_offset_ptr>
void ResctrictSingle(T* coarse_block_input,
                     std_array<T_offset_ptr,27> fine_block_offsets,const T* fine_block_base){
    typedef GatherHelper<block_ysize,block_zsize> T_HELPER;
    ResctrictBlockPair<block_xsize,block_ysize,block_zsize,T_HELPER::xyz_straddle,T_HELPER::xyz_straddle>(coarse_block_input,
                                                                                                          reinterpret_cast<const T *>(fine_block_offsets(0)+(unsigned long)fine_block_base));
    ResctrictBlockPair<block_xsize,block_ysize,block_zsize,T_HELPER::xy_straddle,T_HELPER::xy_straddle>(coarse_block_input,
                                                                                                          reinterpret_cast<const T *>(fine_block_offsets(1)+(unsigned long)fine_block_base));
    ResctrictBlockPair<block_xsize,block_ysize,block_zsize,T_HELPER::xyz_straddle,T_HELPER::xy_straddle>(coarse_block_input,
                                                                                                          reinterpret_cast<const T *>(fine_block_offsets(2)+(unsigned long)fine_block_base));

    ResctrictBlockPair<block_xsize,block_ysize,block_zsize,T_HELPER::xz_straddle,T_HELPER::xz_straddle>(coarse_block_input,
                                                                                                        reinterpret_cast<const T *>(fine_block_offsets(3)+(unsigned long)fine_block_base));
    ResctrictBlockPair<block_xsize,block_ysize,block_zsize,T_HELPER::x_straddle,T_HELPER::x_straddle>(coarse_block_input,
                                                                                                        reinterpret_cast<const T *>(fine_block_offsets(4)+(unsigned long)fine_block_base));
    ResctrictBlockPair<block_xsize,block_ysize,block_zsize,T_HELPER::xz_straddle,T_HELPER::x_straddle>(coarse_block_input,
                                                                                                        reinterpret_cast<const T *>(fine_block_offsets(5)+(unsigned long)fine_block_base));

    ResctrictBlockPair<block_xsize,block_ysize,block_zsize,T_HELPER::xyz_straddle,T_HELPER::xz_straddle>(coarse_block_input,
                                                                                                        reinterpret_cast<const T *>(fine_block_offsets(6)+(unsigned long)fine_block_base));
    ResctrictBlockPair<block_xsize,block_ysize,block_zsize,T_HELPER::xy_straddle,T_HELPER::x_straddle>(coarse_block_input,
                                                                                                        reinterpret_cast<const T *>(fine_block_offsets(7)+(unsigned long)fine_block_base));
    ResctrictBlockPair<block_xsize,block_ysize,block_zsize,T_HELPER::xyz_straddle,T_HELPER::x_straddle>(coarse_block_input,
                                                                                                        reinterpret_cast<const T *>(fine_block_offsets(8)+(unsigned long)fine_block_base));

    ResctrictBlockPair<block_xsize,block_ysize,block_zsize,T_HELPER::yz_straddle,T_HELPER::yz_straddle>(coarse_block_input,
                                                                                                        reinterpret_cast<const T *>(fine_block_offsets(9)+(unsigned long)fine_block_base));
    ResctrictBlockPair<block_xsize,block_ysize,block_zsize,T_HELPER::y_straddle,T_HELPER::y_straddle>(coarse_block_input,
                                                                                                        reinterpret_cast<const T *>(fine_block_offsets(10)+(unsigned long)fine_block_base));
    ResctrictBlockPair<block_xsize,block_ysize,block_zsize,T_HELPER::yz_straddle,T_HELPER::y_straddle>(coarse_block_input,  
                                                                                                       reinterpret_cast<const T *>(fine_block_offsets(11)+(unsigned long)fine_block_base));

    ResctrictBlockPair<block_xsize,block_ysize,block_zsize,T_HELPER::z_straddle,T_HELPER::z_straddle>(coarse_block_input,
                                                                                                       reinterpret_cast<const T *>(fine_block_offsets(12)+(unsigned long)fine_block_base));
    ResctrictBlockPair<block_xsize,block_ysize,block_zsize,T_HELPER::no_straddle,T_HELPER::no_straddle>(coarse_block_input,
                                                                                                       reinterpret_cast<const T *>(fine_block_offsets(13)+(unsigned long)fine_block_base));
    ResctrictBlockPair<block_xsize,block_ysize,block_zsize,T_HELPER::z_straddle,T_HELPER::no_straddle>(coarse_block_input,
                                                                                                       reinterpret_cast<const T *>(fine_block_offsets(14)+(unsigned long)fine_block_base));

    ResctrictBlockPair<block_xsize,block_ysize,block_zsize,T_HELPER::yz_straddle,T_HELPER::z_straddle>(coarse_block_input,
                                                                                                       reinterpret_cast<const T *>(fine_block_offsets(15)+(unsigned long)fine_block_base));
    ResctrictBlockPair<block_xsize,block_ysize,block_zsize,T_HELPER::y_straddle,T_HELPER::no_straddle>(coarse_block_input,
                                                                                                       reinterpret_cast<const T *>(fine_block_offsets(16)+(unsigned long)fine_block_base));
    ResctrictBlockPair<block_xsize,block_ysize,block_zsize,T_HELPER::yz_straddle,T_HELPER::no_straddle>(coarse_block_input,
                                                                                                       reinterpret_cast<const T *>(fine_block_offsets(17)+(unsigned long)fine_block_base));

    ResctrictBlockPair<block_xsize,block_ysize,block_zsize,T_HELPER::xyz_straddle,T_HELPER::yz_straddle>(coarse_block_input,
                                                                                                       reinterpret_cast<const T *>(fine_block_offsets(18)+(unsigned long)fine_block_base));
    ResctrictBlockPair<block_xsize,block_ysize,block_zsize,T_HELPER::xy_straddle,T_HELPER::y_straddle>(coarse_block_input,
                                                                                                       reinterpret_cast<const T *>(fine_block_offsets(19)+(unsigned long)fine_block_base));
    ResctrictBlockPair<block_xsize,block_ysize,block_zsize,T_HELPER::xyz_straddle,T_HELPER::y_straddle>(coarse_block_input,
                                                                                                       reinterpret_cast<const T *>(fine_block_offsets(20)+(unsigned long)fine_block_base));

    ResctrictBlockPair<block_xsize,block_ysize,block_zsize,T_HELPER::xz_straddle,T_HELPER::z_straddle>(coarse_block_input,
                                                                                                       reinterpret_cast<const T *>(fine_block_offsets(21)+(unsigned long)fine_block_base));
    ResctrictBlockPair<block_xsize,block_ysize,block_zsize,T_HELPER::x_straddle,T_HELPER::no_straddle>(coarse_block_input,
                                                                                                       reinterpret_cast<const T *>(fine_block_offsets(22)+(unsigned long)fine_block_base));
    ResctrictBlockPair<block_xsize,block_ysize,block_zsize,T_HELPER::xz_straddle,T_HELPER::no_straddle>(coarse_block_input,
                                                                                                       reinterpret_cast<const T *>(fine_block_offsets(23)+(unsigned long)fine_block_base));

    ResctrictBlockPair<block_xsize,block_ysize,block_zsize,T_HELPER::xyz_straddle,T_HELPER::z_straddle>(coarse_block_input,
                                                                                                       reinterpret_cast<const T *>(fine_block_offsets(24)+(unsigned long)fine_block_base));
    ResctrictBlockPair<block_xsize,block_ysize,block_zsize,T_HELPER::xy_straddle,T_HELPER::no_straddle>(coarse_block_input,
                                                                                                       reinterpret_cast<const T *>(fine_block_offsets(25)+(unsigned long)fine_block_base));
    ResctrictBlockPair<block_xsize,block_ysize,block_zsize,T_HELPER::xyz_straddle,T_HELPER::no_straddle>(coarse_block_input,
                                                                                                       reinterpret_cast<const T *>(fine_block_offsets(26)+(unsigned long)fine_block_base));
}

template <class T, int log2_struct,class T_offset_ptr> void Restriction_Helper_PHI<T,log2_struct,3,T_offset_ptr>::Run_Index_Range(const int index_start,const int index_end)
{
    typedef VECTOR<int,d> T_INDEX;
    #pragma omp parallel for
    for(int index = index_start; index <= index_end;++index){
        ResctrictSingle<T, block_xsize,block_ysize,block_zsize>(reinterpret_cast<T *>((unsigned long)b + offset_coarse[index]),offset_fine_array[index],r);
        // Auto Vectorization?
        for(int i = 0;i < elements_per_block;++i){
            if(!(reinterpret_cast<const unsigned *>((unsigned long)mask + offset_coarse[index])[i] & SPGrid_Solver_Cell_Type_Active)){
                reinterpret_cast<T *>((unsigned long)b + offset_coarse[index])[i] = 0.f;
            }
        }
    }
}
template class Restriction_Helper_PHI<float,4,3,unsigned>;
template class Restriction_Helper_PHI<float,4,3,unsigned long>;

