//#####################################################################
//Removed For Anonymity: Copyright Authors of Submission pap173s1
// Distributed under the FreeBSD license (see license.txt)
//#####################################################################
#include "Prolongation_Helper_Phi.h"
#include <iostream>
#include <Common_Tools/Math_Tools/RANGE_ITERATOR.h>

using namespace SPGrid;
using namespace PhysBAM;

__m512i SPGrid_Cell_Type_Active_v;

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
    {  0,  1,  1,  2,  2,  3,  3,  4,  8,  9,  9, 10, 10, 11, 11, 12},
    { 0xdeadbeef, 0xdeadbeef, 0xdeadbeef, 0xdeadbeef, 0xdeadbeef, 0xdeadbeef, 0xdeadbeef, -4, 0xdeadbeef, 0xdeadbeef, 0xdeadbeef, 0xdeadbeef, 0xdeadbeef, 0xdeadbeef, 0xdeadbeef, 4 },
    { 0xdeadbeef, 0xdeadbeef, 0xdeadbeef, 0xdeadbeef, 0xdeadbeef, 0xdeadbeef, 0xdeadbeef, 0xdeadbeef, -56, -55, -55, -54, -54, -53, -53, -52 },
    { 0xdeadbeef, 0xdeadbeef, 0xdeadbeef, 0xdeadbeef, 0xdeadbeef, 0xdeadbeef, 0xdeadbeef, 0xdeadbeef, 0xdeadbeef, 0xdeadbeef, 0xdeadbeef, 0xdeadbeef, 0xdeadbeef, 0xdeadbeef, 0xdeadbeef, -60 }
};

const uint16_t GatherHelper<8,8>::masks[4][4] = {
    { 0xffff, 0x0000, 0x0000, 0x0000 },
    { 0x7f7f, 0x8080, 0x0000, 0x0000 },
    { 0x00ff, 0x0000, 0xff00, 0x0000 },
    { 0x007f, 0x0080, 0x7f00, 0x8000 }
};

// const int GatherHelper<4,4>::offsets[4][16] __attribute__((aligned(64))) = {
//     {  0,  1,  1,  2,  4,  5,  5,  6,  8,  9,  9, 10, 0xdeadbeef, 0xdeadbeef, 0xdeadbeef, 0xdeadbeef},
//     { 0xdeadbeef, 0xdeadbeef, 0xdeadbeef, -2, 0xdeadbeef, 0xdeadbeef, 0xdeadbeef, 2, 0xdeadbeef, 0xdeadbeef, 0xdeadbeef, 6, 0xdeadbeef, 0xdeadbeef, 0xdeadbeef, 0xdeadbeef },
//     { 0xdeadbeef, 0xdeadbeef, 0xdeadbeef, 0xdeadbeef, 0xdeadbeef, 0xdeadbeef, 0xdeadbeef, 0xdeadbeef, -8, -7, -7, -6, 0xdeadbeef, 0xdeadbeef, 0xdeadbeef, 0xdeadbeef },
//     { 0xdeadbeef, 0xdeadbeef, 0xdeadbeef, 0xdeadbeef, 0xdeadbeef, 0xdeadbeef, 0xdeadbeef, 0xdeadbeef, 0xdeadbeef, 0xdeadbeef, 0xdeadbeef, -10, 0xdeadbeef, 0xdeadbeef, 0xdeadbeef, 0xdeadbeef }
// };

// const uint16_t GatherHelper<4,4>::masks[4][4] = {
//     { 0x0fff, 0x0000, 0x0000, 0x0000 },
//     { 0x0777, 0x0888, 0x0000, 0x0000 },
//     { 0x00ff, 0x0000, 0x0f00, 0x0000 },
//     { 0x0077, 0x0088, 0x0700, 0x0800 }
// };

template<int block_zsize> __m512 Bilinear(const __m512);

template<>
__m512 Bilinear<8>(const __m512 Vin)
{
    __m512 Vone_quarter=_mm512_set1_ps(.25f);
    __m512 Vthree_quarters=_mm512_set1_ps(.75f);
    __m512 Vtmp1=_mm512_mul_ps(Vin,Vthree_quarters);
    Vtmp1=_mm512_fmadd_ps(_mm512_swizzle_ps(Vin,_MM_SWIZ_REG_CDAB),Vone_quarter,Vtmp1);
    __m512 Vtmp2=_mm512_permute4f128_ps(Vtmp1,_MM_PERM_BADC);
    Vtmp1=_mm512_mul_ps(Vtmp1,Vthree_quarters);
    Vtmp1=_mm512_fmadd_ps(Vtmp2,Vone_quarter,Vtmp1);
    return Vtmp1;
}

// template<>
// __m512 Bilinear<4>(const __m512 Vin_raw)
// {
//     __m512 Vin=_mm512_permute4f128_ps(Vin_raw,_MM_PERM_CBBA);
//     __m512 Vone_quarter=_mm512_set1_ps(.25f);
//     __m512 Vthree_quarters=_mm512_set1_ps(.75f);
//     __m512 Vtmp1=_mm512_mul_ps(Vin,Vthree_quarters);
//     Vtmp1=_mm512_fmadd_ps(_mm512_swizzle_ps(Vin,_MM_SWIZ_REG_CDAB),Vone_quarter,Vtmp1);
//     __m512 Vtmp2=_mm512_permute4f128_ps(Vtmp1,_MM_PERM_CDAB);
//     Vtmp1=_mm512_mul_ps(Vtmp1,Vthree_quarters);
//     Vtmp1=_mm512_fmadd_ps(Vtmp2,Vone_quarter,Vtmp1);
//     return Vtmp1;
// }
}

//it is set to be 4 times to compensate the dx*dx terms in the laplace operator
__m512 one_quarter_v=_mm512_set1_ps(.25f*4.0);
__m512 three_quarters_v=_mm512_set1_ps(.75f*4.0);

template<int block_xsize,int block_ysize,int block_zsize,int parity>
void ProlongateSingleTemplatized(float* fine_block_input,
                                 unsigned* mask_input,
                                 float* center_coarse_block_input,
                                 float* z_plus_coarse_block_input,
                                 float* y_plus_coarse_block_input,
                                 float* yz_plus_coarse_block_input,
                                 float* x_plus_coarse_block_input,
                                 float* xz_plus_coarse_block_input,
                                 float* xy_plus_coarse_block_input,
                                 float* xyz_plus_coarse_block_input){
    enum{
        SIMD_width=16,
        elements_per_block=block_xsize*block_ysize*block_zsize,
        block_ystep=SIMD_width/block_zsize
    };

   typedef float (&block_type)[block_xsize][block_ysize][block_zsize];
   typedef unsigned (&mask_type)[block_xsize][block_ysize][block_zsize];
    block_type fine_block=reinterpret_cast<block_type>(fine_block_input[0]);
    mask_type mask_fine=reinterpret_cast<mask_type>(mask_input[0]);
    block_type center_coarse_block =reinterpret_cast<block_type>(center_coarse_block_input[0]);
    block_type z_plus_coarse_block =reinterpret_cast<block_type>(z_plus_coarse_block_input[0]);
    block_type y_plus_coarse_block =reinterpret_cast<block_type>(y_plus_coarse_block_input[0]);
    block_type yz_plus_coarse_block=reinterpret_cast<block_type>(yz_plus_coarse_block_input[0]);
    block_type x_plus_coarse_block =reinterpret_cast<block_type>(x_plus_coarse_block_input[0]);
    block_type xz_plus_coarse_block =reinterpret_cast<block_type>(xz_plus_coarse_block_input[0]);
    block_type xy_plus_coarse_block =reinterpret_cast<block_type>(xy_plus_coarse_block_input[0]);
    block_type xyz_plus_coarse_block=reinterpret_cast<block_type>(xyz_plus_coarse_block_input[0]);

    __m512i center_offsets_v=_mm512_load_epi32(GatherHelper<block_ysize,block_zsize>::offsets[GatherHelper<block_ysize,block_zsize>::no_straddle]);
    __mmask16 center_mask_v=_mm512_int2mask(GatherHelper<block_ysize,block_zsize>::masks[parity&0x3][GatherHelper<block_ysize,block_zsize>::no_straddle]);

    __m512i z_plus_offsets_v=_mm512_load_epi32(GatherHelper<block_ysize,block_zsize>::offsets[GatherHelper<block_ysize,block_zsize>::z_straddle]);
    __mmask16 z_plus_mask_v=_mm512_int2mask(GatherHelper<block_ysize,block_zsize>::masks[parity&0x3][GatherHelper<block_ysize,block_zsize>::z_straddle]);

    __m512i y_plus_offsets_v=_mm512_load_epi32(GatherHelper<block_ysize,block_zsize>::offsets[GatherHelper<block_ysize,block_zsize>::y_straddle]);
    __mmask16 y_plus_mask_v=_mm512_int2mask(GatherHelper<block_ysize,block_zsize>::masks[parity&0x3][GatherHelper<block_ysize,block_zsize>::y_straddle]);

    __m512i yz_plus_offsets_v=_mm512_load_epi32(GatherHelper<block_ysize,block_zsize>::offsets[GatherHelper<block_ysize,block_zsize>::yz_straddle]);
    __mmask16 yz_plus_mask_v=_mm512_int2mask(GatherHelper<block_ysize,block_zsize>::masks[parity&0x3][GatherHelper<block_ysize,block_zsize>::yz_straddle]);
    
    __mmask16 z_special_mask_v=_mm512_int2mask(GatherHelper<block_ysize,block_zsize>::masks[parity&GatherHelper<block_ysize,block_zsize>::z_straddle][GatherHelper<block_ysize,block_zsize>::z_straddle]);

    __mmask16 center_special_mask_v=_mm512_int2mask(GatherHelper<block_ysize,block_zsize>::masks[parity&GatherHelper<block_ysize,block_zsize>::z_straddle][GatherHelper<block_ysize,block_zsize>::no_straddle]);

    for(int x=0;x<block_xsize;x+=2)
    for(int y=0;y<block_ysize;y+=block_ystep){
        __m512i mask_low_v=_mm512_load_epi32(&mask_fine[x][y][0]);
        __mmask16 mask_low_m=_mm512_test_epi32_mask(mask_low_v,SPGrid_Cell_Type_Active_v);

        __m512i mask_high_v=_mm512_load_epi32(&mask_fine[x+1][y][0]);
        __mmask16 mask_high_m=_mm512_test_epi32_mask(mask_high_v,SPGrid_Cell_Type_Active_v);
        if(_mm512_mask2int(mask_high_m) == 0 && _mm512_mask2int(mask_low_m) == 0) continue;
        __m512 u_low_v = _mm512_setzero_ps();
        __m512 u_high_v = _mm512_setzero_ps();
        if(_mm512_mask2int(mask_low_m)) u_low_v = _mm512_load_ps(&fine_block[x][y][0]);
        if(_mm512_mask2int(mask_high_m)) u_high_v = _mm512_load_ps(&fine_block[x+1][y][0]);

        int cx = ((parity&GatherHelper<block_ysize,block_zsize>::x_straddle)?(x+block_xsize)/2:x/2);
        int cy = ((parity&GatherHelper<block_ysize,block_zsize>::y_straddle)?(y+block_ysize)/2:y/2);
        int cz = ((parity&GatherHelper<block_ysize,block_zsize>::z_straddle)?(block_zsize)/2:0);

        __m512 coarse_parents_low_v=_mm512_setzero_ps();
        if(y + block_ystep >= block_ysize){
            coarse_parents_low_v=_mm512_mask_i32gather_ps(coarse_parents_low_v,center_mask_v,center_offsets_v,&center_coarse_block[cx][cy][cz],4);
            coarse_parents_low_v=_mm512_mask_i32gather_ps(coarse_parents_low_v,z_plus_mask_v,z_plus_offsets_v,&z_plus_coarse_block[cx][cy][cz],4);
            coarse_parents_low_v=_mm512_mask_i32gather_ps(coarse_parents_low_v,y_plus_mask_v,y_plus_offsets_v,&y_plus_coarse_block[cx][cy][cz],4);
            coarse_parents_low_v=_mm512_mask_i32gather_ps(coarse_parents_low_v,yz_plus_mask_v,yz_plus_offsets_v,&yz_plus_coarse_block[cx][cy][cz],4);
        }else{
            //the whole thing is still inside the block regarding y direction
            coarse_parents_low_v=_mm512_mask_i32gather_ps(coarse_parents_low_v,center_special_mask_v,center_offsets_v,&center_coarse_block[cx][cy][cz],4);
            coarse_parents_low_v=_mm512_mask_i32gather_ps(coarse_parents_low_v,z_special_mask_v,z_plus_offsets_v,&z_plus_coarse_block[cx][cy][cz],4);
        }
        __m512 coarse_parents_high_v=_mm512_setzero_ps();
        if(cx+1 < block_xsize){
            if(y + block_ystep >= block_ysize){
                coarse_parents_high_v=_mm512_mask_i32gather_ps(coarse_parents_high_v,center_mask_v,center_offsets_v,&center_coarse_block[cx+1][cy][cz],4);
                coarse_parents_high_v=_mm512_mask_i32gather_ps(coarse_parents_high_v,z_plus_mask_v,z_plus_offsets_v,&z_plus_coarse_block[cx+1][cy][cz],4);
                coarse_parents_high_v=_mm512_mask_i32gather_ps(coarse_parents_high_v,y_plus_mask_v,y_plus_offsets_v,&y_plus_coarse_block[cx+1][cy][cz],4);
                coarse_parents_high_v=_mm512_mask_i32gather_ps(coarse_parents_high_v,yz_plus_mask_v,yz_plus_offsets_v,&yz_plus_coarse_block[cx+1][cy][cz],4);            
            }else{
                //the whole thing is still inside the block regarding y direction
                coarse_parents_high_v=_mm512_mask_i32gather_ps(coarse_parents_high_v,center_special_mask_v,center_offsets_v,&center_coarse_block[cx+1][cy][cz],4);
                coarse_parents_high_v=_mm512_mask_i32gather_ps(coarse_parents_high_v,z_special_mask_v,z_plus_offsets_v,&z_plus_coarse_block[cx+1][cy][cz],4);
            }
        }else{
            if(y + block_ystep >= block_ysize){
                coarse_parents_high_v=_mm512_mask_i32gather_ps(coarse_parents_high_v,center_mask_v,center_offsets_v,&x_plus_coarse_block[0][cy][cz],4);
                coarse_parents_high_v=_mm512_mask_i32gather_ps(coarse_parents_high_v,z_plus_mask_v,z_plus_offsets_v,&xz_plus_coarse_block[0][cy][cz],4);
                coarse_parents_high_v=_mm512_mask_i32gather_ps(coarse_parents_high_v,y_plus_mask_v,y_plus_offsets_v,&xy_plus_coarse_block[0][cy][cz],4);
                coarse_parents_high_v=_mm512_mask_i32gather_ps(coarse_parents_high_v,yz_plus_mask_v,yz_plus_offsets_v,&xyz_plus_coarse_block[0][cy][cz],4);              
            }else{
                //the whole thing is still inside the block regarding y direction
                coarse_parents_high_v=_mm512_mask_i32gather_ps(coarse_parents_high_v,center_special_mask_v,center_offsets_v,&x_plus_coarse_block[0][cy][cz],4);
                coarse_parents_high_v=_mm512_mask_i32gather_ps(coarse_parents_high_v,z_special_mask_v,z_plus_offsets_v,&xz_plus_coarse_block[0][cy][cz],4);
            }
        }

        __m512 fine_result_low_v = Bilinear<block_zsize>(coarse_parents_low_v);
        __m512 fine_result_high_v = Bilinear<block_zsize>(coarse_parents_high_v);

        __m512 temp_v=_mm512_mul_ps(fine_result_low_v,one_quarter_v);
        fine_result_low_v=_mm512_mul_ps(fine_result_low_v,three_quarters_v);
        fine_result_low_v=_mm512_fmadd_ps(fine_result_high_v,one_quarter_v,fine_result_low_v);
        fine_result_high_v=_mm512_fmadd_ps(fine_result_high_v,three_quarters_v,temp_v);

        fine_result_low_v = _mm512_mask_add_ps(_mm512_setzero_ps(),mask_low_m,fine_result_low_v,u_low_v);
        fine_result_high_v = _mm512_mask_add_ps(_mm512_setzero_ps(),mask_high_m,fine_result_high_v,u_high_v);

        _mm512_stream_ps(&fine_block[x][y][0],fine_result_low_v);
        _mm512_stream_ps(&fine_block[x+1][y][0],fine_result_high_v);
    }
}

template<int block_xsize,int block_ysize,int block_zsize>
void ProlongateSingle(float* fine_block,
                      unsigned* fine_mask,
                      float* center_coarse_block,
                      float* z_plus_coarse_block,
                      float* y_plus_coarse_block,
                      float* yz_plus_coarse_block,
                      float* x_plus_coarse_block,
                      float* xz_plus_coarse_block,
                      float* xy_plus_coarse_block,
                      float* xyz_plus_coarse_block,
                      int parity){
    switch(parity){
        case 0:
            ProlongateSingleTemplatized<block_xsize,block_ysize,block_zsize,0>(fine_block,fine_mask,center_coarse_block,z_plus_coarse_block,y_plus_coarse_block,yz_plus_coarse_block,
                x_plus_coarse_block,xz_plus_coarse_block,xy_plus_coarse_block,xyz_plus_coarse_block);
            return;
        case 1:
            ProlongateSingleTemplatized<block_xsize,block_ysize,block_zsize,1>(fine_block,fine_mask,center_coarse_block,z_plus_coarse_block,y_plus_coarse_block,yz_plus_coarse_block,
                x_plus_coarse_block,xz_plus_coarse_block,xy_plus_coarse_block,xyz_plus_coarse_block);
            return;
        case 2:
            ProlongateSingleTemplatized<block_xsize,block_ysize,block_zsize,2>(fine_block,fine_mask,center_coarse_block,z_plus_coarse_block,y_plus_coarse_block,yz_plus_coarse_block,
                x_plus_coarse_block,xz_plus_coarse_block,xy_plus_coarse_block,xyz_plus_coarse_block);
            return;
        case 3:
            ProlongateSingleTemplatized<block_xsize,block_ysize,block_zsize,3>(fine_block,fine_mask,center_coarse_block,z_plus_coarse_block,y_plus_coarse_block,yz_plus_coarse_block,
                x_plus_coarse_block,xz_plus_coarse_block,xy_plus_coarse_block,xyz_plus_coarse_block);
            return;
        case 4:
            ProlongateSingleTemplatized<block_xsize,block_ysize,block_zsize,4>(fine_block,fine_mask,center_coarse_block,z_plus_coarse_block,y_plus_coarse_block,yz_plus_coarse_block,
                x_plus_coarse_block,xz_plus_coarse_block,xy_plus_coarse_block,xyz_plus_coarse_block);
            return;
        case 5:
            ProlongateSingleTemplatized<block_xsize,block_ysize,block_zsize,5>(fine_block,fine_mask,center_coarse_block,z_plus_coarse_block,y_plus_coarse_block,yz_plus_coarse_block,
                x_plus_coarse_block,xz_plus_coarse_block,xy_plus_coarse_block,xyz_plus_coarse_block);
            return;
        case 6:
            ProlongateSingleTemplatized<block_xsize,block_ysize,block_zsize,6>(fine_block,fine_mask,center_coarse_block,z_plus_coarse_block,y_plus_coarse_block,yz_plus_coarse_block,
                x_plus_coarse_block,xz_plus_coarse_block,xy_plus_coarse_block,xyz_plus_coarse_block);
            return;
        case 7:
            ProlongateSingleTemplatized<block_xsize,block_ysize,block_zsize,7>(fine_block,fine_mask,center_coarse_block,z_plus_coarse_block,y_plus_coarse_block,yz_plus_coarse_block,
                x_plus_coarse_block,xz_plus_coarse_block,xy_plus_coarse_block,xyz_plus_coarse_block);
            return;
    }
}

template
void ProlongateSingle<4,8,8>(
    float* fine_block_input,
    unsigned* fine_mask_input,
    float* center_coarse_block_input,
    float* z_plus_coarse_block_input,
    float* y_plus_coarse_block_input,
    float* yz_plus_coarse_block_input,
    float* x_plus_coarse_block_input,
    float* xz_plus_coarse_block_input,
    float* xy_plus_coarse_block_input,
    float* xyz_plus_coarse_block_input,
    int parity);

// T_MASK corresponds to the mask for the data (not the mask channel)
template <class T, int log2_struct,class T_offset_ptr> void Prolongation_Helper_PHI<T,log2_struct,3,T_offset_ptr>::Run_Index_Range(const int index_start,const int index_end)
{
    SPGrid_Cell_Type_Active_v = _mm512_set1_epi32(SPGrid_Solver_Cell_Type_Active);
    typedef VECTOR<int,d> T_INDEX;
    #pragma omp parallel for
    for(int index = index_start; index <= index_end;++index){
        for(int parity = 0;parity < 8;++parity){
            if(offset_fine_array[index](parity) == 0xdeadbeef) continue;
            ProlongateSingle<block_xsize,block_ysize,block_zsize>(reinterpret_cast<T*>(offset_fine_array[index](parity) + (unsigned long)u_fine),
                                                                  reinterpret_cast<unsigned *>(offset_fine_array[index](parity) + (unsigned long)mask),
                                                                  reinterpret_cast<T*>(offset_coarse_array[index](0) + (unsigned long)u_coarse),
                                                                  reinterpret_cast<T*>(offset_coarse_array[index](1) + (unsigned long)u_coarse),
                                                                  reinterpret_cast<T*>(offset_coarse_array[index](2) + (unsigned long)u_coarse),
                                                                  reinterpret_cast<T*>(offset_coarse_array[index](3) + (unsigned long)u_coarse),
                                                                  reinterpret_cast<T*>(offset_coarse_array[index](4) + (unsigned long)u_coarse),
                                                                  reinterpret_cast<T*>(offset_coarse_array[index](5) + (unsigned long)u_coarse),
                                                                  reinterpret_cast<T*>(offset_coarse_array[index](6) + (unsigned long)u_coarse),
                                                                  reinterpret_cast<T*>(offset_coarse_array[index](7) + (unsigned long)u_coarse),
                                                                  parity);
        }
    }        
}
template class Prolongation_Helper_PHI<float,4,3,unsigned>;
template class Prolongation_Helper_PHI<float,4,3,unsigned long>;
