//#####################################################################
//Removed For Anonymity: Copyright Authors of Submission pap173s1
// Distributed under the FreeBSD license (see license.txt)
//#####################################################################
#include "Minus_Laplace_Helper_Phi.h"
#include <iostream>
using namespace SPGrid;

//#####################################################################
// Function Run_Index_Range
//#####################################################################
// T_MASK corresponds to the mask for the data (not the mask channel)
template <class T, int log2_struct,class T_offset_ptr,bool accumulative> void Minus_Laplace_Helper_PHI<T,log2_struct,3,T_offset_ptr,accumulative>::Run_Index_Range(const int index_start,const int index_end)
{
    float* tmp = (float*)_mm_malloc(64,64);
    static_assert(elements_per_block%SIMD_LENGTH==0,"elements per block is not multiple of SIMD_LENTH");
    static_assert(sizeof(T)==sizeof(float),"please using float!");
    static_assert(sizeof(unsigned)==sizeof(unsigned int),"please using unsigned int!");
    /////////All the consts.....///////////////////
    enum{
        vectors_per_block=elements_per_block/SIMD_LENGTH
    };
    
    const __m512i zero_v = _mm512_setzero_epi32();

    __mmask16 z_plus_lo_m=_mm512_int2mask(Z_PLUS_MASK_LO);
    __mmask16 z_minus_lo_m=_mm512_int2mask(Z_MINUS_MASK_LO);
    __mmask16 y_plus_lo_m=_mm512_int2mask(Y_PLUS_MASK_LO);
    __mmask16 y_minus_lo_m=_mm512_int2mask(Y_MINUS_MASK_LO);

    typedef const T (&const_block_array)[block_xsize][block_ysize][block_zsize];
    typedef T (&block_array)[block_xsize][block_ysize][block_zsize];
    typedef const unsigned (&const_block_mask)[block_xsize][block_ysize][block_zsize];    

    enum{prefetchl2 = 4,prefetchl1 = 8};
    #pragma omp parallel for
    for(int index=index_start;index<=index_end;index++) { 

        const_block_mask mask_center = reinterpret_cast<const_block_mask>(*reinterpret_cast<const unsigned*>((unsigned long)mask + b[index]));
        const_block_array u_center = reinterpret_cast<const_block_array>(*reinterpret_cast<const T*>((unsigned long)u + b[index]));
        const_block_array u_z_plus = reinterpret_cast<const_block_array>(*reinterpret_cast<const T*>((unsigned long)u + b_z_plus[index]));
        const_block_array u_z_minus = reinterpret_cast<const_block_array>(*reinterpret_cast<const T*>((unsigned long)u + b_z_minus[index]));
        const_block_array u_y_plus = reinterpret_cast<const_block_array>(*reinterpret_cast<const T*>((unsigned long)u + b_y_plus[index]));
        const_block_array u_y_minus = reinterpret_cast<const_block_array>(*reinterpret_cast<const T*>((unsigned long)u + b_y_minus[index]));
        const_block_array u_x_plus = reinterpret_cast<const_block_array>(*reinterpret_cast<const T*>((unsigned long)u + b_x_plus[index]));
        const_block_array u_x_minus = reinterpret_cast<const_block_array>(*reinterpret_cast<const T*>((unsigned long)u + b_x_minus[index]));

        const_block_mask mask_center_l2 = reinterpret_cast<const_block_mask>(*reinterpret_cast<const unsigned*>((unsigned long)mask + b[index+prefetchl2]));
        const_block_array u_center_l2 = reinterpret_cast<const_block_array>(*reinterpret_cast<const T*>((unsigned long)u + b[index+prefetchl2]));
        const_block_array u_z_plus_l2 = reinterpret_cast<const_block_array>(*reinterpret_cast<const T*>((unsigned long)u + b_z_plus[index+prefetchl2]));
        const_block_array u_z_minus_l2 = reinterpret_cast<const_block_array>(*reinterpret_cast<const T*>((unsigned long)u + b_z_minus[index+prefetchl2]));
        const_block_array u_y_plus_l2 = reinterpret_cast<const_block_array>(*reinterpret_cast<const T*>((unsigned long)u + b_y_plus[index+prefetchl2]));
        const_block_array u_y_minus_l2 = reinterpret_cast<const_block_array>(*reinterpret_cast<const T*>((unsigned long)u + b_y_minus[index+prefetchl2]));
        const_block_array u_x_plus_l2 = reinterpret_cast<const_block_array>(*reinterpret_cast<const T*>((unsigned long)u + b_x_plus[index+prefetchl2]));
        const_block_array u_x_minus_l2 = reinterpret_cast<const_block_array>(*reinterpret_cast<const T*>((unsigned long)u + b_x_minus[index+prefetchl2]));

        block_array Lu_center = reinterpret_cast<block_array>(*reinterpret_cast<T*>((unsigned long)Lu + b[index]));

        __m512 tmp_v,z_plus_v,z_minus_v,y_plus_v,y_minus_v,x_plus_v,x_minus_v;
        __m512 center_v;
        __m512i mask_v;

        //#pragma noprefetch
        for(int i=0;i<block_xsize;i++)
        for(int j=0;j<block_ysize;j+=SIMD_width/block_zsize){

            //_mm_prefetch((const char*)&mask_center_l2[i][j][0],_MM_HINT_T1);
            //_mm_prefetch((const char*)&u_center_l2[i][j][0],_MM_HINT_T1);
            //_mm_prefetch((const char*)&u_z_minus_l2[i][j][0],_MM_HINT_T1);
            //_mm_prefetch((const char*)&u_z_plus_l2[i][j][0],_MM_HINT_T1);
            
            center_v = _mm512_load_ps(&u_center[i][j][0]);
            __m512 result_v; 
            if(accumulative)
                result_v=_mm512_sub_ps(_mm512_setzero_ps(),_mm512_load_ps(&Lu_center[i][j][0]));//we have minus sign here because we compute result = Lu - u first. Then negate it.
            else
                result_v=_mm512_setzero_ps();
            
            mask_v = _mm512_load_epi32(&mask_center[i][j][0]);
            __mmask16 mask_m = _mm512_test_epi32_mask(mask_v,SPGrid_Cell_Type_Active_v);
            if(_mm512_mask2int(mask_m)){
                
                tmp_v = _mm512_loadunpackhi_ps(_mm512_setzero_ps(),&u_z_plus[i][j-1][1]+16);
                z_plus_v = _mm512_loadunpacklo_ps(_mm512_setzero_ps(),&u_center[i][j][1]);
                z_plus_v = _mm512_mask_blend_ps(z_plus_lo_m,tmp_v,z_plus_v);
                
                tmp_v = _mm512_loadunpackhi_ps(_mm512_setzero_ps(),&u_center[i][j][-1]+16);
                z_minus_v = _mm512_loadunpacklo_ps(_mm512_setzero_ps(),&u_z_minus[i][j+1][-1]);
                z_minus_v = _mm512_mask_blend_ps(z_minus_lo_m,tmp_v,z_minus_v);
                
                if(yspan_simd > 1){
                    //it is not an aligned read......
                    if(j+SIMD_width/block_zsize>=block_ysize){
                        //_mm_prefetch((const char*)&u_y_plus_l2[i][0][0],_MM_HINT_T1);
                        tmp_v = _mm512_loadunpackhi_ps(_mm512_setzero_ps(),&u_y_plus[i][1-SIMD_width/block_zsize][0]+16);
                    }else{
                        tmp_v = _mm512_loadunpackhi_ps(_mm512_setzero_ps(),&u_center[i][j+1][0]+16);
                    }
                    y_plus_v = _mm512_loadunpacklo_ps(_mm512_setzero_ps(),&u_center[i][j+1][0]);
                    y_plus_v = _mm512_mask_blend_ps(y_plus_lo_m,tmp_v,y_plus_v);            
                }else{
                    //it is actually an aligned read......
                    if(j==block_ysize-1){
                        //_mm_prefetch((const char*)&u_y_plus_l2[i][0][0],_MM_HINT_T1);
                        y_plus_v = _mm512_load_ps(&u_y_plus[i][0][0]);
                    }else{
                        y_plus_v = _mm512_load_ps(&u_center[i][j+1][0]);
                    }
                }
                ////Y MINIUS            
                if(yspan_simd > 1){
                    //it is not an aligned read......
                    tmp_v = _mm512_loadunpackhi_ps(_mm512_setzero_ps(),&u_center[i][j-1][0]+16);         
                    if(j==0){
                        //_mm_prefetch((const char*)&u_y_minus_l2[i][block_ysize-1][0],_MM_HINT_T1);
                        y_minus_v = _mm512_loadunpacklo_ps(_mm512_setzero_ps(),&u_y_minus[i][block_ysize-1][0]);
                    }else{
                         y_minus_v = _mm512_loadunpacklo_ps(_mm512_setzero_ps(),&u_center[i][j-1][0]);
                    }
                    y_minus_v = _mm512_mask_blend_ps(y_minus_lo_m,tmp_v,y_minus_v);      
                }else{
                    //it is actually an aligned read......
                    if(j==0){
                        //_mm_prefetch((const char*)&u_y_minus_l2[i][block_ysize-1][0],_MM_HINT_T1);
                        y_minus_v = _mm512_load_ps(&u_y_minus[i][block_ysize-1][0]);
                    }else{
                        y_minus_v = _mm512_load_ps(&u_center[i][j-1][0]);
                    }
                }

                if(i==block_xsize-1){
                    //_mm_prefetch((const char*)&u_x_plus_l2[0][j][0],_MM_HINT_T1);
                    x_plus_v = _mm512_load_ps(&u_x_plus[0][j][0]);
                }else{
                    x_plus_v = _mm512_load_ps(&u_center[i+1][j][0]);
                }

                if(i==0){
                    //_mm_prefetch((const char*)&u_x_minus_l2[block_xsize-1][j][0],_MM_HINT_T1);
                    x_minus_v = _mm512_load_ps(&u_x_minus[block_xsize-1][j][0]);
                }else{
                    x_minus_v = _mm512_load_ps(&u_center[i-1][j][0]);
                }

                __mmask16 v_mask_as,v_mask_final;
                ////////////////////
                //Z_PLUS
                ////////////////////                    
    
                v_mask_as = _mm512_test_epi32_mask(mask_v,SPGrid_Face_Plus_Z_Active_v);
                v_mask_final = _mm512_kand(v_mask_as,mask_m);
                z_plus_v = _mm512_mask_sub_ps(_mm512_setzero_ps(),v_mask_final,center_v,z_plus_v);
    
                ////////////////////
                //Z_MINUS
                ////////////////////
    
                v_mask_as = _mm512_test_epi32_mask(mask_v,SPGrid_Face_Minus_Z_Active_v);
                v_mask_final = _mm512_kand(v_mask_as,mask_m);
                z_minus_v = _mm512_mask_sub_ps(_mm512_setzero_ps(),v_mask_final,center_v,z_minus_v);
        
                ////////////////////
                //Y_PLUS
                ////////////////////
    
                v_mask_as = _mm512_test_epi32_mask(mask_v,SPGrid_Face_Plus_Y_Active_v);
                v_mask_final = _mm512_kand(v_mask_as,mask_m);
                y_plus_v = _mm512_mask_sub_ps(_mm512_setzero_ps(),v_mask_final,center_v,y_plus_v);
          
                ////////////////////
                //Y_MINUS
                ////////////////////
    
                v_mask_as = _mm512_test_epi32_mask(mask_v,SPGrid_Face_Minus_Y_Active_v);
                v_mask_final = _mm512_kand(v_mask_as,mask_m);
                y_minus_v = _mm512_mask_sub_ps(_mm512_setzero_ps(),v_mask_final,center_v,y_minus_v);
       
                ////////////////////
                //X_PLUS
                ////////////////////
    
                v_mask_as = _mm512_test_epi32_mask(mask_v,SPGrid_Face_Plus_X_Active_v);
                v_mask_final = _mm512_kand(v_mask_as,mask_m);
                x_plus_v = _mm512_mask_sub_ps(_mm512_setzero_ps(),v_mask_final,center_v,x_plus_v);
    
                ////////////////////
                //X_MINUS
                ////////////////////
    
                v_mask_as = _mm512_test_epi32_mask(mask_v,SPGrid_Face_Minus_X_Active_v);
                v_mask_final = _mm512_kand(v_mask_as,mask_m);
                x_minus_v = _mm512_mask_sub_ps(_mm512_setzero_ps(),v_mask_final,center_v,x_minus_v);
        
                ///////////////////
                //Add Togeter
                ///////////////////
                result_v = _mm512_add_ps(result_v,z_plus_v);
                result_v = _mm512_add_ps(result_v,z_minus_v);
                result_v = _mm512_add_ps(result_v,y_plus_v);
                result_v = _mm512_add_ps(result_v,y_minus_v);
                result_v = _mm512_add_ps(result_v,x_plus_v);
                result_v = _mm512_mask_add_ps(_mm512_setzero_ps(),mask_m,result_v,x_minus_v);
                

                result_v =_mm512_sub_ps(_mm512_setzero_ps(),result_v);// negate for the minus laplace

                _mm512_storenrngo_ps(&Lu_center[i][j][0],result_v);
            }
        }
    }
}

template class Minus_Laplace_Helper_PHI<float,4,3,unsigned,true>;
template class Minus_Laplace_Helper_PHI<float,4,3,unsigned,false>;
