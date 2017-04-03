//#####################################################################
//Removed For Anonymity: Copyright Authors of Submission pap173s1
// Distributed under the FreeBSD license (see license.txt)
//#####################################################################
#include "Correction_Helper_Phi.h"
using namespace SPGrid;

namespace{
template<unsigned d> struct BitCount;
template<> struct BitCount<0>  {enum {value=0};};
template<unsigned d> struct BitCount {enum {value=(d&1)+BitCount<(d>>1)>::value};};
float Dinv[64] __attribute__((aligned(64))) = {
                               0.f, 1.f/(float)BitCount< 1>::value, 1.f/(float)BitCount< 2>::value, 1.f/(float)BitCount< 3>::value, 
    1.f/(float)BitCount< 4>::value, 1.f/(float)BitCount< 5>::value, 1.f/(float)BitCount< 6>::value, 1.f/(float)BitCount< 7>::value, 
    1.f/(float)BitCount< 8>::value, 1.f/(float)BitCount< 9>::value, 1.f/(float)BitCount<10>::value, 1.f/(float)BitCount<11>::value, 
    1.f/(float)BitCount<12>::value, 1.f/(float)BitCount<13>::value, 1.f/(float)BitCount<14>::value, 1.f/(float)BitCount<15>::value, 
    1.f/(float)BitCount<16>::value, 1.f/(float)BitCount<17>::value, 1.f/(float)BitCount<18>::value, 1.f/(float)BitCount<19>::value,        
    1.f/(float)BitCount<20>::value, 1.f/(float)BitCount<21>::value, 1.f/(float)BitCount<22>::value, 1.f/(float)BitCount<23>::value, 
    1.f/(float)BitCount<24>::value, 1.f/(float)BitCount<25>::value, 1.f/(float)BitCount<26>::value, 1.f/(float)BitCount<27>::value, 
    1.f/(float)BitCount<28>::value, 1.f/(float)BitCount<29>::value, 1.f/(float)BitCount<30>::value, 1.f/(float)BitCount<31>::value, 
    1.f/(float)BitCount<32>::value, 1.f/(float)BitCount<33>::value, 1.f/(float)BitCount<34>::value, 1.f/(float)BitCount<35>::value, 
    1.f/(float)BitCount<36>::value, 1.f/(float)BitCount<37>::value, 1.f/(float)BitCount<38>::value, 1.f/(float)BitCount<39>::value,
    1.f/(float)BitCount<40>::value, 1.f/(float)BitCount<41>::value, 1.f/(float)BitCount<42>::value, 1.f/(float)BitCount<43>::value, 
    1.f/(float)BitCount<44>::value, 1.f/(float)BitCount<45>::value, 1.f/(float)BitCount<46>::value, 1.f/(float)BitCount<47>::value, 
    1.f/(float)BitCount<48>::value, 1.f/(float)BitCount<49>::value, 1.f/(float)BitCount<50>::value, 1.f/(float)BitCount<51>::value, 
    1.f/(float)BitCount<52>::value, 1.f/(float)BitCount<53>::value, 1.f/(float)BitCount<54>::value, 1.f/(float)BitCount<55>::value, 
    1.f/(float)BitCount<56>::value, 1.f/(float)BitCount<57>::value, 1.f/(float)BitCount<58>::value, 1.f/(float)BitCount<59>::value,
    1.f/(float)BitCount<60>::value, 1.f/(float)BitCount<61>::value, 1.f/(float)BitCount<62>::value, 1.f/(float)BitCount<63>::value
};
}

//#####################################################################
// Function Run_Index_Range
//#####################################################################
// T_MASK corresponds to the mask for the data (not the mask channel)
template <class T, int log2_struct,class T_offset_ptr> void Correction_Helper_PHI<T,log2_struct,3,T_offset_ptr>::Run_Index_Range(const int index_start,const int index_end)
{

    __m512 omega_v=_mm512_set1_ps(omega);
    __m512i SPGrid_Cell_Type_Active_v=_mm512_set1_epi32(SPGrid_Solver_Cell_Type_Active);
    __m512i six_LSB_set_v=_mm512_set1_epi32(0x3f);
    
    #pragma omp parallel for
    for(int index=index_start;index<=index_end;index++)
        for(int i=0;i<T_MASK::elements_per_block;i+=SIMD_width){
              unsigned* mask_ptr = reinterpret_cast<unsigned*>((unsigned long)mask + (index<<12) + i * sizeof(unsigned));
            T* u_ptr = reinterpret_cast<T*>((unsigned long)u + (index<<12) + i * sizeof(T));
            const T* r_ptr = reinterpret_cast<const T*>((unsigned long)r + (index<<12) + i * sizeof(T));
            
            _mm_prefetch((const char*)mask_ptr+(4<<12),_MM_HINT_T1);
            _mm_prefetch((const char*)u_ptr+(4<<12),_MM_HINT_T1);
            _mm_prefetch((const char*)r_ptr+(4<<12),_MM_HINT_T1);

            _mm_prefetch((const char*)mask_ptr+(1<<12),_MM_HINT_T0);
            _mm_prefetch((const char*)u_ptr+(1<<12),_MM_HINT_T0);
            _mm_prefetch((const char*)r_ptr+(1<<12),_MM_HINT_T0);

            __m512i mask_v = _mm512_load_epi32(mask_ptr);
            __mmask16 mask_m = _mm512_test_epi32_mask(mask_v,SPGrid_Cell_Type_Active_v);
 
            if(_mm512_mask2int(mask_m)){
                __m512i lookup_v = _mm512_and_epi32(_mm512_srai_epi32(mask_v,BitLength<SPGrid_Solver_Face_Minus_X_Active>::value-1),six_LSB_set_v);
                __m512 u_new_v = _mm512_load_ps(r_ptr);
                __m512 u_v = _mm512_load_ps(u_ptr);
                __m512 scale_v = _mm512_i32gather_ps(lookup_v,Dinv,4);
                u_new_v = _mm512_mul_ps(u_new_v,scale_v);
                u_new_v = _mm512_mul_ps(u_new_v,omega_v);
                u_new_v = _mm512_add_ps(u_new_v,u_v);
                u_new_v = _mm512_mask_blend_ps(mask_m,u_v,u_new_v);
                _mm512_storenrngo_ps(u_ptr,u_new_v);}}
}
template <class T, int log2_struct,class T_offset_ptr> void Correction_Helper_PHI<T,log2_struct,3,T_offset_ptr>::Run_Boundary_Blocks(const T_offset_ptr* offsets,const int number_of_blocks)
{
    __m512 omega_v=_mm512_set1_ps(omega);
    __m512i SPGrid_Cell_Type_Boundary_v=_mm512_set1_epi32(SPGrid_Solver_Cell_Type_Boundary);
    __m512i six_LSB_set_v=_mm512_set1_epi32(0x3f);

    #pragma omp parallel for
    for(int index=0;index<number_of_blocks;index++){
        _mm_prefetch((const char*)(&offsets[index+32]),_MM_HINT_T1);
        _mm_prefetch((const char*)(&offsets[index+16]),_MM_HINT_T0);

        for(int i=0;i<T_MASK::elements_per_block;i+=SIMD_width){
            
            _mm_prefetch((const char*)mask+offsets[index+4]+i*sizeof(unsigned),_MM_HINT_T1);
            _mm_prefetch((const char*)u+offsets[index+4]+i*sizeof(T),_MM_HINT_T1);
            _mm_prefetch((const char*)r+offsets[index+4]+i*sizeof(T),_MM_HINT_T1);

            /*_mm_prefetch((const char*)mask+offsets[index+1]+i*sizeof(unsigned),_MM_HINT_T0);
            _mm_prefetch((const char*)u+offsets[index+1+i*sizeof(T)],_MM_HINT_T0);
            _mm_prefetch((const char*)r+offsets[index+1+i*sizeof(T)],_MM_HINT_T0);*/

            unsigned* mask_ptr = reinterpret_cast<unsigned*>((unsigned long)mask + offsets[index] + i * sizeof(unsigned));
            T* u_ptr = reinterpret_cast<T*>((unsigned long)u + offsets[index] + i * sizeof(T));
            const T* r_ptr = reinterpret_cast<const T*>((unsigned long)r + offsets[index] + i * sizeof(T));
            
            __m512i mask_v = _mm512_load_epi32(mask_ptr);
            __mmask16 mask_m = _mm512_test_epi32_mask(mask_v,SPGrid_Cell_Type_Boundary_v);
            if(_mm512_mask2int(mask_m)){
                __m512i lookup_v = _mm512_and_epi32(_mm512_srai_epi32(mask_v,BitLength<SPGrid_Solver_Face_Minus_X_Active>::value-1),six_LSB_set_v);
                __m512 scale_v = _mm512_i32gather_ps(lookup_v,Dinv,4);
                __m512 u_new_v = _mm512_load_ps(r_ptr);
                __m512 u_v = _mm512_load_ps(u_ptr);
                u_new_v = _mm512_mul_ps(u_new_v,scale_v);
                u_new_v = _mm512_mul_ps(u_new_v,omega_v);
                u_new_v = _mm512_add_ps(u_new_v,u_v);
                u_new_v = _mm512_mask_blend_ps(mask_m,u_v,u_new_v);
                _mm512_storenrngo_ps(u_ptr,u_new_v);}}
    }
}
template <class T, int log2_struct,class T_offset_ptr> void Correction_Helper_PHI<T,log2_struct,3,T_offset_ptr>::Run_Interior_Blocks(const T_offset_ptr* offsets,const int number_of_blocks)
{
    __m512 omega_v=_mm512_set1_ps(omega);
    __m512i SPGrid_Cell_Type_Boundary_v=_mm512_set1_epi32(SPGrid_Solver_Cell_Type_Boundary);
    __m512i SPGrid_Cell_Type_Active_v=_mm512_set1_epi32(SPGrid_Solver_Cell_Type_Active);
    __m512i six_LSB_set_v=_mm512_set1_epi32(0x3f);

    #pragma omp parallel for
    for(int index=0;index<number_of_blocks;index++){
        _mm_prefetch((const char*)(&offsets[index+32]),_MM_HINT_T1);
        _mm_prefetch((const char*)(&offsets[index+16]),_MM_HINT_T0);

        for(int i=0;i<T_MASK::elements_per_block;i+=SIMD_width){
            
            _mm_prefetch((const char*)mask+offsets[index+4]+i*sizeof(unsigned),_MM_HINT_T1);
            _mm_prefetch((const char*)u+offsets[index+4]+i*sizeof(T),_MM_HINT_T1);
            _mm_prefetch((const char*)r+offsets[index+4]+i*sizeof(T),_MM_HINT_T1);

            /*_mm_prefetch((const char*)mask+offsets[index+1]+i*sizeof(unsigned),_MM_HINT_T0);
            _mm_prefetch((const char*)u+offsets[index+1+i*sizeof(T)],_MM_HINT_T0);
            _mm_prefetch((const char*)r+offsets[index+1+i*sizeof(T)],_MM_HINT_T0);*/

            unsigned* mask_ptr = reinterpret_cast<unsigned*>((unsigned long)mask + offsets[index] + i * sizeof(unsigned));
            T* u_ptr = reinterpret_cast<T*>((unsigned long)u + offsets[index] + i * sizeof(T));
            const T* r_ptr = reinterpret_cast<const T*>((unsigned long)r + offsets[index] + i * sizeof(T));
            
            __m512i mask_v = _mm512_load_epi32(mask_ptr);
            __mmask16 mask_m_active = _mm512_test_epi32_mask(mask_v,SPGrid_Cell_Type_Active_v);
            __mmask16 mask_m_boundary = _mm512_test_epi32_mask(mask_v,SPGrid_Cell_Type_Boundary_v);
            //unmask all boundary cells, this works because all boundary cell must be active cells.
            __mmask16 mask_m = _mm512_kxor(mask_m_active,mask_m_boundary);
            //__mmask16 mask_m = mask_m_active;
            if(_mm512_mask2int(mask_m)){
                __m512i lookup_v = _mm512_and_epi32(_mm512_srai_epi32(mask_v,BitLength<SPGrid_Solver_Face_Minus_X_Active>::value-1),six_LSB_set_v);
                __m512 scale_v = _mm512_i32gather_ps(lookup_v,Dinv,4);//4 here is sizeof(float)
                __m512 u_new_v = _mm512_load_ps(r_ptr);
                __m512 u_v = _mm512_load_ps(u_ptr);
                u_new_v = _mm512_mul_ps(u_new_v,scale_v);
                u_new_v = _mm512_mul_ps(u_new_v,omega_v);
                u_new_v = _mm512_add_ps(u_new_v,u_v);
                u_new_v = _mm512_mask_blend_ps(mask_m,u_v,u_new_v);
                _mm512_storenrngo_ps(u_ptr,u_new_v);}}
    }
}

template class Correction_Helper_PHI<float,4,3,unsigned>;
