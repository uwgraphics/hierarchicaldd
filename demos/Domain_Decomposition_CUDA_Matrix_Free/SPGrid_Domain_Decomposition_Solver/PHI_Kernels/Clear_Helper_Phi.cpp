//#####################################################################
//Removed For Anonymity: Copyright Authors of Submission pap173s1
// Distributed under the FreeBSD license (see license.txt)
//#####################################################################
#include "Clear_Helper_Phi.h"
#include <SPGrid_Fluids/Solvers/Domain_Decomposition/SPGRID_MULTIGRID_FLAGS.h>
using namespace SPGrid;
//#####################################################################
// Function Run_Index_Range
//#####################################################################
// T_MASK corresponds to the mask for the data (not the mask channel)
template <class T, int log2_struct,class T_offset_ptr> void Clear_Helper_PHI<T,log2_struct,3,T_offset_ptr>::Run_Index_Range(const int index_start,const int index_end)
{
    #pragma omp parallel for
    for(int index=index_start;index<=index_end;index++)
        for(int i=0;i<T_MASK::elements_per_block;i+=SIMD_width){
            T* u_ptr = reinterpret_cast<T*>((unsigned long)u + (index<<12) + i * sizeof(T));
            _mm512_storenrngo_ps(u_ptr,_mm512_setzero_ps());}
}
template class Clear_Helper_PHI<float,4,3,unsigned>;

//#####################################################################
// Function Run_Index_Range
//#####################################################################
// T_MASK corresponds to the mask for the data (not the mask channel)
template <class T, int log2_struct,class T_offset_ptr> void Masked_Clear_Helper_PHI<T,log2_struct,3,T_offset_ptr>::Run_Index_Range(const int index_start,const int index_end)
{
    #pragma omp parallel for
    for(int index=index_start;index<=index_end;index++)
        for(int i=0;i<T_MASK::elements_per_block;i+=SIMD_width){
            const unsigned* mask_ptr = reinterpret_cast<const unsigned*>((unsigned long)mask + (unsigned long)b[index] + i * sizeof(unsigned));
            T* u_ptr = reinterpret_cast<T*>((unsigned long)u + (unsigned long)b[index] + i * sizeof(T));
            __m512 u_v = _mm512_load_ps(u_ptr);
            __m512i mask_v = _mm512_load_epi32(mask_ptr);
            __mmask16 mask_m = _mm512_test_epi32_mask(mask_v,clear_mask_v);
            u_v = _mm512_mask_blend_ps(mask_m,u_v,_mm512_setzero_ps());
            _mm512_storenrngo_ps(u_ptr,u_v);
        }
}
template class Masked_Clear_Helper_PHI<float,4,3,unsigned>;
