//#####################################################################
//Removed For Anonymity: Copyright Authors of Submission pap173s1
// Distributed under the FreeBSD license (see license.txt)
//#####################################################################
#ifndef __Clear_Helper_Phi__
#define __Clear_Helper_Phi__
#include <algorithm>
#include <SPGrid/Core/SPGrid_Utilities.h>
#include <SPGrid/Core/SPGrid_Mask.h>
#include <SPGrid_Fluids/Solvers/Domain_Decomposition/SPGRID_MULTIGRID_FLAGS.h>
#include <immintrin.h>
#include <iostream>

namespace SPGrid{

template<class T,int log2_struct, int d,class T_offset_ptr> class Clear_Helper_PHI;
template<class T,int log2_struct, int d,class T_offset_ptr> class Masked_Clear_Helper_PHI;

template<class T,int log2_struct,class T_offset_ptr>
class Clear_Helper_PHI<T,log2_struct,3,T_offset_ptr>
{
    enum{d=3};
    typedef SPGrid_Mask<log2_struct, NextLogTwo<sizeof(T)>::value,d> T_MASK;
    T* const u;         // channel to clear
    const int size;     // number of blocks to process
   
    enum {SIMD_width=16};

public:
    explicit Clear_Helper_PHI(T* const u_input,const int size_input)
        :u(u_input),size(size_input)
    {}
    
    void Run()
    {Run_Index_Range(0,size-1);} 
    
    //#####################################################################
    void Run_Index_Range(const int index_start,const int index_end);
    //#####################################################################
};
template<class T,int log2_struct,class T_offset_ptr>
class Masked_Clear_Helper_PHI<T,log2_struct,3,T_offset_ptr>
{
    enum{d=3};
    typedef SPGrid_Mask<log2_struct, NextLogTwo<sizeof(T)>::value,d> T_MASK;
    T* const u;         // channel to clear
    const unsigned* const mask;
    const T_offset_ptr* const b;   // block offset stream
    const int size;     // number of blocks to process

    const __m512i clear_mask_v;
    enum{SIMD_width=16};

public:
    explicit Masked_Clear_Helper_PHI(T* const u_input,const unsigned* const mask_input,const T_offset_ptr* const b_input,const int size_input,const unsigned clear_flag)
        :u(u_input),mask(mask_input),b(b_input),size(size_input),clear_mask_v(_mm512_set1_epi32(clear_flag))
    {}

    void Run()
    {Run_Index_Range(0,size-1);} 
    
    //#####################################################################
    void Run_Index_Range(const int index_start,const int index_end);
    //#####################################################################

};
}
#endif
