//#####################################################################
//Removed For Anonymity: Copyright Authors of Submission pap173s1
// Distributed under the FreeBSD license (see license.txt)
//#####################################################################
#ifndef __Restriction_Helper_Phi__
#define __Restriction_Helper_Phi__
#include <algorithm>
#include <SPGrid/Core/SPGrid_Utilities.h>
#include <SPGrid/Core/SPGrid_Mask.h>
#include <SPGrid_Fluids/Solvers/Domain_Decomposition/SPGRID_MULTIGRID_FLAGS.h>
#include <immintrin.h>
#include <iostream>

namespace SPGrid{
template<class T,int log2_struct, int d,class T_offset_ptr> class Restriction_Helper_PHI;

#define SIMD_LENGTH 16
template<class T,int log2_struct,class T_offset_ptr>
 class Restriction_Helper_PHI<T,log2_struct,3,T_offset_ptr>
{
    enum{d=3};
    typedef SPGrid_Mask<log2_struct, NextLogTwo<sizeof(T)>::value,d> T_MASK;
    T* const b;         // output stream
    const T* const r;  // first input stream
    const unsigned* const mask;
    const T_offset_ptr* const offset_coarse;   // block offset stream
    const std_array<T_offset_ptr,27>* offset_fine_array;   // block offset stream
    const int size;     // number of blocks to process
   
    const __m512i SPGrid_Cell_Type_Active_v;
    enum {
        block_xsize = 1u << T_MASK::block_xbits,
        block_ysize = 1u << T_MASK::block_ybits,
        block_zsize = 1u << T_MASK::block_zbits,
        elements_per_block=block_xsize*block_ysize*block_zsize,
        SIMD_width=16,
        vector_per_block=elements_per_block/SIMD_width,
    };
    static_assert(block_ysize*block_zsize>=SIMD_width,"Z-aligned face is not a multiple of SIMD width");
    static_assert(block_zsize<=SIMD_width,"Block Z direction exceeds SIMD width");

public:
    explicit Restriction_Helper_PHI(T* const b_input,const T* const r_input,const unsigned* const mask_input
                                    ,const T_offset_ptr*  const offset_coarse_input
                                    ,const std_array<T_offset_ptr,27>*  offset_fine_array_input
                                    ,const int size_input)
        :b(b_input),r(r_input),mask(mask_input),offset_coarse(offset_coarse_input),
        offset_fine_array(offset_fine_array_input),
        size(size_input),
        SPGrid_Cell_Type_Active_v(_mm512_set1_epi32(SPGrid_Solver_Cell_Type_Active))
            {}
    
    void Run()
    {
        Run_Index_Range(0,size-1);
    } 
    
    //#####################################################################
    void Run_Index_Range(const int index_start,const int index_end);
    //#####################################################################
    
};

};
#endif
