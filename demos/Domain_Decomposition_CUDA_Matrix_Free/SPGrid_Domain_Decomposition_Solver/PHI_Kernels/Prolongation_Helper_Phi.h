//#####################################################################
//Removed For Anonymity: Copyright Authors of Submission pap173s1
// Distributed under the FreeBSD license (see license.txt)
//#####################################################################
#ifndef __Prolongation_Helper_Phi__
#define __Prolongation_Helper_Phi__
#include <algorithm>
#include <SPGrid/Core/SPGrid_Utilities.h>
#include <SPGrid/Core/SPGrid_Mask.h>
#include <SPGrid_Fluids/Solvers/Domain_Decomposition/SPGRID_MULTIGRID_FLAGS.h>
#include <immintrin.h>
#include <iostream>
#include <bitset>

namespace SPGrid{
template<class T,int log2_struct, int d,class T_offset_ptr> class Prolongation_Helper_PHI;

#define SIMD_LENGTH 16
template<class T,int log2_struct,class T_offset_ptr>
class Prolongation_Helper_PHI<T,log2_struct,3,T_offset_ptr>
{
    enum{d=3};
    typedef SPGrid_Mask<log2_struct, NextLogTwo<sizeof(T)>::value,d> T_MASK;
    T* const u_fine;         // output stream
    const T* const u_coarse;  // input stream
    const unsigned* const mask;
    std_array<T_offset_ptr,8>* offset_coarse_array;   // block offset stream
    std_array<T_offset_ptr,8>* offset_fine_array;   // block offset stream
    const int size;     // number of blocks to process
    enum {
        block_xsize = 1u << T_MASK::block_xbits,
        block_ysize = 1u << T_MASK::block_ybits,
        block_zsize = 1u << T_MASK::block_zbits,
        elements_per_block=block_xsize*block_ysize*block_zsize,
        SIMD_width=16,
        vector_per_block=elements_per_block/SIMD_width,
        block_0_flag=1,
        block_1_flag=2,
        block_2_flag=4,
        address_reading_flag_bit_shift=2,
        n_fine_parity=2,
        n_coarse_block_offset=2
    };

    static_assert(block_ysize*block_zsize>=SIMD_width,"Z-aligned face is not a multiple of SIMD width");
    static_assert(block_zsize<=SIMD_width,"Block Z direction exceeds SIMD width");

    unsigned GetLinear(int x,int y,int z){
        return ((z%block_zsize) + (y%block_ysize) * block_zsize + (x%block_xsize) * block_ysize * block_zsize)*sizeof(T);
    }
    
 public:
    explicit Prolongation_Helper_PHI(T* const u_fine_input,const T* const u_coarse_input,const unsigned* const mask_input
                                    ,std_array<T_offset_ptr,8>* offset_fine_array_input
                                    ,std_array<T_offset_ptr,8>*  offset_coarse_array_input
                                     ,const int size_input)
        :u_fine(u_fine_input),u_coarse(u_coarse_input),mask(mask_input),offset_coarse_array(offset_coarse_array_input),
        offset_fine_array(offset_fine_array_input),
        size(size_input)
            {}
    
    void Run()
    {Run_Index_Range(0,size-1);} 
    
    //#####################################################################
    void Run_Index_Range(const int index_start,const int index_end);
    //#####################################################################    
 };
};
#endif
