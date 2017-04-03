//#####################################################################
//Removed For Anonymity: Copyright Authors of Submission pap173s1
// Distributed under the FreeBSD license (see license.txt)
//#####################################################################
#ifndef __Prolongation_Helper_CUDA__
#define __Prolongation_Helper_CUDA__
#include <algorithm>
#include <SPGrid/Core/SPGrid_Mask.h>
#include <iostream>
#include <SPGrid_Fluids/Solvers/Domain_Decomposition/SPGRID_MULTIGRID_FLAGS.h>
#include <bitset>
namespace SPGrid{
template<class T,int log2_struct, int d,class T_offset_ptr> class Prolongation_Helper_CUDA;

#define SIMD_LENGTH 16
template<class T,int log2_struct,class T_offset_ptr>
class Prolongation_Helper_CUDA<T,log2_struct,3,T_offset_ptr>
{
    enum{d=3};
    typedef SPGrid_Mask<log2_struct, NextLogTwo<sizeof(T)>::value,d> T_MASK;
    T* const u_fine;         // output stream
    const T* const u_coarse;  // input stream
    const unsigned* const flags;
    T_offset_ptr* offset_coarse_array;   // block offset stream
    T_offset_ptr* offset_fine_array;   // block offset stream
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


    unsigned GetLinear(int x,int y,int z){
        return ((z%block_zsize) + (y%block_ysize) * block_zsize + (x%block_xsize) * block_ysize * block_zsize)*sizeof(T);
    }
    
public:
    explicit Prolongation_Helper_CUDA(T* const u_fine_input,const T* const u_coarse_input,const unsigned* const flags_input
                                      ,T_offset_ptr* offset_fine_array_input
                                      ,T_offset_ptr*  offset_coarse_array_input
                                      ,const int size_input)
        :u_fine(u_fine_input),u_coarse(u_coarse_input),flags(flags_input),offset_coarse_array(offset_coarse_array_input),
         offset_fine_array(offset_fine_array_input),
         size(size_input)
    {
        //STATIC_ASSERT(block_ysize*block_zsize>=SIMD_width,"Z-aligned face is not a multiple of SIMD width");
        //STATIC_ASSERT(block_zsize<=SIMD_width,"Block Z direction exceeds SIMD width");
    }
    
    void Run(cudaStream_t& cuda_stream)
    {Run_Index_Range(0,size-1,cuda_stream);} 
    
    //#####################################################################
    void Run_Index_Range(const unsigned int index_start,const unsigned int index_end,cudaStream_t& cuda_stream);
    //#####################################################################    
 };
};
#endif
