//#####################################################################
//Removed For Anonymity: Copyright Authors of Submission pap173s1
// Distributed under the FreeBSD license (see license.txt)
//#####################################################################
#ifndef __Restriction_Helper_CUDA__
#define __Restriction_Helper_CUDA__
#include <algorithm>
#include <SPGrid/Core/SPGrid_Mask.h>
#include <SPGrid/Core/SPGrid_Utilities.h>
#include <iostream>
#include <SPGrid_Fluids/Solvers/Domain_Decomposition/SPGRID_MULTIGRID_FLAGS.h>

namespace SPGrid{
template<class T,int log2_struct, int d,class T_offset_ptr> class Restriction_Helper_CUDA;
    
#define SIMD_LENGTH 16
template<class T,int log2_struct,class T_offset_ptr>
class Restriction_Helper_CUDA<T,log2_struct,3,T_offset_ptr>
{
    enum{d=3};
    typedef SPGrid_Mask<log2_struct, NextLogTwo<sizeof(T)>::value,d> T_MASK;
    T* const b;         // output stream
    const T* const r;  // first input stream
    const unsigned* const flags;
    T_offset_ptr* const offset_coarse;   // block offset stream
    const T_offset_ptr* const offset_fine_array;   // block offset stream
    const int size;     // number of blocks to process
    
    enum {
        block_xsize = 1u << T_MASK::block_xbits,
        block_ysize = 1u << T_MASK::block_ybits,
        block_zsize = 1u << T_MASK::block_zbits,
        elements_per_block=block_xsize*block_ysize*block_zsize,
        SIMD_width=16,
        vector_per_block=elements_per_block/SIMD_width};
    
public:
    explicit Restriction_Helper_CUDA(T* const b_input,const T* const r_input,const unsigned* const flags_input
                                     ,T_offset_ptr*  const offset_coarse_input
                                     ,const T_offset_ptr*  offset_fine_array_input
                                     ,const int size_input)
        :b(b_input),r(r_input),flags(flags_input),offset_coarse(offset_coarse_input),
         offset_fine_array(offset_fine_array_input),
         size(size_input)
    {}
    
    void Run(cudaStream_t& cuda_stream)
    {
        Run_Index_Range(0,size-1,cuda_stream);
    }
    
    //#####################################################################
    void Run_Index_Range(const unsigned int index_start,const unsigned int index_end,cudaStream_t& cuda_stream);
    //#####################################################################
    
};
    
};
#endif
