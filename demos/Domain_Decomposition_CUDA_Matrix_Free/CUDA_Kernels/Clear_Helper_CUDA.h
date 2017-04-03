//#####################################################################
//Removed For Anonymity: Copyright Authors of Submission pap173s1
// Distributed under the FreeBSD license (see license.txt)
//#####################################################################
#ifndef __Clear_Helper_CUDA__
#define __Clear_Helper_CUDA__
#include <algorithm>
#include <SPGrid/Core/SPGrid_Mask.h>

namespace SPGrid{
template<class T,int log2_struct, int d,class T_offset_ptr> class Clear_Helper_CUDA;
template<class T,int log2_struct, int d,class T_offset_ptr> class Masked_Clear_Helper_CUDA;

template<class T,int log2_struct,class T_offset_ptr>
     class Clear_Helper_CUDA<T,log2_struct,3,T_offset_ptr>
{
    enum{d=3};
    typedef SPGrid_Mask<log2_struct, NextLogTwo<sizeof(T)>::value,d> T_MASK;
    T* const u;         // output stream
    const int size;     // number of blocks to process
    enum {
        prefetch_degree = 0,
        block_xsize = 1u << T_MASK::block_xbits,
        block_ysize = 1u << T_MASK::block_ybits,
        block_zsize = 1u << T_MASK::block_zbits
    };

public:
    explicit Clear_Helper_CUDA(T* const u_input,const int size_input)
        :u(u_input),size(size_input)
    {}

    void Run(cudaStream_t& cuda_stream);
};

template<class T,int log2_struct,class T_offset_ptr>
     class Masked_Clear_Helper_CUDA<T,log2_struct,3,T_offset_ptr>
{
    enum{d=3};
    typedef SPGrid_Mask<log2_struct, NextLogTwo<sizeof(T)>::value,d> T_MASK;
    T* const u;         // output stream
    const unsigned int flag_to_clear;
    const unsigned int* const flags;         // output stream
    const T_offset_ptr* const b;   // block offset stream
    const int size;     // number of blocks to process
    enum {
        prefetch_degree = 0,
        block_xsize = 1u << T_MASK::block_xbits,
        block_ysize = 1u << T_MASK::block_ybits,
        block_zsize = 1u << T_MASK::block_zbits
    };

public:
    explicit Masked_Clear_Helper_CUDA(const unsigned int flag_to_clear_input,const unsigned* const flags_input,T* const u_input,const T_offset_ptr* const b_input,const int size_input)
        :u(u_input),flags(flags_input),flag_to_clear(flag_to_clear_input),b(b_input),size(size_input)
    {}

    void Run(cudaStream_t& cuda_stream);
};
}
#endif
