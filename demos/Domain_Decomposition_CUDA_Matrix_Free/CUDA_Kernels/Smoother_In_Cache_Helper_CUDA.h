//#####################################################################
//Removed For Anonymity: Copyright Authors of Submission pap173s1
// Distributed under the FreeBSD license (see license.txt)
//#####################################################################
#ifndef __Smoother_In_Cache_Helper_CUDA__
#define __Smoother_In_Cache_Helper_CUDA__
#include <iostream>
#include <algorithm>
#include <SPGrid/Core/SPGrid_Mask.h>

namespace SPGrid{
template<class T,int log2_struct, int d,class T_offset_ptr> class Smoother_In_Cache_Helper_CUDA;
template<class T,int log2_struct,class T_offset_ptr>
class Smoother_In_Cache_Helper_CUDA<T,log2_struct,3,T_offset_ptr>
{
    enum{d=3};
    enum{cache_size=96000};//size of the cuda shared memory
    typedef SPGrid_Mask<log2_struct, NextLogTwo<sizeof(T)>::value,d> T_MASK;
    char* data_ptr;
    T_offset_ptr r_offset;
    T_offset_ptr u_offset;
    T_offset_ptr rhs_offset;
    T_offset_ptr mask_offset;
    const T_offset_ptr* const b;   // block offset stream
    const T_offset_ptr* const b_x_minus;   // block offset stream
    const T_offset_ptr* const b_x_plus;   // block offset stream
    const T_offset_ptr* const b_y_minus;   // block offset stream
    const T_offset_ptr* const b_y_plus;   // block offset stream
    const T_offset_ptr* const b_z_minus;   // block offset stream
    const T_offset_ptr* const b_z_plus;   // block offset stream
    const int block_size;     // number of blocks to process
    const T_offset_ptr data_size;     // size of the data buffer
    const int iterations;
    const T omega;
    enum {
        prefetch_degree = 0,
        block_xsize = 1u << T_MASK::block_xbits,
        block_ysize = 1u << T_MASK::block_ybits,
        block_zsize = 1u << T_MASK::block_zbits
    };

public:
    explicit Smoother_In_Cache_Helper_CUDA(char* data_ptr_input,
                                           T_offset_ptr r_offset_input,T_offset_ptr u_offset_input,
                                           T_offset_ptr rhs_offset_input,T_offset_ptr mask_offset_input,
                                           const T_offset_ptr* const b_input,
                                           const T_offset_ptr* const b_x_minus_input,
                                           const T_offset_ptr* const b_x_plus_input,
                                           const T_offset_ptr* const b_y_minus_input,
                                           const T_offset_ptr* const b_y_plus_input,
                                           const T_offset_ptr* const b_z_minus_input,
                                           const T_offset_ptr* const b_z_plus_input,
                                           const int block_size_input,const T_offset_ptr data_size_input,
                                           const T omega_input,const int iterations_input)
    :data_ptr(data_ptr_input),u_offset(u_offset_input),
                            rhs_offset(rhs_offset_input),mask_offset(mask_offset_input),
                            b(b_input),
                            b_x_minus(b_x_minus_input),b_x_plus(b_x_plus_input),
                            b_y_minus(b_y_minus_input),b_y_plus(b_y_plus_input),
                            b_z_minus(b_z_minus_input),b_z_plus(b_z_plus_input),
                            block_size(block_size_input),data_size(data_size_input),
                            omega(omega_input),iterations(iterations_input)
    {
        if(data_size!=((T_offset_ptr)block_size+1)*4096) {std::cerr<<"Number of blocks and buffer size mismatch."<<std::endl;exit(1);}
    }
    void Run_Bottom(cudaStream_t& cuda_stream)
    {Bottom_Smoothing(0,block_size-1,cuda_stream);}
//#####################################################################
    void Bottom_Smoothing(const unsigned int index_start,const unsigned int index_end,cudaStream_t& cuda_stream);
//#####################################################################
};

};
#endif
