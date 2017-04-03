//#####################################################################
//Removed For Anonymity: Copyright Authors of Submission pap173s1
// Distributed under the FreeBSD license (see license.txt)
//#####################################################################
#ifndef __Smoother_Helper_CUDA__
#define __Smoother_Helper_CUDA__
#include <algorithm>
#include <SPGrid/Core/SPGrid_Mask.h>

namespace SPGrid{
template<class T,int log2_struct, int d,class T_offset_ptr> class Smoother_Helper_CUDA;
template<class T,int log2_struct,class T_offset_ptr>
class Smoother_Helper_CUDA<T,log2_struct,3,T_offset_ptr>
{
    enum{d=3};
    typedef SPGrid_Mask<log2_struct, NextLogTwo<sizeof(T)>::value,d> T_MASK;
    T* const r;         // output stream
    const T* const rhs;
    T* const u;
    const unsigned* const mask;
    const T_offset_ptr* const b;   // block offset stream
    const T_offset_ptr* const b_x_minus;   // block offset stream
    const T_offset_ptr* const b_x_plus;   // block offset stream
    const T_offset_ptr* const b_y_minus;   // block offset stream
    const T_offset_ptr* const b_y_plus;   // block offset stream
    const T_offset_ptr* const b_z_minus;   // block offset stream
    const T_offset_ptr* const b_z_plus;   // block offset stream
    const int size;     // number of blocks to process
    const int iterations;
    const T omega;
    enum {
        prefetch_degree = 0,
        block_xsize = 1u << T_MASK::block_xbits,
        block_ysize = 1u << T_MASK::block_ybits,
        block_zsize = 1u << T_MASK::block_zbits
    };

public:
    explicit Smoother_Helper_CUDA(T* const r_input,T* const u_input,const T* const rhs_input,const unsigned* const mask_input,
                                  const T_offset_ptr* const b_input,
                                  const T_offset_ptr* const b_x_minus_input,
                                  const T_offset_ptr* const b_x_plus_input,
                                  const T_offset_ptr* const b_y_minus_input,
                                  const T_offset_ptr* const b_y_plus_input,
                                  const T_offset_ptr* const b_z_minus_input,
                                  const T_offset_ptr* const b_z_plus_input,
                                  const int size_input,const T omega_input,const int iterations_input)
    :r(r_input),u(u_input),rhs(rhs_input),mask(mask_input),
                            b(b_input),
                            b_x_minus(b_x_minus_input),b_x_plus(b_x_plus_input),
                            b_y_minus(b_y_minus_input),b_y_plus(b_y_plus_input),
                            b_z_minus(b_z_minus_input),b_z_plus(b_z_plus_input),
                            size(size_input),omega(omega_input),iterations(iterations_input)
    {}
    void Run_Interior(cudaStream_t& cuda_stream)
    {Interior_Smoothing(0,size-1,cuda_stream);}
    void Run_Boundary(cudaStream_t& cuda_stream)
    {Boundary_Smoothing(0,size-1,cuda_stream);}
    void Run_Bottom(cudaStream_t& cuda_stream)
    {Bottom_Smoothing(0,size-1,cuda_stream);}
//#####################################################################
    void Interior_Smoothing(const unsigned int index_start,const unsigned int index_end,cudaStream_t& cuda_stream);
    void Boundary_Smoothing(const unsigned int index_start,const unsigned int index_end,cudaStream_t& cuda_stream);
    void Bottom_Smoothing(const unsigned int index_start,const unsigned int index_end,cudaStream_t& cuda_stream);
//#####################################################################
};

};
#endif
