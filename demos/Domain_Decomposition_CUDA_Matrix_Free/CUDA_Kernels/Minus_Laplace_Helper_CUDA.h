//#####################################################################
//Removed For Anonymity: Copyright Authors of Submission pap173s1
// Distributed under the FreeBSD license (see license.txt)
//#####################################################################
#ifndef __Minus_Laplace_Helper_CUDA__
#define __Minus_Laplace_Helper_CUDA__
#include <algorithm>
//#include <cuda_runtime_api.h>
//#include <cuda.h>
#include <SPGrid/Core/SPGrid_Mask.h>
#include <SPGrid_Fluids/Solvers/Domain_Decomposition/SPGRID_MULTIGRID_FLAGS.h>


namespace SPGrid{

template<class T,int log2_struct, int d,class T_offset_ptr,bool accumulative=true> class Minus_Laplace_Helper_CUDA;

template<class T,int log2_struct,class T_offset_ptr,bool accumulative>
class Minus_Laplace_Helper_CUDA<T,log2_struct,3,T_offset_ptr,accumulative>
{
    enum{d=3};
    typedef SPGrid_Mask<log2_struct, NextLogTwo<sizeof(T)>::value,d> T_MASK;
    T* const x;         // output stream
    const T* const y;  // first input stream
    const unsigned* const mask;
    const T_offset_ptr* const b;   // block offset stream
    const T_offset_ptr* const b_x_minus;   // block offset stream
    const T_offset_ptr* const b_x_plus;   // block offset stream
    const T_offset_ptr* const b_y_minus;   // block offset stream
    const T_offset_ptr* const b_y_plus;   // block offset stream
    const T_offset_ptr* const b_z_minus;   // block offset stream
    const T_offset_ptr* const b_z_plus;   // block offset stream
    const int size;     // number of blocks to process
    const unsigned flag_to_check;

    enum {
        prefetch_degree = 0,
        block_xsize = 1u << T_MASK::block_xbits,
        block_ysize = 1u << T_MASK::block_ybits,
        block_zsize = 1u << T_MASK::block_zbits
    };

public:
    explicit Minus_Laplace_Helper_CUDA(T* const x_input,const T* const y_input,const unsigned* const mask_input,
                                       const T_offset_ptr* const b_input,
                                       const T_offset_ptr* const b_x_minus_input,
                                       const T_offset_ptr* const b_x_plus_input,
                                       const T_offset_ptr* const b_y_minus_input,
                                       const T_offset_ptr* const b_y_plus_input,
                                       const T_offset_ptr* const b_z_minus_input,
                                       const T_offset_ptr* const b_z_plus_input,
                                       const int size_input,const unsigned flag_to_check_input=SPGrid_Solver_Cell_Type_Active);
    
    void Run(cudaStream_t& cuda_stream)
    { Run_Index_Range(0,size-1,cuda_stream); }

//#####################################################################
    void Run_Index_Range(const unsigned int index_start,const unsigned int index_end,cudaStream_t& cuda_stream);
//#####################################################################

};

};
#endif
