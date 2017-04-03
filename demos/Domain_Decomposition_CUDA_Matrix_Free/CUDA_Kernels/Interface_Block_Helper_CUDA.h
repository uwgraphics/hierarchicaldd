//#####################################################################
//Removed For Anonymity: Copyright Authors of Submission pap173s1
// Distributed under the FreeBSD license (see license.txt)
//#####################################################################
#ifndef __Interface_Block_Helper_CUDA__
#define __Interface_Block_Helper_CUDA__
#include <algorithm>
#include <SPGrid/Core/SPGrid_Utilities.h>
#include <SPGrid/Core/SPGrid_Mask.h>
#include <SPGrid_Fluids/Solvers/Domain_Decomposition/SPGRID_MULTIGRID_FLAGS.h>
#include <cstdio>
#include <iostream>

namespace SPGrid{

template<class T,int log2_struct, int d,class T_offset_ptr> class Interface_Collect_Helper_CUDA;
template<class T,int log2_struct, int d,class T_offset_ptr> class Interface_Distribute_Helper_CUDA;

template<class T,int log2_struct,class T_offset_ptr>
class Interface_Collect_Helper_CUDA<T,log2_struct,3,T_offset_ptr>
{
    enum{d=3};
    typedef SPGrid_Mask<log2_struct, NextLogTwo<sizeof(T)>::value,d> T_MASK;
    T* const interface_data;         // output stream
    const T* const data;             // input stream
    const T_offset_ptr* const b_interface;
    const int interface_block_size;  // number of interface_blocks to process
   
    enum {
        block_xsize = 1u << T_MASK::block_xbits,
        block_ysize = 1u << T_MASK::block_ybits,
        block_zsize = 1u << T_MASK::block_zbits
    };

public:
    explicit Interface_Collect_Helper_CUDA(T* const interface_data_input,const T* const data_input,const T_offset_ptr* const b_interface_input,const int interface_block_size_input)
        :interface_data(interface_data_input),data(data_input),b_interface(b_interface_input),interface_block_size(interface_block_size_input)
    {}
    
    void Run(cudaStream_t& cuda_stream)
    {Run_Index_Range(0,interface_block_size-1,cuda_stream);} 
    
    //#####################################################################
    void Run_Index_Range(const unsigned int index_start,const unsigned int index_end,cudaStream_t& cuda_stream);
    //#####################################################################
};


template<class T,int log2_struct,class T_offset_ptr>
class Interface_Distribute_Helper_CUDA<T,log2_struct,3,T_offset_ptr>
{
    enum{d=3};
    typedef SPGrid_Mask<log2_struct, NextLogTwo<sizeof(T)>::value,d> T_MASK;
    const T* const interface_data;                 // input stream
    T* const data;                                 // output stream
    const unsigned* const interface_flags;         // input stream
    const unsigned* const flags;                   // input stream
    const T_offset_ptr* const b_interface;
    const int interface_block_size;  // number of interface_blocks to process
   
    enum {
        block_xsize = 1u << T_MASK::block_xbits,
        block_ysize = 1u << T_MASK::block_ybits,
        block_zsize = 1u << T_MASK::block_zbits
    };

public:
    explicit Interface_Distribute_Helper_CUDA(const T* const interface_data_input,T* const data_input,
                                              const unsigned* const interface_flags_input,const unsigned* const flags_input,
                                              const T_offset_ptr* const b_interface_input,const int interface_block_size_input)
        :interface_data(interface_data_input),data(data_input),
         interface_flags(interface_flags_input),flags(flags_input),
         b_interface(b_interface_input),interface_block_size(interface_block_size_input)
    {}
    
    void Run(cudaStream_t& cuda_stream)
    {Run_Index_Range(0,interface_block_size-1,cuda_stream);} 
    
    //#####################################################################
    void Run_Index_Range(const unsigned int index_start,const unsigned int index_end,cudaStream_t& cuda_stream);
    //#####################################################################
};

};
#endif
