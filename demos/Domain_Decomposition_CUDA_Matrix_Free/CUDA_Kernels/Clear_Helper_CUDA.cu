//#####################################################################
//Removed For Anonymity: Copyright Authors of Submission pap173s1
// Distributed under the FreeBSD license (see license.txt)
//#####################################################################
#include <iostream>
#include "Clear_Helper_CUDA.h"
#include <SPGrid_Fluids/Solvers/Domain_Decomposition/SPGRID_MULTIGRID_FLAGS.h>
#include <iostream>
using namespace SPGrid;
#define THREADBLOCK 1024
#define PREFETCH 1024

template <class T,int page_size,class T_offset_ptr>
__global__ void Clear_Kernel(T* u){
    unsigned long current_address=((unsigned long)blockIdx.x*page_size+threadIdx.x*sizeof(T));
    *reinterpret_cast<T*>((unsigned long)u + current_address)=0;
}

template <class T, int log2_struct,int xsize,int ysize,int zsize,class T_offset_ptr>
__global__ void Masked_Clear_Kernel(const unsigned int flag_to_clear,const unsigned int* flags,T* u,const T_offset_ptr* offsets,int max_block){
    enum {
        DATABLOCK=xsize*ysize*zsize,
        span=THREADBLOCK/DATABLOCK,
        xstride=ysize*zsize,
        ystride=zsize,
        zstride=1  
    };
    const unsigned int block = threadIdx.x / DATABLOCK;
    const unsigned int entry = threadIdx.x % DATABLOCK;
    typedef SPGrid_Mask<log2_struct, NextLogTwo<sizeof(T)>::value,3> T_MASK;
    __shared__ T_offset_ptr block_index[PREFETCH];
    if(threadIdx.x < PREFETCH)
        if(blockIdx.x * PREFETCH + threadIdx.x < max_block){
            block_index[threadIdx.x] = offsets[blockIdx.x * PREFETCH + threadIdx.x];
        }
    __syncthreads();
    
    for(int i = 0;i < PREFETCH;i += span){
        if (blockIdx.x * PREFETCH + i + block < max_block){
            if(reinterpret_cast<unsigned int*>((unsigned long)flags + (unsigned long)block_index[i + block])[entry]&flag_to_clear){
            reinterpret_cast<T*>((unsigned long)u + (unsigned long)block_index[i + block])[entry]=0;
            }
        }
    }
}
//#####################################################################
// Function Run_Index_Range
//#####################################################################
// T_MASK corresponds to the mask for the data (not the mask channel)
template <class T,int log2_struct,class T_offset_ptr> 
void Clear_Helper_CUDA<T,log2_struct,3,T_offset_ptr>::Run(cudaStream_t& cuda_stream)
{
    if(size==0)return;
    enum{data_per_page=block_xsize*block_ysize*block_zsize,page_size=4096};
    if(page_size!=data_per_page*(1<<log2_struct)) std::cerr<<"Wrong Page Size!"<<std::endl;
    Clear_Kernel<T,page_size,T_offset_ptr>
        <<<size,data_per_page,0,cuda_stream>>>
        (u);
    // cudaDeviceSynchronize(); 
    // cudaError err = cudaGetLastError();
    // if(cudaSuccess != err){
    //     std::cerr << "Error in Clear. Msg: "<< cudaGetErrorString(err) << std::endl;
    //     abort();
    // }
}
//#####################################################################
template class Clear_Helper_CUDA<float,4,3,unsigned int>;

//#####################################################################
// Function Run_Index_Range
//#####################################################################
// T_MASK corresponds to the mask for the data (not the mask channel)
template <class T,int log2_struct,class T_offset_ptr> 
void Masked_Clear_Helper_CUDA<T,log2_struct,3,T_offset_ptr>::Run(cudaStream_t& cuda_stream)
{
    if(size==0)return;
    int number_of_blocks=size;
    int number_of_cuda_blocks = (number_of_blocks%PREFETCH)?(number_of_blocks/PREFETCH+1):(number_of_blocks/PREFETCH);
    if(number_of_cuda_blocks == 0) return;
    Masked_Clear_Kernel<T,log2_struct,block_xsize,block_ysize,block_zsize,T_offset_ptr>
        <<<number_of_cuda_blocks,THREADBLOCK,0,cuda_stream>>>
        (flag_to_clear,flags,u,b,number_of_blocks);
    // cudaDeviceSynchronize(); 
    // cudaError err = cudaGetLastError();
    // if(cudaSuccess != err){
    //     std::cerr << "Error in Masked Clear. Msg: "<< cudaGetErrorString(err) << ". Trying to clear "<<size<<" blocks."<< std::endl;
    //     abort();
    // }
}
//#####################################################################
template class Masked_Clear_Helper_CUDA<float,4,3,unsigned int>;

