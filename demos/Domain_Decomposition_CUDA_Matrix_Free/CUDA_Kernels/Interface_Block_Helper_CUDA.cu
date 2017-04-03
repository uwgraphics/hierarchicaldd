//#####################################################################
//Removed For Anonymity: Copyright Authors of Submission pap173s1
// Distributed under the FreeBSD license (see license.txt)
//#####################################################################
#include "Interface_Block_Helper_CUDA.h"
using namespace SPGrid;

#define THREADBLOCK 1024
#define PREFETCH 1024

using namespace SPGrid;
//#####################################################################
// Kernel Interface_Collecet_Kernel
//#####################################################################
template<class T,int total_number_of_elements_per_block,class T_offset_ptr> 
__global__ void Interface_Collecet_Kernel(T* const interface_data, const T* const data,const T_offset_ptr* const interface_offsets,unsigned int number_of_blocks){
    const int span = THREADBLOCK / total_number_of_elements_per_block;
    const int block = threadIdx.x / total_number_of_elements_per_block;
    const unsigned long entry = threadIdx.x % total_number_of_elements_per_block;
    __shared__ T_offset_ptr block_index[PREFETCH];
    if(threadIdx.x < PREFETCH)
        if(blockIdx.x * PREFETCH + threadIdx.x < number_of_blocks)
            block_index[threadIdx.x] = interface_offsets[blockIdx.x * PREFETCH + threadIdx.x];
    __syncthreads();
    for(int i = 0;i < PREFETCH;i += span){
        unsigned int current_block = blockIdx.x * PREFETCH + i + block;
        if (current_block < number_of_blocks){
            T* interface_data_ptr = &(reinterpret_cast<T*>((unsigned long)interface_data + (unsigned long)current_block * total_number_of_elements_per_block * sizeof(T))[entry]);
            const T* data_ptr = &(reinterpret_cast<const T*>((unsigned long)data + block_index[i+block])[entry]);
            *interface_data_ptr=*data_ptr;
        }
    }
}
//#####################################################################
// Function Run_Index_Range
//#####################################################################
// T_MASK corresponds to the mask for the data (not the mask channel)
template <class T,int log2_struct,class T_offset_ptr> 
void Interface_Collect_Helper_CUDA<T,log2_struct,3,T_offset_ptr>::Run_Index_Range(const unsigned int index_start,
                                                                                  const unsigned int index_end,
                                                                                  cudaStream_t& cuda_stream)
{
    enum{n_channels=(1<<log2_struct)/sizeof(T)};
    int number_of_interface_blocks=index_end-index_start+1;
    int number_of_cuda_blocks=(number_of_interface_blocks%PREFETCH)?(number_of_interface_blocks/PREFETCH+1):(number_of_interface_blocks/PREFETCH);
    if(number_of_cuda_blocks==0)return;
    if(block_xsize*block_ysize*block_zsize*n_channels>THREADBLOCK) std::cerr<<"Not Enough Threads! In function Interface_Collect_Helper."<<std::endl;
    Interface_Collecet_Kernel<T,block_xsize*block_ysize*block_zsize*n_channels/*times the channels here because we are copying everything!*/,T_offset_ptr>
        <<<number_of_cuda_blocks,THREADBLOCK,0,cuda_stream>>>
        (interface_data,data,b_interface,number_of_interface_blocks);
    
    // cudaDeviceSynchronize();
    // cudaError err = cudaGetLastError();
    // if(cudaSuccess != err){
    //     std::cerr << "Error in Interface_Collect_Helper. Msg: "<< cudaGetErrorString(err) << std::endl;
    //     abort();
    // }

}
//#####################################################################
template class Interface_Collect_Helper_CUDA<float,4,3,unsigned int>;
//#####################################################################

//#####################################################################
// Kernel Interface_Distribute_Kernel
//#####################################################################
template<class T,int number_of_elements_per_block,int page_size,class T_offset_ptr> 
__global__ void Interface_Distribute_Kernel(const T* const interface_data,T* const data,const unsigned* const interface_flags,const unsigned* const flags,const T_offset_ptr* const interface_offsets,unsigned int number_of_blocks){
    const int span = THREADBLOCK / number_of_elements_per_block;
    const int block = threadIdx.x / number_of_elements_per_block;
    const unsigned long entry = threadIdx.x % number_of_elements_per_block;
    __shared__ T_offset_ptr block_index[PREFETCH];
    if(threadIdx.x < PREFETCH)
        if(blockIdx.x * PREFETCH + threadIdx.x < number_of_blocks)
            block_index[threadIdx.x] = interface_offsets[blockIdx.x * PREFETCH + threadIdx.x];
    __syncthreads();
    for(int i = 0;i < PREFETCH;i += span){
        unsigned int current_block = blockIdx.x * PREFETCH + i + block;
        if (current_block < number_of_blocks){
            const T* interface_data_ptr = &(reinterpret_cast<const T*>((unsigned long)interface_data + (unsigned long)current_block * page_size)[entry]);
            T* data_ptr = &(reinterpret_cast<T*>((unsigned long)data + block_index[i+block])[entry]);
            const unsigned* interface_flag_ptr = &(reinterpret_cast<const unsigned*>((unsigned long)interface_flags + (unsigned long)current_block * page_size)[entry]);
            unsigned* flag_ptr = &(reinterpret_cast<unsigned*>((unsigned long)flags + block_index[i+block])[entry]);
            //We do not to check the original flag. It is just for debug. Remove for performance in the future
            // something is wrong...do something creazy!
            if((*interface_flag_ptr) & SPGrid_Solver_Cell_Type_Interface){
                if((*flag_ptr)!=(*interface_flag_ptr)) printf("error!\n");
                //printf("aha! %E\n",*interface_data_ptr);
                *data_ptr=*interface_data_ptr;
                //(*flag_ptr)=(*interface_flag_ptr);
            }
        }
    }

}

//#####################################################################
// Function Run_Index_Range
//#####################################################################
// T_MASK corresponds to the mask for the data (not the mask channel)
template <class T,int log2_struct,class T_offset_ptr> 
void Interface_Distribute_Helper_CUDA<T,log2_struct,3,T_offset_ptr>::Run_Index_Range(const unsigned int index_start,
                                                                                     const unsigned int index_end,
                                                                                     cudaStream_t& cuda_stream)
{
    enum{n_channels=(1<<log2_struct)/sizeof(T),page_size=4096};
    int number_of_interface_blocks=index_end-index_start+1;
    int number_of_cuda_blocks=(number_of_interface_blocks%PREFETCH)?(number_of_interface_blocks/PREFETCH+1):(number_of_interface_blocks/PREFETCH);
    if(number_of_cuda_blocks==0)return;
    if(block_xsize*block_ysize*block_zsize*n_channels*sizeof(T)!=page_size)std::cerr<<"Page size wrong!."<<std::endl;
    Interface_Distribute_Kernel<T,block_xsize*block_ysize*block_zsize,page_size,T_offset_ptr>
        <<<number_of_cuda_blocks,THREADBLOCK,0,cuda_stream>>>
        (interface_data,data,interface_flags,flags,b_interface,number_of_interface_blocks);
    // cudaDeviceSynchronize();
    // cudaError err = cudaGetLastError();
    // if(cudaSuccess != err){
    //     std::cerr << "Error in Interface_Distribute_Helper. Msg: "<< cudaGetErrorString(err) << std::endl;
    //     abort();
    // }
}
//#####################################################################
template class Interface_Distribute_Helper_CUDA<float,4,3,unsigned int>;
//#####################################################################
