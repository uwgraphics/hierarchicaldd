//#####################################################################
//Removed For Anonymity: Copyright Authors of Submission pap173s1
// Distributed under the FreeBSD license (see license.txt)
//#####################################################################
#include "Norm_Helper_CUDA.h"
using namespace SPGrid;

#define THREADBLOCK 1024

using namespace SPGrid;
//#####################################################################
// Kernel Norm_Kernel
//#####################################################################
template<class T,int number_of_elements_per_block,int page_size,class T_offset_ptr> 
__global__ void Norm_Kernel(const T* const u,T* buffer,unsigned int number_of_blocks){
    const int span = THREADBLOCK / number_of_elements_per_block;
    const int block = threadIdx.x / number_of_elements_per_block;
    const unsigned long entry = threadIdx.x % number_of_elements_per_block;
    __shared__ double cache[THREADBLOCK];
    unsigned int block_id = blockIdx.x * span + block;
    unsigned long data_offset = ((unsigned long)block_id * page_size + entry * sizeof(T));
    if (block_id < number_of_blocks){
        cache[threadIdx.x] = abs(*(reinterpret_cast<T*>((unsigned long)u + data_offset)));}
    else{cache[threadIdx.x] = 0;}
    __syncthreads();
    int size=THREADBLOCK;
    while(size!=1){
        size=size>>1;
        if(threadIdx.x<size) cache[threadIdx.x] = (cache[threadIdx.x+size]<cache[threadIdx.x])?cache[threadIdx.x]:cache[threadIdx.x+size];
        __syncthreads();}
    buffer[blockIdx.x]=cache[0];
}
//#####################################################################
// Function Run
//#####################################################################
template <class T,int log2_struct,class T_offset_ptr> 
T Norm_Helper_CUDA<T,log2_struct,3,T_offset_ptr>::Run()
{
    enum{number_of_elements_per_block=block_xsize*block_ysize*block_zsize};
    const int span = THREADBLOCK / number_of_elements_per_block;
    unsigned int number_of_cuda_blocks = (size%span)?(size/span+1):(size/span);
    if(number_of_cuda_blocks == 0) return 0;
    T* buffer=NULL;
    if(cudaMalloc((void**)&buffer,number_of_cuda_blocks*sizeof(T))!=cudaSuccess) abort();
    Norm_Kernel<T,number_of_elements_per_block,page_size,T_offset_ptr>
        <<<number_of_cuda_blocks,THREADBLOCK,0>>>
        (u,buffer,size);
    T* host_buffer=NULL;
    if(cudaMallocHost((void**)&host_buffer,number_of_cuda_blocks*sizeof(T))!=cudaSuccess) abort();
    cudaMemcpy(host_buffer,buffer,number_of_cuda_blocks*sizeof(T),cudaMemcpyDeviceToHost);
    T result=0;
    for(int i = 0;i < number_of_cuda_blocks;++i) result=(result>std::abs(host_buffer[i]))?result:std::abs(host_buffer[i]);
    cudaFree(buffer);
    cudaFreeHost(host_buffer);
    return result;
}
//#####################################################################
template class Norm_Helper_CUDA<float,4,3,unsigned int>;
