//#####################################################################
//Removed For Anonymity: Copyright Authors of Submission pap173s1
// Distributed under the FreeBSD license (see license.txt)
//#####################################################################
#include "Prolongation_Helper_CUDA.h"
#include <stdio.h>
using namespace SPGrid;
#define THREADBLOCK 1024
#define PREFETCH 1024
template<class T,int number_of_elements_per_block,class T_offset_ptr,int block_xsize,int block_ysize,int block_zsize> 
__global__ void Prolongation_Kernel_3D(const unsigned* const flags,T* const u_fine, const T* const u_coarse,T_offset_ptr* offset_fine_array,T_offset_ptr* offset_coarse_array,unsigned int number_of_blocks){
    //this kernel assumes we have more threads per block than entry per block
    enum{span = THREADBLOCK / number_of_elements_per_block};
    //if(span == 0) {printf("error, threads per block is too small!");return;}
    const int block = threadIdx.x / number_of_elements_per_block;
    const unsigned long entry = threadIdx.x % number_of_elements_per_block;
    __shared__ int fine_z[THREADBLOCK];
    __shared__ int fine_y[THREADBLOCK];
    __shared__ int fine_x[THREADBLOCK];
    fine_z[threadIdx.x] = entry % block_zsize;
    fine_y[threadIdx.x] = entry / block_zsize % block_ysize;
    fine_x[threadIdx.x] = entry / block_zsize / block_ysize;
    const int cell_parity_z = fine_z[threadIdx.x] % 2;
    const int cell_parity_y = fine_y[threadIdx.x] % 2;
    const int cell_parity_x = fine_x[threadIdx.x] % 2;
   
    __shared__ int coarse_z[THREADBLOCK];
    __shared__ int coarse_y[THREADBLOCK];
    __shared__ int coarse_x[THREADBLOCK];
    __shared__ T_offset_ptr fine_offsets[PREFETCH];
    __shared__ T_offset_ptr coarse_offsets[PREFETCH];
    if(threadIdx.x < PREFETCH)
        if(blockIdx.x * PREFETCH + threadIdx.x < number_of_blocks){
            fine_offsets[threadIdx.x] = offset_fine_array[blockIdx.x * PREFETCH + threadIdx.x];
            coarse_offsets[threadIdx.x] = offset_coarse_array[blockIdx.x * PREFETCH + threadIdx.x];}
    __syncthreads();
    for(int i = 0;i < PREFETCH;i += span){
        const int block_index = i + block;
        const int block_index_coarse = (block_index/8)*8;
        if(block_index + blockIdx.x * PREFETCH  < number_of_blocks){
            if((unsigned long)fine_offsets[block_index] == 0xdeadbeef) continue;
            const unsigned flag = reinterpret_cast<const unsigned*>((unsigned long)flags+(unsigned long)fine_offsets[block_index])[entry];
            if(flag & SPGrid_Solver_Cell_Type_Active){
                coarse_z[threadIdx.x] = (block_index & 0x1)?(fine_z[threadIdx.x]+block_zsize)/2:(fine_z[threadIdx.x])/2;
                coarse_y[threadIdx.x] = (block_index & 0x2)?(fine_y[threadIdx.x]+block_ysize)/2:(fine_y[threadIdx.x])/2;
                coarse_x[threadIdx.x] = (block_index & 0x4)?(fine_x[threadIdx.x]+block_xsize)/2:(fine_x[threadIdx.x])/2;
                T tmp = 0;
                for(int parity_x = 0;parity_x <= 1;++parity_x)
                for(int parity_y = 0;parity_y <= 1;++parity_y)
                for(int parity_z = 0;parity_z <= 1;++parity_z){
                    T weight = 4.f; //weight is multiplied for 4 for the delta x^2 terms in laplace
                    weight *= (cell_parity_z == parity_z)?(3.f/4.f):(1.f/4.f);
                    weight *= (cell_parity_y == parity_y)?(3.f/4.f):(1.f/4.f);
                    weight *= (cell_parity_x == parity_x)?(3.f/4.f):(1.f/4.f);
                    int coarse_z_cell_index = coarse_z[threadIdx.x]+parity_z;
                    int coarse_y_cell_index = coarse_y[threadIdx.x]+parity_y;
                    int coarse_x_cell_index = coarse_x[threadIdx.x]+parity_x;
                    int coarse_block = block_index_coarse;
                    if(coarse_z_cell_index==block_zsize){
                        coarse_z_cell_index = 0;
                        coarse_block += 0x1;}
                    if(coarse_y_cell_index==block_ysize){
                        coarse_y_cell_index = 0;
                        coarse_block += 0x2;}
                    if(coarse_x_cell_index==block_xsize){
                        coarse_x_cell_index = 0;
                        coarse_block += 0x4;}
                    if(coarse_offsets[coarse_block] == 0xdeadbeef) continue;
                    const int coarse_entry = coarse_z_cell_index + coarse_y_cell_index * block_zsize + coarse_x_cell_index * block_zsize * block_ysize;
                    tmp += weight * reinterpret_cast<T*>((unsigned long)u_coarse+(unsigned long)coarse_offsets[coarse_block])[coarse_entry];
                }
                T& u = reinterpret_cast<T*>((unsigned long)u_fine+(unsigned long)fine_offsets[block_index])[entry];
                u += tmp;
            }
        }
    }
}
//#####################################################################
// Function Run_Index_Range
//#####################################################################
// T_MASK corresponds to the mask for the data (not the mask channel)
template <class T,int log2_struct,class T_offset_ptr> 
void Prolongation_Helper_CUDA<T,log2_struct,3,T_offset_ptr>::Run_Index_Range(const unsigned int index_start,
                                                                             const unsigned int index_end,
                                                                             cudaStream_t& cuda_stream)
{
    //cudaEvent_t start,stop;
    //cudaEventCreate(&start);
    //cudaEventCreate(&stop);
    
    int number_of_blocks=(index_end-index_start+1)*8;//times 8 here...I have unflattened the std_array<*,8> here...
    int number_of_cuda_blocks = (number_of_blocks%PREFETCH)?(number_of_blocks/PREFETCH+1):(number_of_blocks/PREFETCH);
    enum{span = THREADBLOCK / block_xsize*block_ysize*block_zsize};
    if(span == 0) {printf("error, threads per block is too small!");abort();}
    if(PREFETCH < 8 || PREFETCH % 8 != 0) {std::cerr << "For Write Dependence Issue, Please Make Sure Prolongation Prefetch Is Multiply Of 8!"<<std::endl;abort();}
    if(number_of_cuda_blocks == 0) return;

    //cudaEventRecord(start);
    Prolongation_Kernel_3D<T,block_xsize*block_ysize*block_zsize,T_offset_ptr,block_xsize,block_ysize,block_zsize>
        <<<number_of_cuda_blocks,THREADBLOCK,0,cuda_stream>>>
        (flags,u_fine,u_coarse,offset_fine_array,offset_coarse_array,number_of_blocks);

    // cudaEventRecord(stop);
    // cudaEventSynchronize(stop);

    // float milliseconds = 0;
    // cudaEventElapsedTime(&milliseconds,start,stop);

    // double data_size = number_of_blocks * 9.0 * (double)4096;

    // std::cout << milliseconds/1000 << " sec for prolongation on " << data_size/2.0/1024.0/1024.0 << " MB of data" << std::endl;
    // std::cout << "bandwidth is: " << data_size/(milliseconds/1000)/1024.0/1024.0/1024.0 << " GB/s" << std::endl;


    // cudaError err = cudaGetLastError();
    // if(cudaSuccess != err){
    //     std::cerr << cudaGetErrorString(err) << std::endl;
    //     abort();
    // }
    // cudaDeviceSynchronize(); 
}

template class Prolongation_Helper_CUDA<float,4,3,unsigned int>;

