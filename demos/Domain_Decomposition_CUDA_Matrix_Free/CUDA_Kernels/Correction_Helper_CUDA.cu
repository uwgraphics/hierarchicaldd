//#####################################################################
//Removed For Anonymity: Copyright Authors of Submission pap173s1
// Distributed under the FreeBSD license (see license.txt)
//#####################################################################
#include "Correction_Helper_CUDA.h"
using namespace SPGrid;

namespace{
template<unsigned d> struct BitCount;
template<> struct BitCount<0>  {enum {value=0};};
template<unsigned d> struct BitCount {enum {value=(d&1)+BitCount<(d>>1)>::value};};
__device__ float Dinv[64] = {  0.f, 1.f/(float)BitCount< 1>::value, 1.f/(float)BitCount< 2>::value, 1.f/(float)BitCount< 3>::value, 
    1.f/(float)BitCount< 4>::value, 1.f/(float)BitCount< 5>::value, 1.f/(float)BitCount< 6>::value, 1.f/(float)BitCount< 7>::value, 
    1.f/(float)BitCount< 8>::value, 1.f/(float)BitCount< 9>::value, 1.f/(float)BitCount<10>::value, 1.f/(float)BitCount<11>::value, 
    1.f/(float)BitCount<12>::value, 1.f/(float)BitCount<13>::value, 1.f/(float)BitCount<14>::value, 1.f/(float)BitCount<15>::value, 
    1.f/(float)BitCount<16>::value, 1.f/(float)BitCount<17>::value, 1.f/(float)BitCount<18>::value, 1.f/(float)BitCount<19>::value,        
    1.f/(float)BitCount<20>::value, 1.f/(float)BitCount<21>::value, 1.f/(float)BitCount<22>::value, 1.f/(float)BitCount<23>::value, 
    1.f/(float)BitCount<24>::value, 1.f/(float)BitCount<25>::value, 1.f/(float)BitCount<26>::value, 1.f/(float)BitCount<27>::value, 
    1.f/(float)BitCount<28>::value, 1.f/(float)BitCount<29>::value, 1.f/(float)BitCount<30>::value, 1.f/(float)BitCount<31>::value, 
    1.f/(float)BitCount<32>::value, 1.f/(float)BitCount<33>::value, 1.f/(float)BitCount<34>::value, 1.f/(float)BitCount<35>::value, 
    1.f/(float)BitCount<36>::value, 1.f/(float)BitCount<37>::value, 1.f/(float)BitCount<38>::value, 1.f/(float)BitCount<39>::value,
    1.f/(float)BitCount<40>::value, 1.f/(float)BitCount<41>::value, 1.f/(float)BitCount<42>::value, 1.f/(float)BitCount<43>::value, 
    1.f/(float)BitCount<44>::value, 1.f/(float)BitCount<45>::value, 1.f/(float)BitCount<46>::value, 1.f/(float)BitCount<47>::value, 
    1.f/(float)BitCount<48>::value, 1.f/(float)BitCount<49>::value, 1.f/(float)BitCount<50>::value, 1.f/(float)BitCount<51>::value, 
    1.f/(float)BitCount<52>::value, 1.f/(float)BitCount<53>::value, 1.f/(float)BitCount<54>::value, 1.f/(float)BitCount<55>::value, 
    1.f/(float)BitCount<56>::value, 1.f/(float)BitCount<57>::value, 1.f/(float)BitCount<58>::value, 1.f/(float)BitCount<59>::value,
    1.f/(float)BitCount<60>::value, 1.f/(float)BitCount<61>::value, 1.f/(float)BitCount<62>::value, 1.f/(float)BitCount<63>::value
};
}
#define THREADBLOCK 256
#define PREFETCH 256

using namespace SPGrid;
//#####################################################################
// Kernel Streaming_Correction_Kernel
//#####################################################################
template<class T,int number_of_elements_per_block,int page_size,class T_offset_ptr> 
__global__ void Streaming_Correction_Kernel(const unsigned* const flags,T* const u, const T* const r,T omega,unsigned int number_of_blocks){
    const int span = THREADBLOCK / number_of_elements_per_block;
    const int block = threadIdx.x / number_of_elements_per_block;
    const unsigned long entry = threadIdx.x % number_of_elements_per_block;
    __shared__ float Dinv_buffered[64];
    if(threadIdx.x < 64)
        Dinv_buffered[threadIdx.x] = Dinv[threadIdx.x];
    __syncthreads();
    for(int i = 0;i < PREFETCH;i += span){
        unsigned int block_id = blockIdx.x * PREFETCH + i + block;
        unsigned long data_offset = ((unsigned long)block_id * page_size + entry * sizeof(T));
        if (block_id < number_of_blocks){
            unsigned int cell_flag = *(reinterpret_cast<unsigned*>((unsigned long)flags + data_offset));
            if(cell_flag & SPGrid_Solver_Cell_Type_Active){
                cell_flag = (cell_flag >> (BitLength<SPGrid_Solver_Face_Minus_X_Active>::value-1))&0x3f;
                *(reinterpret_cast<T*>((unsigned long)u + data_offset)) = *(reinterpret_cast<T*>((unsigned long)u + data_offset)) + 
                    *(reinterpret_cast<T*>((unsigned long)r + data_offset)) * omega * Dinv_buffered[cell_flag];
            }
        }
    }
}
//#####################################################################
// Kernel Blocked_Boundary_Correction_Kernel
//#####################################################################
template<class T,int number_of_elements_per_block,class T_offset_ptr> 
__global__ void Blocked_Boundary_Correction_Kernel(const unsigned* const flags,T* const u, const T* const r,T omega,unsigned int number_of_blocks,const T_offset_ptr* offsets){
    const int span = THREADBLOCK / number_of_elements_per_block;
    const int block = threadIdx.x / number_of_elements_per_block;
    const unsigned long entry = threadIdx.x % number_of_elements_per_block;
    __shared__ float Dinv_buffered[64];
    __shared__ T_offset_ptr block_index[PREFETCH];
    if(threadIdx.x < 64)
        Dinv_buffered[threadIdx.x] = Dinv[threadIdx.x];
    if(threadIdx.x < PREFETCH)
        if(blockIdx.x * PREFETCH + threadIdx.x < number_of_blocks)
            block_index[threadIdx.x] = offsets[blockIdx.x * PREFETCH + threadIdx.x];
    __syncthreads();
    for(int i = 0;i < PREFETCH;i += span){
        if (blockIdx.x * PREFETCH + i + block < number_of_blocks){
            unsigned int cell_flag = reinterpret_cast<unsigned int*>((unsigned long)flags + block_index[i+block])[entry];
            if(cell_flag & SPGrid_Solver_Cell_Type_Boundary){
                cell_flag = (cell_flag >> (BitLength<SPGrid_Solver_Face_Minus_X_Active>::value-1))&0x3f;
                T* u_ptr = &(reinterpret_cast<T*>((unsigned long)u + block_index[i+block])[entry]);
                const T* r_ptr = &(reinterpret_cast<const T*>((unsigned long)r + block_index[i+block])[entry]);
                *(u_ptr) = *(u_ptr) + *(r_ptr)*omega*Dinv_buffered[cell_flag];
            }
        }
    }
}
//#####################################################################
// Kernel Blocked_Interior_Correction_Kernel
//#####################################################################
template<class T,int number_of_elements_per_block,class T_offset_ptr> 
__global__ void Blocked_Interior_Correction_Kernel(const unsigned* const flags,T* const u, const T* const r,T omega,unsigned int number_of_blocks,const T_offset_ptr* offsets){
    const int span = THREADBLOCK / number_of_elements_per_block;
    const int block = threadIdx.x / number_of_elements_per_block;
    const unsigned long entry = threadIdx.x % number_of_elements_per_block;
    __shared__ float Dinv_buffered[64];
    __shared__ T_offset_ptr block_index[PREFETCH];
    if(threadIdx.x < 64)
        Dinv_buffered[threadIdx.x] = Dinv[threadIdx.x];
    if(threadIdx.x < PREFETCH)
        if(blockIdx.x * PREFETCH + threadIdx.x < number_of_blocks)
            block_index[threadIdx.x] = offsets[blockIdx.x * PREFETCH + threadIdx.x];
    __syncthreads();
    for(int i = 0;i < PREFETCH;i += span){
        if (blockIdx.x * PREFETCH + i + block < number_of_blocks){
            unsigned int cell_flag = reinterpret_cast<unsigned int*>((unsigned long)flags + block_index[i+block])[entry];
            if((cell_flag & SPGrid_Solver_Cell_Type_Active) && (!(cell_flag & SPGrid_Solver_Cell_Type_Boundary))){
                cell_flag = (cell_flag >> (BitLength<SPGrid_Solver_Face_Minus_X_Active>::value-1))&0x3f;
                T* u_ptr = &(reinterpret_cast<T*>((unsigned long)u + block_index[i+block])[entry]);
                const T* r_ptr = &(reinterpret_cast<const T*>((unsigned long)r + block_index[i+block])[entry]);
                *(u_ptr) = *(u_ptr) + *(r_ptr)*omega*Dinv_buffered[cell_flag];
            }
        }
    }
}
//#####################################################################
// Function Run_Index_Range
//#####################################################################
// T_MASK corresponds to the mask for the data (not the mask channel)
template <class T,int log2_struct,class T_offset_ptr> 
void Correction_Helper_CUDA<T,log2_struct,3,T_offset_ptr>::Run_Index_Range(const unsigned int index_start,
                                                                           const unsigned int index_end,
                                                                           cudaStream_t& cuda_stream)
{
    enum{page_size=4096};
    if((1<<log2_struct)*block_xsize*block_ysize*block_zsize!=page_size){std::cerr<<"Wrong Page Size!"<<std::endl;abort();}
    int number_of_blocks=index_end-index_start+1;
    int number_of_cuda_blocks = (number_of_blocks%PREFETCH)?(number_of_blocks/PREFETCH+1):(number_of_blocks/PREFETCH);
    if(number_of_cuda_blocks == 0) return;
    Streaming_Correction_Kernel<T,block_xsize*block_ysize*block_zsize,page_size,T_offset_ptr>
        <<<number_of_cuda_blocks,THREADBLOCK,0,cuda_stream>>>
        (mask,u,r,omega,number_of_blocks);
    //cudaError err = cudaDeviceSynchronize(); 
    //if(err!=cudaSuccess) std::cout<<"Something went wrong in correction kernel!"<<std::endl;
}
//#####################################################################
// Function Run_Boundary_Blocks
//#####################################################################
template <class T,int log2_struct,class T_offset_ptr> 
void Correction_Helper_CUDA<T,log2_struct,3,T_offset_ptr>::Run_Boundary_Blocks(const T_offset_ptr* offsets,
                                                                               const unsigned int number_of_blocks,
                                                                               cudaStream_t& cuda_stream)
{
    int number_of_cuda_blocks = (number_of_blocks%PREFETCH)?(number_of_blocks/PREFETCH+1):(number_of_blocks/PREFETCH);
    if(number_of_cuda_blocks == 0) return;
    Blocked_Boundary_Correction_Kernel<T,block_xsize*block_ysize*block_zsize,T_offset_ptr>
        <<<number_of_cuda_blocks,THREADBLOCK,0,cuda_stream>>>
        (mask,u,r,omega,number_of_blocks,offsets);
    //cudaError err = cudaDeviceSynchronize(); 
    //if(err!=cudaSuccess) std::cout<<"Something went wrong in correction boundary kernel!"<<std::endl;
}
//#####################################################################
// Function Run_Interior_Blocks
//#####################################################################
template <class T,int log2_struct,class T_offset_ptr> 
void Correction_Helper_CUDA<T,log2_struct,3,T_offset_ptr>::Run_Interior_Blocks(const T_offset_ptr* offsets,
                                                                               const unsigned int number_of_blocks,
                                                                               cudaStream_t& cuda_stream)
{
    int number_of_cuda_blocks = (number_of_blocks%PREFETCH)?(number_of_blocks/PREFETCH+1):(number_of_blocks/PREFETCH);
    if(number_of_cuda_blocks == 0) return;

    Blocked_Interior_Correction_Kernel<T,block_xsize*block_ysize*block_zsize,T_offset_ptr>
        <<<number_of_cuda_blocks,THREADBLOCK,0,cuda_stream>>>
        (mask,u,r,omega,number_of_blocks,offsets);

    //cudaError err = cudaDeviceSynchronize(); 
    //if(err!=cudaSuccess) std::cout<<"Something went wrong in correction interior kernel!"<<std::endl;
}
//#####################################################################
template class Correction_Helper_CUDA<float,4,3,unsigned int>;
