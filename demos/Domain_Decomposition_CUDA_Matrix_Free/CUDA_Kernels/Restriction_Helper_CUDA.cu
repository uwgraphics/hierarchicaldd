//#####################################################################
//Removed For Anonymity: Copyright Authors of Submission pap173s1
// Distributed under the FreeBSD license (see license.txt)
//#####################################################################
#include "Restriction_Helper_CUDA.h"
#include <stdio.h>
using namespace SPGrid;
#define THREADBLOCK 256
#define PREFETCH 32
template<class T,int number_of_elements_per_block,class T_offset_ptr,int block_xsize,int block_ysize,int block_zsize> 
__global__ void Restriction_Kernel_3D(const unsigned* const flags,T* const b, const T* const r,const T_offset_ptr* const offset_fine_array,T_offset_ptr* const offset_coarse,const unsigned int number_of_blocks){
    //printf("%d\n",(int)offset_fine_array[0]);
    //this kernel assumes we have more threads per block than entry per block
    enum{span = THREADBLOCK / number_of_elements_per_block};
    //if(span == 0) {printf("error, threads per block is too small!");return;}
    const T weight_table_1D[4] = {1.f/8.f,3.f/8.f,3.f/8.f,1.f/8.f};
    const int block = threadIdx.x / number_of_elements_per_block;
    const unsigned long entry = threadIdx.x % number_of_elements_per_block;
    const int coarse_z = entry % block_zsize;
    const int coarse_y = entry / block_zsize % block_ysize;
    const int coarse_x = entry / block_zsize / block_ysize;
    __shared__ T_offset_ptr coarse_offsets[PREFETCH];
    __shared__ T_offset_ptr fine_offsets[PREFETCH*27];
    if(threadIdx.x<PREFETCH){
        for(int i = 0;i < 27;++i){
            unsigned int fine_offset_position=blockIdx.x*27*PREFETCH+PREFETCH*i+threadIdx.x;
            if(fine_offset_position<number_of_blocks*27){
                fine_offsets[PREFETCH*i+threadIdx.x] = offset_fine_array[fine_offset_position];}}        
        if(blockIdx.x*PREFETCH+threadIdx.x<number_of_blocks){
            coarse_offsets[threadIdx.x] = offset_coarse[blockIdx.x*PREFETCH+threadIdx.x];}
    }
    __syncthreads();

    for(int i = 0;i < PREFETCH;i += span){
        const int block_index = i + block;
        const int fine_block_index = block_index * 27;
        if(block_index+blockIdx.x*PREFETCH<number_of_blocks){
            
    
            const unsigned flag = reinterpret_cast<const unsigned*>((unsigned long)flags+(unsigned long)coarse_offsets[block_index])[entry];
            if(flag & SPGrid_Solver_Cell_Type_Active){
                T tmp = 0;
                int fine_x = (coarse_x - 1) * 2;
                for(int i = 0;i < 4;++i,++fine_x){
                    int fine_y = (coarse_y - 1) * 2;
                    for(int j = 0;j < 4;++j,++fine_y){
                        int fine_z = (coarse_z - 1) * 2;
                        for(int k = 0;k < 4;++k,++fine_z){
                            T weight = 1.f;
                            weight *= weight_table_1D[i];
                            weight *= weight_table_1D[j];
                            weight *= weight_table_1D[k];
                            int fine_block = fine_block_index;
                            int fine_entry = 0;
                            if(fine_z < 0){
                                fine_entry = fine_z + block_zsize;
                            }else if(fine_z < block_zsize){
                                fine_entry = fine_z;
                                fine_block += 1;
                            }else{
                                fine_entry = fine_z - block_zsize;
                                fine_block += 2;}
                            
                            if(fine_y < 0){
                                fine_entry += (fine_y + block_ysize) * block_zsize;
                            }else if(fine_y < block_ysize){
                                fine_entry += fine_y * block_zsize;
                                fine_block += 3;
                            }else{
                                fine_entry += (fine_y - block_ysize) * block_zsize;
                                fine_block += 6;}
                            
                            if(fine_x < 0){
                                fine_entry += (fine_x + block_xsize) * block_zsize * block_ysize;
                            }else if(fine_x < block_xsize){
                                fine_entry += fine_x * block_zsize * block_ysize;
                                fine_block += 9;
                            }else{
                                fine_entry += (fine_x - block_xsize) * block_zsize * block_ysize;
                                fine_block += 18;}
                            
                            if(fine_offsets[fine_block] == 0xdeadbeef) continue;
                            tmp+=weight*(reinterpret_cast<const T*>((unsigned long)r+(unsigned long)fine_offsets[fine_block])[fine_entry]);
                            // if(reinterpret_cast<const T*>((unsigned long)r+(unsigned long)fine_offsets[fine_block])[fine_entry] != 0 && 
                            //    reinterpret_cast<const T*>((unsigned long)r+(unsigned long)fine_offsets[fine_block])[fine_entry] != 1){
                            //     printf("%f\n",reinterpret_cast<const T*>((unsigned long)r+(unsigned long)fine_offsets[fine_block])[fine_entry]);
                            //}
                        }
                    }
                }
                reinterpret_cast<T*>((unsigned long)b+(unsigned long)coarse_offsets[block_index])[entry] = tmp;
            }
        }
    }
}
//#####################################################################
// Function Run_Index_Range
//#####################################################################
// T_MASK corresponds to the mask for the data (not the mask channel)
template <class T,int log2_struct,class T_offset_ptr> 
void Restriction_Helper_CUDA<T,log2_struct,3,T_offset_ptr>::Run_Index_Range(const unsigned int index_start,const unsigned int index_end,cudaStream_t& cuda_stream)
{
    int number_of_blocks=index_end-index_start+1;
    int number_of_cuda_blocks = (number_of_blocks%PREFETCH)?(number_of_blocks/PREFETCH+1):(number_of_blocks/PREFETCH);
    enum{span = THREADBLOCK / block_xsize*block_ysize*block_zsize};
    if(span == 0) {printf("error, threads per block is too small!");abort();}
    if(number_of_cuda_blocks == 0) return;

    Restriction_Kernel_3D<T,block_xsize*block_ysize*block_zsize,T_offset_ptr,block_xsize,block_ysize,block_zsize>
        <<<number_of_cuda_blocks,THREADBLOCK,0,cuda_stream>>>
        (flags,b,r,offset_fine_array,offset_coarse,number_of_blocks);

    // cudaError err = cudaGetLastError();
    // cudaDeviceSynchronize(); 
    // if(err!=cudaSuccess) std::cout<<"Something went wrong in restriction kernel!"<<std::endl;
}

template class Restriction_Helper_CUDA<float,4,3,unsigned int>;

