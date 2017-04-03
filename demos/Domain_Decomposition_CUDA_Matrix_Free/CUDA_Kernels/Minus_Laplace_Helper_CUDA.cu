//#####################################################################
//Removed For Anonymity: Copyright Authors of Submission pap173s1
// Distributed under the FreeBSD license (see license.txt)
//#####################################################################
#include <iostream>
#include "Minus_Laplace_Helper_CUDA.h"
using namespace SPGrid;
#define THREADBLOCK 256
#define PREFETCH 256

template <class T, int log2_struct,int xsize,int ysize,int zsize,class T_offset_ptr,bool accumulative>
__global__ void Minus_Laplace_Kernel_3D(const unsigned* masks, T* out,const T* in,
                                        const T_offset_ptr* offsets,
                                        const T_offset_ptr* const b_x_minus,
                                        const T_offset_ptr* const b_x_plus,
                                        const T_offset_ptr* const b_y_minus,
                                        const T_offset_ptr* const b_y_plus,
                                        const T_offset_ptr* const b_z_minus,
                                        const T_offset_ptr* const b_z_plus,
                                        int max_block,const unsigned flag_to_check){
    enum {
        DATABLOCK=xsize*ysize*zsize,
        span=THREADBLOCK/DATABLOCK,
        xstride=ysize*zsize,
        ystride=zsize,
        zstride=1  
    };
    const unsigned int block = threadIdx.x / DATABLOCK;
    const unsigned int entry = threadIdx.x % DATABLOCK;
    const unsigned int z = entry % zsize;
    const unsigned int y = entry / zsize % ysize;
    const unsigned int x = entry / zsize / ysize;
    typedef SPGrid_Mask<log2_struct, NextLogTwo<sizeof(T)>::value,3> T_MASK;

    __shared__ T_offset_ptr block_index[PREFETCH];
    __shared__ T_offset_ptr block_minus_x_index[PREFETCH];
    __shared__ T_offset_ptr block_plus_x_index[PREFETCH];
    __shared__ T_offset_ptr block_minus_y_index[PREFETCH];
    __shared__ T_offset_ptr block_plus_y_index[PREFETCH];
    __shared__ T_offset_ptr block_minus_z_index[PREFETCH];
    __shared__ T_offset_ptr block_plus_z_index[PREFETCH];
    if(threadIdx.x < PREFETCH)
        if(blockIdx.x * PREFETCH + threadIdx.x < max_block){
            block_index[threadIdx.x] = offsets[blockIdx.x * PREFETCH + threadIdx.x];
            block_minus_x_index[threadIdx.x]=b_x_minus[blockIdx.x * PREFETCH + threadIdx.x];
            block_plus_x_index[threadIdx.x]=b_x_plus[blockIdx.x * PREFETCH + threadIdx.x];
            block_minus_y_index[threadIdx.x]=b_y_minus[blockIdx.x * PREFETCH + threadIdx.x];
            block_plus_y_index[threadIdx.x]=b_y_plus[blockIdx.x * PREFETCH + threadIdx.x];
            block_minus_z_index[threadIdx.x]=b_z_minus[blockIdx.x * PREFETCH + threadIdx.x];
            block_plus_z_index[threadIdx.x]=b_z_plus[blockIdx.x * PREFETCH + threadIdx.x];
        }
    __syncthreads();

    for(int i = 0;i < PREFETCH;i += span){
        if (blockIdx.x * PREFETCH + i + block < max_block){
            T* out_base = reinterpret_cast<T*>((unsigned long)out + (unsigned long)block_index[i + block]); 
            unsigned mask_value = reinterpret_cast<unsigned*>((unsigned long)masks + (unsigned long)block_index[i + block])[entry];
            if(mask_value & flag_to_check){             
                T center_value = reinterpret_cast<T*>((unsigned long)in + (unsigned long)block_index[i + block])[entry];
                
                T& x_minus_value = (x==0)
                    ? reinterpret_cast<T*>((unsigned long)in + (unsigned long)block_minus_x_index[i + block])[entry+(xsize-1)*xstride]
                    : reinterpret_cast<T*>((unsigned long)in + (unsigned long)block_index[i + block])[entry-xstride];
                T& x_plus_value = (x==xsize-1)
                    ? reinterpret_cast<T*>((unsigned long)in + (unsigned long)block_plus_x_index[i + block])[entry-(xsize-1)*xstride]
                    : reinterpret_cast<T*>((unsigned long)in + (unsigned long)block_index[i + block])[entry+xstride];
                                
                T& y_minus_value = (y==0)
                    ? reinterpret_cast<T*>((unsigned long)in + (unsigned long)block_minus_y_index[i + block])[entry+(ysize-1)*ystride]
                    : reinterpret_cast<T*>((unsigned long)in + (unsigned long)block_index[i + block])[entry-ystride];
                T& y_plus_value = (y==ysize-1)
                    ? reinterpret_cast<T*>((unsigned long)in + (unsigned long)block_plus_y_index[i + block])[entry-(ysize-1)*ystride]
                    : reinterpret_cast<T*>((unsigned long)in + (unsigned long)block_index[i + block])[entry+ystride];
                    
                T& z_minus_value = (z==0)
                    ? reinterpret_cast<T*>((unsigned long)in + (unsigned long)block_minus_z_index[i + block])[entry+(zsize-1)*zstride]
                    : reinterpret_cast<T*>((unsigned long)in + (unsigned long)block_index[i + block])[entry-zstride];
                T& z_plus_value = (z==zsize-1)
                    ? reinterpret_cast<T*>((unsigned long)in + (unsigned long)block_plus_z_index[i + block])[entry-(zsize-1)*zstride]
                    : reinterpret_cast<T*>((unsigned long)in + (unsigned long)block_index[i + block])[entry+zstride];
                    
                T result=0;
                if(accumulative) result = out_base[entry];
                
                T x_minus = (x_minus_value - center_value);
                T x_plus = (x_plus_value - center_value);
                T y_minus = (y_minus_value - center_value);
                T y_plus = (y_plus_value - center_value);
                T z_minus = (z_minus_value - center_value);
                T z_plus = (z_plus_value - center_value);

                if (mask_value & (SPGrid_Solver_Face_Minus_X_Active))
                    result += x_minus;
                if (mask_value & (SPGrid_Solver_Face_Plus_X_Active))
                    result += x_plus; 

                if (mask_value & (SPGrid_Solver_Face_Minus_Y_Active))
                    result += y_minus;
                if (mask_value & (SPGrid_Solver_Face_Plus_Y_Active))
                    result += y_plus;
 
                if (mask_value & (SPGrid_Solver_Face_Minus_Z_Active))
                    result += z_minus; 
                if (mask_value & (SPGrid_Solver_Face_Plus_Z_Active))
                    result += z_plus; 
                
                //if(result!=0) printf("here,%f\n",result);
                out_base[entry] = result;
            }
        }
    }
}
//#####################################################################
// Constructor 3D
//#####################################################################
template <class T, int log2_struct,class T_offset_ptr,bool accumulative>
Minus_Laplace_Helper_CUDA<T,log2_struct,3,T_offset_ptr,accumulative>::Minus_Laplace_Helper_CUDA(T* const x_input,const T* const y_input,const unsigned* const mask_input,
                                                                                                const T_offset_ptr* const b_input,
                                                                                                const T_offset_ptr* const b_x_minus_input,
                                                                                                const T_offset_ptr* const b_x_plus_input,
                                                                                                const T_offset_ptr* const b_y_minus_input,
                                                                                                const T_offset_ptr* const b_y_plus_input,
                                                                                                const T_offset_ptr* const b_z_minus_input,
                                                                                                const T_offset_ptr* const b_z_plus_input,
                                                                                                const int size_input,const unsigned flag_to_check_input)
:x(x_input),y(y_input),mask(mask_input),
    b(b_input),
    b_x_minus(b_x_minus_input),b_x_plus(b_x_plus_input),
    b_y_minus(b_y_minus_input),b_y_plus(b_y_plus_input),
    b_z_minus(b_z_minus_input),b_z_plus(b_z_plus_input),
    size(size_input),flag_to_check(flag_to_check_input)
{
}
//#####################################################################
// Function Run_Index_Range
//#####################################################################
// T_MASK corresponds to the mask for the data (not the mask channel)
template <class T,int log2_struct,class T_offset_ptr,bool accumulative> 
void Minus_Laplace_Helper_CUDA<T,log2_struct,3,T_offset_ptr,accumulative>::Run_Index_Range(const unsigned int index_start,
                                                                                           const unsigned int index_end,
                                                                                           cudaStream_t& cuda_stream)
{
    int number_of_blocks=index_end-index_start+1;
    int number_of_cuda_blocks = (number_of_blocks%PREFETCH)?(number_of_blocks/PREFETCH+1):(number_of_blocks/PREFETCH);
    if(number_of_cuda_blocks == 0) return;
    Minus_Laplace_Kernel_3D<T,log2_struct,block_xsize,block_ysize,block_zsize,T_offset_ptr,accumulative>
        <<<number_of_cuda_blocks,THREADBLOCK,0,cuda_stream>>>
        (mask,x,y,
         b+index_start,
         b_x_minus+index_start,b_x_plus+index_start,
         b_y_minus+index_start,b_y_plus+index_start,
         b_z_minus+index_start,b_z_plus+index_start,
         number_of_blocks,flag_to_check);

    // cudaDeviceSynchronize();
    // cudaError err = cudaGetLastError();
    // if(cudaSuccess != err){
    //     std::cerr << "Error in Minus Laplace Helper. Msg: "<< cudaGetErrorString(err) << std::endl;
    //     abort();
    // }
}
//#####################################################################
template class Minus_Laplace_Helper_CUDA<float,4,3,unsigned int,true>;
