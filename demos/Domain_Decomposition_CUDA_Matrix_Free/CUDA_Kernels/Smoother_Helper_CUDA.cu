//#####################################################################
//Removed For Anonymity: Copyright Authors of Submission pap173s1
// Distributed under the FreeBSD license (see license.txt)
//#####################################################################
#include <iostream>
#include "Smoother_Helper_CUDA.h"
#include <SPGrid_Fluids/Solvers/Domain_Decomposition/SPGRID_MULTIGRID_FLAGS.h>
using namespace SPGrid;
#define THREADBLOCK 256
#define PREFETCH 256

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

template <class T, int log2_struct,int xsize,int ysize,int zsize,class T_offset_ptr>
__global__ void Boundary_Smoother_Kernel_3D(const unsigned* masks, T* r,T* u,const T* rhs,
                                            const T_offset_ptr* offsets,
                                            const T_offset_ptr* const b_x_minus,
                                            const T_offset_ptr* const b_x_plus,
                                            const T_offset_ptr* const b_y_minus,
                                            const T_offset_ptr* const b_y_plus,
                                            const T_offset_ptr* const b_z_minus,
                                            const T_offset_ptr* const b_z_plus,
                                            const int max_block,const T omega,
                                            const int iterations){
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
    __shared__ float Dinv_buffered[64];

    if(threadIdx.x < 64)
        Dinv_buffered[threadIdx.x] = Dinv[threadIdx.x];
    if(threadIdx.x < PREFETCH)
        if(threadIdx.x < max_block){
            block_index[threadIdx.x] = offsets[threadIdx.x];
            block_minus_x_index[threadIdx.x]=b_x_minus[threadIdx.x];
            block_plus_x_index[threadIdx.x]=b_x_plus[threadIdx.x];
            block_minus_y_index[threadIdx.x]=b_y_minus[threadIdx.x];
            block_plus_y_index[threadIdx.x]=b_y_plus[threadIdx.x];
            block_minus_z_index[threadIdx.x]=b_z_minus[threadIdx.x];
            block_plus_z_index[threadIdx.x]=b_z_plus[threadIdx.x];
        }
    for(int itr=0;itr<iterations;++itr){
        __syncthreads();
        for(int i = 0;i < PREFETCH;i += span){
            if (i + block < max_block){
                T* r_base = reinterpret_cast<T*>((unsigned long)r + (unsigned long)block_index[i + block]); 
                unsigned mask_value = reinterpret_cast<unsigned*>((unsigned long)masks + (unsigned long)block_index[i + block])[entry];
                if(mask_value & SPGrid_Solver_Cell_Type_Active){
                    T center_value = reinterpret_cast<T*>((unsigned long)u + (unsigned long)block_index[i + block])[entry];
                    T b_value = reinterpret_cast<T*>((unsigned long)rhs + (unsigned long)block_index[i + block])[entry];

                    T& x_minus_value = (x==0)
                        ? reinterpret_cast<T*>((unsigned long)u + (unsigned long)block_minus_x_index[i + block])[entry+(xsize-1)*xstride]
                        : reinterpret_cast<T*>((unsigned long)u + (unsigned long)block_index[i + block])[entry-xstride];
                    T& x_plus_value = (x==xsize-1)
                        ? reinterpret_cast<T*>((unsigned long)u + (unsigned long)block_plus_x_index[i + block])[entry-(xsize-1)*xstride]
                        : reinterpret_cast<T*>((unsigned long)u + (unsigned long)block_index[i + block])[entry+xstride];
                                
                    T& y_minus_value = (y==0)
                        ? reinterpret_cast<T*>((unsigned long)u + (unsigned long)block_minus_y_index[i + block])[entry+(ysize-1)*ystride]
                        : reinterpret_cast<T*>((unsigned long)u + (unsigned long)block_index[i + block])[entry-ystride];
                    T& y_plus_value = (y==ysize-1)
                        ? reinterpret_cast<T*>((unsigned long)u + (unsigned long)block_plus_y_index[i + block])[entry-(ysize-1)*ystride]
                        : reinterpret_cast<T*>((unsigned long)u + (unsigned long)block_index[i + block])[entry+ystride];
                    
                    T& z_minus_value = (z==0)
                        ? reinterpret_cast<T*>((unsigned long)u + (unsigned long)block_minus_z_index[i + block])[entry+(zsize-1)*zstride]
                        : reinterpret_cast<T*>((unsigned long)u + (unsigned long)block_index[i + block])[entry-zstride];
                    T& z_plus_value = (z==zsize-1)
                        ? reinterpret_cast<T*>((unsigned long)u + (unsigned long)block_plus_z_index[i + block])[entry-(zsize-1)*zstride]
                        : reinterpret_cast<T*>((unsigned long)u + (unsigned long)block_index[i + block])[entry+zstride];
                    
                    T result=0;
                    T x_minus = (x_minus_value - center_value);
                    T x_plus = (x_plus_value - center_value);
                    T y_minus = (y_minus_value - center_value);
                    T y_plus = (y_plus_value - center_value);
                    T z_minus = (z_minus_value - center_value);
                    T z_plus = (z_plus_value - center_value);

                    if (mask_value & (SPGrid_Solver_Face_Minus_X_Active))
                        result -= x_minus;
                    if (mask_value & (SPGrid_Solver_Face_Plus_X_Active))
                        result -= x_plus; 

                    if (mask_value & (SPGrid_Solver_Face_Minus_Y_Active))
                        result -= y_minus;
                    if (mask_value & (SPGrid_Solver_Face_Plus_Y_Active))
                        result -= y_plus;
 
                    if (mask_value & (SPGrid_Solver_Face_Minus_Z_Active))
                        result -= z_minus; 
                    if (mask_value & (SPGrid_Solver_Face_Plus_Z_Active))
                        result -= z_plus; 
                                   
                    r_base[entry] = b_value - result;
                }
            }
        }
        __syncthreads();
        for(int i = 0;i < PREFETCH;i += span){
            if (i + block < max_block){
                unsigned int cell_flag = reinterpret_cast<unsigned int*>((unsigned long)masks + block_index[i+block])[entry];
                if(cell_flag & SPGrid_Solver_Cell_Type_Boundary){
                    cell_flag = (cell_flag >> (BitLength<SPGrid_Solver_Face_Minus_X_Active>::value-1))&0x3f;
                    T* u_ptr = &(reinterpret_cast<T*>((unsigned long)u + block_index[i+block])[entry]);
                    const T* r_ptr = &(reinterpret_cast<const T*>((unsigned long)r + block_index[i+block])[entry]);
                    *(u_ptr) = *(u_ptr) + *(r_ptr)*omega*Dinv_buffered[cell_flag];
                }
            }
        }
    }
}
template <class T, int log2_struct,int xsize,int ysize,int zsize,class T_offset_ptr>
__global__ void Interior_Smoother_Kernel_3D(const unsigned* masks, T* r,T* u,const T* rhs,
                                            const T_offset_ptr* offsets,
                                            const T_offset_ptr* const b_x_minus,
                                            const T_offset_ptr* const b_x_plus,
                                            const T_offset_ptr* const b_y_minus,
                                            const T_offset_ptr* const b_y_plus,
                                            const T_offset_ptr* const b_z_minus,
                                            const T_offset_ptr* const b_z_plus,
                                            const int max_block,const T omega,
                                            const int iterations){
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
    __shared__ float Dinv_buffered[64];

    if(threadIdx.x < 64)
        Dinv_buffered[threadIdx.x] = Dinv[threadIdx.x];
    if(threadIdx.x < PREFETCH)
        if(threadIdx.x < max_block){
            block_index[threadIdx.x] = offsets[threadIdx.x];
            block_minus_x_index[threadIdx.x]=b_x_minus[threadIdx.x];
            block_plus_x_index[threadIdx.x]=b_x_plus[threadIdx.x];
            block_minus_y_index[threadIdx.x]=b_y_minus[threadIdx.x];
            block_plus_y_index[threadIdx.x]=b_y_plus[threadIdx.x];
            block_minus_z_index[threadIdx.x]=b_z_minus[threadIdx.x];
            block_plus_z_index[threadIdx.x]=b_z_plus[threadIdx.x];
        }
    for(int itr=0;itr<iterations;++itr){
        __syncthreads();
        for(int i = 0;i < PREFETCH;i += span){
            if (i + block < max_block){
                T* r_base = reinterpret_cast<T*>((unsigned long)r + (unsigned long)block_index[i + block]); 
                unsigned mask_value = reinterpret_cast<unsigned*>((unsigned long)masks + (unsigned long)block_index[i + block])[entry];
                if(mask_value & SPGrid_Solver_Cell_Type_Active){
                    T center_value = reinterpret_cast<T*>((unsigned long)u + (unsigned long)block_index[i + block])[entry];
                    T b_value = reinterpret_cast<T*>((unsigned long)rhs + (unsigned long)block_index[i + block])[entry];

                    T& x_minus_value = (x==0)
                        ? reinterpret_cast<T*>((unsigned long)u + (unsigned long)block_minus_x_index[i + block])[entry+(xsize-1)*xstride]
                        : reinterpret_cast<T*>((unsigned long)u + (unsigned long)block_index[i + block])[entry-xstride];
                    T& x_plus_value = (x==xsize-1)
                        ? reinterpret_cast<T*>((unsigned long)u + (unsigned long)block_plus_x_index[i + block])[entry-(xsize-1)*xstride]
                        : reinterpret_cast<T*>((unsigned long)u + (unsigned long)block_index[i + block])[entry+xstride];
                                
                    T& y_minus_value = (y==0)
                        ? reinterpret_cast<T*>((unsigned long)u + (unsigned long)block_minus_y_index[i + block])[entry+(ysize-1)*ystride]
                        : reinterpret_cast<T*>((unsigned long)u + (unsigned long)block_index[i + block])[entry-ystride];
                    T& y_plus_value = (y==ysize-1)
                        ? reinterpret_cast<T*>((unsigned long)u + (unsigned long)block_plus_y_index[i + block])[entry-(ysize-1)*ystride]
                        : reinterpret_cast<T*>((unsigned long)u + (unsigned long)block_index[i + block])[entry+ystride];
                    
                    T& z_minus_value = (z==0)
                        ? reinterpret_cast<T*>((unsigned long)u + (unsigned long)block_minus_z_index[i + block])[entry+(zsize-1)*zstride]
                        : reinterpret_cast<T*>((unsigned long)u + (unsigned long)block_index[i + block])[entry-zstride];
                    T& z_plus_value = (z==zsize-1)
                        ? reinterpret_cast<T*>((unsigned long)u + (unsigned long)block_plus_z_index[i + block])[entry-(zsize-1)*zstride]
                        : reinterpret_cast<T*>((unsigned long)u + (unsigned long)block_index[i + block])[entry+zstride];
                    
                    T result=0;
                    T x_minus = (x_minus_value - center_value);
                    T x_plus = (x_plus_value - center_value);
                    T y_minus = (y_minus_value - center_value);
                    T y_plus = (y_plus_value - center_value);
                    T z_minus = (z_minus_value - center_value);
                    T z_plus = (z_plus_value - center_value);

                    if (mask_value & (SPGrid_Solver_Face_Minus_X_Active))
                        result -= x_minus;
                    if (mask_value & (SPGrid_Solver_Face_Plus_X_Active))
                        result -= x_plus; 

                    if (mask_value & (SPGrid_Solver_Face_Minus_Y_Active))
                        result -= y_minus;
                    if (mask_value & (SPGrid_Solver_Face_Plus_Y_Active))
                        result -= y_plus;
 
                    if (mask_value & (SPGrid_Solver_Face_Minus_Z_Active))
                        result -= z_minus; 
                    if (mask_value & (SPGrid_Solver_Face_Plus_Z_Active))
                        result -= z_plus; 
                                   
                    r_base[entry] = b_value - result;
                }
            }
        }
        __syncthreads();
        for(int i = 0;i < PREFETCH;i += span){
            if (i + block < max_block){
                unsigned int cell_flag = reinterpret_cast<unsigned int*>((unsigned long)masks + block_index[i+block])[entry];
                if((cell_flag & SPGrid_Solver_Cell_Type_Active) && (!(cell_flag & SPGrid_Solver_Cell_Type_Boundary))){
                    cell_flag = (cell_flag >> (BitLength<SPGrid_Solver_Face_Minus_X_Active>::value-1))&0x3f;
                    T* u_ptr = &(reinterpret_cast<T*>((unsigned long)u + block_index[i+block])[entry]);
                    const T* r_ptr = &(reinterpret_cast<const T*>((unsigned long)r + block_index[i+block])[entry]);
                    *(u_ptr) = *(u_ptr) + *(r_ptr)*omega*Dinv_buffered[cell_flag];
                }
            }
        }
    }
}
template <class T, int log2_struct,int xsize,int ysize,int zsize,class T_offset_ptr>
__global__ void Bottom_Smoother_Kernel_3D(const unsigned* masks, T* r,T* u,const T* rhs,
                                          const T_offset_ptr* offsets,
                                          const T_offset_ptr* const b_x_minus,
                                          const T_offset_ptr* const b_x_plus,
                                          const T_offset_ptr* const b_y_minus,
                                          const T_offset_ptr* const b_y_plus,
                                          const T_offset_ptr* const b_z_minus,
                                          const T_offset_ptr* const b_z_plus,
                                          const int max_block,const T omega,
                                          const int iterations){
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
    __shared__ float Dinv_buffered[64];

    if(threadIdx.x < 64)
        Dinv_buffered[threadIdx.x] = Dinv[threadIdx.x];
    if(threadIdx.x < PREFETCH)
        if(threadIdx.x < max_block){
            block_index[threadIdx.x] = offsets[threadIdx.x];
            block_minus_x_index[threadIdx.x]=b_x_minus[threadIdx.x];
            block_plus_x_index[threadIdx.x]=b_x_plus[threadIdx.x];
            block_minus_y_index[threadIdx.x]=b_y_minus[threadIdx.x];
            block_plus_y_index[threadIdx.x]=b_y_plus[threadIdx.x];
            block_minus_z_index[threadIdx.x]=b_z_minus[threadIdx.x];
            block_plus_z_index[threadIdx.x]=b_z_plus[threadIdx.x];
        }
    for(int itr=0;itr<iterations;++itr){
        __syncthreads();
        for(int i = 0;i < PREFETCH;i += span){
            if (i + block < max_block){
                T* r_base = reinterpret_cast<T*>((unsigned long)r + (unsigned long)block_index[i + block]); 
                unsigned mask_value = reinterpret_cast<unsigned*>((unsigned long)masks + (unsigned long)block_index[i + block])[entry];
                if(mask_value & SPGrid_Solver_Cell_Type_Active){
                    T center_value = reinterpret_cast<T*>((unsigned long)u + (unsigned long)block_index[i + block])[entry];
                    T b_value = reinterpret_cast<T*>((unsigned long)rhs + (unsigned long)block_index[i + block])[entry];

                    T& x_minus_value = (x==0)
                        ? reinterpret_cast<T*>((unsigned long)u + (unsigned long)block_minus_x_index[i + block])[entry+(xsize-1)*xstride]
                        : reinterpret_cast<T*>((unsigned long)u + (unsigned long)block_index[i + block])[entry-xstride];
                    T& x_plus_value = (x==xsize-1)
                        ? reinterpret_cast<T*>((unsigned long)u + (unsigned long)block_plus_x_index[i + block])[entry-(xsize-1)*xstride]
                        : reinterpret_cast<T*>((unsigned long)u + (unsigned long)block_index[i + block])[entry+xstride];
                                
                    T& y_minus_value = (y==0)
                        ? reinterpret_cast<T*>((unsigned long)u + (unsigned long)block_minus_y_index[i + block])[entry+(ysize-1)*ystride]
                        : reinterpret_cast<T*>((unsigned long)u + (unsigned long)block_index[i + block])[entry-ystride];
                    T& y_plus_value = (y==ysize-1)
                        ? reinterpret_cast<T*>((unsigned long)u + (unsigned long)block_plus_y_index[i + block])[entry-(ysize-1)*ystride]
                        : reinterpret_cast<T*>((unsigned long)u + (unsigned long)block_index[i + block])[entry+ystride];
                    
                    T& z_minus_value = (z==0)
                        ? reinterpret_cast<T*>((unsigned long)u + (unsigned long)block_minus_z_index[i + block])[entry+(zsize-1)*zstride]
                        : reinterpret_cast<T*>((unsigned long)u + (unsigned long)block_index[i + block])[entry-zstride];
                    T& z_plus_value = (z==zsize-1)
                        ? reinterpret_cast<T*>((unsigned long)u + (unsigned long)block_plus_z_index[i + block])[entry-(zsize-1)*zstride]
                        : reinterpret_cast<T*>((unsigned long)u + (unsigned long)block_index[i + block])[entry+zstride];
                    
                    T result=0;
                    T x_minus = (x_minus_value - center_value);
                    T x_plus = (x_plus_value - center_value);
                    T y_minus = (y_minus_value - center_value);
                    T y_plus = (y_plus_value - center_value);
                    T z_minus = (z_minus_value - center_value);
                    T z_plus = (z_plus_value - center_value);

                    if (mask_value & (SPGrid_Solver_Face_Minus_X_Active))
                        result -= x_minus;
                    if (mask_value & (SPGrid_Solver_Face_Plus_X_Active))
                        result -= x_plus; 

                    if (mask_value & (SPGrid_Solver_Face_Minus_Y_Active))
                        result -= y_minus;
                    if (mask_value & (SPGrid_Solver_Face_Plus_Y_Active))
                        result -= y_plus;
 
                    if (mask_value & (SPGrid_Solver_Face_Minus_Z_Active))
                        result -= z_minus; 
                    if (mask_value & (SPGrid_Solver_Face_Plus_Z_Active))
                        result -= z_plus; 
                                   
                    r_base[entry] = b_value - result;
                }
            }
        }
        __syncthreads();
        for(int i = 0;i < PREFETCH;i += span){
            if (i + block < max_block){
                unsigned int cell_flag = reinterpret_cast<unsigned int*>((unsigned long)masks + block_index[i+block])[entry];
                if(cell_flag & SPGrid_Solver_Cell_Type_Active){
                    cell_flag = (cell_flag >> (BitLength<SPGrid_Solver_Face_Minus_X_Active>::value-1))&0x3f;
                    T* u_ptr = &(reinterpret_cast<T*>((unsigned long)u + block_index[i+block])[entry]);
                    const T* r_ptr = &(reinterpret_cast<const T*>((unsigned long)r + block_index[i+block])[entry]);
                    *(u_ptr) = *(u_ptr) + *(r_ptr)*omega*Dinv_buffered[cell_flag];
                }
            }
        }
    }
}
//#####################################################################
// Function Interior_Smoothing
//#####################################################################
// T_MASK corresponds to the mask for the data (not the mask channel)
template <class T,int log2_struct,class T_offset_ptr> 
void Smoother_Helper_CUDA<T,log2_struct,3,T_offset_ptr>::Interior_Smoothing(const unsigned int index_start,
                                                                            const unsigned int index_end,
                                                                            cudaStream_t& cuda_stream)
{
    int number_of_blocks=index_end-index_start+1;
    int number_of_cuda_blocks = (number_of_blocks%PREFETCH)?(number_of_blocks/PREFETCH+1):(number_of_blocks/PREFETCH);
    if(number_of_cuda_blocks == 0) return;
    if(number_of_cuda_blocks>1) {std::cout<<"The fast smoother kernel only supports single cuda block due to the fact that it requires global synchronization.";exit(1);}

    Interior_Smoother_Kernel_3D<T,log2_struct,block_xsize,block_ysize,block_zsize,T_offset_ptr>
        <<<number_of_cuda_blocks,THREADBLOCK,0,cuda_stream>>>
        (mask,r,u,rhs,
         b+index_start,
         b_x_minus+index_start,b_x_plus+index_start,
         b_y_minus+index_start,b_y_plus+index_start,
         b_z_minus+index_start,b_z_plus+index_start,
         number_of_blocks,omega,iterations);

    // cudaDeviceSynchronize();
    // cudaError err = cudaGetLastError();
    // if(err!=cudaSuccess) {std::cout<<"Something went wrong in residual kernel! Msg: "<< cudaGetErrorString(err)<<std::endl;abort();}
}
//#####################################################################
// Function Boundary_Smoothing
//#####################################################################
// T_MASK corresponds to the mask for the data (not the mask channel)
template <class T,int log2_struct,class T_offset_ptr> 
void Smoother_Helper_CUDA<T,log2_struct,3,T_offset_ptr>::Boundary_Smoothing(const unsigned int index_start,
                                                                            const unsigned int index_end,
                                                                            cudaStream_t& cuda_stream)
{
    int number_of_blocks=index_end-index_start+1;
    int number_of_cuda_blocks = (number_of_blocks%PREFETCH)?(number_of_blocks/PREFETCH+1):(number_of_blocks/PREFETCH);
    if(number_of_cuda_blocks == 0) return;
    if(number_of_cuda_blocks>1) {std::cout<<"The fast smoother kernel only supports single cuda block due to the fact that it requires global synchronization.";exit(1);}

    Boundary_Smoother_Kernel_3D<T,log2_struct,block_xsize,block_ysize,block_zsize,T_offset_ptr>
        <<<number_of_cuda_blocks,THREADBLOCK,0,cuda_stream>>>
        (mask,r,u,rhs,
         b+index_start,
         b_x_minus+index_start,b_x_plus+index_start,
         b_y_minus+index_start,b_y_plus+index_start,
         b_z_minus+index_start,b_z_plus+index_start,
         number_of_blocks,omega,iterations);
    
    // cudaDeviceSynchronize();
    // cudaError err = cudaGetLastError();
    // if(err!=cudaSuccess) {std::cout<<"Something went wrong in residual kernel! Msg: "<< cudaGetErrorString(err)<<std::endl;abort();}
}
//#####################################################################
// Function Bottom_Smoothing
//#####################################################################
// T_MASK corresponds to the mask for the data (not the mask channel)
template <class T,int log2_struct,class T_offset_ptr> 
void Smoother_Helper_CUDA<T,log2_struct,3,T_offset_ptr>::Bottom_Smoothing(const unsigned int index_start,
                                                                          const unsigned int index_end,
                                                                          cudaStream_t& cuda_stream)
{
    int number_of_blocks=index_end-index_start+1;
    int number_of_cuda_blocks = (number_of_blocks%PREFETCH)?(number_of_blocks/PREFETCH+1):(number_of_blocks/PREFETCH);
    if(number_of_cuda_blocks == 0) return;
    if(number_of_cuda_blocks>1) {std::cout<<"The fast smoother kernel only supports single cuda block due to the fact that it requires global synchronization.";exit(1);}

    Bottom_Smoother_Kernel_3D<T,log2_struct,block_xsize,block_ysize,block_zsize,T_offset_ptr>
        <<<number_of_cuda_blocks,THREADBLOCK,0,cuda_stream>>>
        (mask,r,u,rhs,
         b+index_start,
         b_x_minus+index_start,b_x_plus+index_start,
         b_y_minus+index_start,b_y_plus+index_start,
         b_z_minus+index_start,b_z_plus+index_start,
         number_of_blocks,omega,iterations);

    // cudaDeviceSynchronize();
    // cudaError err = cudaGetLastError();
    // if(err!=cudaSuccess) {std::cout<<"Something went wrong in residual kernel! Msg: "<< cudaGetErrorString(err)<<std::endl;abort();}
}
//#####################################################################
template class Smoother_Helper_CUDA<float,4,3,unsigned int>;
