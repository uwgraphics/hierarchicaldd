//#####################################################################
// Copyright 2010, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __Smoother_Helper__
#define __Smoother_Helper__
#include <algorithm>
#include <Threading_Tools/PTHREAD_QUEUE.h>
#include <SPGrid/Core/SPGrid_Utilities.h>
#include <SPGrid/Core/SPGrid_Mask.h>

extern PTHREAD_QUEUE* pthread_queue;

namespace SPGrid{
template<class T,int log2_struct, int d> class Smoother_Helper;

template<class T,int log2_struct>
class Smoother_Helper<T,log2_struct,2>
{
    enum{d=2};
    typedef SPGrid_Mask<log2_struct, NextLogTwo<sizeof(T)>::value,d> T_MASK;
    typedef T T2; // type for calculation (T or double)
    const T* const x;
    const T* const rhs;
    T* const delta;
    const T* const dinv;
    const unsigned* const mask;
    const unsigned long* const b;   // block offset stream
    const int size;     // number of blocks to process
    const T2 laplace_scale_uniform;
    const T2 laplace_scale_nonuniform;
    const T omega;
    const unsigned active_flag_mask;
    
    enum {
        prefetch_degree = 0,
        block_xsize = 1u << T_MASK::block_xbits,
        block_ysize = 1u << T_MASK::block_ybits,
        og_xsize = block_xsize+2,
        og_ysize = block_ysize+2,
        xmin = 1,
        ymin = 1,
        // Inclusive!!! give mins and maxs for actual block within shadow grid
        xmax = og_xsize-2,
        ymax = og_ysize-2,
    };

public:
    explicit Smoother_Helper(const T* const x_input,const T* const rhs_input,T* const delta_input,const T* const dinv_input,const unsigned* const mask_input,
        const unsigned long* const b_input,const int size_input,const T2 laplace_scale_uniform_input,const T2 laplace_scale_nonuniform_input,const T omega_input,const unsigned active_flag_mask_input)
        :x(x_input),rhs(rhs_input),delta(delta_input),dinv(dinv_input),mask(mask_input),b(b_input),size(size_input),laplace_scale_uniform(laplace_scale_uniform_input),
         laplace_scale_nonuniform(laplace_scale_nonuniform_input),omega(omega_input),active_flag_mask(active_flag_mask_input)
    {}
    void Run()
    {Run_Index_Range(0,size-1);}
    static void ComputeShadowGrid(unsigned long* offset_grid_ptr, unsigned* mask_in, unsigned long packed_offset);
//#####################################################################
    void Run_Parallel(const int number_of_partitions);
    void Run_Index_Range(const int index_start,const int index_end);
//#####################################################################
};


template<class T,int log2_struct>
class Smoother_Helper<T,log2_struct,3>
{
    enum{d=3};
    typedef SPGrid_Mask<log2_struct, NextLogTwo<sizeof(T)>::value,d> T_MASK;
    typedef T T2; // type for calculation (T or double)
    const T* const x;
    const T* const rhs;
    T* const delta;
    const T* const dinv;
    const unsigned* const mask;
    const unsigned long* const b;   // block offset stream
    const int size;     // number of blocks to process
    const T2 laplace_scale_uniform;
    const T2 laplace_scale_nonuniform;
    const T omega;
    const unsigned active_flag_mask;
    
    enum {
        prefetch_degree = 0,
        block_xsize = 1u << T_MASK::block_xbits,
        block_ysize = 1u << T_MASK::block_ybits,
        block_zsize = 1u << T_MASK::block_zbits,
        og_xsize = block_xsize+2,
        og_ysize = block_ysize+2,
        og_zsize = block_zsize+2,
        xmin = 1,
        ymin = 1,
        zmin = 1,
        // Inclusive!!! give mins and maxs for actual block within shadow grid
        xmax = og_xsize-2,
        ymax = og_ysize-2,
        zmax = og_zsize-2
    };

public:
    explicit Smoother_Helper(const T* const x_input,const T* const rhs_input,T* const delta_input,const T* const dinv_input,const unsigned* const mask_input,
        const unsigned long* const b_input,const int size_input,const T2 laplace_scale_uniform_input,const T2 laplace_scale_nonuniform_input,const T omega_input,const unsigned active_flag_mask_input)
        :x(x_input),rhs(rhs_input),delta(delta_input),dinv(dinv_input),mask(mask_input),b(b_input),size(size_input),laplace_scale_uniform(laplace_scale_uniform_input),
         laplace_scale_nonuniform(laplace_scale_nonuniform_input),omega(omega_input),active_flag_mask(active_flag_mask_input)
    {}
    void Run()
    {Run_Index_Range(0,size-1);}
    static void ComputeShadowGrid(unsigned long* offset_grid_ptr, unsigned* mask_in, unsigned long packed_offset);
//#####################################################################
    void Run_Parallel(const int number_of_partitions);
    void Run_Index_Range(const int index_start,const int index_end);
//#####################################################################

};

};
#endif
