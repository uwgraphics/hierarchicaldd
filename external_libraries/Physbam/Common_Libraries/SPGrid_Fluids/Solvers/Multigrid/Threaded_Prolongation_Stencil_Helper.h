//#####################################################################
// Copyright 2010, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __Threaded_Prolongation_Stencil_Helper__
#define __Threaded_Prolongation_Stencil_Helper__
#include <algorithm>
#include <Threading_Tools/PTHREAD_QUEUE.h>
#include <SPGrid/Core/SPGrid_Utilities.h>
#include <SPGrid/Core/SPGrid_Mask.h>

extern PTHREAD_QUEUE* pthread_queue;

namespace SPGrid{
template<class T,int log2_struct,int d> class Threaded_Prolongation_Stencil_Helper;
//#####################################################################
// Class Threaded_Prolongation_Stencil_Helper - 2D
//#####################################################################
template<class T,int log2_struct>
class Threaded_Prolongation_Stencil_Helper<T,log2_struct,2>
{
    enum{d=2};
    typedef SPGrid_Mask<log2_struct,NextLogTwo<sizeof(T)>::value,d> T_MASK;
    T* const fine_data;
    const T* const coarse_data;
    const unsigned* const fine_mask;
    const unsigned long* const b;   // block offset stream
    const int size;     // number of blocks to process
    const unsigned active_flag_mask;

    enum{
        block_xsize = 1u << T_MASK::block_xbits,
        block_ysize = 1u << T_MASK::block_ybits,
        og_xsize = block_xsize+2,
        og_ysize = block_ysize+2,
        xmin = 1,
        ymin = 1,
        // Inclusive!!! give mins and maxs for actual block within shadow grid
        xmax = og_xsize-2,
        ymax = og_ysize-2,
        coarse_og_xsize = block_xsize/2+2,
        coarse_og_ysize = block_ysize/2+2
    };
    enum{parents_per_cell=4};
    enum{parity_x_mask=T_MASK::template LinearOffset<1,0>::value,
         parity_y_mask=T_MASK::template LinearOffset<0,1>::value};
//#####################################################################
public:
    explicit Threaded_Prolongation_Stencil_Helper(T* const fine_data_input,const T* const coarse_data_input,const unsigned* const fine_mask_input,const unsigned long* const b_input,
        const int size_input,const unsigned active_flag_mask_input)
        :fine_data(fine_data_input),coarse_data(coarse_data_input),fine_mask(fine_mask_input),b(b_input),size(size_input),active_flag_mask(active_flag_mask_input)
    {}
    void Run()
    {Run_Index_Range(0,size-1);}
    static void ComputeCoarseShadowGrid(unsigned long* offset_grid_ptr,const unsigned long coarse_packed_offset);
//#####################################################################
    void Run_Parallel(const int number_of_partitions);
    void Run_Index_Range(const int index_start,const int index_end);
//#####################################################################
};
//#####################################################################
// Class Threaded_Prolongation_Stencil_Helper - 3D
//#####################################################################
template<class T,int log2_struct>
class Threaded_Prolongation_Stencil_Helper<T,log2_struct,3>
{
    enum{d=3};
    typedef SPGrid_Mask<log2_struct,NextLogTwo<sizeof(T)>::value,d> T_MASK;
    T* const fine_data;
    const T* const coarse_data;
    const unsigned* const fine_mask;
    const unsigned long* const b;   // block offset stream
    const int size;     // number of blocks to process
    const unsigned active_flag_mask;

    enum{
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
        zmax = og_zsize-2,
        coarse_og_xsize = block_xsize/2+2,
        coarse_og_ysize = block_ysize/2+2,
        coarse_og_zsize = block_zsize/2+2
    };
    enum{parents_per_cell=8};
    enum{parity_x_mask=T_MASK::template LinearOffset<1,0,0>::value,
         parity_y_mask=T_MASK::template LinearOffset<0,1,0>::value,
         parity_z_mask=T_MASK::template LinearOffset<0,0,1>::value};
//#####################################################################
public:
    explicit Threaded_Prolongation_Stencil_Helper(T* const fine_data_input,const T* const coarse_data_input,const unsigned* const fine_mask_input,const unsigned long* const b_input,
        const int size_input,const unsigned active_flag_mask_input)
        :fine_data(fine_data_input),coarse_data(coarse_data_input),fine_mask(fine_mask_input),b(b_input),size(size_input),active_flag_mask(active_flag_mask_input)
    {}
    void Run()
    {Run_Index_Range(0,size-1);}
    static void ComputeCoarseShadowGrid(unsigned long* offset_grid_ptr,const unsigned long coarse_packed_offset);
//#####################################################################
    void Run_Parallel(const int number_of_partitions);
    void Run_Index_Range(const int index_start,const int index_end);
//#####################################################################
};
}
#endif
