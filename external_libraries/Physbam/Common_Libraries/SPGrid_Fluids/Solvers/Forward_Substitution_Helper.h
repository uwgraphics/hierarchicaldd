//#####################################################################
// Copyright 2010, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __Forward_Substitution_Helper__
#define __Forward_Substitution_Helper__
#include <algorithm>
#include <Threading_Tools/PTHREAD_QUEUE.h>
#include <SPGrid/Core/SPGrid_Utilities.h>
#include <SPGrid/Core/SPGrid_Mask.h>


extern PTHREAD_QUEUE* pthread_queue;

namespace SPGrid{
template<class T,int log2_struct, int d> class Forward_Substitution_Helper;

template<class T,int log2_struct>
class Forward_Substitution_Helper<T,log2_struct,2>
{
    enum{d=2};
    typedef SPGrid_Mask<log2_struct, NextLogTwo<sizeof(T)>::value,d> T_MASK;

    T* const x;                   // Solution of forward substitution
    const T* const Lx;
    const T* const Ly;
    const unsigned* const mask;   // Flags
    const unsigned long* const b; // Block offset stream
    const int size;               // Number of blocks to process
    
    enum {
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
    explicit Forward_Substitution_Helper(T* const x_input,const T* const Lx_input,const T* const Ly_input,
        const unsigned* const mask_input,const unsigned long* const b_input,const int size_input)
        :x(x_input),Lx(Lx_input),Ly(Ly_input),mask(mask_input),b(b_input),size(size_input)
    {}

    void Run()
    { Run_Index_Range(0,size-1); }

//#####################################################################
    void Run_Parallel(const int number_of_partitions);
    void Run_Index_Range(const int index_start,const int index_end);
//#####################################################################
};


template<class T,int log2_struct>
class Forward_Substitution_Helper<T,log2_struct,3>
{
    enum{d=3};
    typedef SPGrid_Mask<log2_struct, NextLogTwo<sizeof(T)>::value,d> T_MASK;

    T* const x;                   // Solution of forward substitution
    const T* const Lx;
    const T* const Ly;
    const T* const Lz;
    const unsigned* const mask;   // Flags
    const unsigned long* const b; // Block offset stream
    const int size;               // Number of blocks to process
    
    enum {
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
    explicit Forward_Substitution_Helper(T* const x_input,const T* const Lx_input,const T* const Ly_input,const T* const Lz_input,
        const unsigned* const mask_input,const unsigned long* const b_input,const int size_input)
        :x(x_input),Lx(Lx_input),Ly(Ly_input),Lz(Lz_input),mask(mask_input),b(b_input),size(size_input)
    {}

    void Run()
    { Run_Index_Range(0,size-1); }

//#####################################################################
    void Run_Parallel(const int number_of_partitions);
    void Run_Index_Range(const int index_start,const int index_end);
//#####################################################################

};

};
#endif
