//#####################################################################
// Copyright 2010, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __Ghost_Value_Propagate__
#define __Ghost_Value_Propagate__
#include <algorithm>
#include <map>
#include <Threading_Tools/PTHREAD_QUEUE.h>

#include <SPGrid/Core/SPGrid_Utilities.h>
#include <SPGrid/Core/SPGrid_Mask.h>
#include <SPGrid/Core/SPGrid_Allocator.h>
#include <SPGrid/Core/SPGrid_Set.h>


extern PTHREAD_QUEUE* pthread_queue;

namespace SPGrid{

template<class T,class T_STRUCT,int d>
class Ghost_Value_Propagate
{
    unsigned* fine_mask;
    const T* fine_data_ptr;
    const T* coarse_data_ptr;
    const unsigned long* const b;   // block offset stream
    const int size;     // number of blocks to process
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<T>::type Array_type;

public:
    explicit Ghost_Value_Propagate(unsigned* mask_input, const T* fdp, const T* cdp, const unsigned long* const b_input, const int size_input)
        :fine_mask(mask_input),fine_data_ptr(fdp),coarse_data_ptr(cdp),b(b_input),size(size_input)
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
