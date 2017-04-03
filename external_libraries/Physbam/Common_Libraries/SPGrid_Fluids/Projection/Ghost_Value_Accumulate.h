//#####################################################################
// Copyright 2010, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __Ghost_Value_Accumulate__
#define __Ghost_Value_Accumulate__
#include <algorithm>
#include <map>
#include <Threading_Tools/PTHREAD_QUEUE.h>

#include <SPGrid/Core/SPGrid_Utilities.h>
#include <SPGrid/Core/SPGrid_Mask.h>
#include <SPGrid/Core/SPGrid_Allocator.h>
#include <SPGrid/Core/SPGrid_Set.h>


extern PTHREAD_QUEUE* pthread_queue;

namespace SPGrid{

template<class T,class T_STRUCT,int d> class Ghost_Value_Accumulate;

template<class T,class T_STRUCT>
class Ghost_Value_Accumulate<T,T_STRUCT,2>
{
    enum{d=2};
    unsigned* coarse_mask;
    const T* fine_data_ptr;
    const T* coarse_data_ptr;
    const unsigned long* const b;   // block offset stream
    const int size;     // number of blocks to process
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<T>::type Array_type;

public:
    explicit Ghost_Value_Accumulate(unsigned* cmi, const T* fdp, const T* cdp, const unsigned long* const b_input, const int size_input)
        :coarse_mask(cmi),fine_data_ptr(fdp),coarse_data_ptr(cdp),b(b_input),size(size_input)
        {}

    void Run()
    { Run_Index_Range(0,size-1); }

//#####################################################################
    void Run_Parallel(const int number_of_partitions);
    void Run_Index_Range(const int index_start,const int index_end);
//#####################################################################

};

template<class T,class T_STRUCT>
class Ghost_Value_Accumulate<T,T_STRUCT,3>
{
    enum{d=3};
    unsigned* coarse_mask;
    const T* fine_data_ptr;
    const T* coarse_data_ptr;
    const unsigned long* const b;   // block offset stream
    const int size;     // number of blocks to process
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<T>::type Array_type;

public:
    explicit Ghost_Value_Accumulate(unsigned* cmi, const T* fdp, const T* cdp, const unsigned long* const b_input, const int size_input)
        :coarse_mask(cmi),fine_data_ptr(fdp),coarse_data_ptr(cdp),b(b_input),size(size_input)
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
