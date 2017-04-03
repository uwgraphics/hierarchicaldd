//#####################################################################
// Copyright 2010, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __Blocked_Add_Helper__
#define __Blocked_Add_Helper__
namespace SPGrid{

template<class T,int elements_per_block>
class Blocked_Add_Helper
{
    T* const x;         // output stream
    const T* const y1;  // input stream
    const unsigned* const flags;  // flags stream
    const unsigned long* const b;   // block offset stream
    const int size;     // number of blocks to process
    //static const int prefetch_degree = 0;

public:
    explicit Blocked_Add_Helper(T* const x_input,const T* const y_input, const unsigned* const flags_input, const unsigned long* const b_input, const int size_input)
        :x(x_input),y1(y_input),flags(flags_input),b(b_input),size(size_input)
    {
    }

    void Run()
    {
        Run_Index_Range(0,size-1);
    }

//#####################################################################
    void Run_Parallel(const int number_of_partitions);
    void Run_Index_Range(const int index_start,const int index_end);
//#####################################################################
};
}
#endif
