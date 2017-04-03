//#####################################################################
// Copyright 2010, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __Convergence_Norm_Helper__
#define __Convergence_Norm_Helper__
namespace SPGrid{

template<class T,int elements_per_block>
class Convergence_Norm_Helper
{
    const T* const y;  // second stream
    //const unsigned* const flags;  // flags stream
    const unsigned long* const b;   // block offset stream
    const int size;     // number of blocks to process

public:
    explicit Convergence_Norm_Helper(const T* const y_input,/*const unsigned* const flags_input,*/ const unsigned long* const b_input, const int size_input)
        :y(y_input),/*flags(flags_input),*/b(b_input),size(size_input)
    {}

    T Run()
    {
        T result;
        Run_Index_Range(0,size-1,result);
        return result;
    }

//#####################################################################
    T Run_Parallel(const int number_of_partitions);
    void Run_Index_Range(const int index_start,const int index_end,T& result);
//#####################################################################
};
}
#endif
