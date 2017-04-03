//#####################################################################
// Copyright 2010, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __Inner_Product_Helper__
#define __Inner_Product_Helper__
namespace SPGrid{

template<class T,int elements_per_block>
class Inner_Product_Helper
{
    const T* const y1;  // second stream
    const T* const y2;  // second stream
    const unsigned* const flags;  // flags stream
    const unsigned long* const b;   // block offset stream
    const int size;     // number of blocks to process

public:
    explicit Inner_Product_Helper(const T* const y1_input, const T* const y2_input, const unsigned* const flags_input, const unsigned long* const b_input, const int size_input)
        :y1(y1_input),y2(y2_input),flags(flags_input),b(b_input),size(size_input)
    {}

    double Run()
    {
        double result;
        Run_Index_Range(0,size-1,result);
        return result;
    }

//#####################################################################
    double Run_Parallel(const int number_of_partitions);
    void Run_Index_Range(const int index_start,const int index_end,double& result);
//#####################################################################
};
}
#endif
