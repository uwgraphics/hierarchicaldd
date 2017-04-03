//#####################################################################
// Copyright 2010, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __Variable_Beta_Laplace_Helper__
#define __Variable_Beta_Laplace_Helper__
#include <algorithm>
#include "Threading_Tools/PTHREAD_QUEUE.h"

#include "SPGrid/Core/SPGrid_Utilities.h"
#include "SPGrid/Core/SPGrid_Mask.h"
#include <SPGrid_Fluids/Projection/Laplace_Helper.h>

extern PTHREAD_QUEUE* pthread_queue;

namespace SPGrid{
template<class T,int log2_struct, int d> class Laplace_Helper;
template<class T,int log2_struct, int d> class Variable_Beta_Laplace_Helper;

template<class T,int log2_struct>
class Variable_Beta_Laplace_Helper<T,log2_struct,2>:public Laplace_Helper<T,log2_struct,2>
{
    typedef Laplace_Helper<T,log2_struct,2> BASE;
    typedef typename BASE::T_MASK T_MASK;
    using BASE::d;using BASE::x;using BASE::y;using BASE::mask;using BASE::b;
    using BASE::size;using BASE::scale_uniform;using BASE::scale_nonuniform;

    using BASE::prefetch_degree;
    using BASE::og_xsize;using BASE::og_ysize;
    using BASE::xmin;using BASE::ymin;
    using BASE::xmax;using BASE::ymax;

    const T* const rho;  // variable density

public:
    using BASE::Run;
    using BASE::ComputeShadowGrid;
    using BASE::Run_Parallel;

    explicit Variable_Beta_Laplace_Helper(T* const x_input,const T* const y_input,const T* const rho_input,const unsigned* const mask_input,
        const unsigned long* const b_input,const int size_input,const double scale_uniform_in,const double scale_nonuniform_in)
        :BASE(x_input,y_input,mask_input,b_input,size_input,scale_uniform_in,scale_nonuniform_in),rho(rho_input)
    {}
//#####################################################################
    void Run_Index_Range(const int index_start,const int index_end);
//#####################################################################
};

template<class T,int log2_struct>
class Variable_Beta_Laplace_Helper<T,log2_struct,3>:public Laplace_Helper<T,log2_struct,3>
{
    typedef Laplace_Helper<T,log2_struct,3> BASE;
    typedef typename BASE::T_MASK T_MASK;
    using BASE::d;using BASE::x;using BASE::y;using BASE::mask;using BASE::b;
    using BASE::size;using BASE::scale_uniform;using BASE::scale_nonuniform;

    using BASE::prefetch_degree;
    using BASE::og_xsize;using BASE::og_ysize;using BASE::og_zsize;
    using BASE::xmin;using BASE::ymin;using BASE::zmin;
    using BASE::xmax;using BASE::ymax;using BASE::zmax;

    const T* const rho;  // variable density

public:
    using BASE::Run;
    using BASE::ComputeShadowGrid;
    using BASE::Run_Parallel;

    explicit Variable_Beta_Laplace_Helper(T* const x_input,const T* const y_input,const T* const rho_input,const unsigned* const mask_input,
        const unsigned long* const b_input,const int size_input,const double scale_uniform_in,const double scale_nonuniform_in)
        :BASE(x_input,y_input,mask_input,b_input,size_input,scale_uniform_in,scale_nonuniform_in),rho(rho_input)
    {}
//#####################################################################
    void Run_Index_Range(const int index_start,const int index_end);
//#####################################################################
};
}
#endif
