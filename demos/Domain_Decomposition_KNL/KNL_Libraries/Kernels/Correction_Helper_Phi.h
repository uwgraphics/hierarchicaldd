//#####################################################################
//Removed For Anonymity: Copyright Authors of Submission pap173s1
// Distributed under the FreeBSD license (see license.txt)
//#####################################################################
#ifndef __Correction_Helper_Phi__
#define __Correction_Helper_Phi__
#include <algorithm>
#include <SPGrid/Core/SPGrid_Utilities.h>
#include <SPGrid/Core/SPGrid_Mask.h>
#include <SPGrid_Fluids/Solvers/Domain_Decomposition/SPGRID_MULTIGRID_FLAGS.h>
#include <immintrin.h>
#include <iostream>

namespace SPGrid{

template<class T,int log2_struct, int d,class T_offset_ptr> class Correction_Helper_PHI;

template<class T,int log2_struct,class T_offset_ptr>
    class Correction_Helper_PHI<T,log2_struct,3,T_offset_ptr>
{
    enum{d=3};
    typedef SPGrid_Mask<log2_struct, NextLogTwo<sizeof(T)>::value,d> T_MASK;
    T* const u;         // output stream
    const T* const r;  // first input stream
    const unsigned* const mask;
    const int size;     // number of blocks to process
    const T omega;
   
    enum {
        SIMD_width=16,
    };

public:
    explicit Correction_Helper_PHI(T* const u_input,const T* const r_input,const unsigned* const mask_input,const int size_input,T omega_input = 2./3.)
        :r(r_input),u(u_input),mask(mask_input),size(size_input),omega(omega_input)
    {}
    
    void Run()
    {Run_Index_Range(0,size-1);} 
    
    //#####################################################################
    void Run_Index_Range(const int index_start,const int index_end);
    //#####################################################################
    void Run_Boundary_Blocks(const T_offset_ptr* offsets,const int number_of_blocks);
    void Run_Interior_Blocks(const T_offset_ptr* offsets,const int number_of_blocks);
};

};
#endif
