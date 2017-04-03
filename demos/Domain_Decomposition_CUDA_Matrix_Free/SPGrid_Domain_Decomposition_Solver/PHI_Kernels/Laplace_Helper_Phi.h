//#####################################################################
//Removed For Anonymity: Copyright Authors of Submission pap173s1
// Distributed under the FreeBSD license (see license.txt)
//#####################################################################
#ifndef __Laplace_Helper_Phi__
#define __Laplace_Helper_Phi__
#include <algorithm>
#include <SPGrid/Core/SPGrid_Utilities.h>
#include <SPGrid/Core/SPGrid_Mask.h>
#include <SPGrid_Fluids/Solvers/Domain_Decomposition/SPGRID_MULTIGRID_FLAGS.h>
#include <immintrin.h>
#include <iostream>

namespace SPGrid{
template<class T,int log2_struct, int d,class T_offset_ptr> class Laplace_Helper_PHI;

#define SIMD_LENGTH 16
template<class T,int log2_struct,class T_offset_ptr>
class Laplace_Helper_PHI<T,log2_struct,3,T_offset_ptr>
{
    enum{d=3};
    typedef SPGrid_Mask<log2_struct, NextLogTwo<sizeof(T)>::value,d> T_MASK;
    T* const Lu;         // output stream
    const T* const u;  // first input stream
    const unsigned* const mask;
    const T_offset_ptr* const b;   // block offset stream
    const int size;     // number of blocks to process
   
    const __m512i SPGrid_Cell_Type_Active_v;
    const __m512i SPGrid_Face_Minus_X_Active_v;
    const __m512i SPGrid_Face_Plus_X_Active_v;
    
    const __m512i SPGrid_Face_Minus_Y_Active_v;
    const __m512i SPGrid_Face_Plus_Y_Active_v;
    
    const __m512i SPGrid_Face_Minus_Z_Active_v;
    const __m512i SPGrid_Face_Plus_Z_Active_v;
    
    const T_offset_ptr* const b_x_plus;   // block offset stream
    const T_offset_ptr* const b_x_minus;   // block offset stream
    const T_offset_ptr* const b_y_plus;   // block offset stream
    const T_offset_ptr* const b_y_minus;   // block offset stream
    const T_offset_ptr* const b_z_plus;   // block offset stream
    const T_offset_ptr* const b_z_minus;   // block offset stream
    enum {
        block_xsize = 1u << T_MASK::block_xbits,
        block_ysize = 1u << T_MASK::block_ybits,
        block_zsize = 1u << T_MASK::block_zbits,
        elements_per_block=block_xsize*block_ysize*block_zsize,
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
        SIMD_width=16,
        zspan_simd=(SIMD_width<block_zsize)?SIMD_width:block_zsize,
        yspan_simd=(SIMD_width<block_zsize*block_ysize)?SIMD_width/block_zsize:block_ysize,
        xspan_simd=SIMD_width/block_zsize/block_ysize
    };
    enum{
        Z_PLUS_MASK_HI=
        (((0%block_zsize+1)>=zspan_simd)?(1 << 0):0) |
        (((1%block_zsize+1)>=zspan_simd)?(1 << 1):0) |
        (((2%block_zsize+1)>=zspan_simd)?(1 << 2):0) |
        (((3%block_zsize+1)>=zspan_simd)?(1 << 3):0) |
        (((4%block_zsize+1)>=zspan_simd)?(1 << 4):0) |
        (((5%block_zsize+1)>=zspan_simd)?(1 << 5):0) |
        (((6%block_zsize+1)>=zspan_simd)?(1 << 6):0) |
        (((7%block_zsize+1)>=zspan_simd)?(1 << 7):0) |
        (((8%block_zsize+1)>=zspan_simd)?(1 << 8):0) |
        (((9%block_zsize+1)>=zspan_simd)?(1 << 9):0) |
        (((10%block_zsize+1)>=zspan_simd)?(1 << 10):0) |
        (((11%block_zsize+1)>=zspan_simd)?(1 << 11):0) |
        (((12%block_zsize+1)>=zspan_simd)?(1 << 12):0) |
        (((13%block_zsize+1)>=zspan_simd)?(1 << 13):0) |
        (((14%block_zsize+1)>=zspan_simd)?(1 << 14):0) |
        (((15%block_zsize+1)>=zspan_simd)?(1 << 15):0),
        Z_PLUS_MASK_LO=~Z_PLUS_MASK_HI,
        
        Z_MINUS_MASK_LO=
        (((0 % block_zsize - 1) < 0)?(1 << 0):0) |
        (((1 % block_zsize - 1) < 0)?(1 << 1):0) |
        (((2 % block_zsize - 1) < 0)?(1 << 2):0) |
        (((3 % block_zsize - 1) < 0)?(1 << 3):0) |
        (((4 % block_zsize - 1) < 0)?(1 << 4):0) |
        (((5 % block_zsize - 1) < 0)?(1 << 5):0) |
        (((6 % block_zsize - 1) < 0)?(1 << 6):0) |
        (((7 % block_zsize - 1) < 0)?(1 << 7):0) |
        (((8 % block_zsize - 1) < 0)?(1 << 8):0) |
        (((9 % block_zsize - 1) < 0)?(1 << 9):0) |
        (((10 % block_zsize - 1) < 0)?(1 << 10):0) |
        (((11 % block_zsize - 1) < 0)?(1 << 11):0) |
        (((12 % block_zsize - 1) < 0)?(1 << 12):0) |
        (((13 % block_zsize - 1) < 0)?(1 << 13):0) |
        (((14 % block_zsize - 1) < 0)?(1 << 14):0) |
        (((15 % block_zsize - 1) < 0)?(1 << 15):0),
        Z_MINUS_MASK_HI = ~Z_MINUS_MASK_LO,
       
        Y_PLUS_MASK_HI=
        (((0/block_zsize%block_ysize+1)>=yspan_simd)?(1 << 0):0) |
        (((1/block_zsize%block_ysize+1)>=yspan_simd)?(1 << 1):0) |
        (((2/block_zsize%block_ysize+1)>=yspan_simd)?(1 << 2):0) |
        (((3/block_zsize%block_ysize+1)>=yspan_simd)?(1 << 3):0) |
        (((4/block_zsize%block_ysize+1)>=yspan_simd)?(1 << 4):0) |
        (((5/block_zsize%block_ysize+1)>=yspan_simd)?(1 << 5):0) |
        (((6/block_zsize%block_ysize+1)>=yspan_simd)?(1 << 6):0) |
        (((7/block_zsize%block_ysize+1)>=yspan_simd)?(1 << 7):0) |
        (((8/block_zsize%block_ysize+1)>=yspan_simd)?(1 << 8):0) |
        (((9/block_zsize%block_ysize+1)>=yspan_simd)?(1 << 9):0) |
        (((10/block_zsize%block_ysize+1)>=yspan_simd)?(1 << 10):0) |
        (((11/block_zsize%block_ysize+1)>=yspan_simd)?(1 << 11):0) |
        (((12/block_zsize%block_ysize+1)>=yspan_simd)?(1 << 12):0) |
        (((13/block_zsize%block_ysize+1)>=yspan_simd)?(1 << 13):0) |
        (((14/block_zsize%block_ysize+1)>=yspan_simd)?(1 << 14):0) |
        (((15/block_zsize%block_ysize+1)>=yspan_simd)?(1 << 15):0),
        Y_PLUS_MASK_LO=~Y_PLUS_MASK_HI,

        Y_MINUS_MASK_LO=
        (((0/block_zsize%block_ysize-1)<0)?(1 << 0):0) |
        (((1/block_zsize%block_ysize-1)<0)?(1 << 1):0) |
        (((2/block_zsize%block_ysize-1)<0)?(1 << 2):0) |
        (((3/block_zsize%block_ysize-1)<0)?(1 << 3):0) |
        (((4/block_zsize%block_ysize-1)<0)?(1 << 4):0) |
        (((5/block_zsize%block_ysize-1)<0)?(1 << 5):0) |
        (((6/block_zsize%block_ysize-1)<0)?(1 << 6):0) |
        (((7/block_zsize%block_ysize-1)<0)?(1 << 7):0) |
        (((8/block_zsize%block_ysize-1)<0)?(1 << 8):0) |
        (((9/block_zsize%block_ysize-1)<0)?(1 << 9):0) |
        (((10/block_zsize%block_ysize-1)<0)?(1 << 10):0) |
        (((11/block_zsize%block_ysize-1)<0)?(1 << 11):0) |
        (((12/block_zsize%block_ysize-1)<0)?(1 << 12):0) |
        (((13/block_zsize%block_ysize-1)<0)?(1 << 13):0) |
        (((14/block_zsize%block_ysize-1)<0)?(1 << 14):0) |
        (((15/block_zsize%block_ysize-1)<0)?(1 << 15):0),
        Y_MINUS_MASK_HI=~Y_MINUS_MASK_LO,
        
        X_PLUS_MASK_HI=
        ((0/block_zsize/block_ysize%block_xsize+1)>=xspan_simd)?(1<<0):0 |
        ((1/block_zsize/block_ysize%block_xsize+1)>=xspan_simd)?(1<<1):0 |
        ((2/block_zsize/block_ysize%block_xsize+1)>=xspan_simd)?(1<<2):0 |
        ((3/block_zsize/block_ysize%block_xsize+1)>=xspan_simd)?(1<<3):0 |
        ((4/block_zsize/block_ysize%block_xsize+1)>=xspan_simd)?(1<<4):0 |
        ((5/block_zsize/block_ysize%block_xsize+1)>=xspan_simd)?(1<<5):0 |
        ((6/block_zsize/block_ysize%block_xsize+1)>=xspan_simd)?(1<<6):0 |
        ((7/block_zsize/block_ysize%block_xsize+1)>=xspan_simd)?(1<<7):0 |
        ((8/block_zsize/block_ysize%block_xsize+1)>=xspan_simd)?(1<<8):0 |
        ((9/block_zsize/block_ysize%block_xsize+1)>=xspan_simd)?(1<<9):0 |
        ((10/block_zsize/block_ysize%block_xsize+1)>=xspan_simd)?(1<<10):0 |
        ((11/block_zsize/block_ysize%block_xsize+1)>=xspan_simd)?(1<<11):0 |
        ((12/block_zsize/block_ysize%block_xsize+1)>=xspan_simd)?(1<<12):0 |
        ((13/block_zsize/block_ysize%block_xsize+1)>=xspan_simd)?(1<<13):0 |
        ((14/block_zsize/block_ysize%block_xsize+1)>=xspan_simd)?(1<<14):0 |
        ((15/block_zsize/block_ysize%block_xsize+1)>=xspan_simd)?(1<<15):0,
        X_PLUS_MASK_LO=~X_PLUS_MASK_HI,

        X_MINUS_MASK_LO=
        (((0/block_zsize/block_ysize%block_xsize-1)<0)?(1<<0):0) |
        (((1/block_zsize/block_ysize%block_xsize-1)<0)?(1<<1):0) |
        (((2/block_zsize/block_ysize%block_xsize-1)<0)?(1<<2):0) |
        (((3/block_zsize/block_ysize%block_xsize-1)<0)?(1<<3):0) |
        (((4/block_zsize/block_ysize%block_xsize-1)<0)?(1<<4):0) |
        (((5/block_zsize/block_ysize%block_xsize-1)<0)?(1<<5):0) |
        (((6/block_zsize/block_ysize%block_xsize-1)<0)?(1<<6):0) |
        (((7/block_zsize/block_ysize%block_xsize-1)<0)?(1<<7):0) |
        (((8/block_zsize/block_ysize%block_xsize-1)<0)?(1<<8):0) |
        (((9/block_zsize/block_ysize%block_xsize-1)<0)?(1<<9):0) |
        (((10/block_zsize/block_ysize%block_xsize-1)<0)?(1<<10):0) |
        (((11/block_zsize/block_ysize%block_xsize-1)<0)?(1<<11):0) |
        (((12/block_zsize/block_ysize%block_xsize-1)<0)?(1<<12):0) |
        (((13/block_zsize/block_ysize%block_xsize-1)<0)?(1<<13):0) |
        (((14/block_zsize/block_ysize%block_xsize-1)<0)?(1<<14):0) |
        (((15/block_zsize/block_ysize%block_xsize-1)<0)?(1<<15):0),
        X_MINUS_MASK_HI=~X_MINUS_MASK_LO
    };

    static_assert(block_ysize*block_zsize>=SIMD_width,"Z-aligned face is not a multiple of SIMD width");
    static_assert(block_zsize<=SIMD_width,"Block Z direction exceeds SIMD width");

public:
    explicit Laplace_Helper_PHI(T* const Lu_input,const T* const u_input,const unsigned* const mask_input,const T_offset_ptr* const b_input,
                                const T_offset_ptr* const b_x_plus_input,
                                const T_offset_ptr* const b_x_minus_input,
                                const T_offset_ptr* const b_y_plus_input,
                                const T_offset_ptr* const b_y_minus_input,
                                const T_offset_ptr* const b_z_plus_input,
                                const T_offset_ptr* const b_z_minus_input,
                                const int size_input,unsigned active_flag=SPGrid_Solver_Cell_Type_Active)
        :Lu(Lu_input),u(u_input),mask(mask_input),b(b_input),
        b_x_plus(b_x_plus_input),
        b_x_minus(b_x_minus_input),
        b_y_plus(b_y_plus_input),
        b_y_minus(b_y_minus_input),
        b_z_plus(b_z_plus_input),
        b_z_minus(b_z_minus_input),
        size(size_input),
        SPGrid_Cell_Type_Active_v(_mm512_set1_epi32(active_flag)),
        SPGrid_Face_Minus_X_Active_v(_mm512_set1_epi32(SPGrid_Solver_Face_Minus_X_Active)),
        SPGrid_Face_Plus_X_Active_v(_mm512_set1_epi32(SPGrid_Solver_Face_Plus_X_Active)),
            
        SPGrid_Face_Minus_Y_Active_v(_mm512_set1_epi32(SPGrid_Solver_Face_Minus_Y_Active)),
        SPGrid_Face_Plus_Y_Active_v(_mm512_set1_epi32(SPGrid_Solver_Face_Plus_Y_Active)),
        
        SPGrid_Face_Minus_Z_Active_v(_mm512_set1_epi32(SPGrid_Solver_Face_Minus_Z_Active)),
        SPGrid_Face_Plus_Z_Active_v(_mm512_set1_epi32(SPGrid_Solver_Face_Plus_Z_Active))
    {}
    
    
    void Run()
    { 
        Run_Index_Range(0,size-1);            
    } 

//#####################################################################
    void Run_Index_Range(const int index_start,const int index_end);
//#####################################################################

};
};
#endif
