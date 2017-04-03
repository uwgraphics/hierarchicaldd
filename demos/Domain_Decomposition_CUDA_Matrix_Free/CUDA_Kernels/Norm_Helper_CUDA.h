//#####################################################################
//Removed For Anonymity: Copyright Authors of Submission pap173s1
// Distributed under the FreeBSD license (see license.txt)
//#####################################################################
#ifndef __NORM_HELPER_CUDA__
#define __NORM_HELPER_CUDA__
#include <algorithm>
#include <SPGrid/Core/SPGrid_Mask.h>

namespace SPGrid{
template<class T,int log2_struct, int d,class T_offset_ptr> class Norm_Helper_CUDA;

template<class T,int log2_struct,class T_offset_ptr>
class Norm_Helper_CUDA<T,log2_struct,3,T_offset_ptr>
{
    // return max||u||
    enum{d=3};
    typedef SPGrid_Mask<log2_struct, NextLogTwo<sizeof(T)>::value,d> T_MASK;
    const T* const u;         // input stream
    const unsigned int size;
    enum {
        page_size=4096,
        block_xsize = 1u << T_MASK::block_xbits,
        block_ysize = 1u << T_MASK::block_ybits,
        block_zsize = 1u << T_MASK::block_zbits
    };

public:
    explicit Norm_Helper_CUDA(const T* const u_input,const unsigned int size_input)
        :u(u_input),size(size_input)
    {}

    T Run();
};

}
#endif
