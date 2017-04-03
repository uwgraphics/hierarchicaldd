//#####################################################################
//Removed For Anonymity: Copyright Authors of Submission pap173s1
//Removed For Anonymity: Copyright Authors of Submission pap173s1
//#####################################################################
// Class Subdomain_Solver_CUDA
//#####################################################################
#ifndef __Subdomain_Solver_CUDA_H__
#define __Subdomain_Solver_CUDA_H__

#include "../SPGrid_Domain_Decomposition_Solver/SPGRID_DOMAIN_DECOMPOSITION_DATA.h"
#include "../CUDA_Kernels/SPGrid_CUDA_Array_Linearizer.h"
#include <SPGrid/Core/SPGrid_Utilities.h>
//#include "../CUDA_Kernels/SPGrid_V_Cycle_Helper_Phi.h"

namespace Domain_Decomposition{

template<class T,typename T_STRUCT, int d,class T_offset_ptr>
class Subdomain_Solver_CUDA{
public:
    Subdomain_Solver_CUDA(){}
    ~Subdomain_Solver_CUDA(){}        
    // Members
    int levels;
    std::vector<SPGrid::SPGrid_CUDA_Array_Linearizer<T,SPGrid::NextLogTwo<sizeof(T_STRUCT)>::value,d,T_offset_ptr> > linearizer_hierarchy;
    unsigned T_STRUCT::* flags_field;
    T T_STRUCT::* result_field;
    T T_STRUCT::* rhs_field;
    T T_STRUCT::* tmp_field;
};

};
#endif
