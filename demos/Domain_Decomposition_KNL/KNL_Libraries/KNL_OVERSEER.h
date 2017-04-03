//#####################################################################
//Removed For Anonymity: Copyright Authors of Submission pap173s1
//Removed For Anonymity: Copyright Authors of Submission pap173s1
//#####################################################################
// Class KNL_Overseer
//#####################################################################
// The Overseer class managers all the nodes and maintains the maps between subdomains and the accelerators
//#####################################################################
#ifndef __KNL_OVERSEER_H__
#define __KNL_OVERSEER_H__

#include <vector>
#include "../SPGrid_Domain_Decomposition_Solver/SPGrid_Master_Array_Linearizer.h"
#include "../SPGrid_Domain_Decomposition_Solver/SPGrid_V_Cycle_Topology.h"
#include "SPGrid_KNL_Array_Linearizer.h"

namespace Domain_Decomposition{

template <typename T,typename T_STRUCT,int d,typename T_offset_ptr>
class KNL_OVERSEER{
    typedef SPGrid::template SPGrid_Master_Array_Linearizer<T,SPGrid::NextLogTwo<sizeof(T_STRUCT)>::value,d,T_offset_ptr> T_LINEARIZER;
    typedef SPGrid::template SPGrid_KNL_Array_Linearizer<T,SPGrid::NextLogTwo<sizeof(T_STRUCT)>::value,d,T_offset_ptr> T_KNL_LINEARIZER;
    typedef SPGrid::template SPGrid_Allocator<T_STRUCT,d> SPG_Allocator;
    int n_subdomains;
    std::vector<T_KNL_LINEARIZER> knl_linearized_data;
    unsigned T_STRUCT::* knl_flags_field;
    T T_STRUCT::* knl_result_field;
    T T_STRUCT::* knl_rhs_field;
    T T_STRUCT::* knl_tmp_field;
    int n_v_cycles;
    int interior_smoothing;
    int boundary_smoothing;
public:
    KNL_OVERSEER():n_subdomains(0) {};
    ~KNL_OVERSEER() {};
    void Set_Parameters(int n_v_cycles_input,int interior_smoothing_input,int boundary_smoothing_input)
    {
        n_v_cycles=n_v_cycles_input;
        interior_smoothing=interior_smoothing_input;
        boundary_smoothing=boundary_smoothing_input;
    }
    void Set_Channels(unsigned T_STRUCT::* flags_field_in,T T_STRUCT::* result_field_in,T T_STRUCT::* rhs_field_in,T T_STRUCT::* tmp_field_in)
    {
        knl_flags_field=flags_field_in;
        knl_result_field=result_field_in;
        knl_rhs_field=rhs_field_in;
        knl_tmp_field=tmp_field_in;
    }
    //Note with KNL, we do MG one subdomain at a time, so this function call is only for one subdomain initialization
    void Initialize_Subdomain(std::vector<T_LINEARIZER>& linearizer_hierarchy);
    void Copy_Interface_From(std::vector<T_LINEARIZER>& linearizer_hierarchy);
    void Copy_Interface_To(std::vector<T_LINEARIZER>& linearizer_hierarchy);
    void Copy_Data_To(std::vector<T_LINEARIZER>& linearizer_hierarchy);
    void Step_One(std::vector<std::vector<T_LINEARIZER>*>& linearizer_hierarchy,SPG_Allocator& allocator,unsigned T_STRUCT::* flags_field,T T_STRUCT::* ri_field,T T_STRUCT::* sr_field,int number_of_threads=0);
    void Step_Two(std::vector<std::vector<T_LINEARIZER>*>& linearizer_hierarchy,SPG_Allocator& allocator,unsigned T_STRUCT::* flags_field,T T_STRUCT::* z_field,int number_of_threads=0);
};
}
#endif
