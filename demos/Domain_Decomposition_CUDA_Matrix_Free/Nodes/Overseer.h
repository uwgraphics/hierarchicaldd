//#####################################################################
//Removed For Anonymity: Copyright Authors of Submission pap173s1
//Removed For Anonymity: Copyright Authors of Submission pap173s1
//#####################################################################
// Class Overseer
//#####################################################################
// The Overseer class managers all the nodes and maintains the maps between subdomains and the accelerators
//#####################################################################
#ifndef __OVERSEER_H__
#define __OVERSEER_H__

#include <SPGrid/Core/SPGrid_Allocator.h>
#include <SPGrid/Core/SPGrid_Set.h>
#include <vector>
#include "Communicator.h"
#include "../SPGrid_Domain_Decomposition_Solver/SPGrid_Master_Array_Linearizer.h"
#include "../SPGrid_Domain_Decomposition_Solver/SPGrid_V_Cycle_Topology.h"

namespace Domain_Decomposition{

struct Slot{
    int accelerator_index;
    int local_subdomain_index;
};

struct Accelerator{
    Communicator* communicator;
    std::vector<int> subdomain_index;
};

template <typename T,typename T_STRUCT,int d,typename T_offset_ptr>
class Overseer{
    typedef SPGrid::template SPGrid_Master_Array_Linearizer<T,SPGrid::NextLogTwo<sizeof(T_STRUCT)>::value,d,T_offset_ptr> T_LINEARIZER;
    typedef SPGrid::template SPGrid_Allocator<T_STRUCT,d> SPG_Allocator;
    int n_subdomains;
    std::vector<Accelerator> accelerators; 
    std::vector<Slot> subdomain_map; 
    bool MPI_enabled;
    bool SCIF_enabled;
    bool CUDA_enabled;
    int n_mpi_nodes;
    int n_scif_nodes;
    int n_cuda_nodes;
public:
    Overseer()
        :n_subdomains(0),n_mpi_nodes(0),n_scif_nodes(0),n_cuda_nodes(0),
         MPI_enabled(false),SCIF_enabled(false),CUDA_enabled(false){}
    ~Overseer(){for(int i=0;i<accelerators.size();++i) if(accelerators[i].communicator) delete accelerators[i].communicator;}
    void Discover(bool MPI_enabled,bool SCIF_enabled,bool CUDA_enabled);
    void Assign_Subdomains(int n_subdomains);    
    void Initialize_Subdomains(const std::vector<std::vector<T_LINEARIZER>*>& linearizer_hierarchy);
    void Send_Flags(const std::vector<std::vector<T_LINEARIZER>*>& linearizer_hierarchy,unsigned T_STRUCT::* flags_field,int number_of_threads=0);
    void Set_Channels(unsigned T_STRUCT::*& flags_field,T T_STRUCT::*& result_field,T T_STRUCT::*& rhs_field,T T_STRUCT::*& tmp_field);
void Step_One(std::vector<std::vector<T_LINEARIZER>*>& linearizer_hierarchy,SPG_Allocator& allocator,unsigned T_STRUCT::* flags_field,T T_STRUCT::* ri_field,T T_STRUCT::* sr_field,int number_of_threads=0);
    void Step_Two(std::vector<std::vector<T_LINEARIZER>*>& linearizer_hierarchy,SPG_Allocator& allocator,unsigned T_STRUCT::* flags_field,T T_STRUCT::* z_field,int number_of_threads=0);
};
};

#endif
