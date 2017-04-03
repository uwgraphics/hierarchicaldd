//#####################################################################
//Removed For Anonymity: Copyright Authors of Submission pap173s1
//Removed For Anonymity: Copyright Authors of Submission pap173s1
//#####################################################################
// Class KNL_OVERSEER
//#####################################################################
#include "KNL_OVERSEER.h"
#include "../SPGrid_Domain_Decomposition_Solver/SPGRID_DOMAIN_DECOMPOSITION_DATA.h"
#include "../SPGrid_Domain_Decomposition_Solver/SPGrid_Linearized_Data_Copy_Helper.h"
#include "../SPGrid_Domain_Decomposition_Solver/SPGrid_Linearized_Data_Copy_Hashtable_Helper.h"
#include "../SPGrid_Domain_Decomposition_Solver/SPGrid_Master_Array_Linearizer.h"
#include "Kernels/SPGrid_V_Cycle_Helper_Phi.h"
#include "SPGrid_KNL_Array_Linearizer.h"
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>

using namespace PhysBAM;
using namespace SPGrid;
using namespace Domain_Decomposition;

//#####################################################################
// Function Initialize_Subdomain
//#####################################################################
template<typename T,typename T_STRUCT,int d,typename T_offset_ptr>
void KNL_OVERSEER<T,T_STRUCT,d,T_offset_ptr>::Initialize_Subdomain(std::vector<T_LINEARIZER>& linearizer_hierarchy)
{
    knl_linearized_data.clear();
    int nlevels=linearizer_hierarchy.size();
    knl_linearized_data.resize(nlevels);
    for(int level=0;level<nlevels;++level){
        knl_linearized_data[level].number_of_blocks=linearizer_hierarchy[level].number_of_blocks;
        knl_linearized_data[level].number_of_boundary_blocks=linearizer_hierarchy[level].number_of_boundary_blocks;
        knl_linearized_data[level].number_of_interface_blocks=linearizer_hierarchy[level].number_of_interface_blocks;
        knl_linearized_data[level].Allocate();
        #pragma omp parallel for
        for(int i=0;i<knl_linearized_data[level].number_of_blocks;++i){
            knl_linearized_data[level].b[i]=linearizer_hierarchy[level].b[i];
            knl_linearized_data[level].b_x_minus[i]=linearizer_hierarchy[level].b_x_minus[i];
            knl_linearized_data[level].b_x_plus[i]=linearizer_hierarchy[level].b_x_plus[i];
            knl_linearized_data[level].b_y_minus[i]=linearizer_hierarchy[level].b_y_minus[i];
            knl_linearized_data[level].b_y_plus[i]=linearizer_hierarchy[level].b_y_plus[i];
            knl_linearized_data[level].b_z_minus[i]=linearizer_hierarchy[level].b_z_minus[i];
            knl_linearized_data[level].b_z_plus[i]=linearizer_hierarchy[level].b_z_plus[i];}
        
        #pragma omp parallel for
        for(int i=0;i<knl_linearized_data[level].number_of_boundary_blocks;++i){
            knl_linearized_data[level].b_boundary[i]=linearizer_hierarchy[level].b_boundary[i];
            knl_linearized_data[level].b_boundary_x_minus[i]=linearizer_hierarchy[level].b_boundary_x_minus[i];
            knl_linearized_data[level].b_boundary_x_plus[i]=linearizer_hierarchy[level].b_boundary_x_plus[i];
            knl_linearized_data[level].b_boundary_y_minus[i]=linearizer_hierarchy[level].b_boundary_y_minus[i];
            knl_linearized_data[level].b_boundary_y_plus[i]=linearizer_hierarchy[level].b_boundary_y_plus[i];
            knl_linearized_data[level].b_boundary_z_minus[i]=linearizer_hierarchy[level].b_boundary_z_minus[i];
            knl_linearized_data[level].b_boundary_z_plus[i]=linearizer_hierarchy[level].b_boundary_z_plus[i];}

        #pragma omp parallel for
        for(int i=0;i<knl_linearized_data[level].number_of_interface_blocks;++i){
            knl_linearized_data[level].b_interface[i]=linearizer_hierarchy[level].b_interface[i];
            knl_linearized_data[level].b_interface_x_minus[i]=linearizer_hierarchy[level].b_interface_x_minus[i];
            knl_linearized_data[level].b_interface_x_plus[i]=linearizer_hierarchy[level].b_interface_x_plus[i];
            knl_linearized_data[level].b_interface_y_minus[i]=linearizer_hierarchy[level].b_interface_y_minus[i];
            knl_linearized_data[level].b_interface_y_plus[i]=linearizer_hierarchy[level].b_interface_y_plus[i];
            knl_linearized_data[level].b_interface_z_minus[i]=linearizer_hierarchy[level].b_interface_z_minus[i];
            knl_linearized_data[level].b_interface_z_plus[i]=linearizer_hierarchy[level].b_interface_z_plus[i];}
        
        #pragma omp parallel for
        for(T_offset_ptr i=0;i<knl_linearized_data[level].data_buffer_size;++i){
            knl_linearized_data[level].data[i]=linearizer_hierarchy[level].data[i];}

        if(level!=nlevels-1){
            knl_linearized_data[level].number_of_prolongation_blocks=linearizer_hierarchy[level].prolongation_fine_blocks.size();
            knl_linearized_data[level].Allocate_Prolongation_Blocks();
            #pragma omp parallel for
            for(int i=0;i<knl_linearized_data[level].number_of_prolongation_blocks*8;++i){
                knl_linearized_data[level].prolongation_fine_blocks[i]=linearizer_hierarchy[level].prolongation_fine_blocks[i/8](i%8);
                knl_linearized_data[level].prolongation_coarse_blocks[i]=linearizer_hierarchy[level].prolongation_coarse_blocks[i/8](i%8);}}
        if(level!=0){
            knl_linearized_data[level].Allocate_Restriction_Blocks();
            #pragma omp parallel for
            for(int i=0;i<knl_linearized_data[level].number_of_blocks*27;++i){
                knl_linearized_data[level].restriction_fine_blocks[i]=linearizer_hierarchy[level].restriction_fine_blocks[i/27](i%27);}}
    }
}
//#####################################################################
// Function Copy_Interface_From
//#####################################################################
template<typename T,typename T_STRUCT,int d,typename T_offset_ptr>
void KNL_OVERSEER<T,T_STRUCT,d,T_offset_ptr>::Copy_Interface_From(std::vector<T_LINEARIZER>& linearizer_hierarchy)
{
    #pragma omp parallel for
    for(T_offset_ptr i=0;i<knl_linearized_data[0].interface_data_buffer_size;++i){
        knl_linearized_data[0].interface_data[i]=linearizer_hierarchy[0].interface_data[i];}
}
//#####################################################################
// Function Copy_Interface_To
//#####################################################################
template<typename T,typename T_STRUCT,int d,typename T_offset_ptr>
void KNL_OVERSEER<T,T_STRUCT,d,T_offset_ptr>::Copy_Interface_To(std::vector<T_LINEARIZER>& linearizer_hierarchy)
{
    #pragma omp parallel for
    for(T_offset_ptr i=0;i<knl_linearized_data[0].interface_data_buffer_size;++i){
        linearizer_hierarchy[0].interface_data[i]=knl_linearized_data[0].interface_data[i];}
}
//#####################################################################
// Function Copy_Data_To
//#####################################################################
template<typename T,typename T_STRUCT,int d,typename T_offset_ptr>
void KNL_OVERSEER<T,T_STRUCT,d,T_offset_ptr>::Copy_Data_To(std::vector<T_LINEARIZER>& linearizer_hierarchy)
{
    #pragma omp parallel for
    for(T_offset_ptr i=0;i<knl_linearized_data[0].data_buffer_size;++i){
        linearizer_hierarchy[0].data[i]=knl_linearized_data[0].data[i];}
}
//#####################################################################
// Function Step_One
//#####################################################################
template<typename T,typename T_STRUCT,int d,typename T_offset_ptr>
void KNL_OVERSEER<T,T_STRUCT,d,T_offset_ptr>::Step_One(std::vector<std::vector<T_LINEARIZER>*>& linearizer_hierarchies,SPG_Allocator& allocator,unsigned T_STRUCT::* flags_field,T T_STRUCT::* ri_field,T T_STRUCT::* sr_field,int number_of_threads)
{
    LOG::SCOPE scope("Step One");
    // Send_ri
    n_subdomains=linearizer_hierarchies.size();
    {LOG::SCOPE scope("copy to linearizer");
    for(int i=0;i<n_subdomains;++i){
        std::pair<const unsigned long*,unsigned> blocks(&(*linearizer_hierarchies[i])[0].offsets_list[0],(*linearizer_hierarchies[i])[0].offsets_list.size());
        if(blocks.second==0) continue;
        //Do we need to clear the buffer?
        //For this one, we send the interface with the interior, because we need the flags in the interface cells. And the data chennel of the interface will be cleared out at the accelarator
        if(number_of_threads)
            SPGrid_Computations::Threading_Helper<T_STRUCT,d>(allocator,blocks).Run_Parallel(
              SPGrid_Computations::Copy_To_Linearized_Data_With_Interface_Clear<T_STRUCT,T,d,T_offset_ptr>(flags_field,ri_field,(*linearizer_hierarchies[i])[0].data,
                                                                                                           (*linearizer_hierarchies[i])[0].offsets_map),number_of_threads);
        else
            SPGrid_Computations::Copy_To_Linearized_Data_With_Interface_Clear<T_STRUCT,T,d,T_offset_ptr>(allocator,blocks,
                                                                                                         flags_field,ri_field,
                                                                                                         (*linearizer_hierarchies[i])[0].data,
                                                                                                         (*linearizer_hierarchies[i])[0].offsets_map);}
    }
    {LOG::SCOPE scope("Total Subdomain Solve");
    for(int i=0;i<n_subdomains;++i){
        {LOG::SCOPE scope("Initialize Subdomain");
        Initialize_Subdomain(*linearizer_hierarchies[i]);}
        SPGrid_Computations::V_Cycle_Helper_PHI<T_STRUCT,T,d,T_offset_ptr>::clear_interface(knl_flags_field,knl_rhs_field,knl_linearized_data[0]);
        SPGrid_Computations::V_Cycle_Helper_PHI<T_STRUCT,T,d,T_offset_ptr>::clear_interface(knl_flags_field,knl_tmp_field,knl_linearized_data[0]);
        SPGrid_Computations::V_Cycle_Helper_PHI<T_STRUCT,T,d,T_offset_ptr>::clear_interface(knl_flags_field,knl_result_field,knl_linearized_data[0]);
        {LOG::SCOPE scope("Clear u");
        SPGrid_Computations::V_Cycle_Helper_PHI<T_STRUCT,T,d,T_offset_ptr>::clear_u(knl_result_field,knl_linearized_data[0]);}
        {LOG::SCOPE scope("V Cycle");
        for(int j = 0;j < n_v_cycles;++j){
            SPGrid_Computations::V_Cycle_Helper_PHI<T_STRUCT,T,d,T_offset_ptr>::v_cycle(knl_result_field,knl_tmp_field,knl_rhs_field,knl_flags_field,knl_linearized_data,
                                                                                        interior_smoothing,boundary_smoothing);}}
        SPGrid_Computations::V_Cycle_Helper_PHI<T_STRUCT,T,d,T_offset_ptr>::clear_u(knl_tmp_field,knl_linearized_data[0]);
        SPGrid_Computations::V_Cycle_Helper_PHI<T_STRUCT,T,d,T_offset_ptr>::interface_interior_laplace(knl_result_field,knl_tmp_field,knl_flags_field,knl_linearized_data[0]);
        SPGrid_Computations::V_Cycle_Helper_PHI<T_STRUCT,T,d,T_offset_ptr>::collect_interface_blocks(knl_linearized_data[0]);
        {LOG::SCOPE scope("Copy Interface To");
        Copy_Interface_To(*linearizer_hierarchies[i]);}
    }
    }
    {LOG::SCOPE scope("copy from linearizer");
    for(int i=0;i<n_subdomains;++i){
        std::pair<const unsigned long*,unsigned> blocks(&(*linearizer_hierarchies[i])[0].offsets_list_interface[0],(*linearizer_hierarchies[i])[0].offsets_list_interface.size());
        if(blocks.second==0) continue;
        if(number_of_threads)
            SPGrid_Computations::Threading_Helper<T_STRUCT,d>(allocator,blocks).Run_Parallel(SPGrid_Computations::Accumulatively_Interface_Substract_From_Linearized_Data<T_STRUCT,T,d,T_offset_ptr>(flags_field,sr_field,(*linearizer_hierarchies[i])[0].interface_data,(*linearizer_hierarchies[i])[0].offsets_map_interface),number_of_threads);
        else
            SPGrid_Computations::Accumulatively_Interface_Substract_From_Linearized_Data<T_STRUCT,T,d,T_offset_ptr>(allocator,blocks,
                                                                                                                    flags_field,sr_field,
                                                                                                                    (*linearizer_hierarchies[i])[0].interface_data,
                                                                                                                    (*linearizer_hierarchies[i])[0].offsets_map_interface);}
    }
}
//#####################################################################
// Function Step_Two
//#####################################################################
template<typename T,typename T_STRUCT,int d,typename T_offset_ptr>
void KNL_OVERSEER<T,T_STRUCT,d,T_offset_ptr>::Step_Two(std::vector<std::vector<T_LINEARIZER>*>& linearizer_hierarchies,SPG_Allocator& allocator,unsigned T_STRUCT::* flags_field,T T_STRUCT::* z_field,int number_of_threads)
{
    LOG::SCOPE scope("Step Two");
    n_subdomains=linearizer_hierarchies.size();
    {LOG::SCOPE scope("copy to linearizer");
    for(int i=0;i<n_subdomains;++i){
        std::pair<const unsigned long*,unsigned> blocks(&(*linearizer_hierarchies[i])[0].offsets_list_interface[0],(*linearizer_hierarchies[i])[0].offsets_list_interface.size());
        if(number_of_threads)
            SPGrid_Computations::Threading_Helper<T_STRUCT,d>(allocator,blocks).Run_Parallel(SPGrid_Computations::Copy_Interface_To_Linearized_Data<T_STRUCT,T,d,T_offset_ptr>(flags_field,z_field,(*linearizer_hierarchies[i])[0].interface_data,(*linearizer_hierarchies[i])[0].offsets_map_interface),number_of_threads);
        else
            SPGrid_Computations::Copy_Interface_To_Linearized_Data<T_STRUCT,T,d,T_offset_ptr>(allocator,blocks,
                                                                                              flags_field,z_field,
                                                                                              (*linearizer_hierarchies[i])[0].interface_data,
                                                                                              (*linearizer_hierarchies[i])[0].offsets_map_interface);}
    }
    for(int i=0;i<n_subdomains;++i){
        {LOG::SCOPE scope("Initialize Subdomain");
        Initialize_Subdomain(*linearizer_hierarchies[i]);}
        //SPGrid_Computations::V_Cycle_Helper_PHI<T_STRUCT,T,d,T_offset_ptr>::clear_interface(knl_flags_field,knl_rhs_field,knl_linearized_data[0]);
        //SPGrid_Computations::V_Cycle_Helper_PHI<T_STRUCT,T,d,T_offset_ptr>::clear_interface(knl_flags_field,knl_tmp_field,knl_linearized_data[0]);
        //SPGrid_Computations::V_Cycle_Helper_PHI<T_STRUCT,T,d,T_offset_ptr>::clear_interface(knl_flags_field,knl_result_field,knl_linearized_data[0]);
        SPGrid_Computations::V_Cycle_Helper_PHI<T_STRUCT,T,d,T_offset_ptr>::clear_u(knl_result_field,knl_linearized_data[0]);
        {LOG::SCOPE scope("Copy Interface From");
        Copy_Interface_From(*linearizer_hierarchies[i]);}
        SPGrid_Computations::V_Cycle_Helper_PHI<T_STRUCT,T,d,T_offset_ptr>::distribute_interface_blocks_data(knl_flags_field,knl_result_field,knl_linearized_data[0]);
        
        SPGrid_Computations::V_Cycle_Helper_PHI<T_STRUCT,T,d,T_offset_ptr>::interior_interface_accumulative_minus_laplace(knl_result_field,knl_rhs_field,knl_flags_field,knl_linearized_data[0]); 
        SPGrid_Computations::V_Cycle_Helper_PHI<T_STRUCT,T,d,T_offset_ptr>::clear_interface(knl_flags_field,knl_rhs_field,knl_linearized_data[0]);
        SPGrid_Computations::V_Cycle_Helper_PHI<T_STRUCT,T,d,T_offset_ptr>::clear_interface(knl_flags_field,knl_tmp_field,knl_linearized_data[0]);
        SPGrid_Computations::V_Cycle_Helper_PHI<T_STRUCT,T,d,T_offset_ptr>::clear_interface(knl_flags_field,knl_result_field,knl_linearized_data[0]);
        //MAY NOT NEEDED!
        SPGrid_Computations::V_Cycle_Helper_PHI<T_STRUCT,T,d,T_offset_ptr>::clear_u(knl_result_field,knl_linearized_data[0]);
        {LOG::SCOPE scope("V_Cycle");
        for(int j = 0;j < n_v_cycles;++j)
            SPGrid_Computations::V_Cycle_Helper_PHI<T_STRUCT,T,d,T_offset_ptr>::v_cycle(knl_result_field,knl_tmp_field,knl_rhs_field,knl_flags_field,knl_linearized_data,
                                                                                        interior_smoothing,boundary_smoothing);}

        {LOG::SCOPE scope("Copy Data To");
        Copy_Data_To(*linearizer_hierarchies[i]);}

        {LOG::SCOPE scope("Copy Data From Linearizer");
        std::pair<const unsigned long*,unsigned> blocks(&(*linearizer_hierarchies[i])[0].offsets_list[0],(*linearizer_hierarchies[i])[0].offsets_list.size());
        if(number_of_threads)
            SPGrid_Computations::Threading_Helper<T_STRUCT,d>(allocator,blocks).Run_Parallel(SPGrid_Computations::Copy_From_Linearized_Data<T_STRUCT,T,d,T_offset_ptr>(flags_field,z_field,(*linearizer_hierarchies[i])[0].data,(*linearizer_hierarchies[i])[0].offsets_map),number_of_threads);
        else
            SPGrid_Computations::Copy_From_Linearized_Data<T_STRUCT,T,d,T_offset_ptr>(allocator,blocks,
                                                                                      flags_field,z_field,
                                                                                      (*linearizer_hierarchies[i])[0].data,
                                                                                      (*linearizer_hierarchies[i])[0].offsets_map);}}
}
//#####################################################################
template class KNL_OVERSEER<float,SPGRID_DOMAIN_DECOMPOSITION_DATA<float>,3,unsigned int>;
template class KNL_OVERSEER<float,SPGRID_DOMAIN_DECOMPOSITION_DATA<float>,3,unsigned long>;

