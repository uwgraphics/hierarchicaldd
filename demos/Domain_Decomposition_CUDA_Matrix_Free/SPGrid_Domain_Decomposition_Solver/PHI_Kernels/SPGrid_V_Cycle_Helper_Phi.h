//#####################################################################
//Removed For Anonymity: Copyright Authors of Submission pap173s1
//Removed For Anonymity: Copyright Authors of Submission pap173s1
//#####################################################################
// Subroutine SPGrid_Computations::V_Cycle
//#####################################################################
#ifndef __SPGrid_V_Cycle_Helper_h__
#define __SPGrid_V_Cycle_Helper_h__

#include <SPGrid/Core/SPGrid_Mask.h>
#include <SPGrid/Core/SPGrid_Utilities.h>
#include <SPGrid_Fluids/Solvers/Domain_Decomposition/SPGRID_MULTIGRID_FLAGS.h>
#include "SPGrid_PHI_Array_Linearizer.h"
#include "Minus_Laplace_Helper_Phi.h"
#include "Laplace_Helper_Phi.h"
#include "Residual_Helper_Phi.h"
#include "Correction_Helper_Phi.h"
#include "Restriction_Helper_Phi.h"
#include "Prolongation_Helper_Phi.h"
#include "Clear_Helper_Phi.h"
#include <vector>
#include <cstring>

using namespace SPGrid;
namespace SPGrid_Computations{
template<typename T_STRUCT,typename T_DATA,int d,typename T_offset_ptr>
class V_Cycle_Helper_PHI
{
    typedef T_DATA T;
    typedef SPGrid_PHI_Array_Linearizer<T,NextLogTwo<sizeof(T_STRUCT)>::value,d,T_offset_ptr> T_Linearizer;
    typedef SPGrid_Mask<NextLogTwo<sizeof(T_STRUCT)>::value, NextLogTwo<sizeof(T)>::value,d> T_MASK;
    enum{elements_per_block=T_MASK::elements_per_block};
public:
    static inline void collect_interface_blocks(T_Linearizer& linearizer){
        #pragma omp parallel for
        for(int i=0;i<linearizer.number_of_interface_blocks;++i){
            memcpy(reinterpret_cast<void*>((unsigned long)linearizer.interface_data+i*(unsigned long)T_Linearizer::page_size),
                   reinterpret_cast<void*>((unsigned long)linearizer.data+linearizer.b_interface[i]),T_Linearizer::page_size);
        }
    }
    static inline T norm(T T_STRUCT::* data_field,T_Linearizer& linearizer){
        T result = 0;
        for(int i=0;i<linearizer.number_of_blocks;++i){
            T* data = reinterpret_cast<T*>((unsigned long)linearizer.data+i*(unsigned long)T_Linearizer::page_size+elements_per_block*(unsigned long)OffsetOfMember(data_field));
            for(int entry=0;entry<elements_per_block;++entry){
                if(fabs(data[entry])>result) result=fabs(data[entry]);}}
        return result;
    }
    static inline T norm(T T_STRUCT::* data_field,unsigned T_STRUCT::* flags_field,T_Linearizer& linearizer,unsigned& flag){
        T result = 0;
        for(int i=0;i<linearizer.number_of_blocks;++i){
            T* data = reinterpret_cast<T*>((unsigned long)linearizer.data+i*(unsigned long)T_Linearizer::page_size+elements_per_block*(unsigned long)OffsetOfMember(data_field));
            unsigned* flags = reinterpret_cast<unsigned*>((unsigned long)linearizer.data+i*(unsigned long)T_Linearizer::page_size+elements_per_block*(unsigned long)OffsetOfMember(flags_field));
            for(int entry=0;entry<elements_per_block;++entry){
                if(fabs(data[entry])>result) {result=fabs(data[entry]);flag=flags[entry];}}}
        return result;
    }
    static inline bool check_rank(unsigned T_STRUCT::* flags_field,T_Linearizer& linearizer,int rank){
        for(int i=0;i<linearizer.number_of_blocks;++i){
            unsigned* flags = reinterpret_cast<unsigned*>((unsigned long)linearizer.data+i*(unsigned long)T_Linearizer::page_size+elements_per_block*(unsigned long)OffsetOfMember(flags_field));
            for(int entry=0;entry<elements_per_block;++entry){
                if(((flags[entry]&SPGrid_Solver_PartitionID_Mask)!=rank)&&((flags[entry]&SPGrid_Solver_PartitionID_Mask)!=0)) {
                    std::cout<<rank<<","<<flags[entry]<<std::endl;return false;}}}
        return true;
    }
    static inline void set_boundary_one(unsigned T_STRUCT::* flags_field,T T_STRUCT::* data_field,T_Linearizer& linearizer){
        for(int i=0;i<linearizer.number_of_boundary_blocks;++i){
            T* data = reinterpret_cast<T*>((unsigned long)linearizer.data+(unsigned long)linearizer.b_boundary[i]+elements_per_block*(unsigned long)OffsetOfMember(data_field));
            unsigned* flag = reinterpret_cast<unsigned*>((unsigned long)linearizer.data+(unsigned long)linearizer.b_boundary[i]+elements_per_block*(unsigned long)OffsetOfMember(flags_field));
            for(int entry=0;entry<elements_per_block;++entry){
                if(flag[entry]&SPGrid_Solver_Cell_Type_Boundary)data[entry]=1;}}
    }
    static inline bool check_zero(unsigned T_STRUCT::* flags_field,T T_STRUCT::* data_field,T_Linearizer& linearizer,unsigned flags_to_check){
        for(int i=0;i<linearizer.number_of_blocks;++i){
            T* data = reinterpret_cast<T*>((unsigned long)linearizer.data+i*(unsigned long)T_Linearizer::page_size+elements_per_block*(unsigned long)OffsetOfMember(data_field));
            unsigned* flag = reinterpret_cast<unsigned*>((unsigned long)linearizer.data+i*(unsigned long)T_Linearizer::page_size+elements_per_block*(unsigned long)OffsetOfMember(flags_field));
            for(int entry=0;entry<elements_per_block;++entry){
                if((flag==0)&&(data[entry]!=0)) return false;
                if((flag[entry]&flags_to_check)&&(data[entry]!=0)) return false;}}
        return true;
    }
    static inline void distribute_interface_blocks_data(unsigned T_STRUCT::* flags_field,T T_STRUCT::* data_field,T_Linearizer& linearizer){
        #pragma omp parallel for
        for(int i=0;i<linearizer.number_of_interface_blocks;++i){
            unsigned* flags_interface = reinterpret_cast<unsigned*>((unsigned long)linearizer.interface_data+i*(unsigned long)T_Linearizer::page_size+elements_per_block*(unsigned long)OffsetOfMember(flags_field));
            T* data_interface = reinterpret_cast<T*>((unsigned long)linearizer.interface_data+i*(unsigned long)T_Linearizer::page_size+elements_per_block*(unsigned long)OffsetOfMember(data_field));
            unsigned* flags_global = reinterpret_cast<unsigned*>((unsigned long)linearizer.data+linearizer.b_interface[i]+elements_per_block*(unsigned long)OffsetOfMember(flags_field));
            T* data_global = reinterpret_cast<T*>((unsigned long)linearizer.data+linearizer.b_interface[i]+elements_per_block*(unsigned long)OffsetOfMember(data_field));
            for(int entry=0;entry<elements_per_block;++entry){
                //std::cout<<data_interface[entry]<<std::endl;
                if(flags_interface[entry]&SPGrid_Solver_Cell_Type_Interface){
                    if(!(flags_global[entry]&SPGrid_Solver_Cell_Type_Interface)){
                        std::cerr << "Un matched flags... Something is wrong!" << std::endl;
                        abort();}
                    data_global[entry]=data_interface[entry];
                }
            }
        }
    }
    static inline void interface_interior_laplace(T_DATA T_STRUCT::* u_field,T_DATA T_STRUCT::* Lu_field,unsigned T_STRUCT::* flags_field,T_Linearizer& linearizer){
        Laplace_Helper_PHI<T_DATA,NextLogTwo<sizeof(T_STRUCT)>::value,d,T_offset_ptr> 
            laplace_helper_phi(reinterpret_cast<T*>((unsigned long)linearizer.data+elements_per_block*(unsigned long)OffsetOfMember(Lu_field)),
                               reinterpret_cast<T*>((unsigned long)linearizer.data+elements_per_block*(unsigned long)OffsetOfMember(u_field)),
                               reinterpret_cast<unsigned*>((unsigned long)linearizer.data+elements_per_block*(unsigned long)OffsetOfMember(flags_field)),
                               linearizer.b_interface,
                               linearizer.b_interface_x_plus,
                               linearizer.b_interface_x_minus,
                               linearizer.b_interface_y_plus,
                               linearizer.b_interface_y_minus,
                               linearizer.b_interface_z_plus,
                               linearizer.b_interface_z_minus,
                               linearizer.number_of_interface_blocks,SPGrid_Solver_Cell_Type_Interface);
        laplace_helper_phi.Run();
    }
    static inline void interior_interface_accumulative_minus_laplace(T_DATA T_STRUCT::* u_field,T_DATA T_STRUCT::* Lu_field,unsigned T_STRUCT::* flags_field,T_Linearizer& linearizer){
        Minus_Laplace_Helper_PHI<T_DATA,NextLogTwo<sizeof(T_STRUCT)>::value,d,T_offset_ptr,true> 
            Minus_laplace_helper_phi(reinterpret_cast<T*>((unsigned long)linearizer.data+elements_per_block*(unsigned long)OffsetOfMember(Lu_field)),
                                     reinterpret_cast<T*>((unsigned long)linearizer.data+elements_per_block*(unsigned long)OffsetOfMember(u_field)),
                                     reinterpret_cast<unsigned*>((unsigned long)linearizer.data+elements_per_block*(unsigned long)OffsetOfMember(flags_field)),
                                     linearizer.b_boundary,//should I use boundary blocks here?
                                     linearizer.b_boundary_x_plus,//should I use boundary blocks here?
                                     linearizer.b_boundary_x_minus,//should I use boundary blocks here?
                                     linearizer.b_boundary_y_plus,//should I use boundary blocks here?
                                     linearizer.b_boundary_y_minus,//should I use boundary blocks here?
                                     linearizer.b_boundary_z_plus,//should I use boundary blocks here?
                                     linearizer.b_boundary_z_minus,//should I use boundary blocks here?
                                     linearizer.number_of_boundary_blocks,SPGrid_Solver_Cell_Type_Active);
        Minus_laplace_helper_phi.Run();
    }
    static inline void laplace(T_DATA T_STRUCT::* u_field,T_DATA T_STRUCT::* Lu_field,unsigned T_STRUCT::* flags_field,T_Linearizer& linearizer){
        Laplace_Helper_PHI<T_DATA,NextLogTwo<sizeof(T_STRUCT)>::value,d,T_offset_ptr> 
            laplace_helper_phi(reinterpret_cast<T*>((unsigned long)linearizer.data+elements_per_block*(unsigned long)OffsetOfMember(Lu_field)),
                               reinterpret_cast<T*>((unsigned long)linearizer.data+elements_per_block*(unsigned long)OffsetOfMember(u_field)),
                               reinterpret_cast<unsigned*>((unsigned long)linearizer.data+elements_per_block*(unsigned long)OffsetOfMember(flags_field)),
                               linearizer.b,
                               linearizer.b_x_plus,
                               linearizer.b_x_minus,
                               linearizer.b_y_plus,
                               linearizer.b_y_minus,
                               linearizer.b_z_plus,
                               linearizer.b_z_minus,
                               linearizer.number_of_blocks,SPGrid_Solver_Cell_Type_Active);
        laplace_helper_phi.Run();
    }
    static inline void boundary_smoothing(T_DATA T_STRUCT::* u_field,T_DATA T_STRUCT::* r_field,
                                          T_DATA T_STRUCT::* b_field,unsigned T_STRUCT::* flags_field,
                                          T_Linearizer& linearizer){
        Residual_Helper_PHI<T_DATA,NextLogTwo<sizeof(T_STRUCT)>::value,d,T_offset_ptr> 
            residual_helper_phi(reinterpret_cast<T*>((unsigned long)linearizer.data+elements_per_block*(unsigned long)OffsetOfMember(r_field)),
                                reinterpret_cast<T*>((unsigned long)linearizer.data+elements_per_block*(unsigned long)OffsetOfMember(u_field)),
                                reinterpret_cast<T*>((unsigned long)linearizer.data+elements_per_block*(unsigned long)OffsetOfMember(b_field)),
                                reinterpret_cast<unsigned*>((unsigned long)linearizer.data+elements_per_block*(unsigned long)OffsetOfMember(flags_field)),
                                linearizer.b_boundary,
                                linearizer.b_boundary_x_plus,
                                linearizer.b_boundary_x_minus,
                                linearizer.b_boundary_y_plus,
                                linearizer.b_boundary_y_minus,
                                linearizer.b_boundary_z_plus,
                                linearizer.b_boundary_z_minus,
                                linearizer.number_of_boundary_blocks);
        residual_helper_phi.Run();
        Correction_Helper_PHI<T_DATA,NextLogTwo<sizeof(T_STRUCT)>::value,d,T_offset_ptr> 
            correction_helper_phi(reinterpret_cast<T*>((unsigned long)linearizer.data+elements_per_block*(unsigned long)OffsetOfMember(u_field)),
                                  reinterpret_cast<T*>((unsigned long)linearizer.data+elements_per_block*(unsigned long)OffsetOfMember(r_field)),
                                  reinterpret_cast<unsigned*>((unsigned long)linearizer.data+elements_per_block*(unsigned long)OffsetOfMember(flags_field)),
                                  linearizer.number_of_blocks);
        correction_helper_phi.Run_Boundary_Blocks(linearizer.b_boundary,linearizer.number_of_boundary_blocks);
    }
    static inline void interior_smoothing(T_DATA T_STRUCT::* u_field,T_DATA T_STRUCT::* r_field,
                                          T_DATA T_STRUCT::* b_field,unsigned T_STRUCT::* flags_field,
                                          T_Linearizer& linearizer){
        Residual_Helper_PHI<T_DATA,NextLogTwo<sizeof(T_STRUCT)>::value,d,T_offset_ptr> 
            residual_helper_phi(reinterpret_cast<T*>((unsigned long)linearizer.data+elements_per_block*(unsigned long)OffsetOfMember(r_field)),
                                reinterpret_cast<T*>((unsigned long)linearizer.data+elements_per_block*(unsigned long)OffsetOfMember(u_field)),
                                reinterpret_cast<T*>((unsigned long)linearizer.data+elements_per_block*(unsigned long)OffsetOfMember(b_field)),
                                reinterpret_cast<unsigned*>((unsigned long)linearizer.data+elements_per_block*(unsigned long)OffsetOfMember(flags_field)),
                                linearizer.b,
                                linearizer.b_x_plus,
                                linearizer.b_x_minus,
                                linearizer.b_y_plus,
                                linearizer.b_y_minus,
                                linearizer.b_z_plus,
                                linearizer.b_z_minus,
                                linearizer.number_of_blocks);
        residual_helper_phi.Run();
        Correction_Helper_PHI<T_DATA,NextLogTwo<sizeof(T_STRUCT)>::value,d,T_offset_ptr> 
            correction_helper_phi(reinterpret_cast<T*>((unsigned long)linearizer.data+elements_per_block*(unsigned long)OffsetOfMember(u_field)),
                                  reinterpret_cast<T*>((unsigned long)linearizer.data+elements_per_block*(unsigned long)OffsetOfMember(r_field)),
                                  reinterpret_cast<unsigned*>((unsigned long)linearizer.data+elements_per_block*(unsigned long)OffsetOfMember(flags_field)),
                                  linearizer.number_of_blocks);
        correction_helper_phi.Run_Interior_Blocks(linearizer.b,linearizer.number_of_blocks);
    }
    static inline void bottom_smoothing(T_DATA T_STRUCT::* u_field,T_DATA T_STRUCT::* r_field,
                                        T_DATA T_STRUCT::* b_field,unsigned T_STRUCT::* flags_field,
                                        T_Linearizer& linearizer){
        Residual_Helper_PHI<T_DATA,NextLogTwo<sizeof(T_STRUCT)>::value,d,T_offset_ptr> 
            residual_helper_phi(reinterpret_cast<T*>((unsigned long)linearizer.data+elements_per_block*(unsigned long)OffsetOfMember(r_field)),
                                reinterpret_cast<T*>((unsigned long)linearizer.data+elements_per_block*(unsigned long)OffsetOfMember(u_field)),
                                reinterpret_cast<T*>((unsigned long)linearizer.data+elements_per_block*(unsigned long)OffsetOfMember(b_field)),
                                reinterpret_cast<unsigned*>((unsigned long)linearizer.data+elements_per_block*(unsigned long)OffsetOfMember(flags_field)),
                                linearizer.b,
                                linearizer.b_x_plus,
                                linearizer.b_x_minus,
                                linearizer.b_y_plus,
                                linearizer.b_y_minus,
                                linearizer.b_z_plus,
                                linearizer.b_z_minus,
                                linearizer.number_of_blocks);
        residual_helper_phi.Run();
        Correction_Helper_PHI<T_DATA,NextLogTwo<sizeof(T_STRUCT)>::value,d,T_offset_ptr> 
            correction_helper_phi(reinterpret_cast<T*>((unsigned long)linearizer.data+elements_per_block*(unsigned long)OffsetOfMember(u_field)),
                                  reinterpret_cast<T*>((unsigned long)linearizer.data+elements_per_block*(unsigned long)OffsetOfMember(r_field)),
                                  reinterpret_cast<unsigned*>((unsigned long)linearizer.data+elements_per_block*(unsigned long)OffsetOfMember(flags_field)),
                                  linearizer.number_of_blocks);
        correction_helper_phi.Run();
    }
    static inline void compute_residual(T_DATA T_STRUCT::* u_field,T_DATA T_STRUCT::* r_field,
                                        T_DATA T_STRUCT::* b_field,unsigned T_STRUCT::* flags_field,
                                        T_Linearizer& linearizer){
        Residual_Helper_PHI<T_DATA,NextLogTwo<sizeof(T_STRUCT)>::value,d,T_offset_ptr> 
            residual_helper_phi(reinterpret_cast<T*>((unsigned long)linearizer.data+elements_per_block*(unsigned long)OffsetOfMember(r_field)),
                                reinterpret_cast<T*>((unsigned long)linearizer.data+elements_per_block*(unsigned long)OffsetOfMember(u_field)),
                                reinterpret_cast<T*>((unsigned long)linearizer.data+elements_per_block*(unsigned long)OffsetOfMember(b_field)),
                                reinterpret_cast<unsigned*>((unsigned long)linearizer.data+elements_per_block*(unsigned long)OffsetOfMember(flags_field)),
                                linearizer.b,
                                linearizer.b_x_plus,
                                linearizer.b_x_minus,
                                linearizer.b_y_plus,
                                linearizer.b_y_minus,
                                linearizer.b_z_plus,
                                linearizer.b_z_minus,
                                linearizer.number_of_blocks);
        residual_helper_phi.Run();
    }
    static inline void clear_u(T_DATA T_STRUCT::* u_field,T_Linearizer& linearizer){
        Clear_Helper_PHI<T_DATA,NextLogTwo<sizeof(T_STRUCT)>::value,d,T_offset_ptr> 
            clear_helper_phi(reinterpret_cast<T*>((unsigned long)linearizer.data+elements_per_block*(unsigned long)OffsetOfMember(u_field)),
                             linearizer.number_of_blocks);
        clear_helper_phi.Run();
    }
    static inline void clear_interface(unsigned T_STRUCT::* flags_field,T_DATA T_STRUCT::* u_field,T_Linearizer& linearizer){
        Masked_Clear_Helper_PHI<T_DATA,NextLogTwo<sizeof(T_STRUCT)>::value,d,T_offset_ptr> 
            masked_clear_helper_phi(reinterpret_cast<T*>((unsigned long)linearizer.data+elements_per_block*(unsigned long)OffsetOfMember(u_field)),
                                    reinterpret_cast<unsigned*>((unsigned long)linearizer.data+elements_per_block*(unsigned long)OffsetOfMember(flags_field)),
                                    linearizer.b_interface,linearizer.number_of_interface_blocks,SPGrid_Solver_Cell_Type_Interface);
        masked_clear_helper_phi.Run();
    }
    static inline void restriction(T_DATA T_STRUCT::* r_field,T_DATA T_STRUCT::* b_field,unsigned T_STRUCT::* flags_field,
                                   T_Linearizer& coarse_linearizer,T_Linearizer& fine_linearizer){
        Restriction_Helper_PHI<T_DATA,NextLogTwo<sizeof(T_STRUCT)>::value,d,T_offset_ptr> 
            restriction_helper_phi(reinterpret_cast<T*>((unsigned long)coarse_linearizer.data+elements_per_block*(unsigned long)OffsetOfMember(b_field)),
                                   reinterpret_cast<const T*>((unsigned long)fine_linearizer.data+elements_per_block*(unsigned long)OffsetOfMember(r_field)),
                                   reinterpret_cast<const unsigned*>((unsigned long)coarse_linearizer.data+elements_per_block*(unsigned long)OffsetOfMember(flags_field)),
                                   coarse_linearizer.b,
                                   &(coarse_linearizer.restriction_fine_blocks[0]),
                                   coarse_linearizer.number_of_blocks);
        restriction_helper_phi.Run();
    }
    static inline void prolongation(T_DATA T_STRUCT::* u_field,unsigned T_STRUCT::* flags_field,
                                    T_Linearizer& fine_linearizer,T_Linearizer& coarse_linearizer){
        Prolongation_Helper_PHI<T_DATA,NextLogTwo<sizeof(T_STRUCT)>::value,d,T_offset_ptr> 
            prolongation_helper_phi(reinterpret_cast<T*>((unsigned long)fine_linearizer.data+elements_per_block*(unsigned long)OffsetOfMember(u_field)),
                                    reinterpret_cast<T*>((unsigned long)coarse_linearizer.data+elements_per_block*(unsigned long)OffsetOfMember(u_field)),
                                    reinterpret_cast<unsigned*>((unsigned long)fine_linearizer.data+elements_per_block*(unsigned long)OffsetOfMember(flags_field)),
                                    &(fine_linearizer.prolongation_fine_blocks[0]),
                                    &(fine_linearizer.prolongation_coarse_blocks[0]),
                                    fine_linearizer.prolongation_fine_blocks.size());
        prolongation_helper_phi.Run();
    }
    static void v_cycle(T_DATA T_STRUCT::* u_field,
                        T_DATA T_STRUCT::* r_field,
                        T_DATA T_STRUCT::* b_field,
                        unsigned T_STRUCT::* flags_field,
                        std::vector<T_Linearizer>& linearized_topology,
                        int interior_itr = 1,int boundary_itr = 3){
        if(linearized_topology[0].number_of_blocks==0)return;
        const int levels = linearized_topology.size();
        for(int level = 0; level < levels - 1;++level){
            T_Linearizer& fine_linearizer   = linearized_topology[level  ];
            T_Linearizer& coarse_linearizer = linearized_topology[level+1];
            for(int i = 0; i < boundary_itr; ++i) boundary_smoothing(u_field,r_field,b_field,flags_field,fine_linearizer); 
            for(int i = 0; i < interior_itr; ++i) interior_smoothing(u_field,r_field,b_field,flags_field,fine_linearizer); 
            for(int i = 0; i < boundary_itr; ++i) boundary_smoothing(u_field,r_field,b_field,flags_field,fine_linearizer);
            compute_residual(u_field,r_field,b_field,flags_field,fine_linearizer);
            clear_u(u_field,coarse_linearizer);
            // We need to clear b_field too, the restriction operator we wrote does accumulation....
            clear_u(b_field,coarse_linearizer);
            //std::cout<<"residual before restriction(PHI): "<<norm(r_field,fine_linearizer)<<std::endl;
            restriction(r_field,b_field,flags_field,coarse_linearizer,fine_linearizer);
            //std::cout<<"rhs after restriction(PHI): "<<norm(b_field,coarse_linearizer)<<std::endl;
        }
        compute_residual(u_field,r_field,b_field,flags_field,linearized_topology[levels-1]);
        for(int i = 0;i < 200;++i) bottom_smoothing(u_field,r_field,b_field,flags_field,linearized_topology[levels-1]);
        //DEBUG
        compute_residual(u_field,r_field,b_field,flags_field,linearized_topology[levels-1]);
        //std::cout<<"Resudual of the bottom smoothing(PHI): "<<norm(r_field,linearized_topology[levels-1])<<std::endl;
        for(int level = levels-2; level >= 0;--level){
            T_Linearizer& fine_linearizer   = linearized_topology[level  ];
            T_Linearizer& coarse_linearizer = linearized_topology[level+1];
            //std::cout<<"u before prolongation coarse(PHI): "<<norm(u_field,coarse_linearizer)<<std::endl;
            //std::cout<<"u before prolongation fine(PHI): "<<norm(u_field,fine_linearizer)<<std::endl;
            prolongation(u_field,flags_field,fine_linearizer,coarse_linearizer);
            //std::cout<<"u after prolongation(PHI): "<<norm(u_field,fine_linearizer)<<std::endl;
            for(int i = 0; i < boundary_itr; ++i) boundary_smoothing(u_field,r_field,b_field,flags_field,fine_linearizer); 
            for(int i = 0; i < interior_itr; ++i) interior_smoothing(u_field,r_field,b_field,flags_field,fine_linearizer); 
            for(int i = 0; i < boundary_itr; ++i) boundary_smoothing(u_field,r_field,b_field,flags_field,fine_linearizer);
        }
    }
    //////////////////////////////////////////////////////////////////////////////////
    // The Following Codes Are For Debuging Purpose
    //////////////////////////////////////////////////////////////////////////////////
};
//#####################################################################
}
#endif
