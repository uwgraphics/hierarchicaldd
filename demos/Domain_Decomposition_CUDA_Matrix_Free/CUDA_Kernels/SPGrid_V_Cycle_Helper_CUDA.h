//#####################################################################
//Removed For Anonymity: Copyright Authors of Submission pap173s1
//Removed For Anonymity: Copyright Authors of Submission pap173s1
//#####################################################################
// Subroutine SPGrid_Computations::V_Cycle
//#####################################################################
#ifndef __SPGrid_V_Cycle_Helper_CUDA_h__
#define __SPGrid_V_Cycle_Helper_CUDA_h__
#include <SPGrid_Fluids/Solvers/Domain_Decomposition/SPGRID_MULTIGRID_FLAGS.h>
#include <vector>
#include <cstring>
#include <iostream>
#include "Clear_Helper_CUDA.h"
#include "Correction_Helper_CUDA.h"
#include "Interface_Block_Helper_CUDA.h"
#include "Laplace_Helper_CUDA.h"
#include "Minus_Laplace_Helper_CUDA.h"
#include "Residual_Helper_CUDA.h"
#include "Restriction_Helper_CUDA.h"
#include "Prolongation_Helper_CUDA.h"
#include "Smoother_Helper_CUDA.h"
#include "Norm_Helper_CUDA.h"
using namespace SPGrid;
namespace SPGrid_Computations{
template<typename T_STRUCT,typename T_DATA,int d,typename T_offset_ptr>
class V_Cycle_Helper_CUDA
{
    enum{THREADBLOCK_SMOOTHER=256,MORE_SMOOTHING=5};
    typedef T_DATA T;
    typedef SPGrid_CUDA_Array_Linearizer<T,NextLogTwo<sizeof(T_STRUCT)>::value,d,T_offset_ptr> T_Linearizer;
    typedef SPGrid_Mask<NextLogTwo<sizeof(T_STRUCT)>::value, NextLogTwo<sizeof(T)>::value,d> T_MASK;
    enum{elements_per_block=T_MASK::elements_per_block};
public:
    static inline void collect_interface_blocks(T_Linearizer& linearizer,cudaStream_t& cuda_stream){
        Interface_Collect_Helper_CUDA<T_DATA,NextLogTwo<sizeof(T_STRUCT)>::value,d,T_offset_ptr>
            interface_collect_helper_cuda(reinterpret_cast<T*>(linearizer.interface_data),
                                          reinterpret_cast<T*>(linearizer.data),
                                          linearizer.b_interface,
                                          linearizer.number_of_interface_blocks);
        interface_collect_helper_cuda.Run(cuda_stream);
    }
    static inline void distribute_interface_blocks_data(unsigned T_STRUCT::* flags_field,T T_STRUCT::* data_field,T_Linearizer& linearizer,
                                                        cudaStream_t& cuda_stream){
        Interface_Distribute_Helper_CUDA<T_DATA,NextLogTwo<sizeof(T_STRUCT)>::value,d,T_offset_ptr>
            interface_distribute_helper_cuda(reinterpret_cast<T*>((unsigned long)linearizer.interface_data+elements_per_block*(unsigned long)OffsetOfMember(data_field)),
                                             reinterpret_cast<T*>((unsigned long)linearizer.data+elements_per_block*(unsigned long)OffsetOfMember(data_field)),
                                             reinterpret_cast<unsigned*>((unsigned long)linearizer.interface_data+elements_per_block*(unsigned long)OffsetOfMember(flags_field)),
                                             reinterpret_cast<unsigned*>((unsigned long)linearizer.data+elements_per_block*(unsigned long)OffsetOfMember(flags_field)),
                                             linearizer.b_interface,
                                             linearizer.number_of_interface_blocks);
        interface_distribute_helper_cuda.Run(cuda_stream);
    }
    static inline void clear_u(T_DATA T_STRUCT::* u_field,T_Linearizer& linearizer,cudaStream_t& cuda_stream){
        Clear_Helper_CUDA<T_DATA,NextLogTwo<sizeof(T_STRUCT)>::value,d,T_offset_ptr> 
            clear_helper_cuda(reinterpret_cast<T*>((unsigned long)linearizer.data+elements_per_block*(unsigned long)OffsetOfMember(u_field)),
                             linearizer.number_of_blocks);
        clear_helper_cuda.Run(cuda_stream);
    }
    static inline void clear_interface(unsigned T_STRUCT::* flags_field,T_DATA T_STRUCT::* u_field,T_Linearizer& linearizer,cudaStream_t& cuda_stream){
        Masked_Clear_Helper_CUDA<T_DATA,NextLogTwo<sizeof(T_STRUCT)>::value,d,T_offset_ptr> 
            masked_clear_helper_cuda(SPGrid_Solver_Cell_Type_Interface,
                                     reinterpret_cast<unsigned*>((unsigned long)linearizer.data+elements_per_block*(unsigned long)OffsetOfMember(flags_field)),
                                     reinterpret_cast<T*>((unsigned long)linearizer.data+elements_per_block*(unsigned long)OffsetOfMember(u_field)),
                                     linearizer.b_interface,linearizer.number_of_interface_blocks);
        masked_clear_helper_cuda.Run(cuda_stream);
    }
    static inline void interior_interface_accumulative_minus_laplace(T_DATA T_STRUCT::* u_field,T_DATA T_STRUCT::* Lu_field,unsigned T_STRUCT::* flags_field,T_Linearizer& linearizer,cudaStream_t& cuda_stream){
        Minus_Laplace_Helper_CUDA<T_DATA,NextLogTwo<sizeof(T_STRUCT)>::value,d,T_offset_ptr,true> 
            Minus_laplace_helper_cuda(reinterpret_cast<T*>((unsigned long)linearizer.data+elements_per_block*(unsigned long)OffsetOfMember(Lu_field)),
                                      reinterpret_cast<T*>((unsigned long)linearizer.data+elements_per_block*(unsigned long)OffsetOfMember(u_field)),
                                      reinterpret_cast<unsigned*>((unsigned long)linearizer.data+elements_per_block*(unsigned long)OffsetOfMember(flags_field)),
                                      linearizer.b_boundary,//should I use boundary blocks here?
                                      linearizer.b_boundary_x_minus,//should I use boundary blocks here?
                                      linearizer.b_boundary_x_plus,//should I use boundary blocks here?
                                      linearizer.b_boundary_y_minus,//should I use boundary blocks here?
                                      linearizer.b_boundary_y_plus,//should I use boundary blocks here?
                                      linearizer.b_boundary_z_minus,//should I use boundary blocks here?
                                      linearizer.b_boundary_z_plus,//should I use boundary blocks here?
                                      linearizer.number_of_boundary_blocks,SPGrid_Solver_Cell_Type_Active);
        Minus_laplace_helper_cuda.Run(cuda_stream);
    }
    static inline void interface_interior_laplace(T_DATA T_STRUCT::* u_field,T_DATA T_STRUCT::* Lu_field,unsigned T_STRUCT::* flags_field,T_Linearizer& linearizer,cudaStream_t& cuda_stream){
        Laplace_Helper_CUDA<T_DATA,NextLogTwo<sizeof(T_STRUCT)>::value,d,T_offset_ptr> 
            laplace_helper_cuda(reinterpret_cast<T*>((unsigned long)linearizer.data+elements_per_block*(unsigned long)OffsetOfMember(Lu_field)),
                                reinterpret_cast<T*>((unsigned long)linearizer.data+elements_per_block*(unsigned long)OffsetOfMember(u_field)),
                                reinterpret_cast<unsigned*>((unsigned long)linearizer.data+elements_per_block*(unsigned long)OffsetOfMember(flags_field)),
                                linearizer.b_interface,
                                linearizer.b_interface_x_minus,
                                linearizer.b_interface_x_plus,
                                linearizer.b_interface_y_minus,
                                linearizer.b_interface_y_plus,
                                linearizer.b_interface_z_minus,
                                linearizer.b_interface_z_plus,
                                linearizer.number_of_interface_blocks,SPGrid_Solver_Cell_Type_Interface);
        laplace_helper_cuda.Run(cuda_stream);
    }
    static inline void compute_residual(T_DATA T_STRUCT::* u_field,T_DATA T_STRUCT::* r_field,
                                        T_DATA T_STRUCT::* b_field,unsigned T_STRUCT::* flags_field,
                                        T_Linearizer& linearizer,cudaStream_t& cuda_stream){
        Residual_Helper_CUDA<T_DATA,NextLogTwo<sizeof(T_STRUCT)>::value,d,T_offset_ptr> 
            residual_helper_cuda(reinterpret_cast<T*>((unsigned long)linearizer.data+elements_per_block*(unsigned long)OffsetOfMember(r_field)),
                                 reinterpret_cast<T*>((unsigned long)linearizer.data+elements_per_block*(unsigned long)OffsetOfMember(u_field)),
                                 reinterpret_cast<T*>((unsigned long)linearizer.data+elements_per_block*(unsigned long)OffsetOfMember(b_field)),
                                 reinterpret_cast<unsigned*>((unsigned long)linearizer.data+elements_per_block*(unsigned long)OffsetOfMember(flags_field)),
                                 linearizer.b,
                                 linearizer.b_x_minus,
                                 linearizer.b_x_plus,
                                 linearizer.b_y_minus,
                                 linearizer.b_y_plus,
                                 linearizer.b_z_minus,
                                 linearizer.b_z_plus,
                                 linearizer.number_of_blocks);
        residual_helper_cuda.Run(cuda_stream);
    }
    static inline void interior_smoothing(T_DATA T_STRUCT::* u_field,T_DATA T_STRUCT::* r_field,
                                          T_DATA T_STRUCT::* b_field,unsigned T_STRUCT::* flags_field,
                                          T_Linearizer& linearizer,cudaStream_t& cuda_stream,int n_iterations){
        if(linearizer.number_of_blocks<=THREADBLOCK_SMOOTHER){
            //use fast smoother
            Smoother_Helper_CUDA<T_DATA,NextLogTwo<sizeof(T_STRUCT)>::value,d,T_offset_ptr> 
                 smoother_helper_cuda(reinterpret_cast<T*>((unsigned long)linearizer.data+elements_per_block*(unsigned long)OffsetOfMember(r_field)),
                                      reinterpret_cast<T*>((unsigned long)linearizer.data+elements_per_block*(unsigned long)OffsetOfMember(u_field)),
                                      reinterpret_cast<T*>((unsigned long)linearizer.data+elements_per_block*(unsigned long)OffsetOfMember(b_field)),
                                      reinterpret_cast<unsigned*>((unsigned long)linearizer.data+elements_per_block*(unsigned long)OffsetOfMember(flags_field)),
                                      linearizer.b,
                                      linearizer.b_x_minus,
                                      linearizer.b_x_plus,
                                      linearizer.b_y_minus,
                                      linearizer.b_y_plus,
                                      linearizer.b_z_minus,
                                      linearizer.b_z_plus,
                                      linearizer.number_of_blocks,T(2./3.),n_iterations*MORE_SMOOTHING);
            smoother_helper_cuda.Run_Interior(cuda_stream);
        }else{
            for(int i=0;i<n_iterations;++i){
                Residual_Helper_CUDA<T_DATA,NextLogTwo<sizeof(T_STRUCT)>::value,d,T_offset_ptr> 
                    residual_helper_cuda(reinterpret_cast<T*>((unsigned long)linearizer.data+elements_per_block*(unsigned long)OffsetOfMember(r_field)),
                                         reinterpret_cast<T*>((unsigned long)linearizer.data+elements_per_block*(unsigned long)OffsetOfMember(u_field)),
                                         reinterpret_cast<T*>((unsigned long)linearizer.data+elements_per_block*(unsigned long)OffsetOfMember(b_field)),
                                         reinterpret_cast<unsigned*>((unsigned long)linearizer.data+elements_per_block*(unsigned long)OffsetOfMember(flags_field)),
                                         linearizer.b,
                                         linearizer.b_x_minus,
                                         linearizer.b_x_plus,
                                         linearizer.b_y_minus,
                                         linearizer.b_y_plus,
                                         linearizer.b_z_minus,
                                         linearizer.b_z_plus,
                                         linearizer.number_of_blocks);
                residual_helper_cuda.Run(cuda_stream);
                Correction_Helper_CUDA<T_DATA,NextLogTwo<sizeof(T_STRUCT)>::value,d,T_offset_ptr> 
                    correction_helper_cuda(reinterpret_cast<T*>((unsigned long)linearizer.data+elements_per_block*(unsigned long)OffsetOfMember(u_field)),
                                           reinterpret_cast<T*>((unsigned long)linearizer.data+elements_per_block*(unsigned long)OffsetOfMember(r_field)),
                                           reinterpret_cast<unsigned*>((unsigned long)linearizer.data+elements_per_block*(unsigned long)OffsetOfMember(flags_field)),
                                           linearizer.number_of_blocks);
                correction_helper_cuda.Run_Interior_Blocks(linearizer.b,linearizer.number_of_blocks,cuda_stream);}}
    }
    static inline void boundary_smoothing(T_DATA T_STRUCT::* u_field,T_DATA T_STRUCT::* r_field,
                                          T_DATA T_STRUCT::* b_field,unsigned T_STRUCT::* flags_field,
                                          T_Linearizer& linearizer,cudaStream_t& cuda_stream,int n_iterations){
        if(linearizer.number_of_boundary_blocks<=THREADBLOCK_SMOOTHER){
            //use fast smoother
            Smoother_Helper_CUDA<T_DATA,NextLogTwo<sizeof(T_STRUCT)>::value,d,T_offset_ptr> 
                 smoother_helper_cuda(reinterpret_cast<T*>((unsigned long)linearizer.data+elements_per_block*(unsigned long)OffsetOfMember(r_field)),
                                      reinterpret_cast<T*>((unsigned long)linearizer.data+elements_per_block*(unsigned long)OffsetOfMember(u_field)),
                                      reinterpret_cast<T*>((unsigned long)linearizer.data+elements_per_block*(unsigned long)OffsetOfMember(b_field)),
                                      reinterpret_cast<unsigned*>((unsigned long)linearizer.data+elements_per_block*(unsigned long)OffsetOfMember(flags_field)),
                                      linearizer.b_boundary,
                                      linearizer.b_boundary_x_minus,
                                      linearizer.b_boundary_x_plus,
                                      linearizer.b_boundary_y_minus,
                                      linearizer.b_boundary_y_plus,
                                      linearizer.b_boundary_z_minus,
                                      linearizer.b_boundary_z_plus,
                                      linearizer.number_of_boundary_blocks,T(2./3.),n_iterations*MORE_SMOOTHING);
            smoother_helper_cuda.Run_Boundary(cuda_stream);
        }else{
            for(int i=0;i<n_iterations;++i){
                Residual_Helper_CUDA<T_DATA,NextLogTwo<sizeof(T_STRUCT)>::value,d,T_offset_ptr> 
                    residual_helper_cuda(reinterpret_cast<T*>((unsigned long)linearizer.data+elements_per_block*(unsigned long)OffsetOfMember(r_field)),
                                         reinterpret_cast<T*>((unsigned long)linearizer.data+elements_per_block*(unsigned long)OffsetOfMember(u_field)),
                                         reinterpret_cast<T*>((unsigned long)linearizer.data+elements_per_block*(unsigned long)OffsetOfMember(b_field)),
                                         reinterpret_cast<unsigned*>((unsigned long)linearizer.data+elements_per_block*(unsigned long)OffsetOfMember(flags_field)),
                                         linearizer.b_boundary,
                                         linearizer.b_boundary_x_minus,
                                         linearizer.b_boundary_x_plus,
                                         linearizer.b_boundary_y_minus,
                                         linearizer.b_boundary_y_plus,
                                         linearizer.b_boundary_z_minus,
                                         linearizer.b_boundary_z_plus,
                                         linearizer.number_of_boundary_blocks);
                residual_helper_cuda.Run(cuda_stream);
                Correction_Helper_CUDA<T_DATA,NextLogTwo<sizeof(T_STRUCT)>::value,d,T_offset_ptr> 
                    correction_helper_cuda(reinterpret_cast<T*>((unsigned long)linearizer.data+elements_per_block*(unsigned long)OffsetOfMember(u_field)),
                                           reinterpret_cast<T*>((unsigned long)linearizer.data+elements_per_block*(unsigned long)OffsetOfMember(r_field)),
                                           reinterpret_cast<unsigned*>((unsigned long)linearizer.data+elements_per_block*(unsigned long)OffsetOfMember(flags_field)),
                                           linearizer.number_of_blocks);
                correction_helper_cuda.Run_Boundary_Blocks(linearizer.b_boundary,linearizer.number_of_boundary_blocks,cuda_stream);
            }
        }
    }
    static inline void bottom_smoothing(T_DATA T_STRUCT::* u_field,T_DATA T_STRUCT::* r_field,
                                        T_DATA T_STRUCT::* b_field,unsigned T_STRUCT::* flags_field,
                                        T_Linearizer& linearizer,cudaStream_t& cuda_stream,int n_iterations){
        //std::cout<<"BOTTOM NUMBER OF BLOCKS: "<<linearizer.number_of_blocks<<std::endl;
        if(linearizer.number_of_blocks<=THREADBLOCK_SMOOTHER){
            //use fast smoother
            Smoother_Helper_CUDA<T_DATA,NextLogTwo<sizeof(T_STRUCT)>::value,d,T_offset_ptr> 
                 smoother_helper_cuda(reinterpret_cast<T*>((unsigned long)linearizer.data+elements_per_block*(unsigned long)OffsetOfMember(r_field)),
                                      reinterpret_cast<T*>((unsigned long)linearizer.data+elements_per_block*(unsigned long)OffsetOfMember(u_field)),
                                      reinterpret_cast<T*>((unsigned long)linearizer.data+elements_per_block*(unsigned long)OffsetOfMember(b_field)),
                                      reinterpret_cast<unsigned*>((unsigned long)linearizer.data+elements_per_block*(unsigned long)OffsetOfMember(flags_field)),
                                      linearizer.b,
                                      linearizer.b_x_minus,
                                      linearizer.b_x_plus,
                                      linearizer.b_y_minus,
                                      linearizer.b_y_plus,
                                      linearizer.b_z_minus,
                                      linearizer.b_z_plus,
                                      linearizer.number_of_blocks,
                                      T(2./3.),n_iterations);
            smoother_helper_cuda.Run_Bottom(cuda_stream);
        }else{
            for(int i=0;i<n_iterations;++i){
                Residual_Helper_CUDA<T_DATA,NextLogTwo<sizeof(T_STRUCT)>::value,d,T_offset_ptr> 
                    residual_helper_cuda(reinterpret_cast<T*>((unsigned long)linearizer.data+elements_per_block*(unsigned long)OffsetOfMember(r_field)),
                                         reinterpret_cast<T*>((unsigned long)linearizer.data+elements_per_block*(unsigned long)OffsetOfMember(u_field)),
                                         reinterpret_cast<T*>((unsigned long)linearizer.data+elements_per_block*(unsigned long)OffsetOfMember(b_field)),
                                         reinterpret_cast<unsigned*>((unsigned long)linearizer.data+elements_per_block*(unsigned long)OffsetOfMember(flags_field)),
                                         linearizer.b,
                                         linearizer.b_x_minus,
                                         linearizer.b_x_plus,
                                         linearizer.b_y_minus,
                                         linearizer.b_y_plus,
                                         linearizer.b_z_minus,
                                         linearizer.b_z_plus,
                                         linearizer.number_of_blocks);
                //std::cout<<"linearizer.number_of_blocks: "<<linearizer.number_of_blocks<<std::endl;
                residual_helper_cuda.Run(cuda_stream);
                Correction_Helper_CUDA<T_DATA,NextLogTwo<sizeof(T_STRUCT)>::value,d,T_offset_ptr> 
                    correction_helper_cuda(reinterpret_cast<T*>((unsigned long)linearizer.data+elements_per_block*(unsigned long)OffsetOfMember(u_field)),
                                           reinterpret_cast<T*>((unsigned long)linearizer.data+elements_per_block*(unsigned long)OffsetOfMember(r_field)),
                                           reinterpret_cast<unsigned*>((unsigned long)linearizer.data+elements_per_block*(unsigned long)OffsetOfMember(flags_field)),
                                           linearizer.number_of_blocks);
                correction_helper_cuda.Run(cuda_stream);
            }
        }
    }
    static inline void restriction(T_DATA T_STRUCT::* r_field,T_DATA T_STRUCT::* b_field,unsigned T_STRUCT::* flags_field,
                                   T_Linearizer& coarse_linearizer,T_Linearizer& fine_linearizer,cudaStream_t& cuda_stream){
        Restriction_Helper_CUDA<T_DATA,NextLogTwo<sizeof(T_STRUCT)>::value,d,T_offset_ptr> 
            restriction_helper_cuda(reinterpret_cast<T*>((unsigned long)coarse_linearizer.data+elements_per_block*(unsigned long)OffsetOfMember(b_field)),
                                    reinterpret_cast<const T*>((unsigned long)fine_linearizer.data+elements_per_block*(unsigned long)OffsetOfMember(r_field)),
                                    reinterpret_cast<const unsigned*>((unsigned long)coarse_linearizer.data+elements_per_block*(unsigned long)OffsetOfMember(flags_field)),
                                    coarse_linearizer.b,
                                    //Evil stuff here....given that the data lives on GPU and I am trying to flatten a vector of std_array here...
                                    coarse_linearizer.restriction_fine_blocks,
                                    coarse_linearizer.number_of_blocks);
        //std::cout<<"coarse_linearizer.number_of_blocks: "<<coarse_linearizer.number_of_blocks<<std::endl;
        restriction_helper_cuda.Run(cuda_stream);
    }
    static inline void prolongation(T_DATA T_STRUCT::* u_field,unsigned T_STRUCT::* flags_field,
                                    T_Linearizer& fine_linearizer,T_Linearizer& coarse_linearizer,cudaStream_t& cuda_stream){
        Prolongation_Helper_CUDA<T_DATA,NextLogTwo<sizeof(T_STRUCT)>::value,d,T_offset_ptr> 
            prolongation_helper_cuda(reinterpret_cast<T*>((unsigned long)fine_linearizer.data+elements_per_block*(unsigned long)OffsetOfMember(u_field)),
                                     reinterpret_cast<T*>((unsigned long)coarse_linearizer.data+elements_per_block*(unsigned long)OffsetOfMember(u_field)),
                                     reinterpret_cast<unsigned*>((unsigned long)fine_linearizer.data+elements_per_block*(unsigned long)OffsetOfMember(flags_field)),
                                     fine_linearizer.prolongation_fine_blocks,
                                     fine_linearizer.prolongation_coarse_blocks,
                                     fine_linearizer.n_prolongation_meta_blocks);
        prolongation_helper_cuda.Run(cuda_stream);
    }
    static void v_cycle(T_DATA T_STRUCT::* u_field,
                        T_DATA T_STRUCT::* r_field,
                        T_DATA T_STRUCT::* b_field,
                        unsigned T_STRUCT::* flags_field,
                        std::vector<T_Linearizer>& linearized_topology,cudaStream_t& cuda_stream,
                        int interior_itr = 1,int boundary_itr = 3){
        if(linearized_topology[0].number_of_blocks==0)return;
        const int levels = linearized_topology.size();
        for(int level = 0; level < levels - 1;++level){
            T_Linearizer& fine_linearizer   = linearized_topology[level  ];
            T_Linearizer& coarse_linearizer = linearized_topology[level+1];
            boundary_smoothing(u_field,r_field,b_field,flags_field,fine_linearizer,cuda_stream,boundary_itr); 
            interior_smoothing(u_field,r_field,b_field,flags_field,fine_linearizer,cuda_stream,interior_itr); 
            boundary_smoothing(u_field,r_field,b_field,flags_field,fine_linearizer,cuda_stream,boundary_itr);
            compute_residual(u_field,r_field,b_field,flags_field,fine_linearizer,cuda_stream);
            clear_u(u_field,coarse_linearizer,cuda_stream);
            // We need to clear b_field too, the restriction operator we wrote does accumulation....
            clear_u(b_field,coarse_linearizer,cuda_stream);
            restriction(r_field,b_field,flags_field,coarse_linearizer,fine_linearizer,cuda_stream);
        }
        bottom_smoothing(u_field,r_field,b_field,flags_field,linearized_topology[levels-1],cuda_stream,400);
        for(int level = levels-2; level >= 0;--level){
            T_Linearizer& fine_linearizer   = linearized_topology[level  ];
            T_Linearizer& coarse_linearizer = linearized_topology[level+1];
            prolongation(u_field,flags_field,fine_linearizer,coarse_linearizer,cuda_stream);
            boundary_smoothing(u_field,r_field,b_field,flags_field,fine_linearizer,cuda_stream,boundary_itr); 
            interior_smoothing(u_field,r_field,b_field,flags_field,fine_linearizer,cuda_stream,interior_itr); 
            boundary_smoothing(u_field,r_field,b_field,flags_field,fine_linearizer,cuda_stream,boundary_itr);
        }
    }

};
}
#endif
