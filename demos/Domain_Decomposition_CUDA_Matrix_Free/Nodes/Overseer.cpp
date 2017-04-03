//#####################################################################
//Removed For Anonymity: Copyright Authors of Submission pap173s1
//Removed For Anonymity: Copyright Authors of Submission pap173s1
//#####################################################################
// Class Overseer
//#####################################################################
#ifdef __NVCC__
#include <cuda_runtime_api.h>
#include <cuda.h>
#endif
#ifdef MPI_ENABLED
#include <mpi.h>
#endif
#include <iostream>
#include <stdlib.h>
#include "Client_Node_CUDA.h"
#include "Command.h"
#include "Overseer.h"
#include "../SPGrid_Domain_Decomposition_Solver/SPGRID_DOMAIN_DECOMPOSITION_DATA.h"
#include "../SPGrid_Domain_Decomposition_Solver/SPGrid_Linearized_Data_Copy_Helper.h"
#include "../SPGrid_Domain_Decomposition_Solver/SPGrid_Linearized_Data_Copy_Hashtable_Helper.h"
#include "../SPGrid_Domain_Decomposition_Solver/SPGrid_Master_Array_Linearizer.h"

#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>

//DEBUG
#include "../SPGrid_Library/SPGrid_V_Cycle_Helper.h"
using namespace PhysBAM;
using namespace SPGrid;
using namespace Domain_Decomposition;
//#####################################################################
// Function Discover
//#####################################################################
template<typename T,typename T_STRUCT,int d,typename T_offset_ptr>
void Overseer<T,T_STRUCT,d,T_offset_ptr>::Discover(bool MPI_enabled_in,bool SCIF_enabled_in,bool CUDA_enabled_in){
    MPI_enabled = MPI_enabled_in;
    SCIF_enabled = SCIF_enabled_in;
    CUDA_enabled = CUDA_enabled_in;

    for(int i=0;i<accelerators.size();++i) if(accelerators[i].communicator) delete accelerators[i].communicator;
    accelerators.clear();
    int rank;
    if(MPI_enabled){
#ifdef MPI_ENABLED
        MPI_Init(0,NULL);
        MPI_Comm_size(MPI_COMM_WORLD,&n_mpi_nodes);
        MPI_Comm_size(MPI_COMM_WORLD,&n_mpi_nodes);
        MPI_Comm_rank(MPI_COMM_WORLD,&rank);
        if(rank!=0){std::cerr<<"Master is not running on rank 0!"<<std::endl;abort();}
        for(int mpi_node=1;mpi_node<n_mpi_nodes;++mpi_node){
            Accelerator accelerator;
            accelerator.communicator = new Communicator_MPI(mpi_node/*rank*/);
            accelerators.push_back(accelerator);}
#endif
    }
    if(SCIF_enabled){
        //TODO
    }
    if(CUDA_enabled){
#ifdef __NVCC__
        cudaGetDeviceCount(&n_cuda_nodes);
        std::cout<<"Found "<<n_cuda_nodes<<" GPUs"<<std::endl;
        for(int cuda_node=0;cuda_node<n_cuda_nodes;++cuda_node){
            Accelerator accelerator;
            accelerator.communicator = new Communicator_CUDA(cuda_node);
            accelerators.push_back(accelerator);
            new Client_Node_CUDA<T,T_STRUCT,d,T_offset_ptr>(cuda_node,
                                                            dynamic_cast<Communicator_CUDA*>(accelerator.communicator)->command_queue,
                                                            dynamic_cast<Communicator_CUDA*>(accelerator.communicator)->integer_queue,
                                                            dynamic_cast<Communicator_CUDA*>(accelerator.communicator)->data_send_queue,
                                                            dynamic_cast<Communicator_CUDA*>(accelerator.communicator)->data_recv_queue,
                                                            dynamic_cast<Communicator_CUDA*>(accelerator.communicator)->data_recv_confirmation_queue);
        }
#endif
    }
}
//#####################################################################
// Function Assign_Subdomains
//#####################################################################
template<typename T,typename T_STRUCT,int d,typename T_offset_ptr>
void Overseer<T,T_STRUCT,d,T_offset_ptr>::Assign_Subdomains(int n_subdomains){
    std::cout<<"n_subdomains: "<<n_subdomains<<std::endl;
    subdomain_map.resize(n_subdomains);
    int n_accelerators=accelerators.size();
    for(int accelerator=0;accelerator<n_accelerators;++accelerator){
        accelerators[accelerator].subdomain_index.resize(0);}
    // Round Robin
    for(int subdomain=0;subdomain<n_subdomains;++subdomain){
        subdomain_map[subdomain].accelerator_index=subdomain%n_accelerators;
        subdomain_map[subdomain].local_subdomain_index=subdomain/n_accelerators;
        accelerators[subdomain%n_accelerators].subdomain_index.push_back(subdomain);}
}
//#####################################################################
// Function Initialize_Subdomains
//#####################################################################
template<typename T,typename T_STRUCT,int d,typename T_offset_ptr>
void Overseer<T,T_STRUCT,d,T_offset_ptr>::Initialize_Subdomains(const std::vector<std::vector<T_LINEARIZER>*>& linearizer_hierarchies){
    for(int accelerator=0;accelerator<accelerators.size();++accelerator){
        accelerators[accelerator].communicator->Send_Command(INIT);
        accelerators[accelerator].communicator->Send_Int(accelerators[accelerator].subdomain_index.size());
        typedef typename T_LINEARIZER::T_offset_ptr T_offset_ptr;
        for(int subdomain = 0;subdomain < accelerators[accelerator].subdomain_index.size();++subdomain){
            //Use the correct linearizer to send data!
            const std::vector<T_LINEARIZER>& linearizer_hierarchy=*linearizer_hierarchies[accelerators[accelerator].subdomain_index[subdomain]];
            const int levels = linearizer_hierarchy.size();
            accelerators[accelerator].communicator->Send_Int(levels);
            for(int level = 0;level < levels;++level){
                int n_blocks = linearizer_hierarchy[level].number_of_blocks;
                int n_boundary_blocks = linearizer_hierarchy[level].number_of_boundary_blocks;
                int n_interface_blocks = linearizer_hierarchy[level].number_of_interface_blocks;
                accelerators[accelerator].communicator->Send_Int(n_blocks);
                accelerators[accelerator].communicator->Send_Int(n_boundary_blocks);
                accelerators[accelerator].communicator->Send_Int(n_interface_blocks);
                if(n_blocks==0) continue;
                accelerators[accelerator].communicator->Send_Buffer(n_blocks*sizeof(T_offset_ptr),linearizer_hierarchy[level].b);
                accelerators[accelerator].communicator->Send_Buffer(n_blocks*sizeof(T_offset_ptr),linearizer_hierarchy[level].b_x_minus);
                accelerators[accelerator].communicator->Send_Buffer(n_blocks*sizeof(T_offset_ptr),linearizer_hierarchy[level].b_x_plus);
                accelerators[accelerator].communicator->Send_Buffer(n_blocks*sizeof(T_offset_ptr),linearizer_hierarchy[level].b_y_minus);
                accelerators[accelerator].communicator->Send_Buffer(n_blocks*sizeof(T_offset_ptr),linearizer_hierarchy[level].b_y_plus);
                accelerators[accelerator].communicator->Send_Buffer(n_blocks*sizeof(T_offset_ptr),linearizer_hierarchy[level].b_z_minus);
                accelerators[accelerator].communicator->Send_Buffer(n_blocks*sizeof(T_offset_ptr),linearizer_hierarchy[level].b_z_plus);

                accelerators[accelerator].communicator->Send_Buffer(n_boundary_blocks*sizeof(T_offset_ptr),linearizer_hierarchy[level].b_boundary);
                accelerators[accelerator].communicator->Send_Buffer(n_boundary_blocks*sizeof(T_offset_ptr),linearizer_hierarchy[level].b_boundary_x_minus);
                accelerators[accelerator].communicator->Send_Buffer(n_boundary_blocks*sizeof(T_offset_ptr),linearizer_hierarchy[level].b_boundary_x_plus);
                accelerators[accelerator].communicator->Send_Buffer(n_boundary_blocks*sizeof(T_offset_ptr),linearizer_hierarchy[level].b_boundary_y_minus);
                accelerators[accelerator].communicator->Send_Buffer(n_boundary_blocks*sizeof(T_offset_ptr),linearizer_hierarchy[level].b_boundary_y_plus);
                accelerators[accelerator].communicator->Send_Buffer(n_boundary_blocks*sizeof(T_offset_ptr),linearizer_hierarchy[level].b_boundary_z_minus);
                accelerators[accelerator].communicator->Send_Buffer(n_boundary_blocks*sizeof(T_offset_ptr),linearizer_hierarchy[level].b_boundary_z_plus);

                accelerators[accelerator].communicator->Send_Buffer(n_interface_blocks*sizeof(T_offset_ptr),linearizer_hierarchy[level].b_interface);
                accelerators[accelerator].communicator->Send_Buffer(n_interface_blocks*sizeof(T_offset_ptr),linearizer_hierarchy[level].b_interface_x_minus);
                accelerators[accelerator].communicator->Send_Buffer(n_interface_blocks*sizeof(T_offset_ptr),linearizer_hierarchy[level].b_interface_x_plus);
                accelerators[accelerator].communicator->Send_Buffer(n_interface_blocks*sizeof(T_offset_ptr),linearizer_hierarchy[level].b_interface_y_minus);
                accelerators[accelerator].communicator->Send_Buffer(n_interface_blocks*sizeof(T_offset_ptr),linearizer_hierarchy[level].b_interface_y_plus);
                accelerators[accelerator].communicator->Send_Buffer(n_interface_blocks*sizeof(T_offset_ptr),linearizer_hierarchy[level].b_interface_z_minus);
                accelerators[accelerator].communicator->Send_Buffer(n_interface_blocks*sizeof(T_offset_ptr),linearizer_hierarchy[level].b_interface_z_plus);

                int n_meta_blocks = linearizer_hierarchy[level].prolongation_fine_blocks.size();
                if(level != levels-1){
                    accelerators[accelerator].communicator->Send_Int(n_meta_blocks);
                    accelerators[accelerator].communicator->Send_Buffer(n_meta_blocks*8*sizeof(T_offset_ptr),(void*)linearizer_hierarchy[level].prolongation_fine_blocks[0].data);
                    accelerators[accelerator].communicator->Send_Buffer(n_meta_blocks*8*sizeof(T_offset_ptr),(void*)linearizer_hierarchy[level].prolongation_coarse_blocks[0].data);}
                if(level != 0){
                    //JUST A SANITY CHECK...
                    if(n_blocks>=2)
                        PHYSBAM_ASSERT((unsigned long)linearizer_hierarchy[level].restriction_fine_blocks[1].data - 
                                       (unsigned long)linearizer_hierarchy[level].restriction_fine_blocks[0].data == 27*sizeof(T_offset_ptr));
                    accelerators[accelerator].communicator->Send_Buffer(n_blocks*27*sizeof(T_offset_ptr),(void*)linearizer_hierarchy[level].restriction_fine_blocks[0].data);}
            }
        }
    }
}
//#####################################################################
// Function Send_Flags
//#####################################################################
template<typename T,typename T_STRUCT,int d,typename T_offset_ptr>
void Overseer<T,T_STRUCT,d,T_offset_ptr>::Send_Flags(const std::vector<std::vector<T_LINEARIZER>*>& linearizer_hierarchies,unsigned T_STRUCT::* flags_field,int number_of_threads)
{
    n_subdomains=linearizer_hierarchies.size();
    const int nlevels = linearizer_hierarchies[0]->size(); 
    for(int accelerator=0;accelerator<accelerators.size();++accelerator){
        accelerators[accelerator].communicator->Send_Command(SEND_FLAGS);
        for(int subdomain = 0;subdomain < accelerators[accelerator].subdomain_index.size();++subdomain){
            for(int level=0;level<nlevels;++level){
                //Use the correct linearizer to send data!
                const T_LINEARIZER& top=(*linearizer_hierarchies[accelerators[accelerator].subdomain_index[subdomain]])[level];
                accelerators[accelerator].communicator->Send_Buffer(top.data_buffer_size,top.data);}}}
}
//#####################################################################
// Function Set_Channels
//#####################################################################
template<typename T,typename T_STRUCT,int d,typename T_offset_ptr>
void Overseer<T,T_STRUCT,d,T_offset_ptr>::Set_Channels(unsigned T_STRUCT::*& flags_field,T T_STRUCT::*& result_field,T T_STRUCT::*& rhs_field,T T_STRUCT::*& tmp_field){
    for(int accelerator=0;accelerator<accelerators.size();++accelerator){
        accelerators[accelerator].communicator->Send_Command(SET_CHANNELS);
        accelerators[accelerator].communicator->Send_Buffer(sizeof(size_t),(void*)&flags_field);
        accelerators[accelerator].communicator->Send_Buffer(sizeof(size_t),(void*)&result_field);
        accelerators[accelerator].communicator->Send_Buffer(sizeof(size_t),(void*)&rhs_field);
        accelerators[accelerator].communicator->Send_Buffer(sizeof(size_t),(void*)&tmp_field);}
}
//#####################################################################
// Function Step_One
//#####################################################################
template<typename T,typename T_STRUCT,int d,typename T_offset_ptr>
void Overseer<T,T_STRUCT,d,T_offset_ptr>::Step_One(std::vector<std::vector<T_LINEARIZER>*>& linearizer_hierarchies,SPG_Allocator& allocator,unsigned T_STRUCT::* flags_field,T T_STRUCT::* ri_field,T T_STRUCT::* sr_field,int number_of_threads){
    // Send_ri
    n_subdomains = linearizer_hierarchies.size();
    {LOG::SCOPE scope("copy to linearizer");
        for(int i=0;i<n_subdomains;++i){
            std::pair<const unsigned long*,unsigned> blocks(&(*linearizer_hierarchies[i])[0].offsets_list[0],(*linearizer_hierarchies[i])[0].offsets_list.size());
            if(blocks.second==0) continue;
            //Do we need to clear the buffer?
            //For this one, we send the interface with the interior, because we need the flags in the interface cells. And the data chennel of the interface will be cleared out at the accelarator
            if(number_of_threads)
                SPGrid_Computations::Threading_Helper<T_STRUCT,d>(allocator,blocks).Run_Parallel(SPGrid_Computations::Copy_To_Linearized_Data_With_Interface_Clear<T_STRUCT,T,d,T_offset_ptr>(flags_field,ri_field,(*linearizer_hierarchies[i])[0].data,(*linearizer_hierarchies[i])[0].offsets_map),number_of_threads);
            else
                SPGrid_Computations::Copy_To_Linearized_Data_With_Interface_Clear<T_STRUCT,T,d,T_offset_ptr>(allocator,blocks,
                                                                                                             flags_field,ri_field,
                                                                                                             (*linearizer_hierarchies[i])[0].data,
                                                                                                             (*linearizer_hierarchies[i])[0].offsets_map);}
    }
    {LOG::SCOPE scope("Send buffer");
    for(int accelerator=0;accelerator<accelerators.size();++accelerator){
        accelerators[accelerator].communicator->Send_Command(STEP1);
        for(int subdomain = 0;subdomain < accelerators[accelerator].subdomain_index.size();++subdomain){
            const T_LINEARIZER& top=(*linearizer_hierarchies[accelerators[accelerator].subdomain_index[subdomain]])[0];
            accelerators[accelerator].communicator->Send_Buffer(top.data_buffer_size,top.data);}}}
    
    {LOG::SCOPE scope("Recieve buffer");
    for(int accelerator=0;accelerator<accelerators.size();++accelerator){
        for(int subdomain = 0;subdomain < accelerators[accelerator].subdomain_index.size();++subdomain){
            T_LINEARIZER& top=(*linearizer_hierarchies[accelerators[accelerator].subdomain_index[subdomain]])[0];
            accelerators[accelerator].communicator->Recv_Buffer(top.interface_data_buffer_size,top.interface_data);}}}

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
                                                                                                                        (*linearizer_hierarchies[i])[0].offsets_map_interface);
        }
    }
}
//#####################################################################
// Function Step_Two
//#####################################################################
template<typename T,typename T_STRUCT,int d,typename T_offset_ptr>
void Overseer<T,T_STRUCT,d,T_offset_ptr>::Step_Two(std::vector<std::vector<T_LINEARIZER>*>& linearizer_hierarchies,SPG_Allocator& allocator,unsigned T_STRUCT::* flags_field,T T_STRUCT::* z_field,int number_of_threads){
    n_subdomains = linearizer_hierarchies.size();
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
    for(int accelerator=0;accelerator<accelerators.size();++accelerator){
        accelerators[accelerator].communicator->Send_Command(STEP2);
        for(int subdomain = 0;subdomain < accelerators[accelerator].subdomain_index.size();++subdomain){
            const T_LINEARIZER& top=(*linearizer_hierarchies[accelerators[accelerator].subdomain_index[subdomain]])[0];
            accelerators[accelerator].communicator->Send_Buffer(top.interface_data_buffer_size,top.interface_data);}}

    for(int accelerator=0;accelerator<accelerators.size();++accelerator){
        for(int subdomain = 0;subdomain < accelerators[accelerator].subdomain_index.size();++subdomain){
            const T_LINEARIZER& top=(*linearizer_hierarchies[accelerators[accelerator].subdomain_index[subdomain]])[0];
            accelerators[accelerator].communicator->Recv_Buffer(top.data_buffer_size,top.data);}}

    {LOG::SCOPE scope("copy from linearizer");
    for(int i=0;i<n_subdomains;++i){
        std::pair<const unsigned long*,unsigned> blocks(&(*linearizer_hierarchies[i])[0].offsets_list[0],(*linearizer_hierarchies[i])[0].offsets_list.size());
        if(number_of_threads)
            SPGrid_Computations::Threading_Helper<T_STRUCT,d>(allocator,blocks).Run_Parallel(SPGrid_Computations::Copy_From_Linearized_Data<T_STRUCT,T,d,T_offset_ptr>(flags_field,z_field,(*linearizer_hierarchies[i])[0].data,(*linearizer_hierarchies[i])[0].offsets_map),number_of_threads);
        else
            SPGrid_Computations::Copy_From_Linearized_Data<T_STRUCT,T,d,T_offset_ptr>(allocator,blocks,
                                                                                      flags_field,z_field,
                                                                                      (*linearizer_hierarchies[i])[0].data,
                                                                                      (*linearizer_hierarchies[i])[0].offsets_map);}
    }
}
//#####################################################################
template class Overseer<float,SPGRID_DOMAIN_DECOMPOSITION_DATA<float>,3,unsigned int>;
