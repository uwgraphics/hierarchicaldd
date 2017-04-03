//#####################################################################
//Removed For Anonymity: Copyright Authors of Submission pap173s1
//Removed For Anonymity: Copyright Authors of Submission pap173s1
//#####################################################################
// Class Client_Node_MPI
//#####################################################################

#include "Client_Node_MPI.h"
#include "Command.h"
#include "../SPGrid_Domain_Decomposition_Solver/SPGRID_DOMAIN_DECOMPOSITION_DATA.h"
#include "../PHI_Kernels/SPGrid_V_Cycle_Helper_Phi.h"

using namespace SPGrid;
using namespace Domain_Decomposition;
//#####################################################################
// Function Run
//#####################################################################
template<class T,typename T_STRUCT, int d,class T_offset_ptr> void Client_Node_MPI<T,T_STRUCT,d,T_offset_ptr>::
Run(){
    typedef int command_t;
    command_t work;
    MPI_Status status;
    bool running = true;
    while(running){
        MPI_Recv(&work,1,MPI_INT,0,MPI_ANY_TAG,MPI_COMM_WORLD, &status);
        switch(status.MPI_TAG){
        case COMMAND_TAG:{Dispatch(work);}
            break;
        case TERMINATE_TAG:{running = false;}
            break;
        default:{
            std::cout << "Got unexpected input." << std::endl;
            std::cout << "Source: " << status.MPI_SOURCE;
            std::cout << "Tag: " << status.MPI_TAG;
            std::cout << "Data: " << work;}}
    };
}
//#####################################################################
// Function Initialize_Subdomain
//#####################################################################
template<class T,typename T_STRUCT, int d,class T_offset_ptr> void Client_Node_MPI<T,T_STRUCT,d,T_offset_ptr>::
Initialize_Subdomain(int subdomain){
    MPI_Status status;
    int levels;
    MPI_Recv( (void*)&levels, 1, MPI_INT, 0, DATA_TAG, MPI_COMM_WORLD, &status );
    subdomain_solvers[subdomain].linearizer_hierarchy.resize(levels);
    for(int level = 0;level < levels;++level){
        MPI_Recv((void*)&(subdomain_solvers[subdomain].linearizer_hierarchy[level].number_of_blocks), 1, MPI_INT, 0, DATA_TAG, MPI_COMM_WORLD, &status);
        MPI_Recv((void*)&(subdomain_solvers[subdomain].linearizer_hierarchy[level].number_of_boundary_blocks), 1, MPI_INT, 0, DATA_TAG, MPI_COMM_WORLD, &status);
        MPI_Recv((void*)&(subdomain_solvers[subdomain].linearizer_hierarchy[level].number_of_interface_blocks), 1, MPI_INT, 0, DATA_TAG, MPI_COMM_WORLD, &status);
        //allocation
        subdomain_solvers[subdomain].linearizer_hierarchy[level].Allocate();
        if(subdomain_solvers[subdomain].linearizer_hierarchy[level].number_of_blocks==0) continue;
        MPI_Recv((void*)subdomain_solvers[subdomain].linearizer_hierarchy[level].b        ,subdomain_solvers[subdomain].linearizer_hierarchy[level].number_of_blocks*sizeof(T_offset_ptr),MPI_BYTE,0,DATA_TAG,MPI_COMM_WORLD,&status);   
        MPI_Recv((void*)subdomain_solvers[subdomain].linearizer_hierarchy[level].b_x_minus,subdomain_solvers[subdomain].linearizer_hierarchy[level].number_of_blocks*sizeof(T_offset_ptr),MPI_BYTE,0,DATA_TAG,MPI_COMM_WORLD,&status);   
        MPI_Recv((void*)subdomain_solvers[subdomain].linearizer_hierarchy[level].b_x_plus ,subdomain_solvers[subdomain].linearizer_hierarchy[level].number_of_blocks*sizeof(T_offset_ptr),MPI_BYTE,0,DATA_TAG,MPI_COMM_WORLD,&status);   
        MPI_Recv((void*)subdomain_solvers[subdomain].linearizer_hierarchy[level].b_y_minus,subdomain_solvers[subdomain].linearizer_hierarchy[level].number_of_blocks*sizeof(T_offset_ptr),MPI_BYTE,0,DATA_TAG,MPI_COMM_WORLD,&status);   
        MPI_Recv((void*)subdomain_solvers[subdomain].linearizer_hierarchy[level].b_y_plus ,subdomain_solvers[subdomain].linearizer_hierarchy[level].number_of_blocks*sizeof(T_offset_ptr),MPI_BYTE,0,DATA_TAG,MPI_COMM_WORLD,&status);   
        MPI_Recv((void*)subdomain_solvers[subdomain].linearizer_hierarchy[level].b_z_minus,subdomain_solvers[subdomain].linearizer_hierarchy[level].number_of_blocks*sizeof(T_offset_ptr),MPI_BYTE,0,DATA_TAG,MPI_COMM_WORLD,&status);   
        MPI_Recv((void*)subdomain_solvers[subdomain].linearizer_hierarchy[level].b_z_plus ,subdomain_solvers[subdomain].linearizer_hierarchy[level].number_of_blocks*sizeof(T_offset_ptr),MPI_BYTE,0,DATA_TAG,MPI_COMM_WORLD,&status);

        MPI_Recv((void*)subdomain_solvers[subdomain].linearizer_hierarchy[level].b_boundary        ,subdomain_solvers[subdomain].linearizer_hierarchy[level].number_of_boundary_blocks*sizeof(T_offset_ptr),MPI_BYTE,0,DATA_TAG,MPI_COMM_WORLD,&status);   
        MPI_Recv((void*)subdomain_solvers[subdomain].linearizer_hierarchy[level].b_boundary_x_minus,subdomain_solvers[subdomain].linearizer_hierarchy[level].number_of_boundary_blocks*sizeof(T_offset_ptr),MPI_BYTE,0,DATA_TAG,MPI_COMM_WORLD,&status);   
        MPI_Recv((void*)subdomain_solvers[subdomain].linearizer_hierarchy[level].b_boundary_x_plus ,subdomain_solvers[subdomain].linearizer_hierarchy[level].number_of_boundary_blocks*sizeof(T_offset_ptr),MPI_BYTE,0,DATA_TAG,MPI_COMM_WORLD,&status);   
        MPI_Recv((void*)subdomain_solvers[subdomain].linearizer_hierarchy[level].b_boundary_y_minus,subdomain_solvers[subdomain].linearizer_hierarchy[level].number_of_boundary_blocks*sizeof(T_offset_ptr),MPI_BYTE,0,DATA_TAG,MPI_COMM_WORLD,&status);   
        MPI_Recv((void*)subdomain_solvers[subdomain].linearizer_hierarchy[level].b_boundary_y_plus ,subdomain_solvers[subdomain].linearizer_hierarchy[level].number_of_boundary_blocks*sizeof(T_offset_ptr),MPI_BYTE,0,DATA_TAG,MPI_COMM_WORLD,&status);   
        MPI_Recv((void*)subdomain_solvers[subdomain].linearizer_hierarchy[level].b_boundary_z_minus,subdomain_solvers[subdomain].linearizer_hierarchy[level].number_of_boundary_blocks*sizeof(T_offset_ptr),MPI_BYTE,0,DATA_TAG,MPI_COMM_WORLD,&status);   
        MPI_Recv((void*)subdomain_solvers[subdomain].linearizer_hierarchy[level].b_boundary_z_plus ,subdomain_solvers[subdomain].linearizer_hierarchy[level].number_of_boundary_blocks*sizeof(T_offset_ptr),MPI_BYTE,0,DATA_TAG,MPI_COMM_WORLD,&status);   

        MPI_Recv((void*)subdomain_solvers[subdomain].linearizer_hierarchy[level].b_interface        ,subdomain_solvers[subdomain].linearizer_hierarchy[level].number_of_interface_blocks*sizeof(T_offset_ptr),MPI_BYTE,0,DATA_TAG,MPI_COMM_WORLD,&status);   
        MPI_Recv((void*)subdomain_solvers[subdomain].linearizer_hierarchy[level].b_interface_x_minus,subdomain_solvers[subdomain].linearizer_hierarchy[level].number_of_interface_blocks*sizeof(T_offset_ptr),MPI_BYTE,0,DATA_TAG,MPI_COMM_WORLD,&status);   
        MPI_Recv((void*)subdomain_solvers[subdomain].linearizer_hierarchy[level].b_interface_x_plus ,subdomain_solvers[subdomain].linearizer_hierarchy[level].number_of_interface_blocks*sizeof(T_offset_ptr),MPI_BYTE,0,DATA_TAG,MPI_COMM_WORLD,&status);   
        MPI_Recv((void*)subdomain_solvers[subdomain].linearizer_hierarchy[level].b_interface_y_minus,subdomain_solvers[subdomain].linearizer_hierarchy[level].number_of_interface_blocks*sizeof(T_offset_ptr),MPI_BYTE,0,DATA_TAG,MPI_COMM_WORLD,&status);   
        MPI_Recv((void*)subdomain_solvers[subdomain].linearizer_hierarchy[level].b_interface_y_plus ,subdomain_solvers[subdomain].linearizer_hierarchy[level].number_of_interface_blocks*sizeof(T_offset_ptr),MPI_BYTE,0,DATA_TAG,MPI_COMM_WORLD,&status);   
        MPI_Recv((void*)subdomain_solvers[subdomain].linearizer_hierarchy[level].b_interface_z_minus,subdomain_solvers[subdomain].linearizer_hierarchy[level].number_of_interface_blocks*sizeof(T_offset_ptr),MPI_BYTE,0,DATA_TAG,MPI_COMM_WORLD,&status);   
        MPI_Recv((void*)subdomain_solvers[subdomain].linearizer_hierarchy[level].b_interface_z_plus ,subdomain_solvers[subdomain].linearizer_hierarchy[level].number_of_interface_blocks*sizeof(T_offset_ptr),MPI_BYTE,0,DATA_TAG,MPI_COMM_WORLD,&status);   

        if(level != levels-1){
            int n_meta_blocks=0;
            MPI_Recv((void*)&(n_meta_blocks), 1, MPI_INT, 0, DATA_TAG, MPI_COMM_WORLD, &status);
            subdomain_solvers[subdomain].linearizer_hierarchy[level].prolongation_fine_blocks.resize(n_meta_blocks);
            subdomain_solvers[subdomain].linearizer_hierarchy[level].prolongation_coarse_blocks.resize(n_meta_blocks);
            MPI_Recv((void*)subdomain_solvers[subdomain].linearizer_hierarchy[level].prolongation_fine_blocks[0].data,n_meta_blocks*8*sizeof(T_offset_ptr),MPI_BYTE,0,DATA_TAG,MPI_COMM_WORLD,&status);   
            MPI_Recv((void*)subdomain_solvers[subdomain].linearizer_hierarchy[level].prolongation_coarse_blocks[0].data,n_meta_blocks*8*sizeof(T_offset_ptr),MPI_BYTE,0,DATA_TAG,MPI_COMM_WORLD,&status);}  
        if(level != 0){
            subdomain_solvers[subdomain].linearizer_hierarchy[level].restriction_fine_blocks.resize(subdomain_solvers[subdomain].linearizer_hierarchy[level].number_of_blocks);
            MPI_Recv((void*)subdomain_solvers[subdomain].linearizer_hierarchy[level].restriction_fine_blocks[0].data,subdomain_solvers[subdomain].linearizer_hierarchy[level].number_of_blocks*27*sizeof(T_offset_ptr),MPI_BYTE,0,DATA_TAG,MPI_COMM_WORLD,&status);}
    }
}
//#####################################################################
// Function Dispatch
//#####################################################################
template<class T,typename T_STRUCT, int d,class T_offset_ptr> void Client_Node_MPI<T,T_STRUCT,d,T_offset_ptr>::
Dispatch(int command){
    enum{nVcycles=2,BOUNDARY_SMOOTHING=3,INTERIOR_SMOOTHING=1};    
    switch(command){
    case INIT:{
        int n_subdomains=0;
        MPI_Status status;
        MPI_Recv((void*)&n_subdomains,1,MPI_INT,0,DATA_TAG,MPI_COMM_WORLD, &status);
        if(n_subdomains == 0) {std::cout<<"No subdomain is assigned to the accelerator"<<std::endl;abort();}
        subdomain_solvers.resize(n_subdomains);
        std::cout<<n_subdomains<<" subdomains has been assigned to this accelerator"<<std::endl;
        for(int i = 0;i < n_subdomains;++i) Initialize_Subdomain(i);
        std::cout<<"Finished Recieving All The Subdomains"<<std::endl;
    }
        break;
    case SEND_FLAGS:{
        MPI_Status status;
        int n_subdomains=subdomain_solvers.size();
        if(n_subdomains == 0) break;
        const int nlevels = subdomain_solvers[0].linearizer_hierarchy.size();
        for(int i = 0;i < n_subdomains;++i)
            for(int level=0;level<nlevels;++level)
                MPI_Recv((void*)subdomain_solvers[i].linearizer_hierarchy[level].data,
                         subdomain_solvers[i].linearizer_hierarchy[level].data_buffer_size,MPI_BYTE,0,DATA_TAG,MPI_COMM_WORLD,&status);
        std::cout<<"Finished Recieving All The Flags"<<std::endl;        
    }
        break; 
    case SET_CHANNELS:{
        MPI_Status status;
        int n_subdomains=subdomain_solvers.size();
        if(n_subdomains == 0) break;
        MPI_Recv(reinterpret_cast<void*>(&subdomain_solvers[0].flags_field),sizeof(size_t),MPI_BYTE,0,DATA_TAG,MPI_COMM_WORLD,&status);
        MPI_Recv(reinterpret_cast<void*>(&subdomain_solvers[0].result_field),sizeof(size_t),MPI_BYTE,0,DATA_TAG,MPI_COMM_WORLD,&status);
        MPI_Recv(reinterpret_cast<void*>(&subdomain_solvers[0].rhs_field),sizeof(size_t),MPI_BYTE,0,DATA_TAG,MPI_COMM_WORLD,&status);
        MPI_Recv(reinterpret_cast<void*>(&subdomain_solvers[0].tmp_field),sizeof(size_t),MPI_BYTE,0,DATA_TAG,MPI_COMM_WORLD,&status);
        for(int i = 1;i < n_subdomains;++i){
            subdomain_solvers[i].flags_field = subdomain_solvers[0].flags_field;
            subdomain_solvers[i].result_field = subdomain_solvers[0].result_field;
            subdomain_solvers[i].rhs_field = subdomain_solvers[0].rhs_field;
            subdomain_solvers[i].tmp_field = subdomain_solvers[0].tmp_field;}
        std::cout<<"Finished Setting The Channels."<<std::endl;
    }
       break;
    case STEP1:{
        MPI_Status status;
        int n_subdomains=subdomain_solvers.size();
        for(int i = 0;i < n_subdomains;++i){
            MPI_Recv((void*)subdomain_solvers[i].linearizer_hierarchy[0].data,
                     subdomain_solvers[i].linearizer_hierarchy[0].data_buffer_size,MPI_BYTE,0,DATA_TAG,MPI_COMM_WORLD,&status);}
        for(int i = 0;i < n_subdomains;++i){            
            SPGrid_Computations::V_Cycle_Helper_PHI<T_STRUCT,T,d,T_offset_ptr>::clear_interface(subdomain_solvers[i].flags_field,
                                                                                                subdomain_solvers[i].rhs_field,
                                                                                                subdomain_solvers[i].linearizer_hierarchy[0]);
        
            SPGrid_Computations::V_Cycle_Helper_PHI<T_STRUCT,T,d,T_offset_ptr>::clear_u(subdomain_solvers[i].result_field,
                                                                                        subdomain_solvers[i].linearizer_hierarchy[0]);        
            for(int j = 0;j < nVcycles;++j){
                SPGrid_Computations::V_Cycle_Helper_PHI<T_STRUCT,T,d,T_offset_ptr>::v_cycle(subdomain_solvers[i].result_field,
                                                                                            subdomain_solvers[i].tmp_field,
                                                                                            subdomain_solvers[i].rhs_field,
                                                                                            subdomain_solvers[i].flags_field,
                                                                                            subdomain_solvers[i].linearizer_hierarchy,
                                                                                            INTERIOR_SMOOTHING,BOUNDARY_SMOOTHING);}
            // MAY NOT NEEDED!!!!!!!!!!!! We do this to make sure nothing on the interface of the result. 
            SPGrid_Computations::V_Cycle_Helper_PHI<T_STRUCT,T,d,T_offset_ptr>::clear_u(subdomain_solvers[i].tmp_field,
                                                                                        subdomain_solvers[i].linearizer_hierarchy[0]);        
            // NOTE: after this, tmp_field contain gabage content. And the minus equal happens at the host end. This is just result = Lu on the interface.
            // Next step we are going to collect interface from host and store in result. So, all is good! Hopefully.....
            SPGrid_Computations::V_Cycle_Helper_PHI<T_STRUCT,T,d,T_offset_ptr>::interface_interior_laplace(subdomain_solvers[i].result_field,
                                                                                                           subdomain_solvers[i].tmp_field,
                                                                                                           subdomain_solvers[i].flags_field,
                                                                                                           subdomain_solvers[i].linearizer_hierarchy[0]);
            
            SPGrid_Computations::V_Cycle_Helper_PHI<T_STRUCT,T,d,T_offset_ptr>::collect_interface_blocks(subdomain_solvers[i].linearizer_hierarchy[0]);
        }
        for(int i = 0;i < n_subdomains;++i){
            MPI_Send((void*)subdomain_solvers[i].linearizer_hierarchy[0].interface_data,
                     subdomain_solvers[i].linearizer_hierarchy[0].interface_data_buffer_size,MPI_BYTE,0,DATA_TAG,MPI_COMM_WORLD);}

        for(int i = 0;i < n_subdomains;++i){
            SPGrid_Computations::V_Cycle_Helper_PHI<T_STRUCT,T,d,T_offset_ptr>::clear_u(subdomain_solvers[i].result_field,
                                                                                        subdomain_solvers[i].linearizer_hierarchy[0]);}

    }
        break;
    case STEP2:{
        MPI_Status status;
        int n_subdomains=subdomain_solvers.size();
        for(int i = 0;i < n_subdomains;++i){
            MPI_Recv((void*)subdomain_solvers[i].linearizer_hierarchy[0].interface_data,
                     subdomain_solvers[i].linearizer_hierarchy[0].interface_data_buffer_size,MPI_BYTE,0,DATA_TAG,MPI_COMM_WORLD,&status);}

        for(int i = 0;i < n_subdomains;++i){
            SPGrid_Computations::V_Cycle_Helper_PHI<T_STRUCT,T,d,T_offset_ptr>::distribute_interface_blocks_data(subdomain_solvers[i].flags_field,
                                                                                                                 subdomain_solvers[i].result_field,
                                                                                                                 subdomain_solvers[i].linearizer_hierarchy[0]);
        
            SPGrid_Computations::V_Cycle_Helper_PHI<T_STRUCT,T,d,T_offset_ptr>::interior_interface_accumulative_minus_laplace(subdomain_solvers[i].result_field,
                                                                                                                              subdomain_solvers[i].rhs_field,
                                                                                                                              subdomain_solvers[i].flags_field,
                                                                                                                              subdomain_solvers[i].linearizer_hierarchy[0]); 
            //MAY NOT NEEDED!
            SPGrid_Computations::V_Cycle_Helper_PHI<T_STRUCT,T,d,T_offset_ptr>::clear_u(subdomain_solvers[i].result_field,
                                                                                        subdomain_solvers[i].linearizer_hierarchy[0]);
            //there should not be any interface data in any of the data chennels.
            for(int j = 0;j < nVcycles;++j)
                SPGrid_Computations::V_Cycle_Helper_PHI<T_STRUCT,T,d,T_offset_ptr>::v_cycle(subdomain_solvers[i].result_field,
                                                                                            subdomain_solvers[i].tmp_field,
                                                                                            subdomain_solvers[i].rhs_field,
                                                                                            subdomain_solvers[i].flags_field,
                                                                                            subdomain_solvers[i].linearizer_hierarchy,
                                                                                            INTERIOR_SMOOTHING,BOUNDARY_SMOOTHING);

        }
        for(int i = 0;i < n_subdomains;++i){
            MPI_Send((void*)subdomain_solvers[i].linearizer_hierarchy[0].data,
                     subdomain_solvers[i].linearizer_hierarchy[0].data_buffer_size,MPI_BYTE,0,DATA_TAG,MPI_COMM_WORLD);}
    }
        break;
    default: {std::cout << " Unknown command. " << std::endl;}

        
    }
}
//#####################################################################
template class Client_Node_MPI<float,SPGRID_DOMAIN_DECOMPOSITION_DATA<float>,3,unsigned int>;
