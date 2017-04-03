//#####################################################################
//Removed For Anonymity: Copyright Authors of Submission pap173s1
//Removed For Anonymity: Copyright Authors of Submission pap173s1
//#####################################################################
// Class Client_Node_CUDA
//#####################################################################
#include "Client_Node_CUDA.h"
#include "Command.h"
#include "../SPGrid_Domain_Decomposition_Solver/SPGRID_DOMAIN_DECOMPOSITION_DATA.h"
#include "../CUDA_Kernels/SPGrid_V_Cycle_Helper_CUDA.h"
#include <iostream>
#include <unistd.h>
using namespace Domain_Decomposition;
using namespace SPGrid;

//#####################################################################
// Function Run
//#####################################################################
template<class T,typename T_STRUCT, int d,class T_offset_ptr> void Client_Node_CUDA<T,T_STRUCT,d,T_offset_ptr>::
Run(){
    cudaSetDevice(gpu_id);
    while(true){
        int command=command_queue.pop();
        Dispatch(command);
    }
    delete this;
}
//#####################################################################
// Function Initialize_Subdomain
//#####################################################################
template<class T,typename T_STRUCT, int d,class T_offset_ptr> void Client_Node_CUDA<T,T_STRUCT,d,T_offset_ptr>::
Initialize_Subdomain(int subdomain){
    int levels=integer_queue.pop();
    subdomain_solvers[subdomain].linearizer_hierarchy.resize(levels);
    for(int level = 0;level < levels;++level){
        subdomain_solvers[subdomain].linearizer_hierarchy[level].number_of_blocks=integer_queue.pop();
        subdomain_solvers[subdomain].linearizer_hierarchy[level].number_of_boundary_blocks=integer_queue.pop();
        subdomain_solvers[subdomain].linearizer_hierarchy[level].number_of_interface_blocks=integer_queue.pop();
        //std::cout << "number_of_blocks: "<<subdomain_solvers[subdomain].linearizer_hierarchy[level].number_of_blocks<<std::endl;
        //std::cout << "number_of_boundary_blocks: "<<subdomain_solvers[subdomain].linearizer_hierarchy[level].number_of_boundary_blocks<<std::endl;
        //std::cout << "number_of_interface_blocks: "<<subdomain_solvers[subdomain].linearizer_hierarchy[level].number_of_interface_blocks<<std::endl;
        //allocation
        subdomain_solvers[subdomain].linearizer_hierarchy[level].Allocate();
        if(subdomain_solvers[subdomain].linearizer_hierarchy[level].number_of_blocks==0) continue;
        std::pair<size_t,void*> buffer;
        buffer=data_send_queue.pop();cudaMemcpy(subdomain_solvers[subdomain].linearizer_hierarchy[level].b        ,buffer.second,buffer.first,cudaMemcpyHostToDevice);
        buffer=data_send_queue.pop();cudaMemcpy(subdomain_solvers[subdomain].linearizer_hierarchy[level].b_x_minus,buffer.second,buffer.first,cudaMemcpyHostToDevice);
        buffer=data_send_queue.pop();cudaMemcpy(subdomain_solvers[subdomain].linearizer_hierarchy[level].b_x_plus ,buffer.second,buffer.first,cudaMemcpyHostToDevice);
        buffer=data_send_queue.pop();cudaMemcpy(subdomain_solvers[subdomain].linearizer_hierarchy[level].b_y_minus,buffer.second,buffer.first,cudaMemcpyHostToDevice);
        buffer=data_send_queue.pop();cudaMemcpy(subdomain_solvers[subdomain].linearizer_hierarchy[level].b_y_plus ,buffer.second,buffer.first,cudaMemcpyHostToDevice);
        buffer=data_send_queue.pop();cudaMemcpy(subdomain_solvers[subdomain].linearizer_hierarchy[level].b_z_minus,buffer.second,buffer.first,cudaMemcpyHostToDevice);
        buffer=data_send_queue.pop();cudaMemcpy(subdomain_solvers[subdomain].linearizer_hierarchy[level].b_z_plus ,buffer.second,buffer.first,cudaMemcpyHostToDevice);
        if(buffer.first!=subdomain_solvers[subdomain].linearizer_hierarchy[level].number_of_blocks*sizeof(T_offset_ptr)){std::cout<<"Incorrect data size."<<std::endl;abort();}
        
        buffer=data_send_queue.pop();cudaMemcpy(subdomain_solvers[subdomain].linearizer_hierarchy[level].b_boundary        ,buffer.second,buffer.first,cudaMemcpyHostToDevice);
        buffer=data_send_queue.pop();cudaMemcpy(subdomain_solvers[subdomain].linearizer_hierarchy[level].b_boundary_x_minus,buffer.second,buffer.first,cudaMemcpyHostToDevice);
        buffer=data_send_queue.pop();cudaMemcpy(subdomain_solvers[subdomain].linearizer_hierarchy[level].b_boundary_x_plus ,buffer.second,buffer.first,cudaMemcpyHostToDevice);
        buffer=data_send_queue.pop();cudaMemcpy(subdomain_solvers[subdomain].linearizer_hierarchy[level].b_boundary_y_minus,buffer.second,buffer.first,cudaMemcpyHostToDevice);
        buffer=data_send_queue.pop();cudaMemcpy(subdomain_solvers[subdomain].linearizer_hierarchy[level].b_boundary_y_plus ,buffer.second,buffer.first,cudaMemcpyHostToDevice);
        buffer=data_send_queue.pop();cudaMemcpy(subdomain_solvers[subdomain].linearizer_hierarchy[level].b_boundary_z_minus,buffer.second,buffer.first,cudaMemcpyHostToDevice);
        buffer=data_send_queue.pop();cudaMemcpy(subdomain_solvers[subdomain].linearizer_hierarchy[level].b_boundary_z_plus ,buffer.second,buffer.first,cudaMemcpyHostToDevice);
        if(buffer.first!=subdomain_solvers[subdomain].linearizer_hierarchy[level].number_of_boundary_blocks*sizeof(T_offset_ptr)){std::cout<<"Incorrect boundary data size."<<std::endl;abort();}

        buffer=data_send_queue.pop();cudaMemcpy(subdomain_solvers[subdomain].linearizer_hierarchy[level].b_interface        ,buffer.second,buffer.first,cudaMemcpyHostToDevice);
        buffer=data_send_queue.pop();cudaMemcpy(subdomain_solvers[subdomain].linearizer_hierarchy[level].b_interface_x_minus,buffer.second,buffer.first,cudaMemcpyHostToDevice);
        buffer=data_send_queue.pop();cudaMemcpy(subdomain_solvers[subdomain].linearizer_hierarchy[level].b_interface_x_plus ,buffer.second,buffer.first,cudaMemcpyHostToDevice);
        buffer=data_send_queue.pop();cudaMemcpy(subdomain_solvers[subdomain].linearizer_hierarchy[level].b_interface_y_minus,buffer.second,buffer.first,cudaMemcpyHostToDevice);
        buffer=data_send_queue.pop();cudaMemcpy(subdomain_solvers[subdomain].linearizer_hierarchy[level].b_interface_y_plus ,buffer.second,buffer.first,cudaMemcpyHostToDevice);
        buffer=data_send_queue.pop();cudaMemcpy(subdomain_solvers[subdomain].linearizer_hierarchy[level].b_interface_z_minus,buffer.second,buffer.first,cudaMemcpyHostToDevice);
        buffer=data_send_queue.pop();cudaMemcpy(subdomain_solvers[subdomain].linearizer_hierarchy[level].b_interface_z_plus ,buffer.second,buffer.first,cudaMemcpyHostToDevice);
        if(buffer.first!=subdomain_solvers[subdomain].linearizer_hierarchy[level].number_of_interface_blocks*sizeof(T_offset_ptr)){std::cout<<"Incorrect interface data size."<<std::endl;abort();}
        
        cudaError err = cudaGetLastError();
        if(err!=cudaSuccess) std::cout<<"Error In Copy Offsets!!!!"<<std::endl;
        
        if(level != levels-1){
            subdomain_solvers[subdomain].linearizer_hierarchy[level].n_prolongation_meta_blocks=integer_queue.pop();
            const unsigned int& n_meta_blocks = subdomain_solvers[subdomain].linearizer_hierarchy[level].n_prolongation_meta_blocks;
            if(cudaMalloc((void**)&subdomain_solvers[subdomain].linearizer_hierarchy[level].prolongation_fine_blocks  ,n_meta_blocks*8*sizeof(T_offset_ptr))!=cudaSuccess) abort();
            if(cudaMalloc((void**)&subdomain_solvers[subdomain].linearizer_hierarchy[level].prolongation_coarse_blocks,n_meta_blocks*8*sizeof(T_offset_ptr))!=cudaSuccess) abort();
            buffer=data_send_queue.pop();cudaMemcpy(subdomain_solvers[subdomain].linearizer_hierarchy[level].prolongation_fine_blocks,buffer.second,buffer.first,cudaMemcpyHostToDevice);
            if(buffer.first!=subdomain_solvers[subdomain].linearizer_hierarchy[level].n_prolongation_meta_blocks*8*sizeof(T_offset_ptr)){std::cout<<"1 Incorrect prolongation data size."<<std::endl;abort();}
            buffer=data_send_queue.pop();cudaMemcpy(subdomain_solvers[subdomain].linearizer_hierarchy[level].prolongation_coarse_blocks,buffer.second,buffer.first,cudaMemcpyHostToDevice);
            if(buffer.first!=subdomain_solvers[subdomain].linearizer_hierarchy[level].n_prolongation_meta_blocks*8*sizeof(T_offset_ptr)){std::cout<<"2 Incorrect prolongation data size."<<std::endl;abort();}
        }
        if(level != 0){
            if(cudaMalloc((void**)&subdomain_solvers[subdomain].linearizer_hierarchy[level].restriction_fine_blocks,
                          subdomain_solvers[subdomain].linearizer_hierarchy[level].number_of_blocks*27*sizeof(T_offset_ptr))!=cudaSuccess) abort();
            buffer=data_send_queue.pop();cudaMemcpy(subdomain_solvers[subdomain].linearizer_hierarchy[level].restriction_fine_blocks,buffer.second,buffer.first,cudaMemcpyHostToDevice);
            if(buffer.first!=subdomain_solvers[subdomain].linearizer_hierarchy[level].number_of_blocks*27*sizeof(T_offset_ptr)){std::cerr<<"Wrong restriction data size for sending!"<<std::endl;abort();}
        }
    }
}
//#####################################################################
// Function Dispatch
//#####################################################################
template<class T,typename T_STRUCT, int d,class T_offset_ptr> void Client_Node_CUDA<T,T_STRUCT,d,T_offset_ptr>::
Dispatch(int command){    
    enum{nVcycles=5,BOUNDARY_SMOOTHING=5,INTERIOR_SMOOTHING=1};    
    typedef SPGrid_Mask<NextLogTwo<sizeof(T_STRUCT)>::value, NextLogTwo<sizeof(T)>::value,d> T_MASK;
    switch(command){
    case INIT:{
        cudaError err = cudaGetLastError();
        if(err!=cudaSuccess) std::cout<<"Error before Recieving Subdomains"<<std::endl;
        int n_subdomains=integer_queue.pop();
        if(n_subdomains == 0) {std::cout<<"No subdomain is assigned to the accelerator"<<std::endl;break;}
        subdomain_solvers.resize(n_subdomains);
        cuda_streams.resize(n_subdomains);
        for(int i = 0;i < MAX_STREAM;++i) cudaStreamCreate(&streams[i]);        
        for(int i = 0;i < n_subdomains;++i) cuda_streams[i]=streams[i%MAX_STREAM];
        std::cout<<n_subdomains<<" subdomains has been assigned to this accelerator"<<std::endl;
        for(int i = 0;i < n_subdomains;++i) Initialize_Subdomain(i);
        cudaDeviceSynchronize(); 
        err = cudaGetLastError();
        if(err!=cudaSuccess) std::cout<<"Error after Recieving Subdomains"<<std::endl;
        std::cout<<"Finished with subdomain initialization"<<std::endl;
    }
        break;
    case SEND_FLAGS:{
        int n_subdomains=subdomain_solvers.size();
        if(n_subdomains == 0) break;
        const int nlevels = subdomain_solvers[0].linearizer_hierarchy.size();
        for(int i = 0;i < n_subdomains;++i)
            for(int level=0;level<nlevels;++level){
                std::pair<size_t,void*> buffer;
                buffer=data_send_queue.pop();
                if(subdomain_solvers[i].linearizer_hierarchy[level].data_buffer_size!=buffer.first) 
                    std::cerr<<"Incorrect data size. Requested: "<<buffer.first<<". Sending: "<<subdomain_solvers[i].linearizer_hierarchy[level].data_buffer_size<<std::endl;
                cudaMemcpy(subdomain_solvers[i].linearizer_hierarchy[level].data,buffer.second,buffer.first,cudaMemcpyHostToDevice);}
        //std::cout<<"Finished Recieving All The Flags"<<std::endl;
        cudaDeviceSynchronize(); 
        cudaError err = cudaGetLastError();
        if(err!=cudaSuccess) std::cout<<"Error after Recieving Flags"<<std::endl;
    }
        break;
    case SET_CHANNELS:{
        int n_subdomains=subdomain_solvers.size();
        std::pair<size_t,void*> buffer;
        if(n_subdomains == 0){
            buffer=data_send_queue.pop();buffer=data_send_queue.pop();buffer=data_send_queue.pop();buffer=data_send_queue.pop();break;}
        buffer=data_send_queue.pop();subdomain_solvers[0].flags_field = *reinterpret_cast<unsigned T_STRUCT::**>(buffer.second);
        buffer=data_send_queue.pop();subdomain_solvers[0].result_field = *reinterpret_cast<T T_STRUCT::**>(buffer.second);
        buffer=data_send_queue.pop();subdomain_solvers[0].rhs_field = *reinterpret_cast<T T_STRUCT::**>(buffer.second);
        buffer=data_send_queue.pop();subdomain_solvers[0].tmp_field = *reinterpret_cast<T T_STRUCT::**>(buffer.second);
        for(int i = 1;i < n_subdomains;++i){
            subdomain_solvers[i].flags_field = subdomain_solvers[0].flags_field;
            subdomain_solvers[i].result_field = subdomain_solvers[0].result_field;
            subdomain_solvers[i].rhs_field = subdomain_solvers[0].rhs_field;
            subdomain_solvers[i].tmp_field = subdomain_solvers[0].tmp_field;}
        //std::cout<<"Finished Setting The Channels."<<std::endl;
        cudaDeviceSynchronize(); 
        cudaError err = cudaGetLastError();
        if(err!=cudaSuccess) std::cout<<"Error after Setting The Channels"<<std::endl;
    }
        break;
    case STEP1:{
        //std::cout<<" starting STEP1!"<<std::endl;
        int n_subdomains=subdomain_solvers.size();
        if(n_subdomains == 0) break;
        {        
            //cudaDeviceSynchronize();
            //cudaEvent_t start,stop;
            //cudaEventCreate(&start);
            //cudaEventCreate(&stop);
            //cudaEventRecord(start);
            for(int i = 0;i < n_subdomains;++i){
                std::pair<size_t,void*> buffer;
                buffer=data_send_queue.pop();cudaMemcpyAsync(subdomain_solvers[i].linearizer_hierarchy[0].data,buffer.second,buffer.first,cudaMemcpyHostToDevice,cuda_streams[i]);
                if(subdomain_solvers[i].linearizer_hierarchy[0].data_buffer_size!=buffer.first) std::cerr<<"Incorrect data size. Requested: "<<buffer.first<<". Recving: "<<subdomain_solvers[i].linearizer_hierarchy[0].data_buffer_size<<std::endl;
            }
            //cudaEventRecord(stop);
            //cudaEventSynchronize(stop);
            //float milliseconds = 0;
            //cudaEventElapsedTime(&milliseconds,start,stop);            
            //std::cout << milliseconds/1000 << " sec for global data transfer. On " << gpu_id << std::endl;
        }
        
        // cudaDeviceSynchronize();
        // cudaEvent_t start,stop;
        // cudaEventCreate(&start);
        // cudaEventCreate(&stop);
        // cudaEventRecord(start);
        {        
            //cudaDeviceSynchronize();
            //cudaEvent_t start,stop;
            //cudaEventCreate(&start);
            //cudaEventCreate(&stop);
            //cudaEventRecord(start);
            for(int i = 0;i < n_subdomains;++i){
                SPGrid_Computations::V_Cycle_Helper_CUDA<T_STRUCT,T,d,T_offset_ptr>::clear_interface(subdomain_solvers[i].flags_field,
                                                                                                     subdomain_solvers[i].rhs_field,
                                                                                                     subdomain_solvers[i].linearizer_hierarchy[0],
                                                                                                     cuda_streams[i]);}

            for(int i = 0;i < n_subdomains;++i){
                SPGrid_Computations::V_Cycle_Helper_CUDA<T_STRUCT,T,d,T_offset_ptr>::clear_interface(subdomain_solvers[i].flags_field,
                                                                                                     subdomain_solvers[i].tmp_field,
                                                                                                     subdomain_solvers[i].linearizer_hierarchy[0],
                                                                                                     cuda_streams[i]);}

            for(int i = 0;i < n_subdomains;++i){
                SPGrid_Computations::V_Cycle_Helper_CUDA<T_STRUCT,T,d,T_offset_ptr>::clear_interface(subdomain_solvers[i].flags_field,
                                                                                                     subdomain_solvers[i].result_field,
                                                                                                     subdomain_solvers[i].linearizer_hierarchy[0],
                                                                                                     cuda_streams[i]);}
            //cudaEventRecord(stop);
            //cudaEventSynchronize(stop);
            //float milliseconds = 0;
            //cudaEventElapsedTime(&milliseconds,start,stop);            
            //std::cout << milliseconds/1000 << " sec for clear_interface. On " << gpu_id << std::endl;
        }


        {        
            //cudaDeviceSynchronize();
            //cudaEvent_t start,stop;
            //cudaEventCreate(&start);
            //cudaEventCreate(&stop);
            //cudaEventRecord(start);
            for(int i = 0;i < n_subdomains;++i){
                SPGrid_Computations::V_Cycle_Helper_CUDA<T_STRUCT,T,d,T_offset_ptr>::clear_u(subdomain_solvers[i].result_field,
                                                                                             subdomain_solvers[i].linearizer_hierarchy[0],
                                                                                             cuda_streams[i]);}      
            //cudaEventRecord(stop);
            //cudaEventSynchronize(stop);
            //float milliseconds = 0;
            //cudaEventElapsedTime(&milliseconds,start,stop);            
            //std::cout << milliseconds/1000 << " sec for clear_u. On " << gpu_id << std::endl;
        }
        
        {        
            // cudaDeviceSynchronize();
            // cudaEvent_t start,stop;
            // cudaEventCreate(&start);
            // cudaEventCreate(&stop);
            // cudaEventRecord(start);
            for(int i = 0;i < n_subdomains;++i){
                for(int j = 0;j < nVcycles;++j){
                    SPGrid_Computations::V_Cycle_Helper_CUDA<T_STRUCT,T,d,T_offset_ptr>::v_cycle(subdomain_solvers[i].result_field,
                                                                                                 subdomain_solvers[i].tmp_field,
                                                                                                 subdomain_solvers[i].rhs_field,
                                                                                                 subdomain_solvers[i].flags_field,
                                                                                                 subdomain_solvers[i].linearizer_hierarchy,cuda_streams[i],
                                                                                                 INTERIOR_SMOOTHING,BOUNDARY_SMOOTHING);}
            }
            // cudaEventRecord(stop);
            // cudaEventSynchronize(stop);
            // float milliseconds = 0;
            // cudaEventElapsedTime(&milliseconds,start,stop);            
            // std::cout << milliseconds/1000 << " sec for v_cycle. On " << gpu_id << std::endl;
        }
        
        
        {        
            //cudaDeviceSynchronize();
            //cudaEvent_t start,stop;
            //cudaEventCreate(&start);
            //cudaEventCreate(&stop);
            //cudaEventRecord(start);
            for(int i = 0;i < n_subdomains;++i){
                SPGrid_Computations::V_Cycle_Helper_CUDA<T_STRUCT,T,d,T_offset_ptr>::clear_u(subdomain_solvers[i].tmp_field,
                                                                                             subdomain_solvers[i].linearizer_hierarchy[0],cuda_streams[i]);}
            
            //cudaEventRecord(stop);
            //cudaEventSynchronize(stop);
            //float milliseconds = 0;
            //cudaEventElapsedTime(&milliseconds,start,stop);            
            //std::cout << milliseconds/1000 << " sec for clear_u. On " << gpu_id<< std::endl;
        }
        
        {        
            //cudaDeviceSynchronize();
            //cudaEvent_t start,stop;
            //cudaEventCreate(&start);
            //cudaEventCreate(&stop);
            //cudaEventRecord(start);
            for(int i = 0;i < n_subdomains;++i){
                SPGrid_Computations::V_Cycle_Helper_CUDA<T_STRUCT,T,d,T_offset_ptr>::interface_interior_laplace(subdomain_solvers[i].result_field,
                                                                                                                subdomain_solvers[i].tmp_field,
                                                                                                                subdomain_solvers[i].flags_field,
                                                                                                                subdomain_solvers[i].linearizer_hierarchy[0],
                                                                                                                cuda_streams[i]);}
            //cudaEventRecord(stop);
            //cudaEventSynchronize(stop);
            //float milliseconds = 0;
            //cudaEventElapsedTime(&milliseconds,start,stop);            
            //std::cout << milliseconds/1000 << " sec for interface_interior_laplace. On " << gpu_id << std::endl;
        }
        
        {        
            //cudaDeviceSynchronize();
            //cudaEvent_t start,stop;
            //cudaEventCreate(&start);
            //cudaEventCreate(&stop);
            //cudaEventRecord(start);

            for(int i = 0;i < n_subdomains;++i){
                SPGrid_Computations::V_Cycle_Helper_CUDA<T_STRUCT,T,d,T_offset_ptr>::collect_interface_blocks(subdomain_solvers[i].linearizer_hierarchy[0],
                                                                                                              cuda_streams[i]);}
            
            
            //cudaEventRecord(stop);
            //cudaEventSynchronize(stop);
            //float milliseconds = 0;
            //cudaEventElapsedTime(&milliseconds,start,stop);            
            //std::cout << milliseconds/1000 << " sec for collect_interface_blocks. On " << gpu_id << std::endl;
        }
        
        // cudaEventRecord(stop);
        // cudaEventSynchronize(stop);
        // float milliseconds = 0;
        // cudaEventElapsedTime(&milliseconds,start,stop);            
        // std::cout << milliseconds/1000 << " sec for Subdomain Solve. On " << gpu_id << std::endl;
        {        
            //cudaDeviceSynchronize();
            //cudaEvent_t start,stop;
            //cudaEventCreate(&start);
            //cudaEventCreate(&stop);
            //cudaEventRecord(start);

            for(int i = 0;i < n_subdomains;++i){
                std::pair<size_t,void*> buffer;
                buffer=data_recv_queue.pop();
                if(subdomain_solvers[i].linearizer_hierarchy[0].interface_data_buffer_size!=buffer.first) std::cerr<<"Incorrect data size. Requested: "<<buffer.first<<". Sending: "<<subdomain_solvers[i].linearizer_hierarchy[0].interface_data_buffer_size<<std::endl;
                cudaMemcpyAsync(buffer.second,subdomain_solvers[i].linearizer_hierarchy[0].interface_data,buffer.first,cudaMemcpyDeviceToHost,cuda_streams[i]);
                cudaStreamSynchronize(cuda_streams[i]);
                data_recv_confirmation_queue.push(true/*we really don't care about the value....just need to push something into it*/);}
            //cudaEventRecord(stop);
            //cudaEventSynchronize(stop);
            //float milliseconds = 0;
            //cudaEventElapsedTime(&milliseconds,start,stop);            
            //std::cout << milliseconds/1000 << " sec for interface data transfer. On " << gpu_id << std::endl;
        }
        
        for(int i = 0;i < n_subdomains;++i){
            SPGrid_Computations::V_Cycle_Helper_CUDA<T_STRUCT,T,d,T_offset_ptr>::clear_u(subdomain_solvers[i].result_field,
                                                                                         subdomain_solvers[i].linearizer_hierarchy[0],cuda_streams[i]);}
        
        //std::cout<<"Finished With Step1!"<<std::endl;
        //sleep(10);//sleep for 10 sec because we don't want the dead spin occupy CPU...
    }
        break;
    case STEP2:{
        //std::cout<<"Starting Step2!"<<std::endl;
        int n_subdomains=subdomain_solvers.size();
        {
            //cudaDeviceSynchronize();
            //cudaEvent_t start,stop;
            //cudaEventCreate(&start);
            //cudaEventCreate(&stop);
            //cudaEventRecord(start);
            for(int i = 0;i < n_subdomains;++i){
                std::pair<size_t,void*> buffer;
                buffer=data_send_queue.pop();
                cudaMemcpyAsync(subdomain_solvers[i].linearizer_hierarchy[0].interface_data,buffer.second,buffer.first,cudaMemcpyHostToDevice,cuda_streams[i]);
                if(subdomain_solvers[i].linearizer_hierarchy[0].interface_data_buffer_size!=buffer.first) std::cerr<<"Incorrect data size"<<std::endl;
            }
            //cudaEventRecord(stop);
            //cudaEventSynchronize(stop);
            //float milliseconds = 0;
            //cudaEventElapsedTime(&milliseconds,start,stop);            
            //std::cout << milliseconds/1000 << " sec for data transfer." << std::endl;
        }

        //cudaDeviceSynchronize();
        {
            //cudaDeviceSynchronize();
            //cudaEvent_t start,stop;
            //cudaEventCreate(&start);
            //cudaEventCreate(&stop);
            //cudaEventRecord(start);
            // for(int i = 0;i < n_subdomains;++i){
            //     SPGrid_Computations::V_Cycle_Helper_CUDA<T_STRUCT,T,d,T_offset_ptr>::clear_u(subdomain_solvers[i].result_field,
            //                                                                                  subdomain_solvers[i].linearizer_hierarchy[0],cuda_streams[i]);}
            //cudaEventRecord(stop);
            //cudaEventSynchronize(stop);
            //float milliseconds = 0;
            //cudaEventElapsedTime(&milliseconds,start,stop);            
            //std::cout << milliseconds/1000 << " sec for clear_u." << std::endl;
        }
            
        
        {
            //cudaDeviceSynchronize();
            //cudaEvent_t start,stop;
            //cudaEventCreate(&start);
            //cudaEventCreate(&stop);
            //cudaEventRecord(start);
            for(int i = 0;i < n_subdomains;++i){
                SPGrid_Computations::V_Cycle_Helper_CUDA<T_STRUCT,T,d,T_offset_ptr>::
                    distribute_interface_blocks_data(subdomain_solvers[i].flags_field,
                                                     subdomain_solvers[i].result_field,
                                                     subdomain_solvers[i].linearizer_hierarchy[0],cuda_streams[i]);}

            //cudaEventRecord(stop);
            //cudaEventSynchronize(stop);
            //float milliseconds = 0;
            //cudaEventElapsedTime(&milliseconds,start,stop);            
            //std::cout << milliseconds/1000 << " sec for distribute_interface_blocks_data." << std::endl;
        }

        {
            //cudaDeviceSynchronize();
            //cudaEvent_t start,stop;
            //cudaEventCreate(&start);
            //cudaEventCreate(&stop);
            //cudaEventRecord(start);
            for(int i = 0;i < n_subdomains;++i){
                SPGrid_Computations::V_Cycle_Helper_CUDA<T_STRUCT,T,d,T_offset_ptr>::
                    interior_interface_accumulative_minus_laplace(subdomain_solvers[i].result_field,
                                                                  subdomain_solvers[i].rhs_field,
                                                                  subdomain_solvers[i].flags_field,
                                                                  subdomain_solvers[i].linearizer_hierarchy[0],cuda_streams[i]);}
            //cudaEventRecord(stop);
            //cudaEventSynchronize(stop);
            //float milliseconds = 0;
            //cudaEventElapsedTime(&milliseconds,start,stop);            
            //std::cout << milliseconds/1000 << " sec for interior_interface_accumulative_minus_laplace." << std::endl;
        }


        {
            //cudaDeviceSynchronize();
            //cudaEvent_t start,stop;
            //cudaEventCreate(&start);
            //cudaEventCreate(&stop);
            //cudaEventRecord(start);

            // for(int i = 0;i < n_subdomains;++i){
            //     //MAY NOT NEEDED!
            //     SPGrid_Computations::V_Cycle_Helper_CUDA<T_STRUCT,T,d,T_offset_ptr>::clear_u(subdomain_solvers[i].result_field,
            //                                                                                  subdomain_solvers[i].linearizer_hierarchy[0],cuda_streams[i]);}


            //cudaEventRecord(stop);
            //cudaEventSynchronize(stop);
            //float milliseconds = 0;
            //cudaEventElapsedTime(&milliseconds,start,stop);            
            //std::cout << milliseconds/1000 << " sec for clear_u." << std::endl;
        }
        {
            for(int i = 0;i < n_subdomains;++i){
                SPGrid_Computations::V_Cycle_Helper_CUDA<T_STRUCT,T,d,T_offset_ptr>::clear_interface(subdomain_solvers[i].flags_field,
                                                                                                     subdomain_solvers[i].rhs_field,
                                                                                                     subdomain_solvers[i].linearizer_hierarchy[0],
                                                                                                     cuda_streams[i]);}

            for(int i = 0;i < n_subdomains;++i){
                SPGrid_Computations::V_Cycle_Helper_CUDA<T_STRUCT,T,d,T_offset_ptr>::clear_interface(subdomain_solvers[i].flags_field,
                                                                                                     subdomain_solvers[i].tmp_field,
                                                                                                     subdomain_solvers[i].linearizer_hierarchy[0],
                                                                                                     cuda_streams[i]);}

            for(int i = 0;i < n_subdomains;++i){
                SPGrid_Computations::V_Cycle_Helper_CUDA<T_STRUCT,T,d,T_offset_ptr>::clear_interface(subdomain_solvers[i].flags_field,
                                                                                                     subdomain_solvers[i].result_field,
                                                                                                     subdomain_solvers[i].linearizer_hierarchy[0],
                                                                                                     cuda_streams[i]);}
        }
        {
            // cudaDeviceSynchronize();
            // cudaEvent_t start,stop;
            // cudaEventCreate(&start);
            // cudaEventCreate(&stop);
            // cudaEventRecord(start);

            for(int i = 0;i < n_subdomains;++i){
                //there should not be any interface data in any of the data chennels.
                for(int j = 0;j < nVcycles;++j){
                    SPGrid_Computations::V_Cycle_Helper_CUDA<T_STRUCT,T,d,T_offset_ptr>::v_cycle(subdomain_solvers[i].result_field,
                                                                                                 subdomain_solvers[i].tmp_field,
                                                                                                 subdomain_solvers[i].rhs_field,
                                                                                                 subdomain_solvers[i].flags_field,
                                                                                                 subdomain_solvers[i].linearizer_hierarchy,cuda_streams[i],
                                                                                                 INTERIOR_SMOOTHING,BOUNDARY_SMOOTHING);}
            }
            // cudaEventRecord(stop);
            // cudaEventSynchronize(stop);
            // float milliseconds = 0;
            // cudaEventElapsedTime(&milliseconds,start,stop);            
            // std::cout << milliseconds/1000 << " sec for v_cycle @ step 1." << std::endl;
        }

        {
            //cudaDeviceSynchronize();
            //cudaEvent_t start,stop;
            //cudaEventCreate(&start);
            //cudaEventCreate(&stop);
            //cudaEventRecord(start);

            for(int i = 0;i < n_subdomains;++i){
                std::pair<size_t,void*> buffer;
                buffer=data_recv_queue.pop();cudaMemcpyAsync(buffer.second,subdomain_solvers[i].linearizer_hierarchy[0].data,buffer.first,cudaMemcpyDeviceToHost,cuda_streams[i]);
                if(subdomain_solvers[i].linearizer_hierarchy[0].data_buffer_size!=buffer.first) std::cerr<<"Incorrect data size"<<std::endl;
                cudaStreamSynchronize(cuda_streams[i]);
                data_recv_confirmation_queue.push(true/*we really don't care about the value....just need to push something into it*/);}
            //cudaEventRecord(stop);
            //cudaEventSynchronize(stop);
            //float milliseconds = 0;
            //cudaEventElapsedTime(&milliseconds,start,stop);            
            //std::cout << milliseconds/1000 << " sec for interface_data_transfer." << std::endl;
        }

    }
        break;
    case SYNC:{
        cudaDeviceSynchronize();
    }
    default: {std::cout << " Unknown command. " << std::endl;}
    }
}
//#####################################################################
template class Client_Node_CUDA<float,SPGRID_DOMAIN_DECOMPOSITION_DATA<float>,3,unsigned int>;
