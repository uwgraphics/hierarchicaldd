//#####################################################################
//Removed For Anonymity: Copyright Authors of Submission pap173s1
//Removed For Anonymity: Copyright Authors of Submission pap173s1
//#####################################################################
// Class Client_Node_CUDA
//#####################################################################
#ifndef __CLIENT_NODE_CUDA_H__
#define __CLIENT_NODE_CUDA_H__

#include <thread>
#include <vector>
#include <utility>
#include "Command.h"
#include "Subdomain_Solver_CUDA.h"

namespace Domain_Decomposition{
template<class T,typename T_STRUCT, int d,class T_offset_ptr>
class Client_Node_CUDA{
    std::vector< Subdomain_Solver_CUDA<T,T_STRUCT,d,T_offset_ptr> > subdomain_solvers;
    enum{MAX_STREAM=24};
public:
    const int gpu_id;
    Thread_Safe_Queue<int>& command_queue;
    Thread_Safe_Queue<int>& integer_queue;
    Thread_Safe_Queue<std::pair<size_t,void*> >& data_send_queue;
    Thread_Safe_Queue<std::pair<size_t,void*> >& data_recv_queue;
    Thread_Safe_Queue<bool>& data_recv_confirmation_queue;// this is the queue for confirming data has been recieved from GPU
    std::vector<cudaStream_t> cuda_streams;
    std::thread* worker;
    cudaStream_t streams[MAX_STREAM];
    const bool swap;
    Client_Node_CUDA(int gpu_id_in,
                     Thread_Safe_Queue<int>& command_queue_in,
                     Thread_Safe_Queue<int>& integer_queue_in,
                     Thread_Safe_Queue<std::pair<size_t,void*> >& data_send_queue_in,
                     Thread_Safe_Queue<std::pair<size_t,void*> >& data_recv_queue_in,
                     Thread_Safe_Queue<bool>& data_recv_confirmation_queue_in,
                     bool swap_input=false)
        :gpu_id(gpu_id_in),command_queue(command_queue_in),integer_queue(integer_queue_in),data_send_queue(data_send_queue_in),
        data_recv_queue(data_recv_queue_in),data_recv_confirmation_queue(data_recv_confirmation_queue_in),swap(swap_input)
    {worker = new std::thread(&Client_Node_CUDA::Run,this);}
    ~Client_Node_CUDA(){delete worker;}
    void Initialize(){}
    void Run();
    void Dispatch(int command);
    void Initialize_Subdomain(int subdomain);
};
};

#endif
