//#####################################################################
//Removed For Anonymity: Copyright Authors of Submission pap173s1
//Removed For Anonymity: Copyright Authors of Submission pap173s1
//#####################################################################
#ifndef __COMMUNICATOR_H__
#define __COMMUNICATOR_H__

#include <vector>
#ifdef MPIICC
#include <mpi.h>
#endif
#ifdef __NVCC__
#include <utility>
#include "../CUDA_Kernels/SPGrid_CUDA_Array_Linearizer.h"
#endif
#include "Command.h"

namespace Domain_Decomposition{

enum Accelerator_Type {MPI,SCIF,CUDA};
//#####################################################################
// Class Communicator
//#####################################################################
class Communicator{
public:
    const Accelerator_Type type;
    Communicator(Accelerator_Type type_input):type(type_input){};
    ~Communicator(){}; 
    virtual void Send_Command(int number){};
    virtual void Send_Int(int number){};
    virtual void Send_Buffer(size_t size,void* buffer){};
    virtual void Recv_Buffer(size_t size,void* buffer){};
};
#ifdef MPIICC
//#####################################################################
// Class Communicator_MPI
//#####################################################################
class Communicator_MPI:public Communicator{
public:
    const int rank;
    Communicator_MPI(int rank_in):Communicator(MPI),rank(rank_in){};
    ~Communicator_MPI(){};
    virtual void Send_Command(int number){MPI_Send(&number,1,MPI_INT,rank,COMMAND_TAG,MPI_COMM_WORLD);};
    virtual void Send_Int(int number){MPI_Send(&number,1,MPI_INT,rank,DATA_TAG,MPI_COMM_WORLD);};
    virtual void Send_Buffer(size_t size,void* buffer){MPI_Send(buffer,size,MPI_BYTE,rank,DATA_TAG,MPI_COMM_WORLD);};
    virtual void Recv_Buffer(size_t size,void* buffer){MPI_Status status;MPI_Recv(buffer,size,MPI_BYTE,rank,DATA_TAG,MPI_COMM_WORLD,&status);};
};
#endif

#ifdef __NVCC__
//#####################################################################
// Class Communicator_CUDA
//#####################################################################
class Communicator_CUDA:public Communicator{
public:
    Thread_Safe_Queue<int> command_queue;
    Thread_Safe_Queue<int> integer_queue;
    Thread_Safe_Queue<std::pair<size_t,void*> > data_send_queue;
    Thread_Safe_Queue<std::pair<size_t,void*> > data_recv_queue;
    Thread_Safe_Queue<bool> data_recv_confirmation_queue;// this is the queue for confirming data has been recieved from GPU
    int gpu_id;
    Communicator_CUDA(int gpu_id_in):Communicator(CUDA),gpu_id(gpu_id_in){};
    ~Communicator_CUDA(){};
    virtual void Send_Command(int number);
    virtual void Send_Int(int number);
    virtual void Send_Buffer(size_t size,void* buffer);
    virtual void Recv_Buffer(size_t size,void* buffer);
};
#endif

};

#endif
