//#####################################################################
//Removed For Anonymity: Copyright Authors of Submission pap173s1
//Removed For Anonymity: Copyright Authors of Submission pap173s1
//#####################################################################
#include "Communicator.h"

using namespace Domain_Decomposition;

//#####################################################################
// Function Send_Command
//#####################################################################
void Communicator_CUDA::Send_Command(int number){
    command_queue.push(number);
}
//#####################################################################
// Function Send_Int
//#####################################################################
void Communicator_CUDA::Send_Int(int number){
    integer_queue.push(number);
}
//#####################################################################
// Function Send_Buffer
//#####################################################################
void Communicator_CUDA::Send_Buffer(size_t size,void* buffer){
    data_send_queue.push(std::pair<size_t,void*>(size,buffer));
}
//#####################################################################
// Function Recv_Buffer
//#####################################################################
void Communicator_CUDA::Recv_Buffer(size_t size,void* buffer){
    data_recv_queue.push(std::pair<size_t,void*>(size,buffer));//this is actually a request queue.
    data_recv_confirmation_queue.pop();
}
//#####################################################################

