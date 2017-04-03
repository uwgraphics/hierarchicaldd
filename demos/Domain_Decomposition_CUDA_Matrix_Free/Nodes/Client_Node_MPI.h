//#####################################################################
//Removed For Anonymity: Copyright Authors of Submission pap173s1
//Removed For Anonymity: Copyright Authors of Submission pap173s1
//#####################################################################
// Class Client_Node_MPI
//#####################################################################
#ifndef __CLIENT_NODE_MPI_H__
#define __CLIENT_NODE_MPI_H__

#include <mpi.h>
#include <iostream>
#include <vector>
#include "PHI/Subdomain_Solver_PHI.h"

namespace Domain_Decomposition{

template<class T,typename T_STRUCT, int d,class T_offset_ptr>
class Client_Node_MPI{
    std::vector< Subdomain_Solver_PHI<T,T_STRUCT,d,T_offset_ptr> > subdomain_solvers;
    int rank;
public:
    Client_Node_MPI():rank(-1){}
    ~Client_Node_MPI(){}
    void Initialize(){MPI_Init(0,NULL);MPI_Comm_rank(MPI_COMM_WORLD,&rank);}
    void Run();
    void Dispatch(int command);
    void Initialize_Subdomain(int subdomain);

};
};
#endif
