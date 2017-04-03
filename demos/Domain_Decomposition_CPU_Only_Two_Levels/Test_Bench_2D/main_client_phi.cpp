#include "../SPGrid_Domain_Decomposition_Solver/Nodes/Client_Node_SCIF.h"
#include <iostream>
#include <stdlib.h>
#include "../SPGrid_Domain_Decomposition_Solver/SPGRID_DOMAIN_DECOMPOSITION_DATA.h"

using namespace SPGrid;
using namespace Domain_Decomposition;
typedef float T;
typedef SPGRID_DOMAIN_DECOMPOSITION_DATA<T> T_STRUCT;
typedef unsigned T_offset_ptr;
enum{d=3};


int main( int argc, char** argv ){
    std::cout << "Starting CLIENT." << std::endl;
    Domain_Decomposition::Client_Node_SCIF<T,T_STRUCT,d,T_offset_ptr> client_node;
    client_node.Initialize();
    client_node.Run();
    return 0;
}
