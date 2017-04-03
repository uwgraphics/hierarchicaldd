#include "../SPGrid_Domain_Decomposition_Solver/SPGRID_DOMAIN_DECOMPOSITION_DATA.h"
#include "../SPGrid_Domain_Decomposition_Solver/MPI/dd_client_cpu.h"
#include <iostream>
#include <stdlib.h>

using namespace SPGrid;
//using namespace PhysBAM;
namespace PhysBAM{int PhysBAM_number_of_threads=0;}

int main( int argc, char** argv ){
    typedef float T;
    typedef SPGRID_DOMAIN_DECOMPOSITION_DATA<T> T_STRUCT;
    typedef unsigned int T_offset_ptr;

    enum{d=3};
    std::cout << "Starting CLIENT." << std::endl;

    DomainDecompositionMPINode_Client_CPU<T,T_STRUCT,d,T_offset_ptr> client_node( atoi(argv[1]), 1 );

    client_node.InitializeMPI( 0, NULL );
    client_node.Run();
    client_node.TeardownMPI();
}
