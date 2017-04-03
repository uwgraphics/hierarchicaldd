//#####################################################################
//Removed For Anonymity: Copyright Authors of Submission pap173s1
//Removed For Anonymity: Copyright Authors of Submission pap173s1
//#####################################################################
// Class CONJUGATE_GRADIENT
//#####################################################################
#ifndef __CONJUGATE_GRADIENT_DUMPING__
#define __CONJUGATE_GRADIENT_DUMPING__

#include <SPGrid_Fluids/Grids/GRID_HIERARCHY_INITIALIZER.h>
#include <SPGrid_Fluids/Projection/GRID_HIERARCHY_PROJECTION.h>
#include <PhysBAM_Tools/Krylov_Solvers/KRYLOV_SOLVER.h>
#include <SPGrid_Fluids/Simulation/FLUIDS_SIMULATION_DATA.h>
namespace PhysBAM{

// see Golub and Van Loan, 10.2.6, p. 529 for details

// cost per iteration = 1 matrix multiply/project/precondition, 2 inner products, 1 convergence norm, 3 saxpy's
// approximate total flops = 2v + 10n
using namespace SPGrid;
template<class T,int d>
class CONJUGATE_GRADIENT_DUMPING:public KRYLOV_SOLVER<T>
{
    typedef VECTOR<T,d> TV;
    STREAM_TYPE stream_type;
    GRID<TV> grid;
    const std::string output_directory;
    
    typedef KRYLOV_SOLVER<T> BASE;
public:
    using BASE::restart_iterations;using BASE::residual_magnitude_squared;using BASE::iterations_used;using BASE::print_diagnostics;using BASE::print_residuals;
    using BASE::nullspace_measure;using BASE::nullspace_tolerance;using BASE::Solve;

    CONJUGATE_GRADIENT_DUMPING(STREAM_TYPE stream_type_input,const std::string output_directory_input)
        :stream_type(stream_type_input),output_directory(output_directory_input)
    {}

    virtual ~CONJUGATE_GRADIENT_DUMPING();

//#####################################################################
    bool Solve(const KRYLOV_SYSTEM_BASE<T>& system,KRYLOV_VECTOR_BASE<T>& x,const KRYLOV_VECTOR_BASE<T>& b,KRYLOV_VECTOR_BASE<T>& q,KRYLOV_VECTOR_BASE<T>& s,KRYLOV_VECTOR_BASE<T>& r,
               KRYLOV_VECTOR_BASE<T>& k,KRYLOV_VECTOR_BASE<T>& z,const T tolerance,const int min_iterations,const int max_iterations) PHYSBAM_OVERRIDE;
//#####################################################################
};
}
#endif
