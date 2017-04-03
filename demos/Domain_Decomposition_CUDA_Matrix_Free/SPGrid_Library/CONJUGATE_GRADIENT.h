//#####################################################################
//Removed For Anonymity: Copyright Authors of Submission pap173s1
//Removed For Anonymity: Copyright Authors of Submission pap173s1
//#####################################################################
// Class CONJUGATE_GRADIENT
//#####################################################################
#ifndef __CONJUGATE_GRADIENT__
#define __CONJUGATE_GRADIENT__

#include <PhysBAM_Tools/Krylov_Solvers/KRYLOV_SOLVER.h>


namespace PhysBAM{
// see Golub and Van Loan, 10.2.6, p. 529 for details

// cost per iteration = 1 matrix multiply/project/precondition, 2 inner products, 1 convergence norm, 3 saxpy's
// approximate total flops = 2v + 10n
template<class T,class T_STRUCT,int d>
void Write_Output(STREAM_TYPE stream_type,SPGrid_Allocator<T_STRUCT,d>& allocator,SPGrid_Set<typename SPGrid_Allocator<T_STRUCT,d>::Array<unsigned>::type>& set,GRID<VECTOR<T,d> >& grid,std::string output_directory,int frame);

template<class T,typename T_STRUCT,int d>
class CONJUGATE_GRADIENT:public KRYLOV_SOLVER<T>
{
    typedef KRYLOV_SOLVER<T> BASE;
public:
    using BASE::restart_iterations;using BASE::residual_magnitude_squared;using BASE::iterations_used;using BASE::print_diagnostics;using BASE::print_residuals;
    using BASE::nullspace_measure;using BASE::nullspace_tolerance;using BASE::Solve;
    
    STREAM_TYPE stream_type;
    SPGrid_Allocator<T_STRUCT,d>& allocator;
    SPGrid_Set<typename SPGrid_Allocator<T_STRUCT,d>::Array<unsigned>::type>& set;
    GRID<VECTOR<T,d> >& grid;
    std::string output_directory;
    CONJUGATE_GRADIENT(STREAM_TYPE stream_type_in,SPGrid_Allocator<T_STRUCT,d>& allocator_in,
                       SPGrid_Set<typename SPGrid_Allocator<T_STRUCT,d>::Array<unsigned>::type>& set_in,GRID<VECTOR<T,d> >& grid_in,std::string output_directory_in)
        :stream_type(stream_type_in),allocator(allocator_in),set(set_in),grid(grid_in),output_directory(output_directory_in)
    {}
    
    virtual ~CONJUGATE_GRADIENT(){};

//#####################################################################
    bool Solve(const KRYLOV_SYSTEM_BASE<T>& system,KRYLOV_VECTOR_BASE<T>& x,const KRYLOV_VECTOR_BASE<T>& b,KRYLOV_VECTOR_BASE<T>& q,KRYLOV_VECTOR_BASE<T>& s,KRYLOV_VECTOR_BASE<T>& r,
               KRYLOV_VECTOR_BASE<T>& k,KRYLOV_VECTOR_BASE<T>& z,const T tolerance,const int min_iterations,const int max_iterations) PHYSBAM_OVERRIDE{
    // NOTE: you should never try to make copies of VECTOR_T's inside here as they could be indirect.
    static const T small_number=std::numeric_limits<T>::epsilon();
    system.Set_Boundary_Conditions(x);
    T rho_old=(T)FLT_MAX;T convergence_norm=0;
    int iterations;for(iterations=0;;iterations++){
        bool restart=!iterations || (restart_iterations && iterations%restart_iterations==0);
        if(restart){
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
            if(print_residuals) LOG::cout<<"restarting cg"<<std::endl;
#endif
            r=b;system.Multiply(x,q);r-=q;system.Project(r);}
        // stopping conditions
        system.Project_Nullspace(r);
        
        Write_Output<T,T_STRUCT,d>(stream_type,allocator,set,grid,output_directory,iterations);
        
        convergence_norm=system.Convergence_Norm(r);
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
        if(print_residuals) {LOG::cout<<convergence_norm<<std::endl;}
#endif
        if(convergence_norm<=tolerance && (iterations>=min_iterations || convergence_norm<small_number)){
            if(print_diagnostics) LOG::Stat("cg iterations",iterations);if(iterations_used) *iterations_used=iterations;return true;}
        if(iterations==max_iterations) break;
        // actual iteration
        const KRYLOV_VECTOR_BASE<T>& mr=system.Precondition(r,z);
        T rho=(T)system.Inner_Product(mr,r);
        if(restart) s=mr;
        else s.Copy(rho/rho_old,s,mr);
        system.Multiply(s,q);
        system.Project(q);
        T s_dot_q=(T)system.Inner_Product(s,q);
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
        if(s_dot_q<=0) {LOG::cout<<"CG: matrix appears indefinite or singular, s_dot_q/s_dot_s="<<s_dot_q/(T)system.Inner_Product(s,s)<<std::endl;}
#endif
        T alpha=s_dot_q?rho/s_dot_q:(T)FLT_MAX;
        x.Copy(alpha,s,x);
        r.Copy(-alpha,q,r);
        rho_old=rho;}

#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
    if(print_diagnostics) LOG::Stat("cg iterations",iterations);if(iterations_used) *iterations_used=iterations;
    if(print_diagnostics) {LOG::cout<<"cg not converged after "<<max_iterations<<" iterations, error = "<<convergence_norm<<std::endl;}
#endif
    return false;
    
    }
//#####################################################################
};
}
#endif
