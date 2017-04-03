#ifndef __SPGRID_Domain_Decomposition_Preconditioner_H__
#define __SPGRID_Domain_Decomposition_Preconditioner_H__
#include <SPGrid/Tools/SPGrid_Block_Iterator.h>
#include <SPGrid/Tools/SPGrid_Copy.h>
#include "SPGrid_V_Cycle_Topology.h"
#include "SPGrid_Block_Pair_Weaver.h"
#include "SPGrid_Linearized_Data_Copy_Helper.h"
#include "SPGrid_Master_Array_Linearizer.h"
#include "SPGrid_Domain_Divider_Helper.h"

#include <Eigen/Sparse>
#include <Eigen/Eigenvalues>

#include <fstream>

#include "../Interface_Solver/INTERFACE_MATRIX_GENERATOR.h"
#include "../Interface_Solver/SIGMA_MULTIGRID_SOLVER.h"
#include "../Interface_Solver/EIGEN_INTERFACE.h"
#include "../KNL_Libraries/KNL_OVERSEER.h"
#include <PhysBAM_Tools/Krylov_Solvers/CONJUGATE_GRADIENT.h>
namespace SPGrid{
using namespace PhysBAM;
template<typename T,typename T_STRUCT,int d,typename T_offset_ptr>
class SPGrid_Domain_Decomposition_Preconditioner{
    typedef SPGrid_Allocator<T_STRUCT,d> SPG_Allocator;
    typedef typename SPG_Allocator::template Array<T>::type SPG_Data_Array_Type;
    typedef typename SPG_Allocator::template Array<const T>::type SPG_Const_Data_Array_Type;
    typedef typename SPG_Allocator::template Array<unsigned>::type SPG_Flags_Array_Type;
    typedef typename SPG_Allocator::template Array<const unsigned>::type SPG_Const_Flags_Array_Type;
    typedef typename SPG_Allocator::template Array<T>::mask T_MASK;
    typedef SPGrid_Set<SPG_Flags_Array_Type> SPG_Set_Type;
    typedef SPGrid_Master_Array_Linearizer<T,NextLogTwo<sizeof(T_STRUCT)>::value,d,T_offset_ptr> T_Linearizer;
    typedef VECTOR<int,d> T_INDEX;
    typedef VECTOR<T,d> TV;
    typedef Eigen::Matrix<T,Eigen::Dynamic,1> T_EIGEN_VECTOR;

public:

    unsigned T_STRUCT::* flags_field;

    INTERFACE_SOLVER_DATA<T,d> matrices;
    std::vector<std::vector<T_Linearizer>*>& linearizer_hierarchy;
    SPG_Allocator& allocator;
    SPG_Set_Type& set;
    //Preconditioner Parameters
    int levels_sigma;
    int n_v_cycles;
    int n_smoothing_sigma;
    int n_sigma_mg;
    int number_of_threads;
    int interior_smoothing;
    int boundary_smoothing;
    //The following ones are temps used in the preconditioner
    ARRAY<T_EIGEN_VECTOR,T_INDEX> vi;
    Domain_Decomposition::KNL_OVERSEER<T,T_STRUCT,d,T_offset_ptr> overseer;
    SIGMA_MULTIGRID_SOLVER<T,d> sigma_mg;
    INTERFACE_MATRIX_GENERATOR<T,d>& interface_matrix_generator;
    std::vector<T_EIGEN_VECTOR,Eigen::aligned_allocator<std::pair<const int,T_EIGEN_VECTOR> > > r;
    std::vector<T_EIGEN_VECTOR,Eigen::aligned_allocator<std::pair<const int,T_EIGEN_VECTOR> > > u;
    std::vector<T_EIGEN_VECTOR,Eigen::aligned_allocator<std::pair<const int,T_EIGEN_VECTOR> > > b;
    SPGrid_Domain_Decomposition_Preconditioner(SPG_Allocator& allocator_input,SPG_Set_Type& set_input,
                                               INTERFACE_MATRIX_GENERATOR<T,d>& interface_matrix_generator_input,
                                               std::vector<std::vector<T_Linearizer>*>& linearizer_hierarchy_input,
                                               T_INDEX subdomain_size,T_INDEX size,unsigned T_STRUCT::* flags_field_input,
                                               int levels_sigma_input,int n_v_cycles_input,int interior_smoothing_input,int boundary_smoothing_input,
                                               int n_smoothing_sigma_input,int n_sigma_mg_input,int number_of_threads_in)
    //NOTE spgrid size is 2^n+1, we need to minus one to get the power of two cartographer wanted
    :allocator(allocator_input),set(set_input),
        interface_matrix_generator(interface_matrix_generator_input),
        linearizer_hierarchy(linearizer_hierarchy_input),
        flags_field(flags_field_input),
        n_v_cycles(n_v_cycles_input),
        interior_smoothing(interior_smoothing_input),
        boundary_smoothing(boundary_smoothing_input),
        n_sigma_mg(n_sigma_mg_input),
        levels_sigma(levels_sigma_input),
        sigma_mg(matrices,u,r,b,n_smoothing_sigma_input),
        number_of_threads(number_of_threads_in)
    {}
    ~SPGrid_Domain_Decomposition_Preconditioner(){}
    void Initialize(bool MPI_enable,bool SCIF_enable,bool CUDA_enable,bool read_from_file=false){
        LOG::cout<<"data size: "<<set.Get_Blocks().second*(double)4096/1024./1024./1024.<<" GB"<<std::endl;
        // initialize MPI stuff
        if(!read_from_file) interface_matrix_generator.Convert_To_Eigen(matrices);
        sigma_mg.Symbolic_Initialize();
        sigma_mg.Numerical_Initialize();
        overseer.Set_Parameters(n_v_cycles,interior_smoothing,boundary_smoothing);
        //sleep(20);
    }
    void Apply_Preconditioner(T T_STRUCT::* result_field,T T_STRUCT::* rhs_field,T T_STRUCT::* tmp_field){
        //std::cout<<"Apply_Preconditioner!"<<std::endl;
        //LOG::SCOPE scope("Apply_Preconditioner"); 
        //int n_subdomains = linearizer_hierarchy.size();
        SPGrid_Computations::Clear<T_STRUCT,T,d>(allocator,set.Get_Blocks(),result_field);
        //This set channel function trigers an send buffer. if this three variables gets distroyed....we are in huge trouble...
        overseer.Set_Channels(flags_field,result_field,rhs_field,tmp_field);
        // just need the interface blocks to be copied....so, should have done something more efficient
        // And it can be done during better times...
        {LOG::SCOPE scope("Copy");
            SPGrid_Computations::Copy<T_STRUCT,T,d>(allocator,set.Get_Blocks(),rhs_field,tmp_field);} 
        
        {LOG::SCOPE scope("Step ONE");
            overseer.Step_One(linearizer_hierarchy,allocator,flags_field,rhs_field,tmp_field,number_of_threads);}
        
        b[0]=T_EIGEN_VECTOR::Zero(b[0].size());
        r[0]=T_EIGEN_VECTOR::Zero(r[0].size());//may not need to
        {LOG::SCOPE scope("Copy to Eigen");
            EIGEN_INTERFACE<T_STRUCT,T,d>::Copy_To_Eigen_Array(b[0],tmp_field,flags_field,
                                                               allocator,set,
                                                               interface_matrix_generator.interface_map);}

        u[0]=T_EIGEN_VECTOR::Zero(u[0].size());
        sigma_mg.smoother.Compute_Residual(0,u[0],r[0],b[0]);
        for(int j = 0;j < r[0].size();++j){
            r[0](j) = fabs(r[0](j));} 
        //std::cout<<"Residual of sigma in SIGMA_SOLVE: "<<r[0].maxCoeff()<<std::endl;
        for(int i=0;i<n_sigma_mg;++i){
            LOG::SCOPE scope("Sigma_V_Cycle");
            sigma_mg.V_Cycle();
            //sigma_mg.smoother.Compute_Residual(0,u[0],r[0],b[0]);
            //for(int j = 0;j < r[0].size();++j){
                //r[0](j) = fabs(r[0](j));}
            //std::cout<<"Residual of sigma in SIGMA_SOLVE: "<<r[0].maxCoeff()<<std::endl;
        }
        {LOG::SCOPE scope("Copy from Eigen");
            EIGEN_INTERFACE<T_STRUCT,T,d>::Copy_From_Eigen_Array(u[0],result_field,flags_field,
                                                                 allocator,set,
                                                                 interface_matrix_generator.interface_map);}
        {LOG::SCOPE scope("Step Two");
            overseer.Step_Two(linearizer_hierarchy,allocator,flags_field,result_field,number_of_threads);}
    }
};
}
#endif
