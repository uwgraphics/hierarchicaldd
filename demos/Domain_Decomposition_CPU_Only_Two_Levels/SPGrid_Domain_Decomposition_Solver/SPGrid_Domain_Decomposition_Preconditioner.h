#ifndef __SPGRID_Domain_Decomposition_Preconditioner_H__
#define __SPGRID_Domain_Decomposition_Preconditioner_H__
#include <SPGrid/Tools/SPGrid_Block_Iterator.h>
#include <SPGrid/Tools/SPGrid_Copy.h>
#include "SPGrid_V_Cycle_Topology.h"
#include "SPGrid_Block_Pair_Weaver.h"
#include "SPGrid_Domain_Divider_Helper.h"
#include "../SPGrid_Library/SPGrid_V_Cycle_Helper.h"
#include "../SPGrid_Library/SPGrid_Laplace_Legacy.h"
#include <Eigen/Sparse>
#include <Eigen/Eigenvalues>

#include <fstream>

#include "../Interface_Solver/INTERFACE_MATRIX_GENERATOR.h"
#include "../Interface_Solver/SIGMA_MULTIGRID_SOLVER.h"
#include "../Interface_Solver/EIGEN_INTERFACE.h"

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
    typedef VECTOR<int,d> T_INDEX;
    typedef VECTOR<T,d> TV;
    typedef Eigen::Matrix<T,Eigen::Dynamic,1> T_EIGEN_VECTOR;

public:
    SPGrid_Computations::V_Cycle_Topology<T_STRUCT,T,d>* v_cycle_topology;
    // We need a second Multigrid hierarchy because when doing two level DD, we need to enforce non-zero dirichlet condition,
    // That requires them to be moved to the rhs in the interface solve(interface solve always assumes zero-dirichlet), 
    // it means that we have to modify the rhs. In the preconditioner, the rhs is fixed, and can not be touched... 
    SPGrid_Computations::V_Cycle_Topology<T_STRUCT,T,d>* v_cycle_topology_backup;


    unsigned T_STRUCT::* flags_field;

    INTERFACE_SOLVER_DATA<T,d> matrices;
    std::vector<INTERFACE_SOLVER_DATA<T,d>*> matrices_second_level;

    SPG_Allocator& allocator;
    SPG_Set_Type& set;
    int max_subdomain_id;
    //Preconditioner Parameters
    int levels_sigma;
    int n_v_cycles;
    int n_smoothing_sigma;
    int n_sigma_mg;
    int number_of_threads;

    //The following ones are temps used in the preconditioner
    ARRAY<T_EIGEN_VECTOR,T_INDEX> vi;

    SIGMA_MULTIGRID_SOLVER<T,d> sigma_mg;
    std::vector<SIGMA_MULTIGRID_SOLVER<T,d>*> sigma_mg_second_level;
    std::vector<std::vector<T_EIGEN_VECTOR> > r_second_level;
    std::vector<std::vector<T_EIGEN_VECTOR> > u_second_level;
    std::vector<std::vector<T_EIGEN_VECTOR> > b_second_level;

    INTERFACE_MATRIX_GENERATOR<T,d>& interface_matrix_generator;
    std::vector<INTERFACE_MATRIX_GENERATOR<T,d>*>* interface_matrix_generators_second_level;
    bool multilevel;
    std::vector<T_EIGEN_VECTOR> r;
    std::vector<T_EIGEN_VECTOR> u;
    std::vector<T_EIGEN_VECTOR> b;
    std_array<int,d> subdomain_size;
    int interior_smoothing;
    int boundary_smoothing;
    SPGrid_Domain_Decomposition_Preconditioner(SPG_Allocator& allocator_input,SPG_Set_Type& set_input,
                                               INTERFACE_MATRIX_GENERATOR<T,d>& interface_matrix_generator_input,
                                               T_INDEX subdomain_size_input,T_INDEX size,unsigned T_STRUCT::* flags_field_input,
                                               int n_subdomains,int levels_sigma_input,int n_v_cycles_input,
                                               int interior_smoothing_input,int boundary_smoothing_input,
                                               int n_smoothing_sigma_input,int n_sigma_mg_input,int number_of_threads_in,
                                               std::vector<INTERFACE_MATRIX_GENERATOR<T,d>*>* interface_matrix_generators_second_level_input=NULL,
                                               bool multilevel_input=false)
    //NOTE spgrid size is 2^n+1, we need to minus one to get the power of two cartographer wanted
    :allocator(allocator_input),set(set_input),
        interface_matrix_generator(interface_matrix_generator_input),
        flags_field(flags_field_input),
        n_v_cycles(n_v_cycles_input),
        interior_smoothing(interior_smoothing_input),
        boundary_smoothing(boundary_smoothing_input),
        n_sigma_mg(n_sigma_mg_input),max_subdomain_id(n_subdomains),
        levels_sigma(levels_sigma_input),
        sigma_mg(matrices,u,r,b,n_smoothing_sigma_input),
        n_smoothing_sigma(n_smoothing_sigma_input),
        number_of_threads(number_of_threads_in),
        interface_matrix_generators_second_level(interface_matrix_generators_second_level_input),
        multilevel(multilevel_input),
        subdomain_size(subdomain_size_input)
    {
        LOG::cout<<"number of subdomains: "<<n_subdomains<<std::endl;
    }
    ~SPGrid_Domain_Decomposition_Preconditioner(){
        for(int i=0;i<matrices_second_level.size();++i) delete matrices_second_level[i];
        for(int i=0;i<sigma_mg_second_level.size();++i) delete sigma_mg_second_level[i];
    }
    void Initialize(bool MPI_enable,bool SCIF_enable,bool CUDA_enable){
        LOG::cout<<"data size: "<<set.Get_Blocks().second*(double)4096/1024./1024./1024.<<" GB"<<std::endl;
        // initialize MPI stuff
        interface_matrix_generator.Convert_To_Eigen(matrices);
        sigma_mg.Symbolic_Initialize();
        sigma_mg.Numerical_Initialize();
        if(multilevel){
            for(int i=0;i<matrices_second_level.size();++i) delete matrices_second_level[i];
            matrices_second_level.resize(interface_matrix_generators_second_level->size());
            sigma_mg_second_level.resize(interface_matrix_generators_second_level->size());
            r_second_level.resize(interface_matrix_generators_second_level->size());
            u_second_level.resize(interface_matrix_generators_second_level->size());
            b_second_level.resize(interface_matrix_generators_second_level->size());
            for(int i=0;i<matrices_second_level.size();++i){
                matrices_second_level[i]=new INTERFACE_SOLVER_DATA<T,d>();
                (*interface_matrix_generators_second_level)[i]->Convert_To_Eigen(*matrices_second_level[i]);
                sigma_mg_second_level[i]=new SIGMA_MULTIGRID_SOLVER<T,d>(*matrices_second_level[i],u_second_level[i],r_second_level[i],b_second_level[i],n_smoothing_sigma);
                sigma_mg_second_level[i]->Symbolic_Initialize();
                sigma_mg_second_level[i]->Numerical_Initialize();}}
        //sleep(20);
    }
    void Subdomain_DD(T T_STRUCT::* result_field,T T_STRUCT::* rhs_field,T T_STRUCT::* tmp_field){
        for(int i=0;i<sigma_mg_second_level.size();++i){
            b_second_level[i][0]=T_EIGEN_VECTOR::Zero(b_second_level[i][0].size());}
        SPGrid_Computations::Clear<T_STRUCT,T,d>(allocator,set.Get_Blocks(),tmp_field);
        SPGrid_Computations::Laplace_ir_Subdomain<T_STRUCT,T,unsigned,d>(allocator,set.Get_Blocks(),result_field,tmp_field,flags_field,subdomain_size);
        for(int i=0;i<sigma_mg_second_level.size();++i){
            EIGEN_INTERFACE<T_STRUCT,T,d>::Copy_To_Eigen_Array(b_second_level[i][0],tmp_field,flags_field,
                                                               allocator,set,
                                                               (*interface_matrix_generators_second_level)[i]->interface_map);}
        
        for(int i=0;i<n_v_cycles;++i)
            SPGrid_Computations::V_Cycle_Helper<T_STRUCT,T,d>::v_cycle(result_field,tmp_field,rhs_field,flags_field,*v_cycle_topology,number_of_threads,interior_smoothing,boundary_smoothing);
        SPGrid_Computations::Copy<T_STRUCT,T,d>(allocator,set.Get_Blocks(),rhs_field,tmp_field);
        // Note: This A_ri computation will contaminate the tmp_channel of the first level interface....But, we clear them on the top level DD.
        SPGrid_Computations::Laplace_ri<T_STRUCT,T,unsigned,d>(allocator,set.Get_Blocks(),result_field,tmp_field,flags_field);
            
        for(int i=0;i<sigma_mg_second_level.size();++i){
            r_second_level[i][0]=T_EIGEN_VECTOR::Zero(r_second_level[i][0].size());//may not need to
            {LOG::SCOPE scope("Copy to Eigen");
                EIGEN_INTERFACE<T_STRUCT,T,d>::Accumulate_To_Eigen_Array(b_second_level[i][0],tmp_field,flags_field,
                                                                         allocator,set,
                                                                         (*interface_matrix_generators_second_level)[i]->interface_map);}

            u_second_level[i][0]=T_EIGEN_VECTOR::Zero(u_second_level[i][0].size());
            sigma_mg_second_level[i]->smoother.Compute_Residual(0,u_second_level[i][0],r_second_level[i][0],b_second_level[i][0]);
            for(int j = 0;j < r_second_level[i][0].size();++j){
                r_second_level[i][0](j) = fabs(r_second_level[i][0](j));} 
            std::cout<<"Residual of sigma in SIGMA_SOLVE: "<<r_second_level[i][0].maxCoeff()<<std::endl;
            
            for(int itr=0;itr<n_sigma_mg;++itr){
                LOG::SCOPE scope("Sigma_V_Cycle");
                sigma_mg_second_level[i]->V_Cycle();}
            
            sigma_mg_second_level[i]->smoother.Compute_Residual(0,u_second_level[i][0],r_second_level[i][0],b_second_level[i][0]);
            for(int j = 0;j < r_second_level[i][0].size();++j){
                r_second_level[i][0](j) = fabs(r_second_level[i][0](j));} 
            std::cout<<"Residual of sigma after SIGMA_SOLVE: "<<r_second_level[i][0].maxCoeff()<<std::endl;           
            //LOG::cout<<sigma_mg_second_level[i]->solver_data.Prr[0]<<std::endl;
        }
        SPGrid_Computations::Clear<T_STRUCT,T,d>(allocator,set.Get_Blocks(),result_field);
        //Clear the tmp channel here. The restriction is going to read the tmp channel of interface cells which is going to be non-zero, cause MG to have a problem....
        SPGrid_Computations::Clear<T_STRUCT,T,d>(allocator,set.Get_Blocks(),tmp_field);        
        // We cleared the whole result field here....need to put back the interface values....
        {LOG::SCOPE scope("Copy from Eigen");
            EIGEN_INTERFACE<T_STRUCT,T,d>::Copy_From_Eigen_Array(u[0],result_field,flags_field,
                                                                 allocator,set,
                                                                 interface_matrix_generator.interface_map);}
        for(int i=0;i<sigma_mg_second_level.size();++i){
            {LOG::SCOPE scope("Copy from Eigen");
                EIGEN_INTERFACE<T_STRUCT,T,d>::Copy_From_Eigen_Array(u_second_level[i][0],result_field,flags_field,
                                                                     allocator,set,
                                                                     (*interface_matrix_generators_second_level)[i]->interface_map);}}
        for(int i=0;i<n_v_cycles;++i)
            SPGrid_Computations::V_Cycle_Helper<T_STRUCT,T,d>::v_cycle(result_field,tmp_field,rhs_field,flags_field,*v_cycle_topology,number_of_threads,interior_smoothing,boundary_smoothing);
    }
    void Apply_Preconditioner(T T_STRUCT::* result_field,T T_STRUCT::* rhs_field,T T_STRUCT::* tmp_field){
        //LOG::SCOPE scope("Apply_Preconditioner"); 
        u[0]=T_EIGEN_VECTOR::Zero(u[0].size());
        int n_subdomains = max_subdomain_id;
        SPGrid_Computations::Clear<T_STRUCT,T,d>(allocator,set.Get_Blocks(),result_field);
        SPGrid_Computations::Clear<T_STRUCT,T,d>(allocator,set.Get_Blocks(),tmp_field);
        //LOG::cout<<"Runing V cycles: "<<n_v_cycles<<std::endl;
        if(multilevel){
            // Well if it is multilevel, we need to do a DD here....
            Subdomain_DD(result_field,rhs_field,tmp_field);
        }else{
            for(int i=0;i<n_v_cycles;++i)
                SPGrid_Computations::V_Cycle_Helper<T_STRUCT,T,d>::v_cycle(result_field,tmp_field,rhs_field,flags_field,*v_cycle_topology,number_of_threads,interior_smoothing,boundary_smoothing);}
        
        SPGrid_Computations::Copy<T_STRUCT,T,d>(allocator,set.Get_Blocks(),rhs_field,tmp_field);
        SPGrid_Computations::Laplace_ri_Subdomain<T_STRUCT,T,unsigned,d>(allocator,set.Get_Blocks(),result_field,tmp_field,flags_field,subdomain_size);
        //LOG::cout<<"Subdomain_size: "<<subdomain_size<<std::endl;
        b[0]=T_EIGEN_VECTOR::Zero(b[0].size());
        r[0]=T_EIGEN_VECTOR::Zero(r[0].size());//may not need to
        {LOG::SCOPE scope("Copy to Eigen");
            EIGEN_INTERFACE<T_STRUCT,T,d>::Copy_To_Eigen_Array(b[0],tmp_field,flags_field,
                                                               allocator,set,
                                                               interface_matrix_generator.interface_map);}

        u[0]=T_EIGEN_VECTOR::Zero(u[0].size());
        // sigma_mg.smoother.Compute_Residual(0,u[0],r[0],b[0]);
        // for(int j = 0;j < r[0].size();++j){
        //     r[0](j) = fabs(r[0](j));} 
        //std::cout<<"Residual of sigma in SIGMA_SOLVE: "<<r[0].maxCoeff()<<std::endl;
        for(int i=0;i<n_sigma_mg;++i){
            LOG::SCOPE scope("Sigma_V_Cycle");
            sigma_mg.V_Cycle();
            // sigma_mg.smoother.Compute_Residual(0,u[0],r[0],b[0]);
            // for(int j = 0;j < r[0].size();++j){
            //     r[0](j) = fabs(r[0](j));}
            // std::cout<<"Residual of sigma in SIGMA_SOLVE: "<<r[0].maxCoeff()<<std::endl;
        }        
        SPGrid_Computations::Clear<T_STRUCT,T,d>(allocator,set.Get_Blocks(),result_field);
        //Clear the tmp channel here. The restriction is going to read the tmp channel of interface cells which is going to be non-zero, cause MG to have a problem....
        SPGrid_Computations::Clear<T_STRUCT,T,d>(allocator,set.Get_Blocks(),tmp_field);        
        {LOG::SCOPE scope("Copy from Eigen");
            EIGEN_INTERFACE<T_STRUCT,T,d>::Copy_From_Eigen_Array(u[0],result_field,flags_field,
                                                                 allocator,set,
                                                                 interface_matrix_generator.interface_map);}

        if(multilevel){
            // Well if it is multilevel, we need to do a DD here....
            Subdomain_DD(result_field,rhs_field,tmp_field);
        }else{
            for(int i=0;i<n_v_cycles;++i)
                SPGrid_Computations::V_Cycle_Helper<T_STRUCT,T,d>::v_cycle(result_field,tmp_field,rhs_field,flags_field,*v_cycle_topology,number_of_threads,interior_smoothing,boundary_smoothing);}
        
    }
};
}
#endif
