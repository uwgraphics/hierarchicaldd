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
#include <unistd.h>

#include "../Interface_Solver/INTERFACE_MATRIX_GENERATOR.h"
#include "../Interface_Solver/SIGMA_MULTIGRID_SOLVER.h"
#include "../Interface_Solver/EIGEN_INTERFACE.h"
#include "../Nodes/Overseer.h"

#include <PhysBAM_Tools/Krylov_Solvers/CONJUGATE_GRADIENT.h>
//debug 
#include "../SPGrid_Library/SPGrid_Residual_Legacy.h"
#include "../SPGrid_Library/SPGrid_V_Cycle_Helper.h"
#include "SPGrid_Linearized_Data_Copy_Hashtable_Helper.h"
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
    int max_subdomain_id;
    //Preconditioner Parameters
    int levels_sigma;
    int n_v_cycles;
    int n_smoothing_sigma;
    int n_sigma_mg;
    int number_of_threads;

    Domain_Decomposition::Overseer<T,T_STRUCT,d,T_offset_ptr> overseer;    
    //The following ones are temps used in the preconditioner
    ARRAY<T_EIGEN_VECTOR,T_INDEX> vi;

    SIGMA_MULTIGRID_SOLVER<T,d> sigma_mg;
    INTERFACE_MATRIX_GENERATOR<T,d>& interface_matrix_generator;
    std::vector<T_EIGEN_VECTOR> r;
    std::vector<T_EIGEN_VECTOR> u;
    std::vector<T_EIGEN_VECTOR> b;
    SPGrid_Domain_Decomposition_Preconditioner(SPG_Allocator& allocator_input,SPG_Set_Type& set_input,
                                               INTERFACE_MATRIX_GENERATOR<T,d>& interface_matrix_generator_input,
                                               std::vector<std::vector<T_Linearizer>*>& linearizer_hierarchy_input,
                                               T_INDEX subdomain_size,T_INDEX size,unsigned T_STRUCT::* flags_field_input,
                                               int levels_sigma_input,int n_v_cycles_input,
                                               int n_smoothing_sigma_input,int n_sigma_mg_input,int number_of_threads_in)
    //NOTE spgrid size is 2^n+1, we need to minus one to get the power of two cartographer wanted
    :allocator(allocator_input),set(set_input),
        interface_matrix_generator(interface_matrix_generator_input),
        linearizer_hierarchy(linearizer_hierarchy_input),
        flags_field(flags_field_input),
        n_v_cycles(n_v_cycles_input),
        n_sigma_mg(n_sigma_mg_input),max_subdomain_id(linearizer_hierarchy.size()),
        levels_sigma(levels_sigma_input),
        sigma_mg(matrices,u,r,b,n_smoothing_sigma_input),
        number_of_threads(number_of_threads_in)
    {}
    ~SPGrid_Domain_Decomposition_Preconditioner(){}
    void Initialize(bool MPI_enable,bool SCIF_enable,bool CUDA_enable){
        LOG::cout<<"data size: "<<set.Get_Blocks().second*(double)4096/1024./1024./1024.<<" GB"<<std::endl;
        // initialize MPI stuff
        overseer.Discover(false,false,true);
        overseer.Assign_Subdomains(max_subdomain_id);
        overseer.Initialize_Subdomains(linearizer_hierarchy);
        overseer.Send_Flags(linearizer_hierarchy,flags_field,number_of_threads);

        interface_matrix_generator.Convert_To_Eigen(matrices);
        sigma_mg.Symbolic_Initialize();
        sigma_mg.Numerical_Initialize();
        sleep(20);
    }
    void Apply_Preconditioner(T T_STRUCT::* result_field,T T_STRUCT::* rhs_field,T T_STRUCT::* tmp_field){
        for(SPGrid_Block_Iterator<T_MASK> iterator(set.Get_Blocks());iterator.Valid();iterator.Next_Block()){
            unsigned long offset=iterator.Offset();
            T_INDEX base_index=iterator.Index().template Cast<T_INDEX>();
            for(RANGE_ITERATOR<d> iterator(RANGE<T_INDEX>(base_index,base_index+allocator.Block_Size().template Cast<T_INDEX>()-1));
                iterator.Valid();
                iterator.Next(),offset+=sizeof(T)){
                if(set.array(offset)&SPGrid_Solver_Cell_Type_Interface)
                    PHYSBAM_ASSERT(interface_matrix_generator.interface_hash[0].Get_Default(iterator.Index(),-1)!=-1);}}
        

        //std::cout<<"Apply_Preconditioner!"<<std::endl;
        //LOG::SCOPE scope("Apply_Preconditioner"); 
        int n_subdomains = linearizer_hierarchy.size();
        SPGrid_Computations::Clear<T_STRUCT,T,d>(allocator,set.Get_Blocks(),result_field);
        //This set channel function trigers an send buffer. if this three variables gets distroyed....we are in huge trouble...
        overseer.Set_Channels(flags_field,result_field,rhs_field,tmp_field);

        // just need the interface blocks to be copied....so, should have done something more efficient
        // And it can be done during better times...
        {LOG::SCOPE scope("Copy");
            SPGrid_Computations::Copy<T_STRUCT,T,d>(allocator,set.Get_Blocks(),rhs_field,tmp_field);} 
        
        {LOG::SCOPE scope("Step ONE");
            overseer.Step_One(linearizer_hierarchy,allocator,flags_field,rhs_field,tmp_field,number_of_threads);}
        
        //exit(0);

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
        for(int i=0;i<1/*n_sigma_mg*/;++i){
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
