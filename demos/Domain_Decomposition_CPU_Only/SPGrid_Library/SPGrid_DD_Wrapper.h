#ifndef __SPGRID_DD_WRAPPER_H__
#define __SPGRID_DD_WRAPPER_H__
#include <Common_Tools/Grids_Uniform_PDE_Linear/STENCIL_ITERATOR.h>
#include "SPGrid_V_Cycle_Topology.h"
#include "SPGrid_V_Cycle_Helper.h"
#include "SPGrid_Laplace_Legacy.h"
#include <SPGrid/Tools/SPGrid_Copy.h>
#include <SPGrid/Tools/SPGrid_Masked_Set.h>
#include <SPGrid/Tools/SPGrid_Saxpby.h>
#include <Eigen/Sparse>
#include "../Common_Library/Sample_Domains.h"
#include "../Common_Library/Stride_Generator.h"
#include "../Common_Library/Write_Output.h"
#include "../Common_Library/Master_Cartographer.h"
#include "../Common_Library/Galerkin_Operators_Assembler.h"
#include "../Common_Library/Eigen_Interface.h"
#include "../Common_Library/Galerkin_Process_Helper.h"
#include "../Common_Library/Sigma_Multigrid_Solver.h"


using namespace PhysBAM;
using namespace SPGrid;
template<typename T,typename T_STRUCT,int d,typename INDEX>
class SPGrid_DD_Wrapper{
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
    unsigned T_STRUCT::* flags_field;
    T T_STRUCT::* u_v_field;
    T T_STRUCT::* b_v_field;
    T T_STRUCT::* r_v_field;
    SPG_Allocator& allocator;
    SPG_Set_Type& set;
    SPGrid_Computations::V_Cycle_Topology<T_STRUCT,T,d>& v_cycle_topology;
    int n_v_cycles;

    Master_Cartographer<d,INDEX> master_cartographer;
    Galerkin_Operators_Assembler<T,d,INDEX> operators;
    Interface_Solver_Data<T,d,INDEX> matrices;
    std::vector<T_EIGEN_VECTOR> r;
    std::vector<T_EIGEN_VECTOR> u;
    std::vector<T_EIGEN_VECTOR> b;
    Sigma_Multigrid_Solver<T,d,INDEX> sigma_mg;
    int levels_sigma;
    int n_smoothing_sigma;
    int n_sigma_mg;
    int n_interior_smoothing;
    int n_exterior_smoothing;
    //The following ones are temps used in the preconditioner
    ARRAY<T_EIGEN_VECTOR,T_INDEX> vi;
    SPGrid_DD_Wrapper(T_INDEX subdomain_size,
                      SPG_Allocator& allocator_input,SPG_Set_Type& set_input,std_array<int,d> size,
                      unsigned T_STRUCT::* flags_field_input,T T_STRUCT::* u_v_field_input,T T_STRUCT::* b_v_field_input,T T_STRUCT::* r_v_field_input,
                      SPGrid_Computations::V_Cycle_Topology<T_STRUCT,T,d>& v_cycle_topology_input,int levels_sigma_input,
                      int n_v_cycles_input,int n_smoothing_sigma_input,int n_sigma_mg_input,int n_interior_smoothing_input,int n_exterior_smoothing_input)
        //NOTE spgrid size is 2^n+1, we need to minus one to get the power of two cartographer wanted
        :master_cartographer(size.template Cast<T_INDEX>()-1,subdomain_size),
         allocator(allocator_input),set(set_input),flags_field(flags_field_input),
         v_cycle_topology(v_cycle_topology_input),
         u_v_field(u_v_field_input),b_v_field(b_v_field_input),r_v_field(r_v_field_input),
         n_v_cycles(n_v_cycles_input),n_interior_smoothing(n_interior_smoothing_input),n_exterior_smoothing(n_exterior_smoothing_input),
         n_sigma_mg(n_sigma_mg_input),
         operators(master_cartographer),levels_sigma(levels_sigma_input),
         sigma_mg(matrices,operators.prolongation_matrix_interface,u,r,b,n_sigma_mg_input)
    {}
    ~SPGrid_DD_Wrapper(){}
    void Initialize(){
        master_cartographer.Chart(levels_sigma);
        operators.Prebake_Operators();
        Galerkin_Process_Helper<T,d,INDEX>::Resize_Interface_Solver_Data(matrices,operators);
        Galerkin_Process_Helper<T,d,INDEX>::Bake_Interface_Distribution_Matrices(matrices,operators);
        Galerkin_Process_Helper_SPGrid<T_STRUCT,T,d,INDEX>::Galerkin_Process_Full_Laplace_Fast(matrices,operators);
        Galerkin_Process_Helper_SPGrid<T_STRUCT,T,d,INDEX>::Galerkin_Process_Fast(matrices,allocator,set.Get_Blocks(),flags_field,operators);
        //LOG::cout<<matrices.Aii(1)(T_INDEX()+1).rows()<<"x"<<matrices.Aii(1)(T_INDEX()+1).cols()<<std::endl;
        //LOG::cout << matrices.Air(1)(T_INDEX())*matrices.Di(1)(T_INDEX())<<std::endl;
        //LOG::cout<<matrices.Arr(0).rows()<<"x"<<matrices.Arr(0).cols()<<std::endl;
        //LOG::cout << matrices.Arr(1)<<std::endl;
        //exit(1);
        sigma_mg.Symbolic_Initialize();
        sigma_mg.Numerical_Initialize();
    }
    void Update(){
        Galerkin_Process_Helper_SPGrid<T_STRUCT,T,d,INDEX>::Galerkin_Process(matrices,allocator,set.Get_Blocks(),flags_field,operators);
        sigma_mg.Numerical_Initialize();
    }
    ///////////////////////////////////////////////////////////////////////////////////////////
    // Function Apply_DD
    ///////////////////////////////////////////////////////////////////////////////////////////
    void Apply_DD(T T_STRUCT::* x_field,T T_STRUCT::* rhs_field,int number_of_threads = 0){
        LOG::cout<<"number_of_threads: "<<number_of_threads<<std::endl;
        //OK, here we have three channels to work with, please note that here x_field and rhs_field are cg_field not the three channels that are reserved for v_cycle
        LOG::SCOPE scope("Applying DD Precondition");
        // Eigen_Wrapper_SPGrid<T_STRUCT,T,d,INDEX>::Copy_To_Eigen_Array(b[0],u_v_field,flags_field,allocator,set,
        //                                                               master_cartographer.index_maps_1d_to_nd_interface_global[0]);        
        // Eigen_Wrapper_SPGrid<T_STRUCT,T,d,INDEX>::Copy_From_Eigen_Array(b[0],u_v_field,flags_field,allocator,set,
        //                                                               master_cartographer.index_maps_1d_to_nd_interface_global[0]);        
        // SPGrid_Computations::Clear<T_STRUCT,T,d>(allocator,set.Get_Blocks(),x_field);
        // SPGrid_Computations::V_Cycle_Helper<T_STRUCT,T,d>::v_cycle(x_field,r_v_field,rhs_field,flags_field,v_cycle_topology,
        //                                                            number_of_threads);
        SPGrid_Computations::Clear<T_STRUCT,T,d>(allocator,set.Get_Blocks(),x_field);
        SPGrid_Computations::Clear<T_STRUCT,T,d>(allocator,set.Get_Blocks(),u_v_field);
        SPGrid_Computations::Clear<T_STRUCT,T,d>(allocator,set.Get_Blocks(),r_v_field);

        //interface smoothing
        for(int i=0;i<3;++i){
            if(number_of_threads)
                SPGrid_Computations::Threading_Helper<T_STRUCT,d>(allocator,set.Get_Blocks()).Run_Parallel(SPGrid_Computations::Residual<T_STRUCT,T,unsigned,d>(rhs_field,u_v_field,r_v_field,flags_field),number_of_threads);
            else
                SPGrid_Computations::Residual<T_STRUCT,T,unsigned,d>(allocator,set.Get_Blocks(),rhs_field,u_v_field,r_v_field,flags_field);            
            if(number_of_threads)
                SPGrid_Computations::Threading_Helper<T_STRUCT,d>(allocator,set.Get_Blocks()).Run_Parallel(SPGrid_Computations::Correction_Interface<T_STRUCT,T,unsigned,d>(u_v_field,r_v_field,flags_field),number_of_threads);
            else
                SPGrid_Computations::Correction_Interface<T_STRUCT,T,unsigned,d>(allocator,set.Get_Blocks(),u_v_field,r_v_field,flags_field);}

        for(int i=0;i<n_v_cycles;++i){
            LOG::SCOPE scope("v_cycle");
            SPGrid_Computations::V_Cycle_Helper<T_STRUCT,T,d>::v_cycle(u_v_field,r_v_field,rhs_field,flags_field,v_cycle_topology,number_of_threads,
                                                                       n_interior_smoothing,n_exterior_smoothing);}

        {LOG::SCOPE scope("Air");
        if(number_of_threads)
            SPGrid_Computations::Threading_Helper<T_STRUCT,d>(allocator,set.Get_Blocks()).Run_Parallel(SPGrid_Computations::Laplace_ir<T_STRUCT,T,unsigned,d>(u_v_field,r_v_field,flags_field),number_of_threads);
        else
            SPGrid_Computations::Laplace_ir<T_STRUCT,T,unsigned,d>(allocator,set.Get_Blocks(),u_v_field,r_v_field,flags_field);}
        
        b[0]=T_EIGEN_VECTOR::Zero(b[0].size());
        r[0]=T_EIGEN_VECTOR::Zero(r[0].size());
        {LOG::SCOPE scope("Copy to Eigen");
        Eigen_Wrapper_SPGrid<T_STRUCT,T,d,INDEX>::Copy_To_Eigen_Array(b[0],r_v_field,flags_field,allocator,set,
                                                                      master_cartographer.index_maps_1d_to_nd_interface_global[0]);        
        Eigen_Wrapper_SPGrid<T_STRUCT,T,d,INDEX>::Copy_To_Eigen_Array(r[0],rhs_field,flags_field,allocator,set,
                                                                      master_cartographer.index_maps_1d_to_nd_interface_global[0]);}

        b[0]=r[0]-b[0];
        u[0]=T_EIGEN_VECTOR::Zero(u[0].size());
        for(int i=0;i<n_sigma_mg;++i){LOG::SCOPE scope("sigma_mg V_Cycle");sigma_mg.V_Cycle();}

        SPGrid_Computations::Clear<T_STRUCT,T,d>(allocator,set.Get_Blocks(),b_v_field);

        u[0]=-u[0];
        SPGrid_Computations::Clear<T_STRUCT,T,d>(allocator,set.Get_Blocks(),r_v_field);
        Eigen_Wrapper_SPGrid<T_STRUCT,T,d,INDEX>::Copy_From_Eigen_Array(u[0],r_v_field,flags_field,allocator,set,
                                                                        master_cartographer.index_maps_1d_to_nd_interface_global[0]);
        
        {LOG::SCOPE scope("Air");
        if(number_of_threads)
            SPGrid_Computations::Threading_Helper<T_STRUCT,d>(allocator,set.Get_Blocks()).Run_Parallel(SPGrid_Computations::Laplace_ri<T_STRUCT,T,unsigned,d>(r_v_field,b_v_field,flags_field),number_of_threads);
        else
            SPGrid_Computations::Laplace_ri<T_STRUCT,T,unsigned,d>(allocator,set.Get_Blocks(),r_v_field,b_v_field,flags_field);}
        
        SPGrid_Computations::Clear<T_STRUCT,T,d>(allocator,set.Get_Blocks(),r_v_field);
        for(int i=0;i<n_v_cycles;++i){
            LOG::SCOPE scope("v_cycle");
            SPGrid_Computations::V_Cycle_Helper<T_STRUCT,T,d>::v_cycle(x_field,r_v_field,b_v_field,flags_field,v_cycle_topology,number_of_threads,
                                                                       n_interior_smoothing,n_exterior_smoothing);}
        
        SPGrid_Computations::Masked_Saxpby<T_STRUCT,T,unsigned,d>(allocator,set.Get_Blocks(),u_v_field,x_field,T(1.),T(1.),x_field,flags_field,SPGrid_Solver_Cell_Type_Active);

        u[0]=-u[0];
        {LOG::SCOPE scope("Copy from array");
        Eigen_Wrapper_SPGrid<T_STRUCT,T,d,INDEX>::Copy_From_Eigen_Array(u[0],x_field,flags_field,allocator,set,
                                                                        master_cartographer.index_maps_1d_to_nd_interface_global[0]);}

        for(int i=0;i<3;++i){
            if(number_of_threads)
                SPGrid_Computations::Threading_Helper<T_STRUCT,d>(allocator,set.Get_Blocks()).Run_Parallel(SPGrid_Computations::Residual<T_STRUCT,T,unsigned,d>(rhs_field,x_field,r_v_field,flags_field),number_of_threads);
            else
                SPGrid_Computations::Residual<T_STRUCT,T,unsigned,d>(allocator,set.Get_Blocks(),rhs_field,x_field,r_v_field,flags_field);
            if(number_of_threads)
                SPGrid_Computations::Threading_Helper<T_STRUCT,d>(allocator,set.Get_Blocks()).Run_Parallel(SPGrid_Computations::Correction_Interface<T_STRUCT,T,unsigned,d>(x_field,r_v_field,flags_field),number_of_threads);
            else
                SPGrid_Computations::Correction_Interface<T_STRUCT,T,unsigned,d>(allocator,set.Get_Blocks(),x_field,r_v_field,flags_field);}

    }
};
#endif
