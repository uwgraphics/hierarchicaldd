//#####################################################################
//Removed For Anonymity: Copyright Authors of Submission pap173s1
//Removed For Anonymity: Copyright Authors of Submission pap173s1
//#####################################################################
// Subroutine SPGrid_Computations::V_Cycle
//#####################################################################
#ifndef __SPGrid_V_Cycle_Helper_h__
#define __SPGrid_V_Cycle_Helper_h__

#include <SPGrid/Core/SPGrid_Allocator.h>
#include <SPGrid_Fluids/Grids/GRID_TOPOLOGY_HELPER.h>
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <Threading_Tools/PTHREAD_QUEUE.h>
#include <SPGrid/Tools/SPGrid_Threading_Helper.h>
#include <SPGrid/Tools/SPGrid_Block_Iterator.h>
#include <SPGrid/Tools/SPGrid_Clear.h>
#include "SPGrid_Residual_Legacy.h"
#include "SPGrid_Correction_Legacy.h"
#include "SPGrid_Restriction_Legacy.h"
#include "SPGrid_Prolongation_Legacy.h"
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Read_Write/Grids_Uniform/READ_WRITE_GRID.h>
#include <PhysBAM_Tools/Read_Write/Grids_Uniform_Arrays/READ_WRITE_ARRAYS.h>
#include "SPGrid_V_Cycle_Topology.h"
#include <vector>
namespace SPGrid_Computations{

using namespace SPGrid;
using namespace PhysBAM;
template<class T_STRUCT, class T,int d>
class Norm_Helper{
    typedef VECTOR<int,d> T_INDEX;
    typedef SPGrid_Allocator<T_STRUCT,d> SPG_Allocator;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<const unsigned>::type SPG_Flags_Array_Type;
    typedef SPGrid_Set<SPG_Flags_Array_Type> SPG_Set_Type;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<const T>::type Const_Data_Array_Type;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<const T>::mask T_MASK;

public:
    static T L1_Norm(const SPG_Allocator& allocator, const std::pair<const unsigned long*,unsigned>& blocks,T T_STRUCT::* r_field){
        T result = 0.f;
        Const_Data_Array_Type r=allocator.Get_Const_Array(r_field);
        for(SPGrid_Block_Iterator<T_MASK> iterator(blocks);iterator.Valid();iterator.Next_Block()){
            unsigned long offset=iterator.Offset();
            T_INDEX base_index=T_INDEX(iterator.Index()(0),iterator.Index()(1),iterator.Index()(2));
            for(RANGE_ITERATOR<d> iterator(RANGE<T_INDEX>(base_index,base_index+T_INDEX(allocator.Block_Size()(0),allocator.Block_Size()(1),allocator.Block_Size()(2))-1));
                iterator.Valid();
                iterator.Next(),offset+=sizeof(T)){
                if(fabs(r(offset)) > result) result = fabs(r(offset));
            }
        }
        return result;
    }
};

template<class T_STRUCT,class T_DATA,int d>
class V_Cycle_Helper
{
    typedef T_DATA T;
    typedef VECTOR<int,d> T_INDEX;
    typedef VECTOR<T,d> TV;
    typedef SPGrid_Allocator<T_STRUCT,d> SPG_Allocator;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<unsigned>::type SPG_Flags_Array_Type;
    typedef SPGrid_Set<SPG_Flags_Array_Type> SPG_Set_Type;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<const T>::type Const_Data_Array_Type;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<const T>::mask T_MASK;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<T_DATA>::type Data_Array_Type;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<const unsigned>::type Const_Flag_Array_Type;
public:
    static inline void boundary_smoothing(
                                          T_DATA T_STRUCT::* u_field,T_DATA T_STRUCT::* r_field,T_DATA T_STRUCT::* b_field,unsigned T_STRUCT::* flags_field,
                                          SPG_Allocator& allocator, SPG_Set_Type& set,int number_of_threads = 0){
        if(number_of_threads)
            SPGrid_Computations::Threading_Helper<T_STRUCT,d>(allocator,set.Get_Blocks()).Run_Parallel(SPGrid_Computations::Residual<T_STRUCT,T,unsigned,d>(b_field,u_field,r_field,flags_field),number_of_threads);
        else
            SPGrid_Computations::Residual<T_STRUCT,T,unsigned,d>(allocator,set.Get_Blocks(),b_field,u_field,r_field,flags_field);

        if(number_of_threads)
            SPGrid_Computations::Threading_Helper<T_STRUCT,d>(allocator,set.Get_Blocks()).Run_Parallel(SPGrid_Computations::Correction<T_STRUCT,T,unsigned,d>(u_field,r_field,flags_field,false,true),number_of_threads);
        else
            SPGrid_Computations::Correction<T_STRUCT,T,unsigned,d>(allocator,set.Get_Blocks(),u_field,r_field,flags_field,false,true);
    }
    static inline void interior_smoothing(                                          
                                          T_DATA T_STRUCT::* u_field,T_DATA T_STRUCT::* r_field,T_DATA T_STRUCT::* b_field,unsigned T_STRUCT::* flags_field,
                                          SPG_Allocator& allocator, SPG_Set_Type& set,int number_of_threads = 0){
        if(number_of_threads)
            SPGrid_Computations::Threading_Helper<T_STRUCT,d>(allocator,set.Get_Blocks()).Run_Parallel(SPGrid_Computations::Residual<T_STRUCT,T,unsigned,d>(b_field,u_field,r_field,flags_field),number_of_threads);
        else
            SPGrid_Computations::Residual<T_STRUCT,T,unsigned,d>(allocator,set.Get_Blocks(),b_field,u_field,r_field,flags_field);

        if(number_of_threads)
            SPGrid_Computations::Threading_Helper<T_STRUCT,d>(allocator,set.Get_Blocks()).Run_Parallel(SPGrid_Computations::Correction<T_STRUCT,T,unsigned,d>(u_field,r_field,flags_field,true,false),number_of_threads);
        else
            SPGrid_Computations::Correction<T_STRUCT,T,unsigned,d>(allocator,set.Get_Blocks(),u_field,r_field,flags_field,true,false);
    }
    static inline void bottom_smoothing(T_DATA T_STRUCT::* u_field,T_DATA T_STRUCT::* r_field,T_DATA T_STRUCT::* b_field,unsigned T_STRUCT::* flags_field,
                                        SPG_Allocator& allocator, SPG_Set_Type& set,int number_of_threads = 0){
        if(number_of_threads)
            SPGrid_Computations::Threading_Helper<T_STRUCT,d>(allocator,set.Get_Blocks()).Run_Parallel(SPGrid_Computations::Residual<T_STRUCT,T,unsigned,d>(b_field,u_field,r_field,flags_field),number_of_threads);
        else
            SPGrid_Computations::Residual<T_STRUCT,T,unsigned,d>(allocator,set.Get_Blocks(),b_field,u_field,r_field,flags_field);

        if(number_of_threads)
            SPGrid_Computations::Threading_Helper<T_STRUCT,d>(allocator,set.Get_Blocks()).Run_Parallel(SPGrid_Computations::Correction<T_STRUCT,T,unsigned,d>(u_field,r_field,flags_field),number_of_threads);
        else
            SPGrid_Computations::Correction<T_STRUCT,T,unsigned,d>(allocator,set.Get_Blocks(),u_field,r_field,flags_field);
    }
    static inline void global_smoothing(
                                          T_DATA T_STRUCT::* u_field,T_DATA T_STRUCT::* r_field,T_DATA T_STRUCT::* b_field,unsigned T_STRUCT::* flags_field,
                                          SPG_Allocator& allocator, SPG_Set_Type& set,int number_of_threads = 0){
        if(number_of_threads)
            SPGrid_Computations::Threading_Helper<T_STRUCT,d>(allocator,set.Get_Blocks()).Run_Parallel(SPGrid_Computations::Residual_Global<T_STRUCT,T,unsigned,d>(b_field,u_field,r_field,flags_field),number_of_threads);
        else
            SPGrid_Computations::Residual_Global<T_STRUCT,T,unsigned,d>(allocator,set.Get_Blocks(),b_field,u_field,r_field,flags_field);

        if(number_of_threads)
            SPGrid_Computations::Threading_Helper<T_STRUCT,d>(allocator,set.Get_Blocks()).Run_Parallel(SPGrid_Computations::Correction_Global<T_STRUCT,T,unsigned,d>(u_field,r_field,flags_field),number_of_threads);
        else
            SPGrid_Computations::Correction_Global<T_STRUCT,T,unsigned,d>(allocator,set.Get_Blocks(),u_field,r_field,flags_field);
    }
    static inline void compute_residual(T_DATA T_STRUCT::* u_field,T_DATA T_STRUCT::* r_field,T_DATA T_STRUCT::* b_field,unsigned T_STRUCT::* flags_field,
                                        SPG_Allocator& allocator, SPG_Set_Type& set,int number_of_threads = 0){
        if(number_of_threads)
            SPGrid_Computations::Threading_Helper<T_STRUCT,d>(allocator,set.Get_Blocks()).Run_Parallel(SPGrid_Computations::Residual<T_STRUCT,T,unsigned,d>(b_field,u_field,r_field,flags_field),number_of_threads);
        else
            SPGrid_Computations::Residual<T_STRUCT,T,unsigned,d>(allocator,set.Get_Blocks(),b_field,u_field,r_field,flags_field);        
    }
    static inline void compute_residual_global(T_DATA T_STRUCT::* u_field,T_DATA T_STRUCT::* r_field,T_DATA T_STRUCT::* b_field,unsigned T_STRUCT::* flags_field,
                                               SPG_Allocator& allocator, SPG_Set_Type& set,int number_of_threads = 0){
        if(number_of_threads)
            SPGrid_Computations::Threading_Helper<T_STRUCT,d>(allocator,set.Get_Blocks()).Run_Parallel(SPGrid_Computations::Residual_Global<T_STRUCT,T,unsigned,d>(b_field,u_field,r_field,flags_field),number_of_threads);
        else
            SPGrid_Computations::Residual_Global<T_STRUCT,T,unsigned,d>(allocator,set.Get_Blocks(),b_field,u_field,r_field,flags_field);        
    }
    static inline void clear_u(T_DATA T_STRUCT::* u_field,
                               SPG_Allocator& allocator, SPG_Set_Type& set,int number_of_threads = 0){
        if(number_of_threads)
            SPGrid_Computations::Threading_Helper<T_STRUCT,d>(allocator,set.Get_Blocks()).Run_Parallel(SPGrid_Computations::Clear<T_STRUCT,T,d>(u_field),number_of_threads);
        else
            SPGrid_Computations::Clear<T_STRUCT,T,d>(allocator,set.Get_Blocks(),u_field);
    }
    static inline void restriction(T_DATA T_STRUCT::* r_field,T_DATA T_STRUCT::* b_field,unsigned T_STRUCT::* flags_field,
                                   SPG_Allocator& coarse_allocator,SPG_Set_Type& coarse_set,SPG_Allocator& fine_allocator,int number_of_threads = 0){
        if(number_of_threads)
            SPGrid_Computations::Threading_Helper<T_STRUCT,d>(coarse_allocator,coarse_set.Get_Blocks()).Run_Parallel(SPGrid_Computations::Restriction<T_STRUCT,T,unsigned,d>(fine_allocator,b_field,r_field,flags_field),number_of_threads);
        else
            SPGrid_Computations::Restriction<T_STRUCT,T,unsigned,d>(coarse_allocator,coarse_set.Get_Blocks(),fine_allocator,b_field,r_field,flags_field);        
    }
    static inline void prolongation(T_DATA T_STRUCT::* u_field,unsigned T_STRUCT::* flags_field,
                                    SPG_Allocator& fine_allocator,SPG_Set_Type& fine_set,SPG_Allocator& coarse_allocator,int number_of_threads = 0){
        if(number_of_threads)
            SPGrid_Computations::Threading_Helper<T_STRUCT,d>(fine_allocator,fine_set.Get_Blocks()).Run_Parallel(SPGrid_Computations::Prolongation<T_STRUCT,T,unsigned,d>(coarse_allocator,u_field,flags_field),number_of_threads);
        else
            SPGrid_Computations::Prolongation<T_STRUCT,T,unsigned,d>(fine_allocator,fine_set.Get_Blocks(),coarse_allocator,u_field,flags_field);
    }
    static void v_cycle(T_DATA T_STRUCT::* u_field,
                        T_DATA T_STRUCT::* r_field,
                        T_DATA T_STRUCT::* b_field,
                        unsigned T_STRUCT::* flags_field,
                        V_Cycle_Topology<T_STRUCT,T,d>& topology,
                        const int number_of_threads = 0,
                        int interior_itr = 1,int boundary_itr = 3){
        const int levels = topology.levels;
        for(int level = 0; level < levels - 1;++level){
            SPG_Allocator& fine_allocator=topology.Get_Allocator(level);
            SPG_Allocator& coarse_allocator=topology.Get_Allocator(level+1);
            SPG_Set_Type& fine_set=topology.Get_Set(level);
            SPG_Set_Type& coarse_set=topology.Get_Set(level+1);
            for(int i = 0; i < boundary_itr; ++i) boundary_smoothing(u_field,r_field,b_field,flags_field,fine_allocator,fine_set,number_of_threads);
            for(int i = 0; i < interior_itr; ++i) interior_smoothing(u_field,r_field,b_field,flags_field,fine_allocator,fine_set,number_of_threads);
            for(int i = 0; i < boundary_itr; ++i) boundary_smoothing(u_field,r_field,b_field,flags_field,fine_allocator,fine_set,number_of_threads);
            compute_residual(u_field,r_field,b_field,flags_field,fine_allocator, fine_set,number_of_threads);
            clear_u(u_field,coarse_allocator, coarse_set,number_of_threads);
            //std::cout<<"residual before restriction(CPU): "<<Norm_Helper<T_STRUCT,T,d>::L1_Norm(fine_allocator,fine_set.Get_Blocks(),r_field)<<std::endl;
            restriction(r_field,b_field,flags_field,coarse_allocator,coarse_set,fine_allocator,number_of_threads);
            //std::cout<<"rhs after restriction(CPU): "<<Norm_Helper<T_STRUCT,T,d>::L1_Norm(coarse_allocator,coarse_set.Get_Blocks(),b_field)<<std::endl;
        }
        compute_residual(u_field,r_field,b_field,flags_field,topology.Get_Allocator(levels-1),topology.Get_Set(levels-1),number_of_threads);
        for(int i = 0;i < 200;++i) bottom_smoothing(u_field,r_field,b_field,flags_field,topology.Get_Allocator(levels-1),topology.Get_Set(levels-1),number_of_threads);
        compute_residual(u_field,r_field,b_field,flags_field,topology.Get_Allocator(levels-1),topology.Get_Set(levels-1),number_of_threads);
        //std::cout<<"r after bottom solve(CPU): "<<Norm_Helper<T_STRUCT,T,d>::L1_Norm(topology.Get_Allocator(levels-1),topology.Get_Set(levels-1).Get_Blocks(),r_field)<<std::endl;
        for(int level = levels-2; level >= 0;--level){
            SPG_Allocator& fine_allocator=topology.Get_Allocator(level);
            SPG_Allocator& coarse_allocator=topology.Get_Allocator(level+1);
            SPG_Set_Type& fine_set=topology.Get_Set(level);
            SPG_Set_Type& coarse_set=topology.Get_Set(level+1);
            //std::cout<<"u before prolongation coarse(CPU): "<<Norm_Helper<T_STRUCT,T,d>::L1_Norm(coarse_allocator,coarse_set.Get_Blocks(),u_field)<<std::endl;
            //std::cout<<"u before prolongation fine(CPU): "<<Norm_Helper<T_STRUCT,T,d>::L1_Norm(fine_allocator,fine_set.Get_Blocks(),u_field)<<std::endl;
            prolongation(u_field,flags_field,fine_allocator,fine_set,coarse_allocator,number_of_threads);
            //std::cout<<"u after prolongation(CPU): "<<Norm_Helper<T_STRUCT,T,d>::L1_Norm(fine_allocator,fine_set.Get_Blocks(),u_field)<<std::endl;
            for(int i = 0; i < boundary_itr; ++i) boundary_smoothing(u_field,r_field,b_field,flags_field,fine_allocator,fine_set,number_of_threads);
            for(int i = 0; i < interior_itr; ++i) interior_smoothing(u_field,r_field,b_field,flags_field,fine_allocator,fine_set,number_of_threads);
            for(int i = 0; i < boundary_itr; ++i) boundary_smoothing(u_field,r_field,b_field,flags_field,fine_allocator,fine_set,number_of_threads);
        }
    }
    //////////////////////////////////////////////////////////////////////////////////
    // The Following Codes Are For Debuging Purpose
    //////////////////////////////////////////////////////////////////////////////////
};
//#####################################################################
}
#endif
