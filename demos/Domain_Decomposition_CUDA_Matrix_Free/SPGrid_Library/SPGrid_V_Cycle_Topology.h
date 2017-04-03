//#####################################################################
//Removed For Anonymity: Copyright Authors of Submission pap173s1
//Removed For Anonymity: Copyright Authors of Submission pap173s1
//#####################################################################
// Subroutine SPGrid_Computations::V_Cycle_Toplogy
//#####################################################################
#ifndef __SPGrid_V_Cycle_Toplogy_h__
#define __SPGrid_V_Cycle_Toplogy_h__

#include <SPGrid/Core/SPGrid_Allocator.h>
#include <SPGrid_Fluids/Grids/GRID_TOPOLOGY_HELPER.h>
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <Threading_Tools/PTHREAD_QUEUE.h>
#include <SPGrid/Tools/SPGrid_Threading_Helper.h>
#include <SPGrid/Tools/SPGrid_Clear.h>
#include <SPGrid/Tools/SPGrid_Block_Iterator.h>
#include "SPGrid_Flag_Helper.h"
#include "SPGrid_Flag_Boundary_Helper.h"
#include <vector>
namespace SPGrid_Computations{

using namespace SPGrid;
using namespace PhysBAM;
template<class T_STRUCT,class T,int d>
class V_Cycle_Topology
{
    typedef SPGrid_Allocator<T_STRUCT,d> SPG_Allocator;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<const unsigned>::type Const_Flag_Array_Type;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<unsigned>::type SPG_Flags_Array_Type;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<const T>::type Const_Data_Array_Type;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<T>::type Data_Array_Type;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<const T>::mask T_MASK;
    typedef SPGrid_Set<SPG_Flags_Array_Type> SPG_Set_Type;
    typedef VECTOR<int,d> T_INDEX;
    static void Topology_Downsample(SPG_Allocator& fine_allocator, SPG_Set_Type&fine_set,SPG_Allocator& coarse_allocator, SPG_Set_Type&coarse_set, unsigned T_STRUCT::* flag_field){ 
        typedef VECTOR<int,d> T_INDEX;
        SPG_Flags_Array_Type coarse_flags = coarse_allocator.Get_Array(flag_field);
        SPG_Flags_Array_Type fine_flags = fine_allocator.Get_Array(flag_field);
        for(SPGrid_Block_Iterator<T_MASK> fine_iterator(fine_set.Get_Blocks());fine_iterator.Valid();fine_iterator.Next_Block()){
            unsigned long fine_offset=fine_iterator.Offset();
            T_INDEX base_index=fine_iterator.Index().template Cast<T_INDEX>();
            for(RANGE_ITERATOR<d> iterator(RANGE<T_INDEX>(base_index,base_index+fine_allocator.Block_Size().template Cast<T_INDEX>()-1));
                iterator.Valid();
                iterator.Next(),fine_offset+=sizeof(unsigned)){
                T_INDEX fine_index=iterator.Index();
                T_INDEX coarse_index;for(int v=1;v<=d;v++) coarse_index(v)=(fine_index(v)+1)/2;
                unsigned long coarse_offset=T_MASK::DownsampleOffset(fine_offset);
                PHYSBAM_ASSERT(coarse_offset==T_MASK::Linear_Offset(std_array<int,d>(coarse_index)));
                unsigned fine_flag=fine_flags(fine_offset);
                //if the flag is Active on the fine, make sure mark all the coarse cells that it is going to touch during the prolongation
                if(fine_flag & SPGrid_Solver_Cell_Type_Active){
                    T_INDEX coarse_base_index;for(int v=1;v<=d;v++) coarse_base_index(v)=fine_index(v)/2;
                    for(RANGE_ITERATOR<d> coarse_iterator(RANGE<T_INDEX>(coarse_base_index,coarse_base_index+1));coarse_iterator.Valid();coarse_iterator.Next()){
                        coarse_set.MarkPageActive(T_MASK::Linear_Offset(std_array<int,d>(coarse_iterator.Index())));
                    }
                }
                if(!coarse_set.Is_Set(coarse_offset,SPGrid_Solver_Cell_Type_Dirichlet|SPGrid_Solver_Cell_Type_Active|SPGrid_Solver_Cell_Type_Interface)){
                    if(fine_flag & SPGrid_Solver_Cell_Type_Active){
                        coarse_set.Mask(coarse_offset,SPGrid_Solver_Cell_Type_Active);
                    }else if(fine_flag & (SPGrid_Solver_Cell_Type_Dirichlet|SPGrid_Solver_Cell_Type_Interface)){
                        coarse_set.Mask(coarse_offset,fine_flag & (SPGrid_Solver_Cell_Type_Dirichlet));
                    }
                    continue;}
                if(fine_flag & (SPGrid_Solver_Cell_Type_Dirichlet | SPGrid_Solver_Cell_Type_Interface))
                    coarse_flags(coarse_offset)=SPGrid_Solver_Cell_Type_Dirichlet;               
            }
        }
        coarse_set.Refresh_Block_Offsets();
    }
    static void Topology_Up_Sweep(SPG_Allocator& fine_allocator, SPG_Set_Type&fine_set,SPG_Allocator& coarse_allocator, SPG_Set_Type&coarse_set, unsigned T_STRUCT::* flag_field){
        typedef VECTOR<int,d> T_INDEX;
        SPG_Flags_Array_Type coarse_flags = coarse_allocator.Get_Array(flag_field);
        SPG_Flags_Array_Type fine_flags = fine_allocator.Get_Array(flag_field);
        for(SPGrid_Block_Iterator<T_MASK> coarse_iterator(coarse_set.Get_Blocks());coarse_iterator.Valid();coarse_iterator.Next_Block()){
            unsigned long coarse_offset=coarse_iterator.Offset();
            T_INDEX base_index=coarse_iterator.Index().template Cast<T_INDEX>();
            for(RANGE_ITERATOR<d> iterator(RANGE<T_INDEX>(base_index,base_index+coarse_allocator.Block_Size().template Cast<T_INDEX>()-1));
                iterator.Valid();
                iterator.Next(),coarse_offset+=sizeof(unsigned)){
                T_INDEX coarse_index=iterator.Index();
                T_INDEX fine_index;for(int v=1;v<=d;v++) fine_index(v)=coarse_index(v)*2-1;
                unsigned long fine_offset=T_MASK::UpsampleOffset(coarse_offset);
                PHYSBAM_ASSERT(fine_offset==T_MASK::Linear_Offset(std_array<int,d>(fine_index)));
                unsigned coarse_flag=coarse_flags(coarse_offset);
                //if the flag is Active on the fine, make sure mark all the coarse cells that it is going to touch during the prolongation
                if(coarse_flag & SPGrid_Solver_Cell_Type_Active){
                    T_INDEX fine_base_index;for(int v=1;v<=d;v++) fine_base_index(v)=coarse_index(v)*2-2;
                    for(RANGE_ITERATOR<d> fine_iterator(RANGE<T_INDEX>(fine_base_index,fine_base_index+3));fine_iterator.Valid();fine_iterator.Next()){
                        if(!RANGE<T_INDEX>(fine_base_index+1,fine_base_index+2).Lazy_Inside(fine_iterator.Index()))
                            fine_set.MarkPageActive(T_MASK::Linear_Offset(std_array<int,d>(fine_iterator.Index())));
                    }
                }
            }            
        }
        fine_set.Refresh_Block_Offsets();
    }
    std::vector<SPG_Allocator*> allocator_hierarchy;
    std::vector<SPG_Set_Type*> set_hierarchy;
    std::vector<std_array<int,d> > size_hierarchy;
    unsigned T_STRUCT::* flags_field; 
    const bool own_top;
public:
    const int levels;
    V_Cycle_Topology(
            SPG_Allocator& SPGrid_top_allocator_ptr,
            SPG_Set_Type& SPGrid_top_set_ptr,
            std_array<int,d> size_top, int levels_input,unsigned T_STRUCT::* flags_field_input)
        :levels(levels_input),flags_field(flags_field_input),own_top(false)
    {
        //Static_Assert(d == 3);
        allocator_hierarchy.resize(levels);
        set_hierarchy.resize(levels);
        size_hierarchy.resize(levels);
        allocator_hierarchy[0] = &SPGrid_top_allocator_ptr;
        set_hierarchy[0] = &SPGrid_top_set_ptr;
        size_hierarchy[0] = size_top; 
        for(int level = 0; level < levels - 1;++level){
            T_INDEX tmp_index = size_hierarchy[level].template Cast<T_INDEX>();
            size_hierarchy[level + 1] = std_array<int,d>(tmp_index/2 + 1);
            allocator_hierarchy[level + 1] = new SPG_Allocator(size_hierarchy[level + 1]);
            set_hierarchy[level + 1]  = new SPG_Set_Type(allocator_hierarchy[level + 1] -> Get_Array(flags_field));
        }
    }
    V_Cycle_Topology(std_array<int,d> size_top, int levels_input,unsigned T_STRUCT::* flags_field_input)
        :levels(levels_input),flags_field(flags_field_input),own_top(true)
    {
        //Static_Assert(d == 3);
        allocator_hierarchy.resize(levels);
        set_hierarchy.resize(levels);
        size_hierarchy.resize(levels);
        size_hierarchy[0] = size_top;
        allocator_hierarchy[0] = new SPG_Allocator(size_hierarchy[0]);
        set_hierarchy[0]  = new SPG_Set_Type(allocator_hierarchy[0] -> Get_Array(flags_field)); 
        for(int level = 0; level < levels - 1;++level){
            T_INDEX tmp_index = size_hierarchy[level].template Cast<T_INDEX>();
            size_hierarchy[level + 1] = std_array<int,d>(tmp_index/2 + 1);
            allocator_hierarchy[level + 1] = new SPG_Allocator(size_hierarchy[level + 1]);
            set_hierarchy[level + 1]  = new SPG_Set_Type(allocator_hierarchy[level + 1] -> Get_Array(flags_field));
        }
    }
    ~V_Cycle_Topology(){
        if(own_top) {delete allocator_hierarchy[0];delete set_hierarchy[0];}
        for(int level = 1; level < levels;++level){
            delete allocator_hierarchy[level];
            delete set_hierarchy[level];
        }
    }
    void Initialize(int number_of_threads=0){
        for(int level=0;level<levels-1;++level){
            SPG_Allocator& fine_allocator = *(allocator_hierarchy[level]);
            SPG_Allocator& coarse_allocator = *(allocator_hierarchy[level + 1]);
            SPG_Set_Type& fine_set = *(set_hierarchy[level]);
            SPG_Set_Type& coarse_set = *(set_hierarchy[level + 1]);
            {LOG::SCOPE scope("Topology_Downsample");
            Topology_Downsample(fine_allocator,fine_set,coarse_allocator,coarse_set,flags_field);}}
    }
    void Flaging(int number_of_threads=0){
        for(int level = 0; level < levels;++level){
            SPG_Allocator& allocator = *(allocator_hierarchy[level]);
            SPG_Set_Type& set = *(set_hierarchy[level]);
            if(number_of_threads)
                SPGrid_Computations::Threading_Helper<T_STRUCT,d>(allocator,set.Get_Blocks()).Run_Parallel(SPGrid_Computations::Flaging<T_STRUCT,d>(flags_field,set),number_of_threads);
            else
                SPGrid_Computations::Flaging<T_STRUCT,d>(allocator,set.Get_Blocks(),flags_field,set);

            if(number_of_threads)
                SPGrid_Computations::Threading_Helper<T_STRUCT,d>(allocator,set.Get_Blocks()).Run_Parallel(SPGrid_Computations::Flaging_Boundary<T_STRUCT,d>(flags_field,set),number_of_threads);
            else
                SPGrid_Computations::Flaging_Boundary<T_STRUCT,d>(allocator,set.Get_Blocks(),flags_field,set);
        }
    }
    SPG_Allocator& Get_Allocator(int level){
        return *(allocator_hierarchy[level]);
    }
    const SPG_Allocator& Get_Allocator(int level)const{
        return *(allocator_hierarchy[level]);
    }
    SPG_Set_Type& Get_Set(int level){
        return *(set_hierarchy[level]);
    }
    const SPG_Set_Type& Get_Set(int level)const{
        return *(set_hierarchy[level]);
    }
    std_array<unsigned,d> Get_Size(int level)const{
        return size_hierarchy[level];
    } 
    int Levels()const{
        return levels;
    }
};
}
#endif
