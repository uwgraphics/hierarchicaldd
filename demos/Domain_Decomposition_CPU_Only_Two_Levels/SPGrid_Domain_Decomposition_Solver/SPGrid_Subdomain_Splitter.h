//#####################################################################
//Removed For Anonymity: Copyright Authors of Submission pap173s1
//Removed For Anonymity: Copyright Authors of Submission pap173s1
//#####################################################################
// Class SPGrid_Subdomain_Splitter
//#####################################################################
#ifndef __SPGRID_SUBDOMAIN_SPLITTER_H__
#define __SPGRID_SUBDOMAIN_SPLITTER_H__

#include <Common_Tools/Math_Tools/RANGE_ITERATOR.h>
#include <Common_Tools/Math_Tools/HIERARCHICAL_RANGE_ITERATOR.h>
#include <SPGrid/Core/SPGrid_Allocator.h>
#include <SPGrid/Core/SPGrid_Set.h>
#include <SPGrid/Tools/SPGrid_Block_Iterator.h>

#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Read_Write/Grids_Uniform/READ_WRITE_GRID.h>
#include <PhysBAM_Tools/Read_Write/Grids_Uniform_Arrays/READ_WRITE_ARRAYS.h>
#include <PhysBAM_Geometry/Geometry_Particles/REGISTER_GEOMETRY_READ_WRITE.h>
#include <PhysBAM_Geometry/Solids_Geometry/DEFORMABLE_GEOMETRY_COLLECTION.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/FREE_PARTICLES.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/SEGMENTED_CURVE.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <PhysBAM_Geometry/Read_Write/Geometry/READ_WRITE_TRIANGULATED_SURFACE.h>
#include <PhysBAM_Tools/Data_Structures/HASHTABLE.h>
#include <PhysBAM_Tools/Data_Structures/HASHTABLE_ITERATOR.h>

#include "SPGrid_Block_Pair_Weaver.h"
#include "SPGrid_Domain_Divider_Helper.h"
#include "SPGrid_V_Cycle_Topology.h"
#include "SPGrid_Subdomain_Rasterizer.h"
#include "../Interface_Solver/INTERFACE_MATRIX_GENERATOR.h"
#include "../Interface_Solver/ADAPTIVE_SUBDOMAIN_POISSON_DATA.h"

namespace SPGrid{
template<typename T,typename T_STRUCT,int d,typename T_offset_ptr>
class SPGrid_Subdomain_Splitter{
    typedef ADAPTIVE_SUBDOMAIN_POISSON_DATA<T> T_STRUCT_ADAPTATION;
    typedef SPGrid_Allocator<T_STRUCT,d> SPG_Allocator;
    typedef typename SPG_Allocator::template Array<unsigned>::type SPG_Flags_Array_Type;
    typedef typename SPG_Allocator::template Array<T>::type SPG_Data_Array_Type;
    typedef typename SPG_Allocator::template Array<T>::mask T_MASK;
    typedef SPGrid_Set<SPG_Flags_Array_Type> SPG_Set_Type;
    typedef VECTOR<int,d> T_INDEX;
public:
    const T_INDEX domain_size;
    const T_INDEX subdomain_size;
    const SPG_Allocator& allocator;
    SPG_Set_Type& set;
    SPGrid_Subdomain_Splitter(SPG_Allocator& allocator_input,SPG_Set_Type& set_input,T_INDEX domain_size_input,T_INDEX subdomain_size_input)
        :allocator(allocator_input),set(set_input),domain_size(domain_size_input),subdomain_size(subdomain_size_input)
    {}
    int Split(INTERFACE_MATRIX_GENERATOR<T,d>& interface_matrix_generator,const int v_cycle_levels,const int sigma_mg_level,const int number_of_threads=0,
              std::vector<INTERFACE_MATRIX_GENERATOR<T,d>*>* sub_interface_matrix_generators=NULL,bool multilevel=false){
        if(!multilevel){
            std::cout<<"Domain Size: "<<domain_size<<std::endl;
            std::cout<<"SubDomain Size: "<<subdomain_size<<std::endl;
        }
        interface_matrix_generator.Initialize(sigma_mg_level);
        if(sub_interface_matrix_generators) {
            for(int i=0;i<sub_interface_matrix_generators->size();++i)
                delete (*sub_interface_matrix_generators)[i];
            sub_interface_matrix_generators->clear();}
        //This split function analyze the connectivity pattern within each subdomain, split them, linearize them and copy flags into it.
        T_INDEX number_of_subdomains=SPGrid_Computations::Domain_Divider<T_STRUCT,d>::Number_of_Subdomains(domain_size,T_INDEX(),subdomain_size);
        int subdomain_count=0;
        std::vector<T_INDEX> subdomain_id_map;
        for(RANGE_ITERATOR<d> iterator(RANGE<T_INDEX>(T_INDEX(),number_of_subdomains-1));
            iterator.Valid();
            iterator.Next()){
            subdomain_id_map.push_back(iterator.Index());}
        for(int i=0;i<subdomain_id_map.size();++i){
            int subdomain_splited=0;
            RANGE<T_INDEX> subdomain_range=SPGrid_Computations::Domain_Divider<T_STRUCT,d>::Get_Sub_Domain_Range(subdomain_id_map[i],domain_size,T_INDEX(),subdomain_size);
            LOG::cout<<"Processing subdomain "<<subdomain_range<<std::endl;
            std_array<int,d> std_size(domain_size);
            SPG_Allocator subdomain_allocator(std_size);
            SPG_Set_Type subdomain_set(subdomain_allocator.Get_Array(&T_STRUCT::flags));
            RANGE<T_INDEX> hierarchical_range(T_INDEX(),subdomain_allocator.Padded_Size().template Cast<T_INDEX>());
            Subdomain_Rasterizer<T_STRUCT,T,d> rasterizer(set,subdomain_set,subdomain_range.min_corner,subdomain_size);
            for(HIERARCHICAL_RANGE_ITERATOR<d,Subdomain_Rasterizer<T_STRUCT,T,d> > iterator(hierarchical_range,rasterizer);iterator.Valid();iterator.Next());
            subdomain_set.Refresh_Block_Offsets();
#if 0
            PHYSBAM_ASSERT(multilevel==false);
            //Connectivity Analysis
            SPG_Flags_Array_Type flags_ptr=subdomain_allocator.Get_Array(&T_STRUCT::flags);
            while(true){
                std::vector<unsigned long> coloring_queue;
                for(SPGrid_Block_Iterator<T_MASK> iterator(subdomain_set.Get_Blocks());iterator.Valid();iterator.Next()){
                    unsigned long offset=iterator.Offset();
                    if(flags_ptr(offset) & SPGrid_Solver_Cell_Type_Active) {
                        coloring_queue.push_back(offset);
                        break;}}
                if(coloring_queue.empty()) break;

                subdomain_splited++;
                subdomain_count++;
                
                while(!coloring_queue.empty()){
                    unsigned long offset=coloring_queue.back();
                    flags_ptr(offset)=SPGrid_Solver_Cell_Type_Boundary;//Using boundary flag temporarily to mark the connected part 
                    coloring_queue.pop_back();
                    std_array<int,d> coord=T_MASK::LinearToCoord(offset);
                    for(int axis=1;axis<=d;++axis){
                        for(int v=-1;v<=1;v+=2){
                            std_array<int,d> neighbor=coord;
                            neighbor(axis-1)+=v;
                            unsigned long neighbor_offset=T_MASK::Linear_Offset(neighbor);                
                            if(subdomain_set.Is_Set(neighbor_offset,SPGrid_Solver_Cell_Type_Active)){
                                coloring_queue.push_back(neighbor_offset);}}}}
                
                
                SPG_Allocator sub_allocator(std_size);
                SPG_Set_Type sub_set(sub_allocator.Get_Array(&T_STRUCT::flags));
                Subdomain_Masked_Rasterizer<T_STRUCT,T,d> masked_rasterizer(subdomain_set,sub_set,
                                                                            subdomain_range.min_corner,subdomain_size,
                                                                            SPGrid_Solver_Cell_Type_Boundary,subdomain_count);
                for(HIERARCHICAL_RANGE_ITERATOR<d,Subdomain_Masked_Rasterizer<T_STRUCT,T,d> > iterator(hierarchical_range,masked_rasterizer);iterator.Valid();iterator.Next());
                sub_set.Refresh_Block_Offsets();
            
                if(number_of_threads)
                    SPGrid_Computations::Threading_Helper<T_STRUCT,d>(sub_allocator,sub_set.Get_Blocks()).Run_Parallel(SPGrid_Computations::Domain_Divider<T_STRUCT,d>(&T_STRUCT::flags,sub_set,domain_size,T_INDEX(),subdomain_size),number_of_threads);
                else
                    SPGrid_Computations::Domain_Divider<T_STRUCT,d>(sub_allocator,
                                                                    sub_set.Get_Blocks(),
                                                                    &T_STRUCT::flags,sub_set,
                                                                    domain_size,
                                                                    T_INDEX(),subdomain_size);

                for(SPGrid_Block_Iterator<T_MASK> iterator(subdomain_set.Get_Blocks());iterator.Valid();iterator.Next()){
                    unsigned long offset=iterator.Offset();
                    if(flags_ptr(offset) & SPGrid_Solver_Cell_Type_Boundary) {
                        flags_ptr(offset)=0x0;}}
                
                SPGrid_Computations::V_Cycle_Topology<T_STRUCT,T,d>& v_cycle_topology=*(new SPGrid_Computations::V_Cycle_Topology<T_STRUCT,T,d>(sub_allocator,sub_set,domain_size,v_cycle_levels,&T_STRUCT::flags));
                v_cycle_topology.Initialize(number_of_threads);
                v_cycle_topology.Flaging(number_of_threads);
                //LOG::cout<<"subdomain_count:"<<subdomain_count<<std::endl;
                delete &v_cycle_topology;
                std::vector<GRID_HIERARCHY<T_STRUCT_ADAPTATION,T,d>*> hierarchy;
                interface_matrix_generator.template Create_Subdomain_Hierarchy<T_STRUCT_ADAPTATION,T_STRUCT>(hierarchy,sub_set,domain_size,subdomain_range);
                for(int i=0;i<sigma_mg_level;++i) delete hierarchy[i];
            }
            LOG::cout<<"Subdomain splitted: "<<subdomain_splited<<std::endl;            
#else
            SPG_Flags_Array_Type flags_ptr=subdomain_allocator.Get_Array(&T_STRUCT::flags);
            for(SPGrid_Block_Iterator<T_MASK> iterator(subdomain_set.Get_Blocks());iterator.Valid();iterator.Next()){
                unsigned long offset=iterator.Offset();
                if(flags_ptr(offset) & SPGrid_Solver_Cell_Type_Active) {
                    flags_ptr(offset)=SPGrid_Solver_Cell_Type_Boundary;}}
            SPG_Allocator sub_allocator(std_size);
            SPG_Set_Type sub_set(sub_allocator.Get_Array(&T_STRUCT::flags));
            Subdomain_Masked_Rasterizer<T_STRUCT,T,d> masked_rasterizer(subdomain_set,sub_set,
                                                                        subdomain_range.min_corner,subdomain_size,
                                                                        SPGrid_Solver_Cell_Type_Boundary,subdomain_count);
            for(HIERARCHICAL_RANGE_ITERATOR<d,Subdomain_Masked_Rasterizer<T_STRUCT,T,d> > iterator(hierarchical_range,masked_rasterizer);iterator.Valid();iterator.Next());
            sub_set.Refresh_Block_Offsets();
            
            if(number_of_threads)
                SPGrid_Computations::Threading_Helper<T_STRUCT,d>(sub_allocator,sub_set.Get_Blocks()).Run_Parallel(SPGrid_Computations::Domain_Divider<T_STRUCT,d>(&T_STRUCT::flags,sub_set,domain_size,T_INDEX(),subdomain_size),number_of_threads);
            else
                SPGrid_Computations::Domain_Divider<T_STRUCT,d>(sub_allocator,
                                                                sub_set.Get_Blocks(),
                                                                &T_STRUCT::flags,sub_set,
                                                                domain_size,
                                                                T_INDEX(),subdomain_size);

            if(multilevel){
                SPG_Allocator sub_sub_allocator(std_size);
                SPG_Set_Type sub_sub_set(sub_sub_allocator.Get_Array(&T_STRUCT::flags));
                Subdomain_Masked_Rasterizer<T_STRUCT,T,d> masked_rasterizer(subdomain_set,sub_sub_set,
                                                                            subdomain_range.min_corner,subdomain_size,
                                                                            SPGrid_Solver_Cell_Type_Boundary,subdomain_count);
                for(HIERARCHICAL_RANGE_ITERATOR<d,Subdomain_Masked_Rasterizer<T_STRUCT,T,d> > iterator(hierarchical_range,masked_rasterizer);
                    iterator.Valid();iterator.Next());
                sub_sub_set.Refresh_Block_Offsets();
            
                if(number_of_threads)
                    SPGrid_Computations::Threading_Helper<T_STRUCT,d>(sub_sub_allocator,sub_sub_set.Get_Blocks()).Run_Parallel(SPGrid_Computations::Domain_Divider<T_STRUCT,d>(&T_STRUCT::flags,sub_sub_set,domain_size,T_INDEX(),subdomain_size,true),number_of_threads);
                else
                    SPGrid_Computations::Domain_Divider<T_STRUCT,d>(sub_sub_allocator,
                                                                    sub_sub_set.Get_Blocks(),
                                                                    &T_STRUCT::flags,sub_sub_set,
                                                                    domain_size,
                                                                    T_INDEX(),subdomain_size,true);
                // Now subdivide the subdomain again.
                T_INDEX sub_subdomain_size=subdomain_size;
                for(int axis=1;axis<=d;++axis) {PHYSBAM_ASSERT(sub_subdomain_size(axis)%2==0);sub_subdomain_size(axis)/=2;}
                std::cout<<"sub_subdomain_size: "<<sub_subdomain_size<<std::endl;
                SPGrid_Subdomain_Splitter sub_splitter(sub_sub_allocator,sub_sub_set,domain_size,sub_subdomain_size);
                INTERFACE_MATRIX_GENERATOR<T,d>* new_generator=new INTERFACE_MATRIX_GENERATOR<T,d>();
                sub_interface_matrix_generators->push_back(new_generator);
                sub_splitter.Split(*new_generator,v_cycle_levels,sigma_mg_level,number_of_threads,NULL,false);}
            
            for(SPGrid_Block_Iterator<T_MASK> iterator(subdomain_set.Get_Blocks());iterator.Valid();iterator.Next()){
                unsigned long offset=iterator.Offset();
                if(flags_ptr(offset) & SPGrid_Solver_Cell_Type_Boundary){
                    flags_ptr(offset)=0x0;}}
            

            SPGrid_Computations::V_Cycle_Topology<T_STRUCT,T,d>& v_cycle_topology=*(new SPGrid_Computations::V_Cycle_Topology<T_STRUCT,T,d>(sub_allocator,sub_set,
                                                                                                                                            domain_size,v_cycle_levels,
                                                                                                                                            &T_STRUCT::flags));
            v_cycle_topology.Initialize(number_of_threads);
            v_cycle_topology.Flaging(number_of_threads);
            //LOG::cout<<"subdomain_count:"<<subdomain_count<<std::endl;
            subdomain_count++;
            delete &v_cycle_topology;
            std::vector<GRID_HIERARCHY<T_STRUCT_ADAPTATION,T,d>*> hierarchy;
            interface_matrix_generator.template Create_Subdomain_Hierarchy<T_STRUCT_ADAPTATION,T_STRUCT>(hierarchy,sub_set,domain_size,subdomain_range);
            for(int i=0;i<sigma_mg_level;++i) delete hierarchy[i];}
#endif
        return subdomain_count;
    }
};
}
#endif
