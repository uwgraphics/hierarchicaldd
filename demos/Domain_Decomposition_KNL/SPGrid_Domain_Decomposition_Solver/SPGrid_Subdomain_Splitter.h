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
#include "SPGrid_Block_Pair_Weaver.h"
#include "SPGrid_Domain_Divider_Helper.h"
#include "SPGrid_Linearized_Data_Copy_Helper.h"
#include "SPGrid_Master_Array_Linearizer.h"
#include "SPGrid_V_Cycle_Topology.h"
#include "SPGrid_Subdomain_Rasterizer.h"
#include "SPGRID_READ_WRITE.h"
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
    typedef SPGrid_Master_Array_Linearizer<T,NextLogTwo<sizeof(T_STRUCT)>::value,d,T_offset_ptr> T_Linearizer;
public:
    const T_INDEX domain_size;
    const T_INDEX subdomain_size;
    const SPG_Allocator& allocator;
    SPG_Set_Type& set;
    SPGrid_Subdomain_Splitter(SPG_Allocator& allocator_input,SPG_Set_Type& set_input,T_INDEX domain_size_input,T_INDEX subdomain_size_input)
        :allocator(allocator_input),set(set_input),domain_size(domain_size_input),subdomain_size(subdomain_size_input)
    {}
    void Split(T T_STRUCT::* rhs_field,INTERFACE_MATRIX_GENERATOR<T,d>& interface_matrix_generator,
               std::vector<std::vector<T_Linearizer>*>& linearizer_hierarchy,const int v_cycle_levels,const int sigma_mg_level,const int number_of_threads,
               std::string iomode,std::string dump_directory){
        std_array<int,d> std_size(domain_size);
        interface_matrix_generator.Initialize(sigma_mg_level);
        //This split function analyze the connectivity pattern within each subdomain, split them, linearize them and copy flags into it.
        T_INDEX number_of_subdomains=SPGrid_Computations::Domain_Divider<T_STRUCT,d>::Number_of_Subdomains(domain_size,T_INDEX(),subdomain_size);
        int subdomain_count=0;  
        int n_subdomain=0;
        linearizer_hierarchy.clear();
        std::vector<T_INDEX> subdomain_id_map;
        for(RANGE_ITERATOR<d> iterator(RANGE<T_INDEX>(T_INDEX(),number_of_subdomains-1));
            iterator.Valid();
            iterator.Next()){
            subdomain_id_map.push_back(iterator.Index());}

        if(iomode=="input"){
            std::ifstream input;
            input.open(dump_directory+"/number_of_subdomains",std::ios::out|std::ios::binary);        
            input.read(reinterpret_cast<char*>(&n_subdomain),sizeof(n_subdomain));
            input.close();
            for(int subdomain=0;subdomain<n_subdomain;++subdomain){
                SPG_Allocator* allocator;
                SPG_Set_Type* set;
                SPGRID_READ_WRITE<T_STRUCT,T,d>::Read_SPGrid(dump_directory+"/Dumped_data"+STRING_UTILITIES::Value_To_String(subdomain+1),allocator,set,rhs_field);
                SPGrid_Computations::V_Cycle_Topology<T_STRUCT,T,d>& v_cycle_topology=*(new SPGrid_Computations::V_Cycle_Topology<T_STRUCT,T,d>(*allocator,*set,domain_size,v_cycle_levels,&T_STRUCT::flags));
                v_cycle_topology.Initialize(number_of_threads);
                v_cycle_topology.Flaging(number_of_threads);
                //LOG::cout<<"subdomain_count:"<<subdomain_count<<std::endl;
                //Let's linearize the subdomain!
                std::vector<T_Linearizer>& subdomain_linearizer=*(new std::vector<T_Linearizer>());
                linearizer_hierarchy.push_back(&subdomain_linearizer);
                subdomain_linearizer.resize(v_cycle_levels);
                std::vector<std::vector<unsigned long> > subdomain_blocks(v_cycle_levels);
                for(int level=0;level<v_cycle_levels;++level){
                    std::vector<unsigned long> boundary_blocks,interface_blocks;                
                    SPGrid_Computations::Block_Pair_Weaver<T_STRUCT,d>::Weave_Subdomain(subdomain_blocks[level],
                                                                                        v_cycle_topology.Get_Allocator(level),
                                                                                        v_cycle_topology.Get_Set(level).Get_Blocks(),&T_STRUCT::flags);

                    std::pair<unsigned long*,unsigned> sub_blocks(&subdomain_blocks[level][0],subdomain_blocks[level].size());
                
                    SPGrid_Computations::Block_Pair_Weaver<T_STRUCT,d>::Weave_With_Interface_Flags(interface_blocks,
                                                                                                   v_cycle_topology.Get_Allocator(level),
                                                                                                   sub_blocks,&T_STRUCT::flags);
                
                    if(level>0) PHYSBAM_ASSERT(interface_blocks.size()==0);

                    SPGrid_Computations::Block_Pair_Weaver<T_STRUCT,d>::Weave_With_Boundary_Flags(boundary_blocks,
                                                                                                  v_cycle_topology.Get_Allocator(level),
                                                                                                  sub_blocks,&T_STRUCT::flags);
                    
                    //std::cout<<"Subdomain "<<subdomain_count<<" : "<<"boundary blocks "<<boundary_blocks.size()<<std::endl;
                    //std::cout<<"Subdomain "<<subdomain_count<<" : "<<"boundary blocks "<<boundary_blocks.size()<<std::endl;
                    subdomain_linearizer[level].Initialize(std::pair<unsigned long*,unsigned>(&subdomain_blocks[level][0],subdomain_blocks[level].size()));
                    subdomain_linearizer[level].Initialize_Boundary(std::pair<unsigned long*,unsigned>(&boundary_blocks[0],boundary_blocks.size()));
                    subdomain_linearizer[level].Initialize_Interface(std::pair<unsigned long*,unsigned>(&interface_blocks[0],interface_blocks.size()));
                    std::cout<<"interface_blocks.size(): "<<interface_blocks.size()<<std::endl;
                    SPG_Flags_Array_Type flags_ptr=v_cycle_topology.Get_Allocator(level).Get_Array(&T_STRUCT::flags);
                    #pragma omp parallel for
                    for(unsigned int i=0;i<subdomain_blocks[level].size();++i){
                        memcpy(subdomain_linearizer[level].data+(unsigned long)(i)*T_Linearizer::page_size,(char*)&flags_ptr(subdomain_blocks[level][i]),T_Linearizer::page_size);}
                    #pragma omp parallel for
                    for(unsigned int i=0;i<interface_blocks.size();++i){
                        memcpy(subdomain_linearizer[level].interface_data+(unsigned long)(i)*T_Linearizer::page_size,(char*)&flags_ptr(interface_blocks[i]),T_Linearizer::page_size);}
                }
                for(int level=0;level<v_cycle_levels-1;++level){
                    std::pair<unsigned long*,unsigned> fine_blocks(&subdomain_blocks[level][0],subdomain_blocks[level].size());
                    std::pair<unsigned long*,unsigned> coarse_blocks(&subdomain_blocks[level+1][0],subdomain_blocks[level+1].size());
                    subdomain_linearizer[level+1].PopulateRestrictionOffsets(fine_blocks,coarse_blocks,subdomain_linearizer[level]);
                    subdomain_linearizer[level].PopulateProlongationOffsets(fine_blocks,coarse_blocks,subdomain_linearizer[level+1]);
                }
                delete set;
                delete allocator;
                delete &v_cycle_topology;
            }
            int i;
            std::cin>>i;
            return;
        }
        for(int i=0;i<subdomain_id_map.size();++i){
            int subdomain_splited=0;
            RANGE<T_INDEX> subdomain_range=SPGrid_Computations::Domain_Divider<T_STRUCT,d>::Get_Sub_Domain_Range(subdomain_id_map[i],domain_size,T_INDEX(),subdomain_size);
            LOG::cout<<"Processing subdomain "<<subdomain_range<<std::endl;
            SPG_Allocator subdomain_allocator(std_size);
            SPG_Set_Type subdomain_set(subdomain_allocator.Get_Array(&T_STRUCT::flags));
            RANGE<T_INDEX> hierarchical_range(T_INDEX(),subdomain_allocator.Padded_Size().template Cast<T_INDEX>());
            Subdomain_Rasterizer<T_STRUCT,T,d> rasterizer(set,subdomain_set,subdomain_range.min_corner,subdomain_size);
            for(HIERARCHICAL_RANGE_ITERATOR<d,Subdomain_Rasterizer<T_STRUCT,T,d> > iterator(hierarchical_range,rasterizer);iterator.Valid();iterator.Next());
            subdomain_set.Refresh_Block_Offsets();
#if 0
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
                //Let's linearize the subdomain!
                std::vector<T_Linearizer>& subdomain_linearizer=*(new std::vector<T_Linearizer>());

                linearizer_hierarchy.push_back(&subdomain_linearizer);

                subdomain_linearizer.resize(v_cycle_levels);
                std::vector<std::vector<unsigned long> > subdomain_blocks(v_cycle_levels);
                for(int level=0;level<v_cycle_levels;++level){
                    std::vector<unsigned long> boundary_blocks,interface_blocks;
                    SPGrid_Computations::Block_Pair_Weaver<T_STRUCT,d>::Weave_Subdomain(subdomain_blocks[level],
                                                                                        v_cycle_topology.Get_Allocator(level),
                                                                                        v_cycle_topology.Get_Set(level).Get_Blocks(),&T_STRUCT::flags);

                    std::pair<unsigned long*,unsigned> sub_blocks(&subdomain_blocks[level][0],subdomain_blocks[level].size());
                
                    SPGrid_Computations::Block_Pair_Weaver<T_STRUCT,d>::Weave_With_Interface_Flags(interface_blocks,
                                                                                                   v_cycle_topology.Get_Allocator(level),
                                                                                                   sub_blocks,&T_STRUCT::flags);
                
                    SPGrid_Computations::Block_Pair_Weaver<T_STRUCT,d>::Weave_With_Boundary_Flags(boundary_blocks,
                                                                                                  v_cycle_topology.Get_Allocator(level),
                                                                                                  sub_blocks,&T_STRUCT::flags);
                    
                    if(level==0) LOG::cout<<"Top blocks: "<<subdomain_blocks[level].size()<<std::endl;                  
                    //if(level==v_cycle_levels-1) std::cout<<"Bottom blocks: "<<subdomain_blocks[level].size()<<std::endl;
                    //std::cout<<"Subdomain "<<subdomain_count<<" : "<<"boundary blocks "<<boundary_blocks.size()<<std::endl;
                    subdomain_linearizer[level].Initialize(std::pair<unsigned long*,unsigned>(&subdomain_blocks[level][0],subdomain_blocks[level].size()));
                    subdomain_linearizer[level].Initialize_Boundary(std::pair<unsigned long*,unsigned>(&boundary_blocks[0],boundary_blocks.size()));
                    subdomain_linearizer[level].Initialize_Interface(std::pair<unsigned long*,unsigned>(&interface_blocks[0],interface_blocks.size()));
                    SPG_Flags_Array_Type flags_ptr=v_cycle_topology.Get_Allocator(level).Get_Array(&T_STRUCT::flags);
                    for(unsigned int i=0;i<subdomain_blocks[level].size();++i){
                        memcpy(subdomain_linearizer[level].data+(unsigned long)(i)*T_Linearizer::page_size,(char*)&flags_ptr(subdomain_blocks[level][i]),T_Linearizer::page_size);}
                    for(unsigned int i=0;i<interface_blocks.size();++i){
                        memcpy(subdomain_linearizer[level].interface_data+(unsigned long)(i)*T_Linearizer::page_size,(char*)&flags_ptr(interface_blocks[i]),T_Linearizer::page_size);}
                }
                for(int level=0;level<v_cycle_levels-1;++level){
                    std::pair<unsigned long*,unsigned> fine_blocks(&subdomain_blocks[level][0],subdomain_blocks[level].size());
                    std::pair<unsigned long*,unsigned> coarse_blocks(&subdomain_blocks[level+1][0],subdomain_blocks[level+1].size());
                    subdomain_linearizer[level+1].PopulateRestrictionOffsets(fine_blocks,coarse_blocks,subdomain_linearizer[level]);
                    subdomain_linearizer[level].PopulateProlongationOffsets(fine_blocks,coarse_blocks,subdomain_linearizer[level+1]);
                }
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

            n_subdomain++;
            if(iomode=="output")
                SPGRID_READ_WRITE<T_STRUCT,T,d>::Write_SPGrid(dump_directory+"/Dumped_data"+STRING_UTILITIES::Value_To_String(n_subdomain),sub_allocator,sub_set,rhs_field);

            for(SPGrid_Block_Iterator<T_MASK> iterator(subdomain_set.Get_Blocks());iterator.Valid();iterator.Next()){
                unsigned long offset=iterator.Offset();
                if(flags_ptr(offset) & SPGrid_Solver_Cell_Type_Boundary) {
                    flags_ptr(offset)=0x0;}}
            
            SPGrid_Computations::V_Cycle_Topology<T_STRUCT,T,d>& v_cycle_topology=*(new SPGrid_Computations::V_Cycle_Topology<T_STRUCT,T,d>(sub_allocator,sub_set,domain_size,v_cycle_levels,&T_STRUCT::flags));
            v_cycle_topology.Initialize(number_of_threads);
            v_cycle_topology.Flaging(number_of_threads);
            //LOG::cout<<"subdomain_count:"<<subdomain_count<<std::endl;
            //Let's linearize the subdomain!
            std::vector<T_Linearizer>& subdomain_linearizer=*(new std::vector<T_Linearizer>());
            linearizer_hierarchy.push_back(&subdomain_linearizer);
            subdomain_linearizer.resize(v_cycle_levels);
            std::vector<std::vector<unsigned long> > subdomain_blocks(v_cycle_levels);
            for(int level=0;level<v_cycle_levels;++level){
                std::vector<unsigned long> boundary_blocks,interface_blocks;                
                SPGrid_Computations::Block_Pair_Weaver<T_STRUCT,d>::Weave_Subdomain(subdomain_blocks[level],
                                                                                    v_cycle_topology.Get_Allocator(level),
                                                                                    v_cycle_topology.Get_Set(level).Get_Blocks(),&T_STRUCT::flags);

                std::pair<unsigned long*,unsigned> sub_blocks(&subdomain_blocks[level][0],subdomain_blocks[level].size());
                
                SPGrid_Computations::Block_Pair_Weaver<T_STRUCT,d>::Weave_With_Interface_Flags(interface_blocks,
                                                                                               v_cycle_topology.Get_Allocator(level),
                                                                                               sub_blocks,&T_STRUCT::flags);
                
                if(level>0) PHYSBAM_ASSERT(interface_blocks.size()==0);

                SPGrid_Computations::Block_Pair_Weaver<T_STRUCT,d>::Weave_With_Boundary_Flags(boundary_blocks,
                                                                                              v_cycle_topology.Get_Allocator(level),
                                                                                              sub_blocks,&T_STRUCT::flags);
                    
                //std::cout<<"Subdomain "<<subdomain_count<<" : "<<"boundary blocks "<<boundary_blocks.size()<<std::endl;
                //std::cout<<"Subdomain "<<subdomain_count<<" : "<<"boundary blocks "<<boundary_blocks.size()<<std::endl;
                subdomain_linearizer[level].Initialize(std::pair<unsigned long*,unsigned>(&subdomain_blocks[level][0],subdomain_blocks[level].size()));
                subdomain_linearizer[level].Initialize_Boundary(std::pair<unsigned long*,unsigned>(&boundary_blocks[0],boundary_blocks.size()));
                subdomain_linearizer[level].Initialize_Interface(std::pair<unsigned long*,unsigned>(&interface_blocks[0],interface_blocks.size()));
                SPG_Flags_Array_Type flags_ptr=v_cycle_topology.Get_Allocator(level).Get_Array(&T_STRUCT::flags);
                for(unsigned int i=0;i<subdomain_blocks[level].size();++i){
                    memcpy(subdomain_linearizer[level].data+(unsigned long)(i)*T_Linearizer::page_size,(char*)&flags_ptr(subdomain_blocks[level][i]),T_Linearizer::page_size);}
                for(unsigned int i=0;i<interface_blocks.size();++i){
                    memcpy(subdomain_linearizer[level].interface_data+(unsigned long)(i)*T_Linearizer::page_size,(char*)&flags_ptr(interface_blocks[i]),T_Linearizer::page_size);}
            }
            for(int level=0;level<v_cycle_levels-1;++level){
                std::pair<unsigned long*,unsigned> fine_blocks(&subdomain_blocks[level][0],subdomain_blocks[level].size());
                std::pair<unsigned long*,unsigned> coarse_blocks(&subdomain_blocks[level+1][0],subdomain_blocks[level+1].size());
                subdomain_linearizer[level+1].PopulateRestrictionOffsets(fine_blocks,coarse_blocks,subdomain_linearizer[level]);
                subdomain_linearizer[level].PopulateProlongationOffsets(fine_blocks,coarse_blocks,subdomain_linearizer[level+1]);
            }
            delete &v_cycle_topology;
            std::vector<GRID_HIERARCHY<T_STRUCT_ADAPTATION,T,d>*> hierarchy;
            interface_matrix_generator.template Create_Subdomain_Hierarchy<T_STRUCT_ADAPTATION,T_STRUCT>(hierarchy,sub_set,domain_size,subdomain_range);
            for(int i=0;i<sigma_mg_level;++i) delete hierarchy[i];
#endif
        }
    }
};
}
#endif
