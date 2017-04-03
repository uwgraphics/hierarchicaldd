//#####################################################################
//Removed For Anonymity: Copyright Authors of Submission pap173s1
//Removed For Anonymity: Copyright Authors of Submission pap173s1
//#####################################################################
// Class Interface_Matrix_Generator
//#####################################################################
#ifndef __INTERFACE_MATRIX_GENERATOR_H__
#define __INTERFACE_MATRIX_GENERATOR_H__
#include <Common_Tools/Grids_Uniform_PDE_Linear/STENCIL.h>
#include <Common_Tools/Grids_Uniform_PDE_Linear/STENCIL_ITERATOR.h>

#include <SPGrid_Fluids/Grids/GRID_HIERARCHY_ITERATOR.h>
#include <SPGrid_Fluids/Simulation/FLUIDS_SIMULATION_DATA.h>
#include <SPGrid_Fluids/Simulation/FLUIDS_SIMULATION_FLAGS.h>

#include <PhysBAM_Tools/Data_Structures/HASHTABLE.h>
#include <PhysBAM_Tools/Data_Structures/HASHTABLE_ITERATOR.h>

#include "ADAPTIVE_SUBDOMAIN_POISSON_FLAGS.h"
#include "BOUNDARY_REFINEMENT_RASTERIZER.h"
#include "CELL_PROLONGATION_BAKER.h"
#include "INTERFACE_SOLVER_DATA.h"
#include "GALERKIN_PROCESS.h"
#include "PROLONGATION_MATRIX_HELPER.h"
#include "SPGRID_WRITER.h"
#include "STRIDE_GENERATOR.h"

using namespace PhysBAM;

namespace{
template<typename T_MASK,int d> class BLOCK_SIZE_HELPER;
template<typename T_MASK> class BLOCK_SIZE_HELPER<T_MASK,2>{
public:
    static VECTOR<int,2> Get_Block_Size(){
        return VECTOR<int,2>(1<<T_MASK::block_xbits,1<<T_MASK::block_ybits);    
    }
};

template<typename T_MASK> class BLOCK_SIZE_HELPER<T_MASK,3>{
public:
    static VECTOR<int,3> Get_Block_Size(){
        return VECTOR<int,3>(1<<T_MASK::block_xbits,1<<T_MASK::block_ybits,1<<T_MASK::block_zbits);    
    }
};
}

template<class T,int d>
class INTERFACE_MATRIX_GENERATOR{
public:
    typedef VECTOR<int,d> T_INDEX;
    typedef VECTOR<T,d> TV;
    typedef STENCIL<T,1> T_STENCIL;
    int mg_levels;
    std::vector<HASHTABLE<T_INDEX,int> > interface_hash;
    std::vector<T_INDEX> interface_map;
    std::vector<int> interface_dof_counter;
    typedef Eigen::Triplet<T,int> TRIPLET;

    std::vector<std::vector<T_STENCIL*> > Arr;//We use ptr here because we don't know the dimension of Arr until the last subdomain is processed, use ptr to *cheaply* extend the vector.
    std::vector<std::vector<std::vector<T_STENCIL> >*> Air;    
    std::vector<std::vector<std::vector<T_STENCIL> >*> Aii;

    void Embedding_Helper(const VECTOR<int,d>& index,STENCIL<T,d>& parents,const T weight,const VECTOR<ARRAY<int,VECTOR<int,1> >,d>& q_adapted)
    {
        typedef VECTOR<int,d> T_INDEX;
        T_INDEX q_stride_v,min_corner;
        for(int axis=1;axis<=d;++axis){
            q_stride_v(axis) = q_adapted(axis)(index(axis));
            min_corner(axis) = q_adapted(axis).Domain_Indices().min_corner(1);}
        int q_stride=q_stride_v.Min();
        for(int axis=1;axis<=d;++axis){int current_index=index(axis)-min_corner(axis);
            if(current_index%q_stride!=0){
                T_INDEX index1=index,index2=index;
                index1(axis)=(current_index/q_stride)*q_stride+min_corner(axis);
                index2(axis)=(current_index/q_stride+1)*q_stride+min_corner(axis);
                T weight1=weight*(index2(axis)-index(axis))/(index2(axis)-index1(axis));
                T weight2=weight*(index(axis)-index1(axis))/(index2(axis)-index1(axis));
                Embedding_Helper(index1,parents,weight1,q_adapted);Embedding_Helper(index2,parents,weight2,q_adapted);
                return;}}
        // If we got so far, it must be real node
        parents.Get_Or_Insert(index)+=weight;
    }

    ~INTERFACE_MATRIX_GENERATOR(){
        for(int i=0;i<Arr.size();++i)
            for(int j=0;j<Arr[i].size();++j)
                delete Arr[i][j];
        for(int i=0;i<Aii.size();++i)
            delete Aii[i];
        for(int i=0;i<Air.size();++i)
            delete Air[i];
    }
    void Initialize(int mg_levels_input){
        mg_levels=mg_levels_input;
        for(int i=0;i<Arr.size();++i)
            for(int j=0;j<Arr[i].size();++j)
                delete Arr[i][j];
        for(int i=0;i<Aii.size();++i)
            delete Aii[i];
        for(int i=0;i<Air.size();++i)
            delete Air[i];
        Aii.clear();
        Air.clear();
        Arr.resize(mg_levels);
        interface_hash.resize(mg_levels);
        interface_dof_counter.resize(mg_levels);
        for(int level=0;level<mg_levels;++level){
            interface_dof_counter[level]=0;
            interface_hash[level].Remove_All();}
    }
    template<class T_STRUCT_ADAPTATION,class T_STRUCT_SOLVER>
    void Create_Subdomain_Hierarchy(std::vector<GRID_HIERARCHY<T_STRUCT_ADAPTATION,T,d>*>& hierarchy,
                                    const SPGrid_Set<typename SPGrid_Allocator<T_STRUCT_SOLVER,d>::template Array<unsigned>::type>& set,
                                    T_INDEX size,RANGE<T_INDEX> subdomain_range){
        typedef SPGrid_Allocator<T_STRUCT_SOLVER,d> SPG_Allocator;
        typedef typename SPG_Allocator::template Array<unsigned>::type SPG_Flags_Array_Type;
        typedef typename SPG_Allocator::template Array<T>::type SPG_Data_Array_Type;
        typedef typename SPG_Allocator::template Array<T>::mask T_MASK;
        typedef SPGrid_Set<SPG_Flags_Array_Type> SPG_Set_Type;
        typedef SPGrid_Allocator<T_STRUCT_ADAPTATION,d> SPG_Allocator_Adaptation;
        typedef typename SPG_Allocator_Adaptation::template Array<unsigned>::type SPG_Flags_Array_Type_Adaptation;
        typedef typename SPG_Allocator_Adaptation::template Array<int>::type SPG_Index_Array_Type_Adaptation;
        typedef typename SPG_Allocator_Adaptation::template Array<unsigned>::mask T_MASK_ADAPTATION;
        typedef typename SPG_Allocator_Adaptation::template Array<int>::mask T_INDEX_MASK_ADAPTATION;
        // Baking the galerkin matrices for the cells.
        std::vector<std::vector<T_STENCIL> >& Aii_subdomain=*(new std::vector<std::vector<T_STENCIL> >(mg_levels));
        std::vector<std::vector<T_STENCIL> >& Air_subdomain=*(new std::vector<std::vector<T_STENCIL> >(mg_levels));
        Aii.push_back(&Aii_subdomain);
        Air.push_back(&Air_subdomain);
        CELL_PROLONGATION_BAKER<T,d> prolongation_baker;
        prolongation_baker.Bake();
        // std::cout<<prolongation_baker.prolongation_matrices[0]<<std::endl;
        // std::cout<<prolongation_baker.prolongation_matrices[1]<<std::endl;
        // std::cout<<prolongation_baker.prolongation_matrices[2]<<std::endl;
        // std::cout<<prolongation_baker.prolongation_matrices[3]<<std::endl;
        // std::cout<<prolongation_baker.prolongation_matrices[4]<<std::endl;
        // std::cout<<prolongation_baker.prolongation_matrices[5]<<std::endl;
        // std::cout<<prolongation_baker.prolongation_matrices[6]<<std::endl;
        // std::cout<<prolongation_baker.prolongation_matrices[7]<<std::endl;
        typedef ARRAY<int,VECTOR<int,1> > ARRAY_1D;
        VECTOR<ARRAY_1D,d> q_top;
        int levels=STRIDE_GENERATOR<d>::Generate_Adaptive_Stride(subdomain_range,q_top);
        PHYSBAM_ASSERT(mg_levels<=levels);
        ARRAY<VECTOR<ARRAY_1D,d> > q(levels);
        for(int level=1;level<=levels;++level){
            RANGE<T_INDEX> subdomain_range_corrected(subdomain_range);
            for(int axis=1;axis<=d;++axis){
                subdomain_range_corrected.min_corner(axis)>>=(level-1);
                subdomain_range_corrected.max_corner(axis)>>=(level-1);}
            STRIDE_GENERATOR<d>::Generate_Adaptive_Stride(subdomain_range_corrected,q(level));}

   
        //std::cout<<"Size: "<<size<<std::endl;
        //std::cout<<"Subdomain range: "<<subdomain_range<<std::endl;
        //size is 2^n +1. 
        GRID<TV> base_grid(T_INDEX::All_Ones_Vector()*size,RANGE<TV>(TV(),TV(size-1)));
        //std::cout<<levels<<std::endl;
        //std::cout<<base_grid.Numbers_Of_Cells()<<std::endl;
        hierarchy.resize(mg_levels);
        for(int level=1;level<=mg_levels;++level){
            hierarchy[level-1]=new GRID_HIERARCHY<T_STRUCT_ADAPTATION,T,d>(base_grid,levels);
            hierarchy[level-1]->Initialize_Sets();}

        //std::cout<<hierarchy[0]->Allocator(1).Padded_Size()<<std::endl;
#if 1
        {LOG::SCOPE scope("Generating Adaptation Hierarchy");
        SPG_Flags_Array_Type flags_ptr=set.array;
        for(SPGrid_Block_Iterator<T_MASK> block_iterator(set.Get_Blocks());block_iterator.Valid();block_iterator.Next_Block()){
            unsigned long offset=block_iterator.Offset();
            T_INDEX base_index=block_iterator.Index().template Cast<T_INDEX>();
            for(RANGE_ITERATOR<d> iterator(RANGE<T_INDEX>(base_index,base_index+BLOCK_SIZE_HELPER<T_MASK,d>::Get_Block_Size()-1));
                iterator.Valid();
                iterator.Next(),offset+=sizeof(unsigned)){
                const T_INDEX& node_index=iterator.Index();
                const unsigned node_flag=flags_ptr(offset);
                if(node_flag==0) continue;
                for(RANGE_ITERATOR<d> cell_iterator(RANGE<VECTOR<int,d> >(node_index,node_index+1));cell_iterator.Valid();cell_iterator.Next()){
                    const T_INDEX& cell_index=cell_iterator.Index();                
                    if((!subdomain_range.Lazy_Inside(cell_index-1))||(!subdomain_range.Lazy_Inside(cell_index))) continue;//It is a outside cell.
                    unsigned long cell_offset=(T_MASK_ADAPTATION::Linear_Offset(std_array<int,d>(cell_index)));
                    if(hierarchy[0]->Set(1).array(cell_offset)&Uniform_Cell_Traversed) continue;
                    //hierarchy[0]->Set(1).array(cell_offset)|=Uniform_Cell_Traversed;                            
                    hierarchy[0]->Set(1).Mask(cell_offset,Uniform_Cell_Traversed);
                    bool go_next=false;
                    for(int v=1;v<=d;v++)
                        for(RANGE_ITERATOR<d-1> face_iterator(RANGE<VECTOR<int,d-1> >(cell_index.Remove_Index(v)-1,cell_index.Remove_Index(v)));
                            face_iterator.Valid();face_iterator.Next()){
                            if(go_next) continue;
                            std_array<int,d> lower_index(face_iterator.Index().Insert(cell_index(v)-1,v));
                            std_array<int,d> upper_index(face_iterator.Index().Insert(cell_index(v),v));
                            unsigned lower_type=0,upper_type=0;
                            if(set.Is_Set(lower_index,SPGrid_Solver_Cell_Type_Active|SPGrid_Solver_Cell_Type_Interface))
                                lower_type |= Uniform_Node_Type_Active;
                            if(set.Is_Set(lower_index,SPGrid_Solver_Cell_Type_Dirichlet))
                                lower_type |= Uniform_Node_Type_Dirichlet;
                            if(set.Is_Set(upper_index,SPGrid_Solver_Cell_Type_Active|SPGrid_Solver_Cell_Type_Interface))
                                upper_type |= Uniform_Node_Type_Active;
                            if(set.Is_Set(upper_index,SPGrid_Solver_Cell_Type_Dirichlet))
                                upper_type |= Uniform_Node_Type_Dirichlet;

                            if((lower_type&upper_type&Uniform_Node_Type_Active)||
                               ((lower_type|upper_type)&Uniform_Node_Type_Active)&&((lower_type|upper_type)&Uniform_Node_Type_Dirichlet)){
                                T_INDEX cur_index=cell_index;
                                unsigned long cur_offset=cell_offset;
                                RANGE<T_INDEX> cur_range(cur_index-1,cur_index);
                                int cur_level=1;
                                do{
                                    for(int mg_level=1;mg_level<=mg_levels;++mg_level){
                                        if((cur_level==mg_level&&(!subdomain_range.Thickened(-(1<<cur_level)).Contains(cur_range)))||
                                           (cur_level>mg_level&&subdomain_range.Thickened(-(1<<(cur_level-1))).Contains(cur_range)&&
                                            (!subdomain_range.Thickened(-(1<<cur_level)).Contains(cur_range)))){
                                            hierarchy[mg_level-1]->Activate_Cell(cur_level,cur_index,Adaptive_Cell_Type_Interior);
                                        }
                                    }
                                    cur_range.min_corner=(cur_range.min_corner/(1<<cur_level))*(1<<cur_level);
                                    cur_range.max_corner=cur_range.min_corner+(1<<cur_level);
                                    cur_index=(cur_index+1)/2;
                                    ++cur_level;
                                }while(cur_level<=levels);
                                go_next=true;}}}}}
        }
#endif
        for(int level=1;level<=mg_levels;++level){
            T_INDEX adaptation_block_size=hierarchy[level-1]->Allocator(level).Block_Size().template Cast<T_INDEX>();
            //BOUNDARY_REFINEMENT_RASTERIZER<T_STRUCT_ADAPTATION,T_STRUCT_SOLVER,T,d> rasterizer(*hierarchy[level-1],set,subdomain_range,level);
            //{LOG::SCOPE scope("Minimally Refine.");
            //for(GRID_HIERARCHY_ITERATOR<d,BOUNDARY_REFINEMENT_RASTERIZER<T_STRUCT_ADAPTATION,T_STRUCT_SOLVER,T,d> > hierarchy_iterator(hierarchy[level-1]->Grid(levels).Cell_Indices(),levels,rasterizer);hierarchy_iterator.Valid();hierarchy_iterator.Next());}

            hierarchy[level-1]->Update_Block_Offsets();            
            SPG_Index_Array_Type_Adaptation index=hierarchy[level-1]->Allocator(level).Get_Array(&T_STRUCT_ADAPTATION::node_index); 
            if(level==1){
                //If it is the first level, use the original set to check if a interface is a DOF.
                for(int axis=1;axis<=d;++axis)
                    for(int face=0;face<=1;++face)
                        for(RANGE_ITERATOR<d-1> face_iterator(RANGE<VECTOR<int,d-1> >(subdomain_range.min_corner.Remove_Index(axis),
                                                                                      subdomain_range.max_corner.Remove_Index(axis)));
                            face_iterator.Valid();face_iterator.Next()){
                            std_array<int,d> node_index;
                            if(face==0)
                                node_index=std_array<int,d>(face_iterator.Index().Insert(subdomain_range.min_corner(axis),axis));
                            else
                                node_index=std_array<int,d>(face_iterator.Index().Insert(subdomain_range.max_corner(axis),axis));

                            if(set.Is_Set(node_index,SPGrid_Solver_Cell_Type_Interface)){
                                int interface_id=interface_hash[level-1].Get_Default(node_index.template Cast<T_INDEX>(),-1);
                                if(interface_id==-1){
                                    interface_hash[level-1].Insert(node_index.template Cast<T_INDEX>(),interface_dof_counter[level-1]);
                                    interface_id=interface_dof_counter[level-1];
                                    interface_map.push_back(node_index.template Cast<T_INDEX>());
                                    Arr[level-1].push_back(new T_STENCIL());//Append the matrix
                                    interface_dof_counter[level-1]++;}
                                index(node_index)=interface_id;
                                hierarchy[level-1]->Set(level).Mask(node_index,Uniform_Node_Type_Interface|Adaptive_Node_Degree_Marker);
                            }
                            PHYSBAM_ASSERT(!set.Is_Set(node_index,SPGrid_Solver_Cell_Type_Active));
                        }
            }else{
                RANGE<T_INDEX> subdomain_range_corrected(subdomain_range);
                for(int axis=1;axis<=d;++axis){
                    subdomain_range_corrected.min_corner(axis)>>=(level-1);
                    subdomain_range_corrected.max_corner(axis)>>=(level-1);}
                //std::cout<<subdomain_range_corrected<<std::endl;
                //It is not the first level, we have to use Adaptive_Cell_Type_Active to check the interface DOF.
                SPG_Flags_Array_Type_Adaptation flags=hierarchy[level-1]->Allocator(level).Get_Array(&T_STRUCT_ADAPTATION::flags);
                for(SPGrid_Block_Iterator<T_MASK_ADAPTATION> iterator(hierarchy[level-1]->Set(level).Get_Blocks());iterator.Valid();iterator.Next_Block()){
                    unsigned long offset=iterator.Offset();
                    T_INDEX base_index=iterator.Index().template Cast<T_INDEX>();
                    for(RANGE_ITERATOR<d> iterator(RANGE<T_INDEX>(base_index,base_index+adaptation_block_size-1));
                        iterator.Valid();
                        iterator.Next(),offset+=sizeof(unsigned)){                        
                        const T_INDEX& cell_index=iterator.Index();
                        if(flags(offset)&Adaptive_Cell_Type_Interior){
                            for(RANGE_ITERATOR<d> node_iterator(RANGE<T_INDEX>::Unit_Box());
                                node_iterator.Valid();
                                node_iterator.Next()){
                                std_array<int,d> node_index(cell_index+node_iterator.Index()-1);
                                bool is_interface=false;
                                for(int axis=1;axis<=d;++axis){
                                    if(node_index(axis-1)==subdomain_range_corrected.min_corner(axis)||node_index(axis-1)==subdomain_range_corrected.max_corner(axis))
                                        is_interface=true;}
                                if(!is_interface) continue;
                                int interface_id=interface_hash[level-1].Get_Default(node_index.template Cast<T_INDEX>(),-1);
                                if(interface_id==-1){
                                    interface_hash[level-1].Insert(node_index.template Cast<T_INDEX>(),interface_dof_counter[level-1]);
                                    interface_id=interface_dof_counter[level-1];
                                    Arr[level-1].push_back(new T_STENCIL());//Append the matrix
                                    interface_dof_counter[level-1]++;}
                                index(node_index)=interface_id;
                                hierarchy[level-1]->Set(level).Mask(node_index,Uniform_Node_Type_Interface|Adaptive_Node_Degree_Marker);
                            }
                        }
                    }
                }
            }
            hierarchy[level-1]->Update_Block_Offsets();
            //Now deal with the interior DOF
            int node_count=0;
            for(int hierarchy_level=level;hierarchy_level<=levels;++hierarchy_level){
                SPG_Flags_Array_Type_Adaptation flags=hierarchy[level-1]->Allocator(hierarchy_level).Get_Array(&T_STRUCT_ADAPTATION::flags);
                SPG_Index_Array_Type_Adaptation index=hierarchy[level-1]->Allocator(hierarchy_level).Get_Array(&T_STRUCT_ADAPTATION::node_index); 
                for(SPGrid_Block_Iterator<T_MASK_ADAPTATION> iterator(hierarchy[level-1]->Set(hierarchy_level).Get_Blocks());iterator.Valid();iterator.Next_Block()){
                    unsigned long offset=iterator.Offset();
                    T_INDEX base_index=iterator.Index().template Cast<T_INDEX>();
                    for(RANGE_ITERATOR<d> iterator(RANGE<T_INDEX>(base_index,base_index+adaptation_block_size-1));
                        iterator.Valid();
                        iterator.Next(),offset+=sizeof(unsigned)){
                        const unsigned& flag=flags(offset);
                        const T_INDEX& cell_index=iterator.Index();
                        T_INDEX decendent_cell_index=2*cell_index-1;
                        if(flag&Adaptive_Cell_Type_Interior){
                            RANGE<T_INDEX> cell_range(cell_index-1,cell_index);
                            for(int i=hierarchy_level;i>1;--i){
                                cell_range.min_corner=2*cell_range.min_corner;
                                cell_range.max_corner=2*cell_range.max_corner;}
                            for(RANGE_ITERATOR<d> node_iterator(RANGE<T_INDEX>(cell_index-1,cell_index));node_iterator.Valid();node_iterator.Next()){
                                const T_INDEX& node_index=node_iterator.Index();
                                T_INDEX fine_node_index=decendent_cell_index-1+2*(node_index-(cell_index-1));
                                PHYSBAM_ASSERT(fine_node_index==node_index*2);
                                if(hierarchy_level==1&&(!set.Is_Set(std_array<int,d>(node_index),SPGrid_Solver_Cell_Type_Active))) continue;
                                if(hierarchy_level==level&&(hierarchy[level-1]->Set(level).Is_Set(std_array<int,d>(node_index),Uniform_Node_Type_Interface))) continue;
                                T_INDEX q_stride_v,min_corner;
                                for(int axis=1;axis<=d;++axis){
                                    q_stride_v(axis)=q(hierarchy_level)(axis)(node_index(axis));
                                    min_corner(axis)=q(hierarchy_level)(axis).Domain_Indices().min_corner(1);}
                                int q_stride=q_stride_v.Min();
                                bool is_t_junction=false;
                                for(int axis=1;axis<=d;++axis){
                                    int current_index=node_index(axis)-min_corner(axis);
                                    if(current_index%q_stride!=0){
                                        is_t_junction=true;}}
                                //IS ONE LEVEL UP ENOUGH?
                                if(!is_t_junction){
                                    if(hierarchy_level>level){
                                        if(hierarchy[level-1]->Set(hierarchy_level-1).Is_Set(std_array<int,d>(fine_node_index),Adaptive_Node_Degree_Marker)){
                                            index(T_INDEX_MASK_ADAPTATION::Linear_Offset(std_array<int,d>(node_index)))=hierarchy[level-1]->Allocator(hierarchy_level-1).Get_Array(&T_STRUCT_ADAPTATION::node_index)(std_array<int,d>(fine_node_index));
                                            hierarchy[level-1]->Activate_Cell(hierarchy_level,node_index,Adaptive_Node_Degree_Marker|Adaptive_Node_Type_Coarse_Shared);}} 
                                    if(!hierarchy[level-1]->Set(hierarchy_level).Is_Set(std_array<int,d>(node_index),Adaptive_Node_Degree_Marker)){
                                        hierarchy[level-1]->Activate_Cell(hierarchy_level,node_index,Adaptive_Node_Degree_Marker|Adaptive_Node_Type_Active);
                                        index(T_INDEX_MASK_ADAPTATION::Linear_Offset(std_array<int,d>(node_index)))=node_count++;}}
                                else{
                                    hierarchy[level-1]->Activate_Cell(hierarchy_level,node_index,Adaptive_Node_Type_T_Junction);
                                    STENCIL<T,d> parents;
                                    Embedding_Helper(node_index,parents,(T)1.,q(hierarchy_level));
                                    //Activate it's parents.
                                    for(STENCIL_ITERATOR<T,d> stencil_iterator(parents);stencil_iterator.Valid();stencil_iterator.Next()){
                                        const T_INDEX& parent_index=stencil_iterator.Key();
                                        if(!hierarchy[level-1]->Set(hierarchy_level).Is_Set(std_array<int,d>(parent_index),Adaptive_Node_Degree_Marker)){
                                            hierarchy[level-1]->Activate_Cell(hierarchy_level,parent_index,Adaptive_Node_Degree_Marker|Adaptive_Node_Type_Active);
                                            index(T_INDEX_MASK_ADAPTATION::Linear_Offset(std_array<int,d>(parent_index)))=node_count++;}}
                                }
                            }                    
                        }
                    }
                }
            }
            LOG::cout<<"Interior DOF: "<<node_count<<std::endl;
            hierarchy[level-1]->Update_Block_Offsets();
            //Don't clear the bitmaps. The Galerkin function uses Is_Set function
            //hierarchy[level-1]->Clear_Bitmaps();
            
            Aii_subdomain[level-1].resize(node_count);
            Air_subdomain[level-1].resize(node_count);

            typedef float RW;
            STREAM_TYPE stream_type((RW()));
            typedef FLUIDS_SIMULATION_DATA<T> T_STRUCT_VISUALIZATION;
            //SPGRID_WRITER<T_STRUCT_ADAPTATION,T_STRUCT_VISUALIZATION,T,d>::Write_SPGrid(stream_type,base_grid,*hierarchy[level-1],"output",level-1);
            //exit(1);
        }
#if 0
        {LOG::SCOPE scope("Galerin Process");
        static T edge_weight=1.0/(1<<(d-1));
        SPG_Flags_Array_Type flags_ptr=set.array;
        for(SPGrid_Block_Iterator<T_MASK> block_iterator(set.Get_Blocks());block_iterator.Valid();block_iterator.Next_Block()){
            unsigned long offset=block_iterator.Offset();
            T_INDEX base_index=block_iterator.Index().template Cast<T_INDEX>();
            for(RANGE_ITERATOR<d> iterator(RANGE<T_INDEX>(base_index,base_index+BLOCK_SIZE_HELPER<T_MASK,d>::Get_Block_Size()-1));
                iterator.Valid();
                iterator.Next(),offset+=sizeof(unsigned)){
                const T_INDEX& node_index=iterator.Index();
                const unsigned node_flag=flags_ptr(offset);
                if(node_flag==0) continue;
                for(RANGE_ITERATOR<d> cell_iterator(RANGE<VECTOR<int,d> >(node_index,node_index+1));cell_iterator.Valid();cell_iterator.Next()){
                    const T_INDEX& cell_index=cell_iterator.Index();                
                    if((!subdomain_range.Lazy_Inside(cell_index-1))||(!subdomain_range.Lazy_Inside(cell_index))) continue;//It is a outside cell.
                    unsigned long cell_offset=(T_MASK_ADAPTATION::Linear_Offset(std_array<int,d>(cell_index)));
                    if(hierarchy[0]->Set(1).array(cell_offset)&Uniform_Cell_Coarsened) continue;
                    hierarchy[0]->Set(1).array(cell_offset)|=Uniform_Cell_Coarsened;
                    Eigen::Matrix<T,1<<d,1<<d> matrix;
                    matrix.setZero();
                    bool cell_contributing=false;
                    for(int v=1;v<=d;v++){
                        for(RANGE_ITERATOR<d-1> face_iterator(RANGE<VECTOR<int,d-1> >(cell_index.Remove_Index(v)-1,cell_index.Remove_Index(v)));
                            face_iterator.Valid();face_iterator.Next()){
                            std_array<int,d> lower_index(face_iterator.Index().Insert(cell_index(v)-1,v));
                            std_array<int,d> upper_index(face_iterator.Index().Insert(cell_index(v),v));
                            const int lower_id=prolongation_baker.index_map(lower_index.template Cast<T_INDEX>()+1-cell_index);
                            const int upper_id=prolongation_baker.index_map(upper_index.template Cast<T_INDEX>()+1-cell_index);
                            unsigned lower_type=0,upper_type=0;
                            if(set.Is_Set(lower_index,SPGrid_Solver_Cell_Type_Active|SPGrid_Solver_Cell_Type_Interface))
                                lower_type |= Uniform_Node_Type_Active;
                            if(set.Is_Set(lower_index,SPGrid_Solver_Cell_Type_Dirichlet))
                                lower_type |= Uniform_Node_Type_Dirichlet;
                            if(set.Is_Set(upper_index,SPGrid_Solver_Cell_Type_Active|SPGrid_Solver_Cell_Type_Interface))
                                upper_type |= Uniform_Node_Type_Active;
                            if(set.Is_Set(upper_index,SPGrid_Solver_Cell_Type_Dirichlet))
                                upper_type |= Uniform_Node_Type_Dirichlet;
                            
                            if(lower_type&upper_type&Uniform_Node_Type_Active){
                                matrix(lower_id,upper_id)-=edge_weight;
                                matrix(lower_id,lower_id)+=edge_weight;
                                matrix(upper_id,upper_id)+=edge_weight;
                                matrix(upper_id,lower_id)-=edge_weight;
                                cell_contributing=true;}

                            if((lower_type&Uniform_Node_Type_Active)&&
                               (upper_type&Uniform_Node_Type_Dirichlet)){
                                matrix(lower_id,lower_id)+=edge_weight;
                                cell_contributing=true;}
                            
                            if((lower_type&Uniform_Node_Type_Dirichlet)&&
                               (upper_type&Uniform_Node_Type_Active)){
                                matrix(upper_id,upper_id)+=edge_weight;
                                cell_contributing=true;}                            
                        }
                    }
                    if(cell_contributing){
                        T_INDEX cur_index=cell_index;
                        unsigned long cur_offset=cell_offset;
                        RANGE<T_INDEX> cur_range(cur_index-1,cur_index);
                        int cur_level=1;
                        Eigen::Matrix<T,1<<d,1<<d>& cur_matrix=matrix;
                        do{
                            #pragma omp parallel for
                            for(int mg_level=1;mg_level<=mg_levels;++mg_level){
                                if((cur_level==mg_level&&(!subdomain_range.Thickened(-(1<<cur_level)).Contains(cur_range)))||
                                   (cur_level>mg_level&&subdomain_range.Thickened(-(1<<(cur_level-1))).Contains(cur_range)&&
                                    (!subdomain_range.Thickened(-(1<<cur_level)).Contains(cur_range)))){
                                    SPG_Flags_Array_Type_Adaptation flags=hierarchy[mg_level-1]->Allocator(cur_level).Get_Array(&T_STRUCT_ADAPTATION::flags);
                                    SPG_Index_Array_Type_Adaptation index=hierarchy[mg_level-1]->Allocator(cur_level).Get_Array(&T_STRUCT_ADAPTATION::node_index);
                                            
                                    //PHYSBAM_ASSERT(flags(std_array<int,d>(cur_index))&Adaptive_Cell_Type_Interior);
                               
                                    std::vector<STENCIL<T,d> > parents(2<<d);
                                    for(RANGE_ITERATOR<d> node_iterator(RANGE<T_INDEX>::Unit_Box());node_iterator.Valid();node_iterator.Next()){           
                                        T_INDEX node_index=cur_index-1+node_iterator.Index();
                                        Embedding_Helper(node_index,parents[prolongation_baker.index_map(node_iterator.Index())],(T)1.,q(cur_level));}
                                    //We always compute the embedding map. If it is a real DOF, the embedding map would just be 1.0 of itself
                                    for(RANGE_ITERATOR<d> node_1_iterator(RANGE<T_INDEX>::Unit_Box());node_1_iterator.Valid();node_1_iterator.Next()){
                                        for(RANGE_ITERATOR<d> node_2_iterator(RANGE<T_INDEX>::Unit_Box());node_2_iterator.Valid();node_2_iterator.Next()){
                                            const int node_1_local_index=prolongation_baker.index_map(node_1_iterator.Index());
                                            const int node_2_local_index=prolongation_baker.index_map(node_2_iterator.Index());
                                            T matrix_entry=cur_matrix(node_1_local_index,node_2_local_index);
                                            if(matrix_entry==0) continue;
                                            for(STENCIL_ITERATOR<T,d> stencil_1_iterator(parents[node_1_local_index]);
                                                stencil_1_iterator.Valid();stencil_1_iterator.Next()){
                                                for(STENCIL_ITERATOR<T,d> stencil_2_iterator(parents[node_2_local_index]);
                                                    stencil_2_iterator.Valid();stencil_2_iterator.Next()){
                                                    T_INDEX index_1=stencil_1_iterator.Key();
                                                    T_INDEX index_2=stencil_2_iterator.Key();
                                                    unsigned long flag_offset_1=T_MASK_ADAPTATION::Linear_Offset(std_array<int,d>(index_1));
                                                    unsigned long flag_offset_2=T_MASK_ADAPTATION::Linear_Offset(std_array<int,d>(index_2));
                                                    unsigned long index_offset_1=T_INDEX_MASK_ADAPTATION::Linear_Offset(std_array<int,d>(index_1));
                                                    unsigned long index_offset_2=T_INDEX_MASK_ADAPTATION::Linear_Offset(std_array<int,d>(index_2));
                                                    int index_id_1;
                                                    int index_id_2;
                                                    unsigned flag1=flags(flag_offset_1);
                                                    unsigned flag2=flags(flag_offset_2);
                                                    if(flag1&Adaptive_Node_Degree_Marker)
                                                        index_id_1=index(index_offset_1);
                                                    else{
                                                        int current_level=cur_level;
                                                        T_INDEX index=index_1;
                                                        while(!(flag1&Adaptive_Node_Degree_Marker)){
                                                            ++current_level;
                                                            index/=2;
                                                            //We need to go a level coarser
                                                            SPG_Index_Array_Type_Adaptation index_coarse=(hierarchy[mg_level-1]->Allocator(current_level).Get_Array(&T_STRUCT_ADAPTATION::node_index));
                                                            SPG_Flags_Array_Type_Adaptation flags_coarse=hierarchy[mg_level-1]->Allocator(current_level).Get_Array(&T_STRUCT_ADAPTATION::flags);
                                                            index_id_1=index_coarse(T_INDEX_MASK_ADAPTATION::Linear_Offset(std_array<int,d>(index)));
                                                            flag1=flags_coarse(T_INDEX_MASK_ADAPTATION::Linear_Offset(std_array<int,d>(index)));}}

                                                    if(flag2&Adaptive_Node_Degree_Marker)
                                                        index_id_2=index(index_offset_2);
                                                    else{
                                                        int current_level=cur_level;
                                                        T_INDEX index=index_2;
                                                        while(!(flag2&Adaptive_Node_Degree_Marker)){
                                                            ++current_level;
                                                            index/=2;
                                                            //We need to go a level coarser
                                                            SPG_Index_Array_Type_Adaptation index_coarse=(hierarchy[mg_level-1]->Allocator(current_level).Get_Array(&T_STRUCT_ADAPTATION::node_index));
                                                            SPG_Flags_Array_Type_Adaptation flags_coarse=hierarchy[mg_level-1]->Allocator(current_level).Get_Array(&T_STRUCT_ADAPTATION::flags);
                                                            index_id_2=index_coarse(T_INDEX_MASK_ADAPTATION::Linear_Offset(std_array<int,d>(index)));
                                                            flag2=flags_coarse(T_INDEX_MASK_ADAPTATION::Linear_Offset(std_array<int,d>(index))); 
                                                        }                                                
                                                    }
                                                    
                                                    if((flag1&Uniform_Node_Type_Interface)&&(flag2&Uniform_Node_Type_Interface)){
                                                        //PHYSBAM_ASSERT(index_id_2==index(index_offset_2));
                                                        //PHYSBAM_ASSERT(index_id_1==index(index_offset_1));
                                                        Arr[mg_level-1][index_id_1]->Get_Or_Insert(VECTOR<int,1>(index_id_2))+=stencil_1_iterator.Data()*
                                                            stencil_2_iterator.Data()*matrix_entry;}
                                                    if((flag1&(Adaptive_Node_Type_Active|Adaptive_Node_Type_Coarse_Shared))&&
                                                       (flag2&(Adaptive_Node_Type_Active|Adaptive_Node_Type_Coarse_Shared))){
                                                        Aii_subdomain[mg_level-1][index_id_1].Get_Or_Insert(VECTOR<int,1>(index_id_2))+=stencil_1_iterator.Data()*
                                                            stencil_2_iterator.Data()*matrix_entry;}
                                                    if((flag1&(Adaptive_Node_Type_Active|Adaptive_Node_Type_Coarse_Shared))&&
                                                       (flag2&Uniform_Node_Type_Interface))
                                                        Air_subdomain[mg_level-1][index_id_1].Get_Or_Insert(VECTOR<int,1>(index_id_2))+=stencil_1_iterator.Data()*
                                                            stencil_2_iterator.Data()*matrix_entry;
                                            
                                                    PHYSBAM_ASSERT(flag1&Adaptive_Node_Degree_Marker);
                                                    PHYSBAM_ASSERT(flag2&Adaptive_Node_Degree_Marker);
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                            T_INDEX cell_parity;
                            RANGE<T_INDEX> tmp_range=cur_range;
                            for(int axis=1;axis<=d;++axis) cur_range.min_corner(axis)=(cur_range.min_corner(axis)>>cur_level)<<cur_level;
                            for(int axis=1;axis<=d;++axis) cell_parity(axis)=(tmp_range.min_corner(axis)-cur_range.min_corner(axis)==0)?0:1;
                            cur_range.max_corner=cur_range.min_corner+(1<<cur_level);
                            cur_index=(cur_index+1)/2;
                            cur_matrix=prolongation_baker.prolongation_matrices[prolongation_baker.index_map(cell_parity)].transpose()*cur_matrix*prolongation_baker.prolongation_matrices[prolongation_baker.index_map(cell_parity)];
                            ++cur_level;
                        }while(cur_level<=levels);
                    }
                }
            }
        }
        }
#else
        {LOG::SCOPE scope("GALERKIN PROCESS");
        #pragma omp parallel for    
        for(int level=1;level<=mg_levels;++level){
            T_INDEX adaptation_block_size=hierarchy[level-1]->Allocator(level).Block_Size().template Cast<T_INDEX>();
            for(int hierarchy_level=level;hierarchy_level<=levels;++hierarchy_level){
                SPG_Flags_Array_Type_Adaptation flags=hierarchy[level-1]->Allocator(hierarchy_level).Get_Array(&T_STRUCT_ADAPTATION::flags);
                SPG_Index_Array_Type_Adaptation index=hierarchy[level-1]->Allocator(hierarchy_level).Get_Array(&T_STRUCT_ADAPTATION::node_index); 
                for(SPGrid_Block_Iterator<T_MASK_ADAPTATION> iterator(hierarchy[level-1]->Set(hierarchy_level).Get_Blocks());iterator.Valid();iterator.Next_Block()){
                    unsigned long offset=iterator.Offset();
                    T_INDEX base_index=iterator.Index().template Cast<T_INDEX>();
                    for(RANGE_ITERATOR<d> iterator(RANGE<T_INDEX>(base_index,base_index+adaptation_block_size-1));
                        iterator.Valid();
                        iterator.Next(),offset+=sizeof(unsigned)){
                        const unsigned& flag=flags(offset);
                        const T_INDEX& cell_index=iterator.Index();
                        if(flag&Adaptive_Cell_Type_Interior){
                            Eigen::Matrix<T,1<<d,1<<d> cell_matrix=GALERKIN_PROCESS<T_STRUCT_SOLVER,T,d>::Galerkin_Coarsen(set,prolongation_baker.prolongation_matrices,prolongation_baker.index_map,cell_index,hierarchy_level);
                            std::vector<STENCIL<T,d> > parents(2<<d);
                            for(RANGE_ITERATOR<d> node_iterator(RANGE<T_INDEX>::Unit_Box());node_iterator.Valid();node_iterator.Next()){
                                T_INDEX node_index=cell_index-1+node_iterator.Index();
                                Embedding_Helper(node_index,parents[prolongation_baker.index_map(node_iterator.Index())],(T)1.,q(hierarchy_level));}
                            //We always compute the embedding map. If it is a real DOF, the embedding map would just be 1.0 of itself
                            for(RANGE_ITERATOR<d> node_1_iterator(RANGE<T_INDEX>::Unit_Box());node_1_iterator.Valid();node_1_iterator.Next()){
                                for(RANGE_ITERATOR<d> node_2_iterator(RANGE<T_INDEX>::Unit_Box());node_2_iterator.Valid();node_2_iterator.Next()){
                                    const int node_1_local_index=prolongation_baker.index_map(node_1_iterator.Index());
                                    const int node_2_local_index=prolongation_baker.index_map(node_2_iterator.Index());
                                    T matrix_entry=cell_matrix(node_1_local_index,node_2_local_index);
                                    if(matrix_entry==0) continue;
                                    for(STENCIL_ITERATOR<T,d> stencil_1_iterator(parents[node_1_local_index]);
                                        stencil_1_iterator.Valid();stencil_1_iterator.Next()){
                                        for(STENCIL_ITERATOR<T,d> stencil_2_iterator(parents[node_2_local_index]);
                                            stencil_2_iterator.Valid();stencil_2_iterator.Next()){
                                            T_INDEX index_1=stencil_1_iterator.Key();
                                            T_INDEX index_2=stencil_2_iterator.Key();
                                            unsigned long flag_offset_1=T_MASK_ADAPTATION::Linear_Offset(std_array<int,d>(index_1));
                                            unsigned long flag_offset_2=T_MASK_ADAPTATION::Linear_Offset(std_array<int,d>(index_2));
                                            unsigned long index_offset_1=T_INDEX_MASK_ADAPTATION::Linear_Offset(std_array<int,d>(index_1));
                                            unsigned long index_offset_2=T_INDEX_MASK_ADAPTATION::Linear_Offset(std_array<int,d>(index_2));
                                            int index_id_1;
                                            int index_id_2;
                                            unsigned flag1=flags(flag_offset_1);
                                            unsigned flag2=flags(flag_offset_2);
                                            //Sometime at finest level, a node is not a DOF but when coarsened, it is activated. That's why the Adaptive_Node_Degree_Marker is not always set.
                                            //QUESTION: Is going up one level enough? The fact that this T_JUNCTION is a real node in the original matrix, does that automatically make the parents of real in the next level?
                                            if(flag1&Adaptive_Node_Degree_Marker)
                                                index_id_1=index(index_offset_1);
                                            else{
                                                int current_level=hierarchy_level;
                                                T_INDEX index=index_1;
                                                while(!(flag1&Adaptive_Node_Degree_Marker)){
                                                    ++current_level;
                                                    index/=2;
                                                    //We need to go a level coarser
                                                    SPG_Index_Array_Type_Adaptation index_coarse=(hierarchy[level-1]->Allocator(current_level).Get_Array(&T_STRUCT_ADAPTATION::node_index));
                                                    SPG_Flags_Array_Type_Adaptation flags_coarse=hierarchy[level-1]->Allocator(current_level).Get_Array(&T_STRUCT_ADAPTATION::flags);
                                                    index_id_1=index_coarse(T_INDEX_MASK_ADAPTATION::Linear_Offset(std_array<int,d>(index)));
                                                    flag1=flags_coarse(T_INDEX_MASK_ADAPTATION::Linear_Offset(std_array<int,d>(index))); 
                                                }
                                            }

                                            if(flag2&Adaptive_Node_Degree_Marker)
                                                index_id_2=index(index_offset_2);
                                            else{
                                                int current_level=hierarchy_level;
                                                T_INDEX index=index_2;
                                                while(!(flag2&Adaptive_Node_Degree_Marker)){
                                                    ++current_level;
                                                    index/=2;
                                                    //We need to go a level coarser
                                                    SPG_Index_Array_Type_Adaptation index_coarse=(hierarchy[level-1]->Allocator(current_level).Get_Array(&T_STRUCT_ADAPTATION::node_index));
                                                    SPG_Flags_Array_Type_Adaptation flags_coarse=hierarchy[level-1]->Allocator(current_level).Get_Array(&T_STRUCT_ADAPTATION::flags);
                                                    index_id_2=index_coarse(T_INDEX_MASK_ADAPTATION::Linear_Offset(std_array<int,d>(index)));
                                                    flag2=flags_coarse(T_INDEX_MASK_ADAPTATION::Linear_Offset(std_array<int,d>(index))); 
                                                }                                                
                                            }

                                            if((flag1&Uniform_Node_Type_Interface)&&(flag2&Uniform_Node_Type_Interface)){
                                                //PHYSBAM_ASSERT(index_id_2==index(index_offset_2));
                                                //PHYSBAM_ASSERT(index_id_1==index(index_offset_1));
                                                Arr[level-1][index_id_1]->Get_Or_Insert(VECTOR<int,1>(index_id_2))+=stencil_1_iterator.Data()*
                                                    stencil_2_iterator.Data()*matrix_entry;}
                                            if((flag1&(Adaptive_Node_Type_Active|Adaptive_Node_Type_Coarse_Shared))&&
                                               (flag2&(Adaptive_Node_Type_Active|Adaptive_Node_Type_Coarse_Shared))){
                                                Aii_subdomain[level-1][index_id_1].Get_Or_Insert(VECTOR<int,1>(index_id_2))+=stencil_1_iterator.Data()*
                                                    stencil_2_iterator.Data()*matrix_entry;}
                                            if((flag1&(Adaptive_Node_Type_Active|Adaptive_Node_Type_Coarse_Shared))&&
                                               (flag2&Uniform_Node_Type_Interface))
                                                Air_subdomain[level-1][index_id_1].Get_Or_Insert(VECTOR<int,1>(index_id_2))+=stencil_1_iterator.Data()*
                                                    stencil_2_iterator.Data()*matrix_entry;
                                            
                                            PHYSBAM_ASSERT(flag1&Adaptive_Node_Degree_Marker);
                                            PHYSBAM_ASSERT(flag2&Adaptive_Node_Degree_Marker);                                            
                                        }
                                    }                                    
                                }
                            }
                        }
                    }
                }
            }
        }
        }
#endif
        //LOG::cout<<"interface_dof_counter: "<<interface_dof_counter[0]<<std::endl;
        //LOG::cout<<"interface_dof_counter: "<<interface_hash[0].Size()<<std::endl;
    }
    void Convert_To_Eigen(INTERFACE_SOLVER_DATA<T,d>& interface_solver_data){
        const int number_of_subdomains=Aii.size();
        std::cout<<"number_of_subdomains: "<<number_of_subdomains<<std::endl;
        interface_solver_data.Resize(mg_levels,number_of_subdomains);
        //Arr
        for(int level=0;level<mg_levels;++level){
            std::vector<TRIPLET> triplet_list;
            const int interface_dof=interface_dof_counter[level];
            for(int i=0;i<interface_dof;++i){
                for(STENCIL_ITERATOR<T,1> stencil_iterator(*Arr[level][i]);stencil_iterator.Valid();stencil_iterator.Next()){
                    //PHYSBAM_ASSERT(stencil_iterator.Key()(1)<interface_dof);
                    //if(stencil_iterator.Data()==2.5) 
                        //std::cout<<interface_
                    if(stencil_iterator.Data()!=0)
                        triplet_list.push_back(TRIPLET(i,stencil_iterator.Key()(1),stencil_iterator.Data()));}}            
            interface_solver_data.Arr[level].resize(interface_dof,interface_dof);
            interface_solver_data.Arr[level].setFromTriplets(triplet_list.begin(),triplet_list.end());
            interface_solver_data.Arr[level].makeCompressed();}
        
        for(int j=0;j<interface_dof_counter[0];++j)
            PHYSBAM_ASSERT(interface_solver_data.Arr[0].coeff(j,j)!=0);
        //Clear Memory
        for(int i=0;i<Arr.size();++i)
            for(int j=0;j<Arr[i].size();++j)
                delete Arr[i][j];
        Arr.clear();

        std::vector<std::vector<HASHTABLE<int> > > interior_free_dof_hash(number_of_subdomains);        
        //Aii
        for(int i=0;i<number_of_subdomains;++i){
            interior_free_dof_hash[i].resize(mg_levels);
            for(int level=0;level<mg_levels;++level){
                std::vector<TRIPLET> triplet_list;
                const int interior_dof=(*Aii[i])[level].size();
                for(int j=0;j<interior_dof;++j){
                    for(STENCIL_ITERATOR<T,1> stencil_iterator((*Aii[i])[level][j]);stencil_iterator.Valid();stencil_iterator.Next()){
                        if(stencil_iterator.Data()!=0)
                            triplet_list.push_back(TRIPLET(j,stencil_iterator.Key()(1),stencil_iterator.Data()));}}
                interface_solver_data.Aii[level][i].resize(interior_dof,interior_dof);
                interface_solver_data.Aii[level][i].setFromTriplets(triplet_list.begin(),triplet_list.end());
                interface_solver_data.Aii[level][i].makeCompressed();
            
                for(int j=0;j<(*Aii[i])[level].size();++j)
                    if(interface_solver_data.Aii[level][i].coeff(j,j)==0)
                        interior_free_dof_hash[i][level].Insert(j);
            }
        }

        //Clear Memory
        for(int i=0;i<Aii.size();++i)
            delete Aii[i];
        Aii.clear();
        //Air
        for(int i=0;i<number_of_subdomains;++i){
            for(int level=0;level<mg_levels;++level){
                std::vector<TRIPLET> triplet_list;
                const int interface_dof=interface_dof_counter[level];
                const int interior_dof=(*Air[i])[level].size();
                for(int j=0;j<interior_dof;++j){
                    for(STENCIL_ITERATOR<T,1> stencil_iterator((*Air[i])[level][j]);stencil_iterator.Valid();stencil_iterator.Next()){
                        PHYSBAM_ASSERT(stencil_iterator.Key()(1)<interface_dof);
                        if(stencil_iterator.Data()!=0){
                            PHYSBAM_ASSERT(!interior_free_dof_hash[i][level].Contains(j)); 
                            triplet_list.push_back(TRIPLET(j,stencil_iterator.Key()(1),stencil_iterator.Data()));}}}
                interface_solver_data.Air[level][i].resize(interior_dof,interface_dof);
                interface_solver_data.Air[level][i].setFromTriplets(triplet_list.begin(),triplet_list.end());
                interface_solver_data.Air[level][i].makeCompressed();}}
        //Clear Memory
        for(int i=0;i<Air.size();++i)
            delete Air[i];
        Air.clear();
        //Prr
        for(int level=0;level<mg_levels-1;++level){
            PROLONGATION_MATRIX_HELPER<T,d>::Construct_Prolongation_Matrix(interface_solver_data.Prr[level],interface_hash[level],interface_hash[level+1]);};
    }
};
#endif
