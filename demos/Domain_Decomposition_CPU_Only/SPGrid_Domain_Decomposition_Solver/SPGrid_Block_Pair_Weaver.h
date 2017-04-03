//#####################################################################
//Removed For Anonymity: Copyright Authors of Submission pap173s1
//Removed For Anonymity: Copyright Authors of Submission pap173s1
//#####################################################################
// Subroutine SPGrid_Computations::Block_Pair_Weaver
//#####################################################################
#ifndef __SPGrid_Block_Pair_Weaver_h__
#define __SPGrid_Block_Pair_Weaver_h__
#include <Common_Tools/Math_Tools/RANGE_ITERATOR.h>
#include <SPGrid/Core/SPGrid_Allocator.h>
#include <SPGrid_Fluids/Grids/GRID_TOPOLOGY_HELPER.h>
#include <SPGrid_Fluids/Solvers/Domain_Decomposition/SPGRID_MULTIGRID_FLAGS.h>

namespace SPGrid_Computations{

using namespace SPGrid;
using namespace PhysBAM;

template<typename T_STRUCT,int d>
class Block_Pair_Weaver
{
    //This is not particularly thread safe....so....
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<unsigned>::type Flag_Array_Type;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<const unsigned>::type Const_Flag_Array_Type;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<unsigned>::mask T_MASK;
    typedef VECTOR<int,d> T_INDEX;    

public:    
    static void Weave_Subdomain(std::vector<unsigned long>& blocks_output,SPGrid_Allocator<T_STRUCT,d>& allocator,const std::pair<const unsigned long*,unsigned>& blocks_in,unsigned T_STRUCT::* flags_field,int subdomain_id){
        //OK, this function assumes there are enough paddings so that testing neighors of the dirichlet would not reach out to the void.
        static const int number_of_face_neighbors=GRID_TOPOLOGY_HELPER<T_MASK>::faces_per_cell;
        unsigned long face_neighbor_offsets[number_of_face_neighbors];
        GRID_TOPOLOGY_HELPER<T_MASK>::Face_Neighbor_Offsets(face_neighbor_offsets); // order is (-x, x, -y, y, -z, z)

        PHYSBAM_ASSERT(subdomain_id>0);
        blocks_output.clear();
        Const_Flag_Array_Type flags=allocator.Get_Const_Array(flags_field);    
        for(SPGrid_Block_Iterator<T_MASK> iterator(blocks_in);iterator.Valid();iterator.Next_Block()){
            unsigned long offset=iterator.Offset();
            unsigned long page_offset=offset;
            for(RANGE_ITERATOR<d> iterator(RANGE<T_INDEX>(T_INDEX(),allocator.Block_Size().template Cast<T_INDEX>()-1));
                iterator.Valid();
                iterator.Next(),offset+=sizeof(unsigned)){                
                const unsigned& flag=flags(offset);
                if((flags(offset)&SPGrid_Solver_PartitionID_Mask)==subdomain_id){
                    blocks_output.push_back(page_offset);
                    break;
                }else if((flags(offset)&SPGrid_Solver_Cell_Type_Interface)/*||(flags(offset)&SPGrid_Solver_Cell_Type_Dirichlet)*/){
                    //if it is a interface block, check if adjecent cells belongs to that subdomain. If so, push this block to the list.
                    //ok, if it is a dirichlet cell, we will just ignore it. because it is going to be zero dirichlet condition anyway. and the accelerator does not check dirichlet flags anyway.
                    bool belongs=false;
                    for(int face=0;face<number_of_face_neighbors;face++){
                        unsigned long neighbor_offset=T_MASK::Packed_Add(offset,face_neighbor_offsets[face]);
                        if((flags(neighbor_offset)&SPGrid_Solver_PartitionID_Mask)==subdomain_id){
                            belongs=true;
                            break;}}
                    if(belongs){
                        blocks_output.push_back(page_offset);
                        break;}}}}
    }

    static void Weave_Subdomain(std::vector<unsigned long>& blocks_output,SPGrid_Allocator<T_STRUCT,d>& allocator,const std::pair<const unsigned long*,unsigned>& blocks_in,unsigned T_STRUCT::* flags_field){
        //OK, this function assumes there are enough paddings so that testing neighors of the dirichlet would not reach out to the void.
        static const int number_of_face_neighbors=GRID_TOPOLOGY_HELPER<T_MASK>::faces_per_cell;
        unsigned long face_neighbor_offsets[number_of_face_neighbors];
        GRID_TOPOLOGY_HELPER<T_MASK>::Face_Neighbor_Offsets(face_neighbor_offsets); // order is (-x, x, -y, y, -z, z)
        blocks_output.clear();
        Const_Flag_Array_Type flags=allocator.Get_Const_Array(flags_field);    
        for(SPGrid_Block_Iterator<T_MASK> iterator(blocks_in);iterator.Valid();iterator.Next_Block()){
            unsigned long offset=iterator.Offset();
            unsigned long page_offset=offset;
            for(RANGE_ITERATOR<d> iterator(RANGE<T_INDEX>(T_INDEX(),allocator.Block_Size().template Cast<T_INDEX>()-1));
                iterator.Valid();
                iterator.Next(),offset+=sizeof(unsigned)){                
                const unsigned& flag=flags(offset);
                if((flags(offset)&SPGrid_Solver_Cell_Type_Active)||(flags(offset)&SPGrid_Solver_Cell_Type_Interface)){
                    blocks_output.push_back(page_offset);
                    break;}}}
    }

    static void Weave_Subdomain(std::vector<unsigned long>& blocks_output,SPGrid_Allocator<T_STRUCT,d>& allocator,const std::pair<const unsigned long*,unsigned>& blocks_in,unsigned T_STRUCT::* flags_field,RANGE<T_INDEX> subdomain_range/*include interface*/){
        //This function assumes subdomain is a box.....
        static const int number_of_face_neighbors=GRID_TOPOLOGY_HELPER<T_MASK>::faces_per_cell;
        unsigned long face_neighbor_offsets[number_of_face_neighbors];
        GRID_TOPOLOGY_HELPER<T_MASK>::Face_Neighbor_Offsets(face_neighbor_offsets); // order is (-x, x, -y, y, -z, z)
        blocks_output.clear();
        Const_Flag_Array_Type flags=allocator.Get_Const_Array(flags_field);    
        for(SPGrid_Block_Iterator<T_MASK> iterator(blocks_in);iterator.Valid();iterator.Next_Block()){
            unsigned long page_offset=iterator.Offset();
            T_INDEX base_index=iterator.Index().template Cast<T_INDEX>();
            RANGE<T_INDEX> subdomain_range(base_index,base_index+allocator.Block_Size().template Cast<T_INDEX>()-1);
            if(subdomain_range.Intersection(subdomain_range)) blocks_output.push_back(page_offset);}
    }

    static void Weave_With_Interface_Flags(std::vector<unsigned long>& blocks_output,SPGrid_Allocator<T_STRUCT,d>& allocator,const std::pair<const unsigned long*,unsigned>& blocks_in,unsigned T_STRUCT::* flags_field){
        //This returns all blocks that contain interface
        blocks_output.clear();
        Const_Flag_Array_Type flags=allocator.Get_Const_Array(flags_field);    
        for(SPGrid_Block_Iterator<T_MASK> iterator(blocks_in);iterator.Valid();iterator.Next_Block()){
            unsigned long offset=iterator.Offset();
            unsigned long page_offset=offset;
            for(RANGE_ITERATOR<d> iterator(RANGE<T_INDEX>(T_INDEX(),allocator.Block_Size().template Cast<T_INDEX>()-1));
                iterator.Valid();
                iterator.Next(),offset+=sizeof(unsigned)){
                const unsigned& flag=flags(offset);
                if(flags(offset)&SPGrid_Solver_Cell_Type_Interface){
                    blocks_output.push_back(page_offset);
                    break;}}}
    }

    static void Weave_With_Boundary_Flags(std::vector<unsigned long>& blocks_output,SPGrid_Allocator<T_STRUCT,d>& allocator,const std::pair<const unsigned long*,unsigned>& blocks_in,unsigned T_STRUCT::* flags_field,int subdomain_id){
        PHYSBAM_ASSERT(subdomain_id>0);
        blocks_output.clear();
        Const_Flag_Array_Type flags=allocator.Get_Const_Array(flags_field);    
        for(SPGrid_Block_Iterator<T_MASK> iterator(blocks_in);iterator.Valid();iterator.Next_Block()){
            unsigned long offset=iterator.Offset();
            unsigned long page_offset=offset;
            for(RANGE_ITERATOR<d> iterator(RANGE<T_INDEX>(T_INDEX(),allocator.Block_Size().template Cast<T_INDEX>()-1));
                iterator.Valid();
                iterator.Next(),offset+=sizeof(unsigned)){
                const unsigned& flag=flags(offset);
                if(((flags(offset)&SPGrid_Solver_PartitionID_Mask)==subdomain_id)&&(flags(offset)&SPGrid_Solver_Cell_Type_Boundary)){
                    blocks_output.push_back(page_offset);
                    break;}}}
    }
    static void Weave_With_Boundary_Flags(std::vector<unsigned long>& blocks_output,SPGrid_Allocator<T_STRUCT,d>& allocator,const std::pair<const unsigned long*,unsigned>& blocks_in,unsigned T_STRUCT::* flags_field){
        blocks_output.clear();
        Const_Flag_Array_Type flags=allocator.Get_Const_Array(flags_field);    
        for(SPGrid_Block_Iterator<T_MASK> iterator(blocks_in);iterator.Valid();iterator.Next_Block()){
            unsigned long offset=iterator.Offset();
            unsigned long page_offset=offset;
            for(RANGE_ITERATOR<d> iterator(RANGE<T_INDEX>(T_INDEX(),allocator.Block_Size().template Cast<T_INDEX>()-1));
                iterator.Valid();
                iterator.Next(),offset+=sizeof(unsigned)){
                const unsigned& flag=flags(offset);
                if(flags(offset)&SPGrid_Solver_Cell_Type_Boundary){
                    blocks_output.push_back(page_offset);
                    break;}}}
    }
};
}
#endif
