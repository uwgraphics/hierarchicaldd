//#####################################################################
//Removed For Anonymity: Copyright Authors of Submission pap173s1
// Distributed under the FreeBSD license (see license.txt)
//#####################################################################
// Class GRID_HIERARCHY_INITIALIZER
//#####################################################################
#include <Common_Tools/Math_Tools/RANGE_ITERATOR.h>
#include <SPGrid/Tools/SPGrid_Set_Iterator.h>
#include <SPGrid_Fluids/Grids/GRID_HIERARCHY_INITIALIZER.h>
#include <SPGrid_Fluids/Grids/GRID_TOPOLOGY_HELPER.h>
#include <SPGrid_Fluids/Simulation/FLUIDS_SIMULATION_FLAGS.h>
#include "SPGRID_MULTIGRID_DATA.h"

using namespace PhysBAM;
//#####################################################################
// Class Parity_Helper (TO BE DELETED!!!!)
//#####################################################################
namespace{

template<class MASK,int d> class Parity_Helper;

template<class MASK>
class Parity_Helper<MASK,2>
{
public:
    enum {
        GHOST_000 = MASK::template LinearOffset<0,0>::value,
        GHOST_010 = MASK::template LinearOffset<0,1>::value,
        GHOST_100 = MASK::template LinearOffset<1,0>::value,
        GHOST_110 = MASK::template LinearOffset<1,1>::value,
        
        GHOST_001 = MASK::template LinearOffset<1,1>::value,
        GHOST_011 = MASK::template LinearOffset<1,1>::value,
        GHOST_101 = MASK::template LinearOffset<1,1>::value,
        GHOST_111 = MASK::template LinearOffset<1,1>::value,

        X_EDGE = MASK::template LinearOffset<0,1>::value, 
        Y_EDGE = MASK::template LinearOffset<1,0>::value,

        X = 1,
        Y = 2,
        CELL = 0,
    };
};

template<class MASK>
class Parity_Helper<MASK,3>
{
public:
    enum {
        GHOST_000 = MASK::template LinearOffset<0,0,0>::value,
        GHOST_001 = MASK::template LinearOffset<0,0,1>::value,
        GHOST_010 = MASK::template LinearOffset<0,1,0>::value,
        GHOST_011 = MASK::template LinearOffset<0,1,1>::value,
        GHOST_100 = MASK::template LinearOffset<1,0,0>::value,
        GHOST_101 = MASK::template LinearOffset<1,0,1>::value,
        GHOST_110 = MASK::template LinearOffset<1,1,0>::value,
        GHOST_111 = MASK::template LinearOffset<1,1,1>::value,
        
        X_EDGE = MASK::template LinearOffset<0,1,1>::value, 
        Y_EDGE = MASK::template LinearOffset<1,0,1>::value, 
        Z_EDGE = MASK::template LinearOffset<1,1,0>::value,
         
        X_FACE = MASK::template LinearOffset<1,0,0>::value, 
        Y_FACE = MASK::template LinearOffset<0,1,0>::value, 
        Z_FACE = MASK::template LinearOffset<0,0,1>::value, 
        
        X = 1,
        Y = 2,
        Z = 3,
        CELL = 0,
    };
};

}
//#####################################################################
// Function Flag_Ghost_Cells
//#####################################################################
template<class T_STRUCT,class T,int d> void GRID_HIERARCHY_INITIALIZER<T_STRUCT,T,d>::
Flag_Ghost_Cells(Hierarchy_type& hierarchy)
{
    LOG::SCOPE scope("GRID_HIERARCHY_INITIALIZER::Flag_Ghost_Cells");

    typedef Parity_Helper<Flag_array_mask,d> Parity;
    //LOG::cout<<"Tagging ghost cells"<<std::endl;
    static const int number_of_face_neighbors=(d==2)?4:6;
    unsigned long face_neighbor_shifts[number_of_face_neighbors];
    const int levels = hierarchy.Levels();

    for(int axis=0,face_neighbor=0;axis<d;axis++)
        for(int axis_shift=-1;axis_shift<=1;axis_shift+=2){

            std_array<int,d> shift;
            shift.data[axis]=axis_shift;
            face_neighbor_shifts[face_neighbor++]=Flag_array_mask::Linear_Offset(shift);
        }

    for(int level=1;level<levels;level++){

        for(SPGrid_Set_Iterator<Flag_array_type> iterator(hierarchy.Set(level));iterator.Valid();iterator.Next())
        {
            if(iterator.Data() & SPGrid_Cell_Type_Interior)
            {
                for(int face_neighbor=0;face_neighbor<number_of_face_neighbors;face_neighbor++)
                {
                    unsigned long neighbor_offset=iterator.Offset(face_neighbor_shifts[face_neighbor]);
                    unsigned long cur_offset=neighbor_offset;
                    // Now search at this and coarser resolutions for parent
                    int parent_level = 0;
                    for (int cur_level=level;cur_level<=levels;cur_level++)
                    {
                        if(hierarchy.Set(cur_level).Is_Set(cur_offset,(SPGrid_Cell_Type_Interior|SPGrid_Cell_Type_Ghost)))
                        {
                            parent_level = cur_level;
                            break;
                        }
                        cur_offset = Flag_array_mask::DownsampleOffset(cur_offset);
                            
                    }

                    // Fill in ghost values (if parent level==0, nothing to be done!)
                    cur_offset=neighbor_offset;
                    for (int cur_level=level;cur_level<parent_level;cur_level++)
                    {
                        unsigned long parent_offset = Flag_array_mask::DownsampleOffset(cur_offset);
                        hierarchy.Activate_Cell(cur_level,cur_offset,SPGrid_Cell_Type_Ghost);
                            
                        unsigned mask;
                        // Mask bits in parent mask

                        // TODO: Fix, make truly dimension-nonspecific

                        unsigned pmask = cur_offset & Parity::GHOST_111; // Rename to GHOST_ALL_ONES
                        if(pmask == Parity::GHOST_000)
                            mask = SPGrid_Ghost_Child_000;
                        else if(pmask == Parity::GHOST_010)
                            mask = SPGrid_Ghost_Child_010;
                        else if(pmask == Parity::GHOST_100)
                            mask = SPGrid_Ghost_Child_100;
                        else if(pmask == Parity::GHOST_110)
                            mask = SPGrid_Ghost_Child_110;
                        else if(pmask == Parity::GHOST_001)
                        {
                            PHYSBAM_ASSERT(d==3);
                            mask = SPGrid_Ghost_Child_001;
                        }
                        else if(pmask == Parity::GHOST_011)
                        {
                            PHYSBAM_ASSERT(d==3);
                            mask = SPGrid_Ghost_Child_011;
                        }
                        else if(pmask == Parity::GHOST_101)
                        {
                            PHYSBAM_ASSERT(d==3);
                            mask = SPGrid_Ghost_Child_101;
                        }
                        else if(pmask == Parity::GHOST_111)
                        {
                            PHYSBAM_ASSERT(d==3);
                            mask = SPGrid_Ghost_Child_111;
                        }

                        hierarchy.Activate_Cell(cur_level+1,parent_offset,mask);
                        cur_offset = parent_offset;
                    }
                }
            }
        }
    }
    //LOG::cout<<"DONE: Tagging ghost cells"<<std::endl;
}
    
//#####################################################################
// Function Flag_Active_Faces
//#####################################################################
template<class T_STRUCT,class T,int d> void GRID_HIERARCHY_INITIALIZER<T_STRUCT,T,d>::
Flag_Active_Faces(Hierarchy_type& hierarchy)
{
	LOG::SCOPE scope("GRID_HIERARCHY_INITIALIZER::Flag_Active_Faces");

    //LOG::cout<<"Tagging active faces"<<std::endl;  
#if 0
    static const int number_of_face_neighbors=(d==2)?4:6;
    unsigned long face_neighbor_shifts[number_of_face_neighbors];
    const int levels = hierarchy.Levels();

    for(int axis=0,face_neighbor=0;axis<d;axis++)
        for(int axis_shift=-1;axis_shift<=1;axis_shift+=2){

            std_array<int,d> shift;
            shift.data[axis]=axis_shift;
            face_neighbor_shifts[face_neighbor++]=Flag_array_mask::Linear_Offset(shift);
        }

    // Tag edges
    for(int level=1;level<=levels;level++)
        for(SPGrid_Set_Iterator<Flag_array_type> iterator(hierarchy.Set(level));iterator.Valid();iterator.Next())
            if(iterator.Data() & SPGrid_Cell_Type_Interior)
                for(int face_neighbor=0;face_neighbor<number_of_face_neighbors;face_neighbor++)
                {
                    unsigned long neighbor_offset=iterator.Offset(face_neighbor_shifts[face_neighbor]);

                    if(hierarchy.Set(level).Is_Set(neighbor_offset,SPGrid_Cell_Type_Interior|SPGrid_Cell_Type_Dirichlet))
                        iterator.Data() |= (SPGrid_Face_Minus_X_Active<<face_neighbor);
                    else if(hierarchy.Set(level).Is_Set(neighbor_offset,SPGrid_Cell_Type_Ghost))
                        iterator.Data() |= (SPGrid_Face_Minus_X_Scaled|SPGrid_Face_Minus_X_Active)<<face_neighbor;
                }
            else if(iterator.Data() & SPGrid_Cell_Type_Ghost)
                for(int face_neighbor=0;face_neighbor<number_of_face_neighbors;face_neighbor++)
                {
                    unsigned long neighbor_offset=iterator.Offset(face_neighbor_shifts[face_neighbor]);

                    if(hierarchy.Set(level).Is_Set(neighbor_offset,(SPGrid_Cell_Type_Interior|SPGrid_Cell_Type_Dirichlet)))
                        iterator.Data() |= (SPGrid_Face_Minus_X_Scaled|SPGrid_Face_Minus_X_Active)<<face_neighbor;
                }
#else
    const int levels = hierarchy.Levels();

    VECTOR<unsigned long,d> left_cell_offsets;
    for(int axis=1;axis<=d;axis++)
        left_cell_offsets(axis)=Flag_array_mask::Linear_Offset(std_array<int,d>(T_INDEX::Axis_Vector(axis)*(T)(-1.)));

    for(int level=1;level<=levels;level++)
        for(SPGrid_Set_Iterator<Flag_array_type> iterator(hierarchy.Set(level));iterator.Valid();iterator.Next())
            for(int axis=1;axis<=d;axis++){
                unsigned face_valid_mask=GRID_TOPOLOGY_HELPER<Flag_array_mask>::Face_Valid_Mask(axis);
                if(iterator.Data() & face_valid_mask){
                    bool left_cell_not_exterior,right_cell_not_exterior;
                    unsigned long left_cell_offset=Flag_array_mask::Packed_Add(iterator.Offset(),left_cell_offsets(axis));
                    left_cell_not_exterior=((hierarchy.Set(level).Is_Set(left_cell_offset,SPGrid_Cell_Type_Interior)) ||
                                            (hierarchy.Set(level).Is_Set(left_cell_offset,SPGrid_Cell_Type_Dirichlet)) ||
                                            (hierarchy.Set(level).Is_Set(left_cell_offset,SPGrid_Cell_Type_Ghost)));
                    right_cell_not_exterior=((iterator.Data() & SPGrid_Cell_Type_Interior) || 
                                             (iterator.Data() & SPGrid_Cell_Type_Dirichlet) || 
                                             (iterator.Data() & SPGrid_Cell_Type_Ghost));
                    if(left_cell_not_exterior && right_cell_not_exterior)
                        iterator.Data() |= GRID_TOPOLOGY_HELPER<Flag_array_mask>::Face_Active_Mask(axis);
                }}
#endif
    //LOG::cout<<"DONE: Tagging active faces"<<std::endl;
}

//#####################################################################
// Function Flag_Valid_Faces
//#####################################################################
template<class T_STRUCT,class T,int d> void GRID_HIERARCHY_INITIALIZER<T_STRUCT,T,d>::
Flag_Valid_Faces(Hierarchy_type& hierarchy)
{
    LOG::SCOPE scope("GRID_HIERARCHY_INITIALIZER::Flag_Valid_Faces");

    // TODO: Finish this
    //LOG::cout<<"Tagging valid faces"<<std::endl;
    typedef VECTOR<int,d> T_INDEX;
    typedef VECTOR<int,d-1> T_INDEX2;

    /////////////////
    ARRAY<unsigned long,T_INDEX> ghost_child_masks(RANGE<T_INDEX>::Unit_Box());
    int child=0;
    for(RANGE_ITERATOR<d> iterator(RANGE<T_INDEX>::Unit_Box());iterator.Valid();iterator.Next(),child++){
        const T_INDEX& index=iterator.Index();
        if(d==3)
            ghost_child_masks(index) = SPGrid_Ghost_Child_000 << child;
        else
            ghost_child_masks(index) = SPGrid_Ghost_Child_000 << (2*child);}
    /////////////////

    PHYSBAM_ASSERT(child==1<<d);

    static const int number_of_faces=(d==2)?4:6;
    unsigned long face_shifts[number_of_faces];
    unsigned long face_ghost_masks[number_of_faces];
    unsigned long face_neighbor_offsets[number_of_faces];
    const int levels = hierarchy.Levels();

    for(int axis=0,face_number=0;axis<d;axis++)
        for(int axis_shift=0;axis_shift<=1;axis_shift++){
            std_array<int,d> neighbor_offset;
            neighbor_offset(axis)=(axis_shift?1:-1);
            face_neighbor_offsets[face_number]=Flag_array_mask::Linear_Offset(neighbor_offset);
            std_array<int,d> shift;
            shift.data[axis]=axis_shift;
            face_shifts[face_number]=Flag_array_mask::Linear_Offset(shift);
                
            face_ghost_masks[face_number]=(unsigned long)0;

            for(RANGE_ITERATOR<d-1> iterator(RANGE<T_INDEX2>(T_INDEX2(),T_INDEX2::All_Ones_Vector()));iterator.Valid();iterator.Next()){
                const T_INDEX index=T_INDEX::All_Ones_Vector()-iterator.Index().Insert(axis_shift,axis+1);
                face_ghost_masks[face_number] |= ghost_child_masks(index);}
            face_number++;}

    for(int level=1;level<=levels;level++)
        for(SPGrid_Set_Iterator<Flag_array_type> iterator(hierarchy.Set(level));iterator.Valid();iterator.Next())
            if(iterator.Data() & SPGrid_Cell_Type_Interior)
                for(int axis=0,face_number=0;axis<d;axis++)
                    for(int axis_shift=0;axis_shift<=1;axis_shift++){
                        const unsigned long face_offset=Flag_array_mask::Packed_Add(face_shifts[face_number],iterator.Offset());
                        const unsigned long neighbor_offset=Flag_array_mask::Packed_Add(iterator.Offset(),face_neighbor_offsets[face_number]);
                        if(!((iterator.Data() & face_ghost_masks[face_number]) == face_ghost_masks[face_number]))
                            hierarchy.Activate_Cell(level,face_offset,GRID_TOPOLOGY_HELPER<Flag_array_mask>::Face_Valid_Mask(axis+1));
                        else if(hierarchy.Set(level).Is_Set(neighbor_offset,(SPGrid_Cell_Type_Dirichlet | SPGrid_Cell_Type_Interior | SPGrid_Cell_Type_Ghost)))
                            hierarchy.Activate_Cell(level,face_offset,GRID_TOPOLOGY_HELPER<Flag_array_mask>::Face_Valid_Mask(axis+1));
                        face_number++;}                            
}
//#####################################################################
// Function Flag_Active_Nodes
//#####################################################################
template<class T_STRUCT,class T,int d> void GRID_HIERARCHY_INITIALIZER<T_STRUCT,T,d>::
Flag_Active_Nodes(Hierarchy_type& hierarchy)
{
    LOG::SCOPE scope("GRID_HIERARCHY_INITIALIZER::Flag_Active_Nodes");

    static const int nodes_per_cell=GRID_TOPOLOGY_HELPER<Flag_array_mask>::nodes_per_cell;
    unsigned long nodes_of_cell_offsets[nodes_per_cell];
    GRID_TOPOLOGY_HELPER<Flag_array_mask>::Nodes_Of_Cell_Offsets(nodes_of_cell_offsets);
    const unsigned long odd_bits=Flag_array_mask::Linear_Offset(coord_t(T_INDEX::All_Ones_Vector()));        
    const int levels=hierarchy.Levels();

    // Flag active nodes at individual levels
    for(int level=1;level<=levels;level++)
        for(SPGrid_Set_Iterator<Flag_array_type> iterator(hierarchy.Set(level));iterator.Valid();iterator.Next())
            if(iterator.Data() & SPGrid_Cell_Type_Interior)
                for(int node=0;node<nodes_per_cell;node++){
                    unsigned long node_offset=Flag_array_mask::Packed_Add(iterator.Offset(),nodes_of_cell_offsets[node]);
                    hierarchy.Activate_Cell(level,node_offset,SPGrid_Node_Active);}

    // Flag active nodes in intermediate levels
    for(int level=1;level<levels-1;level++)
        for(SPGrid_Set_Iterator<Flag_array_type> iterator(hierarchy.Set(level));iterator.Valid();iterator.Next())
            if(iterator.Data() & SPGrid_Node_Active){

                ARRAY<unsigned long> offsets; // Candidate offsets at in-between levels. First one starts at level+1
                bool found=false;             // Indicates a previously marked active node has been found at a coarser level

                for(int new_level=level+1,new_offset=iterator.Offset();new_level<=levels;new_level++){
                    if((new_offset&odd_bits)!=odd_bits) break; // No co-located node exists at the coarser level
                    new_offset=Flag_array_mask::DownsampleOffset(new_offset);
                    if(hierarchy.Set(new_level).Is_Set(new_offset,SPGrid_Node_Active)) {found=true;break;}
                    offsets.Append(new_offset);}
                   
                if(found)
                    for(int dlevel=1;dlevel<=offsets.m;dlevel++)
                        hierarchy.Activate_Cell(level+dlevel,offsets(dlevel),SPGrid_Node_Active);
            }
}
//#####################################################################
// Function Flag_Shared_Nodes
//#####################################################################
template<class T_STRUCT,class T,int d> void GRID_HIERARCHY_INITIALIZER<T_STRUCT,T,d>::
Flag_Shared_Nodes(Hierarchy_type& hierarchy)
{
    LOG::SCOPE scope("GRID_HIERARCHY_INITIALIZER::Flag_Shared_Nodes");

    const unsigned long odd_bits=Flag_array_mask::Linear_Offset(coord_t(T_INDEX::All_Ones_Vector()));
        
    for(int level=1;level<hierarchy.Levels();level++)
        for(SPGrid_Set_Iterator<Flag_array_type> iterator(hierarchy.Set(level));iterator.Valid();iterator.Next())
            if(iterator.Data() & SPGrid_Node_Active){
                unsigned long offset=iterator.Offset();
                if((offset & odd_bits) == odd_bits){
                    unsigned long coarse_offset=Flag_array_mask::DownsampleOffset(offset);
                    if(hierarchy.Set(level+1).Is_Set(coarse_offset,SPGrid_Node_Active))
                        iterator.Data() |= SPGrid_Node_Coarse_Shared;}}
}
//#####################################################################
// Function Generate_Plus_Minus_Active_Faces
//#####################################################################
template<class T_STRUCT,class T,int d> void GRID_HIERARCHY_INITIALIZER<T_STRUCT,T,d>::
Generate_Plus_Minus_Active_Faces(Hierarchy_type& hierarchy)
{
    LOG::SCOPE scope("GRID_HIERARCHY_INITIALIZER::Generate_Plus_Minus_Active_Faces");

    static const int faces_per_cell=GRID_TOPOLOGY_HELPER<Flag_array_mask>::faces_per_cell;

    VECTOR<unsigned long,d> other_face_offsets;
    for(int axis=1;axis<=d;axis++)
        other_face_offsets(axis)=GRID_TOPOLOGY_HELPER<Flag_array_mask>::Axis_Vector_Offset(axis);

    for(int level=1;level<=hierarchy.Levels();level++)
        for(SPGrid_Set_Iterator<Flag_array_type> iterator(hierarchy.Set(level));iterator.Valid();iterator.Next())
            for(int axis=1;axis<=d;axis++){
                unsigned face_active_mask=GRID_TOPOLOGY_HELPER<Flag_array_mask>::Face_Active_Mask(axis);
                // side 1
                if(iterator.Data() & face_active_mask)
                    iterator.Data() |= GRID_TOPOLOGY_HELPER<Flag_array_mask>::Face_Plus_Minus_Active_Mask(axis,1);
                // side 2
                unsigned long other_face_offset=Flag_array_mask::Packed_Add(iterator.Offset(),other_face_offsets(axis));
                if(hierarchy.Set(level).Is_Set(other_face_offset,face_active_mask))
                    iterator.Data() |= GRID_TOPOLOGY_HELPER<Flag_array_mask>::Face_Plus_Minus_Active_Mask(axis,2);}

    // Scaled stuff
    {unsigned long face_neighbor_offsets[d][2];
    const int levels=hierarchy.Levels();
    for(int axis=0;axis<d;axis++)
        for(int axis_shift=-1;axis_shift<=1;axis_shift+=2){
            const int side=(axis_shift+1)/2;
            std_array<int,d> shift;
            shift(axis)=axis_shift;
            face_neighbor_offsets[axis][side]=Flag_array_mask::Linear_Offset(shift);}

    for(int level=1;level<=levels;level++)
        for(SPGrid_Set_Iterator<Flag_array_type> iterator(hierarchy.Set(level));iterator.Valid();iterator.Next())
            if(iterator.Data() & SPGrid_Cell_Type_Interior)
                for(int axis=1;axis<=d;axis++) for(int side=1;side<=2;side++){
                    const unsigned long neighbor_offset=iterator.Offset(face_neighbor_offsets[axis-1][side-1]);
                    if(hierarchy.Set(level).Is_Set(neighbor_offset,SPGrid_Cell_Type_Ghost)){
                        if(iterator.Data() & GRID_TOPOLOGY_HELPER<Flag_array_mask>::Face_Plus_Minus_Active_Mask(axis,side))
                            iterator.Data() |= GRID_TOPOLOGY_HELPER<Flag_array_mask>::Face_Plus_Minus_Scaled_Mask(axis,side);}}
            else if(iterator.Data() & SPGrid_Cell_Type_Ghost)
                for(int axis=1;axis<=d;axis++) for(int side=1;side<=2;side++){
                    const unsigned long neighbor_offset=iterator.Offset(face_neighbor_offsets[axis-1][side-1]);
                    if(hierarchy.Set(level).Is_Set(neighbor_offset,(SPGrid_Cell_Type_Interior|SPGrid_Cell_Type_Dirichlet))){
                        if(iterator.Data() & GRID_TOPOLOGY_HELPER<Flag_array_mask>::Face_Plus_Minus_Active_Mask(axis,side))
                            iterator.Data() |= GRID_TOPOLOGY_HELPER<Flag_array_mask>::Face_Plus_Minus_Scaled_Mask(axis,side);}}}
}
//#####################################################################
// Function Flag_Active_Cells
//#####################################################################
template<class T_STRUCT,class T,int d> void GRID_HIERARCHY_INITIALIZER<T_STRUCT,T,d>::
Flag_Active_Cells(Hierarchy_type& hierarchy)
{
    LOG::SCOPE scope("GRID_HIERARCHY_INITIALIZER::Flag_Active_Cells");

    for(int level=1;level<=hierarchy.Levels();level++)
        for(SPGrid_Set_Iterator<Flag_array_type> iterator(hierarchy.Set(level));iterator.Valid();iterator.Next())
            if( (iterator.Data() & SPGrid_Cell_Type_Interior) && !(iterator.Data() & SPGrid_Cell_Type_Dirichlet) ){
                bool surrounded_by_neumann=true;
                for(int axis=1;axis<=d;axis++) for(int side=1;side<=2;side++)
                    surrounded_by_neumann = surrounded_by_neumann && !(iterator.Data() & GRID_TOPOLOGY_HELPER<Flag_array_type>::Face_Plus_Minus_Active_Mask(axis,side));
                if(!surrounded_by_neumann) iterator.Data() |= SPGrid_Cell_Type_Active;}
}
//#####################################################################
// Function Flag_T_Junction_Nodes
//#####################################################################
template<class T_STRUCT,class T,int d> void GRID_HIERARCHY_INITIALIZER<T_STRUCT,T,d>::
Flag_T_Junction_Nodes(Hierarchy_type& hierarchy)
{
    LOG::SCOPE scope("GRID_HIERARCHY_INITIALIZER::Flag_T_Junction_Nodes");

    // Note: assumes interior cells, active nodes, ghost cells have been marked
    ARRAY<unsigned long,T_INDEX> nodes_of_cell_offsets(RANGE<T_INDEX>::Unit_Box());
    for(RANGE_ITERATOR<d> node_iterator(RANGE<T_INDEX>::Unit_Box());node_iterator.Valid();node_iterator.Next()){
        const T_INDEX& index=node_iterator.Index();
        nodes_of_cell_offsets(index)=Flag_array_mask::Linear_Offset(std_array<int,d>(index));}
    VECTOR<unsigned long,d> axis_vector_offsets;
    for(int v=1;v<=d;v++) axis_vector_offsets(v)=GRID_TOPOLOGY_HELPER<Flag_array_mask>::Axis_Vector_Offset(v);

    static const unsigned active_node_mask=(SPGrid_Node_Active|SPGrid_Node_Ghost);

    for(int level=1;level<=hierarchy.Levels();level++)
        for(SPGrid_Set_Iterator<Flag_array_type> iterator(hierarchy.Set(level));iterator.Valid();iterator.Next())
            if(iterator.Data()&SPGrid_Cell_Type_Ghost){ // for each ghost cell
                const unsigned long offset=iterator.Offset();
                int coarse_level=level;
                unsigned long coarse_offset=offset;
                T_INDEX parity,index_relative_to_parent;
                while(!hierarchy.Set(coarse_level).Is_Set(coarse_offset,SPGrid_Cell_Type_Interior)){ // find index relative to coarse parent
                    for(int v=1;v<=d;v++) parity(v)=(coarse_offset&axis_vector_offsets(v))?0:1;
                    index_relative_to_parent+=(1<<(coarse_level-level))*parity;
                    coarse_level++;coarse_offset=Flag_array_mask::DownsampleOffset(coarse_offset);}
                const int parent_size=1<<(coarse_level-level);
                for(RANGE_ITERATOR<d> node_iterator(RANGE<T_INDEX>::Unit_Box());node_iterator.Valid();node_iterator.Next()){ // for each node of this ghost cell, see if it is T-Junction
                    const T_INDEX& node_index=node_iterator.Index();
                    const unsigned long node_offset=nodes_of_cell_offsets(node_index);
                    if(!hierarchy.Set(level).Is_Set(Flag_array_mask::Packed_Add(offset,node_offset),active_node_mask)) continue; // if not active node, then not T-Junction (continue)
                    bool is_corner_of_parent=true;
                    const T_INDEX node_index_in_coarse=node_index+index_relative_to_parent;
                    for(int v=1;v<=d;v++) if(node_index_in_coarse(v)!=0&&node_index_in_coarse(v)!=parent_size){is_corner_of_parent=false;break;} // if not corner of coarse parent
                    if(!is_corner_of_parent) iterator.Data(node_offset)|=SPGrid_Node_T_Junction;}}
}
//#####################################################################
// Function Flag_Ghost_Nodes
//#####################################################################
template<class T_STRUCT,class T,int d> void GRID_HIERARCHY_INITIALIZER<T_STRUCT,T,d>::
Flag_Ghost_Nodes(Hierarchy_type& hierarchy)
{
    LOG::SCOPE scope("GRID_HIERARCHY_INITIALIZER::Flag_T_Junction_Nodes");

    // Note: assumes interior cells, ghost cells and scaled faces have been marked

    unsigned long face_neighbor_offsets[d][2];
    for(int axis=0;axis<d;axis++)
        for(int axis_shift=-1;axis_shift<=1;axis_shift+=2){
            const int side=(axis_shift+1)/2;
            std_array<int,d> shift;
            shift(axis)=axis_shift;
            face_neighbor_offsets[axis][side]=Flag_array_mask::Linear_Offset(shift);}
    static const int nodes_per_face=GRID_TOPOLOGY_HELPER<Flag_array_mask>::nodes_per_face;
    unsigned long nodes_of_face_offsets[d][nodes_per_face];
    for(int axis=0;axis<d;axis++) GRID_TOPOLOGY_HELPER<Flag_array_mask>::Nodes_Of_Face_Offsets(nodes_of_face_offsets[axis],axis+1);
    VECTOR<unsigned long,d> axis_vector_offsets;
    for(int v=1;v<=d;v++) axis_vector_offsets(v)=GRID_TOPOLOGY_HELPER<Flag_array_mask>::Axis_Vector_Offset(v);

    const int levels=hierarchy.Levels();
    for(int level=1;level<=levels;level++)
        for(SPGrid_Set_Iterator<Flag_array_type> iterator(hierarchy.Set(level));iterator.Valid();iterator.Next())
            if(iterator.Data()&SPGrid_Cell_Type_Interior)
                for(int axis=0;axis<d;axis++)
                    for(int side=0;side<2;side++){
                        int current_level=level;
                        unsigned long neighbor_offset=iterator.Offset(face_neighbor_offsets[axis][side]);
                        while(hierarchy.Set(current_level).Is_Set(neighbor_offset,SPGrid_Cell_Type_Ghost)){
                            const unsigned long face_offset=(side==0)?Flag_array_mask::Packed_Add(neighbor_offset,axis_vector_offsets(axis+1)):neighbor_offset;
                            for(int node=0;node<nodes_per_face;node++){
                                const unsigned long node_offset=Flag_array_mask::Packed_Add(face_offset,nodes_of_face_offsets[axis][node]);
                                hierarchy.Activate_Cell(current_level,node_offset,SPGrid_Node_Ghost);}
                            neighbor_offset=Flag_array_mask::DownsampleOffset(neighbor_offset);
                            current_level++;}}
}
//#####################################################################
template class GRID_HIERARCHY_INITIALIZER<FLUIDS_SIMULATION_DATA<float>,float,2>;
template class GRID_HIERARCHY_INITIALIZER<FLUIDS_SIMULATION_DATA<float>,float,3>;
