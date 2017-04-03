//#####################################################################
//Removed For Anonymity: Copyright Authors of Submission pap173s1
//Removed For Anonymity: Copyright Authors of Submission pap173s1
//#####################################################################
// Class BOUNDARY_REFINEMENT_RASTERIZER
//#####################################################################
#ifndef __BOUNDARY_REFINEMENT_RASTERIZER__
#define __BOUNDARY_REFINEMENT_RASTERIZER__

#include "SPGrid/Core/SPGrid_Utilities.h"
#include <SPGrid_Fluids/Grids/GRID_HIERARCHY_INITIALIZER.h>
#include <SPGrid_Fluids/Solvers/Domain_Decomposition/SPGRID_MULTIGRID_FLAGS.h>
#include "ADAPTIVE_SUBDOMAIN_POISSON_FLAGS.h"

namespace PhysBAM{
template<class T_STRUCT_ADAPTATION,class T_STRUCT_SOLVER,class T,int d>
class BOUNDARY_REFINEMENT_RASTERIZER
{
    typedef VECTOR<int,d> T_INDEX;
    typedef PAIR<int,T_INDEX> T_CELL;

    typedef SPGrid_Allocator<T_STRUCT_SOLVER,d> SPG_Allocator;
    typedef typename SPG_Allocator::template Array<unsigned>::type SPG_Flags_Array_Type;
    typedef typename SPG_Allocator::template Array<T>::type SPG_Data_Array_Type;
    typedef typename SPG_Allocator::template Array<T>::mask T_MASK;
    typedef SPGrid_Set<SPG_Flags_Array_Type> SPG_Set_Type;

    const SPG_Set_Type& set;
    GRID_HIERARCHY<T_STRUCT_ADAPTATION,T,d>& hierarchy;
    int mg_level;
    RANGE<T_INDEX> subdomain_range;
    ARRAY<unsigned long,T_INDEX> node_offset;
public:
    BOUNDARY_REFINEMENT_RASTERIZER(GRID_HIERARCHY<T_STRUCT_ADAPTATION,T,d>& hierarchy_input,const SPG_Set_Type& set_input,const RANGE<T_INDEX> subdomain_range_input,int mg_level_input)
        :hierarchy(hierarchy_input),set(set_input),mg_level(mg_level_input),subdomain_range(subdomain_range_input)
    {
        node_offset.Resize(RANGE<T_INDEX>::Unit_Box());
        for(RANGE_ITERATOR<d> node_iterator(RANGE<T_INDEX>::Unit_Box());node_iterator.Valid();node_iterator.Next()){
            node_offset(node_iterator.Index())=T_MASK::Linear_Offset(std_array<int,d>(node_iterator.Index()-1));}
    }

    bool Consume(const T_CELL& cell)
    {
        const int level=cell.x;
        const T_INDEX& index=cell.y;
        // If we have hierarchically traversed all levels, and haven't marked this region,
        // then do mark it at level one
        RANGE<T_INDEX> cell_range(index,index);
        for(int i=level;i>1;--i){
            cell_range.min_corner=2*cell_range.min_corner-1;
            cell_range.max_corner=2*cell_range.max_corner;}            
        if(RANGE<T_INDEX>::Intersect(cell_range,subdomain_range).Empty()) return false;
        if(level==mg_level&&(!subdomain_range.Contains(cell_range))) return false;
        if(level==mg_level||subdomain_range.Thickened(-(1<<(level-1))).Contains(cell_range)){
            for(RANGE_ITERATOR<d> iterator(cell_range);iterator.Valid();iterator.Next()){
                const T_INDEX& cell_index=iterator.Index();
                unsigned long cell_offset=T_MASK::Linear_Offset(std_array<int,d>(cell_index));
                for(int v=1;v<=d;v++)
                    for(RANGE_ITERATOR<d-1> face_iterator(RANGE<VECTOR<int,d-1> >(VECTOR<int,d-1>(),VECTOR<int,d-1>()+1));
                        face_iterator.Valid();face_iterator.Next()){
                        T_INDEX lower_index=face_iterator.Index().Insert(0,v);
                        T_INDEX upper_index=face_iterator.Index().Insert(1,v);
                        unsigned long lower_offset=T_MASK::Packed_Add(cell_offset,node_offset(lower_index));
                        unsigned long upper_offset=T_MASK::Packed_Add(cell_offset,node_offset(upper_index));
                        unsigned lower_type=0,upper_type=0;
                        if(set.Is_Set(lower_offset,SPGrid_Solver_Cell_Type_Active)||set.Is_Set(lower_offset,SPGrid_Solver_Cell_Type_Interface))
                            lower_type |= Uniform_Node_Type_Active;
                        if(set.Is_Set(lower_offset,SPGrid_Solver_Cell_Type_Dirichlet))
                            lower_type |= Uniform_Node_Type_Dirichlet;
                        if(set.Is_Set(upper_offset,SPGrid_Solver_Cell_Type_Active)||set.Is_Set(upper_offset,SPGrid_Solver_Cell_Type_Interface))
                            upper_type |= Uniform_Node_Type_Active;
                        if(set.Is_Set(upper_offset,SPGrid_Solver_Cell_Type_Dirichlet))
                            upper_type |= Uniform_Node_Type_Dirichlet;

                        if(lower_type&upper_type&Uniform_Node_Type_Active){ // If both of them are active
                            hierarchy.Activate_Cell(level,index,Adaptive_Cell_Type_Interior);
                            return false;}
                        else if(((lower_type|upper_type)&Uniform_Node_Type_Active) && ((lower_type|upper_type)&Uniform_Node_Type_Dirichlet)){
                            hierarchy.Activate_Cell(level,index,Adaptive_Cell_Type_Interior);
                            return false;}}}
            return false;}
        // Must be on the margin of a level greater than one. Recurse.
        return true;
    }    
};
}
#endif
