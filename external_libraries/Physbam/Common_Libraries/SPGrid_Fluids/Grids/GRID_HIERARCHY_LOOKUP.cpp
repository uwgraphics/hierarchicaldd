//#####################################################################
// Copyright 2013, Raj Setaluri, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <SPGrid_Fluids/Grids/GRID_HIERARCHY_LOOKUP.h>
#include <SPGrid_Fluids/Simulation/FLUIDS_SIMULATION_DATA.h>
#include <SPGrid_Fluids/Simulation/FLUIDS_SIMULATION_FLAGS.h>
using namespace PhysBAM;
//#####################################################################
// Cell_Lookup
//#####################################################################
template<class T_STRUCT,class T,int d> bool GRID_HIERARCHY_LOOKUP<T_STRUCT,T,d>::
Cell_Lookup(const Hierarchy_type& hierarchy,const TV& X,unsigned long &flags_offset,int &level)
{
    TV weights;
    return Cell_Lookup(hierarchy,X,flags_offset,level,weights);
}
//#####################################################################
// Cell_Lookup
//#####################################################################
template<class T_STRUCT,class T,int d> bool GRID_HIERARCHY_LOOKUP<T_STRUCT,T,d>::
Cell_Lookup(const Hierarchy_type& hierarchy,const TV& X,unsigned long &flags_offset,TV_INT &index,int &level,TV& weights)
{
    level=1;
    index=hierarchy.Grid(level).Clamp_To_Cell(X);
    TV node_location=hierarchy.Grid(level).Node(index);    // location of minimum corner node
    weights=(X-node_location)/hierarchy.Grid(level).dX;
    flags_offset=Flag_array_type::MASK::Linear_Offset(std_array<int,d>(index));
    return Cell_Lookup(hierarchy,flags_offset,level,weights);
}
//#####################################################################
// Cell_Lookup
//#####################################################################
template<class T_STRUCT,class T,int d> bool GRID_HIERARCHY_LOOKUP<T_STRUCT,T,d>::
Cell_Lookup(const Hierarchy_type& hierarchy,const TV& X,unsigned long &flags_offset,int &level,TV& weights)
{
    level=1;
    TV_INT cell_index=hierarchy.Grid(level).Clamp_To_Cell(X);
    TV node_location=hierarchy.Grid(level).Node(cell_index);    // location of minimum corner node
    weights=(X-node_location)/hierarchy.Grid(level).dX;
    flags_offset=Flag_array_type::MASK::Linear_Offset(std_array<int,d>(cell_index));
    return Cell_Lookup(hierarchy,flags_offset,level,weights);
}
//#####################################################################
// Cell_Lookup
//#####################################################################
template<class T_STRUCT,class T,int d> bool GRID_HIERARCHY_LOOKUP<T_STRUCT,T,d>::
Cell_Lookup(const Hierarchy_type& hierarchy,unsigned long &flags_offset,int &level,TV& weights)
{
    if(hierarchy.Set(level).Is_Set(flags_offset,SPGrid_Cell_Type_Interior))
        return true;
    if((hierarchy.Set(level).Is_Set(flags_offset,SPGrid_Cell_Type_Ghost)) || !Cell_Lookup_Finer(hierarchy,flags_offset,level,weights))
        return Cell_Lookup_Coarser(hierarchy,flags_offset,level,weights);
    return true;
}
//#####################################################################
// Cell_Lookup_Coarser
//#####################################################################
template<class T_STRUCT,class T,int d> bool GRID_HIERARCHY_LOOKUP<T_STRUCT,T,d>::
Cell_Lookup_Coarser(const Hierarchy_type& hierarchy,unsigned long &flags_offset,int &level,TV& weights)
{
    unsigned long candidate_flags_offset=flags_offset;
    TV candidate_weights=weights;

    for(unsigned candidate_level=level+1;candidate_level<=hierarchy.Levels();candidate_level++){

        // Downsample weights
        candidate_weights*=(T).5f;
        if((flags_offset & Cell_Parity_Helper<Flag_array_mask>::x_mask) == 0)
            candidate_weights(1)+=(T).5;
        if((flags_offset & Cell_Parity_Helper<Flag_array_mask>::y_mask) == 0)
            candidate_weights(2)+=(T).5;
        if(d==3)
            if((flags_offset & Cell_Parity_Helper<Flag_array_mask>::z_mask) == 0)
                candidate_weights(3)+=(T).5;
        
        // Downsample offset
        candidate_flags_offset=Flag_array_mask::DownsampleOffset(candidate_flags_offset);

        // check for this offset being interior
        if(hierarchy.Set(candidate_level).Is_Set(candidate_flags_offset,SPGrid_Cell_Type_Interior)){
            flags_offset=candidate_flags_offset;
            level=candidate_level;
            weights=candidate_weights;
            return true;
        }
    }
    
    return false;

}
//#####################################################################
// Cell_Lookup_Finer
//#####################################################################
template<class T_STRUCT,class T,int d> bool GRID_HIERARCHY_LOOKUP<T_STRUCT,T,d>::
Cell_Lookup_Finer(const Hierarchy_type& hierarchy,unsigned long &flags_offset,int &level,TV& weights)
{
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<unsigned>::type::MASK Flag_array_mask;
    
    unsigned long candidate_flags_offset=flags_offset;
    TV candidate_weights=weights;
    
    for(unsigned candidate_level=level-1;candidate_level>=1;candidate_level--){
        
        // Upsample offset
        candidate_flags_offset=Flag_array_mask::UpsampleOffset(candidate_flags_offset);
        unsigned long child_offset=0UL;

        // Adjust new weights and offset by old weights
        if(candidate_weights(1) >= (T).5){
            child_offset |= Cell_Parity_Helper<Flag_array_mask>::x_mask;
            candidate_weights(1)-=(T).5;
        }
        if(candidate_weights(2) >= (T).5){
            child_offset |= Cell_Parity_Helper<Flag_array_mask>::y_mask;
            candidate_weights(2)-=(T).5;
        }
        if(d==3)
            if(candidate_weights(3) >= (T).5){
            child_offset |= Cell_Parity_Helper<Flag_array_mask>::z_mask;
                candidate_weights(3)-=(T).5;
            }
        candidate_flags_offset=Flag_array_mask::Packed_Add(candidate_flags_offset,child_offset);

        // Upsample weights
        candidate_weights*=(T)2.;
        
        // check for this offset being interior
        if(hierarchy.Set(candidate_level).Is_Set(candidate_flags_offset,SPGrid_Cell_Type_Interior)){
            flags_offset=candidate_flags_offset;
            level=candidate_level;
            weights=candidate_weights;
            return true;
        }        
    }    

    return false;

}
//#####################################################################
// Face_Lookup
//#####################################################################
template<class T_STRUCT,class T,int d> bool GRID_HIERARCHY_LOOKUP<T_STRUCT,T,d>::
Face_Lookup(const Hierarchy_type& hierarchy,unsigned long &flags_offset,int &level,TV2& weights,const unsigned mask,const int axis)
{
    // if exists at this region, return
    if(hierarchy.Set(level).Is_Set(flags_offset,mask)) return true;
    // otherwise look finer
    return Face_Lookup_Finer(hierarchy,flags_offset,level,weights,mask,axis);
}
//#####################################################################
// Face_Lookup_Finer
//#####################################################################
template<class T_STRUCT,class T,int d> bool GRID_HIERARCHY_LOOKUP<T_STRUCT,T,d>::
Face_Lookup_Finer(const Hierarchy_type& hierarchy,unsigned long &flags_offset,int &level,TV2& weights,const unsigned mask,const int axis)
{
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<unsigned>::type::MASK Flag_array_mask;
    
    unsigned long candidate_flags_offset=flags_offset;
    TV2 candidate_weights=weights;
    
    for(unsigned candidate_level=level-1;candidate_level>=1;candidate_level--){
        
        // Upsample offset
        candidate_flags_offset=Flag_array_mask::UpsampleOffset(candidate_flags_offset);
        unsigned long child_offset=0UL;
        // Adjust new weights and offset by old weights
        if(candidate_weights(1) >= (T).5){
            child_offset |= (axis==1)?Cell_Parity_Helper<Flag_array_mask>::y_mask:Cell_Parity_Helper<Flag_array_mask>::x_mask;
            candidate_weights(1)-=(T).5;}
        if(d==3)
        if(candidate_weights(2) >= (T).5){
            child_offset |= (axis==3)?Cell_Parity_Helper<Flag_array_mask>::y_mask:Cell_Parity_Helper<Flag_array_mask>::z_mask;
            candidate_weights(2)-=(T).5;}
        candidate_flags_offset=Flag_array_mask::Packed_Add(candidate_flags_offset,child_offset);

        // Upsample weights
        candidate_weights*=(T)2.;
        
        // check for this offset being interior
        if(hierarchy.Set(candidate_level).Is_Set(candidate_flags_offset,mask)){
            flags_offset=candidate_flags_offset;
            level=candidate_level;
            weights=candidate_weights;
            return true;}
    }    
    return false;
}
//#####################################################################
template class GRID_HIERARCHY_LOOKUP<FLUIDS_SIMULATION_DATA<float>,float,2>;
template class GRID_HIERARCHY_LOOKUP<FLUIDS_SIMULATION_DATA<float>,float,3>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class GRID_HIERARCHY_LOOKUP<FLUIDS_SIMULATION_DATA<double>,double,2>;
template class GRID_HIERARCHY_LOOKUP<FLUIDS_SIMULATION_DATA<double>,double,3>;
#endif
