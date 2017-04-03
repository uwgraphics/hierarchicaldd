//#####################################################################
// Copyright 2013, Raj Setaluri, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Face_Advection_Helper
//#####################################################################
#include <PhysBAM_Tools/Interpolation/LINEAR_INTERPOLATION.h>
#include <SPGrid_Fluids/Advection/Face_Advection_Helper.h>
#include <SPGrid_Fluids/Advection/BACKTRACE.h>
#include <SPGrid_Fluids/Interpolation/GRID_HIERARCHY_INTERPOLATION.h>
#include <SPGrid_Fluids/Simulation/FLUIDS_SIMULATION_DATA.h>
using namespace PhysBAM;
//#####################################################################
// Face_Advection_Helper
//#####################################################################
template<class T_STRUCT,class T_DATA,class T_FLAGS,int d> void Face_Advection_Helper<T_STRUCT,T_DATA,T_FLAGS,d>::
Run(SPGrid_Allocator<T_STRUCT,d>& allocator,const std::pair<const unsigned long*,unsigned>& blocks) const
{
    const T one_over_nodes_per_face=(T)1./(T)nodes_per_face;
    // keep PhysBAM index as well as unsigned offset
    RANGE_ITERATOR<d> range_iterator(RANGE<TV_INT>(TV_INT(),allocator.Block_Size().template Cast<TV_INT>()-1));
    TV_INT base_index;
    Const_flag_array_type flags=allocator.Get_Const_Array(flags_field);
    Const_data_array_type native_face_velocity=allocator.Get_Const_Array(face_velocities(axis));
    Data_array_type result=allocator.Get_Array(result_field); // store result in temp_field
    // iterate over all allocated cells at this level
    for(SPGrid_Block_Iterator<Flag_array_mask> iterator(blocks);iterator.Valid();iterator.Next_Block()){
        range_iterator.Reset();
        base_index=iterator.Index().template Cast<TV_INT>(); // update physbam index
        unsigned long offset=iterator.Offset();
        for(int element=0;element<elements_per_block;element++,range_iterator.Next(),offset+=data_size){
            // if face on this axis is valid
            if(flags(offset)&face_advect_mask){
                const TV_INT index=base_index+range_iterator.Index();
                // grab vector of velocities, one native, others from nodes
                TV velocity=TV();
                velocity(axis)=native_face_velocity(offset);
#if 1
                for(int v=(axis%d)+1;v!=axis;v=(v%d)+1){
                    // average from 2^(d-1) nodes
                    for(int node=0;node<nodes_per_face;node++)
                        velocity(v)+=*reinterpret_cast<T*>(reinterpret_cast<unsigned long>(node_velocity_ptrs(v))+Data_array_type::MASK::Packed_Add(offset,nodes_of_face_offsets[node]));
                    velocity(v)*=one_over_nodes_per_face;}
#else
                for(int node=0;node<nodes_per_face;node++){
                    unsigned long node_offset=Data_array_type::MASK::Packed_Add(offset,nodes_of_face_offsets[node]);
                    for(int v=(axis%d)+1;v!=axis;v=(v%d)+1)
                        velocity(v)+=*reinterpret_cast<T*>(reinterpret_cast<unsigned long>(node_velocity_ptrs(v))+node_offset);}
                for(int v=(axis%d)+1;v!=axis;v=(v%d)+1)
                    velocity(v)*=one_over_nodes_per_face;
#endif
                // call backtrace with information -- get new index and weights
                TV dX=-velocity*dt;
                unsigned long new_offset=offset;
                int new_level=level;
                TV weights=BACKTRACE<T_STRUCT,T,d>::Backtrace(hierarchy,new_level,index,new_offset,intra_cell_dX,dX,rigid_geometry_collection);
                unsigned long other_offset=Flag_array_mask::Packed_Add(new_offset,other_face_offsets(axis));
                T val1=GRID_HIERARCHY_INTERPOLATION<T_STRUCT,T,d>::Face_Interpolation_Helper(hierarchy,nodes_of_face_offsets,axis,new_level,new_offset  ,weights.Remove_Index(axis),flags_field,face_velocities,node_velocities);
                T val2=GRID_HIERARCHY_INTERPOLATION<T_STRUCT,T,d>::Face_Interpolation_Helper(hierarchy,nodes_of_face_offsets,axis,new_level,other_offset,weights.Remove_Index(axis),flags_field,face_velocities,node_velocities);
                result(offset)=LINEAR_INTERPOLATION<T,T>::Linear(val1,val2,weights(axis));}}}
}
//#####################################################################
template class Face_Advection_Helper<FLUIDS_SIMULATION_DATA<float>,float,unsigned,2>;
template class Face_Advection_Helper<FLUIDS_SIMULATION_DATA<float>,float,unsigned,3>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class Face_Advection_Helper<FLUIDS_SIMULATION_DATA<double>,double,unsigned,2>;
template class Face_Advection_Helper<FLUIDS_SIMULATION_DATA<double>,double,unsigned,3>;
#endif
