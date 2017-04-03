//#####################################################################
// Copyright 2013, Raj Setaluri, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Density_Advection_Helper
//#####################################################################
#include <PhysBAM_Tools/Interpolation/LINEAR_INTERPOLATION.h>
#include <SPGrid_Fluids/Advection/Density_Advection_Helper.h>
#include <SPGrid_Fluids/Interpolation/GRID_HIERARCHY_INTERPOLATION.h>
#include <SPGrid_Fluids/Simulation/FLUIDS_SIMULATION_DATA.h>
#include <SPGrid_Fluids/Advection/BACKTRACE.h>
using namespace PhysBAM;
//#####################################################################
// Run
//#####################################################################
template<class T_STRUCT,class T_DATA,class T_FLAGS,int d> void Density_Advection_Helper<T_STRUCT,T_DATA,T_FLAGS,d>::
Run(SPGrid_Allocator<T_STRUCT,d>& allocator,const std::pair<const unsigned long*,unsigned>& blocks) const
{
    // keep PhysBAM index as well as unsigned offset
    RANGE_ITERATOR<d> range_iterator(RANGE<TV_INT>(TV_INT(),allocator.Block_Size().template Cast<TV_INT>()-1));
    TV_INT base_index;
    Const_flag_array_type flags=allocator.Get_Const_Array(flags_field);
    Data_array_type result=allocator.Get_Array(result_field);
    // iterate over all allocated cells at this level
    //PHYSBAM_ASSERT(elements_per_block==(allocator.Block_Size().template Cast<TV_INT>()).Product());
    for(SPGrid_Block_Iterator<Flag_array_mask> iterator(blocks);iterator.Valid();iterator.Next_Block()){
        range_iterator.Reset();
        base_index=iterator.Index().template Cast<TV_INT>(); // update physbam index
        unsigned long offset=iterator.Offset();
        for(int element=0;element<elements_per_block;element++,range_iterator.Next(),offset+=data_size){
        // if this cell is interior
        if(flags(offset)&density_advect_mask){
            const TV_INT index=base_index+range_iterator.Index();
            // Fill out velocity
            TV velocity=TV();
            for(int node=0;node<nodes_per_cell;node++){
                unsigned long node_offset=Data_array_type::MASK::Packed_Add(offset,nodes_of_cell_offsets[node]);
                for(int v=1;v<=d;v++)
                    velocity(v)+=*reinterpret_cast<T*>(reinterpret_cast<unsigned long>(node_velocity_ptrs(v))+node_offset);}
            velocity/=(T)nodes_per_cell;
            // call backtrace with information -- get new index and weights
            TV dX=-velocity*dt;
            unsigned long new_offset=offset;
            int new_level=level;
            TV weights=BACKTRACE<T_STRUCT,T,d>::Backtrace(hierarchy,new_level,index,new_offset,intra_cell_dX,dX,rigid_geometry_collection);
            T val=GRID_HIERARCHY_INTERPOLATION<T_STRUCT,T,d>::Cell_Interpolation_Helper(hierarchy,nodes_of_cell_offsets,new_level,new_offset,weights,density_field,node_density_field);
            result(offset)=val;}}}
}
//#####################################################################
template class Density_Advection_Helper<FLUIDS_SIMULATION_DATA<float>,float,unsigned,2>;
template class Density_Advection_Helper<FLUIDS_SIMULATION_DATA<float>,float,unsigned,3>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class Density_Advection_Helper<FLUIDS_SIMULATION_DATA<double>,double,unsigned,2>;
template class Density_Advection_Helper<FLUIDS_SIMULATION_DATA<double>,double,unsigned,3>;
#endif
