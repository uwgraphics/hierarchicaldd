//#####################################################################
// Copyright 2014, Mridul Aanjaneya, Raj Setaluri.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class GRID_HIERARCHY_VORTICITY
//#####################################################################
#include <SPGrid/Tools/SPGrid_Block_Iterator.h>
#include <SPGrid/Tools/SPGrid_Clear.h>
#include <SPGrid/Tools/SPGrid_Copy.h>
#include <SPGrid/Tools/SPGrid_Threading_Helper.h>
#include <SPGrid_Fluids/Simulation/FLUIDS_SIMULATION_DATA.h>
#include <SPGrid_Fluids/Vorticity/GRID_HIERARCHY_VORTICITY.h>

using namespace PhysBAM;

namespace PhysBAM{
extern int PhysBAM_number_of_threads;
}
//#####################################################################
// Function Compute_Vorticity
//#####################################################################
template<class T_STRUCT,class T,int d> void GRID_HIERARCHY_VORTICITY<T_STRUCT,T,d>::
Compute_Vorticity(T_HIERARCHY& hierarchy,const VECTOR<T T_STRUCT::*,d> node_velocities,unsigned T_STRUCT::* const flags_field,const typename VECTOR<T T_STRUCT::*,d>::SPIN vorticity)
{
    unsigned long nodes_of_cell_offsets[nodes_per_cell];
    GRID_TOPOLOGY_HELPER<typename Flag_array_type::MASK>::Nodes_Of_Cell_Offsets(nodes_of_cell_offsets);
    VECTOR<void*,d> node_velocity_ptrs;void* vorticity_ptr;
    
    for(int level=1;level<=hierarchy.Levels();++level){const TV& one_over_dx=hierarchy.Grid(level).One_Over_DX();
        for(int v=1;v<=d;v++){node_velocity_ptrs(v)=hierarchy.Array(level,node_velocities(v)).Get_Data_Ptr();}
        vorticity_ptr=hierarchy.Array(level,vorticity(1)).Get_Data_Ptr();
        Const_flag_array_type flags=hierarchy.Allocator(level).Get_Const_Array(flags_field);
        for(SPGrid_Block_Iterator<typename Flag_array_type::MASK> iterator(hierarchy.Blocks(level));iterator.Valid();iterator.Next()){
            TV_INT index=iterator.Index().template Cast<TV_INT>();
            unsigned cell_flags=iterator.Data(flags);
            if(cell_flags&(SPGrid_Cell_Type_Interior)){ARRAY<TV> nodal_velocities(nodes_per_cell);
                for(int j=1;j<=nodes_per_cell;++j) for(int k=1;k<=d;++k) nodal_velocities(j)[k]=iterator.template Data<Data_array_type>(node_velocity_ptrs[k],nodes_of_cell_offsets[j-1]);
                T cell_vorticity_x=(T).5*one_over_dx.x*(nodal_velocities(3).y+nodal_velocities(4).y-nodal_velocities(1).y-nodal_velocities(2).y);
                T cell_vorticity_y=(T).5*one_over_dx.y*(nodal_velocities(2).x+nodal_velocities(4).x-nodal_velocities(1).x-nodal_velocities(3).x);
                T cell_vorticity=cell_vorticity_x-cell_vorticity_y;
                iterator.template Data<Data_array_type>(vorticity_ptr)=cell_vorticity;}}}
}
//#####################################################################
template class GRID_HIERARCHY_VORTICITY<FLUIDS_SIMULATION_DATA<float>,float,2>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class GRID_HIERARCHY_VORTICITY<FLUIDS_SIMULATION_DATA<double>,double,2>;
#endif
