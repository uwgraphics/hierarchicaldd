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
    VECTOR<void*,d> node_velocity_ptrs,vorticity_ptrs;

    for(int level=1;level<=hierarchy.Levels();++level){const TV& one_over_dx=hierarchy.Grid(level).One_Over_DX();
        for(int v=1;v<=d;v++){node_velocity_ptrs(v)=hierarchy.Array(level,node_velocities(v)).Get_Data_Ptr();
            vorticity_ptrs(v)=hierarchy.Array(level,vorticity(v)).Get_Data_Ptr();}
        Const_flag_array_type flags=hierarchy.Allocator(level).Get_Const_Array(flags_field);
        for(SPGrid_Block_Iterator<typename Flag_array_type::MASK> iterator(hierarchy.Blocks(level));iterator.Valid();iterator.Next()){
            TV_INT index=iterator.Index().template Cast<TV_INT>();
            unsigned cell_flags=iterator.Data(flags);
            if(cell_flags&(SPGrid_Cell_Type_Interior)){ARRAY<TV> nodal_velocities(nodes_per_cell);
                for(int j=1;j<=nodes_per_cell;++j) for(int k=1;k<=d;++k) nodal_velocities(j)[k]=iterator.template Data<Data_array_type>(node_velocity_ptrs[k],nodes_of_cell_offsets[j-1]);
                typename TV::SPIN cell_vorticity;
                cell_vorticity.x=(T).25*one_over_dx.y*(nodal_velocities(3).z+nodal_velocities(4).z+nodal_velocities(7).z+nodal_velocities(8).z-nodal_velocities(1).z-nodal_velocities(2).z-nodal_velocities(5).z-nodal_velocities(6).z)-(T).25*one_over_dx.z*(nodal_velocities(2).y+nodal_velocities(4).y+nodal_velocities(6).y+nodal_velocities(8).y-nodal_velocities(1).y-nodal_velocities(3).y-nodal_velocities(5).y-nodal_velocities(7).y);
                cell_vorticity.y=(T).25*one_over_dx.z*(nodal_velocities(2).x+nodal_velocities(4).x+nodal_velocities(6).x+nodal_velocities(8).x-nodal_velocities(1).x-nodal_velocities(3).x-nodal_velocities(5).x-nodal_velocities(7).x)-(T).25*one_over_dx.x*(nodal_velocities(5).z+nodal_velocities(6).z+nodal_velocities(7).z+nodal_velocities(8).z-nodal_velocities(1).z-nodal_velocities(2).z-nodal_velocities(3).z-nodal_velocities(4).z);
                cell_vorticity.z=(T).25*one_over_dx.x*(nodal_velocities(5).y+nodal_velocities(6).y+nodal_velocities(7).y+nodal_velocities(8).y-nodal_velocities(1).y-nodal_velocities(2).y-nodal_velocities(3).y-nodal_velocities(4).y)-(T).25*one_over_dx.y*(nodal_velocities(3).x+nodal_velocities(4).x+nodal_velocities(7).x+nodal_velocities(8).x-nodal_velocities(1).x-nodal_velocities(2).x-nodal_velocities(5).x-nodal_velocities(6).x);
                for(int axis=1;axis<=d;++axis) iterator.template Data<Data_array_type>(vorticity_ptrs[axis])=cell_vorticity[axis];}}}
}
//#####################################################################
template class GRID_HIERARCHY_VORTICITY<FLUIDS_SIMULATION_DATA<float>,float,3>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class GRID_HIERARCHY_VORTICITY<FLUIDS_SIMULATION_DATA<double>,double,3>;
#endif
