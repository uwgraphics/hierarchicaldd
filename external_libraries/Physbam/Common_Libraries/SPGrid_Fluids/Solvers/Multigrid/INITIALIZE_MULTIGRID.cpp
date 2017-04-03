//#####################################################################
// Copyright (c) 2014, Raj Setaluri, Mridul Aanjneya
// Distributed under the FreeBSD license (see license.txt)
//#####################################################################
#include <SPGrid_Fluids/Solvers/Multigrid/INITIALIZE_MULTIGRID.h>
#include <SPGrid_Fluids/Simulation/FLUIDS_SIMULATION_DATA.h>
#include <SPGrid_Fluids/Simulation/FLUIDS_SIMULATION_FLAGS.h>
#include <SPGrid/Tools/SPGrid_Block_Iterator.h>
#include <SPGrid_Fluids/Grids/GRID_HIERARCHY_INITIALIZER.h>
#include <Common_Tools/Math_Tools/RANGE_ITERATOR.h>

using namespace PhysBAM;
//#####################################################################
// Initialize_Multigrid_Flags
//#####################################################################
template<class T_STRUCT,class T,int d> void INITIALIZE_MULTIGRID<T_STRUCT,T,d>::
Initialize_Multigrid_Flags(T_HIERARCHY& hierarchy,unsigned T_STRUCT::* flags_channel,unsigned T_STRUCT::* multigrid_flags_channel)
{
    PHYSBAM_NOT_IMPLEMENTED();
}
//#####################################################################
// Initialize_Multigrid_Flags
//#####################################################################
template<class T_STRUCT,class T,int d> void INITIALIZE_MULTIGRID<T_STRUCT,T,d>::
Initialize_Multigrid_Flags(T_HIERARCHY& base_hierarchy,ARRAY<T_HIERARCHY*>& multigrid_hierarchy,const int mg_levels)
{
    const int spgrid_levels=base_hierarchy.Levels();
    GRID<TV> fine_mac_grid(base_hierarchy.Grid(1).Get_MAC_Grid());
    // cleanup
    for(int i=1;i<=multigrid_hierarchy.m;i++) if(multigrid_hierarchy(i)) delete multigrid_hierarchy(i);
    multigrid_hierarchy.Resize(mg_levels);
    multigrid_hierarchy.Fill(0);
    // instance hierarchies
    for(int mg_level=1;mg_level<=mg_levels;mg_level++){
        multigrid_hierarchy(mg_level)=new T_HIERARCHY(fine_mac_grid,max(spgrid_levels,mg_levels));
        multigrid_hierarchy(mg_level)->Initialize_Sets();}
    // copy "interior" information from base hierarchy to mg hierarchies
    for(int spgrid_level=1;spgrid_level<=spgrid_levels;spgrid_level++){
        Const_flag_array_type flags=base_hierarchy.Allocator(spgrid_level).Get_Const_Array(&T_STRUCT::flags);
        for(SPGrid_Block_Iterator<typename Flag_array_type::MASK> iterator(base_hierarchy.Blocks(spgrid_level));iterator.Valid();iterator.Next())
            if(iterator.Data(flags)&SPGrid_Cell_Type_Active){
                unsigned long offset=iterator.Offset();
                for(int mg_level=1;mg_level<=mg_levels;mg_level++){
                    if(mg_level<=spgrid_level) multigrid_hierarchy(mg_level)->Activate_Cell(spgrid_level,offset,SPGrid_Cell_Type_Interior);
                    else{offset=Flag_array_type::MASK::DownsampleOffset(offset);
                        multigrid_hierarchy(mg_level)->Activate_Cell(mg_level,offset,SPGrid_Cell_Type_Interior);}}}}
    // copy "dirichlet" information from base hierarchy to mg hierarchies
    for(int spgrid_level=1;spgrid_level<=spgrid_levels;spgrid_level++){
        Const_flag_array_type flags=base_hierarchy.Allocator(spgrid_level).Get_Const_Array(&T_STRUCT::flags);
        for(SPGrid_Block_Iterator<typename Flag_array_type::MASK> iterator(base_hierarchy.Blocks(spgrid_level));iterator.Valid();iterator.Next())
            if(iterator.Data(flags)&SPGrid_Cell_Type_Dirichlet){
                unsigned long offset=iterator.Offset();
                for(int mg_level=1;mg_level<=mg_levels;mg_level++){
                    if(mg_level<=spgrid_level) multigrid_hierarchy(mg_level)->Activate_Cell(spgrid_level,offset,SPGrid_Cell_Type_Dirichlet);
                    else{offset=Flag_array_type::MASK::DownsampleOffset(offset);
                        multigrid_hierarchy(mg_level)->Activate_Cell(mg_level,offset,SPGrid_Cell_Type_Dirichlet);}}}}
    // update stuff based on cell activation
    for(int mg_level=1;mg_level<=mg_levels;mg_level++){
        GRID_HIERARCHY_INITIALIZER<T_STRUCT,T,d>::Flag_Ghost_Cells(*multigrid_hierarchy(mg_level));
        GRID_HIERARCHY_INITIALIZER<T_STRUCT,T,d>::Flag_Valid_Faces(*multigrid_hierarchy(mg_level));
        GRID_HIERARCHY_INITIALIZER<T_STRUCT,T,d>::Flag_Active_Faces(*multigrid_hierarchy(mg_level));
        GRID_HIERARCHY_INITIALIZER<T_STRUCT,T,d>::Flag_Active_Nodes(*multigrid_hierarchy(mg_level));
        GRID_HIERARCHY_INITIALIZER<T_STRUCT,T,d>::Flag_Shared_Nodes(*multigrid_hierarchy(mg_level));
        GRID_HIERARCHY_INITIALIZER<T_STRUCT,T,d>::Generate_Plus_Minus_Active_Faces(*multigrid_hierarchy(mg_level));
        GRID_HIERARCHY_INITIALIZER<T_STRUCT,T,d>::Flag_Active_Cells(*multigrid_hierarchy(mg_level));
        multigrid_hierarchy(mg_level)->Update_Block_Offsets();}
}
//#####################################################################
// Initialize_Boundary_Flags
//#####################################################################
template<class T_STRUCT,class T,int d> void INITIALIZE_MULTIGRID<T_STRUCT,T,d>::
Initialize_Boundary_Flags(const ARRAY<T_HIERARCHY*>& multigrid_hierarchy,const int boundary_radius,const unsigned mask)
{
    // * debug * -- Clearing all instance of 'mask' being set
    for(int i=1;i<=multigrid_hierarchy.m;i++){
        T_HIERARCHY* hierarchy=multigrid_hierarchy(i);
        for(int level=1;level<=hierarchy->Levels();level++){
            Flag_array_type flags=hierarchy->Array(level,&T_STRUCT::flags);
            for(SPGrid_Block_Iterator<typename Flag_array_type::MASK> iterator(hierarchy->Blocks(level));iterator.Valid();iterator.Next()){
                iterator.Data(flags)&=(~mask);}}}
    // * debug *

    // set up acceleration structure
    ARRAY<unsigned long> neighbor_offsets;
    const TV_INT boundary_radius_vector(TV_INT::Constant_Vector(boundary_radius));
    const TV_INT zero_vector(0*TV_INT::All_Ones_Vector());
    for(RANGE_ITERATOR<d> iterator(RANGE<TV_INT>(-boundary_radius_vector,boundary_radius_vector));iterator.Valid();iterator.Next()){
        const TV_INT& index=iterator.Index();
        if(index!=zero_vector) neighbor_offsets.Append(Flag_array_type::MASK::Linear_Offset(std_array<int,d>(index)));}
    // set mask where it is near boundary
    for(int i=1;i<=multigrid_hierarchy.m;i++){
        T_HIERARCHY* hierarchy=multigrid_hierarchy(i);
        for(int level=1;level<=hierarchy->Levels();level++){
            Flag_array_type flags=hierarchy->Array(level,&T_STRUCT::flags);
            for(SPGrid_Block_Iterator<typename Flag_array_type::MASK> iterator(hierarchy->Blocks(level));iterator.Valid();iterator.Next()){
                if(iterator.Data(flags)&SPGrid_Cell_Type_Interior){ // should be 'SPGrid_Cell_Type_Active' ??
                    const unsigned long offset=iterator.Offset();
                    bool mark_cell;
                    for(int k=1;k<=neighbor_offsets.m;k++){
                        const unsigned long neighbor_offset=Flag_array_type::MASK::Packed_Add(offset,neighbor_offsets(k));
                        mark_cell=(!(hierarchy->Set(level).Is_Set(neighbor_offset,SPGrid_Cell_Type_Interior))); // should be 'SPGrid_Cell_Type_Active' ??
                        if(mark_cell) break;}
                    if(mark_cell) iterator.Data(flags)|=mask;}}}}
}
//#####################################################################
template class INITIALIZE_MULTIGRID<FLUIDS_SIMULATION_DATA<float>,float,2>;
template class INITIALIZE_MULTIGRID<FLUIDS_SIMULATION_DATA<float>,float,3>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class INITIALIZE_MULTIGRID<FLUIDS_SIMULATION_DATA<double>,double,2>;
template class INITIALIZE_MULTIGRID<FLUIDS_SIMULATION_DATA<double>,double,3>;
#endif
