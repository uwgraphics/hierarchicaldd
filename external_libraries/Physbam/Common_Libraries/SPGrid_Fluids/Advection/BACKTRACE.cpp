//#####################################################################
// Copyright 2012-2013, Sean Bauer, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BACKTRACE
//#####################################################################
#include <SPGrid_Fluids/Advection/BACKTRACE.h>

#include <PhysBAM_Geometry/Implicit_Objects/IMPLICIT_OBJECT.h>
#include <SPGrid_Fluids/Simulation/FLUIDS_SIMULATION_DATA.h>
#include <SPGrid_Fluids/Grids/GRID_HIERARCHY_LOOKUP.h>
using namespace PhysBAM;
//#####################################################################
// Helper structs
//#####################################################################
namespace{
template<class T_MASK> inline void Add_Small_Offset(const std_array<int,2> d_index,unsigned long& offset)
{
    for(int i=0;i<d_index(0);i++) offset=T_MASK::Packed_OffsetXdim< 1>(offset);
    for(int i=0;i>d_index(0);i--) offset=T_MASK::Packed_OffsetXdim<-1>(offset);
    for(int i=0;i<d_index(1);i++) offset=T_MASK::Packed_OffsetYdim< 1>(offset);
    for(int i=0;i>d_index(1);i--) offset=T_MASK::Packed_OffsetYdim<-1>(offset);
}
template<class T_MASK> inline void Add_Small_Offset(const std_array<int,3> d_index,unsigned long& offset)
{
    for(int i=0;i<d_index(0);i++) offset=T_MASK::Packed_OffsetXdim< 1>(offset);
    for(int i=0;i>d_index(0);i--) offset=T_MASK::Packed_OffsetXdim<-1>(offset);
    for(int i=0;i<d_index(1);i++) offset=T_MASK::Packed_OffsetYdim< 1>(offset);
    for(int i=0;i>d_index(1);i--) offset=T_MASK::Packed_OffsetYdim<-1>(offset);
    for(int i=0;i<d_index(2);i++) offset=T_MASK::Packed_OffsetZdim< 1>(offset);
    for(int i=0;i>d_index(2);i--) offset=T_MASK::Packed_OffsetZdim<-1>(offset);
}
}
//#####################################################################
// Function Backtrace
//#####################################################################
template<class T_STRUCT, class T,int d> VECTOR<T,d> BACKTRACE<T_STRUCT,T,d>::
Backtrace(const T_HIERARCHY& hierarchy,int& level,const TV_INT& cell_index,
    unsigned long& cell_offset,const TV& intra_cell_dX,const TV& dX,
    const RIGID_GEOMETRY_COLLECTION<TV>* const rigid_geometry_collection)
{
    const int original_level=level;
    const unsigned long original_cell_offset=cell_offset;
    // get raw backtraced location index
    TV d_i=intra_cell_dX+dX*hierarchy.Grid(level).one_over_dX; // delta in cell coord = intra_cell_dX + dX/h  (dX = -vel*dt)
    TV_INT d_index;
    for(int v=1;v<=d;v++) d_index(v)=floor(d_i(v));
    TV weights=d_i-TV(d_index);
    TV_INT backtraced_cell_index=cell_index+d_index;

    // clamp index, weights to grid
    const RANGE<TV_INT> cell_indices=hierarchy.Grid(level).Cell_Indices();
    for(int v=1;v<=d;v++){
        if(backtraced_cell_index(v)<cell_indices.min_corner(v)){
            backtraced_cell_index(v)=cell_indices.min_corner(v);
            weights(v)=T();}
        if(backtraced_cell_index(v)>cell_indices.max_corner(v)){
            backtraced_cell_index(v)=cell_indices.max_corner(v);
            weights(v)=(T)1.;}}

    // Look up where this is
    unsigned long backtraced_cell_offset=Flag_array_mask::Linear_Offset(std_array<int,d>(backtraced_cell_index)); // TODO: see if we can remove this
    // if found, then done
    if(GRID_HIERARCHY_LOOKUP<T_STRUCT,T,d>::Cell_Lookup(hierarchy,backtraced_cell_offset,level,weights)){
        cell_offset=backtraced_cell_offset;
        return weights;}

    // REMEMBER, FOR NO OBJECT, WE SHOULD BE DONE HERE

    // if not found then need to project
    TV X=hierarchy.Grid(level).Node(backtraced_cell_index) + weights*hierarchy.Grid(level).dX;    

    int closest_rigid_geometry_index=0;T phi=std::numeric_limits<T>::max();TV normal;
    for(int i=1;i<=rigid_geometry_collection->particles.array_collection->Size();++i){
        const RIGID_GEOMETRY<TV>& rigid_geometry=rigid_geometry_collection->Rigid_Geometry(i);
        T current_phi;
        TV current_normal=rigid_geometry.Implicit_Geometry_Extended_Normal(X,current_phi);
        if(fabs(current_phi)<fabs(phi)){phi=current_phi;closest_rigid_geometry_index=i;normal=current_normal;}}
    PHYSBAM_ASSERT(closest_rigid_geometry_index>0);
    const TV X_projected=X-normal*phi;

    TV_INT projected_cell_index=hierarchy.Grid(level).Cell(X_projected,0);
    weights=(X_projected-hierarchy.Grid(level).Node(projected_cell_index))*hierarchy.Grid(level).one_over_dX;
    weights=clamp(weights,TV(),TV::All_Ones_Vector());
    backtraced_cell_offset=Flag_array_mask::Linear_Offset(std_array<int,d>(projected_cell_index));
    // now cell_lookup should definitely return true
    bool cell_lookup=GRID_HIERARCHY_LOOKUP<T_STRUCT,T,d>::Cell_Lookup(hierarchy,backtraced_cell_offset,level,weights);

    if(!cell_lookup){
        LOG::cout<<"Backtrace FAILED!!! Originating: "<<cell_index<<", Backtraced: "<<backtraced_cell_index<<", Projected: "<<projected_cell_index<<", Phi: "<<phi<<std::endl;
        // Giving up
        cell_offset=original_cell_offset;
        level=original_level;
        return TV();}

    cell_offset=backtraced_cell_offset;
    return weights;
}
//#####################################################################
template class BACKTRACE<FLUIDS_SIMULATION_DATA<float>,float,2>;
template class BACKTRACE<FLUIDS_SIMULATION_DATA<float>,float,3>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class BACKTRACE<FLUIDS_SIMULATION_DATA<double>,double,2>;
template class BACKTRACE<FLUIDS_SIMULATION_DATA<double>,double,3>;
#endif
