//#####################################################################
// Copyright 2012, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class GRID_HIERARCHY_ITERATOR
//#####################################################################
#include <SPGrid_Fluids/Grids/GRID_HIERARCHY_ITERATOR.h>
#include <Common_Tools/Math_Tools/RANGE_ITERATOR.h>
#include <SPGrid_Fluids/Grids/Rasterizers/HIERARCHICAL_RASTERIZER.h>
#include <SPGrid_Fluids/Simulation/FLUIDS_SIMULATION_DATA.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<int d,class CELL_FUNCTOR> GRID_HIERARCHY_ITERATOR<d,CELL_FUNCTOR>::
GRID_HIERARCHY_ITERATOR(const T_RANGE& range,const int level,CELL_FUNCTOR& functor_input)
    :indicator(VECTOR<unsigned long,d>(range.Edge_Lengths()+1).Product()<<(d*(level-1)),0x100000),functor(functor_input)
{
    for(RANGE_ITERATOR<d> iterator(range);iterator.Valid();iterator.Next())
        stack.Push(T_CELL(level,iterator.Index()));
}
//#####################################################################
// Function Next
//#####################################################################
template<int d,class CELL_FUNCTOR> void GRID_HIERARCHY_ITERATOR<d,CELL_FUNCTOR>::
Next()
{
    T_CELL cell=stack.Pop();
    const int level=cell.x;
    const T_INDEX& index=cell.y;

    bool recurse=functor.Consume(cell);
    if(!recurse) {indicator.Progress(1UL<<(d*(level-1)));return;}

    for(RANGE_ITERATOR<d> iterator(T_RANGE::Unit_Box());iterator.Valid();iterator.Next())
        stack.Push(T_CELL(level-1,2*index-iterator.Index()));
}
//#####################################################################
template class GRID_HIERARCHY_ITERATOR<2,HIERARCHICAL_RASTERIZER<FLUIDS_SIMULATION_DATA<float>,float,2> >;
template class GRID_HIERARCHY_ITERATOR<3,HIERARCHICAL_RASTERIZER<FLUIDS_SIMULATION_DATA<float>,float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class GRID_HIERARCHY_ITERATOR<2,HIERARCHICAL_RASTERIZER<FLUIDS_SIMULATION_DATA<double>,double,2> >;
template class GRID_HIERARCHY_ITERATOR<3,HIERARCHICAL_RASTERIZER<FLUIDS_SIMULATION_DATA<double>,double,3> >;
#endif
