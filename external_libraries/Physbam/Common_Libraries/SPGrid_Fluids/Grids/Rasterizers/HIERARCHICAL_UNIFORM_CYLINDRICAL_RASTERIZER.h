//#####################################################################
// Copyright 2016, Mridul Aanjaneya.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class HIERARCHICAL_UNIFORM_CYLINDRICAL_RASTERIZER
//#####################################################################
#ifndef __HIERARCHICAL_UNIFORM_CYLINDRICAL_RASTERIZER__
#define __HIERARCHICAL_UNIFORM_CYLINDRICAL_RASTERIZER__

#include <Common_Tools/Math_Tools/RANGE_ITERATOR.h>
#include <PhysBAM_Geometry/Basic_Geometry/CYLINDER.h>
#include <SPGrid_Fluids/Grids/Rasterizers/HIERARCHICAL_RASTERIZER.h>
#include <SPGrid_Fluids/Grids/GRID_HIERARCHY.h>
#include <SPGrid_Fluids/Simulation/FLUIDS_SIMULATION_FLAGS.h>

using namespace SPGrid;

namespace PhysBAM{
template<class T_STRUCT,class T,int d>
class HIERARCHICAL_UNIFORM_CYLINDRICAL_RASTERIZER:public HIERARCHICAL_RASTERIZER<T_STRUCT,T,d>
{
    typedef HIERARCHICAL_RASTERIZER<T_STRUCT,T,d> BASE;
    typedef VECTOR<T,d> TV;typedef VECTOR<int,d> T_INDEX;
    typedef RANGE<TV> T_SCALAR_RANGE;typedef PAIR<int,T_INDEX> T_CELL;

public:
    using BASE::hierarchy;
    const ARRAY<CYLINDER<T>*> cylinders;

    HIERARCHICAL_UNIFORM_CYLINDRICAL_RASTERIZER(GRID_HIERARCHY<T_STRUCT,T,d>& hierarchy_input,const ARRAY<CYLINDER<T>*>& cylinders_input)
        :BASE(hierarchy_input),cylinders(cylinders_input)
    {}

    T Phi(const TV& location)
    {
        T phi=(T)1e10;
        for(int i=1;i<=cylinders.Size();++i) phi=min(phi,cylinders(i)->Signed_Distance(location));
        return phi;
    }

    bool Consume(const T_CELL& cell)
    {
        const int level=cell.x;
        const T_INDEX& index=cell.y;
        
        TV center=hierarchy.Grid(level).Center(index);

        const int cell_size=1<<(level-1);
        RANGE<T_INDEX> fine_indices((index-1)*cell_size+1,index*cell_size);
        PHYSBAM_ASSERT(hierarchy.Grid(1).Domain_Indices().Lazy_Inside(fine_indices.min_corner));
        PHYSBAM_ASSERT(hierarchy.Grid(1).Domain_Indices().Lazy_Inside(fine_indices.max_corner));

        //T_INDEX nodes[GRID<TV>::number_of_nodes_per_cell];
        //hierarchy.Grid(level).Nodes_In_Cell_From_Minimum_Corner_Node(index,nodes);
        //T phi[GRID<TV>::number_of_nodes_per_cell];
        //for(int i=0;i<GRID<TV>::number_of_nodes_per_cell;++i) phi[i]=Phi(hierarchy.Grid(level).Node(nodes[i]));
        //bool outside=true;
        //for(int i=0;i<GRID<TV>::number_of_nodes_per_cell;++i) if(phi[i]<root_three*(T).5*hierarchy.Grid(level).dX.Min()){outside=false;break;}
        //if(outside) return false;

        if(level==1 && Phi(center)<=(T)0.) hierarchy.Activate_Cell(level,index,SPGrid_Cell_Type_Interior);
        return (level>1);
    }
//#####################################################################
};
}
#endif
