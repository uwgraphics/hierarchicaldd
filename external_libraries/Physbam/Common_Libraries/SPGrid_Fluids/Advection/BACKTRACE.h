//#####################################################################
// Copyright 2012-2013, Sean Bauer, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BACKTRACE
//#####################################################################
#ifndef __BACKTRACE_h__
#define __BACKTRACE_h__

#include <SPGrid_Fluids/Grids/GRID_HIERARCHY.h>
#include <PhysBAM_Geometry/Solids_Geometry/RIGID_GEOMETRY_COLLECTION.h>
#include <PhysBAM_Geometry/Solids_Geometry/RIGID_GEOMETRY.h>

namespace PhysBAM{

template<class T_STRUCT,class T,int d>
class BACKTRACE
{
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<unsigned>::type Flag_array_type;
    typedef typename Flag_array_type::MASK Flag_array_mask;
    typedef GRID_HIERARCHY<T_STRUCT,T,d> T_HIERARCHY;
    typedef VECTOR<T,d> TV;
    typedef VECTOR<int,d> TV_INT;

//#####################################################################
public:
    static VECTOR<T,d> Backtrace(const T_HIERARCHY& hierarchy,int& level,const VECTOR<int,d>& cell_index,
        unsigned long& cell_offset,const VECTOR<T,d>& intra_cell_dX,const VECTOR<T,d>& dX,
        const RIGID_GEOMETRY_COLLECTION<TV>* const rigid_geometry_collection);
//#####################################################################
};
}
#endif
