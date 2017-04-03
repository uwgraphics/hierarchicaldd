//#####################################################################
// Copyright (c) 2012, Eftychios Sifakis, Sean Bauer
// Distributed under the FreeBSD license (see license.txt)
//#####################################################################
#ifndef __GRID_HIERARCHY_INITIALIZER_h__
#define __GRID_HIERARCHY_INITIALIZER_h__

#include <SPGrid_Fluids/Grids/GRID_HIERARCHY.h>

namespace PhysBAM{
using namespace SPGrid;
//#####################################################################
// Class GRID_HIERARCHY_INITIALIZER
//#####################################################################
template<class T_STRUCT,class T,int d>
class GRID_HIERARCHY_INITIALIZER
{
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<unsigned>::type Flag_array_type;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<unsigned>::mask Flag_array_mask;
    typedef std_array<int,d> coord_t;
    typedef VECTOR<int,d> T_INDEX;
    typedef GRID_HIERARCHY<T_STRUCT,T,d> Hierarchy_type;

// #############################################################################
public:
    static void Flag_Ghost_Cells(Hierarchy_type& hierarchy);
    static void Flag_Edge_Ghost_Cells(Hierarchy_type& hierarchy);
    static void Flag_Active_Faces(Hierarchy_type& hierarchy);
    static void Flag_Valid_Faces(Hierarchy_type& hierarchy);
    static void Flag_Active_Nodes(Hierarchy_type& hierarchy);
    static void Flag_Shared_Nodes(Hierarchy_type& hierarchy);
    static void Generate_Plus_Minus_Active_Faces(Hierarchy_type& hierarchy);
    static void Flag_Active_Cells(Hierarchy_type& hierarchy);
    static void Flag_T_Junction_Nodes(Hierarchy_type& hierarchy);
    static void Flag_Ghost_Nodes(Hierarchy_type& hierarchy);
// #############################################################################    
};
}
#endif
