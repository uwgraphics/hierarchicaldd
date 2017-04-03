//#####################################################################
// Copyright (c) 2012, Eftychios Sifakis, Sean Bauer
// Distributed under the FreeBSD license (see license.txt)
//#####################################################################
#ifndef __READ_HIERARCHY_h__
#define __READ_HIERARCHY_h__

#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <SPGrid_Fluids/Grids/GRID_HIERARCHY.h>

namespace PhysBAM{
using namespace SPGrid;
//#####################################################################
// Class READ_HIERARCHY
//#####################################################################
template<class T_STRUCT,class T,int d>
class READ_HIERARCHY
{
    typedef float RW;
    typedef GRID_HIERARCHY<T_STRUCT,T,d> T_HIERARCHY;
    typedef GRID<VECTOR<T,d> > T_GRID;
    
// #############################################################################
public:
    static T_HIERARCHY& Read_Hierarchy(const std::string& base_dir,const int frame=0);
    static T_HIERARCHY* Read_Hierarchy_Pointer(const std::string& base_dir,const int frame=0);
// #############################################################################
};
}
#endif
