//#####################################################################
// Copyright 2013.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MPI_RENDERING
//#####################################################################
#ifndef __MPI_RENDERING__
#define __MPI_RENDERING__

#include <PhysBAM_Tools/Arrays/ARRAYS_FORWARD.h>
#include <PhysBAM_Tools/Data_Structures/DATA_STRUCTURES_FORWARD.h>
#include <PhysBAM_Tools/Parallel_Computation/MPI_UNIFORM_GRID.h>

namespace PhysBAM{

template<class T_GRID> void Gather_Clipping_Regions_And_Colors(const MPI_UNIFORM_GRID<T_GRID>& mpi_grid,RANGE<VECTOR<int,2> >& clipping_region,ARRAY<ARRAY<RANGE<VECTOR<int,2> > > >& clipping_region_arrays,ARRAY<VECTOR<typename T_GRID::SCALAR,4>,VECTOR<int,2> >& clipping_color,ARRAY<ARRAY<VECTOR<typename T_GRID::SCALAR,4>,VECTOR<int,2> > >& clipping_colors);


}
#endif
