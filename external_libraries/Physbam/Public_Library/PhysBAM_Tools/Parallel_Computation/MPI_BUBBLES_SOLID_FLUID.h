//#####################################################################
// Copyright 2012, Mridul Aanjaneya, Saket Patkar.
// This file is part of PhysBAM whose distribution is governed by the license 
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MPI_BUBBLES_SOLID_FLUID
//#####################################################################
#ifndef __MPI_BUBBLES_SOLID_FLUID__
#define __MPI_BUBBLES_SOLID_FLUID__

#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Utilities/NONCOPYABLE.h>
#include <PhysBAM_Tools/Vectors/VECTOR_ND.h>
#include <PhysBAM_Tools/Parallel_Computation/MPI_UNIFORM_GRID.h>
#include <PhysBAM_Tools/Matrices/SPARSE_MATRIX_FLAT_MXN.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_NODE.h>
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Parallel_Computation/MPI_UTILITIES.h>
#include <PhysBAM_Tools/Parallel_Computation/MPI_PACKAGE.h>

namespace MPI{class Group;class Intracomm;class Request;class Status;class Op;class Comm;}
namespace PhysBAM{

template<class T_GRID>
class MPI_BUBBLES_SOLID_FLUID:public NONCOPYABLE
{
   typedef typename T_GRID::SCALAR T;
   typedef typename T_GRID::NODE_ITERATOR NODE_ITERATOR;typedef typename T_GRID::VECTOR_T TV;typedef typename T_GRID::VECTOR_INT TV_INT;

  public:
    MPI_UNIFORM_GRID<T_GRID>* mpi_grid;
    
  public:

    MPI_BUBBLES_SOLID_FLUID():mpi_grid(0){}
    ~MPI_BUBBLES_SOLID_FLUID(){};

//#####################################################################
    SPARSE_MATRIX_FLAT_MXN<T> Gather_J_Matrix(SPARSE_MATRIX_FLAT_MXN<T>& J);
    VECTOR_ND<T> Gather_Vector(VECTOR_ND<T>& beta);
//#####################################################################
};
}
#endif
