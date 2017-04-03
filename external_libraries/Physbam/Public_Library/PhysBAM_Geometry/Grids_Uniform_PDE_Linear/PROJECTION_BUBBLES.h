//#####################################################################
// Copyright 2012, Mridul Aanjaneya, Saket Patkar.
// This file is part of PhysBAM whose distribution is governed by the license 
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class PROJECTION_BUBBLES
//#####################################################################
#ifndef __PROJECTION_BUBBLES__
#define __PROJECTION_BUBBLES__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Vectors/VECTOR_ND.h>
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Grids_PDE_Linear/LAPLACE.h>
#include <PhysBAM_Tools/Matrices/SPARSE_MATRIX_FLAT_NXN.h>
#include <PhysBAM_Geometry/Parallel_Computation/LAPLACE_UNIFORM_BUBBLE_MPI.h>
namespace PhysBAM{

template<class T_GRID>
class PROJECTION_BUBBLES:public LAPLACE<typename T_GRID::SCALAR>
{
    typedef typename T_GRID::SCALAR T;typedef typename T_GRID::VECTOR_T TV;
    typedef LAPLACE<T> BASE;
  public:
    using BASE::tolerance;

    using BASE::Find_Tolerance;

    ARRAY<T> f_bubble;                                      // right hand side for the bubbles
    ARRAY<T>& density_target;                               // density target per bubble
    ARRAY<int>& cells_per_bubble;    
    ARRAY<ARRAY<T> >& beta_neighbors_per_bubble;            // beta value for each neighbor for each bubble, sum it if the bubble sees it multiple times
    ARRAY<ARRAY<ARRAY<int> > >& neighbor_indices;           // for each bubble, for each incompressible cell, an array of indices of incompressible neighbors
    ARRAY<ARRAY<T> >& f_neighbors_per_bubble;               // right hand side value for each neighbor for each bubble
    LAPLACE_UNIFORM_BUBBLE_MPI<T_GRID>* laplace_bubble_mpi;
    ARRAY<int> bubble_index;                                // matrix indices for all the bubbles
    int number_of_global_bubbles;                           // equivalent to number of regions
    T dt;
    T one_over_dx_squared;                                  // TODO: Assumes uniform grid
    T beta_inc;
    T B;

    PROJECTION_BUBBLES(ARRAY<T>& density_target_input,ARRAY<int>& cells_per_bubble_input,ARRAY<ARRAY<T> >& beta_neighbors_per_bubble_input,ARRAY<ARRAY<ARRAY<int> > >& neighbor_indices_input,
                       ARRAY<ARRAY<T> >& f_neighbors_per_bubble_input,LAPLACE_UNIFORM_BUBBLE_MPI<T_GRID>* laplace_bubble_mpi_input,int number_of_global_bubbles_input,T dt_input,T one_over_dx_squared_input);
    virtual ~PROJECTION_BUBBLES() {}

//#####################################################################
    void Compute_Bubble_RHS(ARRAY<T>& bubble_divergence);
    void Solve(const T time);
    void Find_A(SPARSE_MATRIX_FLAT_NXN<T>& A,VECTOR_ND<T>& b,const int filled_region_cell_count);
    void Find_A_Part_One(ARRAY<int>& row_counts);
    void Find_A_Part_Two(SPARSE_MATRIX_FLAT_NXN<T>& A,VECTOR_ND<T>& b);
//#####################################################################
};
}
#endif
