//#####################################################################
// Copyright 2005-2012, Mridul Aanjaneya, Geoffrey Irving, Michael Lentine, Frank Losasso, Saket Patkar, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class LAPLACE_UNIFORM_BUBBLE_MPI
//#####################################################################
#ifndef __LAPLACE_UNIFORM_BUBBLE_MPI__
#define __LAPLACE_UNIFORM_BUBBLE_MPI__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/GRID_ARRAYS_POLICY_UNIFORM.h>
#include <PhysBAM_Tools/Krylov_Solvers/PCG_SPARSE.h>
#include <PhysBAM_Tools/Math_Tools/RANGE.h>
#include <PhysBAM_Tools/Parallel_Computation/LAPLACE_UNIFORM_MPI.h>
#include <PhysBAM_Tools/Parallel_Computation/MPI_UNIFORM_GRID.h>
#include <PhysBAM_Tools/Parallel_Computation/PTHREAD.h>
#include <PhysBAM_Tools/Parallel_Computation/SPARSE_MATRIX_PARTITION.h>
#include <PhysBAM_Geometry/Parallel_Computation/MPI_BUBBLES.h>
namespace PhysBAM{

class GRAPH;
template<class TV> class GRID;
template<class T_GRID> class LAPLACE_UNIFORM;
template<class T_GRID> class BUBBLES_POISSON;

template<class T_GRID>
class LAPLACE_UNIFORM_BUBBLE_MPI:public LAPLACE_UNIFORM_MPI<T_GRID>
{
    typedef typename T_GRID::VECTOR_T TV;typedef typename TV::SCALAR T;
    typedef typename T_GRID::CELL_ITERATOR CELL_ITERATOR;typedef typename MPI_GRID_POLICY<T_GRID>::MPI_GRID T_MPI_GRID;
    typedef typename T_GRID::VECTOR_INT TV_INT;typedef typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_SCALAR T_ARRAYS_SCALAR;
    typedef typename T_ARRAYS_SCALAR::template REBIND<int>::TYPE T_ARRAYS_INT;typedef typename GRID_ARRAYS_POLICY<T_GRID>::FACE_ARRAYS T_FACE_ARRAYS_SCALAR;
    typedef typename T_GRID::FACE_ITERATOR FACE_ITERATOR;
public:
    typedef LAPLACE_UNIFORM_MPI<T_GRID> BASE;
    using BASE::mpi_grid;using BASE::local_grid;using BASE::filled_region_ranks;using BASE::partitions;using BASE::number_of_regions;using BASE::solve_neumann_regions;using BASE::psi_N;
    using BASE::filled_region_touches_dirichlet;using BASE::local_pcg;using BASE::communicators;

    T_ARRAYS_INT& phase_color_ghost;
    MPI_BUBBLES<TV>* mpi_bubbles;
    ARRAY<ARRAY<int>* > ghost_bubble_indices;           // for each bubble, matrix indices for ghost bubbles cells bordering incompressible cells
    bool sum_int_needs_init;
    int sum_int;
#ifdef USE_PTHREADS
    pthread_mutex_t sum_int_lock,lock;
    pthread_barrier_t barr;
#endif

    LAPLACE_UNIFORM_BUBBLE_MPI(LAPLACE_UNIFORM<T_GRID>& laplace,T_ARRAYS_INT& phase_color_ghost_input,MPI_BUBBLES<TV>* mpi_bubbles_input)
        :LAPLACE_UNIFORM_MPI<T_GRID>(laplace),phase_color_ghost(phase_color_ghost_input),mpi_bubbles(mpi_bubbles_input),sum_int_needs_init(false),sum_int(0)
    {
#ifdef USE_PTHREADS
        pthread_mutex_init(&sum_int_lock,0);pthread_mutex_init(&lock,0);
#endif
    }

    void Find_Matrix_Indices(ARRAY<int,VECTOR<int,1> >& filled_region_cell_count,T_ARRAYS_INT& cell_index_to_matrix_index,ARRAY<ARRAY<TV_INT> >& matrix_index_to_cell_index_array)
    {PHYSBAM_FUNCTION_IS_NOT_DEFINED();}

    void Find_Matrix_Indices_Bubbles(ARRAY<int,VECTOR<int,1> >& filled_region_cell_count,ARRAY<int,VECTOR<int,1> >& incompressible_region_cell_count,T_ARRAYS_INT& cell_index_to_matrix_index,ARRAY<ARRAY<TV_INT> >& matrix_index_to_cell_index_array,BUBBLES_POISSON<T_GRID>* poisson)
    {Find_Matrix_Indices_Bubbles(filled_region_cell_count,incompressible_region_cell_count,cell_index_to_matrix_index,matrix_index_to_cell_index_array,TV(),poisson);}

    int Global_Sum(int input)
    {
#ifdef USE_PTHREADS
        pthread_mutex_lock(&sum_int_lock);
        if(sum_int_needs_init){sum_int=input;sum_int_needs_init=false;}
        else sum_int+=input;
        pthread_mutex_unlock(&sum_int_lock);
        pthread_barrier_wait(&barr);
        sum_int_needs_init=true;
        pthread_barrier_wait(&barr);
        return sum_int;
#else
        return input;
#endif
    }
//#####################################################################
    void Solve(SPARSE_MATRIX_FLAT_NXN<T>& A,VECTOR_ND<T>& x,VECTOR_ND<T>& b,VECTOR_ND<T>& q,VECTOR_ND<T>& s,VECTOR_ND<T>& r,VECTOR_ND<T>& k,VECTOR_ND<T>& z,const T tolerance,const int color,ARRAY<ARRAY<int> >& neighbor_matrix_indices,ARRAY<int>& bubble_index);
//#####################################################################
private:
    void Find_Matrix_Indices_Bubbles(ARRAY<int,VECTOR<int,1> >& filled_region_cell_count,ARRAY<int,VECTOR<int,1> >& incompressible_region_cell_count,T_ARRAYS_INT& cell_index_to_matrix_index,ARRAY<ARRAY<TV_INT> >& matrix_index_to_cell_index_array,const VECTOR<T,1>&,BUBBLES_POISSON<T_GRID>* poisson);
    void Find_Matrix_Indices_Bubbles(ARRAY<int,VECTOR<int,1> >& filled_region_cell_count,ARRAY<int,VECTOR<int,1> >& incompressible_region_cell_count,T_ARRAYS_INT& cell_index_to_matrix_index,ARRAY<ARRAY<TV_INT> >& matrix_index_to_cell_index_array,const VECTOR<T,2>&,BUBBLES_POISSON<T_GRID>* poisson);
    void Find_Matrix_Indices_Bubbles(ARRAY<int,VECTOR<int,1> >& filled_region_cell_count,ARRAY<int,VECTOR<int,1> >& incompressible_region_cell_count,T_ARRAYS_INT& cell_index_to_matrix_index,ARRAY<ARRAY<TV_INT> >& matrix_index_to_cell_index_array,const VECTOR<T,3>&,BUBBLES_POISSON<T_GRID>* poisson);
    void Find_Matrix_Indices_In_Region_Bubbles(const int region_index,const RANGE<TV_INT>& region,ARRAY<int,VECTOR<int,1> >& filled_region_cell_count,ARRAY<int,VECTOR<int,1> >& incompressible_region_cell_count,
                                               T_ARRAYS_INT& cell_index_to_matrix_index,ARRAY<ARRAY<TV_INT> >& matrix_index_to_cell_index_array,BUBBLES_POISSON<T_GRID>* poisson);
    void Find_Boundary_Indices_In_Region(const int side,const RANGE<TV_INT>& region,T_ARRAYS_INT& cell_index_to_matrix_index,BUBBLES_POISSON<T_GRID>* poisson);
//#####################################################################
};
}
#endif
