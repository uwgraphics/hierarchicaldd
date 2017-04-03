//#####################################################################
// Copyright 2012, Mridul Aanjaneya, Saket Patkar.
// This file is part of PhysBAM whose distribution is governed by the license 
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BUBBLES_POISSON
//#####################################################################
#ifndef __BUBBLES_POISSON__
#define __BUBBLES_POISSON__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Grids_PDE_Linear/POISSON.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Geometry/Grids_Uniform_PDE_Linear/POISSON_COLLIDABLE_UNIFORM.h>
#include <PhysBAM_Geometry/Parallel_Computation/LAPLACE_UNIFORM_BUBBLE_MPI.h>
#include <PhysBAM_Geometry/Parallel_Computation/MPI_BUBBLES.h>
namespace PhysBAM{

template<class T_GRID>
class BUBBLES_POISSON:public POISSON_COLLIDABLE_UNIFORM<T_GRID>
{
    typedef POISSON_COLLIDABLE_UNIFORM<T_GRID> BASE;
    typedef typename T_GRID::VECTOR_T TV;typedef typename TV::SCALAR T;
    typedef typename T_GRID::VECTOR_INT TV_INT;typedef typename T_GRID::CELL_ITERATOR CELL_ITERATOR;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_SCALAR T_ARRAYS_SCALAR;typedef typename T_ARRAYS_SCALAR::template REBIND<int>::TYPE T_ARRAYS_INT;
    typedef typename T_ARRAYS_SCALAR::template REBIND<bool>::TYPE T_ARRAYS_BOOL;
    typedef typename T_ARRAYS_SCALAR::template REBIND<HASHTABLE<int> >::TYPE T_ARRAYS_HASH_INT;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::FACE_ARRAYS T_FACE_ARRAYS_SCALAR;typedef typename T_GRID::FACE_ITERATOR FACE_ITERATOR;
    typedef typename T_GRID::INDEX INDEX;
public:
    using LAPLACE_UNIFORM<T_GRID>::Find_Solution_Regions;using POISSON_UNIFORM<T_GRID>::Initialize_Grid;

    using POISSON_UNIFORM<T_GRID>::grid;using POISSON<T>::use_variable_beta;using POISSON<T>::beta_given_on_faces;
    using BASE::second_order_cut_cell_method;using BASE::second_order_cut_cell_threshold;using BASE::dt;using BASE::u_interface;
    using POISSON<T>::use_weighted_divergence;using POISSON<T>::multiphase;
    using POISSON_UNIFORM<T_GRID>::psi_N;using POISSON_UNIFORM<T_GRID>::periodic_boundary;
    using POISSON_UNIFORM<T_GRID>::filled_region_colors;using POISSON_UNIFORM<T_GRID>::beta_face;
    using POISSON_UNIFORM<T_GRID>::filled_region_touches_dirichlet;using POISSON_UNIFORM<T_GRID>::solve_neumann_regions;
    using LAPLACE_UNIFORM<T_GRID>::thread_queue;using LAPLACE_UNIFORM<T_GRID>::pcg_threaded;using LAPLACE_UNIFORM<T_GRID>::pcg;
    using LAPLACE_UNIFORM<T_GRID>::number_of_regions;using LAPLACE_UNIFORM<T_GRID>::mpi_grid;using LAPLACE_UNIFORM<T_GRID>::laplace_mpi;
    using LAPLACE_UNIFORM<T_GRID>::enforce_compatibility;using LAPLACE_UNIFORM<T_GRID>::tolerance;
#ifdef USE_PTHREADS
    using LAPLACE_UNIFORM<T_GRID>::lock;using LAPLACE_UNIFORM<T_GRID>::barr;
#endif
    using LAPLACE_UNIFORM<T_GRID>::u;

    T_ARRAYS_INT bubble_region_colors;                  // different color for each bubble
    T_ARRAYS_INT phase_color_ghost;                     // phase (air,water,bubble) for each color 
    T_ARRAYS_SCALAR cell_density;
    ARRAY<int>& phase_color;
    int number_of_bubble_regions;
    T_FACE_ARRAYS_SCALAR divergence_face_weights;
    LAPLACE_UNIFORM_BUBBLE_MPI<T_GRID>* laplace_bubble_mpi;
    ARRAY<ARRAY<TV_INT> > neighbor_cell_indices;        // for each bubble ordered list of cell indices for neighbors
    ARRAY<ARRAY<int> > neighbor_matrix_indices;         // for each bubble ordered list of matrix indices for neighbors
    ARRAY<ARRAY<TV_INT>* > ghost_bubble_cells;           // for each bubble, cell indices of ghost bubble cells bordering incompressible cells
    ARRAY<int> bubble_index;
    T B;
public:

    BUBBLES_POISSON(const T_GRID& grid_input,T_ARRAYS_SCALAR& u_input,ARRAY<int>& phase_color_input,MPI_BUBBLES<TV>* mpi_bubbles,const bool initialize_grid,const bool multiphase_input,const bool enforce_compatibility_input)
        :BASE(grid_input,u_input,initialize_grid,multiphase_input,enforce_compatibility_input),phase_color(phase_color_input)
    {laplace_bubble_mpi=new LAPLACE_UNIFORM_BUBBLE_MPI<T_GRID>(*this,phase_color_ghost,mpi_bubbles);laplace_mpi=laplace_bubble_mpi;
     B=(T)82646.818923327897;}
    virtual ~BUBBLES_POISSON() {}

    void Save_Bubble_Regions()
    {number_of_bubble_regions=number_of_regions;
    bubble_region_colors=filled_region_colors;}
//#####################################################################
    void Reconstruct_Compressible_Pressure(RANGE<TV_INT>& domain);
    void Solve(const T time,const bool solution_regions_already_computed);
    void Solve_Subregion(ARRAY<INTERVAL<int> >& interior_indices,ARRAY<ARRAY<INTERVAL<int> > >& ghost_indices,ARRAY<TV_INT>& matrix_index_to_cell_index,SPARSE_MATRIX_FLAT_NXN<T>& A,VECTOR_ND<T>& b,const int color,ARRAY<int,TV_INT>* domain_index,ARRAY<ARRAY<int> >& neighbor_matrix_indices,ARRAY<int>& bubble_index);
    void Compute_Matrix_Indices(ARRAY<int,VECTOR<int,1> >& filled_region_cell_count,ARRAY<int,VECTOR<int,1> >& incompressible_region_cell_count,ARRAY<ARRAY<TV_INT> >& matrix_index_to_cell_index_array,T_ARRAYS_INT& cell_index_to_matrix_index);
    void Compute_Matrix_Indices(const RANGE<TV_INT>& domain,ARRAY<int,VECTOR<int,1> >& filled_region_cell_count,ARRAY<int,VECTOR<int,1> >& incompressible_region_cell_count,ARRAY<ARRAY<TV_INT> >& matrix_index_to_cell_index_array,T_ARRAYS_INT& cell_index_to_matrix_index);
    void Find_A_Part_One(RANGE<TV_INT>& domain,T_ARRAYS_INT& cell_index_to_matrix_index,ARRAY<ARRAY<int> >& row_counts);
    void Find_A_Part_Two(RANGE<TV_INT>& domain,ARRAY<SPARSE_MATRIX_FLAT_NXN<T> >& A_array,ARRAY<VECTOR_ND<T> >& b_array,T_ARRAYS_INT& cell_index_to_matrix_index);
};
}
#endif
