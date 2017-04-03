//#####################################################################
// Copyright 2005-2012, Mridul Aanjaneya, Geoffrey Irving, Frank Losasso, Michael Lentine, Saket Patkar, Tamar Shinar, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform_Arrays/FACE_ARRAYS.h>
#include <PhysBAM_Tools/Grids_Uniform_PDE_Linear/LAPLACE_UNIFORM.h>
#include <PhysBAM_Tools/Log/DEBUG_UTILITIES.h>
#include <PhysBAM_Geometry/Grids_Uniform_PDE_Linear/BUBBLES_POISSON.h>
#include <PhysBAM_Geometry/Parallel_Computation/LAPLACE_UNIFORM_BUBBLE_MPI.h>
#include <PhysBAM_Geometry/Parallel_Computation/PCG_SPARSE_BUBBLE_MPI.h>
#ifdef USE_MPI
#include <PhysBAM_Tools/Data_Structures/GRAPH.h>
#include <PhysBAM_Tools/Data_Structures/UNION_FIND.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Parallel_Computation/MPI_PACKAGE.h>
#include <PhysBAM_Tools/Parallel_Computation/MPI_UNIFORM_GRID.h>
#include <PhysBAM_Tools/Parallel_Computation/MPI_UTILITIES.h>
#endif
using namespace PhysBAM;

//#####################################################################
// Function Find_Matrix_Indices
//#####################################################################
template<class T_GRID> void LAPLACE_UNIFORM_BUBBLE_MPI<T_GRID>::
Find_Matrix_Indices_Bubbles(ARRAY<int,VECTOR<int,1> >& filled_region_cell_count,ARRAY<int,VECTOR<int,1> >& incompressible_region_cell_count,T_ARRAYS_INT& cell_index_to_matrix_index,ARRAY<ARRAY<TV_INT> >& matrix_index_to_cell_index_array,const VECTOR<T,1>&,BUBBLES_POISSON<T_GRID>* poisson)
{
    assert(local_grid.Is_MAC_Grid());
    int m=local_grid.counts.x;
    Find_Matrix_Indices_In_Region_Bubbles(0,RANGE<VECTOR<int,1> >(1,m),filled_region_cell_count,incompressible_region_cell_count,cell_index_to_matrix_index,matrix_index_to_cell_index_array,poisson);
    Find_Matrix_Indices_In_Region_Bubbles(1,RANGE<VECTOR<int,1> >(0,0),filled_region_cell_count,incompressible_region_cell_count,cell_index_to_matrix_index,matrix_index_to_cell_index_array,poisson);
    Find_Matrix_Indices_In_Region_Bubbles(2,RANGE<VECTOR<int,1> >(m+1,m+1),filled_region_cell_count,incompressible_region_cell_count,cell_index_to_matrix_index,matrix_index_to_cell_index_array,poisson);
    Find_Boundary_Indices_In_Region(1,RANGE<VECTOR<int,1> >(1,1),cell_index_to_matrix_index,poisson);
    Find_Boundary_Indices_In_Region(2,RANGE<VECTOR<int,1> >(m,m),cell_index_to_matrix_index,poisson);
}
//#####################################################################
// Function Find_Matrix_Indices
//#####################################################################
template<class T_GRID> void LAPLACE_UNIFORM_BUBBLE_MPI<T_GRID>::
Find_Matrix_Indices_Bubbles(ARRAY<int,VECTOR<int,1> >& filled_region_cell_count,ARRAY<int,VECTOR<int,1> >& incompressible_region_cell_count,T_ARRAYS_INT& cell_index_to_matrix_index,ARRAY<ARRAY<TV_INT> >& matrix_index_to_cell_index_array,const VECTOR<T,2>&,BUBBLES_POISSON<T_GRID>* poisson)
{
    assert(local_grid.Is_MAC_Grid());
    int m=local_grid.counts.x,n=local_grid.counts.y;
    Find_Matrix_Indices_In_Region_Bubbles(0,RANGE<VECTOR<int,2> >(1,m,1,n),filled_region_cell_count,incompressible_region_cell_count,cell_index_to_matrix_index,matrix_index_to_cell_index_array,poisson);
    Find_Matrix_Indices_In_Region_Bubbles(1,RANGE<VECTOR<int,2> >(0,0,1,n),filled_region_cell_count,incompressible_region_cell_count,cell_index_to_matrix_index,matrix_index_to_cell_index_array,poisson);
    Find_Matrix_Indices_In_Region_Bubbles(2,RANGE<VECTOR<int,2> >(m+1,m+1,1,n),filled_region_cell_count,incompressible_region_cell_count,cell_index_to_matrix_index,matrix_index_to_cell_index_array,poisson);
    Find_Matrix_Indices_In_Region_Bubbles(3,RANGE<VECTOR<int,2> >(1,m,0,0),filled_region_cell_count,incompressible_region_cell_count,cell_index_to_matrix_index,matrix_index_to_cell_index_array,poisson);
    Find_Matrix_Indices_In_Region_Bubbles(4,RANGE<VECTOR<int,2> >(1,m,n+1,n+1),filled_region_cell_count,incompressible_region_cell_count,cell_index_to_matrix_index,matrix_index_to_cell_index_array,poisson);
    Find_Boundary_Indices_In_Region(1,RANGE<VECTOR<int,2> >(1,1,1,n),cell_index_to_matrix_index,poisson);
    Find_Boundary_Indices_In_Region(2,RANGE<VECTOR<int,2> >(m,m,1,n),cell_index_to_matrix_index,poisson);
    Find_Boundary_Indices_In_Region(3,RANGE<VECTOR<int,2> >(1,m,1,1),cell_index_to_matrix_index,poisson);
    Find_Boundary_Indices_In_Region(4,RANGE<VECTOR<int,2> >(1,m,n,n),cell_index_to_matrix_index,poisson);
}
//#####################################################################
// Function Find_Matrix_Indices
//#####################################################################
template<class T_GRID> void LAPLACE_UNIFORM_BUBBLE_MPI<T_GRID>::
Find_Matrix_Indices_Bubbles(ARRAY<int,VECTOR<int,1> >& filled_region_cell_count,ARRAY<int,VECTOR<int,1> >& incompressible_region_cell_count,T_ARRAYS_INT& cell_index_to_matrix_index,ARRAY<ARRAY<TV_INT> >& matrix_index_to_cell_index_array,const VECTOR<T,3>&,BUBBLES_POISSON<T_GRID>* poisson)
{
    assert(local_grid.Is_MAC_Grid());
    int m=local_grid.counts.x,n=local_grid.counts.y,mn=local_grid.counts.z;
    Find_Matrix_Indices_In_Region_Bubbles(0,RANGE<VECTOR<int,3> >(1,m,1,n,1,mn),filled_region_cell_count,incompressible_region_cell_count,cell_index_to_matrix_index,matrix_index_to_cell_index_array,poisson);
    Find_Matrix_Indices_In_Region_Bubbles(1,RANGE<VECTOR<int,3> >(0,0,1,n,1,mn),filled_region_cell_count,incompressible_region_cell_count,cell_index_to_matrix_index,matrix_index_to_cell_index_array,poisson);
    Find_Matrix_Indices_In_Region_Bubbles(2,RANGE<VECTOR<int,3> >(m+1,m+1,1,n,1,mn),filled_region_cell_count,incompressible_region_cell_count,cell_index_to_matrix_index,matrix_index_to_cell_index_array,poisson);
    Find_Matrix_Indices_In_Region_Bubbles(3,RANGE<VECTOR<int,3> >(1,m,0,0,1,mn),filled_region_cell_count,incompressible_region_cell_count,cell_index_to_matrix_index,matrix_index_to_cell_index_array,poisson);
    Find_Matrix_Indices_In_Region_Bubbles(4,RANGE<VECTOR<int,3> >(1,m,n+1,n+1,1,mn),filled_region_cell_count,incompressible_region_cell_count,cell_index_to_matrix_index,matrix_index_to_cell_index_array,poisson);
    Find_Matrix_Indices_In_Region_Bubbles(5,RANGE<VECTOR<int,3> >(1,m,1,n,0,0),filled_region_cell_count,incompressible_region_cell_count,cell_index_to_matrix_index,matrix_index_to_cell_index_array,poisson);
    Find_Matrix_Indices_In_Region_Bubbles(6,RANGE<VECTOR<int,3> >(1,m,1,n,mn+1,mn+1),filled_region_cell_count,incompressible_region_cell_count,cell_index_to_matrix_index,matrix_index_to_cell_index_array,poisson);
    Find_Boundary_Indices_In_Region(1,RANGE<VECTOR<int,3> >(1,1,1,n,1,mn),cell_index_to_matrix_index,poisson);
    Find_Boundary_Indices_In_Region(2,RANGE<VECTOR<int,3> >(m,m,1,n,1,mn),cell_index_to_matrix_index,poisson);
    Find_Boundary_Indices_In_Region(3,RANGE<VECTOR<int,3> >(1,m,1,1,1,mn),cell_index_to_matrix_index,poisson);
    Find_Boundary_Indices_In_Region(4,RANGE<VECTOR<int,3> >(1,m,n,n,1,mn),cell_index_to_matrix_index,poisson);
    Find_Boundary_Indices_In_Region(5,RANGE<VECTOR<int,3> >(1,m,1,n,1,1),cell_index_to_matrix_index,poisson);
    Find_Boundary_Indices_In_Region(6,RANGE<VECTOR<int,3> >(1,m,1,n,mn,mn),cell_index_to_matrix_index,poisson);
}
//#####################################################################
// Function Find_Matrix_Indices_In_Region
//#####################################################################
template<class T_GRID> void LAPLACE_UNIFORM_BUBBLE_MPI<T_GRID>::
Find_Matrix_Indices_In_Region_Bubbles(const int region_index,const RANGE<TV_INT>& region,ARRAY<int,VECTOR<int,1> >& filled_region_cell_count,ARRAY<int,VECTOR<int,1> >& incompressible_region_cell_count,
                                      T_ARRAYS_INT& cell_index_to_matrix_index,ARRAY<ARRAY<TV_INT> >& matrix_index_to_cell_index_array,BUBBLES_POISSON<T_GRID>* poisson)
{
    if(region_index) for(int color=1;color<=filled_region_ranks.m;color++)partitions(color).ghost_indices(region_index).min_corner=filled_region_cell_count(color)+1;
    else for(int color=1;color<=filled_region_ranks.m;color++)partitions(color).interior_indices.min_corner=filled_region_cell_count(color)+1;
    if(region_index){for(CELL_ITERATOR iterator(local_grid,region);iterator.Valid();iterator.Next()){TV_INT c=iterator.Cell_Index();
        int color=filled_region_colors(c);if(color<1 || (!filled_region_touches_dirichlet(color)&&!solve_neumann_regions)) continue;
        int new_index=++filled_region_cell_count(color);cell_index_to_matrix_index(c)=new_index;
        matrix_index_to_cell_index_array(color)(new_index)=c;}}
    else poisson->Compute_Matrix_Indices(region,filled_region_cell_count,incompressible_region_cell_count,matrix_index_to_cell_index_array,cell_index_to_matrix_index);
    if(region_index) for(int color=1;color<=filled_region_ranks.m;color++)partitions(color).ghost_indices(region_index).max_corner=filled_region_cell_count(color);
    else for(int color=1;color<=filled_region_ranks.m;color++)partitions(color).interior_indices.max_corner=incompressible_region_cell_count(color);
}
//#####################################################################
// Function Find_Boundary_Indices_In_Region
//#####################################################################
template<class T_GRID> void LAPLACE_UNIFORM_BUBBLE_MPI<T_GRID>::
Find_Boundary_Indices_In_Region(const int side,const RANGE<TV_INT>& region,T_ARRAYS_INT& cell_index_to_matrix_index,BUBBLES_POISSON<T_GRID>* poisson)
{
    int axis=(side-1)/2+1,cell_side=side&1;
    RANGE<TV_INT> face_region=region;if(!cell_side) face_region+=TV_INT::Axis_Vector(axis);
    // count boundary indices
    ARRAY<int> counts(partitions.m);
    TV_INT face_offset=cell_side?TV_INT():TV_INT::Axis_Vector(axis);
    for(CELL_ITERATOR iterator(local_grid,region);iterator.Valid();iterator.Next())if(!psi_N.Component(axis)(iterator.Cell_Index()+face_offset)){
        int color=filled_region_colors(iterator.Cell_Index());TV_INT cell=iterator.Cell_Index(),ghost_cell;
        for(int axis=1;axis<=T_GRID::dimension;axis++){TV_INT left_cell=cell-TV_INT::Axis_Vector(axis),right_cell=cell+TV_INT::Axis_Vector(axis);
            if(!local_grid.Inside_Domain(left_cell)){ghost_cell=left_cell;break;}
            else if(!local_grid.Inside_Domain(right_cell)){ghost_cell=right_cell;break;}}
        if(counts.Valid_Index(color) && !(phase_color_ghost(cell)==1 && phase_color_ghost(ghost_cell)==1)) counts(color)++;}
    // fill boundary indices
    for(int color=1;color<=partitions.m;color++)partitions(color).boundary_indices(side).Resize(counts(color));
    ARRAYS_COMPUTATIONS::Fill(counts,0);
    for(CELL_ITERATOR iterator(local_grid,region);iterator.Valid();iterator.Next())if(!psi_N.Component(axis)(iterator.Cell_Index()+face_offset)){
        int color=filled_region_colors(iterator.Cell_Index());TV_INT cell=iterator.Cell_Index(),ghost_cell;
        for(int axis=1;axis<=T_GRID::dimension;axis++){TV_INT left_cell=cell-TV_INT::Axis_Vector(axis),right_cell=cell+TV_INT::Axis_Vector(axis);
            if(!local_grid.Inside_Domain(left_cell)){ghost_cell=left_cell;break;}
            else if(!local_grid.Inside_Domain(right_cell)){ghost_cell=right_cell;break;}}
        if(counts.Valid_Index(color) && !(phase_color_ghost(cell)==1 && phase_color_ghost(ghost_cell)==1)) partitions(color).boundary_indices(side)(++counts(color))=cell_index_to_matrix_index(cell);}
}
//#####################################################################
// Function Solve
//#####################################################################
template<class T_GRID> void LAPLACE_UNIFORM_BUBBLE_MPI<T_GRID>::
Solve(SPARSE_MATRIX_FLAT_NXN<T>& A,VECTOR_ND<T>& x,VECTOR_ND<T>& b,VECTOR_ND<T>& q,VECTOR_ND<T>& s,VECTOR_ND<T>& r,VECTOR_ND<T>& k,VECTOR_ND<T>& z,const T tolerance,const int color,ARRAY<ARRAY<int> >& neighbor_matrix_indices,ARRAY<int>& bubble_index)
{
    if(mpi_bubbles->Fluid_Node() && color>filled_region_ranks.m){local_pcg.Solve(A,x,b,q,s,r,k,z,tolerance);return;}
    else{
#ifdef USE_MPI
        MPI::Intracomm* comm=0;
        if(mpi_bubbles->Fluid_Node()) comm=&(*communicators)(color);

        PCG_SPARSE_BUBBLE_MPI<T_GRID> pcg_mpi(local_pcg,*comm,partitions(color),mpi_bubbles,neighbor_matrix_indices,ghost_bubble_indices,bubble_index);
        pcg_mpi.Parallel_Solve(A,x,b,tolerance);
#endif
    }
}
//#####################################################################
#define P(...) __VA_ARGS__
#define INSTANTIATION_HELPER(T_GRID) \
    template void LAPLACE_UNIFORM_BUBBLE_MPI<T_GRID >::Find_Matrix_Indices_Bubbles(ARRAY<int,VECTOR<int,1> >&,ARRAY<int,VECTOR<int,1> >&,T_ARRAYS_INT&,ARRAY<ARRAY<TV_INT> >&,const TV&,BUBBLES_POISSON<T_GRID>*); \
    template void LAPLACE_UNIFORM_BUBBLE_MPI<T_GRID >::Solve(SPARSE_MATRIX_FLAT_NXN<T>&,VECTOR_ND<T>&,VECTOR_ND<T>&,VECTOR_ND<T>&,VECTOR_ND<T>&,VECTOR_ND<T>&,VECTOR_ND<T>&,VECTOR_ND<T>&,const T,const int,ARRAY<ARRAY<int> >&,ARRAY<int>&);
template LAPLACE_UNIFORM_BUBBLE_MPI<GRID<VECTOR<float,1> > >::LAPLACE_UNIFORM_BUBBLE_MPI(LAPLACE_UNIFORM<P(GRID<VECTOR<float,1> >) >&,T_ARRAYS_INT&,MPI_BUBBLES<TV>*);
template LAPLACE_UNIFORM_BUBBLE_MPI<GRID<VECTOR<float,2> > >::LAPLACE_UNIFORM_BUBBLE_MPI(LAPLACE_UNIFORM<P(GRID<VECTOR<float,2> >) >&,T_ARRAYS_INT&,MPI_BUBBLES<TV>*);
template LAPLACE_UNIFORM_BUBBLE_MPI<GRID<VECTOR<float,3> > >::LAPLACE_UNIFORM_BUBBLE_MPI(LAPLACE_UNIFORM<P(GRID<VECTOR<float,3> >) >&,T_ARRAYS_INT&,MPI_BUBBLES<TV>*);
INSTANTIATION_HELPER(P(GRID<VECTOR<float,1> >));
INSTANTIATION_HELPER(P(GRID<VECTOR<float,2> >));
INSTANTIATION_HELPER(P(GRID<VECTOR<float,3> >));
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template LAPLACE_UNIFORM_BUBBLE_MPI<GRID<VECTOR<double,1> > >::LAPLACE_UNIFORM_BUBBLE_MPI(LAPLACE_UNIFORM<P(GRID<VECTOR<double,1> >) >&,T_ARRAYS_INT&,MPI_BUBBLES<TV>*);
template LAPLACE_UNIFORM_BUBBLE_MPI<GRID<VECTOR<double,2> > >::LAPLACE_UNIFORM_BUBBLE_MPI(LAPLACE_UNIFORM<P(GRID<VECTOR<double,2> >) >&,T_ARRAYS_INT&,MPI_BUBBLES<TV>*);
template LAPLACE_UNIFORM_BUBBLE_MPI<GRID<VECTOR<double,3> > >::LAPLACE_UNIFORM_BUBBLE_MPI(LAPLACE_UNIFORM<P(GRID<VECTOR<double,3> >) >&,T_ARRAYS_INT&,MPI_BUBBLES<TV>*);
INSTANTIATION_HELPER(P(GRID<VECTOR<double,1> >));
INSTANTIATION_HELPER(P(GRID<VECTOR<double,2> >));
INSTANTIATION_HELPER(P(GRID<VECTOR<double,3> >));
#endif
