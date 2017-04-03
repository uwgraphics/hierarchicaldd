//#####################################################################
// Copyright 2012, Mridul Aanjaneya, Saket Patkar.
// This file is part of PhysBAM whose distribution is governed by the license 
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Geometry/Grids_Uniform_PDE_Linear/PROJECTION_BUBBLES.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T_GRID> PROJECTION_BUBBLES<T_GRID>::
PROJECTION_BUBBLES(ARRAY<T>& density_target_input,ARRAY<int>& cells_per_bubble_input,ARRAY<ARRAY<T> >& beta_neighbors_per_bubble_input,ARRAY<ARRAY<ARRAY<int> > >& neighbor_indices_input,
                   ARRAY<ARRAY<T> >& f_neighbors_per_bubble_input,LAPLACE_UNIFORM_BUBBLE_MPI<T_GRID>* laplace_bubble_mpi_input,int number_of_global_bubbles_input,T dt_input,T one_over_dx_squared_input):
    density_target(density_target_input),cells_per_bubble(cells_per_bubble_input),beta_neighbors_per_bubble(beta_neighbors_per_bubble_input),neighbor_indices(neighbor_indices_input),
    f_neighbors_per_bubble(f_neighbors_per_bubble_input),laplace_bubble_mpi(laplace_bubble_mpi_input),number_of_global_bubbles(number_of_global_bubbles_input),dt(dt_input),one_over_dx_squared(one_over_dx_squared_input)
{
    beta_inc=(T)1e-3;B=(T)82646.818923327897;
}
//#####################################################################
// Compute_Bubble_RHS
//#####################################################################
template<class T_GRID> void PROJECTION_BUBBLES<T_GRID>::
Compute_Bubble_RHS(ARRAY<T>& bubble_divergence)
{
    f_bubble.Resize(bubble_divergence.m);
    for(int i=1;i<=bubble_divergence.m;i++) f_bubble(i)=(T)-1.*(cells_per_bubble(i)/dt-bubble_divergence(i));
}
//#####################################################################
// Solve
//#####################################################################
template<class T_GRID> void PROJECTION_BUBBLES<T_GRID>::
Solve(const T time)
{
    // compute a single matrix for all the bubbles
    SPARSE_MATRIX_FLAT_NXN<T> A;VECTOR_ND<T> b;
    
    // compute the number of variables
    int filled_region_cell_count=number_of_global_bubbles;
    for(int bubble=1;bubble<=number_of_global_bubbles;bubble++) filled_region_cell_count+=beta_neighbors_per_bubble(bubble).Size();     // all neighbors + number of bubbles
    
    Find_A(A,b,filled_region_cell_count);

    int number_of_unknowns=filled_region_cell_count;
    A.Negate();b*=(T)-1;
    VECTOR_ND<T> x(number_of_unknowns),q,s,r,k,z;
    Find_Tolerance(b); // needs to happen after b is completely set up
    ARRAY<ARRAY<int> > neighbor_matrix_indices(number_of_global_bubbles);int offset=0;
    for(int i=1;i<=number_of_global_bubbles;i++){
        for(int cell=1;cell<=beta_neighbors_per_bubble(i).Size();cell++)
            neighbor_matrix_indices(i).Append(cell+offset);
        offset+=beta_neighbors_per_bubble(i).Size();}
    laplace_bubble_mpi->Solve(A,x,b,q,s,r,k,z,tolerance,1,neighbor_matrix_indices,bubble_index);
}
//#####################################################################
// Find_A
//#####################################################################
template<class T_GRID> void PROJECTION_BUBBLES<T_GRID>::
Find_A(SPARSE_MATRIX_FLAT_NXN<T>& A,VECTOR_ND<T>& b,const int filled_region_cell_count)
{
    // set number of rows to number of variables
    ARRAY<int> row_counts;row_counts.Resize(filled_region_cell_count,false,false);
    b.Resize(filled_region_cell_count);

    // find length of each row, last rows are for bubbles by convention
    Find_A_Part_One(row_counts);
    A.Set_Row_Lengths(row_counts);
    Find_A_Part_Two(A,b);
}
//#####################################################################
// Find_A_Part_One
//#####################################################################
template<class T_GRID> void PROJECTION_BUBBLES<T_GRID>::
Find_A_Part_One(ARRAY<int>& row_counts)
{
    bubble_index.Resize(number_of_global_bubbles,false,false);int offset=0;
    for(int i=1;i<=number_of_global_bubbles;i++){
        // first resize all the incompressible rows
        for(int cell=1;cell<=beta_neighbors_per_bubble(i).Size();cell++)
            row_counts(cell+offset)=2+neighbor_indices(i)(cell).Size();     // diagonal + one for bubble + incompressible neighbors
        offset+=beta_neighbors_per_bubble(i).Size();}

    laplace_bubble_mpi->partitions.Resize(1);
    laplace_bubble_mpi->partitions(1).interior_indices.min_corner=offset+1;
    // now resize the bubble rows
    for(int i=1;i<=number_of_global_bubbles;i++){
        row_counts(++offset)=beta_neighbors_per_bubble(i).Size()+1;         // bubble sees all incompressible neighbors
        bubble_index(i)=offset;}
    laplace_bubble_mpi->partitions(1).interior_indices.max_corner=offset;
}
//#####################################################################
// Find_A_Part_Two
//#####################################################################
template<class T_GRID> void PROJECTION_BUBBLES<T_GRID>::
Find_A_Part_Two(SPARSE_MATRIX_FLAT_NXN<T>& A,VECTOR_ND<T>& b)
{
    int offset=0;
    for(int i=1;i<=number_of_global_bubbles;i++){T bubble_diagonal=(T)0.;
        for(int cell=1;cell<=beta_neighbors_per_bubble(i).Size();cell++){                   // incompressible cells
            b(cell+offset)=f_neighbors_per_bubble(i)(cell);                                 // RHS for the incompressible cell
            T diagonal=(T)0.;
            for(int j=1;j<=neighbor_indices(i)(cell).Size();j++){                           // incompressible neighbors
                T element=beta_inc*one_over_dx_squared;
                A.Set_Element(cell+offset,neighbor_indices(i)(cell)(j),element);            // off-diagonal incompressible entry
                diagonal-=element;}
            T element=beta_neighbors_per_bubble(i)(cell)*one_over_dx_squared;
            A.Set_Symmetric_Elements(cell+offset,bubble_index(i),element);                  // off-diagonal bubble entry for incompressible cell and vice-versa
            diagonal-=element;
            bubble_diagonal-=element;
            A.Set_Element(cell+offset,cell+offset,diagonal);}                               // diagonal entry for the incompressible cell
        b(bubble_index(i))=f_bubble(i);                                                     // RHS for the bubble cell
        bubble_diagonal-=cells_per_bubble(i)/(B*density_target(i)*dt*dt);
        A.Set_Element(bubble_index(i),bubble_index(i),bubble_diagonal);                     // diagonal bubble entry
        offset+=beta_neighbors_per_bubble(i).Size();}
}
template class PROJECTION_BUBBLES<GRID<VECTOR<float,2> > >;
template class PROJECTION_BUBBLES<GRID<VECTOR<float,3> > >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class PROJECTION_BUBBLES<GRID<VECTOR<double,2> > >;
template class PROJECTION_BUBBLES<GRID<VECTOR<double,3> > >;
#endif
