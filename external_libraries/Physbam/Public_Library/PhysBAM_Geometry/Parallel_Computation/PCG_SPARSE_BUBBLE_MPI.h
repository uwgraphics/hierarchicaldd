//#####################################################################
// Copyright 2012, Mridul Aanjaneya.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class PCG_SPARSE_BUBBLE_MPI
//#####################################################################
#ifndef __PCG_SPARSE_BUBBLE_MPI__
#define __PCG_SPARSE_BUBBLE_MPI__

#ifdef USE_MPI
#include <PhysBAM_Tools/Parallel_Computation/PCG_SPARSE_MPI.h>
#include <PhysBAM_Geometry/Parallel_Computation/MPI_BUBBLES.h>
namespace PhysBAM{

template<class T_GRID>
class PCG_SPARSE_BUBBLE_MPI:public PCG_SPARSE_MPI<T_GRID>
{
    typedef typename T_GRID::VECTOR_T TV;typedef typename TV::SCALAR T;typedef typename T_GRID::VECTOR_INT TV_INT;
  public:
    typedef PCG_SPARSE_MPI<T_GRID> BASE;

    using BASE::pcg;using BASE::partition;

    using BASE::Initialize_Datatypes;

    MPI_BUBBLES<TV>* mpi_bubbles;
    ARRAY<ARRAY<int> >& neighbor_matrix_indices;
    ARRAY<ARRAY<int>* > ghost_bubble_indices;
    ARRAY<int>& bubble_index;

    PCG_SPARSE_BUBBLE_MPI(PCG_SPARSE<T>& pcg_input,MPI::Intracomm& comm_input,SPARSE_MATRIX_PARTITION& partition_input,MPI_BUBBLES<TV>* mpi_bubbles_input,ARRAY<ARRAY<int> >& neighbor_matrix_indices_input,ARRAY<ARRAY<int>* > ghost_bubble_indices_input,ARRAY<int>& bubble_index_input)
        :BASE(pcg_input,comm_input,partition_input),mpi_bubbles(mpi_bubbles_input),neighbor_matrix_indices(neighbor_matrix_indices_input),ghost_bubble_indices(ghost_bubble_indices_input),bubble_index(bubble_index_input)
    {}

    T Global_Max_Bubbles(const T& input)
    {   
        if(mpi_bubbles->number_of_global_bubbles>0) return mpi_bubbles->Global_Max(input);
        return BASE::Global_Max(input);
    }

    T Global_Sum_Bubbles(const T& input)
    {
        if(mpi_bubbles->number_of_global_bubbles>0) return mpi_bubbles->Global_Sum(input);
        return BASE::Global_Sum(input);
    }

    void Fill_Ghost_Cells(VECTOR_ND<T>& v)
    {
        if(mpi_bubbles->Fluid_Node()) BASE::Fill_Ghost_Cells(v);
        if(mpi_bubbles->number_of_global_bubbles>0) mpi_bubbles->Fill_Ghost_Cells(v,neighbor_matrix_indices,bubble_index);

        // fill the ghost bubble indices with the correct bubble pressures
        for(int i=1;i<=ghost_bubble_indices.Size();i++){
            if(ghost_bubble_indices(i) && ghost_bubble_indices(i)->Size()>0){ARRAY<int>& ghost_indices=*ghost_bubble_indices(i);
                for(int j=1;j<=ghost_bubble_indices(i)->Size();j++) v(ghost_indices(j))=v(bubble_index(i));}}
    }
//#####################################################################
    void Parallel_Solve(SPARSE_MATRIX_FLAT_NXN<T>& A,VECTOR_ND<T>& x,VECTOR_ND<T>& b,const T tolerance=1e-7,const bool recompute_preconditioner=true);
//#####################################################################
};
}
#endif
#endif
