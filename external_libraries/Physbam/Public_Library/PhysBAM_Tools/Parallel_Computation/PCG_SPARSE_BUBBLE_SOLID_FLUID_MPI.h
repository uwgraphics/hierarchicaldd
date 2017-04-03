//#####################################################################
// Copyright 2012, Mridul Aanjaneya.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class PCG_SPARSE_BUBBLE_SOLID_FLUID_MPI
//#####################################################################
#ifndef __PCG_SPARSE_BUBBLE_SOLID_FLUID_MPI__
#define __PCG_SPARSE_BUBBLE_SOLID_FLUID_MPI__

#ifdef USE_MPI
#include <PhysBAM_Tools/Parallel_Computation/PCG_SPARSE_MPI.h>
#include <PhysBAM_Tools/Parallel_Computation/MPI_BUBBLES_SOLID_FLUID.h>
namespace PhysBAM{

template<class T_GRID>
class PCG_SPARSE_BUBBLE_SOLID_FLUID_MPI:public PCG_SPARSE_MPI<T_GRID>
{
    typedef typename T_GRID::VECTOR_T TV;typedef typename TV::SCALAR T;typedef typename T_GRID::VECTOR_INT TV_INT;
  public:
    typedef PCG_SPARSE_MPI<T_GRID> BASE;

    using BASE::pcg;using BASE::partition;

    using BASE::Initialize_Datatypes;

    MPI_BUBBLES_SOLID_FLUID<T_GRID>* mpi_bubbles;
    int rigid_start_index;
    int rigid_end_index;

    PCG_SPARSE_BUBBLE_SOLID_FLUID_MPI(PCG_SPARSE<T>& pcg_input,MPI::Intracomm& comm_input,SPARSE_MATRIX_PARTITION& partition_input,MPI_BUBBLES_SOLID_FLUID<T_GRID>* mpi_bubbles_input)
        :BASE(pcg_input,comm_input,partition_input),mpi_bubbles(mpi_bubbles_input)
    {}

    void Fill_Ghost_Cells(VECTOR_ND<T>& v)
    {
        BASE::Fill_Ghost_Cells(v);
        VECTOR_ND<T> lambda(rigid_end_index-rigid_start_index+1);
        for(int i=rigid_start_index;i<=rigid_end_index;i++){
            lambda(i-rigid_start_index+1)=v(i);
        }
        //LOG::cout<<"in pcg "<<lambda<<std::endl;
        lambda=mpi_bubbles->Gather_Vector(lambda);
        //LOG::cout<<"in pcg 1"<<lambda<<std::endl;
        for(int i=1;i<=lambda.Size();i++){
            v(i)=lambda(i);
        }
    }
//#####################################################################
    //void Parallel_Solve(SPARSE_MATRIX_FLAT_NXN<T>& A,VECTOR_ND<T>& x,VECTOR_ND<T>& b,const T tolerance=1e-7,const bool recompute_preconditioner=true);
//#####################################################################
};
}
#endif
#endif
