//#####################################################################
// Copyright 2012, Mridul Aanjaneya, Saket Patkar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MPI_BUBBLES_SOLID_FLUID
//#####################################################################
#include <PhysBAM_Tools/Log/DEBUG_UTILITIES.h>
#include <PhysBAM_Tools/Log/LOG.h>
#ifdef USE_MPI
#include <PhysBAM_Tools/Parallel_Computation/MPI_PACKAGE.h>
#include <PhysBAM_Tools/Parallel_Computation/MPI_UNIFORM_GRID.h>
#endif
#include "MPI_BUBBLES_SOLID_FLUID.h"
using namespace PhysBAM;

#ifdef USE_MPI
//#####################################################################
// Function Gather_J_Matrix
//#####################################################################
template<class T_GRID> SPARSE_MATRIX_FLAT_MXN<typename T_GRID::SCALAR> MPI_BUBBLES_SOLID_FLUID<T_GRID>::
Gather_J_Matrix(SPARSE_MATRIX_FLAT_MXN<T>& J)
{
    int tag=mpi_grid->Get_Unique_Tag(),rank=mpi_grid->rank;
    int buffer_size=MPI_UTILITIES::Pack_Size(J,*(mpi_grid->comm))+1;
    ARRAY<char> buffer(buffer_size);int position=0;
    MPI_UTILITIES::Pack(J,buffer,position,*(mpi_grid->comm));
    ARRAY<MPI::Request> requests;

    for(int i=1;i<=mpi_grid->number_of_processes;i++) if(i-1!=rank){
        requests.Append(mpi_grid->comm->Isend(&(buffer(1)),position,MPI::PACKED,i-1,tag));}

    ARRAY<ARRAY<char> > recv_buffers(mpi_grid->number_of_processes);
    for(int i=1;i<=mpi_grid->number_of_processes;i++) if(i-1!=rank){MPI::Status status;
        mpi_grid->comm->Probe(i-1,tag,status);
        recv_buffers(i).Resize(status.Get_count(MPI::PACKED));
        requests.Append(mpi_grid->comm->Irecv(&recv_buffers(i)(1),recv_buffers(i).m,MPI::PACKED,i-1,tag));}
    MPI_UTILITIES::Wait_All(requests);

    ARRAY<SPARSE_MATRIX_FLAT_MXN<T> > J_row(mpi_grid->number_of_processes);
    for(int i=1;i<=mpi_grid->number_of_processes;i++) if(i-1!=rank){int position=0;
        MPI_UTILITIES::Unpack(J_row(i),recv_buffers(i),position,*(mpi_grid->comm));}
    
    SPARSE_MATRIX_FLAT_MXN<T> result;result.n=0;result.m=0;
    for(int i=1;i<=mpi_grid->number_of_processes;i++){
        if(i!=(rank+1)){
            if(result.m==0) result=J_row(i);   
            else result=result.Concatenate_Rows(J_row(i));}}
    result=result.Concatenate_Rows(J);
    return result;
}

//#####################################################################
// Function Gather_Vector
//#####################################################################
template<class T_GRID> VECTOR_ND<typename T_GRID::SCALAR> MPI_BUBBLES_SOLID_FLUID<T_GRID>::
Gather_Vector(VECTOR_ND<T>& beta)
{
    int tag=mpi_grid->Get_Unique_Tag(),rank=mpi_grid->rank;
    int buffer_size=MPI_UTILITIES::Pack_Size(beta,*(mpi_grid->comm))+1;
    ARRAY<char> buffer(buffer_size);int position=0;
    MPI_UTILITIES::Pack(beta,buffer,position,*(mpi_grid->comm));
    ARRAY<MPI::Request> requests;

    for(int i=1;i<=mpi_grid->number_of_processes;i++) if(i-1!=rank){
        requests.Append(mpi_grid->comm->Isend(&(buffer(1)),position,MPI::PACKED,i-1,tag));}

    ARRAY<ARRAY<char> > recv_buffers(mpi_grid->number_of_processes);
    for(int i=1;i<=mpi_grid->number_of_processes;i++) if(i-1!=rank){MPI::Status status;
        mpi_grid->comm->Probe(i-1,tag,status);
        recv_buffers(i).Resize(status.Get_count(MPI::PACKED));
        requests.Append(mpi_grid->comm->Irecv(&recv_buffers(i)(1),recv_buffers(i).m,MPI::PACKED,i-1,tag));}
    MPI_UTILITIES::Wait_All(requests);

    ARRAY<VECTOR_ND<T> > beta_row(mpi_grid->number_of_processes);
    for(int i=1;i<=mpi_grid->number_of_processes;i++) if(i-1!=rank){int position=0;
        MPI_UTILITIES::Unpack(beta_row(i),recv_buffers(i),position,*(mpi_grid->comm));}

    int size=beta.Size();
    for(int i=1;i<=mpi_grid->number_of_processes;i++){
        if(i!=(rank+1)){size+=beta_row(i).Size();}}
    
    VECTOR_ND<T> result(size);int count=1;
    for(int i=1;i<=mpi_grid->number_of_processes;i++){
        if(i!=(rank+1)){
            for(int j=1;j<=beta_row(i).Size();j++){
                result(count)=beta_row(i)(j);count++;}}}

    for(int j=1;j<=beta.Size();j++){
        result(count)=beta(j);count++;}
    return result;
}
#else
//#####################################################################
template<class T_GRID> SPARSE_MATRIX_FLAT_MXN<typename T_GRID::SCALAR> MPI_BUBBLES_SOLID_FLUID<T_GRID>::
Gather_J_Matrix(SPARSE_MATRIX_FLAT_MXN<T>& J){PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
template<class T_GRID> VECTOR_ND<typename T_GRID::SCALAR> MPI_BUBBLES_SOLID_FLUID<T_GRID>::
Gather_Vector(VECTOR_ND<T>& beta){PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
//#####################################################################
#endif
template class MPI_BUBBLES_SOLID_FLUID<GRID<VECTOR<float,1> > >;
template class MPI_BUBBLES_SOLID_FLUID<GRID<VECTOR<float,2> > >;
template class MPI_BUBBLES_SOLID_FLUID<GRID<VECTOR<float,3> > >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class MPI_BUBBLES_SOLID_FLUID<GRID<VECTOR<double,1> > >;
template class MPI_BUBBLES_SOLID_FLUID<GRID<VECTOR<double,2> > >;
template class MPI_BUBBLES_SOLID_FLUID<GRID<VECTOR<double,3> > >;
#endif
