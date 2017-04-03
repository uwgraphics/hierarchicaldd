//#####################################################################
// Copyright 2013.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/GRID_ARRAYS_POLICY_UNIFORM.h>
#include <PhysBAM_Tools/Log/DEBUG_UTILITIES.h>
#include <PhysBAM_Tools/Parallel_Computation/MPI_UNIFORM_GRID.h>
#ifdef USE_MPI
#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Parallel_Computation/MPI_PACKAGE.h>
#include <PhysBAM_Tools/Parallel_Computation/MPI_UTILITIES.h>
#endif
namespace PhysBAM{

#ifdef USE_MPI

template<class T_GRID> void
Gather_Clipping_Regions_And_Colors(const MPI_UNIFORM_GRID<T_GRID>& mpi_grid,RANGE<VECTOR<int,2> >& clipping_region,ARRAY<ARRAY<RANGE<VECTOR<int,2> > > >& clipping_region_arrays,ARRAY<VECTOR<typename T_GRID::SCALAR,4>,VECTOR<int,2> >& clipping_color,ARRAY<ARRAY<VECTOR<typename T_GRID::SCALAR,4>,VECTOR<int,2> > >& clipping_colors)
{
    ARRAY<MPI_PACKAGE> packages;ARRAY<MPI::Request> requests;
    ////send
    if(mpi_grid.rank!=0){
        ARRAY<RANGE<VECTOR<int,2> > > clipping_region_array;clipping_region_array.Append(clipping_region);
        MPI_PACKAGE package(clipping_region_array);packages.Append(package);requests.Append(package.Isend(*mpi_grid.comm,0,0));}
    ///receive
    if(mpi_grid.rank==0){
        clipping_region_arrays.Resize(mpi_grid.number_of_processes);
        for(int nb_rank=1;nb_rank<mpi_grid.number_of_processes;nb_rank++){
            clipping_region_arrays(nb_rank).Resize(1);
            MPI_PACKAGE package(clipping_region_arrays(nb_rank));packages.Append(package);requests.Append(package.Irecv(*mpi_grid.comm,nb_rank,0));}}
    //finish
    MPI_UTILITIES::Wait_All(requests);MPI_PACKAGE::Free_All(packages);
        
    ////send
    if(mpi_grid.rank!=0){
        MPI_PACKAGE package(clipping_color,clipping_region);
        packages.Append(package);requests.Append(package.Isend(*mpi_grid.comm,0,0));}
    ////receive
    if(mpi_grid.rank==0){
        clipping_colors.Resize(mpi_grid.number_of_processes);
        for(int nb_rank=1;nb_rank<mpi_grid.number_of_processes;nb_rank++){
            RANGE<VECTOR<int,2> > nb_clipping_region=clipping_region_arrays(nb_rank)(1);
            LOG::cout<<"receive: "<<nb_rank<<": "<<nb_clipping_region.min_corner<<", "<<nb_clipping_region.max_corner<<std::endl;
            clipping_colors(nb_rank).Resize(nb_clipping_region);
            MPI_PACKAGE package(clipping_colors(nb_rank),nb_clipping_region);
            packages.Append(package);requests.Append(package.Irecv(*mpi_grid.comm,nb_rank,0));}}
    ////finish
    MPI_UTILITIES::Wait_All(requests);MPI_PACKAGE::Free_All(packages);
}
//#####################################################################

#else

//#####################################################################
template<class T_GRID> void Gather_Clipping_Regions_And_Colors(const MPI_UNIFORM_GRID<T_GRID>& mpi_grid,RANGE<VECTOR<int,2> >& clipping_region,ARRAY<ARRAY<RANGE<VECTOR<int,2> > > >& clipping_region_arrays,ARRAY<VECTOR<typename T_GRID::SCALAR,4>,VECTOR<int,2> >& clipping_color,ARRAY<ARRAY<VECTOR<typename T_GRID::SCALAR,4>,VECTOR<int,2> > >& clipping_colors){PHYSBAM_NOT_IMPLEMENTED();}

//#####################################################################

#endif

//#####################################################################
#define INSTANTIATION_HELPER(T,T_GRID,d) \
    template void Gather_Clipping_Regions_And_Colors(const MPI_UNIFORM_GRID<T_GRID >&,RANGE<VECTOR<int,2> >&,ARRAY<ARRAY<RANGE<VECTOR<int,2> > > >&,ARRAY<VECTOR<T,4>,VECTOR<int,2> >&,ARRAY<ARRAY<VECTOR<T,4>,VECTOR<int,2> > >&);
#define P(...) __VA_ARGS__
INSTANTIATION_HELPER(float,P(GRID<VECTOR<float,1> >),1);
INSTANTIATION_HELPER(float,P(GRID<VECTOR<float,2> >),2);
INSTANTIATION_HELPER(float,P(GRID<VECTOR<float,3> >),3);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
INSTANTIATION_HELPER(double,P(GRID<VECTOR<double,1> >),1);
INSTANTIATION_HELPER(double,P(GRID<VECTOR<double,2> >),2);
INSTANTIATION_HELPER(double,P(GRID<VECTOR<double,3> >),3);
#endif
}
