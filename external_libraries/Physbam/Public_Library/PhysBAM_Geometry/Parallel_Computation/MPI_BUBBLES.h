//#####################################################################
// Copyright 2012, Mridul Aanjaneya.
// This file is part of PhysBAM whose distribution is governed by the license 
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MPI_BUBBLES
//#####################################################################
#ifndef __MPI_BUBBLES__
#define __MPI_BUBBLES__

#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Utilities/NONCOPYABLE.h>
#include <PhysBAM_Tools/Vectors/VECTOR_ND.h>
namespace MPI{class Group;class Intracomm;class Request;class Status;class Op;class Comm;}
namespace PhysBAM{

template<class TV>
class MPI_BUBBLES:public NONCOPYABLE
{
    typedef typename TV::SCALAR T;

  public:
    int rank;       // global rank
    int number_of_processes,number_of_bubble_processes,bubble_node;
    MPI::Intracomm* comm;
    MPI::Group *group,*bubble_group,*fluid_group;
    VECTOR_ND<int> bubble_ranks,fluid_ranks;
    ARRAY<int> bubble_region;       // mask over global regions, 1 only if region is a bubble, 0 otherwise
    int number_of_global_regions,number_of_global_bubbles;
  private:
    mutable int current_tag;
  public:

    MPI_BUBBLES();
    ~MPI_BUBBLES();

    int Number_Of_Processors() const
    {return number_of_processes;}

    int Get_Unique_Tag() const
    {return current_tag=max(32,(current_tag+1)&((1<<15)-1));}

//#####################################################################
    bool Bubble_Node() const;
    bool Fluid_Node() const;
    T Reduce_Min(const T& local_value) const;
    T Global_Max(const T& local_value) const;
    T Global_Sum(const T& local_value) const;
#ifdef USE_MPI
    int Pack_Size(const ARRAY<ARRAY<T> >& data,const MPI::Comm& comm);
    int Pack_Size(const ARRAY<ARRAY<ARRAY<int> > >& data,const MPI::Comm& comm);
    void Pack(const ARRAY<ARRAY<T> >& data,ARRAY_VIEW<char> buffer,int& position,const MPI::Comm& comm);
    void Pack(const ARRAY<ARRAY<ARRAY<int> > >& data,ARRAY_VIEW<char> buffer,int& position,const MPI::Comm& comm);
    void Unpack(ARRAY<ARRAY<T> >& data,ARRAY_VIEW<const char> buffer,int& position,const MPI::Comm& comm);
    void Unpack(ARRAY<ARRAY<ARRAY<int> > >& data,ARRAY_VIEW<const char> buffer,int& position,const MPI::Comm& comm);
#endif
    void Create_Fluid_Comm_For_Bubble_Nodes() const;
    void Populate_Bubble_Regions();
    void Synchronize_Number_Of_Global_Bubbles();
    void Exchange_Data_Before_Projection(ARRAY<int>& global_bubbles_cell_count,ARRAY<ARRAY<T> >& incompressible_beta_sum,ARRAY<ARRAY<ARRAY<int> > >& incompressible_neighbor_indices,
                                         ARRAY<T>& global_bubbles_divergence,ARRAY<T>& density_target,ARRAY<ARRAY<T> >& incompressible_neighbors_rhs);
    void Fill_Ghost_Cells(VECTOR_ND<T>& v,ARRAY<ARRAY<int> >& neighbor_matrix_indices,ARRAY<int>& bubble_index);
//#####################################################################
};
}
#endif
