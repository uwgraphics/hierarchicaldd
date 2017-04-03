//#####################################################################
// Copyright 2012, Mridul Aanjaneya.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MPI_BUBBLES
//#####################################################################
#include <PhysBAM_Tools/Log/DEBUG_UTILITIES.h>
#include <PhysBAM_Tools/Log/LOG.h>
#ifdef USE_MPI
#include <PhysBAM_Tools/Parallel_Computation/MPI_PACKAGE.h>
#include <PhysBAM_Tools/Parallel_Computation/MPI_UNIFORM_GRID.h>
#endif
#include "MPI_BUBBLES.h"
using namespace PhysBAM;

#ifdef USE_MPI

//#####################################################################
// Constructor
//#####################################################################
template<class TV> MPI_BUBBLES<TV>::
MPI_BUBBLES()
    :number_of_bubble_processes(1),bubble_node(0),comm(0),group(0),bubble_group(0),fluid_group(0),fluid_ranks(0),bubble_region(0),number_of_global_regions(0),number_of_global_bubbles(0),current_tag(0)
{
    LOG::SCOPE scope("MPI INITIALIZE","Initializing MPI_BUBBLES");
    number_of_processes=MPI::COMM_WORLD.Get_size();
    LOG::cout<<"number of processes = "<<number_of_processes<<std::endl;

    comm=new MPI::Intracomm;
    *comm=MPI::COMM_WORLD.Dup();
    rank=comm->Get_rank();
    LOG::cout<<"global rank = "<<rank<<std::endl;

    group=new MPI::Group(comm->Get_group());
    bubble_ranks.Resize(number_of_bubble_processes);
    for(int i=1;i<=number_of_bubble_processes;i++) bubble_ranks(i)=i-1;
    bubble_group=new MPI::Group(group->Incl(bubble_ranks.n,&bubble_ranks(1)));
    fluid_ranks.Resize(number_of_processes-number_of_bubble_processes);
    for(int i=1;i<=fluid_ranks.n;i++) fluid_ranks(i)=i+number_of_bubble_processes-1;
    fluid_group=new MPI::Group(group->Incl(fluid_ranks.n,&fluid_ranks(1)));
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> MPI_BUBBLES<TV>::
~MPI_BUBBLES()
{
    if(comm){comm->Free();delete comm;}
    if(group){group->Free();delete group;}
    if(bubble_group){bubble_group->Free();delete bubble_group;}
}
//#####################################################################
// Fluid_Node
//#####################################################################
template<class TV> bool MPI_BUBBLES<TV>::
Fluid_Node() const
{
    return fluid_group->Get_rank()!=MPI::UNDEFINED;
}
//#####################################################################
// Bubble_Node
//#####################################################################
template<class TV> bool MPI_BUBBLES<TV>::
Bubble_Node() const
{
    return bubble_group->Get_rank()!=MPI::UNDEFINED;
}
//#####################################################################
// Function Reduce_Min
//#####################################################################
template<class TV> typename TV::SCALAR MPI_BUBBLES<TV>::
Reduce_Min(const T& local_value) const
{
    T global_value;
    MPI_UTILITIES::Reduce(local_value,global_value,MPI::MIN,*comm);
    return global_value;
}
//#####################################################################
// Function Global_Max
//#####################################################################
template<class TV> typename TV::SCALAR MPI_BUBBLES<TV>::
Global_Max(const T& local_value) const
{
    T global_value;
    MPI_UTILITIES::Reduce(local_value,global_value,MPI::MAX,*comm);
    return global_value;
}
//#####################################################################
// Function Global_Sum
//#####################################################################
template<class TV> typename TV::SCALAR MPI_BUBBLES<TV>::
Global_Sum(const T& local_value) const
{
    T global_value;
    MPI_UTILITIES::Reduce(local_value,global_value,MPI::SUM,*comm);
    return global_value;
}
//#####################################################################
// Function Create_Fluid_Comm_For_Solid_Nodes
//#####################################################################
template<class TV> void MPI_BUBBLES<TV>::
Create_Fluid_Comm_For_Bubble_Nodes() const
{
    MPI::COMM_WORLD.Create(*fluid_group);
}
//#####################################################################
// Function Pack_Size
//#####################################################################
template<class TV> int MPI_BUBBLES<TV>::
Pack_Size(const ARRAY<ARRAY<T> >& data,const MPI::Comm& comm)
{int size=MPI_UTILITIES::Pack_Size<int>(comm);for(int i=1;i<=data.m;i++) size+=MPI_UTILITIES::Pack_Size(data(i),comm);return size;}
//#####################################################################
// Function Pack_Size
//#####################################################################
template<class TV> int MPI_BUBBLES<TV>::
Pack_Size(const ARRAY<ARRAY<ARRAY<int> > >& data,const MPI::Comm& comm)
{int size=MPI_UTILITIES::Pack_Size<int>(comm);for(int i=1;i<=data.m;i++) size+=MPI_UTILITIES::Pack_Size(data(i),comm);return size;}
//#####################################################################
// Function Pack
//#####################################################################
template<class TV> void MPI_BUBBLES<TV>::
Pack(const ARRAY<ARRAY<T> >& data,ARRAY_VIEW<char> buffer,int& position,const MPI::Comm& comm)
{assert(Pack_Size(data,comm)<=buffer.Size()-position);
MPI_UTILITIES::Pack(data.m,buffer,position,comm);
for(int i=1;i<=data.m;i++) MPI_UTILITIES::Pack(data(i),buffer,position,comm);}
//#####################################################################
// Function Pack
//#####################################################################
template<class TV> void MPI_BUBBLES<TV>::
Pack(const ARRAY<ARRAY<ARRAY<int> > >& data,ARRAY_VIEW<char> buffer,int& position,const MPI::Comm& comm)
{assert(Pack_Size(data,comm)<=buffer.Size()-position);
MPI_UTILITIES::Pack(data.m,buffer,position,comm);
for(int i=1;i<=data.m;i++) MPI_UTILITIES::Pack(data(i),buffer,position,comm);}
//#####################################################################
// Function Unpack
//#####################################################################
template<class TV> void MPI_BUBBLES<TV>::
Unpack(ARRAY<ARRAY<T> >& data,ARRAY_VIEW<const char> buffer,int& position,const MPI::Comm& comm)
{int m;MPI_UTILITIES::Unpack(m,buffer,position,comm);data.Resize(m);
for(int i=1;i<=data.m;i++) MPI_UTILITIES::Unpack(data(i),buffer,position,comm);}
//#####################################################################
// Function Unpack
//#####################################################################
template<class TV> void MPI_BUBBLES<TV>::
Unpack(ARRAY<ARRAY<ARRAY<int> > >& data,ARRAY_VIEW<const char> buffer,int& position,const MPI::Comm& comm)
{int m;MPI_UTILITIES::Unpack(m,buffer,position,comm);data.Resize(m);
for(int i=1;i<=data.m;i++) MPI_UTILITIES::Unpack(data(i),buffer,position,comm);}
//#####################################################################
// Function Populate_Bubble_Regions
//#####################################################################
template<class TV> void MPI_BUBBLES<TV>::
Populate_Bubble_Regions()
{
    int tag=Get_Unique_Tag();
    if(Fluid_Node()){int buffer_size=MPI_UTILITIES::Pack_Size(number_of_global_regions,*comm)+1;
        buffer_size+=MPI_UTILITIES::Pack_Size(bubble_region,*comm);
        ARRAY<char> buffer(buffer_size);int position=0;
        MPI_UTILITIES::Pack(number_of_global_regions,buffer,position,*comm);
        MPI_UTILITIES::Pack(bubble_region,buffer,position,*comm);
        comm->Send(&buffer(1),buffer_size,MPI::PACKED,bubble_node,tag);}
    else{ARRAY<ARRAY<int> > fluid_bubble_region(fluid_ranks.n);
        for(int i=1;i<=fluid_ranks.n;i++){MPI::Status status;
            comm->Probe(MPI::ANY_SOURCE,tag,status);int source=status.Get_source();
            ARRAY<char> buffer(status.Get_count(MPI::PACKED));int position=0;
            comm->Recv(&buffer(1),buffer.m,MPI::PACKED,source,tag);
            MPI_UTILITIES::Unpack(number_of_global_regions,buffer,position,*comm);
            MPI_UTILITIES::Unpack(fluid_bubble_region(source),buffer,position,*comm);}

        // compute boolean mask for bubbles over all global regions
        number_of_global_bubbles=0;
        bubble_region.Resize(number_of_global_regions);bubble_region.Fill(0);
        for(int i=1;i<=number_of_global_regions;i++){
            for(int j=1;j<=fluid_ranks.n;j++) if(fluid_bubble_region(j)(i)==1) bubble_region(i)=1;
            if(bubble_region(i)==1) number_of_global_bubbles++;}}

    Synchronize_Number_Of_Global_Bubbles();
}
//#####################################################################
// Function Synchronize_Number_Of_Global_Bubbles
//#####################################################################
template<class TV> void MPI_BUBBLES<TV>::
Synchronize_Number_Of_Global_Bubbles()
{
    int tag=Get_Unique_Tag();
    if(Fluid_Node()){MPI::Status status;
        comm->Probe(MPI::ANY_SOURCE,tag,status);int source=status.Get_source();
        ARRAY<char> buffer(status.Get_count(MPI::PACKED));int position=0;
        comm->Recv(&buffer(1),buffer.m,MPI::PACKED,source,tag);
        MPI_UTILITIES::Unpack(number_of_global_bubbles,buffer,position,*comm);}
    else{for(int i=1;i<=fluid_ranks.n;i++){
            int buffer_size=MPI_UTILITIES::Pack_Size(number_of_global_bubbles,*comm)+1;
            ARRAY<char> buffer(buffer_size);int position=0;
            MPI_UTILITIES::Pack(number_of_global_bubbles,buffer,position,*comm);
            comm->Send(&buffer(1),buffer_size,MPI::PACKED,fluid_ranks(i),tag);}}
}
//#####################################################################
// Function Exchange_Data_Before_Projection
//#####################################################################
template<class TV> void MPI_BUBBLES<TV>::
Exchange_Data_Before_Projection(ARRAY<int>& global_bubbles_cell_count,ARRAY<ARRAY<T> >& incompressible_beta_sum,ARRAY<ARRAY<ARRAY<int> > >& incompressible_neighbor_indices,
                                ARRAY<T>& global_bubbles_divergence,ARRAY<T>& density_target,ARRAY<ARRAY<T> >& incompressible_neighbors_rhs)
{
    int tag=Get_Unique_Tag();
    if(Fluid_Node()){int buffer_size=MPI_UTILITIES::Pack_Size(global_bubbles_cell_count,global_bubbles_divergence,density_target,*comm)+1;
        buffer_size+=Pack_Size(incompressible_neighbors_rhs,*comm);buffer_size+=Pack_Size(incompressible_beta_sum,*comm);
        buffer_size+=Pack_Size(incompressible_neighbor_indices,*comm);
        ARRAY<char> buffer(buffer_size);int position=0;
        MPI_UTILITIES::Pack(global_bubbles_cell_count,global_bubbles_divergence,density_target,buffer,position,*comm);
        Pack(incompressible_neighbors_rhs,buffer,position,*comm);Pack(incompressible_beta_sum,buffer,position,*comm);
        Pack(incompressible_neighbor_indices,buffer,position,*comm);
        comm->Send(&buffer(1),buffer_size,MPI::PACKED,bubble_node,tag);}
    else{ARRAY<ARRAY<int> > global_bubbles_cell_counts;global_bubbles_cell_counts.Resize(fluid_ranks.n);
        ARRAY<ARRAY<T> > global_bubbles_divergences,global_bubbles_density_target;
        global_bubbles_density_target.Resize(fluid_ranks.n);global_bubbles_divergences.Resize(fluid_ranks.n);
        ARRAY<ARRAY<ARRAY<T> > > global_bubbles_incompressible_neighbors_rhs;global_bubbles_incompressible_neighbors_rhs.Resize(fluid_ranks.n);
        ARRAY<ARRAY<ARRAY<T> > > global_bubbles_incompressible_beta_sum;global_bubbles_incompressible_beta_sum.Resize(fluid_ranks.n);
        ARRAY<ARRAY<ARRAY<ARRAY<int> > > > global_bubbles_incompressible_neighbor_indices;global_bubbles_incompressible_neighbor_indices.Resize(fluid_ranks.n);
        for(int i=1;i<=fluid_ranks.n;i++){MPI::Status status;
            comm->Probe(MPI::ANY_SOURCE,tag,status);
            int source=status.Get_source();
            ARRAY<char> buffer(status.Get_count(MPI::PACKED));int position=0;
            comm->Recv(&buffer(1),buffer.m,MPI::PACKED,source,tag);
            MPI_UTILITIES::Unpack(global_bubbles_cell_counts(source),global_bubbles_divergences(source),global_bubbles_density_target(source),buffer,position,*comm);
            Unpack(global_bubbles_incompressible_neighbors_rhs(source),buffer,position,*comm);
            Unpack(global_bubbles_incompressible_beta_sum(source),buffer,position,*comm);
            Unpack(global_bubbles_incompressible_neighbor_indices(source),buffer,position,*comm);}
        
        // reassemble stuff
        global_bubbles_cell_count.Resize(number_of_global_bubbles);global_bubbles_cell_count.Fill(0);
        global_bubbles_divergence.Resize(number_of_global_bubbles);global_bubbles_divergence.Fill((T)0.);
        incompressible_neighbors_rhs.Resize(number_of_global_bubbles);incompressible_beta_sum.Resize(number_of_global_bubbles);
        incompressible_neighbor_indices.Resize(number_of_global_bubbles);density_target.Resize(number_of_global_bubbles);int cnt=0,previous=0;
        for(int j=1;j<=number_of_global_regions;j++) if(bubble_region(j)==1){cnt++;
            for(int i=1;i<=fluid_ranks.n;i++){global_bubbles_cell_count(cnt)+=global_bubbles_cell_counts(i)(j);
                global_bubbles_divergence(cnt)+=global_bubbles_divergences(i)(j);
                if(i==1) density_target(cnt)=global_bubbles_density_target(i)(j);
                for(int k=1;k<=global_bubbles_incompressible_neighbors_rhs(i)(j).m;k++)
                    incompressible_neighbors_rhs(cnt).Append(global_bubbles_incompressible_neighbors_rhs(i)(j)(k));
                for(int k=1;k<=global_bubbles_incompressible_beta_sum(i)(j).m;k++)
                    incompressible_beta_sum(cnt).Append(global_bubbles_incompressible_beta_sum(i)(j)(k));
                for(int k=1;k<=global_bubbles_incompressible_neighbor_indices(i)(j).m;k++){
                    ARRAY<int> indices(global_bubbles_incompressible_neighbor_indices(i)(j)(k).m);
                    for(int t=1;t<=global_bubbles_incompressible_neighbor_indices(i)(j)(k).m;t++) indices(t)=global_bubbles_incompressible_neighbor_indices(i)(j)(k)(t)+previous;
                    incompressible_neighbor_indices(cnt).Append(indices);}
                previous+=global_bubbles_incompressible_neighbor_indices(i)(j).m;}}}
}
//#####################################################################
// Function Fill_Ghost_Cells
//#####################################################################
template<class TV> void MPI_BUBBLES<TV>::
Fill_Ghost_Cells(VECTOR_ND<T>& v,ARRAY<ARRAY<int> >& neighbor_matrix_indices,ARRAY<int>& bubble_index)
{
    int tag=Get_Unique_Tag();
    if(Fluid_Node()){ARRAY<ARRAY<T> > neighbor_incompressible_data(neighbor_matrix_indices.m);
        for(int i=1;i<=neighbor_matrix_indices.m;i++){
            for(int j=1;j<=neighbor_matrix_indices(i).m;j++){
                int index=neighbor_matrix_indices(i)(j);
                neighbor_incompressible_data(i).Append(v(index));}}

        // send the incompressible pressures
        int buffer_size=Pack_Size(neighbor_incompressible_data,*comm)+1;
        ARRAY<char> send_buffer(buffer_size);int position=0;
        Pack(neighbor_incompressible_data,send_buffer,position,*comm);
        comm->Send(&send_buffer(1),buffer_size,MPI::PACKED,bubble_node,tag);
        
        // receive and set the bubble pressures
        ARRAY<T> bubble_pressures;MPI::Status status;
        comm->Probe(MPI::ANY_SOURCE,tag,status);int source=status.Get_source();
        ARRAY<char> recv_buffer(status.Get_count(MPI::PACKED));position=0;
        comm->Recv(&recv_buffer(1),recv_buffer.m,MPI::PACKED,source,tag);
        MPI_UTILITIES::Unpack(bubble_pressures,recv_buffer,position,*comm);

        int cnt=0;
        for(int i=1;i<=number_of_global_regions;i++) if(bubble_region(i)==1){cnt++;
            int index=bubble_index(i);v(index)=bubble_pressures(cnt);}}
    else{ARRAY<T> bubble_pressures;
        for(int i=1;i<=bubble_index.m;i++) bubble_pressures.Append(v(bubble_index(i)));

        // send the bubble pressures
        for(int i=1;i<=fluid_ranks.n;i++){
            int buffer_size=MPI_UTILITIES::Pack_Size(bubble_pressures,*comm)+1;
            ARRAY<char> buffer(buffer_size);int position=0;
            MPI_UTILITIES::Pack(bubble_pressures,buffer,position,*comm);
            comm->Send(&buffer(1),buffer_size,MPI::PACKED,fluid_ranks(i),tag);}

        // receive the incompressible pressures
        ARRAY<ARRAY<ARRAY<T> > > global_incompressible_pressures(fluid_ranks.n);
        for(int i=1;i<=fluid_ranks.n;i++){MPI::Status status;
            comm->Probe(MPI::ANY_SOURCE,tag,status);
            int source=status.Get_source();
            ARRAY<char> buffer(status.Get_count(MPI::PACKED));int position=0;
            comm->Recv(&buffer(1),buffer.m,MPI::PACKED,source,tag);
            Unpack(global_incompressible_pressures(source),buffer,position,*comm);}

        int offset=0;
        for(int j=1;j<=number_of_global_regions;j++){
            if(bubble_region(j)==1){
                for(int i=1;i<=fluid_ranks.n;i++){
                    for(int k=1;k<=global_incompressible_pressures(i)(j).m;k++){
                        int index=k+offset;
                        v(index)=global_incompressible_pressures(i)(j)(k);}
                offset+=global_incompressible_pressures(i)(j).Size();}}}}
}
#else
//#####################################################################
template<class TV> MPI_BUBBLES<TV>::MPI_BUBBLES(){PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
template<class TV> MPI_BUBBLES<TV>::~MPI_BUBBLES(){}
template<class TV> bool MPI_BUBBLES<TV>::Fluid_Node() const
{PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
template<class TV> bool MPI_BUBBLES<TV>::Bubble_Node() const
{PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
template<class TV> typename TV::SCALAR MPI_BUBBLES<TV>::Reduce_Min(const T& local_value) const
{PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
template<class TV> typename TV::SCALAR MPI_BUBBLES<TV>::Global_Max(const T& local_value) const
{PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
template<class TV> typename TV::SCALAR MPI_BUBBLES<TV>::Global_Sum(const T& local_value) const
{PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
template<class TV> void MPI_BUBBLES<TV>::Create_Fluid_Comm_For_Bubble_Nodes() const
{PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
template<class TV> void MPI_BUBBLES<TV>::Populate_Bubble_Regions()
{PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
template<class TV> void MPI_BUBBLES<TV>::Synchronize_Number_Of_Global_Bubbles()
{PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
template<class TV> void MPI_BUBBLES<TV>::Exchange_Data_Before_Projection(ARRAY<int>& global_bubbles_cell_count,ARRAY<ARRAY<T> >& incompressible_beta_sum,ARRAY<ARRAY<ARRAY<int> > >& incompressible_neighbor_indices,
                                                                         ARRAY<T>& global_bubbles_divergence,ARRAY<T>& density_target,ARRAY<ARRAY<T> >& incompressible_neighbors_rhs)
{PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
template<class TV> void MPI_BUBBLES<TV>::Fill_Ghost_Cells(VECTOR_ND<T>& v,ARRAY<ARRAY<int> >& neighbor_matrix_indices_input,ARRAY<int>& bubble_index_input)
{PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
//#####################################################################
#endif
template class MPI_BUBBLES<VECTOR<float,1> >;
template class MPI_BUBBLES<VECTOR<float,2> >;
template class MPI_BUBBLES<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class MPI_BUBBLES<VECTOR<double,1> >;
template class MPI_BUBBLES<VECTOR<double,2> >;
template class MPI_BUBBLES<VECTOR<double,3> >;
#endif
