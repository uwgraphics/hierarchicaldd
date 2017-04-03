//#####################################################################
// Copyright 2013, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#include "Ghost_Value_Accumulate.h"
#include <SPGrid_Fluids/Simulation/FLUIDS_SIMULATION_DATA.h>
#include <SPGrid_Fluids/Simulation/FLUIDS_SIMULATION_FLAGS.h>

#include "stdio.h"
#include <map>
#include <utility>

using namespace SPGrid;
using namespace PhysBAM;

template<class T,class T_STRUCT, int d>
struct Ghost_Value_Thread_Accumulate:public PTHREAD_QUEUE::TASK
{
    Ghost_Value_Accumulate<T,T_STRUCT, d>* const obj;
    const int index_start,index_end;
    Ghost_Value_Thread_Accumulate(Ghost_Value_Accumulate<T,T_STRUCT, d>* const obj_input,const int index_start_input,const int index_end_input)
        :obj(obj_input),index_start(index_start_input),index_end(index_end_input) {}
    void Run(){obj->Run_Index_Range(index_start,index_end);}
};

//#####################################################################
// Function Run_Parallel
//#####################################################################
template<class T, class T_STRUCT> void Ghost_Value_Accumulate<T,T_STRUCT,2>::Run_Parallel(const int number_of_partitions)
{
    for(int partition=0;partition<number_of_partitions;partition++){
        // Calculate indicies of current partition
        int first_index_of_partition=(size/number_of_partitions)*(partition)+std::min(size%number_of_partitions,partition);
        int last_index_of_partition=(size/number_of_partitions)*(partition+1)+std::min(size%number_of_partitions,partition+1)-1;
        // Create helper object
        Ghost_Value_Thread_Accumulate<T,T_STRUCT,2>* task=new Ghost_Value_Thread_Accumulate<T,T_STRUCT,2>(this,first_index_of_partition,last_index_of_partition);
        // Enqueue
        pthread_queue->Queue(task);
    }
    
    // Wait for all tasks to complete
    pthread_queue->Wait();
}

//#####################################################################
// Function Run_Index_Range
//#####################################################################
template <class T, class T_STRUCT> void Ghost_Value_Accumulate<T,T_STRUCT,2>::Run_Index_Range(const int index_start,const int index_end)
{   
    // Iterating through all block indices
    for(int index=index_start;index<=index_end;index++) 
    {
        unsigned* coarse_mask_in  = reinterpret_cast<unsigned*>((unsigned long)coarse_mask + b[index]);
        T* coarse_data = reinterpret_cast<T*>((unsigned long)coarse_data_ptr + b[index]);

        unsigned long coarse_offset = b[index];
        for(unsigned int i=0;i<Array_type::MASK::elements_per_block;i++)
        {
            unsigned mask = coarse_mask_in[i];
            if ( mask & (SPGrid_Cell_Type_Ghost|SPGrid_Cell_Type_Active) )
            {
                unsigned long fine_offset = Array_type::MASK::UpsampleOffset(coarse_offset);
                
                if(mask & SPGrid_Ghost_Child_000)
                {
                    unsigned long offset = Array_type::MASK::template Packed_Offset<1,1>(fine_offset);
                    coarse_data[i] += *reinterpret_cast<T*>((unsigned long)fine_data_ptr + offset);
                }
                if(mask & SPGrid_Ghost_Child_010)
                {
                    unsigned long offset = Array_type::MASK::template Packed_Offset<1,0>(fine_offset);
                    coarse_data[i] += *reinterpret_cast<T*>((unsigned long)fine_data_ptr + offset);
                }
                if(mask & SPGrid_Ghost_Child_100)
                {
                    unsigned long offset = Array_type::MASK::template Packed_Offset<0,1>(fine_offset);
                    coarse_data[i] += *reinterpret_cast<T*>((unsigned long)fine_data_ptr + offset);
                }
                if(mask & SPGrid_Ghost_Child_110)
                {
                    unsigned long offset = Array_type::MASK::template Packed_Offset<0,0>(fine_offset);
                    coarse_data[i] += *reinterpret_cast<T*>((unsigned long)fine_data_ptr + offset);
                }
                
            }
            coarse_offset+=sizeof(T);
        }
    }
}

//#####################################################################
// Function Run_Parallel
//#####################################################################
template<class T, class T_STRUCT> void Ghost_Value_Accumulate<T,T_STRUCT,3>::Run_Parallel(const int number_of_partitions)
{
    for(int partition=0;partition<number_of_partitions;partition++){
        // Calculate indicies of current partition
        int first_index_of_partition=(size/number_of_partitions)*(partition)+std::min(size%number_of_partitions,partition);
        int last_index_of_partition=(size/number_of_partitions)*(partition+1)+std::min(size%number_of_partitions,partition+1)-1;
        // Create helper object
        Ghost_Value_Thread_Accumulate<T,T_STRUCT,3>* task=new Ghost_Value_Thread_Accumulate<T,T_STRUCT,3>(this,first_index_of_partition,last_index_of_partition);
        // Enqueue
        pthread_queue->Queue(task);
    }
    
    // Wait for all tasks to complete
    pthread_queue->Wait();
}

//#####################################################################
// Function Run_Index_Range
//#####################################################################
template <class T, class T_STRUCT> void Ghost_Value_Accumulate<T,T_STRUCT,3>::Run_Index_Range(const int index_start,const int index_end)
{   
    // Iterating through all block indices
    for(int index=index_start;index<=index_end;index++) 
    {
        unsigned* coarse_mask_in  = reinterpret_cast<unsigned*>((unsigned long)coarse_mask + b[index]);
        T* coarse_data = reinterpret_cast<T*>((unsigned long)coarse_data_ptr + b[index]);

        unsigned long coarse_offset = b[index];
        for(unsigned int i=0;i<Array_type::MASK::elements_per_block;i++)
        {
            unsigned mask = coarse_mask_in[i];
            if ( mask & (SPGrid_Cell_Type_Ghost|SPGrid_Cell_Type_Active) )
            {
                unsigned long fine_offset = Array_type::MASK::UpsampleOffset(coarse_offset);
                
                if(mask & SPGrid_Ghost_Child_000)
                {
                    unsigned long offset = Array_type::MASK::template Packed_Offset<1,1,1>(fine_offset);
                    coarse_data[i] += *reinterpret_cast<T*>((unsigned long)fine_data_ptr + offset);
                }
                if(mask & SPGrid_Ghost_Child_001)
                {
                    unsigned long offset = Array_type::MASK::template Packed_Offset<1,1,0>(fine_offset);
                    coarse_data[i] += *reinterpret_cast<T*>((unsigned long)fine_data_ptr + offset);
                }
                if(mask & SPGrid_Ghost_Child_010)
                {
                    unsigned long offset = Array_type::MASK::template Packed_Offset<1,0,1>(fine_offset);
                    coarse_data[i] += *reinterpret_cast<T*>((unsigned long)fine_data_ptr + offset);
                }
                if(mask & SPGrid_Ghost_Child_011)
                {
                    unsigned long offset = Array_type::MASK::template Packed_Offset<1,0,0>(fine_offset);
                    coarse_data[i] += *reinterpret_cast<T*>((unsigned long)fine_data_ptr + offset);
                }
                if(mask & SPGrid_Ghost_Child_100)
                {
                    unsigned long offset = Array_type::MASK::template Packed_Offset<0,1,1>(fine_offset);
                    coarse_data[i] += *reinterpret_cast<T*>((unsigned long)fine_data_ptr + offset);
                }
                if(mask & SPGrid_Ghost_Child_101)
                {
                    unsigned long offset = Array_type::MASK::template Packed_Offset<0,1,0>(fine_offset);
                    coarse_data[i] += *reinterpret_cast<T*>((unsigned long)fine_data_ptr + offset);
                }
                if(mask & SPGrid_Ghost_Child_110)
                {
                    unsigned long offset = Array_type::MASK::template Packed_Offset<0,0,1>(fine_offset);
                    coarse_data[i] += *reinterpret_cast<T*>((unsigned long)fine_data_ptr + offset);
                }
                if(mask & SPGrid_Ghost_Child_111)
                {
                    unsigned long offset = Array_type::MASK::template Packed_Offset<0,0,0>(fine_offset);
                    coarse_data[i] += *reinterpret_cast<T*>((unsigned long)fine_data_ptr + offset);
                }
                
            }
            coarse_offset+=sizeof(T);
        }
    }
}


template class Ghost_Value_Accumulate<unsigned,FLUIDS_SIMULATION_DATA<float>,2>;
template class Ghost_Value_Accumulate<unsigned,FLUIDS_SIMULATION_DATA<float>,3>;
template class Ghost_Value_Accumulate<float,FLUIDS_SIMULATION_DATA<float>,2>;
template class Ghost_Value_Accumulate<float,FLUIDS_SIMULATION_DATA<float>,3>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class Ghost_Value_Accumulate<unsigned,FLUIDS_SIMULATION_DATA<double>,2>;
template class Ghost_Value_Accumulate<unsigned,FLUIDS_SIMULATION_DATA<double>,3>;
template class Ghost_Value_Accumulate<double,FLUIDS_SIMULATION_DATA<double>,2>;
template class Ghost_Value_Accumulate<double,FLUIDS_SIMULATION_DATA<double>,3>;
#endif
