//#####################################################################
// Copyright 2010, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <SPGrid_Fluids/Solvers/Inner_Product_Helper.h>

#include <Threading_Tools/PTHREAD_QUEUE.h>
#include <PhysBAM_Tools/Math_Tools/maxabs.h>
#include <SPGrid_Fluids/Simulation/FLUIDS_SIMULATION_FLAGS.h>

using namespace SPGrid;
using namespace PhysBAM;

extern PTHREAD_QUEUE* pthread_queue;

//#####################################################################
// Function Run_Parallel
//#####################################################################
namespace{
template<class T, int elements_per_block>
struct Inner_Product_Helper_Thread_Helper:public PTHREAD_QUEUE::TASK
{
    Inner_Product_Helper<T,elements_per_block>* const obj;
    const int index_start,index_end;
    double& result;

    Inner_Product_Helper_Thread_Helper(Inner_Product_Helper<T,elements_per_block>* const obj_input,const int index_start_input,const int index_end_input,double& result_input)
        :obj(obj_input),index_start(index_start_input),index_end(index_end_input),result(result_input) {}
    void Run(){obj->Run_Index_Range(index_start,index_end,result);}
};
}

template<class T, int elements_per_block> double Inner_Product_Helper<T,elements_per_block>::
Run_Parallel(const int number_of_partitions)
{
    double partial_results[number_of_partitions];

    for(int partition=0;partition<number_of_partitions;partition++){
        // Calculate indicies of current partition
        int first_index_of_partition=(size/number_of_partitions)*(partition)+std::min(size%number_of_partitions,partition);
        int last_index_of_partition=(size/number_of_partitions)*(partition+1)+std::min(size%number_of_partitions,partition+1)-1;
        // Create helper object
        Inner_Product_Helper_Thread_Helper<T,elements_per_block>* task=new Inner_Product_Helper_Thread_Helper<T,elements_per_block>(this,first_index_of_partition,last_index_of_partition,partial_results[partition]);
        // Enqueue
        pthread_queue->Queue(task);
    }

    // Wait for all tasks to complete
    pthread_queue->Wait();

    double result = 0;
    for(int i=0;i<number_of_partitions;i++)
        result += partial_results[i];

    return result;
}
//#####################################################################
// Function Run_Index_Range
//#####################################################################
template<class T, int elements_per_block> 
void Inner_Product_Helper<T,elements_per_block>::
Run_Index_Range(const int index_start,const int index_end,double& result)
{   
    int num_elements = elements_per_block;
    double temp_result = 0;
    
    // Iterating through all block indices
    for(int index=index_start;index<=index_end;index++) 
    {
      T* y1_in  = reinterpret_cast<T*>((unsigned long)y1 + b[index]);
      T* y2_in  = reinterpret_cast<T*>((unsigned long)y2 + b[index]);
      unsigned* flags_in  = reinterpret_cast<unsigned*>((unsigned long)flags + b[index]);

      for(int i=0;i<num_elements;i++){
          if(flags_in[i] & SPGrid_Cell_Type_Interior)
            temp_result += y1_in[i]*y2_in[i];
      }
    }

    result = temp_result;
}
//#####################################################################
template class Inner_Product_Helper<float,16>;
template class Inner_Product_Helper<float,32>;
template class Inner_Product_Helper<float,64>;
template class Inner_Product_Helper<float,128>;
template class Inner_Product_Helper<float,256>;
template class Inner_Product_Helper<float,512>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class Inner_Product_Helper<double,16>;
template class Inner_Product_Helper<double,32>;
template class Inner_Product_Helper<double,64>;
template class Inner_Product_Helper<double,128>;
template class Inner_Product_Helper<double,256>;
template class Inner_Product_Helper<double,512>;
#endif
