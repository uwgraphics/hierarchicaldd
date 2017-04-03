//#####################################################################
// Copyright 2013, Sean Bauer, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include "Smoother_Helper.h"
#include <SPGrid_Fluids/Simulation/FLUIDS_SIMULATION_FLAGS.h>
#include <SPGrid_Fluids/Simulation/FLUIDS_SIMULATION_DATA.h>

using namespace SPGrid;
using namespace PhysBAM;

template<class T,int log2_struct,int d>
struct Smoother_Thread_Helper:public PTHREAD_QUEUE::TASK
{
    Smoother_Helper<T,log2_struct,d>* const obj;
    const int index_start,index_end;
    Smoother_Thread_Helper(Smoother_Helper<T,log2_struct,d>* const obj_input,const int index_start_input,const int index_end_input)
        :obj(obj_input),index_start(index_start_input),index_end(index_end_input) {}
    void Run(){obj->Run_Index_Range(index_start,index_end);}
};

//#####################################################################
// Function Run_Parallel
//#####################################################################
template<class T, int log2_struct> void Smoother_Helper<T,log2_struct,2>::Run_Parallel(const int number_of_partitions)
{
    for(int partition=0;partition<number_of_partitions;partition++){
        // Calculate indicies of current partition
        int first_index_of_partition=(size/number_of_partitions)*(partition)+std::min(size%number_of_partitions,partition);
        int last_index_of_partition=(size/number_of_partitions)*(partition+1)+std::min(size%number_of_partitions,partition+1)-1;
        // Create helper object
        Smoother_Thread_Helper<T,log2_struct,d>* task=new Smoother_Thread_Helper<T,log2_struct,d>(this,first_index_of_partition,last_index_of_partition);
        // Enqueue
        pthread_queue->Queue(task);
    }
    // Wait for all tasks to complete
    pthread_queue->Wait();
}

template <class T, int log2_struct>
void Smoother_Helper<T,log2_struct,2>::ComputeShadowGrid( unsigned long* offset_grid_ptr, 
                                                         unsigned* mask_in, 
                                                         unsigned long packed_offset)
{
    typedef unsigned long (&offset_grid_type)[og_xsize][og_ysize];
    offset_grid_type o_grid = reinterpret_cast<offset_grid_type>(*offset_grid_ptr);
    
    unsigned long simple_offset = 0;
    // Fill in simple offsets
    for (int i = xmin; i<=xmax; i++)
    for (int j = ymin; j<=ymax; j++)
    {
      o_grid[i][j] = packed_offset + simple_offset;  // Can do simple addition here since addresses are within block
      simple_offset += sizeof(T);
    }

    // First let's do the starting points
    o_grid[xmin-1][ymin] =  T_MASK::template Packed_OffsetXdim<-1>(o_grid[xmin][ymin]);
    o_grid[xmax+1][ymin] =  T_MASK::template Packed_OffsetXdim< 1>(o_grid[xmax][ymin]);

    o_grid[xmin][ymin-1] =  T_MASK::template Packed_OffsetYdim<-1>(o_grid[xmin][ymin]);
    o_grid[xmin][ymax+1] =  T_MASK::template Packed_OffsetYdim< 1>(o_grid[xmin][ymax]);

    // Fill in edge offsets (cube faces, but not edges will be correct after this)
    // This is ok for 6 neighbors, but one more pass will be needed for kernels that use edges
    {
      // Left and Right face
      for (int i=xmin-1; i<=xmax+1; i+= (xmax-xmin)+2)
      {
        simple_offset = o_grid[i][ymin];
        for (int j=ymin; j<=ymax; j++)
        {
          o_grid[i][j] = simple_offset;
          simple_offset += sizeof(T);  // Simple addition (going through neighboring block in same manner)
        }
      }
    }
    
    {
      // Top and bottom face
      for (int j = ymin-1; j<=ymax+1; j+= (ymax-ymin)+2)
      {
        simple_offset = o_grid[xmin][j];
        for (int i=xmin; i<=xmax; i++)
        {
          o_grid[i][j] = simple_offset;
          simple_offset += sizeof(T) * (block_ysize);
        }
      }
    }

}
//#####################################################################
// Function Run_Index_Range
//#####################################################################
// T_MASK corresponds to the mask for the data (not the mask channel)
template <class T, int log2_struct> void Smoother_Helper<T,log2_struct,2>::Run_Index_Range(const int index_start,const int index_end)
{
    // Compute shadow grid of linear offsets
    unsigned long* offset_grid_ptr = (unsigned long*)malloc( (og_xsize) * (og_ysize) * sizeof(unsigned long));
    typedef unsigned long (&offset_grid_type)[og_xsize][og_ysize];
    offset_grid_type o_grid = reinterpret_cast<offset_grid_type>(*offset_grid_ptr);
    
    // Iterating through all block indices
    for(int index=index_start;index<=index_end;index++) {
      T* x_ptr     = reinterpret_cast<T*>((unsigned long)x + b[index]);
      T* rhs_ptr   = reinterpret_cast<T*>((unsigned long)rhs + b[index]);
      T* delta_ptr = reinterpret_cast<T*>((unsigned long)delta + b[index]);
      T* dinv_ptr  = reinterpret_cast<T*>((unsigned long)dinv + b[index]);

      unsigned* mask_in  = reinterpret_cast<unsigned*>((unsigned long)mask + b[index]);
      unsigned long packed_offset = (unsigned long)x_ptr - (unsigned long)x;
      unsigned long x_base_addr = reinterpret_cast<unsigned long>(x);
      unsigned long mask_base_addr = reinterpret_cast<unsigned long>(mask);

      ComputeShadowGrid(offset_grid_ptr, mask_in, packed_offset);

      int cur_index = 0;
      // Actually process elements
      for(int i=xmin;i<=xmax;i++)
      for(int j=ymin;j<=ymax;j++)
      {
        unsigned mask_value = mask_in[cur_index];

        if (mask_value & active_flag_mask){
            T2 result=(T2)0.;
            const T& dinv_i = dinv_ptr[cur_index];

            result = rhs_ptr[cur_index]*dinv_i-x_ptr[cur_index];

            // - x
            if((*reinterpret_cast<unsigned*>(mask_base_addr + o_grid[i-1][j])) & SPGrid_Cell_Type_Active)
                result -= laplace_scale_uniform*(*reinterpret_cast<T*>(x_base_addr + o_grid[i-1][j]))*dinv_i;
            else if((*reinterpret_cast<unsigned*>(mask_base_addr + o_grid[i-1][j])) & SPGrid_Cell_Type_Ghost)
                result -= laplace_scale_nonuniform*(*reinterpret_cast<T*>(x_base_addr + o_grid[i-1][j]))*dinv_i;

            // + x
            if((*reinterpret_cast<unsigned*>(mask_base_addr + o_grid[i+1][j])) & SPGrid_Cell_Type_Active)
                result -= laplace_scale_uniform*(*reinterpret_cast<T*>(x_base_addr + o_grid[i+1][j]))*dinv_i;
            else if((*reinterpret_cast<unsigned*>(mask_base_addr + o_grid[i+1][j])) & SPGrid_Cell_Type_Ghost)
                result -= laplace_scale_nonuniform*(*reinterpret_cast<T*>(x_base_addr + o_grid[i+1][j]))*dinv_i;

            // - y
            if((*reinterpret_cast<unsigned*>(mask_base_addr + o_grid[i][j-1])) & SPGrid_Cell_Type_Active)
                result -= laplace_scale_uniform*(*reinterpret_cast<T*>(x_base_addr + o_grid[i][j-1]))*dinv_i;
            else if((*reinterpret_cast<unsigned*>(mask_base_addr + o_grid[i][j-1])) & SPGrid_Cell_Type_Ghost)
                result -= laplace_scale_nonuniform*(*reinterpret_cast<T*>(x_base_addr + o_grid[i][j-1]))*dinv_i;

            // + y
            if((*reinterpret_cast<unsigned*>(mask_base_addr + o_grid[i][j+1])) & SPGrid_Cell_Type_Active)
                result -= laplace_scale_uniform*(*reinterpret_cast<T*>(x_base_addr + o_grid[i][j+1]))*dinv_i;
            else if((*reinterpret_cast<unsigned*>(mask_base_addr + o_grid[i][j+1])) & SPGrid_Cell_Type_Ghost)
                result -= laplace_scale_nonuniform*(*reinterpret_cast<T*>(x_base_addr + o_grid[i][j+1]))*dinv_i;

            delta_ptr[cur_index] = (T)(result*(T2)omega);
        }
        cur_index++;
      }
    }
    free(offset_grid_ptr);
}

//#####################################################################
// Function Run_Parallel
//#####################################################################
template<class T, int log2_struct> void Smoother_Helper<T,log2_struct,3>::Run_Parallel(const int number_of_partitions)
{
    for(int partition=0;partition<number_of_partitions;partition++){
        // Calculate indicies of current partition
        int first_index_of_partition=(size/number_of_partitions)*(partition)+std::min(size%number_of_partitions,partition);
        int last_index_of_partition=(size/number_of_partitions)*(partition+1)+std::min(size%number_of_partitions,partition+1)-1;
        // Create helper object
        Smoother_Thread_Helper<T,log2_struct,d>* task=new Smoother_Thread_Helper<T,log2_struct,d>(this,first_index_of_partition,last_index_of_partition);
        // Enqueue
        pthread_queue->Queue(task);
    }
    // Wait for all tasks to complete
    pthread_queue->Wait();
}


template <class T, int log2_struct>
void Smoother_Helper<T,log2_struct,3>::ComputeShadowGrid( unsigned long* offset_grid_ptr, 
                                                                            unsigned* mask_in, 
                                                                            unsigned long packed_offset)
{
    typedef unsigned long (&offset_grid_type)[og_xsize][og_ysize][og_zsize];
    offset_grid_type o_grid = reinterpret_cast<offset_grid_type>(*offset_grid_ptr);
    
    unsigned long simple_offset = 0;
    // Fill in simple offsets
    for (int i = xmin; i<=xmax; i++)
    for (int j = ymin; j<=ymax; j++)
    for (int k = zmin; k<=zmax; k++)
    {
      o_grid[i][j][k] = packed_offset + simple_offset;  // Can do simple addition here since addresses are within block
      simple_offset += sizeof(T);
    }

    // First let's do the starting points
    o_grid[xmin-1][ymin][zmin] =  T_MASK::template Packed_OffsetXdim<-1>(o_grid[xmin][ymin][zmin]);
    o_grid[xmax+1][ymin][zmin] =  T_MASK::template Packed_OffsetXdim< 1>(o_grid[xmax][ymin][zmin]);

    o_grid[xmin][ymin][zmin-1] =  T_MASK::template Packed_OffsetZdim<-1>(o_grid[xmin][ymin][zmin]);
    o_grid[xmin][ymin][zmax+1] =  T_MASK::template Packed_OffsetZdim< 1>(o_grid[xmin][ymin][zmax]);

    o_grid[xmin][ymin-1][zmin] =  T_MASK::template Packed_OffsetYdim<-1>(o_grid[xmin][ymin][zmin]);
    o_grid[xmin][ymax+1][zmin] =  T_MASK::template Packed_OffsetYdim< 1>(o_grid[xmin][ymax][zmin]);

    // Fill in edge offsets (cube faces, but not edges will be correct after this)
    // This is ok for 6 neighbors, but one more pass will be needed for kernels that use edges
    {
      // Left and Right face
      for (int i=xmin-1; i<=xmax+1; i+= (xmax-xmin)+2)
      {
        simple_offset = o_grid[i][ymin][zmin];
        for (int j=ymin; j<=ymax; j++)
        for (int k=zmin; k<=zmax; k++)
        {
          o_grid[i][j][k] = simple_offset;
          simple_offset += sizeof(T);  // Simple addition (going through neighboring block in same manner)
        }
      }
    }

    {
      // Front and Back face
      for (int k=zmin-1; k<=zmax+1; k+= (zmax-zmin)+2)
      {
        simple_offset = o_grid[xmin][ymin][k];
        for (int i=xmin; i<=xmax; i++)
        for (int j=ymin; j<=ymax; j++)
        {
          o_grid[i][j][k] = simple_offset;
          simple_offset += block_zsize*sizeof(T);  
        }
      }
    }
    
    {
      // Top and bottom face
      for (int j=ymin-1; j<=ymax+1; j+= (ymax-ymin)+2)
      {
        simple_offset = o_grid[xmin][j][zmin];
        for (int i=xmin; i<=xmax; i++)
        {
          for (int k=zmin; k<=zmax; k++)
          {
            o_grid[i][j][k] = simple_offset;
            simple_offset += sizeof(T);  
          }
          simple_offset += sizeof(T) * (block_ysize-1) * (block_zsize);
        }
      }
    }

}

//#####################################################################
// Function Run_Index_Range
//#####################################################################
// T_MASK corresponds to the mask for the data (not the mask channel)
template <class T, int log2_struct> void Smoother_Helper<T,log2_struct,3>::Run_Index_Range(const int index_start,const int index_end)
{
    // Compute shadow grid of linear offsets
    unsigned long* offset_grid_ptr = (unsigned long*)malloc( (og_xsize) * (og_ysize) * (og_zsize) * sizeof(unsigned long));
    typedef unsigned long (&offset_grid_type)[og_xsize][og_ysize][og_zsize];
    offset_grid_type o_grid = reinterpret_cast<offset_grid_type>(*offset_grid_ptr);
    
    // Iterating through all block indices
    for(int index=index_start;index<=index_end;index++) {
      T* x_ptr     = reinterpret_cast<T*>((unsigned long)x + b[index]);
      T* rhs_ptr   = reinterpret_cast<T*>((unsigned long)rhs + b[index]);
      T* delta_ptr = reinterpret_cast<T*>((unsigned long)delta + b[index]);
      T* dinv_ptr  = reinterpret_cast<T*>((unsigned long)dinv + b[index]);

      unsigned* mask_in  = reinterpret_cast<unsigned*>((unsigned long)mask + b[index]);
      unsigned long packed_offset = (unsigned long)x_ptr - (unsigned long)x;
      unsigned long x_base_addr = reinterpret_cast<unsigned long>(x);
      unsigned long mask_base_addr = reinterpret_cast<unsigned long>(mask);

      ComputeShadowGrid(offset_grid_ptr, mask_in, packed_offset);

      int cur_index = 0;
      // Actually process elements
      for(int i=xmin;i<=xmax;i++)
      for(int j=ymin;j<=ymax;j++)
      for(int k=zmin;k<=zmax;k++)
      {
        unsigned mask_value = mask_in[cur_index];

        if (mask_value & active_flag_mask){
            T result=(T)0.;
            const T dinv_i = dinv_ptr[cur_index];

            result = rhs_ptr[cur_index]*dinv_i-x_ptr[cur_index];

            // - x
            if((*reinterpret_cast<unsigned*>(mask_base_addr + o_grid[i-1][j][k])) & SPGrid_Cell_Type_Active)
                result -= laplace_scale_uniform*(*reinterpret_cast<T*>(x_base_addr + o_grid[i-1][j][k]))*dinv_i;
            else if((*reinterpret_cast<unsigned*>(mask_base_addr + o_grid[i-1][j][k])) & SPGrid_Cell_Type_Ghost)
                result -= laplace_scale_nonuniform*(*reinterpret_cast<T*>(x_base_addr + o_grid[i-1][j][k]))*dinv_i;

            // + x
            if((*reinterpret_cast<unsigned*>(mask_base_addr + o_grid[i+1][j][k])) & SPGrid_Cell_Type_Active)
                result -= laplace_scale_uniform*(*reinterpret_cast<T*>(x_base_addr + o_grid[i+1][j][k]))*dinv_i;
            else if((*reinterpret_cast<unsigned*>(mask_base_addr + o_grid[i+1][j][k])) & SPGrid_Cell_Type_Ghost)
                result -= laplace_scale_nonuniform*(*reinterpret_cast<T*>(x_base_addr + o_grid[i+1][j][k]))*dinv_i;

            // - y
            if((*reinterpret_cast<unsigned*>(mask_base_addr + o_grid[i][j-1][k])) & SPGrid_Cell_Type_Active)
                result -= laplace_scale_uniform*(*reinterpret_cast<T*>(x_base_addr + o_grid[i][j-1][k]))*dinv_i;
            else if((*reinterpret_cast<unsigned*>(mask_base_addr + o_grid[i][j-1][k])) & SPGrid_Cell_Type_Ghost)
                result -= laplace_scale_nonuniform*(*reinterpret_cast<T*>(x_base_addr + o_grid[i][j-1][k]))*dinv_i;

            // + y
            if((*reinterpret_cast<unsigned*>(mask_base_addr + o_grid[i][j+1][k])) & SPGrid_Cell_Type_Active)
                result -= laplace_scale_uniform*(*reinterpret_cast<T*>(x_base_addr + o_grid[i][j+1][k]))*dinv_i;
            else if((*reinterpret_cast<unsigned*>(mask_base_addr + o_grid[i][j+1][k])) & SPGrid_Cell_Type_Ghost)
                result -= laplace_scale_nonuniform*(*reinterpret_cast<T*>(x_base_addr + o_grid[i][j+1][k]))*dinv_i;

            // - z
            if((*reinterpret_cast<unsigned*>(mask_base_addr + o_grid[i][j][k-1])) & SPGrid_Cell_Type_Active)
                result -= laplace_scale_uniform*(*reinterpret_cast<T*>(x_base_addr + o_grid[i][j][k-1]))*dinv_i;
            else if((*reinterpret_cast<unsigned*>(mask_base_addr + o_grid[i][j][k-1])) & SPGrid_Cell_Type_Ghost)
                result -= laplace_scale_nonuniform*(*reinterpret_cast<T*>(x_base_addr + o_grid[i][j][k-1]))*dinv_i;

            // + z
            if((*reinterpret_cast<unsigned*>(mask_base_addr + o_grid[i][j][k+1])) & SPGrid_Cell_Type_Active)
                result -= laplace_scale_uniform*(*reinterpret_cast<T*>(x_base_addr + o_grid[i][j][k+1]))*dinv_i;
            else if((*reinterpret_cast<unsigned*>(mask_base_addr + o_grid[i][j][k+1])) & SPGrid_Cell_Type_Ghost)
                result -= laplace_scale_nonuniform*(*reinterpret_cast<T*>(x_base_addr + o_grid[i][j][k+1]))*dinv_i;

            result*=omega;
            delta_ptr[cur_index] = result;
        }
        cur_index++;
      }
    }
    free(offset_grid_ptr);
}

template class Smoother_Helper<float,5,2>;
template class Smoother_Helper<float,6,2>;

template class Smoother_Helper<float,5,3>;
template class Smoother_Helper<float,6,3>;
//template class Smoother_Helper<float,5,2>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class Smoother_Helper<double,5,2>;
template class Smoother_Helper<double,6,2>;

template class Smoother_Helper<double,5,3>;
template class Smoother_Helper<double,6,3>;
#endif
