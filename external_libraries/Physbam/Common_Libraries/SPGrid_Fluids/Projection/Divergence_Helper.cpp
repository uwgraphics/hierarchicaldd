//#####################################################################
// Copyright 2013, Sean Bauer, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include "Divergence_Helper.h"
#include "Laplace_Helper.h"
#include <SPGrid_Fluids/Simulation/FLUIDS_SIMULATION_DATA.h>
#include <SPGrid_Fluids/Simulation/FLUIDS_SIMULATION_FLAGS.h>

using namespace SPGrid;
using namespace PhysBAM;

template<class T,int log2_struct,int d>
struct Divergence_Thread_Helper:public PTHREAD_QUEUE::TASK
{
    Divergence_Helper<T,log2_struct,d>* const obj;
    const int index_start,index_end;
    Divergence_Thread_Helper(Divergence_Helper<T,log2_struct,d>* const obj_input,const int index_start_input,const int index_end_input)
        :obj(obj_input),index_start(index_start_input),index_end(index_end_input) {}
    void Run(){obj->Run_Index_Range(index_start,index_end);}
};

//#####################################################################
// Function Run_Parallel
//#####################################################################
template<class T, int log2_struct> void Divergence_Helper<T,log2_struct,2>::Run_Parallel(const int number_of_partitions)
{
    
    ///* DEBUG */ std::cout<<"Divergence_Helper::Run_Parallel"<<std::endl;
    
    for(int partition=0;partition<number_of_partitions;partition++){
        // Calculate indicies of current partition
        int first_index_of_partition=(size/number_of_partitions)*(partition)+std::min(size%number_of_partitions,partition);
        int last_index_of_partition=(size/number_of_partitions)*(partition+1)+std::min(size%number_of_partitions,partition+1)-1;
        // Create helper object
        Divergence_Thread_Helper<T,log2_struct,d>* task=new Divergence_Thread_Helper<T,log2_struct,d>(this,first_index_of_partition,last_index_of_partition);
        // Enqueue
        pthread_queue->Queue(task);
    }
    // Wait for all tasks to complete
    pthread_queue->Wait();
}

//#####################################################################
// Function Run_Index_Range
//#####################################################################
// T_MASK corresponds to the mask for the data (not the mask channel)
template <class T, int log2_struct> void Divergence_Helper<T,log2_struct,2>::Run_Index_Range(const int index_start,const int index_end)
{  
   // Compute shadow grid of linear offsets
    unsigned long* offset_grid_ptr = (unsigned long*)malloc( (og_xsize) * (og_ysize) * sizeof(unsigned long));
    typedef unsigned long (&offset_grid_type)[og_xsize][og_ysize];
    offset_grid_type o_grid = reinterpret_cast<offset_grid_type>(*offset_grid_ptr);
    
    // Iterating through all block indices
    for(int index=index_start;index<=index_end;index++) {
      T* x_face_input = reinterpret_cast<T*>((unsigned long)x_faces + b[index]);
      T* y_face_input = reinterpret_cast<T*>((unsigned long)y_faces + b[index]);
      T* cell_data_out  = reinterpret_cast<T*>((unsigned long)cell_data + b[index]);

      unsigned* mask_in  = reinterpret_cast<unsigned*>((unsigned long)mask + b[index]);
      unsigned long packed_offset = b[index];
      unsigned long x_face_base_addr = reinterpret_cast<unsigned long>(x_faces);
      unsigned long y_face_base_addr = reinterpret_cast<unsigned long>(y_faces);

      Laplace_Helper<T,log2_struct,d>::ComputeShadowGrid(offset_grid_ptr, mask_in, packed_offset);

      int cur_index = 0;
      // Actually process elements
      for(int i=xmin;i<=xmax;i++)
      for(int j=ymin;j<=ymax;j++)
      {
        unsigned mask_value = mask_in[cur_index];

        if ( mask_value & (SPGrid_Cell_Type_Interior|SPGrid_Cell_Type_Ghost)) {
            double result=(double)0;
            const double x_face_input_current = (double)(x_face_input[cur_index]);
            const double y_face_input_current = (double)(y_face_input[cur_index]);
#if 0
            if (mask_value & SPGrid_Face_Minus_X_Active)
                result -= scale * x_face_input[cur_index];
            if (mask_value & SPGrid_Face_Plus_X_Active)
                result += scale * (*reinterpret_cast<T*>(x_face_base_addr + o_grid[i+1][j]));

            if (mask_value & SPGrid_Face_Minus_Y_Active)
                result -= scale * y_face_input[cur_index];
            if (mask_value & SPGrid_Face_Plus_Y_Active)
                result += scale * (*reinterpret_cast<T*>(y_face_base_addr + o_grid[i][j+1]));
#else
            unsigned long mask_base_addr = reinterpret_cast<unsigned long>(mask);
            unsigned mask_value_minus_x = mask_in[cur_index];
            unsigned mask_value_minus_y = mask_in[cur_index];
            // careful!! -- TODO: fix
            unsigned mask_value_plus_x  = (*reinterpret_cast<unsigned*>(mask_base_addr + o_grid[i+1][j]));
            unsigned mask_value_plus_y  = (*reinterpret_cast<unsigned*>(mask_base_addr + o_grid[i][j+1]));
            // careful!! -- TODO: fix
            const double x_face_input_current_plus = (double)(*reinterpret_cast<T*>(x_face_base_addr + o_grid[i+1][j]));
            const double y_face_input_current_plus = (double)(*reinterpret_cast<T*>(y_face_base_addr + o_grid[i][j+1]));
            
            // if (mask_value_minus_x & SPGrid_Face_Type_X_Valid)
            //     result -= scale * x_face_input[cur_index];
            // if (mask_value_plus_x  & SPGrid_Face_Type_X_Valid)
            //     result += scale * (*reinterpret_cast<T*>(x_face_base_addr + o_grid[i+1][j]));
                                                        
            // if (mask_value_minus_y & SPGrid_Face_Type_Y_Valid)
            //     result -= scale * y_face_input[cur_index];
            // if (mask_value_plus_y  & SPGrid_Face_Type_Y_Valid)
            //     result += scale * (*reinterpret_cast<T*>(y_face_base_addr + o_grid[i][j+1]));                                                        

            if (mask_value_minus_x & SPGrid_Face_Type_X_Valid)
                result -= scale * x_face_input_current;
            if (mask_value_plus_x  & SPGrid_Face_Type_X_Valid)
                result += scale * x_face_input_current_plus;

            if (mask_value_minus_y & SPGrid_Face_Type_Y_Valid)
                result -= scale * y_face_input_current;
            if (mask_value_plus_y  & SPGrid_Face_Type_Y_Valid)
                result += scale * y_face_input_current_plus;
#endif
            
            cell_data_out[cur_index] = result;
        }
        cur_index++;
      }
    }
    free(offset_grid_ptr);
}

//#####################################################################
// Function Run_Parallel
//#####################################################################
template<class T, int log2_struct> void Divergence_Helper<T,log2_struct,3>::Run_Parallel(const int number_of_partitions)
{
    for(int partition=0;partition<number_of_partitions;partition++){
        // Calculate indicies of current partition
        int first_index_of_partition=(size/number_of_partitions)*(partition)+std::min(size%number_of_partitions,partition);
        int last_index_of_partition=(size/number_of_partitions)*(partition+1)+std::min(size%number_of_partitions,partition+1)-1;
        // Create helper object
        Divergence_Thread_Helper<T,log2_struct,d>* task=new Divergence_Thread_Helper<T,log2_struct,d>(this,first_index_of_partition,last_index_of_partition);
        // Enqueue
        pthread_queue->Queue(task);
    }
    // Wait for all tasks to complete
    pthread_queue->Wait();
}

//#####################################################################
// Function Run_Index_Range
//#####################################################################
// T_MASK corresponds to the mask for the data (not the mask channel)
template <class T, int log2_struct> void Divergence_Helper<T,log2_struct,3>::Run_Index_Range(const int index_start,const int index_end)
{  

    // Compute shadow grid of linear offsets
    unsigned long* offset_grid_ptr = (unsigned long*)malloc( (og_xsize) * (og_ysize) * (og_zsize) * sizeof(unsigned long));
    typedef unsigned long (&offset_grid_type)[og_xsize][og_ysize][og_zsize];
    offset_grid_type o_grid = reinterpret_cast<offset_grid_type>(*offset_grid_ptr);
    
    // Iterating through all block indices
    for(int index=index_start;index<=index_end;index++) {
      T* x_face_input = reinterpret_cast<T*>((unsigned long)x_faces + b[index]);
      T* y_face_input = reinterpret_cast<T*>((unsigned long)y_faces + b[index]);
      T* z_face_input = reinterpret_cast<T*>((unsigned long)z_faces + b[index]);
      T* cell_data_out  = reinterpret_cast<T*>((unsigned long)cell_data + b[index]);

      unsigned* mask_in  = reinterpret_cast<unsigned*>((unsigned long)mask + b[index]);
      unsigned long packed_offset = b[index];
      unsigned long x_face_base_addr = reinterpret_cast<unsigned long>(x_faces);
      unsigned long y_face_base_addr = reinterpret_cast<unsigned long>(y_faces);
      unsigned long z_face_base_addr = reinterpret_cast<unsigned long>(z_faces);

      Laplace_Helper<T,log2_struct,d>::ComputeShadowGrid(offset_grid_ptr, mask_in, packed_offset);

      int cur_index = 0;
      // Actually process elements
      for(int i=xmin;i<=xmax;i++)
      for(int j=ymin;j<=ymax;j++)
      for(int k=zmin;k<=zmax;k++)
      {
        unsigned mask_value = mask_in[cur_index];

        if ( mask_value & (SPGrid_Cell_Type_Interior|SPGrid_Cell_Type_Ghost)) {
            T result = 0;

            const double x_face_input_current = (double)(x_face_input[cur_index]);
            const double y_face_input_current = (double)(y_face_input[cur_index]);
            const double z_face_input_current = (double)(z_face_input[cur_index]);
#if 0
            if (mask_value & SPGrid_Face_Minus_X_Active)
                result -= scale * x_face_input[cur_index];
            if (mask_value & SPGrid_Face_Plus_X_Active)
                result += scale * (*reinterpret_cast<T*>(x_face_base_addr + o_grid[i+1][j][k]));

            if (mask_value & SPGrid_Face_Minus_Y_Active)
                result -= scale * y_face_input[cur_index];
            if (mask_value & SPGrid_Face_Plus_Y_Active)
                result += scale * (*reinterpret_cast<T*>(y_face_base_addr + o_grid[i][j+1][k]));

            if (mask_value & SPGrid_Face_Minus_Z_Active)
                result -= scale * z_face_input[cur_index];
            if (mask_value & SPGrid_Face_Plus_Z_Active)
                result += scale * (*reinterpret_cast<T*>(z_face_base_addr + o_grid[i][j][k+1]));
#else
            unsigned long mask_base_addr = reinterpret_cast<unsigned long>(mask);
            unsigned mask_value_minus_x = mask_in[cur_index];
            unsigned mask_value_minus_y = mask_in[cur_index];
            unsigned mask_value_minus_z = mask_in[cur_index];
            // careful!! -- TODO: fix
            unsigned mask_value_plus_x  = (*reinterpret_cast<unsigned*>(mask_base_addr + o_grid[i+1][j][k]));
            unsigned mask_value_plus_y  = (*reinterpret_cast<unsigned*>(mask_base_addr + o_grid[i][j+1][k]));
            unsigned mask_value_plus_z  = (*reinterpret_cast<unsigned*>(mask_base_addr + o_grid[i][j][k+1]));
            // careful!! -- TODO: fix
            const double x_face_input_current_plus = (double)(*reinterpret_cast<T*>(x_face_base_addr + o_grid[i+1][j][k]));
            const double y_face_input_current_plus = (double)(*reinterpret_cast<T*>(y_face_base_addr + o_grid[i][j+1][k]));
            const double z_face_input_current_plus = (double)(*reinterpret_cast<T*>(z_face_base_addr + o_grid[i][j][k+1]));

            if (mask_value_minus_x & SPGrid_Face_Type_X_Valid)
                result -= scale * x_face_input_current;
            if (mask_value_plus_x  & SPGrid_Face_Type_X_Valid)
                result += scale * x_face_input_current_plus;

            if (mask_value_minus_y & SPGrid_Face_Type_Y_Valid)
                result -= scale * y_face_input_current;
            if (mask_value_plus_y  & SPGrid_Face_Type_Y_Valid)
                result += scale * y_face_input_current_plus;

            if (mask_value_minus_z & SPGrid_Face_Type_Z_Valid)
                result -= scale * z_face_input_current;
            if (mask_value_plus_z  & SPGrid_Face_Type_Z_Valid)
                result += scale * z_face_input_current_plus;
#endif

            cell_data_out[cur_index] = result;
              
        }
        cur_index++;
      }
    }
    free(offset_grid_ptr);
}

template class Divergence_Helper<float,5,2>;
template class Divergence_Helper<float,6,2>;

template class Divergence_Helper<float,5,3>;
template class Divergence_Helper<float,6,3>;
//template class Divergence_Helper<float,5,2>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class Divergence_Helper<double,5,2>;
template class Divergence_Helper<double,6,2>;

template class Divergence_Helper<double,5,3>;
template class Divergence_Helper<double,6,3>;
#endif
