//#####################################################################
// Copyright 2013, Sean Bauer, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <SPGrid_Fluids/Solvers/Forward_Substitution_Helper.h>

#include <SPGrid_Fluids/Projection/Laplace_Helper.h>
#include <SPGrid_Fluids/Simulation/FLUIDS_SIMULATION_FLAGS.h>
#include <SPGrid_Fluids/Simulation/FLUIDS_SIMULATION_DATA.h>

using namespace SPGrid;
using namespace PhysBAM;

template<class T,int log2_struct,int d>
struct Forward_Substitution_Thread_Helper:public PTHREAD_QUEUE::TASK
{
    Forward_Substitution_Helper<T,log2_struct,d>* const obj;
    const int index_start,index_end;
    Forward_Substitution_Thread_Helper(Forward_Substitution_Helper<T,log2_struct,d>* const obj_input,const int index_start_input,const int index_end_input)
        :obj(obj_input),index_start(index_start_input),index_end(index_end_input) {}
    void Run(){obj->Run_Index_Range(index_start,index_end);}
};

//#####################################################################
// Function Run_Parallel
//#####################################################################
template<class T, int log2_struct> void Forward_Substitution_Helper<T,log2_struct,2>::Run_Parallel(const int number_of_partitions)
{
    for(int partition=0;partition<number_of_partitions;partition++){
        // Calculate indicies of current partition
        int first_index_of_partition=(size/number_of_partitions)*(partition)+std::min(size%number_of_partitions,partition);
        int last_index_of_partition=(size/number_of_partitions)*(partition+1)+std::min(size%number_of_partitions,partition+1)-1;
        // Create helper object
        Forward_Substitution_Thread_Helper<T,log2_struct,d>* task=new Forward_Substitution_Thread_Helper<T,log2_struct,d>(this,first_index_of_partition,last_index_of_partition);
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
template <class T, int log2_struct> void Forward_Substitution_Helper<T,log2_struct,2>::Run_Index_Range(const int index_start,const int index_end)
{  
    Static_Assert(sizeof(T)==sizeof(unsigned));

    // Compute shadow grid of linear offsets
    unsigned long* offset_grid_ptr = (unsigned long*)malloc( (og_xsize) * (og_ysize) * sizeof(unsigned long));
    typedef unsigned long (&offset_grid_type)[og_xsize][og_ysize];
    offset_grid_type o_grid = reinterpret_cast<offset_grid_type>(*offset_grid_ptr);
    
    // Iterating through all block indices
    for(int index=index_start;index<=index_end;index++){
        T* output = reinterpret_cast<T*>((unsigned long)x + b[index]);

        unsigned* mask_in  = reinterpret_cast<unsigned*>((unsigned long)mask + b[index]);
        unsigned long packed_offset = (unsigned long)output - (unsigned long)x;

        unsigned long mask_base_addr = reinterpret_cast<unsigned long>(mask);
        unsigned long x_base_addr = reinterpret_cast<unsigned long>(x);
        unsigned long Lx_base_addr = reinterpret_cast<unsigned long>(Lx);
        unsigned long Ly_base_addr = reinterpret_cast<unsigned long>(Ly);

        Laplace_Helper<T,log2_struct,2>::ComputeShadowGrid(offset_grid_ptr, mask_in, packed_offset);

        int cur_index = 0;
        // Actually process elements
        for(int i=xmin;i<=xmax;i++)
        for(int j=ymin;j<=ymax;j++){

            const unsigned mask_value=mask_in[cur_index];

            if(mask_value & SPGrid_Cell_Type_Active){

                const T x_value=output[cur_index];

                if(mask_value & SPGrid_Face_Minus_X_Active){
                    bool use_neighbor=o_grid[i-1][j]>o_grid[i][j];
                    if(!use_neighbor){
                        const unsigned mask_other_value=(*reinterpret_cast<unsigned*>(mask_base_addr + o_grid[i-1][j]));
                        if(mask_other_value & SPGrid_Cell_Type_Ghost) use_neighbor=true;}
                    if(use_neighbor){
                        const T coefficient=*reinterpret_cast<T*>(Lx_base_addr + o_grid[i][j]);
                        if(coefficient){
                            T& x_other_value=(*reinterpret_cast<T*>(x_base_addr + o_grid[i-1][j]));
                            x_other_value -= coefficient*x_value;}}}

                if(mask_value & SPGrid_Face_Plus_X_Active){
                    bool use_neighbor=o_grid[i+1][j]>o_grid[i][j];
                    if(!use_neighbor){
                        const unsigned mask_other_value=(*reinterpret_cast<unsigned*>(mask_base_addr + o_grid[i+1][j]));
                        if(mask_other_value & SPGrid_Cell_Type_Ghost) use_neighbor=true;}
                    if(use_neighbor){
                        const T coefficient=*reinterpret_cast<T*>(Lx_base_addr + o_grid[i+1][j]);
                        if(coefficient){
                            T& x_other_value=(*reinterpret_cast<T*>(x_base_addr + o_grid[i+1][j]));
                            x_other_value -= coefficient*x_value;}}}

                if(mask_value & SPGrid_Face_Minus_Y_Active){
                    bool use_neighbor=o_grid[i][j-1]>o_grid[i][j];
                    if(!use_neighbor){
                        const unsigned mask_other_value=(*reinterpret_cast<unsigned*>(mask_base_addr + o_grid[i][j-1]));
                        if(mask_other_value & SPGrid_Cell_Type_Ghost) use_neighbor=true;}
                    if(use_neighbor){
                        const T coefficient=*reinterpret_cast<T*>(Ly_base_addr + o_grid[i][j]);
                        if(coefficient){
                            T& x_other_value=(*reinterpret_cast<T*>(x_base_addr + o_grid[i][j-1]));
                            x_other_value -= coefficient*x_value;}}}

                if(mask_value & SPGrid_Face_Plus_Y_Active){
                    bool use_neighbor=o_grid[i][j+1]>o_grid[i][j];
                    if(!use_neighbor){
                        const unsigned mask_other_value=(*reinterpret_cast<unsigned*>(mask_base_addr + o_grid[i][j+1]));
                        if(mask_other_value & SPGrid_Cell_Type_Ghost) use_neighbor=true;}
                    if(use_neighbor){
                        const T coefficient=*reinterpret_cast<T*>(Ly_base_addr + o_grid[i][j+1]);
                        if(coefficient){
                            T& x_other_value=(*reinterpret_cast<T*>(x_base_addr + o_grid[i][j+1]));
                            x_other_value -= coefficient*x_value;}}}

            }
            cur_index++;
        }
    }

    free(offset_grid_ptr);
}
//#####################################################################
// Function Run_Parallel
//#####################################################################
template<class T, int log2_struct> void Forward_Substitution_Helper<T,log2_struct,3>::Run_Parallel(const int number_of_partitions)
{
    for(int partition=0;partition<number_of_partitions;partition++){
        // Calculate indicies of current partition
        int first_index_of_partition=(size/number_of_partitions)*(partition)+std::min(size%number_of_partitions,partition);
        int last_index_of_partition=(size/number_of_partitions)*(partition+1)+std::min(size%number_of_partitions,partition+1)-1;
        // Create helper object
        Forward_Substitution_Thread_Helper<T,log2_struct,d>* task=new Forward_Substitution_Thread_Helper<T,log2_struct,d>(this,first_index_of_partition,last_index_of_partition);
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
template <class T, int log2_struct> void Forward_Substitution_Helper<T,log2_struct,3>::Run_Index_Range(const int index_start,const int index_end)
{  
    Static_Assert(sizeof(T)==sizeof(unsigned));

    // Compute shadow grid of linear offsets
    unsigned long* offset_grid_ptr = (unsigned long*)malloc( (og_xsize) * (og_ysize) * (og_zsize) * sizeof(unsigned long));
    typedef unsigned long (&offset_grid_type)[og_xsize][og_ysize][og_zsize];
    offset_grid_type o_grid = reinterpret_cast<offset_grid_type>(*offset_grid_ptr);
    
    // Iterating through all block indices
    for(int index=index_start;index<=index_end;index++){
        T* output = reinterpret_cast<T*>((unsigned long)x + b[index]);

        unsigned* mask_in  = reinterpret_cast<unsigned*>((unsigned long)mask + b[index]);
        unsigned long packed_offset = (unsigned long)output - (unsigned long)x;

        unsigned long mask_base_addr = reinterpret_cast<unsigned long>(mask);
        unsigned long x_base_addr = reinterpret_cast<unsigned long>(x);
        unsigned long Lx_base_addr = reinterpret_cast<unsigned long>(Lx);
        unsigned long Ly_base_addr = reinterpret_cast<unsigned long>(Ly);
        unsigned long Lz_base_addr = reinterpret_cast<unsigned long>(Lz);

        Laplace_Helper<T,log2_struct,3>::ComputeShadowGrid(offset_grid_ptr, mask_in, packed_offset);

        int cur_index = 0;
        // Actually process elements
        for(int i=xmin;i<=xmax;i++)
        for(int j=ymin;j<=ymax;j++)
        for(int k=zmin;k<=zmax;k++){

            const unsigned mask_value=mask_in[cur_index];

            if(mask_value & SPGrid_Cell_Type_Active){

                const T x_value=output[cur_index];

                if(mask_value & SPGrid_Face_Minus_X_Active){
                    bool use_neighbor=o_grid[i-1][j][k]>o_grid[i][j][k];
                    if(!use_neighbor){
                        const unsigned mask_other_value=(*reinterpret_cast<unsigned*>(mask_base_addr + o_grid[i-1][j][k]));
                        if(mask_other_value & SPGrid_Cell_Type_Ghost) use_neighbor=true;}
                    if(use_neighbor){
                        const T coefficient=*reinterpret_cast<T*>(Lx_base_addr + o_grid[i][j][k]);
                        if(coefficient){
                            T& x_other_value=(*reinterpret_cast<T*>(x_base_addr + o_grid[i-1][j][k]));
                            x_other_value -= coefficient*x_value;}}}

                if(mask_value & SPGrid_Face_Plus_X_Active){
                    bool use_neighbor=o_grid[i+1][j][k]>o_grid[i][j][k];
                    if(!use_neighbor){
                        const unsigned mask_other_value=(*reinterpret_cast<unsigned*>(mask_base_addr + o_grid[i+1][j][k]));
                        if(mask_other_value & SPGrid_Cell_Type_Ghost) use_neighbor=true;}
                    if(use_neighbor){
                        const T coefficient=*reinterpret_cast<T*>(Lx_base_addr + o_grid[i+1][j][k]);
                        if(coefficient){
                            T& x_other_value=(*reinterpret_cast<T*>(x_base_addr + o_grid[i+1][j][k]));
                            x_other_value -= coefficient*x_value;}}}

                if(mask_value & SPGrid_Face_Minus_Y_Active){
                    bool use_neighbor=o_grid[i][j-1][k]>o_grid[i][j][k];
                    if(!use_neighbor){
                        const unsigned mask_other_value=(*reinterpret_cast<unsigned*>(mask_base_addr + o_grid[i][j-1][k]));
                        if(mask_other_value & SPGrid_Cell_Type_Ghost) use_neighbor=true;}
                    if(use_neighbor){
                        const T coefficient=*reinterpret_cast<T*>(Ly_base_addr + o_grid[i][j][k]);
                        if(coefficient){
                            T& x_other_value=(*reinterpret_cast<T*>(x_base_addr + o_grid[i][j-1][k]));
                            x_other_value -= coefficient*x_value;}}}

                if(mask_value & SPGrid_Face_Plus_Y_Active){
                    bool use_neighbor=o_grid[i][j+1][k]>o_grid[i][j][k];
                    if(!use_neighbor){
                        const unsigned mask_other_value=(*reinterpret_cast<unsigned*>(mask_base_addr + o_grid[i][j+1][k]));
                        if(mask_other_value & SPGrid_Cell_Type_Ghost) use_neighbor=true;}
                    if(use_neighbor){
                        const T coefficient=*reinterpret_cast<T*>(Ly_base_addr + o_grid[i][j+1][k]);
                        if(coefficient){
                            T& x_other_value=(*reinterpret_cast<T*>(x_base_addr + o_grid[i][j+1][k]));
                            x_other_value -= coefficient*x_value;}}}

                if(mask_value & SPGrid_Face_Minus_Z_Active){
                    bool use_neighbor=o_grid[i][j][k-1]>o_grid[i][j][k];
                    if(!use_neighbor){
                        const unsigned mask_other_value=(*reinterpret_cast<unsigned*>(mask_base_addr + o_grid[i][j][k-1]));
                        if(mask_other_value & SPGrid_Cell_Type_Ghost) use_neighbor=true;}
                    if(use_neighbor){
                        const T coefficient=*reinterpret_cast<T*>(Lz_base_addr + o_grid[i][j][k]);
                        if(coefficient){
                            T& x_other_value=(*reinterpret_cast<T*>(x_base_addr + o_grid[i][j][k-1]));
                            x_other_value -= coefficient*x_value;}}}

                if(mask_value & SPGrid_Face_Plus_Z_Active){
                    bool use_neighbor=o_grid[i][j][k+1]>o_grid[i][j][k];
                    if(!use_neighbor){
                        const unsigned mask_other_value=(*reinterpret_cast<unsigned*>(mask_base_addr + o_grid[i][j][k+1]));
                        if(mask_other_value & SPGrid_Cell_Type_Ghost) use_neighbor=true;}
                    if(use_neighbor){
                        const T coefficient=*reinterpret_cast<T*>(Lz_base_addr + o_grid[i][j][k+1]);
                        if(coefficient){
                            T& x_other_value=(*reinterpret_cast<T*>(x_base_addr + o_grid[i][j][k+1]));
                            x_other_value -= coefficient*x_value;}}}

            }
            cur_index++;
        }
    }

    free(offset_grid_ptr);
}
//#####################################################################
template class Forward_Substitution_Helper<float,5,2>;
template class Forward_Substitution_Helper<float,6,2>;
template class Forward_Substitution_Helper<float,5,3>;
template class Forward_Substitution_Helper<float,6,3>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class Forward_Substitution_Helper<double,5,2>;
template class Forward_Substitution_Helper<double,6,2>;
template class Forward_Substitution_Helper<double,5,3>;
template class Forward_Substitution_Helper<double,6,3>;
#endif
