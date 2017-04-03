//#####################################################################
// Copyright 2013, Sean Bauer, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include "Threaded_Restriction_Stencil_Helper.h"
#include <SPGrid_Fluids/Simulation/FLUIDS_SIMULATION_FLAGS.h>
#include <SPGrid_Fluids/Simulation/FLUIDS_SIMULATION_DATA.h>

using namespace SPGrid;
using namespace PhysBAM;
//#####################################################################
// Threaded_Restriction_Stencil_Thread_Helper
//#####################################################################
template<class T,int log2_struct,int d>
struct Threaded_Restriction_Stencil_Thread_Helper:public PTHREAD_QUEUE::TASK
{
    Threaded_Restriction_Stencil_Helper<T,log2_struct,d>* const obj;
    const int index_start,index_end;
    Threaded_Restriction_Stencil_Thread_Helper(Threaded_Restriction_Stencil_Helper<T,log2_struct,d>* const obj_input,const int index_start_input,const int index_end_input)
        :obj(obj_input),index_start(index_start_input),index_end(index_end_input) {}
    void Run(){obj->Run_Index_Range(index_start,index_end);}
};
//#####################################################################
// Function Run_Parallel - 2D
//#####################################################################
template<class T,int log2_struct> void Threaded_Restriction_Stencil_Helper<T,log2_struct,2>::
Run_Parallel(const int number_of_partitions)
{
    for(int partition=0;partition<number_of_partitions;partition++){
        // Calculate indicies of current partition
        int first_index_of_partition=(size/number_of_partitions)*(partition)+std::min(size%number_of_partitions,partition);
        int last_index_of_partition=(size/number_of_partitions)*(partition+1)+std::min(size%number_of_partitions,partition+1)-1;
        // Create helper object
        Threaded_Restriction_Stencil_Thread_Helper<T,log2_struct,d>* task=new Threaded_Restriction_Stencil_Thread_Helper<T,log2_struct,d>(this,first_index_of_partition,last_index_of_partition);
        // Enqueue
        pthread_queue->Queue(task);}
    // Wait for all tasks to complete
    pthread_queue->Wait();
}
//#####################################################################
// Function ComputeFineShadowGrid  - 2D
//#####################################################################
template<class T,int log2_struct> void Threaded_Restriction_Stencil_Helper<T,log2_struct,2>::
ComputeFineShadowGrid(unsigned long* offset_grid_ptr,const unsigned long fine_packed_offset)
{
    typedef unsigned long (&offset_grid_type)[fine_og_xsize][fine_og_ysize];
    offset_grid_type o_grid=reinterpret_cast<offset_grid_type>(*offset_grid_ptr);
    unsigned long current_offset=fine_packed_offset;
    for(int i=0;i<fine_og_xsize;i++){o_grid[i][0]=current_offset;
        current_offset=T_MASK::template Packed_OffsetXdim<1>(current_offset);}
    for(int i=0;i<fine_og_xsize;i++){current_offset=o_grid[i][0];
    for(int j=1;j<fine_og_ysize;j++){current_offset=T_MASK::template Packed_OffsetYdim<1>(current_offset);
        o_grid[i][j]=current_offset;}}
}
//#####################################################################
// Function Run_Index_Range  - 2D
//#####################################################################
template<class T,int log2_struct> void Threaded_Restriction_Stencil_Helper<T,log2_struct,2>::
Run_Index_Range(const int index_start,const int index_end)
{
    // Shadow grid declaration
    unsigned long* offset_grid_ptr=(unsigned long*)malloc((fine_og_xsize)*(fine_og_ysize)*sizeof(unsigned long));
    typedef unsigned long (&offset_grid_type)[fine_og_xsize][fine_og_ysize];
    offset_grid_type o_grid=reinterpret_cast<offset_grid_type>(*offset_grid_ptr);
    const unsigned long fine_data_base_addr=reinterpret_cast<unsigned long>(fine_data);
    // Iterating through all block indices
    for(int index=index_start;index<=index_end;index++){
        T* const coarse_data_ptr=reinterpret_cast<T*>((unsigned long)coarse_data+b[index]);
        const T* const fine_data_ptr=reinterpret_cast<T*>((unsigned long)fine_data+b[index]);
        const unsigned* const coarse_mask_ptr=reinterpret_cast<unsigned*>((unsigned long)coarse_mask+b[index]);
        const unsigned long packed_offset=(unsigned long)coarse_data_ptr-(unsigned long)coarse_data;
        unsigned long fine_packed_offset=T_MASK::UpsampleOffset(packed_offset);
        fine_packed_offset=T_MASK::Packed_Offset<-1,-1>(fine_packed_offset);
        ComputeFineShadowGrid(offset_grid_ptr,fine_packed_offset);
        // Iterate within block
        int cur_index=0;
        for(int i=xmin;i<=xmax;i++)
        for(int j=ymin;j<=ymax;j++){
            const unsigned& coarse_mask_value=coarse_mask_ptr[cur_index];
            if(coarse_mask_value&active_flag_mask){
                T value=(T)0.;
                const int ii(2*(i-1)),jj(2*(j-1));
                
                value+=(T)1*(*reinterpret_cast<T*>(fine_data_base_addr+o_grid[ii+0][jj+0]));
                value+=(T)3*(*reinterpret_cast<T*>(fine_data_base_addr+o_grid[ii+0][jj+1]));
                value+=(T)3*(*reinterpret_cast<T*>(fine_data_base_addr+o_grid[ii+0][jj+2]));
                value+=(T)1*(*reinterpret_cast<T*>(fine_data_base_addr+o_grid[ii+0][jj+3]));
                
                value+=(T)3*(*reinterpret_cast<T*>(fine_data_base_addr+o_grid[ii+1][jj+0]));
                value+=(T)9*(*reinterpret_cast<T*>(fine_data_base_addr+o_grid[ii+1][jj+1]));
                value+=(T)9*(*reinterpret_cast<T*>(fine_data_base_addr+o_grid[ii+1][jj+2]));
                value+=(T)3*(*reinterpret_cast<T*>(fine_data_base_addr+o_grid[ii+1][jj+3]));
                
                value+=(T)3*(*reinterpret_cast<T*>(fine_data_base_addr+o_grid[ii+2][jj+0]));
                value+=(T)9*(*reinterpret_cast<T*>(fine_data_base_addr+o_grid[ii+2][jj+1]));
                value+=(T)9*(*reinterpret_cast<T*>(fine_data_base_addr+o_grid[ii+2][jj+2]));
                value+=(T)3*(*reinterpret_cast<T*>(fine_data_base_addr+o_grid[ii+2][jj+3]));
                
                value+=(T)1*(*reinterpret_cast<T*>(fine_data_base_addr+o_grid[ii+3][jj+0]));
                value+=(T)3*(*reinterpret_cast<T*>(fine_data_base_addr+o_grid[ii+3][jj+1]));
                value+=(T)3*(*reinterpret_cast<T*>(fine_data_base_addr+o_grid[ii+3][jj+2]));
                value+=(T)1*(*reinterpret_cast<T*>(fine_data_base_addr+o_grid[ii+3][jj+3]));
                
                coarse_data_ptr[cur_index]+=(value*scale);}
            cur_index++;}}
    free(offset_grid_ptr);
}
//#####################################################################
// Function Run_Parallel - 3D
//#####################################################################
template<class T,int log2_struct> void Threaded_Restriction_Stencil_Helper<T,log2_struct,3>::Run_Parallel(const int number_of_partitions)
{
    for(int partition=0;partition<number_of_partitions;partition++){
        // Calculate indicies of current partition
        int first_index_of_partition=(size/number_of_partitions)*(partition)+std::min(size%number_of_partitions,partition);
        int last_index_of_partition=(size/number_of_partitions)*(partition+1)+std::min(size%number_of_partitions,partition+1)-1;
        // Create helper object
        Threaded_Restriction_Stencil_Thread_Helper<T,log2_struct,d>* task=new Threaded_Restriction_Stencil_Thread_Helper<T,log2_struct,d>(this,first_index_of_partition,last_index_of_partition);
        // Enqueue
        pthread_queue->Queue(task);}
    // Wait for all tasks to complete
    pthread_queue->Wait();
}
//#####################################################################
// Function ComputeShadowGrid - 3D
//#####################################################################
template<class T,int log2_struct> void Threaded_Restriction_Stencil_Helper<T,log2_struct,3>::
ComputeFineShadowGrid(unsigned long* offset_grid_ptr,const unsigned long fine_packed_offset)
{
    typedef unsigned long (&offset_grid_type)[fine_og_xsize][fine_og_ysize][fine_og_zsize];
    offset_grid_type o_grid=reinterpret_cast<offset_grid_type>(*offset_grid_ptr);
    unsigned long current_offset=fine_packed_offset;
    for(int i=0;i<fine_og_xsize;i++){o_grid[i][0][0]=current_offset;
        current_offset=T_MASK::template Packed_OffsetXdim<1>(current_offset);}
    for(int i=0;i<fine_og_xsize;i++){current_offset=o_grid[i][0][0];
    for(int j=1;j<fine_og_ysize;j++){current_offset=T_MASK::template Packed_OffsetYdim<1>(current_offset);
        o_grid[i][j][0]=current_offset;}}
    for(int i=0;i<fine_og_xsize;i++)
    for(int j=0;j<fine_og_ysize;j++){current_offset=o_grid[i][j][0];
    for(int k=1;k<fine_og_zsize;k++){current_offset=T_MASK::template Packed_OffsetZdim<1>(current_offset);
        o_grid[i][j][k]=current_offset;}}
}
//#####################################################################
// Function Run_Index_Range - 3D
//#####################################################################
template<class T,int log2_struct> void Threaded_Restriction_Stencil_Helper<T,log2_struct,3>::Run_Index_Range(const int index_start,const int index_end)
{
    // Shadow grid declaration
    unsigned long* offset_grid_ptr=(unsigned long*)malloc((fine_og_xsize)*(fine_og_ysize)*(fine_og_zsize)*sizeof(unsigned long));
    typedef unsigned long (&offset_grid_type)[fine_og_xsize][fine_og_ysize][fine_og_zsize];
    offset_grid_type o_grid=reinterpret_cast<offset_grid_type>(*offset_grid_ptr);
    const unsigned long fine_data_base_addr=reinterpret_cast<unsigned long>(fine_data);
    // Iterating through all block indices
    for(int index=index_start;index<=index_end;index++){
        T* const coarse_data_ptr=reinterpret_cast<T*>((unsigned long)coarse_data+b[index]);
        const T* const fine_data_ptr=reinterpret_cast<T*>((unsigned long)fine_data+b[index]);
        const unsigned* const coarse_mask_ptr=reinterpret_cast<unsigned*>((unsigned long)coarse_mask+b[index]);
        const unsigned long packed_offset=(unsigned long)coarse_data_ptr-(unsigned long)coarse_data;
        unsigned long fine_packed_offset=T_MASK::UpsampleOffset(packed_offset);
        fine_packed_offset=T_MASK::Packed_Offset<-1,-1,-1>(fine_packed_offset);
        ComputeFineShadowGrid(offset_grid_ptr,fine_packed_offset);
        // Iterate within block
        int cur_index=0;
        for(int i=xmin;i<=xmax;i++)
        for(int j=ymin;j<=ymax;j++)
        for(int k=zmin;k<=zmax;k++){
            const unsigned& coarse_mask_value=coarse_mask_ptr[cur_index];
            if(coarse_mask_value&active_flag_mask){
                T value=(T)0.;
                const int ii(2*(i-1)),jj(2*(j-1)),kk(2*(k-1));
                
                value+=(T) 1.*(*reinterpret_cast<T*>(fine_data_base_addr+o_grid[ii+0][jj+0][kk+0]));
                value+=(T) 3.*(*reinterpret_cast<T*>(fine_data_base_addr+o_grid[ii+0][jj+0][kk+1]));
                value+=(T) 3.*(*reinterpret_cast<T*>(fine_data_base_addr+o_grid[ii+0][jj+0][kk+2]));
                value+=(T) 1.*(*reinterpret_cast<T*>(fine_data_base_addr+o_grid[ii+0][jj+0][kk+3]));
                value+=(T) 3.*(*reinterpret_cast<T*>(fine_data_base_addr+o_grid[ii+0][jj+1][kk+0]));
                value+=(T) 9.*(*reinterpret_cast<T*>(fine_data_base_addr+o_grid[ii+0][jj+1][kk+1]));
                value+=(T) 9.*(*reinterpret_cast<T*>(fine_data_base_addr+o_grid[ii+0][jj+1][kk+2]));
                value+=(T) 3.*(*reinterpret_cast<T*>(fine_data_base_addr+o_grid[ii+0][jj+1][kk+3]));
                value+=(T) 3.*(*reinterpret_cast<T*>(fine_data_base_addr+o_grid[ii+0][jj+2][kk+0]));
                value+=(T) 9.*(*reinterpret_cast<T*>(fine_data_base_addr+o_grid[ii+0][jj+2][kk+1]));
                value+=(T) 9.*(*reinterpret_cast<T*>(fine_data_base_addr+o_grid[ii+0][jj+2][kk+2]));
                value+=(T) 3.*(*reinterpret_cast<T*>(fine_data_base_addr+o_grid[ii+0][jj+2][kk+3]));
                value+=(T) 1.*(*reinterpret_cast<T*>(fine_data_base_addr+o_grid[ii+0][jj+3][kk+0]));
                value+=(T) 3.*(*reinterpret_cast<T*>(fine_data_base_addr+o_grid[ii+0][jj+3][kk+1]));
                value+=(T) 3.*(*reinterpret_cast<T*>(fine_data_base_addr+o_grid[ii+0][jj+3][kk+2]));
                value+=(T) 1.*(*reinterpret_cast<T*>(fine_data_base_addr+o_grid[ii+0][jj+3][kk+3]));
                
                value+=(T) 3.*(*reinterpret_cast<T*>(fine_data_base_addr+o_grid[ii+1][jj+0][kk+0]));
                value+=(T) 9.*(*reinterpret_cast<T*>(fine_data_base_addr+o_grid[ii+1][jj+0][kk+1]));
                value+=(T) 9.*(*reinterpret_cast<T*>(fine_data_base_addr+o_grid[ii+1][jj+0][kk+2]));
                value+=(T) 3.*(*reinterpret_cast<T*>(fine_data_base_addr+o_grid[ii+1][jj+0][kk+3]));
                value+=(T) 9.*(*reinterpret_cast<T*>(fine_data_base_addr+o_grid[ii+1][jj+1][kk+0]));
                value+=(T)27.*(*reinterpret_cast<T*>(fine_data_base_addr+o_grid[ii+1][jj+1][kk+1]));
                value+=(T)27.*(*reinterpret_cast<T*>(fine_data_base_addr+o_grid[ii+1][jj+1][kk+2]));
                value+=(T) 9.*(*reinterpret_cast<T*>(fine_data_base_addr+o_grid[ii+1][jj+1][kk+3]));
                value+=(T) 9.*(*reinterpret_cast<T*>(fine_data_base_addr+o_grid[ii+1][jj+2][kk+0]));
                value+=(T)27.*(*reinterpret_cast<T*>(fine_data_base_addr+o_grid[ii+1][jj+2][kk+1]));
                value+=(T)27.*(*reinterpret_cast<T*>(fine_data_base_addr+o_grid[ii+1][jj+2][kk+2]));
                value+=(T) 9.*(*reinterpret_cast<T*>(fine_data_base_addr+o_grid[ii+1][jj+2][kk+3]));
                value+=(T) 3.*(*reinterpret_cast<T*>(fine_data_base_addr+o_grid[ii+1][jj+3][kk+0]));
                value+=(T) 9.*(*reinterpret_cast<T*>(fine_data_base_addr+o_grid[ii+1][jj+3][kk+1]));
                value+=(T) 9.*(*reinterpret_cast<T*>(fine_data_base_addr+o_grid[ii+1][jj+3][kk+2]));
                value+=(T) 3.*(*reinterpret_cast<T*>(fine_data_base_addr+o_grid[ii+1][jj+3][kk+3]));
                
                value+=(T) 3.*(*reinterpret_cast<T*>(fine_data_base_addr+o_grid[ii+2][jj+0][kk+0]));
                value+=(T) 9.*(*reinterpret_cast<T*>(fine_data_base_addr+o_grid[ii+2][jj+0][kk+1]));
                value+=(T) 9.*(*reinterpret_cast<T*>(fine_data_base_addr+o_grid[ii+2][jj+0][kk+2]));
                value+=(T) 3.*(*reinterpret_cast<T*>(fine_data_base_addr+o_grid[ii+2][jj+0][kk+3]));
                value+=(T) 9.*(*reinterpret_cast<T*>(fine_data_base_addr+o_grid[ii+2][jj+1][kk+0]));
                value+=(T)27.*(*reinterpret_cast<T*>(fine_data_base_addr+o_grid[ii+2][jj+1][kk+1]));
                value+=(T)27.*(*reinterpret_cast<T*>(fine_data_base_addr+o_grid[ii+2][jj+1][kk+2]));
                value+=(T) 9.*(*reinterpret_cast<T*>(fine_data_base_addr+o_grid[ii+2][jj+1][kk+3]));
                value+=(T) 9.*(*reinterpret_cast<T*>(fine_data_base_addr+o_grid[ii+2][jj+2][kk+0]));
                value+=(T)27.*(*reinterpret_cast<T*>(fine_data_base_addr+o_grid[ii+2][jj+2][kk+1]));
                value+=(T)27.*(*reinterpret_cast<T*>(fine_data_base_addr+o_grid[ii+2][jj+2][kk+2]));
                value+=(T) 9.*(*reinterpret_cast<T*>(fine_data_base_addr+o_grid[ii+2][jj+2][kk+3]));
                value+=(T) 3.*(*reinterpret_cast<T*>(fine_data_base_addr+o_grid[ii+2][jj+3][kk+0]));
                value+=(T) 9.*(*reinterpret_cast<T*>(fine_data_base_addr+o_grid[ii+2][jj+3][kk+1]));
                value+=(T) 9.*(*reinterpret_cast<T*>(fine_data_base_addr+o_grid[ii+2][jj+3][kk+2]));
                value+=(T) 3.*(*reinterpret_cast<T*>(fine_data_base_addr+o_grid[ii+2][jj+3][kk+3]));
                
                value+=(T) 1.*(*reinterpret_cast<T*>(fine_data_base_addr+o_grid[ii+3][jj+0][kk+0]));
                value+=(T) 3.*(*reinterpret_cast<T*>(fine_data_base_addr+o_grid[ii+3][jj+0][kk+1]));
                value+=(T) 3.*(*reinterpret_cast<T*>(fine_data_base_addr+o_grid[ii+3][jj+0][kk+2]));
                value+=(T) 1.*(*reinterpret_cast<T*>(fine_data_base_addr+o_grid[ii+3][jj+0][kk+3]));
                value+=(T) 3.*(*reinterpret_cast<T*>(fine_data_base_addr+o_grid[ii+3][jj+1][kk+0]));
                value+=(T) 9.*(*reinterpret_cast<T*>(fine_data_base_addr+o_grid[ii+3][jj+1][kk+1]));
                value+=(T) 9.*(*reinterpret_cast<T*>(fine_data_base_addr+o_grid[ii+3][jj+1][kk+2]));
                value+=(T) 3.*(*reinterpret_cast<T*>(fine_data_base_addr+o_grid[ii+3][jj+1][kk+3]));
                value+=(T) 3.*(*reinterpret_cast<T*>(fine_data_base_addr+o_grid[ii+3][jj+2][kk+0]));
                value+=(T) 9.*(*reinterpret_cast<T*>(fine_data_base_addr+o_grid[ii+3][jj+2][kk+1]));
                value+=(T) 9.*(*reinterpret_cast<T*>(fine_data_base_addr+o_grid[ii+3][jj+2][kk+2]));
                value+=(T) 3.*(*reinterpret_cast<T*>(fine_data_base_addr+o_grid[ii+3][jj+2][kk+3]));
                value+=(T) 1.*(*reinterpret_cast<T*>(fine_data_base_addr+o_grid[ii+3][jj+3][kk+0]));
                value+=(T) 3.*(*reinterpret_cast<T*>(fine_data_base_addr+o_grid[ii+3][jj+3][kk+1]));
                value+=(T) 3.*(*reinterpret_cast<T*>(fine_data_base_addr+o_grid[ii+3][jj+3][kk+2]));
                value+=(T) 1.*(*reinterpret_cast<T*>(fine_data_base_addr+o_grid[ii+3][jj+3][kk+3]));
                
                coarse_data_ptr[cur_index]+=(value*scale);}
            cur_index++;}}
    free(offset_grid_ptr);
}
//#####################################################################
template class Threaded_Restriction_Stencil_Helper<float,5,2>;
template class Threaded_Restriction_Stencil_Helper<float,6,2>;
template class Threaded_Restriction_Stencil_Helper<float,5,3>;
template class Threaded_Restriction_Stencil_Helper<float,6,3>;
//template class Threaded_Restriction_Stencil_Helper<float,5,2>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class Threaded_Restriction_Stencil_Helper<double,5,2>;
template class Threaded_Restriction_Stencil_Helper<double,6,2>;
template class Threaded_Restriction_Stencil_Helper<double,5,3>;
template class Threaded_Restriction_Stencil_Helper<double,6,3>;
#endif
