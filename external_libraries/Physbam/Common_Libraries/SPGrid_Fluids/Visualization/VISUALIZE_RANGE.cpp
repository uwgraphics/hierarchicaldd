//#####################################################################
// Copyright (c) 2012, Eftychios Sifakis, Sean Bauer
// Distributed under the FreeBSD license (see license.txt)
//#####################################################################
#include <SPGrid_Fluids/Visualization/VISUALIZE_RANGE.h>

#include <SPGrid/Tools/SPGrid_Clear.h>
#include <SPGrid/Tools/SPGrid_Block_Iterator.h>
#include <SPGrid_Fluids/Read_Write/READ_HIERARCHY.h>
#include <PhysBAM_Tools/Read_Write/Grids_Uniform/READ_WRITE_GRID.h>
#include <PhysBAM_Tools/Read_Write/Grids_Uniform_Arrays/READ_WRITE_ARRAYS.h>
#include <PhysBAM_Tools/Read_Write/Grids_Uniform_Arrays/READ_WRITE_FACE_ARRAYS.h>
#include <SPGrid_Fluids/Simulation/FLUIDS_SIMULATION_DATA.h>
#include <SPGrid_Fluids/Simulation/FLUIDS_SIMULATION_FLAGS.h>
#include <SPGrid_Fluids/Grids/GRID_TOPOLOGY_HELPER.h>

using namespace PhysBAM;
//#####################################################################
// Visualize_Range
//#####################################################################
template<class T_STRUCT,class T,int d> void VISUALIZE_RANGE<T_STRUCT,T,d>::
Visualize_Range(const RANGE<TV>& range,const ARRAY<int>& frames,const std::string input_dir,const std::string output_dir)
{
    // Setup
    RW rw=RW(); STREAM_TYPE stream_type(rw);
    T T_STRUCT::* density_channel =&T_STRUCT::ch3;
    VECTOR<T T_STRUCT::*,d> face_velocity_channels;
    face_velocity_channels(1)=&T_STRUCT::ch0;
    face_velocity_channels(2)=&T_STRUCT::ch1;
    if(d==3) face_velocity_channels(3)=&T_STRUCT::ch2;    

    // Read Hierarchy
    T_HIERARCHY hierarchy=READ_HIERARCHY<T_STRUCT,T,d>::Read_Hierarchy(input_dir);
    const int levels=hierarchy.Levels();
    GRID<TV> coarsest_grid=hierarchy.Grid(levels).Get_MAC_Grid();
    GRID<TV> finest_grid=hierarchy.Grid(1).Get_MAC_Grid();
    
    // workaround
    const TV epsilon=(T).5*finest_grid.dX*T();

    // Adjust range
    RANGE<TV> working_range(range);
    ARRAY<RANGE<T_INDEX> > working_cell_domains(levels);
     
    for(int v=1;v<=d;v++){
        working_range.min_corner(v)=max(working_range.min_corner(v),coarsest_grid.Domain().min_corner(v));
        working_range.max_corner(v)=min(working_range.max_corner(v),coarsest_grid.Domain().max_corner(v));}

    for(int level=1;level<=levels;level++)
        working_cell_domains(level)=RANGE<T_INDEX>(hierarchy.Grid(level).Cell(working_range.min_corner+epsilon,0),
                                                   hierarchy.Grid(level).Cell(working_range.max_corner-epsilon+hierarchy.Grid(level).dX,0));

    working_range.min_corner=coarsest_grid.Node(coarsest_grid.Cell(working_range.min_corner+epsilon,0));
    working_range.max_corner=coarsest_grid.Node(coarsest_grid.Cell(working_range.max_corner-epsilon+coarsest_grid.dX,0));
    LOG::cout<<"Working range: "<<working_range<<std::endl;
    
    // iterate over frames
    int output_frame=0;
    for(int i=1;i<=frames.m;i++){        
        const int input_frame=frames(i);
        std::string input_f=STRING_UTILITIES::string_sprintf("%d",input_frame);

        // Clear channels
        for(int level=1;level<=levels;level++){
            SPGrid_Computations::Clear<T_STRUCT,T,d>(hierarchy.Allocator(level),hierarchy.Blocks(level),density_channel);
            for(int axis=1;axis<=d;axis++)
                SPGrid_Computations::Clear<T_STRUCT,T,d>(hierarchy.Allocator(level),hierarchy.Blocks(level),face_velocity_channels(axis));}
        
        // Read channels
        hierarchy.Read_Data_Channel(input_dir+"/"+input_f+"/spgrid_density",density_channel);
        hierarchy.Read_Data_Channel(input_dir+"/"+input_f+"/spgrid_u",face_velocity_channels(1));
        hierarchy.Read_Data_Channel(input_dir+"/"+input_f+"/spgrid_v",face_velocity_channels(2));
        if(d==3) hierarchy.Read_Data_Channel(input_dir+"/"+input_f+"/spgrid_w",face_velocity_channels(3));
        

        // Visualize
        Visualize_Range_Helper(hierarchy,working_range,working_cell_domains,output_frame,output_dir);
        
        output_frame++;
    }
}
//#####################################################################
// Visualize_Range_Helper
//#####################################################################
template<class T_STRUCT,class T,int d> void VISUALIZE_RANGE<T_STRUCT,T,d>::
Visualize_Range_Helper(T_HIERARCHY& hierarchy,const RANGE<TV>& working_range,const ARRAY<RANGE<T_INDEX> >& working_cell_domains,const int output_frame,const std::string output_dir)
{
    // setup
    RW rw=RW(); STREAM_TYPE stream_type(rw);
    const int levels=hierarchy.Levels();    
    T T_STRUCT::* density_channel =&T_STRUCT::ch3;
    VECTOR<T T_STRUCT::*,d> face_velocity_channels;
    face_velocity_channels(1)=&T_STRUCT::ch0;
    face_velocity_channels(2)=&T_STRUCT::ch1;
    if(d==3) face_velocity_channels(3)=&T_STRUCT::ch2;
    
    // create top level directory and one per level inside
    if(output_frame==0){
        FILE_UTILITIES::Create_Directory(output_dir);
        for(int level=1;level<=levels;level++){
            std::string l_str=STRING_UTILITIES::string_sprintf("%d",level);
            FILE_UTILITIES::Create_Directory(output_dir+"/"+l_str);
            FILE_UTILITIES::Create_Directory(output_dir+"/"+l_str+"/common");}}
    
    // iterate over levels
    for(int level=1;level<=levels;level++){
        
        VECTOR<void*,d> face_velocity_data_ptrs;
        for(int v=1;v<=d;v++)
            face_velocity_data_ptrs(v)=hierarchy.Array(level,face_velocity_channels(v)).Get_Data_Ptr();

        T_GRID g(working_cell_domains(level).max_corner-working_cell_domains(level).min_corner,working_range);
        T_GRID mac_grid=g.Get_MAC_Grid();

        Const_flag_array_type flags=hierarchy.Allocator(level).Get_Const_Array(&T_STRUCT::flags);
        Const_data_array_type spgrid_density=hierarchy.Allocator(level).Get_Const_Array(density_channel);

        T_ARRAYS_SCALAR density(mac_grid.Domain_Indices(3));
        T_FACE_ARRAYS face_velocities(mac_grid,1);

        T_INDEX shift=working_cell_domains(level).min_corner;

        for(SPGrid_Block_Iterator<typename Flag_array_type::MASK> iterator(hierarchy.Blocks(level));iterator.Valid();iterator.Next()){

            T_INDEX index=iterator.Index().template Cast<T_INDEX>();
            
            if(working_cell_domains(level).Lazy_Inside(index)){
                
                unsigned cell_flags=iterator.Data(flags);
                
                // density and pressure
                if(cell_flags&(SPGrid_Cell_Type_Interior|SPGrid_Cell_Type_Ghost))
                    density(index-shift)=iterator.Data(spgrid_density);
                
                // face velocities
                for(int axis=1;axis<=d;axis++)
                    if(cell_flags&GRID_TOPOLOGY_HELPER<typename Flag_array_type::MASK>::Face_Valid_Mask(axis))
                        face_velocities(axis,index-shift)=iterator.template Data<Data_array_type>(face_velocity_data_ptrs(axis));
            }
        }

        std::string l_str=STRING_UTILITIES::string_sprintf("%d",level);
        std::string output_frame_str=STRING_UTILITIES::string_sprintf("%d",output_frame);
    
        if(output_frame==0)
            FILE_UTILITIES::Write_To_File(stream_type,output_dir+"/"+l_str+"/common/grid",mac_grid);
    
        std::string current_dir=output_dir+"/"+l_str+"/"+output_frame_str;
        FILE_UTILITIES::Create_Directory(current_dir);
        FILE_UTILITIES::Write_To_File(stream_type,current_dir+"/density",density);
        FILE_UTILITIES::Write_To_File(stream_type,current_dir+"/mac_velocities",face_velocities);
        FILE_UTILITIES::Write_To_File(stream_type,current_dir+"/grid",mac_grid);
    }
}
//#####################################################################
template class VISUALIZE_RANGE<FLUIDS_SIMULATION_DATA<float>,float,2>;
template class VISUALIZE_RANGE<FLUIDS_SIMULATION_DATA<float>,float,3>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class VISUALIZE_RANGE<FLUIDS_SIMULATION_DATA<double>,double,2>;
template class VISUALIZE_RANGE<FLUIDS_SIMULATION_DATA<double>,double,3>;
#endif
