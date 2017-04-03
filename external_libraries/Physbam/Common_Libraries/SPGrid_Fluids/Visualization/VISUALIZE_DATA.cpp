//#####################################################################
// Copyright (c) 2012, Eftychios Sifakis, Sean Bauer
// Distributed under the FreeBSD license (see license.txt)
//#####################################################################
#include <SPGrid_Fluids/Visualization/VISUALIZE_DATA.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/FACE_ARRAYS.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/GRID_ARRAYS_POLICY_UNIFORM.h>
#include <PhysBAM_Tools/Read_Write/Grids_Uniform/READ_WRITE_GRID.h>
#include <PhysBAM_Tools/Read_Write/Grids_Uniform_Arrays/READ_WRITE_ARRAYS.h>
#include <PhysBAM_Tools/Read_Write/Grids_Uniform_Arrays/READ_WRITE_FACE_ARRAYS.h>

#include <SPGrid/Core/SPGrid_Utilities.h>
#include <SPGrid/Core/SPGrid_Geometry.h>
#include <SPGrid/Core/SPGrid_Allocator.h>
#include <SPGrid/Core/SPGrid_Set.h>
#include <SPGrid/Tools/SPGrid_Block_Iterator.h>
#include <SPGrid_Fluids/Simulation/FLUIDS_SIMULATION_FLAGS.h>
#include <SPGrid_Fluids/Simulation/FLUIDS_SIMULATION_DATA.h>
#include <SPGrid_Fluids/Grids/GRID_TOPOLOGY_HELPER.h>
using namespace PhysBAM;
//#####################################################################
// Visualize_Flags
//#####################################################################
template<class T_STRUCT,class T,int d> void VISUALIZE_DATA<T_STRUCT,T,d>::
Visualize_Data(T_HIERARCHY& hierarchy,const int frame,const std::string directory)
{
    typedef float RW;
    RW rw=RW(); STREAM_TYPE stream_type(rw);
    
    // create top level directory and one per level inside
    if(frame==0){
        FILE_UTILITIES::Create_Directory(directory);
        for(int level=1;level<=hierarchy.Levels();level++){
            std::string l_str=STRING_UTILITIES::string_sprintf("%d",level);
            FILE_UTILITIES::Create_Directory(directory+"/"+l_str);
            FILE_UTILITIES::Create_Directory(directory+"/"+l_str+"/common");}}
    
    for(int level=1;level<=hierarchy.Levels();level++)
        Write_Output(stream_type,hierarchy,frame,directory,level);
}
template<class T_STRUCT,class T,int d> void VISUALIZE_DATA<T_STRUCT,T,d>::
Write_Output(STREAM_TYPE stream_type,T_HIERARCHY& hierarchy,const int frame,const std::string directory,const int level)
{
    typedef ARRAY<T,T_INDEX> T_ARRAYS_SCALAR;
    typedef ARRAY<T,FACE_INDEX<d> > T_FACE_ARRAYS;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::FACE_ARRAYS::template REBIND<bool>::TYPE T_FACE_ARRAYS_BOOL;
    typedef ARRAY<bool,T_INDEX> T_ARRAYS_BOOL;
    
    VECTOR<T T_STRUCT::*,d> face_velocity_channels;
    face_velocity_channels(1)=&T_STRUCT::ch0;
    face_velocity_channels(2)=&T_STRUCT::ch1;
    if(d==3) face_velocity_channels(3)=&T_STRUCT::ch2;

    VECTOR<void*,d> face_velocity_data_ptrs;
        for(int v=1;v<=d;v++)
            face_velocity_data_ptrs(v)=hierarchy.Array(level,face_velocity_channels(v)).Get_Data_Ptr();

    T_GRID mac_grid=hierarchy.Grid(level).Get_MAC_Grid();    

    Const_flag_array_type flags=hierarchy.Allocator(level).Get_Const_Array(&T_STRUCT::flags);
    Const_data_array_type spgrid_density=hierarchy.Allocator(level).Get_Const_Array(&T_STRUCT::ch3);
    Const_data_array_type spgrid_pressure=hierarchy.Allocator(level).Get_Const_Array(&T_STRUCT::ch9);

    T_ARRAYS_SCALAR density(mac_grid.Domain_Indices(3));
    T_ARRAYS_SCALAR pressure(mac_grid.Domain_Indices(1));
    T_FACE_ARRAYS face_velocities(mac_grid,1);

    for(SPGrid_Block_Iterator<typename Flag_array_type::MASK> iterator(hierarchy.Blocks(level));iterator.Valid();iterator.Next()){

        T_INDEX index=iterator.Index().template Cast<T_INDEX>();
        unsigned cell_flags=iterator.Data(flags);

        // density and pressure
        if(cell_flags&(SPGrid_Cell_Type_Interior|SPGrid_Cell_Type_Ghost)){
            density(index)=iterator.Data(spgrid_density);
            pressure(index)=iterator.Data(spgrid_pressure);}

        // face velocities
        for(int axis=1;axis<=d;axis++)
            if(cell_flags&GRID_TOPOLOGY_HELPER<typename Flag_array_type::MASK>::Face_Valid_Mask(axis))
                face_velocities(axis,index)=iterator.template Data<Data_array_type>(face_velocity_data_ptrs(axis));
    }

    std::string l_str=STRING_UTILITIES::string_sprintf("%d",level);
    std::string frame_str=STRING_UTILITIES::string_sprintf("%d",frame);
    
    if(frame==0)
        FILE_UTILITIES::Write_To_File(stream_type,directory+"/"+l_str+"/common/grid",mac_grid);
    
    std::string current_dir=directory+"/"+l_str+"/"+frame_str;
    FILE_UTILITIES::Create_Directory(current_dir);
    FILE_UTILITIES::Write_To_File(stream_type,current_dir+"/density",density);
    FILE_UTILITIES::Write_To_File(stream_type,current_dir+"/mac_velocities",face_velocities);
    FILE_UTILITIES::Write_To_File(stream_type,current_dir+"/grid",mac_grid);
    FILE_UTILITIES::Write_To_File(stream_type,current_dir+"/pressure",pressure);
}
//#####################################################################
template class VISUALIZE_DATA<FLUIDS_SIMULATION_DATA<float>,float,2>;
template class VISUALIZE_DATA<FLUIDS_SIMULATION_DATA<float>,float,3>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class VISUALIZE_DATA<FLUIDS_SIMULATION_DATA<double>,double,2>;
template class VISUALIZE_DATA<FLUIDS_SIMULATION_DATA<double>,double,3>;
#endif
