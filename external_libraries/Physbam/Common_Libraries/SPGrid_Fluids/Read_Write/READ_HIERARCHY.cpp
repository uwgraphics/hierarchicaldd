//#####################################################################
// Copyright (c) 2014, Raj Setaluri, Mridul Aanjneya
// Distributed under the FreeBSD license (see license.txt)
//#####################################################################
#include <SPGrid_Fluids/Read_Write/READ_HIERARCHY.h>
#include <PhysBAM_Tools/Read_Write/Grids_Uniform/READ_WRITE_GRID.h>
#include <SPGrid_Fluids/Simulation/FLUIDS_SIMULATION_DATA.h>

using namespace PhysBAM;
//#####################################################################
// Read_Hierarchy
//#####################################################################
template<class T_STRUCT,class T,int d> GRID_HIERARCHY<T_STRUCT,T,d>& READ_HIERARCHY<T_STRUCT,T,d>::
Read_Hierarchy(const std::string& base_dir,const int frame)
{
    RW rw=RW(); STREAM_TYPE stream_type(rw);
    T_HIERARCHY* hierarchy;T_GRID fine_mac_grid;
    std::string f=STRING_UTILITIES::string_sprintf("%d",frame);
    int levels=0;std::string directory_filename=base_dir+"/"+f;
    FILE_UTILITIES::template Read_From_Text_File<int>(STRING_UTILITIES::string_sprintf("%s/levels",directory_filename.c_str()),levels);
    FILE_UTILITIES::Read_From_File(stream_type,base_dir+"/common/fine_grid",fine_mac_grid);
    hierarchy=new GRID_HIERARCHY<T_STRUCT,T,d>(fine_mac_grid,levels);
    hierarchy->Initialize_Sets();
    hierarchy->Read_Block_Offsets(base_dir+"/"+f+"/block_offsets");
    hierarchy->Read_Flags_Channel(base_dir+"/"+f+"/flags");
    return *hierarchy;
}
//#####################################################################
// Read_Hierarchy_Pointer
//#####################################################################
template<class T_STRUCT,class T,int d> GRID_HIERARCHY<T_STRUCT,T,d>* READ_HIERARCHY<T_STRUCT,T,d>::
Read_Hierarchy_Pointer(const std::string& base_dir,const int frame)
{
    RW rw=RW(); STREAM_TYPE stream_type(rw);
    T_HIERARCHY* hierarchy;T_GRID fine_mac_grid;
    std::string f=STRING_UTILITIES::string_sprintf("%d",frame);
    int levels=0;std::string directory_filename=base_dir+"/"+f;
    FILE_UTILITIES::template Read_From_Text_File<int>(STRING_UTILITIES::string_sprintf("%s/levels",directory_filename.c_str()),levels);
    FILE_UTILITIES::Read_From_File(stream_type,base_dir+"/common/fine_grid",fine_mac_grid);
    hierarchy=new GRID_HIERARCHY<T_STRUCT,T,d>(fine_mac_grid,levels);
    hierarchy->Initialize_Sets();
    hierarchy->Read_Block_Offsets(base_dir+"/"+f+"/block_offsets");
    hierarchy->Read_Flags_Channel(base_dir+"/"+f+"/flags");
    return hierarchy;
}
//#####################################################################
template class READ_HIERARCHY<FLUIDS_SIMULATION_DATA<float>,float,2>;
template class READ_HIERARCHY<FLUIDS_SIMULATION_DATA<float>,float,3>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class READ_HIERARCHY<FLUIDS_SIMULATION_DATA<double>,double,2>;
template class READ_HIERARCHY<FLUIDS_SIMULATION_DATA<double>,double,3>;
#endif
