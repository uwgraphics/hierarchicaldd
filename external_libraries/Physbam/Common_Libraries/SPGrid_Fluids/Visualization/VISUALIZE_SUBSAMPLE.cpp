//#####################################################################
// Copyright (c) 2012, Eftychios Sifakis, Sean Bauer
// Distributed under the FreeBSD license (see license.txt)
//#####################################################################
#include <SPGrid_Fluids/Visualization/VISUALIZE_SUBSAMPLE.h>

#include <SPGrid/Tools/SPGrid_Block_Iterator.h>
#include <SPGrid_Fluids/Grids/GRID_TOPOLOGY_HELPER.h>
#include <SPGrid_Fluids/Read_Write/READ_HIERARCHY.h>
#include <PhysBAM_Tools/Read_Write/Grids_Uniform/READ_WRITE_GRID.h>
#include <PhysBAM_Tools/Read_Write/Grids_Uniform_Arrays/READ_WRITE_ARRAYS.h>
#include <PhysBAM_Tools/Read_Write/Grids_Uniform_Arrays/READ_WRITE_FACE_ARRAYS.h>
#include <Common_Tools/Math_Tools/RANGE_ITERATOR.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <SPGrid_Fluids/Simulation/FLUIDS_SIMULATION_DATA.h>

using namespace PhysBAM;
//#####################################################################
// Visualize_Subsample
//#####################################################################
template<class T_STRUCT,class T,int d> void VISUALIZE_SUBSAMPLE<T_STRUCT,T,d>::
Visualize_Subsample(const int visualization_level,const ARRAY<int>& frames,const std::string input_dir,const std::string output_dir)
{
    // Setup
    RW rw=RW(); STREAM_TYPE stream_type(rw);
    T T_STRUCT::* density_channel=&T_STRUCT::ch3;

    // Read Hierarchy
    T_HIERARCHY hierarchy=READ_HIERARCHY<T_STRUCT,T,d>::Read_Hierarchy(input_dir);

    // Make sure level is valid
    PHYSBAM_ASSERT(visualization_level >=1 && visualization_level <= hierarchy.Levels());

    // Working grid
    T_GRID mac_grid=hierarchy.Grid(visualization_level).Get_MAC_Grid();
    
    // workaround
    TV epsilon=(T).5*hierarchy.Grid(1).dX;

    // Setup output
    ARRAY<T,T_INDEX> density;
    density.Resize(mac_grid.Domain_Indices(3));

    // for subsasmpling
    T_ARRAYS_SCALAR density_weights(density);
    
    // Create dirs
    FILE_UTILITIES::Create_Directory(output_dir);
    FILE_UTILITIES::Create_Directory(output_dir+"/common");

    // iterate over frames
    int output_frame=0;
    for(int i=1;i<=frames.m;i++){
        
        density.Fill(T());
        density_weights.Fill(T());
        
        const int input_frame=frames(i);
        std::string input_f=STRING_UTILITIES::string_sprintf("%d",input_frame);
        std::string output_f=STRING_UTILITIES::string_sprintf("%d",output_frame);
        FILE_UTILITIES::Create_Directory(output_dir+"/"+output_f);        
        hierarchy.Read_Data_Channel(input_dir+"/"+input_f+"/spgrid_density",density_channel);
        
        for(int level=1;level<=hierarchy.Levels();level++){
            Const_data_array_type spgrid_density=hierarchy.Allocator(level).Get_Const_Array(density_channel);
            Const_flag_array_type flags=hierarchy.Allocator(level).Get_Const_Array(&T_STRUCT::flags);    
            RANGE_ITERATOR<d> range_iterator(RANGE<T_INDEX>(T_INDEX(),hierarchy.Allocator(level).Block_Size().template Cast<T_INDEX>()-1));
            T_INDEX base_index;
            for(SPGrid_Block_Iterator<typename Data_array_type::MASK> iterator(hierarchy.Blocks(level));iterator.Valid();iterator.Next()){
                if((iterator.Offset() & 0xfffUL) == 0){
                    range_iterator.Reset();
                    base_index=iterator.Index().template Cast<T_INDEX>();}
                const T_INDEX index=base_index+range_iterator.Index();
                // for finer levels, average to this one
                if(level < visualization_level){
                    T_INDEX coarse_index=hierarchy.Grid(visualization_level).Cell(hierarchy.Grid(level).Node(index)+epsilon,0);
                    if(iterator.Data(flags) & SPGrid_Cell_Type_Interior)
                        if(mac_grid.Cell_Indices(1).Lazy_Inside(coarse_index)){
                            density(coarse_index)+=iterator.Data(spgrid_density);
                            density_weights(coarse_index)++;}

                // for coarser (and same) levels, copy to all cells
                } else {
                    T_INDEX base_fine_index=hierarchy.Grid(visualization_level).Cell(hierarchy.Grid(level).Node(index)+epsilon,0);
                    T_INDEX max_fine_index=hierarchy.Grid(visualization_level).Cell(hierarchy.Grid(level).Node(index+1)+epsilon,0);
                    RANGE<T_INDEX> fine_cells(base_fine_index,max_fine_index-1);
                    for(RANGE_ITERATOR<d> fine_cell_iterator(fine_cells);fine_cell_iterator.Valid();fine_cell_iterator.Next()){
                        const T_INDEX& fine_index=fine_cell_iterator.Index();
                        if(mac_grid.Cell_Indices(1).Lazy_Inside(fine_index))
                            if(iterator.Data(flags) & SPGrid_Cell_Type_Interior)
                                density(fine_index)=iterator.Data(spgrid_density);}
                }
                range_iterator.Next();
            }
        }
        
        // Divide by weights for averaged values
        for(typename GRID<TV>::CELL_ITERATOR iterator(mac_grid);iterator.Valid();iterator.Next())
            if(density_weights(iterator.Cell_Index()))
                density(iterator.Cell_Index())/=density_weights(iterator.Cell_Index());
        
        FILE_UTILITIES::Write_To_File(stream_type,output_dir+"/"+output_f+"/grid",mac_grid);
        FILE_UTILITIES::Write_To_File(stream_type,output_dir+"/common/grid",mac_grid);
        FILE_UTILITIES::Write_To_File(stream_type,output_dir+"/"+output_f+"/density",density);
        
        output_frame++;
    }
}
//#####################################################################
template class VISUALIZE_SUBSAMPLE<FLUIDS_SIMULATION_DATA<float>,float,2>;
template class VISUALIZE_SUBSAMPLE<FLUIDS_SIMULATION_DATA<float>,float,3>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class VISUALIZE_SUBSAMPLE<FLUIDS_SIMULATION_DATA<double>,double,2>;
template class VISUALIZE_SUBSAMPLE<FLUIDS_SIMULATION_DATA<double>,double,3>;
#endif
