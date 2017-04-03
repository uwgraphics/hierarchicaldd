//#####################################################################
// Copyright 2013, Raj Setaluri, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Subroutine SPGrid_Print_Hierarchy::Print_Hierarchy
//#####################################################################
#ifndef __SPGrid_Print_Hierarchy_h__
#define __SPGrid_Print_Hierarchy_h__

#include <SPGrid_Fluids/Grids/GRID_HIERARCHY.h>
#include <SPGrid_Fluids/Read_Write/READ_HIERARCHY.h>
#include <SPGrid/Tools/SPGrid_Block_Iterator.h>
#include <SPGrid_Fluids/Simulation/FLUIDS_SIMULATION_FLAGS.h>

namespace SPGrid_Print_Hierarchy{

using namespace SPGrid;
using namespace PhysBAM;

template<class T_STRUCT,class T,int d> void
Read_And_Print_Hierarchy(const std::string& base_dir)
{
    typedef GRID_HIERARCHY<T_STRUCT,T,d> T_HIERARCHY;    
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<const unsigned>::type Const_flag_array_type;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<unsigned>::type Flag_array_type;

    T_HIERARCHY hierarchy=READ_HIERARCHY<T_STRUCT,T,d>::Read_Hierarchy(base_dir);
    Print_Hierarchy(hierarchy);
}

template<class T_STRUCT,class T,int d> void
Print_Hierarchy(GRID_HIERARCHY<T_STRUCT,T,d>& hierarchy)
{
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<const unsigned>::type Const_flag_array_type;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<unsigned>::type Flag_array_type;

    const int levels=hierarchy.Levels();
    ARRAY<int> number_of_interior_cells(levels);number_of_interior_cells.Fill(0);
    ARRAY<int> number_of_active_cells(levels);number_of_active_cells.Fill(0);

    for(int level=1;level<=levels;level++){
        Const_flag_array_type flags=hierarchy.Allocator(level).Get_Const_Array(&T_STRUCT::flags);
        for(SPGrid_Block_Iterator<typename Flag_array_type::MASK> iterator(hierarchy.Blocks(level));iterator.Valid();iterator.Next()){
            if(iterator.Data(flags) & SPGrid_Cell_Type_Interior)
                number_of_interior_cells(level)++;
            if(iterator.Data(flags) & SPGrid_Cell_Type_Active)
                number_of_active_cells(level)++;}}

    int total_interior_cells=0; for(int level=1;level<=levels;level++) total_interior_cells+=number_of_interior_cells(level);
    int total_active_cells=0; for(int level=1;level<=levels;level++) total_active_cells+=number_of_active_cells(level);
    
    {
    LOG::SCOPE scope("Hierarchy");
    LOG::cout<<"Number of levels               : "<<STRING_UTILITIES::string_sprintf("%d",levels)<<std::endl;
    LOG::cout<<"Total number of Interior Cells : "<<total_interior_cells<<std::endl;
    LOG::cout<<"Total number of Active Cells   : "<<total_active_cells<<std::endl;

    for(int level=1;level<=levels;level++){
        LOG::SCOPE scope("Level : "+STRING_UTILITIES::string_sprintf("%d",level));
        LOG::cout<<"Grid Domain      : "<<hierarchy.Grid(level).Domain()<<std::endl;
        LOG::cout<<"Grid Cell Counts : "<<hierarchy.Grid(level).Numbers_Of_Cells()<<std::endl;
        LOG::cout<<"Interior Cells   : "<<number_of_interior_cells(level)<<std::endl;
        LOG::cout<<"Active Cells     : "<<number_of_active_cells(level)<<std::endl;
        LOG::cout<<"Number of Blocks : "<<hierarchy.Blocks(level).second<<std::endl;}
    }
    
    LOG::cout<<std::endl;
}
//#####################################################################
}
#endif
