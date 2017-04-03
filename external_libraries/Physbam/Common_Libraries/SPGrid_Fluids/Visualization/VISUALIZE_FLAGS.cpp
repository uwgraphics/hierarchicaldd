//#####################################################################
// Copyright (c) 2012, Eftychios Sifakis, Sean Bauer
// Distributed under the FreeBSD license (see license.txt)
//#####################################################################
#include <SPGrid_Fluids/Visualization/VISUALIZE_FLAGS.h>
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

#define DEBUG

using namespace PhysBAM;
//#####################################################################
// Visualize_Flags
//#####################################################################
template<class T_STRUCT,class T,int d> void VISUALIZE_FLAGS<T_STRUCT,T,d>::
Visualize_Flags(T_HIERARCHY& hierarchy,const std::string directory)
{
    typedef float RW;
    RW rw=RW(); STREAM_TYPE stream_type(rw);
    for(int level=1;level<=hierarchy.Levels();level++)
        Write_Output(stream_type,hierarchy,directory,level);
}
template<class T_STRUCT,class T,int d> void VISUALIZE_FLAGS<T_STRUCT,T,d>::
Write_Output(STREAM_TYPE stream_type,T_HIERARCHY& hierarchy,const std::string directory,const int level)
{
    typedef ARRAY<T,T_INDEX> T_ARRAYS_SCALAR;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::FACE_ARRAYS::template REBIND<bool>::TYPE T_FACE_ARRAYS;
    typedef ARRAY<bool,T_INDEX> T_ARRAYS_BOOL;

    Flag_array_type flags=hierarchy.Allocator(level).Get_Array(&T_STRUCT::flags);

    T_GRID mac_grid=hierarchy.Grid(level).Get_MAC_Grid();

    T_ARRAYS_SCALAR cell_type(mac_grid.Domain_Indices(3));
    T_FACE_ARRAYS active_faces(mac_grid,1);
   
    T_FACE_ARRAYS valid_faces(mac_grid,1);

    T_ARRAYS_BOOL dirichlet_cells(mac_grid.Domain_Indices(3)); dirichlet_cells.Fill(false);
    T_ARRAYS_BOOL active_cells(mac_grid.Domain_Indices(3)); active_cells.Fill(false);
    
    T_ARRAYS_BOOL active_nodes(mac_grid.Domain_Indices(3));
    T_ARRAYS_BOOL coarse_shared_nodes(mac_grid.Domain_Indices(3));

    T_FACE_ARRAYS active_faces_pm(mac_grid,1);
#ifdef DEBUG
    T_FACE_ARRAYS active_faces_p(mac_grid,1);
    T_FACE_ARRAYS active_faces_m(mac_grid,1);
#endif

    T_FACE_ARRAYS scaled_faces_pm(mac_grid,1);
#ifdef DEBUG
    T_FACE_ARRAYS scaled_faces_p(mac_grid,1);
    T_FACE_ARRAYS scaled_faces_m(mac_grid,1);
#endif

    for(SPGrid_Block_Iterator<typename Flag_array_type::MASK> iterator(hierarchy.Blocks(level));iterator.Valid();iterator.Next()){

        T_INDEX index=iterator.Index().template Cast<T_INDEX>();
        unsigned cell_flags=iterator.Data(flags);
        
        // Interior, Ghost, Exterior Cells
        if(mac_grid.Domain_Indices(3).Lazy_Inside(index)){
            if(cell_flags & SPGrid_Cell_Type_Interior)
                cell_type(index)=(T).85;
            else if(cell_flags & SPGrid_Cell_Type_Ghost)
                cell_type(index)=(T).55;
            else
                cell_type(index)=(T).4;}
        
#if 1
        // Face P/M Active
        if(cell_flags & SPGrid_Face_Minus_X_Active)
            active_faces_pm(1,index)=true;
        if(cell_flags & SPGrid_Face_Plus_X_Active)
            active_faces_pm(1,index+T_INDEX::Axis_Vector(1))=true;

        if(cell_flags & SPGrid_Face_Minus_Y_Active)
            active_faces_pm(2,index)=true;
        if(cell_flags & SPGrid_Face_Plus_Y_Active)
            active_faces_pm(2,index+T_INDEX::Axis_Vector(2))=true;

        if(cell_flags & SPGrid_Face_Minus_Z_Active)
            active_faces_pm(3,index)=true;
        if(cell_flags & SPGrid_Face_Plus_Z_Active)
            active_faces_pm(3,index+T_INDEX::Axis_Vector(3))=true;
#ifdef DEBUG
        if(cell_flags & SPGrid_Face_Minus_X_Active)
            active_faces_m(1,index)=true;
        if(cell_flags & SPGrid_Face_Plus_X_Active)
            active_faces_p(1,index)=true;

        if(cell_flags & SPGrid_Face_Minus_Y_Active)
            active_faces_m(2,index)=true;
        if(cell_flags & SPGrid_Face_Plus_Y_Active)
            active_faces_p(2,index)=true;

        if(cell_flags & SPGrid_Face_Minus_Z_Active)
            active_faces_m(3,index)=true;
        if(cell_flags & SPGrid_Face_Plus_Z_Active)
            active_faces_p(3,index)=true;
#endif

        // Face P/M Scaled
        if(cell_flags & SPGrid_Face_Minus_X_Scaled)
            scaled_faces_pm(1,index)=true;
        if(cell_flags & SPGrid_Face_Plus_X_Scaled)
            scaled_faces_pm(1,index+T_INDEX::Axis_Vector(1))=true;

        if(cell_flags & SPGrid_Face_Minus_Y_Scaled)
            scaled_faces_pm(2,index)=true;
        if(cell_flags & SPGrid_Face_Plus_Y_Scaled)
            scaled_faces_pm(2,index+T_INDEX::Axis_Vector(2))=true;

        if(cell_flags & SPGrid_Face_Minus_Z_Scaled)
            scaled_faces_pm(3,index)=true;
        if(cell_flags & SPGrid_Face_Plus_Z_Scaled)
            scaled_faces_pm(3,index+T_INDEX::Axis_Vector(3))=true;
#ifdef DEBUG
        if(cell_flags & SPGrid_Face_Minus_X_Scaled)
            scaled_faces_m(1,index)=true;
        if(cell_flags & SPGrid_Face_Plus_X_Scaled)
            scaled_faces_p(1,index)=true;

        if(cell_flags & SPGrid_Face_Minus_Y_Scaled)
            scaled_faces_m(2,index)=true;
        if(cell_flags & SPGrid_Face_Plus_Y_Scaled)
            scaled_faces_p(2,index)=true;

        if(cell_flags & SPGrid_Face_Minus_Z_Scaled)
            scaled_faces_m(3,index)=true;
        if(cell_flags & SPGrid_Face_Plus_Z_Scaled)
            scaled_faces_p(3,index)=true;
#endif


        // Face Active
        if(cell_flags & SPGrid_Face_Type_X_Active)
            active_faces(1,index)=true;
        
        if(cell_flags & SPGrid_Face_Type_Y_Active)
            active_faces(2,index)=true;

        if(cell_flags & SPGrid_Face_Type_Z_Active)
            active_faces(3,index)=true;

        // Face Valid
        if(cell_flags & SPGrid_Face_Type_X_Valid)
            valid_faces(1,index)=true;

        if(cell_flags & SPGrid_Face_Type_Y_Valid)
            valid_faces(2,index)=true;

        if(cell_flags & SPGrid_Face_Type_Z_Valid)
            valid_faces(3,index)=true;
#endif        
        // Dirichlet Cells
        if(cell_flags & SPGrid_Cell_Type_Dirichlet)
            dirichlet_cells(index)=true;

        // Active Cells
        if(cell_flags & SPGrid_Cell_Type_Active)
            active_cells(index)=true;

        // Active Nodes
        if(cell_flags & SPGrid_Node_Active)
            active_nodes(index)=true;

        // Coarse Shared Nodes
        if(cell_flags & SPGrid_Node_Coarse_Shared)
            coarse_shared_nodes(index)=true;

    }
    
#ifdef DEBUG    
    // *** BEGIN DEBUG *** 
    VECTOR<RANGE<T_INDEX>,d> Face_Indices(mac_grid.Face_Indices(0));
    for(SPGrid_Block_Iterator<typename Flag_array_type::MASK> iterator(hierarchy.Blocks(level));iterator.Valid();iterator.Next()){        
        T_INDEX index=iterator.Index().template Cast<T_INDEX>();
        unsigned cell_flags=iterator.Data(flags);
        
        for(int axis=1;axis<=d;axis++)
            if(Face_Indices(axis).Lazy_Inside(index)){
                PHYSBAM_ASSERT(active_faces_p(axis,index)==active_faces_m(axis,index+T_INDEX::Axis_Vector(axis)));
                PHYSBAM_ASSERT(scaled_faces_p(axis,index)==scaled_faces_m(axis,index+T_INDEX::Axis_Vector(axis)));}
    }
    // *** END DEBUG *** 
#endif

    std::string l_str=STRING_UTILITIES::string_sprintf("%d",level);
    FILE_UTILITIES::Create_Directory(directory+"/"+l_str+"/common");
    FILE_UTILITIES::Write_To_File(stream_type,directory+"/"+l_str+"/common/grid",mac_grid);
    FILE_UTILITIES::Create_Directory(directory+"/"+l_str+"/0");
    FILE_UTILITIES::Write_To_File(stream_type,directory+"/"+l_str+"/0/density",cell_type);
    FILE_UTILITIES::Write_To_File(stream_type,directory+"/"+l_str+"/0/psi_N",active_faces);
    FILE_UTILITIES::Write_To_File(stream_type,directory+"/"+l_str+"/0/psi_D",dirichlet_cells);
    FILE_UTILITIES::Create_Directory(directory+"/"+l_str+"/1");
    FILE_UTILITIES::Write_To_File(stream_type,directory+"/"+l_str+"/1/density",cell_type);
    FILE_UTILITIES::Write_To_File(stream_type,directory+"/"+l_str+"/1/psi_N",active_faces_pm);
    FILE_UTILITIES::Write_To_File(stream_type,directory+"/"+l_str+"/1/psi_D",dirichlet_cells);
    FILE_UTILITIES::Create_Directory(directory+"/"+l_str+"/2");
    FILE_UTILITIES::Write_To_File(stream_type,directory+"/"+l_str+"/2/density",cell_type);
    FILE_UTILITIES::Write_To_File(stream_type,directory+"/"+l_str+"/2/psi_N",valid_faces);
    FILE_UTILITIES::Write_To_File(stream_type,directory+"/"+l_str+"/2/psi_D",dirichlet_cells);
    FILE_UTILITIES::Create_Directory(directory+"/"+l_str+"/3");
    FILE_UTILITIES::Write_To_File(stream_type,directory+"/"+l_str+"/3/density",cell_type);
    FILE_UTILITIES::Write_To_File(stream_type,directory+"/"+l_str+"/3/psi_N",valid_faces);
    FILE_UTILITIES::Write_To_File(stream_type,directory+"/"+l_str+"/3/psi_D",active_cells);
    FILE_UTILITIES::Create_Directory(directory+"/"+l_str+"/4");
    FILE_UTILITIES::Write_To_File(stream_type,directory+"/"+l_str+"/4/density",cell_type);
    FILE_UTILITIES::Write_To_File(stream_type,directory+"/"+l_str+"/4/psi_N",valid_faces);
    FILE_UTILITIES::Write_To_File(stream_type,directory+"/"+l_str+"/4/psi_D",active_nodes);
    FILE_UTILITIES::Create_Directory(directory+"/"+l_str+"/5");
    FILE_UTILITIES::Write_To_File(stream_type,directory+"/"+l_str+"/5/density",cell_type);
    FILE_UTILITIES::Write_To_File(stream_type,directory+"/"+l_str+"/5/psi_N",valid_faces);
    FILE_UTILITIES::Write_To_File(stream_type,directory+"/"+l_str+"/5/psi_D",coarse_shared_nodes);
    FILE_UTILITIES::Create_Directory(directory+"/"+l_str+"/6");
    FILE_UTILITIES::Write_To_File(stream_type,directory+"/"+l_str+"/6/density",cell_type);
    FILE_UTILITIES::Write_To_File(stream_type,directory+"/"+l_str+"/6/psi_N",scaled_faces_pm);
    FILE_UTILITIES::Write_To_File(stream_type,directory+"/"+l_str+"/6/psi_D",active_cells);
}
//#####################################################################
template class VISUALIZE_FLAGS<FLUIDS_SIMULATION_DATA<float>,float,2>;
template class VISUALIZE_FLAGS<FLUIDS_SIMULATION_DATA<float>,float,3>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class VISUALIZE_FLAGS<FLUIDS_SIMULATION_DATA<double>,double,2>;
template class VISUALIZE_FLAGS<FLUIDS_SIMULATION_DATA<double>,double,3>;
#endif
