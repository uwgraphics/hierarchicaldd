//#####################################################################
//Removed For Anonymity: Copyright Authors of Submission pap173s1
//Removed For Anonymity: Copyright Authors of Submission pap173s1
//#####################################################################
// Class SPGrid_Writer
//#####################################################################
#ifndef __SPGRID_WRITER_H__
#define __SPGRID_WRITER_H__

#include <PhysBAM_Tools/Read_Write/Grids_Uniform/READ_WRITE_GRID.h>

#include <SPGrid_Fluids/Simulation/FLUIDS_SIMULATION_DATA.h>
#include <SPGrid_Fluids/Simulation/FLUIDS_SIMULATION_FLAGS.h>
#include <SPGrid_Fluids/Grids/GRID_HIERARCHY.h>


namespace SPGrid{
using namespace PhysBAM;
template<typename T_STRUCT,typename T_STRUCT_VISUALIZATION,typename T,int d>
class SPGRID_WRITER{
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<T>::mask T_MASK_1;
    typedef typename SPGrid_Allocator<T_STRUCT_VISUALIZATION,d>::template Array<T>::mask T_MASK_2;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<const T>::type SPG_Array_Type_1;
    typedef typename SPGrid_Allocator<T_STRUCT_VISUALIZATION,d>::template Array<T>::type SPG_Array_Type_2;

    static void Convert_SPGrid(SPGrid_Allocator<T_STRUCT,d>& allocator1,SPGrid_Allocator<T_STRUCT_VISUALIZATION,d>& allocator2,
                               SPGrid_Set<typename SPGrid_Allocator<T_STRUCT,d>::template Array<unsigned>::type>& set1,
                               SPGrid_Set<typename SPGrid_Allocator<T_STRUCT_VISUALIZATION,d>::template Array<unsigned>::type>& set2,
                               T T_STRUCT::* ch1,T T_STRUCT_VISUALIZATION::* ch2){
        LOG::SCOPE scope("converting");
        //From 1 to 2.
        SPG_Array_Type_1 Array1=allocator1.Get_Const_Array(ch1);
        SPG_Array_Type_2 Array2=allocator2.Get_Array(ch2);
        for(SPGrid_Block_Iterator<T_MASK_1> iterator(set1.Get_Blocks());iterator.Valid();iterator.Next()){
            Array2(T_MASK_2::Linear_Offset(iterator.Index()))=Array1(iterator.Offset());}
    }
    static void Initialize_SPGrid(SPGrid_Allocator<T_STRUCT,d>& allocator1,SPGrid_Allocator<T_STRUCT_VISUALIZATION,d>& allocator2,
                                  SPGrid_Set<typename SPGrid_Allocator<T_STRUCT,d>::template Array<unsigned>::type>& set1,
                                  SPGrid_Set<typename SPGrid_Allocator<T_STRUCT_VISUALIZATION,d>::template Array<unsigned>::type>& set2){
        LOG::SCOPE scope("Initializing");
        //From 1 to 2.
        for(SPGrid_Block_Iterator<T_MASK_1> iterator(set1.Get_Blocks());iterator.Valid();iterator.Next()){
            if(set1.array(iterator.Offset())&Uniform_Node_Type_Interface/*Adaptive_Node_Type_Active*/)
                set2.Mask(T_MASK_2::Linear_Offset(iterator.Index()+1),SPGrid_Node_Active);
            if(set1.array(iterator.Offset())&Adaptive_Node_Type_T_Junction)
                set2.Mask(T_MASK_2::Linear_Offset(iterator.Index()+1),SPGrid_Node_T_Junction);
            if(set1.array(iterator.Offset())&Adaptive_Node_Type_Coarse_Shared)
                set2.Mask(T_MASK_2::Linear_Offset(iterator.Index()+1),SPGrid_Node_Coarse_Shared);
            if(set1.array(iterator.Offset())&Uniform_Node_Type_Dirichlet){
                set2.Mask(T_MASK_2::Linear_Offset(iterator.Index()),SPGrid_Cell_Type_Dirichlet);}
            if(set1.array(iterator.Offset())&(Adaptive_Cell_Type_Interior)/*Uniform_Node_Type_Dirichlet*/){
                set2.Mask(T_MASK_2::Linear_Offset(iterator.Index()),SPGrid_Cell_Type_Interior|SPGrid_Cell_Type_Active);}
        }
        set2.Refresh_Block_Offsets();
    }

public:
    static void Write_SPGrid(STREAM_TYPE stream_type,GRID<VECTOR<T,d> >& grid,GRID_HIERARCHY<T_STRUCT,T,d>& hierarchy,const std::string output_directory,const int frame)
    {
        GRID_HIERARCHY<T_STRUCT_VISUALIZATION,T,d> hierarchy_visualization(grid,hierarchy.Levels());
        hierarchy_visualization.Initialize_Grids();
        hierarchy_visualization.Initialize_Allocators();
        hierarchy_visualization.Initialize_Sets();
        for(int level=1;level<=hierarchy.Levels();level++){
            Initialize_SPGrid(hierarchy.Allocator(level),hierarchy_visualization.Allocator(level),
                              hierarchy.Set(level),hierarchy_visualization.Set(level));}

        std::string f=STRING_UTILITIES::string_sprintf("%d",frame);
        FILE_UTILITIES::Create_Directory(output_directory);
        FILE_UTILITIES::Create_Directory(output_directory+"/common");
        FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/common/fine_grid",grid);
        FILE_UTILITIES::Create_Directory(output_directory+"/"+f);
        FILE_UTILITIES::Write_To_Text_File(output_directory+"/"+f+"/levels",STRING_UTILITIES::string_sprintf("%d",hierarchy_visualization.Levels()));
        hierarchy_visualization.Write_Flags_Channel(output_directory+"/"+f+"/flags");
        hierarchy_visualization.Write_Block_Offsets(output_directory+"/"+f+"/block_offsets");
        //SPGrid_hierarchy.Write_Data_Channel(output_directory+"/"+f+"/spgrid_density",density_field);
    }
};

}

#endif
