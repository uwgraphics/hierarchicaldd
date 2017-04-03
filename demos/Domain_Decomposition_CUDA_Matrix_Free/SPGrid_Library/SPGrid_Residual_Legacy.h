//#####################################################################
//Removed For Anonymity: Copyright Authors of Submission pap173s1
//Removed For Anonymity: Copyright Authors of Submission pap173s1
//#####################################################################
// Subroutine SPGrid_Computations::Residual
//#####################################################################
#ifndef __SPGrid_Residual_h__
#define __SPGrid_Residual_h__

#include <SPGrid/Core/SPGrid_Allocator.h>
#include <SPGrid_Fluids/Grids/GRID_TOPOLOGY_HELPER.h>

namespace SPGrid_Computations{

using namespace SPGrid;

template<class T_STRUCT,class T_DATA,class T_FLAGS,int d>
class Residual
{
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<const T_DATA>::type Const_data_array_type;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<T_DATA>::type Data_array_type;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<const T_FLAGS>::type Const_flag_array_type;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<T_FLAGS>::mask Flag_array_mask;

    T_DATA T_STRUCT::* u_field;
    T_DATA T_STRUCT::* r_field;
    T_DATA T_STRUCT::* b_field;
    T_FLAGS T_STRUCT::* flags_field;

public:
    Residual(
        T_DATA T_STRUCT::* b_field_input,T_DATA T_STRUCT::* u_field_input,T_DATA T_STRUCT::* r_field_input,T_FLAGS T_STRUCT::* flags_field_input)
        :b_field(b_field_input),u_field(u_field_input),r_field(r_field_input),flags_field(flags_field_input)
    {}
    
     Residual(SPGrid_Allocator<T_STRUCT,d>& allocator,const std::pair<const unsigned long*,unsigned>& blocks,
              T_DATA T_STRUCT::* b_field_input,T_DATA T_STRUCT::* u_field_input,T_DATA T_STRUCT::* r_field_input,T_FLAGS T_STRUCT::* flags_field_input)
         :b_field(b_field_input),u_field(u_field_input),r_field(r_field_input),flags_field(flags_field_input)
    {Run(allocator,blocks);}
    
    void Run(SPGrid_Allocator<T_STRUCT,d>& allocator,const std::pair<const unsigned long*,unsigned>& blocks) const
    {
        static const int number_of_face_neighbors=GRID_TOPOLOGY_HELPER<Flag_array_mask>::faces_per_cell;
        unsigned long face_neighbor_offsets[number_of_face_neighbors];
        GRID_TOPOLOGY_HELPER<Flag_array_mask>::Face_Neighbor_Offsets(face_neighbor_offsets); // order is (-x, x, -y, y, -z, z)

        Const_data_array_type u=allocator.Get_Const_Array(u_field);
        Const_data_array_type b=allocator.Get_Const_Array(b_field);
        Data_array_type r=allocator.Get_Array(r_field);
        Const_flag_array_type flags=allocator.Get_Const_Array(flags_field);    
        
        for(SPGrid_Block_Iterator<Flag_array_mask> iterator(blocks);iterator.Valid();iterator.Next()){
            unsigned flag = iterator.Data(flags);
            if(flag & SPGrid_Solver_Cell_Type_Active){
                T_DATA cell_value=iterator.Data(u),result=0;
                for(int face=0;face<number_of_face_neighbors;face++){
                    unsigned long offset = face_neighbor_offsets[face];                    
                    if(flag & (SPGrid_Solver_Face_Minus_X_Active<<face)){
                        T_DATA neighbor_value=iterator.Data(u,offset);
                        result -= (neighbor_value - cell_value);}}
                iterator.Data(r) = iterator.Data(b) - result;}}
    }
};

template<class T_STRUCT,class T_DATA,class T_FLAGS,int d>
class Residual_Global
{
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<const T_DATA>::type Const_data_array_type;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<T_DATA>::type Data_array_type;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<const T_FLAGS>::type Const_flag_array_type;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<T_FLAGS>::mask Flag_array_mask;

    T_DATA T_STRUCT::* u_field;
    T_DATA T_STRUCT::* r_field;
    T_DATA T_STRUCT::* b_field;
    T_FLAGS T_STRUCT::* flags_field;

public:
    Residual_Global(
        T_DATA T_STRUCT::* b_field_input,T_DATA T_STRUCT::* u_field_input,T_DATA T_STRUCT::* r_field_input,T_FLAGS T_STRUCT::* flags_field_input)
        :b_field(b_field_input),u_field(u_field_input),r_field(r_field_input),flags_field(flags_field_input)
    {}
    
     Residual_Global(SPGrid_Allocator<T_STRUCT,d>& allocator,const std::pair<const unsigned long*,unsigned>& blocks,
              T_DATA T_STRUCT::* b_field_input,T_DATA T_STRUCT::* u_field_input,T_DATA T_STRUCT::* r_field_input,T_FLAGS T_STRUCT::* flags_field_input)
         :b_field(b_field_input),u_field(u_field_input),r_field(r_field_input),flags_field(flags_field_input)
    {Run(allocator,blocks);}
    
    void Run(SPGrid_Allocator<T_STRUCT,d>& allocator,const std::pair<const unsigned long*,unsigned>& blocks) const
    {
        static const int number_of_face_neighbors=GRID_TOPOLOGY_HELPER<Flag_array_mask>::faces_per_cell;
        unsigned long face_neighbor_offsets[number_of_face_neighbors];
        GRID_TOPOLOGY_HELPER<Flag_array_mask>::Face_Neighbor_Offsets(face_neighbor_offsets); // order is (-x, x, -y, y, -z, z)

        Const_data_array_type u=allocator.Get_Const_Array(u_field);
        Const_data_array_type b=allocator.Get_Const_Array(b_field);
        Data_array_type r=allocator.Get_Array(r_field);
        Const_flag_array_type flags=allocator.Get_Const_Array(flags_field);    
        
        for(SPGrid_Block_Iterator<Flag_array_mask> iterator(blocks);iterator.Valid();iterator.Next()){
            unsigned flag = iterator.Data(flags);
            if(flag & (SPGrid_Solver_Cell_Type_Active | SPGrid_Solver_Cell_Type_Interface)){             
                T_DATA cell_value=iterator.Data(u),result=0;
                for(int face=0;face<number_of_face_neighbors;face++){
                    unsigned long offset = face_neighbor_offsets[face];                    
                    if(flag & (SPGrid_Solver_Face_Minus_X_Active<<face)){
                        T_DATA neighbor_value=iterator.Data(u,offset);
                        result -= (neighbor_value - cell_value);}}
                iterator.Data(r) = iterator.Data(b) - result;}}
    }
};
//#####################################################################
}
#endif
