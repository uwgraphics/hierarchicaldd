//#####################################################################
//Removed For Anonymity: Copyright Authors of Submission pap173s1
//Removed For Anonymity: Copyright Authors of Submission pap173s1
//#####################################################################
// Subroutine SPGrid_Computations::Correction
//#####################################################################
#ifndef __SPGrid_Correction_h__
#define __SPGrid_Correction_h__

#include <SPGrid/Core/SPGrid_Allocator.h>
#include <SPGrid_Fluids/Grids/GRID_TOPOLOGY_HELPER.h>

#include <iostream>

namespace SPGrid_Computations{

using namespace SPGrid;

template<class T_STRUCT,class T_DATA,class T_FLAGS,int d>
class Correction
{
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<const T_DATA>::type Const_data_array_type;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<T_DATA>::type Data_array_type;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<const T_FLAGS>::type Const_flag_array_type;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<T_FLAGS>::mask Flag_array_mask;

    T_DATA T_STRUCT::* u_field;
    T_DATA T_STRUCT::* r_field;
    T_FLAGS T_STRUCT::* flags_field;
    const bool write_interior;
    const bool write_boundary;
    const T_DATA omega;
public:
    Correction(T_DATA T_STRUCT::* u_field_input,T_DATA T_STRUCT::* r_field_input,T_FLAGS T_STRUCT::* flags_field_input,
        bool write_interior_input=true,bool write_boundary_input=true,T_DATA omega_input=2./3.)
        :u_field(u_field_input),r_field(r_field_input),flags_field(flags_field_input),
        write_interior(write_interior_input),write_boundary(write_boundary_input),omega(omega_input)
    {}
    
     Correction(SPGrid_Allocator<T_STRUCT,d>& allocator,const std::pair<const unsigned long*,unsigned>& blocks,
         T_DATA T_STRUCT::* u_field_input,T_DATA T_STRUCT::* r_field_input,T_FLAGS T_STRUCT::* flags_field_input,
         bool write_interior_input=true,bool write_boundary_input=true,T_DATA omega_input=2./3.)
         :u_field(u_field_input),r_field(r_field_input),flags_field(flags_field_input),
         write_interior(write_interior_input),write_boundary(write_boundary_input),omega(omega_input)
    {Run(allocator,blocks);}
    
    void Run(SPGrid_Allocator<T_STRUCT,d>& allocator,const std::pair<const unsigned long*,unsigned>& blocks) const
    {
        static const int number_of_face_neighbors=GRID_TOPOLOGY_HELPER<Flag_array_mask>::faces_per_cell;
        unsigned long face_neighbor_offsets[number_of_face_neighbors];
        GRID_TOPOLOGY_HELPER<Flag_array_mask>::Face_Neighbor_Offsets(face_neighbor_offsets); // order is (-x, x, -y, y, -z, z)

        Const_data_array_type r=allocator.Get_Const_Array(r_field);
        Data_array_type u=allocator.Get_Array(u_field);
        Const_flag_array_type flags=allocator.Get_Const_Array(flags_field);    
        
        for(SPGrid_Block_Iterator<Flag_array_mask> iterator(blocks);iterator.Valid();iterator.Next()){
            unsigned flag = iterator.Data(flags);

            switch(flag & (SPGrid_Solver_Cell_Type_Active|SPGrid_Solver_Cell_Type_Boundary)){
                case SPGrid_Solver_Cell_Type_Active:
                    if(!write_interior) continue;
                    break;
                case SPGrid_Solver_Cell_Type_Active|SPGrid_Solver_Cell_Type_Boundary:
                    if(!write_boundary) continue;
                    break;
                case SPGrid_Solver_Cell_Type_Boundary:
                    PHYSBAM_FATAL_ERROR();
                default:
                    continue;}
            
            int face_count=0;
            for(int face=0;face<number_of_face_neighbors;face++)
                if(flag & (SPGrid_Solver_Face_Minus_X_Active<<face)) face_count++;
            PHYSBAM_ASSERT(face_count > 0);
            iterator.Data(u) += (omega/(T_DATA)face_count)*iterator.Data(r);
        }
    }
};



template<class T_STRUCT,class T_DATA,class T_FLAGS,int d>
class Correction_Interface
{
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<const T_DATA>::type Const_data_array_type;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<T_DATA>::type Data_array_type;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<const T_FLAGS>::type Const_flag_array_type;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<T_FLAGS>::mask Flag_array_mask;

    T_DATA T_STRUCT::* u_field;
    T_DATA T_STRUCT::* r_field;
    T_FLAGS T_STRUCT::* flags_field;
    const T_DATA omega;
public:
    Correction_Interface(T_DATA T_STRUCT::* u_field_input,T_DATA T_STRUCT::* r_field_input,T_FLAGS T_STRUCT::* flags_field_input,T_DATA omega_input=2./3.)
        :u_field(u_field_input),r_field(r_field_input),flags_field(flags_field_input),omega(omega_input)
    {}
    
    Correction_Interface(SPGrid_Allocator<T_STRUCT,d>& allocator,const std::pair<const unsigned long*,unsigned>& blocks,
                         T_DATA T_STRUCT::* u_field_input,T_DATA T_STRUCT::* r_field_input,T_FLAGS T_STRUCT::* flags_field_input,T_DATA omega_input=2./3.)
        :u_field(u_field_input),r_field(r_field_input),flags_field(flags_field_input),omega(omega_input)
    {Run(allocator,blocks);}
    
    void Run(SPGrid_Allocator<T_STRUCT,d>& allocator,const std::pair<const unsigned long*,unsigned>& blocks) const
    {
        static const int number_of_face_neighbors=GRID_TOPOLOGY_HELPER<Flag_array_mask>::faces_per_cell;
        unsigned long face_neighbor_offsets[number_of_face_neighbors];
        GRID_TOPOLOGY_HELPER<Flag_array_mask>::Face_Neighbor_Offsets(face_neighbor_offsets); // order is (-x, x, -y, y, -z, z)

        Const_data_array_type r=allocator.Get_Const_Array(r_field);
        Data_array_type u=allocator.Get_Array(u_field);
        Const_flag_array_type flags=allocator.Get_Const_Array(flags_field);    
        
        for(SPGrid_Block_Iterator<Flag_array_mask> iterator(blocks);iterator.Valid();iterator.Next()){
            unsigned flag = iterator.Data(flags);

            if(flag&(SPGrid_Solver_Cell_Type_Boundary|SPGrid_Solver_Cell_Type_Interface)){
                int face_count=0;
                for(int face=0;face<number_of_face_neighbors;face++)
                    if(flag&(SPGrid_Solver_Face_Minus_X_Active<<face)) face_count++;
                PHYSBAM_ASSERT(face_count > 0);
                iterator.Data(u) += (omega/(T_DATA)face_count)*iterator.Data(r);
            }
        }
    }
};
//#####################################################################
}
#endif
