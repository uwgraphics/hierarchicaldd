//#####################################################################
//Removed For Anonymity: Copyright Authors of Submission pap173s1
//Removed For Anonymity: Copyright Authors of Submission pap173s1
//#####################################################################
// Subroutine SPGrid_Computations::Laplace
//#####################################################################
#ifndef __SPGrid_Laplace_h__
#define __SPGrid_Laplace_h__

#include <SPGrid/Core/SPGrid_Allocator.h>
#include <SPGrid_Fluids/Grids/GRID_TOPOLOGY_HELPER.h>
#include <SPGrid/Tools/SPGrid_Block_Iterator.h>
#include <Common_Tools/Math_Tools/RANGE_ITERATOR.h>
#include <SPGrid_Fluids/Solvers/Domain_Decomposition/SPGRID_MULTIGRID_FLAGS.h>

namespace SPGrid_Computations{

using namespace PhysBAM;
using namespace SPGrid;

template<class T_STRUCT,class T_DATA,class T_FLAGS,int d>
class Laplace
{
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<const T_DATA>::type Const_data_array_type;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<T_DATA>::type Data_array_type;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<const T_FLAGS>::type Const_flag_array_type;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<T_FLAGS>::mask Flag_array_mask;

    T_DATA T_STRUCT::* u_field;
    T_DATA T_STRUCT::* Lu_field;
    T_FLAGS T_STRUCT::* flags_field;

public:
    Laplace(T_DATA T_STRUCT::* u_field_input,T_DATA T_STRUCT::* Lu_field_input,T_FLAGS T_STRUCT::* flags_field_input)
        :u_field(u_field_input),Lu_field(Lu_field_input),flags_field(flags_field_input)
    {}
    
    Laplace(SPGrid_Allocator<T_STRUCT,d>& allocator,const std::pair<const unsigned long*,unsigned>& blocks,
            T_DATA T_STRUCT::* u_field_input,T_DATA T_STRUCT::* Lu_field_input,T_FLAGS T_STRUCT::* flags_field_input)
        :u_field(u_field_input),Lu_field(Lu_field_input),flags_field(flags_field_input)
    {Run(allocator,blocks);}
    
    void Run(SPGrid_Allocator<T_STRUCT,d>& allocator,const std::pair<const unsigned long*,unsigned>& blocks) const
    {
        static const int number_of_face_neighbors=GRID_TOPOLOGY_HELPER<Flag_array_mask>::faces_per_cell;
        unsigned long face_neighbor_offsets[number_of_face_neighbors];
        GRID_TOPOLOGY_HELPER<Flag_array_mask>::Face_Neighbor_Offsets(face_neighbor_offsets); // order is (-x, x, -y, y, -z, z)

        Const_data_array_type u=allocator.Get_Const_Array(u_field);
        Data_array_type Lu=allocator.Get_Array(Lu_field);
        Const_flag_array_type flags=allocator.Get_Const_Array(flags_field);    
        
        for (SPGrid_Block_Iterator<Flag_array_mask> iterator(blocks);iterator.Valid();iterator.Next())
            {
            unsigned flag = iterator.Data(flags);
            if(flag & (SPGrid_Solver_Cell_Type_Active | SPGrid_Solver_Cell_Type_Interface)){             
                double cell_value=(double)(iterator.Data(u));
                double result=(double)0.;
                for(int face=0;face<number_of_face_neighbors;face++){
                    unsigned long offset = face_neighbor_offsets[face];
                    if(flag & (SPGrid_Solver_Face_Minus_X_Active<<face)){
                        double neighbor_value=(double)(iterator.Data(u,offset));
                        result -= (neighbor_value - cell_value);}
                }
                iterator.Data(Lu) = (T_DATA)result;
            }
        }
    }
};
template<class T_STRUCT,class T_DATA,class T_FLAGS,int d>
class Laplace_ri
{
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<const T_DATA>::type Const_data_array_type;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<T_DATA>::type Data_array_type;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<const T_FLAGS>::type Const_flag_array_type;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<T_FLAGS>::mask Flag_array_mask;

    T_DATA T_STRUCT::* u_field;
    T_DATA T_STRUCT::* Lu_field;
    T_FLAGS T_STRUCT::* flags_field;

public:
    Laplace_ri(T_DATA T_STRUCT::* u_field_input,T_DATA T_STRUCT::* Lu_field_input,T_FLAGS T_STRUCT::* flags_field_input)
        :u_field(u_field_input),Lu_field(Lu_field_input),flags_field(flags_field_input)
    {}
    
    Laplace_ri(SPGrid_Allocator<T_STRUCT,d>& allocator,const std::pair<const unsigned long*,unsigned>& blocks,
            T_DATA T_STRUCT::* u_field_input,T_DATA T_STRUCT::* Lu_field_input,T_FLAGS T_STRUCT::* flags_field_input)
        :u_field(u_field_input),Lu_field(Lu_field_input),flags_field(flags_field_input)
    {Run(allocator,blocks);}
    
    void Run(SPGrid_Allocator<T_STRUCT,d>& allocator,const std::pair<const unsigned long*,unsigned>& blocks) const
    {
        static const int number_of_face_neighbors=GRID_TOPOLOGY_HELPER<Flag_array_mask>::faces_per_cell;
        unsigned long face_neighbor_offsets[number_of_face_neighbors];
        GRID_TOPOLOGY_HELPER<Flag_array_mask>::Face_Neighbor_Offsets(face_neighbor_offsets); // order is (-x, x, -y, y, -z, z)

        Const_data_array_type u=allocator.Get_Const_Array(u_field);
        Data_array_type Lu=allocator.Get_Array(Lu_field);
        Const_flag_array_type flags=allocator.Get_Const_Array(flags_field);    
        
        for(SPGrid_Block_Iterator<Flag_array_mask> iterator(blocks);iterator.Valid();iterator.Next()){
            unsigned flag = iterator.Data(flags);
            if(flag & SPGrid_Solver_Cell_Type_Interface){             
                PHYSBAM_ASSERT(!(flag & SPGrid_Solver_Cell_Type_Active));
                double cell_value=(double)(iterator.Data(u));
                double result=(double)0.;
                for(int face=0;face<number_of_face_neighbors;face++){
                    unsigned long offset=face_neighbor_offsets[face];                    
                    if(flag & (SPGrid_Solver_Face_Minus_X_Active<<face)){
                        unsigned neighbor_flag=iterator.Data(flags,offset);
                        if(neighbor_flag&SPGrid_Solver_Cell_Type_Active){             
                            PHYSBAM_ASSERT(!(neighbor_flag & SPGrid_Solver_Cell_Type_Interface));
                            double neighbor_value=(double)(iterator.Data(u,offset));
                            result-=neighbor_value;}}}
                iterator.Data(Lu)-=(T_DATA)result;}}
    }
};
template<class T_STRUCT,class T_DATA,class T_FLAGS,int d>
class Laplace_ri_Subdomain
{
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<const T_DATA>::type Const_data_array_type;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<T_DATA>::type Data_array_type;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<const T_FLAGS>::type Const_flag_array_type;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<T_FLAGS>::mask Flag_array_mask;

    T_DATA T_STRUCT::* u_field;
    T_DATA T_STRUCT::* Lu_field;
    T_FLAGS T_STRUCT::* flags_field;
    std_array<int,d> subdomain_size;

public:
    Laplace_ri_Subdomain(T_DATA T_STRUCT::* u_field_input,T_DATA T_STRUCT::* Lu_field_input,T_FLAGS T_STRUCT::* flags_field_input,
                         std_array<int,d>& subdomain_size_input)
        :u_field(u_field_input),Lu_field(Lu_field_input),flags_field(flags_field_input),subdomain_size(subdomain_size_input)
    {}
    
    Laplace_ri_Subdomain(SPGrid_Allocator<T_STRUCT,d>& allocator,const std::pair<const unsigned long*,unsigned>& blocks,
                         T_DATA T_STRUCT::* u_field_input,T_DATA T_STRUCT::* Lu_field_input,T_FLAGS T_STRUCT::* flags_field_input,
                         std_array<int,d>& subdomain_size_input)
        :u_field(u_field_input),Lu_field(Lu_field_input),flags_field(flags_field_input),subdomain_size(subdomain_size_input)
    {Run(allocator,blocks);}
    
    static bool Is_Interface(const std_array<int,d>& index,const std_array<int,d>& subdomain_size)
    {
        for(int axis=0;axis<d;++axis){if(index(axis)%subdomain_size(axis)==0) return true;}
        return false;
    }

    void Run(SPGrid_Allocator<T_STRUCT,d>& allocator,const std::pair<const unsigned long*,unsigned>& blocks) const
    {
        static const int number_of_face_neighbors=GRID_TOPOLOGY_HELPER<Flag_array_mask>::faces_per_cell;
        unsigned long face_neighbor_offsets[number_of_face_neighbors];
        GRID_TOPOLOGY_HELPER<Flag_array_mask>::Face_Neighbor_Offsets(face_neighbor_offsets); // order is (-x, x, -y, y, -z, z)

        Const_data_array_type u=allocator.Get_Const_Array(u_field);
        Data_array_type Lu=allocator.Get_Array(Lu_field);
        Const_flag_array_type flags=allocator.Get_Const_Array(flags_field);    
        
        for(SPGrid_Block_Iterator<Flag_array_mask> iterator(blocks);iterator.Valid();iterator.Next()){
            unsigned flag = iterator.Data(flags);
            if(flag&(SPGrid_Solver_Cell_Type_Interface|SPGrid_Solver_Cell_Type_Active)){             
                std_array<int,d> index=iterator.Index();
                if(!Is_Interface(index,subdomain_size)) continue;
                double cell_value=(double)(iterator.Data(u));
                double result=(double)0.;
                for(int face=0;face<number_of_face_neighbors;face++){
                    unsigned long offset=face_neighbor_offsets[face];                    
                    if(flag&(SPGrid_Solver_Face_Minus_X_Active<<face)){
                        unsigned neighbor_flag=iterator.Data(flags,offset);
                        if(neighbor_flag&(SPGrid_Solver_Cell_Type_Interface|SPGrid_Solver_Cell_Type_Active)){             
                            std_array<int,d> neighbor_index=Flag_array_mask::LinearToCoord(Flag_array_mask::Packed_Add(iterator.Offset(),offset));
                            if(Is_Interface(neighbor_index,subdomain_size)) continue;
                            double neighbor_value=(double)(iterator.Data(u,offset));
                            result-=neighbor_value;}}}
                iterator.Data(Lu)-=(T_DATA)result;}}
    }
};

template<class T_STRUCT,class T_DATA,class T_FLAGS,int d>
class Laplace_ir
{
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<const T_DATA>::type Const_data_array_type;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<T_DATA>::type Data_array_type;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<const T_FLAGS>::type Const_flag_array_type;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<T_FLAGS>::mask Flag_array_mask;

    T_DATA T_STRUCT::* u_field;
    T_DATA T_STRUCT::* Lu_field;
    T_FLAGS T_STRUCT::* flags_field;

public:
    Laplace_ir(T_DATA T_STRUCT::* u_field_input,T_DATA T_STRUCT::* Lu_field_input,T_FLAGS T_STRUCT::* flags_field_input)
        :u_field(u_field_input),Lu_field(Lu_field_input),flags_field(flags_field_input)
    {}
    
    Laplace_ir(SPGrid_Allocator<T_STRUCT,d>& allocator,const std::pair<const unsigned long*,unsigned>& blocks,
            T_DATA T_STRUCT::* u_field_input,T_DATA T_STRUCT::* Lu_field_input,T_FLAGS T_STRUCT::* flags_field_input)
        :u_field(u_field_input),Lu_field(Lu_field_input),flags_field(flags_field_input)
    {Run(allocator,blocks);}
    
    void Run(SPGrid_Allocator<T_STRUCT,d>& allocator,const std::pair<const unsigned long*,unsigned>& blocks) const
    {
        static const int number_of_face_neighbors=GRID_TOPOLOGY_HELPER<Flag_array_mask>::faces_per_cell;
        unsigned long face_neighbor_offsets[number_of_face_neighbors];
        GRID_TOPOLOGY_HELPER<Flag_array_mask>::Face_Neighbor_Offsets(face_neighbor_offsets); // order is (-x, x, -y, y, -z, z)

        Const_data_array_type u=allocator.Get_Const_Array(u_field);
        Data_array_type Lu=allocator.Get_Array(Lu_field);
        Const_flag_array_type flags=allocator.Get_Const_Array(flags_field);    
        
        for(SPGrid_Block_Iterator<Flag_array_mask> iterator(blocks);iterator.Valid();iterator.Next()){
            unsigned flag = iterator.Data(flags);
            if(flag & SPGrid_Solver_Cell_Type_Active){             
                PHYSBAM_ASSERT(!(flag & SPGrid_Solver_Cell_Type_Interface));
                double cell_value=(double)(iterator.Data(u));
                double result=(double)0.;
                for(int face=0;face<number_of_face_neighbors;face++){
                    unsigned long offset = face_neighbor_offsets[face];                    
                    if(flag & (SPGrid_Solver_Face_Minus_X_Active<<face)){
                        unsigned neighbor_flag = iterator.Data(flags,offset);
                        if(neighbor_flag & SPGrid_Solver_Cell_Type_Interface){             
                            PHYSBAM_ASSERT(!(neighbor_flag & SPGrid_Solver_Cell_Type_Active));
                            double neighbor_value=(double)(iterator.Data(u,offset));
                            result -= neighbor_value;}}}
                // Note: this is inverted
                iterator.Data(Lu) = -(T_DATA)result;}}
    }
};

template<class T_STRUCT,class T_DATA,class T_FLAGS,int d>
class Laplace_ir_Subdomain
{
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<const T_DATA>::type Const_data_array_type;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<T_DATA>::type Data_array_type;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<const T_FLAGS>::type Const_flag_array_type;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<T_FLAGS>::mask Flag_array_mask;

    T_DATA T_STRUCT::* u_field;
    T_DATA T_STRUCT::* Lu_field;
    T_FLAGS T_STRUCT::* flags_field;
    std_array<int,d> subdomain_size;

public:
    Laplace_ir_Subdomain(T_DATA T_STRUCT::* u_field_input,T_DATA T_STRUCT::* Lu_field_input,T_FLAGS T_STRUCT::* flags_field_input,std_array<int,d>& subdomain_size_input)
        :u_field(u_field_input),Lu_field(Lu_field_input),flags_field(flags_field_input),subdomain_size(subdomain_size_input)
    {}
    
    Laplace_ir_Subdomain(SPGrid_Allocator<T_STRUCT,d>& allocator,const std::pair<const unsigned long*,unsigned>& blocks,
                         T_DATA T_STRUCT::* u_field_input,T_DATA T_STRUCT::* Lu_field_input,T_FLAGS T_STRUCT::* flags_field_input,
                         std_array<int,d>& subdomain_size_input)
        :u_field(u_field_input),Lu_field(Lu_field_input),flags_field(flags_field_input),subdomain_size(subdomain_size_input)
    {Run(allocator,blocks);}
    
    static bool Is_Interface(const std_array<int,d>& index,const std_array<int,d>& subdomain_size)
    {
        for(int axis=0;axis<d;++axis){if(index(axis)%subdomain_size(axis)==0) return true;}
        return false;
    }

    void Run(SPGrid_Allocator<T_STRUCT,d>& allocator,const std::pair<const unsigned long*,unsigned>& blocks) const
    {
        static const int number_of_face_neighbors=GRID_TOPOLOGY_HELPER<Flag_array_mask>::faces_per_cell;
        unsigned long face_neighbor_offsets[number_of_face_neighbors];
        GRID_TOPOLOGY_HELPER<Flag_array_mask>::Face_Neighbor_Offsets(face_neighbor_offsets); // order is (-x, x, -y, y, -z, z)

        Const_data_array_type u=allocator.Get_Const_Array(u_field);
        Data_array_type Lu=allocator.Get_Array(Lu_field);
        Const_flag_array_type flags=allocator.Get_Const_Array(flags_field);    
        
        for(SPGrid_Block_Iterator<Flag_array_mask> iterator(blocks);iterator.Valid();iterator.Next()){
            unsigned flag = iterator.Data(flags);
            if(flag&(SPGrid_Solver_Cell_Type_Active|SPGrid_Solver_Cell_Type_Interface)){             
                std_array<int,d> index=iterator.Index();
                if(Is_Interface(index,subdomain_size)) continue;
                double cell_value=(double)(iterator.Data(u));
                double result=(double)0.;
                for(int face=0;face<number_of_face_neighbors;face++){
                    unsigned long offset = face_neighbor_offsets[face];                    
                    if(flag&(SPGrid_Solver_Face_Minus_X_Active<<face)){
                        unsigned neighbor_flag = iterator.Data(flags,offset);
                        if(neighbor_flag&(SPGrid_Solver_Cell_Type_Interface|SPGrid_Solver_Cell_Type_Active)){
                            std_array<int,d> neighbor_index=Flag_array_mask::LinearToCoord(Flag_array_mask::Packed_Add(iterator.Offset(),offset));
                            if(!Is_Interface(neighbor_index,subdomain_size)) continue;
                            double neighbor_value=(double)(iterator.Data(u,offset));
                            result -= neighbor_value;}}}
                // Note: this is inverted
                iterator.Data(Lu) = -(T_DATA)result;}}
    }
};
//#####################################################################
}
#endif
