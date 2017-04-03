//#####################################################################
//Removed For Anonymity: Copyright Authors of Submission pap173s1
//Removed For Anonymity: Copyright Authors of Submission pap173s1
//#####################################################################
// Subroutine SPGrid_Computations::Flaging_Boundary
//#####################################################################
#ifndef __SPGrid_Flaging_Boundary_h__
#define __SPGrid_Flaging_Boundary_h__

#include <Common_Tools/Math_Tools/RANGE_ITERATOR.h>
#include <SPGrid/Core/SPGrid_Allocator.h>
#include <SPGrid_Fluids/Solvers/Domain_Decomposition/SPGRID_MULTIGRID_FLAGS.h>
#include <algorithm>
namespace SPGrid_Computations{

using namespace SPGrid;
using namespace PhysBAM;

template<class T_STRUCT,int d> class Flaging_Boundary;
template<class T_STRUCT>
class Flaging_Boundary<T_STRUCT,3>
{
    enum{d=3};
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<unsigned>::type Flag_Array_Type;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<unsigned>::mask T_MASK;
    typedef VECTOR<int,d> T_INDEX;

    unsigned T_STRUCT::* flags_field;
    SPGrid_Set<Flag_Array_Type>& set;

public:
    Flaging_Boundary(unsigned T_STRUCT::* flags_field_input,SPGrid_Set<Flag_Array_Type>& set_input)
        :flags_field(flags_field_input),set(set_input)
    {}
    
   Flaging_Boundary(SPGrid_Allocator<T_STRUCT,d>& allocator,const std::pair<const unsigned long*,unsigned>& blocks,unsigned T_STRUCT::* flags_field_input,SPGrid_Set<Flag_Array_Type>& set_input)
       :flags_field(flags_field_input),set(set_input)
    {Run(allocator,blocks);}
    
    void Run(SPGrid_Allocator<T_STRUCT,d>& allocator,const std::pair<const unsigned long*,unsigned>& blocks) const
    {
        Flag_Array_Type flags=allocator.Get_Array(flags_field);    
        enum{
            block_xsize = 1u << T_MASK::block_xbits,
            block_ysize = 1u << T_MASK::block_ybits,
            block_zsize = 1u << T_MASK::block_zbits  
        };
        T_INDEX block_size=T_INDEX(block_xsize,block_ysize,block_zsize);
        typedef unsigned (&block_flag)[block_xsize][block_ysize][block_zsize]; 
        for(SPGrid_Block_Iterator<T_MASK> iterator(blocks);iterator.Valid();iterator.Next_Block()){
            T_INDEX base_index = iterator.Index().template Cast<T_INDEX>();
            ARRAY<unsigned,T_INDEX> flag_buffer(RANGE<T_INDEX>(T_INDEX()-4,block_size+4));
            ARRAY<unsigned,T_INDEX> distance_buffer(RANGE<T_INDEX>(T_INDEX()-4,block_size+4));
            unsigned long block_base_offsets[3][3][3];
            for(int i = -1;i < 2;++i)
            for(int j = -1;j < 2;++j)
            for(int k = -1;k < 2;++k){
                PHYSBAM_ASSERT(T_MASK::Linear_Offset(std_array<unsigned int,d>(base_index(1) * i * block_xsize,base_index(2) * j * block_ysize,base_index(3) * k * block_zsize))%4096 == 0);
                std_array<int,d> coord(base_index(1) +  i * block_xsize,base_index(2) + j * block_ysize,base_index(3) + k * block_zsize);
                unsigned long block_offset = T_MASK::Linear_Offset(coord);
                if(set.IsPageActive(block_offset))
                    block_base_offsets[i + 1][j + 1][k + 1] = block_offset;
                else
                    block_base_offsets[i + 1][j + 1][k + 1] = 0xdeadbeef;
            }
            for(RANGE_ITERATOR<d> range_iterator(RANGE<T_INDEX>(T_INDEX()-4,block_size+4));
                range_iterator.Valid();
                range_iterator.Next()){
                T_INDEX index = range_iterator.Index();
                //first determin which block are we pick
                const int block_i = (index.x < 0)?0:((index.x >= block_xsize)?2:1);
                const int block_j = (index.y < 0)?0:((index.y >= block_ysize)?2:1);
                const int block_k = (index.z < 0)?0:((index.z >= block_zsize)?2:1);
                if(block_base_offsets[block_i][block_j][block_k] == 0xdeadbeef) continue;
                const int cell_x = (index.x + block_xsize) % block_xsize;
                const int cell_y = (index.y + block_ysize) % block_ysize;
                const int cell_z = (index.z + block_zsize) % block_zsize;
                flag_buffer(index) = reinterpret_cast<block_flag>(flags(block_base_offsets[block_i][block_j][block_k]))[cell_x][cell_y][cell_z];
            }
            //init the distance field
            for(RANGE_ITERATOR<d> range_iterator(RANGE<T_INDEX>(T_INDEX()-4, block_size+4));
                range_iterator.Valid();
                range_iterator.Next()){
                T_INDEX index = range_iterator.Index();
                if((flag_buffer(index) == 0x0u) || (flag_buffer(index) & (SPGrid_Solver_Cell_Type_Interface | SPGrid_Solver_Cell_Type_Dirichlet)))
                    distance_buffer(index) = 0;
                else
                    distance_buffer(index) = 0xffff;
            }
            for(int dis = 0;dis < 3;++dis){
                for(RANGE_ITERATOR<d> range_iterator(RANGE<T_INDEX>(T_INDEX() - 3, block_size + 3));
                    range_iterator.Valid();
                    range_iterator.Next()){
                    T_INDEX index = range_iterator.Index();
                    if(distance_buffer(index) != 0){
                        for(int v = 1;v <= d;++v) distance_buffer(index) = std::min<unsigned int>(distance_buffer(index),distance_buffer(index+T_INDEX::Axis_Vector(v))+1);
                        for(int v = 1;v <= d;++v) distance_buffer(index) = std::min<unsigned int>(distance_buffer(index),distance_buffer(index-T_INDEX::Axis_Vector(v))+1);
                    }
                }
            }
            unsigned long flag_offset = iterator.Offset();
            for(RANGE_ITERATOR<d> range_iterator(RANGE<T_INDEX>(T_INDEX(), block_size-1));
                range_iterator.Valid();
                range_iterator.Next(),flag_offset+=sizeof(unsigned)){
                if((flags(flag_offset) & SPGrid_Solver_Cell_Type_Active) && (distance_buffer(range_iterator.Index()) <= 3))
                   flags(flag_offset) |= SPGrid_Solver_Cell_Type_Boundary;
            }
        }
    } 
};

template<class T_STRUCT>
class Flaging_Boundary<T_STRUCT,2>
{
    enum{d=2};
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<unsigned>::type Flag_Array_Type;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<unsigned>::mask T_MASK;
    typedef VECTOR<int,d> T_INDEX;

    unsigned T_STRUCT::* flags_field;
    SPGrid_Set<Flag_Array_Type>& set;

public:
    Flaging_Boundary(unsigned T_STRUCT::* flags_field_input,SPGrid_Set<Flag_Array_Type>& set_input)
        :flags_field(flags_field_input),set(set_input)
    {}
    
   Flaging_Boundary(SPGrid_Allocator<T_STRUCT,d>& allocator,const std::pair<const unsigned long*,unsigned>& blocks,unsigned T_STRUCT::* flags_field_input,SPGrid_Set<Flag_Array_Type>& set_input)
       :flags_field(flags_field_input),set(set_input)
    {Run(allocator,blocks);}
    
    void Run(SPGrid_Allocator<T_STRUCT,d>& allocator,const std::pair<const unsigned long*,unsigned>& blocks) const
    {
        Flag_Array_Type flags=allocator.Get_Array(flags_field);    
        enum{
            block_xsize = 1u << T_MASK::block_xbits,
            block_ysize = 1u << T_MASK::block_ybits
        };
        T_INDEX block_size=T_INDEX(block_xsize,block_ysize);
        typedef unsigned (&block_flag)[block_xsize][block_ysize];
        for(SPGrid_Block_Iterator<T_MASK> iterator(blocks);iterator.Valid();iterator.Next_Block()){
            T_INDEX base_index = iterator.Index().template Cast<T_INDEX>();
            ARRAY<unsigned,T_INDEX> flag_buffer(RANGE<T_INDEX>(T_INDEX() - 4,block_size + 4));
            ARRAY<unsigned,T_INDEX> distance_buffer(RANGE<T_INDEX>(T_INDEX() - 4,block_size + 4));
            unsigned long block_base_offsets[3][3];
            for(int i = -1;i < 2;++i)
            for(int j = -1;j < 2;++j){
                PHYSBAM_ASSERT(T_MASK::Linear_Offset(std_array<unsigned int,d>(base_index(1) * i * block_xsize,base_index(2) * j * block_ysize))%4096 == 0);
                std_array<int,d> coord(base_index(1) +  i * block_xsize,base_index(2) + j * block_ysize);
                unsigned long block_offset = T_MASK::Linear_Offset(coord);
                if(set.IsPageActive(block_offset))
                    block_base_offsets[i + 1][j + 1] = block_offset;
                else
                    block_base_offsets[i + 1][j + 1] = 0xdeadbeef;
            }
            for(RANGE_ITERATOR<d> range_iterator(RANGE<T_INDEX>(T_INDEX() - 3, block_size + 3));
                range_iterator.Valid();
                range_iterator.Next()){
                T_INDEX index = range_iterator.Index();
                //first determin which block are we pick
                const int block_i = (index.x < 0)?0:((index.x >= block_xsize)?2:1);
                const int block_j = (index.y < 0)?0:((index.y >= block_ysize)?2:1);
                if(block_base_offsets[block_i][block_j] == 0xdeadbeef) continue;
                const int cell_x = (index.x + block_xsize) % block_xsize;
                const int cell_y = (index.y + block_ysize) % block_ysize;
                flag_buffer(index) = reinterpret_cast<block_flag>(flags(block_base_offsets[block_i][block_j]))[cell_x][cell_y];
            }
            //init the distance field
            for(RANGE_ITERATOR<d> range_iterator(RANGE<T_INDEX>(T_INDEX() - 4, block_size + 4));
                range_iterator.Valid();
                range_iterator.Next()){
                T_INDEX index = range_iterator.Index();
                if((flag_buffer(index) == 0x0u) || (flag_buffer(index) & (SPGrid_Solver_Cell_Type_Interface | SPGrid_Solver_Cell_Type_Dirichlet)))
                    distance_buffer(index) = 0;
                else
                    distance_buffer(index) = 0xffff;
            }
            for(int dis = 0;dis < 3;++dis){
                for(RANGE_ITERATOR<d> range_iterator(RANGE<T_INDEX>(T_INDEX() - 3, block_size + 3));
                    range_iterator.Valid();
                    range_iterator.Next()){
                    T_INDEX index = range_iterator.Index();
                    if(distance_buffer(index) != 0){
                        for(int v = 1;v <= d;++v) distance_buffer(index) = std::min<unsigned int>(distance_buffer(index),distance_buffer(index+T_INDEX::Axis_Vector(v))+1);
                        for(int v = 1;v <= d;++v) distance_buffer(index) = std::min<unsigned int>(distance_buffer(index),distance_buffer(index-T_INDEX::Axis_Vector(v))+1);
                    }
                }
            }
            unsigned long flag_offset = iterator.Offset();
            for(RANGE_ITERATOR<d> range_iterator(RANGE<T_INDEX>(T_INDEX(), block_size-1));
                range_iterator.Valid();
                range_iterator.Next(),flag_offset+=sizeof(unsigned)){
                if((flags(flag_offset) & SPGrid_Solver_Cell_Type_Active) && (distance_buffer(range_iterator.Index()) <= 3))
                   flags(flag_offset) |= SPGrid_Solver_Cell_Type_Boundary;
            }
        }
    } 
};
//#####################################################################
}
#endif
