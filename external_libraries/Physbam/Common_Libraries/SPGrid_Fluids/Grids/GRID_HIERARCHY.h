//#####################################################################
// Copyright (c) 2012, Eftychios Sifakis, Sean Bauer
// Distributed under the FreeBSD license (see license.txt)
//#####################################################################
#ifndef __GRID_HIERARCHY_h__
#define __GRID_HIERARCHY_h__

#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <SPGrid/Core/SPGrid_Allocator.h>
#include <SPGrid/Core/SPGrid_Set.h>
#include <SPGrid/Tools/SPGrid_Block_Iterator.h>

namespace PhysBAM{
using namespace SPGrid;
//#####################################################################
// Class GRID_HIERARCHY
//#####################################################################
template<class T_STRUCT,class T,int d>
class GRID_HIERARCHY
{
    typedef SPGrid_Allocator<T_STRUCT,d> Allocator_type;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<unsigned>::type Flag_array_type;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<const unsigned>::type Const_flag_array_type;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<T>::type Data_array_type;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<const T>::type Const_data_array_type;
    typedef SPGrid_Set<Flag_array_type> Set_type;
    typedef std_array<unsigned int,d> coord_t;

    typedef GRID<VECTOR<T,d> > T_GRID;
    typedef VECTOR<int,d> T_INDEX;
    
    const T_GRID& base_grid;
    const int levels;

    ARRAY<T_GRID> grids;
    ARRAY<Allocator_type*> allocators;
    ARRAY<Set_type*> sets;
    ARRAY<std::vector<std::pair<const unsigned long*,unsigned> > > red_blocks;
    ARRAY<std::vector<std::pair<const unsigned long*,unsigned> > > black_blocks;
    ARRAY<std::vector<unsigned long> > copy_of_blocks;
    ARRAY<std::vector<unsigned long> > boundary_blocks;

public:

    const T_GRID& Grid(const int level) const
    {return grids(level);}

    // Access to allocators
    const Allocator_type& Allocator(const int level) const
    {return *allocators(level);}

    Allocator_type& Allocator(const int level)
    {return *allocators(level);}
    
    // Access to sets (more often use blocks)
    const Set_type& Set(const int level) const
    {return *sets(level);}

    Set_type& Set(const int level)
    {return *sets(level);}

    // Access to red blocks
    const std::vector<std::pair<const unsigned long*,unsigned> >& Red_Partition(const int level) const
    {return red_blocks(level);}

    // Access to black blocks
    const std::vector<std::pair<const unsigned long*,unsigned> >& Black_Partition(const int level) const
    {return black_blocks(level);}

    // Access to blocks (for optimized iteration)
    std::pair<const unsigned long*,unsigned> Blocks(const int level) const
    {return sets(level)->Get_Blocks();}

    // Access to boundary blocks (for optimized iteration)
    std::pair<const unsigned long*,unsigned> Boundary_Blocks(const int level) const
    {if((boundary_blocks(level)).size()) return std::pair<const unsigned long*,unsigned>(&((boundary_blocks(level))[0]),(boundary_blocks(level)).size());
    else return std::pair<const unsigned long*,unsigned>((const unsigned long *)0,0);}

    // Access to arrays - also the same as Allocator(level).Get_Array(field)
    Data_array_type Array(const int level,T T_STRUCT::* field)
    {return allocators(level)->Get_Array(field);}
    
    Const_data_array_type Array(const int level,const T T_STRUCT::* field) const
    {return allocators(level)->Get_Const_Array(field);}
    
    Flag_array_type Array(const int level,unsigned T_STRUCT::* field)
    {return allocators(level)->Get_Array(field);}
    
    Const_flag_array_type Array(const int level,const unsigned T_STRUCT::* field) const
    {return allocators(level)->Get_Const_Array(field);}
    
    int Levels() const
    {return levels;}
    
    void Activate_Cell(const int level,const T_INDEX& cell_index,unsigned long mask)
    {sets(level)->Mask(coord_t(cell_index),mask);}
    
    void Activate_Cell(const int level,unsigned long offset,unsigned long mask)
    {sets(level)->Mask(offset,mask);}

// #############################################################################
public:
    GRID_HIERARCHY(const T_GRID& base_grid_input,const int levels_input);    
    ~GRID_HIERARCHY();    
    void Initialize_Grids();
    void Initialize_Allocators();
    void Initialize_Sets();
    void Initialize_Red_Black_Partition(const int number_of_partitions);
    void Initialize_Boundary_Blocks();
    void Update_Block_Offsets();
    void Clear_Bitmaps();
    void Print_Grids(std::ostream& output);
    int Hash(T T_STRUCT::* field);
    void Write_Data_Channel(const std::string& filename,const T T_STRUCT::* field) const;
    void Read_Data_Channel(const std::string& filename,T T_STRUCT::* field);
    void Write_Flags_Channel(const std::string& filename) const;
    void Read_Flags_Channel(const std::string& filename);
    void Write_Block_Offsets(const std::string& filename) const;
    void Read_Block_Offsets(const std::string& filename) const;
// #############################################################################
};
}
#endif
