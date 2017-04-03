#ifndef __SPGRID_SUBDOMAIN_RASTERIZER_H__
#define __SPGRID_SUBDOMAIN_RASTERIZER_H__
#include <SPGrid/Core/SPGrid_Allocator.h>
#include <SPGrid/Core/SPGrid_Set.h>
#include <SPGrid/Tools/SPGrid_Block_Iterator.h>
namespace SPGrid{
template<class T_STRUCT, class T,int d>
class Subdomain_Rasterizer
{
    typedef SPGrid_Allocator<T_STRUCT,d> SPG_Allocator;
    typedef typename SPG_Allocator::Array<unsigned>::type SPG_Flags_Array_Type;
    typedef typename SPG_Allocator::Array<const unsigned>::type SPG_Const_Flags_Array_Type;
    typedef SPGrid_Set<SPG_Flags_Array_Type> SPG_Set_Type;

    typedef VECTOR<int,d> T_INDEX;
    typedef VECTOR<T,d> TV;

    SPG_Set_Type& set;
    SPG_Set_Type& reference_set;
    SPG_Flags_Array_Type& flags;
    SPG_Flags_Array_Type& reference_flags;
    T_INDEX block_size;
    const T_INDEX subdomain_origin;
    const T_INDEX subdomain_size;
    const RANGE<T_INDEX> subdomain;

public:
    Subdomain_Rasterizer(SPG_Set_Type& reference_set_input,SPG_Set_Type& set_input,const T_INDEX& subdomain_origin_input,const T_INDEX& subdomain_size_input)
        :reference_set(reference_set_input),reference_flags(reference_set.array),
         set(set_input),flags(set.array),
         subdomain_origin(subdomain_origin_input),subdomain_size(subdomain_size_input),
         subdomain(subdomain_origin,subdomain_origin+subdomain_size)
    {
        block_size=flags.geometry.Block_Size().template Cast<T_INDEX>();
    }

    bool Consume(const RANGE<VECTOR<int,d> >& range){
        // Make sure we did not descend to a sub-block size
        PHYSBAM_ASSERT(range.Edge_Lengths().All_Greater_Equal(block_size));

        RANGE<T_INDEX> adjusted_range(range);adjusted_range.max_corner-=1; // Adjusted range is inclusive of max corner
        // If the block is fully exterior, do nothing
        if(RANGE<T_INDEX>::Intersect(adjusted_range,subdomain).Empty()) return false;

        if(range.Edge_Lengths()!=block_size) return true;
        // Recurse until we are down to the size of a single block
        unsigned long offset=SPG_Flags_Array_Type::MASK::Linear_Offset(std_array<int,d>(adjusted_range.min_corner));
        if(!reference_set.IsPageActive(offset)) return false;
        //std::cout<<"here?"<<std::endl;
        for(RANGE_ITERATOR<d> iterator(adjusted_range);iterator.Valid();iterator.Next(),offset+=sizeof(unsigned)){
            //set.Mask(offset,SPGrid_Solver_Cell_Type_Active);
            if(subdomain.Lazy_Inside(iterator.Index()))
                set.Mask(offset,reference_flags(offset));
        }
        return false;
    }
};

template<class T_STRUCT, class T,int d>
class Subdomain_Masked_Rasterizer
{
    //THE WATER WILL BE phi<0
    typedef SPGrid_Allocator<T_STRUCT,d> SPG_Allocator;
    typedef typename SPG_Allocator::Array<unsigned>::type SPG_Flags_Array_Type;
    typedef SPGrid_Set<SPG_Flags_Array_Type> SPG_Set_Type;

    typedef VECTOR<int,d> T_INDEX;
    typedef VECTOR<T,d> TV;

    SPG_Set_Type& set;
    SPG_Set_Type& reference_set;
    SPG_Flags_Array_Type& flags;
    SPG_Flags_Array_Type& reference_flags;
    T_INDEX block_size;
    const T_INDEX subdomain_origin;
    const T_INDEX subdomain_size;
    const RANGE<T_INDEX> subdomain;
    const unsigned active_flag;
    const int subdomain_id;
public:
    Subdomain_Masked_Rasterizer(SPG_Set_Type& reference_set_input,SPG_Set_Type& set_input,const T_INDEX& subdomain_origin_input,const T_INDEX& subdomain_size_input,unsigned active_flag_input,int subdomain_id_input)
        :reference_set(reference_set_input),reference_flags(reference_set.array),
         set(set_input),flags(set.array),
         subdomain_origin(subdomain_origin_input),subdomain_size(subdomain_size_input),
         subdomain(subdomain_origin,subdomain_origin+subdomain_size),
         active_flag(active_flag_input),
        subdomain_id(subdomain_id_input)
    {
        block_size=flags.geometry.Block_Size().template Cast<T_INDEX>();
    }

    bool Consume(const RANGE<VECTOR<int,d> >& range){
        // Make sure we did not descend to a sub-block size
        PHYSBAM_ASSERT(range.Edge_Lengths().All_Greater_Equal(block_size));

        // If the block is fully exterior, do nothing
        RANGE<T_INDEX> adjusted_range(range);adjusted_range.max_corner-=1; // Adjusted range is inclusive of max corner
        if(RANGE<T_INDEX>::Intersect(adjusted_range,subdomain).Empty()) return false;

        // Recurse until we are down to the size of a single block
        if(range.Edge_Lengths()!=block_size) return true;

        unsigned long base_offset=SPG_Flags_Array_Type::MASK::Linear_Offset(std_array<int,d>(adjusted_range.min_corner));
        if(!reference_set.IsPageActive(base_offset)) return false;

        RANGE<T_INDEX> expanded_range=adjusted_range.Thickened(1);
        ARRAY<unsigned,T_INDEX> flags_local(expanded_range);
        flags_local.Fill(0x0);
        for(RANGE_ITERATOR<d> iterator(expanded_range);iterator.Valid();iterator.Next()){
            T_INDEX cell_index=iterator.Index();
            unsigned long cell_offset=SPG_Flags_Array_Type::MASK::Linear_Offset(std_array<int,d>(cell_index));
            if(subdomain.Lazy_Inside(cell_index)&&reference_set.IsPageActive(cell_offset)){
                flags_local(cell_index)=reference_flags(cell_offset);
                if(flags_local(cell_index)&active_flag) {flags_local(cell_index)=SPGrid_Solver_Cell_Type_Active/*|subdomain_id*/;}
                else if(flags_local(cell_index)&SPGrid_Solver_Cell_Type_Active){flags_local(cell_index)=0x0;}}}   

        bool active_or_dirichlet_found=false;
        for(RANGE_ITERATOR<d> iterator(adjusted_range);iterator.Valid();iterator.Next()){
            T_INDEX cell_index=iterator.Index();
            if(flags_local(cell_index) & SPGrid_Solver_Cell_Type_Active) active_or_dirichlet_found=true;
            if ((flags_local(cell_index) & SPGrid_Solver_Cell_Type_Dirichlet)||(flags_local(cell_index))){
                bool active_neighbor_found=false;
                for(int v=1;v<=d;v++)
                    if(flags_local(cell_index+T_INDEX::Axis_Vector(v)) & SPGrid_Solver_Cell_Type_Active ||
                       flags_local(cell_index-T_INDEX::Axis_Vector(v)) & SPGrid_Solver_Cell_Type_Active)
                        active_neighbor_found=true;
                if(active_neighbor_found) active_or_dirichlet_found=true;
                else flags_local(cell_index) = 0x0u;
            }
        }

        if(active_or_dirichlet_found){
            set.MarkPageActive(SPG_Flags_Array_Type::MASK::Linear_Offset(std_array<int,d>(range.min_corner)));
            unsigned* flags_ptr=&flags(std_array<int,d>(adjusted_range.min_corner));
            for(RANGE_ITERATOR<d> iterator(adjusted_range);iterator.Valid();iterator.Next(),flags_ptr++)
                *flags_ptr=flags_local(iterator.Index());
        }
        return false;
    }
};
}
#endif
