//#####################################################################
// Copyright (c) 2013, Sean Bauer, Eftychios Sifakis
// Distributed under the FreeBSD license (see license.txt)
//#####################################################################
#ifndef __HIERARCHY_PRECONDITIONER_h__
#define __HIERARCHY_PRECONDITIONER_h__

#include "SPGrid/Core/SPGrid_Utilities.h"
#include "SPGrid/Core/SPGrid_Allocator.h"
#include "SPGrid/Core/SPGrid_Set.h"
#include <SPGrid_Fluids/Projection/Ghost_Value_Accumulate.h>
#include <SPGrid_Fluids/Projection/Ghost_Value_Propagate.h>
#include <SPGrid/Tools/SPGrid_Block_Iterator.h>
#include <SPGrid/Tools/SPGrid_Masked_Set.h>
#include <SPGrid_Fluids/Solvers/HIERARCHY_NEIGHBOR_ITERATOR.h>
#include <SPGrid_Fluids/Solvers/TRAILING_NEIGHBOR_ITERATOR.h>
#include <SPGrid_Fluids/Solvers/TRAILING_NGN_ITERATOR.h>
#include <SPGrid_Fluids/Solvers/Forward_Substitution_Helper.h>
#include <SPGrid_Fluids/Solvers/Backward_Substitution_Helper.h>
#include <SPGrid_Fluids/Grids/GRID_TOPOLOGY_HELPER.h>

#include <SPGrid/Tools/SPGrid_Clear.h>
#include <SPGrid/Tools/SPGrid_Masked_Copy_Or_Clear.h>
#include <SPGrid/Tools/SPGrid_Threading_Helper.h>

//#define TIMING

using namespace SPGrid;

namespace PhysBAM{

extern int PhysBAM_number_of_threads;

template<class T_STRUCT,class T,int d> class HIERARCHY_PRECONDITIONER;

//#####################################################################
// Class HIERARCHY_PRECONDITIONER
//#####################################################################
template<class T_STRUCT, class T>
class HIERARCHY_PRECONDITIONER<T_STRUCT,T,2>
{
    enum {d=2};
    typedef GRID_HIERARCHY<T_STRUCT,T,d> Hierarchy_type;
    typedef SPGrid_Allocator<T_STRUCT,d> Allocator_type;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<unsigned>::type Flag_array_type;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<T>::type Data_array_type;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<const unsigned>::type Const_flag_array_type;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<const T>::type Const_data_array_type;
    typedef SPGrid_Set<Flag_array_type> Set_type;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<unsigned>::mask MASK;
public:
    static void In_Place_Incomplete_Cholesky_Factorization(
        Hierarchy_type& hierarchy,
        T T_STRUCT::* diag_field,
        VECTOR<T T_STRUCT::*,2> L_channels,
        VECTOR<T T_STRUCT::*,2> U_channels)
    {
        In_Place_Incomplete_Cholesky_Factorization(hierarchy,diag_field,
                                                   L_channels(1),L_channels(2),
                                                   U_channels(1),U_channels(2));
    }

    static void In_Place_Incomplete_Cholesky_Factorization(
        Hierarchy_type& hierarchy,
        T T_STRUCT::* diag_field,
        T T_STRUCT::* Lx_field,
        T T_STRUCT::* Ly_field,
        T T_STRUCT::* Ux_field,
        T T_STRUCT::* Uy_field
        )
    {
        int levels = hierarchy.Levels();
        int total_cells = 0;
        int cell_count = 0;
        
        for(int level=1;level<=levels;level++)
        {
            Flag_array_type flags = hierarchy.Array(level,&T_STRUCT::flags); 
            for(SPGrid_Block_Iterator<typename Data_array_type::MASK> iterator(hierarchy.Blocks(level));iterator.Valid();iterator.Next())
                if(iterator.Data(flags) & SPGrid_Cell_Type_Active)
                    total_cells++;
        }
        
        bool modified = false;
        T modified_coefficient = (T).97;
        T zero_tolerance   = (T)1e-8;
        T zero_replacement = (T)1e-8;
        
        for(int level=1;level<=levels;level++)
        {
            Flag_array_type flags = hierarchy.Array(level,&T_STRUCT::flags); 

            for(SPGrid_Block_Iterator<typename Data_array_type::MASK> iterator(hierarchy.Blocks(level));iterator.Valid();iterator.Next())
                if(iterator.Data(flags) & SPGrid_Cell_Type_Active)
                {
                    cell_count++;
                    T sum = 0;
                    CELL_ID i_cid(level,iterator.Offset());

                    HIERARCHY_NEIGHBOR_ITERATOR<T_STRUCT,T,d> itrI(true,
                        hierarchy,
                        i_cid,
                        &T_STRUCT::flags,
                        diag_field, // Use diag_field so we can retrieve via Data() call
                        diag_field,
                        Lx_field,
                        Ly_field,
                        Ux_field,
                        Uy_field
                    );

                    for(;itrI.Valid() && itrI.NCID() < i_cid; itrI.Next())
                    {
                        CELL_ID k_cid = itrI.NCID();

                        TRAILING_NGN_ITERATOR<T_STRUCT,T,d> itrK(true,
                            hierarchy,
                            k_cid,
                            &T_STRUCT::flags,
                            diag_field, // Use diag_field so we can easily retrieve via Data() call
                            diag_field,
                            Ux_field,
                            Uy_field
                        );

                        //itrK.Set_To_Diagonal(); this is already true
                        T& diag_k = itrK.Coefficient();
                        itrI.Coefficient() *= diag_k;

                        HIERARCHY_NEIGHBOR_ITERATOR<T_STRUCT,T,d> itrJ(itrI);
                        itrJ.Next();

                        for(itrK.Next();itrK.Valid();itrK.Next())
                        {
                            T dot_product_term = itrI.Coefficient() * itrK.Coefficient();
                            while(itrJ.Valid() && itrJ.NCID() < itrK.NCID()) itrJ.Next();
                            if(itrJ.Valid() && itrJ.NCID()==itrK.NCID()) itrJ.Coefficient() -= dot_product_term;
                            else if(modified) sum += dot_product_term;
                        }
                    }


                    T& diag_i = itrI.Coefficient();
                    T denominator = diag_i - modified_coefficient*sum;
                    if(cell_count==total_cells && denominator <= zero_tolerance) diag_i = 1/zero_replacement;
                    else diag_i = 1/denominator;
                }
        }
    }

    static void Clear_Cross_Partition_Coefficients(
        Hierarchy_type& hierarchy,
        const ARRAY<int>& partitions,
        VECTOR<T T_STRUCT::*,2> L_channels,
        T T_STRUCT::* temp_channel)
    {
        typedef std::pair<const unsigned long*,unsigned> T_BLOCK;

        const int levels=hierarchy.Levels();
        PHYSBAM_ASSERT(partitions.m==levels);

        for(int level=1;level<=levels;level++)
            if(partitions(level)>1){

                // Create a balanced partitioning of blocks
                const unsigned long* block_offsets=hierarchy.Blocks(level).first;
                const int size=hierarchy.Blocks(level).second;
                int number_of_partitions=partitions(level);
                ARRAY<T_BLOCK> blocks_of_partition;               

                for(int partition=0;partition<number_of_partitions;partition++){
                    int first_index_of_partition=(size/number_of_partitions)*(partition)+std::min(size%number_of_partitions,partition);
                    int last_index_of_partition=(size/number_of_partitions)*(partition+1)+std::min(size%number_of_partitions,partition+1)-1;
                    int block_size=(size/number_of_partitions)+((partition<size%number_of_partitions)?1:0);
                    PHYSBAM_ASSERT(block_size==(last_index_of_partition-first_index_of_partition+1));
                    blocks_of_partition.Append(T_BLOCK(block_offsets+first_index_of_partition,block_size));}

                // Clear temporary channel
                SPGrid_Computations::Clear<T_STRUCT,T,d>(hierarchy.Allocator(level),hierarchy.Blocks(level),temp_channel);

                // Tag active/ghost cells according to their respective partition
                for(int partition=1;partition<=number_of_partitions;partition++)
                    SPGrid_Computations::Masked_Set(
                        hierarchy.Allocator(level),
                        blocks_of_partition(partition),
                        temp_channel,
                        &T_STRUCT::flags,
                        unsigned(SPGrid_Cell_Type_Active|SPGrid_Cell_Type_Ghost),
                        T(partition));

                // Clear active faces on partition interfaces
                Const_flag_array_type flags=hierarchy.Allocator(level).Get_Const_Array(&T_STRUCT::flags);
                Const_data_array_type partitions=hierarchy.Allocator(level).Get_Const_Array(temp_channel);
                for(int axis=1;axis<=d;axis++){                    
                    Data_array_type coefficients=hierarchy.Allocator(level).Get_Array(L_channels(axis));
                    const unsigned long Negative_Axis_Vector_Offset=GRID_TOPOLOGY_HELPER<typename Flag_array_type::MASK>::Negative_Axis_Vector_Offset(axis);
                    const unsigned face_active_mask=GRID_TOPOLOGY_HELPER<typename Flag_array_type::MASK>::Face_Active_Mask(axis);
                    for(SPGrid_Block_Iterator<typename Data_array_type::MASK> iterator(hierarchy.Blocks(level));iterator.Valid();iterator.Next())
                        if(iterator.Data(flags) & face_active_mask)
                            if(iterator.Data(partitions) && iterator.Data(partitions,Negative_Axis_Vector_Offset) && (iterator.Data(partitions) != iterator.Data(partitions,Negative_Axis_Vector_Offset)))
                                iterator.Data(coefficients)=(T)0.;}
            }
    }

    static void Solve_Forward_Substitution(
        Hierarchy_type& hierarchy,
        T T_STRUCT::* b_field, 
        T T_STRUCT::* temp_field,
        T T_STRUCT::* diag_field,
        VECTOR<T T_STRUCT::*,2> L_channels,
        const ARRAY<int>& subdivision_partitions
        )
    {
        Solve_Forward_Substitution(hierarchy,b_field,temp_field,diag_field,L_channels(1),L_channels(2),subdivision_partitions);
    }

    static void Solve_Forward_Substitution(
        Hierarchy_type& hierarchy,
        T T_STRUCT::* b_field, 
        T T_STRUCT::* temp_field,
        T T_STRUCT::* diag_field,
        T T_STRUCT::* Lx_field,
        T T_STRUCT::* Ly_field,
        const ARRAY<int>& subdivision_partitions
        )
    {
#ifdef TIMING
        LOG::SCOPE scope("HIERARCHY_PRECONDITIONER::Solve_Forward_Substitution");
#endif
        int levels = hierarchy.Levels();
        // Copy r into temp, zero ghost values
        for(int level=1; level<=levels;level++)
        {
#ifdef TIMING
            LOG::SCOPE scope("Masked_Copy_Or_Clear");
#endif
            
            if(PhysBAM_number_of_threads)
                SPGrid_Computations::Threading_Helper<T_STRUCT,d>(hierarchy.Allocator(level),hierarchy.Blocks(level)).Run_Parallel(
                    SPGrid_Computations::Masked_Copy_Or_Clear<T_STRUCT,T,unsigned,d>(temp_field,b_field,&T_STRUCT::flags,
                        (unsigned)SPGrid_Cell_Type_Active,(unsigned)SPGrid_Cell_Type_Ghost),
                    PhysBAM_number_of_threads);
            else
                SPGrid_Computations::Masked_Copy_Or_Clear<T_STRUCT,T,unsigned,d>(hierarchy.Allocator(level),hierarchy.Blocks(level),
                    temp_field,b_field,&T_STRUCT::flags,(unsigned)SPGrid_Cell_Type_Active,(unsigned)SPGrid_Cell_Type_Ghost);
        }
        
        // Perform forward sub
        for(int level=1;level<=levels;level++)
        {
#ifdef TIMING
            LOG::SCOPE scope("Forward Sub");
#endif
            Flag_array_type flags = hierarchy.Array(level,&T_STRUCT::flags); 
            Data_array_type data_array = hierarchy.Array(level,temp_field);

#if 0
            for(SPGrid_Block_Iterator<typename Data_array_type::MASK> iterator(hierarchy.Blocks(level));iterator.Valid();iterator.Next())
                if(iterator.Data(flags) & SPGrid_Cell_Type_Active)
                {
                    TRAILING_NEIGHBOR_ITERATOR<T_STRUCT,T,d> nitr(
                        hierarchy,
                        CELL_ID(level,iterator.Offset()),
                        &T_STRUCT::flags,
                        temp_field,
                        diag_field,
                        Lx_field,
                        Ly_field
                    );

                    // Trailing neighbors only include nid's larger than current
                    for(;nitr.Valid();nitr.Next())
                        nitr.Data() -= nitr.Coefficient()*iterator.Data(data_array);
                }
#else
            Forward_Substitution_Helper<T,NextLogTwo<sizeof(T_STRUCT)>::value,d> helper(
                (T*)hierarchy.Array(level,temp_field).Get_Data_Ptr(),
                (T*)hierarchy.Array(level,Lx_field).Get_Data_Ptr(),
                (T*)hierarchy.Array(level,Ly_field).Get_Data_Ptr(),
                (unsigned*)hierarchy.Array(level,&T_STRUCT::flags).Get_Data_Ptr(),
                hierarchy.Blocks(level).first,
                hierarchy.Blocks(level).second);
            if(PhysBAM_number_of_threads)
                helper.Run_Parallel(subdivision_partitions(level));
            else
                helper.Run();
#endif                
                // Accumulate ghosts to higher level (to properly update data entries)
                if(level<levels)
                {
#ifdef SERIAL
                    SPGRID_LAPLACE<T_STRUCT,T,2>::AccumulateGhostValues(
                        hierarchy.Set(level),
                        hierarchy.Allocator(level).Get_Array(temp_field),
                        hierarchy.Allocator(level+1).Get_Array(temp_field)
                        );
#else
                    Ghost_Value_Accumulate<T,T_STRUCT,d>
                        helper( (unsigned*)hierarchy.Set(level+1).array.Get_Data_Ptr(),
                                (T*)hierarchy.Array(level, temp_field).Get_Data_Ptr(),
                                (T*)hierarchy.Array(level+1, temp_field).Get_Data_Ptr(),
                                hierarchy.Blocks(level+1).first,
                                hierarchy.Blocks(level+1).second);
                    if(PhysBAM_number_of_threads)
                        helper.Run_Parallel(PhysBAM_number_of_threads);
                    else
                        helper.Run();
#endif
                }
        }
    }
    
    static void Solve_Backward_Substitution(
        Hierarchy_type& hierarchy,
        T T_STRUCT::* b_field, 
        T T_STRUCT::* x_field,
        T T_STRUCT::* diag_field,
        VECTOR<T T_STRUCT::*,2> L_channels,
        const ARRAY<int>& subdivision_partitions
        )
    {
        Solve_Backward_Substitution(hierarchy,b_field,x_field,diag_field,L_channels(1),L_channels(2),subdivision_partitions);
    }
    
    static void Solve_Backward_Substitution(
        Hierarchy_type& hierarchy,
        T T_STRUCT::* b_field, 
        T T_STRUCT::* x_field,
        T T_STRUCT::* diag_field,
        T T_STRUCT::* Lx_field,
        T T_STRUCT::* Ly_field,
        const ARRAY<int>& subdivision_partitions
        )
    {
#ifdef TIMING
        LOG::SCOPE scope("HIERARCHY_PRECONDITIONER::Solve_Backward_Substitution");
#endif
        for(int level=hierarchy.Levels();level>=1;level--)
        {
#ifdef TIMING
            LOG::SCOPE scope("Backward Sub");
#endif

#if 0
            Flag_array_type flags = hierarchy.Array(level,&T_STRUCT::flags); 
            Data_array_type b_array = hierarchy.Array(level,b_field);
            Data_array_type x_array = hierarchy.Array(level,x_field);
            Data_array_type diag_array = hierarchy.Array(level,diag_field);

            for(SPGrid_Reverse_Block_Iterator<typename Data_array_type::MASK> iterator(hierarchy.Blocks(level));iterator.Valid();iterator.Next())
                if(iterator.Data(flags) & SPGrid_Cell_Type_Active)
                {
                    iterator.Data(x_array) = iterator.Data(b_array)*iterator.Data(diag_array);
                    TRAILING_NEIGHBOR_ITERATOR<T_STRUCT,T,d> nitr(
                        hierarchy,
                        CELL_ID(level,iterator.Offset()),
                        &T_STRUCT::flags,
                        x_field,
                        diag_field,
                        Lx_field,
                        Ly_field
                    );
                    // Trailing neighbors only include nid's larger than current
                    for(;nitr.Valid();nitr.Next())
                        iterator.Data(x_array) -= nitr.Coefficient()*nitr.Data();
                }
#else
            Backward_Substitution_Helper<T,NextLogTwo<sizeof(T_STRUCT)>::value,d> helper(
                (T*)hierarchy.Array(level,x_field).Get_Data_Ptr(),
                (T*)hierarchy.Array(level,diag_field).Get_Data_Ptr(),
                (T*)hierarchy.Array(level,Lx_field).Get_Data_Ptr(),
                (T*)hierarchy.Array(level,Ly_field).Get_Data_Ptr(),
                (unsigned*)hierarchy.Array(level,&T_STRUCT::flags).Get_Data_Ptr(),
                hierarchy.Blocks(level).first,
                hierarchy.Blocks(level).second);
            if(PhysBAM_number_of_threads)
                helper.Run_Parallel(subdivision_partitions(level));
            else
                helper.Run();
#endif
                
            if(level>1)
            {
                Ghost_Value_Propagate<T,T_STRUCT,d>
                    helper( (unsigned*)hierarchy.Set(level-1).array.Get_Data_Ptr(),
                            (T*)hierarchy.Array(level-1, x_field).Get_Data_Ptr(),
                            (T*)hierarchy.Array(level, x_field).Get_Data_Ptr(),
                            hierarchy.Blocks(level-1).first,
                            hierarchy.Blocks(level-1).second);
                if(PhysBAM_number_of_threads)
                    helper.Run_Parallel(PhysBAM_number_of_threads);
                else
                    helper.Run();
            }
        }
    }

};


template<class T_STRUCT, class T>
class HIERARCHY_PRECONDITIONER<T_STRUCT,T,3>
{
    enum {d=3};
    typedef GRID_HIERARCHY<T_STRUCT,T,d> Hierarchy_type;
    typedef SPGrid_Allocator<T_STRUCT,d> Allocator_type;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<unsigned>::type Flag_array_type;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<T>::type Data_array_type;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<const unsigned>::type Const_flag_array_type;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<const T>::type Const_data_array_type;
    typedef SPGrid_Set<Flag_array_type> Set_type;
    
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<unsigned>::mask MASK;
public:
    static void In_Place_Incomplete_Cholesky_Factorization(
        Hierarchy_type& hierarchy,
        T T_STRUCT::* diag_field,
        VECTOR<T T_STRUCT::*,3> L_channels,
        VECTOR<T T_STRUCT::*,3> U_channels)
    {
        In_Place_Incomplete_Cholesky_Factorization(hierarchy,diag_field,
                                                   L_channels(1),L_channels(2),L_channels(3),
                                                   U_channels(1),U_channels(2),U_channels(3));
    }

    static void In_Place_Incomplete_Cholesky_Factorization(
        Hierarchy_type& hierarchy,
        T T_STRUCT::* diag_field,
        T T_STRUCT::* Lx_field,
        T T_STRUCT::* Ly_field,
        T T_STRUCT::* Lz_field,
        T T_STRUCT::* Ux_field,
        T T_STRUCT::* Uy_field,
        T T_STRUCT::* Uz_field
        )
    {
        int levels = hierarchy.Levels();
        int total_cells = 0;
        int cell_count = 0;
        
        for(int level=1;level<=levels;level++)
        {
            Flag_array_type flags = hierarchy.Array(level,&T_STRUCT::flags); 
            for(SPGrid_Block_Iterator<typename Data_array_type::MASK> iterator(hierarchy.Blocks(level));iterator.Valid();iterator.Next())
                if(iterator.Data(flags) & SPGrid_Cell_Type_Active)
                    total_cells++;
        }
        
        bool modified = false;
        T modified_coefficient = (T).97;
        T zero_tolerance   = (T)1e-8;
        T zero_replacement = (T)1e-8;
        
        for(int level=1;level<=levels;level++)
        {
            Flag_array_type flags = hierarchy.Array(level,&T_STRUCT::flags); 

            for(SPGrid_Block_Iterator<typename Data_array_type::MASK> iterator(hierarchy.Blocks(level));iterator.Valid();iterator.Next())
                if(iterator.Data(flags) & SPGrid_Cell_Type_Active)
                {
                    cell_count++;
                    T sum = 0;
                    CELL_ID i_cid(level,iterator.Offset());
                    
                    HIERARCHY_NEIGHBOR_ITERATOR<T_STRUCT,T,d> itrI(true,
                        hierarchy,
                        i_cid,
                        &T_STRUCT::flags,
                        diag_field, // Use diag_field so we can retrieve via Data() call
                        diag_field,
                        Lx_field,
                        Ly_field,
                        Lz_field,
                        Ux_field,
                        Uy_field,
                        Uz_field
                    );

                    for(;itrI.Valid() && itrI.NCID() < i_cid; itrI.Next())
                    {
                        CELL_ID k_cid = itrI.NCID();

                        TRAILING_NGN_ITERATOR<T_STRUCT,T,d> itrK(true,
                            hierarchy,
                            k_cid,
                            &T_STRUCT::flags,
                            diag_field, // Use diag_field so we can easily retrieve via Data() call
                            diag_field,
                            Ux_field,
                            Uy_field,
                            Uz_field
                        );

                        //itrK.Set_To_Diagonal(); this is already true
                        T& diag_k = itrK.Coefficient();
                        itrI.Coefficient() *= diag_k;

                        HIERARCHY_NEIGHBOR_ITERATOR<T_STRUCT,T,d> itrJ(itrI);
                        itrJ.Next();

                        for(itrK.Next();itrK.Valid();itrK.Next())
                        {
                            T dot_product_term = itrI.Coefficient() * itrK.Coefficient();
                            while(itrJ.Valid() && itrJ.NCID() < itrK.NCID())
                                itrJ.Next();
                            if(itrJ.Valid() && itrJ.NCID()==itrK.NCID())
                                itrJ.Coefficient() -= dot_product_term;
                            else if(modified) sum += dot_product_term;
                        }
                    }


                    T& diag_i = itrI.Coefficient();
                    T denominator = diag_i - modified_coefficient*sum;
                    if(cell_count==total_cells && denominator <= zero_tolerance) diag_i = 1/zero_replacement;
                    else diag_i = 1/denominator;
                }
        }
    }

    static void Clear_Cross_Partition_Coefficients(
        Hierarchy_type& hierarchy,
        const ARRAY<int>& partitions,
        VECTOR<T T_STRUCT::*,3> L_channels,
        T T_STRUCT::* temp_channel)
    {
        typedef std::pair<const unsigned long*,unsigned> T_BLOCK;

        const int levels=hierarchy.Levels();
        PHYSBAM_ASSERT(partitions.m==levels);

        for(int level=1;level<=levels;level++)
            if(partitions(level)>1){

                // Create a balanced partitioning of blocks
                const unsigned long* block_offsets=hierarchy.Blocks(level).first;
                const int size=hierarchy.Blocks(level).second;
                int number_of_partitions=partitions(level);
                ARRAY<T_BLOCK> blocks_of_partition;

                for(int partition=0;partition<number_of_partitions;partition++){
                    int first_index_of_partition=(size/number_of_partitions)*(partition)+std::min(size%number_of_partitions,partition);
                    int last_index_of_partition=(size/number_of_partitions)*(partition+1)+std::min(size%number_of_partitions,partition+1)-1;
                    int block_size=(size/number_of_partitions)+((partition<size%number_of_partitions)?1:0);
                    PHYSBAM_ASSERT(block_size==(last_index_of_partition-first_index_of_partition+1));
                    blocks_of_partition.Append(T_BLOCK(block_offsets+first_index_of_partition,block_size));}

                // Clear temporary channel
                SPGrid_Computations::Clear<T_STRUCT,T,d>(hierarchy.Allocator(level),hierarchy.Blocks(level),temp_channel);

                // Tag active/ghost cells according to their respective partition
                for(int partition=1;partition<=number_of_partitions;partition++)
                    SPGrid_Computations::Masked_Set(
                        hierarchy.Allocator(level),
                        blocks_of_partition(partition),
                        temp_channel,
                        &T_STRUCT::flags,
                        unsigned(SPGrid_Cell_Type_Active|SPGrid_Cell_Type_Ghost),
                        T(partition));

                // Clear active faces on partition interfaces
                Const_flag_array_type flags=hierarchy.Allocator(level).Get_Const_Array(&T_STRUCT::flags);
                Const_data_array_type partitions=hierarchy.Allocator(level).Get_Const_Array(temp_channel);
                for(int axis=1;axis<=d;axis++){                    
                    Data_array_type coefficients=hierarchy.Allocator(level).Get_Array(L_channels(axis));
                    const unsigned long Negative_Axis_Vector_Offset=GRID_TOPOLOGY_HELPER<typename Flag_array_type::MASK>::Negative_Axis_Vector_Offset(axis);
                    const unsigned face_active_mask=GRID_TOPOLOGY_HELPER<typename Flag_array_type::MASK>::Face_Active_Mask(axis);
                    for(SPGrid_Block_Iterator<typename Data_array_type::MASK> iterator(hierarchy.Blocks(level));iterator.Valid();iterator.Next())
                        if(iterator.Data(flags) & face_active_mask)
                            if(iterator.Data(partitions) && iterator.Data(partitions,Negative_Axis_Vector_Offset) && (iterator.Data(partitions) != iterator.Data(partitions,Negative_Axis_Vector_Offset)))
                                iterator.Data(coefficients)=(T)0.;}
            }
    }

    static void Solve_Forward_Substitution(
        Hierarchy_type& hierarchy,
        T T_STRUCT::* b_field, 
        T T_STRUCT::* temp_field,
        T T_STRUCT::* diag_field,
        VECTOR<T T_STRUCT::*,3> L_channels,
        const ARRAY<int>& subdivision_partitions
        )
    {
        Solve_Forward_Substitution(hierarchy,b_field,temp_field,diag_field,L_channels(1),L_channels(2),L_channels(3),subdivision_partitions);
    }

    static void Solve_Forward_Substitution(
        Hierarchy_type& hierarchy,
        T T_STRUCT::* b_field, 
        T T_STRUCT::* temp_field,
        T T_STRUCT::* diag_field,
        T T_STRUCT::* Lx_field,
        T T_STRUCT::* Ly_field,
        T T_STRUCT::* Lz_field,
        const ARRAY<int>& subdivision_partitions
        )
    {
#ifdef TIMING
        LOG::SCOPE scope("HIERARCHY_PRECONDITIONER::Solve_Forward_Substitution");
#endif
        int levels = hierarchy.Levels();
        // Copy r into temp, zero ghost values
        for(int level=1; level<=levels;level++)
        {
#ifdef TIMING
            LOG::SCOPE scope("Masked_Copy_Or_Clear");
#endif     
            if(PhysBAM_number_of_threads)
                SPGrid_Computations::Threading_Helper<T_STRUCT,d>(hierarchy.Allocator(level),hierarchy.Blocks(level)).Run_Parallel(
                    SPGrid_Computations::Masked_Copy_Or_Clear<T_STRUCT,T,unsigned,d>(temp_field,b_field,&T_STRUCT::flags,
                        (unsigned)SPGrid_Cell_Type_Active,(unsigned)SPGrid_Cell_Type_Ghost),
                    PhysBAM_number_of_threads);
            else
                SPGrid_Computations::Masked_Copy_Or_Clear<T_STRUCT,T,unsigned,d>(hierarchy.Allocator(level),hierarchy.Blocks(level),
                    temp_field,b_field,&T_STRUCT::flags,(unsigned)SPGrid_Cell_Type_Active,(unsigned)SPGrid_Cell_Type_Ghost);
        }
        
        // Perform forward sub
        for(int level=1;level<=levels;level++)
        {
#ifdef TIMING
            LOG::SCOPE scope("Forward Sub");
#endif
            Flag_array_type flags = hierarchy.Array(level,&T_STRUCT::flags); 
            Data_array_type data_array = hierarchy.Array(level,temp_field);
#if 0
            for(SPGrid_Block_Iterator<typename Data_array_type::MASK> iterator(hierarchy.Blocks(level));iterator.Valid();iterator.Next())
                if(iterator.Data(flags) & SPGrid_Cell_Type_Active)
                {
                    TRAILING_NEIGHBOR_ITERATOR<T_STRUCT,T,d> nitr(
                        hierarchy,
                        CELL_ID(level,iterator.Offset()),
                        &T_STRUCT::flags,
                        temp_field,
                        diag_field,
                        Lx_field,
                        Ly_field,
                        Lz_field
                    );

                    // Trailing neighbors only include nid's larger than current
                    for(;nitr.Valid();nitr.Next())
                        nitr.Data() -= nitr.Coefficient()*iterator.Data(data_array);
                }
#else
            Forward_Substitution_Helper<T,NextLogTwo<sizeof(T_STRUCT)>::value,d> helper(
                (T*)hierarchy.Array(level,temp_field).Get_Data_Ptr(),
                (T*)hierarchy.Array(level,Lx_field).Get_Data_Ptr(),
                (T*)hierarchy.Array(level,Ly_field).Get_Data_Ptr(),
                (T*)hierarchy.Array(level,Lz_field).Get_Data_Ptr(),
                (unsigned*)hierarchy.Array(level,&T_STRUCT::flags).Get_Data_Ptr(),
                hierarchy.Blocks(level).first,
                hierarchy.Blocks(level).second);
            if(PhysBAM_number_of_threads)
                helper.Run_Parallel(subdivision_partitions(level));
            else
                helper.Run();
#endif
                // Accumulate ghosts to higher level (to properly update data entries)
                if(level<levels)
                {
#ifdef SERIAL
                    SPGRID_LAPLACE<T_STRUCT,T,3>::AccumulateGhostValues(
                        hierarchy.Set(level),
                        hierarchy.Allocator(level).Get_Array(temp_field),
                        hierarchy.Allocator(level+1).Get_Array(temp_field)
                        );
#else
                    Ghost_Value_Accumulate<T,T_STRUCT,d>
                        helper( (unsigned*)hierarchy.Set(level+1).array.Get_Data_Ptr(),
                                (T*)hierarchy.Array(level, temp_field).Get_Data_Ptr(),
                                (T*)hierarchy.Array(level+1, temp_field).Get_Data_Ptr(),
                                hierarchy.Blocks(level+1).first,
                                hierarchy.Blocks(level+1).second);
                    if(PhysBAM_number_of_threads)
                        helper.Run_Parallel(PhysBAM_number_of_threads);
                    else
                        helper.Run();
#endif
                }
        }
    }

    static void Solve_Backward_Substitution(
        Hierarchy_type& hierarchy,
        T T_STRUCT::* b_field, 
        T T_STRUCT::* x_field,
        T T_STRUCT::* diag_field,
        VECTOR<T T_STRUCT::*,3> L_channels,
        const ARRAY<int>& subdivision_partitions
        )
    {
        Solve_Backward_Substitution(hierarchy,b_field,x_field,diag_field,L_channels(1),L_channels(2),L_channels(3),subdivision_partitions);
    }
    
    static void Solve_Backward_Substitution(
        Hierarchy_type& hierarchy,
        T T_STRUCT::* b_field, 
        T T_STRUCT::* x_field,
        T T_STRUCT::* diag_field,
        T T_STRUCT::* Lx_field,
        T T_STRUCT::* Ly_field,
        T T_STRUCT::* Lz_field,
        const ARRAY<int>& subdivision_partitions
        )
    {
#ifdef TIMING
        LOG::SCOPE scope("HIERARCHY_PRECONDITIONER::Solve_Backward_Substitution");
#endif
        int count = 0;
        for(int level=hierarchy.Levels();level>=1;level--)
        {
#ifdef TIMING
            LOG::SCOPE scope("Backward Sub");
#endif

#if 0
            Flag_array_type flags = hierarchy.Array(level,&T_STRUCT::flags); 
            Data_array_type b_array = hierarchy.Array(level,b_field);
            Data_array_type x_array = hierarchy.Array(level,x_field);
            Data_array_type diag_array = hierarchy.Array(level,diag_field);

            for(SPGrid_Reverse_Block_Iterator<typename Data_array_type::MASK> iterator(hierarchy.Blocks(level));iterator.Valid();iterator.Next())
                if(iterator.Data(flags) & SPGrid_Cell_Type_Active)
                {
                    iterator.Data(x_array) = iterator.Data(b_array)*iterator.Data(diag_array);
                    TRAILING_NEIGHBOR_ITERATOR<T_STRUCT,T,d> nitr(
                        hierarchy,
                        CELL_ID(level,iterator.Offset()),
                        &T_STRUCT::flags,
                        x_field,
                        diag_field,
                        Lx_field,
                        Ly_field,
                        Lz_field
                    );
                    // Trailing neighbors only include nid's larger than current
                    for(;nitr.Valid();nitr.Next())
                        iterator.Data(x_array) -= nitr.Coefficient()*nitr.Data();
                }
#else
            Backward_Substitution_Helper<T,NextLogTwo<sizeof(T_STRUCT)>::value,d> helper(
                (T*)hierarchy.Array(level,x_field).Get_Data_Ptr(),
                (T*)hierarchy.Array(level,diag_field).Get_Data_Ptr(),
                (T*)hierarchy.Array(level,Lx_field).Get_Data_Ptr(),
                (T*)hierarchy.Array(level,Ly_field).Get_Data_Ptr(),
                (T*)hierarchy.Array(level,Lz_field).Get_Data_Ptr(),
                (unsigned*)hierarchy.Array(level,&T_STRUCT::flags).Get_Data_Ptr(),
                hierarchy.Blocks(level).first,
                hierarchy.Blocks(level).second);
            if(PhysBAM_number_of_threads)
                helper.Run_Parallel(subdivision_partitions(level));
            else
                helper.Run();
#endif       
            if(level>1)
            {
                Ghost_Value_Propagate<T,T_STRUCT,d>
                    helper( (unsigned*)hierarchy.Set(level-1).array.Get_Data_Ptr(),
                            (T*)hierarchy.Array(level-1, x_field).Get_Data_Ptr(),
                            (T*)hierarchy.Array(level, x_field).Get_Data_Ptr(),
                            hierarchy.Blocks(level-1).first,
                            hierarchy.Blocks(level-1).second);
                if(PhysBAM_number_of_threads)
                    helper.Run_Parallel(PhysBAM_number_of_threads);
                else
                    helper.Run();
            }
        }
    }

};
}
#endif
