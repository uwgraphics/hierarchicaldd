//#####################################################################
// Copyright 2013, Raj Setaluri, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <SPGrid_Fluids/Solvers/BUILD_MATRIX.h>
#include <SPGrid_Fluids/Simulation/FLUIDS_SIMULATION_DATA.h>
#include <SPGrid_Fluids/Simulation/FLUIDS_SIMULATION_FLAGS.h>
#include <SPGrid_Fluids/Projection/GRID_HIERARCHY_PROJECTION.h>
#include <SPGrid_Fluids/Grids/GRID_TOPOLOGY_HELPER.h>
#include <SPGrid/Tools/SPGrid_Clear.h>
#include <SPGrid/Tools/SPGrid_Copy.h>
#include <PhysBAM_Tools/Random_Numbers/RANDOM_NUMBERS.h>
#include <SPGrid_Fluids/Solvers/HIERARCHY_NEIGHBOR_ITERATOR.h>
#include <SPGrid_Fluids/Solvers/HIERARCHY_PRECONDITIONER.h>
using namespace PhysBAM;
//#####################################################################
// Build_Matrix
//#####################################################################
template<class T_STRUCT,class T,int d> void BUILD_MATRIX<T_STRUCT,T,d>::
Build_Matrix(Hierarchy_type& hierarchy,T T_STRUCT::* diag_channel,VECTOR<T T_STRUCT::*,d> L_channels,VECTOR<T T_STRUCT::*,d> U_channels,unsigned T_STRUCT::* flags_channel)
{
    const int levels=hierarchy.Levels();
    
    float laplace_scale_uniform[levels];
    float laplace_scale_nonuniform[levels];
    
    VECTOR<unsigned long,d> negative_axis_offsets;
    for(int v=1;v<=d;v++)
        negative_axis_offsets(v)=GRID_TOPOLOGY_HELPER<typename Flag_array_type::MASK>::Negative_Axis_Vector_Offset(v);

    // preparation
    for(int level=1;level<=levels;level++){
        SPGrid_Computations::Clear<T_STRUCT,T,d>(hierarchy.Allocator(level),hierarchy.Blocks(level),diag_channel);
        for(int v=1;v<=d;v++)
            SPGrid_Computations::Clear<T_STRUCT,T,d,2>(hierarchy.Allocator(level),hierarchy.Blocks(level),U_channels(v),L_channels(v));
        laplace_scale_uniform[level-1]=GRID_HIERARCHY_PROJECTION<T_STRUCT,T,d>::Laplace_Scale_Uniform(hierarchy,level);
        laplace_scale_nonuniform[level-1]=GRID_HIERARCHY_PROJECTION<T_STRUCT,T,d>::Laplace_Scale_Nonuniform(hierarchy,level);}

    // actual iterations -- only build L
    for(int level=1;level<=levels;level++){
        Const_flag_array_type flags=hierarchy.Allocator(level).Get_Const_Array(flags_channel);
        Data_array_type diag=hierarchy.Allocator(level).Get_Array(diag_channel);
        VECTOR<void*,d> L_data_ptrs,U_data_ptrs;
        for(int v=1;v<=d;v++){
            L_data_ptrs(v)=hierarchy.Array(level,L_channels(v)).Get_Data_Ptr();
            U_data_ptrs(v)=hierarchy.Array(level,U_channels(v)).Get_Data_Ptr();}

        for(SPGrid_Block_Iterator<typename Data_array_type::MASK> iterator(hierarchy.Blocks(level));iterator.Valid();iterator.Next()){
            // TODO: check these flags
            const unsigned flag=iterator.Data(flags);
            for(int axis=1;axis<=d;axis++){
                if(flag & GRID_TOPOLOGY_HELPER<typename Flag_array_type::MASK>::Face_Minus_Scaled_Mask(axis)){
                    iterator.template Data<Data_array_type>(L_data_ptrs(axis))=laplace_scale_nonuniform[level-1];
                    iterator.Data(diag) += (T)-1.*laplace_scale_nonuniform[level-1];
                    iterator.Data(diag,negative_axis_offsets(axis)) += (T)-1.*laplace_scale_nonuniform[level-1];
                } else if(flag & GRID_TOPOLOGY_HELPER<typename Flag_array_type::MASK>::Face_Active_Mask(axis)) {
                    iterator.template Data<Data_array_type>(L_data_ptrs(axis))=laplace_scale_uniform[level-1];
                    iterator.Data(diag) += (T)-1.*laplace_scale_uniform[level-1];
                    iterator.Data(diag,negative_axis_offsets(axis)) += (T)-1.*laplace_scale_uniform[level-1];
                }
            }
        }
    }
    
    // Accumulate diagonal ghosts
    GRID_HIERARCHY_PROJECTION<T_STRUCT,T,d>::Accumulate_Ghost_Values(hierarchy,flags_channel,diag_channel);

    // Clear out dirichlet entries
    for(int level=1;level<=levels;level++){
        Const_flag_array_type flags=hierarchy.Allocator(level).Get_Const_Array(flags_channel);
        VECTOR<void*,d> L_data_ptrs,U_data_ptrs;
        for(int v=1;v<=d;v++){
            L_data_ptrs(v)=hierarchy.Array(level,L_channels(v)).Get_Data_Ptr();
            U_data_ptrs(v)=hierarchy.Array(level,U_channels(v)).Get_Data_Ptr();}
        for(SPGrid_Block_Iterator<typename Data_array_type::MASK> iterator(hierarchy.Blocks(level));iterator.Valid();iterator.Next()){
            if(iterator.Data(flags) & SPGrid_Cell_Type_Interior){
                HIERARCHY_NEIGHBOR_ITERATOR<T_STRUCT,T,d> neighbor_iterator(false,hierarchy,CELL_ID(level,iterator.Offset()),flags_channel,diag_channel,diag_channel,L_channels,L_channels); // both as L channels
                for(;neighbor_iterator.Valid();neighbor_iterator.Next()){
                    if(((hierarchy.Array(neighbor_iterator.Neighbor().level,flags_channel)(neighbor_iterator.Neighbor().offset)) & SPGrid_Cell_Type_Dirichlet))
                        neighbor_iterator.Coefficient()=(T)0.;
                }
            }
        }
    }

    // Copy L to U
    for(int level=1;level<=levels;level++)
        for(int v=1;v<=d;v++)
            SPGrid_Computations::Copy<T_STRUCT,T,d>(hierarchy.Allocator(level),hierarchy.Blocks(level),L_channels(v),U_channels(v));
 
    // Clear ghosts and Dirichlet on diagonal
    for(int level=1;level<=levels;level++)
        SPGrid_Computations::Masked_Clear(hierarchy.Allocator(level),hierarchy.Blocks(level),diag_channel,flags_channel,(unsigned)(SPGrid_Cell_Type_Ghost|SPGrid_Cell_Type_Dirichlet));
}
//#####################################################################
// Build_Variable_Beta_Matrix
//#####################################################################
template<class T_STRUCT,class T,int d> void BUILD_MATRIX<T_STRUCT,T,d>::
Build_Variable_Beta_Matrix(Hierarchy_type& hierarchy,T T_STRUCT::* diag_channel,VECTOR<T T_STRUCT::*,d> L_channels,VECTOR<T T_STRUCT::*,d> U_channels,T T_STRUCT::* rho_channel,unsigned T_STRUCT::* flags_channel)
{
    const int levels=hierarchy.Levels();
    
    float laplace_scale_uniform[levels];
    float laplace_scale_nonuniform[levels];
    
    VECTOR<unsigned long,d> negative_axis_offsets;
    for(int v=1;v<=d;v++)
        negative_axis_offsets(v)=GRID_TOPOLOGY_HELPER<typename Flag_array_type::MASK>::Negative_Axis_Vector_Offset(v);

    // preparation
    for(int level=1;level<=levels;level++){
        SPGrid_Computations::Clear<T_STRUCT,T,d>(hierarchy.Allocator(level),hierarchy.Blocks(level),diag_channel);
        for(int v=1;v<=d;v++)
            SPGrid_Computations::Clear<T_STRUCT,T,d,2>(hierarchy.Allocator(level),hierarchy.Blocks(level),U_channels(v),L_channels(v));
        laplace_scale_uniform[level-1]=GRID_HIERARCHY_PROJECTION<T_STRUCT,T,d>::Laplace_Scale_Uniform(hierarchy,level);
        laplace_scale_nonuniform[level-1]=GRID_HIERARCHY_PROJECTION<T_STRUCT,T,d>::Laplace_Scale_Nonuniform(hierarchy,level);}

    // TODO: populate rho in ghost cells

    // NOTE: currently will only work for uniform grids!!

    // actual iterations -- only build L
    for(int level=1;level<=levels;level++){
        Const_flag_array_type flags=hierarchy.Allocator(level).Get_Const_Array(flags_channel);
        Data_array_type diag=hierarchy.Allocator(level).Get_Array(diag_channel);
        Const_data_array_type rho=hierarchy.Allocator(level).Get_Const_Array(rho_channel);
        VECTOR<void*,d> L_data_ptrs,U_data_ptrs;
        for(int v=1;v<=d;v++){
            L_data_ptrs(v)=hierarchy.Array(level,L_channels(v)).Get_Data_Ptr();
            U_data_ptrs(v)=hierarchy.Array(level,U_channels(v)).Get_Data_Ptr();}

        for(SPGrid_Block_Iterator<typename Data_array_type::MASK> iterator(hierarchy.Blocks(level));iterator.Valid();iterator.Next()){
            // TODO: check these flags
            const unsigned flag=iterator.Data(flags);
            for(int axis=1;axis<=d;axis++){
                if(flag & GRID_TOPOLOGY_HELPER<typename Flag_array_type::MASK>::Face_Minus_Scaled_Mask(axis)){ // NOTE: assuming rho in ghost cells is correct!
                    const double rho_right=(double)iterator.Data(rho);
                    const double rho_left=(double)iterator.Data(rho,negative_axis_offsets(axis));
                    const double rho_face=(double).5*(rho_right+rho_left);
                    const double one_over_rho_face=(double)1./rho_face;
                    iterator.template Data<Data_array_type>(L_data_ptrs(axis))=one_over_rho_face*laplace_scale_nonuniform[level-1];
                    iterator.Data(diag) += (T)-1.*one_over_rho_face*laplace_scale_nonuniform[level-1];
                    iterator.Data(diag,negative_axis_offsets(axis)) += (T)-1.*one_over_rho_face*laplace_scale_nonuniform[level-1];
                } else if(flag & GRID_TOPOLOGY_HELPER<typename Flag_array_type::MASK>::Face_Active_Mask(axis)) {
                    const double rho_right=(double)iterator.Data(rho);
                    const double rho_left=(double)iterator.Data(rho,negative_axis_offsets(axis));
                    const double rho_face=(double).5*(rho_right+rho_left);
                    const double one_over_rho_face=(double)1./rho_face;
                    iterator.template Data<Data_array_type>(L_data_ptrs(axis))=one_over_rho_face*laplace_scale_uniform[level-1];
                    iterator.Data(diag) += (T)-1.*one_over_rho_face*laplace_scale_uniform[level-1];
                    iterator.Data(diag,negative_axis_offsets(axis)) += (T)-1.*one_over_rho_face*laplace_scale_uniform[level-1];
                }
            }
        }
    }
    
    // Accumulate diagonal ghosts
    GRID_HIERARCHY_PROJECTION<T_STRUCT,T,d>::Accumulate_Ghost_Values(hierarchy,flags_channel,diag_channel);

    // Clear out dirichlet entries
    for(int level=1;level<=levels;level++){
        Const_flag_array_type flags=hierarchy.Allocator(level).Get_Const_Array(flags_channel);
        VECTOR<void*,d> L_data_ptrs,U_data_ptrs;
        for(int v=1;v<=d;v++){
            L_data_ptrs(v)=hierarchy.Array(level,L_channels(v)).Get_Data_Ptr();
            U_data_ptrs(v)=hierarchy.Array(level,U_channels(v)).Get_Data_Ptr();}
        for(SPGrid_Block_Iterator<typename Data_array_type::MASK> iterator(hierarchy.Blocks(level));iterator.Valid();iterator.Next()){
            if(iterator.Data(flags) & SPGrid_Cell_Type_Interior){
                HIERARCHY_NEIGHBOR_ITERATOR<T_STRUCT,T,d> neighbor_iterator(false,hierarchy,CELL_ID(level,iterator.Offset()),flags_channel,diag_channel,diag_channel,L_channels,L_channels); // both as L channels
                for(;neighbor_iterator.Valid();neighbor_iterator.Next()){
                    if(((hierarchy.Array(neighbor_iterator.Neighbor().level,flags_channel)(neighbor_iterator.Neighbor().offset)) & SPGrid_Cell_Type_Dirichlet))
                        neighbor_iterator.Coefficient()=(T)0.;
                }
            }
        }
    }

    // Copy L to U
    for(int level=1;level<=levels;level++)
        for(int v=1;v<=d;v++)
            SPGrid_Computations::Copy<T_STRUCT,T,d>(hierarchy.Allocator(level),hierarchy.Blocks(level),L_channels(v),U_channels(v));
 
    // Clear ghosts and Dirichlet on diagonal
    for(int level=1;level<=levels;level++)
        SPGrid_Computations::Masked_Clear(hierarchy.Allocator(level),hierarchy.Blocks(level),diag_channel,flags_channel,(unsigned)(SPGrid_Cell_Type_Ghost|SPGrid_Cell_Type_Dirichlet));
}
//#####################################################################
// Explicit_Matrix_Multiply
//#####################################################################
template<class T_STRUCT,class T,int d> void BUILD_MATRIX<T_STRUCT,T,d>::
Explicit_Matrix_Multiply(Hierarchy_type& hierarchy,T T_STRUCT::* diag_channel,VECTOR<T T_STRUCT::*,d> L_channels,VECTOR<T T_STRUCT::*,d> U_channels,
                         unsigned T_STRUCT::* flags_channel,T T_STRUCT::* x_channel,T T_STRUCT::* y_channel)
{
    const int levels=hierarchy.Levels();

    // clear y
    for(int level=1;level<=levels;level++)
        SPGrid_Computations::Clear<T_STRUCT,T,d>(hierarchy.Allocator(level),hierarchy.Blocks(level),y_channel);
    
    // Propagate Ghosts
    GRID_HIERARCHY_PROJECTION<T_STRUCT,T,d>::Propagate_Ghost_Values(hierarchy,flags_channel,x_channel);
    
    // Mat mul per level
    for(int level=1;level<=levels;level++){
        Const_flag_array_type flags=hierarchy.Allocator(level).Get_Const_Array(flags_channel);
        Data_array_type x=hierarchy.Allocator(level).Get_Array(x_channel); 
        Const_data_array_type diag=hierarchy.Allocator(level).Get_Const_Array(diag_channel);
        VECTOR<void*,d> L_data_ptrs,U_data_ptrs;
        for(int v=1;v<=d;v++){
            L_data_ptrs(v)=hierarchy.Array(level,L_channels(v)).Get_Data_Ptr();
            U_data_ptrs(v)=hierarchy.Array(level,U_channels(v)).Get_Data_Ptr();}
        Data_array_type y=hierarchy.Allocator(level).Get_Array(y_channel);
        for(SPGrid_Block_Iterator<typename Data_array_type::MASK> iterator(hierarchy.Blocks(level));iterator.Valid();iterator.Next()){
            const unsigned flag=iterator.Data(flags);
            for(int v=1;v<=d;v++){
                std_array<int,d> index_offset;
                index_offset(v-1)=-1;
                unsigned long offset=Flag_array_mask::Linear_Offset(index_offset);
                if(flag & GRID_TOPOLOGY_HELPER<typename Flag_array_type::MASK>::Face_Active_Mask(v)){
                    iterator.Data(y)+=iterator.template Data<Data_array_type>(L_data_ptrs(v))*iterator.template Data<Data_array_type>(x,offset);
                    iterator.template Data<Data_array_type>(y,offset)+=iterator.template Data<Data_array_type>(L_data_ptrs(v))*iterator.Data(x);}}
            if(flag & SPGrid_Cell_Type_Active)
                iterator.Data(y)+=iterator.Data(diag)*iterator.Data(x);
        }
    }
    
    // Accumulate Ghosts
    GRID_HIERARCHY_PROJECTION<T_STRUCT,T,d>::Accumulate_Ghost_Values(hierarchy,flags_channel,y_channel);

    // Clear ghost and dirichlet cells of result
    for(int level=1;level<=levels;level++){
        SPGrid_Computations::Masked_Clear(hierarchy.Allocator(level),hierarchy.Blocks(level),y_channel,flags_channel,(unsigned)SPGrid_Cell_Type_Dirichlet);
        SPGrid_Computations::Masked_Clear(hierarchy.Allocator(level),hierarchy.Blocks(level),y_channel,flags_channel,(unsigned)SPGrid_Cell_Type_Ghost);}
}
//#####################################################################
// Verify_Matrix
//#####################################################################
template<class T_STRUCT,class T,int d> void BUILD_MATRIX<T_STRUCT,T,d>::
Verify_Matrix(Hierarchy_type& hierarchy,T T_STRUCT::* diag_channel,VECTOR<T T_STRUCT::*,d> L_channels,VECTOR<T T_STRUCT::*,d> U_channels,T T_STRUCT::* u_channel,
              T T_STRUCT::* temp1_channel,T T_STRUCT::* temp2_channel,unsigned T_STRUCT::* flags_channel)
{
    // VERIFY THAT ALL THE FOLLOWING ARE EQUIVALENT:
    // -- SPGRID_LAPLACE
    // -- EXPLICIT_MATRIX_MULTIPLY

    const int levels=hierarchy.Levels();

    // clear u and temp1,temp2 (results)
    for(int level=1;level<=levels;level++){
        SPGrid_Computations::Clear<T_STRUCT,T,d>(hierarchy.Allocator(level),hierarchy.Blocks(level),u_channel);
        SPGrid_Computations::Clear<T_STRUCT,T,d,2>(hierarchy.Allocator(level),hierarchy.Blocks(level),temp1_channel,temp2_channel);}

    // Randomize u
    const T random_number_range=(T)1.;
    RANDOM_NUMBERS<T> random_numbers;random_numbers.Set_Seed(1);
    for(int level=1;level<=levels;level++){
        Data_array_type u=hierarchy.Array(level,u_channel);
        for(SPGrid_Block_Iterator<typename Data_array_type::MASK> iterator(hierarchy.Blocks(level));iterator.Valid();iterator.Next())
            iterator.Data(u)=random_numbers.Get_Uniform_Number((T)-1.*random_number_range,random_number_range);
    }

    // Copy u --> temp_2
    for(int level=1;level<=levels;level++)
        SPGrid_Computations::Copy<T_STRUCT,T,d>(hierarchy.Allocator(level),hierarchy.Blocks(level),u_channel,temp2_channel);

    // Clear Dirichlet cells in u
    for(int level=1;level<=levels;level++){
        Const_flag_array_type flags=hierarchy.Allocator(level).Get_Const_Array(flags_channel);
        Data_array_type u=hierarchy.Allocator(level).Get_Array(u_channel);
        for(SPGrid_Block_Iterator<typename Data_array_type::MASK> iterator(hierarchy.Blocks(level));iterator.Valid();iterator.Next())
            if(!(iterator.Data(flags) & SPGrid_Cell_Type_Active))
                iterator.Data(u)=T();}

    // Compute Laplace(u) --> temp1
    GRID_HIERARCHY_PROJECTION<T_STRUCT,T,d>::Compute_Laplacian(hierarchy,flags_channel,u_channel,temp1_channel);
    
    // Copy temp_2 --> u
    for(int level=1;level<=levels;level++)
        SPGrid_Computations::Copy<T_STRUCT,T,d>(hierarchy.Allocator(level),hierarchy.Blocks(level),temp2_channel,u_channel);

    // Multiply Matrix : Laplace(u) --> temp_2
    Explicit_Matrix_Multiply(hierarchy,diag_channel,L_channels,U_channels,flags_channel,u_channel,temp2_channel);

    // Test it
    for(int level=1;level<=levels;level++){
        Data_array_type result1=hierarchy.Array(level,temp1_channel);
        Data_array_type result2=hierarchy.Array(level,temp2_channel);
        Const_flag_array_type flags=hierarchy.Allocator(level).Get_Const_Array(flags_channel);
        for(SPGrid_Block_Iterator<typename Data_array_type::MASK> iterator(hierarchy.Blocks(level));iterator.Valid();iterator.Next())
            if(!(iterator.Data(flags)&SPGrid_Cell_Type_Ghost)){
                T rel_abs_diff=abs(iterator.Data(result1)-iterator.Data(result2))/iterator.Data(result1);
                //PHYSBAM_ASSERT( iterator.Data(result1) == iterator.Data(result2) );
                //if(! (iterator.Data(result1) == iterator.Data(result2)) ){
                if(rel_abs_diff > 1e-5){
                    LOG::cout<<"level="<<level<<std::endl;
                    LOG::cout<<"index="<<iterator.Index()<<std::endl;
                    LOG::cout<<"explicit="<<iterator.Data(result2)<<std::endl;
                    LOG::cout<<"implicit="<<iterator.Data(result1)<<std::endl;
                    LOG::cout<<"relative diff="<<rel_abs_diff<<std::endl;
                    LOG::cout<<"-------------------"<<std::endl;}}
    }
}
//#####################################################################
// Generate_Conversion_Structures
//#####################################################################
template<class T_STRUCT,class T,int d> void BUILD_MATRIX<T_STRUCT,T,d>::
Generate_Conversion_Structures(HASHTABLE<CID,int>& active_dof_hash,ARRAY<CID>& active_dof_array,Hierarchy_type& hierarchy,unsigned T_STRUCT::* flags_channel)
{
    const int levels=hierarchy.Levels();
    for(int level=1;level<=levels;level++){
        Const_flag_array_type flags=hierarchy.Allocator(level).Get_Const_Array(flags_channel);
        for(SPGrid_Block_Iterator<typename Data_array_type::MASK> iterator(hierarchy.Blocks(level));iterator.Valid();iterator.Next())
            if(iterator.Data(flags) & SPGrid_Cell_Type_Active){
                // add to hash table
                int number = active_dof_hash.Size()+1;
                unsigned long offset = iterator.Offset();
                active_dof_hash.Insert(CID(level,offset),number);
                // add to array
                active_dof_array.Append(CID(level,offset));
                // check for consistency
                PHYSBAM_ASSERT(active_dof_array.m==active_dof_hash.Size());
            }
    }
    PHYSBAM_ASSERT(active_dof_array.m==active_dof_hash.Size());
}
//#####################################################################
// SPGrid_Channel_To_Vector
//#####################################################################
template<class T_STRUCT,class T,int d> void BUILD_MATRIX<T_STRUCT,T,d>::
SPGrid_Channel_To_Vector(ARRAY<CID>& active_dof_array,Hierarchy_type& hierarchy,T T_STRUCT::* data_channel,VECTOR_ND<T>& data_vec)
{
    const int levels=hierarchy.Levels();

    data_vec.Resize(active_dof_array.m);
    data_vec.Fill(T());
    for(int i=1;i<=active_dof_array.m;i++)
        data_vec(i)=hierarchy.Array(active_dof_array(i).x,data_channel)(active_dof_array(i).y);
}
//#####################################################################
// Vector_To_SPGrid_Channel
//#####################################################################
template<class T_STRUCT,class T,int d> void BUILD_MATRIX<T_STRUCT,T,d>::
Vector_To_SPGrid_Channel(ARRAY<CID>& active_dof_array,Hierarchy_type& hierarchy,T T_STRUCT::* data_channel,VECTOR_ND<T>& data_vec)
{
    PHYSBAM_ASSERT(active_dof_array.m==data_vec.Size());

    const int levels=hierarchy.Levels();

    for(int level=1;level<=levels;level++)
        SPGrid_Computations::Clear<T_STRUCT,T,d>(hierarchy.Allocator(level),hierarchy.Blocks(level),data_channel);
    
    for(int i=1;i<=active_dof_array.m;i++)
        hierarchy.Array(active_dof_array(i).x,data_channel)(active_dof_array(i).y)=data_vec(i);
}
//#####################################################################
// SPGrid_Channels_To_Matrix
//#####################################################################
template<class T_STRUCT,class T,int d> void BUILD_MATRIX<T_STRUCT,T,d>::
SPGrid_Channels_To_Matrix(Hierarchy_type& hierarchy,T T_STRUCT::* diag_channel,
    VECTOR<T T_STRUCT::*,d> L_channels,VECTOR<T T_STRUCT::*,d> U_channels,
    SPARSE_MATRIX_FLAT_NXN<T>& matrix,HASHTABLE<CID,int>& active_dof_hash,
    unsigned T_STRUCT::* flags_channel)
{
    typedef NEIGHBOR_STRUCT<T> NEIGHBOR;
    const int levels=hierarchy.Levels();

    ARRAY<int> row_lengths;
    row_lengths.Resize(active_dof_hash.Size());
    for(int level=1;level<=levels;level++){
        Const_flag_array_type flags=hierarchy.Allocator(level).Get_Const_Array(flags_channel);
        for(SPGrid_Block_Iterator<typename Data_array_type::MASK> iterator(hierarchy.Blocks(level));iterator.Valid();iterator.Next())
            if(iterator.Data(flags) & SPGrid_Cell_Type_Active){
                int my_id=active_dof_hash.Get(CID(level,iterator.Offset()));
                HIERARCHY_NEIGHBOR_ITERATOR<T_STRUCT,T,d> neighbor_iterator(true,hierarchy,CELL_ID(level,iterator.Offset()),flags_channel,diag_channel,diag_channel,L_channels,U_channels);
                row_lengths(my_id)=neighbor_iterator.Size();}}
    matrix.Set_Row_Lengths(row_lengths);

    for(int level=1;level<=levels;level++){
        Const_flag_array_type flags=hierarchy.Allocator(level).Get_Const_Array(flags_channel);        
        for(SPGrid_Block_Iterator<typename Data_array_type::MASK> iterator(hierarchy.Blocks(level));iterator.Valid();iterator.Next())
            if(iterator.Data(flags) & SPGrid_Cell_Type_Active){
                HIERARCHY_NEIGHBOR_ITERATOR<T_STRUCT,T,d> neighbor_iterator(true,hierarchy,CELL_ID(level,iterator.Offset()),flags_channel,diag_channel,diag_channel,L_channels,U_channels); // use separate L and U
                int my_id = active_dof_hash.Get(CID(level,iterator.Offset()));
                for(;neighbor_iterator.Valid();neighbor_iterator.Next()){
                    NEIGHBOR neighbor=neighbor_iterator.Neighbor();
                    if((hierarchy.Array(neighbor.level,flags_channel)(neighbor.offset)) & SPGrid_Cell_Type_Active){
                        int neighbor_id=active_dof_hash.Get(CID(neighbor.level,neighbor.offset));
                        matrix.Set_Element(my_id,neighbor_id,neighbor_iterator.Coefficient());}}
            }
    }
}
//#####################################################################
// Verify_Factorization
//#####################################################################
template<class T_STRUCT,class T,int d> void BUILD_MATRIX<T_STRUCT,T,d>::
Verify_Factorization(Hierarchy_type& hierarchy,T T_STRUCT::* diag_channel,VECTOR<T T_STRUCT::*,d> L_channels,VECTOR<T T_STRUCT::*,d> U_channels,unsigned T_STRUCT::* flags_channel,
                     T T_STRUCT::* temp0_channel,T T_STRUCT::* temp1_channel,T T_STRUCT::* temp2_channel,T T_STRUCT::* temp3_channel)
{
#if 0
    // ##########################
    // Verify IChol
    // ##########################

    const int levels=hierarchy.Levels();

    // Build PhysBAM SPARSE_NXN
    HASHTABLE<CID,int> active_dof_hash;
    ARRAY<CID> active_dof_array;
    SPARSE_MATRIX_FLAT_NXN<T> matrix;
    Generate_Conversion_Structures(active_dof_hash,active_dof_array,hierarchy,flags_channel);
    SPGrid_Channels_To_Matrix(hierarchy,diag_channel,L_channels,U_channels,matrix,active_dof_hash,flags_channel);

    PHYSBAM_ASSERT(active_dof_hash.Size() == matrix.n);

    // PhysBAM In-Place Incomplete Cholesky
    SPARSE_MATRIX_FLAT_NXN<T> physbam_factorized_matrix=matrix;
    physbam_factorized_matrix.In_Place_Incomplete_Cholesky_Factorization(false,(T).97,(T)1e-8,(T)1e-8);

#if 0
    {
        const SPARSE_MATRIX_FLAT_NXN<T>& A=physbam_factorized_matrix;
        const int n=A.n;
        const ARRAY<int>& offsets=A.offsets;
        for(int i=1;i<=n;i++)
          for(int index=offsets(i);index<A.diagonal_index(i);index++){
              const int j=A.A(index).j;
              const T Lij=A.A(index).a;
              const T Uji=A(j,i);
              const T Dinv_ii=A(i,i);
              const T Dinv_jj=A(j,j);
              if(fabs(Lij-Dinv_ii*Uji)>1e-2){
                  LOG::cout<<"i = "<<i<<std::endl;
                  LOG::cout<<"j = "<<j<<std::endl;
                  LOG::cout<<"Lij = "<<Lij<<std::endl;
                  LOG::cout<<"Uji = "<<Uji<<std::endl;
                  LOG::cout<<"Dinv_ii = "<<Dinv_ii<<std::endl;
                  LOG::cout<<"Dinv_jj = "<<Dinv_jj<<std::endl;
                  LOG::cout<<"fabs(Lij-Dinv_ii*Uji) = "<<fabs(Lij-Dinv_ii*Uji)<<std::endl;
                  PHYSBAM_FATAL_ERROR();
                  
              }
          }
      LOG::cout<<"Succeeded"<<std::endl;
      exit(0);        
    }
#endif

    // SPGrid In-Place Incomplete Cholesky
    HIERARCHY_PRECONDITIONER<T_STRUCT,T,d>::In_Place_Incomplete_Cholesky_Factorization(hierarchy,diag_channel,L_channels,U_channels);

    // Copy factorized version over to PhysBAM matrix
    SPARSE_MATRIX_FLAT_NXN<T> factorized_spgrid_to_physbam_matrix(matrix);
    factorized_spgrid_to_physbam_matrix*=(T)0.; // clear
    SPGrid_Channels_To_Matrix(hierarchy,diag_channel,L_channels,U_channels,factorized_spgrid_to_physbam_matrix,active_dof_hash,flags_channel); // use L & U channels

    // Check they are same
    PHYSBAM_ASSERT(factorized_spgrid_to_physbam_matrix.n == physbam_factorized_matrix.n);
    int num_discrepancies=0;
    const int matrix_n=factorized_spgrid_to_physbam_matrix.n;
    for(int i=1;i<=matrix_n;i++)
        for(int j=1;j<=matrix_n;j++){
            PHYSBAM_ASSERT( physbam_factorized_matrix.Element_Present(i,j) == factorized_spgrid_to_physbam_matrix.Element_Present(i,j) );
            if(physbam_factorized_matrix.Element_Present(i,j)){
                if( physbam_factorized_matrix(i,j) != factorized_spgrid_to_physbam_matrix(i,j)){
                    num_discrepancies++;
                    LOG::cout<<"Discrepancy In Factorization!!"<<std::endl;
                    LOG::cout<<"Element "<<i<<", "<<j<<std::endl;
                    LOG::cout<<"Element "<<Flag_array_type::MASK::LinearToCoord(active_dof_array(i).y)<<", "<<Flag_array_type::MASK::LinearToCoord(active_dof_array(j).y)<<std::endl;
                    LOG::cout<<"physbam_factorized_matrix("          <<i<<")("<<j<<")="<<physbam_factorized_matrix(i,j)          <<std::endl;
                    LOG::cout<<"factorized_spgrid_to_physbam_matrix("<<i<<")("<<j<<")="<<factorized_spgrid_to_physbam_matrix(i,j)<<std::endl;
                    LOG::cout<<"matrix("                             <<i<<")("<<j<<")="<<matrix(i,j)                              <<std::endl;
                    LOG::cout<<"------------------------------"<<std::endl;}}}
    LOG::cout<<"num_discrepancies="<<num_discrepancies<<std::endl;
    if(num_discrepancies > 0){LOG::cout<<"exiting due to Incomplete Cholesky Factorization Errors!"<<std::endl;exit(0);}

    // ##########################
    // Verify Back/Forward Sub
    // ##########################

    // Clear x,b
    for(int level=1;level<=levels;level++){
        SPGrid_Computations::Clear<T_STRUCT,T,d,2>(hierarchy.Allocator(level),hierarchy.Blocks(level),temp0_channel,temp1_channel);
        SPGrid_Computations::Clear<T_STRUCT,T,d,2>(hierarchy.Allocator(level),hierarchy.Blocks(level),temp2_channel,temp3_channel);}

    // b <-- random
    const T random_number_range=(T)1.;
    RANDOM_NUMBERS<T> random_numbers;random_numbers.Set_Seed(1);
    for(int level=1;level<=levels;level++){
        Data_array_type b=hierarchy.Array(level,temp1_channel);
        for(SPGrid_Block_Iterator<typename Data_array_type::MASK> iterator(hierarchy.Blocks(level));iterator.Valid();iterator.Next())
            iterator.Data(b)=random_numbers.Get_Uniform_Number((T)-1.*random_number_range,random_number_range);}

    // Copy b to b_vector (VECTOR_ND)    
    VECTOR_ND<T> b_vector;
    SPGrid_Channel_To_Vector(active_dof_array,hierarchy,temp1_channel,b_vector);
    VECTOR_ND<T> z_vector(b_vector); z_vector*=(T)0.;
    VECTOR_ND<T> x_vector(b_vector); x_vector*=(T)0.;
    physbam_factorized_matrix.Solve_Forward_Substitution (b_vector,z_vector,true);
    physbam_factorized_matrix.Solve_Backward_Substitution(z_vector,x_vector,false,true);

    // SPGrid version
    HIERARCHY_PRECONDITIONER<T_STRUCT,T,d>::Solve_Forward_Substitution(hierarchy,temp1_channel,temp2_channel,diag_channel,L_channels);
    for(int level=1;level<=levels;level++) 
        SPGrid_Computations::Masked_Clear(hierarchy.Allocator(level),hierarchy.Blocks(level),temp2_channel,flags_channel,(unsigned)(SPGrid_Cell_Type_Ghost|SPGrid_Cell_Type_Dirichlet));
    HIERARCHY_PRECONDITIONER<T_STRUCT,T,d>::Solve_Backward_Substitution(hierarchy,temp2_channel,temp0_channel,diag_channel,L_channels);

    // Copy over PhysBAM result to SPGrid (temp3)
    Vector_To_SPGrid_Channel(active_dof_array,hierarchy,temp3_channel,x_vector);

    // Compare
    const T tolerance=(T)1e-8;
    for(int level=1;level<=levels;level++){
        Const_flag_array_type flags=hierarchy.Allocator(level).Get_Const_Array(flags_channel);
        Const_data_array_type x_array=hierarchy.Allocator(level).Get_Const_Array(temp0_channel);
        for(SPGrid_Block_Iterator<typename Data_array_type::MASK> iterator(hierarchy.Blocks(level));iterator.Valid();iterator.Next())
            if(iterator.Data(flags) & SPGrid_Cell_Type_Active){
                CID cid(level,iterator.Offset());
                if( abs(iterator.Data(x_array) - x_vector(active_dof_hash.Get(cid))) > tolerance){
                    LOG::cout<<"Error in subtitution!"<<std::endl;
                }
            }
    }
    
    LOG::cout<<"Done checking factorization!"<<std::endl;
#endif
}
//#####################################################################
template class BUILD_MATRIX<FLUIDS_SIMULATION_DATA<float>,float,2>;
template class BUILD_MATRIX<FLUIDS_SIMULATION_DATA<float>,float,3>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class BUILD_MATRIX<FLUIDS_SIMULATION_DATA<double>,double,2>;
template class BUILD_MATRIX<FLUIDS_SIMULATION_DATA<double>,double,3>;
#endif
