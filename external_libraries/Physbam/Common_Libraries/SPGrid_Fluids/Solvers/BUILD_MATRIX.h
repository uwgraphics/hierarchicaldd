//#####################################################################
// Copyright 2013, Raj Setaluri, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __BUILD_MATRIX_h__
#define __BUILD_MATRIX_h__

#include <SPGrid/Core/SPGrid_Mask.h>
#include <SPGrid_Fluids/Grids/GRID_HIERARCHY.h>
#include <SPGrid_Fluids/Solvers/CELL_ID.h>
#include <PhysBAM_Tools/Data_Structures/PAIR.h>
#include <PhysBAM_Tools/Data_Structures/HASHTABLE.h>
#include <PhysBAM_Tools/Data_Structures/DATA_STRUCTURES_FORWARD.h>
#include <PhysBAM_Tools/Vectors/VECTOR_ND.h>
#include <PhysBAM_Tools/Matrices/SPARSE_MATRIX_FLAT_NXN.h>

namespace PhysBAM{
inline HASH Hash_Reduce(const PAIR<int,unsigned long>& key)
{STATIC_ASSERT(sizeof(unsigned long)==2*sizeof(int));
union {unsigned long ul;int i[2];} raw;raw.ul=key.y;
return HASH(key.x,raw.i[0],raw.i[1]);}
}

namespace PhysBAM{
using namespace SPGrid;
template<class T_STRUCT,class T,int d>
//#####################################################################
// Class BUILD_MATRIX
//#####################################################################
class BUILD_MATRIX
{
    typedef GRID_HIERARCHY<T_STRUCT,T,d> Hierarchy_type;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<unsigned>::type Flag_array_type;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<const unsigned>::type Const_flag_array_type;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<T>::type Data_array_type;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<const T>::type Const_data_array_type;
    typedef typename Flag_array_type::MASK Flag_array_mask;
    typedef PAIR<int,unsigned long> CID;

// #############################################################################
public:
    static void Build_Matrix(Hierarchy_type& hierarchy,T T_STRUCT::* diag_channel,VECTOR<T T_STRUCT::*,d> L_channels,VECTOR<T T_STRUCT::*,d> U_channels,unsigned T_STRUCT::* flags_channel);
    static void Build_Variable_Beta_Matrix(Hierarchy_type& hierarchy,T T_STRUCT::* diag_channel,VECTOR<T T_STRUCT::*,d> L_channels,VECTOR<T T_STRUCT::*,d> U_channels,T T_STRUCT::* rho_channel,unsigned T_STRUCT::* flags_channel);
    static void Explicit_Matrix_Multiply(Hierarchy_type& hierarchy,T T_STRUCT::* diag_channel,VECTOR<T T_STRUCT::*,d> L_channels,VECTOR<T T_STRUCT::*,d> U_channels,
                                         unsigned T_STRUCT::* flags_channel,T T_STRUCT::* x_channel,T T_STRUCT::* y_channel);
    static void Verify_Matrix(Hierarchy_type& hierarchy,T T_STRUCT::* diag_channel,VECTOR<T T_STRUCT::*,d> L_channels,VECTOR<T T_STRUCT::*,d> U_channels,
                              T T_STRUCT::* u_channel,T T_STRUCT::* temp1_channel,T T_STRUCT::* temp2_channel,unsigned T_STRUCT::* flags_channel);
    static void Generate_Conversion_Structures(HASHTABLE<CID,int>& active_dof_hash,ARRAY<CID>& active_dof_array,Hierarchy_type& hierarchy,unsigned T_STRUCT::* flags_channel);
    static void SPGrid_Channel_To_Vector(ARRAY<CID>& active_dof_array,Hierarchy_type& hierarchy,T T_STRUCT::* data_channel,VECTOR_ND<T>& data_vec);
    static void Vector_To_SPGrid_Channel(ARRAY<CID>& active_dof_array,Hierarchy_type& hierarchy,T T_STRUCT::* data_channel,VECTOR_ND<T>& data_vec);
    static void SPGrid_Channels_To_Matrix(Hierarchy_type& hierarchy,T T_STRUCT::* diag_channel,
                                          VECTOR<T T_STRUCT::*,d> L_channels,VECTOR<T T_STRUCT::*,d> U_channels,
                                          SPARSE_MATRIX_FLAT_NXN<T>& matrix,HASHTABLE<CID,int>& active_dof_hash,
                                          unsigned T_STRUCT::* flags_channel);
    static void Verify_Factorization(Hierarchy_type& hierarchy,T T_STRUCT::* diag_channel,VECTOR<T T_STRUCT::*,d> L_channels,VECTOR<T T_STRUCT::*,d> U_channels,unsigned T_STRUCT::* flags_channel,
                                     T T_STRUCT::* temp0_channel,T T_STRUCT::* temp1_channel,T T_STRUCT::* temp2_channel,T T_STRUCT::* temp3_channel);
// #############################################################################
};
}
#endif
