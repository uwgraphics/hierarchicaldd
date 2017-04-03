//#####################################################################
// Copyright 2011, Sean Bauer, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class CG_VECTOR
//#####################################################################
#ifndef __CG_VECTOR__
#define __CG_VECTOR__

#include <PhysBAM_Tools/Krylov_Solvers/KRYLOV_VECTOR_BASE.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Vectors/VECTOR_2D.h>
#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>
#include <SPGrid_Fluids/Simulation/FLUIDS_SIMULATION_FLAGS.h>
#include <SPGrid_Fluids/Grids/GRID_HIERARCHY.h>

//using namespace SPGrid;

namespace PhysBAM{

template<class T_STRUCT,class T,int d> class CG_SYSTEM;

template<class T_STRUCT,class T,int d>
class CG_VECTOR:public KRYLOV_VECTOR_BASE<T>
{
    typedef KRYLOV_VECTOR_BASE<T> BASE;

    typedef GRID_HIERARCHY<T_STRUCT,T,d> Hierarchy_type;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<T>::type Data_array_type;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<const T>::type Const_data_array_type;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<const unsigned>::type Const_flag_array_type;

    Hierarchy_type& hierarchy;
    T T_STRUCT::* field;  // should check to make sure these types are the same

public:
    CG_VECTOR(Hierarchy_type& hierarchy_input, T T_STRUCT::* field_input)
        :hierarchy(hierarchy_input),field(field_input) {}

    static const CG_VECTOR& Cg_Vector(const BASE& base)
    {return dynamic_cast<const CG_VECTOR&>(base);}

    static const Hierarchy_type& Hierarchy(const BASE& base)
    {return dynamic_cast<const CG_VECTOR&>(base).hierarchy;}

    static Hierarchy_type& Hierarchy(BASE& base)
    {return dynamic_cast<CG_VECTOR&>(base).hierarchy;}

    friend class CG_SYSTEM<T_STRUCT,T,d>;

//#####################################################################
    BASE& operator+=(const BASE& bv);
    BASE& operator-=(const BASE& bv);
    BASE& operator*=(const T a);
    void Copy(const T c,const BASE& bv);
    void Copy(const T c1,const BASE& bv1,const BASE& bv2);
    int Raw_Size() const;
    T& Raw_Get(int i);
//#####################################################################   
};
}
#endif
