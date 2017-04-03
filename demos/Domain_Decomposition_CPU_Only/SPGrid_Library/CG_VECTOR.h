//#####################################################################
//Removed For Anonymity: Copyright Authors of Submission pap173s1
//Removed For Anonymity: Copyright Authors of Submission pap173s1
//#####################################################################
// Class CG_VECTOR
//#####################################################################
#ifndef __CG_VECTOR__
#define __CG_VECTOR__

#include <PhysBAM_Tools/Krylov_Solvers/KRYLOV_VECTOR_BASE.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Vectors/VECTOR_2D.h>
#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>
#include <SPGrid/Core/SPGrid_Allocator.h>
#include <SPGrid/Core/SPGrid_Set.h>
#include <SPGrid/Tools/SPGrid_Block_Iterator.h>
#include <SPGrid_Fluids/Simulation/FLUIDS_SIMULATION_DATA.h>
#include "SPGRID_MULTIGRID_FLAGS.h"


using namespace SPGrid;

namespace PhysBAM{

template<class T_STRUCT,class T,int d,typename INDEX=int> class CG_SYSTEM;

template<class T_STRUCT,class T,int d>
class CG_VECTOR:public KRYLOV_VECTOR_BASE<T>
{
    typedef KRYLOV_VECTOR_BASE<T> BASE;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<T>::type Data_array_type;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<const T>::type Const_data_array_type;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<const unsigned>::type Const_flag_array_type;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<unsigned>::type Flag_array_type;
    typedef SPGrid_Allocator<T_STRUCT,d> SPG_Allocator;
    typedef SPGrid_Set<Flag_array_type> SPG_Set_Type;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<const unsigned>::mask T_MASK;
    typedef VECTOR<int,d> T_INDEX;

    SPG_Allocator& allocator;
    const SPG_Set_Type& set;

 public:
    T T_STRUCT::* field;  // should check to make sure these types are the same

    CG_VECTOR(SPG_Allocator& allocator_input,SPG_Set_Type& set_input,T T_STRUCT::* field_input)
        :allocator(allocator_input),set(set_input),field(field_input) {}
    
    static const CG_VECTOR& Cg_Vector(const BASE& base)
    {return dynamic_cast<const CG_VECTOR&>(base);}

    static const SPG_Allocator& Allocator(const BASE& base)
    {return dynamic_cast<const CG_VECTOR&>(base).allocator;}

    static SPG_Allocator& Allocator(BASE& base)
    {return dynamic_cast<CG_VECTOR&>(base).allocator;}

    static const SPG_Set_Type& Set(const BASE& base)
    {return dynamic_cast<const CG_VECTOR&>(base).set;}

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
