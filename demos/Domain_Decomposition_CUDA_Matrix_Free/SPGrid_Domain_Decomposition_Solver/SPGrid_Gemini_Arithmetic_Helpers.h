//#####################################################################
//Removed For Anonymity: Copyright Authors of Submission pap173s1
//Removed For Anonymity: Copyright Authors of Submission pap173s1
//#####################################################################
#ifndef __SPGRID_GEMINI_ARITHMETIC_HELPER_H__
#define __SPGRID_GEMINI_ARITHMETIC_HELPER_H__

#include <Common_Tools/Math_Tools/RANGE_ITERATOR.h>
#include <SPGrid/Core/SPGrid_Allocator.h>
#include <SPGrid/Tools/SPGrid_Block_Iterator.h>
#include <SPGrid_Fluids/Grids/GRID_TOPOLOGY_HELPER.h>
#include <SPGrid_Fluids/Simulation/FLUIDS_SIMULATION_FLAGS.h>
#include <SPGrid_Fluids/Solvers/Domain_Decomposition/SPGRID_MULTIGRID_FLAGS.h>
#include "SPGrid_Gemini.h"

namespace SPGrid_Computations{

using namespace PhysBAM;
using namespace SPGrid;

template<class T_STRUCT,class T_DATA,int d> class SPGrid_Plus_Equal;
template<class T_STRUCT,class T_DATA,int d> class SPGrid_Minus_Equal;
template<class T_STRUCT,class T_DATA,int d> class SPGrid_Scale;
template<class T_STRUCT,class T_DATA,int d> class SPGrid_Multiply_Equal;
template<class T_STRUCT,class T_DATA,int d> class SPGrid_Copy;
template<class T_STRUCT,class T_DATA,int d> class SPGrid_Copy2;
template<class T_STRUCT,class T_DATA,class T_FLAG,int d> class SPGrid_Laplace;
template<class T_STRUCT,class T_DATA,class T_FLAG,int d> class SPGrid_Clear;

template<class T_STRUCT,class T_DATA,int d>
class SPGrid_Plus_Equal
{
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<const T_DATA>::type Const_Data_Array_Type;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<T_DATA>::type Data_Array_Type;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<T_DATA>::mask T_MASK;

    const SPGrid_Allocator<T_STRUCT,d>& allocator_in;
    T_DATA T_STRUCT::*const in_field;
    T_DATA T_STRUCT::* out_field;

public:
    SPGrid_Plus_Equal(const SPGrid_Allocator<T_STRUCT,d>& allocator_in_input,T_DATA T_STRUCT::*const in_field_input,T_DATA T_STRUCT::* out_field_input)
        :allocator_in(allocator_in_input),in_field(in_field_input),out_field(out_field_input)
    {}
    
    SPGrid_Plus_Equal(SPGrid_Allocator<T_STRUCT,d>& allocator,const std::pair<const unsigned long*,unsigned>& blocks,
                      const SPGrid_Allocator<T_STRUCT,d>& allocator_in_input,T_DATA T_STRUCT::*const in_field_input,T_DATA T_STRUCT::* out_field_input)
        :allocator_in(allocator_in_input),in_field(in_field_input),out_field(out_field_input)
    {Run(allocator,blocks);}
    
    void Run(SPGrid_Allocator<T_STRUCT,d>& allocator,const std::pair<const unsigned long*,unsigned>& blocks)const{
        Data_Array_Type data_out=allocator.Get_Array(out_field);
        Const_Data_Array_Type data_in=allocator_in.Get_Const_Array(in_field);
        for(SPGrid_Block_Iterator<T_MASK> iterator(blocks);iterator.Valid();iterator.Next()){
            iterator.Data(data_out)+=iterator.Data(data_in);}
    }
};

template<class T_STRUCT,class T_DATA,int d>
class SPGrid_Minus_Equal
{
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<const T_DATA>::type Const_Data_Array_Type;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<T_DATA>::type Data_Array_Type;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<T_DATA>::mask T_MASK;

    const SPGrid_Allocator<T_STRUCT,d>& allocator_in;
    T_DATA T_STRUCT::*const in_field;
    T_DATA T_STRUCT::* out_field;

public:
    SPGrid_Minus_Equal(const SPGrid_Allocator<T_STRUCT,d>& allocator_in_input,T_DATA T_STRUCT::*const in_field_input,T_DATA T_STRUCT::* out_field_input)
        :allocator_in(allocator_in_input),in_field(in_field_input),out_field(out_field_input)
    {}
    
    SPGrid_Minus_Equal(SPGrid_Allocator<T_STRUCT,d>& allocator,const std::pair<const unsigned long*,unsigned>& blocks,
                      const SPGrid_Allocator<T_STRUCT,d>& allocator_in_input,T_DATA T_STRUCT::*const in_field_input,T_DATA T_STRUCT::* out_field_input)
        :allocator_in(allocator_in_input),in_field(in_field_input),out_field(out_field_input)
    {Run(allocator,blocks);}
    
    void Run(SPGrid_Allocator<T_STRUCT,d>& allocator,const std::pair<const unsigned long*,unsigned>& blocks)const{
        Data_Array_Type data_out=allocator.Get_Array(out_field);
        Const_Data_Array_Type data_in=allocator_in.Get_Const_Array(in_field);
        for(SPGrid_Block_Iterator<T_MASK> iterator(blocks);iterator.Valid();iterator.Next()){
            iterator.Data(data_out)-=iterator.Data(data_in);}
    }
};

template<class T_STRUCT,class T_DATA,int d>
class SPGrid_Scale
{
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<const T_DATA>::type Const_Data_Array_Type;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<T_DATA>::type Data_Array_Type;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<T_DATA>::mask T_MASK;

    T_DATA c;
    T_DATA T_STRUCT::* out_field;

public:
    SPGrid_Scale(T_DATA T_STRUCT::* out_field_input,T_DATA c_in)
        :out_field(out_field_input),c(c_in)
    {}
    
    SPGrid_Scale(SPGrid_Allocator<T_STRUCT,d>& allocator,const std::pair<const unsigned long*,unsigned>& blocks,
                 T_DATA T_STRUCT::* out_field_input,T_DATA c_in)
        :out_field(out_field_input),c(c_in)
    {Run(allocator,blocks);}
    
    void Run(SPGrid_Allocator<T_STRUCT,d>& allocator,const std::pair<const unsigned long*,unsigned>& blocks)const{
        Data_Array_Type data_out=allocator.Get_Array(out_field);
        for(SPGrid_Block_Iterator<T_MASK> iterator(blocks);iterator.Valid();iterator.Next()){
            iterator.Data(data_out)*=c;}
    }
};

template<class T_STRUCT,class T_DATA,int d>
class SPGrid_Multiply_Equal
{
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<const T_DATA>::type Const_Data_Array_Type;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<T_DATA>::type Data_Array_Type;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<T_DATA>::mask T_MASK;

    const SPGrid_Allocator<T_STRUCT,d>& allocator_in;
    T_DATA T_STRUCT::*const in_field;
    T_DATA T_STRUCT::* out_field;

public:
    SPGrid_Multiply_Equal(const SPGrid_Allocator<T_STRUCT,d>& allocator_in_input,T_DATA T_STRUCT::*const in_field_input,T_DATA T_STRUCT::* out_field_input)
        :allocator_in(allocator_in_input),in_field(in_field_input),out_field(out_field_input)
    {}
    
    SPGrid_Multiply_Equal(SPGrid_Allocator<T_STRUCT,d>& allocator,const std::pair<const unsigned long*,unsigned>& blocks,
                          const SPGrid_Allocator<T_STRUCT,d>& allocator_in_input,T_DATA T_STRUCT::*const in_field_input,T_DATA T_STRUCT::* out_field_input)
        :allocator_in(allocator_in_input),in_field(in_field_input),out_field(out_field_input)
    {Run(allocator,blocks);}
    
    void Run(SPGrid_Allocator<T_STRUCT,d>& allocator,const std::pair<const unsigned long*,unsigned>& blocks)const{
        Data_Array_Type data_out=allocator.Get_Array(out_field);
        Const_Data_Array_Type data_in=allocator_in.Get_Const_Array(in_field);
        for(SPGrid_Block_Iterator<T_MASK> iterator(blocks);iterator.Valid();iterator.Next()){
            iterator.Data(data_out)*=iterator.Data(data_in);}
    }
};

template<class T_STRUCT,class T_DATA,int d>
class SPGrid_Copy
{
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<const T_DATA>::type Const_Data_Array_Type;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<T_DATA>::type Data_Array_Type;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<T_DATA>::mask T_MASK;

    const SPGrid_Allocator<T_STRUCT,d>& allocator_in;
    T_DATA c;
    T_DATA T_STRUCT::*const in_field;
    T_DATA T_STRUCT::* out_field;

public:
    SPGrid_Copy(const SPGrid_Allocator<T_STRUCT,d>& allocator_in_input,T_DATA T_STRUCT::*const in_field_input,T_DATA T_STRUCT::* out_field_input,T_DATA c_in)
        :allocator_in(allocator_in_input),in_field(in_field_input),out_field(out_field_input),c(c_in)
    {}
    
    SPGrid_Copy(SPGrid_Allocator<T_STRUCT,d>& allocator,const std::pair<const unsigned long*,unsigned>& blocks,
                const SPGrid_Allocator<T_STRUCT,d>& allocator_in_input,T_DATA T_STRUCT::*const in_field_input,T_DATA T_STRUCT::* out_field_input,T_DATA c_in)
        :allocator_in(allocator_in_input),in_field(in_field_input),out_field(out_field_input),c(c_in)
    {Run(allocator,blocks);}
    
    void Run(SPGrid_Allocator<T_STRUCT,d>& allocator,const std::pair<const unsigned long*,unsigned>& blocks)const{
        Data_Array_Type data_out=allocator.Get_Array(out_field);
        Const_Data_Array_Type data_in=allocator_in.Get_Const_Array(in_field);
        for(SPGrid_Block_Iterator<T_MASK> iterator(blocks);iterator.Valid();iterator.Next()){
            iterator.Data(data_out)=iterator.Data(data_in)*c;}
    }
};

template<class T_STRUCT,class T_DATA,int d>
class SPGrid_Copy2
{
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<const T_DATA>::type Const_Data_Array_Type;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<T_DATA>::type Data_Array_Type;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<T_DATA>::mask T_MASK;

    const SPGrid_Allocator<T_STRUCT,d>& allocator_in_1;
    const SPGrid_Allocator<T_STRUCT,d>& allocator_in_2;
    T_DATA c;
    T_DATA T_STRUCT::*const in_field_1;
    T_DATA T_STRUCT::*const in_field_2;
    T_DATA T_STRUCT::* out_field;

public:
    SPGrid_Copy2(const SPGrid_Allocator<T_STRUCT,d>& allocator_in_1_input,const SPGrid_Allocator<T_STRUCT,d>& allocator_in_2_input,
                 T_DATA T_STRUCT::*const in_field_1_input,T_DATA T_STRUCT::*const in_field_2_input,T_DATA T_STRUCT::* out_field_input,T_DATA c_in)
        :allocator_in_1(allocator_in_1_input),allocator_in_2(allocator_in_2_input),in_field_1(in_field_1_input),in_field_2(in_field_2_input),out_field(out_field_input),c(c_in)
    {}
    
    SPGrid_Copy2(SPGrid_Allocator<T_STRUCT,d>& allocator,const std::pair<const unsigned long*,unsigned>& blocks,
                 const SPGrid_Allocator<T_STRUCT,d>& allocator_in_1_input,const SPGrid_Allocator<T_STRUCT,d>& allocator_in_2_input,
                 T_DATA T_STRUCT::*const in_field_1_input,T_DATA T_STRUCT::*const in_field_2_input,T_DATA T_STRUCT::* out_field_input,T_DATA c_in)
        :allocator_in_1(allocator_in_1_input),allocator_in_2(allocator_in_2_input),in_field_1(in_field_1_input),in_field_2(in_field_2_input),out_field(out_field_input),c(c_in)
    {Run(allocator,blocks);}
    
    void Run(SPGrid_Allocator<T_STRUCT,d>& allocator,const std::pair<const unsigned long*,unsigned>& blocks)const{
        Data_Array_Type data_out=allocator.Get_Array(out_field);
        Const_Data_Array_Type data_in_1=allocator_in_1.Get_Const_Array(in_field_1);
        Const_Data_Array_Type data_in_2=allocator_in_2.Get_Const_Array(in_field_2);
        for(SPGrid_Block_Iterator<T_MASK> iterator(blocks);iterator.Valid();iterator.Next()){
            iterator.Data(data_out)=iterator.Data(data_in_1)*c+iterator.Data(data_in_2);}
    }
};

template<class T_STRUCT,class T_DATA,class T_FLAGS,int d>
class SPGrid_Laplace
{
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<const T_DATA>::type Const_Data_Array_Type;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<T_DATA>::type Data_Array_Type;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<const T_FLAGS>::type Const_Flag_Array_Type;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<T_FLAGS>::mask T_MASK;

    const SPGrid_Allocator<T_STRUCT,d>& allocator_u;
    const SPGrid_Allocator<T_STRUCT,d>& allocator_flags;
    T_FLAGS active_flag;
    T_DATA T_STRUCT::* u_field;
    T_DATA T_STRUCT::* Lu_field;
    T_FLAGS T_STRUCT::* flags_field;

public:
    SPGrid_Laplace(const SPGrid_Allocator<T_STRUCT,d>& allocator_u_input,const SPGrid_Allocator<T_STRUCT,d>& allocator_flags_input,
                   T_DATA T_STRUCT::* u_field_input,T_DATA T_STRUCT::* Lu_field_input,T_FLAGS T_STRUCT::* flags_field_input,
                   T_FLAGS active_flag_input)
        :allocator_u(allocator_u_input),allocator_flags(allocator_flags_input),
         u_field(u_field_input),Lu_field(Lu_field_input),flags_field(flags_field_input),active_flag(active_flag_input)
    {}
    
    SPGrid_Laplace(SPGrid_Allocator<T_STRUCT,d>& allocator,const std::pair<const unsigned long*,unsigned>& blocks,
                   const SPGrid_Allocator<T_STRUCT,d>& allocator_u_input,const SPGrid_Allocator<T_STRUCT,d>& allocator_flags_input,
                   T_DATA T_STRUCT::* u_field_input,T_DATA T_STRUCT::* Lu_field_input,T_FLAGS T_STRUCT::* flags_field_input,
                   T_FLAGS active_flag_input)
        :allocator_u(allocator_u_input),allocator_flags(allocator_flags_input),
         u_field(u_field_input),Lu_field(Lu_field_input),flags_field(flags_field_input),active_flag(active_flag_input)
    {Run(allocator,blocks);}
    
    void Run(SPGrid_Allocator<T_STRUCT,d>& allocator,const std::pair<const unsigned long*,unsigned>& blocks) const
    {
        static const int number_of_face_neighbors=GRID_TOPOLOGY_HELPER<T_MASK>::faces_per_cell;
        unsigned long face_neighbor_offsets[number_of_face_neighbors];
        GRID_TOPOLOGY_HELPER<T_MASK>::Face_Neighbor_Offsets(face_neighbor_offsets); // order is (-x, x, -y, y, -z, z)

        Const_Data_Array_Type u=allocator_u.Get_Const_Array(u_field);
        Data_Array_Type Lu=allocator.Get_Array(Lu_field);
        Const_Flag_Array_Type flags=allocator_flags.Get_Const_Array(flags_field);    
        
        for (SPGrid_Block_Iterator<T_MASK> iterator(blocks);iterator.Valid();iterator.Next()){
            unsigned flag = iterator.Data(flags);
            if(flag&active_flag){
                double cell_value=(double)(iterator.Data(u));
                double result=(double)0.;
                for(int face=0;face<number_of_face_neighbors;face++){
                    unsigned long offset=face_neighbor_offsets[face];                    
                    if(flag&(SPGrid_Solver_Face_Minus_X_Active<<face)){
                        double neighbor_value=(double)(iterator.Data(u,offset));
                        result-=(neighbor_value-cell_value);}}
                iterator.Data(Lu) = (T_DATA)result;}}
    }
};

}
#endif
