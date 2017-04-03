//#####################################################################
//Removed For Anonymity: Copyright Authors of Submission pap173s1
//Removed For Anonymity: Copyright Authors of Submission pap173s1
//#####################################################################
#ifndef __SPGRID_GEMINI_H__
#define __SPGRID_GEMINI_H__

#include <SPGrid/Core/SPGrid_Allocator.h>
#include <SPGrid_Fluids/Grids/GRID_TOPOLOGY_HELPER.h>
#include <SPGrid/Tools/SPGrid_Block_Iterator.h>
#include <Common_Tools/Math_Tools/RANGE_ITERATOR.h>

namespace SPGrid{
template<class T_STRUCT,int d>
class SPGrid_Gemini{
public:
    //Assumes that T_STRUCT_FIRST holds the flags data. The flags channel and ch3 are unioned...so don't assess them in the wrong allocator!
    typedef SPGrid_Allocator<T_STRUCT,d> SPG_Allocator;
    typedef typename SPG_Allocator::template Array<unsigned>::type SPG_Gemini_Flags_Array_Type;
    typedef SPGrid_Set<SPG_Gemini_Flags_Array_Type> SPG_Gemini_Set_Type;

    SPG_Allocator* allocator_first;
    SPG_Allocator* allocator_second;
    SPG_Gemini_Set_Type* set;
    SPGrid_Gemini()
        :allocator_first(NULL),allocator_second(NULL),set(NULL){}
    ~SPGrid_Gemini(){
        if(allocator_first) delete allocator_first;
        if(allocator_second) delete allocator_second;
        if(set) delete set;
    }
    SPG_Allocator& First_Allocator(){
        return *allocator_first;
    }
    const SPG_Allocator& First_Allocator()const{
        return *allocator_first;
    }
    SPG_Allocator& Second_Allocator(){
        return *allocator_second;
    }
    const SPG_Allocator& Second_Allocator()const{
        return *allocator_second;
    }
    SPG_Gemini_Set_Type& Set(){
        return *set;
    }
    const SPG_Gemini_Set_Type& Set()const{
        return *set;
    }
};
}
#endif
