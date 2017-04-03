//#####################################################################
// Copyright 2015, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license 
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class STENCIL_ITERATOR
//#####################################################################
#ifndef __STENCIL_ITERATOR__
#define __STENCIL_ITERATOR__

#include <Common_Tools/Grids_Uniform_PDE_Linear/STENCIL.h>

namespace PhysBAM{
template<class T,int d>
class STENCIL_ITERATOR
{
    typedef typename IF<IS_CONST<T>::value,const STENCIL<typename REMOVE_CONST<T>::TYPE,d>,STENCIL<T,d> >::TYPE T_STENCIL;

    T_STENCIL& stencil;
    int current_index;
public:    

    STENCIL_ITERATOR(T_STENCIL& stencil_input)
        :stencil(stencil_input)
    {Reset();}

    void Reset()
    {current_index=1;}

    const VECTOR<int,d>& Key() const
    {return stencil.entries(current_index).x;}

    T& Data()
    {return stencil.entries(current_index).y;}

    const T& Data() const
    {return stencil.entries(current_index).y;}

    bool Valid() const
    {return current_index<=stencil.Size();}

    void Next()
    {assert(Valid());current_index++;}
//#####################################################################
};
}
#endif
