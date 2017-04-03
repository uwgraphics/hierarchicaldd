//#####################################################################
// Copyright 2012-2013, Sean Bauer, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class Parity_Helper
//#####################################################################
#ifndef __Parity_Helper__
#define __Parity_Helper__

namespace PhysBAM{

template<class MASK,int d> class Parity_Helper;

template<class MASK>
class Parity_Helper<MASK,2>
{
public:
    enum {
        GHOST_000 = MASK::template LinearOffset<0,0>::value,
        GHOST_010 = MASK::template LinearOffset<0,1>::value,
        GHOST_100 = MASK::template LinearOffset<1,0>::value,
        GHOST_110 = MASK::template LinearOffset<1,1>::value,
        
        GHOST_001 = MASK::template LinearOffset<1,1>::value,
        GHOST_011 = MASK::template LinearOffset<1,1>::value,
        GHOST_101 = MASK::template LinearOffset<1,1>::value,
        GHOST_111 = MASK::template LinearOffset<1,1>::value,

        X_EDGE = MASK::template LinearOffset<0,1>::value, 
        Y_EDGE = MASK::template LinearOffset<1,0>::value,

        X = 1,
        Y = 2,
        CELL = 0,
    };
};

template<class MASK>
class Parity_Helper<MASK,3>
{
public:
    enum {
        GHOST_000 = MASK::template LinearOffset<0,0,0>::value,
        GHOST_001 = MASK::template LinearOffset<0,0,1>::value,
        GHOST_010 = MASK::template LinearOffset<0,1,0>::value,
        GHOST_011 = MASK::template LinearOffset<0,1,1>::value,
        GHOST_100 = MASK::template LinearOffset<1,0,0>::value,
        GHOST_101 = MASK::template LinearOffset<1,0,1>::value,
        GHOST_110 = MASK::template LinearOffset<1,1,0>::value,
        GHOST_111 = MASK::template LinearOffset<1,1,1>::value,
        
        X_EDGE = MASK::template LinearOffset<0,1,1>::value, 
        Y_EDGE = MASK::template LinearOffset<1,0,1>::value, 
        Z_EDGE = MASK::template LinearOffset<1,1,0>::value,
         
        X_FACE = MASK::template LinearOffset<1,0,0>::value, 
        Y_FACE = MASK::template LinearOffset<0,1,0>::value, 
        Z_FACE = MASK::template LinearOffset<0,0,1>::value, 
        
        X = 1,
        Y = 2,
        Z = 3,
        CELL = 0,
    };
};
}

#endif
