//#####################################################################
// Copyright 2012-2013, Sean Bauer, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class NEIGHBOR_STRUCT
//#####################################################################
#ifndef __NEIGHBOR_STRUCT__
#define __NEIGHBOR_STRUCT__

//#include "GRID_HIERARCHY.h"
#include "CELL_ID.h"
#include <stddef.h>

namespace PhysBAM{

template<class T>
struct NEIGHBOR_STRUCT
{
    int level;
    unsigned long offset;
    T* data;
    T* coefficient;

    NEIGHBOR_STRUCT(int l, unsigned long o, T* d, T* c) : level(l), offset(o), data(d), coefficient(c) {}
    NEIGHBOR_STRUCT(int l, unsigned long o) : level(l), offset(o), data(NULL), coefficient(NULL) {}
    NEIGHBOR_STRUCT(CELL_ID& cid, T* d, T* c) : level(cid.level), offset(cid.offset), data(d), coefficient(c) {}
    NEIGHBOR_STRUCT(CELL_ID& cid) : level(cid.level), offset(cid.offset), data(NULL), coefficient(NULL) {}
    NEIGHBOR_STRUCT() : level(0), offset(0), data(NULL), coefficient(NULL) {}

    // Might be best to define a single comparison operator on CELL_ID and simply have cellids be the basic form of identification
    bool operator<(const NEIGHBOR_STRUCT& n) const
    {
        if(level==n.level) return offset<n.offset;
        return level<n.level;
    }

    CELL_ID Get_CID()
    {
        return CELL_ID(level,offset);
    }
};

}
#endif
