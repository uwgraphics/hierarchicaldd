//#####################################################################
// Copyright 2012-2013, Sean Bauer, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class CELL_ID
//#####################################################################
#ifndef __CELL_ID__
#define __CELL_ID__

namespace PhysBAM{

struct CELL_ID
{
    int level;
    unsigned long offset;

    CELL_ID(int l, unsigned long o) : level(l), offset(o) {}
    CELL_ID() : level(0), offset(0) {}

    bool operator<(CELL_ID cid) const {
        if(level==cid.level) return offset<cid.offset;
        return level<cid.level;
    }
   
    bool operator==(CELL_ID cid) const {
        if(level==cid.level) return offset==cid.offset;
        return false;
    }
};
}
#endif
