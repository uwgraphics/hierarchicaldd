//#####################################################################
// Copyright 2013, Raj Setaluri, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __SPGRID_MULTIGRID_DATA_H__
#define __SPGRID_MULTIGRID_DATA_H__

namespace SPGrid{

template<class T>
struct SPGRID_MULTIGRID_DATA
{
    unsigned flags;
    T ch0;  // u_x
    T ch1;  // u_y
    T ch2;  // u_z (3D only)
};

}
#endif
