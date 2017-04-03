//#####################################################################
// Copyright 2013, Raj Setaluri, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __FLUIDS_SIMULATION_DATA_H__
#define __FLUIDS_SIMULATION_DATA_H__

namespace PhysBAM{

template<class T>
struct FLUIDS_SIMULATION_DATA
{
    unsigned flags;
    T ch0;  // u_x
    T ch1;  // u_y
    T ch2;  // u_z (3D only)
    T ch3;  // density
    T ch4;  // un_x
    T ch5;  // un_y
    T ch6;  // un_z
    T ch7;  // weights (averaging)
    T ch8;
    T ch9;
    T ch10;
    T ch11;
    T ch12;
    T ch13;
    T ch14;
};

}
#endif
