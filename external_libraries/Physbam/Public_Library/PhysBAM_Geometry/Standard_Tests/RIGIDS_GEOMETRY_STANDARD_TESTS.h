//#####################################################################
// Copyright 2009.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RIGIDS_GEOMETRY_STANDARD_TESTS
//#####################################################################
#ifndef __RIGIDS_GEOMETRY_STANDARD_TESTS__
#define __RIGIDS_GEOMETRY_STANDARD_TESTS__

#include <PhysBAM_Tools/Data_Structures/DATA_STRUCTURES_FORWARD.h>
#include <PhysBAM_Tools/Vectors/VECTOR_FORWARD.h>
#include <string>

namespace PhysBAM{

template<class TV> class EXAMPLE;
template<class TV> class RIGID_GEOMETRY;
template<class TV> class RIGID_GEOMETRY_COLLECTION;

template<class TV>
class RIGIDS_GEOMETRY_STANDARD_TESTS
{
    typedef typename TV::SCALAR T;
public:
    EXAMPLE<TV>& example;
    RIGID_GEOMETRY_COLLECTION<TV>& rigid_geometry_collection;

    RIGIDS_GEOMETRY_STANDARD_TESTS(EXAMPLE<TV>& example_input,RIGID_GEOMETRY_COLLECTION<TV>& rigid_body_collection_input);

//#####################################################################
    RIGID_GEOMETRY<TV>& Add_Rigid_Body(const std::string& name,const T scaling_factor,const T friction,const bool read_implicit=true,const bool always_read_object=false);
    RIGID_GEOMETRY<TV>& Add_Ground(const T friction=(T).3,const T height=0,const T scale=(T)1);
    RIGID_GEOMETRY<TV>& Add_Analytic_Box(const VECTOR<T,1>& scaling_factor);
    RIGID_GEOMETRY<TV>& Add_Analytic_Box(const VECTOR<T,2>& scaling_factor,int segments_per_side=1);
    RIGID_GEOMETRY<TV>& Add_Analytic_Box(const VECTOR<T,3>& scaling_factor);
    RIGID_GEOMETRY<TV>& Add_Analytic_Ellipse(const VECTOR<T,2> radii,const T density,int levels=4);
//#####################################################################
};
}
#endif
