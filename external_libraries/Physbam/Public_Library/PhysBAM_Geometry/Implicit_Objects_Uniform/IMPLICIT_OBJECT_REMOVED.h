//#####################################################################
// Copyright 2007, Michael Lentine, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class IMPLICIT_OBJECT_REMOVED
//#####################################################################
#ifndef __IMPLICIT_OBJECT_REMOVED__
#define __IMPLICIT_OBJECT_REMOVED__

#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/GRID_ARRAYS_POLICY_UNIFORM.h>
#include <PhysBAM_Tools/Matrices/MATRIX_3X3.h>
#include <PhysBAM_Geometry/Basic_Geometry/BASIC_GEOMETRY_POLICY.h>
#include <PhysBAM_Geometry/Basic_Geometry/ORIENTED_BOX.h>
#include <PhysBAM_Geometry/Basic_Geometry/RAY.h>
#include <PhysBAM_Geometry/Implicit_Objects/IMPLICIT_OBJECT.h>

namespace PhysBAM{

template<class TV>
class IMPLICIT_OBJECT_REMOVED:public IMPLICIT_OBJECT<TV>
{
    typedef typename TV::SCALAR T;
    typedef typename GRID_ARRAYS_POLICY<GRID<TV> >::ARRAYS_SCALAR T_ARRAYS_SCALAR;
    typedef typename BASIC_GEOMETRY_POLICY<TV>::ORIENTED_BOX T_ORIENTED_BOX;
    typedef typename MATRIX_POLICY<TV>::SYMMETRIC_MATRIX T_SYMMETRIC_MATRIX;
    enum WORKAROUND {d=TV::m};
public:
    typedef IMPLICIT_OBJECT<TV> BASE;
    using BASE::box;

    IMPLICIT_OBJECT<TV> *implicit_object1,*implicit_object2;
    bool owns_implicit_object1,owns_implicit_object2;
    T alpha;
    
    IMPLICIT_OBJECT_REMOVED(IMPLICIT_OBJECT<TV>* implicit_object1_input,bool owns_implicit_object1_input,IMPLICIT_OBJECT<TV>* implicit_object2_input, bool owns_implicit_object2_input)
        :implicit_object1(implicit_object1_input),implicit_object2(implicit_object2_input),owns_implicit_object1(owns_implicit_object1_input),owns_implicit_object2(owns_implicit_object2_input)
    {}

    virtual ~IMPLICIT_OBJECT_REMOVED()
    {if(owns_implicit_object1) delete implicit_object1;if(owns_implicit_object2) delete implicit_object2;}

    void Set_Weights(T alpha_input)
    {alpha=alpha_input;}

    void Update_Box() PHYSBAM_OVERRIDE
    {implicit_object1->Update_Box();implicit_object2->Update_Box();box=implicit_object1->box;}

    void Update_Minimum_Cell_Size(const int maximum_depth=0) PHYSBAM_OVERRIDE
    {implicit_object1->Update_Minimum_Cell_Size(maximum_depth);implicit_object2->Update_Minimum_Cell_Size(maximum_depth);}

    // TODO: box is not used.  Is it still needed?
    T Minimum_Cell_Size_Within_Box(const RANGE<TV>& box) const PHYSBAM_OVERRIDE
    {return min(implicit_object1->Minimum_Cell_Size_Within_Box(box),implicit_object2->Minimum_Cell_Size_Within_Box(box));} 

    T operator()(const TV& location) const PHYSBAM_OVERRIDE
    {return max((*implicit_object1)(location),-(*implicit_object2)(location));}

    T Signed_Distance(const TV& location) const PHYSBAM_OVERRIDE
    {return (*this)(location);} // to make this class compatible with the geometry classes

    T Extended_Phi(const TV& location) const PHYSBAM_OVERRIDE
    {return Extended_Value(location);}

    T Phi_Secondary(const TV& location) const PHYSBAM_OVERRIDE
    {return max(implicit_object1->Phi_Secondary(location),-implicit_object2->Phi_Secondary(location));}

    T Phi_With_Index(const TV& location,int& index) const
    {if((*implicit_object1)(location)>-(*implicit_object2)(location)){index=1;return (*implicit_object1)(location);}
    else{index=2;return -(*implicit_object2)(location);}}

    TV Normal(const TV& location,const int aggregate=-1) const PHYSBAM_OVERRIDE
    {assert((aggregate >= 1 && aggregate <= 2*GRID<TV>::dimension) || aggregate == -1);
    if(aggregate != -1) return box.Normal(aggregate);
    else{
        int index=0;Phi_With_Index(location,index);return index==1?implicit_object1->Normal(location):-implicit_object2->Normal(location);}}

    TV Extended_Normal(const TV& location,const int aggregate=-1) const PHYSBAM_OVERRIDE
    {assert((aggregate >= 1 && aggregate <= 2*GRID<TV>::dimension) || aggregate == -1);
    if(aggregate != -1) return box.Normal(aggregate);
    else{
        int index=0;Phi_With_Index(location,index);return index==1?implicit_object1->Extended_Normal(location):-implicit_object2->Extended_Normal(location);}}

    void Compute_Normals() PHYSBAM_OVERRIDE
    {implicit_object1->Compute_Normals();implicit_object2->Compute_Normals();}

    void Compute_Cell_Minimum_And_Maximum(const bool recompute_if_exists=true) PHYSBAM_OVERRIDE
    {implicit_object1->Compute_Cell_Minimum_And_Maximum(recompute_if_exists);implicit_object2->Compute_Cell_Minimum_And_Maximum(recompute_if_exists);}
    
    void Rescale(const T scaling_factor) PHYSBAM_OVERRIDE
    {implicit_object1->Rescale(scaling_factor);implicit_object2->Rescale(scaling_factor);}

    void Inflate(const T inflation_distance) PHYSBAM_OVERRIDE
    {implicit_object1->Inflate(inflation_distance);implicit_object2->Inflate(inflation_distance);}

    bool Lazy_Inside(const TV& location,const T contour_value=0) const PHYSBAM_OVERRIDE
    {return box.Inside(location,-contour_value) && (*this)(location)<=contour_value;}

    bool Lazy_Inside_And_Value(const TV& location,T& phi,const T contour_value=0) const PHYSBAM_OVERRIDE
    {if(box.Inside(location,-contour_value)){phi=(*this)(location);if(phi<=contour_value) return true;} return false;}

    bool Lazy_Inside_Extended_Levelset(const TV& location,const T contour_value=0) const PHYSBAM_OVERRIDE
    {return box.Inside(location,-contour_value) && Extended_Value(location)<=contour_value;}

    bool Lazy_Inside_Extended_Levelset_And_Value(const TV& location,T& phi_value,const T contour_value=0) const PHYSBAM_OVERRIDE
    {if(box.Inside(location,-contour_value)){phi_value=Extended_Phi(location);if(phi_value<=contour_value) return true;} return false;}

    bool Lazy_Outside(const TV& location,const T contour_value=0) const PHYSBAM_OVERRIDE
    {return !Lazy_Inside(location,contour_value);}

    bool Lazy_Outside_Extended_Levelset(const TV& location,const T contour_value=0) const PHYSBAM_OVERRIDE
    {return !Lazy_Inside_Extended_Levelset(location,contour_value);}

    bool Lazy_Outside_Extended_Levelset_And_Value(const TV& location,T& phi_value,const T contour_value=0) const PHYSBAM_OVERRIDE
    {return !Lazy_Inside_Extended_Levelset_And_Value(location,phi_value,contour_value);}

    T Min_Phi() const PHYSBAM_OVERRIDE
    {PHYSBAM_NOT_IMPLEMENTED();}

    T Minimum_Cell_Size() const PHYSBAM_OVERRIDE
    {return min(implicit_object1->Minimum_Cell_Size(),implicit_object2->Minimum_Cell_Size());}

    T Value(const TV& location) const
    {return (*this)(location);}

    T Extended_Value(const TV& location) const
    {return Extended_Phi(location);}

    bool Intersection(RAY<TV>& ray,const T thickness) const PHYSBAM_OVERRIDE
    {PHYSBAM_NOT_IMPLEMENTED();}

    static std::string Static_Name() 
    {return STRING_UTILITIES::string_sprintf("LEVELSET_IMPLICIT_OBJECT<T,VECTOR<T,%d> >",TV::dimension);}

    static std::string Static_Extension()
    {return TV::dimension==2?"phi2d":"phi";}

    virtual std::string Name() const PHYSBAM_OVERRIDE {return Static_Name();}
    virtual std::string Extension() const PHYSBAM_OVERRIDE {return Static_Extension();}

    static IMPLICIT_OBJECT_REMOVED<TV>* Create()
    {PHYSBAM_FATAL_ERROR();return 0;}

//#####################################################################
};
}
#endif

