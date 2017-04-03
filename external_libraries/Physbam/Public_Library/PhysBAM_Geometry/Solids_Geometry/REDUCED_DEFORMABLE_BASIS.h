//#####################################################################
// Copyright 2012-2013, Rahul Sheth.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class REDUCED_DEFORMABLE_BASIS
//#####################################################################
#ifndef __REDUCED_DEFORMABLE_BASIS__
#define __REDUCED_DEFORMABLE_BASIS__

#include <PhysBAM_Tools/Matrices/MATRIX_MXN.h>

namespace PhysBAM{

template<class TV>
class REDUCED_DEFORMABLE_BASIS
{
    typedef typename TV::SCALAR T;
  public:
    MATRIX_MXN<T> reduced_basis;

    REDUCED_DEFORMABLE_BASIS() {};
    REDUCED_DEFORMABLE_BASIS(const MATRIX_MXN<T>& m_in):reduced_basis(m_in) {};
    ~REDUCED_DEFORMABLE_BASIS() {};

    virtual const MATRIX_MXN<T>& Get_Reduced_Basis() const{return reduced_basis;};
    virtual void Set_Reduced_Basis(const MATRIX_MXN<T>& m_in){reduced_basis=m_in;};

    template<class T_VECTOR> const LAZY_MATRIX_VECTOR_PRODUCT<MATRIX_MXN<T>,T_VECTOR> operator*(const T_VECTOR& vector_in) const
    {return reduced_basis.Lazy_Multiply(vector_in);};

    //can't return a reference, so don't use this for assignment
    //TODO: use a vector view type class to make this cleaner
    TV operator()(const int particle_index,const int column_index) const {
        TV temp_vec;
        for(int i=1;i<=TV::dimension;i++) {
            temp_vec(i)=reduced_basis(particle_index*TV::dimension-(TV::dimension-i),column_index);}
        return temp_vec;}
    void Set_Particle_Vector(const TV& vector_in,const int particle_index,const int column_index) {
        for(int i=1;i<=TV::dimension;i++) {
            reduced_basis(particle_index*TV::dimension-(TV::dimension-i),column_index)=vector_in(i);}}
//##################################################################### 
};
}
#endif
