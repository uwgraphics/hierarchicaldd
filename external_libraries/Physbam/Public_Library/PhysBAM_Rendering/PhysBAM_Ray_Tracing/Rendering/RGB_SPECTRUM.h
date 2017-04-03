//#####################################################################
// Copyright 2013, Ed Quigley
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RGB_SPECTRUM
//#####################################################################
#ifndef __RGB_SPECTRUM__
#define __RGB_SPECTRUM__

#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering/CIE_XYZ.h>
#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Matrices/MATRIX.h>
#include <PhysBAM_Tools/Matrices/MATRIX_3X3.h>
#include <PhysBAM_Tools/Matrices/MATRIX_MXN.h>
#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>
namespace PhysBAM{

template<class T>
class RGB_SPECTRUM
{
    typedef VECTOR<T,3> TV;
public:
    T w1,w2,w3;
    T w1_inverse,w2_inverse,w3_inverse;
    const T wmin,wmax;
    T lc1,lc2,lc3;
    ARRAY<T> x;
    TV color;
    CIE_XYZ<T> cie;
    MATRIX<T,3> coefficients_inverse;
    const int max_samples;
    int samples;
    T (RGB_SPECTRUM<T>::*basis_function)(int,T) const;

    RGB_SPECTRUM(int samples_);
//#####################################################################
    ARRAY<T> RGB_To_Spectrum(const VECTOR<T,3>& rgb) const;
    VECTOR<T,3> Spectrum_To_RGB(const ARRAY<T>& spectrum) const;
    T Local_Saturation(int i,int j);
    void Compute_Weights();
    T F_Gaussian(int i,T lambda) const;
    T F_Fourier(int i,T lambda) const;
    T F_Delta(int i,T lambda) const;
    T F_Box(int i,T lambda) const;
//#####################################################################
};
}
#endif
