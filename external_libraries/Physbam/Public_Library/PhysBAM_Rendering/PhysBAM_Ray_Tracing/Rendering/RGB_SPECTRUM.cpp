//#####################################################################
// Copyright 2013, Ed Quigley
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering/RGB_SPECTRUM.h>
using namespace PhysBAM;

template<class T>
RGB_SPECTRUM<T>::RGB_SPECTRUM(int samples_)
    : wmin((T)10e-9),wmax((T)150e-9),lc1((T)680e-9),lc2((T)550e-9),lc3((T)420e-9),color((T)0,(T)0,(T)0),/*k(680),*/max_samples(81),samples(samples_)
{
    basis_function=&RGB_SPECTRUM<T>::F_Fourier;
    Compute_Weights();
    ARRAY<ARRAY<T> > t;
    for(int i=1;i<=3;i++){
        t.Append(ARRAY<T>());
        ARRAY<T,VECTOR<int,1> >& color_match=(i==1?cie.X_spectrum:(i==2?cie.Y_spectrum:cie.Z_spectrum));
        for(int j=1;j<=3;j++){
            t(i).Append((T)0);
            for(int m=1;m<=samples;m++){
                T wavelength=cie.starting_wavelength+(T)(m-1)*(cie.ending_wavelength-cie.starting_wavelength)/(samples-1);
                int index=(int)(m*max_samples/samples);
                t(i)(j)+=color_match(index)*(*this.*basis_function)(j,wavelength);
            }
        }
    }
    for(int i=1;i<=3;i++){
        VECTOR<T,3> xyz(t(1)(i),t(2)(i),t(3)(i));
        VECTOR<T,3> rgb=cie.XYZ_To_RGB(xyz);
        t(1)(i)=rgb(1);t(2)(i)=rgb(2);t(3)(i)=rgb(3);
    }
    MATRIX<T,3> coefficients(t(1)(1),t(2)(1),t(3)(1),t(1)(2),t(2)(2),t(3)(2),t(1)(3),t(2)(3),t(3)(3));
    coefficients_inverse=coefficients.Inverse();
}

template<class T>
ARRAY<T> RGB_SPECTRUM<T>::RGB_To_Spectrum(const VECTOR<T,3>& rgb) const
{
    VECTOR<T,3> xyz=cie.RGB_To_XYZ(rgb);
    MATRIX_MXN<T> colors(1,3);
    colors.Set_Row(1,xyz);
    MATRIX_MXN<T> result=coefficients_inverse*colors;
    TV x(result(1,1),result(2,1),result(3,1));
    ARRAY<T> spectrum;
    for(int i=1;i<=samples;i++){
        spectrum.Append((T)0);
        T wavelength=cie.starting_wavelength+(T)(i-1)*(cie.ending_wavelength-cie.starting_wavelength)/(samples-1);
        for(int j=1;j<=3;j++){
            spectrum(i)+=x(j)*(*this.*basis_function)(j,wavelength);
        }
    }
    return spectrum;
}

template<class T>
VECTOR<T,3> RGB_SPECTRUM<T>::Spectrum_To_RGB(const ARRAY<T>& spectrum) const
{
    TV xyz;
    for(int i=1;i<=3;i++){
        xyz(i)=(T)0;
        const ARRAY<T,VECTOR<int,1> >& color_match=(i==1?cie.X_spectrum:(i==2?cie.Y_spectrum:cie.Z_spectrum));
        for(int j=1;j<=samples;j++){
            int index=(int)(j*max_samples/samples);
            xyz(i)+=spectrum(j)*color_match(index);
        }
    }
    return cie.XYZ_To_RGB(xyz);
}

template<class T>
T RGB_SPECTRUM<T>::Local_Saturation(int i,int j)
{
    if(color(i)+color(j)==(T)0)return (T)0;
    return abs(color(i)-color(j))/(color(i)+color(j));
}

template<class T>
void RGB_SPECTRUM<T>::Compute_Weights()
{
    w1=Local_Saturation(1,2)*wmin+((T)1-Local_Saturation(1,2))*wmax;
    w2=Local_Saturation(3,2)*wmin+((T)1-Local_Saturation(3,2))*wmax;
    w3=min(w1,w2);
    w1_inverse=(w1==(T)0?(T)0:(T)1./w1);
    w2_inverse=(w2==(T)0?(T)0:(T)1./w2);
    w3_inverse=(w3==(T)0?(T)0:(T)1./w3);
}

template<class T>
T RGB_SPECTRUM<T>::F_Gaussian(int i,T lambda) const
{
    T w_inverse=(i==1?w1_inverse:(i==2?w2_inverse:w3_inverse));
    T lc=(i==1?lc1:(i==2?lc2:lc3));
    return (T)exp(-1*log(2)*(2*(lambda-lc)*w_inverse)*(2*(lambda-lc)*w_inverse));
}

template<class T>
T RGB_SPECTRUM<T>::F_Fourier(int i,T lambda) const
{
    static const T lambda_min=(T)380e-9,lambda_max=(T)780e-9;
    static const T lambda_max_minus_min_inverse=(T)1/(lambda_max-lambda_min);
    T return_value;
    switch(i){
    case 1:return_value=1;break;
    case 2:return_value=(T)cos(two_pi*(lambda-lambda_min)*lambda_max_minus_min_inverse);break;
    case 3:default:return_value=(T)sin(two_pi*(lambda-lambda_min)*lambda_max_minus_min_inverse);break;
    }
    return return_value;
}

template<class T>
T RGB_SPECTRUM<T>::F_Delta(int i,T lambda) const
{
    T l=(i==1?(T)590e-9:(i==2?(T)560e-9:(T)440e-9));
    const T epsilon=(T)0.5e-9;
    T return_value;
    if(abs(l-lambda)<epsilon)return_value=(T)1;
    else return_value=(T)0;
    return return_value;
}

template<class T>
T RGB_SPECTRUM<T>::F_Box(int i,T lambda) const
{
    ARRAY<T> l;
    for(int j=1;j<=4;j++)l.Append((T)(300e-9+j*100e-9));
    T return_value;
    if(lambda>l(i)&&lambda<l(i+1))return_value=(T)1;
    else return_value=(T)0;
    return return_value;
}

template class RGB_SPECTRUM<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class RGB_SPECTRUM<double>;
#endif
