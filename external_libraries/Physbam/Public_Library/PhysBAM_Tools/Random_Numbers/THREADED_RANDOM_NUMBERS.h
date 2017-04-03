//#####################################################################
// Copyright 2013, Bo Zhu
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RANDOM_NUMBERS
//#####################################################################
#ifndef __THREADED_RANDOM_NUMBERS__
#define __THREADED_RANDOM_NUMBERS__
#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Random_Numbers/RANDOM_NUMBERS.h>
#include <time.h>
namespace PhysBAM{
template<class T,class GENERATOR=MT19937<T> >
class THREADED_RANDOM_NUMBERS:public NONCOPYABLE
{
public:
    int thread_number;
    ARRAY<RANDOM_NUMBERS<T,GENERATOR>* > random_numbers;
public:
    T Get_Number(int tid){return random_numbers(tid)->Get_Number();}

    explicit THREADED_RANDOM_NUMBERS(const int thread_number_input=1):thread_number(0) {Initialize_Threads(thread_number_input);}
    virtual ~THREADED_RANDOM_NUMBERS(){}

    void Initialize_Threads(int thread_number_input,int seed=time(0))
    {Clear_Threads();thread_number=thread_number_input;random_numbers.Resize(thread_number);for(int i=1;i<=thread_number;i++)random_numbers(i)=new RANDOM_NUMBERS<T,GENERATOR>(seed);}

    void Clear_Threads()
    {for(int i=1;i<=thread_number;i++)delete random_numbers(i);thread_number=0;}

    void Set_Seeds(const ARRAY<unsigned int> seeds_input)
    {for(int i=1;i<=thread_number;i++)Set_Seed(seeds_input(i),i);}

    template<class T2,class T_ARRAY,class ID> void Fill_Uniform(ARRAY_BASE<T2,T_ARRAY,ID>& array,const T a,const T b,const int tid)
    {random_numbers(tid)->Fill_Uniform(array,a,b);}

    //#####################################################################
    void Set_Seed(const unsigned int seed_input,const int tid){random_numbers(tid)->Set_Seed(seed_input);}
    int Get_Uniform_Integer(const int a,const int b,const int tid){return random_numbers(tid)->Get_Uniform_Integer(a,b);}
    T Get_Uniform_Number(const T a,const T b,const int tid){return random_numbers(tid)->Get_Uniform_Number(a,b);}
    template<int d> VECTOR<T,d> Get_Uniform_Vector(const VECTOR<T,d>& v0,const VECTOR<T,d>& v1,const int tid){return random_numbers(tid)->template Get_Uniform_Vector<d>(v0,v1);}
    template<int d> VECTOR<T,d> Get_Uniform_Vector(const T a,const T b,const int tid){return random_numbers(tid)->template Get_Uniform_Vector<d>(a,b);}
    template<class TV> TV Get_Uniform_Vector(const RANGE<TV>& box,const int tid){return random_numbers(tid)->template Get_Uniform_Vector<TV>(box);}
    void Fill_Uniform(T& x,const T a,const T b,const int tid){random_numbers(tid)->Fill_Uniform(x,a,b);}
    template<class T_VECTOR> void Fill_Uniform(VECTOR_BASE<T,T_VECTOR>& v,const T a,const T b,const int tid){random_numbers(tid)->Fill_Uniform(v,a,b);}
    template<class T_MATRIX> void Fill_Uniform(MATRIX_BASE<T,T_MATRIX>& m,const T a,const T b,const int tid){random_numbers(tid)->Fill_Uniform(m,a,b);}
    template<int d> void Fill_Uniform(DIAGONAL_MATRIX<T,d>& m,const T a,const T b,const int tid){random_numbers(tid)->template Fill_Uniform<d>(m,a,b);}
    template<int d> void Fill_Uniform(SYMMETRIC_MATRIX<T,d>& m,const T a,const T b,const int tid){random_numbers(tid)->template Fill_Uniform<d>(m,a,b);}
    template<int d> void Fill_Uniform(UPPER_TRIANGULAR_MATRIX<T,d>& m,const T a,const T b,const int tid){random_numbers(tid)->template Fill_Uniform<d>(m,a,b);}
    template<class TV> void Fill_Uniform(TWIST<TV>& m,const T a,const T b,const int tid){random_numbers(tid)->template Fill_Uniform<TV>(m,a,b);}
    T Get_Gaussian(const int tid){return random_numbers(tid)->Get_Gaussian(tid);}
    template<class TV> TV Get_Vector_In_Unit_Sphere(const int tid){return random_numbers(tid)->template Get_Vector_In_Unit_Sphere<TV>();}
    template<class TV> TV Get_Direction(const int tid){return random_numbers(tid)->template Get_Direction<TV>();}
    template<class TV> ROTATION<TV> Get_Rotation(const int tid){return random_numbers(tid)->template Get_Rotation<TV>();}
    template<class TV> FRAME<TV> Get_Frame(const TV& v0,const TV& v1,const int tid){return random_numbers(tid)->template Get_Frame<TV>(v0,v1);}
    template<class TV> TWIST<TV> Get_Twist(const T& a,const int tid){return random_numbers(tid)->template Get_Twist<TV>(a);}
    //#####################################################################
};
}
#endif