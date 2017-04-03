//#####################################################################
//Removed For Anonymity: Copyright Authors of Submission pap173s1
// This file is part of PhysBAM whose distribution is governed by the license 
//Removed For Anonymity: Copyright Authors of Submission pap173s1
//#####################################################################
#ifndef __STRIDE_GENERATOR_H__
#define __STRIDE_GENERATOR_H__

#include <Common_Tools/Math_Tools/RANGE_ITERATOR.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Vectors/VECTOR_ND.h>
#include <PhysBAM_Tools/Data_Structures/STACK.h>
#include <PhysBAM_Tools/Data_Structures/PAIR.h>
#include <limits>
using namespace PhysBAM;

template<int d>
class STRIDE_GENERATOR{
    typedef VECTOR<int,d> T_INDEX;
    typedef ARRAY<int,VECTOR<int,1> > ARRAY_1D;

public:
    static int Generate_Adaptive_Stride(const RANGE<T_INDEX>& range,VECTOR<ARRAY_1D,d>& q){
        int adaptation_level=std::numeric_limits<int>::max();
        T_INDEX counts = range.max_corner - range.min_corner;
        for(int axis = 1;axis <= d;++axis){
            int axis_adaptation_level=0;
            q(axis).Resize(RANGE<VECTOR<int,1> >(range.min_corner(axis),range.max_corner(axis)));
            ARRAY<PAIR<int,int> > adaptation;
            adaptation.Append(PAIR<int,int>(2,1));
            int current_size=1,domain_size=2;
            while(domain_size<counts(axis)){int quotient=min(2,(counts(axis)-domain_size)/current_size);
                adaptation.Append(PAIR<int,int>(quotient,current_size));domain_size+=quotient*current_size;
                if((counts(axis)-domain_size)%(2*current_size)==0) {current_size*=2;axis_adaptation_level++;}}

            int min_corner = range.min_corner(axis);
            STACK<int> s;int cnt=1;
            for(int i=1;i<=adaptation.m;++i){for(int k=1;k<=adaptation(i).y;++k) {q(axis)(min_corner + cnt - 1)=adaptation(i).y;++cnt;};
                if(adaptation(i).x==1) {q(axis)(min_corner + cnt - 1) = adaptation(i).y;++cnt;}
                else{if(i==adaptation.m) {q(axis)(min_corner + cnt - 1) = adaptation(i).y;++cnt;}
                    s.Push(adaptation(i).y);}}
            while(!s.Empty()){int current_size=s.Pop();
                for(int i=1;i<=current_size;++i) {q(axis)(min_corner + cnt - 1) = current_size;++cnt;}}
            int min = q(axis).Domain_Indices().min_corner(1);
            int max = q(axis).Domain_Indices().max_corner(1);
            
            //LOG::cout << "q[" << min << "," << max  <<"]on axis " << axis << ":[";
            //for(int i = min;i < max;++i) {LOG::cout << q(axis)(i) << " ";}
            //LOG::cout << q(axis)(max) << "]" << std::endl;

            if(axis_adaptation_level<adaptation_level) adaptation_level = axis_adaptation_level;
        }
        return adaptation_level;
    }
};
#endif
