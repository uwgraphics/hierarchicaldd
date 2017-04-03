//#####################################################################
//Removed For Anonymity: Copyright Authors of Submission pap173s1
//Removed For Anonymity: Copyright Authors of Submission pap173s1
//#####################################################################
#ifndef __SPGRID_DOMAIN_DECOMPOSITION_DATA_H__
#define __SPGRID_DOMAIN_DECOMPOSITION_DATA_H__

namespace SPGrid{

template<class T>
struct SPGRID_DOMAIN_DECOMPOSITION_DATA
{
    union{
        unsigned flags;
        T ch3;
    };
    T ch0;  
    T ch1;
    T ch2;
};

}
#endif
