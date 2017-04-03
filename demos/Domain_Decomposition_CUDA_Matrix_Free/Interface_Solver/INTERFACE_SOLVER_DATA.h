//#####################################################################
//Removed For Anonymity: Copyright Authors of Submission pap173s1
// This file is part of PhysBAM whose distribution is governed by the license 
//Removed For Anonymity: Copyright Authors of Submission pap173s1
//#####################################################################
#ifndef __Interface_Solver_Data_h__
#define __Interface_Solver_Data_h__

#include <Eigen/SparseCore>
#include <Eigen/SparseCholesky>
#include <array>
#include <vector>
#include "mkl_pardiso.h"
#include "mkl_types.h"

//#define USING_EIGEN

template<typename T,int d>
class INTERFACE_SOLVER_DATA{
    //In Eigen 3.0 dense matrices are always indexed by int of size of pointers
    typedef Eigen::Matrix<T,Eigen::Dynamic,1> T_EIGEN_VECTOR;
public:

    std::vector<std::vector<std::array<void*,64> > > pardiso_pt_i; // One 64-element array of void*'s per hierarchy level, to pass to pardiso_64 (for subdomains)
    std::vector<std::vector<std::array<MKL_INT,64> > > pardiso_iparm_i;

    std::vector<std::array<void*,64> > pardiso_pt_r; // One 64-element array of void*'s per hierarchy level, to pass to pardiso_64 (for interface)
    std::vector<std::array<MKL_INT,64> > pardiso_iparm_r;

    int maxfct_i;

    std::vector<std::vector<Eigen::SparseMatrix<T,Eigen::ColMajor,int> > > Aii;
    std::vector<std::vector<T* > >           Aii_a;  // Aii(level)(index) is A_{ii}^{\ast 2^{level-1}h}
    std::vector<std::vector<MKL_INT > >      Aii_n;  // Per-subdomain, per-level, adaptively coarsened (and locally indexed) interior operator matrix    
    std::vector<std::vector<MKL_INT* > >     Aii_ia; // Aii_a contains the matrix data in CSR format (varies per subdomain)
    std::vector<std::vector<MKL_INT* > >     Aii_ja; // Aii_n is the matrix size, Aii_ia is row offsets, Aii_ja column indices (do not vary per subdomain)

    mutable std::vector<std::vector<T_EIGEN_VECTOR > > Aii_vx;//vectors used in subdomain solve 
    mutable std::vector<std::vector<T_EIGEN_VECTOR > > Aii_vb;//vectors used in subdomain solve 

    std::vector<Eigen::SparseMatrix<T,Eigen::ColMajor,int> > Arr;  // Arr(level) is A_{ii}^{\ast 2^{level-1}h}
    std::vector<T*> Arr_a;
    std::vector<MKL_INT>      Arr_n;
    std::vector<MKL_INT*>     Arr_ia;
    std::vector<MKL_INT*>     Arr_ja;

    std::vector<std::vector<Eigen::SparseMatrix<T,Eigen::RowMajor,int> > > Air;                // Air, \gamma is indexed locally
    std::vector<Eigen::SparseMatrix<T,Eigen::RowMajor,int> > Prr;
    //####################################################################################3
    void Resize(const int levels,const int number_of_subdomains){
        Aii.resize(levels);
        Aii_a.resize(levels);
        Aii_n.resize(levels);
        Aii_ia.resize(levels);
        Aii_ja.resize(levels);
        Aii_vx.resize(levels);
        Aii_vb.resize(levels);
        Air.resize(levels);

        Arr.resize(levels);
        Arr_a.resize(levels);
        Arr_n.resize(levels);
        Arr_ia.resize(levels);
        Arr_ja.resize(levels);

        Prr.resize(levels-1);
    
        pardiso_pt_i.resize(levels);
        pardiso_pt_r.resize(levels);
        pardiso_iparm_i.resize(levels);
        pardiso_iparm_r.resize(levels);

        for(int level=0;level<levels;++level){
            Aii[level].resize(number_of_subdomains);
            Aii_n[level].resize(number_of_subdomains);
            Aii_a[level].resize(number_of_subdomains);
            Aii_ia[level].resize(number_of_subdomains);
            Aii_ja[level].resize(number_of_subdomains);
            Air[level].resize(number_of_subdomains);
            Aii_vx[level].resize(number_of_subdomains);
            Aii_vb[level].resize(number_of_subdomains);

            pardiso_pt_i[level].resize(number_of_subdomains);
            pardiso_iparm_i[level].resize(number_of_subdomains);
            
            Arr_a[level]=NULL;
            Arr_ia[level]=NULL;
            Arr_ja[level]=NULL;
            
            for(int i=0;i<number_of_subdomains;++i){
                Aii_a[level][i]=NULL;            
                Aii_ia[level][i]=NULL;
                Aii_ja[level][i]=NULL;}}
    }
    int Get_Number_Of_Levels()const{
        return Aii_a.size();
    }
    int Get_Number_Of_Subdomains()const{
        return Aii_a[0].size();
    }
};
#endif
