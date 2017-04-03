//#####################################################################
//Removed For Anonymity: Copyright Authors of Submission pap173s1
// This file is part of PhysBAM whose distribution is governed by the license 
//Removed For Anonymity: Copyright Authors of Submission pap173s1
//#####################################################################
#include "SIGMA_SMOOTHER.h"
#include <Eigen/SparseCore>
#include <Eigen/SparseCholesky>
#include <chrono>
#include <fstream>
#include <sstream>
using namespace PhysBAM;

namespace {
    class Timer_Helper
    {
    public:
        Timer_Helper() : beg_(clock_::now()) {}
        void reset() { beg_ = clock_::now(); }
        double elapsed() const { 
            return std::chrono::duration_cast<second_>
                (clock_::now() - beg_).count(); }

    private:
        typedef std::chrono::high_resolution_clock clock_;
        typedef std::chrono::duration<double, std::ratio<1> > second_;
        std::chrono::time_point<clock_> beg_;
    };
};

//#####################################################################
// Function Check_Symmetry (DEBUG)
//#####################################################################
template<typename T> void
Check_Symmetry(const Eigen::SparseMatrix<T,Eigen::ColMajor,int>& matrix,T tolorence=1e-5){
    int counter=0;
    for(int k=0;k<matrix.outerSize();++k){
        for(typename Eigen::SparseMatrix<T,Eigen::ColMajor,int>::InnerIterator itr(matrix,k);itr;++itr){
            PHYSBAM_ASSERT(itr.value()!=0);
            if(fabs(itr.value()-matrix.coeff(itr.col(),itr.row()))>tolorence)
                std::cout<<matrix<<std::endl;;
            PHYSBAM_ASSERT(fabs(itr.value()-matrix.coeff(itr.col(),itr.row()))<=tolorence);
            ++counter;}}
    PHYSBAM_ASSERT(counter==matrix.nonZeros());
}
//#####################################################################
// Function Eigen_To_Pardiso
//#####################################################################
template<typename T> void
Eigen_To_Pardiso(const Eigen::SparseMatrix<T,Eigen::ColMajor,int>& matrix,MKL_INT& n,T*& a,MKL_INT*& ia,MKL_INT*& ja,T epsilon=1e-5){
    //We assumes that the matrix is SPD
    //Pardiso wants the (U)pper triangle of a Row major matrix,
    //So we give it the (L)ower triangle of a Col major matrix.
    n=matrix.rows();
    PHYSBAM_ASSERT(matrix.cols()==n);
    PHYSBAM_ASSERT(matrix.isCompressed());
    if( a!=NULL) delete[] a;
    if(ia!=NULL) delete[]ia;
    if(ja!=NULL) delete[]ja;
    
    MKL_INT nnz=matrix.nonZeros();
    
    ia=new MKL_INT[n+1];
    
    MKL_INT counter=0;
    int free_dof=0;
    int k_start=-1;
    for(int k=0;k<matrix.outerSize();++k){
        for(typename Eigen::SparseMatrix<T,Eigen::ColMajor,int>::InnerIterator itr(matrix,k);itr;++itr){
            if(k!=k_start){
                //if(k!=k_start+1) std::cout<<k<<" vs "<<k_start<<std::endl;
                for(int i=k_start;i<k-1;++i) {++counter;++free_dof;}//Those are the DOF has zero values every where
                k_start=k;}
            if(itr.row()>=itr.col()){++counter;}
            if(itr.row()==itr.col()){PHYSBAM_ASSERT(itr.value()!=0);}}}
    for(int i=k_start;i<matrix.outerSize()-1;++i) {++counter;++free_dof;}//Those are the DOF has zero values every where

    MKL_INT nnz_lower=(nnz+n+free_dof)/2;

    ja=new MKL_INT[counter];
    a =new T[counter];

    counter=0;
    k_start=-1;
    for(int k=0;k<matrix.outerSize();++k){
        for(typename Eigen::SparseMatrix<T,Eigen::ColMajor,int>::InnerIterator itr(matrix,k);itr;++itr){
            if(k!=k_start){
                for(int i=k_start;i<k-1;++i){
                    //Those are the DOF has zero values every where
                    //std::cout<<"Errr....Here?"<<std::endl;
                    ia[i+1]=counter;
                    ja[counter]=i+1;
                    a[counter] =epsilon;
                    ++counter;}
                ia[k]=counter;
                k_start=k;}
            if(itr.row()>=itr.col()){
                ja[counter]=itr.row();
                a[counter] =itr.value();
                ++counter;}
            if(itr.row()==itr.col()){
                a[counter-1]+=epsilon;}}}
    for(int i=k_start;i<matrix.outerSize()-1;++i){
        ia[i+1]=counter;
        ja[counter]=i+1;
        a[counter] =epsilon;
        ++counter;}
    
    ia[n]=counter;
    PHYSBAM_ASSERT(counter==nnz_lower);
}
//#####################################################################
// Function Symbolic_Factorization
//#####################################################################
template<typename T,int d> void SIGMA_SMOOTHER<T,d>::
Symbolic_Factorization(){
    LOG::SCOPE scope("Smoother Symbolic Factorization");
    //First... we need to prob out the lower triangle of the matrices....
    for(int level=0;level<solver_data.Get_Number_Of_Levels();++level){
        std::cout<<"@LEVEL: "<<level<<std::endl;
        //############################################
        // Interface
        //############################################
        //Check_Symmetry<T>(solver_data.Arr[level]);
        Eigen_To_Pardiso<T>(solver_data.Arr[level],
                            solver_data.Arr_n[level],solver_data.Arr_a[level],
                            solver_data.Arr_ia[level],solver_data.Arr_ja[level]);
        
        for(int i=0;i<64;i++){solver_data.pardiso_iparm_r[level][i]=0;}
        solver_data.pardiso_iparm_r[level][0] = 1;         /* No solver default */
        solver_data.pardiso_iparm_r[level][1] = 3;         /* Fill-in reordering from METIS */
        solver_data.pardiso_iparm_r[level][3] = 0;         /* No iterative-direct algorithm */
        solver_data.pardiso_iparm_r[level][4] = 0;         /* No user fill-in reducing permutation */
        solver_data.pardiso_iparm_r[level][5] = 0;         /* Write solution into x */
        solver_data.pardiso_iparm_r[level][6] = 0;         /* Not in use */
        solver_data.pardiso_iparm_r[level][7] = 0;         /* Max numbers of iterative refinement steps */
        solver_data.pardiso_iparm_r[level][8] = 0;         /* Not in use */
        solver_data.pardiso_iparm_r[level][9] = 13;        /* Perturb the pivot elements with 1E-13 */
        solver_data.pardiso_iparm_r[level][10] = 1;        /* Use nonsymmetric permutation and scaling MPS */
        solver_data.pardiso_iparm_r[level][11] = 0;        /* Not in use */
        solver_data.pardiso_iparm_r[level][12] = 0;        /* Maximum weighted matching algorithm is switched-off (default for symmetric). Try iparm[12] = 1 in case of inappropriate accuracy */
        solver_data.pardiso_iparm_r[level][13] = 0;        /* Output: Number of perturbed pivots */
        solver_data.pardiso_iparm_r[level][14] = 0;        /* Not in use */
        solver_data.pardiso_iparm_r[level][15] = 0;        /* Not in use */
        solver_data.pardiso_iparm_r[level][16] = 0;        /* Not in use */
        solver_data.pardiso_iparm_r[level][17] = -1;       /* Output: Number of nonzeros in the factor LU */
        solver_data.pardiso_iparm_r[level][18] = -1;       /* Output: Mflops for LU factorization */
        solver_data.pardiso_iparm_r[level][19] = 0;        /* Output: Numbers of CG Iterations */
        solver_data.pardiso_iparm_r[level][23] = 1;        /* Two-level factorization*/
        solver_data.pardiso_iparm_r[level][26] = 1;        /* Check matrix for errors */
        solver_data.pardiso_iparm_r[level][27] = std::is_same<T,float>::value;        /* float or double precision */   
        solver_data.pardiso_iparm_r[level][34] = 1;        /* Zero based indexing */   
        for(int i=0;i<64;i++){solver_data.pardiso_pt_r[level][i]=0;}
        MKL_INT maxfct = 1;           /* Maximum number of numerical factorizations. */
        MKL_INT mnum = 1;             /* Which factorization to use. */
        MKL_INT msglvl = 0;           /* Print statistical information in file */
        MKL_INT error = 0;            /* Initialize error flag */
        MKL_INT nrhs=1;
        MKL_INT mtype=2;
        MKL_INT idum;
        T ddum;        
        /* -------------------------------------------------------------------- */
        /* .. Reordering and Symbolic Factorization. This step also allocates */
        /* all memory that is necessary for the factorization. */
        /* -------------------------------------------------------------------- */
        MKL_INT phase = 11;
        PARDISO(&solver_data.pardiso_pt_r[level][0],&maxfct,&mnum,&mtype,&phase,
                &solver_data.Arr_n[level],solver_data.Arr_a[level],solver_data.Arr_ia[level],solver_data.Arr_ja[level],
                &idum,&nrhs,&solver_data.pardiso_iparm_r[level][0],&msglvl,&ddum,&ddum, &error);
        if(error!=0){
            printf("\nERROR during symbolic factorization: %d",error);
            exit(1);}
        printf("\nNumber of nonzeros in factors = %d",solver_data.pardiso_iparm_r[level][17]);
        printf("\nNumber of factorization MFLOPS = %d",solver_data.pardiso_iparm_r[level][18]);
        
        //############################################
        // Interior
        //############################################
        for(int subdomain_index=0;subdomain_index<solver_data.Get_Number_Of_Subdomains();++subdomain_index){
            std::cout<<"@Subdomain: "<<subdomain_index<<std::endl;
            for(int i=0;i<64;i++){solver_data.pardiso_iparm_i[level][subdomain_index][i]=0;}
            solver_data.pardiso_iparm_i[level][subdomain_index][0] = 1;         /* No solver default */
            solver_data.pardiso_iparm_i[level][subdomain_index][1] = 3;         /* Fill-in reordering from METIS */
            solver_data.pardiso_iparm_i[level][subdomain_index][3] = 0;         /* No iterative-direct algorithm */
            solver_data.pardiso_iparm_i[level][subdomain_index][4] = 0;         /* No user fill-in reducing permutation */
            solver_data.pardiso_iparm_i[level][subdomain_index][5] = 0;         /* Write solution into x */
            solver_data.pardiso_iparm_i[level][subdomain_index][6] = 0;         /* Not in use */
            solver_data.pardiso_iparm_i[level][subdomain_index][7] = 0;         /* Max numbers of iterative refinement steps */
            solver_data.pardiso_iparm_i[level][subdomain_index][8] = 0;         /* Not in use */
            solver_data.pardiso_iparm_i[level][subdomain_index][9] = 13;        /* Perturb the pivot elements with 1E-13 */
            solver_data.pardiso_iparm_i[level][subdomain_index][10] = 1;        /* Use nonsymmetric permutation and scaling MPS */
            solver_data.pardiso_iparm_i[level][subdomain_index][11] = 0;        /* Not in use */
            solver_data.pardiso_iparm_i[level][subdomain_index][12] = 0;        /* Maximum weighted matching algorithm is switched-off (default for symmetric). Try iparm[12] = 1 in case of inappropriate accuracy */
            solver_data.pardiso_iparm_i[level][subdomain_index][13] = 0;        /* Output: Number of perturbed pivots */
            solver_data.pardiso_iparm_i[level][subdomain_index][14] = 0;        /* Not in use */
            solver_data.pardiso_iparm_i[level][subdomain_index][15] = 0;        /* Not in use */
            solver_data.pardiso_iparm_i[level][subdomain_index][16] = 0;        /* Not in use */
            solver_data.pardiso_iparm_i[level][subdomain_index][17] = -1;       /* Output: Number of nonzeros in the factor LU */
            solver_data.pardiso_iparm_i[level][subdomain_index][18] = -1;       /* Output: Mflops for LU factorization */
            solver_data.pardiso_iparm_i[level][subdomain_index][19] = 0;        /* Output: Numbers of CG Iterations */
            solver_data.pardiso_iparm_i[level][subdomain_index][23] = 1;        /* Two-level factorization*/
            solver_data.pardiso_iparm_i[level][subdomain_index][26] = 1;        /* Check matrix for errors */
            solver_data.pardiso_iparm_i[level][subdomain_index][27] = std::is_same<T,float>::value;        /* float or double precision */   
            solver_data.pardiso_iparm_i[level][subdomain_index][34] = 1;        /* Zero based indexing */   
            for(int i=0;i<64;i++){solver_data.pardiso_pt_i[level][subdomain_index][i]=0;}

            //We just need the first one for factorization
            //Check_Symmetry<T>(solver_data.Aii[level][subdomain_index]);
            Eigen_To_Pardiso<T>(solver_data.Aii[level][subdomain_index],solver_data.Aii_n[level][subdomain_index],
                                      solver_data.Aii_a[level][subdomain_index],solver_data.Aii_ia[level][subdomain_index],
                                      solver_data.Aii_ja[level][subdomain_index]);
            solver_data.Aii[level][subdomain_index].resize(1,1);
            solver_data.Aii[level][subdomain_index].data().squeeze(); 
            solver_data.Aii_vx[level][subdomain_index].resize(solver_data.Aii_n[level][subdomain_index]);
            solver_data.Aii_vb[level][subdomain_index].resize(solver_data.Aii_n[level][subdomain_index]);

            maxfct = 1;                              /* Maximum number of numerical factorizations. */
            mnum = 1;                                /* Which factorization to use. */
            /* -------------------------------------------------------------------- */
            /* .. Reordering and Symbolic Factorization. This step also allocates */
            /* all memory that is necessary for the factorization. */
            /* -------------------------------------------------------------------- */
            phase = 11;
            //std::cout<<"1."<<std::endl;
            PARDISO(&solver_data.pardiso_pt_i[level][subdomain_index][0],&maxfct,&mnum,&mtype,&phase,
                    &solver_data.Aii_n[level][subdomain_index],solver_data.Aii_a[level][subdomain_index],
                    solver_data.Aii_ia[level][subdomain_index],solver_data.Aii_ja[level][subdomain_index],
                    &idum,&nrhs,&solver_data.pardiso_iparm_i[level][subdomain_index][0],&msglvl,&ddum,&ddum,&error);
            if(error!=0){
                printf("\nERROR during symbolic factorization: %d",error);
                exit(1);}
            //std::cout<<"2."<<std::endl;
        }
        printf("\nNumber of nonzeros in factors = %d",solver_data.pardiso_iparm_i[level][17]);
        printf("\nNumber of factorization MFLOPS = %d",solver_data.pardiso_iparm_i[level][18]);
    }
    symbolic_factorized=true;
}
//#####################################################################
// Function Numerical_Factorization
//#####################################################################
template<typename T,int d> void SIGMA_SMOOTHER<T,d>::
Numerical_Factorization(){
    LOG::SCOPE scope("Smoother Numerical Factorization");
    PHYSBAM_ASSERT(symbolic_factorized);
    //First... we need to prob out the lower triangle of the matrices....
    for(int level=0;level<solver_data.Get_Number_Of_Levels();++level){
        Timer_Helper timer;
        //############################################
        // Interface
        //############################################
        MKL_INT maxfct = 1;           /* Maximum number of numerical factorizations. */
        MKL_INT mnum = 1;             /* Which factorization to use. */
        MKL_INT msglvl = 0;           /* Print statistical information in file */
        MKL_INT error = 0;            /* Initialize error flag */
        MKL_INT nrhs=1;
        MKL_INT mtype=2;
        MKL_INT idum;
        T ddum;                
        /* -------------------------------------------------------------------- */
        /* .. Numerical factorization. */
        /* -------------------------------------------------------------------- */
        MKL_INT phase=22;
        PARDISO(&solver_data.pardiso_pt_r[level][0],&maxfct,&mnum,&mtype,&phase,
                &solver_data.Arr_n[level],solver_data.Arr_a[level],solver_data.Arr_ia[level],solver_data.Arr_ja[level],
                &idum,&nrhs,&solver_data.pardiso_iparm_r[level][0],&msglvl,&ddum,&ddum,&error);
        if(error!=0){
            printf("\nERROR during numerical factorization of interface: %d",error);
            exit(2);}
        std::cout<<"Numerical factorization of interface takes: "<<timer.elapsed()<<" secs"<<std::endl;
        timer.reset();
        //############################################
        // Interior
        //############################################
        for(int subdomain_index=0;subdomain_index<solver_data.Get_Number_Of_Subdomains();++subdomain_index){
            /* -------------------------------------------------------------------- */
            /* .. Numerical factorization. */
            /* -------------------------------------------------------------------- */
            phase=22;
            PARDISO(&solver_data.pardiso_pt_i[level][subdomain_index][0],&maxfct,&mnum,&mtype,&phase,
                    &solver_data.Aii_n[level][subdomain_index],solver_data.Aii_a[level][subdomain_index],
                    solver_data.Aii_ia[level][subdomain_index],solver_data.Aii_ja[level][subdomain_index],
                    &idum,&nrhs,&solver_data.pardiso_iparm_i[level][subdomain_index][0],&msglvl,&ddum,&ddum,&error);
            if(error!=0){
                printf("\nERROR during numerical factorization in interface at level %d: %d",level,error);
                exit(2);}}
        std::cout<<"Numerical factorization of interior takes: "<<timer.elapsed()<<" secs"<<std::endl;
    }
}
//#####################################################################
// Function Smooth
//#####################################################################
template<typename T,int d> void SIGMA_SMOOTHER<T,d>::
Smooth(const int level,T_EIGEN_VECTOR& u,T_EIGEN_VECTOR& tmp,const T_EIGEN_VECTOR& b)const{
    //std::stringstream ss;
    //ss<<level;
    //LOG::SCOPE scope(std::string("Smooth @ ") + ss.str());
    //PHYSBAM_ASSERT((b-solver_data.Arr_dof_distribution_matrix[level]*(solver_data.Arr_dof_distribution_matrix[level].transpose()*b)).maxCoeff()==0);
    tmp=b;
    MKL_INT maxfct=1;                   /* Maximum number of numerical factorizations. */
    MKL_INT mnum=1;                     /* Which factorization to use. */
    MKL_INT msglvl=0;                   /* Print statistical information in file */
    MKL_INT error=0;                    /* Initialize error flag */
    MKL_INT nrhs=1;
    MKL_INT mtype=2;
    MKL_INT idum;
    MKL_INT phase=33;
    for(int subdomain_index=0;subdomain_index<solver_data.Get_Number_Of_Subdomains();++subdomain_index){
        if(solver_data.Aii_n[level][subdomain_index]==0) continue;
        solver_data.Aii_vb[level][subdomain_index]=solver_data.Air[level][subdomain_index]*u;
        PARDISO(&solver_data.pardiso_pt_i[level][subdomain_index][0],&maxfct,&mnum,&mtype,&phase,
                &solver_data.Aii_n[level][subdomain_index],solver_data.Aii_a[level][subdomain_index],
                solver_data.Aii_ia[level][subdomain_index],solver_data.Aii_ja[level][subdomain_index],
                &idum,&nrhs,&solver_data.pardiso_iparm_i[level][subdomain_index][0],&msglvl,
                solver_data.Aii_vb[level][subdomain_index].data(),solver_data.Aii_vx[level][subdomain_index].data(),&error);
        if(error!=0){printf("\nERROR during interior solve: %d",error);exit(3);}
        tmp+=solver_data.Air[level][subdomain_index].transpose()*solver_data.Aii_vx[level][subdomain_index];
    }
    PARDISO(&solver_data.pardiso_pt_r[level][0],&maxfct,&mnum,&mtype,&phase,
            &solver_data.Arr_n[level],solver_data.Arr_a[level],
            solver_data.Arr_ia[level],solver_data.Arr_ja[level],
            &idum,&nrhs,&solver_data.pardiso_iparm_r[level][0],&msglvl,
            tmp.data(),u.data(),&error);
    if(error!=0){printf("\nERROR during interface solve: %d",error);exit(3);}
}
//#####################################################################
// Function Compute_Residual
//#####################################################################
template<typename T,int d> void SIGMA_SMOOTHER<T,d>::
Compute_Residual(const int level,const T_EIGEN_VECTOR& u,T_EIGEN_VECTOR& r,const T_EIGEN_VECTOR& b)const{
    LOG::SCOPE scope("Compute Residual");
    r=b;
    MKL_INT maxfct=1;                   /* Maximum number of numerical factorizations. */
    MKL_INT mnum=1;                     /* Which factorization to use. */
    MKL_INT msglvl=0;                   /* Print statistical information in file */
    MKL_INT error=0;                    /* Initialize error flag */
    MKL_INT nrhs=1;
    MKL_INT mtype=2;
    MKL_INT idum;
    MKL_INT phase=33;
    for(int subdomain_index=0;subdomain_index<solver_data.Get_Number_Of_Subdomains();++subdomain_index){
        if(solver_data.Aii_n[level][subdomain_index]==0) continue;
        //T_EIGEN_VECTOR tmp=solver_data.Air[level][subdomain_index]*u;
        //solver_data.Aii_vb[level][subdomain_index]=tmp;
        solver_data.Aii_vb[level][subdomain_index]=solver_data.Air[level][subdomain_index]*u;
        PARDISO(&solver_data.pardiso_pt_i[level][subdomain_index][0],&maxfct,&mnum,&mtype,&phase,
                &solver_data.Aii_n[level][subdomain_index],solver_data.Aii_a[level][subdomain_index],
                solver_data.Aii_ia[level][subdomain_index],solver_data.Aii_ja[level][subdomain_index],
                &idum,&nrhs,&solver_data.pardiso_iparm_i[level][subdomain_index][0],&msglvl,
                solver_data.Aii_vb[level][subdomain_index].data(),solver_data.Aii_vx[level][subdomain_index].data(),&error);
        T_EIGEN_VECTOR tmp2=solver_data.Air[level][subdomain_index].transpose()*solver_data.Aii_vx[level][subdomain_index];
        r+=solver_data.Air[level][subdomain_index].transpose()*solver_data.Aii_vx[level][subdomain_index];
        if(error!=0){printf("\nERROR during interior solve: %d",error);exit(3);}
    }
    r-=solver_data.Arr[level]*u;
}
//#####################################################################
// Function Multiply
//#####################################################################
template<typename T,int d> void SIGMA_SMOOTHER<T,d>::
Multiply(const int level,const T_EIGEN_VECTOR& u,T_EIGEN_VECTOR& result)const{
    result=T_EIGEN_VECTOR::Zero(result.size());
    MKL_INT maxfct=1;                   /* Maximum number of numerical factorizations. */
    MKL_INT mnum=1;                     /* Which factorization to use. */
    MKL_INT msglvl=0;                   /* Print statistical information in file */
    MKL_INT error=0;                    /* Initialize error flag */
    MKL_INT nrhs=1;
    MKL_INT mtype=2;
    MKL_INT idum;
    MKL_INT phase=33;
    for(int subdomain_index=0;subdomain_index<solver_data.Get_Number_Of_Subdomains();++subdomain_index){
        if(solver_data.Aii_n[level][subdomain_index]==0) continue;
        solver_data.Aii_vb[level][subdomain_index]=solver_data.Air[level][subdomain_index]*u;
        PARDISO(&solver_data.pardiso_pt_i[level][subdomain_index][0],&maxfct,&mnum,&mtype,&phase,
            &solver_data.Aii_n[level][subdomain_index],solver_data.Aii_a[level][subdomain_index],
            solver_data.Aii_ia[level][subdomain_index],solver_data.Aii_ja[level][subdomain_index],
            &idum,&nrhs,&solver_data.pardiso_iparm_i[level][subdomain_index][0],&msglvl,
            solver_data.Aii_vb[level][subdomain_index].data(),solver_data.Aii_vx[level][subdomain_index].data(),&error);
        result-=solver_data.Air[level][subdomain_index].transpose()*solver_data.Aii_vx[level][subdomain_index];
        if(error!=0){printf("\nERROR during interior solve: %d",error);exit(3);}
    }
    result+=solver_data.Arr[level]*u;
}
//#####################################################################
template class SIGMA_SMOOTHER<float,2>;
template class SIGMA_SMOOTHER<float,3>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class SIGMA_SMOOTHER<double,2>;
template class SIGMA_SMOOTHER<double,3>;
#endif

