//#####################################################################
//Removed For Anonymity: Copyright Authors of Submission pap173s1
// This file is part of PhysBAM whose distribution is governed by the license 
//Removed For Anonymity: Copyright Authors of Submission pap173s1
//#####################################################################
#ifndef __Interface_Solver_Data_h__
#define __Interface_Solver_Data_h__

#include <Eigen/SparseCore>
#include <Eigen/SparseCholesky>
#include <unsupported/Eigen/SparseExtra>
#include <array>
#include <fstream>
#include <vector>
#include "mkl_pardiso.h"
#include "mkl_types.h"

//#define USING_EIGEN


template <typename T1, int whatever, typename IND>
void Serialize(const std::string file_name,Eigen::SparseMatrix<T1, whatever, IND>& m) {
    typedef Eigen::Triplet<int> Trip;
    std::vector<Trip> res;
    int sz = m.nonZeros();
    m.makeCompressed();

    std::fstream writeFile;
    writeFile.open(file_name, std::ios::binary | std::ios::out);

    if(writeFile.is_open())
        {
            IND rows, cols, nnzs, outS, innS;
            rows = m.rows()     ;
            cols = m.cols()     ;
            nnzs = m.nonZeros() ;
            outS = m.outerSize();
            innS = m.innerSize();

            writeFile.write((const char *)&(rows), sizeof(IND));
            writeFile.write((const char *)&(cols), sizeof(IND));
            writeFile.write((const char *)&(nnzs), sizeof(IND));
            writeFile.write((const char *)&(outS), sizeof(IND));
            writeFile.write((const char *)&(innS), sizeof(IND));

            writeFile.write((const char *)(m.valuePtr()),       sizeof(T1 ) * m.nonZeros());
            writeFile.write((const char *)(m.outerIndexPtr()),  sizeof(IND) * m.outerSize());
            writeFile.write((const char *)(m.innerIndexPtr()),  sizeof(IND) * m.nonZeros());

            writeFile.close();
        }
}

template <typename T1, int whatever, typename IND>
void Deserialize(const std::string file_name,Eigen::SparseMatrix<T1, whatever, IND>& m) {
    std::fstream readFile;
    readFile.open(file_name, std::ios::binary | std::ios::in);
    if(readFile.is_open())
        {
            IND rows, cols, nnz, inSz, outSz;
            readFile.read((char*)&rows , sizeof(IND));
            readFile.read((char*)&cols , sizeof(IND));
            readFile.read((char*)&nnz  , sizeof(IND));
            readFile.read((char*)&inSz , sizeof(IND));
            readFile.read((char*)&outSz, sizeof(IND));

            m.resize(rows, cols);
            m.makeCompressed();
            m.resizeNonZeros(nnz);

            readFile.read((char*)(m.valuePtr())     , sizeof(T1 ) * nnz  );
            readFile.read((char*)(m.outerIndexPtr()), sizeof(IND) * outSz);
            readFile.read((char*)(m.innerIndexPtr()), sizeof(IND) * nnz );

            m.finalize();
            readFile.close();

        } // file is open
}

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
    void Write_To_File(std::string directory){
        std::ofstream output;
        const int nlevels=Aii.size();
        const int nsubdomains=Aii[0].size();
        output.open(directory+"/levels",std::ios::out|std::ios::binary);        
        output.write(reinterpret_cast<const char*>(&nlevels),sizeof(nlevels));
        output.close();
        output.open(directory+"/number_of_subdomains",std::ios::out|std::ios::binary);
        output.write(reinterpret_cast<const char*>(&nsubdomains),sizeof(nsubdomains));
        output.close();
        for(int level=0;level<nlevels;++level){
            for(int subdomain=0;subdomain<nsubdomains;++subdomain){
                saveMarket(Aii[level][subdomain],directory+"/Aii_"+std::to_string((long long)level)+"_"+std::to_string((long long)subdomain));
                saveMarket(Air[level][subdomain],directory+"/Air_"+std::to_string((long long)level)+"_"+std::to_string((long long)subdomain));}
                //Serialize<T,Eigen::ColMajor,int>(directory+"/Aii_"+std::to_string((long long)level)+"_"+std::to_string((long long)subdomain),Aii[level][subdomain]);
                //Serialize<T,Eigen::RowMajor,int>(directory+"/Air_"+std::to_string((long long)level)+"_"+std::to_string((long long)subdomain),Air[level][subdomain]);}
            saveMarket(Arr[level],directory+"/Arr_"+std::to_string((long long)level));
            //Serialize<T,Eigen::ColMajor,int>(directory+"/Arr_"+std::to_string((long long)level),Arr[level]);
            if(level!=nlevels-1)
                saveMarket(Prr[level],directory+"/Prr_"+std::to_string((long long)level));}
                //Serialize<T,Eigen::RowMajor,int>(directory+"/Prr_"+std::to_string((long long)level),Prr[level]);}
    }
    void Read_From_File(std::string directory){
        std::ifstream input;
        int nlevels;
        int nsubdomains;
        input.open(directory+"/levels",std::ios::in|std::ios::binary);        
        input.read(reinterpret_cast<char*>(&nlevels),sizeof(nlevels));
        input.close();
        std::cout<<"sizeof(nlevels): "<<sizeof(nlevels)<<std::endl;
        std::cout<<"Reading : "<<directory+"/levels"<<std::endl;
        std::cout<<"Read levels: "<<nlevels<<std::endl;
        input.open(directory+"/number_of_subdomains",std::ios::in|std::ios::binary);
        input.read(reinterpret_cast<char*>(&nsubdomains),sizeof(nsubdomains));
        input.close();
        Resize(nlevels,nsubdomains);
        for(int level=0;level<nlevels;++level){
            for(int subdomain=0;subdomain<nsubdomains;++subdomain){
                loadMarket(Aii[level][subdomain],directory+"/Aii_"+std::to_string((long long)level)+"_"+std::to_string((long long)subdomain));
                loadMarket(Air[level][subdomain],directory+"/Air_"+std::to_string((long long)level)+"_"+std::to_string((long long)subdomain));}
                //Deserialize<T,Eigen::ColMajor,int>(directory+"/Aii_"+std::to_string((long long)level)+"_"+std::to_string((long long)subdomain),Aii[level][subdomain]);
                //Deserialize<T,Eigen::RowMajor,int>(directory+"/Air_"+std::to_string((long long)level)+"_"+std::to_string((long long)subdomain),Air[level][subdomain]);}
            loadMarket(Arr[level],directory+"/Arr_"+std::to_string((long long)level));
            //Deserialize<T,Eigen::ColMajor,int>(directory+"/Arr_"+std::to_string((long long)level),Arr[level]);
            if(level!=nlevels-1) 
                loadMarket(Prr[level],directory+"/Prr_"+std::to_string((long long)level));}
                //Deserialize<T,Eigen::RowMajor,int>(directory+"/Prr_"+std::to_string((long long)level),Prr[level]);}
    }
};
#endif
