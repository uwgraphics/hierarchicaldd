//#####################################################################
//Removed For Anonymity: Copyright Authors of Submission pap173s1
// This file is part of PhysBAM whose distribution is governed by the license 
//Removed For Anonymity: Copyright Authors of Submission pap173s1
//#####################################################################
// Class SPGRID_EIGEN_WRAPPER
//#####################################################################
#ifndef __SPGRID_EIGEN_WRAPPER_H__
#define __SPGRID_EIGEN_WRAPPER_H__
#include <vector>
#include <Eigen/Sparse>
#include <Common_Tools/Math_Tools/RANGE_ITERATOR.h>
#include <PhysBAM_Tools/Data_Structures/HASHTABLE.h>
#include <PhysBAM_Tools/Data_Structures/HASHTABLE_ITERATOR.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <SPGrid/Core/SPGrid_Allocator.h>
#include <SPGrid/Tools/SPGrid_Block_Iterator.h>
#include <SPGrid_Fluids/Grids/GRID_TOPOLOGY_HELPER.h>
#include "SPGRID_MULTIGRID_FLAGS.h"
#include "../Common_Library/Eigen_Matrix_Wrapper.h"
using namespace PhysBAM;

template<typename T,int d>
class SPGrid_Eigen_Wrapper{
    typedef VECTOR<int,d> T_INDEX;
    typedef Eigen::Matrix<T,Eigen::Dynamic,1> T_EIGEN_VECTOR;
public:
    template<typename T_STRUCT>
    static void Construct_From_SPGrid(Eigen::SparseMatrix<T>& matrix,SPGrid_Allocator<T_STRUCT,d>& allocator,
                                      const std::pair<const unsigned long*,unsigned>& blocks,unsigned T_STRUCT::* flags_field,
                                      const ARRAY<int,VECTOR<int,d> >& index_map_nd_to_1d){
        typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<const unsigned>::type Const_flag_array_type;
        typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<unsigned>::mask Flag_array_mask;
        //Please resize the matrix before passing in this function
        static const int number_of_face_neighbors=GRID_TOPOLOGY_HELPER<Flag_array_mask>::faces_per_cell;
        unsigned long face_neighbor_offsets[number_of_face_neighbors];
        GRID_TOPOLOGY_HELPER<Flag_array_mask>::Face_Neighbor_Offsets(face_neighbor_offsets); // order is (-x, x, -y, y, -z, z)
        std_array<T_INDEX,d*2> index_offsets;
        for(int axis=1;axis<=d;++axis){
            index_offsets((axis-1)*2)  =-T_INDEX::Axis_Vector(axis);
            index_offsets((axis-1)*2+1)= T_INDEX::Axis_Vector(axis);}
        Const_flag_array_type flags=allocator.Get_Const_Array(flags_field);    
        typedef Eigen::Triplet<T> TRIPLET;
        const RANGE<T_INDEX>& range=index_map_nd_to_1d.Domain_Indices();
        std::vector<TRIPLET> triplet_list;        
        for(SPGrid_Block_Iterator<Flag_array_mask> iterator(blocks);iterator.Valid();iterator.Next()){
            unsigned flag=iterator.Data(flags);
            T diagnal=0;
            T_INDEX index=Flag_array_mask::LinearToCoord(iterator.Offset()).template Cast<T_INDEX>();
            if((!index_map_nd_to_1d.Domain_Indices().Lazy_Inside(index))||index_map_nd_to_1d(index)==-1)continue;
            if(flag&(SPGrid_Solver_Cell_Type_Active|SPGrid_Solver_Cell_Type_Interface)){             
                for(int face=0;face<number_of_face_neighbors;face++){
                    unsigned long offset=face_neighbor_offsets[face];                    
                    if(flag&(SPGrid_Solver_Face_Minus_X_Active<<face)){
                        ++diagnal;                        
                        unsigned neighbor_flag=(iterator.Data(flags,offset));
                        if(neighbor_flag&(SPGrid_Solver_Cell_Type_Active|SPGrid_Solver_Cell_Type_Interface)){
                            T_INDEX neighbor_index=index_offsets(face)+index;
                            if((!index_map_nd_to_1d.Domain_Indices().Lazy_Inside(neighbor_index))||index_map_nd_to_1d(neighbor_index)==-1)continue;
                            triplet_list.push_back(TRIPLET(index_map_nd_to_1d(index),index_map_nd_to_1d(neighbor_index),-1));}}}
                triplet_list.push_back(TRIPLET(index_map_nd_to_1d(index),index_map_nd_to_1d(index),diagnal));}}
        matrix.setFromTriplets(triplet_list.begin(),triplet_list.end());
        //------------------------------ Should we compress here? ----------------------------------//
        matrix.makeCompressed();
    };
    template<typename T_STRUCT>
    static void Construct_From_SPGrid_With_Correction(Eigen::SparseMatrix<T>& matrix,SPGrid_Allocator<T_STRUCT,d>& allocator,
                                                      const std::pair<const unsigned long*,unsigned>& blocks,unsigned T_STRUCT::* flags_field,
                                                      const ARRAY<int,VECTOR<int,d> >& index_map_nd_to_1d,T epsilon=1e-4){
        typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<const unsigned>::type Const_flag_array_type;
        typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<unsigned>::mask Flag_array_mask;
        typedef Eigen::Triplet<T> TRIPLET;
        typedef STENCIL<T,d> T_STENCIL;
        typedef ARRAY<T_STENCIL,T_INDEX> T_MATRIX;
        //Please resize the matrix before passing in this function
        //And this function adding epsilon edges on the GLOBAL laplace matrix!!
        //And this function should always returns PD matrices
        //This function is really inefficient right now......
        static const int number_of_face_neighbors=GRID_TOPOLOGY_HELPER<Flag_array_mask>::faces_per_cell;
        unsigned long face_neighbor_offsets[number_of_face_neighbors];
        GRID_TOPOLOGY_HELPER<Flag_array_mask>::Face_Neighbor_Offsets(face_neighbor_offsets); // order is (-x, x, -y, y, -z, z)
        std_array<T_INDEX,d*2> index_offsets;
        for(int axis=1;axis<=d;++axis){
            index_offsets((axis-1)*2)  =-T_INDEX::Axis_Vector(axis);
            index_offsets((axis-1)*2+1)= T_INDEX::Axis_Vector(axis);}
        Const_flag_array_type flags=allocator.Get_Const_Array(flags_field);
        const RANGE<T_INDEX>& range=index_map_nd_to_1d.Domain_Indices();
        LOG::cout<<matrix.rows()<<" "<<matrix.cols()<<std::endl;
        const T weight=T(1<<(d-1));

        T_MATRIX system_tmp(range);
        //First fill the matrix up with the spgrid block iterator
        //This loop take care of all the edges that has a dof(active or interface) node.
        for(SPGrid_Block_Iterator<Flag_array_mask> iterator(blocks);iterator.Valid();iterator.Next()){
            unsigned flag=iterator.Data(flags);
            T diagnal=0;
            T_INDEX index=Flag_array_mask::LinearToCoord(iterator.Offset()).template Cast<T_INDEX>();
            if((!index_map_nd_to_1d.Domain_Indices().Lazy_Inside(index))||index_map_nd_to_1d(index)==-1)continue;
            if(flag&(SPGrid_Solver_Cell_Type_Active|SPGrid_Solver_Cell_Type_Interface)){             
                for(int face=0;face<number_of_face_neighbors;face++){
                    unsigned long offset=face_neighbor_offsets[face];                    
                    if(flag&(SPGrid_Solver_Face_Minus_X_Active<<face)){
                        ++diagnal;                        
                        unsigned neighbor_flag=(iterator.Data(flags,offset));
                        if(neighbor_flag&(SPGrid_Solver_Cell_Type_Active|SPGrid_Solver_Cell_Type_Interface)){
                            T_INDEX neighbor_index=index_offsets(face)+index;
                            if((!index_map_nd_to_1d.Domain_Indices().Lazy_Inside(neighbor_index))||index_map_nd_to_1d(neighbor_index)==-1)continue;
                            system_tmp(index).Insert(neighbor_index,-1);}}}
                system_tmp(index).Insert(index,diagnal);}}
        Eigen_Wrapper<T,d>::Construct_From_Stencil_With_Correction(matrix,system_tmp,index_map_nd_to_1d);
    };
    template<typename T_STRUCT>
    static void Construct_From_SPGrid(Eigen::SparseMatrix<T>& matrix,SPGrid_Allocator<T_STRUCT,d>& allocator,
                                      const std::pair<const unsigned long*,unsigned>& blocks,unsigned T_STRUCT::* flags_field,
                                      HASHTABLE<VECTOR<int,d>,int>& index_map_nd_to_1d){
        typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<const unsigned>::type Const_flag_array_type;
        typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<unsigned>::mask Flag_array_mask;
        //Please resize the matrix before passing in this function
        static const int number_of_face_neighbors=GRID_TOPOLOGY_HELPER<Flag_array_mask>::faces_per_cell;
        unsigned long face_neighbor_offsets[number_of_face_neighbors];
        GRID_TOPOLOGY_HELPER<Flag_array_mask>::Face_Neighbor_Offsets(face_neighbor_offsets); // order is (-x, x, -y, y, -z, z)
        std_array<T_INDEX,d*2> index_offsets;
        for(int axis=1;axis<=d;++axis){
            index_offsets((axis-1)*2)  =-T_INDEX::Axis_Vector(axis);
            index_offsets((axis-1)*2+1)= T_INDEX::Axis_Vector(axis);}
        Const_flag_array_type flags=allocator.Get_Const_Array(flags_field);    
        typedef Eigen::Triplet<T> TRIPLET;
        std::vector<TRIPLET> triplet_list;
        for(SPGrid_Block_Iterator<Flag_array_mask> iterator(blocks);iterator.Valid();iterator.Next()){
            unsigned flag=iterator.Data(flags);
            T diagnal=0;
            T_INDEX index=Flag_array_mask::LinearToCoord(iterator.Offset()).template Cast<T_INDEX>();
            int index_1d=index_map_nd_to_1d.Get_Default(index,-1);
            if(index_1d==-1)continue;
            PHYSBAM_ASSERT(!(flag&SPGrid_Solver_Cell_Type_Active));
            if(flag&(SPGrid_Solver_Cell_Type_Active|SPGrid_Solver_Cell_Type_Interface)){             
                for(int face=0;face<number_of_face_neighbors;face++){
                    unsigned long offset=face_neighbor_offsets[face];                    
                    if(flag&(SPGrid_Solver_Face_Minus_X_Active<<face)){
                        ++diagnal;       
                        unsigned neighbor_flag=(iterator.Data(flags,offset));
                        if(neighbor_flag&(SPGrid_Solver_Cell_Type_Active|SPGrid_Solver_Cell_Type_Interface)){
                            T_INDEX neighbor_index=index_offsets(face)+index;
                            int neighbor_index_1d=index_map_nd_to_1d.Get_Default(neighbor_index,-1);
                            if(neighbor_index_1d==-1)continue;
                            triplet_list.push_back(TRIPLET(index_1d,neighbor_index_1d,-1));}}}
                triplet_list.push_back(TRIPLET(index_1d,index_1d,diagnal));}}
        matrix.setFromTriplets(triplet_list.begin(),triplet_list.end());
        //------------------------------ Should we compress here? ----------------------------------//
        matrix.makeCompressed();
    };
    template<typename T_STRUCT>
    static void Construct_From_SPGrid(Eigen::SparseMatrix<T>& matrix,SPGrid_Allocator<T_STRUCT,d>& allocator,
                                       const std::pair<const unsigned long*,unsigned>& blocks,unsigned T_STRUCT::* flags_field,
                                       const ARRAY<int,VECTOR<int,d> >& index_map1,const HASHTABLE<VECTOR<int,d>,int>& index_map2){
        //Please resize the matrix before passing in this function
        //This function is used to generate the Air.
        typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<const unsigned>::type Const_flag_array_type;
        typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<unsigned>::mask Flag_array_mask;
        //Please resize the matrix before passing in this function
        static const int number_of_face_neighbors=GRID_TOPOLOGY_HELPER<Flag_array_mask>::faces_per_cell;
        unsigned long face_neighbor_offsets[number_of_face_neighbors];
        GRID_TOPOLOGY_HELPER<Flag_array_mask>::Face_Neighbor_Offsets(face_neighbor_offsets); // order is (-x, x, -y, y, -z, z)
        std_array<T_INDEX,d*2> index_offsets;
        for(int axis=1;axis<=d;++axis){
            index_offsets((axis-1)*2)  =-T_INDEX::Axis_Vector(axis);
            index_offsets((axis-1)*2+1)= T_INDEX::Axis_Vector(axis);}
        Const_flag_array_type flags=allocator.Get_Const_Array(flags_field);    
        typedef Eigen::Triplet<T> TRIPLET;
        std::vector<TRIPLET> triplet_list;        
        for(SPGrid_Block_Iterator<Flag_array_mask> iterator(blocks);iterator.Valid();iterator.Next()){
            unsigned flag=iterator.Data(flags);
            T diagnal=0;
            T_INDEX index=Flag_array_mask::LinearToCoord(iterator.Offset()).template Cast<T_INDEX>();
            if((!index_map1.Domain_Indices().Lazy_Inside(index))||index_map1(index)==-1)continue;
            int index_1d=index_map1(index);
            if(index_1d==-1)continue;
            if(flag&(SPGrid_Solver_Cell_Type_Active|SPGrid_Solver_Cell_Type_Interface)){             
                for(int face=0;face<number_of_face_neighbors;face++){
                    unsigned long offset=face_neighbor_offsets[face];                    
                    if(flag&(SPGrid_Solver_Face_Minus_X_Active<<face)){
                        ++diagnal;                        
                        unsigned neighbor_flag=(iterator.Data(flags,offset));
                        if(neighbor_flag&(SPGrid_Solver_Cell_Type_Active|SPGrid_Solver_Cell_Type_Interface)){
                            T_INDEX neighbor_index=index_offsets(face)+index;
                            int neighbor_index_1d=index_map2.Get_Default(neighbor_index,-1);
                            if(neighbor_index_1d==-1)continue;
                            triplet_list.push_back(TRIPLET(index_1d,neighbor_index_1d,-1));}}}
                /*PHYSBAM_ASSERT(diagnal==1||diagnal==0);*/}}
        matrix.setFromTriplets(triplet_list.begin(),triplet_list.end());
        //------------------------------ Should we compress here? ----------------------------------//
        matrix.makeCompressed();
    };
    template<typename T_STRUCT>
    static void Copy_To_Eigen_Array(T_EIGEN_VECTOR& array_out,SPGrid_Allocator<T_STRUCT,d>& allocator,
                                    const std::pair<const unsigned long*,unsigned>& blocks,T T_STRUCT::* u_field,
                                    const std::vector<VECTOR<int,d> >& index_map_1d_to_nd){
        //----------------------------------CAUTION! THIS FUNCTION WILL TOUCH UNMAPPED BLOCKS---------------------------//
        typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<const unsigned>::type Const_flag_array_type;
        typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<const T>::type Const_data_array_type;
        typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<unsigned>::mask Flag_array_mask;
        Const_data_array_type u=allocator.Get_Const_Array(u_field);
        for(unsigned int i = 0;i < index_map_1d_to_nd.size();++i)
            array_out(i) = u(Flag_array_mask::Linear_Offset(std_array<int,d>(index_map_1d_to_nd[i])));
    }
    template<typename T_STRUCT>
    static void Copy_From_Eigen_Array(const T_EIGEN_VECTOR& array_in,SPGrid_Allocator<T_STRUCT,d>& allocator,
                                      const std::pair<const unsigned long*,unsigned>& blocks,T T_STRUCT::* u_field,unsigned T_STRUCT::* flag_field,
                                      const std::vector<VECTOR<int,d> >& index_map_1d_to_nd){
        //----------------------------------CAUTION! THIS FUNCTION WILL TOUCH UNMAPPED BLOCKS---------------------------//
        typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<const unsigned>::type Const_flag_array_type;
        typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<T>::type Data_array_type;
        typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<unsigned>::mask Flag_array_mask;
        Data_array_type u=allocator.Get_Array(u_field);
        Const_flag_array_type flags=allocator.Get_Const_Array(flag_field);
        for(unsigned int i=0;i<index_map_1d_to_nd.size();++i){
            unsigned long offset=Flag_array_mask::Linear_Offset(std_array<int,d>(index_map_1d_to_nd[i]));
            if(flags(offset)&(SPGrid_Solver_Cell_Type_Active|SPGrid_Solver_Cell_Type_Interface))
                u(offset)=array_in(i);}
    }
    template<typename T_STRUCT>
    static void Add_From_Eigen_Array(const T_EIGEN_VECTOR& array_in,SPGrid_Allocator<T_STRUCT,d>& allocator,
                                     const std::pair<const unsigned long*,unsigned>& blocks,T T_STRUCT::* u_field,unsigned T_STRUCT::* flag_field,
                                     const std::vector<VECTOR<int,d> >& index_map_1d_to_nd){
        //----------------------------------CAUTION! THIS FUNCTION WILL TOUCH UNMAPPED BLOCKS---------------------------//
        typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<const unsigned>::type Const_flag_array_type;
        typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<T>::type Data_array_type;
        typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<unsigned>::mask Flag_array_mask;
        Data_array_type u=allocator.Get_Array(u_field);
        Const_flag_array_type flags=allocator.Get_Const_Array(flag_field);
        for(unsigned int i=0;i<index_map_1d_to_nd.size();++i){
            unsigned long offset=Flag_array_mask::Linear_Offset(std_array<int,d>(index_map_1d_to_nd[i]));
            if(flags(offset)&(SPGrid_Solver_Cell_Type_Active|SPGrid_Solver_Cell_Type_Interface))
                u(offset)+=array_in(i);}
    }
    static void Make_Positive_Definite(Eigen::SparseMatrix<T>& matrix){
        //basically this function put ones on the zero diagnal entries
        PHYSBAM_ASSERT(matrix.rows() == matrix.cols());
        PHYSBAM_ASSERT(matrix.rows() == matrix.outerSize());
        for(int k = 0;k < matrix.outerSize();++k)
                if(matrix.coeff(k,k) == 0) matrix.insert(k,k) = 1;
        matrix.makeCompressed();
    }
    static bool Matrix_Integrity_Check(const Eigen::SparseMatrix<T>& matrix){
        for(int k = 0;k<matrix.outerSize();++k)            
            if(matrix.coeff(k,k)==0)//this node is not a dof. make sure nothing in that row (or column) has any values
                for(typename Eigen::SparseMatrix<T>::InnerIterator itr(matrix,k);itr;++itr)
                    if(itr.value()!=0) {LOG::cout<<itr.row()<<" "<<itr.col()<<": "<<itr.value()<<std::endl;return false;}
        return true;
    }
    template<typename T_STRUCT>
    static void SPGrid_Linear_Iterating_Cartographer(const std::pair<const unsigned long*,unsigned>& blocks,const RANGE<T_INDEX>& range,ARRAY<int,VECTOR<int,d> >& index_map_nd_to_1d,std::vector<VECTOR<int,d> >& index_map_1d_to_nd){
        typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<const unsigned>::type Const_flag_array_type;
        typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<unsigned>::mask Flag_array_mask;
        typedef VECTOR<int,d> T_INDEX;
        index_map_nd_to_1d.Resize(range);
        index_map_nd_to_1d.Fill(-1);
        int dof_counter=0;//note that it is start from 0;
        for(SPGrid_Block_Iterator<Flag_array_mask> iterator(blocks);iterator.Valid();iterator.Next()){
            T_INDEX index=Flag_array_mask::LinearToCoord(iterator.Offset()).template Cast<T_INDEX>();
            if(!range.Lazy_Inside(index)) continue;
            index_map_nd_to_1d(index)=dof_counter++;
            index_map_1d_to_nd.push_back(index);}
        PHYSBAM_ASSERT(dof_counter==index_map_1d_to_nd.size());
    }
};

#endif
