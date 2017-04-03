//#####################################################################
//Removed For Anonymity: Copyright Authors of Submission pap173s1
// Distributed under the FreeBSD license (see license.txt)
//#####################################################################
#ifndef __SPGRID_CUDA_ARRAY_LINEARIZER__
#define __SPGRID_CUDA_ARRAY_LINEARIZER__
#include <cuda_runtime_api.h>
#include <cuda.h>
#include <vector>
#include <SPGrid/Core/SPGrid_Allocator.h>
namespace SPGrid{
    template<class T,int log2_struct, int d,class T_offset_ptr> class SPGrid_CUDA_Array_Linearizer;
    ////////////////////////////////////////////////////////////////////////
    //Basicly this class maps cpu mmaped spgrid onto a dense array on GPU //
    ////////////////////////////////////////////////////////////////////////
    template<class T,int log2_struct,class T_offset_ptr>
    class SPGrid_CUDA_Array_Linearizer<T,log2_struct,3,T_offset_ptr>{
    public:
        typedef T_offset_ptr T_offset_ptr;
        enum{d=3};
        typedef SPGrid_Mask<log2_struct, NextLogTwo<sizeof(T)>::value,d> T_MASK;
        T_offset_ptr* b;   // block offset stream
        T_offset_ptr* b_x_minus;   // block offset stream
        T_offset_ptr* b_x_plus;   // block offset stream
        T_offset_ptr* b_y_minus;   // block offset stream
        T_offset_ptr* b_y_plus;   // block offset stream
        T_offset_ptr* b_z_minus;   // block offset stream
        T_offset_ptr* b_z_plus;   // block offset stream

        T_offset_ptr* b_boundary;           // block offset stream
        T_offset_ptr* b_boundary_x_minus;   // block offset stream
        T_offset_ptr* b_boundary_x_plus;    // block offset stream
        T_offset_ptr* b_boundary_y_minus;   // block offset stream
        T_offset_ptr* b_boundary_y_plus;    // block offset stream
        T_offset_ptr* b_boundary_z_minus;   // block offset stream
        T_offset_ptr* b_boundary_z_plus;    // block offset stream

        T_offset_ptr* b_interface;           // block offset stream
        T_offset_ptr* b_interface_x_minus;   // block offset stream
        T_offset_ptr* b_interface_x_plus;    // block offset stream
        T_offset_ptr* b_interface_y_minus;   // block offset stream
        T_offset_ptr* b_interface_y_plus;    // block offset stream
        T_offset_ptr* b_interface_z_minus;   // block offset stream
        T_offset_ptr* b_interface_z_plus;    // block offset stream

        T_offset_ptr* prolongation_fine_blocks;
        T_offset_ptr* prolongation_coarse_blocks;
        T_offset_ptr* restriction_fine_blocks;
        
        char* data;//this contains the linearized data!
        T_offset_ptr data_buffer_size;
        unsigned int number_of_blocks;
        unsigned int number_of_boundary_blocks;

        char* interface_data;//this is the buffer for transfer blocks containing interface data;
        unsigned int interface_data_buffer_size;
        unsigned int number_of_interface_blocks;
        unsigned int n_prolongation_meta_blocks;
        enum {
            page_size = 4096u,
            block_base_alignment = page_size*8u,
            block_base_parity_mask = block_base_alignment-1,
            block_xsize = 1u << T_MASK::block_xbits,
            block_ysize = 1u << T_MASK::block_ybits,
            block_zsize = 1u << T_MASK::block_zbits,
            elements_per_block=block_xsize*block_ysize*block_zsize,
        };
        ////////////////////////////////////////////////////
        //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        //Hashtable can not hash unsigned long. It is converted to int during hashing process!!!!
        //TODO: fix it!
        ////////////////////////////////////////////////////
        SPGrid_CUDA_Array_Linearizer()
            :b(NULL),b_x_minus(NULL),b_x_plus(NULL),b_y_minus(NULL),
             b_y_plus(NULL),b_z_minus(NULL),b_z_plus(NULL),
             b_boundary(NULL),b_boundary_x_minus(NULL),b_boundary_x_plus(NULL),b_boundary_y_minus(NULL),
             b_boundary_y_plus(NULL),b_boundary_z_minus(NULL),b_boundary_z_plus(NULL),
             b_interface(NULL),b_interface_x_minus(NULL),b_interface_x_plus(NULL),b_interface_y_minus(NULL),
             b_interface_y_plus(NULL),b_interface_z_minus(NULL),b_interface_z_plus(NULL),
             data(NULL),interface_data(NULL)
        {}
        void Allocate(){
            if(cudaMalloc((void**)&b        ,number_of_blocks*sizeof(T_offset_ptr))!=cudaSuccess) abort();
            if(cudaMalloc((void**)&b_x_minus,number_of_blocks*sizeof(T_offset_ptr))!=cudaSuccess) abort();
            if(cudaMalloc((void**)&b_x_plus ,number_of_blocks*sizeof(T_offset_ptr))!=cudaSuccess) abort();
            if(cudaMalloc((void**)&b_y_minus,number_of_blocks*sizeof(T_offset_ptr))!=cudaSuccess) abort();
            if(cudaMalloc((void**)&b_y_plus ,number_of_blocks*sizeof(T_offset_ptr))!=cudaSuccess) abort();
            if(cudaMalloc((void**)&b_z_minus,number_of_blocks*sizeof(T_offset_ptr))!=cudaSuccess) abort();
            if(cudaMalloc((void**)&b_z_plus ,number_of_blocks*sizeof(T_offset_ptr))!=cudaSuccess) abort();

            if(cudaMalloc((void**)&b_boundary        ,number_of_boundary_blocks*sizeof(T_offset_ptr))!=cudaSuccess) abort();
            if(cudaMalloc((void**)&b_boundary_x_minus,number_of_boundary_blocks*sizeof(T_offset_ptr))!=cudaSuccess) abort();
            if(cudaMalloc((void**)&b_boundary_x_plus ,number_of_boundary_blocks*sizeof(T_offset_ptr))!=cudaSuccess) abort();
            if(cudaMalloc((void**)&b_boundary_y_minus,number_of_boundary_blocks*sizeof(T_offset_ptr))!=cudaSuccess) abort();
            if(cudaMalloc((void**)&b_boundary_y_plus ,number_of_boundary_blocks*sizeof(T_offset_ptr))!=cudaSuccess) abort();
            if(cudaMalloc((void**)&b_boundary_z_minus,number_of_boundary_blocks*sizeof(T_offset_ptr))!=cudaSuccess) abort();
            if(cudaMalloc((void**)&b_boundary_z_plus ,number_of_boundary_blocks*sizeof(T_offset_ptr))!=cudaSuccess) abort();

            if(cudaMalloc((void**)&b_interface        ,number_of_interface_blocks*sizeof(T_offset_ptr))!=cudaSuccess) abort();
            if(cudaMalloc((void**)&b_interface_x_minus,number_of_interface_blocks*sizeof(T_offset_ptr))!=cudaSuccess) abort();
            if(cudaMalloc((void**)&b_interface_x_plus ,number_of_interface_blocks*sizeof(T_offset_ptr))!=cudaSuccess) abort();
            if(cudaMalloc((void**)&b_interface_y_minus,number_of_interface_blocks*sizeof(T_offset_ptr))!=cudaSuccess) abort();
            if(cudaMalloc((void**)&b_interface_y_plus ,number_of_interface_blocks*sizeof(T_offset_ptr))!=cudaSuccess) abort();
            if(cudaMalloc((void**)&b_interface_z_minus,number_of_interface_blocks*sizeof(T_offset_ptr))!=cudaSuccess) abort();
            if(cudaMalloc((void**)&b_interface_z_plus ,number_of_interface_blocks*sizeof(T_offset_ptr))!=cudaSuccess) abort();
            
            data_buffer_size=(number_of_blocks+1) * page_size;//plus one for the deadbeef block
            if(cudaMalloc((void**)&data,data_buffer_size)!=cudaSuccess) abort();
            //cudaMemset(reinterpret_cast<void*>((unsigned long)data + number_of_blocks * page_size),0,page_size);
            //std::cout<<"allocated deadbeef: "<<number_of_blocks*page_size<<std::endl;

            interface_data_buffer_size=number_of_interface_blocks * page_size;
            if(cudaMalloc((void**)&interface_data,interface_data_buffer_size)!=cudaSuccess) abort();
        }
        ~SPGrid_CUDA_Array_Linearizer(){
            if(b) cudaFree(b); 
            if(b_x_minus) cudaFree(b_x_minus);
            if(b_x_plus) cudaFree(b_x_plus);
            if(b_y_minus) cudaFree(b_y_minus);
            if(b_y_plus) cudaFree(b_y_plus);
            if(b_z_minus) cudaFree(b_z_minus);
            if(b_z_plus) cudaFree(b_z_plus);

            if(b_boundary) cudaFree(b_boundary); 
            if(b_boundary_x_minus) cudaFree(b_boundary_x_minus);
            if(b_boundary_x_plus) cudaFree(b_boundary_x_plus);
            if(b_boundary_y_minus) cudaFree(b_boundary_y_minus);
            if(b_boundary_y_plus) cudaFree(b_boundary_y_plus);
            if(b_boundary_z_minus) cudaFree(b_boundary_z_minus);
            if(b_boundary_z_plus) cudaFree(b_boundary_z_plus);

            if(data) cudaFree(data);
            if(interface_data) cudaFree(interface_data);
        }
    };  
};
#endif
