//#####################################################################
//Removed For Anonymity: Copyright Authors of Submission pap173s1
// Distributed under the FreeBSD license (see license.txt)
//#####################################################################
#ifndef __SPGRID_PHI_ARRAY_LINEARIZER__
#define __SPGRID_PHI_ARRAY_LINEARIZER__
#include <SPGrid/Core/SPGrid_Allocator.h>
#include <vector>
namespace SPGrid{
    template<class T,int log2_struct, int d,class T_offset_ptr> class SPGrid_PHI_Array_Linearizer;
    /////////////////////////////////////////////////////////////////
    //Basicly this class maps cpu mmaped spgrid onto a dense array on phi //
    /////////////////////////////////////////////////////////////////
    template<class T,int log2_struct,class T_offset_ptr>
    class SPGrid_PHI_Array_Linearizer<T,log2_struct,3,T_offset_ptr>{
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
        unsigned int number_of_prolongation_blocks;

        char* data;//this contains the linearized data!
        T_offset_ptr data_buffer_size;
        unsigned int number_of_blocks;
        unsigned int number_of_boundary_blocks;

        char* interface_data;//this is the buffer for transfer blocks containing interface data;
        unsigned int interface_data_buffer_size;
        unsigned int number_of_interface_blocks;
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
        SPGrid_PHI_Array_Linearizer()
            :b(NULL),b_x_minus(NULL),b_x_plus(NULL),b_y_minus(NULL),
             b_y_plus(NULL),b_z_minus(NULL),b_z_plus(NULL),
             b_boundary(NULL),b_boundary_x_minus(NULL),b_boundary_x_plus(NULL),b_boundary_y_minus(NULL),
             b_boundary_y_plus(NULL),b_boundary_z_minus(NULL),b_boundary_z_plus(NULL),
             b_interface(NULL),b_interface_x_minus(NULL),b_interface_x_plus(NULL),b_interface_y_minus(NULL),
             b_interface_y_plus(NULL),b_interface_z_minus(NULL),b_interface_z_plus(NULL),
             data(NULL),interface_data(NULL),prolongation_fine_blocks(NULL),prolongation_coarse_blocks(NULL),restriction_fine_blocks(NULL)
        {}
        void Allocate(){
            if(posix_memalign((void**)&b        ,page_size,(number_of_blocks*sizeof(T_offset_ptr)/page_size+1)*page_size)) abort();
            if(posix_memalign((void**)&b_x_minus,page_size,(number_of_blocks*sizeof(T_offset_ptr)/page_size+1)*page_size)) abort();
            if(posix_memalign((void**)&b_x_plus ,page_size,(number_of_blocks*sizeof(T_offset_ptr)/page_size+1)*page_size)) abort();
            if(posix_memalign((void**)&b_y_minus,page_size,(number_of_blocks*sizeof(T_offset_ptr)/page_size+1)*page_size)) abort();
            if(posix_memalign((void**)&b_y_plus ,page_size,(number_of_blocks*sizeof(T_offset_ptr)/page_size+1)*page_size)) abort();
            if(posix_memalign((void**)&b_z_minus,page_size,(number_of_blocks*sizeof(T_offset_ptr)/page_size+1)*page_size)) abort();
            if(posix_memalign((void**)&b_z_plus ,page_size,(number_of_blocks*sizeof(T_offset_ptr)/page_size+1)*page_size)) abort();

            if(posix_memalign((void**)&b_boundary        ,page_size,(number_of_boundary_blocks*sizeof(T_offset_ptr)/page_size+1)*page_size))
                abort();
            if(posix_memalign((void**)&b_boundary_x_minus,page_size,(number_of_boundary_blocks*sizeof(T_offset_ptr)/page_size+1)*page_size)) 
                abort();
            if(posix_memalign((void**)&b_boundary_x_plus ,page_size,(number_of_boundary_blocks*sizeof(T_offset_ptr)/page_size+1)*page_size)) 
                abort();
            if(posix_memalign((void**)&b_boundary_y_minus,page_size,(number_of_boundary_blocks*sizeof(T_offset_ptr)/page_size+1)*page_size)) 
                abort();
            if(posix_memalign((void**)&b_boundary_y_plus ,page_size,(number_of_boundary_blocks*sizeof(T_offset_ptr)/page_size+1)*page_size)) 
                abort();
            if(posix_memalign((void**)&b_boundary_z_minus,page_size,(number_of_boundary_blocks*sizeof(T_offset_ptr)/page_size+1)*page_size)) 
                abort();
            if(posix_memalign((void**)&b_boundary_z_plus ,page_size,(number_of_boundary_blocks*sizeof(T_offset_ptr)/page_size+1)*page_size)) 
                abort();

            if(posix_memalign((void**)&b_interface        ,page_size,(number_of_interface_blocks*sizeof(T_offset_ptr)/page_size+1)*page_size))
                abort();
            if(posix_memalign((void**)&b_interface_x_minus,page_size,(number_of_interface_blocks*sizeof(T_offset_ptr)/page_size+1)*page_size))
                abort();
            if(posix_memalign((void**)&b_interface_x_plus ,page_size,(number_of_interface_blocks*sizeof(T_offset_ptr)/page_size+1)*page_size))
                abort();
            if(posix_memalign((void**)&b_interface_y_minus,page_size,(number_of_interface_blocks*sizeof(T_offset_ptr)/page_size+1)*page_size))
                abort();
            if(posix_memalign((void**)&b_interface_y_plus ,page_size,(number_of_interface_blocks*sizeof(T_offset_ptr)/page_size+1)*page_size))
                abort();
            if(posix_memalign((void**)&b_interface_z_minus,page_size,(number_of_interface_blocks*sizeof(T_offset_ptr)/page_size+1)*page_size))
                abort();
            if(posix_memalign((void**)&b_interface_z_plus ,page_size,(number_of_interface_blocks*sizeof(T_offset_ptr)/page_size+1)*page_size))
                abort();
            
            data_buffer_size=(number_of_blocks+1) * page_size;//plus one for the deadbeef block
            if(posix_memalign((void**)&data,page_size,data_buffer_size)) abort();
            //std::cout<<"allocated deadbeef: "<<number_of_blocks*page_size<<std::endl;

            interface_data_buffer_size=(number_of_interface_blocks+1)*page_size;
            if(posix_memalign((void**)&interface_data,page_size,interface_data_buffer_size)) abort();
        }
        bool Validate_Dead_Beef(){
            T* data_ptr=reinterpret_cast<T*>((unsigned long)data + number_of_blocks * page_size);
            for(int i = 0;i < elements_per_block*4/*4 for 4 channels*/;++i){
                if(data_ptr[i]!=0) {std::cerr<<"Non zero values at deadbeef block. Something went wrong. Aborting...."<<std::endl;abort();}
            }
        }
        ~SPGrid_PHI_Array_Linearizer(){
            if(b) free(b); 
            if(b_x_minus) free(b_x_minus);
            if(b_x_plus) free(b_x_plus);
            if(b_y_minus) free(b_y_minus);
            if(b_y_plus) free(b_y_plus);
            if(b_z_minus) free(b_z_minus);
            if(b_z_plus) free(b_z_plus);

            if(b_boundary) free(b_boundary); 
            if(b_boundary_x_minus) free(b_boundary_x_minus);
            if(b_boundary_x_plus) free(b_boundary_x_plus);
            if(b_boundary_y_minus) free(b_boundary_y_minus);
            if(b_boundary_y_plus) free(b_boundary_y_plus);
            if(b_boundary_z_minus) free(b_boundary_z_minus);
            if(b_boundary_z_plus) free(b_boundary_z_plus);

            if(data) free(data);
            if(interface_data) free(interface_data);
        }
    private:
        template<int v=8>
        inline void Clear(std_array<T_offset_ptr,v>& array)const{for(int i = 0;i < v;++i) array(i) = 0xdeadbeef;}
        inline void CorrectingParity(int& parity)const{
            switch(T_MASK::lob){
            case 0:
                break;
            case 1:
                parity = (parity >> 1) | ((parity & 1) << 2);
                break;
            case 2:
                parity = ((parity & 3) << 2) | (parity >> 2);
                break;
            default:
                abort();
            } 
        }
        void CheckParity(int parity,unsigned long current_offset,unsigned long base_offset){
            if(current_offset == 0xdeadbeef) return;
            //std::cout << parity << std::endl;
            switch(parity){
            case 0:
                if(current_offset != base_offset){std::cerr << "Wrong Parity" << std::endl;abort();}
                break;
            case 1:
                if(current_offset !=T_MASK::Packed_Offset<0                  ,0                  ,block_zsize>(base_offset)){std::cerr << "Wrong Parity" << std::endl;abort();};
                break;
            case 2:
                if(current_offset !=T_MASK::Packed_Offset<0                  ,block_ysize,0                 >(base_offset)){std::cerr << "Wrong Parity" << std::endl;abort();};
                break;
            case 3:
                if(current_offset !=T_MASK::Packed_Offset<0                  ,block_ysize,block_zsize>(base_offset)){std::cerr << "Wrong Parity" << std::endl;abort();};
                break;
            case 4:
                if(current_offset !=T_MASK::Packed_Offset<block_xsize,0                 ,0                  >(base_offset)){std::cerr << "Wrong Parity" << std::endl;abort();};
                break;
            case 5:
                if(current_offset !=T_MASK::Packed_Offset<block_xsize,0                  ,block_zsize>(base_offset)){std::cerr << "Wrong Parity" << std::endl;abort();};
                break;
            case 6:
                if(current_offset !=T_MASK::Packed_Offset<block_xsize,block_ysize,0                 >(base_offset)){std::cerr << "Wrong Parity" << std::endl;abort();};
                break;
            case 7:
                if(current_offset !=T_MASK::Packed_Offset<block_xsize,block_ysize,block_zsize>(base_offset)){std::cerr << "Wrong Parity" << std::endl;abort();};
                break;
            default:
                break;
            }
        }
    };  
};
#endif
