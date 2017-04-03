//#####################################################################
//Removed For Anonymity: Copyright Authors of Submission pap173s1
// Distributed under the FreeBSD license (see license.txt)
//#####################################################################
#ifndef __SPGRID_MASTER_ARRAY_LINEARIZER_H__
#define __SPGRID_MASTER_ARRAY_LINEARIZER_H__
#include <PhysBAM_Tools/Data_Structures/HASHTABLE.h>
#include <PhysBAM_Tools/Data_Structures/HASHTABLE_ITERATOR.h>
#include <SPGrid/Core/SPGrid_Allocator.h>
#include <SPGrid_Fluids/Solvers/Domain_Decomposition/SPGRID_MULTIGRID_FLAGS.h>
#include <vector>
#include <limits>
#include <cuda_runtime.h>
namespace SPGrid{
    template<class T,int log2_struct, int d,class T_offset_ptr> class SPGrid_Master_Array_Linearizer;
    ////////////////////////////////////////////////////////////////////////////////
    //Basicly this class maps cpu mmaped spgrid onto a dense array on accelerator //
    ////////////////////////////////////////////////////////////////////////////////
    template<class T,int log2_struct,class T_offset_ptr>
    class SPGrid_Master_Array_Linearizer<T,log2_struct,3,T_offset_ptr>{
    public:
        typedef T_offset_ptr T_offset_ptr;
        enum{d=3};
        typedef SPGrid_Mask<log2_struct, NextLogTwo<sizeof(T)>::value,d> T_MASK;
        T_offset_ptr* b;           // block offset stream
        T_offset_ptr* b_x_minus;   // block offset stream
        T_offset_ptr* b_x_plus;    // block offset stream
        T_offset_ptr* b_y_minus;   // block offset stream
        T_offset_ptr* b_y_plus;    // block offset stream
        T_offset_ptr* b_z_minus;   // block offset stream
        T_offset_ptr* b_z_plus;    // block offset stream

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

        T_offset_ptr deadbeef_block;
        std::vector<std_array<T_offset_ptr,8> >  prolongation_fine_blocks;
        std::vector<std_array<T_offset_ptr,8> >  prolongation_coarse_blocks;
        std::vector<std_array<T_offset_ptr,27> > restriction_fine_blocks;
        mutable char* data;//this contains the linearized data!
        unsigned int data_buffer_size;
        unsigned int number_of_blocks;
        unsigned int number_of_boundary_blocks;

        mutable char* interface_data;//this is the buffer for transfer blocks containing interface data;
        unsigned int interface_data_buffer_size;
        unsigned int number_of_interface_blocks;
        int subdomain_id;
        enum {
            page_size = 4096u,
            block_base_alignment = page_size*8u,
            block_base_parity_mask = block_base_alignment-1,
            block_xsize = 1u << T_MASK::block_xbits,
            block_ysize = 1u << T_MASK::block_ybits,
            block_zsize = 1u << T_MASK::block_zbits
        };
        ////////////////////////////////////////////////////
        //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        //Hashtable can not hash unsigned long. It is converted to int during hashing process!!!!
        //TODO: fix it!
        ////////////////////////////////////////////////////
        PhysBAM::HASHTABLE<unsigned long,T_offset_ptr> offsets_map;
        std::vector<unsigned long> offsets_list;
        PhysBAM::HASHTABLE<unsigned long,T_offset_ptr> offsets_map_interface;
        std::vector<unsigned long> offsets_list_interface;
        SPGrid_Master_Array_Linearizer()
            :b(NULL),b_x_minus(NULL),b_x_plus(NULL),b_y_minus(NULL),
             b_y_plus(NULL),b_z_minus(NULL),b_z_plus(NULL),
             b_boundary(NULL),b_boundary_x_minus(NULL),b_boundary_x_plus(NULL),b_boundary_y_minus(NULL),
             b_boundary_y_plus(NULL),b_boundary_z_minus(NULL),b_boundary_z_plus(NULL),
             b_interface(NULL),b_interface_x_minus(NULL),b_interface_x_plus(NULL),b_interface_y_minus(NULL),
             b_interface_y_plus(NULL),b_interface_z_minus(NULL),b_interface_z_plus(NULL),
             data(NULL),interface_data(NULL)
        {}
        SPGrid_Master_Array_Linearizer(const std::pair<const unsigned long*,unsigned>& blocks)
            :b(NULL),b_x_minus(NULL),b_x_plus(NULL),b_y_minus(NULL),
             b_y_plus(NULL),b_z_minus(NULL),b_z_plus(NULL),
             b_boundary(NULL),b_boundary_x_minus(NULL),b_boundary_x_plus(NULL),b_boundary_y_minus(NULL),
             b_boundary_y_plus(NULL),b_boundary_z_minus(NULL),b_boundary_z_plus(NULL),
             b_interface(NULL),b_interface_x_minus(NULL),b_interface_x_plus(NULL),b_interface_y_minus(NULL),
             b_interface_y_plus(NULL),b_interface_z_minus(NULL),b_interface_z_plus(NULL),
             data(NULL),interface_data(NULL)
        {
            Initialize(blocks);
        }
        void Initialize(const std::pair<const unsigned long*,unsigned>& blocks){
            number_of_blocks = blocks.second;
            if((unsigned long)(number_of_blocks+1) * page_size > std::template numeric_limits<T_offset_ptr>::max()){
                std::cerr<<"Allocating more than the pointer type."<<std::endl;abort();}
            b         = new T_offset_ptr[number_of_blocks];
            b_x_minus = new T_offset_ptr[number_of_blocks];
            b_x_plus  = new T_offset_ptr[number_of_blocks];
            b_y_minus = new T_offset_ptr[number_of_blocks];
            b_y_plus  = new T_offset_ptr[number_of_blocks];
            b_z_minus = new T_offset_ptr[number_of_blocks];
            b_z_plus  = new T_offset_ptr[number_of_blocks];
            offsets_map.Initialize_New_Table(number_of_blocks);
            offsets_list.resize(number_of_blocks);
            //we have plus one here for the deadbeef block at the end
            data_buffer_size=(number_of_blocks+1)*page_size;

            //if(posix_memalign((void**)&data,64,data_buffer_size)) abort();
            cudaMallocHost(&data,data_buffer_size);

            //populate the hashtable first
            for(int i = 0;i < number_of_blocks;++i){
                b[i] = i * page_size;
                offsets_map.Insert(blocks.first[i],b[i]);
                offsets_list[i]=blocks.first[i];}

            //initialize the deadbeef block
            deadbeef_block = number_of_blocks * page_size;
            std::memset(reinterpret_cast<void*>(reinterpret_cast<unsigned long>(data)+(unsigned long)deadbeef_block),0,page_size);

            //Compute the adjecent offsets
            for(int i = 0;i < number_of_blocks;++i){
                unsigned long offset=blocks.first[i];
                
                unsigned long x_plus_offset =T_MASK::Packed_OffsetXdim<+block_xsize>(offset);
                unsigned long x_minus_offset=T_MASK::Packed_OffsetXdim<-block_xsize>(offset);
                unsigned long y_plus_offset =T_MASK::Packed_OffsetYdim<+block_ysize>(offset);
                unsigned long y_minus_offset=T_MASK::Packed_OffsetYdim<-block_ysize>(offset);
                unsigned long z_plus_offset =T_MASK::Packed_OffsetZdim<+block_zsize>(offset);
                unsigned long z_minus_offset=T_MASK::Packed_OffsetZdim<-block_zsize>(offset);
                
                T_offset_ptr offset_tmp;
                if(offsets_map.Get(offset,offset_tmp)){
                    b[i] = offset_tmp;
                }else{
                    std::cerr << "Interior offset mismatched" << std::endl;
                    abort();}
                if(offsets_map.Get(x_minus_offset,offset_tmp)){
                    b_x_minus[i] = offset_tmp;
                }else{
                    b_x_minus[i] = deadbeef_block;}

                if(offsets_map.Get(x_plus_offset,offset_tmp)){
                    b_x_plus[i] = offset_tmp;
                }else{
                    b_x_plus[i] = deadbeef_block;}

                if(offsets_map.Get(y_minus_offset,offset_tmp)){
                    b_y_minus[i] = offset_tmp;
                }else{
                    b_y_minus[i] = deadbeef_block;}

                if(offsets_map.Get(y_plus_offset,offset_tmp)){
                    b_y_plus[i] = offset_tmp;
                }else{
                    b_y_plus[i] = deadbeef_block;}

                if(offsets_map.Get(z_minus_offset,offset_tmp)){
                    b_z_minus[i] = offset_tmp;
                }else{
                    b_z_minus[i] = deadbeef_block;}

                if(offsets_map.Get(z_plus_offset,offset_tmp)){
                    b_z_plus[i] = offset_tmp;
                }else{
                    b_z_plus[i] = deadbeef_block;}
            }
        }
        void Initialize_Interface(const std::pair<const unsigned long*,unsigned>& blocks){
            number_of_interface_blocks = blocks.second;

            b_interface         = new T_offset_ptr[number_of_interface_blocks];
            b_interface_x_minus = new T_offset_ptr[number_of_interface_blocks];
            b_interface_x_plus  = new T_offset_ptr[number_of_interface_blocks];
            b_interface_y_minus = new T_offset_ptr[number_of_interface_blocks];
            b_interface_y_plus  = new T_offset_ptr[number_of_interface_blocks];
            b_interface_z_minus = new T_offset_ptr[number_of_interface_blocks];
            b_interface_z_plus  = new T_offset_ptr[number_of_interface_blocks];
            offsets_map_interface.Initialize_New_Table(number_of_interface_blocks);
            offsets_list_interface.resize(number_of_interface_blocks);
            interface_data_buffer_size=number_of_interface_blocks*page_size;

            //if(posix_memalign((void**)&interface_data,64,interface_data_buffer_size)) abort();
            cudaMallocHost(&interface_data,data_buffer_size);

            //populate the hashtable first
            for(int i = 0;i < number_of_interface_blocks;++i){
                b_interface[i] = offsets_map.Get(blocks.first[i]);
                offsets_map_interface.Insert(blocks.first[i],i*page_size);
                offsets_list_interface[i]=blocks.first[i];}

            //Compute the adjecent offsets
            for(int i = 0;i < number_of_interface_blocks;++i){
                unsigned long offset=blocks.first[i];
                
                unsigned long x_plus_offset =T_MASK::Packed_OffsetXdim<+block_xsize>(offset);
                unsigned long x_minus_offset=T_MASK::Packed_OffsetXdim<-block_xsize>(offset);
                unsigned long y_plus_offset =T_MASK::Packed_OffsetYdim<+block_ysize>(offset);
                unsigned long y_minus_offset=T_MASK::Packed_OffsetYdim<-block_ysize>(offset);
                unsigned long z_plus_offset =T_MASK::Packed_OffsetZdim<+block_zsize>(offset);
                unsigned long z_minus_offset=T_MASK::Packed_OffsetZdim<-block_zsize>(offset);
                
                T_offset_ptr offset_tmp;
                if(offsets_map.Get(offset,offset_tmp)){
                    b_interface[i] = offset_tmp;
                }else{
                    std::cerr << "Interface offset mismatched" << std::endl;
                    abort();}
                if(offsets_map.Get(x_minus_offset,offset_tmp)){
                    b_interface_x_minus[i] = offset_tmp;
                }else{
                    b_interface_x_minus[i] = deadbeef_block;}

                if(offsets_map.Get(x_plus_offset,offset_tmp)){
                    b_interface_x_plus[i] = offset_tmp;
                }else{
                    b_interface_x_plus[i] = deadbeef_block;}

                if(offsets_map.Get(y_minus_offset,offset_tmp)){
                    b_interface_y_minus[i] = offset_tmp;
                }else{
                    b_interface_y_minus[i] = deadbeef_block;}

                if(offsets_map.Get(y_plus_offset,offset_tmp)){
                    b_interface_y_plus[i] = offset_tmp;
                }else{
                    b_interface_y_plus[i] = deadbeef_block;}

                if(offsets_map.Get(z_minus_offset,offset_tmp)){
                    b_interface_z_minus[i] = offset_tmp;
                }else{
                    b_interface_z_minus[i] = deadbeef_block;}

                if(offsets_map.Get(z_plus_offset,offset_tmp)){
                    b_interface_z_plus[i] = offset_tmp;
                }else{
                    b_interface_z_plus[i] = deadbeef_block;}
            }
        }
        void Initialize_Boundary(const std::pair<const unsigned long*,unsigned>& blocks){
            number_of_boundary_blocks = blocks.second;

            b_boundary         = new T_offset_ptr[number_of_boundary_blocks];
            b_boundary_x_minus = new T_offset_ptr[number_of_boundary_blocks];
            b_boundary_x_plus  = new T_offset_ptr[number_of_boundary_blocks];
            b_boundary_y_minus = new T_offset_ptr[number_of_boundary_blocks];
            b_boundary_y_plus  = new T_offset_ptr[number_of_boundary_blocks];
            b_boundary_z_minus = new T_offset_ptr[number_of_boundary_blocks];
            b_boundary_z_plus  = new T_offset_ptr[number_of_boundary_blocks];

            //Compute the adjecent offsets
            for(int i = 0;i < number_of_boundary_blocks;++i){
                unsigned long offset=blocks.first[i];
                
                unsigned long x_plus_offset =T_MASK::Packed_OffsetXdim<+block_xsize>(offset);
                unsigned long x_minus_offset=T_MASK::Packed_OffsetXdim<-block_xsize>(offset);
                unsigned long y_plus_offset =T_MASK::Packed_OffsetYdim<+block_ysize>(offset);
                unsigned long y_minus_offset=T_MASK::Packed_OffsetYdim<-block_ysize>(offset);
                unsigned long z_plus_offset =T_MASK::Packed_OffsetZdim<+block_zsize>(offset);
                unsigned long z_minus_offset=T_MASK::Packed_OffsetZdim<-block_zsize>(offset);
                
                T_offset_ptr offset_tmp;
                if(offsets_map.Get(offset,offset_tmp)){
                    b_boundary[i] = offset_tmp;
                }else{
                    std::cerr << "Boundary offset mismatched" << std::endl;
                    abort();}
                if(offsets_map.Get(x_minus_offset,offset_tmp)){
                    b_boundary_x_minus[i] = offset_tmp;
                }else{
                    b_boundary_x_minus[i] = deadbeef_block;}

                if(offsets_map.Get(x_plus_offset,offset_tmp)){
                    b_boundary_x_plus[i] = offset_tmp;
                }else{
                    b_boundary_x_plus[i] = deadbeef_block;}

                if(offsets_map.Get(y_minus_offset,offset_tmp)){
                    b_boundary_y_minus[i] = offset_tmp;
                }else{
                    b_boundary_y_minus[i] = deadbeef_block;}

                if(offsets_map.Get(y_plus_offset,offset_tmp)){
                    b_boundary_y_plus[i] = offset_tmp;
                }else{
                    b_boundary_y_plus[i] = deadbeef_block;}

                if(offsets_map.Get(z_minus_offset,offset_tmp)){
                    b_boundary_z_minus[i] = offset_tmp;
                }else{
                    b_boundary_z_minus[i] = deadbeef_block;}

                if(offsets_map.Get(z_plus_offset,offset_tmp)){
                    b_boundary_z_plus[i] = offset_tmp;
                }else{
                    b_boundary_z_plus[i] = deadbeef_block;}
            }
        }
        void PopulateRestrictionOffsets(const std::pair<const unsigned long*,unsigned>& fine_blocks,const std::pair<const unsigned long*,unsigned>& coarse_blocks,
                                        const SPGrid_Master_Array_Linearizer& fine_linearizer){
            if(coarse_blocks.second != number_of_blocks) abort();
            if(fine_blocks.second != fine_linearizer.number_of_blocks) abort();
            if(fine_blocks.second==0&&coarse_blocks.second==0) return;
            const int coarse_n_blocks = number_of_blocks;
            const int fine_n_blocks = fine_linearizer.number_of_blocks;
            restriction_fine_blocks.clear();
            for(int i = 0;i < coarse_n_blocks;++i){
                unsigned long fine_base = T_MASK::template Packed_Offset<1,1,1>(T_MASK::UpsampleOffset(coarse_blocks.first[i]));
                if(fine_base % page_size) {std::cerr << "Upsampling does not much the page." << std::endl;abort();}                        
                std_array<T_offset_ptr,27> fine_tmp;
                Clear<27>(fine_tmp);
                std_array<unsigned long,27> key_tmp;
                key_tmp(0)  = T_MASK::Packed_Offset<-block_xsize,-block_ysize,-block_zsize>(fine_base);
                key_tmp(1)  = T_MASK::Packed_Offset<-block_xsize,-block_ysize,0           >(fine_base);
                key_tmp(2)  = T_MASK::Packed_Offset<-block_xsize,-block_ysize,+block_zsize>(fine_base);

                key_tmp(3)  = T_MASK::Packed_Offset<-block_xsize,0           ,-block_zsize>(fine_base);
                key_tmp(4)  = T_MASK::Packed_Offset<-block_xsize,0           ,0           >(fine_base);
                key_tmp(5)  = T_MASK::Packed_Offset<-block_xsize,0           ,+block_zsize>(fine_base);

                key_tmp(6)  = T_MASK::Packed_Offset<-block_xsize,+block_ysize,-block_zsize>(fine_base);
                key_tmp(7)  = T_MASK::Packed_Offset<-block_xsize,+block_ysize,0           >(fine_base);
                key_tmp(8)  = T_MASK::Packed_Offset<-block_xsize,+block_ysize,+block_zsize>(fine_base);

                key_tmp(9)  = T_MASK::Packed_Offset<0           ,-block_ysize,-block_zsize>(fine_base);
                key_tmp(10) = T_MASK::Packed_Offset<0           ,-block_ysize,0           >(fine_base);
                key_tmp(11) = T_MASK::Packed_Offset<0           ,-block_ysize,+block_zsize>(fine_base);

                key_tmp(12) = T_MASK::Packed_Offset<0           ,0           ,-block_zsize>(fine_base);
                key_tmp(13) = T_MASK::Packed_Offset<0           ,0           ,0           >(fine_base);
                key_tmp(14) = T_MASK::Packed_Offset<0           ,0           ,+block_zsize>(fine_base);

                key_tmp(15) = T_MASK::Packed_Offset<0           ,+block_ysize,-block_zsize>(fine_base);
                key_tmp(16) = T_MASK::Packed_Offset<0           ,+block_ysize,0           >(fine_base);
                key_tmp(17) = T_MASK::Packed_Offset<0           ,+block_ysize,+block_zsize>(fine_base);

                key_tmp(18) = T_MASK::Packed_Offset<+block_xsize,-block_ysize,-block_zsize>(fine_base);
                key_tmp(19) = T_MASK::Packed_Offset<+block_xsize,-block_ysize,0           >(fine_base);
                key_tmp(20) = T_MASK::Packed_Offset<+block_xsize,-block_ysize,+block_zsize>(fine_base);

                key_tmp(21) = T_MASK::Packed_Offset<+block_xsize,0           ,-block_zsize>(fine_base);
                key_tmp(22) = T_MASK::Packed_Offset<+block_xsize,0           ,0           >(fine_base);
                key_tmp(23) = T_MASK::Packed_Offset<+block_xsize,0           ,+block_zsize>(fine_base);

                key_tmp(24) = T_MASK::Packed_Offset<+block_xsize,+block_ysize,-block_zsize>(fine_base);
                key_tmp(25) = T_MASK::Packed_Offset<+block_xsize,+block_ysize,0           >(fine_base);
                key_tmp(26) = T_MASK::Packed_Offset<+block_xsize,+block_ysize,+block_zsize>(fine_base);

                T_offset_ptr offset_tmp;
                for(int p = 0;p < 27;++p){
                    if(fine_linearizer.offsets_map.Get(key_tmp(p),offset_tmp)){
                        fine_tmp(p) = offset_tmp;
                    }else{
                        fine_tmp(p) = fine_linearizer.deadbeef_block;}}

                restriction_fine_blocks.push_back(fine_tmp);
            }
        }
        void PopulateProlongationOffsets(const std::pair<const unsigned long*,unsigned>& fine_blocks,const std::pair<const unsigned long*,unsigned>& coarse_blocks,
                                         const SPGrid_Master_Array_Linearizer& coarse_linearizer){
            prolongation_fine_blocks.clear();
            prolongation_coarse_blocks.clear();
            if(fine_blocks.second != number_of_blocks) abort();
            if(coarse_blocks.second != coarse_linearizer.number_of_blocks) abort();
            if(fine_blocks.second==0&&coarse_blocks.second==0) return;
            std_array<T_offset_ptr,8> fine_tmp;
            unsigned long fine_base=fine_blocks.first[0] & ~((unsigned long)block_base_parity_mask);; 
            Clear(fine_tmp);
            const int fine_n_blocks = number_of_blocks;
            const int coarse_n_blocks = coarse_linearizer.number_of_blocks;
            for(int i = 0;i < fine_n_blocks;++i){
                T_offset_ptr fine_offset_linearized;
                if(i >= 1) if(fine_blocks.first[i] <= fine_blocks.first[i-1]) abort();
                if(!offsets_map.Get(fine_blocks.first[i],fine_offset_linearized)) {std::cerr << "missing offset encountered." << std::endl;abort();}                
                if(fine_offset_linearized != b[i]) {std::cerr << "mismatched lookup table." << std::endl;abort();}                
                // "Left over bits" - lob=0 means page address starts with z-bit, lob=1 is x-bit, lob=2 is y-bit
                int parity = (fine_blocks.first[i] & block_base_parity_mask) / page_size;
                //Use the correction to reorder the 8 fine blocks into z,y,x order. Instead of by the order of the lob.
                CorrectingParity(parity);
                //std::cout << parity << ","<< T_MASK::LinearToCoord(fine_blocks.first[i])<<std::endl;
                fine_tmp(parity) = fine_offset_linearized;

                CheckParity(parity,fine_blocks.first[i],fine_base);
                if((i == fine_n_blocks - 1) || fine_base != (fine_blocks.first[i + 1] & ~((unsigned long)block_base_parity_mask))){                    
                    prolongation_fine_blocks.push_back(fine_tmp);
                    const unsigned long coarse_base = T_MASK::DownsampleOffset(fine_base);
                    if(coarse_base % page_size) {std::cerr << "Downsampling does not much the page." << std::endl;abort();}                        
                    std_array<T_offset_ptr,8> coarse_tmp;
                    Clear(coarse_tmp);
                    std_array<unsigned long,8> key_tmp;
                    key_tmp(0) = coarse_base;
                    key_tmp(1) =T_MASK::Packed_Offset<0          ,0          ,block_zsize>(coarse_base);
                    key_tmp(2) =T_MASK::Packed_Offset<0          ,block_ysize,0          >(coarse_base);
                    key_tmp(3) =T_MASK::Packed_Offset<0          ,block_ysize,block_zsize>(coarse_base);
                    key_tmp(4) =T_MASK::Packed_Offset<block_xsize,0          ,0          >(coarse_base);
                    key_tmp(5) =T_MASK::Packed_Offset<block_xsize,0          ,block_zsize>(coarse_base);
                    key_tmp(6) =T_MASK::Packed_Offset<block_xsize,block_ysize,0          >(coarse_base);
                    key_tmp(7) =T_MASK::Packed_Offset<block_xsize,block_ysize,block_zsize>(coarse_base);
                    T_offset_ptr offset_tmp;
                    for(int p = 0;p < 8;++p){
                        if(coarse_linearizer.offsets_map.Get(key_tmp(p),offset_tmp)){
                            coarse_tmp(p) = offset_tmp;
                        }else{
                            coarse_tmp(p) = coarse_linearizer.deadbeef_block;}}

                    prolongation_coarse_blocks.push_back(coarse_tmp);
                    if(i < fine_n_blocks - 1)
                        fine_base = fine_blocks.first[i + 1] & ~((unsigned long)block_base_parity_mask);
                    Clear(fine_tmp);
                }
            }
        }
        ~SPGrid_Master_Array_Linearizer(){
            if(b) delete[] b; 
            if(b_x_minus) delete[] b_x_minus;
            if(b_x_plus) delete[] b_x_plus;
            if(b_y_minus) delete[] b_y_minus;
            if(b_y_plus) delete[] b_y_plus;
            if(b_z_minus) delete[] b_z_minus;
            if(b_z_plus) delete[] b_z_plus;

            if(b_interface) delete[] b_interface; 
            if(b_interface_x_minus) delete[] b_interface_x_minus;
            if(b_interface_x_plus) delete[] b_interface_x_plus;
            if(b_interface_y_minus) delete[] b_interface_y_minus;
            if(b_interface_y_plus) delete[] b_interface_y_plus;
            if(b_interface_z_minus) delete[] b_interface_z_minus;
            if(b_interface_z_plus) delete[] b_interface_z_plus;

            if(b_boundary) delete[] b_boundary; 
            if(b_boundary_x_minus) delete[] b_boundary_x_minus;
            if(b_boundary_x_plus) delete[] b_boundary_x_plus;
            if(b_boundary_y_minus) delete[] b_boundary_y_minus;
            if(b_boundary_y_plus) delete[] b_boundary_y_plus;
            if(b_boundary_z_minus) delete[] b_boundary_z_minus;
            if(b_boundary_z_plus) delete[] b_boundary_z_plus;

            if(data) cudaFreeHost(data);
            if(interface_data) cudaFreeHost(interface_data);

        }
        void Clear_Buffer()const{
            memset(data,0,data_buffer_size);
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
                abort();} 
        }
        void CheckParity(int parity,unsigned long current_offset,unsigned long base_offset){
            if(current_offset == 0xdeadbeef) return;
            //std::cout << parity << std::endl;
            switch(parity){
            case 0:
                if(current_offset != base_offset){std::cerr << "Wrong Parity" << std::endl;abort();}
                break;
            case 1:
                if(current_offset !=T_MASK::Packed_Offset<0          ,0          ,block_zsize>(base_offset)){std::cerr << "Wrong Parity" << std::endl;abort();};
                break;
            case 2:
                if(current_offset !=T_MASK::Packed_Offset<0          ,block_ysize,0          >(base_offset)){std::cerr << "Wrong Parity" << std::endl;abort();};
                break;
            case 3:
                if(current_offset !=T_MASK::Packed_Offset<0          ,block_ysize,block_zsize>(base_offset)){std::cerr << "Wrong Parity" << std::endl;abort();};
                break;
            case 4:
                if(current_offset !=T_MASK::Packed_Offset<block_xsize,0          ,0          >(base_offset)){std::cerr << "Wrong Parity" << std::endl;abort();};
                break;
            case 5:
                if(current_offset !=T_MASK::Packed_Offset<block_xsize,0          ,block_zsize>(base_offset)){std::cerr << "Wrong Parity" << std::endl;abort();};
                break;
            case 6:
                if(current_offset !=T_MASK::Packed_Offset<block_xsize,block_ysize,0          >(base_offset)){std::cerr << "Wrong Parity" << std::endl;abort();};
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
