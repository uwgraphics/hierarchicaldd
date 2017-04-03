//#####################################################################
// Copyright (c) 2012-2013, Sean Bauer, Eftychios Sifakis.
// Distributed under the FreeBSD license (see license.txt)
//#####################################################################
// Utility classes/functions
//#####################################################################
#include <iostream>
#include <stdexcept>
#include <stdio.h>
#include <errno.h>
#include <sys/mman.h>
#include <unistd.h>

#include <SPGrid/Core/SPGrid_Utilities.h>

#ifdef HASWELL
#include <xmmintrin.h>
#endif

namespace SPGrid{

//#####################################################################
// Functions Bit_Spread/Bit_Pack
//#####################################################################

#ifdef HASWELL
unsigned long Bit_Spread(const unsigned int data,const unsigned long mask)
{unsigned long uldata=data;return _pdep_u64(uldata,mask);}
unsigned long Bit_Spread(const int data,const unsigned long mask)
{union{ signed long sldata; unsigned long uldata; };sldata=data;return _pdep_u64(uldata,mask);}
#endif

int Bit_Pack(const unsigned long data, const unsigned long mask)
{
    union{ signed long slresult; unsigned long ulresult; };    
#ifdef HASWELL
    ulresult=_pext_u64(data,mask);
#else
    unsigned long uldata=data; int count=0; ulresult=0;
    for(unsigned long bit=1;bit;bit<<=1)
        if(bit & mask) ulresult |= (uldata & bit)>>count; else count++;
#endif
    return (int)slresult;
}

//#####################################################################
// Function Fatal_Error
//#####################################################################
void Fatal_Error(const char* function,const char* file,unsigned int line)
{
    Fatal_Error(function,file,line,"Fatal error");
}
void Fatal_Error(const char* function,const char* file,unsigned int line,const char* message)
{
    Fatal_Error(function,file,line,std::string(message));
}
void Fatal_Error(const char* function,const char* file,unsigned int line,const std::string& message)
{
    static char buffer[2048];
    sprintf(buffer,"%s:%s:%d: %s",file,function,line,message.c_str());
    std::string error(buffer);
    std::cout<<std::flush;std::cerr<<"\n";
    std::cerr<<"\n*** ERROR: "<<error<<'\n'<<std::endl;
    throw std::runtime_error(error);
}

//#####################################################################
// Function Check_Compliance
//#####################################################################
void Check_Compliance()
{
    if(sysconf(_SC_PAGESIZE)!=4096) FATAL_ERROR("Page size different than 4KB detected");
    if(sizeof(unsigned long)!=8) FATAL_ERROR("unsigned long is not 64-bit integer");
    if(sizeof(size_t)!=8) FATAL_ERROR("size_t is not 64-bit long");
    if(sizeof(void*)!=8) FATAL_ERROR("void* is not 64-bit long");
    typedef enum {dummy=0xffffffffffffffffUL} Long_Enum;
    if(sizeof(Long_Enum)!=8) FATAL_ERROR("Missing support for 64-bit enums");
}

//#####################################################################
// Function Raw_Allocate
//#####################################################################
void* Raw_Allocate(const size_t size)
{
    void *ptr=mmap(NULL,size,PROT_READ|PROT_WRITE,MAP_PRIVATE|MAP_ANONYMOUS|MAP_NORESERVE,-1,0);
    if(ptr==MAP_FAILED) FATAL_ERROR("Failed to allocate "+Value_To_String(size)+" bytes");
    if(0xfffUL&(unsigned long)ptr) FATAL_ERROR("Allocated pointer value "+Value_To_String(ptr)+" is not page-aligned");
    return ptr;
}

//#####################################################################
// Function Raw_Deallocate
//#####################################################################
void Raw_Deallocate(void* data, const size_t size)
{
    if(munmap(data,size)!=0) FATAL_ERROR("Failed to deallocate "+Value_To_String(size)+" bytes");
}

//#####################################################################
// Function Check_Address_Resident
//#####################################################################
void Check_Address_Resident(const void* addr)
{
    void* page_addr=reinterpret_cast<void*>(reinterpret_cast<unsigned long>(addr)&0xfffffffffffff000UL);
    unsigned char status;
    if(mincore(page_addr,4096,&status)==-1)
        switch(errno){
            case ENOMEM: FATAL_ERROR("In Check_Address_Resident() : Input address "+Value_To_String(addr)+" has not been mapped");
            default: FATAL_ERROR(" In Check_Address_Resident() : mincore() failed with errno="+Value_To_String(errno));}
    if(!status) FATAL_ERROR("In Check_Address_Resident() : Input address "+Value_To_String(addr)+" is not resident in physical memory");
}
//#####################################################################
// Function Validate_Memory_Use
//#####################################################################
// WARNING: This validation function is very conservative. It can fail (without truly indicating a problem) when using transparent hugepages
void Validate_Memory_Use(unsigned long number_of_pages,void *data_ptr,unsigned long *page_mask_array)
{
    if(number_of_pages%64) FATAL_ERROR("In Validate_Memory_Use(): Number of mapped pages must be divisible by 64");
    for(unsigned long base_page=0;base_page<number_of_pages;base_page+=64){
        unsigned char vec[64];
        if(mincore(data_ptr,262144,vec)==-1)
            switch(errno){
                case ENOMEM: FATAL_ERROR("In Validate_Memory_Use() : Input address "+Value_To_String(data_ptr)+" has not been mapped");
                default: FATAL_ERROR(" In Validate_Memory_Use() : mincore() failed with errno="+Value_To_String(errno));}
        unsigned long page_bit=1UL;
        for(int i=0;i<64;i++,page_bit<<=1){
            bool mincore_set = vec[i];
            bool page_mask_set = (*page_mask_array & page_bit);
            if(mincore_set && !page_mask_set){
                for(int qword=0;qword<512;qword++)
                    if(*reinterpret_cast<unsigned long*>(reinterpret_cast<unsigned long>(data_ptr)+i*4096UL+qword*8UL))
                        FATAL_ERROR("In Validate_Memory_Use(): Page "+Value_To_String(base_page+i)+" is resident (and nonzero) in physical memory, but not marked in SPGrid_Set");}
            else if(!mincore_set && page_mask_set)
                FATAL_ERROR("In Validate_Memory_Use(): Page "+Value_To_String(base_page+i)+" is not resident in physical memory, but marked in SPGrid_Set");}
        data_ptr=reinterpret_cast<void*>(reinterpret_cast<unsigned long>(data_ptr)+262144UL);
        page_mask_array++;
    }
}
//#####################################################################
}

