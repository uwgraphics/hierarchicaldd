//#####################################################################
//Removed For Anonymity: Copyright Authors of Submission pap173s1
//Removed For Anonymity: Copyright Authors of Submission pap173s1
//#####################################################################
// Subroutine SPGrid_Computations::Domain_Divider
//#####################################################################
#ifndef __SPGrid_Domain_Divider_h__
#define __SPGrid_Domain_Divider_h__

#include <Common_Tools/Math_Tools/RANGE_ITERATOR.h>
#include <SPGrid/Core/SPGrid_Allocator.h>
#include "SPGRID_MULTIGRID_FLAGS.h"

namespace SPGrid_Computations{

using namespace SPGrid;
using namespace PhysBAM;

template<typename T_STRUCT,int d>
class Domain_Divider
{
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<unsigned>::type Flag_Array_Type;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<unsigned>::mask T_MASK;
    typedef VECTOR<int,d> T_INDEX;
    
    unsigned T_STRUCT::* flags_field;
    SPGrid_Set<Flag_Array_Type>& set;
    const T_INDEX origin,sub_domain_size,domain_size;
    
public:
    static T_INDEX Number_of_Subdomains(T_INDEX domain_max_corner,T_INDEX origin,T_INDEX sub_domain_size){
        T_INDEX n_subdomains=(domain_max_corner-origin-1)/sub_domain_size;
        for(int axis=1;axis<=d;++axis) if((domain_max_corner(axis)-origin(axis)-1)%sub_domain_size(axis))++n_subdomains(axis);
        return n_subdomains;
    }
    static RANGE<T_INDEX> Get_Sub_Domain_Range(T_INDEX subdomain_index,T_INDEX domain_size,T_INDEX origin,T_INDEX sub_domain_size){
        T_INDEX base = origin + subdomain_index * sub_domain_size;
        T_INDEX top = origin + (subdomain_index + 1) * sub_domain_size;// + T_INDEX::All_Ones_Vector();
        return RANGE<T_INDEX>(base,top);
    }
    Domain_Divider(unsigned T_STRUCT::* flags_field_input,SPGrid_Set<Flag_Array_Type>& set_input,T_INDEX domain_size_input,T_INDEX origin_input,T_INDEX sub_domain_size_input)
        :flags_field(flags_field_input),set(set_input),origin(origin_input),sub_domain_size(sub_domain_size_input),domain_size(domain_size_input)
    {}    
    Domain_Divider(SPGrid_Allocator<T_STRUCT,d>& allocator,const std::pair<const unsigned long*,unsigned>& blocks,unsigned T_STRUCT::* flags_field_input,SPGrid_Set<Flag_Array_Type>& set_input,
                   T_INDEX domain_size_input,T_INDEX origin_input,T_INDEX sub_domain_size_input)
        :flags_field(flags_field_input),set(set_input),origin(origin_input),sub_domain_size(sub_domain_size_input),domain_size(domain_size_input)
    {Run(allocator,blocks);}
    
    void Run(SPGrid_Allocator<T_STRUCT,d>& allocator,const std::pair<const unsigned long*,unsigned>& blocks) const
    {
        Flag_Array_Type flags=allocator.Get_Array(flags_field);    
        T_INDEX number_of_subdomains = Number_of_Subdomains(domain_size,origin,sub_domain_size);
        for(SPGrid_Block_Iterator<T_MASK> iterator(blocks);iterator.Valid();iterator.Next_Block()){
            unsigned long offset=iterator.Offset();
            T_INDEX base_index=iterator.Index().template Cast<T_INDEX>();
            for(RANGE_ITERATOR<d> iterator(RANGE<T_INDEX>(base_index,base_index+allocator.Block_Size().template Cast<T_INDEX>()-1));
                iterator.Valid();
                iterator.Next(),offset+=sizeof(unsigned)){
                unsigned& flag=flags(offset);
                T_INDEX index=iterator.Index()-origin;
                if(flag & SPGrid_Solver_Cell_Type_Active){
                    bool is_interface=false;
                    for(int axis=1;axis<=d;++axis) if(index(axis) % sub_domain_size(axis)==0) is_interface=true;
                    if(is_interface){
                        flag-=SPGrid_Solver_Cell_Type_Active;//Mark this cell no longer Active, but rather Interface
                        flag|=SPGrid_Solver_Cell_Type_Interface;
                        int id=(index(1)/sub_domain_size(1));
                        for(int axis=1;axis<=d;++axis) id=id*number_of_subdomains(axis)+(index(axis)/sub_domain_size(axis));
                        PHYSBAM_ASSERT(id==id&SPGrid_Solver_PartitionID_Mask);
                        flag|=id;}}}}
    }     
};
    //#####################################################################
}
#endif
