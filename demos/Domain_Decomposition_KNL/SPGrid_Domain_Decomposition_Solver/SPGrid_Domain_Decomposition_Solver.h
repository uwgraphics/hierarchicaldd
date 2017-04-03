//#####################################################################
//Removed For Anonymity: Copyright Authors of Submission pap173s1
//Removed For Anonymity: Copyright Authors of Submission pap173s1
//#####################################################################
// Class SPGrid_Domain_Decomposition_Solver
//#####################################################################
#ifndef __SPGRID_DOMAIN_DECOMPOSITION_SOLVER_H__
#define __SPGRID_DOMAIN_DECOMPOSITION_SOLVER_H__

#include <SPGrid/Core/SPGrid_Allocator.h>
#include <SPGrid/Core/SPGrid_Set.h>

#include "../Interface_Solver/INTERFACE_MATRIX_GENERATOR.h"
#include "CG_SYSTEM.h"
#include "CG_VECTOR.h"
#include "SPGrid_Convertor.h"
#include "SPGrid_Domain_Decomposition_Preconditioner.h"
#include "SPGrid_Gemini.h"
#include "SPGrid_Subdomain_Splitter.h"
#include "SPGRID_READ_WRITE.h"

namespace SPGrid{
template<typename T,typename T_STRUCT_4,typename T_STRUCT_16,int d,typename T_offset_ptr>
class SPGrid_Domain_Decomposition_Solver{
public:
    template<typename T_FIELD>
    struct Field{
        T_FIELD T_STRUCT_4::* field;
        SPGrid_Allocator<T_STRUCT_4,d>* allocator;
    };
    // All the field for the DDPCG
    Field<T> x_field,b_field,q_field,s_field,z_field,r_field,k_field;
    Field<unsigned> flags_field;
    typedef PhysBAM::VECTOR<int,d> T_INDEX;
    int number_of_threads;
    SPGrid_Domain_Decomposition_Preconditioner<T,T_STRUCT_4,d,T_offset_ptr>* dd_preconditioner;
    typedef SPGrid_Master_Array_Linearizer<T,NextLogTwo<sizeof(T_STRUCT_4)>::value,d,T_offset_ptr> T_Linearizer;
    std::vector<std::vector<T_Linearizer>*> linearizer_hierarchy;
    SPGrid_Gemini<T_STRUCT_4,d> spgrid_gemini;
    INTERFACE_MATRIX_GENERATOR<T,d> interface_matrix_generator;
    SPGrid_Domain_Decomposition_Solver(int number_of_threads_input=0):number_of_threads(number_of_threads_input),dd_preconditioner(0){}
    ~SPGrid_Domain_Decomposition_Solver(){
        if(dd_preconditioner) delete dd_preconditioner;
        for(int i=0;i<linearizer_hierarchy.size();++i) delete linearizer_hierarchy[i];};
    //#####################################################################
    // Function Clear_Non_Active_Values
    //#####################################################################
    template<typename T_STRUCT>
    static void Clear_Non_Active_Values(SPGrid_Allocator<T_STRUCT,d>& flag_allocator,SPGrid_Allocator<T_STRUCT,d>& data_allocator,SPGrid_Set<typename SPGrid_Allocator<T_STRUCT,d>::template Array<unsigned>::type>& set,T T_STRUCT::* data_field,unsigned T_STRUCT::* flags_field){
        typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<const unsigned>::type Const_Flag_Array_Type;
        typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<T>::type Data_Array_Type;
        typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<unsigned>::mask T_MASK;
        Data_Array_Type data_array = data_allocator.Get_Array(data_field);
        Const_Flag_Array_Type flags_array = flag_allocator.Get_Const_Array(flags_field);
        
        for(SPGrid_Block_Iterator<T_MASK> iterator(set.Get_Blocks());iterator.Valid();iterator.Next_Block()){
            unsigned long offset=iterator.Offset();
            T_INDEX base_index=iterator.Index().template Cast<T_INDEX>();
            for(RANGE_ITERATOR<d> iterator(RANGE<T_INDEX>(base_index,base_index+flag_allocator.Block_Size().template Cast<T_INDEX>()-1));
                iterator.Valid();
                iterator.Next(),offset+=sizeof(unsigned)){
                const unsigned& flag = flags_array(offset);
                T& data = data_array(offset);
                if(!(flag & (SPGrid_Solver_Cell_Type_Active | SPGrid_Solver_Cell_Type_Interface))) {
                    //if(data!=T(0)) std::cerr << "Found non-zero in non active cells!"<<std::endl;
                    data = T(0);
                }else{
                    //data = T(1);
                }
            }}
    }
    //#####################################################################
    // Function Convert
    //#####################################################################
    void Read_SPGrid(const std::string& base_dir){
        SPGRID_READ_WRITE<T_STRUCT_4,T,d>::Read_SPGrid(base_dir,spgrid_gemini.allocator_second,spgrid_gemini.allocator_first,spgrid_gemini.set,&T_STRUCT_4::ch1);
        x_field.allocator = spgrid_gemini.allocator_second;x_field.field=&T_STRUCT_4::ch0;
        b_field.allocator = spgrid_gemini.allocator_second;b_field.field=&T_STRUCT_4::ch1;
        q_field.allocator = spgrid_gemini.allocator_second;q_field.field=&T_STRUCT_4::ch2;
        s_field.allocator = spgrid_gemini.allocator_second;s_field.field=&T_STRUCT_4::ch3;
        z_field.allocator = spgrid_gemini.allocator_first;z_field.field=&T_STRUCT_4::ch0;
        r_field.allocator = spgrid_gemini.allocator_first;r_field.field=&T_STRUCT_4::ch1;
        k_field.allocator = spgrid_gemini.allocator_first;k_field.field=&T_STRUCT_4::ch2;
        flags_field.allocator = spgrid_gemini.allocator_first;flags_field.field=&T_STRUCT_4::flags;
    }
    //#####################################################################
    // Function Convert
    //#####################################################################
    void Convert(SPGrid_Allocator<T_STRUCT_16,d>& allocator_16,SPGrid_Set<typename SPGrid_Allocator<T_STRUCT_16,d>::template Array<unsigned>::type>& set_16,T T_STRUCT_16::* b_field_16,unsigned T_STRUCT_16::* flags_field_16){
        SPGrid_Computations::SPGrid_Creator_16_to_4<T_STRUCT_16,T_STRUCT_4,unsigned,d>::Create(spgrid_gemini,allocator_16,set_16);        
        x_field.allocator = spgrid_gemini.allocator_second;x_field.field=&T_STRUCT_4::ch0;
        b_field.allocator = spgrid_gemini.allocator_second;b_field.field=&T_STRUCT_4::ch1;
        q_field.allocator = spgrid_gemini.allocator_second;q_field.field=&T_STRUCT_4::ch2;
        s_field.allocator = spgrid_gemini.allocator_second;s_field.field=&T_STRUCT_4::ch3;
        z_field.allocator = spgrid_gemini.allocator_first;z_field.field=&T_STRUCT_4::ch0;
        r_field.allocator = spgrid_gemini.allocator_first;r_field.field=&T_STRUCT_4::ch1;
        k_field.allocator = spgrid_gemini.allocator_first;k_field.field=&T_STRUCT_4::ch2;
        flags_field.allocator = spgrid_gemini.allocator_first;flags_field.field=&T_STRUCT_4::flags;        
        if(number_of_threads)
            SPGrid_Computations::Threading_Helper<T_STRUCT_4,d>(*(b_field.allocator),spgrid_gemini.Set().Get_Blocks()).Run_Parallel(SPGrid_Computations::SPGrid_Convertor_16_To_4<T_STRUCT_16,T_STRUCT_4,T,d>(allocator_16,set_16,b_field_16,b_field.field),number_of_threads);
        else
            SPGrid_Computations::SPGrid_Convertor_16_To_4<T_STRUCT_16,T_STRUCT_4,T,d>(*(b_field.allocator),spgrid_gemini.Set().Get_Blocks(),allocator_16,set_16,b_field_16,b_field.field);
        if(number_of_threads)
            SPGrid_Computations::Threading_Helper<T_STRUCT_4,d>(*(flags_field.allocator),spgrid_gemini.Set().Get_Blocks()).Run_Parallel(SPGrid_Computations::SPGrid_Flag_Convertor_16_To_4<T_STRUCT_16,T_STRUCT_4,d>(allocator_16,set_16,flags_field_16,flags_field.field),number_of_threads);
        else
            SPGrid_Computations::SPGrid_Flag_Convertor_16_To_4<T_STRUCT_16,T_STRUCT_4,d>(*(flags_field.allocator),spgrid_gemini.Set().Get_Blocks(),allocator_16,set_16,flags_field_16,flags_field.field);    
    }
    void Copy(SPGrid_Allocator<T_STRUCT_16,d>& allocator_16,SPGrid_Set<typename SPGrid_Allocator<T_STRUCT_16,d>::template Array<unsigned>::type>& set_16,T T_STRUCT_4::* data_field_4,T T_STRUCT_16::* data_field_16){
        if(number_of_threads)
            SPGrid_Computations::Threading_Helper<T_STRUCT_4,d>(spgrid_gemini.First_Allocator(),spgrid_gemini.Set().Get_Blocks()).Run_Parallel(SPGrid_Computations::SPGrid_Convertor_16_To_4<T_STRUCT_16,T_STRUCT_4,T,d>(allocator_16,set_16,data_field_16,data_field_4),number_of_threads);
        else
            SPGrid_Computations::SPGrid_Convertor_16_To_4<T_STRUCT_16,T_STRUCT_4,T,d>(spgrid_gemini.First_Allocator(),spgrid_gemini.Set().Get_Blocks(),allocator_16,set_16,data_field_16,data_field_4);
    }
    void Copy_Back(SPGrid_Allocator<T_STRUCT_16,d>& allocator_16,SPGrid_Set<typename SPGrid_Allocator<T_STRUCT_16,d>::template Array<unsigned>::type>& set_16,T T_STRUCT_4::* data_field_4,T T_STRUCT_16::* data_field_16){
        if(number_of_threads)
            SPGrid_Computations::Threading_Helper<T_STRUCT_4,d>(spgrid_gemini.First_Allocator(),spgrid_gemini.Set().Get_Blocks()).Run_Parallel(SPGrid_Computations::SPGrid_Convertor_4_To_16<T_STRUCT_16,T_STRUCT_4,T,d>(allocator_16,set_16,data_field_16,data_field_4),number_of_threads);
        else
            SPGrid_Computations::SPGrid_Convertor_4_To_16<T_STRUCT_16,T_STRUCT_4,T,d>(spgrid_gemini.First_Allocator(),spgrid_gemini.Set().Get_Blocks(),allocator_16,set_16,data_field_16,data_field_4);
    }
    //#####################################################################
    // Function Initialize_Preconditioner
    //#####################################################################
    void Initialize_Preconditioner(int v_cycle_levels,T_INDEX size,T_INDEX subdomain_size,int sigma_levels,int nVcycles,int interior_smoothing,int boundary_smoothing,int sigma_smoothing,int n_sigma_mg,std::string iomode="output",std::string dump_directory="matrices_dump"){
        SPGrid_Subdomain_Splitter<T,T_STRUCT_4,d,T_offset_ptr> splitter(*(flags_field.allocator),spgrid_gemini.Set(),size,subdomain_size);
        splitter.Split(b_field.field,interface_matrix_generator,linearizer_hierarchy,v_cycle_levels,sigma_levels,number_of_threads,iomode,dump_directory);

        if(number_of_threads)
            SPGrid_Computations::Threading_Helper<T_STRUCT_4,d>(*(flags_field.allocator),spgrid_gemini.Set().Get_Blocks()).Run_Parallel(SPGrid_Computations::Domain_Divider<T_STRUCT_4,d>(&T_STRUCT_4::flags,spgrid_gemini.Set(),size,T_INDEX(),subdomain_size),number_of_threads);
        else
            SPGrid_Computations::Domain_Divider<T_STRUCT_4,d>(*(flags_field.allocator),
                                                              spgrid_gemini.Set().Get_Blocks(),
                                                              &T_STRUCT_4::flags,spgrid_gemini.Set(),
                                                              size,
                                                              T_INDEX(),subdomain_size);

        if(number_of_threads)
            SPGrid_Computations::Threading_Helper<T_STRUCT_4,d>(*(flags_field.allocator),spgrid_gemini.Set().Get_Blocks()).Run_Parallel(SPGrid_Computations::Flaging<T_STRUCT_4,d>(flags_field.field,spgrid_gemini.Set()),number_of_threads);
        else
            SPGrid_Computations::Flaging<T_STRUCT_4,d>(*(flags_field.allocator),spgrid_gemini.Set().Get_Blocks(),flags_field.field,spgrid_gemini.Set());
        
        dd_preconditioner = new SPGrid_Domain_Decomposition_Preconditioner<T,T_STRUCT_4,d,T_offset_ptr>(*(flags_field.allocator),spgrid_gemini.Set(),
                                                                                                        interface_matrix_generator,linearizer_hierarchy,
                                                                                                        subdomain_size,size,flags_field.field,sigma_levels,
                                                                                                        nVcycles,interior_smoothing,boundary_smoothing,
                                                                                                        sigma_smoothing,n_sigma_mg,number_of_threads);
        
        if(iomode=="input"){
            dd_preconditioner->matrices.Read_From_File(dump_directory);
            dd_preconditioner->interface_matrix_generator.Read_Interface_Map(dump_directory);}

        dd_preconditioner->Initialize(true,false,false,iomode=="input");

        if(iomode=="output"){
            FILE_UTILITIES::Create_Directory(dump_directory);
            dd_preconditioner->matrices.Write_To_File(dump_directory);
            dd_preconditioner->interface_matrix_generator.Write_Interface_Map(dump_directory);}
    }
    //#####################################################################
    // Function Solve
    //#####################################################################
    void Solve(){
        LOG::cout<<"Start Solving"<<std::endl;
        PhysBAM::CG_VECTOR<T_STRUCT_4,T,d> x_V(*x_field.allocator,*spgrid_gemini.set,x_field.field,number_of_threads);
        PhysBAM::CG_VECTOR<T_STRUCT_4,T,d> b_V(*b_field.allocator,*spgrid_gemini.set,b_field.field,number_of_threads);
        PhysBAM::CG_VECTOR<T_STRUCT_4,T,d> q_V(*q_field.allocator,*spgrid_gemini.set,q_field.field,number_of_threads);
        PhysBAM::CG_VECTOR<T_STRUCT_4,T,d> s_V(*s_field.allocator,*spgrid_gemini.set,s_field.field,number_of_threads);
        PhysBAM::CG_VECTOR<T_STRUCT_4,T,d> z_V(*z_field.allocator,*spgrid_gemini.set,z_field.field,number_of_threads);
        PhysBAM::CG_VECTOR<T_STRUCT_4,T,d> r_V(*r_field.allocator,*spgrid_gemini.set,r_field.field,number_of_threads);
        PhysBAM::CG_VECTOR<T_STRUCT_4,T,d> k_V(*k_field.allocator,*spgrid_gemini.set,k_field.field,number_of_threads);

        PhysBAM::CG_SYSTEM<T_STRUCT_4,T,d,T_offset_ptr> cg_system(spgrid_gemini.First_Allocator(),&T_STRUCT_4::flags,*dd_preconditioner,&T_STRUCT_4::ch2,number_of_threads);
        cg_system.use_preconditioner=true;

        PhysBAM::CONJUGATE_GRADIENT<T> cg;
        cg.print_residuals=true;
        cg.print_diagnostics=true;
        {LOG::SCOPE scope("solve");
        cg.Solve(cg_system,
                 x_V,
                 b_V,
                 q_V,
                 s_V,
                 r_V,
                 k_V,
                 z_V,
                 1e-5,0,5000);
        }
    }
};
};
#endif
