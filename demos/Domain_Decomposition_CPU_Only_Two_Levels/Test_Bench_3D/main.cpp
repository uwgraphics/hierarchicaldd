#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <SPGrid/Core/SPGrid_Allocator.h>
#include <SPGrid/Core/SPGrid_Set.h>
#include <SPGrid/Tools/SPGrid_Block_Iterator.h>
#include <SPGrid/Tools/SPGrid_Clear.h>
#include <SPGrid_Fluids/Simulation/FLUIDS_SIMULATION_DATA.h>
#include <Common_Tools/Math_Tools/RANGE_ITERATOR.h>
#include <Threading_Tools/PTHREAD_QUEUE.h>
#include <SPGrid/Tools/SPGrid_Threading_Helper.h>
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Krylov_Solvers/CONJUGATE_GRADIENT.h>
#include <PhysBAM_Tools/Read_Write/Grids_Uniform/READ_WRITE_GRID.h>
#include <PhysBAM_Tools/Read_Write/Grids_Uniform_Arrays/READ_WRITE_ARRAYS.h>
#include <SPGrid_Fluids/Grids/Rasterizers/HIERARCHICAL_RASTERIZER.h>
#include <SPGrid_Fluids/Grids/GRID_HIERARCHY_INITIALIZER.h>
#include <SPGrid_Fluids/Projection/GRID_HIERARCHY_PROJECTION.h>
#include <PhysBAM_Geometry/Geometry_Particles/REGISTER_GEOMETRY_READ_WRITE.h>
#include <PhysBAM_Geometry/Solids_Geometry/DEFORMABLE_GEOMETRY_COLLECTION.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/FREE_PARTICLES.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/SEGMENTED_CURVE.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <PhysBAM_Geometry/Read_Write/Geometry/READ_WRITE_TRIANGULATED_SURFACE.h>
#include <PhysBAM_Tools/Data_Structures/HASHTABLE.h>
#include <PhysBAM_Tools/Data_Structures/HASHTABLE_ITERATOR.h>
#include "../SPGrid_Domain_Decomposition_Solver/CG_SYSTEM.h"
#include "../SPGrid_Domain_Decomposition_Solver/CG_VECTOR.h"
#include "../SPGrid_Domain_Decomposition_Solver/SPGrid_Block_Pair_Weaver.h"
#include "../SPGrid_Domain_Decomposition_Solver/SPGrid_Convertor.h"
#include "../SPGrid_Domain_Decomposition_Solver/SPGRID_DOMAIN_DECOMPOSITION_DATA.h"
#include "../SPGrid_Domain_Decomposition_Solver/SPGrid_Domain_Divider_Helper.h"
#include "../SPGrid_Domain_Decomposition_Solver/SPGrid_Domain_Decomposition_Preconditioner.h"
#include "../SPGrid_Domain_Decomposition_Solver/SPGrid_Flag_Helper.h"
#include "../SPGrid_Domain_Decomposition_Solver/SPGrid_Gemini.h"
#include "../SPGrid_Domain_Decomposition_Solver/SPGrid_Linearized_Data_Copy_Helper.h"
#include "../SPGrid_Domain_Decomposition_Solver/SPGrid_Master_Array_Linearizer.h"
#include "../SPGrid_Domain_Decomposition_Solver/SPGrid_V_Cycle_Topology.h"
#include "../SPGrid_Domain_Decomposition_Solver/SPGrid_Domain_Decomposition_Solver.h"
#include <vector>
#include <omp.h>
using namespace PhysBAM;
using namespace SPGrid;
namespace PhysBAM{int PhysBAM_number_of_threads=0;}

typedef enum{EXTERIOR=0,DIRICHLET,INTERIOR} CELL_TYPE;


template<typename T_STRUCT,typename T,int d>
void Read_SPGrid(STREAM_TYPE stream_type,GRID<VECTOR<T,d> >& grid,GRID_HIERARCHY<T_STRUCT,T,d>*& SPGrid_hierarchy,T T_STRUCT::* velocity_u_field,T T_STRUCT::* velocity_v_field,T T_STRUCT::* velocity_w_field,const std::string input_directory,const int frame)
{
    std::string f=STRING_UTILITIES::string_sprintf("%d",frame);
    int levels=0;
    FILE_UTILITIES::Read_From_File(stream_type,input_directory+"common/fine_grid",grid);
    LOG::cout<<"Grid Read: "<<grid<<std::endl;
    FILE_UTILITIES::template Read_From_Text_File<int>(input_directory+"/"+f+"/levels",levels);
    //PHYSBAM_ASSERT(levels==1);
    if(SPGrid_hierarchy) delete SPGrid_hierarchy;
    SPGrid_hierarchy=new GRID_HIERARCHY<T_STRUCT,T,d>(grid,levels);
    SPGrid_hierarchy->Initialize_Sets();
    SPGrid_hierarchy->Read_Block_Offsets(input_directory+"/"+f+"/block_offsets");
    SPGrid_hierarchy->Read_Flags_Channel(input_directory+"/"+f+"/flags"); 
    SPGrid_hierarchy->Read_Data_Channel(input_directory+"/"+f+"/spgrid_u",velocity_u_field);
    SPGrid_hierarchy->Read_Data_Channel(input_directory+"/"+f+"/spgrid_v",velocity_v_field);
    if(d == 3) SPGrid_hierarchy->Read_Data_Channel(input_directory+"/"+f+"/spgrid_w",velocity_w_field);
}

template<typename T_STRUCT,typename T,int d>
void Read_SPGrid(STREAM_TYPE stream_type,GRID<VECTOR<T,d> >& grid,GRID_HIERARCHY<T_STRUCT,T,d>*& SPGrid_hierarchy,T T_STRUCT::* b_field,const std::string input_directory,const int frame)
{
    std::string f=STRING_UTILITIES::string_sprintf("%d",frame);
    int levels=0;
    FILE_UTILITIES::Read_From_File(stream_type,input_directory+"common/fine_grid",grid);
    FILE_UTILITIES::template Read_From_Text_File<int>(input_directory+"/"+f+"/levels",levels);
    //PHYSBAM_ASSERT(levels==1);
    if(SPGrid_hierarchy) delete SPGrid_hierarchy;
    SPGrid_hierarchy=new GRID_HIERARCHY<T_STRUCT,T,d>(grid,levels);
    SPGrid_hierarchy->Initialize_Sets();
    SPGrid_hierarchy->Read_Block_Offsets(input_directory+"/"+f+"/block_offsets");
    SPGrid_hierarchy->Read_Flags_Channel(input_directory+"/"+f+"/flags"); 
    //SPGrid_hierarchy->Read_Data_Channel(input_directory+"/"+f+"/b_field",b_field);
    SPGrid_hierarchy->Read_Data_Channel(input_directory+"/"+f+"/residual",b_field);
}

template<typename T_STRUCT,typename T,int d>
void Write_SPGrid(STREAM_TYPE stream_type,GRID<VECTOR<T,d> >& grid,GRID_HIERARCHY<T_STRUCT,T,d>*& SPGrid_hierarchy,T T_STRUCT::* density_field,const std::string output_directory,const int frame)
{
    std::string f=STRING_UTILITIES::string_sprintf("%d",frame);
    FILE_UTILITIES::Create_Directory(output_directory);
    FILE_UTILITIES::Create_Directory(output_directory+"/common");
    FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/common/fine_grid",grid);
    FILE_UTILITIES::Create_Directory(output_directory+"/"+f);
    FILE_UTILITIES::Write_To_Text_File(output_directory+"/"+f+"/levels",STRING_UTILITIES::string_sprintf("%d",SPGrid_hierarchy->Levels()));
    SPGrid_hierarchy->Write_Flags_Channel(output_directory+"/"+f+"/flags");
    SPGrid_hierarchy->Write_Block_Offsets(output_directory+"/"+f+"/block_offsets");
    SPGrid_hierarchy->Write_Data_Channel(output_directory+"/"+f+"/spgrid_density",density_field);
}
template<class T_STRUCT, class T,int d>
class Norm_Helper{
    typedef VECTOR<int,d> T_INDEX;
    typedef SPGrid_Allocator<T_STRUCT,d> SPG_Allocator;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<const unsigned>::type SPG_Flags_Array_Type;
    typedef SPGrid_Set<SPG_Flags_Array_Type> SPG_Set_Type;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<const T>::type Const_data_array_type;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<const T>::mask T_MASK;

public:
    static T L1_Norm(const SPG_Allocator& allocator, const std::pair<const unsigned long*,unsigned>& blocks,T T_STRUCT::* r_field){
        T result = 0.f;
        Const_data_array_type r=allocator.Get_Const_Array(r_field);
        for(SPGrid_Block_Iterator<T_MASK> iterator(blocks);iterator.Valid();iterator.Next_Block()){
            unsigned long offset=iterator.Offset();
            T_INDEX base_index=iterator.Index().template Cast<T_INDEX>();
            for(RANGE_ITERATOR<d> iterator(RANGE<T_INDEX>(base_index,base_index+allocator.Block_Size().template Cast<T_INDEX>()-1));
                iterator.Valid();
                iterator.Next(),offset+=sizeof(T)){
                if(fabs(r(offset)) > result) result = fabs(r(offset));}}
        return result;
    }
};
template<class T,class T_STRUCT,int d>
void Write_Output(STREAM_TYPE stream_type,SPGrid_Set<typename SPGrid_Allocator<T_STRUCT,d>::template Array<unsigned>::type>& set,
                  GRID<VECTOR<T,d> >& grid,std::string output_directory,int frame,
                  SPGrid_Allocator<T_STRUCT,d>& flags_allocator,unsigned T_STRUCT::* flags_field,
                  SPGrid_Allocator<T_STRUCT,d>& u_allocator,T T_STRUCT::* u_field,
                  SPGrid_Allocator<T_STRUCT,d>& b_allocator,T T_STRUCT::* b_field,
                  SPGrid_Allocator<T_STRUCT,d>& r_allocator,T T_STRUCT::* r_field)
{
    typedef VECTOR<int,d> T_INDEX;
    typedef SPGrid_Allocator<T_STRUCT,d> SPG_Allocator;
    typedef typename SPG_Allocator::Array<const unsigned>::type SPG_Const_Flags_Array_Type;
    typedef typename SPG_Allocator::Array<const T>::type SPG_Const_Data_Array_Type;
    typedef typename SPG_Allocator::Array<T>::mask T_MASK;

    RANGE<T_INDEX> range=grid.Cell_Indices();
    ARRAY<T,T_INDEX> u_field_array(range);
    ARRAY<T,T_INDEX> b_field_array(range);
    ARRAY<T,T_INDEX> r_field_array(range);
    ARRAY<T,T_INDEX> face_count_array(range);
    ARRAY<bool,T_INDEX> psi_A(range);
    ARRAY<bool,T_INDEX> psi_D(range);
    ARRAY<bool,T_INDEX> psi_B(range);
    ARRAY<bool,T_INDEX> psi_I(range);
    ARRAY<bool,T_INDEX> psi_E(range);
    psi_E.Fill(true);

    SPG_Const_Flags_Array_Type flags = flags_allocator.Get_Const_Array(flags_field);
    SPG_Const_Data_Array_Type u_array = u_allocator.Get_Const_Array(u_field);
    SPG_Const_Data_Array_Type b_array = u_allocator.Get_Const_Array(b_field);
    SPG_Const_Data_Array_Type r_array = u_allocator.Get_Const_Array(r_field);

    for(SPGrid_Block_Iterator<T_MASK> iterator(set.Get_Blocks());iterator.Valid();iterator.Next_Block()){
        unsigned long offset=iterator.Offset();
        T_INDEX base_index=iterator.Index().template Cast<T_INDEX>();

        for(RANGE_ITERATOR<d> iterator(RANGE<T_INDEX>(base_index,base_index+flags_allocator.Block_Size().template Cast<T_INDEX>()-1));
            iterator.Valid();
            iterator.Next(),offset+=sizeof(unsigned)){
            T_INDEX index=iterator.Index()+1;//plus one here is because grid array starts from 1, SPGrid starts at 0
            unsigned flag = flags(offset);
            T u = u_array(offset);
            T r = r_array(offset);
            T b = b_array(offset);
            
            if(flag) PHYSBAM_ASSERT(range.Lazy_Inside(index));
            if(!flag) continue;
            if(range.Lazy_Inside(index)) psi_E(index) = false;
            if(flag & SPGrid_Solver_Cell_Type_Active)    psi_A(index) = true;
            if(flag & SPGrid_Solver_Cell_Type_Dirichlet) psi_D(index) = true;
            if(flag & SPGrid_Solver_Cell_Type_Boundary)  psi_B(index) = true;
            if(flag & SPGrid_Solver_Cell_Type_Interface) psi_I(index) = true;

            if(flag & (SPGrid_Solver_Cell_Type_Active | SPGrid_Solver_Cell_Type_Interface) ){
                for(int i = 0;i < 6;++i)
                    if(flag & (SPGrid_Solver_Face_Minus_X_Active << i)){
                        face_count_array(index) += 1.f/6.f;}}
            if(flag){
                u_field_array(index) = u;
                b_field_array(index) = b;
                r_field_array(index) = fabs(r);
            }else{
                PHYSBAM_ASSERT(u == 0.f);
                PHYSBAM_ASSERT(r == 0.f);}
        }
    }
    std::string f=STRING_UTILITIES::string_sprintf("%d",frame);
    FILE_UTILITIES::Create_Directory(output_directory);
    FILE_UTILITIES::Create_Directory(output_directory+"/common");
    FILE_UTILITIES::Create_Directory(output_directory+"/"+f);
    FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/common/grid",grid);
    FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/u_field",u_field_array);
    FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/b_field",b_field_array);
    FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/r_field",r_field_array);
    FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/face_count_field",face_count_array);
    FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/psi_D",psi_D);
    FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/psi_I",psi_I);
    FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/psi_A",psi_A);
    FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/psi_B",psi_B);
    FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/psi_E",psi_E);
}
int main(int argc,char *argv[])
{
    enum{d=3};
    typedef float T;
    typedef int INDEX;
    typedef FLUIDS_SIMULATION_DATA<T> T_STRUCT_SIMULATION;
    //typedef SPGRID_MULTIGRID_DATA<T> T_STRUCT;
    unsigned T_STRUCT_SIMULATION::* flags_field=&T_STRUCT_SIMULATION::flags;
    T T_STRUCT_SIMULATION::* velocity_u_field=&T_STRUCT_SIMULATION::ch0;
    T T_STRUCT_SIMULATION::* velocity_v_field=&T_STRUCT_SIMULATION::ch1;
    T T_STRUCT_SIMULATION::* velocity_w_field=&T_STRUCT_SIMULATION::ch2;
    T T_STRUCT_SIMULATION::* b_field=&T_STRUCT_SIMULATION::ch3;
    T T_STRUCT_SIMULATION::* u_field=&T_STRUCT_SIMULATION::ch4;
    T T_STRUCT_SIMULATION::* r_field=&T_STRUCT_SIMULATION::ch5;
    T T_STRUCT_SIMULATION::* tmp_field=&T_STRUCT_SIMULATION::ch6;

    typedef SPGrid_Allocator<T_STRUCT_SIMULATION,d> SPG_Allocator_Simulation;
    typedef SPG_Allocator_Simulation::Array<T>::type SPG_Simulation_Data_Array_Type;
    typedef SPG_Allocator_Simulation::Array<unsigned>::type SPG_Simulation_Flags_Array_Type;
    typedef SPG_Allocator_Simulation::Array<const unsigned>::type SPG_Const_Simulation_Flags_Array_Type;
    typedef SPG_Allocator_Simulation::Array<T>::mask T_SIMULATION_MASK;
    typedef SPGrid_Set<SPG_Simulation_Flags_Array_Type> SPG_Simulation_Set_Type;

    enum{simulation_block_xsize=1<<T_SIMULATION_MASK::block_xbits,
         simulation_block_ysize=1<<T_SIMULATION_MASK::block_ybits,
         simulation_block_zsize=1<<T_SIMULATION_MASK::block_zbits};

    typedef VECTOR<int,d> T_INDEX;
    typedef VECTOR<T,d> TV;

    LOG::Initialize_Logging();

    Initialize_Geometry_Particle();
    Initialize_Read_Write_Structures();

    PARSE_ARGS parse_args;
    parse_args.Add_Vector_3D_Argument("-size",VECTOR<double,3>(257,257,257),"n n n","Size of 3D grid");
    parse_args.Add_Integer_Argument("-threads",0,"n","Number of pThreads (0 means serial)");
    parse_args.Add_Integer_Argument("-vlevels",3,"n","Number of V cycle levels");
    parse_args.Add_Integer_Argument("-sigma_levels",3,"n","Number of sigma MG levels");
    parse_args.Add_Integer_Argument("-sigma_smoothing",3,"n","Number of sigma smoothing iteration");
    parse_args.Add_Integer_Argument("-sigma_mg",1,"n","Number of sigma MG iteration");
    parse_args.Add_Integer_Argument("-interior_smoothing",1,"n","Number of interior smoothing iteration");
    parse_args.Add_Integer_Argument("-boundary_smoothing",3,"n","Number of boundary smoothing iteration");
    parse_args.Add_Integer_Argument("-nVcycles",1,"n","Number of V cycles in DD");
    parse_args.Add_Integer_Argument("-frame",0,"n","The frame we are loading");
    parse_args.Add_Vector_3D_Argument("-origin",VECTOR<double,3>(0,0,0),"n n n","Origin of The domain");
    parse_args.Add_Vector_3D_Argument("-sub_size",VECTOR<double,3>(64,64,64),"n n n","Size of the sub-domains");
    parse_args.Add_String_Argument("-d","","Directory");
    parse_args.Add_String_Argument("-o","output","Output directory");
    parse_args.Add_String_Argument("-mode","velocity","velocity, rhs, or no_convertion");
    parse_args.Add_String_Argument("-multilevel","no","yes or no");
    parse_args.Parse(argc,argv);
    
    std::string input_directory=parse_args.Get_String_Value("-d");
    std::string output_directory=parse_args.Get_String_Value("-o");
    std::string mode=parse_args.Get_String_Value("-mode");
    bool multilevel=(parse_args.Get_String_Value("-multilevel")=="yes");
    int frame=parse_args.Get_Integer_Value("-frame");

    std_array<int,d> size;

    for(int v=1;v<=d;v++) size(v-1)=static_cast<unsigned int>(parse_args.Get_Vector_3D_Value("-size")(v));
    LOG::cout<<"Size of 3D grid          : "<<size<<std::endl;

    T_INDEX subdomain_size;
    for(int v=1;v<=d;v++) subdomain_size(v)=static_cast<unsigned int>(parse_args.Get_Vector_3D_Value("-sub_size")(v));
    LOG::cout<<"Size of subdomain        : "<<subdomain_size<<std::endl;

    T_INDEX origin;
    for(int v=1;v<=d;v++) origin(v)=static_cast<unsigned int>(parse_args.Get_Vector_3D_Value("-origin")(v));
    LOG::cout<<"Origin                   : "<<origin<<std::endl;

    int number_of_threads=parse_args.Get_Integer_Value("-threads");
    if(number_of_threads){
        pthread_queue=new PTHREAD_QUEUE(number_of_threads);
        LOG::cout<<"pThreads             : "<<number_of_threads<<std::endl;}
    else
        LOG::cout<<"pThreads             : None (serial)"<<std::endl;

    if(number_of_threads)
        omp_set_num_threads(number_of_threads);

    int v_cycle_levels=parse_args.Get_Integer_Value("-vlevels");
    LOG::cout<<"V cycle levels           : "<<v_cycle_levels<<std::endl;

    int sigma_levels=parse_args.Get_Integer_Value("-sigma_levels");
    LOG::cout<<"Sigma MG levels          : "<<sigma_levels<<std::endl;

    int sigma_smoothing=parse_args.Get_Integer_Value("-sigma_smoothing");
    LOG::cout<<"Sigma smoothing iteration: "<<sigma_smoothing<<std::endl;

    int n_sigma_mg=parse_args.Get_Integer_Value("-sigma_mg");
    LOG::cout<<"Sigma MG iteration: "<<n_sigma_mg<<std::endl;

    int interior_smoothing=parse_args.Get_Integer_Value("-interior_smoothing");
    LOG::cout<<"Interior smoothing iteration: "<<interior_smoothing<<std::endl;

    int boundary_smoothing=parse_args.Get_Integer_Value("-boundary_smoothing");
    LOG::cout<<"boundary smoothing iteration: "<<boundary_smoothing<<std::endl;

    int nVcycles=parse_args.Get_Integer_Value("-nVcycles");
    LOG::cout<<"V cycle iterations in DD : "<<nVcycles<<std::endl;

    typedef float RW;
    STREAM_TYPE stream_type((RW()));
    GRID<TV> grid;
    GRID_HIERARCHY<T_STRUCT_SIMULATION,T,d>* SPGrid_hierarchy=NULL;
    if(mode=="no_convertion"){
        typedef SPGRID_DOMAIN_DECOMPOSITION_DATA<T> T_STRUCT_4;
        typedef unsigned int T_offset_ptr;
    
        SPGrid_Domain_Decomposition_Solver<T,T_STRUCT_4,T_STRUCT_SIMULATION,d,T_offset_ptr> solver(number_of_threads);
        solver.Read_SPGrid(input_directory);
        LOG::cout<<"Finished Reading SPGrid"<<std::endl;
        T norm_rhs=Norm_Helper<T_STRUCT_4,T,d>::L1_Norm(*(solver.spgrid_gemini.allocator_second),(*(solver.spgrid_gemini.set)).Get_Blocks(),&T_STRUCT_4::ch1);
        LOG::cout<<"Norm of the RHS: "<<norm_rhs<<std::endl;
        solver.Initialize_Preconditioner(v_cycle_levels,size.Cast<T_INDEX>(),subdomain_size,sigma_levels,nVcycles,interior_smoothing,boundary_smoothing,sigma_smoothing,n_sigma_mg);
        solver.Solve();
    }else{
        if(mode=="velocity"){
            Read_SPGrid<T_STRUCT_SIMULATION,T,d>(stream_type,grid,SPGrid_hierarchy,velocity_u_field,velocity_v_field,velocity_w_field,input_directory,frame);
            {LOG::SCOPE scope("Compute Divergence");
                GRID_HIERARCHY_PROJECTION<T_STRUCT_SIMULATION,T,d>::Compute_Divergence(*SPGrid_hierarchy,flags_field,velocity_u_field,velocity_v_field,velocity_w_field,b_field);}}
        else if(mode=="rhs"){
            Read_SPGrid<T_STRUCT_SIMULATION,T,d>(stream_type,grid,SPGrid_hierarchy,b_field,input_directory,frame);
        }
        //Write_SPGrid(stream_type,grid,SPGrid_hierarchy,b_field,output_directory,0);

        SPG_Allocator_Simulation& allocator_simulation=SPGrid_hierarchy->Allocator(1);
        SPG_Simulation_Set_Type& set_simulation=SPGrid_hierarchy->Set(1);

        LOG::cout<<"Padded size: "<<allocator_simulation.Padded_Size()<<std::endl;
    
        T norm_rhs=Norm_Helper<T_STRUCT_SIMULATION,T,d>::L1_Norm(allocator_simulation,set_simulation.Get_Blocks(),b_field);
        LOG::cout<<"Norm of the RHS: "<<norm_rhs<<std::endl;

        typedef SPGRID_DOMAIN_DECOMPOSITION_DATA<T> T_STRUCT_4;
        typedef unsigned int T_offset_ptr;
    
        SPGrid_Domain_Decomposition_Solver<T,T_STRUCT_4,T_STRUCT_SIMULATION,d,T_offset_ptr> solver(number_of_threads);
        solver.Convert(allocator_simulation,set_simulation,b_field,flags_field);
    
        delete SPGrid_hierarchy;
        
        solver.Initialize_Preconditioner(v_cycle_levels,size.Cast<T_INDEX>(),subdomain_size,sigma_levels,nVcycles,interior_smoothing,boundary_smoothing,sigma_smoothing,n_sigma_mg);
        solver.Solve();
    }
    LOG::Finish_Logging();
    return 0;
}
