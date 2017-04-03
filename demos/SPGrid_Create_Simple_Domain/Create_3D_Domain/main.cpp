#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <SPGrid/Core/SPGrid_Allocator.h>
#include <SPGrid/Core/SPGrid_Set.h>
#include <SPGrid/Tools/SPGrid_Block_Iterator.h>
#include <SPGrid_Fluids/Solvers/Domain_Decomposition/SPGRID_MULTIGRID_FLAGS.h>
#include <SPGrid_Fluids/Solvers/Domain_Decomposition/SPGRID_MULTIGRID_DATA.h>
#include <SPGrid_Fluids/Simulation/FLUIDS_SIMULATION_DATA.h>
#include "../HIERARCHICAL_RANGE_ITERATOR.h"
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Read_Write/Grids_Uniform/READ_WRITE_GRID.h>
#include <PhysBAM_Tools/Read_Write/Grids_Uniform_Arrays/READ_WRITE_ARRAYS.h>
#include <SPGrid_Fluids/Grids/Rasterizers/HIERARCHICAL_RASTERIZER.h>
#include <SPGrid_Fluids/Grids/GRID_HIERARCHY_INITIALIZER.h>
#include <PhysBAM_Tools/Krylov_Solvers/CONJUGATE_GRADIENT.h>
#include <PhysBAM_Geometry/Geometry_Particles/REGISTER_GEOMETRY_READ_WRITE.h>
#include <PhysBAM_Geometry/Solids_Geometry/DEFORMABLE_GEOMETRY_COLLECTION.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/FREE_PARTICLES.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/SEGMENTED_CURVE.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <PhysBAM_Geometry/Read_Write/Geometry/READ_WRITE_TRIANGULATED_SURFACE.h>
#include <PhysBAM_Tools/Data_Structures/HASHTABLE.h>
#include <PhysBAM_Tools/Data_Structures/HASHTABLE_ITERATOR.h>
#include "../SPGRID_READ_WRITE.h"
#include "../SPGrid_Flag_Helper.h"

using namespace PhysBAM;
using namespace SPGrid;
namespace PhysBAM{
int PhysBAM_number_of_threads=0;
}
typedef float T;
//typedef FLUIDS_SIMULATION_DATA<T> T_STRUCT;
typedef SPGRID_MULTIGRID_DATA<T> T_STRUCT;
unsigned T_STRUCT::* flags_field=&T_STRUCT::flags;
float T_STRUCT::* u_field=&T_STRUCT::ch0;
float T_STRUCT::* b_field=&T_STRUCT::ch1;
float T_STRUCT::* r_field=&T_STRUCT::ch2;

template<class T>
T PseudoRandom(T& x)
{
    x+=(pi-3.f);
    if(x>1.f) x-=1.f;
    return x;
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
                if(fabs(r(offset)) > result) result = fabs(r(offset));
            }
        }
        return result;
    }
};
template<class T_STRUCT, class T,int d>
class Rasterizer
{
    typedef SPGrid_Allocator<T_STRUCT,d> SPG_Allocator;
    typedef typename SPG_Allocator::Array<unsigned>::type SPG_Flags_Array_Type;
    typedef SPGrid_Set<SPG_Flags_Array_Type> SPG_Set_Type;

    typedef VECTOR<int,d> T_INDEX;
    typedef VECTOR<T,d> TV;

    SPG_Set_Type& set;
    SPG_Flags_Array_Type& flags;

    const RANGE<T_INDEX> active_domain,interior_domain;
    T_INDEX block_size;

    GRID<TV> grid;

public:
    Rasterizer(SPG_Set_Type& set_input,const RANGE<T_INDEX>& active_domain_input,const int boundary_width=1)
        :set(set_input),flags(set.array),active_domain(active_domain_input),interior_domain(active_domain.Thickened(-boundary_width))
    {
        block_size=flags.geometry.Block_Size().template Cast<T_INDEX>();
        PHYSBAM_ASSERT(active_domain.min_corner==T_INDEX());
        
        grid=GRID<TV>(active_domain.max_corner-1,RANGE<TV>::Unit_Box());
    }

    bool Consume(const RANGE<VECTOR<int,d> >& range){

        // Make sure we did not descend to a sub-block size
        PHYSBAM_ASSERT(range.Edge_Lengths().All_Greater_Equal(block_size));

        // Recurse until we are down to the size of a single block
        if(range.Edge_Lengths()!=block_size) return true;

        // If the block is fully exterior, do nothing
        RANGE<T_INDEX> adjusted_range(range);adjusted_range.max_corner-=1; // Adjusted range is inclusive of max corner
        if(RANGE<T_INDEX>::Intersect(adjusted_range,active_domain).Empty()) return false;

        RANGE<T_INDEX> expanded_range=adjusted_range.Thickened(1);
        ARRAY<unsigned,T_INDEX> flags_local(expanded_range);

        for(RANGE_ITERATOR<d> iterator(expanded_range);iterator.Valid();iterator.Next()){
            T_INDEX cell_index=iterator.Index();

            if(interior_domain.Lazy_Outside(cell_index)) { flags_local(cell_index) = SPGrid_Solver_Cell_Type_Dirichlet; continue; } // This is an exterior cell -- do nothing

            TV center_X=grid.Center(cell_index);
            if((center_X-.5f).Magnitude_Squared() <= .09f) continue; // Neumann bubbles are exterior no matter what
            
            flags_local(cell_index) = SPGrid_Solver_Cell_Type_Active;
        }

        bool active_or_dirichlet_found=false;

        for(RANGE_ITERATOR<d> iterator(adjusted_range);iterator.Valid();iterator.Next()){
            T_INDEX cell_index=iterator.Index();

            if(flags_local(cell_index) == SPGrid_Solver_Cell_Type_Active) active_or_dirichlet_found=true;
            if(flags_local(cell_index) == SPGrid_Solver_Cell_Type_Dirichlet){
                bool active_neighbor_found=false;
                for(int v=1;v<=d;v++)
                    if(flags_local(cell_index+T_INDEX::Axis_Vector(v)) == SPGrid_Solver_Cell_Type_Active ||
                       flags_local(cell_index-T_INDEX::Axis_Vector(v)) == SPGrid_Solver_Cell_Type_Active)
                        active_neighbor_found=true;
                if(active_neighbor_found) active_or_dirichlet_found=true;
                else flags_local(cell_index) = 0x0u;
            }
        }

        if(active_or_dirichlet_found){
            set.MarkPageActive(SPG_Flags_Array_Type::MASK::Linear_Offset(std_array<int,d>(range.min_corner)));
            unsigned* flags_ptr=&flags(std_array<int,d>(adjusted_range.min_corner));
            for(RANGE_ITERATOR<d> iterator(adjusted_range);iterator.Valid();iterator.Next(),flags_ptr++)
                *flags_ptr=flags_local(iterator.Index());
        }

        return false;
    }
};

template<class T,class T_STRUCT,int d>
void Visualize_Heightfield(STREAM_TYPE stream_type,SPGrid_Allocator<T_STRUCT,d>& allocator,SPGrid_Set<typename SPGrid_Allocator<T_STRUCT,d>::Array<unsigned>::type>& set,
                           GRID<VECTOR<T,d> >& grid,T T_STRUCT::* data_channel,unsigned T_STRUCT::* flag_channel,const std::string output_directory,const int frame,const T scale=1){
    typedef VECTOR<int,d> T_INDEX;
    typedef VECTOR<T,d> TV;
    typedef SPGrid_Allocator<T_STRUCT,d> SPG_Allocator;
    typedef typename SPG_Allocator::Array<const unsigned>::type Const_Flags_Array_Type;
    typedef typename SPG_Allocator::Array<const T>::type Const_Data_Array_Type;
    typedef typename SPG_Allocator::Array<T>::mask T_MASK;
    
    typedef VECTOR<T,3> TV_3;
    PHYSBAM_ASSERT(sizeof(T) == sizeof(unsigned));
    Initialize_Geometry_Particle();Initialize_Read_Write_Structures();

    GEOMETRY_PARTICLES<TV_3> particles;

    TRIANGULATED_SURFACE<T>* value_surface = TRIANGULATED_SURFACE<T>::Create(particles);
    
    HASHTABLE<T_INDEX, int> cell_hash; 
    Const_Flags_Array_Type flags=allocator.Get_Const_Array(flag_channel);
    Const_Data_Array_Type data=allocator.Get_Const_Array(data_channel);
    int counter = 1;
    for(SPGrid_Block_Iterator<T_MASK> iterator(set.Get_Blocks());iterator.Valid();iterator.Next_Block()){
        unsigned long offset=iterator.Offset();
        T_INDEX base_index=iterator.Index().template Cast<T_INDEX>();
        for(RANGE_ITERATOR<d> iterator(RANGE<T_INDEX>(base_index,base_index+allocator.Block_Size().template Cast<T_INDEX>()-1));
            iterator.Valid();
            iterator.Next(),offset+=sizeof(unsigned int)){
            if(flags(offset) & (SPGrid_Solver_Cell_Type_Active | SPGrid_Solver_Cell_Type_Interface)){
                cell_hash.Insert(iterator.Index(),counter++);
            }
        }
    }
    particles.array_collection->Add_Elements(cell_hash.Size());
    // Create segment lattice for interior cells
    for(HASHTABLE_ITERATOR<T_INDEX, int> iterator(cell_hash);iterator.Valid();iterator.Next()){
        unsigned long offset = T_MASK::Linear_Offset(std_array<int,d>(iterator.Key()));
        const TV X=grid.Center(iterator.Key());
        particles.X(iterator.Data())=X.Insert((T)scale*data(offset),3);
    }
    for(HASHTABLE_ITERATOR<T_INDEX, int> iterator(cell_hash);iterator.Valid();iterator.Next()){
        unsigned long offset = T_MASK::Linear_Offset(std_array<int,d>(iterator.Key()));
        if(flags(offset) & (SPGrid_Solver_Cell_Type_Active | SPGrid_Solver_Cell_Type_Interface)){
            int x_plus=cell_hash.Get_Default(iterator.Key() + T_INDEX::Axis_Vector(1),-1);
            int y_plus=cell_hash.Get_Default(iterator.Key() + T_INDEX::Axis_Vector(2),-1);
            int xy_plus=cell_hash.Get_Default(iterator.Key() + T_INDEX::Axis_Vector(1) + T_INDEX::Axis_Vector(2),-1);
            if(x_plus != -1 && xy_plus != -1)
                value_surface->mesh.elements.Append(VECTOR<int,3>(iterator.Data(),x_plus,xy_plus));
            if(y_plus != -1 && xy_plus != -1)
                value_surface->mesh.elements.Append(VECTOR<int,3>(xy_plus,y_plus,iterator.Data()));
            if(xy_plus == -1 && (x_plus != -1 && y_plus != -1))
                value_surface->mesh.elements.Append(VECTOR<int,3>(iterator.Data(),x_plus,y_plus));
        }
    }
    value_surface->Update_Number_Nodes();
    DEFORMABLE_GEOMETRY_COLLECTION<TV_3> collection(particles);
    collection.Add_Structure(value_surface);
    FILE_UTILITIES::Create_Directory(output_directory+"/"+STRING_UTILITIES::Value_To_String(frame));
    collection.Write(stream_type,output_directory,frame,frame,true);
    FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+STRING_UTILITIES::Value_To_String(frame)+"/heightfield.tri",*value_surface);
    //delete value_surface;
}

int main(int argc,char *argv[])
{
    enum{d=3};
    typedef SPGrid_Allocator<T_STRUCT,d> SPG_Allocator;
    typedef SPG_Allocator::Array<T>::type SPG_Data_Array_Type;
    typedef SPG_Allocator::Array<unsigned>::type SPG_Flags_Array_Type;
    typedef SPG_Allocator::Array<const unsigned>::type SPG_Const_Flags_Array_Type;
    typedef SPG_Allocator::Array<T>::mask T_MASK;
    typedef SPGrid_Set<SPG_Flags_Array_Type> SPG_Set_Type;

    enum{
        block_xsize=1<<SPG_Flags_Array_Type::MASK::block_xbits,
        block_ysize=1<<SPG_Flags_Array_Type::MASK::block_ybits,
        block_zsize=1<<SPG_Flags_Array_Type::MASK::block_zbits,
    };

    typedef VECTOR<int,d> T_INDEX;
    LOG::Initialize_Logging();

    PARSE_ARGS parse_args;
    if(d==2){
        parse_args.Add_Vector_2D_Argument("-size",VECTOR<double,2>(256,256),"n n","Size of 2D grid");
        parse_args.Add_Vector_2D_Argument("-position_delta",VECTOR<double,2>(64,128),"n n","Position of the delta function");
    }else if(d==3){
        parse_args.Add_Vector_3D_Argument("-size",VECTOR<double,3>(256,256,256),"n n n","Size of 3D grid");
        parse_args.Add_Vector_3D_Argument("-position_delta",VECTOR<double,3>(64,128,128),"n n n","Position of the delta function");}
    parse_args.Add_String_Argument("-d","output","Directory");
    parse_args.Parse(argc,argv);

    std::string output_directory=parse_args.Get_String_Value("-d");
    
    std_array<int,d> size;
    if(d==2){
        for(int v=1;v<=d;v++) size(v-1)=static_cast<unsigned int>(parse_args.Get_Vector_2D_Value("-size")(v));    
        LOG::cout<<"Size of 2D grid         : "<<size<<std::endl;
    }else if(d==3){
        for(int v=1;v<=d;v++) size(v-1)=static_cast<unsigned int>(parse_args.Get_Vector_3D_Value("-size")(v));    
        LOG::cout<<"Size of 3D grid         : "<<size<<std::endl;}
    
    T dx=1;

    T_INDEX delta_position; 
    if(d==2)
        for(int v=1;v<=d;v++) delta_position(v)=static_cast<unsigned int>(parse_args.Get_Vector_2D_Value("-position_delta")(v));
    else if(d==3)
        for(int v=1;v<=d;v++) delta_position(v)=static_cast<unsigned int>(parse_args.Get_Vector_3D_Value("-position_delta")(v));
    LOG::cout<<"Position of the delta function: "<<delta_position<<std::endl;

    typedef VECTOR<T,d> TV;
    GRID<TV> grid(size.Cast<T_INDEX>(),RANGE<TV>(TV(),TV::All_Ones_Vector()));

    SPG_Allocator allocator(size);
    T_INDEX padded_size(allocator.xsize_padded,allocator.ysize_padded,allocator.zsize_padded);
    LOG::cout<<"Padded (allocated) size : "<<padded_size<<std::endl;

    SPG_Flags_Array_Type flags=allocator.Get_Array(flags_field);
    SPG_Data_Array_Type u_array=allocator.Get_Array(u_field);
    SPG_Data_Array_Type b_array=allocator.Get_Array(b_field);
    SPG_Data_Array_Type r_array=allocator.Get_Array(r_field);
    SPG_Set_Type set(flags);
    RANGE<T_INDEX> active_range(T_INDEX(),size.Cast<T_INDEX>()-1);
    RANGE<T_INDEX> hierarchical_range(T_INDEX(),padded_size);
    Rasterizer<T_STRUCT,T,d> rasterizer(set,active_range,1);
    {LOG::SCOPE scope("Rasterization");
        for(HIERARCHICAL_RANGE_ITERATOR<d,Rasterizer<T_STRUCT,T,d> > iterator(hierarchical_range,rasterizer);iterator.Valid();iterator.Next());}
    set.Refresh_Block_Offsets();

    SPG_Data_Array_Type b=allocator.Get_Array(b_field);
    LOG::cout<<"Delta Position: "<<delta_position<<std::endl;
    PHYSBAM_ASSERT(flags(T_MASK::Linear_Offset(std_array<int,d>(delta_position)))!=0);
    b(T_MASK::Linear_Offset(std_array<int,d>(delta_position))) = 1/(dx*dx) * 4;
   
    typedef float RW;
    STREAM_TYPE stream_type((RW()));
    
    //T output_scale=100;

    //FILE_UTILITIES::Create_Directory("output_height_field");
    //FILE_UTILITIES::Create_Directory("output_height_field/common");
    //FILE_UTILITIES::Write_To_File(stream_type,"output_height_field/common/grid",grid);
    //Visualize_Heightfield(stream_type,allocator,set,grid,b_field,flags_field,"output_height_field",0,1.f/dx*output_scale);

    SPGrid_Computations::Flaging<T_STRUCT,d>(allocator,set.Get_Blocks(),flags_field,set);
    LOG::cout<<"Data size: "<<(double)(set.Get_Blocks().second)*sizeof(T)*4*set.Get_Blocks().second/1024./1024./1024.<<" GB"<<std::endl;
    FILE_UTILITIES::Create_Directory(output_directory);
    SPGRID_READ_WRITE<T_STRUCT,T,d>::Write_SPGrid(output_directory,allocator,set,b_field);

    LOG::Finish_Logging();
    return 0;
}
