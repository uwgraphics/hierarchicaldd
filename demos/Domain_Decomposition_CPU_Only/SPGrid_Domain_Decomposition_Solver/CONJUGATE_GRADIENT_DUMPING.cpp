//#####################################################################
//Removed For Anonymity: Copyright Authors of Submission pap173s1
//Removed For Anonymity: Copyright Authors of Submission pap173s1
//#####################################################################
// Class CONJUGATE_GRADIENT
//#####################################################################
#include "SPGRID_DOMAIN_DECOMPOSITION_DATA.h"
#include "CONJUGATE_GRADIENT_DUMPING.h"
#include "CG_VECTOR.h"
#include "SPGrid_Convertor.h"
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Read_Write/Grids_Uniform/READ_WRITE_GRID.h>
#include <PhysBAM_Tools/Read_Write/Grids_Uniform_Arrays/READ_WRITE_ARRAYS.h>
#include <PhysBAM_Tools/Krylov_Solvers/KRYLOV_VECTOR_BASE.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Geometry/Geometry_Particles/REGISTER_GEOMETRY_READ_WRITE.h>
#include <PhysBAM_Geometry/Solids_Geometry/DEFORMABLE_GEOMETRY_COLLECTION.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/FREE_PARTICLES.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/SEGMENTED_CURVE.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <PhysBAM_Geometry/Read_Write/Geometry/READ_WRITE_TRIANGULATED_SURFACE.h>
#include <PhysBAM_Tools/Data_Structures/HASHTABLE.h>
#include <PhysBAM_Tools/Data_Structures/HASHTABLE_ITERATOR.h>
#include <SPGrid_Fluids/Grids/Rasterizers/HIERARCHICAL_RASTERIZER.h>
#include <SPGrid_Fluids/Grids/GRID_HIERARCHY_INITIALIZER.h>
#include <SPGrid_Fluids/Projection/GRID_HIERARCHY_PROJECTION.h>
#include <cfloat>
#include <limits>
using namespace PhysBAM;
using namespace SPGrid;

template<class T,class T_STRUCT,int d>
class Visualize_Heightfield_Local;

template<class T,class T_STRUCT>
class Visualize_Heightfield_Local<T,T_STRUCT,3>{
public:
    static void Visualize(STREAM_TYPE stream_type,const SPGrid_Allocator<T_STRUCT,3>& allocator,
                          const SPGrid_Allocator<T_STRUCT,3>& flag_allocator,
                          const SPGrid_Set<typename SPGrid_Allocator<T_STRUCT,3>::Array<unsigned>::type>& set,
                          const GRID<VECTOR<T,3> >& grid,T T_STRUCT::* data_channel,unsigned T_STRUCT::* flag_channel,
                          const std::string output_directory,const int frame,const T scale=1){}
};

template<class T,class T_STRUCT>
class Visualize_Heightfield_Local<T,T_STRUCT,2>{
public: 
    static void Visualize(STREAM_TYPE stream_type,const SPGrid_Allocator<T_STRUCT,2>& allocator,
                          const SPGrid_Allocator<T_STRUCT,2>& flag_allocator,
                          const SPGrid_Set<typename SPGrid_Allocator<T_STRUCT,2>::Array<unsigned>::type>& set,
                          const GRID<VECTOR<T,2> >& grid,T T_STRUCT::* data_channel,unsigned T_STRUCT::* flag_channel,
                          const std::string output_directory,const int frame,const T scale=1){
        enum{d=2};
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
        Const_Flags_Array_Type flags=flag_allocator.Get_Const_Array(flag_channel);
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
};

template<class T,int d>
class Initialize_Grid;

template<class T>
class Initialize_Grid<T,2>{
public:
    static void Initialize(GRID<VECTOR<T,2> >& grid,T x_size,T y_size)
    {
        typedef VECTOR<T,2> TV;
        grid.Initialize(x_size,y_size,RANGE<TV>(TV(),TV::All_Ones_Vector()),true);
    }
};

template<class T>
class Initialize_Grid<T,3>{
public:
    static void Initialize(GRID<VECTOR<T,3> >& grid,T x_size,T y_size)
    {
    }
};

//#####################################################################
// Destructor
//#####################################################################
template<class T,int d> CONJUGATE_GRADIENT_DUMPING<T,d>::
~CONJUGATE_GRADIENT_DUMPING()
{}
//#####################################################################
// Function Solve
//#####################################################################
template<class T,int d> bool CONJUGATE_GRADIENT_DUMPING<T,d>::
Solve(const KRYLOV_SYSTEM_BASE<T>& system,KRYLOV_VECTOR_BASE<T>& x,const KRYLOV_VECTOR_BASE<T>& b,KRYLOV_VECTOR_BASE<T>& q,KRYLOV_VECTOR_BASE<T>& s,KRYLOV_VECTOR_BASE<T>& r,
      KRYLOV_VECTOR_BASE<T>& k,KRYLOV_VECTOR_BASE<T>& z,const T tolerance,const int min_iterations,const int max_iterations)
{
    Initialize_Grid<T,d>::Initialize(grid,CG_VECTOR<SPGRID_DOMAIN_DECOMPOSITION_DATA<T>,T,d>::Allocator(r).xsize,
                                     CG_VECTOR<SPGRID_DOMAIN_DECOMPOSITION_DATA<T>,T,d>::Allocator(r).ysize);

    // NOTE: you should never try to make copies of VECTOR_T's inside here as they could be indirect.
    static const T small_number=std::numeric_limits<T>::epsilon();
    system.Set_Boundary_Conditions(x);
    T rho_old=(T)FLT_MAX;T convergence_norm=0;
    int iterations;for(iterations=0;;iterations++){
        bool restart=!iterations || (restart_iterations && iterations%restart_iterations==0);
        if(restart){
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
            if(print_residuals) LOG::cout<<"restarting cg"<<std::endl;
#endif
            r=b;system.Multiply(x,q);r-=q;system.Project(r);}
        // stopping conditions
        system.Project_Nullspace(r);
        convergence_norm=system.Convergence_Norm(r);
        
        T scale=1;
        
        Visualize_Heightfield_Local<T,SPGRID_DOMAIN_DECOMPOSITION_DATA<T>,d>::
            Visualize(stream_type,CG_VECTOR<SPGRID_DOMAIN_DECOMPOSITION_DATA<T>,T,d>::Allocator(r),
                      CG_VECTOR<SPGRID_DOMAIN_DECOMPOSITION_DATA<T>,T,d>::Allocator(r),
                      CG_VECTOR<SPGRID_DOMAIN_DECOMPOSITION_DATA<T>,T,d>::Set(r),
                      grid,CG_VECTOR<SPGRID_DOMAIN_DECOMPOSITION_DATA<T>,T,d>::Cg_Vector(r).field,
                      &SPGRID_DOMAIN_DECOMPOSITION_DATA<T>::flags,
                      output_directory+"_residual",iterations,scale);

        Visualize_Heightfield_Local<T,SPGRID_DOMAIN_DECOMPOSITION_DATA<T>,d>::
            Visualize(stream_type,CG_VECTOR<SPGRID_DOMAIN_DECOMPOSITION_DATA<T>,T,d>::Allocator(x),
                      CG_VECTOR<SPGRID_DOMAIN_DECOMPOSITION_DATA<T>,T,d>::Allocator(r),
                      CG_VECTOR<SPGRID_DOMAIN_DECOMPOSITION_DATA<T>,T,d>::Set(x),
                      grid,CG_VECTOR<SPGRID_DOMAIN_DECOMPOSITION_DATA<T>,T,d>::Cg_Vector(x).field,
                      &SPGRID_DOMAIN_DECOMPOSITION_DATA<T>::flags,
                      output_directory+"_result",iterations,scale);
        
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
        if(print_residuals) {LOG::cout<<convergence_norm<<std::endl;}
#endif
        if(convergence_norm<=tolerance && (iterations>=min_iterations || convergence_norm<small_number)){
            if(print_diagnostics) LOG::Stat("cg iterations",iterations);if(iterations_used) *iterations_used=iterations;return true;}
        if(iterations==max_iterations) break;
        // actual iteration
        const KRYLOV_VECTOR_BASE<T>& mr=system.Precondition(r,z);
        T rho=(T)system.Inner_Product(mr,r);
        if(restart) s=mr;
        else s.Copy(rho/rho_old,s,mr);
        system.Multiply(s,q);
        system.Project(q);
        T s_dot_q=(T)system.Inner_Product(s,q);
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
        if(s_dot_q<=0) {LOG::cout<<"CG: matrix appears indefinite or singular, s_dot_q/s_dot_s="<<s_dot_q/(T)system.Inner_Product(s,s)<<std::endl;}
#endif
        T alpha=s_dot_q?rho/s_dot_q:(T)FLT_MAX;
        x.Copy(alpha,s,x);
        r.Copy(-alpha,q,r);
        rho_old=rho;}

#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
    if(print_diagnostics) LOG::Stat("cg iterations",iterations);if(iterations_used) *iterations_used=iterations;
    if(print_diagnostics) {LOG::cout<<"cg not converged after "<<max_iterations<<" iterations, error = "<<convergence_norm<<std::endl;}
#endif
    return false;
}
//#####################################################################
template class CONJUGATE_GRADIENT_DUMPING<float,2>;
template class CONJUGATE_GRADIENT_DUMPING<float,3>;
