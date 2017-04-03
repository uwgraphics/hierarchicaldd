//#####################################################################
// Copyright 2002-2007, Doug Enright, Ron Fedkiw, Eran Guendelman, Geoffrey Irving, Sergey Koltakov, Neil Molino, Andrew Selle, Eftychios Sifakis, Rachel Weinstein.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class CHIMERA_PARTICLE_LEVELSET_IMPLICIT_OBJECT
//#####################################################################
#ifndef __CHIMERA_PARTICLE_LEVELSET_IMPLICIT_OBJECT__
#define __CHIMERA_PARTICLE_LEVELSET_IMPLICIT_OBJECT__

#include <PhysBAM_Tools/Matrices/FRAME.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_1D.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_2D.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_3D.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_POLICY_UNIFORM.h>
#include <PhysBAM_Geometry/Implicit_Objects/IMPLICIT_OBJECT.h>
#include <PhysBAM_Geometry/Implicit_Objects_Uniform/CHIMERA_LEVELSET_IMPLICIT_OBJECT.h>
#include <PhysBAM_Geometry/Basic_Geometry_Intersections/RAY_BOX_INTERSECTION.h>
#include <PhysBAM_Geometry/Level_Sets/IMPLICIT_OBJECT_ON_A_RAY.h>
#include <PhysBAM_Geometry/Level_Sets/LEVELSET_UTILITIES.h>
#include <PhysBAM_Dynamics/Particles/PARTICLE_LEVELSET_REMOVED_PARTICLES.h>
#include <PhysBAM_Dynamics/Level_Sets/REMOVED_PARTICLES_BLENDER_3D.h>
#include <PhysBAM_Tools/Nonlinear_Equations/ITERATIVE_SOLVER.h>
#include <PhysBAM_Geometry/Basic_Geometry_Intersections/RAY_PLANE_INTERSECTION.h>
#include <PhysBAM_Geometry/Basic_Geometry/PLANE.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_NODE.h>
#include <PhysBAM_Tools/Math_Tools/clamp.h>
#include <PhysBAM_Tools/Math_Tools/cbrt.h>
#include <PhysBAM_Tools/Polynomials/CUBIC.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TOPOLOGY_BASED_SIMPLEX_POLICY.h>
#include <PhysBAM_Geometry/Topology/TOPOLOGY_POLICY.h>
#include <PhysBAM_Geometry/Geometry_Particles/GEOMETRY_PARTICLES.h>
#include <PhysBAM_Geometry/Spatial_Acceleration/BOX_HIERARCHY.h>
namespace PhysBAM{

//template<class TV> class CHIMERA_PARTICLE_LEVELSET_IMPLICIT_OBJECT;
//
//template<class TV>
//class IMPLICIT_OBJECT_ON_A_RAY<CHIMERA_PARTICLE_LEVELSET_IMPLICIT_OBJECT<TV> >:public NONLINEAR_FUNCTION<typename TV::SCALAR(typename TV::SCALAR)>
//{
//    typedef typename TV::SCALAR T;
//public:
//    const CHIMERA_PARTICLE_LEVELSET_IMPLICIT_OBJECT<TV>& implicit_object;
//    RAY<TV>& ray;
//    //ARRAY<TRIPLE<TV,TV,T> > candidate_intersect_particles;
//        
//    IMPLICIT_OBJECT_ON_A_RAY(const CHIMERA_PARTICLE_LEVELSET_IMPLICIT_OBJECT<TV>& implicit_object_input,RAY<TV>& ray_input)
//        :implicit_object(implicit_object_input),ray(ray_input)
//    {
//        //for(int i=1;i<=implicit_object.grids.Size();i++){
//        //    int grid_index=implicit_object.grid_order(i);
//        //    for(UNIFORM_GRID_ITERATOR_NODE<TV> it(*implicit_object.grids(grid_index));it.Valid();it.Next())
//        //        if(PARTICLE_LEVELSET_REMOVED_PARTICLES<VECTOR<T,3> >* particles=(*implicit_object.negative_particle_data(grid_index))(it.Node_Index())){
//        //            //TODO
//        //        }
//        //}
//    }
//
//    T operator()(const T x) const PHYSBAM_OVERRIDE
//    {return implicit_object(ray.Point(x));}
//
////#####################################################################
//};   

template<class TV> class GRID;

template<class TV>
class CHIMERA_PARTICLE_LEVELSET_IMPLICIT_OBJECT:public CHIMERA_LEVELSET_IMPLICIT_OBJECT<TV>
{
    typedef typename TV::SCALAR T;
    typedef typename TV::template REBIND<int>::TYPE TV_INT;
    //typedef typename GRID_ARRAYS_POLICY<GRID<TV> >::ARRAYS_SCALAR T_ARRAYS_SCALAR;typedef typename T_ARRAYS_SCALAR::template REBIND<TV>::TYPE T_ARRAYS_VECTOR;
    typedef typename LEVELSET_POLICY<GRID<TV> >::LEVELSET T_LEVELSET;
    //typedef typename MATRIX_POLICY<TV>::SYMMETRIC_MATRIX T_SYMMETRIC_MATRIX;
    //enum WORKAROUND {d=TV::m};
    //typedef VECTOR<T,d-1> T_PRINCIPAL_CURVATURES;
public:
    typedef CHIMERA_LEVELSET_IMPLICIT_OBJECT<TV> BASE;
    using BASE::box;using BASE::minimum_cell_size;

    using BASE::grids;
    using BASE::grid_frames;
    using BASE::levelset_data;
    using BASE::levelsets;

    //using BASE::number_of_ghost_cells;
    using BASE::blend_size;
    using BASE::use_secondary_interpolation;
    using BASE::extend_levelset_to_infinity;
    using BASE::water_level;
    
    ARRAY<ARRAY<PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>*,TV_INT>*> positive_particle_data;
    ARRAY<ARRAY<PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>*,TV_INT>*> negative_particle_data;

    std::string particle_processing_mode;
    T particle_scale;
    T blending_parameter;
    CUBIC<T> kernel;

    T average_particle_radius;

    const bool use_ellipsoid;
    const T dt;

    int particle_number_of_ghost_cells;

    const T integration_step_scale;
    const T levelset_kernel_scale;
    const T particle_maximum_radius;
    const bool particle_use_velocity_scaling;
    const T particle_elongation_cap;

    T root;

    ARRAY<ARRAY<ARRAY<PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>*,TV_INT>*>*> particle_data_array;
    ARRAY<ARRAY<bool,TV_INT> > nodes_covered_by_finer_grids;
    GRID<TV> background_grid;
    ARRAY<int,TV_INT> particle_count_touching_cell;

    BOX_HIERARCHY<TV> particle_hierarchy;
    BOX_HIERARCHY<TV> thickened_particle_hierarchy;
    ARRAY<TV> particle_positions;
    ARRAY<TV> particle_velocities;
    ARRAY<T> particle_radii;
    ARRAY<T> particle_radii_long;
    ARRAY<TV> particle_axis_directions;

    using BASE::interpolation_domains;
    using BASE::grid_order;
    using BASE::interpolation;
public:

    CHIMERA_PARTICLE_LEVELSET_IMPLICIT_OBJECT(const ARRAY<GRID<TV>*>& grids_input,const ARRAY<FRAME<TV>*>& grid_frames_input,const ARRAY<ARRAY<T,TV_INT>*>& levelset_data_input,T blendsz,int levelset_ghost,
                                              const ARRAY<ARRAY<PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>*,TV_INT>*>& positive_particle_data_input,
                                              const ARRAY<ARRAY<PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>*,TV_INT>*>& negative_particle_data_input,std::string pmode,T pscale,T pblend,T inte_step_scale,T lset_kern_scale,T p_max_r,bool p_use_vel_scaling,T p_elong_cap,T p_elong_fps)
    :CHIMERA_LEVELSET_IMPLICIT_OBJECT<TV>(grids_input,grid_frames_input,levelset_data_input,blendsz,levelset_ghost),
     positive_particle_data(positive_particle_data_input),negative_particle_data(negative_particle_data_input),particle_processing_mode(pmode),particle_scale(pscale),
     blending_parameter(pblend),kernel(2,-3,0,1),use_ellipsoid(true),dt((T)1/p_elong_fps),particle_number_of_ghost_cells(levelset_ghost),integration_step_scale(inte_step_scale),levelset_kernel_scale(lset_kern_scale),particle_maximum_radius(p_max_r),particle_use_velocity_scaling(p_use_vel_scaling),particle_elongation_cap(p_elong_cap)
    {
        kernel.c0-=blending_parameter;
        kernel.Compute_Roots_Noniterative_In_Interval(0,1);
        root=kernel.root1;
        LOG::cout<<"root1 = "<<kernel.root1<<std::endl;
        if(kernel.roots<1){LOG::cerr<<"Error: kernel.roots=="<<kernel.roots<<std::endl;PHYSBAM_FATAL_ERROR();}
        kernel.c0=1;

        int n_cells=0;
        for(int grid_index=1;grid_index<=grids.Size();grid_index++)
            n_cells+=grids(grid_index)->Counts().Product();
        background_grid=GRID<TV>::Create_Even_Sized_Grid_Given_Cell_Size(box,cbrt(box.Size()/n_cells),false);
        particle_count_touching_cell.Resize(background_grid.Domain_Indices(),true,false,0);

        T sum_particle_radii(0);
        int sum_particle_number=0;
        LOG::cout<<"particle-number-of-ghost-cells="<<particle_number_of_ghost_cells<<std::endl;
        for(int grid_index=1;grid_index<=grids.Size();grid_index++){
            T maximum_particle_radius(-FLT_MAX);
            int n_particle=0;
            int histogram[10]={0,0,0,0,0,0,0,0,0,0};
            int current_grid_particle_number_of_ghost_cells=negative_particle_data(grid_index)->Number_Of_Ghost_Cells();
            if(current_grid_particle_number_of_ghost_cells!=particle_number_of_ghost_cells) LOG::cout<<"Warning: particle-number-of-ghost-cells from array "<<grid_index<<" ="<<current_grid_particle_number_of_ghost_cells<<std::endl;
            for(UNIFORM_GRID_ITERATOR_NODE<TV> it(*grids(grid_index),current_grid_particle_number_of_ghost_cells);it.Valid();it.Next())
                if(PARTICLE_LEVELSET_REMOVED_PARTICLES<VECTOR<T,3> >* particles=(*negative_particle_data(grid_index))(it.Node_Index()))
                    //TODO TODO TODO move this after building hierarchy
                    for(int p=1;p<=particles->array_collection->Size();p++){
                        TV world_X=grid_frames(grid_index)->operator*(particles->X(p));
                        T radius_from_array=particles->radius(p);
                        T radius_long=particle_scale*radius_from_array;TV axis_direction(1,0,0);
                        Get_Elongated_Influence_Radius(particles->V(p),radius_long,axis_direction);
                        T influence_radius=particle_scale*radius_from_array/root;
                        T influence_radius_long=radius_long/root;
                        RANGE<TV> bbox=Ellipsoid_Bounding_Box(world_X,influence_radius,influence_radius_long,axis_direction);
                        for(UNIFORM_GRID_ITERATOR_CELL<TV> itb(background_grid,background_grid.Clamp_To_Cell(bbox,0));itb.Valid();itb.Next())
                            particle_count_touching_cell(itb.Cell_Index())+=1;
                        if(radius_from_array>maximum_particle_radius) maximum_particle_radius=radius_from_array;
                        sum_particle_radii+=particles->radius(p);
                        ++sum_particle_number;
                        ++n_particle;
                        int bin=min(9,int(particles->radius(p)*particle_scale/(root*grids(grid_index)->dX(1))*10));  // pradius / root < .5 dx        pradius / root = (bin / 10) dx
                        histogram[bin]+=1;
                    }
            for(UNIFORM_GRID_ITERATOR_FACE<TV> it(*grids(grid_index));it.Valid();it.Next())
                if(levelset_data(grid_index)->operator()(it.First_Cell_Index())*levelset_data(grid_index)->operator()(it.Second_Cell_Index())<=0){
                    TV world_X=grid_frames(grid_index)->operator*(it.Location());
                    RANGE<TV> bbox(world_X-(T).5*grids(grid_index)->dX(1),world_X+(T).5*grids(grid_index)->dX(1));
                    for(UNIFORM_GRID_ITERATOR_CELL<TV> itb(background_grid,background_grid.Clamp_To_Cell(bbox,0));itb.Valid();itb.Next())
                        particle_count_touching_cell(itb.Cell_Index())+=1;
                }

            LOG::cout<<"grid "<<grid_index<<" particle statistics: n_particles="<<n_particle<<" max_r="<<maximum_particle_radius<<" dx="<<grids(grid_index)->dX(1)<<std::endl;
            LOG::cout<<"r >";
            for(int bin=0;bin<10;++bin) LOG::cout<<"\t"<<(T)bin/(T)10*grids(grid_index)->dX(1)*root;
            LOG::cout<<std::endl;
            LOG::cout<<"r > (of dx)";
            for(int bin=0;bin<10;++bin) LOG::cout<<"\t"<<(T)bin/(T)10*root;
            LOG::cout<<std::endl;
            LOG::cout<<"count";
            for(int bin=0;bin<10;++bin) LOG::cout<<"\t"<<histogram[bin];
            LOG::cout<<std::endl;
            LOG::cout<<"percent";
            for(int bin=0;bin<10;++bin) LOG::cout<<"\t"<<(T)histogram[bin]/n_particle*T(100)<<"%";
            LOG::cout<<std::endl;
        }
        if(sum_particle_number)
            average_particle_radius=(sum_particle_radii/sum_particle_number/root);
        else
            average_particle_radius=grids(grid_order(1))->dX(1);
        LOG::cout<<"AVERAGE PARTICLE RADIUS="<<average_particle_radius*root<<std::endl;
        if(particle_maximum_radius>1e-7) LOG::cout<<"Particle radii capped at "<<particle_maximum_radius<<std::endl;
        if(particle_use_velocity_scaling) LOG::cout<<"Using velocity-proportional elongation"<<std::endl;
        else LOG::cout<<"Using fixed elongation=3"<<std::endl;

        particle_data_array.Append(&positive_particle_data);
        particle_data_array.Append(&negative_particle_data);

        Build_Bounding_Box_Hierarchy_For_Particles(grids_input,grid_frames_input,negative_particle_data_input);
    }

    void Get_Elongated_Influence_Radius(const TV particle_velocity,T& long_radius,TV& axis_direction) const
    {
        if(use_ellipsoid){
            T velocity_magnitude_squared=particle_velocity.Magnitude_Squared();
            if(velocity_magnitude_squared>1e-8){
                T particle_speed=sqrt(velocity_magnitude_squared);
                if(particle_use_velocity_scaling){T max_elong=particle_elongation_cap*long_radius;long_radius+=(T).5*dt*particle_speed;long_radius=min(long_radius,max_elong);}
                else long_radius*=3;
                axis_direction=particle_velocity/particle_speed;
            }
        }
    }
    RANGE<TV> Ellipsoid_Bounding_Box(const TV& center,const T radius,const T radius_long, const TV& axis_direction)
    {
        if(radius_long/radius<(T)1.0001) return RANGE<TV>(center).Thickened(radius);
        TV y=axis_direction.Unit_Orthogonal_Vector();
        TV z=TV::Cross_Product(axis_direction,y);
        TV long_edge=(T)2*radius_long*axis_direction;
        TV short_edge_y=(T)2*radius*y;
        TV short_edge_z=(T)2*radius*z;
        TV corner=center-(T).5*(long_edge+short_edge_y+short_edge_z);
        ORIENTED_BOX<TV> oriented_box(corner,long_edge,short_edge_y,short_edge_z);
        return RANGE<TV>::Intersect(RANGE<TV>(center).Thickened(radius_long),oriented_box.Axis_Aligned_Bounding_Box());
    }


    void Build_Bounding_Box_Hierarchy_For_Particles(const ARRAY<GRID<TV>*>& grids,const ARRAY<FRAME<TV>*>& grid_frames,const ARRAY<ARRAY<PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>*,TV_INT>*>& particle_data)
    {
        ARRAY<RANGE<TV> > particle_ranges;
        ARRAY<RANGE<TV> > thickened_particle_ranges;

        //LOG::cout<<"Elongation=";
        for(int grid_index=1;grid_index<=particle_data.Size();grid_index++){
            for(UNIFORM_GRID_ITERATOR_NODE<TV> iterator(*grids(grid_index),particle_data(grid_index)->Number_Of_Ghost_Cells());iterator.Valid();iterator.Next())
                if(PARTICLE_LEVELSET_REMOVED_PARTICLES<VECTOR<T,3> >* particles=(*particle_data(grid_index))(iterator.Node_Index())){
                    for(int p=1;p<=particles->array_collection->Size();p++){
                        TV location=(*grid_frames(grid_index)) * particles->X(p);
                        TV velocity=grid_frames(grid_index)->Rotation().Rotate(particles->V(p));
                        T radius_from_array=particles->radius(p);
                        if(particle_maximum_radius>1e-7){radius_from_array=min(radius_from_array,particle_maximum_radius);} // this caps particle radii by a global constant
                        T radius_long=particle_scale*radius_from_array;TV axis_direction(1,0,0);
                        Get_Elongated_Influence_Radius(velocity,radius_long,axis_direction);
                        T influence_radius=particle_scale*radius_from_array/root;
                        T influence_radius_long=radius_long/root;
                        //LOG::cout<<" "<<radius_long/(particle_scale*radius_from_array);
                        RANGE<TV> range=Ellipsoid_Bounding_Box(location,influence_radius,influence_radius_long,axis_direction);
                        //LOG::cout<<"a ratio="<<influence_radius_long/influence_radius<<" axis="<<axis_direction<<" range#="<<(range-location)/influence_radius<<std::endl;
                        particle_positions.Append(location);
                        particle_velocities.Append(velocity);
                        particle_radii.Append(influence_radius);
                        particle_radii_long.Append(influence_radius_long);
                        particle_axis_directions.Append(axis_direction);
                        particle_ranges.Append(range);
                        int finest_enclosing_index=Find_Finest_Grid(location);
                        T local_cell_size=grids(finest_enclosing_index)->dX(1);
                        range=range.Thickened(local_cell_size);
                        thickened_particle_ranges.Append(range);
                    }
                }
        }
        //LOG::cout<<std::endl;
        LOG::cout<<"Number of particles is "<<particle_ranges.Size()<<std::endl;
        particle_hierarchy.Clean_Memory();thickened_particle_hierarchy.Clean_Memory();
        if(particle_ranges.Size()){
            particle_hierarchy.Set_Leaf_Boxes(particle_ranges,true);
            thickened_particle_hierarchy.Set_Leaf_Boxes(thickened_particle_ranges,true);}
    }

    T Integration_Step_World_Position(/*const T phi,*/TV location) const
    {
        //return BASE::Integration_Step(phi);

        int finest_enclosing_index=Find_Finest_Grid(location);
        T local_cell_size=grids(finest_enclosing_index)->dX(1);
        TV location_local=grid_frames(finest_enclosing_index)->Inverse_Times(location);

        T distance=BASE::operator()(location);
        T step_size;
        if(distance > 3*local_cell_size){
            /*RANGE<TV> local_location_range(location_local-(T).5*distance,location_local+(T).5*distance);
            //for every grid finer than source grid of distance, check if range intersects with its domain
            if(!grids(finest_enclosing_index)->domain.Lazy_Intersection(local_location_range))
                step_size=(T).5*distance; // not robust - need to 
            else*/ 
                step_size=local_cell_size;//(T).5*distance;    
        }
        else if(distance > local_cell_size) step_size=(T).25*distance;
        else step_size=(T).1*local_cell_size; 
       
        ARRAY<int> neighbor_particle_indices;
        thickened_particle_hierarchy.Intersection_List(location,neighbor_particle_indices);
            
        for(int i=1;i<=neighbor_particle_indices.Size();i++){
            T particle_radius=particle_radii_long(neighbor_particle_indices(i));
            T distance=(location-particle_positions(neighbor_particle_indices(i))).Magnitude();
            T candidate_step_size;
            if(distance>particle_radius*10) candidate_step_size=(T)9.*particle_radius;
            else if(distance>particle_radius*5) candidate_step_size=(T)4.*particle_radius;
            else if(distance>particle_radius*3) candidate_step_size=(T)2.*particle_radius;
            else if(distance>particle_radius*1.5) candidate_step_size=(T).5*particle_radius;
            else  candidate_step_size=(T).2*particle_radius;
            step_size=std::min(step_size,candidate_step_size);
        }
        
        return step_size;
    }

    int Find_Finest_Grid(TV location) const
    {
        int finest_enclosing_index=grid_order(grids.m);
        for(int i=1;i<=grids.m;i++){
            int grid_index=grid_order(i);
            TV location_local=grid_frames(grid_index)->Inverse_Times(location);
            if(grids(grid_index)->domain.Inside(location_local,-grids(1)->dX(1))){finest_enclosing_index=grid_index;break;}}
        return finest_enclosing_index;
    }
    T Local_DX(TV location) const
    {
        int finest_enclosing_index=grid_order(grids.m);
        for(int i=1;i<=grids.m;i++){
            int grid_index=grid_order(i);
            TV location_local=grid_frames(grid_index)->Inverse_Times(location);
            if(grids(grid_index)->domain.Inside(location_local,-grids(grid_index)->dX(1))){finest_enclosing_index=grid_index;break;}
        }
        return grids(finest_enclosing_index)->dX(1);
    }

    void Advance_To_Next_Background_Cell(const RAY<TV>& ray,TV_INT& background_index,T& t) const
    {
        T min_dist_to_cross_cell=FLT_MAX;
        for(int d=1;d<=3;++d)
            if(abs(ray.direction(d))>(T)1e-6){
                T dist=background_grid.domain.min_corner(d)+(background_index(d)+(ray.direction(d)>0?0:-1))*background_grid.dX(d)-ray.Point(t)(d);
                while(dist/ray.direction(d)<(T)1e-6){
                    //LOG::cout<<"dist="<<dist<<" ray.direction(d)="<<ray.direction(d)<<std::endl;
                    background_index(d)+=(ray.direction(d)>0?1:-1);dist+=(ray.direction(d)>0?1:-1)*background_grid.dX(d);
                }
                min_dist_to_cross_cell=min(min_dist_to_cross_cell,dist/ray.direction(d));}
        //LOG::cout<<"direction="<<ray.direction<<" min_dist_to_cross_cell="<<min_dist_to_cross_cell<<" dx="<<background_grid.dX<<std::endl;
        t+=min_dist_to_cross_cell+background_grid.dX.Min()*(T)1e-6;
    }

    bool Intersection(RAY<TV>& ray,const T thickness) const
    {
        bool save_semi_infinite=ray.semi_infinite;T save_t_max=ray.t_max;int save_aggregate=ray.aggregate_id;

        T t_start,t_end;
        bool exit_intersection=false;int exit_aggregate=0;T exit_t_max=0;
        int intersect_box=INTERSECTION::Intersects(ray,box,thickness);
        int outside_box=box.Outside(ray.endpoint,thickness);

        RAY<TV> plane_ray(ray.endpoint,ray.direction);plane_ray.semi_infinite=save_semi_infinite;plane_ray.t_max=save_t_max;plane_ray.aggregate_id=save_aggregate;
        TV ground_normal=TV(0,1,0);
        PLANE<T> water_plane(ground_normal,TV(0,water_level,0));
        bool intersect_plane=(extend_levelset_to_infinity && INTERSECTION::Intersects(plane_ray,water_plane,thickness));
        if(outside_box && intersect_plane && plane_ray.t_max<ray.t_max){
            ray.semi_infinite=plane_ray.semi_infinite;ray.t_max=plane_ray.t_max;ray.aggregate_id=-1;return true;}
 
        if(outside_box && !intersect_box) return false; // missed the box
        else if(outside_box){ // intersected the box from the outside
            TV point=ray.Point(ray.t_max); // intersection point with the box
            point=box.Thickened(-4*thickness).Clamp(point); // moves the point inside the box
            t_start=ray.t_max; // intersection with the box
            ray.semi_infinite=save_semi_infinite;ray.t_max=save_t_max;ray.aggregate_id=save_aggregate; // box intersection doesn't count
            RAY<TV> new_ray(point,ray.direction);INTERSECTION::Intersects(new_ray,box,thickness);
            if(ray.semi_infinite) t_end=t_start+new_ray.t_max;
            else t_end=min(t_start+new_ray.t_max,ray.t_max);
            exit_intersection=true;exit_t_max=t_end;exit_aggregate=new_ray.aggregate_id;} // save for exiting rays
        else if(!intersect_box){t_start=0;t_end=ray.t_max;} // intersected some object inside the box
        else{ // intersects the box from inside
            t_start=0;t_end=ray.t_max;
            exit_intersection=true;exit_t_max=ray.t_max;exit_aggregate=ray.aggregate_id; // save for exiting rays
            ray.semi_infinite=save_semi_infinite;ray.t_max=save_t_max;ray.aggregate_id=save_aggregate;} // box intersection doesn't count

        if(!use_secondary_interpolation){
            // set up marching
            IMPLICIT_OBJECT_ON_A_RAY<CHIMERA_PARTICLE_LEVELSET_IMPLICIT_OBJECT> implicit_surface_on_a_ray(*this,ray);
            ITERATIVE_SOLVER<T> iterative_solver;iterative_solver.tolerance=Iterative_Solver_Tolerance<T>()*thickness;
            // start marching
            T t1=t_start+thickness;
            T phi1=(*this)(ray.Point(t1));
            while(true){
                TV_INT background_index=background_grid.Clamp_To_Cell(ray.Point(t1));
                if(particle_count_touching_cell(background_index)==0){
                    while(particle_count_touching_cell(background_index)==0){
                        //LOG::cout<<"t1="<<t1<<" background_index="<<background_index<<std::endl;
                        T old_t1=t1;
                        Advance_To_Next_Background_Cell(ray,background_index,t1);
                        if(t1>t_end){t1=old_t1;break;}
                        background_index=background_grid.Clamp_To_Cell(ray.Point(t1));
                    }
                    t1-=background_grid.dX.Min()*(T)1e-6;
                    phi1=(*this)(ray.Point(t1));
                }
                T t2=t1+Integration_Step_World_Position(/*phi1,*/ray.Point(t1));
                if(t2>t_end) break;
                T phi2=(*this)(ray.Point(t2));
                if(LEVELSET_UTILITIES<T>::Interface(phi1,phi2)){ray.semi_infinite=false;ray.t_max=iterative_solver.Bisection_Secant_Root(implicit_surface_on_a_ray,t1,t2);ray.aggregate_id=-1;
                    return true;}
                else{t1=t2;phi1=phi2;}
            }
            T t2=t_end;T phi2=(*this)(ray.Point(t2));
            if(LEVELSET_UTILITIES<T>::Interface(phi1,phi2)){ray.semi_infinite=false;ray.t_max=iterative_solver.Bisection_Secant_Root(implicit_surface_on_a_ray,t1,t2);ray.aggregate_id=-1;return true;}
            else if(exit_intersection){
                return BASE::Infinite_Levelset_Intersection(ray,exit_t_max,thickness);
            }
            else return false;}
        else{ // use_secondary_interpolation
            PHYSBAM_NOT_IMPLEMENTED();}
    }
    
    void Get_Neighbor_Particles(const TV& location,T box_radius,ARRAY<TRIPLE<TV,TV,T> >& particle_info)
    {
        ARRAY<int> neighbor_particle_indices;
        RANGE<TV> range(location);range=range.Thickened(box_radius);
        particle_hierarchy.Intersection_List(range,neighbor_particle_indices);

        for(int i=1;i<=neighbor_particle_indices.Size();i++){
            particle_info.Append(Tuple(
                particle_positions(neighbor_particle_indices(i)),
                particle_velocities(neighbor_particle_indices(i)),
                particle_radii(neighbor_particle_indices(i))
            ));
        }
    }


    T operator()(const TV& location) const PHYSBAM_OVERRIDE
    {
        //return BASE::operator()(location);

        T water_distance=BASE::operator()(location);
        T particle_phi(0);

        ARRAY<int> neighbor_particle_indices;
        particle_hierarchy.Intersection_List(location,neighbor_particle_indices);

        for(int i=1;i<=neighbor_particle_indices.Size();i++){
            T radius_x=particle_radii_long(neighbor_particle_indices(i));
            T radius=particle_radii(neighbor_particle_indices(i));
            TV major_axis=particle_axis_directions(neighbor_particle_indices(i));
            TV dX=location-particle_positions(neighbor_particle_indices(i));
            T dist=sqrt(dX.Magnitude_Squared()+(sqr(radius/radius_x)-(T)1)*sqr(TV::Dot_Product(dX,major_axis)));
            if(dist<=radius) particle_phi+=T(-1)*kernel(dist/radius);
        }

        T water_phi(0);
        T r=water_distance/(average_particle_radius*levelset_kernel_scale)+root;
        if(r<(T)1){
            if(r<(T)0) water_phi=-1;
            else water_phi=-kernel(r);
        }
        return /*minimum_cell_size* */(blending_parameter+particle_phi+water_phi);
    }
    TV Normal(const TV& location,const int aggregate) const PHYSBAM_OVERRIDE
    {
        //return BASE::Normal(location,aggregate);

        if(aggregate!=-1) return grid_frames(grid_order(grids.m))->Rotation().Rotate(box.Normal(aggregate));
        T water_distance=BASE::operator()(location);
        TV particle_normal;

        ARRAY<int> neighbor_particle_indices;
        particle_hierarchy.Intersection_List(location,neighbor_particle_indices);

        for(int i=1;i<=neighbor_particle_indices.Size();i++){
            T radius_x=particle_radii_long(neighbor_particle_indices(i));
            T radius=particle_radii(neighbor_particle_indices(i));
            TV major_axis=particle_axis_directions(neighbor_particle_indices(i));
            TV dX=location-particle_positions(neighbor_particle_indices(i));
            T dot=TV::Dot_Product(dX,major_axis);
            T dist=sqrt(dX.Magnitude_Squared()+(sqr(radius/radius_x)-(T)1)*sqr(dot));
            TV normal=dX+(sqr(radius/radius_x)-(T)1)*dot*major_axis;
            if(dist<=radius) particle_normal+=(T(-1)*kernel.Prime(dist/radius)/radius/dist)*normal;
        }

        TV water_normal;
        T r=water_distance/(average_particle_radius*levelset_kernel_scale)+root;
        if(r>(T)0 && r<(T)1) water_normal=(-kernel.Prime(r)/(average_particle_radius*levelset_kernel_scale))*BASE::Normal(location,aggregate);
        return (particle_normal+water_normal).Normalized();
    }
};
}
#endif
