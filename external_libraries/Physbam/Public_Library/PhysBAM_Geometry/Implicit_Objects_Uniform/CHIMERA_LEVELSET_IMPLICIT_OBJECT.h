//#####################################################################
// Copyright 2002-2007, Doug Enright, Ron Fedkiw, Eran Guendelman, Geoffrey Irving, Sergey Koltakov, Neil Molino, Andrew Selle, Eftychios Sifakis, Rachel Weinstein.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class CHIMERA_LEVELSET_IMPLICIT_OBJECT
//#####################################################################
#ifndef __CHIMERA_LEVELSET_IMPLICIT_OBJECT__
#define __CHIMERA_LEVELSET_IMPLICIT_OBJECT__

#include <PhysBAM_Tools/Matrices/FRAME.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_1D.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_2D.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_3D.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_POLICY_UNIFORM.h>
#include <PhysBAM_Geometry/Implicit_Objects/IMPLICIT_OBJECT.h>
#include <PhysBAM_Geometry/Implicit_Objects_Uniform/LEVELSET_IMPLICIT_OBJECT.h>
#include <PhysBAM_Geometry/Basic_Geometry_Intersections/RAY_BOX_INTERSECTION.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Geometry/Level_Sets/IMPLICIT_OBJECT_ON_A_RAY.h>
#include <PhysBAM_Geometry/Level_Sets/LEVELSET_UTILITIES.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/FAST_MARCHING_METHOD_UNIFORM.h>
#include <PhysBAM_Tools/Nonlinear_Equations/ITERATIVE_SOLVER.h>
#include <PhysBAM_Geometry/Basic_Geometry_Intersections/RAY_PLANE_INTERSECTION.h>
#include <PhysBAM_Geometry/Basic_Geometry/PLANE.h>
#include <PhysBAM_Geometry/Basic_Geometry/ORIENTED_BOX.h>

namespace PhysBAM{

template<class TV> class GRID;

//#####################################################################
// Function Iterative_Solver_Tolerance
//#####################################################################
namespace{
template<class T> inline T Iterative_Solver_Tolerance(){STATIC_ASSERT((T)false);}
template<> inline float Iterative_Solver_Tolerance(){return (float).01;}
template<> inline double Iterative_Solver_Tolerance(){return .001;}
}

template<class TV>
class CHIMERA_LEVELSET_IMPLICIT_OBJECT:public LEVELSET_IMPLICIT_OBJECT<TV>
{
    typedef typename TV::SCALAR T;
    typedef typename TV::template REBIND<int>::TYPE TV_INT;
    //typedef typename GRID_ARRAYS_POLICY<GRID<TV> >::ARRAYS_SCALAR T_ARRAYS_SCALAR;typedef typename T_ARRAYS_SCALAR::template REBIND<TV>::TYPE T_ARRAYS_VECTOR;
    typedef typename LEVELSET_POLICY<GRID<TV> >::LEVELSET T_LEVELSET;
    //typedef typename MATRIX_POLICY<TV>::SYMMETRIC_MATRIX T_SYMMETRIC_MATRIX;
    //enum WORKAROUND {d=TV::m};
    //typedef VECTOR<T,d-1> T_PRINCIPAL_CURVATURES;
public:
    typedef LEVELSET_IMPLICIT_OBJECT<TV> BASE;
    using BASE::box;using BASE::levelset;using BASE::minimum_cell_size;
    using IMPLICIT_OBJECT<TV>::use_secondary_interpolation;

    ARRAY<GRID<TV>*> grids;
    ARRAY<FRAME<TV>*> grid_frames;
    ARRAY<ARRAY<T,TV_INT>*> levelset_data;
    ARRAY<T_LEVELSET*> levelsets;

    int number_of_ghost_cells;
    int levelset_number_of_ghost_cells;
    T blend_size;
    bool extend_levelset_to_infinity;
    T water_level;

    BOX<TV> old_box;
protected:
    ARRAY<BOX<TV> > interpolation_domains;
    ARRAY<ORIENTED_BOX<TV> > world_space_domains;
    ARRAY<int> grid_order;
    LINEAR_INTERPOLATION_UNIFORM<GRID<TV>,T> interpolation;
public:

    CHIMERA_LEVELSET_IMPLICIT_OBJECT(const ARRAY<GRID<TV>*>& grids_input,const ARRAY<FRAME<TV>*>& grid_frames_input,const ARRAY<ARRAY<T,TV_INT>*>& levelset_data_input,T blendsz,int levelset_ghost_cells)
    :LEVELSET_IMPLICIT_OBJECT<TV>(*grids_input(1),*levelset_data_input(1)),grids(grids_input),grid_frames(grid_frames_input),levelset_data(levelset_data_input),number_of_ghost_cells(2),levelset_number_of_ghost_cells(levelset_ghost_cells),blend_size(blendsz),extend_levelset_to_infinity(false),water_level((T)0)
    {
        LOG::cerr<<"levelset number of ghost cells"<<levelset_number_of_ghost_cells<<std::endl;
        LOG::cerr<<"chimera blend_size="<<blend_size<<std::endl;
        interpolation_domains.Resize(grids.m);
        world_space_domains.Resize(grids.m);
        ARRAY<PAIR<T,int> > grid_sizes(grids.m);

        for(int grid_index=1;grid_index<=grids.m;grid_index++){
            interpolation_domains(grid_index)=BOX<TV>(Interpolation_Domain(grid_index,(T)-1));
            world_space_domains(grid_index)=ORIENTED_BOX<TV>(grids(grid_index)->Domain(),*grid_frames(grid_index));
            T dx=grids(grid_index)->DX().Magnitude();
            grid_sizes(grid_index)=Tuple(dx,grid_index);
            
            LOG::cout << "grid " << grid_index << " dx=" << grids(grid_index)->DX()(1) << " grid frame="<< *grid_frames(grid_index) << " interpolation domain=" << interpolation_domains(grid_index) << std::endl;
        }
        
        Stable_Sort(grid_sizes);
        grid_order.Resize(grids.m);
        for(int i=1;i<=grids.m;i++)
            grid_order(i)=grid_sizes(i).y;
        LOG::cout << "sorted grids " << grid_order << "  grids.Size()=" << grids.Size() <<" = " << grids.m << std::endl;

        minimum_cell_size=grids(grid_order(1))->dX(1);
        box=grids(grid_order(grids.m))->Domain();
        LOG::cout<<"box="<<box<<" min_cell_size="<<minimum_cell_size<<std::endl;

        for(int i=1;i<=grids.Size();++i){ //fill ghost cells (only for smoothing)
            int grid_index=grid_order(i);
            for(UNIFORM_GRID_ITERATOR_CELL<TV> it(*grids(grid_index),levelset_number_of_ghost_cells);it.Valid();it.Next()){
                if(!grids(grid_index)->Inside_Domain(it.Cell_Index())){
                    TV world_location=grid_frames(grid_index)->operator*(it.Location());
                    for(int j=1;j<=grids.Size();++j){
                        if(j==i) continue;
                        int finer_grid_index=grid_order(j);
                        TV location_local=grid_frames(finer_grid_index)->Inverse_Times(world_location);
                        if(Interpolation_Domain(finer_grid_index).Lazy_Inside(location_local)){
                            levelset_data(grid_index)->operator()(it.Cell_Index())=interpolation.Clamped_To_Array(*grids(finer_grid_index),*levelset_data(finer_grid_index),location_local);
                            break;}
                    }
                }
            }
        }
        for(int i=1;i<=grids.m;i++) levelsets.Append(new T_LEVELSET(*grids(i),*levelset_data(i)));
        for(int i=1;i<=levelsets.m;++i){
            LOG::cout<<"levelset "<<i<<" domain indices="<<levelset_data(i)->Domain_Indices()<<std::endl;
            levelset.Fast_Marching_Method(FLT_MAX,0);
            //levelset.Fast_Marching_Method((T)8*grids(i)->dX(1),0);
        }
        old_box=box;
    }
    void Smooth(T radius)
    {
        for(int i=1;i<=grids.Size();++i){
            int grid_index=grid_order(i);
            //T factor=grids(grid_index)->dX(1)/(radius*grids(grid_order(1))->dX(1));
            T factor=(T)1/radius;
            if(factor>(T)2.5) continue;
            //RANGE<TV_INT> cell_range=levelset_data(grid_index)->Domain_Indices();
            ARRAY<T,TV_INT> phi_ghost(*levelset_data(grid_index));
            //for(UNIFORM_GRID_ITERATOR_CELL<TV> it(*grids(grid_index),levelset_number_of_ghost_cells);it.Valid();it.Next()){
            //    phi_ghost(it.Cell_Index())=levelset_data(grid_index)->operator()(it.Cell_Index());
            //}
            //for(UNIFORM_GRID_ITERATOR_CELL<TV> it(*grids(grid_index),3);it.Valid();it.Next()){
            //    if(!grids(grid_index)->Inside_Domain(it.Cell_Index())){
            //        LOG::cout<<"grid "<<grid_index<<" order "<<i<<" ghost cell at "<<it.Cell_Index()<<" = "<<levelset_data(grid_index)->operator()(it.Cell_Index())<<std::endl;
            //    }
            //}
            for(UNIFORM_GRID_ITERATOR_CELL<TV> it(*grids(grid_index));it.Valid();it.Next()){
                T accum=0.0;
                T normalization_weight=0.0;
                for(int i=-2;i<=2;i++){
                    for(int j=-2;j<=2;j++){
                        for(int k=-2;k<=2;k++){
                            TV_INT delta(i,j,k);
                            T s=factor*sqrt((T)delta.Magnitude_Squared());
                            T cell_weight;
                            if(0.0<=s && s<0.5){
                                cell_weight=((s*s*s*s)/4.0)-((5.0*s*s)/8.0)+(115.0/192.0);}
                            else if(0.5<=s && s<1.5){
                                cell_weight=-((s*s*s*s)/6.0)+((5.0*s*s*s)/6.0)-((5.0*s*s)/4.0)+((5.0*s)/24.0)+(55.0/96.0);}
                            else if(1.5<=s && s<2.5){
                                cell_weight=(((2.5-s)*(2.5-s)*(2.5-s)*(2.5-s))/24.0);}
                            else{
                                cell_weight=0.0;}
                            if(cell_weight!=0.0){
                                normalization_weight+=cell_weight;
                                accum+=cell_weight*phi_ghost(it.Cell_Index()+delta);}}}}
                levelset_data(grid_index)->operator()(it.Cell_Index())=accum/normalization_weight;
            }
        }
        //for(int i=1;i<=levelsets.m;++i){
        //    LOG::cout<<"levelset "<<i<<" domain indices="<<levelset_data(i)->Domain_Indices()<<std::endl;
        //    FAST_MARCHING_METHOD_UNIFORM<GRID<TV> > fmm(*levelsets(i),levelset_number_of_ghost_cells,0);
        //    fmm.Fast_Marching_Method(*levelset_data(i),0,0,false,false,true);
        //}
    }
    T Integration_Step(const T phi) const PHYSBAM_OVERRIDE
    {
        return BASE::Integration_Step(phi);
    }

    void Update_Box()
    {
        old_box=box;//.Thickened((T)-1.5*grids(grid_order(grids.m))->dX(1));
        interpolation_domains(grid_order(grids.m))=old_box.Thickened((T)-1.5*grids(grid_order(grids.m))->dX(1));
    }

    inline T Local_Phi(int grid_index,const TV& location_local) const
    {
        //return levelsets(grid_index)->Phi((extend_levelset_to_infinity && grid_index==grid_order(grids.m))?old_box.Clamp(location_local):location_local);
       return levelsets(grid_index)->Phi(location_local);
    }
    inline TV Local_Normal(int grid_index,const TV& location) const
    {
        const TV& dx=grids(grid_index)->dX;
        return TV(
            (Local_Phi(grid_index,TV(location.x+dx.x,location.y,location.z))-Local_Phi(grid_index,TV(location.x-dx.x,location.y,location.z)))/(2*dx.x),
            (Local_Phi(grid_index,TV(location.x,location.y+dx.y,location.z))-Local_Phi(grid_index,TV(location.x,location.y-dx.y,location.z)))/(2*dx.y),
            (Local_Phi(grid_index,TV(location.x,location.y,location.z+dx.z))-Local_Phi(grid_index,TV(location.x,location.y,location.z-dx.z)))/(2*dx.z)).Normalized();
    }

    T Signed_Distance(const TV& location) const PHYSBAM_OVERRIDE
    {
        PHYSBAM_FATAL_ERROR("signed distance function should not be called!");
        T abs_value(FLT_MAX);
        T value(FLT_MAX);
        for(int i=1;i<=grids.Size();i++){
            int grid_index=grid_order(i);
            TV location_local=grid_frames(grid_index)->Inverse_Times(location);
            T signed_distance=interpolation_domains(grid_index).Signed_Distance(location_local);
            if(signed_distance<0){
                T dist=levelsets(grid_index)->Phi(location_local);
                if(abs(dist)<abs_value){value=dist;abs_value=abs(dist);}
            }
        }
        return value;
    }
    T operator()(const TV& location) const PHYSBAM_OVERRIDE
    {
        T value(0);
        T weight_complement=1,weight;
        for(int i=1;i<=grids.Size();i++){
            int grid_index=grid_order(i);
            TV location_local=grid_frames(grid_index)->Inverse_Times(location);
            T signed_distance=interpolation_domains(grid_index).Signed_Distance(location_local);
            if(signed_distance<0 || (!extend_levelset_to_infinity && i==grids.Size())){
                weight=blend_size<1e-5?1:min((T)1,-signed_distance/(blend_size*grids(grid_index)->dX(1)));
                value+=weight_complement*weight*Local_Phi(grid_index,location_local);
                weight_complement*=(1-weight);
                if(weight_complement<.0001) break;}}
        if(extend_levelset_to_infinity) value+=weight_complement*(location.y-water_level);
        
        return value;
    }
    TV Normal(const TV& location,const int aggregate) const PHYSBAM_OVERRIDE
    {
        if(aggregate!=-1) return box.Normal(aggregate);
        TV sum=TV();
        T weight_complement=1,weight;
        for(int i=1;i<=grids.Size();i++){
            int grid_index=grid_order(i);
            TV location_local=grid_frames(grid_index)->Inverse_Times(location);
            T signed_distance=interpolation_domains(grid_index).Signed_Distance(location_local);
            if(signed_distance<0 || (!extend_levelset_to_infinity && i==grids.Size())){
                weight=blend_size<1e-5?1:min((T)1,-signed_distance/(blend_size*grids(grid_index)->dX(1)));
                sum+=weight_complement*weight*grid_frames(grid_index)->Rotation().Rotate(Local_Normal(grid_index,location_local));
                weight_complement*=(1-weight);
                if(weight_complement<.0001) break;}}
        if(extend_levelset_to_infinity) sum+=weight_complement*TV(0,1,0);
        
        return sum.Normalized();
    }

    RANGE<TV> Interpolation_Domain(int grid_index,T n_ghost_cells=0) const PHYSBAM_OVERRIDE
    {
        T thicken_cells=n_ghost_cells-(T).5;
        RANGE<TV> domain=grids(grid_index)->domain;
        TV dx=grids(grid_index)->dX;
        domain.Change_Size(thicken_cells*dx);
        return domain;
    }

    bool Inside(const TV& location,const T thickness_over_two) const
    {
        return (extend_levelset_to_infinity || box.Inside(grid_frames(1)->Inverse_Times(location),thickness_over_two)) && (*this)(location)<=-thickness_over_two;
    }

    bool Outside(const TV& location,const T thickness_over_two) const
    {
        return (!extend_levelset_to_infinity && box.Outside(grid_frames(1)->Inverse_Times(location),thickness_over_two)) || (*this)(location)>=thickness_over_two;
    }
    
    inline bool Infinite_Levelset_Intersection(RAY<TV>& ray,const T exit_t_max,const T thickness) const
    {
        bool save_semi_infinite=ray.semi_infinite;T save_t_max=ray.t_max;int save_aggregate=ray.aggregate_id;
        TV ground_normal=TV(0,1,0);
        PLANE<T> water_plane(ground_normal,TV(0,water_level,0));
        if(INTERSECTION::Intersects(ray,water_plane,thickness) && ray.t_max>exit_t_max){
            ray.aggregate_id=-1;
            return true;
        }else{
            ray.semi_infinite=save_semi_infinite;ray.t_max=save_t_max;ray.aggregate_id=save_aggregate;
            return false;
        }
    }

    inline int Finest_Grid_Order(const TV& location) const
    {
        int finest_grid_order=grid_order(grids.m);
        for(int i=1;i<=grids.m;++i){
            int grid_index=grid_order(i);
            if(world_space_domains(grid_index).Lazy_Inside(location)){finest_grid_order=i;break;}
        }
        return finest_grid_order;
    }
    inline T Integration_Step_World_Location(const RAY<TV>& ray,const T phi,const T t_ray,const T thickness) const
    {
        return BASE::Integration_Step(phi);

        T distance=abs(phi);
        int start_finest_grid_order=Finest_Grid_Order(ray.Point(t_ray));
        T local_cell_size=grids(grid_order(start_finest_grid_order))->dX(1);
        T trial_distance=0;
        if(distance > 3*local_cell_size) trial_distance=(T).5*distance;
        else if(distance > local_cell_size) trial_distance=(T).25*distance;
        else trial_distance=(T).1*local_cell_size;
        //int end_finest_grid_order=Finest_Grid_Order(ray.Point(t_ray+trial_distance));
        //if(end_finest_grid_order<start_finest_grid_order){
        //    // object space
        //    RAY<TV> box_ray(grid_frames(grid_order(end_finest_grid_order))->Inverse_Times(ray.Point(trial_distance)),
        //                    grid_frames(grid_order(end_finest_grid_order))->Rotation().Inverse_Rotate(ray.direction));
        //    box_ray.semi_infinite=ray.semi_infinite;box_ray.t_max=ray.t_max-trial_distance;box_ray.aggregate_id=ray.aggregate_id;
        //    if(box_ray.t_max>0){
        //        if(INTERSECTION::Intersects(box_ray,grids(grid_order(end_finest_grid_order))->domain,thickness)){return box_ray.t_max+4*thickness;}
        //    }
        //}
        return trial_distance;
    }
    inline void Compute_Box_Intersections(const RAY<TV>& ws_ray, ARRAY<RAY<TV> >& object_space_rays, ARRAY<INTERVAL<T> >& box_intersections, ARRAY<INTERVAL<T> >& interior_box_intersections, const T t_start, const T t_end) const
    {
        for(int i=1;i<=grids.m;++i){
            int grid_index=grid_order(i);
            RAY<TV> ray(grid_frames(grid_index)->Inverse_Times(ws_ray.Point(t_start)),grid_frames(grid_index)->Rotation().Inverse_Rotate(ws_ray.direction));
            ray.semi_infinite=ws_ray.semi_infinite;ray.t_max=t_end-t_start;ray.aggregate_id=ws_ray.aggregate_id;
            object_space_rays(i)=ray;

            BOX<TV> box=interpolation_domains(grid_index);
            T box_start_t=FLT_MAX,box_end_t=-FLT_MAX;
            INTERSECTION::Get_Intersection_Range(ray,box,box_start_t,box_end_t);
            box_intersections(i)=INTERVAL<T>(box_start_t,box_end_t);

            BOX<TV> interior_box=interpolation_domains(grid_index).Thickened(-blend_size*grids(grid_index)->dX(1));
            T interior_box_start_t=FLT_MAX,interior_box_end_t=-FLT_MAX;
            INTERSECTION::Get_Intersection_Range(ray,interior_box,interior_box_start_t,interior_box_end_t);
            interior_box_intersections(i)=INTERVAL<T>(interior_box_start_t,interior_box_end_t);
        }
    }
    T Source_Term_On_Ray(const TV& location,const T t,const ARRAY<RAY<TV> >& object_space_rays,const ARRAY<INTERVAL<T> >& box_intersections,const ARRAY<INTERVAL<T> >& interior_box_intersections) const PHYSBAM_OVERRIDE
    {
        int finest_order=grids.Size();
        int interior_finest_order=-1;

        for(int i=1;i<grids.Size();i++) if(box_intersections(i).Lazy_Inside(t)) {finest_order=i;break;}
        for(int i=1;i<=grids.Size();i++) if(interior_box_intersections(i).Lazy_Inside(t)) {interior_finest_order=i;break;}

        if(finest_order!=interior_finest_order) return (*this)(location);

        int grid_index=grid_order(finest_order);
        TV location_local=object_space_rays(finest_order).Point(t);
        return Local_Phi(grid_index,location_local);
    }
    bool Intersection(RAY<TV>& ray,const T thickness) const
    {
        bool save_semi_infinite=ray.semi_infinite;T save_t_max=ray.t_max;int save_aggregate=ray.aggregate_id;

        T t_start,t_end;
        bool exit_intersection=false;int exit_aggregate=0;T exit_t_max=0;
        int intersect_box=INTERSECTION::Intersects(ray,old_box,thickness);
        int outside_box=old_box.Outside(ray.endpoint,thickness);

        if(extend_levelset_to_infinity){
            RAY<TV> plane_ray(ray.endpoint,ray.direction);plane_ray.semi_infinite=save_semi_infinite;plane_ray.t_max=save_t_max;plane_ray.aggregate_id=save_aggregate;
            TV ground_normal=TV(0,1,0);
            PLANE<T> water_plane(ground_normal,TV(0,water_level,0));
            bool intersect_plane=INTERSECTION::Intersects(plane_ray,water_plane,thickness);
            if(outside_box && intersect_plane && (!intersect_box || plane_ray.t_max<ray.t_max)){
                ray.semi_infinite=plane_ray.semi_infinite;ray.t_max=plane_ray.t_max;ray.aggregate_id=-1;return true;}
        }
 
        if(outside_box && !intersect_box) return false; // missed the box
        else if(outside_box){ // intersected the box from the outside
            TV point=ray.Point(ray.t_max); // intersection point with the box
            point=old_box.Thickened(-4*thickness).Clamp(point); // moves the point inside the box
            if(!extend_levelset_to_infinity && (*this)(point) <= 0){return true;}// level set is on the edge of the box
            t_start=ray.t_max; // intersection with the box
            ray.semi_infinite=save_semi_infinite;ray.t_max=save_t_max;ray.aggregate_id=save_aggregate; // box intersection doesn't count
            RAY<TV> new_ray(point,ray.direction);INTERSECTION::Intersects(new_ray,box,thickness);
            if(ray.semi_infinite) t_end=t_start+new_ray.t_max;
            else t_end=min(t_start+new_ray.t_max,ray.t_max);
            if(extend_levelset_to_infinity){exit_intersection=true;exit_t_max=t_end;exit_aggregate=new_ray.aggregate_id;}} // save for exiting rays
        else if(!intersect_box){t_start=0;t_end=ray.t_max;} // intersected some object inside the box
        else{ // intersects the box from inside
            t_start=0;t_end=ray.t_max;
            exit_intersection=true;exit_t_max=ray.t_max;exit_aggregate=ray.aggregate_id; // save for exiting rays
            ray.semi_infinite=save_semi_infinite;ray.t_max=save_t_max;ray.aggregate_id=save_aggregate;} // box intersection doesn't count
    
        if(!use_secondary_interpolation){
            ARRAY<RAY<TV> > object_space_rays(grids.m);
            ARRAY<INTERVAL<T> > box_intersections(grids.m);
            ARRAY<INTERVAL<T> > interior_box_intersections(grids.m);
            Compute_Box_Intersections(ray,object_space_rays,box_intersections,interior_box_intersections,t_start,t_end);
            // set up marching
            IMPLICIT_OBJECT_ON_A_RAY<CHIMERA_LEVELSET_IMPLICIT_OBJECT> implicit_surface_on_a_ray(*this,ray);
            ITERATIVE_SOLVER<T> iterative_solver;iterative_solver.tolerance=Iterative_Solver_Tolerance<T>()*thickness;
            // start marching
            T t1=t_start+thickness,phi1=(*this)(ray.Point(t1)),t2=t1+Integration_Step_World_Location(ray,phi1,t1,thickness);
            // march through the line segment
            int n_steps=0;
            while(t2 <= t_end){
                ++n_steps;
                //T phi2=(*this)(ray.Point(t2));
                T phi2=Source_Term_On_Ray(ray.Point(t2),t2-t_start,object_space_rays,box_intersections,interior_box_intersections);
                //std::cout<<"lqiu debug ray point "<<ray.Point(t2)<<std::endl;
                if(LEVELSET_UTILITIES<T>::Interface(phi1,phi2)){ray.semi_infinite=false;ray.t_max=iterative_solver.Bisection_Secant_Root(implicit_surface_on_a_ray,t1,t2);ray.aggregate_id=-1;
                //LOG::cout<<"n_steps = "<<n_steps<<std::endl;
                return true;}
                else{t1=t2;phi1=phi2;t2=t1+Integration_Step_World_Location(ray,phi1,t1,thickness);}}
            //LOG::cout<<"n_steps = "<<n_steps+1<<std::endl;
            // check the last piece of the line segment
            t2=t_end;T phi2=(*this)(ray.Point(t_end));
            if(LEVELSET_UTILITIES<T>::Interface(phi1,phi2)){ray.semi_infinite=false;ray.t_max=iterative_solver.Bisection_Secant_Root(implicit_surface_on_a_ray,t1,t2);ray.aggregate_id=-1;return true;}
            else if(exit_intersection){
                if(extend_levelset_to_infinity) return Infinite_Levelset_Intersection(ray,exit_t_max,thickness);
                else if(phi2 <= 0){ray.semi_infinite=false;ray.t_max=exit_t_max;ray.aggregate_id=exit_aggregate;return true;} // exiting ray
                else return false;
            }
            else return false;
        }
        else{
            PHYSBAM_NOT_IMPLEMENTED();
        }
    }
};
}
#endif
