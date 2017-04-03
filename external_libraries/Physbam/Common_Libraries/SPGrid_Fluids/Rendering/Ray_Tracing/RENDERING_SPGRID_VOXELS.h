//#####################################################################
// Copyright (c) 2014, Mridul Aanjaneya.
// Distributed under the FreeBSD license (see license.txt)
//#####################################################################
// Class RENDERING_SPGRID_VOXELS
//##################################################################### 
#ifndef __RENDERING_SPGRID_VOXELS__
#define __RENDERING_SPGRID_VOXELS__
#include <PhysBAM_Geometry/Basic_Geometry_Intersections/RAY_PLANE_INTERSECTION.h>
#include <PhysBAM_Geometry/Implicit_Objects/IMPLICIT_OBJECT.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering/RENDER_WORLD.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering_Objects/RENDERING_VOXELS.h>
#include <PhysBAM_Tools/Data_Structures/HASHTABLE.h>
#include <PhysBAM_Tools/Grids_Uniform_Interpolation/INTERPOLATION_UNIFORM.h>
#include <PhysBAM_Tools/Grids_Uniform_Interpolation/LINEAR_INTERPOLATION_UNIFORM.h>
#include <SPGrid/Tools/SPGrid_Clear.h>
#include <SPGrid_Fluids/Advection/BACKTRACE.h>
#include <SPGrid_Fluids/Grids/GRID_HIERARCHY.h>
#include <SPGrid_Fluids/Grids/GRID_HIERARCHY_LOOKUP.h>
#include <SPGrid_Fluids/Grids/GRID_TOPOLOGY_HELPER.h>
#include <SPGrid_Fluids/Interpolation/GRID_HIERARCHY_AVERAGING.h>
#include <SPGrid_Fluids/Interpolation/GRID_HIERARCHY_INTERPOLATION.h>
#include <SPGrid_Fluids/Simulation/FLUIDS_SIMULATION_DATA.h>
namespace PhysBAM{

template<class T>
class RENDERING_SPGRID_VOXELS:public RENDERING_VOXELS<T>
{
public:
    using RENDERING_VOXELS<T>::box;using RENDERING_VOXELS<T>::small_number;using RENDERING_VOXELS<T>::precompute_single_scattering;
    typedef VECTOR<T,3> TV;typedef VECTOR<int,3> TV_INT;
    typedef FLUIDS_SIMULATION_DATA<T> T_STRUCT;typedef GRID<TV> T_GRID;
    typedef typename SPGrid_Allocator<T_STRUCT,3>::template Array<unsigned>::type Flag_array_type;
    typedef typename SPGrid_Allocator<T_STRUCT,3>::template Array<const unsigned>::type Const_flag_array_type;
    typedef typename SPGrid_Allocator<T_STRUCT,3>::template Array<T>::type Data_array_type;
    typedef typename SPGrid_Allocator<T_STRUCT,3>::template Array<const T>::type Const_data_array_type;
    GRID<TV>& fine_mac_grid;
    GRID_HIERARCHY<T_STRUCT,T,3>& hierarchy;
    T T_STRUCT::* density_channel;
    T T_STRUCT::* node_density_channel;
    ARRAY<ARRAY<VECTOR<T,3> ,VECTOR<int,3> >*> precomputed_light;
    ARRAY<ARRAY<bool,VECTOR<int,3> >*> precomputed_light_valid;
    ARRAY<VECTOR<int,3> > map_from_accessor_index_to_my_index;
    T volumetric_step;
    INTERPOLATION_UNIFORM<GRID<TV>,VECTOR<T,3> >* voxel_light_interpolation;
    LINEAR_INTERPOLATION_UNIFORM<GRID<TV>,VECTOR<T,3> > default_voxel_light_interpolation;
    int number_of_smoothing_steps;
    GRID<TV> light_grid;
    enum{nodes_per_cell=GRID_TOPOLOGY_HELPER<typename Flag_array_type::MASK>::nodes_per_cell};
    unsigned long nodes_of_cell_offsets[nodes_per_cell];

    RENDERING_SPGRID_VOXELS(GRID<TV> &fine_mac_grid_input,GRID_HIERARCHY<T_STRUCT,T,3>& hierarchy_input,T T_STRUCT::* density_channel_input,const T volumetric_step)
        :fine_mac_grid(fine_mac_grid_input),hierarchy(hierarchy_input),density_channel(density_channel_input),volumetric_step(volumetric_step),number_of_smoothing_steps(0)
    {
        box=fine_mac_grid.domain;
        node_density_channel=&T_STRUCT::ch8;
        T T_STRUCT::* temp_channel=&T_STRUCT::ch7;
        GRID_HIERARCHY_AVERAGING<T_STRUCT,T,3>::Average_Cell_Density_To_Nodes(hierarchy,density_channel,node_density_channel,&T_STRUCT::flags,temp_channel);
        voxel_light_interpolation=&default_voxel_light_interpolation;
        light_grid.Initialize(hierarchy.Grid(hierarchy.Levels()).counts,hierarchy.Grid(1).domain,true);
        GRID_TOPOLOGY_HELPER<typename Flag_array_type::MASK>::Nodes_Of_Cell_Offsets(nodes_of_cell_offsets);
    }

    void Get_Hits_With_Plane(const PLANE<T>& plane,RAY<VECTOR<T,3> > ray,T &min_hit,T &second_min_hit) const
    {if(INTERSECTION::Intersects(ray,plane) && ray.t_max<min_hit){second_min_hit=min_hit;min_hit=ray.t_max;}}
        
    T Volumetric_Integration_Step(const RAY<VECTOR<T,3> > &ray,const T xi) const PHYSBAM_OVERRIDE
    {return xi*volumetric_step;}

    T Source_Term(const int source_term_index,const VECTOR<T,3>& location) const PHYSBAM_OVERRIDE
    {unsigned long offset;int level;TV weights;
    bool ret=GRID_HIERARCHY_LOOKUP<T_STRUCT,T,3>::Cell_Lookup(hierarchy,location,offset,level,weights);
    if(!ret) return (T)0.;
    return GRID_HIERARCHY_INTERPOLATION<T_STRUCT,T,3>::Cell_Interpolation_Helper(hierarchy,const_cast<unsigned long *>(nodes_of_cell_offsets),level,offset,weights,density_channel,node_density_channel);}

    void Get_Node_Locations(ARRAY<VECTOR<T,3> >& locations) PHYSBAM_OVERRIDE
    {locations.Resize(light_grid.counts.Product());map_from_accessor_index_to_my_index.Resize(locations.m);int index=1;
    for(int i=1;i<=light_grid.counts.x;i++)for(int j=1;j<=light_grid.counts.y;j++)for(int ij=1;ij<=light_grid.counts.z;ij++){
        map_from_accessor_index_to_my_index(index)=VECTOR<int,3>(i,j,ij);
        locations(index)=light_grid.X(i,j,ij);index++;}}

    bool Use_Precomputed_Light_Data(const VECTOR<T,3>& location,const int light_index) const PHYSBAM_OVERRIDE
    {if(!precompute_single_scattering)return false;
    VECTOR<int,3> index=INTERPOLATION_UNIFORM<GRID<TV>,VECTOR<T,3> >::Clamped_Index_End_Minus_One(light_grid,*precomputed_light(light_index),location);
    int i=index.x,j=index.y,ij=index.z;
    if(light_grid.Outside(location))return false;
    bool i_j_ij=(*precomputed_light_valid(light_index))(i,j,ij),i_j_ij1=(*precomputed_light_valid(light_index))(i,j,ij+1),i_j1_ij=(*precomputed_light_valid(light_index))(i,j+1,ij),
        i_j1_ij1=(*precomputed_light_valid(light_index))(i,j+1,ij+1),i1_j_ij=(*precomputed_light_valid(light_index))(i+1,j,ij),i1_j_ij1=(*precomputed_light_valid(light_index))(i+1,j,ij+1),
        i1_j1_ij=(*precomputed_light_valid(light_index))(i+1,j+1,ij),i1_j1_ij1=(*precomputed_light_valid(light_index))(i+1,j+1,ij+1);
        return i_j_ij&&i_j_ij1&&i_j1_ij&&i_j1_ij1&&i1_j_ij&&i1_j_ij1&&i1_j1_ij&&i1_j1_ij1;}

    void Set_Precomputed_Light_Data(const int location_index,const int light_index,const VECTOR<T,3>& light_value) PHYSBAM_OVERRIDE
    {(*precomputed_light(light_index))(map_from_accessor_index_to_my_index(location_index))=light_value;}

    void Set_Precomputed_Light_Valid(const int location_index,const int light_index,const bool value) PHYSBAM_OVERRIDE
    {VECTOR<int,3> index=map_from_accessor_index_to_my_index(location_index);
    (*precomputed_light_valid(light_index))(index)=value;}

    VECTOR<T,3> Precomputed_Light_Data(const VECTOR<T,3>& location,const int light) const PHYSBAM_OVERRIDE
    {return voxel_light_interpolation->Clamped_To_Array(light_grid,*precomputed_light(light),location);}

    void Set_Custom_Light_Interpolation(INTERPOLATION_UNIFORM<GRID<TV>,VECTOR<T,3> >* interpolation)
    {voxel_light_interpolation=interpolation;}

protected:
    void Prepare_For_Precomputation(RENDER_WORLD<T>& world) PHYSBAM_OVERRIDE
    {precomputed_light.Resize(world.Lights().m);precomputed_light_valid.Resize(world.Lights().m);
    for(int i=1;i<=precomputed_light.m;i++)precomputed_light(i)=new ARRAY<VECTOR<T,3> ,VECTOR<int,3> >(1,light_grid.counts.x,1,light_grid.counts.y,1,light_grid.counts.z);
    for(int i=1;i<=precomputed_light.m;i++){precomputed_light_valid(i)=new ARRAY<bool,VECTOR<int,3> >(1,light_grid.counts.x,1,light_grid.counts.y,1,light_grid.counts.z);precomputed_light_valid(i)->Fill(false);}}

//#####################################################################
    void Postprocess_Light_Field() PHYSBAM_OVERRIDE;
//#####################################################################
};
}
#endif
