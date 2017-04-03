//#####################################################################
// Copyright 2014.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __RENDERING_SPOTLIGHT_TEXTURED__
#define __RENDERING_SPOTLIGHT_TEXTURED__

#include <PhysBAM_Tools/Matrices/MATRIX_3X3.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering/RENDER_WORLD.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering_Lights/RENDERING_LIGHT.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_NODE.h>
#include <PhysBAM_Tools/Grids_Uniform_Boundaries/BOUNDARY_UNIFORM.h>
#include <PhysBAM_Tools/Grids_Uniform_Interpolation/LINEAR_INTERPOLATION_UNIFORM.h>
#include <PhysBAM_Tools/Images/IMAGE.h>
#include <cstdio>
namespace PhysBAM{

template<class T>
class RENDERING_SPOTLIGHT_TEXTURED:public RENDERING_LIGHT<T>
{
public:
    using RENDERING_LIGHT<T>::world;using RENDERING_LIGHT<T>::position;using RENDERING_LIGHT<T>::color;using RENDERING_LIGHT<T>::brightness;using RENDERING_LIGHT<T>::global_photon_random;
    using RENDERING_LIGHT<T>::caustic_photon_random;using RENDERING_LIGHT<T>::volume_photon_random;using RENDERING_LIGHT<T>::sample_points_random;
    using RENDERING_LIGHT<T>::supports_global_photon_mapping;using RENDERING_LIGHT<T>::supports_caustic_photon_mapping;using RENDERING_LIGHT<T>::supports_volume_photon_mapping;

    T cone_angle,penumbra_angle;
    T cos_cone_angle,cos_penumbra_angle;
    T power_to_intensity;
    VECTOR<T,3> spot_direction;
    MATRIX<T,3> local_vector_to_world_vector;
    GRID<VECTOR<T,2> > grid;
    ARRAY<VECTOR<T,3> ,VECTOR<int,2> > pixels;
    LINEAR_INTERPOLATION_UNIFORM<GRID<VECTOR<T,2> >,VECTOR<T,3> > linear_interpolation;
    
    RENDERING_SPOTLIGHT_TEXTURED(const VECTOR<T,3>& position_input,const VECTOR<T,3>& color_input,const T brightness_input,const VECTOR<T,3>& direction_input,const T cone_angle_input,
        const T penumbra_angle_input,RENDER_WORLD<T>& world_input,const bool supports_global_photons,const bool supports_caustic_photons,const bool supports_volume_photons,
        const bool photon_source_only,const std::string texture_filename)
        :RENDERING_LIGHT<T>(position_input,color_input,brightness_input,world_input,supports_global_photons,supports_caustic_photons,supports_volume_photons,photon_source_only),
        cone_angle(cone_angle_input),penumbra_angle(penumbra_angle_input),spot_direction(direction_input)
    {
        spot_direction.Normalize();
        cos_cone_angle=cos(cone_angle);cos_penumbra_angle=cos(penumbra_angle_input);
        power_to_intensity=T(1)/(2*T(pi)*(1-T(.5)*(cos_cone_angle+cos_penumbra_angle))); // the inverse of the solid angle covered by a cone halfway between penumbra begin and full width
        // setup matrix to transform the canonical hemisphere directions (probably should change this to a quaternion)
        VECTOR<T,3> temporary_vector=spot_direction.Orthogonal_Vector();
        VECTOR<T,3> u=VECTOR<T,3>::Cross_Product(temporary_vector,spot_direction);u.Normalize();
        VECTOR<T,3> v=VECTOR<T,3>::Cross_Product(spot_direction,u);
        local_vector_to_world_vector=MATRIX<T,3>(u,v,spot_direction);

        IMAGE<T>::Read(texture_filename,pixels);
        int ghost_cells=3;
        grid.Initialize(pixels.counts.x,pixels.counts.y,RANGE<VECTOR<T,2> >::Unit_Box(),true);pixels.Resize(grid.Domain_Indices(ghost_cells));
        BOUNDARY_UNIFORM<GRID<VECTOR<T,2> >,VECTOR<T,3> >().Fill_Ghost_Cells(grid,pixels,pixels,0,0,ghost_cells);
    }

    void Sample_Points(const VECTOR<T,3>& surface_position,const VECTOR<T,3>& surface_normal,ARRAY<RAY<VECTOR<T,3> > >& sample_array) const PHYSBAM_OVERRIDE
    {sample_array.Resize(1);sample_array(1)=RAY<VECTOR<T,3> >(SEGMENT_3D<T>(surface_position,position));}

    T Spotlight_Falloff(const VECTOR<T,3>& direction) const
    {T cos_theta=VECTOR<T,3>::Dot_Product(direction,spot_direction);
    if(cos_theta<cos_cone_angle)return 0;
    else if(cos_theta>cos_penumbra_angle)return 1;
    else {T alpha=(cos_theta-cos_cone_angle)/(cos_penumbra_angle-cos_cone_angle);return alpha*alpha*alpha*alpha;}}

    T Texture_Falloff(const VECTOR<T,3>& direction) const
    {
        VECTOR<T,3> u_direction=VECTOR<T,3>::Cross_Product(VECTOR<T,3>(1,0,0),spot_direction).Normalized();
        VECTOR<T,3> v_direction=VECTOR<T,3>::Cross_Product(spot_direction,u_direction).Normalized();
        T u=VECTOR<T,3>::Dot_Product(direction,u_direction);
        T v=VECTOR<T,3>::Dot_Product(direction,v_direction);
        PHYSBAM_ASSERT(cos_cone_angle>0);
        T sin_cone_angle=sqrt(1-sqr(cos_cone_angle));
        VECTOR<T,2> texcoord=VECTOR<T,2>(u,v)/sin_cone_angle*(T).5+(T).5;
        return linear_interpolation.Clamped_To_Array(grid,pixels,texcoord)(1);
    }

    VECTOR<T,3> Emitted_Light(const RENDERING_RAY<T>& ray) const PHYSBAM_OVERRIDE
    {return brightness*color*power_to_intensity*Texture_Falloff(-ray.ray.direction)*Spotlight_Falloff(-ray.ray.direction)/sqr(ray.ray.t_max);}

    int Emit_Photons(RENDERING_RAY<T>& parent_ray,PHOTON_MAP<T>& photon_map,const typename PHOTON_MAP<T>::PHOTON_MAP_TYPE type,const int tid=1) const  PHYSBAM_OVERRIDE
    {int photons_emitted=0;
    while(photon_map.Light_Emission_Quota_Remains()){
        printf("photons_emitted %d\n",photons_emitted);
        printf("photons_stored %d\n",photon_map.photons_stored);
        RANDOM_NUMBERS<T>* random=0;
        if(type==PHOTON_MAP<T>::GLOBAL_PHOTON_MAP) random=&global_photon_random;
        else if(type==PHOTON_MAP<T>::CAUSTIC_PHOTON_MAP) random=&caustic_photon_random;
        else if(type==PHOTON_MAP<T>::VOLUME_PHOTON_MAP) random=&volume_photon_random;
        world.random.Set_Seed((int)abs(random->Get_Uniform_Number((T)0,(T)1000000)),tid);
        photons_emitted++;
        // Sample spotlight cone uniformly (i.e. theta=arccos(1-xi_1+xi_1*cos_cone_angle), phi=2*pi*xi_2)
        T xi_1=world.random.Get_Uniform_Number((T)0,(T)1,tid),xi_2=world.random.Get_Uniform_Number((T)0,(T)1,tid);
        T phi=2*T(pi)*xi_2,theta=acos(1-xi_1+xi_1*cos_cone_angle);
        VECTOR<T,3> local_coordinate_direction(cos(phi)*sin(theta),sin(phi)*sin(theta),cos(theta));
        VECTOR<T,3> direction=local_vector_to_world_vector*local_coordinate_direction;
        // shoot
        RENDERING_RAY<T> photon_ray(RAY<VECTOR<T,3> >(position,direction),1,parent_ray.current_object);
        world.Cast_Photon(photon_ray,parent_ray,color*brightness,type,0,0);
    }
    return photons_emitted;}

//#####################################################################
};
}
#endif
