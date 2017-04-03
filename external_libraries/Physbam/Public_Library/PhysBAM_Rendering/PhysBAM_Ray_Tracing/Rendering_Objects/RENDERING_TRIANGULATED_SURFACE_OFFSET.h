//#####################################################################
// Copyright 2013, Yue Yu.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RENDERING_TRIANGULATED_SURFACE_OFFSET
//#####################################################################
#ifndef __RENDERING_TRIANGULATED_SURFACE_OFFSET__
#define __RENDERING_TRIANGULATED_SURFACE_OFFSET__


#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering_Objects/RENDERING_TRIANGULATED_SURFACE.h>
namespace PhysBAM{

template<class T>
class RENDERING_TRIANGULATED_SURFACE_OFFSET:public RENDERING_TRIANGULATED_SURFACE<T>
{
    typedef VECTOR<T,3> TV;
public:
    using RENDERING_TRIANGULATED_SURFACE<T>::triangulated_surface;

    const ARRAY_VIEW<T>& height;
    const ARRAY_VIEW<TV>& normal;
    T minimum_renderable_height;
    T thickness_scale;

    RENDERING_TRIANGULATED_SURFACE_OFFSET(const ARRAY_VIEW<T>& height_input,const ARRAY_VIEW<TV>& normal_input,T minimum_renderable_height_input,T thickness_scale_input,TRIANGULATED_SURFACE<T>& triangulated_surface_input,const int triangles_per_hierarchy_group=0)
        :RENDERING_TRIANGULATED_SURFACE<T>(triangulated_surface_input,triangles_per_hierarchy_group),
        height(height_input),normal(normal_input),minimum_renderable_height(minimum_renderable_height_input),thickness_scale(thickness_scale_input)
    {
        if(triangulated_surface.particles.X.m!=height.m||height.m!=normal.m){
            std::stringstream ss;ss<<"triangulated_surface and height/normal do not have matching #particles! triangulated_surface.particles.X.m="<<triangulated_surface.particles.X.m<<" height.m="<<height.m<<" normal.m="<<normal.m;
            PHYSBAM_FATAL_ERROR(ss.str());}
        for(int j=1;j<=triangulated_surface.particles.X.m;++j){
            triangulated_surface.particles.X(j)+=normal(j)*(height(j)>minimum_renderable_height?thickness_scale*height(j):-minimum_renderable_height);
        }
        triangulated_surface.Update_Bounding_Box();triangulated_surface.Initialize_Hierarchy(true,triangles_per_hierarchy_group);
        if(triangulated_surface.use_vertex_normals)triangulated_surface.Update_Vertex_Normals();
    }

    virtual ~RENDERING_TRIANGULATED_SURFACE_OFFSET()
    {}

//#####################################################################
};   
}
#endif
