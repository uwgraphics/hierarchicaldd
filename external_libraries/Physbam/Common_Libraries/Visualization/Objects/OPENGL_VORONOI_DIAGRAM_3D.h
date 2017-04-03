//#####################################################################
// Copyright 2015, Mridul Aanjaneya.
// This file is part of PhysBAM whose distribution is governed by the license 
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_VORONOI_DIAGRAM_3D
//#####################################################################
#ifndef __OPENGL_VORONOI_DIAGRAM_3D__
#define __OPENGL_VORONOI_DIAGRAM_3D__

#include <Common_Geometry/Topology_Based_Geometry/VORONOI_DIAGRAM.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_COLOR.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_OBJECT.h>

namespace PhysBAM{
template<class T>
class OPENGL_VORONOI_DIAGRAM_3D: public OPENGL_OBJECT
{
  public:
    typedef VECTOR<T,3> TV;

    VORONOI_DIAGRAM<TV> voronoi_diagram;
    OPENGL_COLOR color;
    bool draw;
    int current_frame;
    RANGE<VECTOR<float,3> > domain;

    OPENGL_VORONOI_DIAGRAM_3D(const OPENGL_COLOR& color_input=OPENGL_COLOR::Cyan())
        :color(color_input),draw(true),current_frame(-1)
    {}

//#####################################################################
    virtual RANGE<VECTOR<float,3> > Bounding_Box() const PHYSBAM_OVERRIDE;
    void Display(const int in_color=1) const PHYSBAM_OVERRIDE;
    void Update_Domain();
//#####################################################################
};
}
#endif
