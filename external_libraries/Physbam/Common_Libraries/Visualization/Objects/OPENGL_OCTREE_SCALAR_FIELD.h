//#####################################################################
// Copyright 2015, Mridul Aanjaneya.
// This file is part of PhysBAM whose distribution is governed by the license 
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_OCTREE_SCALAR_FIELD
//##################################################################### 
#ifndef __OPENGL_OCTREE_SCALAR_FIELD__
#define __OPENGL_OCTREE_SCALAR_FIELD__

#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_CALLBACK.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_COLOR_MAP.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_OBJECT.h>
#include <PhysBAM_Tools/Data_Structures/HASHTABLE.h>
#include <Visualization/Objects/OPENGL_OCTREE.h>

namespace PhysBAM{
template<class T,class T2=T>
class OPENGL_OCTREE_SCALAR_FIELD: public OPENGL_OBJECT
{
    typedef VECTOR<T,3> TV;typedef VECTOR<int,3> T_INDEX;
public:
    HASHTABLE<int,T> data;
    OPENGL_OCTREE<T>& opengl_octree;

    ARRAY<OPENGL_COLOR_MAP<T2>*> color_maps; // all owned by us
    int current_color_map;
private:
    bool scale_range;
    T2 scale_range_min,scale_range_dx;

public:
    OPENGL_OCTREE_SCALAR_FIELD(OPENGL_OCTREE<T>& opengl_octree_input,OPENGL_COLOR_MAP<T2>* color_map_input);
    ~OPENGL_OCTREE_SCALAR_FIELD();

    void Set_Scale_Range(const T2 range_min,const T2 range_max);
    void Reset_Scale_Range();
    T2 Pre_Map_Value(const T2 value) const;

    void Display(const int in_color=1) const PHYSBAM_OVERRIDE;
    virtual RANGE<VECTOR<float,3> > Bounding_Box() const PHYSBAM_OVERRIDE;

    void Print_Selection_Info(std::ostream& output_stream,OPENGL_SELECTION* current_selection) const PHYSBAM_OVERRIDE;
    void Toggle_Color_Map();

private:
    void Initialize_Color_Maps(OPENGL_COLOR_MAP<T2>* color_map_input);
};
}
#endif
