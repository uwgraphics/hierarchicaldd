//#####################################################################
// Copyright 2015, Mridul Aanjaneya.
// This file is part of PhysBAM whose distribution is governed by the license 
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_COLOR_RAMP.h>
#include <Visualization/Objects/OPENGL_QUADTREE_SCALAR_FIELD.h>
using namespace PhysBAM;
//#####################################################################
// Function Constructor
//#####################################################################
template<class T,class T2> OPENGL_QUADTREE_SCALAR_FIELD<T,T2>::
OPENGL_QUADTREE_SCALAR_FIELD(OPENGL_QUADTREE<T>& opengl_quadtree_input,OPENGL_COLOR_MAP<T2>* color_map_input)
    :opengl_quadtree(opengl_quadtree_input),current_color_map(1),scale_range(false)
{
    PHYSBAM_ASSERT(color_map_input);
    Initialize_Color_Maps(color_map_input);
}
//#####################################################################
// Function Destructor
//#####################################################################
template<class T,class T2> OPENGL_QUADTREE_SCALAR_FIELD<T,T2>::
~OPENGL_QUADTREE_SCALAR_FIELD()
{
    color_maps.Delete_Pointers_And_Clean_Memory();
}
//#####################################################################
// Function Initialize_Color_Maps
//#####################################################################
template<> void OPENGL_QUADTREE_SCALAR_FIELD<float,bool>::
Initialize_Color_Maps(OPENGL_COLOR_MAP<bool>* color_map_input)
{
    color_maps.Append(color_map_input);
}
//#####################################################################
// Function Initialize_Color_Maps
//#####################################################################
template<> void OPENGL_QUADTREE_SCALAR_FIELD<double,bool>::
Initialize_Color_Maps(OPENGL_COLOR_MAP<bool>* color_map_input)
{
    color_maps.Append(color_map_input);
}
//#####################################################################
// Function Initialize_Color_Maps
//#####################################################################
template<> void OPENGL_QUADTREE_SCALAR_FIELD<float,int>::
Initialize_Color_Maps(OPENGL_COLOR_MAP<int>* color_map_input)
{
    color_maps.Append(color_map_input);
}
//#####################################################################
// Function Initialize_Color_Maps
//#####################################################################
template<> void OPENGL_QUADTREE_SCALAR_FIELD<double,int>::
Initialize_Color_Maps(OPENGL_COLOR_MAP<int>* color_map_input)
{
    color_maps.Append(color_map_input);
}
//#####################################################################
// Function Initialize_Color_Maps
//#####################################################################
template<class T,class T2> void OPENGL_QUADTREE_SCALAR_FIELD<T,T2>::
Initialize_Color_Maps(OPENGL_COLOR_MAP<T2>* color_map_input)
{
    color_maps.Append(color_map_input);
    color_maps.Append(OPENGL_COLOR_RAMP<T2>::Matlab_Jet(0,1));
    color_maps.Append(OPENGL_COLOR_RAMP<T2>::Matlab_Hot(0,1));
}
//#####################################################################
// Function Set_Scale_Range
//#####################################################################
template<> void OPENGL_QUADTREE_SCALAR_FIELD<float,bool>::
Set_Scale_Range(const bool range_min,const bool range_max)
{PHYSBAM_FATAL_ERROR();}
//#####################################################################
// Function Set_Scale_Range
//#####################################################################
template<> void OPENGL_QUADTREE_SCALAR_FIELD<double,bool>::
Set_Scale_Range(const bool range_min,const bool range_max)
{PHYSBAM_FATAL_ERROR();}
//#####################################################################
// Function Set_Scale_Range
//#####################################################################
template<> void OPENGL_QUADTREE_SCALAR_FIELD<float,int>::
Set_Scale_Range(const int range_min,const int range_max)
{PHYSBAM_FATAL_ERROR();}
//#####################################################################
// Function Set_Scale_Range
//#####################################################################
template<> void OPENGL_QUADTREE_SCALAR_FIELD<double,int>::
Set_Scale_Range(const int range_min,const int range_max)
{PHYSBAM_FATAL_ERROR();}
//#####################################################################
// Function Set_Scale_Range
//#####################################################################
template<class T,class T2> void OPENGL_QUADTREE_SCALAR_FIELD<T,T2>::
Set_Scale_Range(const T2 range_min,const T2 range_max)
{
    scale_range=true;
    scale_range_min=range_min;
    T2 range_length=(range_max-range_min);
    scale_range_dx=range_length>1e-10?(T2)1/range_length:(T2)0;
}
//#####################################################################
// Function Reset_Scale_Range
//#####################################################################
template<class T,class T2> void OPENGL_QUADTREE_SCALAR_FIELD<T,T2>::
Reset_Scale_Range()
{
    scale_range=false;
}
//#####################################################################
// Function Pre_Map_Value
//#####################################################################
template<> bool OPENGL_QUADTREE_SCALAR_FIELD<float,bool>::
Pre_Map_Value(const bool value) const
{
    return value;
}
//#####################################################################
// Function Pre_Map_Value
//#####################################################################
template<> bool OPENGL_QUADTREE_SCALAR_FIELD<double,bool>::
Pre_Map_Value(const bool value) const
{
    return value;
}
//#####################################################################
// Function Pre_Map_Value
//#####################################################################
template<class T,class T2> T2 OPENGL_QUADTREE_SCALAR_FIELD<T,T2>::
Pre_Map_Value(const T2 value) const
{
    if(!scale_range) return value;
    else return (value-scale_range_min)*scale_range_dx; 
}
//#####################################################################
// Function Display
//#####################################################################
template<class T,class T2> void OPENGL_QUADTREE_SCALAR_FIELD<T,T2>::
Display(const int in_color) const
{ 
    glPushAttrib(GL_ENABLE_BIT|GL_DEPTH_BUFFER_BIT|GL_CURRENT_BIT);
    glDisable(GL_LIGHTING);
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_BLEND);
    glDepthMask(GL_FALSE);

    OpenGL_Begin(GL_QUADS);
    for(int i=1;i<=opengl_quadtree.tree.nodes.Size();++i){QUADTREE_NODE<T> *node=opengl_quadtree.tree.nodes(i);
        if(node->north_east_index==0){TV min_corner=node->center-node->dx_over_two,max_corner=node->center+node->dx_over_two;
            color_maps(current_color_map)->Lookup(data.Get(i)).Send_To_GL_Pipeline();
            OpenGL_Vertex(min_corner);
            OpenGL_Vertex(min_corner+TV::Axis_Vector(1)*(T)2.*node->dx_over_two(1));
            OpenGL_Vertex(max_corner);
            OpenGL_Vertex(min_corner+TV::Axis_Vector(2)*(T)2.*node->dx_over_two(2));}}
    OpenGL_End();
    glPopAttrib();
}
//#####################################################################
// Function Bounding_Box
//#####################################################################
template<class T,class T2> RANGE<VECTOR<float,3> > OPENGL_QUADTREE_SCALAR_FIELD<T,T2>::
Bounding_Box() const
{
    // May not be the exact bounds, but close enough...
    return World_Space_Box(RANGE<VECTOR<float,3> >(opengl_quadtree.domain.min_corner.x,opengl_quadtree.domain.max_corner.x,opengl_quadtree.domain.min_corner.y,opengl_quadtree.domain.max_corner.y,0,0));
}
//#####################################################################
// Function Print_Selection_Info
//#####################################################################
template<class T,class T2> void OPENGL_QUADTREE_SCALAR_FIELD<T,T2>::
Print_Selection_Info(std::ostream& output_stream,OPENGL_SELECTION* current_selection) const
{
    if(current_selection && current_selection->type==OPENGL_SELECTION::GRID_CELL_2D){
        int index=((OPENGL_SELECTION_QUADTREE_CELL<T>*)current_selection)->index;
        output_stream<<data.Get(index)<<std::endl;}
}
//#####################################################################
// Function Toggle_Color_Map
//#####################################################################
template<class T,class T2> void OPENGL_QUADTREE_SCALAR_FIELD<T,T2>::
Toggle_Color_Map()
{
    current_color_map=current_color_map%color_maps.m+1;
}
//#####################################################################
template class OPENGL_QUADTREE_SCALAR_FIELD<float,int>;
template class OPENGL_QUADTREE_SCALAR_FIELD<float,bool>;
template class OPENGL_QUADTREE_SCALAR_FIELD<float,float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class OPENGL_QUADTREE_SCALAR_FIELD<double,int>;
template class OPENGL_QUADTREE_SCALAR_FIELD<double,bool>;
template class OPENGL_QUADTREE_SCALAR_FIELD<double,float>;
#endif
