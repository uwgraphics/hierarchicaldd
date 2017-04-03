//#####################################################################
// Copyright (c) 2013, Eftychios Sifakis, Raj Setaluri.
// Distributed under the FreeBSD license (see license.txt)
//#####################################################################
#ifndef __GRID_TOPOLOGY_HELPER_h__
#define __GRID_TOPOLOGY_HELPER_h__

#include <PhysBAM_Tools/Log/DEBUG_UTILITIES.h>
#include <Common_Tools/Math_Tools/RANGE_ITERATOR.h>
#include <SPGrid/Core/SPGrid_Mask.h>
#include <SPGrid_Fluids/Simulation/FLUIDS_SIMULATION_FLAGS.h>

namespace PhysBAM{
using namespace SPGrid;

// Axis_Vector_Offset_Helper
namespace{
template<class T_MASK> struct Axis_Vector_Offset_Helper;        
template<int log2_struct,int log2_field>
struct Axis_Vector_Offset_Helper<SPGrid_Mask<log2_struct,log2_field,2> >
{
    enum{
        x_mask=SPGrid_Mask<log2_struct,log2_field,2>::template LinearOffset<1,0>::value,
        y_mask=SPGrid_Mask<log2_struct,log2_field,2>::template LinearOffset<0,1>::value,
        z_mask=0UL
    };
};    
template<int log2_struct,int log2_field>
struct Axis_Vector_Offset_Helper<SPGrid_Mask<log2_struct,log2_field,3> >
{
    enum{
        x_mask=SPGrid_Mask<log2_struct,log2_field,3>::template LinearOffset<1,0,0>::value,
        y_mask=SPGrid_Mask<log2_struct,log2_field,3>::template LinearOffset<0,1,0>::value,
        z_mask=SPGrid_Mask<log2_struct,log2_field,3>::template LinearOffset<0,0,1>::value
    };
};
}
// Negative_Axis_Vector_Offset_Helper
namespace{
template<class T_MASK> struct Negative_Axis_Vector_Offset_Helper;
template<int log2_struct,int log2_field>
struct Negative_Axis_Vector_Offset_Helper<SPGrid_Mask<log2_struct,log2_field,2> >
{
    enum{
        x_mask=SPGrid_Mask<log2_struct,log2_field,2>::template LinearOffset<-1,0>::value,
        y_mask=SPGrid_Mask<log2_struct,log2_field,2>::template LinearOffset<0,-1>::value,
        z_mask=0UL
    };
};    
template<int log2_struct,int log2_field>
struct Negative_Axis_Vector_Offset_Helper<SPGrid_Mask<log2_struct,log2_field,3> >
{
    enum{
        x_mask=SPGrid_Mask<log2_struct,log2_field,3>::template LinearOffset<-1,0,0>::value,
        y_mask=SPGrid_Mask<log2_struct,log2_field,3>::template LinearOffset<0,-1,0>::value,
        z_mask=SPGrid_Mask<log2_struct,log2_field,3>::template LinearOffset<0,0,-1>::value
    };
};
}
// Shadow Grid Helper
namespace{
template<class T_MASK> struct Shadow_Grid_Stencil_Offsets_Helper;        
template<int log2_struct,int log2_field>
struct Shadow_Grid_Stencil_Offsets_Helper<SPGrid_Mask<log2_struct,log2_field,2> >
{
    enum{d=2};
    typedef SPGrid_Mask<log2_struct,log2_field,d> T_MASK;
    enum{og_xsize = (1u << T_MASK::block_xbits)+2,og_ysize = (1u << T_MASK::block_ybits)+2};
    enum{x_shift=og_ysize,y_shift=1};
    static int Offset(std_array<int,d> vector){return x_shift*vector(0)+y_shift*vector(1);}
};
template<int log2_struct,int log2_field>
struct Shadow_Grid_Stencil_Offsets_Helper<SPGrid_Mask<log2_struct,log2_field,3> >
{
    enum{d=3};
    typedef SPGrid_Mask<log2_struct,log2_field,d> T_MASK;
    enum{og_xsize = (1u << T_MASK::block_xbits)+2,og_ysize = (1u << T_MASK::block_ybits)+2,og_zsize = (1u << T_MASK::block_zbits)+2};
    enum{x_shift=og_ysize*og_zsize,y_shift=og_zsize,z_shift=1};
    static int Offset(std_array<int,d> vector){return x_shift*vector(0)+y_shift*vector(1)+z_shift*vector(2);}
};    
}
//#####################################################################
// Class GRID_TOPOLOGY_HELPER
//#####################################################################
template<class T_MASK>
struct GRID_TOPOLOGY_HELPER
{
    enum {d=T_MASK::dim};
    STATIC_ASSERT(d==2 || d==3);

    typedef VECTOR<int,d> T_INDEX;
    enum {nodes_per_face=(d==2)?2:4};
    enum {nodes_per_cell=(d==2)?4:8};
    enum {faces_per_cell=(d==2)?4:6};

    static void Nodes_Of_Face_Offsets(unsigned long nodes_of_face_offsets[nodes_per_face],int axis)
    {
        PHYSBAM_ASSERT(1<=axis && axis<=d);
        int node=0;
        for(RANGE_ITERATOR<d> iterator(RANGE<T_INDEX>(T_INDEX(),T_INDEX::All_Ones_Vector()-T_INDEX::Axis_Vector(axis)));iterator.Valid();iterator.Next())
            nodes_of_face_offsets[node++]=T_MASK::Linear_Offset(std_array<int,d>(iterator.Index()));
    }

    static void Nodes_Of_Cell_Offsets(unsigned long nodes_of_cell_offsets[nodes_per_cell])
    {
        int node=0;
        for(RANGE_ITERATOR<d> iterator(RANGE<T_INDEX>(T_INDEX(),T_INDEX::All_Ones_Vector()));iterator.Valid();iterator.Next())
            nodes_of_cell_offsets[node++]=T_MASK::Linear_Offset(std_array<int,d>(iterator.Index()));
    }

    static void Nodes_Of_Face_Shadow_Grid_Offsets(int nodes_of_face_offsets[nodes_per_face],int axis)
    {
        PHYSBAM_ASSERT(1<=axis && axis<=d);
        int node=0;
        for(RANGE_ITERATOR<d> iterator(RANGE<T_INDEX>(T_INDEX(),T_INDEX::All_Ones_Vector()-T_INDEX::Axis_Vector(axis)));iterator.Valid();iterator.Next())
            nodes_of_face_offsets[node++]=Shadow_Grid_Stencil_Offsets_Helper<T_MASK>::Offset(std_array<int,d>(iterator.Index()));
    }

    static void Nodes_Of_Cell_Shadow_Grid_Offsets(int nodes_of_cell_offsets[nodes_per_cell])
    {
        int node=0;
        for(RANGE_ITERATOR<d> iterator(RANGE<T_INDEX>(T_INDEX(),T_INDEX::All_Ones_Vector()));iterator.Valid();iterator.Next())
            nodes_of_cell_offsets[node++]=Shadow_Grid_Stencil_Offsets_Helper<T_MASK>::Offset(std_array<int,d>(iterator.Index()));
    }

    static void Face_Neighbor_Offsets(unsigned long face_neighbor_offsets[faces_per_cell])
    {
        int face=0;
        for(int axis=0;axis<d;axis++)
            for(int side=-1;side<=1;side+=2){
                std_array<int,d> shift;
                shift(axis)=side;
                face_neighbor_offsets[face]=T_MASK::Linear_Offset(shift);
                face++;}
    }

    static unsigned Face_Valid_Mask(int axis)
    {
        PHYSBAM_ASSERT(1<=axis && axis<=d);
        if(axis==1) return SPGrid_Face_Type_X_Valid;
        if(axis==2) return SPGrid_Face_Type_Y_Valid;
        return SPGrid_Face_Type_Z_Valid;
    }

    static unsigned Face_Active_Mask(int axis)
    {
        PHYSBAM_ASSERT(1<=axis && axis<=d);
        if(axis==1) return SPGrid_Face_Type_X_Active;
        if(axis==2) return SPGrid_Face_Type_Y_Active;
        return SPGrid_Face_Type_Z_Active;
    }

    static unsigned Face_Plus_Minus_Active_Mask(int axis,int side)
    {
        PHYSBAM_ASSERT(1<=axis && axis<=d && 1<=side && side<=2);
        if(side==1){
            if(axis==1) return SPGrid_Face_Minus_X_Active;
            if(axis==2) return SPGrid_Face_Minus_Y_Active;
            return SPGrid_Face_Minus_Z_Active;}
        else{
            if(axis==1) return SPGrid_Face_Plus_X_Active;
            if(axis==2) return SPGrid_Face_Plus_Y_Active;
            return SPGrid_Face_Plus_Z_Active;}
    }

    static unsigned Face_Minus_Scaled_Mask(int axis)
    {
        PHYSBAM_ASSERT(1<=axis && axis<=d);
        if(axis==1) return SPGrid_Face_Minus_X_Scaled;
        if(axis==2) return SPGrid_Face_Minus_Y_Scaled;
        return SPGrid_Face_Minus_Z_Scaled;
    }

    static unsigned Face_Plus_Minus_Scaled_Mask(int axis,int side)
    {
        PHYSBAM_ASSERT(1<=axis && axis<=d && 1<=side && side<=2);
        if(side==1){
            if(axis==1) return SPGrid_Face_Minus_X_Scaled;
            if(axis==2) return SPGrid_Face_Minus_Y_Scaled;
            return SPGrid_Face_Minus_Z_Scaled;}
        else{
            if(axis==1) return SPGrid_Face_Plus_X_Scaled;
            if(axis==2) return SPGrid_Face_Plus_Y_Scaled;
            return SPGrid_Face_Plus_Z_Scaled;}
    }

    static unsigned long Axis_Vector_Offset(int axis)
    {
        PHYSBAM_ASSERT(1<=axis && axis<=d);
        if(axis==1) return Axis_Vector_Offset_Helper<T_MASK>::x_mask;
        if(axis==2) return Axis_Vector_Offset_Helper<T_MASK>::y_mask;
        return Axis_Vector_Offset_Helper<T_MASK>::z_mask;
    }
    
    static unsigned long Negative_Axis_Vector_Offset(int axis)
    {
        PHYSBAM_ASSERT(1<=axis && axis<=d);
        if(axis==1) return Negative_Axis_Vector_Offset_Helper<T_MASK>::x_mask;
        if(axis==2) return Negative_Axis_Vector_Offset_Helper<T_MASK>::y_mask;
        return Negative_Axis_Vector_Offset_Helper<T_MASK>::z_mask;
    }

};

}
#endif
