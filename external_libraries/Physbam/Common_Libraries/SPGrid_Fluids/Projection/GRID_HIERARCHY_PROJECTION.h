//#####################################################################
// Copyright 2012-2013, Sean Bauer, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class GRID_HIERARCHY_PROJECTION
//#####################################################################
#ifndef __GRID_HIERARCHY_PROJECTION__
#define __GRID_HIERARCHY_PROJECTION__

#include <SPGrid_Fluids/Grids/GRID_HIERARCHY.h>
#include <SPGrid_Fluids/Simulation/FLUIDS_SIMULATION_DATA.h>

#include "Divergence_Helper.h"
#include "Gradient_Helper.h"

using namespace PhysBAM;

namespace SPGrid{
template<class T_STRUCT, class T,int d>
class GRID_HIERARCHY_PROJECTION
{
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<unsigned>::type Flag_array_type;
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<T>::type Data_array_type;
    typedef typename Flag_array_type::MASK Flag_array_mask;
    typedef GRID_HIERARCHY<T_STRUCT,T,d> Hierarchy_type;

    typedef std_array<T,d> Grid_Coord;
    typedef std_array<int,d> Cell_Coord;

//#define EXTERNAL_SCALING ( ( (double)hierarchy.Grid(1).one_over_dX.x ) * ( (double)hierarchy.Grid(1).one_over_dX.x ) )
#define EXTERNAL_SCALING ( (d==2) ? ( ( (double)hierarchy.Grid(1).one_over_dX.x ) * ( (double)hierarchy.Grid(1).one_over_dX.x ) ) : ( ( (double)hierarchy.Grid(1).one_over_dX.x ) * ( (double)hierarchy.Grid(1).one_over_dX.x ) * ( (double)hierarchy.Grid(1).one_over_dX.x ) ) )

public:
    // new functions
#define OLD_SCALING
    
    // -( h^(d-2) ) OR -( 1 / h^d )
    static double Laplace_Scale_Uniform(Hierarchy_type& hierarchy,const int level)
    {        
#ifdef OLD_SCALING
        return EXTERNAL_SCALING * ( (d==2) ? (double)(-1.) : (double)( (double)(-1.) * (double)(hierarchy.Grid(level).dX.x) ) );
#else
        return EXTERNAL_SCALING * ( (d==2) ? 
            (double)( (double)(-1.) * (double)(hierarchy.Grid(level).one_over_dX.x) * (double)(hierarchy.Grid(level).one_over_dX.x) ) :
            (double)( (double)(-1.) * (double)(hierarchy.Grid(level).one_over_dX.x) * (double)(hierarchy.Grid(level).one_over_dX.x) * (double)(hierarchy.Grid(level).one_over_dX.x) ) );
#endif
    }

    // 2/3 * uniform
    static double Laplace_Scale_Nonuniform(Hierarchy_type& hierarchy,const int level)
    {return (double) ( (double)(2./3.) * Laplace_Scale_Uniform(hierarchy,level) );}

    // 1/h
    static double Gradient_Scale_Uniform(Hierarchy_type& hierarchy,const int level)
    {return (double)(hierarchy.Grid(level).one_over_dX.x);}

    // 2/3 * uniform
    static double Gradient_Scale_Nonuniform(Hierarchy_type& hierarchy,const int level)
    {return (double) ( (double)(2./3.) * Gradient_Scale_Uniform(hierarchy,level) );}

    // -( h^(d-1) ) OR -( 1 / h^(d-1) )
    static double Divergence_Scale_Uniform(Hierarchy_type& hierarchy,const int level)
    {
#ifdef OLD_SCALING
        return EXTERNAL_SCALING * ( (d==2) ? (double)(-1.) * (double)(hierarchy.Grid(level).dX.x) : (double)( (double)(-1.) * (double)(hierarchy.Grid(level).dX.x) * (double)(hierarchy.Grid(level).dX.x) ) );
#else
        return EXTERNAL_SCALING * ( (d==2) ? 
            (double)( (double)(-1.) * (double)(hierarchy.Grid(level).one_over_dX.x) ):
            (double)( (double)(-1.) * (double)(hierarchy.Grid(level).one_over_dX.x) * (double)(hierarchy.Grid(level).one_over_dX.x) ) );
#endif        
    }

    static void Propagate_Ghost_Values(
        Hierarchy_type& hierarchy,
        unsigned T_STRUCT::* flags_field, 
        T T_STRUCT::* v_field); 
    
    static void Accumulate_Ghost_Values(
        Hierarchy_type& hierarchy,
        unsigned T_STRUCT::* flags_field, 
        T T_STRUCT::* v_field);

    static void Compute_Laplacian(
        Hierarchy_type& hierarchy,
        unsigned T_STRUCT::* flags_field,
        T T_STRUCT::* u_field,
        T T_STRUCT::* result_field);

    static void Compute_Variable_Beta_Laplacian( // variable beta
        Hierarchy_type& hierarchy,
        unsigned T_STRUCT::* flags_field,
        T T_STRUCT::* u_field,
        T T_STRUCT::* result_field,
        T T_STRUCT::* variable_beta_field);

    // TODO: compress into one
    
    // 2D
    static void Compute_Divergence(
        Hierarchy_type& hierarchy,
        unsigned T_STRUCT::* flags_field, 
        T T_STRUCT::* d_field,
        T T_STRUCT::* u_field,
        T T_STRUCT::* v_field);
    
    // 3D
    static void Compute_Divergence(
        Hierarchy_type& hierarchy,
        unsigned T_STRUCT::* flags_field, 
        T T_STRUCT::* d_field,
        T T_STRUCT::* u_field,
        T T_STRUCT::* v_field,
        T T_STRUCT::* w_field);

    // 2D
    static void Compute_Gradient(
        Hierarchy_type& hierarchy,
        unsigned T_STRUCT::* flags_field, 
        T T_STRUCT::* gu_field,
        T T_STRUCT::* gv_field,
        T T_STRUCT::* d_field);
    
    // 3D
    static void Compute_Gradient(
        Hierarchy_type& hierarchy,
        unsigned T_STRUCT::* flags_field, 
        T T_STRUCT::* gu_field,
        T T_STRUCT::* gv_field,
        T T_STRUCT::* gw_field,
        T T_STRUCT::* d_field);

//#####################################################################
};
}
#endif
