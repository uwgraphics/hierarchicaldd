//#####################################################################
// Copyright 2013, Sean Bauer, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include "Variable_Beta_Laplace_Helper.h"
#include <SPGrid_Fluids/Simulation/FLUIDS_SIMULATION_FLAGS.h>
#include <SPGrid_Fluids/Simulation/FLUIDS_SIMULATION_DATA.h>

using namespace SPGrid;
using namespace PhysBAM;

//#####################################################################
// Function Run_Index_Range
//#####################################################################
// T_MASK corresponds to the mask for the data (not the mask channel)
template <class T, int log2_struct> void Variable_Beta_Laplace_Helper<T,log2_struct,2>::Run_Index_Range(const int index_start,const int index_end)
{  
    // Compute shadow grid of linear offsets
    unsigned long* offset_grid_ptr = (unsigned long*)malloc( (og_xsize) * (og_ysize) * sizeof(unsigned long));
    typedef unsigned long (&offset_grid_type)[og_xsize][og_ysize];
    offset_grid_type o_grid = reinterpret_cast<offset_grid_type>(*offset_grid_ptr);
    
    // Iterating through all block indices
    for(int index=index_start;index<=index_end;index++) {
      T* output = reinterpret_cast<T*>((unsigned long)x + b[index]);
      T* y1_in  = reinterpret_cast<T*>((unsigned long)y + b[index]);
      T* rho_in  = reinterpret_cast<T*>((unsigned long)rho + b[index]);

      unsigned* mask_in  = reinterpret_cast<unsigned*>((unsigned long)mask + b[index]);
      unsigned long packed_offset = (unsigned long)y1_in - (unsigned long)y;
      unsigned long y_base_addr = reinterpret_cast<unsigned long>(y);
      unsigned long rho_base_addr = reinterpret_cast<unsigned long>(rho);

      ComputeShadowGrid(offset_grid_ptr, mask_in, packed_offset);

      int cur_index = 0;
      // Actually process elements
      for(int i=xmin;i<=xmax;i++)
      for(int j=ymin;j<=ymax;j++)
      {
        unsigned mask_value = mask_in[cur_index];

        if ( mask_value & (SPGrid_Cell_Type_Active|SPGrid_Cell_Type_Ghost)) {
            double result=(double)0.;

            if (mask_value & SPGrid_Face_Minus_X_Scaled){
                const double neighbor_rho=(double)(*reinterpret_cast<T*>(rho_base_addr + o_grid[i-1][j]));
                const double beta=(double)2./(rho_in[cur_index]+neighbor_rho);
                result += beta * scale_nonuniform * ((double)(*reinterpret_cast<T*>(y_base_addr + o_grid[i-1][j])) - y1_in[cur_index]);}
            else if (mask_value & SPGrid_Face_Minus_X_Active){
                const double neighbor_rho=(double)(*reinterpret_cast<T*>(rho_base_addr + o_grid[i-1][j]));
                const double beta=(double)2./(rho_in[cur_index]+neighbor_rho);
                result += beta * scale_uniform    * ((double)(*reinterpret_cast<T*>(y_base_addr + o_grid[i-1][j])) - y1_in[cur_index]);}

            if (mask_value & SPGrid_Face_Plus_X_Scaled){
                const double neighbor_rho=(double)(*reinterpret_cast<T*>(rho_base_addr + o_grid[i+1][j]));
                const double beta=(double)2./(rho_in[cur_index]+neighbor_rho);
                result += beta * scale_nonuniform * ((double)(*reinterpret_cast<T*>(y_base_addr + o_grid[i+1][j])) - y1_in[cur_index]);}
            else if (mask_value & SPGrid_Face_Plus_X_Active){
                const double neighbor_rho=(double)(*reinterpret_cast<T*>(rho_base_addr + o_grid[i+1][j]));
                const double beta=(double)2./(rho_in[cur_index]+neighbor_rho);
                result += beta * scale_uniform    * ((double)(*reinterpret_cast<T*>(y_base_addr + o_grid[i+1][j])) - y1_in[cur_index]);}

            
            if (mask_value & SPGrid_Face_Minus_Y_Scaled){
                const double neighbor_rho=(double)(*reinterpret_cast<T*>(rho_base_addr + o_grid[i][j-1]));
                const double beta=(double)2./(rho_in[cur_index]+neighbor_rho);
                result += beta * scale_nonuniform * ((double)(*reinterpret_cast<T*>(y_base_addr + o_grid[i][j-1])) - y1_in[cur_index]);}
            else if (mask_value & SPGrid_Face_Minus_Y_Active){                                
                const double neighbor_rho=(double)(*reinterpret_cast<T*>(rho_base_addr + o_grid[i][j-1]));
                const double beta=(double)2./(rho_in[cur_index]+neighbor_rho);
                result += beta * scale_uniform    * ((double)(*reinterpret_cast<T*>(y_base_addr + o_grid[i][j-1])) - y1_in[cur_index]);}

            if (mask_value & SPGrid_Face_Plus_Y_Scaled){                                      
                const double neighbor_rho=(double)(*reinterpret_cast<T*>(rho_base_addr + o_grid[i][j+1]));
                const double beta=(double)2./(rho_in[cur_index]+neighbor_rho);
                result += beta * scale_nonuniform * ((double)(*reinterpret_cast<T*>(y_base_addr + o_grid[i][j+1])) - y1_in[cur_index]);}
            else if (mask_value & SPGrid_Face_Plus_Y_Active){                                 
                const double neighbor_rho=(double)(*reinterpret_cast<T*>(rho_base_addr + o_grid[i][j+1]));
                const double beta=(double)2./(rho_in[cur_index]+neighbor_rho);
                result += beta * scale_uniform    * ((double)(*reinterpret_cast<T*>(y_base_addr + o_grid[i][j+1])) - y1_in[cur_index]);}
            
            output[cur_index] = (T)result;
        }
        cur_index++;
      }
    }
    free(offset_grid_ptr);
}
//#####################################################################
// Function Run_Index_Range
//#####################################################################
// T_MASK corresponds to the mask for the data (not the mask channel)
template <class T, int log2_struct> void Variable_Beta_Laplace_Helper<T,log2_struct,3>::Run_Index_Range(const int index_start,const int index_end)
{  
    // Compute shadow grid of linear offsets
    unsigned long* offset_grid_ptr = (unsigned long*)malloc( (og_xsize) * (og_ysize) * (og_zsize) * sizeof(unsigned long));
    typedef unsigned long (&offset_grid_type)[og_xsize][og_ysize][og_zsize];
    offset_grid_type o_grid = reinterpret_cast<offset_grid_type>(*offset_grid_ptr);
    
    // Iterating through all block indices
    for(int index=index_start;index<=index_end;index++) {
      T* output = reinterpret_cast<T*>((unsigned long)x + b[index]);
      T* y1_in  = reinterpret_cast<T*>((unsigned long)y + b[index]);
      T* rho_in  = reinterpret_cast<T*>((unsigned long)rho + b[index]);

      unsigned* mask_in  = reinterpret_cast<unsigned*>((unsigned long)mask + b[index]);
      unsigned long packed_offset = (unsigned long)y1_in - (unsigned long)y;
      unsigned long y_base_addr = reinterpret_cast<unsigned long>(y);
      unsigned long rho_base_addr = reinterpret_cast<unsigned long>(rho);

      ComputeShadowGrid(offset_grid_ptr, mask_in, packed_offset);

      int cur_index = 0;
      // Actually process elements
      for(int i=xmin;i<=xmax;i++)
      for(int j=ymin;j<=ymax;j++)
      for(int k=zmin;k<=zmax;k++)
      {
        unsigned mask_value = mask_in[cur_index];

        if ( mask_value & (SPGrid_Cell_Type_Active|SPGrid_Cell_Type_Ghost)) {
            double result=(double)0.;

            if (mask_value & SPGrid_Face_Minus_X_Scaled){
                const double neighbor_rho=(double)(*reinterpret_cast<T*>(rho_base_addr + o_grid[i-1][j][k]));
                const double beta=(double)2./(rho_in[cur_index]+neighbor_rho);
                result += beta * scale_nonuniform * ((double)(*reinterpret_cast<T*>(y_base_addr + o_grid[i-1][j][k])) - y1_in[cur_index]);}
            else if (mask_value & SPGrid_Face_Minus_X_Active){
                const double neighbor_rho=(double)(*reinterpret_cast<T*>(rho_base_addr + o_grid[i-1][j][k]));
                const double beta=(double)2./(rho_in[cur_index]+neighbor_rho);
                result += beta * scale_uniform    * ((double)(*reinterpret_cast<T*>(y_base_addr + o_grid[i-1][j][k])) - y1_in[cur_index]);}

            if (mask_value & SPGrid_Face_Plus_X_Scaled){
                const double neighbor_rho=(double)(*reinterpret_cast<T*>(rho_base_addr + o_grid[i+1][j][k]));
                const double beta=(double)2./(rho_in[cur_index]+neighbor_rho);
                result += beta * scale_nonuniform * ((double)(*reinterpret_cast<T*>(y_base_addr + o_grid[i+1][j][k])) - y1_in[cur_index]);}
            else if (mask_value & SPGrid_Face_Plus_X_Active){
                const double neighbor_rho=(double)(*reinterpret_cast<T*>(rho_base_addr + o_grid[i+1][j][k]));
                const double beta=(double)2./(rho_in[cur_index]+neighbor_rho);
                result += beta * scale_uniform    * ((double)(*reinterpret_cast<T*>(y_base_addr + o_grid[i+1][j][k])) - y1_in[cur_index]);}

            
            if (mask_value & SPGrid_Face_Minus_Y_Scaled){
                const double neighbor_rho=(double)(*reinterpret_cast<T*>(rho_base_addr + o_grid[i][j-1][k]));
                const double beta=(double)2./(rho_in[cur_index]+neighbor_rho);
                result += beta * scale_nonuniform * ((double)(*reinterpret_cast<T*>(y_base_addr + o_grid[i][j-1][k])) - y1_in[cur_index]);}
            else if (mask_value & SPGrid_Face_Minus_Y_Active){                                
                const double neighbor_rho=(double)(*reinterpret_cast<T*>(rho_base_addr + o_grid[i][j-1][k]));
                const double beta=(double)2./(rho_in[cur_index]+neighbor_rho);
                result += beta * scale_uniform    * ((double)(*reinterpret_cast<T*>(y_base_addr + o_grid[i][j-1][k])) - y1_in[cur_index]);}

            if (mask_value & SPGrid_Face_Plus_Y_Scaled){                                      
                const double neighbor_rho=(double)(*reinterpret_cast<T*>(rho_base_addr + o_grid[i][j+1][k]));
                const double beta=(double)2./(rho_in[cur_index]+neighbor_rho);
                result += beta * scale_nonuniform * ((double)(*reinterpret_cast<T*>(y_base_addr + o_grid[i][j+1][k])) - y1_in[cur_index]);}
            else if (mask_value & SPGrid_Face_Plus_Y_Active){                                 
                const double neighbor_rho=(double)(*reinterpret_cast<T*>(rho_base_addr + o_grid[i][j+1][k]));
                const double beta=(double)2./(rho_in[cur_index]+neighbor_rho);
                result += beta * scale_uniform    * ((double)(*reinterpret_cast<T*>(y_base_addr + o_grid[i][j+1][k])) - y1_in[cur_index]);}
            
            
            if (mask_value & SPGrid_Face_Minus_Z_Scaled){
                const double neighbor_rho=(double)(*reinterpret_cast<T*>(rho_base_addr + o_grid[i][j][k-1]));
                const double beta=(double)2./(rho_in[cur_index]+neighbor_rho);
                result += beta * scale_nonuniform * ((double)(*reinterpret_cast<T*>(y_base_addr + o_grid[i][j][k-1])) - y1_in[cur_index]);}
            else if (mask_value & SPGrid_Face_Minus_Z_Active){                                   
                const double neighbor_rho=(double)(*reinterpret_cast<T*>(rho_base_addr + o_grid[i][j][k-1]));
                const double beta=(double)2./(rho_in[cur_index]+neighbor_rho);
                result += beta * scale_uniform    * ((double)(*reinterpret_cast<T*>(y_base_addr + o_grid[i][j][k-1])) - y1_in[cur_index]);}
            
            if (mask_value & SPGrid_Face_Plus_Z_Scaled){                                         
                const double neighbor_rho=(double)(*reinterpret_cast<T*>(rho_base_addr + o_grid[i][j][k+1]));
                const double beta=(double)2./(rho_in[cur_index]+neighbor_rho);
                result += beta * scale_nonuniform * ((double)(*reinterpret_cast<T*>(y_base_addr + o_grid[i][j][k+1])) - y1_in[cur_index]);}
            else if (mask_value & SPGrid_Face_Plus_Z_Active){                                    
                const double neighbor_rho=(double)(*reinterpret_cast<T*>(rho_base_addr + o_grid[i][j][k+1]));
                const double beta=(double)2./(rho_in[cur_index]+neighbor_rho);
                result += beta * scale_uniform    * ((double)(*reinterpret_cast<T*>(y_base_addr + o_grid[i][j][k+1])) - y1_in[cur_index]);}
            
            output[cur_index] = (T)result;
        }
        cur_index++;
      }
    }
    free(offset_grid_ptr);
}

template class Variable_Beta_Laplace_Helper<float,5,2>;
template class Variable_Beta_Laplace_Helper<float,6,2>;

template class Variable_Beta_Laplace_Helper<float,5,3>;
template class Variable_Beta_Laplace_Helper<float,6,3>;
//template class Variable_Beta_Laplace_Helper<float,5,2>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class Variable_Beta_Laplace_Helper<double,5,2>;
template class Variable_Beta_Laplace_Helper<double,6,2>;

template class Variable_Beta_Laplace_Helper<double,5,3>;
template class Variable_Beta_Laplace_Helper<double,6,3>;
#endif
