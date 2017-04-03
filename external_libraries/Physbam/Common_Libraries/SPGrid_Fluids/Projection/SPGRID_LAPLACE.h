//#####################################################################
// Copyright 2012-2013, Sean Bauer, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SPGRID_LAPLACE
//#####################################################################
#ifndef __SPGRID_LAPLACE__
#define __SPGRID_LAPLACE__

#include <SPGrid/Core/SPGrid_Set.h>
#include <SPGrid_Fluids/Simulation/FLUIDS_SIMULATION_FLAGS.h>
#include <SPGrid_Fluids/Grids/GRID_TOPOLOGY_HELPER.h>

using namespace PhysBAM;

namespace SPGrid{

template<class Data_array_type,class Flag_array_type>
class SPGRID_LAPLACE
{
    typedef typename Flag_array_type::MASK Flag_array_mask;
    typedef typename Data_array_type::DATA T;
    enum{ d=Data_array_type::dim};

public:

    static void Compute(const std::pair<const unsigned long*,unsigned>& blocks,const Data_array_type u,Data_array_type Lu,Flag_array_type flags,const double scale_uniform,const double scale_nonuniform)
    {
        static const int number_of_face_neighbors=GRID_TOPOLOGY_HELPER<Flag_array_mask>::faces_per_cell;
        unsigned long face_neighbor_offsets[number_of_face_neighbors];
        GRID_TOPOLOGY_HELPER<Flag_array_mask>::Face_Neighbor_Offsets(face_neighbor_offsets); // order is (-x, x, -y, y, -z, z)

        for (SPGrid_Block_Iterator<Flag_array_mask> iterator(blocks);iterator.Valid();iterator.Next())
        {
            unsigned flag = iterator.Data(flags);
            if(flag & (SPGrid_Cell_Type_Active|SPGrid_Cell_Type_Ghost)){             
                double cell_value=(double)(iterator.Data(u));
                double result=(double)0.;
                for(int face=0;face<number_of_face_neighbors;face++){
                    const int side=face%2+1;
                    const int axis=face/2+1;
                    unsigned long offset = face_neighbor_offsets[face];                    
                    if (flag & GRID_TOPOLOGY_HELPER<Flag_array_mask>::Face_Scaled_Mask(axis,side)){
                        double neighbor_value=(double)(iterator.Data(u,offset));
                        result += scale_nonuniform*(neighbor_value - cell_value);}
                    else if(flag & GRID_TOPOLOGY_HELPER<Flag_array_mask>::Face_Active_Mask(axis,side)){
                        double neighbor_value=(double)(iterator.Data(u,offset));
                        result += scale_uniform*(neighbor_value - cell_value);}
                }
                iterator.Data(Lu) = (T)result;
            }
        }
    }
//#####################################################################
};
}
#endif
