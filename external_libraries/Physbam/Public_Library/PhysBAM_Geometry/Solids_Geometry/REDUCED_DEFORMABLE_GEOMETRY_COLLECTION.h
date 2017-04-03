//#####################################################################
// Copyright 2013, Rahul Sheth.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class REDUCED_DEFORMABLE_GEOMETRY_COLLECTION
//#####################################################################
#ifndef __REDUCED_DEFORMABLE_GEOMETRY_COLLECTION__
#define __REDUCED_DEFORMABLE_GEOMETRY_COLLECTION__

#include <PhysBAM_Geometry/Geometry_Particles/REDUCED_DEFORMABLE_GEOMETRY_PARTICLES.h>
#include <PhysBAM_Geometry/Solids_Geometry/REDUCED_DEFORMABLE_GEOMETRY_STATE.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/STRUCTURE.h>

namespace PhysBAM{

template<class TV>
class REDUCED_DEFORMABLE_GEOMETRY_COLLECTION : public NONCOPYABLE
{
    typedef typename TV::SCALAR T;
public:
    //"rigid" array collection - 1 particle per body
    REDUCED_DEFORMABLE_GEOMETRY_PARTICLES<TV>* reduced_deformable_geometry_particles;
    //"deformable" array collection - n particles per body
    GEOMETRY_PARTICLES<TV>* rest_geometry_particles;
    bool own_particles;

    REDUCED_DEFORMABLE_GEOMETRY_COLLECTION();
    REDUCED_DEFORMABLE_GEOMETRY_COLLECTION(REDUCED_DEFORMABLE_GEOMETRY_PARTICLES<TV>* reduced_deformable_geometry_particles_in,GEOMETRY_PARTICLES<TV>* geometry_particles_in);
    virtual ~REDUCED_DEFORMABLE_GEOMETRY_COLLECTION();

//#####################################################################
    //returns new "rigid" particle index
    int Add_Body(REDUCED_DEFORMABLE_GEOMETRY_STATE<TV>& state_in,int index_in=0);

#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
    STRUCTURE<TV>* Find_Or_Read_Structure(const STREAM_TYPE stream_type,const std::string& filename,const T scaling_factor);
    void Read(const STREAM_TYPE stream_type,const std::string& prefix,const int frame,const bool include_static_variables);
    void Write(const STREAM_TYPE stream_type,const std::string& prefix,const int frame,const bool include_static_variables,const bool write_basis=true) const;
private:
    void Read_Static_Variables(const STREAM_TYPE stream_type,const std::string& prefix);
    void Read_Dynamic_Variables(const STREAM_TYPE stream_type,const std::string& prefix,const int frame);
    void Write_Static_Variables(const STREAM_TYPE stream_type,const std::string& prefix) const;
    void Write_Basis(const STREAM_TYPE stream_type,const std::string& prefix) const;
    void Write_Dynamic_Variables(const STREAM_TYPE stream_type,const std::string& prefix,const int frame) const;
#endif
//#####################################################################
};

}
#endif
