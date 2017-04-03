//#####################################################################
// Copyright 2003-2009, Zhaosheng Bao, Ronald Fedkiw, Geoffrey Irving, Igor Neverov, Craig Schroeder, Andrew Selle, Tamar Shinar, Eftychios Sifakis, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class DEFORMABLE_GEOMETRY_COLLECTION
//#####################################################################
#ifndef __DEFORMABLE_GEOMETRY_COLLECTION__
#define __DEFORMABLE_GEOMETRY_COLLECTION__

#include <PhysBAM_Tools/Utilities/Find_Type.h>
#include <PhysBAM_Geometry/Geometry_Particles/GEOMETRY_PARTICLES.h>
namespace PhysBAM{

template<class TV> class STRUCTURE;

template<class TV>
class DEFORMABLE_GEOMETRY_COLLECTION:public NONCOPYABLE
{
    typedef typename TV::SCALAR T;
public:
    GEOMETRY_PARTICLES<TV>& particles;
    ARRAY<STRUCTURE<TV>*> structures;
    bool own_structures;

    DEFORMABLE_GEOMETRY_COLLECTION(GEOMETRY_PARTICLES<TV>& particles_input,bool own_structures_input=true);
    virtual ~DEFORMABLE_GEOMETRY_COLLECTION();

    int Add_Structure(STRUCTURE<TV>* structure)
    {return structures.Append(structure);}

    template<class T_STRUCTURE> T_STRUCTURE
    Find_Structure(const int index=1)
    {return Find_Type<T_STRUCTURE>(structures,index);}

    template<class T_STRUCTURE> const T_STRUCTURE
    Find_Structure(const int index=1) const
    {return Find_Type<T_STRUCTURE>(structures,index);}

#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
    void Read(const STREAM_TYPE stream_type,const std::string& prefix,const int frame,const int static_frame,const bool include_static_variables,const std::string& postfix=std::string())
    {Read(stream_type,prefix,prefix,frame,static_frame,include_static_variables,postfix);}

    void Write(const STREAM_TYPE stream_type,const std::string& prefix,const int frame,const int static_frame,const bool include_static_variables,const std::string& postfix=std::string()) const
    {Write(stream_type,prefix,prefix,frame,static_frame,include_static_variables,postfix);}

//#####################################################################
    void Read(const STREAM_TYPE,const std::string& prefix,const std::string& static_prefix,const int frame,const int static_frame,const bool include_static_variables,const std::string& postfix=std::string());
    void Write(const STREAM_TYPE,const std::string& prefix,const std::string& static_prefix,const int frame,const int static_frame,const bool include_static_variables,const std::string& postfix=std::string()) const;
private:
    void Read_Static_Variables(const STREAM_TYPE stream_type,const std::string& prefix,const int frame,const std::string& postfix=std::string());
    void Read_Dynamic_Variables(const STREAM_TYPE stream_type,const std::string& prefix,const int frame,const std::string& postfix=std::string());
    void Write_Static_Variables(const STREAM_TYPE stream_type,const std::string& prefix,const int frame,const std::string& postfix=std::string()) const;
    void Write_Dynamic_Variables(const STREAM_TYPE stream_type,const std::string& prefix,const int frame,const std::string& postfix=std::string()) const;
//#####################################################################
#endif
};
}
#endif
