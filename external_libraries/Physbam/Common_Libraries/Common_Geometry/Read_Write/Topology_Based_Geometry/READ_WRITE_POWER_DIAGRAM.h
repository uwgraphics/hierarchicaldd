//#####################################################################
// Copyright 2015, Mridul Aanjaneya.
// This file is part of PhysBAM whose distribution is governed by the license 
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class READ_WRITE_POWER_DIAGRAM
//#####################################################################
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
#ifndef __READ_WRITE_POWER_DIAGRAM__
#define __READ_WRITE_POWER_DIAGRAM__

#include <Common_Geometry/Topology_Based_Geometry/POWER_DIAGRAM.h>
#include <PhysBAM_Tools/Read_Write/Arrays/READ_WRITE_ARRAY.h>
#include <PhysBAM_Tools/Read_Write/Utilities/READ_WRITE.h>

namespace PhysBAM{
template<class RW,class TV>
class Read_Write<POWER_DIAGRAM<TV>,RW>
{
  public:
    static void Read(std::istream& input,POWER_DIAGRAM<TV>& object)
    {
        object.Clean_Memory();
        // read center
        Read_Binary<RW>(input,object.center);

        // read centroid
        Read_Binary<RW>(input,object.centroid);

        // read face vertices
        int face_vertices_size;
        Read_Binary<RW>(input,face_vertices_size);
        object.face_vertices.Resize(face_vertices_size);
        for(int i=1;i<=face_vertices_size;++i){int current_face_vertices_size;
            Read_Binary<RW>(input,current_face_vertices_size);
            object.face_vertices(i).Resize(current_face_vertices_size);
            for(int j=1;j<=current_face_vertices_size;++j) Read_Binary<RW>(input,object.face_vertices(i)(j));}

        // read neighbors
        int neighbors_size;
        Read_Binary<RW>(input,neighbors_size);
        object.neighbors.Resize(neighbors_size);
        for(int i=1;i<=neighbors_size;++i) Read_Binary<RW>(input,object.neighbors(i));
    }

    static void Write(std::ostream& output,const POWER_DIAGRAM<TV>& object)
    {
        // write center
        Write_Binary<RW>(output,object.center);

        // write centroid
        Write_Binary<RW>(output,object.centroid);

        // write face vertices
        Write_Binary<RW>(output,object.face_vertices.Size());
        for(int i=1;i<=object.face_vertices.Size();++i){
            Write_Binary<RW>(output,object.face_vertices(i).Size());
            for(int j=1;j<=object.face_vertices(i).Size();++j)
                Write_Binary<RW>(output,object.face_vertices(i)(j));}

        // write neighbors
        Write_Binary<RW>(output,object.neighbors.Size());
        for(int i=1;i<=object.neighbors.Size();++i)
            Write_Binary<RW>(output,object.neighbors(i));
    }
};
}
#endif
#endif
