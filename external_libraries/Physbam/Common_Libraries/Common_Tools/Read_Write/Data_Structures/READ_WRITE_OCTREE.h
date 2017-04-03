//#####################################################################
// Copyright 2015, Mridul Aanjaneya, Ming Gao.
// This file is part of PhysBAM whose distribution is governed by the license 
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class READ_WRITE_OCTREE
//#####################################################################
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
#ifndef __READ_WRITE_OCTREE__
#define __READ_WRITE_OCTREE__

#include <Common_Tools/Data_Structures/OCTREE.h>
#include <PhysBAM_Tools/Read_Write/Arrays/READ_WRITE_ARRAY.h>
#include <PhysBAM_Tools/Read_Write/Utilities/READ_WRITE.h>
#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>

namespace PhysBAM{
template<class RW,class T>
class Read_Write<OCTREE<T>,RW>
{
    typedef VECTOR<T,3> TV;
  public:
    static void Read(std::istream& input,OCTREE<T>& object)
    {
        object.Clean_Memory();
        Read_Binary<RW>(input,object.domain);
        int number_of_nodes=0;
        Read_Binary<RW>(input,number_of_nodes);
        for(int i=1;i<=number_of_nodes;++i){TV center,dx_over_two;
            Read_Binary<RW>(input,center);
            Read_Binary<RW>(input,dx_over_two);
            ARRAY<int> children(8);
            for(int j=1;j<=8;++j) Read_Binary<RW>(input,children(j));
            OCTREE_NODE<T> *node=new OCTREE_NODE<T>(center,dx_over_two,children);
            object.nodes.Append(node);}
    }

    static void Write(std::ostream& output,const OCTREE<T>& object)
    {
        Write_Binary<RW>(output,object.domain);
        Write_Binary<RW>(output,object.nodes.Size());
        for(int i=1;i<=object.nodes.Size();++i){
            Write_Binary<RW>(output,object.nodes(i)->center);
            Write_Binary<RW>(output,object.nodes(i)->dx_over_two);
            for(int j=1;j<=8;++j) Write_Binary<RW>(output,object.nodes(i)->children(j));}
    }
};
}
#endif
#endif
