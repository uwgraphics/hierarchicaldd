//#####################################################################
// Copyright 2015, Mridul Aanjaneya, Ming Gao.
// This file is part of PhysBAM whose distribution is governed by the license 
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class READ_WRITE_QUADTREE
//#####################################################################
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
#ifndef __READ_WRITE_QUADTREE__
#define __READ_WRITE_QUADTREE__

#include <Common_Tools/Data_Structures/QUADTREE.h>
#include <PhysBAM_Tools/Read_Write/Arrays/READ_WRITE_ARRAY.h>
#include <PhysBAM_Tools/Read_Write/Utilities/READ_WRITE.h>
#include <PhysBAM_Tools/Vectors/VECTOR_2D.h>

namespace PhysBAM{
template<class RW,class T>
class Read_Write<QUADTREE<T>,RW>
{
    typedef VECTOR<T,2> TV;
  public:
    static void Read(std::istream& input,QUADTREE<T>& object)
    {
        object.Clean_Memory();
        Read_Binary<RW>(input,object.domain);
        int number_of_nodes=0;
        Read_Binary<RW>(input,number_of_nodes);
        for(int i=1;i<=number_of_nodes;++i){TV center,dx_over_two;
            Read_Binary<RW>(input,center);
            Read_Binary<RW>(input,dx_over_two);
            int north_east_index,north_west_index,south_east_index,south_west_index;
            int index,parent_index;
            Read_Binary<RW>(input,north_east_index);
            Read_Binary<RW>(input,north_west_index);
            Read_Binary<RW>(input,south_east_index);
            Read_Binary<RW>(input,south_west_index);
            Read_Binary<RW>(input,parent_index);
            Read_Binary<RW>(input,index);
            QUADTREE_NODE<T> *node=new QUADTREE_NODE<T>(center,dx_over_two,north_east_index,north_west_index,south_east_index,south_west_index,parent_index,index);
            object.nodes.Append(node);}
    }

    static void Write(std::ostream& output,const QUADTREE<T>& object)
    {
        Write_Binary<RW>(output,object.domain);
        Write_Binary<RW>(output,object.nodes.Size());
        for(int i=1;i<=object.nodes.Size();++i){
            Write_Binary<RW>(output,object.nodes(i)->center);
            Write_Binary<RW>(output,object.nodes(i)->dx_over_two);
            Write_Binary<RW>(output,object.nodes(i)->north_east_index);
            Write_Binary<RW>(output,object.nodes(i)->north_west_index);
            Write_Binary<RW>(output,object.nodes(i)->south_east_index);
            Write_Binary<RW>(output,object.nodes(i)->south_west_index);
            Write_Binary<RW>(output,object.nodes(i)->parent_index);
            Write_Binary<RW>(output,object.nodes(i)->index);}
    }
};
}
#endif
#endif
