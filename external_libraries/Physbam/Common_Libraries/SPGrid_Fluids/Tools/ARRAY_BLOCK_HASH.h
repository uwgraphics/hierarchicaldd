//#####################################################################
// Copyright 2012-2013, Sean Bauer, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ARRAY_BLOCK_HASH
//#####################################################################
#ifndef __ARRAY_BLOCK_HASH_h__
#define __ARRAY_BLOCK_HASH_h__

#include <PhysBAM_Tools/Math_Tools/Hash.h>
#include <SPGrid/Tools/SPGrid_Block_Iterator.h>

namespace PhysBAM{
using namespace SPGrid;

template<class T_ARRAY>
int ARRAY_BLOCK_HASH(T_ARRAY array,const std::pair<const unsigned long*,unsigned>& blocks)
{
    typedef typename T_ARRAY::DATA T;
    static const unsigned int elements_per_block = T_ARRAY::MASK::elements_per_block;
    const unsigned& size=blocks.second;

    ARRAY<int> hashes;
    if(size%2==0) hashes.Append(HASH().value);

    for(SPGrid_Block_Iterator<typename T_ARRAY::MASK> iterator(blocks);iterator.Valid();iterator.Next_Block()){
        PhysBAM::ARRAY_VIEW<T> block_array(elements_per_block,&iterator.Data(array));
        hashes.Append(Hash(block_array));
        if(hashes.m==3){
            hashes(1)=HASH(hashes(1),hashes(2),hashes(3)).value;
            hashes.Resize(1);}}

    PHYSBAM_ASSERT(hashes.m==1);
    return hashes(1);
}

}
#endif
