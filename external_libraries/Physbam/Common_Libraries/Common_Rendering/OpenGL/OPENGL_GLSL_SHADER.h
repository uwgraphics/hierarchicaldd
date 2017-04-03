//#####################################################################
// Copyright 2015. Nathan Mitchell
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
#ifndef __OPENGL_GLSL_SHADER_H__
#define __OPENGL_GLSL_SHADER_H__

#include <string>
#include <PhysBAM_Tools/Data_Structures/HASHTABLE.h>
#include <stdint.h>

namespace PhysBAM {

    class OPENGL_GLSL_SHADER{
    public:
        const char** vertex_shader;
        const char** fragment_shader;
        PhysBAM::HASHTABLE<std::string, int> attr_bindings;
        PhysBAM::HASHTABLE<std::string, int> unif_bindings;
        
        uint32_t _vshader, _fshader, _program;
        bool _program_ready;

        bool CreateProgram();
        void UseProgram() const;
        void DisableProgram() const;

    };

}
#endif
