#define GL_GLEXT_PROTOTYPES 1

#include <Common_Rendering/OpenGL/OPENGL_GLSL_SHADER.h>
#include <PhysBAM_Tools/Data_Structures/HASHTABLE_ITERATOR.h>
#include <PhysBAM_Tools/Log/LOG.h>

#ifdef _WIN32
#include <windows.h>
#endif

#ifndef __APPLE__
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>
#include <GL/glext.h>
#else
#include <GLUT/glut.h>
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#include <OpenGL/glext.h>
#endif

using namespace PhysBAM;


bool OPENGL_GLSL_SHADER::
CreateProgram() {
    GLint compiled, linked;
    
    _vshader = glCreateShader( GL_VERTEX_SHADER );
    glShaderSource( _vshader, 1, vertex_shader, NULL );
    glCompileShader( _vshader );
    glGetShaderiv( _vshader, GL_COMPILE_STATUS, &compiled );
    if( !compiled ){
        GLint length;
        GLchar* log;
        glGetShaderiv( _vshader, GL_INFO_LOG_LENGTH, &length );
        log = (GLchar*)( malloc( length ) );
        glGetShaderInfoLog( _vshader, length, &length, log );
        LOG::cerr << "Vertex shader compile log = " << log << std::endl;
        free( log );
        _program_ready = false;
        return _program_ready;
    }

    _fshader = glCreateShader( GL_FRAGMENT_SHADER );
    glShaderSource( _fshader, 1, fragment_shader, NULL );
    glCompileShader( _fshader );
    glGetShaderiv( _fshader, GL_COMPILE_STATUS, &compiled );
    if( !compiled ){
        GLint length;
        GLchar* log;
        glGetShaderiv( _fshader, GL_INFO_LOG_LENGTH, &length );
        log = (GLchar*)( malloc( length ) );
        glGetShaderInfoLog( _fshader, length, &length, log );
        LOG::cerr << "Fragment shader compile log = " << log << std::endl;
        free( log );
        _program_ready = false;
        return _program_ready;
    }

    _program = glCreateProgram();
    glAttachShader( _program, _vshader );
    glAttachShader( _program, _fshader );

    glLinkProgram(_program);
                     
    glGetProgramiv( _program, GL_LINK_STATUS, &linked );
    if( !linked ){
        GLint length;
        GLchar* log;
        glGetProgramiv( _program, GL_INFO_LOG_LENGTH, &length );
        log = (GLchar*)( malloc( length ) );
        glGetProgramInfoLog( _program, length, &length, log );
        LOG::cerr << "Shader Link log = " << log << std::endl;
        free( log );
        _program_ready = false;
        return _program_ready;
    }
        
    for( HASHTABLE_ITERATOR<std::string,int> biter( attr_bindings ); biter.Valid(); biter.Next() ){
        int location = glGetAttribLocation( _program,  biter.Key().c_str() );
        if( location != -1 )
            biter.Data() = location;
        else{
            biter.Data() = -1;
            LOG::cout << "Could not find a location for custom attribute '" << biter.Key() << "'." << std::endl;
        }
    }


    for( HASHTABLE_ITERATOR<std::string,int> biter( unif_bindings ); biter.Valid(); biter.Next() ){
        int location = glGetUniformLocation( _program,  biter.Key().c_str() );
        if( location != -1 )
            biter.Data() = location;
        else{
            biter.Data() = -1;
            LOG::cout << "Could not find a location for uniform '" << biter.Key() << "'." << std::endl;
        }
    }

    _program_ready = true;
    return _program_ready;
}

void OPENGL_GLSL_SHADER::
UseProgram() const {
    PHYSBAM_ASSERT( _program_ready );
    glUseProgram( _program );
}

void OPENGL_GLSL_SHADER::
DisableProgram() const {
    glUseProgram( 0 );
}
