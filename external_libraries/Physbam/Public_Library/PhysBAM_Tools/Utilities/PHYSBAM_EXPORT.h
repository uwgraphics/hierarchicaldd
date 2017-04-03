//#####################################################################
// Copyright 2013, Rahul Sheth.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __PHYSBAM_EXPORT__
#define __PHYSBAM_EXPORT__

#if defined(WIN32) && defined (PHYSBAM_BUILD_SHARED_LIBS)
#if defined(PHYSBAM_EXPORT_TO_DLL)
    #define PHYSBAM_EXPORT __declspec(dllexport)
#else
    #define PHYSBAM_EXPORT __declspec(dllimport)
#endif
#else
    #define PHYSBAM_EXPORT
#endif

#endif
