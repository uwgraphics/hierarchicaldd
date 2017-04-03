//#####################################################################
// Copyright 2013, Rahul Sheth.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class REDUCED_DEFORMABLE_GEOMETRY_STATE
//#####################################################################
#ifndef __REDUCED_DEFORMABLE_GEOMETRY_STATE__
#define __REDUCED_DEFORMABLE_GEOMETRY_STATE__

#include <PhysBAM_Tools/Matrices/FRAME.h>
#include <PhysBAM_Tools/Vectors/TWIST.h>
#include <PhysBAM_Tools/Vectors/VECTOR_ND.h>

namespace PhysBAM {

//TODO: remove me!
template<class TV> class STANDARD_TESTS_2D;

template<class TV>
class REDUCED_DEFORMABLE_GEOMETRY_STATE
{
    //TODO: remove me!
    template<class TV2> friend class STANDARD_TESTS_2D;
  private:
    typedef typename TV::SCALAR T;
    typedef typename TV::SPIN T_SPIN;
  protected:
    int reduced_dimension;
    //R and t
    FRAME<TV> frame;
    //R_dot and t_dot
    TWIST<TV> twist;
    //q and q_dot
    VECTOR_ND<T> reduced_displacements;
    VECTOR_ND<T> reduced_velocities;
  public:
    REDUCED_DEFORMABLE_GEOMETRY_STATE(const int reduced_dimension_in=0,const FRAME<TV> frame_in=FRAME<TV>(),const TWIST<TV> twist_in=TWIST<TV>())
        :frame(frame_in),twist(twist_in)
    {
        Set_Reduced_Dimension(reduced_dimension_in);
    };
    template<class TV2> explicit REDUCED_DEFORMABLE_GEOMETRY_STATE(const REDUCED_DEFORMABLE_GEOMETRY_STATE<TV2>& state_in)
        :frame(state_in.frame),twist(state_in.twist),reduced_dimension(state_in.reduced_dimension),reduced_displacements(state_in.reduced_displacements),reduced_velocities(state_in.reduced_velocities)
    {};
    ~REDUCED_DEFORMABLE_GEOMETRY_STATE(){};

    //Get and Set Methods
    void Get_Rotation(ROTATION<TV>& rotation_out) const
    {rotation_out=frame.r;};
    void Set_Rotation(const ROTATION<TV>& rotation_in)
    {frame.r=rotation_in;};
    void Get_Angular_Velocity(T_SPIN& angular_out) const
    {angular_out=twist.angular;};
    void Set_Angular_Velocity(const T_SPIN& angular_in)
    {twist.angular=angular_in;};
    void Get_Translation(TV& translation_out) const
    {translation_out=frame.t;};
    void Set_Translation(const TV& translation_in)
    {frame.t=translation_in;};
    void Get_Translational_Velocity(TV& translational_velocity_out) const
    {translational_velocity_out=twist.linear;};
    void Set_Translational_Velocity(const TV& translational_velocity_in)
    {twist.linear=translational_velocity_in;};
    void Get_Reduced_Displacements(VECTOR_ND<T>& reduced_displacements_out) const
    {reduced_displacements_out=reduced_displacements;};
    void Set_Reduced_Displacements(const VECTOR_ND<T>& reduced_displacements_in)
    {reduced_displacements=reduced_displacements_in;};
    void Get_Reduced_Velocities(VECTOR_ND<T>& reduced_velocities_out) const
    {reduced_velocities_out=reduced_velocities;};
    void Set_Reduced_Velocities(const VECTOR_ND<T>& reduced_velocities_in)
    {reduced_velocities=reduced_velocities_in;};
    void Set_Reduced_Dimension(const int reduced_dimension_in)
    {reduced_dimension=reduced_dimension_in;reduced_displacements.Resize(reduced_dimension);reduced_velocities.Resize(reduced_dimension);};
    const int Get_Reduced_Dimension() const
    {return reduced_dimension;};

    //Frame Transformation Methods
    TV Object_Space_Point(const TV& world_space_point) const
    {return frame.Inverse_Times(world_space_point);}
    TV Object_Space_Vector(const TV& world_space_vector) const
    {return frame.r.Inverse_Rotate(world_space_vector);}
    TV World_Space_Point(const TV& object_space_point) const
    {return frame*object_space_point;}
    TV World_Space_Vector(const TV& object_space_vector) const
    {return frame.r.Rotate(object_space_vector);}
//#####################################################################
};

}
#endif
