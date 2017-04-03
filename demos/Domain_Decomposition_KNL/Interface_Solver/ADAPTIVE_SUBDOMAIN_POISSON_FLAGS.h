//#####################################################################
//Removed For Anonymity: Copyright Authors of Submission pap173s1
//Removed For Anonymity: Copyright Authors of Submission pap173s1
//#####################################################################
// Class ADAPTIVE_SUBDOMAIN_POISSON_FLAGS
//#####################################################################
#ifndef __ADAPTIVE_SUBDOMAIN_POISSON_FLAGS__
#define __ADAPTIVE_SUBDOMAIN_POISSON_FLAGS__

namespace PhysBAM{
enum {
    Uniform_Node_Type_Active             = 0x00000010u,
    Uniform_Node_Type_Dirichlet          = 0x00000020u,
    Uniform_Node_Type_Interface          = 0x00000040u,
    Adaptive_Cell_Type_Minimally_Refined = 0x00000080u,
    Uniform_Cell_Traversed               = 0x01000000u,
    Uniform_Cell_Coarsened               = 0x02000000u,
    Adaptive_Cell_Type_Interior          = 0x00000001u,
    Adaptive_Node_Degree_Marker          = 0x10000000u,//This is simply used to mark if a node is a ADAPTED INTERIOR dof at that level.
    Adaptive_Node_Type_T_Junction        = 0x00000100u,//This is simply used to mark if a node is a ADAPTED INTERIOR dof at that level.
    Adaptive_Node_Type_Coarse_Shared     = 0x00000200u,//This is simply used to mark if a node is a ADAPTED INTERIOR dof at that level.
    Adaptive_Node_Type_Active            = 0x00000400u,//This is simply used to mark if a node is a ADAPTED INTERIOR dof at that level.
    // Uniform_Cell_Type_Dirichlet = 0x00000002u, // Set by user
    // SPGrid_Cell_Type_Ghost     = 0x00000004u, // Generated automatically
    // SPGrid_Cell_Type_Active    = 0x00000008u, // Active = Interior but not Dirichlet, also not fully su
};
}
#endif
