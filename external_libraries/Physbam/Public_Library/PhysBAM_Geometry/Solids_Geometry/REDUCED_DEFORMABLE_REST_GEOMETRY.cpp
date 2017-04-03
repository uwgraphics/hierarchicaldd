//#####################################################################
// Copyright 2013, Rahul Sheth.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Geometry/Implicit_Objects/IMPLICIT_OBJECT_TRANSFORMED.h>
#include <PhysBAM_Geometry/Solids_Geometry/REDUCED_DEFORMABLE_REST_GEOMETRY.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/POINT_SIMPLICES_1D.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/SEGMENTED_CURVE_2D.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TETRAHEDRALIZED_VOLUME.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_AREA.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <PhysBAM_Tools/Utilities/DEBUG_CAST.h>

namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class TV> REDUCED_DEFORMABLE_REST_GEOMETRY<TV>::
REDUCED_DEFORMABLE_REST_GEOMETRY(REDUCED_DEFORMABLE_GEOMETRY_COLLECTION<TV>* const reduced_deformable_geometry_collection_in,int parent_particle_index_in)
    :reduced_deformable_geometry_collection(reduced_deformable_geometry_collection_in),rest_geometry_particles(reduced_deformable_geometry_collection_in->rest_geometry_particles),
     simplicial_object(NULL),implicit_object(NULL),rest_particle_indices(0),own_structures(true)
{
    if(parent_particle_index_in==-1) return; //for child classes to handle initial setup
    if(parent_particle_index_in) parent_particle_index=parent_particle_index_in;
    else parent_particle_index=reduced_deformable_geometry_collection->reduced_deformable_geometry_particles->array_collection->Add_Element();
    assert(!reduced_deformable_geometry_collection->reduced_deformable_geometry_particles->rest_geometry(parent_particle_index));
    reduced_deformable_geometry_collection->reduced_deformable_geometry_particles->rest_geometry(parent_particle_index)=this;
}
//#####################################################################
// "Copy" Constructor - For Transferring to Child Class
//#####################################################################
template<class TV> REDUCED_DEFORMABLE_REST_GEOMETRY<TV>::
REDUCED_DEFORMABLE_REST_GEOMETRY(REDUCED_DEFORMABLE_REST_GEOMETRY<TV>* rest_geometry_in)
    :reduced_deformable_geometry_collection(rest_geometry_in->reduced_deformable_geometry_collection),rest_geometry_particles(rest_geometry_in->rest_geometry_particles),
     simplicial_object(rest_geometry_in->simplicial_object),implicit_object(rest_geometry_in->implicit_object),parent_particle_index(rest_geometry_in->parent_particle_index),
     structures(rest_geometry_in->structures),rest_particle_indices(rest_geometry_in->rest_particle_indices),own_structures(true)
{
    rest_geometry_in->own_structures=false;
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> REDUCED_DEFORMABLE_REST_GEOMETRY<TV>::
~REDUCED_DEFORMABLE_REST_GEOMETRY()
{
    reduced_deformable_geometry_collection->reduced_deformable_geometry_particles->rest_geometry(parent_particle_index)=NULL;
    if(own_structures) {
        if(simplicial_object) delete simplicial_object;
        if(implicit_object) delete implicit_object; 
        structures.Delete_Pointers_And_Clean_Memory();}
}
//#####################################################################
// Function Add_Structure (1D)
//#####################################################################
template<class TV> template<int d> 
typename ENABLE_IF<d==1>::TYPE REDUCED_DEFORMABLE_REST_GEOMETRY<TV>::Add_Structure(STRUCTURE<VECTOR<T,d> >& structure)
{
    if(POINT_SIMPLICES_1D<T>* point_simplices=dynamic_cast<POINT_SIMPLICES_1D<T>*>(&structure)){
        if(this->simplicial_object) this->Remove_Structure(this->simplicial_object);
        this->simplicial_object=point_simplices;}
    else if(IMPLICIT_OBJECT<TV>* implicit_object_input=dynamic_cast<IMPLICIT_OBJECT<TV>*>(&structure)){
        if(this->implicit_object) this->Remove_Structure(this->implicit_object);
        this->implicit_object=implicit_object_input;}
    this->structures.Append(&structure);
}
//#####################################################################
// Function Add_Structure (2D)
//#####################################################################
template<class TV> template<int d> 
typename ENABLE_IF<d==2>::TYPE REDUCED_DEFORMABLE_REST_GEOMETRY<TV>::Add_Structure(STRUCTURE<VECTOR<T,d> >& structure)
{
    if(SEGMENTED_CURVE_2D<T>* segmented_curve=dynamic_cast<SEGMENTED_CURVE_2D<T>*>(&structure)){
        if(this->simplicial_object) this->Remove_Structure(this->simplicial_object);
        this->simplicial_object=segmented_curve;}
    else if(IMPLICIT_OBJECT<TV>* implicit_object_input=dynamic_cast<IMPLICIT_OBJECT<TV>*>(&structure)){
        if(this->implicit_object) this->Remove_Structure(this->implicit_object);
        this->implicit_object=implicit_object_input;}
    else if(dynamic_cast<TRIANGULATED_AREA<T>*>(&structure)){
        if(TRIANGULATED_AREA<T>* old_triangulated_area=this->template Find_Structure<TRIANGULATED_AREA<T>*>()) this->Remove_Structure(old_triangulated_area);}
    this->structures.Append(&structure);
}
//#####################################################################
// Function Add_Structure (3D)
//#####################################################################
template<class TV> template<int d> 
typename ENABLE_IF<d==3>::TYPE REDUCED_DEFORMABLE_REST_GEOMETRY<TV>::Add_Structure(STRUCTURE<VECTOR<T,d> >& structure)
{
    if(TRIANGULATED_SURFACE<T>* triangulated_surface_input=dynamic_cast<TRIANGULATED_SURFACE<T>*>(&structure)){
        if(this->simplicial_object) this->Remove_Structure(this->simplicial_object);
        this->simplicial_object=triangulated_surface_input;}
    else if(IMPLICIT_OBJECT_TRANSFORMED<TV,FRAME<TV> >* implicit_object_input=dynamic_cast<IMPLICIT_OBJECT_TRANSFORMED<TV,FRAME<TV> >*>(&structure)){
        if(this->implicit_object) this->Remove_Structure(this->implicit_object);
        this->implicit_object=implicit_object_input->object_space_implicit_object;}
    else if(IMPLICIT_OBJECT<TV>* implicit_object_input=dynamic_cast<IMPLICIT_OBJECT<TV>*>(&structure)){
        if(this->implicit_object) this->Remove_Structure(this->implicit_object);
        this->implicit_object=implicit_object_input;}
    else if(dynamic_cast<TETRAHEDRALIZED_VOLUME<T>*>(&structure)){
        if(TETRAHEDRALIZED_VOLUME<T>* old_tetrahedralized_volume=this->template Find_Structure<TETRAHEDRALIZED_VOLUME<T>*>()) this->Remove_Structure(old_tetrahedralized_volume);}
    this->structures.Append(&structure);
}
//#####################################################################
// Function Remove_Structure
//#####################################################################
template<class TV> void REDUCED_DEFORMABLE_REST_GEOMETRY<TV>::
Remove_Structure(STRUCTURE<TV>* structure)
{
    structures.Remove_Index_Lazy(structures.Find(structure));
}
//#####################################################################
// Update_Rest_Particle_Indices
//#####################################################################
template<class TV> void REDUCED_DEFORMABLE_REST_GEOMETRY<TV>::
Update_Rest_Particle_Indices()
{
    ARRAY<int> particle_marks(rest_geometry_particles->array_collection->Size(),false);
    particle_marks.Compact();
    ARRAYS_COMPUTATIONS::Fill(particle_marks,0);
    assert(structures(1)!=NULL); //TODO: don't assume useful rest structure is first
    structures(1)->Mark_Nodes_Referenced(particle_marks,1);
    rest_particle_indices.Resize(0);
    for(int j=1;j<=particle_marks.Size();j++){
        if(particle_marks(j)==1){rest_particle_indices.Append(j);}}
}
//#####################################################################
#define INSTANTIATION_HELPER(T,d) \
    template class REDUCED_DEFORMABLE_REST_GEOMETRY<VECTOR<T,d> >; \
    template void REDUCED_DEFORMABLE_REST_GEOMETRY<VECTOR<T,d> >::Add_Structure<d>(STRUCTURE<VECTOR<T,d> >& structure);
INSTANTIATION_HELPER(float,1);
INSTANTIATION_HELPER(float,2);
INSTANTIATION_HELPER(float,3);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
INSTANTIATION_HELPER(double,1);
INSTANTIATION_HELPER(double,2);
INSTANTIATION_HELPER(double,3);
#endif

}
