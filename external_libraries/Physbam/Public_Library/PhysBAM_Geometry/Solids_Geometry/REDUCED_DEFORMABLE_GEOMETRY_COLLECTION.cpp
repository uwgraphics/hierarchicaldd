//#####################################################################
// Copyright 2013, Rahul Sheth.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Geometry/Solids_Geometry/REDUCED_DEFORMABLE_GEOMETRY_COLLECTION.h>
#include <PhysBAM_Geometry/Solids_Geometry/REDUCED_DEFORMABLE_REST_GEOMETRY.h>
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
#include <PhysBAM_Tools/Read_Write/Point_Clouds/READ_WRITE_POINT_CLOUD.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Geometry/Read_Write/Geometry/READ_WRITE_STRUCTURE.h>
#endif
namespace PhysBAM {
//#####################################################################
// Constructor (owns reduced_deformable_geometry_particles)
//#####################################################################
template<class TV> REDUCED_DEFORMABLE_GEOMETRY_COLLECTION<TV>::
REDUCED_DEFORMABLE_GEOMETRY_COLLECTION()
{
    own_particles=true;
    reduced_deformable_geometry_particles=new REDUCED_DEFORMABLE_GEOMETRY_PARTICLES<TV>();
    rest_geometry_particles=new GEOMETRY_PARTICLES<TV>();
};
//#####################################################################
// Constructor (doesn't own reduced_deformable_geometry_particles)
//#####################################################################
template<class TV> REDUCED_DEFORMABLE_GEOMETRY_COLLECTION<TV>::
REDUCED_DEFORMABLE_GEOMETRY_COLLECTION(REDUCED_DEFORMABLE_GEOMETRY_PARTICLES<TV>* reduced_deformable_geometry_particles_in,GEOMETRY_PARTICLES<TV>* geometry_particles_in)
{
    own_particles=false;
    reduced_deformable_geometry_particles=reduced_deformable_geometry_particles_in;
    rest_geometry_particles=geometry_particles_in;
};
//#####################################################################
// Destructor
//#####################################################################
template<class TV> REDUCED_DEFORMABLE_GEOMETRY_COLLECTION<TV>::
~REDUCED_DEFORMABLE_GEOMETRY_COLLECTION()
{
    if(own_particles) {delete reduced_deformable_geometry_particles;delete rest_geometry_particles;}
};
//#####################################################################
// Function Add_Body
//#####################################################################
template<class TV> int REDUCED_DEFORMABLE_GEOMETRY_COLLECTION<TV>::
Add_Body(REDUCED_DEFORMABLE_GEOMETRY_STATE<TV>& state_in,int index_in)
{
    //set up rest geometry first, which also creates the reduced_deformable_geometry_particle at index
    REDUCED_DEFORMABLE_REST_GEOMETRY<TV>* rdrg=new REDUCED_DEFORMABLE_REST_GEOMETRY<TV>(this,index_in);
    int index=rdrg->parent_particle_index;
    reduced_deformable_geometry_particles->Set_State(state_in,index);
    return index;
}
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
//#####################################################################
// Function Find_Or_Read_Structure
//#####################################################################
template<class TV> STRUCTURE<TV> * REDUCED_DEFORMABLE_GEOMETRY_COLLECTION<TV>::
Find_Or_Read_Structure(const STREAM_TYPE stream_type,const std::string& filename,const T scaling_factor)
{
    STRUCTURE<TV>* structure=0;
    if(!FILE_UTILITIES::File_Exists(filename)){return structure;};
    if(!stream_type.use_doubles) {structure=Read_Write<STRUCTURE<TV>,float>::Create_From_File(filename);}
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
    else {structure=Read_Write<STRUCTURE<TV>,double>::Create_From_File(filename);}
#endif
    if(scaling_factor!=1){structure->Rescale(scaling_factor);}
    return structure;
}
//#####################################################################
// Function Read
//#####################################################################
template<class TV> void REDUCED_DEFORMABLE_GEOMETRY_COLLECTION<TV>::
Read(const STREAM_TYPE stream_type,const std::string& prefix,const int frame,const bool include_static_variables)
{
    Read_Dynamic_Variables(stream_type,prefix,frame);
    if(include_static_variables) Read_Static_Variables(stream_type,prefix);
}
//#####################################################################
// Function Write
//#####################################################################
template<class TV> void REDUCED_DEFORMABLE_GEOMETRY_COLLECTION<TV>::
Write(const STREAM_TYPE stream_type,const std::string& prefix,const int frame,const bool include_static_variables,const bool write_basis) const
{
    Write_Dynamic_Variables(stream_type,prefix,frame);
    if(include_static_variables) Write_Static_Variables(stream_type,prefix);
    if(write_basis) Write_Basis(stream_type,prefix);
}
//#####################################################################
// Function Read_Static_Variables
//#####################################################################
template<class TV> void REDUCED_DEFORMABLE_GEOMETRY_COLLECTION<TV>::
Read_Static_Variables(const STREAM_TYPE stream_type,const std::string& prefix)
{
    FILE_UTILITIES::Read_From_File(stream_type,prefix+"/common/reduced_deformable_rest_particles",*rest_geometry_particles);
    std::istream* input_raw=FILE_UTILITIES::Safe_Open_Input(prefix+"/common/reduced_deformable_rest_geometry");
    std::istream* input_basis_raw=FILE_UTILITIES::Safe_Open_Input(prefix+"/common/reduced_deformable_bases");
    TYPED_ISTREAM input(*input_raw,stream_type);
    TYPED_ISTREAM input_basis(*input_basis_raw,stream_type);
    int number_of_particles,number_of_structures;
    Read_Binary(input,number_of_particles);
    for(int j=1;j<=number_of_particles;j++) {
        Read_Binary(input,number_of_structures);
        if(!reduced_deformable_geometry_particles->rest_geometry(j)) {
            reduced_deformable_geometry_particles->rest_geometry(j)=new REDUCED_DEFORMABLE_REST_GEOMETRY<TV>(this,j);}
        if(!reduced_deformable_geometry_particles->reduced_basis(j)) {
            reduced_deformable_geometry_particles->reduced_basis(j)=new REDUCED_DEFORMABLE_BASIS<TV>();}
        int m=0,n=0;
        Read_Binary(input_basis,m,n);
        reduced_deformable_geometry_particles->reduced_basis(j)->reduced_basis.Resize(m,n);
        ARRAY<STRUCTURE<TV>*>& structure_array=reduced_deformable_geometry_particles->rest_geometry(j)->structures;
        if(!stream_type.use_doubles) {
            Read_Binary_Array<float>(*input_basis_raw,reduced_deformable_geometry_particles->reduced_basis(j)->reduced_basis.x,m*n);}
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
        else {
            Read_Binary_Array<double>(*input_basis_raw,reduced_deformable_geometry_particles->reduced_basis(j)->reduced_basis.x,m*n);}
#endif
        if(!structure_array.m){
            structure_array.Resize(number_of_structures);
            if(!stream_type.use_doubles) {
                for(int k=1;k<=structure_array.m;k++) {
                    structure_array(k)=Read_Write<STRUCTURE<TV>,float>::Create_Structure(*input_raw,*rest_geometry_particles); }}
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
            else {
                for(int k=1;k<=structure_array.m;k++) {
                    structure_array(k)=Read_Write<STRUCTURE<TV>,double>::Create_Structure(*input_raw,*rest_geometry_particles); }}
#endif
        }
        else if(structure_array.m<=number_of_structures){
            int old_number_of_structures=structure_array.m;structure_array.Resize(number_of_structures);
            if(!stream_type.use_doubles){
                for(int k=1;k<=old_number_of_structures;k++){
                    Read_Write<STRUCTURE<TV>,float>::Read_Structure(*input_raw,*structure_array(k));}
                for(int k=old_number_of_structures+1;k<=number_of_structures;k++){
                    structure_array(k)=Read_Write<STRUCTURE<TV>,float>::Create_Structure(*input_raw,*rest_geometry_particles);}}
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
            else{
                for(int k=1;k<=old_number_of_structures;k++){
                    Read_Write<STRUCTURE<TV>,double>::Read_Structure(*input_raw,*structure_array(k));}
                for(int k=old_number_of_structures+1;k<=number_of_structures;k++){
                    structure_array(k)=Read_Write<STRUCTURE<TV>,double>::Create_Structure(*input_raw,*rest_geometry_particles);}}
#endif
        }
        else{
            std::stringstream ss;ss<<"Current number of structures in reduced deformable "<<j<<" ("<<structure_array.m<<") is greater than number in file ("<<number_of_structures<<").";LOG::filecout(ss.str());
            PHYSBAM_FATAL_ERROR();}
        reduced_deformable_geometry_particles->rest_geometry(j)->Update_Rest_Particle_Indices();
    }
    delete input_raw;
    delete input_basis_raw;
}
//#####################################################################
// Function Read_Dynamic_Variables
//#####################################################################
template<class TV> void REDUCED_DEFORMABLE_GEOMETRY_COLLECTION<TV>::
Read_Dynamic_Variables(const STREAM_TYPE stream_type,const std::string& prefix,const int frame)
{
    ARRAY<REDUCED_DEFORMABLE_REST_GEOMETRY<TV>*> geometry(reduced_deformable_geometry_particles->array_collection->Size());
    ARRAY<REDUCED_DEFORMABLE_BASIS<TV>*> bases(reduced_deformable_geometry_particles->array_collection->Size());
    ARRAYS_COMPUTATIONS::Copy_With_Offset(reduced_deformable_geometry_particles->rest_geometry,geometry,0);
    ARRAYS_COMPUTATIONS::Copy_With_Offset(reduced_deformable_geometry_particles->reduced_basis,bases,0);
    FILE_UTILITIES::Read_From_File(stream_type,prefix+"/"+FILE_UTILITIES::Number_To_String(frame)+"/reduced_deformable_particles",*reduced_deformable_geometry_particles);
    ARRAYS_COMPUTATIONS::Copy_With_Offset(geometry,reduced_deformable_geometry_particles->rest_geometry,0);
    ARRAYS_COMPUTATIONS::Copy_With_Offset(bases,reduced_deformable_geometry_particles->reduced_basis,0);
}
//#####################################################################
// Function Write_Static_Variables
//#####################################################################
template<class TV> void REDUCED_DEFORMABLE_GEOMETRY_COLLECTION<TV>::
Write_Static_Variables(const STREAM_TYPE stream_type,const std::string& prefix) const
{
    std::ostream* output_raw=FILE_UTILITIES::Safe_Open_Output(prefix+"/common/reduced_deformable_rest_geometry");
    TYPED_OSTREAM output(*output_raw,stream_type);
    int number_of_particles=reduced_deformable_geometry_particles->array_collection->Size();
    Write_Binary(output,number_of_particles);
    for(int j=1;j<=number_of_particles;j++) {
        ARRAY<STRUCTURE<TV>*>& structure_array=reduced_deformable_geometry_particles->rest_geometry(j)->structures;
        Write_Binary(output,structure_array.m);
        if(!stream_type.use_doubles) {
            for(int k=1;k<=structure_array.m;k++) {
                Read_Write<STRUCTURE<TV>,float>::Write_Structure(*output_raw,*structure_array(k)); }}
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
        else {
            for(int k=1;k<=reduced_deformable_geometry_particles->rest_geometry(j)->structures.m;k++) {
                Read_Write<STRUCTURE<TV>,double>::Write_Structure(*output_raw,*structure_array(k)); }}
#endif
    }
    delete output_raw;
    FILE_UTILITIES::Write_To_File(stream_type,prefix+"/common/reduced_deformable_rest_particles",*rest_geometry_particles);
}
//#####################################################################
// Function Write_Basis
//#####################################################################
template<class TV> void REDUCED_DEFORMABLE_GEOMETRY_COLLECTION<TV>::
Write_Basis(const STREAM_TYPE stream_type,const std::string& prefix) const
{
    std::ostream* output_basis_raw=FILE_UTILITIES::Safe_Open_Output(prefix+"/common/reduced_deformable_bases");
    TYPED_OSTREAM output_basis(*output_basis_raw,stream_type);
    int number_of_particles=reduced_deformable_geometry_particles->array_collection->Size();
    for(int j=1;j<=number_of_particles;j++) {
        Write_Binary(output_basis,reduced_deformable_geometry_particles->reduced_basis(j)->reduced_basis.Rows(),
                     reduced_deformable_geometry_particles->reduced_basis(j)->reduced_basis.Columns());
        MATRIX_MXN<T> scaled_basis=reduced_deformable_geometry_particles->scaling_factors(j)*sqrt((T)reduced_deformable_geometry_particles->reduced_basis(j)->reduced_basis.Rows())*reduced_deformable_geometry_particles->reduced_basis(j)->reduced_basis;
        if(!stream_type.use_doubles) {
            Write_Binary_Array<float>(*output_basis_raw,scaled_basis.x,reduced_deformable_geometry_particles->reduced_basis(j)->reduced_basis.Rows()*reduced_deformable_geometry_particles->reduced_basis(j)->reduced_basis.Columns());}
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
        else {
            Write_Binary_Array<double>(*output_basis_raw,scaled_basis.x,reduced_deformable_geometry_particles->reduced_basis(j)->reduced_basis.Rows()*reduced_deformable_geometry_particles->reduced_basis(j)->reduced_basis.Columns());}
#endif
    }
    delete output_basis_raw;
}
//#####################################################################
// Function Write_Dynamic_Variables
//#####################################################################
template<class TV> void REDUCED_DEFORMABLE_GEOMETRY_COLLECTION<TV>::
Write_Dynamic_Variables(const STREAM_TYPE stream_type,const std::string& prefix,const int frame) const
{
    FILE_UTILITIES::Write_To_File(stream_type,prefix+"/"+FILE_UTILITIES::Number_To_String(frame)+"/reduced_deformable_particles",*reduced_deformable_geometry_particles);
}
#endif
//#####################################################################
#define INSTANTIATION_HELPER(T) \
    template class REDUCED_DEFORMABLE_GEOMETRY_COLLECTION<VECTOR<T,1> >; \
    template class REDUCED_DEFORMABLE_GEOMETRY_COLLECTION<VECTOR<T,2> >; \
    template class REDUCED_DEFORMABLE_GEOMETRY_COLLECTION<VECTOR<T,3> >;
INSTANTIATION_HELPER(float)
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
INSTANTIATION_HELPER(double)
#endif
}
