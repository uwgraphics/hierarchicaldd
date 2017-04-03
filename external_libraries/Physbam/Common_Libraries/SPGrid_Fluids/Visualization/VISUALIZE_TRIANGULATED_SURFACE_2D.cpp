//#####################################################################
// Copyright (c) 2014, Raj Setaluri, Mridul Aanjneya
// Distributed under the FreeBSD license (see license.txt)
//#####################################################################
#include <SPGrid_Fluids/Visualization/VISUALIZE_TRIANGULATED_SURFACE.h>
#include <SPGrid_Fluids/Simulation/FLUIDS_SIMULATION_DATA.h>
#include <Common_Tools/Math_Tools/RANGE_ITERATOR.h>
#include <PhysBAM_Geometry/Geometry_Particles/REGISTER_GEOMETRY_READ_WRITE.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <PhysBAM_Tools/Data_Structures/HASHTABLE_ITERATOR.h>
#include <SPGrid/Tools/SPGrid_Block_Iterator.h>
#include <SPGrid_Fluids/Simulation/FLUIDS_SIMULATION_FLAGS.h>
#include <PhysBAM_Geometry/Read_Write/Geometry/READ_WRITE_TRIANGULATED_SURFACE.h>

using namespace PhysBAM;
//#####################################################################
// Visualize_Triangulated_Surface
//#####################################################################
template<class T_STRUCT,class T,int d> void VISUALIZE_TRIANGULATED_SURFACE<T_STRUCT,T,d>::
Visualize_Triangulated_Surface(STREAM_TYPE stream_type,T_HIERARCHY& hierarchy,const std::string output_directory,const int axis,const T h,const bool flat,const T scale)
{
    typedef VECTOR<T,3> TV_3;
    typedef VECTOR<int,3> TV_INT_3;
    ARRAY<TV_INT> unit_box_di;
    ARRAY<VECTOR<int,2> > connected_cell_edges;ARRAY<VECTOR<int,3> > triangles;
    unit_box_di.Append(TV_INT(0,0)); // mannually add these to maintain convention for triangle normals
    unit_box_di.Append(TV_INT(0,1));
    unit_box_di.Append(TV_INT(1,0));
    unit_box_di.Append(TV_INT(1,1));
    connected_cell_edges.Append(VECTOR<int,2>(1,3));
    connected_cell_edges.Append(VECTOR<int,2>(3,4));
    connected_cell_edges.Append(VECTOR<int,2>(4,2));
    connected_cell_edges.Append(VECTOR<int,2>(2,1));
    triangles.Append(VECTOR<int,3>(1,3,2));
    triangles.Append(VECTOR<int,3>(3,4,2));

    Initialize_Geometry_Particle();
    Initialize_Read_Write_Structures();
    const int levels=hierarchy.Levels();

    GEOMETRY_PARTICLES<TV_3> particles;
    TRIANGULATED_SURFACE<T>* total_surface;
    ARRAY<TRIANGULATED_SURFACE<T>* > interior_cell_surface(levels);
    total_surface=TRIANGULATED_SURFACE<T>::Create(particles);
    for(int level=1;level<=levels;level++) interior_cell_surface(level)=TRIANGULATED_SURFACE<T>::Create(particles);
    // for flat checkered tris
    TRIANGULATED_SURFACE<T>* total_surface_alternate;
    ARRAY<TRIANGULATED_SURFACE<T>* > interior_cell_surface_alternate(levels);
    total_surface_alternate=TRIANGULATED_SURFACE<T>::Create(particles);
    for(int level=1;level<=levels;level++) interior_cell_surface_alternate(level)=TRIANGULATED_SURFACE<T>::Create(particles);

    // cell hash
    HASHTABLE<TV_INT,CELL_INFO> cell_hash;
    for(int level=1;level<=levels;level++){const int cell_size=1<<(level-1);
        Const_flag_array_type flags=hierarchy.Allocator(level).Get_Const_Array(&T_STRUCT::flags);
        for(SPGrid_Block_Iterator<Flag_array_mask> iterator(hierarchy.Blocks(level));iterator.Valid();iterator.Next())
            if(iterator.Data(flags)&SPGrid_Cell_Type_Interior){
                const TV_INT native_index=iterator.Index().template Cast<TV_INT>();
                const TV_INT refined_index=(native_index-1)*2*cell_size+cell_size+1; // refined_grid based resolution
                CELL_INFO info(cell_hash.Size()+1,level,native_index);
                cell_hash.Insert(refined_index,info);}}
    // node hash
    typedef CELL_INFO NODE_INFO;
    HASHTABLE<TV_INT,NODE_INFO> node_hash; // maps nodes id, 4 neighbors
    HASHTABLE<TV,NODE_INFO> active_node_hash;
    int particle_offset=0;
    for(HASHTABLE_ITERATOR<TV_INT,CELL_INFO> iterator(cell_hash);iterator.Valid();iterator.Next()){    
        const int level=iterator.Data().y;
        const TV_INT native_index=iterator.Data().z;
        const int cell_size=1<<(level-1);
        const unsigned& flag=hierarchy.Array(level,&T_STRUCT::flags)(std_array<int,d>(native_index));
        if(flag&SPGrid_Cell_Type_Interior){
            ARRAY<TV_INT> native_indices;
            for(int i=1;i<=unit_box_di.m;i++) native_indices.Append(native_index+unit_box_di(i));
            const TV_INT base_fine_index=(native_index-1)*2*cell_size+1;
            for(int i=1;i<=unit_box_di.m;i++){const TV_INT fine_index=base_fine_index+(2*cell_size)*unit_box_di(i);
                node_hash.Get_Or_Insert(fine_index,NODE_INFO(node_hash.Size()+1,level,native_indices(i)));}}
        if(flag&SPGrid_Node_Active){const TV_INT fine_index=(native_index-1)*2*cell_size+1;
            active_node_hash.Get_Or_Insert(TV(fine_index),NODE_INFO(active_node_hash.Size()+1,level,native_index));}}
    // Create cell centered particles
    if(!flat){particles.array_collection->Add_Elements(cell_hash.Size());
        for(HASHTABLE_ITERATOR<TV_INT,CELL_INFO> iterator(cell_hash);iterator.Valid();iterator.Next()){
            const int p=particle_offset+iterator.Data().x;
            const int level=iterator.Data().y;
            const TV_INT& index=iterator.Data().z;
            TV_3 X=hierarchy.Grid(level).Center(index).Insert(h,axis);
            X(axis)+=scale*hierarchy.Grid(level).dX.x*(T).5;
            particles.X(p)=X;}}
    // Create node particles
    particle_offset=particles.array_collection->Size();
    particles.array_collection->Add_Elements(node_hash.Size());
    for(HASHTABLE_ITERATOR<TV_INT,NODE_INFO> iterator(node_hash);iterator.Valid();iterator.Next()){    
        const int p=particle_offset+iterator.Data().x;
        const int level=iterator.Data().y;
        const TV_INT& index=iterator.Data().z;
        particles.X(p)=hierarchy.Grid(level).Node(index).Insert(h,axis);}
    // Create triangulated surface for interior cells
    for(HASHTABLE_ITERATOR<TV_INT, CELL_INFO> iterator(cell_hash);iterator.Valid();iterator.Next()){
        const int cell_center_index=iterator.Data().x;
        const int level=iterator.Data().y;
        const TV_INT native_index=iterator.Data().z;
        const int cell_size=1<<(level-1);
        const unsigned& flag=hierarchy.Array(level,&T_STRUCT::flags)(std_array<unsigned int,d>(native_index));
        if(flag&SPGrid_Cell_Type_Interior){
            const TV_INT base_fine_index=(native_index-1)*2*cell_size+1;
            ARRAY<int> p;
            for(int i=1;i<=unit_box_di.m;i++){const TV_INT fine_index=base_fine_index+(2*cell_size)*unit_box_di(i);
                p.Append(node_hash.Get(fine_index).x+particle_offset);}
            if(!flat){
                for(int i=1;i<=connected_cell_edges.m;i++){const VECTOR<int,2>& edge=connected_cell_edges(i);
                    interior_cell_surface(level)->mesh.elements.Append(VECTOR<int,3>(p(edge(1)),p(edge(2)),cell_center_index));
                    total_surface->mesh.elements.Append(VECTOR<int,3>(p(edge(1)),p(edge(2)),cell_center_index));}}
            else{
                if(native_index.Sum()%2){
                    for(int i=1;i<=triangles.m;i++){const VECTOR<int,3>& t=triangles(i);
                        interior_cell_surface(level)->mesh.elements.Append(VECTOR<int,3>(p(t(1)),p(t(2)),p(t(3))));
                        total_surface->mesh.elements.Append(VECTOR<int,3>(p(t(1)),p(t(2)),p(t(3))));}}
                else{
                    for(int i=1;i<=triangles.m;i++){const VECTOR<int,3>& t=triangles(i);
                        interior_cell_surface_alternate(level)->mesh.elements.Append(VECTOR<int,3>(p(t(1)),p(t(2)),p(t(3))));
                        total_surface_alternate->mesh.elements.Append(VECTOR<int,3>(p(t(1)),p(t(2)),p(t(3))));}}
            }
        }
    }
    // Write output
    for(int level=1;level<=levels;level++) interior_cell_surface(level)->Update_Number_Nodes();
    total_surface->Update_Number_Nodes();
    FILE_UTILITIES::Create_Directory(output_directory);
    for(int level=1;level<=levels;level++){const std::string l=STRING_UTILITIES::string_sprintf("%d",level);
        FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/grid_"+l+".tri",*interior_cell_surface(level));}
    FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/grid.tri",*total_surface);
    if(flat){
        for(int level=1;level<=levels;level++) interior_cell_surface_alternate(level)->Update_Number_Nodes();
        total_surface_alternate->Update_Number_Nodes();
        for(int level=1;level<=levels;level++){const std::string l=STRING_UTILITIES::string_sprintf("%d",level);
            FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/grid_"+l+"_alternate.tri",*interior_cell_surface_alternate(level));}
        FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/grid_alternate.tri",*total_surface_alternate);}
}
//#####################################################################
template class VISUALIZE_TRIANGULATED_SURFACE<FLUIDS_SIMULATION_DATA<float>,float,2>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class VISUALIZE_TRIANGULATED_SURFACE<FLUIDS_SIMULATION_DATA<double>,double,2>;
#endif
