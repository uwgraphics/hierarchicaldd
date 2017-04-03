//#####################################################################
// Copyright (c) 2014, Raj Setaluri, Mridul Aanjneya
// Distributed under the FreeBSD license (see license.txt)
//#####################################################################
#include <SPGrid_Fluids/Visualization/VISUALIZE_HEIGHTFIELD.h>
#include <SPGrid/Tools/SPGrid_Block_Iterator.h>
#include <SPGrid_Fluids/Simulation/FLUIDS_SIMULATION_DATA.h>
#include <SPGrid_Fluids/Simulation/FLUIDS_SIMULATION_FLAGS.h>
#include <PhysBAM_Tools/Read_Write/Grids_Uniform_Arrays/READ_WRITE_ARRAYS.h>
#include <PhysBAM_Tools/Data_Structures/HASHTABLE_ITERATOR.h>
#include <PhysBAM_Tools/Data_Structures/TRIPLE.h>
#include <PhysBAM_Geometry/Geometry_Particles/REGISTER_GEOMETRY_READ_WRITE.h>
#include <PhysBAM_Geometry/Solids_Geometry/DEFORMABLE_GEOMETRY_COLLECTION.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/FREE_PARTICLES.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/SEGMENTED_CURVE.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <PhysBAM_Geometry/Read_Write/Geometry/READ_WRITE_TRIANGULATED_SURFACE.h>

using namespace PhysBAM;
//#####################################################################
// Visualize_Heightfield
//#####################################################################
template<class T_STRUCT,class T,int d> void VISUALIZE_HEIGHTFIELD<T_STRUCT,T,d>::
Visualize_Heightfield(STREAM_TYPE stream_type,T_HIERARCHY& hierarchy,T T_STRUCT::* data_channel,const std::string output_directory,const int frame,const T scale)
{
    typedef VECTOR<T,3> TV_3;

    Initialize_Geometry_Particle();Initialize_Read_Write_Structures();

    GEOMETRY_PARTICLES<TV_3> particles;

    ARRAY<FREE_PARTICLES<TV_3>*> interior_cell_centers(hierarchy.Levels());
    ARRAY<FREE_PARTICLES<TV_3>*> dirichlet_cell_centers(hierarchy.Levels());
    ARRAY<FREE_PARTICLES<TV_3>*> ghost_cell_centers(hierarchy.Levels());
    ARRAY<FREE_PARTICLES<TV_3>*> node_particles(hierarchy.Levels());
    ARRAY<SEGMENTED_CURVE<TV_3>*> interior_cell_curve(hierarchy.Levels());
    TRIANGULATED_SURFACE<T>* value_surface;
    SEGMENTED_CURVE<TV_3>* dual_graph_curve;
    
    for(int level=1;level<=hierarchy.Levels();level++){
        interior_cell_centers(level) =FREE_PARTICLES<TV_3>::Create(particles);
        dirichlet_cell_centers(level)=FREE_PARTICLES<TV_3>::Create(particles);
        ghost_cell_centers(level)    =FREE_PARTICLES<TV_3>::Create(particles);
        node_particles(level)        =FREE_PARTICLES<TV_3>::Create(particles);
        interior_cell_curve(level)   =SEGMENTED_CURVE<TV_3>::Create(particles);}
    value_surface=TRIANGULATED_SURFACE<T>::Create(particles);
    dual_graph_curve=SEGMENTED_CURVE<TV_3>::Create(particles);

    // Create cell hash and map to non-ghost parents
    typedef TRIPLE<int,int,TV_INT> CELL_INFO;
    HASHTABLE<TV_INT, CELL_INFO> cell_hash; 
    ARRAY<int> non_ghost_parent_of_cell;
    for(int level=hierarchy.Levels();level>=1;level--){
        Flag_array_type flags=hierarchy.Allocator(level).Get_Array(&T_STRUCT::flags);
        const int cell_size=1<<(level-1);
        for(SPGrid_Block_Iterator<typename Flag_array_type::MASK> iterator(hierarchy.Blocks(level));iterator.Valid();iterator.Next()){
            if(iterator.Data(flags) & (SPGrid_Cell_Type_Interior|SPGrid_Cell_Type_Ghost|SPGrid_Cell_Type_Dirichlet)){
                // Iterate over all active cells, add to cell hash
                TV_INT native_index_2d=iterator.Index().template Cast<TV_INT>();
                TV_INT refined_index_2d=(native_index_2d-1)*2*cell_size + cell_size+1;  // refined_grid based resolution                
                // Cell hashing
                CELL_INFO info(cell_hash.Size()+1,level,native_index_2d);
                cell_hash.Insert(refined_index_2d,info);
                // Find parents of each cell...
                if(iterator.Data(flags) & SPGrid_Cell_Type_Ghost){
                    // T_RAW_INDEX cur_coord(native_index_2d);
                    TV_INT current_native_index_2d=native_index_2d;
                    int parent_level;
                    for(parent_level=level+1; parent_level<=hierarchy.Levels();parent_level++){
                        current_native_index_2d=(current_native_index_2d-1)/2+1;
                        if(hierarchy.Set(parent_level).Is_Set(std_array<int,d>(current_native_index_2d), SPGrid_Cell_Type_Interior)) break;}
                    if(parent_level<=hierarchy.Levels()){
                        TV_INT current_refined_index_2d=(current_native_index_2d-1)*(1<<parent_level)+(1<<(parent_level-1))+1;
                        non_ghost_parent_of_cell.Append(cell_hash.Get(current_refined_index_2d).x);}
                    else non_ghost_parent_of_cell.Append(0);}
                else non_ghost_parent_of_cell.Append(info.x);}}}

    PHYSBAM_ASSERT(non_ghost_parent_of_cell.Size()==cell_hash.Size());
    
    // Fill in node_hash
    typedef CELL_INFO NODE_INFO;
    HASHTABLE<TV_INT,NODE_INFO> node_hash; // maps nodes id, 4 neighbors
    for(HASHTABLE_ITERATOR<TV_INT, CELL_INFO> iterator(cell_hash);iterator.Valid();iterator.Next()){    
        const int level=iterator.Data().y;
        const TV_INT native_index_2d=iterator.Data().z;
        const int cell_size=1<<(level-1);
        Flag_array_type flags=hierarchy.Allocator(level).Get_Array(&T_STRUCT::flags);
        unsigned cell_type=flags(std_array<unsigned int,d>(native_index_2d));
        if(cell_type & (SPGrid_Cell_Type_Interior|SPGrid_Cell_Type_Dirichlet)){
            // Node hashing @ doubly refined resolution
            TV_INT native_index_00=native_index_2d;
            TV_INT native_index_10=native_index_00+TV_INT(1,0);
            TV_INT native_index_11=native_index_00+TV_INT(1,1);
            TV_INT native_index_01=native_index_00+TV_INT(0,1);                
            TV_INT index_00=(native_index_2d-1)*2*cell_size+1;
            TV_INT index_10=index_00+TV_INT(2*cell_size,0);
            TV_INT index_11=index_00+TV_INT(2*cell_size,2*cell_size);
            TV_INT index_01=index_00+TV_INT(0,2*cell_size);
            node_hash.Get_Or_Insert(index_00,NODE_INFO(node_hash.Size()+1, level, native_index_00));
            node_hash.Get_Or_Insert(index_01,NODE_INFO(node_hash.Size()+1, level, native_index_01));
            node_hash.Get_Or_Insert(index_11,NODE_INFO(node_hash.Size()+1, level, native_index_11));
            node_hash.Get_Or_Insert(index_10,NODE_INFO(node_hash.Size()+1, level, native_index_10));}}

    // Fill in node_neighbors    
    typedef VECTOR<int,4> NODE_NEIGHBORS;
    ARRAY<NODE_NEIGHBORS > node_neighbors(node_hash.Size());
    for(HASHTABLE_ITERATOR<TV_INT, CELL_INFO> iterator(cell_hash);iterator.Valid();iterator.Next()){    
        const int cell_id=iterator.Data().x;
        const int level=iterator.Data().y;
        const TV_INT native_index_2d=iterator.Data().z;
        const int cell_size=1<<(level-1);
        Flag_array_type flags=hierarchy.Allocator(level).Get_Array(&T_STRUCT::flags);
        unsigned cell_type=flags(std_array<unsigned int,d>(native_index_2d));
        if(cell_type & (SPGrid_Cell_Type_Interior|SPGrid_Cell_Type_Dirichlet|SPGrid_Cell_Type_Ghost)){
            // Node hashing @ doubly refined resolution
            TV_INT index_00=(native_index_2d-1)*2*cell_size+1;
            TV_INT index_01=index_00+TV_INT(0,2*cell_size);
            TV_INT index_10=index_00+TV_INT(2*cell_size,0);
            TV_INT index_11=index_00+TV_INT(2*cell_size,2*cell_size);
            int n_00=node_hash.Get_Default(index_00,NODE_INFO(0,0,TV_INT())).x;
            int n_01=node_hash.Get_Default(index_01,NODE_INFO(0,0,TV_INT())).x;
            int n_10=node_hash.Get_Default(index_10,NODE_INFO(0,0,TV_INT())).x;
            int n_11=node_hash.Get_Default(index_11,NODE_INFO(0,0,TV_INT())).x;
            if(n_00) node_neighbors(n_00)(3)=non_ghost_parent_of_cell(cell_id);
            if(n_01) node_neighbors(n_01)(2)=non_ghost_parent_of_cell(cell_id);
            if(n_10) node_neighbors(n_10)(4)=non_ghost_parent_of_cell(cell_id);
            if(n_11) node_neighbors(n_11)(1)=non_ghost_parent_of_cell(cell_id);}}

    // Cell-centered particles on domain plane
    {int particle_offset=particles.array_collection->Size();
    particles.array_collection->Add_Elements(cell_hash.Size());    
    // Plot particles at cell centers
    for(HASHTABLE_ITERATOR<TV_INT, CELL_INFO> iterator(cell_hash);iterator.Valid();iterator.Next()){    
        const int cell_id=iterator.Data().x;
        const int p=cell_id + particle_offset;
        const int level=iterator.Data().y;
        const TV_INT native_index_2d=iterator.Data().z;
        Flag_array_type flags=hierarchy.Allocator(level).Get_Array(&T_STRUCT::flags);
        const unsigned cell_type=flags(std_array<unsigned int,d>(native_index_2d));        
        const TV X=hierarchy.Grid(level).Center(native_index_2d);
        particles.X(p)=X.Insert(0,2);        
        if(cell_type & SPGrid_Cell_Type_Interior)
            interior_cell_centers(level)->nodes.Append(p);
        else if(cell_type & SPGrid_Cell_Type_Dirichlet)
            dirichlet_cell_centers(level)->nodes.Append(p);
        else if(cell_type & SPGrid_Cell_Type_Ghost)
            ghost_cell_centers(level)->nodes.Append(p);}}

    // Cell centered particles on value plot
    {const int particle_offset=particles.array_collection->Size();
    particles.array_collection->Add_Elements(cell_hash.Size());
    
    // Plot particles at cell centers
    for(HASHTABLE_ITERATOR<TV_INT, CELL_INFO> iterator(cell_hash);iterator.Valid();iterator.Next()){    
        const int cell_id=iterator.Data().x;
        const int p=cell_id + particle_offset;
        const int level=iterator.Data().y;
        const TV_INT native_index_2d=iterator.Data().z;        
        const TV X=hierarchy.Grid(level).Center(native_index_2d);
        Data_array_type values=hierarchy.Allocator(level).Get_Array(data_channel);
        particles.X(p)=X.Insert((T)scale*values(std_array<unsigned int,d>(native_index_2d)),2);}    
    for(int node_id=1;node_id<=node_neighbors.m;node_id++){
        ARRAY<int> neighbor_particles;
        for(int i=1;i<=4;i++){
            const int cell_id=node_neighbors(node_id)(i);
            if(cell_id) neighbor_particles.Append_Unique(cell_id+particle_offset);}
        switch(neighbor_particles.m){
        case 3:
            value_surface->mesh.elements.Append(VECTOR<int,3>(neighbor_particles));
            break;
        case 4:
            value_surface->mesh.elements.Append(VECTOR<int,3>(neighbor_particles(1),neighbor_particles(2),neighbor_particles(3)));
            value_surface->mesh.elements.Append(VECTOR<int,3>(neighbor_particles(1),neighbor_particles(3),neighbor_particles(4)));
            break;
        default:;}
        if(neighbor_particles.m<3) continue;
        for(int i=1;i<=neighbor_particles.m;i++) dual_graph_curve->mesh.elements.Append(TV_INT(neighbor_particles(i),neighbor_particles(i%neighbor_particles.m+1)));}}

    // Node particles on domain plane
    {int particle_offset=particles.array_collection->Size();
    particles.array_collection->Add_Elements(node_hash.Size());    
    // Create node particles
    for(HASHTABLE_ITERATOR<TV_INT, NODE_INFO> iterator(node_hash);iterator.Valid();iterator.Next()){
        const int node_id=iterator.Data().x;
        const int p=node_id + particle_offset;
        const int level=iterator.Data().y;
        const TV_INT native_index_2d=iterator.Data().z;        
        const TV X=hierarchy.Grid(level).Node(native_index_2d);
        particles.X(p)=X.Insert(0,2);        
        node_particles(level)->nodes.Append(p);}
    // Create segment lattice for interior cells    
    for(HASHTABLE_ITERATOR<TV_INT, CELL_INFO> iterator(cell_hash);iterator.Valid();iterator.Next()){
        const int level=iterator.Data().y;
        const TV_INT native_index_2d=iterator.Data().z;
        const int cell_size=1<<(level-1);        
        Flag_array_type flags=hierarchy.Allocator(level).Get_Array(&T_STRUCT::flags);
        const unsigned cell_type=flags(std_array<unsigned int,d>(native_index_2d));        
        if(cell_type & (SPGrid_Cell_Type_Interior|SPGrid_Cell_Type_Dirichlet)){
            TV_INT refined_index_00=(native_index_2d-1)*2*cell_size+1;
            TV_INT refined_index_10=refined_index_00+TV_INT(2*cell_size,0);
            TV_INT refined_index_11=refined_index_00+TV_INT(2*cell_size,2*cell_size);
            TV_INT refined_index_01=refined_index_00+TV_INT(0,2*cell_size);            
            int p_00=node_hash.Get(refined_index_00).x+particle_offset;
            int p_01=node_hash.Get(refined_index_01).x+particle_offset;
            int p_11=node_hash.Get(refined_index_11).x+particle_offset;
            int p_10=node_hash.Get(refined_index_10).x+particle_offset;            
            interior_cell_curve(level)->mesh.elements.Append(TV_INT(p_00,p_01));
            interior_cell_curve(level)->mesh.elements.Append(TV_INT(p_01,p_11));
            interior_cell_curve(level)->mesh.elements.Append(TV_INT(p_11,p_10));
            interior_cell_curve(level)->mesh.elements.Append(TV_INT(p_10,p_00));}}}

    for(int level=1;level<=hierarchy.Levels();level++){
        interior_cell_centers(level)->Update_Number_Nodes();
        dirichlet_cell_centers(level)->Update_Number_Nodes();
        ghost_cell_centers(level)->Update_Number_Nodes();
        node_particles(level)->Update_Number_Nodes();
        interior_cell_curve(level)->Update_Number_Nodes();}
    dual_graph_curve->Update_Number_Nodes();
    value_surface->Update_Number_Nodes();

    DEFORMABLE_GEOMETRY_COLLECTION<TV_3> collection(particles);

    for(int level=1;level<=hierarchy.Levels();level++){
        // collection.Add_Structure(interior_cell_centers(level));
        delete interior_cell_centers(level);

        // collection.Add_Structure(dirichlet_cell_centers(level));
        delete dirichlet_cell_centers(level);

        // collection.Add_Structure(ghost_cell_centers(level));
        delete ghost_cell_centers(level);

        // collection.Add_Structure(node_particles(level));
        delete node_particles(level);

        // collection.Add_Structure(interior_cell_curve(level));
        delete interior_cell_curve(level);}

    collection.Add_Structure(dual_graph_curve);
    // delete dual_graph_curve;

    collection.Add_Structure(value_surface);
    // delete value_surface;

    // Write PhysBAM-style output
    FILE_UTILITIES::Create_Directory(output_directory+"/"+STRING_UTILITIES::Value_To_String(frame));
    collection.Write(stream_type,output_directory,frame,frame,true);
    FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+STRING_UTILITIES::Value_To_String(frame)+"/heightfield.tri",*value_surface);

    // generate uv coordinates
    {typedef float RW;
    const int backward_compatible=0;
    ARRAY<VECTOR<T,2> > texture_coordinates;
    ARRAY<VECTOR<int,3> > triangle_texture_coordinates;
    for(int p=1;p<=value_surface->particles.array_collection->Size();p++)
        texture_coordinates.Append(VECTOR<T,2>(value_surface->particles.X(p).x,value_surface->particles.X(p).z));
    for(int i=1;i<=value_surface->mesh.elements.m;++i){
        int index1,index2,index3;value_surface->mesh.elements(i).Get(index1,index2,index3);
        triangle_texture_coordinates.Append(VECTOR<int,3>(index1,index2,index3));}
    FILE_UTILITIES::Write_To_File<RW>(output_directory+"/"+STRING_UTILITIES::Value_To_String(frame)+"/heightfield.uv",texture_coordinates,backward_compatible,triangle_texture_coordinates);}
}
//#####################################################################
template class VISUALIZE_HEIGHTFIELD<FLUIDS_SIMULATION_DATA<float>,float,2>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class VISUALIZE_HEIGHTFIELD<FLUIDS_SIMULATION_DATA<double>,double,2>;
#endif
