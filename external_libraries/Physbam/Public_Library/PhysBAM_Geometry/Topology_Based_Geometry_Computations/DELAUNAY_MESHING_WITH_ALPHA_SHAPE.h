//#####################################################################
// Copyright 2012, Wen Zheng.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __DELAUNAY_MESHING_WITH_ALPHA_SHAPE__
#define __DELAUNAY_MESHING_WITH_ALPHA_SHAPE__

#include <PhysBAM_Tools/Arrays/INDIRECT_ARRAY.h>
#include <PhysBAM_Tools/Arrays_Computations/ARRAY_COPY.h>
#include <PhysBAM_Tools/Arrays_Computations/SORT.h>
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_NODE.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_EDGE.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Data_Structures/QUEUE.h>
#include <PhysBAM_Tools/Data_Structures/HASHTABLE.h>
#include <PhysBAM_Tools/Log/DEBUG_SUBSTEPS.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_AREA.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TETRAHEDRALIZED_VOLUME.h>
#include <PhysBAM_Geometry/Topology/SEGMENT_MESH.h>
#include <PhysBAM_Geometry/Topology/TRIANGLE_MESH.h>
#include <PhysBAM_Geometry/Basic_Geometry/BOX.h>
#include <PhysBAM_Geometry/Basic_Geometry/TRIANGLE_2D.h>
#include <PhysBAM_Geometry/Basic_Geometry/TRIANGLE_3D.h>
#include <PhysBAM_Geometry/Basic_Geometry/SEGMENT_2D.h>
#include <PhysBAM_Geometry/Basic_Geometry/TETRAHEDRON.h>
#include <PhysBAM_Geometry/Basic_Geometry_Intersections/SEGMENT_2D_SEGMENT_2D_INTERSECTION.h>
#include <PhysBAM_Geometry/Basic_Geometry_Intersections/TRIANGLE_3D_TRIANGLE_3D_INTERSECTION.h>
#include <PhysBAM_Geometry/Basic_Geometry_Intersections/SEGMENT_3D_TRIANGLE_3D_INTERSECTION.h>
#include <PhysBAM_Geometry/Basic_Geometry_Intersections/BOX_SEGMENT_2D_INTERSECTION.h>
#include <PhysBAM_Geometry/Basic_Geometry_Intersections/BOX_TRIANGLE_3D_INTERSECTION.h>
#include <PhysBAM_Tools/Random_Numbers/RANDOM_NUMBERS.h>
#include <PhysBAM_Tools/Arrays_Computations/SORT.h>
#include <PhysBAM_Tools/Data_Structures/QUEUE.h>
#include <PhysBAM_Tools/Parallel_Computation/MPI_UNIFORM_GRID.h>
#include <PhysBAM_Tools/Parallel_Computation/MPI_PARTICLES.h>
#include <PhysBAM_Geometry/Geometry_Particles/GEOMETRY_PARTICLES.h>

namespace ADAPTIVE_EXACT_PREDICATES{
double insphere(double *pa,double *pb,double *pc,double *pd,double *pe);
double incircle(double *pa,double *pb,double *pc,double *pe);
double orient3d(double *pa,double *pb,double *pc,double *pd);
double orient2d(double *pa,double *pb,double *pc);
void exactinit();
};

namespace PhysBAM
{
template<class TV> class INSIDE_DOMAIN_PREDICATOR
{
public:
    virtual bool operator()(const TV& x,const typename TV::SCALAR thickness=0) const {return true;}
};

template<class TV> class DELAUNAY_MESHING_WITH_ALPHA_SHAPE;

// 2D version
template<class T> 
class DELAUNAY_MESHING_WITH_ALPHA_SHAPE<VECTOR<T,2> >
{
    typedef VECTOR<T,2> TV;
    typedef VECTOR<int,2> TV_INT;
    typedef TRIANGLE_2D<T> T_ELEMENT;

public:
    const GRID<TV>& grid;
    GEOMETRY_PARTICLES<TV>& particles;
    ARRAY<ARRAY<int>,TV_INT> cell_particles;
    ARRAY<bool> particle_finished;
    ARRAY<int> particle_counts;
    QUEUE<TV_INT> frontier;
    HASHTABLE<TV_INT,int> face_map;
    int ghost_cells,int_radius;
    T grid_radius,radius,radius_squared;
    const ARRAY<int>* valid_indices;
    const ARRAY<bool,TV_INT>* blocked_cells;
    const SEGMENT_MESH* initial_faces;
    ARRAY<ARRAY<TV_INT>,TV_INT> cell_to_faces;
    bool use_alpha_shape,use_manifold_mesh;
    MPI_PARTICLES<GRID<TV> >* mpi_particles,*mpi_particles2;
    bool clean_memory_after_meshing;

    DELAUNAY_MESHING_WITH_ALPHA_SHAPE(const GRID<TV>& grid_input,GEOMETRY_PARTICLES<TV>& particles_input)
        :grid(grid_input),particles(particles_input),frontier(1),grid_radius((T)0.88),valid_indices(0),blocked_cells(0),initial_faces(0),use_alpha_shape(true),use_manifold_mesh(true),mpi_particles(0),mpi_particles2(0),clean_memory_after_meshing(false)
    {
        ADAPTIVE_EXACT_PREDICATES::exactinit();
    }

    ~DELAUNAY_MESHING_WITH_ALPHA_SHAPE(){}

    void Initialize_Grid(int ghost_cells_input=3)
    {
        ghost_cells=ghost_cells_input;
        cell_particles.Resize(grid.Domain_Indices(ghost_cells));
        cell_to_faces.Resize(grid.Domain_Indices(ghost_cells));
    }

    void Construct_Mesh(TRIANGLE_MESH& mesh)
    {
        LOG::SCOPE scope("Construct Mesh");
        Initialize_Meshing_Parameters();
        Build_Grid_Indexing();
        mesh.elements.Clean_Memory();
        Clear_Auxilary_Data();
        while(Initialize_Front_With_A_Single_Triangle(mesh)) Advancing_Front(mesh);
        if(blocked_cells) Clear_Intersection_With_Blocked_Cells(grid,*blocked_cells,mesh);
    }

    void Construct_Constrained_Mesh(TRIANGLE_MESH& mesh,SEGMENT_MESH& boundary_mesh,bool skip_dangenous=false,INSIDE_DOMAIN_PREDICATOR<TV>* domain=0,INSIDE_DOMAIN_PREDICATOR<TV>* mpi_local_domain=0)
    {
        LOG::SCOPE scope("Construct Constrained Mesh");
        // index boundary elements
        LOG::Time("Index boundary circumspheres");
        ARRAY<ARRAY<int>,TV_INT> cell_elements(grid.Domain_Indices(ghost_cells));
        ARRAY<TV> element_circumcener(boundary_mesh.elements.m);
        ARRAY<T> element_circumradius(boundary_mesh.elements.m);
        ARRAY<bool> is_boundary_nodes(particles.X.m);
        ARRAYS_COMPUTATIONS::Fill(is_boundary_nodes,false);
        for(int i=1;i<=boundary_mesh.elements.m;i++){
            const TV_INT& e=boundary_mesh.elements(i);
            for(int j=1;j<=TV::m;j++) is_boundary_nodes(e(j))=true;
            element_circumcener(i)=(T)0.5*(particles.X(e.x)+particles.X(e.y));
            element_circumradius(i)=(element_circumcener(i)-particles.X(e.x)).Magnitude();
            TV_INT cell=Clamp_Cell(element_circumcener(i));
            int range=(int)ceil(element_circumradius(i)/grid.dX.Max());
            RANGE<TV_INT> cell_range(Clamp_Cell(cell-range),Clamp_Cell(cell+range));
            for(UNIFORM_GRID_ITERATOR_CELL<TV> iter(grid,cell_range);iter.Valid();iter.Next()) cell_elements(iter.Cell_Index()).Append(i);}
        // filter out dangerous particles
        LOG::Time("Filter out dangerous particles");
        ARRAY<int> full_indices,safe_indices,dangerous_indices;
        if(valid_indices) full_indices=*valid_indices;
        else{full_indices.Resize(particles.X.m);for(int i=1;i<=particles.X.m;i++) full_indices(i)=i;}
        for(int i=1;i<=full_indices.m;i++){
            int p=full_indices(i);
            if(is_boundary_nodes(p)) safe_indices.Append(p);
            else{
                TV_INT cell=Clamp_Cell(particles.X(p));
                const ARRAY<int>& list=cell_elements(cell);
                bool dangerous=false;
                for(int j=1;j<=list.m;j++)if((particles.X(p)-element_circumcener(list(j))).Magnitude()<element_circumradius(list(j))+1e-8*grid.dX.Min()){dangerous=true;break;}
                if(dangerous) dangerous_indices.Append(p);
                else safe_indices.Append(p);}}
        // meshing without dangerous particles
        Set_Valid_Indices(&safe_indices);
        Construct_Mesh(mesh);
        PHYSBAM_DEBUG_WRITE_SUBSTEP("after constructing mesh for step 2",0,2);
        if(domain){
            ARRAY<int> deletion_list;
            for(int i=1;i<=mesh.elements.m;i++){
                const VECTOR<int,TV::m+1>& element=mesh.elements(i);
                TV barycenter=T_ELEMENT::Center(particles.X.Subset(element));
                if(!(*domain)(barycenter) || (mpi_local_domain && !(*mpi_local_domain)(barycenter))){deletion_list.Append(i);continue;}
                TV circumcenter=T_ELEMENT::Circumcenter(particles.X.Subset(element));
                if(!(*domain)(circumcenter))deletion_list.Append(i);}
            Sort(deletion_list);
            int last=mesh.elements.m;
            for(int k=deletion_list.m;k>=1;k--) mesh.elements(deletion_list(k))=mesh.elements(last--);
            mesh.elements.Resize(last);}
        PHYSBAM_DEBUG_WRITE_SUBSTEP("after deleting outside elements for step 2",0,2);
        if(skip_dangenous) return;
        // add in dangerous particles
        LOG::Time("Index all elements");
        Set_Valid_Indices(&full_indices);
        // bin all elements for point location and construct face-to-element reference
        for(int i=1;i<=cell_elements.array.m;i++) cell_elements.array(i).Resize(0);
        HASHTABLE<TV_INT,VECTOR<int,2> > face_to_elements;
        for(int i=1;i<=mesh.elements.m;i++){
            const VECTOR<int,TV::m+1>& e=mesh.elements(i);
            RANGE<TV> bounding_box=RANGE<TV>::Bounding_Box(particles.X.Subset(e));
            RANGE<TV_INT> cell_range(Clamp_Cell(bounding_box.min_corner),Clamp_Cell(bounding_box.max_corner));
            for(UNIFORM_GRID_ITERATOR_CELL<TV> iter(grid,cell_range);iter.Valid();iter.Next()) cell_elements(iter.Cell_Index()).Append(i);
            VECTOR<TV_INT,TV::m+1> tet_faces=Get_Segments_In_Triangle(e);
            for(int j=1;j<=TV::m+1;j++){
                TV_INT face=tet_faces(j);
                bool sign=Order_Face_With_Sign(face);
                int index=sign?1:2;
                if(face_to_elements.Contains(face)){
                    VECTOR<int,2> value=face_to_elements.Get(face);
                    PHYSBAM_ASSERT(value(index)==0);
                    value(index)=i;
                    face_to_elements.Set(face,value);}
                else{
                    VECTOR<int,2> value;value(index)=i;
                    face_to_elements.Set(face,value);}}}
        // incrementally insert particles
        LOG::Time("Add in dangerous particles");
        ARRAY<int> deletion_list;
        ARRAY<bool> element_deleted(mesh.elements.m);
        ARRAYS_COMPUTATIONS::Fill(element_deleted,false);
        for(int i=1;i<=dangerous_indices.m;i++){
            int p=dangerous_indices(i);
            // locate
            TV_INT cell=Clamp_Cell(particles.X(p));
            const ARRAY<int>& list=cell_elements(cell);
            int initial_element_index=0;
            for(int j=1;j<=list.m;j++)if(!element_deleted(list(j))){
                TRIANGLE_2D<T> element(particles.X.Subset(mesh.elements(list(j))));
                if(element.Inside(particles.X(p))){
                    initial_element_index=list(j);
                    break;}}
            if(initial_element_index==0) continue;// the point is outside
            // expand (depth-first search)
            ARRAY<int> stack;
            ARRAY<TV_INT> delayed_face_list;
            ARRAY<TV_INT> face_list;
            stack.Append(initial_element_index);
            element_deleted(initial_element_index)=true;
            deletion_list.Append(initial_element_index);
            while(stack.m>0){
                int element_index=stack.Pop();
                const VECTOR<int,TV::m+1>& current_element=mesh.elements(element_index);
                VECTOR<TV_INT,TV::m+1> tet_faces=Get_Segments_In_Triangle(current_element);
                for(int j=1;j<=TV::m+1;j++){
                    TV_INT face=tet_faces(j);
                    bool sign=Order_Face_With_Sign(face);
                    int face_index=sign?1:2;
                    VECTOR<int,2> value=face_to_elements.Get(face);
                    PHYSBAM_ASSERT(value(face_index)!=0);
                    delayed_face_list.Append(tet_faces(j));
                    // search nearby elements
                    int other_face_index=sign?2:1;
                    int other_element_index=value(other_face_index);
                    if(other_element_index==0){// boundary face
                        face_list.Append(tet_faces(j));}
                    else if(!element_deleted(other_element_index)){
                        const VECTOR<int,TV::m+1>& e=mesh.elements(other_element_index);
                        if(StarShapeTest(e,face,p) && InCircle(e(1),e(2),e(3),p)){
                            element_deleted(other_element_index)=true;
                            deletion_list.Append(other_element_index);
                            stack.Append(other_element_index);}
                        else{
                            face_list.Append(tet_faces(j));}}}}
            // update reference data structure
            for(int j=1;j<=delayed_face_list.m;j++){
                TV_INT face=delayed_face_list(j);
                bool sign=Order_Face_With_Sign(face);
                int face_index=sign?1:2;
                VECTOR<int,2> value=face_to_elements.Get(face);
                value(face_index)=0;
                face_to_elements.Set(face,value);}
            // fill
            for(int j=1;j<=face_list.m;j++){
                TV_INT face=face_list(j);
                VECTOR<int,TV::m+1> new_element(face(1),face(2),p);
                mesh.elements.Append(new_element);
                element_deleted.Append(false);
                // add to bin
                RANGE<TV> bounding_box=RANGE<TV>::Bounding_Box(particles.X.Subset(new_element));
                RANGE<TV_INT> cell_range(Clamp_Cell(bounding_box.min_corner),Clamp_Cell(bounding_box.max_corner));
                for(UNIFORM_GRID_ITERATOR_CELL<TV> iter(grid,cell_range);iter.Valid();iter.Next()) cell_elements(iter.Cell_Index()).Append(mesh.elements.m);
                // update face-to-elements 
                VECTOR<TV_INT,TV::m+1> tet_faces=Get_Segments_In_Triangle(new_element);
                for(int k=1;k<=TV::m+1;k++){
                    TV_INT new_face=tet_faces(k);
                    bool sign=Order_Face_With_Sign(new_face);
                    int index=sign?1:2;
                    if(face_to_elements.Contains(new_face)){
                        VECTOR<int,2> value=face_to_elements.Get(new_face);
                        //PHYSBAM_ASSERT(value(index)==0);
                        value(index)=mesh.elements.m;
                        face_to_elements.Set(new_face,value);}
                    else{
                        VECTOR<int,2> value;value(index)=mesh.elements.m;
                        face_to_elements.Set(new_face,value);}}}}
        Sort(deletion_list);
        int last=mesh.elements.m;
        for(int k=deletion_list.m;k>=1;k--) mesh.elements(deletion_list(k))=mesh.elements(last--);
        mesh.elements.Resize(mesh.elements.m-deletion_list.m);
        PHYSBAM_DEBUG_WRITE_SUBSTEP("after adding dangerous particles for step 2",0,2);
    }

    void Set_Valid_Indices(const ARRAY<int>* valid_indices_input)
    {valid_indices=valid_indices_input;}

    void Set_Blocked_Cells(const ARRAY<bool,TV_INT>* blocked_cells_input)
    {blocked_cells=blocked_cells_input;}

    void Set_Grid_Radius(const T r=0.88)
    {grid_radius=r;}

    void Clear_Intersection_With_Blocked_Cells(const GRID<TV>& grid,const ARRAY<bool,TV_INT>& blocked_cells,TRIANGLE_MESH& mesh)
    {
        ARRAY<int> deletion_list;
        for(int i=1;i<=mesh.elements.m;i++){
            const VECTOR<int,3>& element=mesh.elements(i);
            bool intersected=false;
            for(int j=1;j<=3;j++){TV_INT face(element(j),element(j%3+1));
                SEGMENT_2D<T> segment(particles.X(face.x),particles.X(face.y));
                RANGE<TV_INT> local_range=RANGE<TV_INT>::Bounding_Box(grid.Cell(segment.x1,ghost_cells),grid.Cell(segment.x2,ghost_cells)).Thickened(1);
                for(UNIFORM_GRID_ITERATOR_CELL<TV> iter(grid,local_range);iter.Valid();iter.Next())if(grid.Inside_Domain(iter.Cell_Index(),ghost_cells) && blocked_cells(iter.Cell_Index())){
                    TV x=iter.Location();
                    RANGE<TV> box(x-grid.dX*(T)(0.5-1e-3),x+grid.dX*(T)(0.5-1e-3));
                    if(INTERSECTION::Intersects(box,segment)){intersected=true;break;}}
                if(intersected) break;}
            if(intersected) deletion_list.Append(i);}
        Sort(deletion_list);
        int last=mesh.elements.m;
        for(int k=deletion_list.m;k>=1;k--) mesh.elements(deletion_list(k))=mesh.elements(last--);
        mesh.elements.Resize(mesh.elements.m-deletion_list.m);
    }

protected:
    bool StarShapeTest(const VECTOR<int,TV::m+1>& e,const VECTOR<int,TV::m>& f,const int p)
    {
        for(int i=1;i<=TV::m+1;i++){
            TV_INT face(e(i),e(i%(TV::m+1)+1));
            if(Ordered_Face(face)==f) continue;
            if(ADAPTIVE_EXACT_PREDICATES::orient2d(VECTOR<double,2>(particles.X(face(1))).begin(),VECTOR<double,2>(particles.X(face(2))).begin(),VECTOR<double,2>(particles.X(p)).begin())<=0) return false;}
        return true;
    }

    void Initialize_Meshing_Parameters()
    {
        radius=grid.dX.Max()*grid_radius;
        radius_squared=sqr(radius);
        int_radius=(int)ceil(radius/grid.dX.Max()+.5);
    }

    bool InCircle(int ia, int ib, int ic, int ip)
    {
        double result=ADAPTIVE_EXACT_PREDICATES::incircle(VECTOR<double,2>(particles.X(ia)).begin(),VECTOR<double,2>(particles.X(ib)).begin(),VECTOR<double,2>(particles.X(ic)).begin(),VECTOR<double,2>(particles.X(ip)).begin());
        if(result==0){
            // symbolic perturbation
            ARRAY<int> local_ids(4),global_ids(4),sort_keys(4);
            local_ids(1)=1;local_ids(2)=2;local_ids(3)=3;local_ids(4)=4;
            global_ids(1)=ia;global_ids(2)=ib;global_ids(3)=ic;global_ids(4)=ip;
            if(mpi_particles){sort_keys(1)=particles.id(ia);sort_keys(2)=particles.id(ib);sort_keys(3)=particles.id(ic);sort_keys(4)=particles.id(ip);}
            else sort_keys=global_ids;
            Sort(local_ids,INDIRECT_COMPARE<ARRAY<int> >(sort_keys));
            for(int i=1;i<=4;i++){
                int skip=local_ids(i);
                ARRAY<VECTOR<double,2> > data(3);
                for(int local_id=1,count=0;local_id<=4;local_id++){
                    if(local_id!=skip){
                        data(++count)=VECTOR<double,2>(particles.X(global_ids(local_id)));}}
                double perturbation=ADAPTIVE_EXACT_PREDICATES::orient2d(data(1).begin(),data(2).begin(),data(3).begin());
                if(perturbation==0) continue;
                if(skip%2==0) perturbation=-perturbation;
                return perturbation>0;}
            PHYSBAM_FATAL_ERROR("InCircle test fails");
            return false;}
        return result>0;
    }

    void Clear_Auxilary_Data()
    {
        particle_finished.Resize(particles.X.m);ARRAYS_COMPUTATIONS::Fill(particle_finished,false);
        particle_counts.Resize(particles.X.m);ARRAYS_COMPUTATIONS::Fill(particle_counts,0);
        frontier.Remove_All();
        face_map.Clean_Memory();
        for(int i=1;i<=cell_to_faces.array.m;i++) cell_to_faces.array(i).Resize(0);
    }
    
    VECTOR<TV_INT,TV::m+1> Get_Segments_In_Triangle(VECTOR<int,TV::m+1> e)
    {
        VECTOR<TV_INT,TV::m+1> segs;
        segs(1)=TV_INT(e(1),e(2));
        segs(2)=TV_INT(e(2),e(3));
        segs(3)=TV_INT(e(3),e(1));
        return segs;
    }

    bool Initialize_Front_With_A_Single_Triangle(TRIANGLE_MESH& mesh)
    {
        VECTOR<int,3> initial_triangle=Find_Initial_Triangle();
        if(initial_triangle==VECTOR<int,3>()) return false;
        mesh.elements.Append(initial_triangle);
        for(int i=1;i<=3;i++){
            TV_INT initial_face(initial_triangle(i%3+1),initial_triangle(i));
            Insert_Face(initial_face);}
        return true;
    }

    void Advancing_Front(TRIANGLE_MESH& mesh)
    {
        while(!frontier.Empty()){
            TV_INT face=frontier.Dequeue();
            int type=face_map.Get(Ordered_Face(face));
            if(type==0) continue;
            int p=Find_Empty_Circumcircle_Point(face,type);
            if(p>0){
                VECTOR<int,3> new_triangle=VECTOR<int,3>(face.x,face.y,p);
                TV_INT new_face1(new_triangle.x,p),new_face2(p,new_triangle.y);
                if(Is_Face_Valid(new_face1) && Is_Face_Valid(new_face2)){
                    Insert_Face(new_face1);
                    Insert_Face(new_face2);
                    mesh.elements.Append(new_triangle);}
                else p=0;}
            Remove_Face(face);
            if(use_manifold_mesh){
                if(particle_counts(face.x)==0) particle_finished(face.x)=true;
                if(particle_counts(face.y)==0) particle_finished(face.y)=true;
                if(p>0 && particle_counts(p)==0) particle_finished(p)=true;}
            else{
                if(p>0){
                    if(particle_counts(face.x)==0) particle_finished(face.x)=true;
                    if(particle_counts(face.y)==0) particle_finished(face.y)=true;
                    if(particle_counts(p)==0) particle_finished(p)=true;}}}
    }

    TV_INT Clamp_Cell(const TV& p)
    {
        return grid.Domain_Indices().Clamp(grid.Cell(p,ghost_cells));
    }

    TV_INT Clamp_Cell(const TV_INT& p)
    {
        return grid.Domain_Indices().Clamp(p);
    }

    void Build_Grid_Indexing()
    {
        for(int i=1;i<=cell_particles.array.m;i++) cell_particles.array(i).Resize(0);
        if(!valid_indices) for(int i=1;i<=particles.X.m;i++){
            TV_INT cell=Clamp_Cell(particles.X(i));
            if(!grid.Inside_Domain(cell,ghost_cells)) continue;
            cell=Clamp_Cell(particles.X(i));
            cell_particles(cell).Append(i);}
        else for(int j=1;j<=valid_indices->m;j++){int i=(*valid_indices)(j);
            TV_INT cell=Clamp_Cell(particles.X(i));
            if(!grid.Inside_Domain(cell,ghost_cells)) continue;
            cell=Clamp_Cell(particles.X(i));
            cell_particles(cell).Append(i);}
    }

    void Remove_Face(TV_INT face)
    {
        face_map.Set(Ordered_Face(face),0);
        particle_counts(face.x)--;
        particle_counts(face.y)--;
    }

    void Insert_Face(TV_INT face,int type=1)
    {
        TV_INT ordered_face=Ordered_Face(face);
        int* value=face_map.Get_Pointer(ordered_face);
        if(!value){
            face_map.Insert(ordered_face,type);
            if(type==2) Add_Rasterized_Face(ordered_face);
            frontier.Safe_Enqueue(face);
            particle_counts(face.x)++;
            particle_counts(face.y)++;}
        else{
            if(*value!=0) Remove_Face(face);
            else PHYSBAM_FATAL_ERROR("This should not happen");}
    }

    bool Is_Face_Valid(TV_INT face,int type=1)
    {
        TV_INT ordered_face=Ordered_Face(face);
        int* value=face_map.Get_Pointer(ordered_face);
        if(!value) return true;
        else if(*value!=0) return true;
        else return false;
    }

    void Add_Rasterized_Face(TV_INT face)
    {
        RANGE<TV> box=RANGE<TV>::Bounding_Box(particles.X(face.x),particles.X(face.y));
        RANGE<TV_INT> cell_range(Clamp_Cell(box.min_corner),Clamp_Cell(box.max_corner));
        for(UNIFORM_GRID_ITERATOR_CELL<TV> iter(grid,cell_range);iter.Valid();iter.Next()) cell_to_faces(iter.Cell_Index()).Append(Ordered_Face(face));
    }

    void Remove_Rasterized_Face(TV_INT face)
    {
        RANGE<TV> box=RANGE<TV>::Bounding_Box(particles.X(face.x),particles.X(face.y));
        RANGE<TV_INT> cell_range(Clamp_Cell(box.min_corner),Clamp_Cell(box.max_corner));
        for(UNIFORM_GRID_ITERATOR_CELL<TV> iter(grid,cell_range);iter.Valid();iter.Next()){
            ARRAY<TV_INT>& list=cell_to_faces(iter.Cell_Index());
            for(int i=1;i<=list.m;i++) if(list(i)==Ordered_Face(face)){list.Remove_Index_Lazy(i);i--;}}
    }

    bool Face_Collides_With_Front(TV_INT face)
    {
        SEGMENT_2D<T> segment1(particles.X(face.x),particles.X(face.y));
        RANGE<TV> box=RANGE<TV>::Bounding_Box(particles.X(face.x),particles.X(face.y));
        RANGE<TV_INT> cell_range(Clamp_Cell(box.min_corner),Clamp_Cell(box.max_corner));
        for(UNIFORM_GRID_ITERATOR_CELL<TV> iter(grid,cell_range);iter.Valid();iter.Next()){
            ARRAY<TV_INT>& list=cell_to_faces(iter.Cell_Index());
            for(int i=1;i<=list.m;i++){
                const TV_INT& face2=list(i);
                if(face.x==face2.x || face.x==face2.y || face.y==face2.x || face.y==face2.y) continue;
                SEGMENT_2D<T> segment2(particles.X(face2.x),particles.X(face2.y));
                if(INTERSECTION::Intersects(segment1,segment2)) return true;}}
        return false;
    }

    bool Element_Collides_With_Front(VECTOR<int,3> element)
    {
        TRIANGLE_2D<T> triangle(particles.X(element.x),particles.X(element.y),particles.X(element.z));
        RANGE<TV> box=RANGE<TV>::Bounding_Box(particles.X(element.x),particles.X(element.y),particles.X(element.z));
        RANGE<TV_INT> cell_range(Clamp_Cell(box.min_corner),Clamp_Cell(box.max_corner));
        for(UNIFORM_GRID_ITERATOR_CELL<TV> iter(grid,cell_range);iter.Valid();iter.Next()){
            ARRAY<int>& list=cell_particles(iter.Cell_Index());
            for(int i=1;i<=list.m;i++){
                int particle_index=list(i);
                if(element.x==particle_index || element.y==particle_index || element.z==particle_index) continue;
                if(triangle.Inside(particles.X(particle_index))) return true;}}
        return Face_Collides_With_Front(TV_INT(element.x,element.z)) || Face_Collides_With_Front(TV_INT(element.y,element.z));
    }

    TV_INT Ordered_Face(TV_INT face){return TV_INT(min(face.x,face.y),max(face.x,face.y));};

    bool Order_Face_With_Sign(TV_INT& f){
        bool sign=true;
        if(f(1)>f(2)){int tmp=f(1);f(1)=f(2);f(2)=tmp;sign=!sign;}
        return sign;
    }

    bool Triangle_In_Mesh(const int p0,const int p1,const int p2)
    {
        TV_INT f0=Ordered_Face(TV_INT(p0,p1));
        TV_INT f1=Ordered_Face(TV_INT(p1,p2));
        TV_INT f2=Ordered_Face(TV_INT(p2,p0));
        return (face_map.Get_Pointer(f0)||face_map.Get_Pointer(f1)||face_map.Get_Pointer(f2));
    }

    VECTOR<int,3> Find_Initial_Triangle()
    {
        for(UNIFORM_GRID_ITERATOR_CELL<TV> iter(grid,ghost_cells);iter.Valid();iter.Next()){
            TV_INT cell=iter.Cell_Index();
            for(int i=1;i<=cell_particles(cell).m;i++){int p0=cell_particles(cell)(i);
                if(particle_finished(p0)) continue;
                int p1=Find_Nearest_Point(p0);
                if(p1>0){
                    int p2=Find_Empty_Circumcircle_Point(TV_INT(p0,p1));
                    if(p2>0 && Is_Face_Valid(TV_INT(p0,p1)) && Is_Face_Valid(TV_INT(p1,p2)) && Is_Face_Valid(TV_INT(p0,p2)) &&
                        (use_manifold_mesh||!Triangle_In_Mesh(p0,p1,p2)))
                        return VECTOR<int,3>(p0,p1,p2);
                    p2=Find_Empty_Circumcircle_Point(TV_INT(p1,p0));
                    if(p2>0 && Is_Face_Valid(TV_INT(p0,p1)) && Is_Face_Valid(TV_INT(p1,p2)) && Is_Face_Valid(TV_INT(p0,p2)) &&
                        (use_manifold_mesh||!Triangle_In_Mesh(p0,p1,p2))) 
                        return VECTOR<int,3>(p1,p0,p2);}}}
        return VECTOR<int,3>();
    }

    int Find_Empty_Circumcircle_Point(const TV_INT& face,int type=1)
    {
        TV x1=particles.X(face.x),x2=particles.X(face.y);
        TV center=(x1+x2)*(T)0.5;
        TV face_normal=(x2-x1).Orthogonal_Vector().Normalized();
        center+=face_normal*sqrt(sqr(radius)-(x1-center).Magnitude_Squared());
        TV_INT cell=Clamp_Cell(center);
        RANGE<TV_INT> range(cell-int_radius,cell+int_radius);
        int minimum_p=0;
        for(UNIFORM_GRID_ITERATOR_CELL<TV> iter(grid,range);iter.Valid();iter.Next())if(grid.Inside_Domain(iter.Cell_Index(),ghost_cells)){
            for(int i=1;i<=cell_particles(iter.Cell_Index()).m;i++){int p=cell_particles(iter.Cell_Index())(i);
                if(p==face.x || p==face.y) continue;
                if((particles.X(p)-center).Magnitude_Squared()>radius_squared*1.1) continue;
                if(ADAPTIVE_EXACT_PREDICATES::orient2d(VECTOR<double,2>(x1).begin(),VECTOR<double,2>(x2).begin(),VECTOR<double,2>(particles.X(p)).begin())<=0) continue;
                if(minimum_p!=0) if(!InCircle(face.x,face.y,minimum_p,p)) continue;
                minimum_p=p;}}
        if(minimum_p!=0){
            if(particle_finished(minimum_p)) return 0;
            if(use_alpha_shape && type==1 && (particles.X(minimum_p)-center).Magnitude_Squared()>radius_squared) return 0;}
        return minimum_p;
    }

    int Find_Nearest_Point(const int index)
    {
        TV x=particles.X(index);
        TV_INT cell=Clamp_Cell(x);
        RANGE<TV_INT> range(cell-int_radius,cell+int_radius);
        int minimum_p=0;T minimum_dd=FLT_MAX;
        for(UNIFORM_GRID_ITERATOR_CELL<TV> iter(grid,range);iter.Valid();iter.Next())if(grid.Inside_Domain(iter.Cell_Index(),ghost_cells)){
            const ARRAY<int>& list=cell_particles(iter.Cell_Index());
            for(int i=1;i<=list.m;i++){int p=list(i);
                if(p==index || particle_finished(p)) continue;
                TV y=particles.X(p);
                T dd=(x-y).Magnitude();
                if(dd<minimum_dd && dd<2*radius){
                    if(Face_Collides_With_Front(TV_INT(index,p))) continue;
                    minimum_dd=dd;
                    minimum_p=p;}}}
        return minimum_p;
    }
};

// 3D version
template<class T> 
class DELAUNAY_MESHING_WITH_ALPHA_SHAPE<VECTOR<T,3> >
{
    typedef VECTOR<T,3> TV;
    typedef VECTOR<int,3> TV_INT;
    typedef TETRAHEDRON<T> T_ELEMENT;

public:
    const GRID<TV>& grid;
    GEOMETRY_PARTICLES<TV>& particles;
    ARRAY<ARRAY<int>,TV_INT> cell_particles;
    ARRAY<bool> particle_finished;
    ARRAY<int> particle_counts;
    QUEUE<TV_INT> frontier;
    HASHTABLE<TV_INT,int> face_map;
    T grid_radius,grid_radius_squared,grid_radius_ratio;
    int grid_r_steps_in_bkgrid,ghost_cells;
    const ARRAY<int>* valid_indices;
    ARRAY<bool> dangerous_indices;
    const ARRAY<bool,TV_INT>* blocked_cells;
    ARRAY<ARRAY<TV_INT>,TV_INT> cell_to_faces;
    ARRAY<bool> particle_on_boundary;
    bool use_alpha_shape;bool use_manifold_mesh;
    MPI_PARTICLES<GRID<TV> >* mpi_particles,* mpi_particles2;
    bool clean_memory_after_meshing;

    DELAUNAY_MESHING_WITH_ALPHA_SHAPE(const GRID<TV>& grid_input,GEOMETRY_PARTICLES<TV>& particles_input)
        :grid(grid_input),particles(particles_input),frontier(1),grid_radius_ratio((T).88),valid_indices(0),blocked_cells(0),use_alpha_shape(true),use_manifold_mesh(true),mpi_particles(0),mpi_particles2(0),clean_memory_after_meshing(false)
    {
        ADAPTIVE_EXACT_PREDICATES::exactinit();
    }

    ~DELAUNAY_MESHING_WITH_ALPHA_SHAPE(){}

    void Initialize_Grid(int ghost_cells_input=3)
    {
        ghost_cells=ghost_cells_input;
        cell_particles.Resize(grid.Domain_Indices(ghost_cells));
        cell_to_faces.Resize(grid.Domain_Indices(ghost_cells));
    }

    void Construct_Mesh(TETRAHEDRON_MESH& mesh)
    {
        LOG::SCOPE scope("Construct Mesh");
        if(clean_memory_after_meshing) Initialize_Grid(ghost_cells);
        Initialize_Meshing_Parameters();
        Build_Grid_Indexing();
        mesh.elements.Clean_Memory();
        Clear_Auxilary_Data();
        while(Initialize_Front_With_A_Single_Tetrahedron(mesh)) Advancing_Front(mesh);
        if(blocked_cells) Clear_Intersection_With_Blocked_Cells(grid,*blocked_cells,mesh);
        if(clean_memory_after_meshing) Clean_Memory();
    }

    void Construct_Constrained_Mesh(TETRAHEDRON_MESH& mesh,TRIANGLE_MESH& boundary_mesh,bool skip_dangenous=false,INSIDE_DOMAIN_PREDICATOR<TV>* domain=0,INSIDE_DOMAIN_PREDICATOR<TV>* mpi_local_domain=0)
    {
        LOG::SCOPE scope("Construct Constrained Mesh");
        // index boundary elements
        LOG::Time("Index boundary circumspheres");
        boundary_mesh.Initialize_Segment_Mesh();
        ARRAY<ARRAY<int>,TV_INT> cell_elements(grid.Domain_Indices(ghost_cells));
        ARRAY<TV> element_circumcener(boundary_mesh.elements.m+boundary_mesh.segment_mesh->elements.m);
        ARRAY<T> element_circumradius(boundary_mesh.elements.m+boundary_mesh.segment_mesh->elements.m);
        ARRAY<bool> is_boundary_nodes(particles.X.m);
        ARRAYS_COMPUTATIONS::Fill(is_boundary_nodes,false);
        for(int i=1;i<=boundary_mesh.elements.m;i++){
            const TV_INT& e=boundary_mesh.elements(i);
            for(int j=1;j<=TV::m;j++) is_boundary_nodes(e(j))=true;
            element_circumcener(i)=TRIANGLE_3D<T>::Circumcenter(particles.X(e.x),particles.X(e.y),particles.X(e.z));
            element_circumradius(i)=(element_circumcener(i)-particles.X(e.x)).Magnitude();
            if(element_circumradius(i)<grid.dX.Max()*(T).5){
                TV_INT cell=Clamp_Cell(element_circumcener(i));
                int range=(int)ceil(element_circumradius(i)/grid.dX.Max());
                RANGE<TV_INT> cell_range(Clamp_Cell(cell-range),Clamp_Cell(cell+range));
                for(UNIFORM_GRID_ITERATOR_CELL<TV> iter(grid,cell_range);iter.Valid();iter.Next()) cell_elements(iter.Cell_Index()).Append(i);}}
        for(int i=1;i<=boundary_mesh.segment_mesh->elements.m;i++){
            const VECTOR<int,TV::m-1>& e=boundary_mesh.segment_mesh->elements(i);
            int index=boundary_mesh.elements.m+i;
            element_circumcener(index)=(particles.X(e.x)+particles.X(e.y))*0.5;
            element_circumradius(index)=(element_circumcener(index)-particles.X(e.x)).Magnitude();
            TV_INT cell=Clamp_Cell(element_circumcener(index));
            int range=(int)ceil(element_circumradius(index)/grid.dX.Max());
            RANGE<TV_INT> cell_range(Clamp_Cell(cell-range),Clamp_Cell(cell+range));
            for(UNIFORM_GRID_ITERATOR_CELL<TV> iter(grid,cell_range);iter.Valid();iter.Next()) cell_elements(iter.Cell_Index()).Append(index);}
        boundary_mesh.Delete_Auxiliary_Structures();
        // filter out dangerous particles
        LOG::Time("Filter out dangerous particles");
        ARRAY<int> full_indices,safe_indices,dangerous_indices;
        if(valid_indices) full_indices=*valid_indices;
        else{full_indices.Resize(particles.X.m);for(int i=1;i<=particles.X.m;i++) full_indices(i)=i;}
        ARRAY<bool> dangerous(particles.X.m);ARRAYS_COMPUTATIONS::Fill(dangerous,false);
        for(int i=1;i<=full_indices.m;i++){
            int p=full_indices(i);
            if(is_boundary_nodes(p)) continue;
            TV_INT cell=Clamp_Cell(particles.X(p));
            const ARRAY<int>& list=cell_elements(cell);
            for(int j=1;j<=list.m;j++)if((particles.X(p)-element_circumcener(list(j))).Magnitude()<element_circumradius(list(j))+1e-8*grid.dX.Min()){
                dangerous(p)=true;break;}}
        if(mpi_particles || mpi_particles2){
            ARRAY_VIEW<bool> dangerous_view(dangerous);
            if(mpi_particles) mpi_particles->Exchange_Boundary_Particle_Values(dangerous_view);
            if(mpi_particles2) mpi_particles2->Exchange_Boundary_Particle_Values(dangerous_view);}
        for(int i=1;i<=full_indices.m;i++){
            int p=full_indices(i);
            if(dangerous(p)) dangerous_indices.Append(p);
            else safe_indices.Append(p);}
        dangerous.Clean_Memory();
        is_boundary_nodes.Clean_Memory();
        element_circumcener.Clean_Memory();
        element_circumradius.Clean_Memory();
        // meshing without dangerous particles
        Set_Valid_Indices(&safe_indices);
        Construct_Mesh(mesh);
        PHYSBAM_DEBUG_WRITE_SUBSTEP("after constructing mesh for step 2",0,2);
        if(domain){
            ARRAY<int> deletion_list;
            for(int i=1;i<=mesh.elements.m;i++){
                const VECTOR<int,TV::m+1>& element=mesh.elements(i);
                TV barycenter=T_ELEMENT::Center(particles.X.Subset(element));
                if(!(*domain)(barycenter) || (mpi_local_domain && !(*mpi_local_domain)(barycenter))){deletion_list.Append(i);continue;}
                TV circumcenter=T_ELEMENT::Circumcenter(particles.X.Subset(element));
                if(!(*domain)(circumcenter))deletion_list.Append(i);}
            Sort(deletion_list);
            int last=mesh.elements.m;
            for(int k=deletion_list.m;k>=1;k--) mesh.elements(deletion_list(k))=mesh.elements(last--);
            mesh.elements.Resize(last);}
        PHYSBAM_DEBUG_WRITE_SUBSTEP("after deleting outside elements for step 2",0,2);
        if(skip_dangenous) return;
        // add in dangerous particles
        LOG::Time("Index all elements");
        Set_Valid_Indices(&full_indices);
        // bin all elements for point location and construct face-to-element reference
        for(int i=1;i<=cell_elements.array.m;i++) cell_elements.array(i).Resize(0);
        HASHTABLE<TV_INT,VECTOR<int,2> > face_to_elements;
        for(int i=1;i<=mesh.elements.m;i++){
            const VECTOR<int,TV::m+1>& e=mesh.elements(i);
            RANGE<TV> bounding_box=RANGE<TV>::Bounding_Box(particles.X.Subset(e));
            RANGE<TV_INT> cell_range(Clamp_Cell(bounding_box.min_corner),Clamp_Cell(bounding_box.max_corner));
            for(UNIFORM_GRID_ITERATOR_CELL<TV> iter(grid,cell_range);iter.Valid();iter.Next()) cell_elements(iter.Cell_Index()).Append(i);
            VECTOR<TV_INT,TV::m+1> tet_faces=Get_Triangles_In_Tetrahedron(e);
            for(int j=1;j<=TV::m+1;j++){
                TV_INT face=tet_faces(j);
                bool sign=Order_Face_With_Sign(face);
                int index=sign?1:2;
                if(face_to_elements.Contains(face)){
                    VECTOR<int,2> value=face_to_elements.Get(face);
                    PHYSBAM_ASSERT(value(index)==0);
                    value(index)=i;
                    face_to_elements.Set(face,value);}
                else{
                    VECTOR<int,2> value;value(index)=i;
                    face_to_elements.Set(face,value);}}}
        // incrementally insert particles
        LOG::Time("Add in dangerous particles");
        ARRAY<int> deletion_list;
        ARRAY<bool> element_deleted(mesh.elements.m);
        ARRAYS_COMPUTATIONS::Fill(element_deleted,false);
        for(int i=1;i<=dangerous_indices.m;i++){
            int p=dangerous_indices(i);
            // locate
            TV_INT cell=Clamp_Cell(particles.X(p));
            const ARRAY<int>& list=cell_elements(cell);
            int initial_element_index=0;
            for(int j=1;j<=list.m;j++)if(!element_deleted(list(j))){
                TETRAHEDRON<T> element(particles.X.Subset(mesh.elements(list(j))));
                if(element.Inside(particles.X(p))){
                    initial_element_index=list(j);
                    break;}}
            if(initial_element_index==0) continue;// the point is outside
            // expand (depth-first search)
            ARRAY<int> stack;
            ARRAY<TV_INT> delayed_face_list;
            ARRAY<TV_INT> face_list;
            stack.Append(initial_element_index);
            element_deleted(initial_element_index)=true;
            deletion_list.Append(initial_element_index);
            while(stack.m>0){
                int element_index=stack.Pop();
                const VECTOR<int,TV::m+1>& current_element=mesh.elements(element_index);
                VECTOR<TV_INT,TV::m+1> tet_faces=Get_Triangles_In_Tetrahedron(current_element);
                for(int j=1;j<=TV::m+1;j++){
                    TV_INT face=tet_faces(j);
                    bool sign=Order_Face_With_Sign(face);
                    int face_index=sign?1:2;
                    VECTOR<int,2> value=face_to_elements.Get(face);
                    PHYSBAM_ASSERT(value(face_index)!=0);
                    delayed_face_list.Append(tet_faces(j));
                    // search nearby elements
                    int other_face_index=sign?2:1;
                    int other_element_index=value(other_face_index);
                    if(other_element_index==0){// boundary face
                        face_list.Append(tet_faces(j));}
                    else if(!element_deleted(other_element_index)){
                        const VECTOR<int,TV::m+1>& e=mesh.elements(other_element_index);
                        if(StarShapeTest(e,face,p) && InSphere(e(1),e(2),e(3),e(4),p)){
                            element_deleted(other_element_index)=true;
                            deletion_list.Append(other_element_index);
                            stack.Append(other_element_index);}
                        else{
                            face_list.Append(tet_faces(j));}}}}
            // update reference data structure
            for(int j=1;j<=delayed_face_list.m;j++){
                TV_INT face=delayed_face_list(j);
                bool sign=Order_Face_With_Sign(face);
                int face_index=sign?1:2;
                VECTOR<int,2> value=face_to_elements.Get(face);
                value(face_index)=0;
                face_to_elements.Set(face,value);}
            // fill
            for(int j=1;j<=face_list.m;j++){
                TV_INT face=face_list(j);
                VECTOR<int,TV::m+1> new_element(face(1),face(2),face(3),p);
                mesh.elements.Append(new_element);
                element_deleted.Append(false);
                // add to bin
                RANGE<TV> bounding_box=RANGE<TV>::Bounding_Box(particles.X.Subset(new_element));
                RANGE<TV_INT> cell_range(Clamp_Cell(bounding_box.min_corner),Clamp_Cell(bounding_box.max_corner));
                for(UNIFORM_GRID_ITERATOR_CELL<TV> iter(grid,cell_range);iter.Valid();iter.Next()) cell_elements(iter.Cell_Index()).Append(mesh.elements.m);
                // update face-to-elements 
                VECTOR<TV_INT,TV::m+1> tet_faces=Get_Triangles_In_Tetrahedron(new_element);
                for(int k=1;k<=TV::m+1;k++){
                    TV_INT new_face=tet_faces(k);
                    bool sign=Order_Face_With_Sign(new_face);
                    int index=sign?1:2;
                    if(face_to_elements.Contains(new_face)){
                        VECTOR<int,2> value=face_to_elements.Get(new_face);
                        //PHYSBAM_ASSERT(value(index)==0);
                        value(index)=mesh.elements.m;
                        face_to_elements.Set(new_face,value);}
                    else{
                        VECTOR<int,2> value;value(index)=mesh.elements.m;
                        face_to_elements.Set(new_face,value);}}}}
        Sort(deletion_list);
        int last=mesh.elements.m;
        for(int k=deletion_list.m;k>=1;k--) mesh.elements(deletion_list(k))=mesh.elements(last--);
        mesh.elements.Resize(mesh.elements.m-deletion_list.m);
        PHYSBAM_DEBUG_WRITE_SUBSTEP("after adding dangerous particles for step 2",0,2);
    }

    void Set_Grid_Radius(const T r=0.88)
    {grid_radius_ratio=r;}
    
    void Set_Valid_Indices(const ARRAY<int>* valid_indices_input)
    {valid_indices=valid_indices_input;}
    
    void Set_Blocked_Cells(const ARRAY<bool,TV_INT>* blocked_cells_input)
    {blocked_cells=blocked_cells_input;}

    void Clear_Intersection_With_Blocked_Cells(const GRID<TV>& grid,const ARRAY<bool,TV_INT>& blocked_cells,TETRAHEDRON_MESH& mesh)
    {
        ARRAY<int> deletion_list;
        for(int i=1;i<=mesh.elements.m;i++){
            const VECTOR<int,4>& element=mesh.elements(i);
            VECTOR<TV_INT,4> faces=Get_Triangles_In_Tetrahedron(element);
            bool intersected=false;
            for(int j=1;j<=4;j++){TV_INT face=faces(j);
                TRIANGLE_3D<T> triangle(particles.X(face.x),particles.X(face.y),particles.X(face.z));
                RANGE<TV> bounding_box=triangle.Bounding_Box();
                RANGE<TV_INT> local_range(grid.Cell(bounding_box.min_corner,ghost_cells)-1,grid.Cell(bounding_box.max_corner,ghost_cells)+1);
                for(UNIFORM_GRID_ITERATOR_CELL<TV> iter(grid,local_range);iter.Valid();iter.Next())if(grid.Inside_Domain(iter.Cell_Index(),ghost_cells) && blocked_cells(iter.Cell_Index())){
                    TV x=iter.Location();
                    RANGE<TV> box(x-grid.dX*(T)(0.5-1e-3),x+grid.dX*(T)(0.5-1e-3));
                    if(INTERSECTION::Intersects(box,triangle)){intersected=true;break;};}
                if(intersected) break;}
            if(intersected) deletion_list.Append(i);}
        Sort(deletion_list);
        int last=mesh.elements.m;
        for(int k=deletion_list.m;k>=1;k--) mesh.elements(deletion_list(k))=mesh.elements(last--);
        mesh.elements.Resize(mesh.elements.m-deletion_list.m);
    }

protected:
    bool StarShapeTest(const VECTOR<int,TV::m+1>& e,const VECTOR<int,TV::m>& f,const int p)
    {
        VECTOR<VECTOR<int,TV::m>,TV::m+1> faces=Get_Triangles_In_Tetrahedron(e);
        for(int i=1;i<=TV::m+1;i++){
            TV_INT face=faces(i);
            if(Ordered_Face(face)==f) continue;
            if(ADAPTIVE_EXACT_PREDICATES::orient3d(VECTOR<double,3>(particles.X(face(1))).begin(),VECTOR<double,3>(particles.X(face(2))).begin(),VECTOR<double,3>(particles.X(face(3))).begin(),VECTOR<double,3>(particles.X(p)).begin())>=0) return false;}
        return true;
    }

    VECTOR<VECTOR<int,3>,4> Get_Triangles_In_Tetrahedron(VECTOR<int,4> tetra_index,bool test_signed=false)
    {
        VECTOR<TV_INT,4> tetra_tris;
        if(!test_signed || TETRAHEDRON<T>::Signed_Size(particles.X.Subset(tetra_index))>0){
            tetra_tris(1)=TV_INT(tetra_index(1),tetra_index(2),tetra_index(3));
            tetra_tris(2)=TV_INT(tetra_index(1),tetra_index(4),tetra_index(2));
            tetra_tris(3)=TV_INT(tetra_index(1),tetra_index(3),tetra_index(4));
            tetra_tris(4)=TV_INT(tetra_index(2),tetra_index(4),tetra_index(3));}
        else{
            tetra_tris(1)=TV_INT(tetra_index(1),tetra_index(3),tetra_index(2));
            tetra_tris(2)=TV_INT(tetra_index(1),tetra_index(2),tetra_index(4));
            tetra_tris(3)=TV_INT(tetra_index(1),tetra_index(4),tetra_index(3));
            tetra_tris(4)=TV_INT(tetra_index(2),tetra_index(3),tetra_index(4));}
        return tetra_tris;
    }

    void Initialize_Meshing_Parameters()
    {
        grid_radius=grid.dX.Max()*grid_radius_ratio;
        grid_radius_squared=grid_radius*grid_radius;
        grid_r_steps_in_bkgrid=(int)ceil(grid_radius/grid.dX.Max()+(T).5);
    }

    bool InSphere(int ia, int ib, int ic, int id, int ip)
    {
        double result=ADAPTIVE_EXACT_PREDICATES::insphere(VECTOR<double,3>(particles.X(ia)).begin(),VECTOR<double,3>(particles.X(ib)).begin(),VECTOR<double,3>(particles.X(ic)).begin(),VECTOR<double,3>(particles.X(id)).begin(),VECTOR<double,3>(particles.X(ip)).begin());
        if(result==0){
            // symbolic perturbation
            ARRAY<int> local_ids(5),global_ids(5),sort_keys(5);
            local_ids(1)=1;local_ids(2)=2;local_ids(3)=3;local_ids(4)=4;local_ids(5)=5;
            global_ids(1)=ia;global_ids(2)=ib;global_ids(3)=ic;global_ids(4)=id;global_ids(5)=ip;
            if(mpi_particles){sort_keys(1)=particles.id(ia);sort_keys(2)=particles.id(ib);sort_keys(3)=particles.id(ic);sort_keys(4)=particles.id(id);sort_keys(5)=particles.id(ip);}
            else sort_keys=global_ids;
            Sort(local_ids,INDIRECT_COMPARE<ARRAY<int> >(sort_keys));
            for(int i=1;i<=5;i++){
                int skip=local_ids(i);
                ARRAY<VECTOR<double,3> > data(4);
                for(int local_id=1,count=0;local_id<=5;local_id++){
                    if(local_id!=skip){
                        data(++count)=VECTOR<double,3>(particles.X(global_ids(local_id)));}}
                double perturbation=ADAPTIVE_EXACT_PREDICATES::orient3d(data(1).begin(),data(2).begin(),data(3).begin(),data(4).begin());
                if(perturbation==0) continue;
                if(skip%2==0) perturbation=-perturbation;
                return perturbation<0;}
            PHYSBAM_FATAL_ERROR("InSphere fails");
            return false;}
        return result<0;
    }

    void Clear_Auxilary_Data()
    {
        particle_finished.Resize(particles.X.m);ARRAYS_COMPUTATIONS::Fill(particle_finished,false);
        particle_counts.Resize(particles.X.m);ARRAYS_COMPUTATIONS::Fill(particle_counts,0);
        face_map.Clean_Memory();
        frontier.Remove_All();frontier.array.Resize(1000);
        particle_on_boundary.Resize(particles.X.m);ARRAYS_COMPUTATIONS::Fill(particle_on_boundary,false);
        for(int i=1;i<=cell_to_faces.array.m;i++) cell_to_faces.array(i).Resize(0);
    }

    void Clean_Memory()
    {
        particle_finished.Clean_Memory();
        particle_counts.Clean_Memory();
        face_map.Clean_Memory();
        frontier.Remove_All();frontier.array.Clean_Memory();
        particle_on_boundary.Clean_Memory();
        cell_to_faces.Clean_Memory();
    }

    bool Initialize_Front_With_A_Single_Tetrahedron(TETRAHEDRON_MESH& mesh)
    {
        VECTOR<int,4> initial_tetrahedron=Find_Initial_Tetrahedron();
        if(initial_tetrahedron==VECTOR<int,4>()) return false;
        VECTOR<TV_INT,4> tetra_tris=Get_Triangles_In_Tetrahedron(initial_tetrahedron);
        for(int i=1;i<=4;i++){
            TV_INT initial_face(tetra_tris(i)(1),tetra_tris(i)(3),tetra_tris(i)(2));    ////inversed sequence
            Insert_Face(initial_face);}
        mesh.elements.Append(initial_tetrahedron);
        return true;
    }

    void Advancing_Front(TETRAHEDRON_MESH& mesh)
    {
        while(!frontier.Empty()){
            TV_INT face=frontier.Dequeue();
            int type=face_map.Get(Ordered_Face(face));
            if(type==0) continue;
            int p=Find_Empty_Circumsphere_Point_With_Grid_Radius(face,type);
            if(p>0){
                VECTOR<int,4> new_tet=VECTOR<int,4>(face(1),face(2),face(3),p);
                TV_INT new_face1(new_tet(1),new_tet(2),new_tet(4)),////inversed sequence
                    new_face2(new_tet(1),new_tet(4),new_tet(3)),////inversed sequence
                    new_face3(new_tet(2),new_tet(3),new_tet(4));////inversed sequence
                if(Is_Face_Valid(new_face1) && Is_Face_Valid(new_face2) && Is_Face_Valid(new_face3)){
                    Insert_Face(new_face1);
                    Insert_Face(new_face2);
                    Insert_Face(new_face3);
                    mesh.elements.Append(new_tet);}
                else p=0;}
            Remove_Face(face);
            for(int i=1;i<=3;i++) if(particle_counts(face(i))==0) particle_finished(face(i))=true;
            if(p>0 && particle_counts(p)==0) particle_finished(p)=true;}
    }

    TV_INT Clamp_Cell(const TV& p)
    {
        return grid.Domain_Indices().Clamp(grid.Cell(p,ghost_cells));
    }

    TV_INT Clamp_Cell(const TV_INT& p)
    {
        return grid.Domain_Indices().Clamp(p);
    }

    void Build_Grid_Indexing()
    {
        for(int i=1;i<=cell_particles.array.m;i++) cell_particles.array(i).Resize(0);
        if(!valid_indices) for(int i=1;i<=particles.X.m;i++){
            TV_INT cell=Clamp_Cell(particles.X(i));
            if(!grid.Inside_Domain(cell,ghost_cells)) continue;
            cell=Clamp_Cell(particles.X(i));
            cell_particles(cell).Append(i);}
        else for(int j=1;j<=valid_indices->m;j++){int i=(*valid_indices)(j);
            TV_INT cell=Clamp_Cell(particles.X(i));
            if(!grid.Inside_Domain(cell,ghost_cells)) continue;
            cell=Clamp_Cell(particles.X(i));
            cell_particles(cell).Append(i);}
    }

    void Remove_Face(TV_INT face)
    {
        face_map.Set(Ordered_Face(face),0);
        for(int i=1;i<=3;i++) particle_counts(face(i))--;
    }

    void Insert_Face(TV_INT face,int type=1)
    {
        TV_INT ordered_face=Ordered_Face(face);
        int* value=face_map.Get_Pointer(ordered_face);
        if(!value){
            face_map.Insert(ordered_face,type);
            frontier.Safe_Enqueue(face);
            particle_counts(face.x)++;
            particle_counts(face.y)++;
            particle_counts(face.z)++;}
        else{
            if(*value!=0) Remove_Face(face);
            else PHYSBAM_FATAL_ERROR("This should not happen.");}
    }

    bool Is_Face_Valid(TV_INT face,int type=1)
    {
        TV_INT ordered_face=Ordered_Face(face);
        int* value=face_map.Get_Pointer(ordered_face);
        if(!value) return true;
        else if(*value!=0) return true;
        else return false;
    }

    TV_INT Ordered_Face(TV_INT face){
        TV_INT f=face;
        if(f(1)>f(2)){int tmp=f(1);f(1)=f(2);f(2)=tmp;}
        if(f(2)>f(3)){int tmp=f(2);f(2)=f(3);f(3)=tmp;}
        if(f(1)>f(2)){int tmp=f(1);f(1)=f(2);f(2)=tmp;}
        return f;
    }

    bool Order_Face_With_Sign(TV_INT& f){
        bool sign=true;
        if(f(1)>f(2)){int tmp=f(1);f(1)=f(2);f(2)=tmp;sign=!sign;}
        if(f(2)>f(3)){int tmp=f(2);f(2)=f(3);f(3)=tmp;sign=!sign;}
        if(f(1)>f(2)){int tmp=f(1);f(1)=f(2);f(2)=tmp;sign=!sign;}
        return sign;
    }

    TV_INT Inverted_Face(TV_INT face){
        TV_INT f=face;
        int tmp=f(1);f(1)=f(2);f(2)=tmp;
        return f;
    }

    void Add_Rasterized_Face(TV_INT face)
    {
        RANGE<TV> box=RANGE<TV>::Bounding_Box(particles.X(face.x),particles.X(face.y),particles.X(face.z));
        RANGE<TV_INT> cell_range(Clamp_Cell(box.min_corner),Clamp_Cell(box.max_corner));
        for(UNIFORM_GRID_ITERATOR_CELL<TV> iter(grid,cell_range);iter.Valid();iter.Next()) cell_to_faces(iter.Cell_Index()).Append(Ordered_Face(face));
    }

    void Remove_Rasterized_Face(TV_INT face)
    {
        RANGE<TV> box=RANGE<TV>::Bounding_Box(particles.X(face.x),particles.X(face.y));
        RANGE<TV_INT> cell_range(Clamp_Cell(box.min_corner),Clamp_Cell(box.max_corner));
        for(UNIFORM_GRID_ITERATOR_CELL<TV> iter(grid,cell_range);iter.Valid();iter.Next()){
            ARRAY<TV_INT>& list=cell_to_faces(iter.Cell_Index());
            for(int i=1;i<=list.m;i++) if(list(i)==Ordered_Face(face)){list.Remove_Index_Lazy(i);i--;}}
    }

    bool Segment_Collides_With_Front(VECTOR<int,2> seg)
    {
        HASHTABLE<TV_INT> visited_faces;
        SEGMENT_3D<T> segment(particles.X(seg.x),particles.X(seg.y));
        RANGE<TV> box=RANGE<TV>::Bounding_Box(particles.X(seg.x),particles.X(seg.y));
        RANGE<TV_INT> cell_range(Clamp_Cell(box.min_corner),Clamp_Cell(box.max_corner));
        for(UNIFORM_GRID_ITERATOR_CELL<TV> iter(grid,cell_range);iter.Valid();iter.Next()){
            ARRAY<TV_INT>& list=cell_to_faces(iter.Cell_Index());
            for(int i=1;i<=list.m;i++){
                const TV_INT& face2=list(i);
                if(visited_faces.Contains(face2)) continue;
                visited_faces.Set(face2);
                if(seg.x==face2.x || seg.x==face2.y  || seg.x==face2.z || seg.y==face2.x || seg.y==face2.y || seg.y==face2.z) continue;
                TRIANGLE_3D<T> triangle2(particles.X(face2.x),particles.X(face2.y),particles.X(face2.z));
                if(INTERSECTION::Intersects(segment,triangle2)) return true;}}
        return false;
    }

    bool Face_Collides_With_Front(TV_INT face)
    {
        HASHTABLE<TV_INT> visited_faces;
        TRIANGLE_3D<T> triangle1(particles.X(face.x),particles.X(face.y),particles.X(face.z));
        RANGE<TV> box=RANGE<TV>::Bounding_Box(particles.X(face.x),particles.X(face.y),particles.X(face.z));
        RANGE<TV_INT> cell_range(Clamp_Cell(box.min_corner),Clamp_Cell(box.max_corner));
        for(UNIFORM_GRID_ITERATOR_CELL<TV> iter(grid,cell_range);iter.Valid();iter.Next()){
            ARRAY<TV_INT>& list=cell_to_faces(iter.Cell_Index());
            for(int i=1;i<=list.m;i++){
                const TV_INT& face2=list(i);
                if(visited_faces.Contains(face2)) continue;
                visited_faces.Set(face2);
                if(face.x==face2.x || face.x==face2.y  || face.x==face2.z || face.y==face2.x || face.y==face2.y || face.y==face2.z || face.z==face2.x || face.z==face2.y || face.z==face2.z) continue;
                TRIANGLE_3D<T> triangle2(particles.X(face2.x),particles.X(face2.y),particles.X(face2.z));
                if(INTERSECTION::Intersects(triangle1,triangle2)) return true;}}
        return false;
    }

    bool Element_Collides_With_Front(VECTOR<int,4> element)
    {
        TETRAHEDRON<T> tetrahedron(particles.X(element(1)),particles.X(element(2)),particles.X(element(3)),particles.X(element(4)));
        RANGE<TV> box=RANGE<TV>::Bounding_Box(particles.X(element(1)),particles.X(element(2)),particles.X(element(3)),particles.X(element(4)));
        RANGE<TV_INT> cell_range(Clamp_Cell(box.min_corner),Clamp_Cell(box.max_corner));
        for(UNIFORM_GRID_ITERATOR_CELL<TV> iter(grid,cell_range);iter.Valid();iter.Next()){
            ARRAY<int>& list=cell_particles(iter.Cell_Index());
            for(int i=1;i<=list.m;i++){
                int particle_index=list(i);
                if(!particle_on_boundary(particle_index)) continue;
                if(element(1)==particle_index || element(2)==particle_index || element(3)==particle_index || element(4)==particle_index) continue;
                if(tetrahedron.Inside(particles.X(particle_index))) return true;}}
        return Face_Collides_With_Front(TV_INT(element(1),element(2),element(4))) || Face_Collides_With_Front(TV_INT(element(2),element(3),element(4))) || Face_Collides_With_Front(TV_INT(element(3),element(1),element(4)));
    }

    VECTOR<int,4> Find_Initial_Tetrahedron()
    {
        for(UNIFORM_GRID_ITERATOR_CELL<TV> iter(grid,ghost_cells);iter.Valid();iter.Next()){
            TV_INT cell=iter.Cell_Index();
            for(int i=1;i<=cell_particles(cell).m;i++){int p0=cell_particles(cell)(i);
                if(particle_finished(p0))continue;
                int p1=Find_Nearest_Point(p0);
                if(p1>0){
                    int p2=Find_Smallest_Circumcircle_Area_Point(p0,p1);
                    if(p2>0){
                        int p3=Find_Smallest_Circumsphere_Volume_Point(p0,p1,p2);
                        if(p3>0 && Is_Face_Valid(TV_INT(p0,p1,p2)) && Is_Face_Valid(TV_INT(p0,p2,p3)) && Is_Face_Valid(TV_INT(p0,p1,p3)) && Is_Face_Valid(TV_INT(p1,p2,p3)))
                            return VECTOR<int,4>(p0,p1,p2,p3);
                        p3=Find_Smallest_Circumsphere_Volume_Point(p1,p0,p2);
                        if(p3>0 && Is_Face_Valid(TV_INT(p0,p1,p2)) && Is_Face_Valid(TV_INT(p0,p2,p3)) && Is_Face_Valid(TV_INT(p0,p1,p3)) && Is_Face_Valid(TV_INT(p1,p2,p3)))
                            return VECTOR<int,4>(p1,p0,p2,p3);}}}}
        return VECTOR<int,4>();
    }

    T Circumcircle_Radius(TV x1,TV x2,TV x3)
    {
        TV x12=x2-x1,x13=x3-x1,x23=x3-x2;
        return x23.Magnitude()*x12.Magnitude()*x13.Magnitude()*(T)0.5/TV::Cross_Product(x12,x13).Magnitude();
    }

    void Process_Face(TV xx1,TV xx2,TV xx3,T max_radius,/*result*/TV &face_normal,/*result*/TV &circumcircle_center,/*result*/TV &max_circumsphere_center)
    {
        T a1=TV::Angle_Between(xx2-xx1,xx3-xx1);T a2=TV::Angle_Between(xx1-xx2,xx3-xx2);T a3=TV::Angle_Between(xx2-xx3,xx1-xx3);
        TV x1,x2,x3;
        if(a1<=a2&&a1<=a3){x3=xx1;x1=xx2;x2=xx3;}
        else if(a2<=a1&&a2<=a3){x3=xx2;x1=xx3;x2=xx1;}
        else{x3=xx3;x1=xx1;x2=xx2;}
        T circumcircle_r=Circumcircle_Radius(x1,x2,x3);
        face_normal=TV::Cross_Product(x2-x1,x3-x1).Normalized();
        circumcircle_center=(x1+x2)*0.5+TV::Cross_Product(face_normal,x2-x1).Normalized()*sqrt(abs(circumcircle_r*circumcircle_r-(x2-x1).Magnitude_Squared()*0.25));
        max_circumsphere_center=circumcircle_center+face_normal*sqrt(max_radius*max_radius-circumcircle_r*circumcircle_r);
    }

    static TV Get_Max_Circumsphere_Center(TV x1,TV x2,TV x3,T max_radius)
    {
        TV circumcircle_center=TRIANGLE_3D<T>::Circumcenter(x1,x2,x3);
        T circumcircle_r=(x1-circumcircle_center).Magnitude();
        TV face_normal=(TV::Cross_Product(x2-x1,x3-x1)).Normalized();
        return circumcircle_center+face_normal*sqrt(abs(max_radius*max_radius-circumcircle_r*circumcircle_r));
    }

    int Find_Empty_Circumsphere_Point_With_Grid_Radius(TV_INT face,int type=1)
    {
        TV x1=particles.X(face.x),x2=particles.X(face.y),x3=particles.X(face.z);
        TV max_circumsphere_center=Get_Max_Circumsphere_Center(x1,x2,x3,grid_radius);
        TV_INT bk_cell=Clamp_Cell(max_circumsphere_center);
        RANGE<TV_INT> range(bk_cell-grid_r_steps_in_bkgrid,bk_cell+grid_r_steps_in_bkgrid);
        int minimum_p=0;
        for(UNIFORM_GRID_ITERATOR_CELL<TV> iter(grid,range);iter.Valid();iter.Next())if(grid.Inside_Domain(iter.Cell_Index(),ghost_cells)){
            for(int i=1;i<=cell_particles(iter.Cell_Index()).m;i++){int p=cell_particles(iter.Cell_Index())(i);
                if(p==face.x || p==face.y || p==face.z) continue;
                if((particles.X(p)-max_circumsphere_center).Magnitude_Squared()>grid_radius_squared*1.1) continue;
                if(ADAPTIVE_EXACT_PREDICATES::orient3d(VECTOR<double,3>(x1).begin(),VECTOR<double,3>(x2).begin(),VECTOR<double,3>(x3).begin(),VECTOR<double,3>(particles.X(p)).begin())>=0) continue;
                if(minimum_p!=0) if(!InSphere(face.x,face.y,face.z,minimum_p,p)) continue;
                minimum_p=p;}}
        if(minimum_p!=0){
            if(particle_finished(minimum_p)) return 0;
            if(use_alpha_shape && type==1 && (particles.X(minimum_p)-max_circumsphere_center).Magnitude_Squared()>grid_radius_squared) return 0;}
        return minimum_p;
    }

    int Find_Nearest_Point(const int index)
    {
        TV x=particles.X(index);
        TV_INT bk_cell=Clamp_Cell(x);
        RANGE<TV_INT> range(bk_cell-grid_r_steps_in_bkgrid,bk_cell+grid_r_steps_in_bkgrid);
        int minimum_p=0;T minimum_dd=FLT_MAX;
        for(UNIFORM_GRID_ITERATOR_CELL<TV> iter(grid,range);iter.Valid();iter.Next())if(grid.Inside_Domain(iter.Cell_Index(),ghost_cells)){
            const ARRAY<int>& list=cell_particles(iter.Cell_Index());
            for(int i=1;i<=list.m;i++){int p=list(i);
                if(p==index || particle_finished(p)) continue;
                TV y=particles.X(p);
                T dd=(x-y).Magnitude();
                if(dd<minimum_dd && dd<2*grid_radius){
                    minimum_dd=dd;
                    minimum_p=p;}}}
        return minimum_p;
    }

    int Find_Smallest_Circumcircle_Area_Point(int i1,int i2)
    {
        TV x1=particles.X(i1),x2=particles.X(i2);
        RANGE<TV> bbox=RANGE<TV>::Bounding_Box(x1,x2);
        TV_INT cell1=Clamp_Cell(bbox.min_corner),cell2=Clamp_Cell(bbox.max_corner);
        RANGE<TV_INT> range(cell1-grid_r_steps_in_bkgrid,cell2+grid_r_steps_in_bkgrid);
        int min_p=0;T min_dd=FLT_MAX;
        for(UNIFORM_GRID_ITERATOR_CELL<TV> itr(grid,range);itr.Valid();itr.Next())
            if(grid.Inside_Domain(itr.Cell_Index(),ghost_cells)){
                const ARRAY<int>& list=cell_particles(itr.Cell_Index());
                for(int i=1;i<=list.m;i++){
                    int p=list(i);
                    if(p==i1||p==i2||particle_finished(p)) continue;
                    TV x3=particles.X(p);
                    T dd=Circumcircle_Radius(x1,x2,x3);
                    if(dd<min_dd && dd<2*grid_radius){
                        min_dd=dd;
                        min_p=p;}}}
        return min_p;
    }

    int Find_Smallest_Circumsphere_Volume_Point(int i1,int i2,int i3)
    {
        return Find_Empty_Circumsphere_Point_With_Grid_Radius(TV_INT(i1,i2,i3));
    }
};

}
#endif
