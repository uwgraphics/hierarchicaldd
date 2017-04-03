//#####################################################################
// Copyright 2012, Mridul Aanjaneya, Saket Patkar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Geometry/Grids_Uniform_PDE_Linear/BUBBLES_POISSON.h>
using namespace PhysBAM;
//#####################################################################
// Function Solve
//#####################################################################
template<class T_GRID> void BUBBLES_POISSON<T_GRID>::
Solve(const T time,const bool solution_regions_already_computed)
{
    if(!solution_regions_already_computed) Find_Solution_Regions();
#ifdef USE_PTHREADS
    if(thread_queue){
        pthread_mutex_init(&lock,0);
        pthread_barrier_init(&barr,0,thread_queue->Number_Of_Threads());
        pcg_threaded->maximum_iterations=pcg.maximum_iterations;
        pcg_threaded->number_of_threads=thread_queue->Number_Of_Threads();pthread_barrier_init(&pcg_threaded->barr,0,pcg_threaded->number_of_threads);}
#endif
    // matrix to cell index is currently resized to 1 row per region
    ARRAY<ARRAY<TV_INT> > matrix_index_to_cell_index_array(number_of_regions);T_ARRAYS_INT cell_index_to_matrix_index(grid.Domain_Indices(1));
    ARRAY<int,VECTOR<int,1> > filled_region_cell_count(-1,number_of_regions),incompressible_region_cell_count(-1,number_of_regions);
    // Array is resized to number of regions
    ARRAY<SPARSE_MATRIX_FLAT_NXN<T> > A_array(number_of_regions);ARRAY<VECTOR_ND<T> > b_array(number_of_regions);
    ARRAY<int> belongs_to_region(number_of_bubble_regions);     // the region that the bubble belongs to
    // count number of cells in each color
    for(CELL_ITERATOR iterator(grid,1);iterator.Valid();iterator.Next()){TV_INT cell=iterator.Cell_Index();
        if(grid.Inside_Domain(cell) && phase_color_ghost(cell)==1) belongs_to_region(bubble_region_colors(cell))=filled_region_colors(cell);
        else filled_region_cell_count(filled_region_colors(cell))++;}
    for(int color=1;color<=number_of_bubble_regions;color++){
        if(phase_color(color)==1)filled_region_cell_count(belongs_to_region(color))++;}
    // for regions that we would solve resize matrix_index_to_cell_index_array to number of cells in that region
    for(int color=1;color<=number_of_regions;color++) if(filled_region_touches_dirichlet(color)||solve_neumann_regions){
        matrix_index_to_cell_index_array(color).Resize(filled_region_cell_count(color));}
    filled_region_cell_count.Fill(0);   // reusing this array in order to make the indirection arrays
    incompressible_region_cell_count.Fill(0);
    DOMAIN_ITERATOR_THREADED_ALPHA<LAPLACE_UNIFORM<T_GRID>,TV> threaded_iterator(grid.Domain_Indices(1),thread_queue,1,1,2,1);
    ARRAY<int,TV_INT> domain_index(grid.Domain_Indices(1),false);
    // identify which domain each cell belongs to (when threading breaks the grid into many domains)
    for(int i=1;i<=threaded_iterator.domains.m;i++){
        RANGE<TV_INT> interior_domain(threaded_iterator.domains(i));interior_domain.max_corner-=TV_INT::All_Ones_Vector();interior_domain.min_corner+=TV_INT::All_Ones_Vector();
        for(CELL_ITERATOR iterator(grid,interior_domain);iterator.Valid();iterator.Next()) domain_index(iterator.Cell_Index())=i;}
    ARRAY<ARRAY<INTERVAL<int> > > interior_indices(number_of_regions);
    ARRAY<ARRAY<ARRAY<INTERVAL<int> > > > ghost_indices(number_of_regions);
    for(int color=1;color<=number_of_regions;color++){
        interior_indices(color).Resize(threaded_iterator.number_of_domains);ghost_indices(color).Resize(threaded_iterator.number_of_domains);
        for(int i=1;i<=threaded_iterator.domains.m;i++) ghost_indices(color)(i).Resize(2*TV::dimension);}
    // fill the matrix index to cell index & cell index to matrix index matrices    
    if(!mpi_grid && !thread_queue) Compute_Matrix_Indices(grid.Domain_Indices(1),filled_region_cell_count,incompressible_region_cell_count,matrix_index_to_cell_index_array,cell_index_to_matrix_index);
    else if(!thread_queue) laplace_bubble_mpi->Find_Matrix_Indices_Bubbles(filled_region_cell_count,incompressible_region_cell_count,cell_index_to_matrix_index,matrix_index_to_cell_index_array,this);

    // fill the ghost bubble indices
    laplace_bubble_mpi->ghost_bubble_indices.Delete_Pointers_And_Clean_Memory();
    laplace_bubble_mpi->ghost_bubble_indices.Resize(ghost_bubble_cells.Size());
    for(int i=1;i<=ghost_bubble_cells.Size();i++){
        if(ghost_bubble_cells(i) && ghost_bubble_cells(i)->Size()>0){ARRAY<TV_INT>& ghost_cells=*ghost_bubble_cells(i);
            laplace_bubble_mpi->ghost_bubble_indices(i)=new ARRAY<int>(ghost_bubble_cells(i)->Size());
            ARRAY<int>& ghost_indices=*laplace_bubble_mpi->ghost_bubble_indices(i);
            for(int j=1;j<=ghost_bubble_cells(i)->Size();j++) ghost_indices(j)=cell_index_to_matrix_index(ghost_cells(j));}}

    RANGE<TV_INT> domain=grid.Domain_Indices(1);
    // create the actual matrix
    Find_A(domain,A_array,b_array,filled_region_cell_count,cell_index_to_matrix_index);

    // add the compressible cell entries in the matrix.
    for(CELL_ITERATOR iterator(grid,1);iterator.Valid();iterator.Next()){INDEX cell_index=iterator.Cell_Index();
        int color=filled_region_colors(cell_index);
        if(color!=-1 && (filled_region_touches_dirichlet(color)||solve_neumann_regions)){
            int matrix_index=cell_index_to_matrix_index(cell_index);
            SPARSE_MATRIX_FLAT_NXN<T>& A=A_array(filled_region_colors(cell_index));
            if(phase_color_ghost(cell_index)==1) A(matrix_index,matrix_index)-=(T)1./(B*cell_density(cell_index)*dt*dt);}}

    for(int color=1;color<=number_of_regions;color++) if(filled_region_cell_count(color)>0 && (filled_region_touches_dirichlet(color)||solve_neumann_regions)){
        pcg.Enforce_Compatibility(!filled_region_touches_dirichlet(color)&&enforce_compatibility);
        Solve_Subregion(interior_indices(color),ghost_indices(color),matrix_index_to_cell_index_array(color),A_array(color),b_array(color),color,&domain_index,neighbor_matrix_indices,bubble_index);}
    if(!solve_neumann_regions) for(CELL_ITERATOR iterator(grid,1);iterator.Valid();iterator.Next()){
        int filled_region_color=filled_region_colors(iterator.Cell_Index());if(filled_region_color>0 && !filled_region_touches_dirichlet(filled_region_color)) u(iterator.Cell_Index())=0;}

    Reconstruct_Compressible_Pressure(domain);
}
//#####################################################################
// Function Solve_Subregion
//#####################################################################
template<class T_GRID> void BUBBLES_POISSON<T_GRID>::
Solve_Subregion(ARRAY<INTERVAL<int> >& interior_indices,ARRAY<ARRAY<INTERVAL<int> > >& ghost_indices,ARRAY<TV_INT>& matrix_index_to_cell_index,SPARSE_MATRIX_FLAT_NXN<T>& A,VECTOR_ND<T>& b,const int color,ARRAY<int,TV_INT>* domain_index,ARRAY<ARRAY<int> >& neighbor_matrix_indices,ARRAY<int>& bubble_index)
{
    int number_of_unknowns=matrix_index_to_cell_index.m;
    A.Negate();b*=(T)-1;
    VECTOR_ND<T> x(number_of_unknowns),q,s,r,k,z;
    for(int i=1;i<=number_of_unknowns;i++) x(i)=u(matrix_index_to_cell_index(i));
    Find_Tolerance(b); // needs to happen after b is completely set up
    if(pcg.show_results){std::stringstream ss;ss<< "solving " << number_of_unknowns << " cells to tolerance " << tolerance << std::endl;LOG::filecout(ss.str());}
    DOMAIN_ITERATOR_THREADED_ALPHA<PCG_SPARSE_THREADED<TV>,TV> threaded_iterator(grid.Domain_Indices(1),thread_queue,1,1,2,1);
    static const int min_unknowns_for_threading=100;
    bool use_threaded_solve=thread_queue&&number_of_unknowns>=min_unknowns_for_threading;
    if(!mpi_grid){
        if(use_threaded_solve){pcg_threaded->p.Resize(A.n,false);pcg_threaded->temp.Resize(A.n,false);}
        if(use_threaded_solve) threaded_iterator.template Run<const ARRAY<int,TV_INT>&,const ARRAY<INTERVAL<int> >&,const ARRAY<ARRAY<INTERVAL<int> > >&,SPARSE_MATRIX_FLAT_NXN<T>&,VECTOR_ND<T>&,VECTOR_ND<T>&,T>(*pcg_threaded,&PCG_SPARSE_THREADED<TV>::Solve,*domain_index,interior_indices,ghost_indices,A,x,b,tolerance);
        //if(use_threaded_solve) pcg_threaded->Solve_In_Parts(threaded_iterator,domain_index,interior_indices,ghost_indices,A,x,b,tolerance);
        //if(use_threaded_solve) pcg_threaded->Solve_In_Parts(A,x,b,tolerance);
        else pcg.Solve(A,x,b,q,s,r,k,z,tolerance);}
    else{
        if(use_threaded_solve) DOMAIN_ITERATOR_THREADED_ALPHA<LAPLACE_MPI<T_GRID>,TV>(grid.Domain_Indices(1),thread_queue,1,1,2,1).template Run<const ARRAY<int,TV_INT>&,ARRAY<INTERVAL<int> >&,ARRAY<ARRAY<INTERVAL<int> > >&,SPARSE_MATRIX_FLAT_NXN<T>&,VECTOR_ND<T>&,VECTOR_ND<T>&,VECTOR_ND<T>&,VECTOR_ND<T>&,VECTOR_ND<T>&,VECTOR_ND<T>&,VECTOR_ND<T>&,T,int,int>(*laplace_mpi,&LAPLACE_MPI<T_GRID>::Solve_Threaded,*domain_index,interior_indices,ghost_indices,A,x,b,q,s,r,k,z,tolerance,color,1);
        else laplace_bubble_mpi->Solve(A,x,b,q,s,r,k,z,tolerance,color,neighbor_matrix_indices,bubble_index);}
    for(int i=1;i<=number_of_unknowns;i++){TV_INT cell_index=matrix_index_to_cell_index(i);u(cell_index)=x(i);}
}
//#####################################################################
// Function Reconstruct_Compressible_Pressure
//#####################################################################
template<class T_GRID> void BUBBLES_POISSON<T_GRID>::
Reconstruct_Compressible_Pressure(RANGE<TV_INT>& domain)
{
    ARRAY<T> first_bubble_pressure(number_of_bubble_regions);
    ARRAY<bool> region_encountered(number_of_bubble_regions);region_encountered.Fill(false);
    
    for(CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){TV_INT cell=iterator.Cell_Index();
        int color=filled_region_colors(cell),bubble_color=bubble_region_colors(cell);
        if(color>0&&(filled_region_touches_dirichlet(color)||solve_neumann_regions)){
            if(phase_color_ghost(cell)==1){
                if(!region_encountered(bubble_color)){first_bubble_pressure(bubble_color)=u(cell);region_encountered(bubble_color)=true;}
                else u(cell)=first_bubble_pressure(bubble_color);}}}
}
//#####################################################################
// Function Compute_Matrix_Indices
//#####################################################################
template<class T_GRID> void BUBBLES_POISSON<T_GRID>::
Compute_Matrix_Indices(ARRAY<int,VECTOR<int,1> >& filled_region_cell_count,ARRAY<int,VECTOR<int,1> >& incompressible_region_cell_count,ARRAY<ARRAY<TV_INT> >& matrix_index_to_cell_index_array,T_ARRAYS_INT& cell_index_to_matrix_index)
{
    RANGE<TV_INT> domain(grid.Domain_Indices(1));
    Compute_Matrix_Indices(domain,filled_region_cell_count,incompressible_region_cell_count,matrix_index_to_cell_index_array,cell_index_to_matrix_index);
}
//#####################################################################
// Function Compute_Matrix_Indices
//#####################################################################
template<class T_GRID> void BUBBLES_POISSON<T_GRID>::
Compute_Matrix_Indices(const RANGE<TV_INT>& domain,ARRAY<int,VECTOR<int,1> >& filled_region_cell_count,ARRAY<int,VECTOR<int,1> >& incompressible_region_cell_count,ARRAY<ARRAY<TV_INT> >& matrix_index_to_cell_index_array,T_ARRAYS_INT& cell_index_to_matrix_index)
{
    ARRAY<TV_INT> first_bubble_cell(number_of_bubble_regions);ARRAY<int> belongs_to_region(number_of_bubble_regions);
    ARRAY<bool> region_encountered(number_of_bubble_regions);region_encountered.Fill(false);
   
    for(CELL_ITERATOR iterator(grid,domain);iterator.Valid();iterator.Next()){
        TV_INT cell=iterator.Cell_Index();int color=filled_region_colors(cell);
        if(color>0&&(filled_region_touches_dirichlet(color)||solve_neumann_regions)){
            if(grid.Inside_Domain(cell) && phase_color_ghost(cell)==1){int bubble_color=bubble_region_colors(cell);
                belongs_to_region(bubble_color)=filled_region_colors(cell);
                if(!region_encountered(bubble_color)){first_bubble_cell(bubble_color)=cell;region_encountered(bubble_color)=true;}}
            else{cell_index_to_matrix_index(cell)=++filled_region_cell_count(color);
                matrix_index_to_cell_index_array(color)(filled_region_cell_count(color))=cell;}}}

    for(int i=1;i<=number_of_regions;i++) incompressible_region_cell_count(i)=filled_region_cell_count(i);

    for(int color=number_of_bubble_regions;color>=1;color--){
        if(phase_color(color)==1 && region_encountered(color)){filled_region_cell_count(belongs_to_region(color))++;
            if(laplace_bubble_mpi->mpi_bubbles && color>laplace_bubble_mpi->mpi_bubbles->number_of_global_regions) incompressible_region_cell_count(belongs_to_region(color))++;
            matrix_index_to_cell_index_array(belongs_to_region(color))(filled_region_cell_count(belongs_to_region(color)))=first_bubble_cell(color);
            cell_index_to_matrix_index(first_bubble_cell(color))=filled_region_cell_count(belongs_to_region(color));}}

    for(CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){TV_INT cell=iterator.Cell_Index();
        int color=filled_region_colors(cell),bubble_color=bubble_region_colors(cell);
        if(color>0&&(filled_region_touches_dirichlet(color)||solve_neumann_regions)){
            if(phase_color_ghost(cell)==1 && region_encountered(bubble_color)){
                cell_index_to_matrix_index(cell)=cell_index_to_matrix_index(first_bubble_cell(bubble_color));}}}

    // fill the matrix indices for all global bubbles
    if(mpi_grid){bubble_index.Resize(number_of_bubble_regions);
        for(int color=1;color<=number_of_bubble_regions;color++)
            if(region_encountered(color) && phase_color(color)==1) bubble_index(color)=cell_index_to_matrix_index(first_bubble_cell(color));}

    // fill the neighbor matrix indices
    for(int i=1;i<=neighbor_cell_indices.Size();i++){
        for(int j=1;j<=neighbor_cell_indices(i).Size();j++){
            neighbor_matrix_indices(i)(j)=cell_index_to_matrix_index(neighbor_cell_indices(i)(j));}}
}
//#####################################################################
// Find_A_Part_One
//#####################################################################
template<class T_GRID> void BUBBLES_POISSON<T_GRID>::
Find_A_Part_One(RANGE<TV_INT>& domain,T_ARRAYS_INT& cell_index_to_matrix_index,ARRAY<ARRAY<int> >& row_counts)
{
    // compute neighbor cells for individual bubbles
    ARRAY<int> incompressible_neighbors(number_of_bubble_regions);incompressible_neighbors.Fill(0);
    T_ARRAYS_HASH_INT visited(grid.Domain_Indices(1));

    for(CELL_ITERATOR iterator(grid,1);iterator.Valid();iterator.Next()){TV_INT cell=iterator.Cell_Index();
        if(grid.Inside_Domain(cell)){int color=filled_region_colors(cell);      // all collapsed bubble cells are inside the grid
            for(int axis=1;axis<=T_GRID::dimension;axis++){TV_INT offset;offset[axis]=1;
                TV_INT left_cell=cell-offset,right_cell=cell+offset;int bubble_color=bubble_region_colors(cell);
                bool is_bubble=(phase_color_ghost(cell)==1);

                bool left_cell_ghost=!grid.Inside_Domain(left_cell),is_left_water=(phase_color_ghost(left_cell)==0);
                if(is_bubble && ((!left_cell_ghost && is_left_water)||left_cell_ghost) && !visited(left_cell).Contains(bubble_color)){
                    if((filled_region_colors(left_cell)==color || (!grid.Inside_Domain(left_cell) && periodic_boundary[axis])) && !psi_N.Component(axis)(cell)){ 
                        incompressible_neighbors(bubble_color)++;visited(left_cell).Insert(bubble_color);}}

                bool right_cell_ghost=!grid.Inside_Domain(right_cell),is_right_water=(phase_color_ghost(right_cell)==0);
                if(is_bubble && ((!right_cell_ghost && is_right_water)||right_cell_ghost) && !visited(right_cell).Contains(bubble_color)){
                    if((filled_region_colors(right_cell)==color || (!grid.Inside_Domain(right_cell) && periodic_boundary[axis])) && !psi_N.Component(axis)(cell+offset)){ 
                        incompressible_neighbors(bubble_color)++;visited(right_cell).Insert(bubble_color);}}}}}

    // TODO: this should be rewritten in terms of faces cause this got really hacky with MPI
    for(CELL_ITERATOR iterator(grid,domain);iterator.Valid();iterator.Next()){TV_INT cell_index=iterator.Cell_Index();
        int color=filled_region_colors(cell_index);assert(color!=0);

        if(!grid.Inside_Domain(cell_index) && color!=-1 && (filled_region_touches_dirichlet(color)||solve_neumann_regions)){int row_count=1;
            for(int axis=1;axis<=T_GRID::dimension;axis++){TV_INT offset;offset[axis]=1;
                if(((filled_region_colors.Valid_Index(cell_index-offset) && filled_region_colors(cell_index-offset)==color) ||
                    (grid.Domain_Indices().Lazy_Outside(cell_index-offset) && periodic_boundary[axis])) && !psi_N.Component(axis)(cell_index)) row_count++;
                if(((filled_region_colors.Valid_Index(cell_index+offset) && filled_region_colors(cell_index+offset)==color) ||
                    (grid.Domain_Indices().Lazy_Outside(cell_index+offset) && periodic_boundary[axis])) && !psi_N.Component(axis)(cell_index+offset)) row_count++;}
            row_counts(color)(cell_index_to_matrix_index(cell_index))=row_count;}
        else{int bubble_color=bubble_region_colors(cell_index);
            if(phase_color_ghost(cell_index)==0){
                if(color!=-1 && (filled_region_touches_dirichlet(color)||solve_neumann_regions)){
                    //counter for valid neighbouring cells for each cell, can be dim*2+1 at max
                    int row_count=1;HASHTABLE<int> bubble_colors_visited;
                    for(int axis=1;axis<=T_GRID::dimension;axis++){TV_INT offset;offset[axis]=1;
                        if(((filled_region_colors.Valid_Index(cell_index-offset) && filled_region_colors(cell_index-offset)==color) ||
                            (grid.Domain_Indices().Lazy_Outside(cell_index-offset) && periodic_boundary[axis])) && !psi_N.Component(axis)(cell_index)){
                            int bubble_color=bubble_region_colors(cell_index-offset);

                            if(!grid.Inside_Domain(cell_index-offset) || phase_color_ghost(cell_index-offset)==0 || bubble_color<=0) row_count++;
                            else if(phase_color_ghost(cell_index-offset)==1 && !bubble_colors_visited.Contains(bubble_color)){row_count++;bubble_colors_visited.Insert(bubble_color);}}
                        if(((filled_region_colors.Valid_Index(cell_index+offset) && filled_region_colors(cell_index+offset)==color) ||
                            (grid.Domain_Indices().Lazy_Outside(cell_index+offset) && periodic_boundary[axis])) && !psi_N.Component(axis)(cell_index+offset)){
                            int bubble_color=bubble_region_colors(cell_index+offset);

                            if(!grid.Inside_Domain(cell_index+offset) || phase_color_ghost(cell_index+offset)==0 || bubble_color<=0) row_count++;
                            else if(phase_color_ghost(cell_index+offset)==1 && !bubble_colors_visited.Contains(bubble_color)){row_count++;bubble_colors_visited.Insert(bubble_color);}}}
                    //set the row size of the row the cell corresponds to 
                    row_counts(color)(cell_index_to_matrix_index(cell_index))=row_count;}}
            else if(color!=-1 && (filled_region_touches_dirichlet(color)||solve_neumann_regions) && phase_color_ghost(cell_index)==1)
                row_counts(color)(cell_index_to_matrix_index(cell_index))=incompressible_neighbors(bubble_color)+1;}}
}
//#####################################################################
// Function Find_A_Part_Two
//#####################################################################
template<class T_GRID> void BUBBLES_POISSON<T_GRID>::
Find_A_Part_Two(RANGE<TV_INT>& domain,ARRAY<SPARSE_MATRIX_FLAT_NXN<T> >& A_array,ARRAY<VECTOR_ND<T> >& b_array,T_ARRAYS_INT& cell_index_to_matrix_index)
{
    TV one_over_dx2=Inverse(grid.dX*grid.dX);TV_INT grid_counts=grid.counts;

    if(use_weighted_divergence)
        for(CELL_ITERATOR iterator(grid,domain);iterator.Valid();iterator.Next()){
            int color=filled_region_colors(iterator.Cell_Index());
            if(color!=-1 && (filled_region_touches_dirichlet(color)||solve_neumann_regions)){const TV_INT& cell_index=iterator.Cell_Index();
                //for every valid cell
                int matrix_index=cell_index_to_matrix_index(cell_index);
                SPARSE_MATRIX_FLAT_NXN<T>& A=A_array(filled_region_colors(cell_index));VECTOR_ND<T>& b=b_array(filled_region_colors(cell_index));b(matrix_index)+=f(cell_index);
                T diagonal=0;
                //for all neighbours
                for(int axis=1;axis<=T_GRID::dimension;axis++){TV_INT offset=TV_INT::Axis_Vector(axis);
                    if(filled_region_colors.Valid_Index(cell_index-offset)){
                        if(!psi_N.Component(axis)(cell_index)){
                            T element=divergence_face_weights.Component(axis)(cell_index)*beta_face.Component(axis)(cell_index)*one_over_dx2[axis];
                            diagonal-=element;
                            if(grid.Domain_Indices().Lazy_Outside(cell_index-offset) && periodic_boundary[axis]){
                                TV_INT periodic_offset_cell=cell_index-offset;
                                int axis_periodic_cell=1+wrap(periodic_offset_cell[axis]-1,grid_counts[axis]);
                                periodic_offset_cell[axis]=axis_periodic_cell;
                                A.Add_Element(matrix_index,cell_index_to_matrix_index(periodic_offset_cell),element);}
                            else if(psi_D(cell_index-offset)) b(matrix_index)-=element*u(cell_index-offset);
                            else A.Add_Element(matrix_index,cell_index_to_matrix_index(cell_index-offset),element);}}
                    if(filled_region_colors.Valid_Index(cell_index+offset)){
                        if(!psi_N.Component(axis)(cell_index+offset)){
                            T element=divergence_face_weights.Component(axis)(cell_index+offset)*beta_face.Component(axis)(cell_index+offset)*one_over_dx2[axis];
                            diagonal-=element;
                            if(grid.Domain_Indices().Lazy_Outside(cell_index+offset) && periodic_boundary[axis]){
                                TV_INT periodic_offset_cell=cell_index+offset;
                                int axis_periodic_cell=1+wrap(periodic_offset_cell[axis]-1,grid_counts[axis]);
                                periodic_offset_cell[axis]=axis_periodic_cell;
                                A.Add_Element(matrix_index,cell_index_to_matrix_index(periodic_offset_cell),element);}
                            else if(psi_D(cell_index+offset)) b(matrix_index)-=element*u(cell_index+offset);
                            else A.Add_Element(matrix_index,cell_index_to_matrix_index(cell_index+offset),element);}}}
                A.Add_Element(matrix_index,matrix_index,diagonal);}} // set diagonal and right hand side
    else
        for(CELL_ITERATOR iterator(grid,domain);iterator.Valid();iterator.Next()){
            int color=filled_region_colors(iterator.Cell_Index());
            if(color!=-1 && (filled_region_touches_dirichlet(color)||solve_neumann_regions)){const TV_INT& cell_index=iterator.Cell_Index();
                int matrix_index=cell_index_to_matrix_index(cell_index);
                SPARSE_MATRIX_FLAT_NXN<T>& A=A_array(filled_region_colors(cell_index));VECTOR_ND<T>& b=b_array(filled_region_colors(cell_index));b(matrix_index)+=f(cell_index);
                T diagonal=0;
                for(int axis=1;axis<=T_GRID::dimension;axis++){TV_INT offset=TV_INT::Axis_Vector(axis);
                    if(filled_region_colors.Valid_Index(cell_index-offset) && (filled_region_colors(cell_index-offset)!=-1 || psi_D(cell_index-offset))){
                        if(!psi_N.Component(axis)(cell_index)){T element=beta_face.Component(axis)(cell_index)*one_over_dx2[axis];
                            diagonal-=element;
                            if(grid.Domain_Indices().Lazy_Outside(cell_index-offset) && periodic_boundary[axis]){
                                TV_INT periodic_offset_cell=cell_index-offset;
                                int axis_periodic_cell=1+wrap(periodic_offset_cell[axis]-1,grid_counts[axis]);
                                periodic_offset_cell[axis]=axis_periodic_cell;
                                A.Add_Element(matrix_index,cell_index_to_matrix_index(periodic_offset_cell),element);}
                            else if(psi_D(cell_index-offset)) b(matrix_index)-=element*u(cell_index-offset);
                            else A.Add_Element(matrix_index,cell_index_to_matrix_index(cell_index-offset),element);}}
                    if(filled_region_colors.Valid_Index(cell_index+offset) && (filled_region_colors(cell_index+offset)!=-1 || psi_D(cell_index+offset))){
                        if(!psi_N.Component(axis)(cell_index+offset)){T element=beta_face.Component(axis)(cell_index+offset)*one_over_dx2[axis];
                            diagonal-=element;
                            if(grid.Domain_Indices().Lazy_Outside(cell_index+offset) && periodic_boundary[axis]){
                                TV_INT periodic_offset_cell=cell_index+offset;
                                int axis_periodic_cell=1+wrap(periodic_offset_cell[axis]-1,grid_counts[axis]);
                                periodic_offset_cell[axis]=axis_periodic_cell;
                                A.Add_Element(matrix_index,cell_index_to_matrix_index(periodic_offset_cell),element);}
                            else if(psi_D(cell_index+offset)) b(matrix_index)-=element*u(cell_index+offset);
                            else A.Add_Element(matrix_index,cell_index_to_matrix_index(cell_index+offset),element);}}}
                A.Add_Element(matrix_index,matrix_index,diagonal);}} // set diagonal and right hand side
    if(second_order_cut_cell_method) BASE::Apply_Second_Order_Cut_Cell_Method(domain,A_array,b_array,cell_index_to_matrix_index);
}
template class BUBBLES_POISSON<GRID<VECTOR<float,1> > >;
template class BUBBLES_POISSON<GRID<VECTOR<float,2> > >;
template class BUBBLES_POISSON<GRID<VECTOR<float,3> > >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class BUBBLES_POISSON<GRID<VECTOR<double,1> > >;
template class BUBBLES_POISSON<GRID<VECTOR<double,2> > >;
template class BUBBLES_POISSON<GRID<VECTOR<double,3> > >;
#endif
