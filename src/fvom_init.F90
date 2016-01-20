!> @author
!> Sergei Danilov
!> @brief
!> fvom_init.F90
!! A set of routines that does domain partitioning for a distributed setup.
!! Cell-vertex finite-volume discretization
!! routines used here are the stripped down versions of the main setup.
!! Do no mix them up!

!=============================================================================
!> @brief
!> Main driver routine for initialization
program MAIN
  use o_PARAM
  use o_MESH
  use g_PARSUP
  use g_CONFIG
  use g_rotate_grid
  implicit none
  character(len=100)   :: nmlfile  !> name of configuration namelist file
  
  nmlfile ='../config/namelist.config'
  open (20,file=nmlfile)
  read (20,NML=paths)         ! We need MeshPath
  read (20,NML=geometry)      ! We need cyclic_length and cartesian
  close (20)
  cyclic_length=cyclic_length*rad 
  call set_mesh_transform_matrix  !(rotated grid) 
  call par_init
  call read_mesh_ini
  call test_tri_ini
  call find_edges_ini
  call find_elem_neighbors_ini
  call find_levels

  if (mype==0) call stiff_mat_ini  ! Only the master will call Metis and needs the matrix
  call set_par_support_ini
  call save_dist_mesh
  call par_ex
end program MAIN
!=============================================================================
!> @brief
!> Reads mesh files 
subroutine read_mesh_ini
USE o_MESH
USE o_PARAM
USE g_PARSUP
use g_CONFIG
use g_rotate_grid
!
IMPLICIT NONE
!
INTEGER           :: nq, nodes(3)
INTEGER           :: n1,n2,n3
INTEGER           :: n, nz, exit_flag
REAL(kind=WP)     :: x1, x2
INTEGER		  :: tag
! ===================
! Surface mesh
! ===================
  open (20,file=trim(meshpath)//'nod2d.out', status='old')
  open (21,file=trim(meshpath)//'elem2d.out', status='old')
  READ(20,*) nod2D
  ALLOCATE(coord_nod2D(2,nod2D))
    
  do n=1,nod2D
     read(20,*) nq, x1, x2, tag
     if (force_rotation) then
        call g2r(x1*rad, x2*rad, x1, x2)
        x1=x1/rad
        x2=x2/rad
     end if      
     coord_nod2D(1,nq)=x1*rad
     coord_nod2D(2,nq)=x2*rad
  end do
  CLOSE(20) 
      
  READ(21,*)  elem2D    
  ALLOCATE(elem2D_nodes(3,elem2D))
  
  do n=1,elem2D
     read(21,*) n1,n2,n3
     elem2D_nodes(1,n)=n1
     elem2D_nodes(2,n)=n2
     elem2D_nodes(3,n)=n3
  end do
  !
  CLOSE(21)
   	   	 
write(*,*) '========================='
write(*,*) 'Mesh is read'
write(*,*) '========================='
END SUBROUTINE read_mesh_ini
!=======================================================================
!> @brief 
!> Check the order of nodes in triangles; correct it if necessary to make
!! it same sense (clockwise) 
SUBROUTINE test_tri_ini
USE o_MESH
USE o_PARAM
USE g_CONFIG
IMPLICIT NONE
real(kind=WP)   ::  a(2), b(2), c(2),  r
integer         ::  n, nx, elnodes(3)

   
   DO n=1, elem2D
      elnodes=elem2D_nodes(:,n)
	  
          a=coord_nod2D(:,elnodes(1))
	  b=coord_nod2D(:,elnodes(2))-a
	  c=coord_nod2D(:,elnodes(3))-a
          
	  if(b(1)>cyclic_length/2.) b(1)=b(1)-cyclic_length
          if(b(1)<-cyclic_length/2.) b(1)=b(1)+cyclic_length
	  if(c(1)>cyclic_length/2.) c(1)=c(1)-cyclic_length
          if(c(1)<-cyclic_length/2.) c(1)=c(1)+cyclic_length
	  
	    
	  r=b(1)*c(2)-b(2)*c(1)
	  if (r>0) then
	  ! Vector b is to right of c
	  ! Exchange second and third nodes:
	  
	  nx=elnodes(2)
	  elnodes(2)=elnodes(3)
	  elnodes(3)=nx
	  elem2D_nodes(:,n)=elnodes
      end if
   END DO

END SUBROUTINE  test_tri_ini
!=========================================================================
!> @brief
!> Finds edges. Creates 3 files: edgenum.out, edges.out, edge_tri.out
SUBROUTINE find_edges_ini
USE o_MESH
USE o_PARAM
USE g_PARSUP
USE g_CONFIG
IMPLICIT NONE
integer, allocatable, dimension(:)    :: aux, aux1, aux2
integer, allocatable, dimension(:,:)  :: auxne, auxnn
integer                               :: counter, n, k,q,q1
integer                               :: elem,elem1, elems(2)
integer                               :: elnodes(3), ed(2),flag, eledges(3)
integer                               :: temp(100), node 
real(kind=WP)                         :: xc(2), xe(2), ax(3), amin
! ====================
! (a) find edges. To make the procedure fast 
! one needs neighbourhood arrays
! ====================

allocate(aux(nod2d), aux1(nod2D), aux2(nod2D))
aux=0
aux1=0
aux2=0
DO n=1,elem2D         
   elnodes=elem2D_nodes(:,n)
   aux(elnodes)=aux(elnodes)+1
END DO
   k=maxval(aux)     ! maximum number of neighbour elements
                     ! (there can be k+1 neighbor nodes)
allocate(auxne(k,nod2D),auxnn(k+1, nod2D))
aux=0
DO n=1,elem2D
   elnodes=elem2D_nodes(:,n)
   DO q=1,3         
   	  aux(elnodes(q))=aux(elnodes(q))+1
      auxne(aux(elnodes(q)),elnodes(q))=n
   END DO
END DO				 ! neighbor elements are found 
! count neighbour nodes				     
DO n=1, nod2D
     counter=0
	 DO k=1, aux(n)
	    elem=auxne(k,n)
		elnodes=elem2D_nodes(:,elem)
		DO q=1,3
		if(elnodes(q)==n) CYCLE 
		if(aux1(elnodes(q)).ne.1) then
		counter=counter+1
		aux1(elnodes(q))=1
		temp(counter)=elnodes(q)
		end if
		END DO
	 END DO
	 aux2(n)=counter
	 aux1(temp(1:counter))=0
	 auxnn(1:counter,n)=temp(1:counter)
END DO              ! neighbor nodes are found (excluding n) 
! ====================
! (b) Find edges and triangles containing them.
!     Write information to auxiliary file
! ====================
if (mype==0) then
 open(10, file='edges.out')
end if 
 counter=0

DO n=1,nod2D
      ! ==================== 
      ! form edges with n by cycling over neighbor
      ! nodes (if edges are not accounted yet). 
      ! New edges are added only if neighbor>n  
      ! ====================	  
   DO q=1,aux2(n)
      node=auxnn(q,n)
	  if(node<n) CYCLE
	  counter=counter+1   ! new edge (n,node)
	  ! find triangles containing n and node
	  flag=0
	  DO k=1, aux(n)
	  elem=auxne(k,n)
	  elnodes=elem2D_nodes(:,elem)
	    DO q1=1,3
		   if (elnodes(q1)==node) then
		   flag=flag+1
		   elems(flag)=elem
		   EXIT
		   end if
		END DO
	  END DO
	  if(mype==0)then
	  if(flag==2) write(10,'(4I10)') n, node, elems
          if(flag==1) write(10,'(4I10)') n, node, elems(1), -999
	  ! flag=1 if second triangle is absent (boundary edge)
	  if ((flag>2).or.(flag==0)) write(*,*) 'flag'
	  end if
   END DO
END DO
 if(mype==0) close(10)
 deallocate(auxnn, auxne)
 deallocate(aux2, aux1, aux)
 call MPI_Barrier(MPI_COMM_WORLD, MPIERR)   ! Wait untill PE=0 does its job
! ====================
! (c) Read list of edges and tri containing them from file 
! ====================
open(10, file='edges.out')
allocate(edges(2,counter))
allocate(edge_tri(2,counter))
edge2D=counter
 counter=0
DO n=1,edge2D
   read(10,*) ed,elem,elem1
   if(elem1>0) then
   counter=counter+1
   edges(:,counter)=ed
   edge_tri(1,counter)=elem
   edge_tri(2,counter)=elem1
   end if
END DO
rewind(10)
edge2D_in=counter
DO n=1,edge2D
   read(10,*) ed,elem,elem1
   if(elem1==-999) then
   counter=counter+1
   edges(:,counter)=ed
   edge_tri(1,counter)=elem
   edge_tri(2,counter)=elem1
   end if
END DO
! Edges from edge2D_in+1 to edge2D lie on the horizontal boundary
! The rest (1:edge2D_in) are internal edges
 close(10)
! ====================
! (d) the list of elements on both sides of edge e, edge_tri(:,e)
! should be ordered so that the first is a triangle which is to the left of 
! the edge (vector pointing from the first to the second node of the edge. 
! If the edge is on the boundary, there is only the first triangle.  
! ====================

DO n=1, edge2D
   ed=edges(:,n)
   if (edge_tri(1,n)<=0) then
         elem=edge_tri(1,n)
         elem1=edge_tri(2,n)
         edge_tri(1,n)=elem1
         edge_tri(2,n)=elem
   endif
   call elem_center(edge_tri(1,n), xc(1), xc(2))
   xc=xc-coord_nod2D(:,ed(1))
   xe=coord_nod2D(:,ed(2))-coord_nod2D(:,ed(1))
     if(xe(1)>=cyclic_length/2.) xe(1)=xe(1)-cyclic_length
     if(xe(1)<-cyclic_length/2.) xe(1)=xe(1)+cyclic_length
     if(xc(1)>=cyclic_length/2.) xc(1)=xc(1)-cyclic_length
     if(xc(1)<-cyclic_length/2.) xc(1)=xc(1)+cyclic_length

    if(xc(1)*xe(2)-xc(2)*xe(1)>0) then
	 ! Vector drawn to the center of the first triangle is to the right
	 ! of the edge vector. Triangles have to be exchanged:
	 elem=edge_tri(1,n)
	 elem1=edge_tri(2,n)
	 if(elem1>0) then  !TODO
	   edge_tri(1,n)=elem1
	   edge_tri(2,n)=elem
	 else
           if (elem<=0) write(*,*) '2 neighbouring elems are on the ground.'
	   elem=edges(2,n)         ! change the order of nodes
	   edges(2,n)=edges(1,n)
	   edges(1,n)=elem
   	 end if
    end if
END DO
 

	
   

! ====================
! (e) We need an array inverse to edge_tri listing edges
! of a given triangle 
! ====================
allocate(elem_edges(3,elem2D))
allocate(aux(elem2D))
aux=0
DO n=1, edge2D
   DO k=1,2
      q=edge_tri(k,n)   ! triangle number
	  if (q>0) then
	  aux(q)=aux(q)+1
	  elem_edges(aux(q),q)=n
	  end if
   END DO
END DO
deallocate(aux)
!> The edges in this list should be ordered so that they
!! are listed in the same rotation sense as nodes.
DO elem=1,elem2D
   elnodes=elem2D_nodes(:,elem)
   eledges=elem_edges(:,elem)
   DO q=1,3
      DO k=1,3
         if((edges(1,eledges(k)).ne.elnodes(q)).and. &
            (edges(2,eledges(k)).ne.elnodes(q))) then
           elem_edges(q,elem)=eledges(k)
	   exit
         end if
      END DO
   END DO
END DO
!> The edge and elem lists agree in the sense that edge 1 does not
!! contain node 1 and so on
 if(mype==0) then
 open(11, file=trim(meshpath)//'edgenum.out')
 write(11,*) edge2D
 write(11,*) edge2D_in
 close(11)
 open(10, file=trim(meshpath)//'edges.out')
 open(12, file=trim(meshpath)//'edge_tri.out')
 do n=1,edge2D
 write(10,*) edges(:,n)
 write(12,*) edge_tri(:,n)
 end do
 close(10)
 close(12)
 end if
END SUBROUTINE find_edges_ini
!===================================================================
!> @brief
!> Finds elemental and nodal levels.
!> Does some thresholding: if (depth>zbar(4)) x=zbar(4)
!> Fixes rough topography, by converting some oceans cells to ground cell(reflected by changing levels arrays)
!> Creates 2 files: elvls.out, nlvls.out
subroutine find_levels
use g_config
use o_mesh
use g_parsup
implicit none
INTEGER           :: nodes(3),elems(3)
integer :: elem, eledges(3), elem1, j, n, node, enum,count1,count2,exit_flag,i,nz,fileID=111
real*8 :: x,dmean
integer :: thers_lev=5
character*80 :: file_name
ALLOCATE(depth(nod2D))
        
file_name=trim(meshpath)//'aux3d.out'
open(fileID, file=file_name)
read(fileID,*) nl          ! the number of levels 
allocate(zbar(nl))         ! their standard depths
read(fileID,*) zbar
if(zbar(2)>0) zbar=-zbar   ! zbar is negative 
allocate(Z(nl-1))
Z=zbar(1:nl-1)+zbar(2:nl)  ! mid-depths of cells
Z=0.5_WP*Z
DO n=1,nod2D
         read(fileID,*) x
         if (x>0) x=-x
         if (x>zbar(thers_lev)) x=zbar(thers_lev) !TODO KK threshholding for depth
         depth(n)=x
END DO
close(fileID)
        
if(depth(2)>0) depth=-depth  ! depth is negative


allocate(nlevels(elem2D))
allocate(nlevels_nod2D(nod2D))


! ===================
! Compute the number of levels
! ===================
   DO n=1, elem2D
      nodes=elem2D_nodes(:,n)
      dmean=minval(depth(nodes))
      exit_flag=0
          DO nz=1,nl-1
                 if(Z(nz)<dmean) then
                        exit_flag=1
                    nlevels(n)=nz
                    exit
                 end if
          END DO
          if((exit_flag==0).and.(dmean<0)) nlevels(n)=nl
          if(dmean>=0) nlevels(n)=thers_lev
          if(nlevels(n)<thers_lev) nlevels(n)=thers_lev
    END DO


  do nz=4,nl
    exit_flag=0
    count1=0
    do while((exit_flag==0).and.(count1<100))
        exit_flag=1
        count1=count1+1
        do n=1,elem2D
         !find ocean cell
         if (nlevels(n)>=nz) then
          count2=0
          elems=elem_neighbors(:,n)
          do i=1,3
            if (elems(i)>1) then
            if (nlevels(elems(i))>=nz) then
                !count neighbours
                count2=count2+1
            endif
            endif
          enddo
          if (count2<2) then
                !if cell is "bad" convert to bottom cell
                nlevels(n)=nz-1
                !force recheck for all current ocean cells
                exit_flag=0
          endif
         endif
        enddo
    enddo
 enddo

  nlevels_nod2D=0
  DO n=1,elem2D
       DO j=1,3
       node=elem2D_nodes(j,n)
       if(nlevels_nod2D(node)<nlevels(n)) then
       nlevels_nod2D(node)=nlevels(n)
       end if
       END DO
  END DO

if (mype==0) then
   file_name=trim(meshpath)//'elvls.out'
   open(fileID, file=file_name)
   do n=1,elem2D
      write(fileID,*) nlevels(n)
   enddo
   close(fileID)

   file_name=trim(meshpath)//'nlvls.out'
   open(fileID, file=file_name)
   do n=1,nod2D
      write(fileID,*) nlevels_nod2D(n)
   enddo
   close(fileID)

   write(*,*) '========================='
   write(*,*) 'Mesh is read : ', 'nod2D=', nod2D,' elem2D=', elem2D, ' nl=', nl
   write(*,*) 'Min/max depth on mype: ',  -zbar(minval(nlevels)),-zbar(maxval(nlevels))
   write(*,*) '3D tracer nodes on mype ', sum(nlevels_nod2d)-(elem2D)
   write(*,*) '========================='
endif


end subroutine find_levels
!===================================================================

subroutine edge_center(n1, n2, x, y)
USE o_MESH
USE g_CONFIG
!
! Returns coordinates of edge center in x and y
! 
implicit none
integer       :: n1, n2   ! nodes of the edge
real(kind=WP) :: x, y, a(2), b(2)

a=coord_nod2D(:,n1)
b=coord_nod2D(:,n2)
if(a(1)-b(1)>cyclic_length/2.0) a(1)=a(1)-cyclic_length
if(a(1)-b(1)<-cyclic_length/2.0) b(1)=b(1)-cyclic_length
x=0.5_WP*(a(1)+b(1))
y=0.5_WP*(a(2)+b(2))
end subroutine edge_center
!====================================================================
subroutine elem_center(elem, x, y)
!
! Returns coordinates of elem center in x and y
!
USE o_MESH
USE g_CONFIG
implicit none
integer       :: elem, elnodes(3), k    
real(kind=WP) :: x, y, ax(3), amin

   elnodes=elem2D_nodes(:,elem)
   ax=coord_nod2D(1, elnodes)
   amin=minval(ax)
   DO k=1,3
   if(ax(k)-amin>cyclic_length/2.0) ax(k)=ax(k)-cyclic_length
   END DO
   x=sum(ax)/3.0_WP
   y=sum(coord_nod2D(2,elnodes))/3.0_WP
   
end subroutine elem_center
!=======================================================================
SUBROUTINE find_elem_neighbors_ini
! For each element three its element neighbors are found
USE o_MESH
USE g_PARSUP
implicit none
integer    :: elem, eledges(3), elem1, j, n, elnodes(3)
allocate(elem_neighbors(3,elem2D))
elem_neighbors=0
DO elem=1,elem2D
   
   eledges=elem_edges(:,elem)
   DO j=1,3
   elem1=edge_tri(1,eledges(j))
   if(elem1==elem) elem1=edge_tri(2,eledges(j))
   elem_neighbors(j,elem)=elem1
   END DO
   
END DO
 ! Among elem_neighbors there can be negative numbers. These correspond to 
 ! boundary elements for which neighbours are absent. However, an element 
 ! should have at least two valid neighbors

 ! Test that there are at least two neighbors at the surface:

DO elem=1,elem2D
   elem1=0
   DO j=1,3
   if(elem_neighbors(j,elem)>0) elem1=elem1+1
   END DO
   if (elem1<2) then
   write(*,*) 'find_elem_neighbors_ini:Insufficient number of neighbors ',elem
   write(*,*) 'find_elem_neighbors_ini:Elem neighbors ',elem_neighbors(:,elem)
   STOP
   end if
END DO    

 ! The rotation sense: corresponds to edges, and edges correspond 
 ! to nodes

 ! =============
 ! To facilitate computations the neibourhood
 ! information is assembled
 ! =============	 
 allocate(nod_in_elem2D_num(nod2D))
 nod_in_elem2D_num=0
 do n=1,elem2D
    elnodes=elem2D_nodes(:,n)
    nod_in_elem2D_num(elnodes)=nod_in_elem2D_num(elnodes)+1
 end do
 allocate(nod_in_elem2D(maxval(nod_in_elem2D_num),nod2D))
 nod_in_elem2D=0
 
 nod_in_elem2D_num=0
 do n=1,elem2D   
    elnodes=elem2D_nodes(:,n)
    nod_in_elem2D_num(elnodes)=nod_in_elem2D_num(elnodes)+1
    do j=1, 3
       nod_in_elem2D(nod_in_elem2D_num(elnodes(j)),elnodes(j))=n
    end do
 end do  
END SUBROUTINE find_elem_neighbors_ini
!===================================================================
! Stiffness matrix for the elevation
subroutine stiff_mat_ini
  use o_MESH
  
  !
  implicit none
  integer                :: n, n1, n2, q, ed
  integer                :: n_num(nod2D)
  !
  ! a) Allocate, initialize fields
! Only the master does this work (and allocates the memory)
  ssh_stiff%dim = nod2D   
  ssh_stiff%nza = 2*edge2D ! Two matrix entries for each edge
  allocate(ssh_stiff%rowptr(nod2D+1))
  allocate(ssh_stiff%colind(2*edge2D))

  ssh_stiff%colind(:) = -1
  n_num(:)=0

  ! b) Neighbourhood information
  !     For Metis, the diagonal (node connected to itself) is not needed 
  
  Do ed=1, edge2D
     ! n_num contains the number of neighbors
     n_num(edges(1:2,ed)) = n_num(edges(1:2,ed)) + 1
  END DO 

  ssh_stiff%rowptr(1)=1	 
  do n=1,nod2D
     ssh_stiff%rowptr(n+1) = ssh_stiff%rowptr(n)+n_num(n)
  end do

  ! c)
  ! Loop through the edges again, collect the connection information
  ! in the compact storage ssh_stiff.
  ! Life is easy, because every pair (n1,n2) is given exactly once by
  ! the corresponding edge.
  do ed=1,edge2D
     n1 = edges(1,ed)
     n2 = edges(2,ed)

     ! insert (n1,n2)
     do q = ssh_stiff%rowptr(n1), ssh_stiff%rowptr(n1+1)-1
        if (ssh_stiff%colind(q) < 0) then
           ssh_stiff%colind(q) = n2
           exit
        end if
     enddo
     ! insert (n2,n1)
     do q = ssh_stiff%rowptr(n2), ssh_stiff%rowptr(n2+1)-1
        if (ssh_stiff%colind(q) < 0) then
           ssh_stiff%colind(q) = n1
           exit
        end if
     enddo
  end do

end subroutine stiff_mat_ini

