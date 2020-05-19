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
  use MOD_MESH
  use o_MESH
  use g_PARSUP
  use g_CONFIG
  use g_rotate_grid
  
  implicit none

interface
   subroutine read_mesh_ini(mesh)
     use mod_mesh
     type(t_mesh), intent(inout)  , target :: mesh
   end subroutine read_mesh_ini
end interface

interface
   subroutine test_tri_ini(mesh)
     use mod_mesh
     type(t_mesh), intent(inout)  , target :: mesh
   end subroutine test_tri_ini
end interface
interface
   subroutine find_edges_ini(mesh)
     use mod_mesh
     type(t_mesh), intent(inout)  , target :: mesh
   end subroutine find_edges_ini
end interface
interface
   subroutine find_elem_neighbors_ini(mesh)
     use mod_mesh
     type(t_mesh), intent(inout)  , target :: mesh
   end subroutine find_elem_neighbors_ini
end interface
interface
   subroutine find_levels(mesh)
     use mod_mesh
     type(t_mesh), intent(inout)  , target :: mesh
   end subroutine find_levels
end interface
interface
   subroutine stiff_mat_ini(mesh)
     use mod_mesh
     type(t_mesh), intent(inout)  , target :: mesh
   end subroutine stiff_mat_ini
end interface
interface
   subroutine set_par_support_ini(mesh)
     use mod_mesh
     type(t_mesh), intent(inout)  , target :: mesh
   end subroutine set_par_support_ini
end interface
interface
   subroutine communication_ini(mesh)
     use mod_mesh
     type(t_mesh), intent(inout)  , target :: mesh
   end subroutine communication_ini
end interface

  character(len=1000)         :: nmlfile  !> name of configuration namelist file
  integer                     :: start_t, interm_t, finish_t, rate_t
  type(t_mesh), target, save  :: mesh

  call system_clock(start_t, rate_t)
  interm_t = start_t
  
  nmlfile ='namelist.config'
  open (20,file=nmlfile)
  read (20,NML=paths)         ! We need MeshPath
  read (20,NML=geometry)      ! We need cyclic_length and cartesian
  read (20,NML=machine)       ! We need partitioning hierarchy
  close (20)
  cyclic_length=cyclic_length*rad 
  call set_mesh_transform_matrix  !(rotated grid) 
  call read_mesh_ini(mesh)
  call system_clock(finish_t)
  print '("**** Reading initial mesh time = ",f12.3," seconds. ****")', &
       real(finish_t-interm_t)/real(rate_t)
  interm_t = finish_t
  call test_tri_ini(mesh)
  call find_edges_ini(mesh)
  call find_elem_neighbors_ini(mesh)
  call find_levels(mesh)

! NR Some arrays are not needed for partitioning, after setting up the grid
  deallocate(mesh%coord_nod2D)
  deallocate(mesh%edge_tri)
  deallocate(mesh%zbar,mesh%Z,mesh%nlevels,mesh%depth)
  call system_clock(finish_t)
  print '("**** Checking and initializing mesh time = ",f12.3," seconds. ****")', &
       real(finish_t-interm_t)/real(rate_t)
  interm_t = finish_t

  call stiff_mat_ini(mesh)
  call set_par_support_ini(mesh)
  call system_clock(finish_t)
  print '("**** Partitioning time = ",f12.3," seconds. ****")', &
       real(finish_t-interm_t)/real(rate_t)
  interm_t = finish_t
  call communication_ini(mesh)
  call system_clock(finish_t)
  print '("**** Storing partitioned mesh time = ",f12.3," seconds. ****")', &
       real(finish_t-interm_t)/real(rate_t)
  print '("**** Total time = ",f12.3," seconds. ****")', &
       real(finish_t-start_t)/real(rate_t)
end program MAIN
!=============================================================================
!> @brief
!> Reads mesh files 
subroutine read_mesh_ini(mesh)
USE MOD_MESH
USE o_PARAM
USE g_PARSUP
use g_CONFIG
use g_rotate_grid
!
IMPLICIT NONE
!
type(t_mesh), intent(inout), target :: mesh
INTEGER                             :: nq
INTEGER                             :: n1,n2,n3
INTEGER                             :: n, nz, exit_flag
REAL(kind=WP)                       :: x1, x2
INTEGER	                            :: tag
INTEGER, allocatable                :: elem_data(:)
INTEGER                             :: i_error
#include "associate_mesh_ini.h"
! ===================
! Surface mesh
! ===================
  open (20,file=trim(meshpath)//'nod2d.out', status='old')
  open (21,file=trim(meshpath)//'elem2d.out', status='old')
  READ(20,*) mesh%nod2D
  ALLOCATE(mesh%coord_nod2D(2,mesh%nod2D)) 
  coord_nod2D => mesh%coord_nod2D !required after the allocation, otherwise the pointer remains undefined
    
  do n=1, mesh%nod2D
     read(20,*) nq, x1, x2, tag
     if (force_rotation) then
        call g2r(x1*rad, x2*rad, x1, x2)
        x1=x1/rad
        x2=x2/rad
     end if      
     mesh%coord_nod2D(1,nq)=x1*rad
     mesh%coord_nod2D(2,nq)=x2*rad
  end do
  CLOSE(20) 
      
  READ(21,*)  mesh%elem2D    
  ALLOCATE(mesh%elem2D_nodes(4,mesh%elem2D))
  elem2D_nodes => mesh%elem2D_nodes !required after the allocation, otherwise the pointer remains undefined
  ALLOCATE(elem_data(4*mesh%elem2D))
  elem_data(:)=-1
  
  ! meshes with quads have 4 columns, but TsunAWI grids may be
  ! purely triangular, with 3 columns each. Test, how many
  ! columns there are!  
  read(21,*,iostat=i_error) elem_data(1:4*mesh%elem2D)
  if (i_error == 0) then      ! There is a fourth column => quad or mixed mesh (not working yet!)
     mesh%elem2D_nodes = reshape(elem_data, shape(mesh%elem2D_nodes))
  else     ! No fourth column => triangles only
     mesh%elem2D_nodes(1:3,:) = reshape(elem_data, shape(mesh%elem2D_nodes(1:3,:)))
     mesh%elem2D_nodes(4,:)   = mesh%elem2D_nodes(1,:)
  end if
     
  deallocate(elem_data)
!!$  do n=1,elem2D
!!$     read(21,*) n1,n2,n3
!!$     elem2D_nodes(1,n)=n1
!!$     elem2D_nodes(2,n)=n2
!!$     elem2D_nodes(3,n)=n3
!!$  end do
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
SUBROUTINE test_tri_ini(mesh)
USE MOD_MESH
USE o_PARAM
USE g_CONFIG
IMPLICIT NONE

real(kind=WP)                       ::  a(2), b(2), c(2),  r
integer                             ::  n, nx, elnodes(3)
type(t_mesh), intent(inout), target :: mesh
#include "associate_mesh_ini.h"
   
   DO n=1, elem2D
      elnodes=elem2D_nodes(1:3,n)
	  
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
	  elem2D_nodes(1:3,n)=elnodes
      end if
   END DO
END SUBROUTINE  test_tri_ini
!=========================================================================
!> @brief
!> Finds edges. Creates 3 files: edgenum.out, edges.out, edge_tri.out
SUBROUTINE find_edges_ini(mesh)
USE MOD_MESH
USE o_MESH
USE o_PARAM
USE g_PARSUP
USE g_CONFIG
IMPLICIT NONE

interface
   subroutine elem_center(elem, x, y, mesh)
     USE MOD_MESH
     USE g_CONFIG
     integer, intent(in)        :: elem
     real(kind=WP), intent(out) :: x, y
     type(t_mesh), intent(in), target   :: mesh
   end subroutine elem_center
end interface

integer, allocatable                  :: aux1(:), ne_num(:), ne_pos(:,:)
integer                               :: counter, counter_in, n, k, q
integer                               :: elem, elem1, elems(2), q1, q2
integer                               :: elnodes(4), ed(2), flag, eledges(4)
integer                               :: temp(100), node 
real(kind=WP)                         :: xc(2), xe(2), ax(3), amin
type(t_mesh), intent(inout), target :: mesh
#include "associate_mesh_ini.h"
! ====================
! (a) find edges. To make the procedure fast 
! one needs neighbourhood arrays
! ====================

allocate(ne_num(nod2d), ne_pos(MAX_ADJACENT, nod2D), nn_num(nod2D))
ne_num=0
DO n=1,elem2D
    elnodes=elem2D_nodes(:,n)
    q1=4
    if(elnodes(1)==elnodes(4)) q1=3
    DO q=1,q1
       ne_num(elnodes(q))=ne_num(elnodes(q))+1
       if (ne_num(elnodes(q)) > MAX_ADJACENT ) then
          print *,'Parameter in o_MESH from ocean_modules.F90, too small.'
          print *,'Recompile with larger value for MAX_ADJACENT.'
          stop
       else
          ne_pos(ne_num(elnodes(q)),elnodes(q))=n
       endif
    END Do
END DO                 ! neighbor elements are found 

! count neighbour nodes				     
! In quads we should count the nodes that are 
! connected by edges!
allocate(aux1(nod2D))                   
aux1=0

DO n=1, nod2D
    counter=0
    DO k=1, ne_num(n)
        elem=ne_pos(k,n)
        elnodes=elem2D_nodes(:,elem)        
        if(elnodes(1)==elnodes(4)) then
            DO q=1,3
                if(elnodes(q)==n) CYCLE 
                if(aux1(elnodes(q)).ne.1) then
                    counter=counter+1
                    aux1(elnodes(q))=1
                    temp(counter)=elnodes(q)
                end if
            END DO
        else
        ! Find the position of n in elnodes:
            if(elnodes(1)==n .or. elnodes(3)==n) then
                ed(1)=elnodes(2)
                ed(2)=elnodes(4)
            else
                ed(1)=elnodes(1)
                ed(2)=elnodes(3)
            end if
            DO q=1,2
                if(aux1(ed(q)).ne.1) then
                    counter=counter+1
                    aux1(ed(q))=1
                    temp(counter)=ed(q)
                end if
            END DO
        end if
    END DO
    nn_num(n)=counter         
    aux1(temp(1:counter))=0
END DO              
 
allocate(nn_pos(maxval(nn_num)+1,nod2D))
nn_pos = -1
aux1=0
 
DO n=1, nod2D
    counter=0
    DO k=1, ne_num(n)
        elem=ne_pos(k,n)
        elnodes=elem2D_nodes(:,elem)        
        if(elnodes(1)==elnodes(4)) then 
            DO q=1,3
                if(elnodes(q)==n) CYCLE 
                if(aux1(elnodes(q)).ne.1) then
                    counter=counter+1
                    aux1(elnodes(q))=1
                    temp(counter)=elnodes(q)
                end if
            END DO
        else
            ! Find the position of n in elnodes:
            if(elnodes(1)==n .or. elnodes(3)==n) then
                ed(1)=elnodes(2)
                ed(2)=elnodes(4)
            else
                ed(1)=elnodes(1)
                ed(2)=elnodes(3)
            end if
            DO q=1,2
                if(aux1(ed(q)).ne.1) then
                    counter=counter+1
                    aux1(ed(q))=1
                    temp(counter)=ed(q)
                end if
            END DO
        end if
    END DO
    nn_num(n)=counter+1
    aux1(temp(1:counter))=0
    nn_pos(1,n)=n 
    nn_pos(2:counter+1,n)=temp(1:counter)
END DO              
deallocate(aux1)
! neighboring nodes are found. First in the list is the node itself 

! ====================
! (b) Find edges and elements containing them.
!     Write information to auxiliary file
! ====================
! open(10, file='edges.out')

 ! Count edges: 
 ! ==================== 
 ! form edges with n by cycling over neighbor
 ! nodes (if edges are not accounted yet). 
 ! New edges are added only if neighbor>n  
 ! ====================
 counter = 0
 DO n=1,nod2D
    counter = counter + count(nn_pos(2:nn_num(n),n)>nn_pos(1,n))
 end do
 edge2D=counter

 allocate(mesh%edges   (2, edge2D))
 allocate(mesh%edge_tri(2, edge2D))
 edges    => mesh%edges    !required after the allocation, otherwise the pointer remains undefined
 edge_tri => mesh%edge_tri !required after the allocation, otherwise the pointer remains undefined
 counter_in=0 

 DO n=1,nod2D
    DO q=2,nn_num(n)
       node=nn_pos(q,n)
       if(node<n) CYCLE
       counter=counter+1   ! new edge (n,node)
       ! find elements containing n and node
       flag=0
       DO k=1, ne_num(n)
          elem=ne_pos(k,n)
          elnodes=elem2D_nodes(:,elem)
          q2=4
          if(elnodes(1)==elnodes(4)) q2=3
          DO q1=1,q2
             if (elnodes(q1)==node) then
                flag=flag+1
                elems(flag)=elem
                EXIT
             end if
          END DO
       END DO
       if(flag==2) then
          counter_in=counter_in+1
          edges(1,counter_in)=n
          edges(2,counter_in)=node
          edge_tri(:,counter_in)=elems
!          write(10,'(4I10)') n, node, elems
!       else if (flag==1) then
!          write(10,'(4I10)') n, node, elems(1), -999
!       else
!          write(*,*) 'flag'
       end if
    END DO
 END DO
 edge2D_in=counter_in

 ! Repeat to collect boundary edges:   
 counter=0
 DO n=1,nod2D
    DO q=2,nn_num(n)
       node=nn_pos(q,n)
       if(node<n) CYCLE
       counter=counter+1   ! new edge (n,node)
       ! find triangles containing n and node
       flag=0
       DO k=1, ne_num(n)
          elem=ne_pos(k,n)
          elnodes=elem2D_nodes(:,elem)
          q2=4
          if(elnodes(1)==elnodes(4)) q2=3
          DO q1=1,q2
             if (elnodes(q1)==node) then
                flag=flag+1
                elems(flag)=elem
                EXIT
             end if
          END DO
       END DO
       if(flag==1) then
          counter_in=counter_in+1
          edges(1,counter_in)=n
          edges(2,counter_in)=node
          elems(2)=-999
          edge_tri(:,counter_in)=elems
       end if
    END DO
 END DO
 !close(10)
 ! Edges from edge2D_in+1 to edge2D lie on the horizontal boundary
 ! The rest (1:edge2D_in) are internal edges

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
    call elem_center(edge_tri(1,n), xc(1), xc(2), mesh)
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
 allocate(mesh%elem_edges(4,elem2D))
 elem_edges => mesh%elem_edges !required after the allocation, otherwise the pointer remains undefined
 allocate(aux1(elem2D))
 aux1=0
 DO n=1, edge2D
    DO k=1,2
       q=edge_tri(k,n)   ! triangle number
       if (q>0) then
	  aux1(q)=aux1(q)+1
	  elem_edges(aux1(q),q)=n
       end if
    END DO
 END DO
 deallocate(aux1)
 ! We order the edges in this list so that they
 ! are listed in the same rotation sense as nodes.
 ! First is the edge formed by elnodes(1:2), and so on
 DO elem=1,elem2D
    elnodes=elem2D_nodes(:,elem)
    q1=4
    if(elnodes(1)==elnodes(4)) q1=3
    eledges=elem_edges(:,elem)
    DO q=1,q1-1
       DO k=1,q1
          if(((edges(1,eledges(k))==elnodes(q)).and. &
               (edges(2,eledges(k))==elnodes(q+1))).or. &
               ((edges(1,eledges(k))==elnodes(q+1)).and. &
               (edges(2,eledges(k))==elnodes(q)))) then
             elem_edges(q,elem)=eledges(k)
             exit
          end if
       END DO
    END DO
    DO k=1,q1
       if(((edges(1,eledges(k))==elnodes(q1)).and. &
            (edges(2,eledges(k))==elnodes(1))).or. &
            ((edges(1,eledges(k))==elnodes(1)).and. &
            (edges(2,eledges(k))==elnodes(q1)))) then
          elem_edges(q1,elem)=eledges(k)
          exit
       end if
    END DO
    if(q1==3) elem_edges(4,elem)=elem_edges(1,elem)
 END DO

 !> The edge and elem lists agree in the sense that edge 1 does not
 !! contain node 1 and so on
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
 deallocate(ne_num, ne_pos)
END SUBROUTINE find_edges_ini
!=========================================================================
!> @brief
!> Finds elemental and nodal levels.
!> Does some thresholding: if (depth>zbar(4)) x=zbar(4)
!> Fixes rough topography, by converting some oceans cells to ground cell(reflected by changing levels arrays)
!> Creates 2 files: elvls.out, nlvls.out
subroutine find_levels(mesh)
use g_config
use mod_mesh
use g_parsup
implicit none
INTEGER :: nodes(3), elems(3), eledges(3)
integer :: elem, elem1, j, n, q, node, enum,count1,count2,exit_flag,i,nz,fileID=111
real(kind=WP) :: x,dmean
integer :: thers_lev=5
character*200 :: file_name
type(t_mesh), intent(inout), target :: mesh
#include "associate_mesh_ini.h"


ALLOCATE(mesh%depth(nod2D))
depth => mesh%depth !required after the allocation, otherwise the pointer remains undefined
file_name=trim(meshpath)//'aux3d.out'
open(fileID, file=file_name)
read(fileID,*) nl          ! the number of levels 
allocate(mesh%zbar(nl))         ! their standard depths
zbar => mesh%zbar !required after the allocation, otherwise the pointer remains undefined
read(fileID,*) zbar
if(zbar(2)>0) zbar=-zbar   ! zbar is negative 
allocate(mesh%Z(nl-1))
Z => mesh%Z !required after the allocation, otherwise the pointer remains undefined
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


allocate(mesh%nlevels(elem2D))
nlevels => mesh%nlevels             !required after the allocation, otherwise the pointer remains undefined
allocate(mesh%nlevels_nod2D(nod2D))
nlevels_nod2D => mesh%nlevels_nod2D !required after the allocation, otherwise the pointer remains undefined

! ===================
! Compute the number of levels
! ===================
   DO n=1, elem2D
      nodes=elem2D_nodes(1:3,n)
      
      ! depth of element is  shallowest depth of sorounding vertices
      !dmean=maxval(depth(nodes))
      
      ! depth of element is deepest depth of sorounding vertices
      !dmean=minval(depth(nodes))
      
      ! depth of element is  mean depth of sorounding vertices
      dmean=sum(depth(nodes))/3.0
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
       do while((exit_flag==0).and.(count1<1000))
          exit_flag=1
          count1=count1+1
          do n=1,elem2D
             q = merge(3,4,elem2D_nodes(1,n) == elem2D_nodes(4,n))
             !find ocean cell
             if (nlevels(n)>=nz) then
                count2=0
                elems=elem_neighbors(1:3,n)
                do i=1,q
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
     q = merge(3,4,elem2D_nodes(1,n) == elem2D_nodes(4,n))
     DO j=1,q
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

subroutine edge_center(n1, n2, x, y, mesh)
USE MOD_MESH
USE g_CONFIG
!
! Returns coordinates of edge center in x and y
! 
implicit none
integer       :: n1, n2   ! nodes of the edge
real(kind=WP) :: x, y, a(2), b(2)
type(t_mesh), intent(inout), target :: mesh
#include "associate_mesh_ini.h"

a=coord_nod2D(:,n1)
b=coord_nod2D(:,n2)
if(a(1)-b(1)>cyclic_length/2.0) a(1)=a(1)-cyclic_length
if(a(1)-b(1)<-cyclic_length/2.0) b(1)=b(1)-cyclic_length
x=0.5_WP*(a(1)+b(1))
y=0.5_WP*(a(2)+b(2))
end subroutine edge_center
!====================================================================
subroutine elem_center(elem, x, y, mesh)
!
! Returns coordinates of elem center in x and y
!
USE MOD_MESH
USE g_CONFIG
implicit none
integer, intent(in)  :: elem
integer              ::  elnodes(3), k    
real(kind=WP), intent(out) :: x, y
real(kind=WP)        ::  ax(3), amin
type(t_mesh), intent(in), target :: mesh
#include "associate_mesh_ini.h"

   elnodes=elem2D_nodes(1:3,elem)
   ax=coord_nod2D(1, elnodes)
   amin=minval(ax)
   DO k=1,3
   if(ax(k)-amin>cyclic_length/2.0) ax(k)=ax(k)-cyclic_length
   END DO
   x=sum(ax)/3.0_WP
   y=sum(coord_nod2D(2,elnodes))/3.0_WP
   
end subroutine elem_center
!=======================================================================
SUBROUTINE find_elem_neighbors_ini(mesh)
! For each element three its element neighbors are found
USE MOD_MESH
USE g_PARSUP
implicit none
integer    :: elem, eledges(3), elem1, j, n, elnodes(3)
type(t_mesh), intent(inout), target :: mesh
#include "associate_mesh_ini.h"

allocate(mesh%elem_neighbors(4,elem2D))
elem_neighbors => mesh%elem_neighbors !required after the allocation, otherwise the pointer remains undefined
elem_neighbors=0
DO elem=1,elem2D
   
   eledges=elem_edges(1:3,elem)
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
    if (mype==0) then 
    write(*,*) '____________________________________________________________________'
    write(*,*) ' ERROR: The mesh you want to partitioning contains triangles that'
    write(*,*) '        have just one neighbor, this was OK for FESOM1.4 but not'
    write(*,*) '        for FESOM2.0.                                           '
    write(*,*) '                                                                '
    write(*,*) '          #########################################             '
    write(*,*) '          ################### o ###################             '
    write(*,*) '          ################# ./|\. #################             '
    write(*,*) '              Land    ### ./|||||\. ###    Land                 '
    write(*,*) '          ############## /|||||||||\ ##############             '
    write(*,*) '          --o-----------o-----------o-----------o--             '
    write(*,*) '          ./ \.       ./ \.       ./ \.       ./ \.             '
    write(*,*) '               \.   ./     \.   ./     \.   ./                  '
    write(*,*) '                 \ /         \ /         \ /                    '
    write(*,*) '            ------o-----------o-----------o-------              '
    write(*,*) '                ./ \.       ./ \.       ./ \.                   '
    write(*,*) '                                                                '
    write(*,*) '        Take a programm of your choice (Python, Matlab ...) and '
    write(*,*) '        eliminate these triangles and the corresponding         '
    write(*,*) '        unconnected vertice and try to re-partitioning again    '
    write(*,*) '____________________________________________________________________'
    end if 
   STOP
   end if
END DO    

 ! The rotation sense: corresponds to edges, and edges correspond 
 ! to nodes

 ! =============
 ! To facilitate computations the neibourhood
 ! information is assembled
 ! =============	 
 allocate(mesh%nod_in_elem2D_num(nod2D))
 nod_in_elem2D_num => mesh%nod_in_elem2D_num !required after the allocation, otherwise the pointer remains undefined
 nod_in_elem2D_num=0
 do n=1,elem2D
    elnodes=elem2D_nodes(1:3,n)
    nod_in_elem2D_num(elnodes)=nod_in_elem2D_num(elnodes)+1
 end do
 allocate(mesh%nod_in_elem2D(maxval(nod_in_elem2D_num),nod2D))
 nod_in_elem2D => mesh%nod_in_elem2D
 nod_in_elem2D=0
 
 nod_in_elem2D_num=0
 do n=1,elem2D   
    elnodes=elem2D_nodes(1:3,n)
    nod_in_elem2D_num(elnodes)=nod_in_elem2D_num(elnodes)+1
    do j=1, 3
       nod_in_elem2D(nod_in_elem2D_num(elnodes(j)),elnodes(j))=n
    end do
 end do  
END SUBROUTINE find_elem_neighbors_ini
!===================================================================
! Stiffness matrix for the elevation
subroutine stiff_mat_ini(mesh)
  use MOD_MESH
  
  !
  implicit none
  integer                :: i, j, n, q, el, elem_nodes_max, nod(4)
  integer, allocatable   :: num_ne(:), ne(:,:)
  !
  type(t_mesh), intent(inout), target :: mesh
#include "associate_mesh_ini.h"

  ssh_stiff%dim = nod2D   
  allocate(mesh%ssh_stiff%rowptr(nod2D+1))
  ssh_stiff => mesh%ssh_stiff !required after the allocation, otherwise the pointer remains undefined

  allocate(num_ne(nod2D), ne(MAX_ADJACENT,nod2D))
  num_ne(:)           = 0
  ne(:,:)             = -1
 
  ! Check the maximum number of nodes in an element (FESOM triangular meshes = 3, Hybrid meshes = 4)
  elem_nodes_max = size(elem2D_nodes, 1)

  ! Fill node adjacency info
  ! all nodes in an element are adjacent in the sense of being halo nodes
  ! (also the opposite nodes of a quad: no edge, but the indirect connection
  ! should be taken into account by metis domain decomposition)
  do el=1,elem2D
     nod(1:elem_nodes_max) = elem2D_nodes(1:elem_nodes_max,el)  ! Fortran-numbering
     q = elem_nodes_max
     if (nod(1) == nod(elem_nodes_max)) q = q-1  ! triangle
    
     do i=2,q
        do j=1,i-1
           if (all(ne(:,nod(i)) /= nod(j))) then
              num_ne(nod(i)) = num_ne(nod(i)) + 1
              num_ne(nod(j)) = num_ne(nod(j)) + 1

              if (max(num_ne(nod(i)), num_ne(nod(j))) > MAX_ADJACENT ) then
                 print *,'Parameter in o_MESH from ocean_modules.F90, too small.'
                 print *,'Recompile with larger value for MAX_ADJACENT.'
                 stop
              else
                 ne(num_ne(nod(i)), nod(i)) = nod(j)
                 ne(num_ne(nod(j)), nod(j)) = nod(i)
              endif
           endif
        end do
     end do
  end do

! copy adjacency matrix to CSR-format
  ssh_stiff%rowptr(1) = 1
  do n=1,nod2D
     ssh_stiff%rowptr(n+1) = ssh_stiff%rowptr(n)+num_ne(n)
  end do

  allocate(ssh_stiff%colind(ssh_stiff%rowptr(nod2D+1)-1))  
  ssh_stiff => mesh%ssh_stiff !required after the allocation, otherwise the pointer remains undefined

  !required after the allocation, otherwise the pointer remains undefined
  do n=1,nod2D
     ssh_stiff%colind(ssh_stiff%rowptr(n):ssh_stiff%rowptr(n+1)-1) = ne(1:num_ne(n),n)
  end do

  deallocate(num_ne, ne)

end subroutine stiff_mat_ini

!===================================================================
! Setup of communication arrays
subroutine communication_ini(mesh)
  use MOD_MESH
  USE g_CONFIG
  USE g_PARSUP
  use omp_lib
  implicit none

  integer        :: n
  character*10   :: npes_string
  character*200  :: dist_mesh_dir
  LOGICAL        :: L_EXISTS
  type(t_mesh), intent(inout), target :: mesh
#include "associate_mesh_ini.h"

  ! Create the distributed mesh subdirectory
  write(npes_string,"(I10)") npes
  dist_mesh_dir=trim(meshpath)//'dist_'//trim(ADJUSTL(npes_string))//'/'
  INQUIRE(file=trim(dist_mesh_dir), EXIST=L_EXISTS)
  if (.not. L_EXISTS) call system('mkdir '//trim(dist_mesh_dir))

#ifdef OMP_MAX_THREADS
!$OMP PARALLEL NUM_THREADS(OMP_MAX_THREADS)
  if (OMP_GET_THREAD_NUM() == 0) then
     write(*,*) 'Setting up communication arrays using ', OMP_GET_NUM_THREADS(), ' threads'
  endif
#else
!$OMP PARALLEL NUM_THREADS(1)
  write(*,*) 'Setting up communication arrays using 1 thread (serially)'
#endif
  
!$OMP DO
  do n = 0, npes-1
     mype = n ! mype is threadprivate and must not be iterator
     call communication_nodn(mesh)
     call communication_elemn(mesh)
     call save_dist_mesh(mesh)         ! Write out communication file com_infoxxxxx.out
  end do
!$OMP END DO
!$OMP END PARALLEL

  deallocate(mesh%elem_neighbors)
  deallocate(mesh%elem_edges)
  deallocate(part)
  write(*,*) 'Communication arrays have been set up'   
end subroutine communication_ini
!=================================================================
subroutine set_par_support_ini(mesh)
  use g_PARSUP
  use iso_c_binding, only: idx_t=>C_INT32_T
  use MOD_MESH
  use g_config
  implicit none

interface 
   subroutine check_partitioning(mesh)
     use MOD_MESH
     type(t_mesh), intent(inout)  , target :: mesh
   end subroutine check_partitioning
end interface

  integer         :: n, j, k, nini, nend, ierr
  integer(idx_t)  :: np(10)
  type(t_mesh), intent(inout), target :: mesh

  interface 
     subroutine partit(n,ptr,adj,wgt,np,part) bind(C)
       use iso_c_binding, only: idx_t=>C_INT32_T
       integer(idx_t), intent(in)  :: n, ptr(*), adj(*), wgt(*), np(*)
       integer(idx_t), intent(out) :: part(*)
     end subroutine partit
  end interface
#include "associate_mesh_ini.h"

  ! Construct partitioning vector
  if (n_levels<1 .OR. n_levels>10) then
     print *,'Number of hierarchic partition levels is out of range [1-10]! Aborting...'
     call MPI_ABORT( MPI_COMM_FESOM, 1 )
  end if

  np(:) = n_part(:)               ! Number of partitions on each hierarchy level
  if (n_part(1) == 0) then        ! Backward compatibility case: Take the number of 
     np(1) = npes                 ! partitions from the number of MPI processes
     n_levels = 1
  end if
  if (n_levels < 10) then         ! 0 is an indicator of the last hierarchy level
     np(n_levels+1) = 0 
  end if

  allocate(part(nod2D))
  part=0

  npes = PRODUCT(np(1:n_levels))
  if(npes<2) then
     print *,'Total number of parallel partitions is less than one! Aborting...'
     stop
  end if
  
  write(*,*) 'Calling partit for npes=', np
  call partit(ssh_stiff%dim, ssh_stiff%rowptr, ssh_stiff%colind, &
       nlevels_nod2D, np, part)

  call check_partitioning(mesh)

  write(*,*) 'Partitioning is done.'

! The stiffness matrix is no longer needed. 
  deallocate(mesh%ssh_stiff%rowptr)
  deallocate(mesh%ssh_stiff%colind)
        
  !NR No longer needed - last use was as weight for partitioning
  deallocate(mesh%nlevels_nod2D)
end subroutine set_par_support_ini
!=======================================================================
subroutine check_partitioning(mesh)

  ! In general, METIS 5 has several advantages compared to METIS 4, e.g.,
  !   * neighbouring tasks get neighbouring partitions (important for multicore computers!)
  !   * lower maximum of weights per partition (better load balancing)
  !   * lower memory demand
  !
  ! BUT: there might be outliers, single nodes connected to their partition by
  !      only one edge or even completely isolated. This spoils everything :-(
  !
  ! This routine checks for isolated nodes and moves them to an adjacent partition,
  ! trying not to spoil the load balance.

  use MOD_MESH
  use g_PARSUP
  integer :: i, j, k, n, n_iso, n_iter, is, ie, kmax, np
  integer :: nod_per_partition(2,0:npes-1)
  integer :: max_nod_per_part(2), min_nod_per_part(2)
  integer :: average_nod_per_part(2), node_neighb_part(100)
  logical :: already_counted, found_part

  integer :: max_adjacent_nodes
  integer, allocatable :: ne_part(:), ne_part_num(:), ne_part_load(:,:)
  type(t_mesh), intent(inout), target :: mesh
#include "associate_mesh_ini.h"

  ! Check load balancing
  do i=0,npes-1
     nod_per_partition(1,i) = count(part(:) == i)
     nod_per_partition(2,i) = sum(nlevels_nod2D,part(:) == i)
  enddo

  min_nod_per_part(1) = minval( nod_per_partition(1,:))
  min_nod_per_part(2) = minval( nod_per_partition(2,:))

  max_nod_per_part(1) = maxval( nod_per_partition(1,:))
  max_nod_per_part(2) = maxval( nod_per_partition(2,:))

  average_nod_per_part(1) = nod2D / npes
  average_nod_per_part(2) = sum(nlevels_nod2D(:)) / npes

  ! Now check for isolated nodes (connect by one or even no edge to other
  ! nodes of its partition) and repair, if possible

  max_adjacent_nodes = maxval(ssh_stiff%rowptr(2:nod2D+1) - ssh_stiff%rowptr(1:nod2D))
  allocate(ne_part(max_adjacent_nodes), ne_part_num(max_adjacent_nodes), &
       ne_part_load(2,max_adjacent_nodes))

  print *,' '
  print *,'Check for isolated nodes ========'
  n_iso = 0
  checkloop: do n=1,nod2D
     is = ssh_stiff%rowptr(n)
     ie = ssh_stiff%rowptr(n+1) -1

     node_neighb_part(1:ie-is) = part(ssh_stiff%colind(is:ie))
     if (count(node_neighb_part(1:ie-is) == part(n)) <= 1) then

        n_iso = n_iso+1
        print *,'Isolated node',n, 'in partition', part(n)
        print *,'Neighbouring nodes are in partitions',  node_neighb_part(1:ie-is)

        ! count the adjacent nodes of the other PEs

        np=1
        ne_part(1) = node_neighb_part(1)
        ne_part_num(1) = 1
        ne_part_load(1,1) = nod_per_partition(1,ne_part(1)) + 1
        ne_part_load(2,1) = nod_per_partition(2,ne_part(1)) + nlevels_nod2D(n)

        do i=1,ie-is
           if (node_neighb_part(i)==part(n)) cycle
           already_counted = .false.
           do k=1,np
              if (node_neighb_part(i) == ne_part(k)) then
                 ne_part_num(k) = ne_part_num(k) + 1
                 already_counted = .true.
                 exit
              endif
           enddo
           if (.not. already_counted) then
              np = np+1
              ne_part(np) = node_neighb_part(i)
              ne_part_num(np) = 1
              ne_part_load(1,np) = nod_per_partition(1,ne_part(np)) + 1
              ne_part_load(2,np) = nod_per_partition(2,ne_part(np)) + nlevels_nod2D(n)
           endif
        enddo

        ! Now, check for two things: The load balance, and if 
        ! there is more than one node of that partition.
        ! Otherwise, it would become isolated again.

        ! Best choice would be the partition with most adjacent nodes (edgecut!)
        ! Choose, if it does not decrease the load balance. 
        !        (There might be two partitions with the same number of adjacent
        !         nodes. Don't care about this here)

        kmax = maxloc(ne_part_num(1:np),1)

        if (ne_part_num(kmax) <= 1) then
           ! No chance - this is probably a boundary
           ! node that has only two neighbors.
           cycle checkloop
        endif

        if  (ne_part_load(1,kmax) <= max_nod_per_part(1) .and. &
             ne_part_load(2,kmax) <= max_nod_per_part(2) ) then
           k = kmax
        else
           ! Don't make it too compicated. Reject partitions that have only one
           ! adjacent node. Take the next not violating the load balance.
           found_part = .false.
           do k=1,np
              if (ne_part_num(k)==1 .or. k==kmax) cycle

              if  (ne_part_load(1,k) <= max_nod_per_part(1) .and. &
                   ne_part_load(2,k) <= max_nod_per_part(2) ) then

                 found_part = .true.
                 exit
              endif
           enddo

           if (.not. found_part) then
              ! Ok, don't think to much. Simply go for minimized edge cut.
              k = kmax
           endif
        endif

        ! Adjust the load balancing

        nod_per_partition(1,ne_part(k)) = nod_per_partition(1,ne_part(k)) + 1 
        nod_per_partition(2,ne_part(k)) = nod_per_partition(2,ne_part(k)) + nlevels_nod2D(n)
        nod_per_partition(1,part(n))    = nod_per_partition(1,part(n)) - 1
        nod_per_partition(2,part(n))    = nod_per_partition(2,part(n)) - nlevels_nod2D(n)

        ! And, finally, move nod n to other partition        
        part(n) = ne_part(k)
        print *,'Node',n,'is moved to part',part(n)
     endif
  enddo checkloop

  deallocate(ne_part, ne_part_num, ne_part_load)

  print *,'=== LOAD BALANCING ==='
  print *,'2D nodes: min, aver, max per part',min_nod_per_part(1), &
       average_nod_per_part(1),max_nod_per_part(1)

  write(*,"('2D nodes: percent min, aver, max ',f8.3,'%, 100%, ',f8.3,'%')") &
       100.*real(min_nod_per_part(1)) / real(average_nod_per_part(1)), &
       100.*real(max_nod_per_part(1)) / real(average_nod_per_part(1))

  print *,'3D nodes: Min, aver, max per part',min_nod_per_part(2), &
       average_nod_per_part(2),max_nod_per_part(2)
  write(*,"('3D nodes: percent min, aver, max ',f8.3,'%, 100%, ',f8.3,'%')") &
       100.*real(min_nod_per_part(2)) / real(average_nod_per_part(2)), &
       100.*real(max_nod_per_part(2)) / real(average_nod_per_part(2))

end subroutine check_partitioning


