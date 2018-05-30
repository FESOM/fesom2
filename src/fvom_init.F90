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
  character(len=1000)   :: nmlfile  !> name of configuration namelist file
  integer               :: start_t, interm_t, finish_t, rate_t
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
!!$  call par_init
  call read_mesh_ini
  call system_clock(finish_t)
  print '("**** Reading initial mesh time = ",f12.3," seconds. ****")', &
       real(finish_t-interm_t)/real(rate_t)
  interm_t = finish_t
  call test_tri_ini
  call find_edges_ini
  call find_elem_neighbors_ini
  call find_levels

! NR Some arrays are not needed for partitioning, after setting up the grid
  deallocate(coord_nod2D)
  deallocate(edge_tri)
  deallocate(zbar,Z,nlevels,depth)
  call system_clock(finish_t)
  print '("**** Checking and initializing mesh time = ",f12.3," seconds. ****")', &
       real(finish_t-interm_t)/real(rate_t)
  interm_t = finish_t

  call stiff_mat_ini
  call set_par_support_ini
  call system_clock(finish_t)
  print '("**** Partitioning time = ",f12.3," seconds. ****")', &
       real(finish_t-interm_t)/real(rate_t)
  interm_t = finish_t
  call communication_ini
  call system_clock(finish_t)
  print '("**** Storing partitioned mesh time = ",f12.3," seconds. ****")', &
       real(finish_t-interm_t)/real(rate_t)
  print '("**** Total time = ",f12.3," seconds. ****")', &
       real(finish_t-start_t)/real(rate_t)
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
INTEGER               :: nq
INTEGER               :: n1,n2,n3
INTEGER               :: n, nz, exit_flag
REAL(kind=WP)         :: x1, x2
INTEGER	              :: tag
INTEGER, allocatable  :: elem_data(:)
INTEGER               :: i_error
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
  ALLOCATE(elem2D_nodes(4,elem2D))
  ALLOCATE(elem_data(4*elem2D))
  elem_data(:)=-1
  
  ! meshes with quads have 4 columns, but TsunAWI grids may be
  ! purely triangular, with 3 columns each. Test, how many
  ! columns there are!  
  read(21,*,iostat=i_error) elem_data(1:4*elem2D)
  if (i_error == 0) then      ! There is a fourth column => quad or mixed mesh (not working yet!)
     elem2D_nodes = reshape(elem_data, shape(elem2D_nodes))
  else     ! No fourth column => triangles only
     elem2D_nodes(1:3,:) = reshape(elem_data, shape(elem2D_nodes(1:3,:)))
     elem2D_nodes(4,:) = elem2D_nodes(1,:)
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
SUBROUTINE test_tri_ini
USE o_MESH
USE o_PARAM
USE g_CONFIG
IMPLICIT NONE
real(kind=WP)   ::  a(2), b(2), c(2),  r
integer         ::  n, nx, elnodes(3)

   
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
SUBROUTINE find_edges_ini
USE o_MESH
USE o_PARAM
USE g_PARSUP
USE g_CONFIG
IMPLICIT NONE
integer, allocatable                  :: aux1(:), ne_num(:), ne_pos(:,:)
!!$integer, allocatable, dimension(:,:)  :: auxne, auxnn
integer                               :: counter, counter_in, n, k, q
integer                               :: elem, elem1, elems(2), q1, q2
integer                               :: elnodes(4), ed(2), flag, eledges(4)
integer                               :: temp(100), node 
real(kind=WP)                         :: xc(2), xe(2), ax(3), amin
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
 open(10, file='edges.out')

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

 allocate(edges(2,edge2D), edge_tri(2, edge2D))
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
          write(10,'(4I10)') n, node, elems
       else if (flag==1) then
          write(10,'(4I10)') n, node, elems(1), -999
       else
          write(*,*) 'flag'
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
 close(10)
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
 allocate(elem_edges(4,elem2D))
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
subroutine find_levels
use g_config
use o_mesh
use g_parsup
implicit none
INTEGER :: nodes(3), elems(3), eledges(3)
integer :: elem, elem1, j, n, q, node, enum,count1,count2,exit_flag,i,nz,fileID=111
real*8 :: x,dmean
integer :: thers_lev=5
character*200 :: file_name
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
      nodes=elem2D_nodes(1:3,n)
      dmean=maxval(depth(nodes))
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
SUBROUTINE find_elem_neighbors_ini
! For each element three its element neighbors are found
USE o_MESH
USE g_PARSUP
implicit none
integer    :: elem, eledges(3), elem1, j, n, elnodes(3)
allocate(elem_neighbors(4,elem2D))
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
    elnodes=elem2D_nodes(1:3,n)
    nod_in_elem2D_num(elnodes)=nod_in_elem2D_num(elnodes)+1
 end do
 allocate(nod_in_elem2D(maxval(nod_in_elem2D_num),nod2D))
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
subroutine stiff_mat_ini
  use o_MESH
  
  !
  implicit none
  integer                :: i, j, n, q, el, elem_nodes_max, nod(4)
  integer, allocatable   :: num_ne(:), ne(:,:)
  !

  ssh_stiff%dim = nod2D   
  allocate(ssh_stiff%rowptr(nod2D+1))

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
  do n=1,nod2D
     ssh_stiff%colind(ssh_stiff%rowptr(n):ssh_stiff%rowptr(n+1)-1) = ne(1:num_ne(n),n)
  end do

  deallocate(num_ne, ne)

end subroutine stiff_mat_ini

!===================================================================
! Setup of communication arrays
subroutine communication_ini
  use o_MESH
  USE g_CONFIG
  USE g_PARSUP
  use omp_lib
  implicit none

  integer        :: n
  character*10   :: npes_string
  character*200  :: dist_mesh_dir
  LOGICAL        :: L_EXISTS

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
     call communication_nodn
     call communication_elemn

     call save_dist_mesh         ! Write out communication file com_infoxxxxx.out
  end do
!$OMP END DO
!$OMP END PARALLEL

  deallocate(elem_neighbors)
  deallocate(elem_edges)
  deallocate(part)
  write(*,*) 'Communication arrays have been set up'   
end subroutine communication_ini
