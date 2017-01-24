! Driving routine. The distributed mesh information and mesh proper 
! are read from files.
! Auxiliary arrays with mesh information are assembled.
! At the beginning of each routine I list arrays it initializes.
! Array sizes vary (sometimes we need only myDim, yet sometimes more)! 
! S. Danilov, 2012
SUBROUTINE mesh_setup
USE g_parsup
USE g_ROTATE_grid
IMPLICIT NONE
      call set_mesh_transform_matrix  !(rotated grid)
      call read_mesh
      call set_par_support
      call find_levels
      call test_tri
      call load_edges
      call find_neighbors
      call mesh_areas
      call mesh_auxiliary_arrays
END SUBROUTINE mesh_setup
!======================================================================
! Reads distributed mesh
SUBROUTINE read_mesh
USE o_PARAM
USE g_CONFIG
USE o_MESH
USE o_ARRAYS
USE g_PARSUP
USE g_rotate_grid 
IMPLICIT NONE

 Integer        :: n, m, fileID, ind, nini, nend, n1, n2, n3, n4
 integer        :: vert_nodes(100)
 real(kind=WP)  :: x, y
 character*10   :: mype_string,npes_string
 character*200   :: file_name
 character*200   :: dist_mesh_dir
 integer, allocatable :: mapping(:)
 
 write(mype_string,'(i5.5)') mype  
 write(npes_string,"(I10)") npes
 dist_mesh_dir=trim(meshpath)//'dist_'//trim(ADJUSTL(npes_string))//'/'
 
 
         !=======================
	 ! rank partitioning vectors
	 !=======================
         file_name=trim(dist_mesh_dir)//'rpart.out' 
	 fileID=10+mype
	 open(fileID, file=trim(file_name)) 
         allocate(part(npes+1))
	 
	 read(fileID,*) n
	 if (n.ne.npes) then
	  write(*,*) 'NPES does not coincide with that of the mesh'
	  call par_ex(1)
	  STOP
	 end if
	 part(1)=1
	 read(fileID,*) part(2:npes+1)
	 DO n=2, npes+1
	 part(n)=part(n-1)+part(n)
         END DO
	 close(fileID)
	 if (mype==0) write(*,*) mype,'rpart is read'
	 !===========================
	 ! Lists of nodes and elements 
	 ! in global indexing. Not everything
	 ! is needed
	 !===========================
 
         file_name=trim(dist_mesh_dir)//'my_list'//trim(mype_string)//'.out'  
         fileID=10+mype  
    
         open(fileID, file=trim(file_name))
         read(fileID,*) n
	 
	 read(fileID,*) myDim_nod2D
	 read(fileID,*) eDim_nod2D
	 allocate(myList_nod2D(myDim_nod2D+eDim_nod2D)) 	 
	 read(fileID,*) myList_nod2D
	 
	 read(fileID,*) myDim_elem2D
	 read(fileID,*) eDim_elem2D
	 read(fileID,*) eXDim_elem2D
	 allocate(myList_elem2D(myDim_elem2D+eDim_elem2D+eXDim_elem2D))
	 read(fileID,*) myList_elem2D
	
	 read(fileID,*) myDim_edge2D
	 read(fileID,*) eDim_edge2D
         allocate(myList_edge2D(myDim_edge2D+eDim_edge2D))
	 read(fileID,*) myList_edge2D ! m

	 close(fileID)
	 if (mype==0) write(*,*) 'myLists are read'
         !==============================
         ! Allocate mapping array
         !==============================
	 nod2D=part(npes+1)-1
	 file_name=trim(meshpath)//'elem2d.out'
	 open(fileID, file=file_name)
	 read(fileID,*) elem2D
	 close(fileID)
	 Allocate(mapping(elem2D))  
         mapping=0 
	                       
	 DO n=1, myDim_nod2D+eDim_nod2D
	 mapping(myList_nod2D(n))=n
	 END DO
	 
	 !==============================
         ! read 2d node data
	 !==============================
	 
	 ALLOCATE(coord_nod2D(2,myDim_nod2D+eDim_nod2D))
	 
	 file_name=trim(meshpath)//'nod2d.out' 
	 open(fileID, file=file_name)
	 read(fileID,*) n      ! nod2D, we know it already
	                      
         DO n=1,nod2D
            read(fileID,*) m, x, y
            if (mapping(n)>0) then
               if (force_rotation) then
                  call g2r(x*rad, y*rad, x, y)
                  x=x/rad
                  y=y/rad
               end if
               coord_nod2D(1,mapping(n))=x*rad
               coord_nod2D(2,mapping(n))=y*rad
	 end if
	 END DO
         close(fileID)
	
	 !==============================
         ! read depth data
	 !==============================
	 
	 ALLOCATE(depth(myDim_nod2D+eDim_nod2D))
	 
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
         if (x>zbar(5)) x=zbar(5) !TODO KK threshholding for depth
	 if (mapping(n)>0) then
	  depth(mapping(n))=x
	 end if
	 END DO
	 mapping(1:nod2D)=0
	 close(fileID)
	 
         if(depth(2)>0) depth=-depth  ! depth is negative
	 
	 !==============================
         ! read 2d elem data
	 !==============================
	 file_name=trim(meshpath)//'elem2d.out' 
	 open(fileID, file=file_name)
	 
	 ALLOCATE(elem2D_nodes(3, myDim_elem2D))
         DO n=1, myDim_elem2D
	 mapping(myList_elem2D(n))=n
	 END DO
	 read(fileID,*) elem2d    
	 DO n=1,elem2D
	 read(fileID,*) n1, n2, n3
	 if (mapping(n)>0) then
	 elem2D_nodes(1,mapping(n))=n1
	 elem2D_nodes(2,mapping(n))=n2
	 elem2D_nodes(3,mapping(n))=n3
	 end if
	 END DO
	 close(fileID)
	 ! nodes in elem2d are in global numbering. convert to local:
 
	 mapping(1:elem2D)=0
	 DO n=1, myDim_nod2D+eDim_nod2D
	 mapping(myList_nod2D(n))=n
	 END DO
	 DO n=1, myDim_elem2D
	    DO m=1,3
               n1=elem2D_nodes(m,n)	 
	       elem2D_nodes(m,n)=mapping(n1)	 
	    END DO
	 END DO   
	 mapping(1:nod2D)=0
	 
	 if (mype==0) write(*,*) 'elements are read' 
	 
	 ! ==============================
         ! Communication information
         ! ==============================
         file_name=trim(dist_mesh_dir)//'com_info'//trim(mype_string)//'.out'  
         fileID=10+mype  
         open(fileID, file=file_name)
	 read(fileID,*)  n
	 read(fileID,*) com_nod2D%rPEnum
	 ALLOCATE(com_nod2D%rPE(com_nod2D%rPEnum))
	 read(fileID,*) com_nod2D%rPE
	 ALLOCATE(com_nod2D%rptr(com_nod2D%rPEnum+1))
	 read(fileID,*) com_nod2D%rptr
	 ALLOCATE(com_nod2D%rlist(eDim_nod2D))
	 read(fileID,*) com_nod2D%rlist
	 
	 read(fileID,*) com_nod2D%sPEnum
	 ALLOCATE(com_nod2D%sPE(com_nod2D%sPEnum))
	 read(fileID,*) com_nod2D%sPE
	 ALLOCATE(com_nod2D%sptr(com_nod2D%sPEnum+1))
	 read(fileID,*) com_nod2D%sptr
	 n=com_nod2D%sptr(com_nod2D%sPEnum+1)-1
	 ALLOCATE(com_nod2D%slist(n))
	 read(fileID,*) com_nod2D%slist
	 
	 read(fileID,*) com_elem2D%rPEnum
	 ALLOCATE(com_elem2D%rPE(com_elem2D%rPEnum))
	 read(fileID,*) com_elem2D%rPE
	 ALLOCATE(com_elem2D%rptr(com_elem2D%rPEnum+1))
	 read(fileID,*) com_elem2D%rptr
	 ALLOCATE(com_elem2D%rlist(eDim_elem2D))
	 read(fileID,*) com_elem2D%rlist
	 
	 read(fileID,*) com_elem2D%sPEnum
	 ALLOCATE(com_elem2D%sPE(com_elem2D%sPEnum))
	 read(fileID,*) com_elem2D%sPE
	 ALLOCATE(com_elem2D%sptr(com_elem2D%sPEnum+1))
	 read(fileID,*) com_elem2D%sptr
	 n=com_elem2D%sptr(com_elem2D%sPEnum+1)-1
	 ALLOCATE(com_elem2D%slist(n))
	 read(fileID,*) com_elem2D%slist
	 
	 read(fileID,*) com_elem2D_full%rPEnum
	 ALLOCATE(com_elem2D_full%rPE(com_elem2D_full%rPEnum))
	 read(fileID,*) com_elem2D_full%rPE
	 ALLOCATE(com_elem2D_full%rptr(com_elem2D_full%rPEnum+1))
	 read(fileID,*) com_elem2D_full%rptr
	 ALLOCATE(com_elem2D_full%rlist(eDim_elem2D+eXDim_elem2D))
	 read(fileID,*) com_elem2D_full%rlist
	 
	 read(fileID,*) com_elem2D_full%sPEnum
	 ALLOCATE(com_elem2D_full%sPE(com_elem2D_full%sPEnum))
	 read(fileID,*) com_elem2D_full%sPE
	 ALLOCATE(com_elem2D_full%sptr(com_elem2D_full%sPEnum+1))
	 read(fileID,*) com_elem2D_full%sptr
	 n=com_elem2D_full%sptr(com_elem2D_full%sPEnum+1)-1
	 ALLOCATE(com_elem2D_full%slist(n))
	 read(fileID,*) com_elem2D_full%slist

	 read(fileID,*) com_edge2D%rPEnum
	 ALLOCATE(com_edge2D%rPE(com_edge2D%rPEnum))
	 read(fileID,*) com_edge2D%rPE
	 ALLOCATE(com_edge2D%rptr(com_edge2D%rPEnum+1))
	 read(fileID,*) com_edge2D%rptr
	 ALLOCATE(com_edge2D%rlist(eDim_edge2D))
	 read(fileID,*) com_edge2D%rlist
	 
	 read(fileID,*) com_edge2D%sPEnum
	 ALLOCATE(com_edge2D%sPE(com_edge2D%sPEnum))
	 read(fileID,*) com_edge2D%sPE
	 ALLOCATE(com_edge2D%sptr(com_edge2D%sPEnum+1))
	 read(fileID,*) com_edge2D%sptr
	 n=com_edge2D%sptr(com_edge2D%sPEnum+1)-1
	 ALLOCATE(com_edge2D%slist(n))
	 read(fileID,*) com_edge2D%slist
	 close(fileID)
	 if (mype==0) write(*,*) 'communication arrays are read'
	 deallocate(mapping)
 END subroutine  read_mesh
!============================================================ 
subroutine find_levels
USE o_MESH
USE o_PARAM
USE g_PARSUP
  use g_config
!
IMPLICIT NONE
!
character*200 :: file_name
integer :: n,fileID=111,x
integer, allocatable                  :: mapping(:)
allocate(nlevels(myDim_elem2D+eDim_elem2D+eXDim_elem2D))

file_name=trim(meshpath)//'elvls.out'
open(fileID, file=file_name)
allocate(mapping(elem2D))
mapping=0 
DO n=1,myDim_elem2D+eDim_elem2D+eXDim_elem2D
      mapping(myList_elem2D(n))=n 
END DO
DO n=1,elem2D
read(fileID,*) x
if (mapping(n)>0) then
      nlevels(mapping(n)) = x
endif
ENDDO
deallocate(mapping)
close(fileID)

allocate(nlevels_nod2D(myDim_nod2D+eDim_nod2D))
file_name=trim(meshpath)//'nlvls.out'
open(fileID, file=file_name)
allocate(mapping(nod2D))
mapping=0 
DO n=1,myDim_nod2D+eDim_nod2D
      mapping(myList_nod2D(n))=n
END DO
DO n=1,nod2D
read(fileID,*) x
if (mapping(n)>0) then
      nlevels_nod2D(mapping(n)) = x
endif
ENDDO
deallocate(mapping)
close(fileID)

if (mype==0) then
write(*,*) '========================='
write(*,*) 'Mesh is read : ', 'nod2D=', nod2D,' elem2D=', elem2D, ' nl=', nl
write(*,*) 'Min/max depth on mype: ', mype, -zbar(minval(nlevels)),-zbar(maxval(nlevels))
write(*,*) '3D tracer nodes on mype ', mype, sum(nlevels_nod2d)-(myDim_elem2D+eDim_elem2D)
write(*,*) 'Further info on mype (1)', mype, minval(nlevels),maxval(nlevels)
write(*,*) 'Further info on mype (2)', mype, minval(nlevels_nod2d),maxval(nlevels_nod2d)
write(*,*) '========================='
endif
END SUBROUTINE find_levels
!===========================================================================
SUBROUTINE test_tri
USE o_MESH
USE o_PARAM
USE g_PARSUP
USE g_CONFIG
IMPLICIT NONE
! Check the order of nodes in triangles; correct it if necessary to make
! it same sense (clockwise) 
real(kind=WP)   ::  a(2), b(2), c(2),  r
integer         ::  n, nx, elnodes(3)

   
   DO n=1, myDim_elem2D
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

END SUBROUTINE  test_tri
!=========================================================================
SUBROUTINE load_edges
USE o_MESH
USE o_PARAM
USE g_PARSUP
USE g_CONFIG
IMPLICIT NONE
integer                               :: counter, n, k,q
integer                               :: elems(2), elem
integer                               :: elnodes(3), ed(2), eledges(3)
integer, allocatable                  :: mapping(:), aux(:)         
integer              :: n1, n2, m
! Edge array is already available (we computed it in the init phase)
! 
! (a) Read list of edges and tri containing them from file 
!
open(11, file=trim(meshpath)//'edgenum.out')
 read(11,*) edge2D
 read(11,*) edge2D_in
 close(11) 
open(10, file=trim(meshpath)//'edges.out')
open(12, file=trim(meshpath)//'edge_tri.out')
allocate(edges(2,myDim_edge2D+eDim_edge2D))
allocate(edge_tri(2,myDim_edge2D+eDim_edge2D))
allocate(mapping(edge2D))
 mapping=0
   DO n=1,myDim_edge2D+eDim_edge2D 
      mapping(myList_edge2D(n))=n 
   END DO
   
 DO n=1,edge2D
   read(10,*) ed
   read(12,*) elems
   if(mapping(n)>0) then
   edges(:,mapping(n))=ed
   edge_tri(:,mapping(n))=elems
   end if
 END DO
   close(10)
   close(12) 
 
   ! =========
   ! Local numbers for nodes
   ! =========
   mapping=0 
   DO n=1,myDim_nod2D+eDim_nod2D 
     mapping(myList_nod2D(n))=n 
   END DO
   DO n=1,myDim_edge2D+eDim_edge2D
      ed=edges(:,n)
      edges(1,n)=mapping(ed(1))
      edges(2,n)=mapping(ed(2))
   END DO
   ! =========
   ! Local numbers for elements
   ! =========
   mapping(1:nod2D)=0 
   DO n=1,myDim_elem2D+eDim_elem2D
     mapping(myList_elem2D(n))=n 
   END DO
   DO n=1,myDim_edge2D+eDim_edge2D
      ed=edge_tri(:,n)
      edge_tri(1,n)=mapping(ed(1))   ! Neighbor elements may appear
      if(ed(2)>0) then               ! with eDim edges
      edge_tri(2,n)=mapping(ed(2))   
      else
      edge_tri(2,n)=(ed(2))
      end if
   END DO
 
 deallocate(mapping)   
! Now the only way to check whether an edge is on boundary is 
! through myList_edge2D(n):  myList_edge2D(n)>edge2D_in == boundary edge

! (b) We need an array inverse to edge_tri listing edges
! of a given triangle 
allocate(elem_edges(3,myDim_elem2D))
allocate(aux(myDim_elem2D))
aux=0
DO n=1, myDim_edge2D+eDim_edge2D
   DO k=1,2
      q=edge_tri(k,n)   ! triangle number
	  if((q>0).and.(q<=myDim_elem2D)) then
	  aux(q)=aux(q)+1
	  elem_edges(aux(q),q)=n
	  end if
   END DO
END DO
deallocate(aux)
! The edges in this list should be ordered so that they
! are listed in the same rotation sense as nodes.
DO elem=1,myDim_elem2D
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
! The edge and elem lists agree in the sense that edge1 does not
! contain node 1 and so on
END SUBROUTINE load_edges
!===========================================================================
SUBROUTINE find_neighbors
! For each element three its element neighbors are found
! For each node the elements containing it are found
! Allocated are:
! elem_neighbors(3,myDim_elem2D)
! nod_in_elem2D_num(myDim_nod2D)
! nod_in_elem2D(:, myDim_nod2D)
! 

USE o_PARAM
USE o_MESH
USE g_PARSUP
  use g_comm_auto
implicit none
integer               :: elem, eledges(3), elem1, j, n, node, enum,elems(3),count1,count2,exit_flag,i,nz
integer, allocatable  :: temp_i(:)
integer               :: mymax(npes), rmax(npes)    
 ! =============
 ! elem neighbors == those that share edges
 ! =============
   allocate(elem_neighbors(3,myDim_elem2D))
   elem_neighbors=0
  
DO elem=1,myDim_elem2D
   eledges=elem_edges(:,elem)
   DO j=1,3
   elem1=edge_tri(1,eledges(j))
   if(elem1==elem) elem1=edge_tri(2,eledges(j))
   elem_neighbors(j,elem)=elem1
   END DO
END DO
 ! =============
 ! Node neighbourhood
 ! == elements that contain node n
 ! We need eDim neighborhood too for MUSCL advection. 
 ! And we already have the place allocated for all 
 ! these neighbor elements: it is eDim_elem2D+eXDim_elem2D
 ! =============	 
 allocate(nod_in_elem2D_num(myDim_nod2D+eDim_nod2D)) 
 nod_in_elem2D_num=0
 do n=1,myDim_elem2D
    do j=1,3
    node=elem2D_nodes(j,n)
    if (node>myDim_nod2D) cycle
    nod_in_elem2D_num(node)=nod_in_elem2D_num(node)+1
    end do
 end do
 mymax=0
 rmax=0
 mymax(mype+1)=maxval(nod_in_elem2D_num(1:myDim_nod2D))
 call MPI_AllREDUCE( mymax, rmax, &
       npes, MPI_INTEGER,MPI_SUM, &
       MPI_COMM_WORLD, MPIerr)
 
 allocate(nod_in_elem2D(maxval(rmax),myDim_nod2D+eDim_nod2D))
 nod_in_elem2D=0
 
 nod_in_elem2D_num=0
 do n=1,myDim_elem2D   
    do j=1,3
    node=elem2D_nodes(j,n)
    if (node>myDim_nod2D) cycle 
    nod_in_elem2D_num(node)=nod_in_elem2D_num(node)+1
    nod_in_elem2D(nod_in_elem2D_num(node),node)=n
    end do
 end do


 call exchange_nod(nod_in_elem2D_num)
 allocate (temp_i(myDim_nod2D+eDim_nod2D))
 temp_i=0
 DO n=1, maxval(rmax)
       ! Exchange global element numbers
       do j=1,myDim_nod2D
         if (nod_in_elem2D(n,j)>0) temp_i(j)=myList_elem2D(nod_in_elem2D(n,j))
       enddo
       call exchange_nod(temp_i)
       nod_in_elem2D(n,:)=temp_i
 END DO
 deallocate(temp_i)
 ! Substitute back local element numbers
 allocate(temp_i(elem2D))
 temp_i=0
 Do n=1, myDim_elem2D+eDim_elem2D+eXDim_elem2D
    temp_i(myList_elem2D(n))=n
 END DO
 DO n=1, myDim_nod2D+eDim_nod2D     
    DO j=1, nod_in_elem2D_num(n)
       nod_in_elem2D(j,n)=temp_i(nod_in_elem2D(j,n))
    END DO
 END DO
 deallocate(temp_i)       
 ! Among elem_neighbors there can be negative numbers. These correspond to 
 ! boundary elements for which neighbours are absent. However, an element 
 ! should have at least two valid neighbors
 ! ============
 ! Test that there are at least two neighbors at the surface:
 ! ============ 
DO elem=1,myDim_elem2D
   elem1=0
   DO j=1,3
   if(elem_neighbors(j,elem)>0) elem1=elem1+1
   END DO
   if (elem1<2) then
    write(*,*) 'Insufficient number of neighbors ', myList_elem2D(elem)
    call par_ex(1)
    STOP
   end if
END DO    

END SUBROUTINE find_neighbors
!==========================================================================
subroutine edge_center(n1, n2, x, y)
USE o_MESH
USE o_PARAM
USE g_CONFIG 
implicit none
integer      :: n1, n2   ! nodes of the edge
real(kind=WP) :: x, y, a(2), b(2)

a=coord_nod2D(:,n1)
b=coord_nod2D(:,n2)
if(a(1)-b(1)>cyclic_length/2.0) a(1)=a(1)-cyclic_length
if(a(1)-b(1)<-cyclic_length/2.0) b(1)=b(1)-cyclic_length
x=0.5_WP*(a(1)+b(1))
y=0.5_WP*(a(2)+b(2))
end subroutine edge_center
!==========================================================================
subroutine elem_center(elem, x, y)
USE o_MESH
USE o_PARAM
USE g_CONFIG  
implicit none
integer      :: elem, elnodes(3), k    
real(kind=WP) :: x, y, ax(3), amin

   elnodes=elem2D_nodes(:,elem)
   ax=coord_nod2D(1, elnodes)
   amin=minval(ax)
   DO k=1,3
   if(ax(k)-amin>=cyclic_length/2.0) ax(k)=ax(k)-cyclic_length
   if(ax(k)-amin<-cyclic_length/2.0) ax(k)=ax(k)+cyclic_length
   END DO
   x=sum(ax)/3.0_WP
   y=sum(coord_nod2D(2,elnodes))/3.0_WP
   
end subroutine elem_center
!==========================================================================
SUBROUTINE mesh_areas
USE o_MESH
USE o_PARAM
USE g_PARSUP
USE g_ROTATE_GRID
  use g_comm_auto
IMPLICIT NONE
! Collects auxilliary information on the mesh
! Allocated and filled in are:
! elem_area(myDim_elem2D)
! area(nl, myDim_nod2D)


integer                                   :: n,j,q, elnodes(3), ed(2), elem, nz
real(kind=8)	                          :: a(2), b(2), ax, ay, lon, lat, vol
real(kind=WP), allocatable,dimension(:)   :: work_array
	 
 allocate(elem_area(myDim_elem2D+eDim_elem2D+eXDim_elem2D))
 !allocate(elem_area(myDim_elem2D))
 allocate(area(nl,myDim_nod2d+eDim_nod2D))   !! Extra size just for simplicity
                                             !! in some further routines
 allocate(area_inv(nl,myDim_nod2d+eDim_nod2D)) 
 allocate(mesh_resolution(myDim_nod2d+eDim_nod2D))
 ! ============
 ! The areas of triangles:
 ! ============
 DO n=1, myDim_elem2D
 !DO n=1, myDim_elem2D+eDim_elem2D+eXDim_elem2D
    elnodes=elem2D_nodes(:,n)
    ay=sum(coord_nod2D(2,elnodes))/3.0_WP
    ay=cos(ay)
    if (cartesian) ay=1.0_WP
    a=coord_nod2D(:,elnodes(2))-coord_nod2D(:,elnodes(1))
    b=coord_nod2D(:,elnodes(3))-coord_nod2D(:,elnodes(1))
    if(a(1)>cyclic_length/2.) a(1)=a(1)-cyclic_length
    if(a(1)<-cyclic_length/2.) a(1)=a(1)+cyclic_length
    if(b(1)>cyclic_length/2.) b(1)=b(1)-cyclic_length
    if(b(1)<-cyclic_length/2.) b(1)=b(1)+cyclic_length
    a(1)=a(1)*ay
    b(1)=b(1)*ay
    elem_area(n)=0.5_WP*abs(a(1)*b(2)-b(1)*a(2))
 END DO
 call exchange_elem(elem_area)
 ! =============
 ! Scalar element 
 ! areas at different levels (there can be partly land)
 ! =============
 
 area=0.0_WP
 DO n=1, myDim_nod2D
    DO j=1,nod_in_elem2D_num(n)
       elem=nod_in_elem2D(j,n)
       DO nz=1,nlevels(elem)-1
       area(nz,n)=area(nz,n)+elem_area(elem)/3.0_WP
       END DO
    END DO
 END DO
 
 ! Only areas through which there is exchange are counted

 ! ===========
 ! Update to proper dimension
 ! ===========
 elem_area=elem_area*r_earth*r_earth
 area=area*r_earth*r_earth
 
 call exchange_nod(area)

do n=1,myDim_nod2d+eDim_nod2D
   do nz=1,nl
      if (area(nz,n) > 0._WP) then
         area_inv(nz,n) = 1._WP/area(nz,n)
      else
         area_inv(nz,n) = 0._WP
      end if
   end do
end do
 ! coordinates are in radians, edge_dxdy are in meters,
 ! and areas are in m^2
 

 allocate(work_array(myDim_nod2D))
 mesh_resolution=sqrt(area(1, :)/pi)*2._WP
 DO q=1, 3 !apply mass matrix N times to smooth the field
    DO n=1, myDim_nod2D
       vol=0._WP
       work_array(n)=0._WP
       DO j=1, nod_in_elem2D_num(n)
          elem=nod_in_elem2D(j, n)
          elnodes=elem2D_nodes(:,elem)
          work_array(n)=work_array(n)+sum(mesh_resolution(elnodes))/3._WP*elem_area(elem)
          vol=vol+elem_area(elem)
       END DO
       work_array(n)=work_array(n)/vol
    END DO
    DO n=1,myDim_nod2D
       mesh_resolution(n)=work_array(n)
    ENDDO
    call exchange_nod(mesh_resolution)
 END DO
 deallocate(work_array)

 vol=0.0_WP
 do n=1, myDim_nod2D
    vol=vol+area(1, n)
 end do
 ocean_area=0.0
 call MPI_AllREDUCE(vol, ocean_area, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
       MPI_COMM_WORLD, MPIerr)
if (mype==0) then
 write(*,*)  mype, 'Mesh statistics:'
 write(*,*)  mype, 'maxArea ',maxval(elem_area), '   MinArea ', minval(elem_area)
 write(*,*)  mype, 'maxScArea ',maxval(area(1,:)), &
            '   MinScArea ', minval(area(1,:))
 write(*,*)  mype, 'Edges:    ', edge2D, ' internal ', edge2D_in
 if (mype==0) then
    write(*,*) 'Total ocean area is: ', ocean_area, ' m^2'
 end if
endif
END SUBROUTINE mesh_areas

!===================================================================

SUBROUTINE mesh_auxiliary_arrays
! Collects auxiliary information needed to speed up computations 
! of gradients, div. This also makes implementation of cyclicity 
! much more straightforward
! Allocated and filled in are:
! edge_dxdy(2,myDim_edge2D+eDim_edge2D)
! edge_cross_dxdy(4,myDim_edge2D)
! gradient_node(6,myDim_elem2D)
! gradient_elem(6,myDim_elem2d)
! metric_factor(myDim_elem2D+eDim_elem2D)
! elem_cos(myDim_elem2D+eDim_elem2D)
! coriolis(myDim_elem2D)

USE o_MESH
USE o_PARAM
USE g_PARSUP
USE o_ARRAYS
USE g_ROTATE_grid
use g_comm_auto
IMPLICIT NONE

integer              :: n,j,q, elnodes(3), ed(2), elem, el(2), elnodes_(3)
real(kind=WP)	     :: a(2), b(2), ax, ay, dfactor, lon, lat
real(kind=WP)	     :: deltaX31, deltaX21, deltaY31, deltaY21
real(kind=WP)         :: x(3), y(3), cxx, cxy, cyy, d
real(kind=WP), allocatable :: center_x(:), center_y(:), temp(:) 


!real*8,allocatable :: arr2Dglobal(:,:) 
 
 allocate(edge_dxdy(2,myDim_edge2D+eDim_edge2D))
 allocate(edge_cross_dxdy(4,myDim_edge2D+eDim_edge2D))
 allocate(gradient_sca(6,myDim_elem2D))	 
 allocate(gradient_vec(6,myDim_elem2D))
 allocate(metric_factor(myDim_elem2D+eDim_elem2D+eXDim_elem2D))
 allocate(elem_cos(myDim_elem2D+eDim_elem2D+eXDim_elem2D))
 allocate(coriolis(myDim_elem2D))
 allocate(coriolis_node(myDim_nod2D+eDim_nod2D))
 allocate(geo_coord_nod2D(2,myDim_nod2D+eDim_nod2D))
 
 allocate(center_x(myDim_elem2D+eDim_elem2D+eXDim_elem2D))
 allocate(center_y(myDim_elem2D+eDim_elem2D+eXDim_elem2D)) 
 

 ! ============
 ! coriolis
 ! ============
 DO n=1,myDim_nod2D+eDim_nod2D 
 call r2g(lon, lat, coord_nod2D(1,n), coord_nod2D(2,n))
 coriolis_node(n)=2*omega*sin(lat)	 
 END DO
 DO n=1,myDim_nod2D+eDim_nod2D 
 call r2g(lon, lat, coord_nod2D(1,n), coord_nod2D(2,n))
 geo_coord_nod2D(1,n)=lon
 geo_coord_nod2D(2,n)=lat	 
 END DO
 
 DO n=1,myDim_elem2D 
 call elem_center(n, ax, ay)
 call r2g(lon, lat, ax, ay)
 coriolis(n)=2*omega*sin(lat)	 
 END DO
 
 
 if(fplane) then 
 coriolis=2*omega*0.71
 end if
 
 ! ============
 ! cos on elements + metric factor (tan/R_earth) 
 ! ============
 DO n=1,myDim_elem2D
 call elem_center(n, ax, ay)
 center_x(n)=ax
 center_y(n)=ay
 elem_cos(n)=cos(ay)
 metric_factor=tan(ay)/r_earth
 END DO

 call exchange_elem(metric_factor)
 call exchange_elem(elem_cos)
 call exchange_elem(center_x)
 call exchange_elem(center_y)  
 if (cartesian) then
 elem_cos=1.0_WP
 metric_factor=0.0_WP
 end if
 
 ! ===========
 ! Distances along the edge
 ! We need them in radian measure!
 ! ===========
 DO n=1, myDim_edge2D+eDim_edge2D
 ed=edges(:,n)
 a=coord_nod2D(:,ed(2))-coord_nod2D(:, ed(1))
 if(a(1)>cyclic_length/2) a(1)=a(1)-cyclic_length
 if(a(1)<-cyclic_length/2) a(1)=a(1)+cyclic_length
      !a(1)=a(1)*aux_cos_edge(n)
      !a=a*r_earth
 edge_dxdy(:,n)=a
 END DO

 ! ===========
 ! Cross-distances for the edge
 ! They are in physical measure!
 ! ===========
 DO n=1, myDim_edge2D
 ed=edges(:,n)
 el=edge_tri(:,n)
 
 call elem_center(el(1), b(1), b(2))
 call edge_center(ed(1), ed(2), a(1), a(2))
 b=b-a
 
 if(b(1)>cyclic_length/2)  b(1)=b(1)-cyclic_length
 if(b(1)<-cyclic_length/2) b(1)=b(1)+cyclic_length
 
 b(1)=b(1)*elem_cos(el(1))
 b=b*r_earth
 edge_cross_dxdy(1:2,n)=b(1:2)
 
 if(el(2)>0) then
 call elem_center(el(2), b(1), b(2))
 b=b-a
 if(b(1)>cyclic_length/2) b(1)=b(1)-cyclic_length
 if(b(1)<-cyclic_length/2) b(1)=b(1)+cyclic_length
 
 b(1)=b(1)*elem_cos(el(2))
 b=b*r_earth
 edge_cross_dxdy(3:4,n)=b(1:2)
 else
 edge_cross_dxdy(3:4,n)=0.0
 end if
 END DO
    allocate(temp(myDim_edge2D+eDim_edge2D))
    do n=1,4  
       temp(1:myDim_edge2D)=edge_cross_dxdy(n,1:myDim_edge2D)
       call exchange_edge(temp)
       edge_cross_dxdy(n,:)=temp
    end do       

 ! ==========================
 ! Derivatives of scalar quantities
 ! ==========================

DO elem=1, myDim_elem2D
   elnodes=elem2D_nodes(:,elem)
   
   deltaX31=coord_nod2D(1,elnodes(3))-coord_nod2D(1,elnodes(1))
   if(deltaX31>cyclic_length/2) deltaX31=deltaX31-cyclic_length
   if(deltaX31<-cyclic_length/2) deltaX31=deltaX31+cyclic_length
   deltaX31=elem_cos(elem)*deltaX31
   
   deltaX21=coord_nod2D(1,elnodes(2))-coord_nod2D(1,elnodes(1))
   if(deltaX21>cyclic_length/2) deltaX21=deltaX21-cyclic_length
   if(deltaX21<-cyclic_length/2) deltaX21=deltaX21+cyclic_length
   deltaX21=elem_cos(elem)*deltaX21
   
   deltaY31=coord_nod2D(2,elnodes(3))-coord_nod2D(2,elnodes(1))
   deltaY21=coord_nod2D(2,elnodes(2))-coord_nod2D(2,elnodes(1))
   
   dfactor=-0.5_8*r_earth/elem_area(elem)
   gradient_sca(1,elem)=(-deltaY31+deltaY21)*dfactor
   gradient_sca(2,elem)=deltaY31*dfactor
   gradient_sca(3,elem)=-deltaY21*dfactor
   
   gradient_sca(4,elem)=(deltaX31-deltaX21)*dfactor
   gradient_sca(5,elem)=-deltaX31*dfactor
   gradient_sca(6,elem)=deltaX21*dfactor
END DO

 ! ==========================
 ! Derivatives of vector quantities
 ! Least squares interpolation is used
 ! ==========================

   DO elem=1,myDim_elem2D
             !elnodes=elem2D_nodes(:,elem)
      a(1)=center_x(elem)
      a(2)=center_y(elem)
      DO j=1,3
      el(1)=elem_neighbors(j,elem)
      if (el(1)>0) then
             !elnodes_=elem2D_nodes(:,el(1))
      b(1)=center_x(el(1))
      b(2)=center_y(el(1))
      x(j)=b(1)-a(1)
      if(x(j)>cyclic_length/2) x(j)=x(j)-cyclic_length
      if(x(j)<-cyclic_length/2) x(j)=x(j)+cyclic_length
      y(j)=b(2)-a(2)
      else
      ! Virtual element center is taken
      ed=edges(:,elem_edges(j,elem))
      call edge_center(ed(1), ed(2), b(1), b(2))
      x(j)=(b(1)-a(1))
      if(x(j)>cyclic_length/2)   x(j)=x(j)-cyclic_length
      if(x(j)<-cyclic_length/2)  x(j)=x(j)+cyclic_length
      x(j)=2*x(j)
      y(j)=2*(b(2)-a(2))
      end if
      END DO
      x=x*elem_cos(elem)*r_earth
      y=y*r_earth
      cxx=sum(x**2)
      cxy=sum(x*y)
      cyy=sum(y**2)
      d=cxy*cxy-cxx*cyy
	  ! coefficients to compute gradients of velocity
      gradient_vec(1:3,elem)=(cxy*y-cyy*x)/d
      gradient_vec(4:6,elem)=(cxy*x-cxx*y)/d
    END DO
    deallocate(center_y, center_x)

END SUBROUTINE mesh_auxiliary_arrays

!===================================================================

SUBROUTINE check_mesh_consistency
USE o_MESH
USE o_PARAM
USE g_PARSUP
USE g_ROTATE_GRID
  use g_comm_auto
IMPLICIT NONE
! Collects auxilliary information on the mesh
! Allocated and filled in are:
! elem_area(myDim_elem2D)
! area(nl, myDim_nod2D)


integer              :: nz, n, elem , elnodes(3)
real(kind=8)	     :: vol_n(nl), vol_e(nl), aux(nl)

   vol_n=0.
   vol_e=0.

   aux=0.
   do n=1, myDim_nod2D
      do nz=1, nlevels_nod2D(n)-1
         aux(nz)=aux(nz)+area(nz, n)
      end do
   end do
   call MPI_AllREDUCE(aux, vol_n, nl, MPI_DOUBLE_PRECISION, MPI_SUM, &
       MPI_COMM_WORLD, MPIerr)

   aux=0.
   do elem=1, myDim_elem2D
      elnodes=elem2D_nodes(:, elem)
      if (elnodes(1) > myDim_nod2D) CYCLE
      do nz=1, nlevels(elem)         
         aux(nz)=aux(nz)+elem_area(elem)
      end do
   end do
   call MPI_AllREDUCE(aux, vol_e, nl, MPI_DOUBLE_PRECISION, MPI_SUM, &
       MPI_COMM_WORLD, MPIerr)

if (mype==0) then
write(*,*) '***start level area_test***'
do nz=1, nl
   write(*,*) vol_n(nz), vol_e(nz)
end do
write(*,*) '***end level area_test***'
end if

!call par_ex
!stop
END SUBROUTINE check_mesh_consistency
