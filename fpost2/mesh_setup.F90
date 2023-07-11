subroutine ocean_mesh_setup
  use o_param
  use o_elements
  use o_mesh
  use g_config
  use g_rotate_grid
  implicit none

  domain_length=domain_length*rad
  if(rotated_grid) call calculate_rotate_matrix

  call read_2Dmesh
  call read_3Dmesh
  call mesh_scaling                         ! long., lat. are transf. into rad
  write(*,*) 'mesh_scaling: DONE'
  call standard_element_definition_2D       ! Basis functions and scalar
  write(*,*) 'standard_element_definition_2D: DONE'
  call basisfunctions_2D            
  write(*,*) 'basisfunctions_2D: DONE'
  call build_nghbr_arrays                   ! Builds arrays nod_in_elem2D
  write(*,*) 'build_nghbr_arrays: DONE'
  call test_tri
  call load_edges
  call find_neighbors
  call nodal_areas
  call mesh_auxiliary_arrays
end subroutine ocean_mesh_setup
!
!---------------------------------------------------------------------------
!
subroutine mesh_scaling  
  !
  ! Transforms degrees in rad in coord_nod2D(2,nod2D)
  ! Constructs the arrays cos_elem2D(elem2D)
  ! Constructs num_layers_below_nod2D, 
  ! and does transform to rad in coord_nod3D

  use g_config
  use o_param
  use o_MESH
  use o_ELEMENTS
  use g_rotate_grid

  implicit none

  integer         :: i,j,ind2,ind3
  integer         :: n,  node, nodeup
  real(kind=8)    :: lon,lat,rlon,rlat

  ! =======================
  !  Lon and lat to radians
  ! =======================
  coord_nod2D(1,:)=coord_nod2D(1,:)*rad
  coord_nod2D(2,:)=coord_nod2D(2,:)*rad

  ! =======================
  ! Mean cos on 2D elements
  ! This sets spherical geometry!
  ! =======================
  allocate(cos_elem2D(elem2D))
  do i=1, elem2D  
     cos_elem2D(i)=cos(sum(coord_nod2D(2,elem2D_nodes(:,i)))/3.0) !sum(cos(coord_nod2D(2,elem2D_nodes(:,i))))/3.0
  end do
  if(cartesian) cos_elem2D=1.0
end subroutine mesh_scaling
!
!=====================================================================
!
subroutine standard_element_definition_2D
  !
  !    - BASISFUNCTIONS ON 2D STANDARD ELEMENTS 
  !         stdbafunc(1)=1-x-y  stdbafunc(2)=x  stdbafunc(3)=y
  !
  use o_ELEMENTS
  use o_mesh
  !
  implicit none
  !
  integer :: i,j
  !
  Vol2D =  1.0_8/2.0_8               ! Vol2D = < 1,1>
  !        
  sProd_2Di = 1.0_8/6.0_8            ! <stdbafunc(i),1.>
  !
  allocate(sProd_2Dij(3,3))
  sProd_2Dij=1.0_8/24.0_8            ! <stdbafunc(i),stdbafunc(j)>
  do j=1,3
     sProd_2Dij(j,j)=1.0_8/12.0_8 
  end do

  ! Scalar products are only required as sProd_2D/Vol2D:
  sProd_2Dij=sProd_2Dij/Vol2D
  sProd_2Di=sProd_2Di/Vol2D

  !   derivative_stdbafu_x(i,j) = d(Fi(j))/dx(i) on the standard element 
  allocate(derivative_stdbafu_x_2D(2,3))
  !
  derivative_stdbafu_x_2D= 0.
  derivative_stdbafu_x_2D(:,1)= -1.
  derivative_stdbafu_x_2D(1,2)= 1.
  derivative_stdbafu_x_2D(2,3)= 1.
  !
end subroutine standard_element_definition_2D
!
!=====================================================================
! 
subroutine basisfunctions_2D
  use o_ELEMENTS
  use o_MESH
  implicit none
  !
  real(kind=8)                           :: DET2D
  real(kind=8), dimension(2,3)           :: derivative_locbafu_x_2D
  real(kind=8), dimension(2,2)           :: jacobian2D, jacobian2D_inv
  integer                                :: elem, i


  allocate(bafux_2d(3, elem2d), bafuy_2d(3, elem2d))
  allocate(voltriangle(elem2d))

  bafux_2d = 0.0
  bafuy_2d = 0.0
  voltriangle = 0.0

  do elem=1,elem2d
     call local_element_def_2D(elem, DET2D, derivative_locbafu_x_2D)
     do i=1,3
        bafux_2d(i,elem) = derivative_locbafu_x_2D(1,i)
        bafuy_2d(i,elem) = derivative_locbafu_x_2D(2,i)
     enddo
     voltriangle(elem) = abs(DET2D) * Vol2D
  enddo
end subroutine basisfunctions_2D
!
!=====================================================================
! 
subroutine local_element_def_2D(element, DET, derivative_locbafu_x_2D)

  use o_ELEMENTS
  use o_MESH
  use o_param
  use g_config
  implicit none
  !
  integer, intent(IN)                        :: element
  real(kind=8), dimension(2,2)               :: jacobian2D
  real(kind=8), dimension(2,2)               :: jacobian2D_inv
  real(kind=8), intent(OUT)                  :: DET
  real(kind=8), dimension(2,3), intent(OUT)  :: derivative_locbafu_x_2D
  !
  real(kind=8), dimension(2,3)               :: local_cart
  real(kind=8), dimension(3,2)               :: der_transp
  integer                                    :: i, node
  real(kind=8)                               :: meancos
  !
  meancos=cos_elem2D(element)
  do i=1,3
     node=elem2D_nodes(i,element)
     !
     !  scaled cartesian coordinates
     !
     local_cart(1,i)=coord_nod2D(1,node) 
     local_cart(2,i)=r_earth * coord_nod2D(2,node) 
  end do
  !
  !  jacobian
  !
  do i=1,2
     jacobian2D(:,i)= local_cart(:,i+1)-local_cart(:,1)
     if (jacobian2D(1,i)> domain_length/2.0) jacobian2D(1,i)=jacobian2D(1,i)-domain_length
     if (jacobian2D(1,i)<-domain_length/2.0) jacobian2D(1,i)=jacobian2D(1,i)+domain_length
  end do
  jacobian2D(1,:)=jacobian2D(1,:)*meancos *r_earth
  !
  !  inverse of jacobian
  !
  call matrix_inverse_2x2(jacobian2D, jacobian2D_inv, DET)
  !
  der_transp=matmul(transpose(derivative_stdbafu_x_2D), jacobian2D_inv)
  derivative_locbafu_x_2D=transpose(der_transp)
  !
  !
end subroutine local_element_def_2D

!
!=======================================================================
!
subroutine  matrix_inverse_2x2 (A, AINV, DET)
  !
  !
  implicit none
  !
  real(kind=8), dimension(2,2), intent(IN)  :: A
  real(kind=8), dimension(2,2), intent(OUT) :: AINV
  real(kind=8), intent(OUT)                 :: DET
  !
  integer                                   :: i,j
  !
  DET  = A(1,1)*A(2,2) - A(1,2)*A(2,1)
  if ( DET .eq. 0.0 )  then
     do j=1,2
        write(*,*) (A(i,j),i=1,2)
     end do
     stop 'SINGULAR 2X2 MATRIX'
  else
     AINV(1,1) =  A(2,2)/DET
     AINV(1,2) = -A(1,2)/DET
     AINV(2,1) = -A(2,1)/DET
     AINV(2,2) =  A(1,1)/DET
  endif
end subroutine matrix_inverse_2x2
!
!=======================================================================
!
subroutine matrix_inverse_3x3(A, AINV, DET)
  !
  implicit none
  !
  real(kind=8), dimension(3,3), intent(IN)  :: A
  real(kind=8), dimension(3,3), intent(OUT) :: AINV
  real(kind=8), intent(OUT)                 :: DET
  !
  integer                                   :: i,j
  !
  AINV(1,1) =  A(2,2)*A(3,3) - A(3,2)*A(2,3)
  AINV(2,1) = -A(2,1)*A(3,3) + A(3,1)*A(2,3)
  AINV(3,1) =  A(2,1)*A(3,2) - A(3,1)*A(2,2)
  AINV(1,2) = -A(1,2)*A(3,3) + A(3,2)*A(1,3)
  AINV(2,2) =  A(1,1)*A(3,3) - A(3,1)*A(1,3)
  AINV(3,2) = -A(1,1)*A(3,2) + A(3,1)*A(1,2)
  AINV(1,3) =  A(1,2)*A(2,3) - A(2,2)*A(1,3)
  AINV(2,3) = -A(1,1)*A(2,3) + A(2,1)*A(1,3)
  AINV(3,3) =  A(1,1)*A(2,2) - A(2,1)*A(1,2)
  DET = A(1,1)*AINV(1,1) + A(1,2)*AINV(2,1) + A(1,3)*AINV(3,1)
  !
  if ( DET .eq. 0.0 )  then
     do j=1,3
        write(*,*) (A(i,j),i=1,3)
     end do
     stop 'SINGULAR 3X3 MATRIX'
  else
     AINV = AINV/DET
  endif
  !
end subroutine matrix_inverse_3x3
!
! =======================================================================
!
subroutine build_nghbr_arrays
  !
  ! Assembles additional arrays which list for each node the elements 
  ! containing the node and node neighbours

  use o_MESH
  use o_ELEMENTS
  implicit none

  integer                            :: j,k,m,n,a,tr(3),tet(4), counter, el,ml(1)
  integer, allocatable, dimension(:) :: ind
  integer, dimension(100)            :: AUX=0

  !--------------- 2D mesh:
  ! Builds nod_in_elem2D
     
  allocate(ind(nod2D))
  ind=0
  do j=1,elem2D
     tr=elem2D_nodes(:,j)
     ind(tr)=ind(tr)+1
  end do
  allocate(nod_in_elem2D(nod2D))
  nod_in_elem2D%nmb=ind(1:nod2D)    
  do j=1,nod2D   
     allocate(nod_in_elem2D(j)%addresses(ind(j)))
  end do
  ind=0
  do j=1,elem2D   
     tr=elem2D_nodes(:,j)
     ind(tr)=ind(tr)+1
     do k=1,3
        if(tr(k)<=nod2D) then 
           nod_in_elem2D(tr(k))%addresses(ind(tr(k)))=j
        end if
     end do
  end do
  ! the list of elements is ordered, and no sorting is needed
  ! Builds nghbr_nod2D
  allocate(nghbr_nod2D(nod2D))
  ind=0
  do j=1, nod2D
     counter=0
     do m=1,nod_in_elem2D(j)%nmb
        el=nod_in_elem2D(j)%addresses(m)
        do k=1, 3
           a=elem2D_nodes(k,el)       
           if (ind(a)==0) then  
              ind(a)=1 
              counter=counter+1         
              aux(counter)=a
           end if
        end do
     end do
     nghbr_nod2D(j)%nmb=counter
     allocate(nghbr_nod2D(j)%addresses(counter))

     ! we need to sort array aux(1:counter)
     do m=counter,1,-1
        ml=maxloc(aux(1:counter))
        a=ml(1)
        nghbr_nod2D(j)%addresses(m)=aux(a)
        ind(aux(a))=0
        aux(a)=-999
     end do
  end do
deallocate(ind)
end subroutine build_nghbr_arrays
!
! =======================================================================
!
SUBROUTINE test_tri
use o_MESH
use o_ELEMENTS
use o_param
use g_config
IMPLICIT NONE
! Check the order of nodes in triangles; correct it if necessary to make
! it same sense (clockwise) 
real(kind=8)   ::  a(2), b(2), c(2),  r
integer         ::  n, nx, elnodes(3)

   
   DO n=1, elem2D
      elnodes=elem2D_nodes(:,n)
	  
          a=coord_nod2D(:,elnodes(1))
	  b=coord_nod2D(:,elnodes(2))-a
	  c=coord_nod2D(:,elnodes(3))-a
          
	  if(b(1)>domain_length/2.) b(1)=b(1)-domain_length
          if(b(1)<-domain_length/2.) b(1)=b(1)+domain_length
	  if(c(1)>domain_length/2.) c(1)=c(1)-domain_length
          if(c(1)<-domain_length/2.) c(1)=c(1)+domain_length
	  
	    
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
!
! =======================================================================
!
SUBROUTINE load_edges
use o_MESH
use o_ELEMENTS
use o_param
use g_config
IMPLICIT NONE
integer                               :: counter, n, k,q
integer                               :: elems(2), elem
integer                               :: elnodes(3), ed(2), eledges(3)
integer, allocatable                  :: aux(:)         
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
allocate(edges(2,edge2D))
allocate(edge_tri(2,edge2D))

 DO n=1,edge2D
   read(10,*) ed
   read(12,*) elems
   edges(:,n)=ed
   edge_tri(:,n)=elems
 END DO
 close(10)
 close(12) 
 ! =========
 ! Local numbers for nodes
 ! =========
 DO n=1, edge2D
    ed=edges(:,n)
    edges(1,n)=ed(1)
    edges(2,n)=ed(2)
 END DO
 ! =========
 ! Local numbers for elements
 ! =========
 DO n=1, edge2D
    ed=edge_tri(:,n)
    edge_tri(1,n)=ed(1)  ! Neighbor elements may appear
    if(ed(2)>0) then     ! with eDim edges
    edge_tri(2,n)=ed(2)
    else
    edge_tri(2,n)=(ed(2))
    end if
 END DO
 
! Now the only way to check whether an edge is on boundary is 
! through myList_edge2D(n):  myList_edge2D(n)>edge2D_in == boundary edge

! (b) We need an array inverse to edge_tri listing edges
! of a given triangle 
allocate(elem_edges(3,elem2D))
allocate(aux(elem2D))
aux=0
DO n=1, edge2D
   DO k=1,2
      q=edge_tri(k,n)   ! triangle number
      if ((q>0).and.(q<=elem2D)) then
	 aux(q)=aux(q)+1
	 elem_edges(aux(q),q)=n
      end if
   END DO
END DO
deallocate(aux)
! The edges in this list should be ordered so that they
! are listed in the same rotation sense as nodes.
DO elem=1, elem2D
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
!
! =======================================================================
!
SUBROUTINE find_neighbors
! For each element three its element neighbors are found
! For each node the elements containing it are found
! Allocated are:
! elem_neighbors(3,elem2D)
use o_MESH
use o_ELEMENTS
use o_param
use g_config
use o_elements
implicit none
integer               :: elem, eledges(3), elem1, j, n, node, enum,elems(3),count1,count2,exit_flag,i,nz
! =============
! elem neighbors == those that share edges
! =============
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
 ! ============
 ! Test that there are at least two neighbors at the surface:
 ! ============ 
DO elem=1,elem2D
   elem1=0
   DO j=1,3
   if(elem_neighbors(j,elem)>0) elem1=elem1+1
   END DO
   if (elem1<2) then
    write(*,*) 'Insufficient number of neighbors ', elem
    STOP
   end if
END DO    
END SUBROUTINE find_neighbors
!
! =======================================================================
!
SUBROUTINE mesh_auxiliary_arrays
! Collects auxiliary information needed to speed up computations 
! of gradients, div. This also makes implementation of cyclicity 
! much more straightforward
! Allocated and filled in are:
! edge_dxdy(2,edge2D)
! edge_cross_dxdy(4,myDim_edge2D)
! gradient_node(6,elem2D)
! gradient_elem(6,elem2D)
! metric_factor(elem2D+eDim_elem2D)
! elem_cos(elem2D+eDim_elem2D)

use o_MESH
use o_ELEMENTS
use o_param
use g_config
use o_elements
use g_rotate_grid
IMPLICIT NONE

integer              :: n,j,q, elnodes(3), ed(2), elem, el(2), elnodes_(3)
real(kind=8)	     :: a(2), b(2), ax, ay, dfactor, lon, lat
real(kind=8)	     :: deltaX31, deltaX21, deltaY31, deltaY21
real(kind=8)         :: x(3), y(3), cxx, cxy, cyy, d
real(kind=8), allocatable :: center_x(:), center_y(:)


!real*8,allocatable :: arr2Dglobal(:,:) 
write(*,*) 'mesh_auxiliary_arrays starts:'   
 allocate(edge_dxdy(2,edge2D))
 allocate(edge_cross_dxdy(4,edge2D))
 allocate(gradient_sca(6,elem2D))	 
 allocate(gradient_vec(6,elem2D))
 allocate(metric_factor(elem2D))
 allocate(elem_cos(elem2D))
 
 allocate(center_x(elem2D))
 allocate(center_y(elem2D)) 
 
write(*,*) 'allocation complete...'  
! ============
! cos on elements
! ============
 DO n=1, elem2D
    elnodes=elem2D_nodes(:,n)
    ay=sum(coord_nod2D(2, elnodes))/3.0
    elem_cos(n)=cos(ay)
 END DO
 ! ===========
 ! Distances along the edge
 ! We need them in radian measure!
 ! ===========
 DO n=1, edge2D
 ed=edges(:,n)
 a=coord_nod2D(:,ed(2))-coord_nod2D(:, ed(1))
 if(a(1)>domain_length/2) a(1)=a(1)-domain_length
 if(a(1)<-domain_length/2) a(1)=a(1)+domain_length
      !a(1)=a(1)*aux_cos_edge(n)
      !a=a*r_earth
 edge_dxdy(:,n)=a
 END DO
write(*,*) 'edge_dydy complete...'   
 ! ===========
 ! Cross-distances for the edge
 ! They are in physical measure!
 ! ===========
 DO n=1, edge2D
 ed=edges(:,n)
 el=edge_tri(:,n)
 
 call elem_center(el(1), b(1), b(2))
 call edge_center(ed(1), ed(2), a(1), a(2))
 b=b-a
 
 if(b(1)>domain_length/2)  b(1)=b(1)-domain_length
 if(b(1)<-domain_length/2) b(1)=b(1)+domain_length
 
 b(1)=b(1)*elem_cos(el(1))
 b=b*r_earth
 edge_cross_dxdy(1:2,n)=b(1:2)
 
 if(el(2)>0) then
 call elem_center(el(2), b(1), b(2))
 b=b-a
 if(b(1)> domain_length/2) b(1)=b(1)-domain_length
 if(b(1)<-domain_length/2) b(1)=b(1)+domain_length
 
 b(1)=b(1)*elem_cos(el(2))
 b=b*r_earth
 edge_cross_dxdy(3:4,n)=b(1:2)
 else
 edge_cross_dxdy(3:4,n)=0.0
 end if
 END DO
write(*,*) 'edge_cross_dydy complete...'   
 ! ==========================
 ! Derivatives of scalar quantities
 ! ==========================

DO elem=1, elem2D
   elnodes=elem2D_nodes(:,elem)
   
   deltaX31=coord_nod2D(1,elnodes(3))-coord_nod2D(1,elnodes(1))
   if(deltaX31>domain_length/2) deltaX31=deltaX31-domain_length
   if(deltaX31<-domain_length/2) deltaX31=deltaX31+domain_length
   deltaX31=elem_cos(elem)*deltaX31
   
   deltaX21=coord_nod2D(1,elnodes(2))-coord_nod2D(1,elnodes(1))
   if(deltaX21>domain_length/2) deltaX21=deltaX21-domain_length
   if(deltaX21<-domain_length/2) deltaX21=deltaX21+domain_length
   deltaX21=elem_cos(elem)*deltaX21
   
   deltaY31=coord_nod2D(2,elnodes(3))-coord_nod2D(2,elnodes(1))
   deltaY21=coord_nod2D(2,elnodes(2))-coord_nod2D(2,elnodes(1))
   
   dfactor=-0.5_8*r_earth/voltriangle(elem)
   gradient_sca(1,elem)=(-deltaY31+deltaY21)*dfactor
   gradient_sca(2,elem)=deltaY31*dfactor
   gradient_sca(3,elem)=-deltaY21*dfactor
   
   gradient_sca(4,elem)=(deltaX31-deltaX21)*dfactor
   gradient_sca(5,elem)=-deltaX31*dfactor
   gradient_sca(6,elem)=deltaX21*dfactor
END DO
write(*,*) 'Derivatives of scalar quantities complete...'   
 ! ==========================
 ! Derivatives of vector quantities
 ! Least squares interpolation is used
 ! ==========================

   DO elem=1, elem2D
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
      if(x(j)>domain_length/2) x(j)=x(j)-domain_length
      if(x(j)<-domain_length/2) x(j)=x(j)+domain_length
      y(j)=b(2)-a(2)
      else
      ! Virtual element center is taken
      ed=edges(:,elem_edges(j,elem))
      call edge_center(ed(1), ed(2), b(1), b(2))
      x(j)=(b(1)-a(1))
      if(x(j)>domain_length/2)   x(j)=x(j)-domain_length
      if(x(j)<-domain_length/2)  x(j)=x(j)+domain_length
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
write(*,*) 'Derivatives of vector quantities complete...'   
    deallocate(center_y, center_x)
write(*,*) 'mesh_auxiliary_arrays DONE!'
END SUBROUTINE mesh_auxiliary_arrays
!
! =======================================================================
!
subroutine edge_center(n1, n2, x, y)
USE o_MESH
use o_elements
USE o_PARAM
USE g_CONFIG 
implicit none
integer      :: n1, n2   ! nodes of the edge
real(kind=8) :: x, y, a(2), b(2)

a=coord_nod2D(:,n1)
b=coord_nod2D(:,n2)
if(a(1)-b(1)>domain_length/2.0) a(1)=a(1)-domain_length
if(a(1)-b(1)<-domain_length/2.0) b(1)=b(1)-domain_length
x=0.5*(a(1)+b(1))
y=0.5*(a(2)+b(2))
end subroutine edge_center
!
! =======================================================================
!
subroutine elem_center(elem, x, y)
USE o_MESH
use o_elements
USE o_PARAM
USE g_CONFIG  
implicit none
integer      :: elem, elnodes(3), k    
real(kind=8) :: x, y, ax(3), amin

   elnodes=elem2D_nodes(:,elem)
   ax=coord_nod2D(1, elnodes)
   amin=minval(ax)
   DO k=1,3
   if(ax(k)-amin>=domain_length/2.0) ax(k)=ax(k)-domain_length
   if(ax(k)-amin<-domain_length/2.0) ax(k)=ax(k)+domain_length
   END DO
   x=sum(ax)/3.0
   y=sum(coord_nod2D(2,elnodes))/3.0
   
end subroutine elem_center
!
! =======================================================================
!
SUBROUTINE nodal_areas
use o_MESH
use o_ELEMENTS
use o_param
use g_config
use o_elements
use g_rotate_grid
IMPLICIT NONE
! Collects auxilliary information on the mesh
! Allocated and filled in are:
! elem_area(myDim_elem2D)
! area(nl, myDim_nod2D)

integer                                   :: n,j,q, elnodes(3), ed(2), elem, nz
real(kind=8)	                          :: a(2), b(2), ax, ay, lon, lat, vol
real(kind=8), allocatable,dimension(:)   :: work_array
	 
 allocate(area(nl_1, nod2D))   !! Extra size just for simplicity
                             !! in some further routines
 ! =============
 ! Scalar element 
 ! areas at different levels (there can be partly land)
 ! =============
 
 area=0.0
 DO n=1, nod2D
    DO j=1,nod_in_elem2D(n)%nmb
       elem=nod_in_elem2D(n)%addresses(j)
       DO nz=1, elvls(elem)-1
          area(nz,n)=area(nz,n)+voltriangle(elem)/3.0
       END DO
    END DO
 END DO
END SUBROUTINE nodal_areas
