module iceberg_element
 use MOD_PARTIT
 USE MOD_MESH
 USE MOD_DYN
 USE MOD_PARSUP
 use iceberg_params
 use iceberg_thermodynamics
 use iceberg_ocean_coupling
! use iceberg_dynamics
! use iceberg_step

 implicit none

 public ::   mean_gradient
 public ::   nodal_average
 public ::   FEM_eval
 public ::   FEM_3eval
 public ::   iceberg_elem4all
 public ::   find_new_iceberg_elem
 public ::   global2local
 public ::   com_integer
 public ::   matrix_inverse_2x2

 contains

subroutine mean_gradient(mesh, partit, dynamics, elem, lon_rad, lat_rad, nablaeta)
 
 integer, intent(IN) :: elem
 real, 	  intent(IN) :: lon_rad, lat_rad
 real, dimension(2), intent(OUT) :: nablaeta

!LA 2023-03-07
real(kind=WP), dimension(:), pointer    :: eta_n_ib
 integer :: m, node
 real, dimension(3) :: gradientx, gradienty
 logical :: notmynode

type(t_mesh), intent(in) , target :: mesh
type(t_partit), intent(inout), target :: partit
type(t_dyn)   , intent(inout), target :: dynamics
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"
 
 nablaeta = 0.
 gradientx= 0.
 gradienty= 0.
 notmynode= .false.
 
 do m = 1, 3
  node = elem2D_nodes(m,elem)
  call nodal_average(mesh, partit, dynamics, node, gradientx(m), gradienty(m), notmynode) 
  if (notmynode) exit
 end do
 
 if (notmynode) then 	!do no smoothing on gradient this time; change this
 			!so its not dependent on the processor distribution!
  
  !nablaeta(1) = sum( ssh(elem2D_nodes(:,elem)) * bafux_2D(:,elem) )
  !nablaeta(2) = sum( ssh(elem2D_nodes(:,elem)) * bafuy_2D(:,elem) )

!LA 2023-03-07
eta_n_ib => dynamics%eta_n_ib(:)
! kh 18.03.21 use eta_n_ib buffered values here
  nablaeta(1) = sum( eta_n_ib(elem2D_nodes(:,elem)) * gradient_sca(1:3, elem)) 
  nablaeta(2) = sum( eta_n_ib(elem2D_nodes(:,elem)) * gradient_sca(4:6, elem)) 
 else
  call FEM_3eval(mesh, partit, nablaeta(1),nablaeta(2),lon_rad,lat_rad,gradientx,gradienty,elem)
 end if
  
end subroutine mean_gradient


 !***************************************************************************************************************************
 !***************************************************************************************************************************

subroutine nodal_average(mesh, partit, dynamics, local_idx, gradientx, gradienty, notmynode)  
 use MOD_PARTIT
 USE MOD_MESH
 USE MOD_DYN

 implicit none
 
 integer, intent(IN) :: local_idx
 real,    intent(OUT):: gradientx, gradienty
 logical, intent(OUT):: notmynode
 
 integer :: k, node, idx_elem, elem
 real :: area_, patch

!LA 2023-03-07
real(kind=WP), dimension(:), pointer    :: eta_n_ib

type(t_mesh), intent(in) , target :: mesh
type(t_partit), intent(inout), target :: partit
type(t_dyn)   , intent(inout), target :: dynamics
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"
 
 area_ = 0.
 patch = 0.
 gradientx = 0.
 gradienty = 0.
 notmynode = .false.
 
 if(local_idx > myDim_nod2D) then
  notmynode = .true.
  return !not my node
 end if 

 do node=1,myDim_nod2D
   do k=1, nod_in_elem2D_num(node) 
     elem  = nod_in_elem2D(k,node) 

     !do idx_elem = 1,  nod_in_elem2D(local_idx)%nmb
     !elem = nod_in_elem2D(local_idx)%addresses(idx_elem)
     !area_ = voltriangle(elem)
     area_ = elem_area(elem)
     patch= patch + area_
     
     !gradientx = gradientx + area * sum( ssh(elem2D_nodes(:,elem)) * bafux_2D(:,elem) )
     !gradienty = gradienty + area * sum( ssh(elem2D_nodes(:,elem)) * bafuy_2D(:,elem) )

!LA 2023-03-07
eta_n_ib => dynamics%eta_n_ib(:)
! kh 18.03.21 use eta_n_ib buffered values here
     gradientx = gradientx + area_ * sum( eta_n_ib(elem2D_nodes(:,elem)) * gradient_sca(1:3, elem)) 
     gradienty = gradienty + area_ * sum( eta_n_ib(elem2D_nodes(:,elem)) * gradient_sca(4:6, elem)) 
   end do
 end do
 
 gradientx = gradientx / patch
 gradienty = gradienty / patch
 
end subroutine nodal_average


 !***************************************************************************************************************************
 !***************************************************************************************************************************

!=================================================================
! evaluates a given FEM vectorfield at location lon, lat (radiant)
! OUT: u_at_ib, v_at_ib (u and v component at location of iceberg)
! IN : lon, lat (position of iceberg in radiant)
!      field_u, field_v (u and v components of vector field)
!      elem	(the LOCAL element the iceberg lies in)
!=================================================================
subroutine FEM_eval(mesh, partit, u_at_ib,v_at_ib,lon,lat,field_u,field_v,elem)
  use MOD_PARTIT !for myDim_nod2D, eDim_nod2D
  use o_param !for rad
  USE MOD_MESH

  implicit none
  real, intent(in)			:: lon, lat 
  integer, intent(in)			:: elem
  real, intent(out) 			:: u_at_ib
  real, intent(out) 			:: v_at_ib
  real, DIMENSION(3)			:: phi, values_u, values_v
  real					:: lon_deg, lat_deg
  real, dimension(2)             :: coords_tmp

type(t_mesh), intent(in) , target :: mesh
type(t_partit), intent(inout), target :: partit
  real, dimension(partit%myDim_nod2D+partit%eDim_nod2D), intent(in)  	:: field_u
  real, dimension(partit%myDim_nod2D+partit%eDim_nod2D), intent(in)  	:: field_v
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"
  
  !convert to deg
  lon_deg = lon/rad
  lat_deg = lat/rad
  coords_tmp = [lon_deg, lat_deg]
  
  !values of the 3 local basisfunctions at the 
  !position 'coords'
  call locbafu_2D(mesh, partit, phi,elem,coords_tmp)

  values_u = field_u(elem2D_nodes(:,elem))
  values_v = field_v(elem2D_nodes(:,elem))

  u_at_ib = sum( values_u(:) * phi(:))
  v_at_ib = sum( values_v(:) * phi(:))

  !correct small errors
  if (u_at_ib < minval(values_u, 1)) u_at_ib=minval(values_u, 1)
  if (v_at_ib < minval(values_v, 1)) v_at_ib=minval(values_v, 1)
  if (u_at_ib > maxval(values_u, 1)) u_at_ib=maxval(values_u, 1)
  if (v_at_ib > maxval(values_v, 1)) v_at_ib=maxval(values_v, 1)

  !from E6F:
  !interp = sum( values(:) * phi(:))
  !if (interp < minval(values, 1)) interp=minval(values, 1)
  !if (interp > maxval(values, 1)) interp=maxval(values, 1)

end subroutine FEM_eval


 !***************************************************************************************************************************
 !***************************************************************************************************************************

!=================================================================
! evaluates a given FEM vectorfield at location lon, lat (radiant)
! OUT: u_at_ib, v_at_ib (u and v component at location of iceberg)
! IN : lon, lat (position of iceberg in radiant)
!      field_u, field_v (u and v components of vector field)
!      elem	(the LOCAL element the iceberg lies in)
!=================================================================
subroutine FEM_eval_old(mesh, partit, u_at_ib,v_at_ib,lon,lat,field_u,field_v,elem)
  use o_param
  use g_clock
  use g_forcing_arrays
  use g_rotate_grid
  
  !use iceberg module
  !use iceberg_params
  
  implicit none
  real, intent(in)			:: lon, lat 
  integer				:: elem
  real, intent(out) 			:: u_at_ib
  real, intent(out) 			:: v_at_ib
  
  real					:: x, y, x1,x2,x3,y1,y2,y3
  real					:: T1_u, T1_v, T2_u, T2_v, T3_u, T3_v
  real, dimension(2,2) 			:: inv_matrix
  real, dimension(2)			:: alphabeta
  real					:: maxlon, minlon, maxlat, minlat

type(t_mesh), intent(in) , target :: mesh
type(t_partit), intent(inout), target :: partit
  real, dimension(partit%myDim_nod2D+partit%eDim_nod2D), intent(in)  	:: field_u
  real, dimension(partit%myDim_nod2D+partit%eDim_nod2D), intent(in)  	:: field_v
!type(t_ice),    intent(inout), target :: ice
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"
  
  !location of iceberg
  x = lon
  y = lat
   
  
  !coords of the 3 nodes
  x1 = coord_nod2D(1,elem2D_nodes(1,elem))
  y1 = coord_nod2D(2,elem2D_nodes(1,elem))
  x2 = coord_nod2D(1,elem2D_nodes(2,elem))
  y2 = coord_nod2D(2,elem2D_nodes(2,elem))  
  x3 = coord_nod2D(1,elem2D_nodes(3,elem))
  y3 = coord_nod2D(2,elem2D_nodes(3,elem))
  
  
  !check whether iceberg position is in the element
  maxlon = maxval( (/ x1, x2, x3 /) )
  minlon = minval( (/ x1, x2, x3 /) )
  maxlat = maxval( (/ y1, y2, y3 /) )
  minlat = minval( (/ y1, y2, y3 /) )
  if( (x > maxlon) .OR. (x < minlon) ) then
    write(*,*) 'FEM_eval error: iceberg lon ', x, ' outside element!'
    write(*,*) 'maxlon:', maxlon, ' minlon:', minlon
    call par_ex (partit%MPI_COMM_FESOM, partit%mype)
    stop
  else if( (y > maxlat) .OR. (y < minlat)) then
    write(*,*) 'FEM_eval error: iceberg lat', y, ' outside element!'
    write(*,*) 'maxlat:', maxlat, ' minlat:', minlat
    call par_ex (partit%MPI_COMM_FESOM, partit%mype)
    stop
  else
    !everything okay
  end if 
  
  
  !distances wrt node 1
  x2 = x2 - x1
  y2 = y2 - y1 
  x3 = x3 - x1
  y3 = y3 - y1
  x  = x - x1
  y  = y - y1 
  
  !nodal values Ti_u, Ti_v
  T1_u = field_u(elem2D_nodes(1,elem))
  T1_v = field_v(elem2D_nodes(1,elem))
  T2_u = field_u(elem2D_nodes(2,elem))
  T2_v = field_v(elem2D_nodes(2,elem))  
  T3_u = field_u(elem2D_nodes(3,elem))
  T3_v = field_v(elem2D_nodes(3,elem))
  
  !differences wrt node 1
  T2_u = T2_u - T1_u
  T2_v = T2_v - T1_v  
  T3_u = T3_u - T1_u
  T3_v = T3_v - T1_v
  
  !determine alpha and beta for u velocity
  inv_matrix(1,1) = y2
  inv_matrix(1,2) = -y3
  inv_matrix(2,1) = -x2
  inv_matrix(2,2) = x3
  inv_matrix = (1./(x3*y2 - x2*y3)) * inv_matrix
  alphabeta = MATMUL(inv_matrix, (/ T3_u, T2_u /))
  
  u_at_ib = T1_u + alphabeta(1)*x + alphabeta(2)*y
  
  !same for v velocity
  alphabeta = MATMUL(inv_matrix, (/ T3_v, T2_v /))
  
  v_at_ib = T1_v + alphabeta(1)*x + alphabeta(2)*y 
end subroutine FEM_eval_old


 !***************************************************************************************************************************
 !***************************************************************************************************************************

!=================================================================
! interpolates ocean velocity to ib's location (lon, lat) [radiant]
! OUT: u_at_ib, v_at_ib (u and v component at location of iceberg)
! IN : lon, lat (position of iceberg in radiant)
!      ocean_u, ocean_v (3 nodal values for u and v component where
!                        ocean_u(m) is value of 
!                        elem2D_nodes(m,iceberg_elem), m=1,2,3 )
!      elem	(the LOCAL element the iceberg lies in)
!=================================================================
subroutine FEM_3eval(mesh, partit, u_at_ib,v_at_ib,lon,lat,ocean_u,ocean_v,elem)
  use MOD_PARTIT !for myDim_nod2D, eDim_nod2D
  use o_param !for rad
  use MOD_MESH

  implicit none
  real, intent(in)			:: lon, lat 
  real, dimension(3), intent(in)  	:: ocean_u
  real, dimension(3), intent(in)  	:: ocean_v
  integer, intent(in)			:: elem
  real, intent(out) 			:: u_at_ib
  real, intent(out) 			:: v_at_ib
  real, DIMENSION(3)			:: phi
  real					:: lon_deg, lat_deg
  real, dimension(2)             :: coords_tmp
  
type(t_mesh), intent(in) , target :: mesh
type(t_partit), intent(inout), target :: partit
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"
  
  !convert to deg
  lon_deg = lon/rad
  lat_deg = lat/rad
  
  !values of the 3 local basisfunctions at the 
  !position 'coords'
  coords_tmp = [lon_deg, lat_deg]
  call locbafu_2D(mesh, partit, phi,elem, coords_tmp)
  
  u_at_ib = sum( ocean_u(:) * phi(:))
  v_at_ib = sum( ocean_v(:) * phi(:))

  !correct small errors
  if (u_at_ib < minval(ocean_u, 1)) u_at_ib=minval(ocean_u, 1)
  if (v_at_ib < minval(ocean_v, 1)) v_at_ib=minval(ocean_v, 1)
  if (u_at_ib > maxval(ocean_u, 1)) u_at_ib=maxval(ocean_u, 1)
  if (v_at_ib > maxval(ocean_v, 1)) v_at_ib=maxval(ocean_v, 1)

  !from E6F:
  !interp = sum( values(:) * phi(:))
  !if (interp < minval(values, 1)) interp=minval(values, 1)
  !if (interp > maxval(values, 1)) interp=maxval(values, 1)
end subroutine FEM_3eval

 !***************************************************************************************************************************
 !***************************************************************************************************************************

subroutine iceberg_elem4all(mesh, partit, elem, lon_deg, lat_deg) 
 USE MOD_MESH
 use MOD_PARTIT		!for myDim_nod2D, myList_elem2D
 
 implicit none
 
 integer, intent(INOUT) :: elem
 real, intent(IN) :: lon_deg, lat_deg
 logical :: i_have_element, reject_tmp

type(t_mesh), intent(in) , target :: mesh
type(t_partit), intent(inout), target :: partit
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"
 
  call point_in_triangle(mesh, partit, elem, (/lon_deg, lat_deg/)) !all PEs search here
  i_have_element= (elem .ne. 0) !up to 3 PEs .true.
  
  if(i_have_element) then
   i_have_element= elem2D_nodes(1,elem) <= myDim_nod2D !1 PE still .true.
   if (use_cavity) then 
      !reject_tmp = all( (mesh%cavity_depth(elem2D_nodes(:,elem))/=0.0) .OR. (mesh%bc_index_nod2D(elem2D_nodes(:,elem))==0.0) )
      reject_tmp = any(mesh%cavity_depth(elem2D_nodes(:,elem))/=0.0) .OR. all(mesh%bc_index_nod2D(elem2D_nodes(:,elem))==0.0)
      if(reject_tmp) then
      !if( reject_elem(mesh, partit, elem) ) then
       elem=0 !reject element
       i_have_element=.false.
       write(*,*) 'elem4all: iceberg found in shelf region: elem = 0'
      else 
       elem=myList_elem2D(elem) !global now
      end if 
   else
    elem=myList_elem2D(elem) !global now
   endif 
  end if
  call com_integer(partit, i_have_element,elem) !max 1 PE sends element here; 
end subroutine iceberg_elem4all


 !***************************************************************************************************************************

subroutine find_new_iceberg_elem(mesh, partit, old_iceberg_elem, pt, left_mype)
  use o_param

  implicit none
  
  logical                :: reject_tmp
  INTEGER, INTENT(INOUT) :: old_iceberg_elem
  REAL, DIMENSION(2), INTENT(IN) :: pt
  real, INTENT(OUT) :: left_mype
  
  INTEGER :: m, n2, idx_elem_containing_n2, elem_containing_n2, ibelem_tmp
  REAL, DIMENSION(3) :: werte2D

type(t_mesh), intent(in) , target :: mesh
type(t_partit), intent(inout), target :: partit
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"
  

ibelem_tmp=old_iceberg_elem
left_mype=0.0

!for each node of the old iceberg element...
do m=1, 3
 n2=elem2D_nodes(m,old_iceberg_elem)
 
 if(n2 > myDim_nod2D) cycle !n2 is not my node, so i cannot access all elements around it
 
 !...and for each element containing this node (so we get all neighbour elements)...
 !do idx_elem_containing_n2 = 1,  nod_in_elem2D(n2)%nmb
 do idx_elem_containing_n2 = 1,  nod_in_elem2D_num(n2)
 
  !elem_containing_n2 = nod_in_elem2D(n2)%addresses(idx_elem_containing_n2)
  elem_containing_n2 = nod_in_elem2D(idx_elem_containing_n2,n2) 
    
  call locbafu_2D(mesh, partit, werte2D, elem_containing_n2, pt)
   
  if (ALL(werte2D <= 1.+ 1.0e-07) .AND. ALL(werte2D >= 0.0- 1.0e-07) ) then
   old_iceberg_elem=elem_containing_n2
   if (use_cavity) then 
      !if( reject_elem(mesh, partit, old_iceberg_elem) ) then
      !reject_tmp = all( (mesh%cavity_depth(elem2D_nodes(:,ibelem_tmp))/=0.0) .OR. (mesh%bc_index_nod2D(elem2D_nodes(:,ibelem_tmp))==0.0) )
      reject_tmp = any(mesh%cavity_depth(elem2D_nodes(:,ibelem_tmp))/=0.0) .OR. all(mesh%bc_index_nod2D(elem2D_nodes(:,ibelem_tmp))==0.0)
      if(reject_tmp) then
         left_mype=1.0
         write(*,*) 'iceberg found in shelf region: left_mype = 1'
         old_iceberg_elem=ibelem_tmp
      end if
   endif
   RETURN 
  end if
 end do
end do  

!no element found (including old element!)
!or the iceberg is way too fast..
left_mype=1.0
end subroutine find_new_iceberg_elem


 !***************************************************************************************************************************
 !***************************************************************************************************************************

SUBROUTINE point_in_triangle(mesh, partit, el2D, pt)
  !  returns triangle containing the point pt
  !  lon, lat in deg
  
  USE o_param
  !use o_mesh
  USE MOD_MESH
  USE MOD_PARTIT !for myDim_elem2D
  
  IMPLICIT NONE
  
  REAL, DIMENSION(2), INTENT(IN) :: pt
  INTEGER, INTENT(OUT) :: el2D
  
  INTEGER :: i, k, l
  REAL, DIMENSION(3) :: werte2D
  REAL               :: xdiff, ydiff

type(t_mesh), intent(in) , target :: mesh
type(t_partit), intent(inout), target :: partit
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"

  el2D=0
  !DO l=1,elem2D
  DO l=1,partit%myDim_elem2D
 
     call locbafu_2D(mesh, partit, werte2D, l, pt)
     
     if (ALL(werte2D <= 1.+ 1.0e-07) .AND. ALL(werte2D >= 0.0- 1.0e-07) ) then
        el2D=l
        !print *,'"point_in_triangle": Das 2D-Element ist ',l, werte2D
        EXIT
     end if
  END DO
END SUBROUTINE point_in_triangle


 !***************************************************************************************************************************
 !***************************************************************************************************************************

!coords in deg, elem is LOCAL element; 
!values of the 3 local basisfunctions at the 
!position 'coords' are returned
SUBROUTINE locbafu_2D(mesh, partit, values, elem, coords)
  !use o_mesh
  USE MOD_MESH
  use MOD_PARTIT
  USE o_param

  IMPLICIT NONE
  
  REAL, DIMENSION(3), INTENT(OUT) :: values
  INTEGER, INTENT(IN) :: elem
  REAL, DIMENSION(2), INTENT(IN) :: coords
  
  INTEGER :: i,j
  INTEGER :: node
  INTEGER, DIMENSION(3) :: local_nodes
  REAL, DIMENSION(2,3) :: local_coords, local_cart_coords
  REAL, DIMENSION(2,2) :: TRANS, TRANS_inv
  REAL :: DET
  REAL, DIMENSION(2) :: x, x_cart, stdel_coords
  REAL, DIMENSION(2) :: vec

type(t_mesh), intent(in) , target :: mesh
type(t_partit), intent(inout), target :: partit
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"
 
  x(1)=coords(1)*rad
  x(2)=coords(2)*rad
  do i=1,3
     node=mesh%elem2D_nodes(i,elem)
     local_nodes(i)=node
     local_coords(:,i)=mesh%coord_nod2D(:,node)
  end do

  DO i=1, 2
     if (local_coords(1,i+1)-local_coords(1,1) > 180.*rad) then 
         local_coords(1,i+1)=local_coords(1,i+1)-360.*rad
     end if
     if (local_coords(1,i+1)-local_coords(1,1) <-180.*rad) then 
         local_coords(1,i+1)=local_coords(1,i+1)+360.*rad
     end if
  END DO  

  if (x(1)-local_coords(1,1) > 180.*rad) then 
      x(1)=x(1)-360.*rad
  end if
  if (x(1)-local_coords(1,1) <-180.*rad) then 
      x(1)=x(1)+360.*rad
  end if
  !  cartesian coordinates
  x_cart(1) = r_earth * COS(x(2)) * x(1)
  x_cart(2) = r_earth * x(2)  

  do i=1,3
     local_cart_coords(1,i) = r_earth * COS(local_coords(2,i)) * local_coords(1,i)
     local_cart_coords(2,i) = r_earth * local_coords(2,i)
  end do  
  DO i=1, 2
     TRANS(:,i) = local_cart_coords(:,i+1)-local_cart_coords(:,1)     
  END DO  
  call matrix_inverse_2x2(TRANS, TRANS_inv, DET, elem, coords)  

  vec=x_cart-local_cart_coords(:,1)  
  stdel_coords = MATMUL(TRANS_inv, vec)
  call stdbafu_2D(values,  stdel_coords,1)
END SUBROUTINE locbafu_2D


 !**************************************************************************************************************************
 !***************************************************************************************************************************

SUBROUTINE stdbafu_2D(values,x,m) !(stdbafu,x,m)
  
  IMPLICIT NONE
  
  INTEGER, INTENT(IN) :: m
  REAL, DIMENSION(3,m), INTENT(OUT) :: values
  REAL, DIMENSION(2,m), INTENT(IN) :: x  
  INTEGER :: k

  do k=1,m
     values(1,k)=1.-x(1,k)-x(2,k)
     values(2,k)=   x(1,k)
     values(3,k)=          x(2,k)
  end do
END SUBROUTINE stdbafu_2D


 !***************************************************************************************************************************
 !***************************************************************************************************************************

subroutine global2local(mesh, partit, aux, tmp)
 use MOD_PARTIT !for myDim_elem2D, myList_elem2D
 !use o_mesh
 USE MOD_MESH
 implicit none
 
 integer, intent(in):: tmp
 integer, dimension(tmp), intent(out):: aux
 integer:: n

type(t_mesh), intent(in) , target :: mesh
type(t_partit), intent(inout), target :: partit
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"

 aux = 0
 do n = 1, myDim_elem2D
  aux(myList_elem2D(n)) = n
 end do
end subroutine global2local


 !***************************************************************************************************************************
 !***************************************************************************************************************************

subroutine com_integer(partit, i_have_element, iceberg_element)
 use MOD_PARTIT !for npes
 implicit none
 
 logical, intent(in):: i_have_element
 integer, intent(inout):: iceberg_element
 
 integer:: status(MPI_STATUS_SIZE)
 integer:: req
 logical:: completed
type(t_partit), intent(inout), target :: partit
!#include "associate_part_def.h"
!#include "associate_part_ass.h"

! kh 10.02.22
 if(i_have_element) then
!$omp critical
     call MPI_IAllreduce(MPI_IN_PLACE, iceberg_element, 1, MPI_INTEGER, MPI_SUM, partit%MPI_COMM_FESOM_IB, req, partit%MPIERR_IB)
!$omp end critical
 else
!$omp critical
     call MPI_IAllreduce(0,            iceberg_element, 1, MPI_INTEGER, MPI_SUM, partit%MPI_COMM_FESOM_IB, req, partit%MPIERR_IB)
!$omp end critical
 end if

 completed = .false.
 do while (.not. completed)
!$omp critical
     CALL MPI_TEST(req, completed, status, partit%MPIERR_IB)
!$omp end critical
 end do

 end subroutine com_integer



 !***************************************************************************************************************************
 !***************************************************************************************************************************
 
! ! kh 10.02.21
! subroutine com_values_old_dont_use(partit, i_have_element, arr, iceberg_element)
! use MOD_PARTIT !for npes
! implicit none
! 
! logical, intent(in)   	:: i_have_element
! real,    intent(inout)	:: arr(15)
! integer, intent(inout)	:: iceberg_element
! 
! logical:: he_has_element
! real   :: arr_r(15)
! integer:: i, sender, status(MPI_STATUS_SIZE)
!!type(t_partit), intent(inout), target :: partit
!!#include "associate_part_def.h"
!!#include "associate_part_ass.h"
!
! if (mype==0) then
!    do i=1, npes-1
!       CALL MPI_RECV(he_has_element, 1, MPI_LOGICAL, MPI_ANY_SOURCE, 0, MPI_COMM_FESOM, status, MPIerr )		      
!       sender = status(MPI_SOURCE)
!       if (he_has_element) then
!          CALL MPI_RECV(arr_r, 15, MPI_DOUBLE_PRECISION, sender, 1, MPI_COMM_FESOM, status, MPIerr )
!	  CALL MPI_RECV(iceberg_element, 1, MPI_DOUBLE_PRECISION, sender, 2, MPI_COMM_FESOM, status, MPIerr )
!	  arr=arr_r
!       end if
!    end do
! else
!       CALL MPI_SEND(i_have_element, 1, MPI_LOGICAL, 0, 0, MPI_COMM_FESOM, MPIerr )
!       if (i_have_element) then
!          CALL MPI_SEND(arr, 15, MPI_DOUBLE_PRECISION,0, 1, MPI_COMM_FESOM, MPIerr )
!	  CALL MPI_SEND(iceberg_element, 1, MPI_INTEGER,0, 2, MPI_COMM_FESOM, MPIerr )
!       end if 
! end if
! 
!! if (mype==0) then
!!    do i=1, npes-1
!!       CALL MPI_RECV(he_has_element, 1, MPI_LOGICAL, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, status, MPIerr )		      
!!       sender = status(MPI_SOURCE)
!!       if (he_has_element) then
!!          CALL MPI_RECV(arr_r, 15, MPI_DOUBLE_PRECISION, sender, 1, MPI_COMM_WORLD, status, MPIerr )
!!	  CALL MPI_RECV(iceberg_element, 1, MPI_DOUBLE_PRECISION, sender, 2, MPI_COMM_WORLD, status, MPIerr )
!!	  arr=arr_r
!!       end if
!!    end do
!! else
!!       CALL MPI_SEND(i_have_element, 1, MPI_LOGICAL, 0, 0, MPI_COMM_WORLD, MPIerr )
!!       if (i_have_element) then
!!          CALL MPI_SEND(arr, 15, MPI_DOUBLE_PRECISION,0, 1, MPI_COMM_WORLD, MPIerr )
!!	  CALL MPI_SEND(iceberg_element, 1, MPI_INTEGER,0, 2, MPI_COMM_WORLD, MPIerr )
!!       end if 
!! end if
! 
! ! *** PROC 0 SENDS ICEBERG ELEMENT TO ALL OTHERS ***
! !
! !1. buffer - Startadresse des Datenpuffers
! !2. count - Anzahl der Elemente im Puffer (integer)
! !3. datatype - Datentyp der Pufferelemente (handle)
! !4. root - Wurzelproze�; der, welcher sendet (integer)
! !5. comm - Kommunikator (handle) 
! CALL MPI_BCAST(arr, 15, MPI_DOUBLE_PRECISION, 0, MPI_COMM_FESOM, MPIerr)
! CALL MPI_BCAST(iceberg_element, 1, MPI_INTEGER, 0, MPI_COMM_FESOM, MPIerr)
! !CALL MPI_BCAST(arr, 15, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, MPIerr)
! !CALL MPI_BCAST(iceberg_element, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, MPIerr)
! 
! ! kh 10.02.21
! end subroutine com_values_old_dont_use


 !***************************************************************************************************************************
 !***************************************************************************************************************************
 
! !==============================================================================
!! routine for visualizing the distribution of the processors involved
!!
!!   Thomas Rackow, 30.07.2010
!!==============================================================================
!SUBROUTINE processor_distr(mesh, partit)
!  !use o_mesh
!  USE MOD_MESH
!  use o_param
!!  use i_therm_param
!!  use i_param
!!  use i_arrays
!  use MOD_PARTIT
!
!! kh 18.03.21 not really used here
!! use o_arrays
!
!  use g_clock
!  use g_forcing_arrays
!  use g_rotate_grid
!
!  IMPLICIT NONE
!  
!  character			:: mype_char*3
!  INTEGER			:: m, row
!
!type(t_mesh), intent(in) , target :: mesh
!type(t_partit), intent(inout), target :: partit
!#include "associate_part_def.h"
!#include "associate_mesh_def.h"
!#include "associate_part_ass.h"
!#include "associate_mesh_ass.h"
!  
!  write(mype_char,*) mype
!  !left-adjust the string..
!  mype_char = adjustl(mype_char)
!
!  DO m=1, myDim_nod2d
!  row=myList_nod2d(m)
!  open(unit=mype+66, file='/work/ab0046/a270046/results/ICB02processor' // trim(mype_char) // '.dat', position='append')
!  !		global local    lon		   lat		PE	
!  write(mype+66,*) row, m, coord_nod2D(1,m), coord_nod2D(2,m), mype
!  close(mype+66)
!  END DO
!END SUBROUTINE processor_distr 


 !***************************************************************************************************************************
 !***************************************************************************************************************************

!!==============================================================================
!! routine for visualizing the amplitude of the main 4 tidal constituents
!!
!!   Thomas Rackow, 18.12.2010
!!==============================================================================
!SUBROUTINE eides_distr(partit)
!!  use o_mesh
!  use o_param
! ! use i_therm_param
! ! use i_param
! ! use i_arrays
!  use MOD_PARTIT
!! kh 18.03.21 not really used here
!! use o_arrays
!
!  use g_clock
!  use g_forcing_arrays
!  use g_rotate_grid
!  
!  IMPLICIT NONE
!  
!  character			:: mype_char*3
!  INTEGER			:: m, row
!type(t_partit), intent(inout), target :: partit
!#include "associate_part_def.h"
!#include "associate_part_ass.h"
!  
!  write(mype_char,*) mype
!  !left-adjust the string..
!  mype_char = adjustl(mype_char)
!
!  !DO m=1, myDim_nod2d
!  !row=myList_nod2d(m)
!  !open(unit=mype+66, file='/work/ab0046/a270046/results/TIDESprocessor' // trim(mype_char) // '.dat', position='append')
!  !!		global local    M2		  S2		  K1		   O1	
!  !write(mype+66,*) row, m, tide_z_amp(m,1), tide_z_amp(m,2), tide_z_amp(m,3), tide_z_amp(m,4)
!  !close(mype+66)
!  !END DO
!END SUBROUTINE tides_distr

!LA from oce_mesh_setup ofr iceberg coupling
subroutine  matrix_inverse_2x2 (A, AINV, DET, elem, coords)
  !
  ! Coded by Sergey Danilov
  ! Reviewed by Qiang Wang
  !-------------------------------------------------------------
  
  implicit none
 
  integer                                   :: elem
  REAL, DIMENSION(2), INTENT(IN) :: coords
  
  real(kind=8), dimension(2,2), intent(IN)  :: A
  real(kind=8), dimension(2,2), intent(OUT) :: AINV
  real(kind=8), intent(OUT)                 :: DET
  integer                                   :: i,j
  
  DET  = A(1,1)*A(2,2) - A(1,2)*A(2,1)
  if ( DET .eq. 0.0 )  then
     do j=1,2
        write(*,*) (A(i,j),i=1,2)
     end do
     write(*,*) " * elem: ", elem, ", coords: ", coords
     stop 'SINGULAR 2X2 MATRIX'
  else
     AINV(1,1) =  A(2,2)/DET
     AINV(1,2) = -A(1,2)/DET
     AINV(2,1) = -A(2,1)/DET
     AINV(2,2) =  A(1,1)/DET
  endif
end subroutine matrix_inverse_2x2
end module iceberg_element
