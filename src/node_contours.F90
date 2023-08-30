subroutine node_contours(my_x_corners, my_y_corners, partit, mesh)
    USE MOD_MESH
    USE MOD_PARTIT
    USE MOD_PARSUP
    USE o_PARAM
    use g_comm_auto
    use o_ARRAYS
    use g_rotate_grid

    IMPLICIT NONE
    type(t_mesh),   intent(inout), target :: mesh
    type(t_partit), intent(inout), target :: partit
    real(kind=WP), allocatable, intent(inout) :: my_x_corners(:,:)     ! longitude node corners
    real(kind=WP), allocatable, intent(inout) :: my_y_corners(:,:)     ! latitude node corners    
    integer                               :: bEdge_left, bEdge_right
    integer,              dimension(2)    :: belem_left, belem_right
    integer                               :: edge_left, edge_right
    integer                               :: n, ee, elem, nn, el(2), flag, nn1, nn2
    integer, allocatable, dimension(:)    :: nedges, nelems, nedges1, nelems1, nedges2, nelems2
    real(kind=WP)              :: this_x_coord, this_y_coord
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"

ALLOCATE(my_x_corners(myDim_nod2D, 25)) !maxval(nod_in_elem2D_num, 1)*2+2))
ALLOCATE(my_y_corners(myDim_nod2D, 25)) !maxval(nod_in_elem2D_num, 1)*2+2))

do n=1, myDim_nod2D
    ! find the type of node: internal or at boundary
    bEdge_left =0
    belem_left =0
    bEdge_right=0
    belem_right=0

    do ee=1, nod_in_elem2D_num(n)
       elem=nod_in_elem2D(ee,n)
       if (elem2D_nodes(1,elem)==n) then
          edge_left=elem_edges(3,elem)
          edge_right=elem_edges(2,elem)
       elseif (elem2D_nodes(2,elem)==n) then
          edge_left=elem_edges(1,elem)
          edge_right=elem_edges(3,elem)
       else
          edge_left=elem_edges(2,elem)
          edge_right=elem_edges(1,elem)
       end if
       if (myList_edge2D(edge_left)>edge2D_in) then
          bEdge_left=bEdge_left+1
          belem_left(bEdge_left)=elem
       end if
       if (myList_edge2D(edge_right)>edge2D_in) then
          bEdge_right=bEdge_right+1
          belem_right(bEdge_right)=elem
       end if
    end do

! now we have three cases
   if (bEdge_left==0) then      ! inner contour
      elem=nod_in_elem2D(1, n)  ! we can start from any
      allocate(nedges(nod_in_elem2D_num(n)))
      nedges=0
      allocate(nelems(nod_in_elem2D_num(n)))
      nelems=0
      !!!!!!! inner_node_contour
#include "node_contour_inner.h"
      do nn=1, nod_in_elem2D_num(n)
         call edge_center(edges(1, nedges(nn)), edges(2, nedges(nn)), my_x_corners(n, (nn-1)*2+1), my_y_corners(n, (nn-1)*2+1), mesh)
         call elem_center(nelems(nn), my_x_corners(n, (nn-1)*2+2), my_y_corners(n, (nn-1)*2+2), mesh)
      end do
      do nn=nod_in_elem2D_num(n)+1, size(my_x_corners, 2)
         my_x_corners(n, nn)=my_x_corners(n, nod_in_elem2D_num(n))
         my_y_corners(n, nn)=my_y_corners(n, nod_in_elem2D_num(n))
      end do
      deallocate(nedges, nelems)
   end if

   if (bEdge_left==1) then ! standard boundary node
      elem=belem_left(1)
      allocate(nedges(nod_in_elem2D_num(n)+1))
      nedges=0
      allocate(nelems(nod_in_elem2D_num(n)))
      nelems=0
      !!!!!!!boundary_node_contour
#include "node_contour_boundary.h"
      do nn=1, nod_in_elem2D_num(n)
         call edge_center(edges(1, nedges(nn)), edges(2, nedges(nn)), my_x_corners(n, (nn-1)*2+1), my_y_corners(n, (nn-1)*2+1), mesh)
         call elem_center(nelems(nn), my_x_corners(n, (nn-1)*2+2), my_y_corners(n, (nn-1)*2+2), mesh)
      end do
      nn=nod_in_elem2D_num(n)+1
      call edge_center(edges(1, nedges(nn)), edges(2, nedges(nn)), my_x_corners(n, (nn-1)*2+1), my_y_corners(n, (nn-1)*2+1), mesh)
      nn=nod_in_elem2D_num(n)+2
      my_x_corners(n, nn)=coord_nod2D(1,n)
      my_y_corners(n, nn)=coord_nod2D(2,n)
      do nn=nod_in_elem2D_num(n)+3, size(my_x_corners, 2)
         my_x_corners(n, nn)=my_x_corners(n, nod_in_elem2D_num(n)+2)
         my_y_corners(n, nn)=my_y_corners(n, nod_in_elem2D_num(n)+2)
      end do
      !!!!!!!
      deallocate(nedges, nelems)
   end if

   if (bEdge_left==2) then  ! strange boundary node
       elem=belem_left(1)
       allocate(nedges (nod_in_elem2D_num(n)+1))
       allocate(nedges1(nod_in_elem2D_num(n)+1))
       nedges =0
       nedges1=0
       allocate(nelems (nod_in_elem2D_num(n)))
       allocate(nelems1(nod_in_elem2D_num(n)))
       nelems=0
       nelems1=0
       !!!!!!!boundary_node_contour
#include "node_contour_boundary.h"
       where (nedges>0)
             nedges1=nedges
       end where
       where (nelems>0)
             nelems1=nelems
       end where
       nn1=nn
       do nn=1, nn1
          call edge_center(edges(1, nedges1(nn)), edges(2, nedges1(nn)), my_x_corners(n, (nn-1)*2+1), my_y_corners(n, (nn-1)*2+1), mesh)
          call elem_center(nelems1(nn), my_x_corners(n, (nn-1)*2+2), my_y_corners(n, (nn-1)*2+2), mesh)
       end do
       nn=nn1+1
       call edge_center(edges(1, nedges1(nn)), edges(2, nedges1(nn)), my_x_corners(n, (nn-1)*2+1), my_y_corners(n, (nn-1)*2+1), mesh)
       nn=nn1+2
       my_x_corners(n, nn)=coord_nod2D(1,n)
       my_y_corners(n, nn)=coord_nod2D(2,n)
       !!!!!!!
       elem=belem_left(2)
       allocate(nedges2(nod_in_elem2D_num(n)+1))
       nedges =0
       nedges2=0
       allocate(nelems2(nod_in_elem2D_num(n)))
       nelems =0
       nelems2=0
       !!!!!!!boundary_node_contour
#include "node_contour_boundary.h"
       where (nedges>0)
            nedges2=nedges
       end where
       where (nelems>0)
             nelems2=nelems
       end where
       nn2=nn
       do nn=nn1+3, nn1+nn2+2
          call edge_center(edges(1, nedges2(nn)), edges(2, nedges2(nn)), my_x_corners(n, (nn-1)*2+1), my_y_corners(n, (nn-1)*2+1), mesh)
          call elem_center(nelems2(nn), my_x_corners(n, (nn-1)*2+2), my_y_corners(n, (nn-1)*2+2), mesh)
       end do
       nn=nn1+nn2+3
       call edge_center(edges(1, nedges2(nn)), edges(2, nedges2(nn)), my_x_corners(n, (nn-1)*2+1), my_y_corners(n, (nn-1)*2+1), mesh)
       nn=nn1+nn2+4
       my_x_corners(n, nn)=coord_nod2D(1,n)
       my_y_corners(n, nn)=coord_nod2D(2,n)
       do nn=nn1+nn2+5, size(my_x_corners, 2)
          my_x_corners(n, nn)=my_x_corners(n, nn1+nn2+4)
          my_y_corners(n, nn)=my_y_corners(n, nn1+nn2+4)
       end do
       !!!!!!!
       deallocate(nedges, nelems, nedges1, nelems1, nedges2, nelems2)
   end if
end do
do n=1, myDim_nod2D
   do nn=1, size(my_x_corners, 2)
      this_x_coord=my_x_corners(n, nn)
      this_y_coord=my_y_corners(n, nn)
      call r2g(my_x_corners(n, nn), my_y_corners(n, nn), this_x_coord, this_y_coord)
   end do
end do
my_x_corners=my_x_corners/rad
my_y_corners=my_y_corners/rad
end subroutine node_contours