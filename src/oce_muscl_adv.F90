module find_up_downwind_triangles_interface
  interface
    subroutine find_up_downwind_triangles(twork, partit, mesh)
      use MOD_MESH
      USE MOD_PARTIT
      USE MOD_PARSUP
      use MOD_TRACER
      type(t_mesh),        intent(in)  ,  target :: mesh
      type(t_partit),      intent(inout), target :: partit
      type(t_tracer_work), intent(inout), target :: twork
    end subroutine find_up_downwind_triangles
  end interface
end module find_up_downwind_triangles_interface

! A set of routines to implement MUSCL-type of advection
! For description, see Abalakin, I., Dervieux, A., Kozubskaya, T., 2002. A
! vertex-centered high-order MUSCL scheme applying to linearized Euler acoustics.
! INRIA, Rapport de recherche 4459.
!
!  The advection routine is solve_tracer_muscl
!  muscl_adv_init initializes several arrays needed for the algorithm
!  
!  The algorithm works with the concept of upwind and downwind triangles to a given
!  edge. This introduces additional halo communications.
!  sergey.danilov@awi.de 2012
!
!Contains:
!	muscl_adv_init
!	find_up_downwind_triangles
!	fill_up_dn_grad
!	adv_tracer_muscl
subroutine muscl_adv_init(twork, partit, mesh)
    use MOD_MESH
    USE MOD_PARTIT
    USE MOD_PARSUP
    use MOD_TRACER
    use o_ARRAYS
    use o_PARAM
    use g_comm_auto
    use g_config
    use find_up_downwind_triangles_interface
    IMPLICIT NONE
    integer     :: n, k, n1, n2

    type(t_mesh),        intent(inout), target :: mesh
    type(t_partit),      intent(inout), target :: partit
    type(t_tracer_work), intent(inout), target :: twork

#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"

    !___________________________________________________________________________
    ! find upwind and downwind triangle for each local edge 
    call find_up_downwind_triangles(twork, partit, mesh)
    
    !___________________________________________________________________________
    nn_size=0
    k=0
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(n)
!$OMP DO REDUCTION(max: k)
    do n=1, myDim_nod2D
        ! get number of  neighbouring nodes from sparse stiffness matrix
        ! stiffnes matrix filled up in subroutine init_stiff_mat_ale
        ! --> SSH_stiff%rowptr... compressed row index of sparse matrix
        ! --> SSH_stiff%values... array with values at row column location, has legth nod2d+1
        ! --> SSH_stiff%rowptr(n) ... gives index location in SSH_stiff%values where the 
        !                             next value switches to a new row
        ! --> SSH_stiff%rowptr(n+1)-SSH_stiff%rowptr(n) gives maximum number of 
        !     neighbouring nodes within a single row of the sparse matrix
        k=max(k, SSH_stiff%rowptr(n+1)-SSH_stiff%rowptr(n))
    end do
!$OMP END DO    
!$OMP END PARALLEL
    nn_size=k
    !___________________________________________________________________________
    allocate(mesh%nn_num(myDim_nod2D), mesh%nn_pos(nn_size,myDim_nod2D))
    nn_num(1:myDim_nod2D)            => mesh%nn_num(:)
    nn_pos(1:nn_size, 1:myDim_nod2D) => mesh%nn_pos(:,:)
    ! These are the same arrays that we also use in quadratic reconstruction
    !MOVE IT TO SOMEWHERE ELSE
!$OMP PARALLEL DO
    do n=1, myDim_nod2d
        ! number of neigbouring nodes to node n
        nn_num(n)=1
        ! local position of neigbouring nodes
        nn_pos(1,n)=n
    end do   
!$OMP END PARALLEL DO
    !___________________________________________________________________________
    allocate(twork%nboundary_lay(myDim_nod2D+eDim_nod2D)) !node n becomes a boundary node after layer twork%nboundary_lay(n)
    twork%nboundary_lay=nl-1
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(n, k, n1, n2)
!$OMP DO
    do n=1, myDim_edge2D
        ! n1 and n2 are local indices 
        n1=edges(1,n)
        n2=edges(2,n)

#if defined(__openmp_reproducible)
!$OMP ORDERED
#endif

        ! ... if(n1<=myDim_nod2D) --> because dont use extended nodes
        if(n1<=myDim_nod2D) then
#if defined(_OPENMP)  && !defined(__openmp_reproducible)
           call omp_set_lock(partit%plock(n1))
#endif
            nn_pos(nn_num(n1)+1,n1)=n2
            nn_num(n1)=nn_num(n1)+1
#if defined(_OPENMP)  && !defined(__openmp_reproducible)
           call omp_unset_lock(partit%plock(n1))
#endif
        end if
        
        ! ... if(n2<=myDim_nod2D) --> because dont use extended nodes
        if(n2<=myDim_nod2D) then
#if defined(_OPENMP)  && !defined(__openmp_reproducible)
           call omp_set_lock(partit%plock(n2))
#endif
            nn_pos(nn_num(n2)+1,n2)=n1
            nn_num(n2)=nn_num(n2)+1
#if defined(_OPENMP)  && !defined(__openmp_reproducible)
           call omp_unset_lock(partit%plock(n2))
#endif
        end if
        
        if (any(edge_tri(:,n)<=0)) then
            ! this edge nodes is already at the surface at the boundary ...
            ! later here ...sign(1, twork%nboundary_lay(enodes(1))-nz) for nz=1 must be negativ
            ! thats why here twork%nboundary_lay(edges(:,n))=0
            twork%nboundary_lay(edges(:,n))=0
        else
            ! this edge nodes become boundary edge with increasing depth due to bottom topography
            ! at the depth twork%nboundary_lay the edge (edgepoints) still has two valid ocean triangles
            ! below that depth, edge becomes boundary edge
#if defined(_OPENMP)  && !defined(__openmp_reproducible)
           call omp_set_lock  (partit%plock(edges(1,n)))
#endif
            twork%nboundary_lay(edges(1,n))=min(twork%nboundary_lay(edges(1,n)), minval(nlevels(edge_tri(:,n)))-1)
#if defined(_OPENMP)  && !defined(__openmp_reproducible)
           call omp_unset_lock(partit%plock(edges(1,n)))
           call omp_set_lock  (partit%plock(edges(2,n)))
#endif
            twork%nboundary_lay(edges(2,n))=min(twork%nboundary_lay(edges(2,n)), minval(nlevels(edge_tri(:,n)))-1)
#if defined(_OPENMP)  && !defined(__openmp_reproducible)
           call omp_unset_lock(partit%plock(edges(2,n)))
#endif
        end if

#if defined(__openmp_reproducible)
!$OMP END ORDERED
#endif
    end do
!$OMP END DO
!$OMP END PARALLEL
end SUBROUTINE muscl_adv_init
!
!
!_______________________________________________________________________________
SUBROUTINE find_up_downwind_triangles(twork, partit, mesh)
USE MOD_MESH
USE MOD_PARTIT
USE MOD_PARSUP
USE MOD_TRACER
USE o_ARRAYS
USE o_PARAM
USE g_CONFIG
use g_comm_auto
IMPLICIT NONE
integer                    :: n, k, ednodes(2), elem, el
real(kind=WP)              :: x(2),b(2), c(2), cr, bx, by, xx, xy, ab, ax
real(kind=WP), allocatable :: coord_elem(:, :,:), temp(:)
integer, allocatable       :: temp_i(:), e_nodes(:,:)

type(t_mesh),        intent(in)   , target :: mesh
type(t_partit),      intent(inout), target :: partit
type(t_tracer_work), intent(inout), target :: twork
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"

allocate(twork%edge_up_dn_tri(2,myDim_edge2D))
allocate(twork%edge_up_dn_grad(4,nl-1,myDim_edge2D))
twork%edge_up_dn_tri=0
! =====
! In order that this procedure works, we need to know nodes and their coordinates 
! on the extended set of elements (not only my, but myDim+eDim+eXDim) 
! =====
allocate(coord_elem(2, 3, myDim_elem2D+eDim_elem2D+eXDim_elem2D))
allocate(temp(myDim_elem2D+eDim_elem2D+eXDim_elem2D))
   DO n=1,3
        DO k=1,2
!$OMP PARALLEL
!$OMP DO
           DO el=1,myDim_elem2D
              temp(el)=coord_nod2D(k,elem2D_nodes(n,el))
           END DO
!$OMP END DO
!$OMP MASTER
	   call exchange_elem(temp, partit)
!$OMP END MASTER
!$OMP BARRIER
!$OMP DO
           DO el=1, myDim_elem2D+eDim_elem2D+eXDim_elem2D
              coord_elem(k,n,el)=temp(el)
           END DO
!$OMP END DO
!$OMP END PARALLEL
	END DO
   END DO
deallocate(temp)

allocate(e_nodes(3,myDim_elem2D+eDim_elem2D+eXDim_elem2D))
allocate(temp_i(myDim_elem2D+eDim_elem2D+eXDim_elem2D))
    DO n=1,3
!$OMP PARALLEL
!$OMP DO
       do el=1,myDim_elem2D
          temp_i(el)=myList_nod2D(elem2D_nodes(n,el))
       end do
!$OMP END DO
!$OMP MASTER
       call exchange_elem(temp_i, partit)
!$OMP END MASTER
!$OMP BARRIER
!$OMP DO
       DO el=1, myDim_elem2D+eDim_elem2D+eXDim_elem2D
          e_nodes(n, el)=temp_i(el)
       END DO
!$OMP END DO
!$OMP END PARALLEL
    END DO   
deallocate(temp_i)

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(n, k, ednodes, elem, el, x,b, c, cr, bx, by, xx, xy, ab, ax)
!$OMP DO
DO n=1, myDim_edge2d
   ednodes=edges(:,n) 
   x=coord_nod2D(:,ednodes(2))-coord_nod2D(:,ednodes(1))
      	 if(x(1)>cyclic_length/2._WP) x(1)=x(1)-cyclic_length
         if(x(1)<-cyclic_length/2._WP) x(1)=x(1)+cyclic_length
	
   ! Find upwind (in the sense of x) triangle, i. e. 
   ! find which triangle contains -x:
   x=-x
   DO k=1,nod_in_elem2D_num(ednodes(1))
      elem=nod_in_elem2D(k,ednodes(1))
      if(e_nodes(1,elem)==myList_nod2D(ednodes(1))) then
	 b=coord_elem(:,2,elem)-coord_elem(:,1,elem)
	 c=coord_elem(:,3,elem)-coord_elem(:,1,elem)
      elseif(e_nodes(2,elem)==myList_nod2D(ednodes(1))) then
	 b=coord_elem(:,1,elem)-coord_elem(:,2,elem)
	 c=coord_elem(:,3,elem)-coord_elem(:,2,elem)
      else	 
	 b=coord_elem(:,1,elem)-coord_elem(:,3,elem)
	 c=coord_elem(:,2,elem)-coord_elem(:,3,elem)
      end if
      	 if(b(1)>cyclic_length/2._WP) b(1)=b(1)-cyclic_length
         if(b(1)<-cyclic_length/2._WP) b(1)=b(1)+cyclic_length
	 if(c(1)>cyclic_length/2._WP) c(1)=c(1)-cyclic_length
         if(c(1)<-cyclic_length/2._WP) c(1)=c(1)+cyclic_length
      ! the vector x has to be between b and c
      ! Decompose b and x into parts along c and along (-cy,cx), i.e.
      ! 90 degree counterclockwise
      cr=sum(c*c)
      bx=sum(b*c)/cr
      by=(-b(1)*c(2)+b(2)*c(1))/cr
      xx=sum(x*c)/cr
      xy=(-x(1)*c(2)+x(2)*c(1))/cr
      ab=atan2(by,bx)
      ax=atan2(xy,xx)
      ! Since b and c are the sides of triangle, |ab|<pi, and atan2 should 
      ! be what is needed
      if((ab>0.0_WP).and.(ax>0.0_WP).and.(ax<ab)) then
      twork%edge_up_dn_tri(1,n)=elem
      cycle
      endif
      if((ab<0.0_WP).and.(ax<0.0_WP).and.(ax>ab)) then
      twork%edge_up_dn_tri(1,n)=elem
      cycle
      endif
      if((ab==ax).or.(ax==0.0_WP)) then
      twork%edge_up_dn_tri(1,n)=elem
      cycle
      endif
   END DO
   ! Find downwind element
   x=-x
   DO k=1,nod_in_elem2D_num(ednodes(2))
      elem=nod_in_elem2D(k,ednodes(2))
      if(e_nodes(1,elem)==myList_nod2D(ednodes(2))) then
      	 b=coord_elem(:,2,elem)-coord_elem(:,1,elem)
	 c=coord_elem(:,3,elem)-coord_elem(:,1,elem)
      elseif(e_nodes(2, elem)==myList_nod2D(ednodes(2))) then
	 b=coord_elem(:,1,elem)-coord_elem(:,2,elem)
	 c=coord_elem(:,3,elem)-coord_elem(:,2,elem)
      else	 
	 b=coord_elem(:,1,elem)-coord_elem(:,3,elem)
	 c=coord_elem(:,2,elem)-coord_elem(:,3,elem)
      end if
      	 if(b(1)>cyclic_length/2.) b(1)=b(1)-cyclic_length
         if(b(1)<-cyclic_length/2.) b(1)=b(1)+cyclic_length
	 if(c(1)>cyclic_length/2.) c(1)=c(1)-cyclic_length
         if(c(1)<-cyclic_length/2.) c(1)=c(1)+cyclic_length
      ! the vector x has to be between b and c
      ! Decompose b and x into parts along c and along (-cy,cx), i.e.
      ! 90 degree counterclockwise
      cr=sum(c*c)
      bx=sum(b*c)/cr
      by=(-b(1)*c(2)+b(2)*c(1))/cr
      xx=sum(x*c)/cr
      xy=(-x(1)*c(2)+x(2)*c(1))/cr
      ab=atan2(by,bx)
      ax=atan2(xy,xx)
      ! Since b and c are the sides of triangle, |ab|<pi, and atan2 should 
      ! be what is needed
      if((ab>0.0_WP).and.(ax>0.0_WP).and.(ax<ab)) then
      twork%edge_up_dn_tri(2,n)=elem
      cycle
      endif
      if((ab<0.0_WP).and.(ax<0.0_WP).and.(ax>ab)) then
      twork%edge_up_dn_tri(2,n)=elem
      cycle
      endif
      if((ab==ax).or.(ax==0.0)) then
      twork%edge_up_dn_tri(2,n)=elem
      cycle
      endif
   END DO
END DO
!$OMP END DO
!$OMP END PARALLEL

! For edges touching the boundary --- up or downwind elements may be absent.  
! We return to the standard Miura at nodes that
! belong to such edges. Same at the depth.
! Count the number of 'good' edges:
!k=0
!DO n=1, myDim_edge2D
!   if((twork%edge_up_dn_tri(1,n).ne.0).and.(twork%edge_up_dn_tri(2,n).ne.0)) k=k+1
!END DO

!$OMP PARALLEL DO
DO n=1, myDim_edge2D
   twork%edge_up_dn_grad(:, :, n)=0.0_WP
END DO
!$OMP END PARALLEL DO
deallocate(e_nodes, coord_elem)
end SUBROUTINE find_up_downwind_triangles
!
!
!_______________________________________________________________________________
SUBROUTINE fill_up_dn_grad(twork, partit, mesh)
! ttx, tty  elemental gradient of tracer 
USE o_PARAM
USE MOD_MESH
USE MOD_PARTIT
USE MOD_PARSUP
USE MOD_TRACER
USE o_ARRAYS
IMPLICIT NONE
integer                  :: edge, n, nz, elem, k, ednodes(2), nzmin, nzmax
real(kind=WP)            :: tvol, tx, ty
type(t_mesh),        intent(in),    target :: mesh
type(t_partit),      intent(inout), target :: partit
type(t_tracer_work), intent(inout), target :: twork
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"
	!___________________________________________________________________________
	! loop over edge segments
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(edge, n, nz, elem, k, ednodes, nzmin, nzmax, tvol, tx, ty)
!$OMP DO
	DO edge=1, myDim_edge2D
		ednodes=edges(:,edge)
		!_______________________________________________________________________
		! case when edge has upwind and downwind triangle on the surface
		if((twork%edge_up_dn_tri(1,edge).ne.0.0_WP).and.(twork%edge_up_dn_tri(2,edge).ne.0.0_WP)) then
			nzmin = maxval(ulevels_nod2D_max(ednodes))
			nzmax = minval(nlevels_nod2D_min(ednodes))
			
			!___________________________________________________________________
			! loop over not shared depth levels of edge node 1 (ednodes(1))
			DO nz=ulevels_nod2D(ednodes(1)), nzmin-1
				tvol=0.0_WP
				tx=0.0_WP
				ty=0.0_WP
				! loop over number triangles that share the nodeedge points ednodes(1)
				! --> calculate mean gradient at ednodes(1) over the sorounding 
				!     triangle gradients
				DO k=1, nod_in_elem2D_num(ednodes(1))
					elem=nod_in_elem2D(k,ednodes(1))
					!!PS if(nlevels(elem)-1 < nz) cycle
					if(nlevels(elem)-1<nz .or. nz<ulevels(elem)) cycle
					tvol=tvol+elem_area(elem)
					tx=tx+tr_xy(1,nz,elem)*elem_area(elem)
					ty=ty+tr_xy(2,nz,elem)*elem_area(elem)
				END DO
				twork%edge_up_dn_grad(1,nz,edge)=tx/tvol
				twork%edge_up_dn_grad(3,nz,edge)=ty/tvol
			END DO
			
			!___________________________________________________________________
			! loop over not shared depth levels of edge node 2 (ednodes(2))
			DO nz=ulevels_nod2D(ednodes(2)),nzmin-1
				tvol=0.0_WP
				tx=0.0_WP
				ty=0.0_WP
				! loop over number triangles that share the nodeedge points ednodes(2)
				! --> calculate mean gradient at ednodes(2) over the sorounding 
				!     triangle gradients
				DO k=1, nod_in_elem2D_num(ednodes(2))
					elem=nod_in_elem2D(k,ednodes(2))
					!!PS if(nlevels(elem)-1 < nz) cycle
					if(nlevels(elem)-1<nz .or. nz<ulevels(elem)) cycle
					tvol=tvol+elem_area(elem)
					tx=tx+tr_xy(1,nz,elem)*elem_area(elem)
					ty=ty+tr_xy(2,nz,elem)*elem_area(elem)
				END DO
				twork%edge_up_dn_grad(2,nz,edge)=tx/tvol
				twork%edge_up_dn_grad(4,nz,edge)=ty/tvol
			END DO
			
			!___________________________________________________________________
			! loop over shared depth levels
			!!PS DO nz=1, minval(nlevels_nod2D_min(ednodes))-1
			DO nz=nzmin, nzmax-1
				! tracer gradx for upwind and downwind tri
				twork%edge_up_dn_grad(1:2,nz,edge)=tr_xy(1,nz,twork%edge_up_dn_tri(:,edge))
				! tracer grady for upwind and downwind tri
				twork%edge_up_dn_grad(3:4,nz,edge)=tr_xy(2,nz,twork%edge_up_dn_tri(:,edge))
			END DO
			
			!___________________________________________________________________
			! loop over not shared depth levels of edge node 1 (ednodes(1))
			!!PS DO nz=minval(nlevels_nod2D_min(ednodes)),nlevels_nod2D(ednodes(1))-1
			DO nz=nzmax, nlevels_nod2D(ednodes(1))-1
				tvol=0.0_WP
				tx=0.0_WP
				ty=0.0_WP
				! loop over number triangles that share the nodeedge points ednodes(1)
				! --> calculate mean gradient at ednodes(1) over the sorounding 
				!     triangle gradients
				DO k=1, nod_in_elem2D_num(ednodes(1))
					elem=nod_in_elem2D(k,ednodes(1))
					!!PS if(nlevels(elem)-1 < nz) cycle
					if(nlevels(elem)-1<nz .or. nz<ulevels(elem)) cycle
					tvol=tvol+elem_area(elem)
					tx=tx+tr_xy(1,nz,elem)*elem_area(elem)
					ty=ty+tr_xy(2,nz,elem)*elem_area(elem)
				END DO
				twork%edge_up_dn_grad(1,nz,edge)=tx/tvol
				twork%edge_up_dn_grad(3,nz,edge)=ty/tvol
			END DO
			!___________________________________________________________________
			! loop over not shared depth levels of edge node 2 (ednodes(2))
			!!PS DO nz=minval(nlevels_nod2D_min(ednodes)),nlevels_nod2D(ednodes(2))-1
			DO nz=nzmax, nlevels_nod2D(ednodes(2))-1
				tvol=0.0_WP
				tx=0.0_WP
				ty=0.0_WP
				! loop over number triangles that share the nodeedge points ednodes(2)
				! --> calculate mean gradient at ednodes(2) over the sorounding 
				!     triangle gradients
				DO k=1, nod_in_elem2D_num(ednodes(2))
					elem=nod_in_elem2D(k,ednodes(2))
					!!PS if(nlevels(elem)-1 < nz) cycle
					if(nlevels(elem)-1<nz .or. nz<ulevels(elem)) cycle
					tvol=tvol+elem_area(elem)
					tx=tx+tr_xy(1,nz,elem)*elem_area(elem)
					ty=ty+tr_xy(2,nz,elem)*elem_area(elem)
				END DO
				twork%edge_up_dn_grad(2,nz,edge)=tx/tvol
				twork%edge_up_dn_grad(4,nz,edge)=ty/tvol
			END DO
		!_______________________________________________________________________
		! case when edge either upwind or downwind triangle on the surface
		! --> surface boundary edge
		else
			! Only linear reconstruction part
			nzmin = ulevels_nod2D(ednodes(1))
			nzmax = nlevels_nod2D(ednodes(1))
			!!PS DO nz=1,nlevels_nod2D(ednodes(1))-1
			DO nz=nzmin,nzmax-1
				tvol=0.0_WP
				tx=0.0_WP
				ty=0.0_WP
				DO k=1, nod_in_elem2D_num(ednodes(1))
					elem=nod_in_elem2D(k,ednodes(1))
					!!PS if(nlevels(elem)-1 < nz) cycle
					if(nlevels(elem)-1 < nz .or. nz<ulevels(elem) ) cycle
					tvol=tvol+elem_area(elem)
					tx=tx+tr_xy(1,nz,elem)*elem_area(elem)
					ty=ty+tr_xy(2,nz,elem)*elem_area(elem)
				END DO
				twork%edge_up_dn_grad(1,nz,edge)=tx/tvol
				twork%edge_up_dn_grad(3,nz,edge)=ty/tvol
			END DO
			nzmin = ulevels_nod2D(ednodes(2))
			nzmax = nlevels_nod2D(ednodes(2))
			!!PS DO nz=1,nlevels_nod2D(ednodes(2))-1
			DO nz=nzmin,nzmax-1
				tvol=0.0_WP
				tx=0.0_WP
				ty=0.0_WP
				DO k=1, nod_in_elem2D_num(ednodes(2))
					elem=nod_in_elem2D(k,ednodes(2))
					!!PS if(nlevels(elem)-1 < nz) cycle
					if(nlevels(elem)-1 < nz .or. nz<ulevels(elem) ) cycle
					tvol=tvol+elem_area(elem)
					tx=tx+tr_xy(1,nz,elem)*elem_area(elem)
					ty=ty+tr_xy(2,nz,elem)*elem_area(elem)
				END DO
				twork%edge_up_dn_grad(2,nz,edge)=tx/tvol
				twork%edge_up_dn_grad(4,nz,edge)=ty/tvol
			END DO
		end if  
	END DO
!$OMP END DO
!$OMP END PARALLEL
END SUBROUTINE fill_up_dn_grad
