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
SUBROUTINE muscl_adv_init
USE o_MESH
USE o_ARRAYS
USE o_PARAM
USE g_PARSUP
  use g_comm_auto
use g_config
IMPLICIT NONE
integer     :: n, k, n1, n2, n_num

 call find_up_downwind_triangles
 n_num=0
  DO n=1, myDim_nod2D
  k=SSH_stiff%rowptr(n+1)-SSH_stiff%rowptr(n)
  if(k>n_num) n_num=k
  END DO
  allocate(nlevels_nod2D_min(myDim_nod2D+eDim_nod2D))
  allocate(nn_num(myDim_nod2D), nn_pos(n_num,myDim_nod2D))
                    !! These are the same arrays that we also use in quadratic
		    !! reconstruction
  DO n=1,myDim_nod2d
     nn_num(n)=1
     nn_pos(1,n)=n
  end do   
  Do n=1, myDim_edge2D
     n1=edges(1,n)
     n2=edges(2,n)
     if(n1<=myDim_nod2D) then
     nn_pos(nn_num(n1)+1,n1)=n2
     nn_num(n1)=nn_num(n1)+1
     end if
     if(n2<=myDim_nod2D) then
     nn_pos(nn_num(n2)+1,n2)=n1
     nn_num(n2)=nn_num(n2)+1
     end if
  END DO  
  DO n=1,myDim_nod2D
    nlevels_nod2D_min(n)=minval(nlevels_nod2D(nn_pos(1:nn_num(n),n))) 
  END DO
  call exchange_nod(nlevels_nod2D_min)
end SUBROUTINE muscl_adv_init
!=======================================================================
SUBROUTINE find_up_downwind_triangles
USE o_MESH
USE o_ARRAYS
USE o_PARAM
USE g_PARSUP
USE g_CONFIG
use g_comm_auto
IMPLICIT NONE
integer                    :: n, k, ednodes(2), elem
real(kind=WP)              :: x(2),b(2), c(2), cr, bx, by, xx, xy, ab, ax
real(kind=WP), allocatable :: coord_elem(:, :,:), temp(:)
integer, allocatable       :: temp_i(:), e_nodes(:,:)

allocate(edge_up_dn_tri(2,myDim_edge2D))
allocate(edge_up_dn_grad(4,nl-1,myDim_edge2D))
edge_up_dn_tri=0
! =====
! In order that this procedure works, we need to know nodes and their coordinates 
! on the extended set of elements (not only my, but myDim+eDim+eXDim) 
! =====
allocate(coord_elem(2,3,myDim_elem2D+eDim_elem2D+eXDim_elem2D))
allocate(temp(myDim_elem2D+eDim_elem2D+eXDim_elem2D))
   DO n=1,3
        DO k=1,2
           temp(1:myDim_elem2D)=coord_nod2D(k,elem2D_nodes(n,:))
	   call exchange_elem(temp)
	   coord_elem(k,n,:)=temp(:)
	END DO
   END DO
deallocate(temp)

allocate(e_nodes(3,myDim_elem2D+eDim_elem2D+eXDim_elem2D))
allocate(temp_i(myDim_elem2D+eDim_elem2D+eXDim_elem2D))
    DO n=1,3
       temp_i(1:myDim_elem2D)=myList_nod2D(elem2D_nodes(n,:))
       call exchange_elem(temp_i)
       e_nodes(n,:)=temp_i(:)
    END DO   
deallocate(temp_i)


DO n=1, myDim_edge2d
   ednodes=edges(:,n) 
   x=coord_nod2D(:,ednodes(2))-coord_nod2D(:,ednodes(1))
      	 if(x(1)>cyclic_length/2.) x(1)=x(1)-cyclic_length
         if(x(1)<-cyclic_length/2.) x(1)=x(1)+cyclic_length
	
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
      if((ab>0).and.(ax>0).and.(ax<ab)) then
      edge_up_dn_tri(1,n)=elem
      cycle
      endif
      if((ab<0).and.(ax<0).and.(ax>ab)) then
      edge_up_dn_tri(1,n)=elem
      cycle
      endif
      if((ab==ax).or.(ax==0.0)) then
      edge_up_dn_tri(1,n)=elem
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
      if((ab>0).and.(ax>0).and.(ax<ab)) then
      edge_up_dn_tri(2,n)=elem
      cycle
      endif
      if((ab<0).and.(ax<0).and.(ax>ab)) then
      edge_up_dn_tri(2,n)=elem
      cycle
      endif
      if((ab==ax).or.(ax==0.0)) then
      edge_up_dn_tri(2,n)=elem
      cycle
      endif
   END DO
END DO
! For edges touching the boundary --- up or downwind elements may be absent.  
! We return to the standard Miura at nodes that
! belong to such edges. Same at the depth. 

! Count the number of 'good' edges:
k=0 
DO n=1,myDim_edge2D
   if((edge_up_dn_tri(1,n).ne.0).and.(edge_up_dn_tri(2,n).ne.0)) k=k+1
END DO

deallocate(e_nodes, coord_elem)


edge_up_dn_grad=0.0
end SUBROUTINE find_up_downwind_triangles
!=======================================================================
SUBROUTINE fill_up_dn_grad

! ttx, tty  elemental gradient of tracer 
USE o_PARAM
USE o_MESH
USE o_ARRAYS
USE g_PARSUP
IMPLICIT NONE
integer        :: n, nz, elem, k, edge, ednodes(2)
real(kind=8)   :: tvol, tx, ty

  DO edge=1,myDim_edge2D
     ednodes=edges(:,edge)
    if((edge_up_dn_tri(1,edge).ne.0).and.(edge_up_dn_tri(2,edge).ne.0)) then
    
     DO nz=1, minval(nlevels_nod2D_min(ednodes))-1
        edge_up_dn_grad(1:2,nz,edge)=tt_xy_stored(1,nz,edge_up_dn_tri(:,edge))
	edge_up_dn_grad(3:4,nz,edge)=tt_xy_stored(2,nz,edge_up_dn_tri(:,edge))
     END DO
     DO nz=minval(nlevels_nod2D_min(ednodes)),nlevels_nod2D(ednodes(1))-1
        tvol=0.0
	tx=0.0
	ty=0.0
	DO k=1, nod_in_elem2D_num(ednodes(1))
           elem=nod_in_elem2D(k,ednodes(1))
           if(nlevels(elem)-1<nz) cycle
           tvol=tvol+elem_area(elem)
           tx=tx+tt_xy_stored(1,nz,elem)*elem_area(elem)
	   ty=ty+tt_xy_stored(2,nz,elem)*elem_area(elem)
        END DO
	edge_up_dn_grad(1,nz,edge)=tx/tvol
	edge_up_dn_grad(3,nz,edge)=ty/tvol
     END DO
     DO nz=minval(nlevels_nod2D_min(ednodes)),nlevels_nod2D(ednodes(2))-1
        tvol=0.0
	tx=0.0
	ty=0.0
	DO k=1, nod_in_elem2D_num(ednodes(2))
           elem=nod_in_elem2D(k,ednodes(2))
           if(nlevels(elem)-1<nz) cycle
           tvol=tvol+elem_area(elem)
           tx=tx+tt_xy_stored(1,nz,elem)*elem_area(elem)
	   ty=ty+tt_xy_stored(2,nz,elem)*elem_area(elem)
        END DO
	edge_up_dn_grad(2,nz,edge)=tx/tvol
	edge_up_dn_grad(4,nz,edge)=ty/tvol
     END DO
    else
    
      ! Only linear reconstruction part
     DO nz=1,nlevels_nod2D(ednodes(1))-1
        tvol=0.0
	tx=0.0
	ty=0.0
	DO k=1, nod_in_elem2D_num(ednodes(1))
           elem=nod_in_elem2D(k,ednodes(1))
           if(nlevels(elem)-1<nz) cycle
           tvol=tvol+elem_area(elem)
           tx=tx+tt_xy_stored(1,nz,elem)*elem_area(elem)
	   ty=ty+tt_xy_stored(2,nz,elem)*elem_area(elem)
        END DO
	edge_up_dn_grad(1,nz,edge)=tx/tvol
	edge_up_dn_grad(3,nz,edge)=ty/tvol
     END DO
     DO nz=1,nlevels_nod2D(ednodes(2))-1
        tvol=0.0
	tx=0.0
	ty=0.0
	DO k=1, nod_in_elem2D_num(ednodes(2))
           elem=nod_in_elem2D(k,ednodes(2))
           if(nlevels(elem)-1<nz) cycle
           tvol=tvol+elem_area(elem)
           tx=tx+tt_xy_stored(1,nz,elem)*elem_area(elem)
	   ty=ty+tt_xy_stored(2,nz,elem)*elem_area(elem)
        END DO
	edge_up_dn_grad(2,nz,edge)=tx/tvol
	edge_up_dn_grad(4,nz,edge)=ty/tvol
     END DO
    end if  
   END DO 
   
END SUBROUTINE fill_up_dn_grad
!===========================================================================
! It is assumed that velocity is at n+1/2, hence only tracer field 
! is AB2 interpolated to n+1/2. 
SUBROUTINE adv_tracer_muscl(ttf, dttf, ttfold)
USE o_MESH
USE o_ARRAYS
USE o_PARAM
USE g_PARSUP
USE g_CONFIG
use g_comm_auto
IMPLICIT NONE
 integer      :: el(2), enodes(2), n, nz, edge
 integer      :: nl1, nl2,tr_num
 real(kind=8) :: c1, c2, deltaX1, deltaY1, deltaX2, deltaY2, flux=0.0 
 real(kind=8) :: tvert(nl), a, b, c, d, da, db, dg
 real(kind=8) :: Tx, Ty, Tmean, rdata=0.0
 real(kind=8) :: ttf(nl-1, myDim_nod2D+eDim_nod2D), dttf(nl-1, myDim_nod2D+eDim_nod2D)
 real(kind=8) :: ttfold(nl-1, myDim_nod2D+eDim_nod2D)

! Clean the rhs
ttrhs=0d0  
! Horizontal advection
  DO edge=1, myDim_edge2D
   enodes=edges(:,edge)   
   el=edge_tri(:,edge)
   c1=0.0
   c2=0.0
   nl1=nlevels(el(1))-1
   deltaX1=edge_cross_dxdy(1,edge)
   deltaY1=edge_cross_dxdy(2,edge)
   nl2=0
   a=r_earth*elem_cos(el(1))
   if(el(2)>0) then
   deltaX2=edge_cross_dxdy(3,edge)
   deltaY2=edge_cross_dxdy(4,edge)
   nl2=nlevels(el(2))-1
   b=r_earth*elem_cos(el(2))
   end if     
   
   ! ============
   ! First segment
   ! ============
   DO nz=1, nl1
   ! ============
   ! MUSCL type reconstruction
   ! ============
   if(UV(2,nz,el(1))*deltaX1- UV(1,nz,el(1))*deltaY1>0) then   
      Tmean=ttfold(nz, enodes(2))- &
      (2.0_8*(ttfold(nz, enodes(2))-ttfold(nz,enodes(1)))+ &
      edge_dxdy(1,edge)*a*edge_up_dn_grad(2,nz,edge)+ &
      edge_dxdy(2,edge)*r_earth*edge_up_dn_grad(4,nz,edge))/6.0_8    
   else
      Tmean=ttfold(nz, enodes(1))+ &
      (2.0_8*(ttfold(nz, enodes(2))-ttfold(nz,enodes(1)))+ &
      edge_dxdy(1,edge)*a*edge_up_dn_grad(1,nz,edge)+ &
      edge_dxdy(2,edge)*r_earth*edge_up_dn_grad(3,nz,edge))/6.0_8    
   end if
   c1=UV(2,nz,el(1))*Tmean*deltaX1- UV(1,nz,el(1))*Tmean*deltaY1
   ttrhs(nz,enodes(1))=ttrhs(nz,enodes(1))+c1
   ttrhs(nz,enodes(2))=ttrhs(nz,enodes(2))-c1
   END DO
   ! ============
   ! Second segment
   ! ============
   if(el(2)>0)then
   DO nz=1, nl2
   ! ============
   ! Linear upwind reconstruction
   ! ============
   if(UV(2,nz,el(2))*deltaX2- UV(1,nz,el(2))*deltaY2<0) then   
      Tmean=ttfold(nz, enodes(2))- &
      (2.0_8*(ttfold(nz, enodes(2))-ttfold(nz,enodes(1)))+ &
      edge_dxdy(1,edge)*b*edge_up_dn_grad(2,nz,edge)+ &
      edge_dxdy(2,edge)*r_earth*edge_up_dn_grad(4,nz,edge))/6.0_8    
   else
      Tmean=ttfold(nz, enodes(1))+ &
      (2.0_8*(ttfold(nz, enodes(2))-ttfold(nz,enodes(1)))+ &
      edge_dxdy(1,edge)*b*edge_up_dn_grad(1,nz,edge)+ &
      edge_dxdy(2,edge)*r_earth*edge_up_dn_grad(3,nz,edge))/6.0_8     
   end if
   c2=-UV(2,nz,el(2))*Tmean*deltaX2+ UV(1,nz,el(2))*Tmean*deltaY2
   ttrhs(nz,enodes(1))=ttrhs(nz,enodes(1))+c2
   ttrhs(nz,enodes(2))=ttrhs(nz,enodes(2))-c2
   END DO
   end if
   
 END DO

! ===================
! Vertical advection
! ===================
 
 DO n=1, myDim_nod2D                  !! P (d) n=1, nod2D

  ! ===========
  ! Fluxes in the column
  ! ===========
 
   tvert(1)= -Wvel(1,n)*ttfold(1,n)*area(1,n)
		    
  ! Bottom conditions	  
   tvert(nlevels_nod2D(n))=0.	    
  
   DO nz=2, nlevels_nod2D(n)-1
      ! ============
      ! QUICK upwind (3rd order)
      ! ============
      if(Wvel(nz,n)>0) then
        if(nz==nlevels_nod2D(n)-1) then
	  Tmean=0.5_8*(ttfold(nz-1,n)+ttfold(nz,n))  ! or replace this with 
	                                             ! the first order 
						     ! upwind  tttfold(nz,n)
	else
	a=Z(nz-1)-zbar(nz)
	b=zbar(nz)-Z(nz)
	c=zbar(nz)-Z(nz+1)
	dg=a*b/(c+a)/(b-c)
	db=-a*c/(b+a)/(b-c)
	da=1.0-dg-db
	Tmean=ttfold(nz-1,n)*da+ttfold(nz,n)*db+ttfold(nz+1,n)*dg
	end if
      end if

      if(Wvel(nz,n)<0) then
        if(nz==2) then
	  Tmean=0.5_8*(ttfold(nz-1,n)+ttfold(nz,n))        ! or ttfold(nz-1,n)
	else  
	a=zbar(nz)-Z(nz)
	b=Z(nz-1)-zbar(nz)
	c=Z(nz-2)-zbar(nz)
	dg=a*b/(c+a)/(b-c)
	db=-a*c/(b+a)/(b-c)
	da=1.0-dg-db
	Tmean=ttfold(nz,n)*da+ttfold(nz-1,n)*db+ttfold(nz-2,n)*dg
	end if
      end if
      tvert(nz)= -Tmean*Wvel(nz,n)*area(nz,n)
   END DO
 
   DO nz=1,nlevels_nod2D(n)-1
      ttrhs(nz,n)=(ttrhs(nz,n)+ &
                    (tvert(nz)-tvert(nz+1))/(zbar(nz)-zbar(nz+1)))
   END DO
 END DO

! =================
! Update ttfold (to be used on the next time level)
! and compute new ttf
! =================          
  
  DO n=1, myDim_nod2D+eDim_nod2D      !! P (h)  n=1, nod2D
                                      !! n=myList_nod2D(m)
     DO nz=1,nlevels_nod2D(n)-1
        dttf(nz,n)=dttf(nz,n)+ttrhs(nz,n)*dt/area(nz,n)
     END DO
  END DO
end subroutine adv_tracer_muscl

!===========================================================================
