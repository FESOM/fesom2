!=======================================================================
! Advection based on the gradient reconstruction.
! The last argument 'num_ord<1' defines the share of the
! 4th order centered contribution, and (1-num_ord) is done with 3rd order upwind.
! Dissipation comes only from the first part. num_ord=0.75--0.85 is 
! recommended if stable. 
! It is assumed that velocity is at n+1/2, where n is time step, hence only tracer field 
! is AB2 interpolated to n+1/2. 
SUBROUTINE adv_tracer_muscl34(dttf, ttfold, num_ord)
USE o_MESH
USE o_ARRAYS
USE o_PARAM
USE g_PARSUP
USE g_CONFIG
use g_comm_auto
IMPLICIT NONE
 integer      :: el(2), enodes(2), n, nz, edge
 integer      :: nl1, nl2,n2
 real(kind=8) :: c1, a, deltaX1, deltaY1, deltaX2, deltaY2, vflux=0.0 
 real(kind=8) :: Tmean, Tmean1, Tmean2
 real(kind=WP):: tvert(nl), b, c, da, db, dg
 real(kind=8) :: ttfold(nl-1, myDim_nod2D+eDim_nod2D), dttf(nl-1, myDim_nod2D+eDim_nod2D)
 real(kind=8) :: num_ord

! Clean the rhs
  ttrhs=0d0  
! =================
! Horizontal advection
! =================
  DO edge=1, myDim_edge2D
   enodes=edges(:,edge)   
   el=edge_tri(:,edge)
   c1=0.0
   nl1=nlevels(el(1))-1
   deltaX1=edge_cross_dxdy(1,edge)
   deltaY1=edge_cross_dxdy(2,edge)
   nl2=0
   a=r_earth*elem_cos(el(1))
   if(el(2)>0) then
   deltaX2=edge_cross_dxdy(3,edge)
   deltaY2=edge_cross_dxdy(4,edge)
   nl2=nlevels(el(2))-1
   a=0.5_8*(a+r_earth*elem_cos(el(2)))
   end if     
   n2=min(nl1,nl2)
   ! ============
   ! Both segments 
   ! ============
   DO nz=1, n2
   ! ============
   ! MUSCL-type reconstruction
   ! ============
   Tmean2=ttfold(nz, enodes(2))- &
      (2.0_8*(ttfold(nz, enodes(2))-ttfold(nz,enodes(1)))+ &
      edge_dxdy(1,edge)*a*edge_up_dn_grad(2,nz,edge)+ &
      edge_dxdy(2,edge)*r_earth*edge_up_dn_grad(4,nz,edge))/6.0_8    
  
   Tmean1=ttfold(nz, enodes(1))+ &
      (2.0_8*(ttfold(nz, enodes(2))-ttfold(nz,enodes(1)))+ &
      edge_dxdy(1,edge)*a*edge_up_dn_grad(1,nz,edge)+ &
      edge_dxdy(2,edge)*r_earth*edge_up_dn_grad(3,nz,edge))/6.0_8    
   ! ------
   ! volume flux across the segments
   ! ------
   vflux=-UV(2,nz,el(1))*deltaX1 + UV(1,nz,el(1))*deltaY1 &
         +UV(2,nz,el(2))*deltaX2 - UV(1,nz,el(2))*deltaY2
   ! ------
   ! tracer flux upwind
   ! ------
   c1=(vflux+abs(vflux))*Tmean1+(vflux-abs(vflux))*Tmean2
   ! ------
   ! combined with centered
   ! ------
   c1=-0.5_8*((1.0_8-num_ord)*c1+vflux*num_ord*(Tmean1+Tmean2))
   ttrhs(nz,enodes(1))= ttrhs(nz,enodes(1))+c1
   ttrhs(nz,enodes(2))= ttrhs(nz,enodes(2))-c1
   END DO
   ! ============
   ! remaining segments on the left or on the right
   ! ============
   if(nl1>nl2) then 
   DO nz=1+n2,nl1
   Tmean2=ttfold(nz, enodes(2))- &
      (2.0_8*(ttfold(nz, enodes(2))-ttfold(nz,enodes(1)))+ &
      edge_dxdy(1,edge)*a*edge_up_dn_grad(2,nz,edge)+ &
      edge_dxdy(2,edge)*r_earth*edge_up_dn_grad(4,nz,edge))/6.0_8    
  
   Tmean1=ttfold(nz, enodes(1))+ &
      (2.0_8*(ttfold(nz, enodes(2))-ttfold(nz,enodes(1)))+ &
      edge_dxdy(1,edge)*a*edge_up_dn_grad(1,nz,edge)+ &
      edge_dxdy(2,edge)*r_earth*edge_up_dn_grad(3,nz,edge))/6.0_8    
   ! ------
   ! volume flux across the segments
   ! ------
   vflux=-UV(2,nz,el(1))*deltaX1 + UV(1,nz,el(1))*deltaY1 
   ! ------
   ! tracer flux upwind
   ! ------
   c1=(vflux+abs(vflux))*Tmean1+(vflux-abs(vflux))*Tmean2
   ! ------
   ! combined with centered
   ! ------
   c1=-0.5_8*((1.0_8-num_ord)*c1+vflux*num_ord*(Tmean1+Tmean2))
   ttrhs(nz,enodes(1))= ttrhs(nz,enodes(1))+c1
   ttrhs(nz,enodes(2))= ttrhs(nz,enodes(2))-c1
   END DO
   else
   DO nz=n2+1,nl2
   Tmean2=ttfold(nz, enodes(2))- &
      (2.0_8*(ttfold(nz, enodes(2))-ttfold(nz,enodes(1)))+ &
      edge_dxdy(1,edge)*a*edge_up_dn_grad(2,nz,edge)+ &
      edge_dxdy(2,edge)*r_earth*edge_up_dn_grad(4,nz,edge))/6.0_8    
  
   Tmean1=ttfold(nz, enodes(1))+ &
      (2.0_8*(ttfold(nz, enodes(2))-ttfold(nz,enodes(1)))+ &
      edge_dxdy(1,edge)*a*edge_up_dn_grad(1,nz,edge)+ &
      edge_dxdy(2,edge)*r_earth*edge_up_dn_grad(3,nz,edge))/6.0_8    
   ! ------
   ! volume flux across the segments
   ! ------
   vflux=UV(2,nz,el(2))*deltaX2 - UV(1,nz,el(2))*deltaY2
   ! ------
   ! tracer flux upwind
   ! ------
   c1=(vflux+abs(vflux))*Tmean1+(vflux-abs(vflux))*Tmean2
   ! ------
   ! combined with centered
   ! ------
   c1=-0.5_8*((1.0_8-num_ord)*c1+vflux*num_ord*(Tmean1+Tmean2))
   ttrhs(nz,enodes(1))= ttrhs(nz,enodes(1))+c1
   ttrhs(nz,enodes(2))= ttrhs(nz,enodes(2))-c1
   END DO
   end if
 END DO

! ===================
! Vertical advection
! ===================
 DO n=1, myDim_nod2D
  ! ===========
  ! Fluxes in the column
  ! ===========
 
   tvert(1)= -Wvel_e(1,n)*ttfold(1,n)*area(1,n)
		    
  ! Bottom conditions	  
   tvert(nlevels_nod2D(n))=0.	    
  
   DO nz=2, nlevels_nod2D(n)-1
      ! ============
      ! QUICK upwind (3rd order)
      ! ============
      if(Wvel_e(nz,n)>0) then
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

      if(Wvel_e(nz,n)<0) then
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
      tvert(nz)= -Tmean*Wvel_e(nz,n)*area(nz,n)
   END DO
 
   DO nz=1,nlevels_nod2D(n)-1
      ttrhs(nz,n)=(ttrhs(nz,n)+ &
                    (tvert(nz)-tvert(nz+1))/(zbar(nz)-zbar(nz+1)))
   END DO
 END DO
! =================
! Update dttf (will be also used on the next time level)
! and compute new ttf
! =================            
  DO n=1, myDim_nod2D+eDim_nod2D
     DO nz=1,nlevels_nod2D(n)-1
        dttf(nz,n)=dttf(nz,n)+ttrhs(nz,n)*dt/area(nz,n)
     END DO
  END DO
end subroutine adv_tracer_muscl34

!===========================================================================
