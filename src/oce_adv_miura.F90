SUBROUTINE adv_tracer_miura(ttf, dttf, tr_num,tt_xy,tt_xynodes)
use g_parsup
use o_PARAM
use o_arrays
use o_mesh
use g_config
IMPLICIT NONE
real(kind=WP) :: tsv(nl,myDim_nod2D+eDIm_nod2D),tvol 
real(kind=WP) :: ttf(nl-1, myDim_nod2D+eDim_nod2D), dttf(nl-1, myDim_nod2D+eDim_nod2D)
integer       :: tr_num
integer       :: el(2), enodes(2), n, nz, edge
integer       :: nl1, nl2
real(kind=WP) :: tvert(nl), c1, c2, deltaX1, deltaY1, deltaX2, deltaY2
real(kind=WP) :: a, b, c, d, da, db, dg, wm, dbm, dgm
real(kind=WP) :: Tx, Ty, Tmean
real*8      :: tt_xy(2,nl-1,myDim_elem2D+eDim_elem2D+eXDim_elem2D)
real*8      :: tt_xynodes(2,nl-1,myDim_nod2D+eDim_nod2D)
!real*8 t1,t2 !delete me
!t1=MPI_Wtime()
! =================
! Fill the rhs with zeros
! =================
 ttrhs=0d0
! =================
! Horizontal advection
! =================
 DO edge=1, myDim_edge2D    !! m=1, myDim_edge2D       
                            !! edge=myList_edge2D(m)
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
   ! the first column
   ! ============
   DO nz=1, nl1
   ! ============
   ! Miura upwind implementation (second order, linear reconstruction)
   ! ============
   if(UV(2,nz,el(1))*deltaX1- UV(1,nz,el(1))*deltaY1>0) then   
      Tmean=ttf(nz, enodes(2))+0.5*(-edge_dxdy(1,edge)*a+ deltaX1- &
                              UV(1,nz,el(1))*dt)*tt_xynodes(1,nz,enodes(2))+&
                               0.5*(-edge_dxdy(2,edge)*r_earth+ deltaY1- &
			      UV(2,nz,el(1))*dt)*tt_xynodes(2,nz,enodes(2))    
   else
      Tmean=ttf(nz, enodes(1))+0.5*(edge_dxdy(1,edge)*a+ deltaX1 - &
                              UV(1,nz,el(1))*dt)*tt_xynodes(1,nz,enodes(1))+&
                               0.5*(edge_dxdy(2,edge)*r_earth+ deltaY1 - &
			      UV(2,nz,el(1))*dt)*tt_xynodes(2,nz,enodes(1))    
   end if
   c1=(UV(2,nz,el(1))*Tmean)*deltaX1- (UV(1,nz,el(1))*Tmean)*deltaY1
   ttrhs(nz,enodes(1))=ttrhs(nz,enodes(1)) + c1
   ttrhs(nz,enodes(2))=ttrhs(nz,enodes(2)) - c1
   END DO
   if(el(2)>0)then
   ! ============
   ! the second column
   ! ============
   DO nz=1, nl2
   if(UV(2,nz,el(2))*deltaX2- UV(1,nz,el(2))*deltaY2<0) then   
      Tmean=ttf(nz, enodes(2))+0.5*(-edge_dxdy(1,edge)*b+ deltaX2- &
                                 UV(1,nz,el(2))*dt)*tt_xynodes(1,nz,enodes(2))+&
                               0.5*(-edge_dxdy(2,edge)*r_earth+ deltaY2- &
			         UV(2,nz,el(2))*dt)*tt_xynodes(2,nz,enodes(2))    
   else
      Tmean=ttf(nz, enodes(1))+0.5*(edge_dxdy(1,edge)*b+ deltaX2- &
                                 UV(1,nz,el(2))*dt)*tt_xynodes(1,nz,enodes(1))+&
                               0.5*(edge_dxdy(2,edge)*r_earth+ deltaY2- &
			         UV(2,nz,el(2))*dt)*tt_xynodes(2,nz,enodes(1))    
   end if
   c2=-(UV(2,nz,el(2))*Tmean)*deltaX2 + (UV(1,nz,el(2))*Tmean)*deltaY2
   ttrhs(nz,enodes(1))=ttrhs(nz,enodes(1))+c2
   ttrhs(nz,enodes(2))=ttrhs(nz,enodes(2))-c2
   END DO
   endif
 END DO
! ===================
! Vertical advection
! ===================
   DO n=1,myDim_nod2D+eDim_nod2D
     DO nz=2,nlevels_nod2D(n)-1
        tvol=Z(nz-1)-Z(nz)
        tsv(nz,n)=(ttf(nz-1,n)-ttf(nz,n))/tvol
     END DO
        tsv(1,n)=0.0_WP
	nz=nlevels_nod2D(n)
	tsv(nz,n)=0.0_WP
  END DO

 DO n=1, myDim_nod2D !! m=1, myDim_nod2D
                     !! n=myList_nod2D(m)
  ! Surface forcing 
  tvert(1) = -Wvel(1,n)*ttf(1,n)*area(1,n)
  ! Bottom conditions	  
  tvert(nlevels_nod2D(n))=0.	    
  
  DO nz=2, nlevels_nod2D(n)-1
      ! ============
      ! QUICK upwind (3rd order)
      ! ============

      if(Wvel(nz,n)>0) then
        if (nz==nlevels_nod2D(n)-1) then
	  Tmean=0.5_8*(ttf(nz-1,n)+ttf(nz,n))  ! or replace this with first
	                                       ! order upwind  ttf(nz,n)
	else
	a=Z(nz-1)-zbar(nz)
	b=zbar(nz)-Z(nz)
	c=zbar(nz)-Z(nz+1)
	wm=Wvel(nz,n)*dt*.5
	dgm=(b-a-wm)*wm
	dbm=(a-c+wm)*wm
	dg=(a*b+dgm)/(c+a)/(b-c)
	db=(-a*c+dbm)/(b+a)/(b-c)
	da=1.0-dg-db
	
	Tmean=ttf(nz-1,n)*da+ttf(nz,n)*db+ttf(nz+1,n)*dg
	end if
      end if

      if(Wvel(nz,n)<0) then
        if(nz==2) then
	  Tmean=0.5_8*(ttf(nz-1,n)+ttf(nz,n))  
	else  
	  a=zbar(nz)-Z(nz)
	  b=Z(nz-1)-zbar(nz)
	  c=Z(nz-2)-zbar(nz)
	  wm=-Wvel(nz,n)*dt*.5
	  dgm=(b-a-wm)*wm
	  dbm=(a-c+wm)*wm
	  dg=(a*b+dgm)/(c+a)/(b-c)
	  db=(-a*c+dgm)/(b+a)/(b-c)
	  da=1.0-dg-db
	  Tmean=ttf(nz,n)*da+ttf(nz-1,n)*db+ttf(nz-2,n)*dg
	end if
      end if
      tvert(nz)= -Tmean*Wvel(nz,n)*area(nz,n)
   END DO
!Update rhs 
   DO nz=1,nlevels_nod2D(n)-1
      ttrhs(nz,n)=(ttrhs(nz,n)+ &
                    (tvert(nz)-tvert(nz+1))/(zbar(nz)-zbar(nz+1)))
   END DO
END DO
!Update sol
DO n=1, myDim_nod2D
     DO nz=1,nlevels_nod2D(n)-1
        dttf(nz,n)=dttf(nz,n)+ttrhs(nz,n)*dt/area(nz,n)
     END DO
END DO
end subroutine adv_tracer_miura
