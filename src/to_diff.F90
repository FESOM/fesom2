

!Contains:
!	adv_tracer_fct
!	fct_init
!	fct_hor
!	fct_ver
!	fct_muscl_solve_LH
! ==========================================================================
SUBROUTINE adv_tracer_fct(ttf,dttf,ttfold,num_ord)
USE o_ARRAYS
USE o_MESH
USE g_PARSUP
IMPLICIT NONE
real(kind=WP)  :: ttf(nl-1, myDim_nod2D+eDim_nod2D),dttf(nl-1, myDim_nod2D+eDim_nod2D)
real(kind=WP)  :: ttfold(nl-1, myDim_nod2D+eDim_nod2D)
real(kind=WP)  :: num_ord
 call fct_muscl_LH(ttf, ttfold,num_ord)
 call fct(ttf, dttf)

END SUBROUTINE adv_tracer_fct
! ==========================================================================
SUBROUTINE fct_init
USE o_MESH
USE o_ARRAYS
USE o_PARAM
USE g_PARSUP
IMPLICIT NONE
integer      :: my_size

my_size=myDim_nod2D+eDim_nod2D
allocate(fct_aec(nl-1,myDim_edge2D))   ! Antidiffusive edge contributions
allocate(fct_LO(nl-1, my_size))        ! Low-order solution 
allocate(fct_aec_ver(nl, myDim_nod2D)) ! antidiffusive 
                                       ! vertical fluxes
allocate(fct_ttf_max(nl-1, my_size),fct_ttf_min(nl-1, my_size))
allocate(fct_plus(nl-1, my_size),fct_minus(nl-1, my_size))
! Initialize with zeros: 
 fct_aec=0.0
 fct_LO=0.0
 
 fct_aec_ver=0.0
 fct_ttf_max=0.0
 fct_ttf_min=0.0
 fct_plus=0.0
 fct_minus=0.0
 
 write(*,*) 'FCT is initialized'
 
END SUBROUTINE fct_init
!===========================================================================
SUBROUTINE fct(ttf, dttf)
!
! 3D Flux Corrected Transport scheme
! Limits antidiffusive fluxes==the difference in flux HO-LO
! LO ==Low-order  (first-order upwind)
! HO ==High-order (3rd/4th order gradient reconstruction method)
! Adds limited fluxes to the LO solution   
USE o_MESH
USE o_ARRAYS
USE o_PARAM
USE g_PARSUP

IMPLICIT NONE
integer       :: n, nz, k, elem, enodes(3), num, el(2), nl1, nl2, edge
real(kind=WP) :: flux, ae,tvert_max(nl-1),tvert_min(nl-1) 
real(kind=WP) :: ttf(nl-1, myDim_nod2D+eDim_nod2D)
real(kind=WP) :: dttf(nl-1, myDim_nod2D+eDim_nod2D)
real*8        :: flux_eps=1e-16
real*8        :: bignumber=1e3
integer       :: vlimit=1
! ----------------------------------
! ttf is the tracer field on step n
! dttf is the increment 
! vlimit sets the version of limiting, see below
! ----------------------------------
! ==========
! a1. max, min between old solution and updated low-order solution 
! ==========
DO n=1,myDim_nod2D + edim_nod2d
   DO nz=1, nlevels_nod2D(n)-1 
      fct_ttf_max(nz,n)=max(fct_LO(nz,n), ttf(nz,n))
      fct_ttf_min(nz,n)=min(fct_LO(nz,n), ttf(nz,n))
   END DO
END DO       

! ==========
! a2. Admissible increments on elements
!     (only layers below the first and above the last layer)
! ==========
DO elem=1, myDim_elem2D
   enodes=elem2D_nodes(:,elem)
   DO nz=1, nlevels(elem)-1
      U_rhs(nz,elem)=maxval(fct_ttf_max(nz,enodes))
      V_rhs(nz,elem)=minval(fct_ttf_min(nz,enodes))
   END DO
if (nlevels(elem)<=nl-1) then
   DO nz=nlevels(elem),nl-1
      U_rhs(nz,elem)=-bignumber
      V_rhs(nz,elem)= bignumber
   END DO
endif
END DO
! ==========
! a3. Bounds on clusters and admissible increments
! ==========
if(vlimit==1) then
!Horizontal
DO n=1, myDim_nod2D
   DO nz=1,nlevels_nod2D(n)-1
      tvert_max(nz)= maxval(U_rhs(nz,nod_in_elem2D(1:nod_in_elem2D_num(n),n)))
      tvert_min(nz)= minval(V_rhs(nz,nod_in_elem2D(1:nod_in_elem2D_num(n),n)))
   END DO
   !fct_ttf_max(:,n)=tvert_max
   !fct_ttf_min(:,n)=tvert_min 
   fct_ttf_max(1,n)=tvert_max(1)-fct_LO(1,n)
   fct_ttf_min(1,n)=tvert_min(1)-fct_LO(1,n)
   DO nz=2,nlevels_nod2D(n)-2  
   fct_ttf_max(nz,n)=maxval(tvert_max(nz-1:nz+1))-fct_LO(nz,n)
   fct_ttf_min(nz,n)=minval(tvert_min(nz-1:nz+1))-fct_LO(nz,n)
   END DO
   nz=nlevels_nod2D(n)-1
   fct_ttf_max(nz,n)=tvert_max(nz)-fct_LO(nz,n)
   fct_ttf_min(nz,n)=tvert_min(nz)-fct_LO(nz,n)  
END DO
!Vertical1: In this version we look at the bounds on the clusters
!           above and below, which leaves wide bounds because typically 
!           vertical gradients are larger.  
end if

if(vlimit==2) then
DO n=1, myDim_nod2D
   DO nz=1,nlevels_nod2D(n)-1
      tvert_max(nz)= maxval(U_rhs(nz,nod_in_elem2D(1:nod_in_elem2D_num(n),n)))
      tvert_min(nz)= minval(V_rhs(nz,nod_in_elem2D(1:nod_in_elem2D_num(n),n)))
   END DO
   DO nz=2, nlevels_nod2D(n)-2
      tvert_max(nz)=max(tvert_max(nz),maxval(fct_ttf_max(nz-1:nz+1,n)))
      tvert_min(nz)=min(tvert_min(nz),minval(fct_ttf_max(nz-1:nz+1,n)))
   END DO
   DO nz=1,nlevels_nod2D(n)-1
       fct_ttf_max(nz,n)=tvert_max(nz)-fct_LO(nz,n)
       fct_ttf_min(nz,n)=tvert_min(nz)-fct_LO(nz,n)  
   END DO
END DO
!Vertical2: Similar to the version above, but the vertical bounds are more local  
END IF

if(vlimit==3) then
DO n=1, myDim_nod2D
   DO nz=1,nlevels_nod2D(n)-1
      tvert_max(nz)= maxval(U_rhs(nz,nod_in_elem2D(1:nod_in_elem2D_num(n),n)))
      tvert_min(nz)= minval(V_rhs(nz,nod_in_elem2D(1:nod_in_elem2D_num(n),n)))
   END DO
   DO nz=2, nlevels_nod2D(n)-2
      tvert_max(nz)=min(tvert_max(nz),maxval(fct_ttf_max(nz-1:nz+1,n)))
      tvert_min(nz)=max(tvert_min(nz),minval(fct_ttf_max(nz-1:nz+1,n)))
   END DO
   DO nz=1,nlevels_nod2D(n)-1
       fct_ttf_max(nz,n)=tvert_max(nz)-fct_LO(nz,n)
       fct_ttf_min(nz,n)=tvert_min(nz)-fct_LO(nz,n)  
   END DO
END DO
!Vertical3: Vertical bounds are taken into account only if they are narrower than the
!           horizontal ones  
END IF

! ==========
! b1. Split positive and negative antidiffusive contributions
! ==========
DO n=1, myDim_nod2D
   DO nz=1,nlevels_nod2D(n)-1
      fct_plus(nz,n)=0.
      fct_minus(nz,n)=0.
   END DO
END DO
!Vertical
DO n=1, myDim_nod2D
   DO nz=1,nlevels_nod2D(n)-1
      fct_plus(nz,n)=fct_plus(nz,n)+ &
                    (max(0.0_8,fct_aec_ver(nz,n))+max(0.0_8,-fct_aec_ver(nz+1,n))) &
                     /(zbar(nz)-zbar(nz+1))
      fct_minus(nz,n)=fct_minus(nz,n)+ &
                     (min(0.0_8,fct_aec_ver(nz,n))+min(0.0_8,-fct_aec_ver(nz+1,n))) &
                     /(zbar(nz)-zbar(nz+1))
   END DO
END DO
!Horizontal
DO edge=1, myDim_edge2D
   enodes(1:2)=edges(:,edge)   
   el=edge_tri(:,edge)
   nl1=nlevels(el(1))-1
   nl2=0
   if(el(2)>0) then
     nl2=nlevels(el(2))-1
   end if   
   DO nz=1, max(nl1,nl2)
      fct_plus (nz,enodes(1))=fct_plus (nz,enodes(1)) + max(0.0, fct_aec(nz,edge))
      fct_minus(nz,enodes(1))=fct_minus(nz,enodes(1)) + min(0.0, fct_aec(nz,edge))  
      fct_plus (nz,enodes(2))=fct_plus (nz,enodes(2)) + max(0.0,-fct_aec(nz,edge))
      fct_minus(nz,enodes(2))=fct_minus(nz,enodes(2)) + min(0.0,-fct_aec(nz,edge)) 
   END DO   
END DO 
! ==========
! b2. Limiting factors
! ==========
DO n=1,myDim_nod2D
   DO nz=1,nlevels_nod2D(n)-1
     flux=fct_plus(nz,n)*dt/area(nz,n)+flux_eps
     fct_plus(nz,n)=min(1.0,fct_ttf_max(nz,n)/flux)
     flux=fct_minus(nz,n)*dt/area(nz,n)-flux_eps
     fct_minus(nz,n)=min(1.0,fct_ttf_min(nz,n)/flux)
   END DO
END DO 
! fct_minus and fct_plus must be known to neighbouring PE
  call exchange_nod3D(fct_plus)
  call exchange_nod3D(fct_minus)	   
! 
!========================	 
! b3. Limiting   
!========================	 
!Vertical
DO n=1, myDim_nod2D
   nz=1
   ae=1.0_8
   flux=fct_aec_ver(nz,n)
   if(flux>=0.0) then 
      ae=min(ae,fct_plus(nz,n))
   else
      ae=min(ae,fct_minus(nz,n))
   end if
   fct_aec_ver(nz,n)=ae*fct_aec_ver(nz,n)    
   DO nz=2,nlevels_nod2D(n)-1
      ae=1.0
      flux=fct_aec_ver(nz,n)
        if(flux>=0.) then 
	   ae=min(ae,fct_minus(nz-1,n))
	   ae=min(ae,fct_plus(nz,n))
	else
	   ae=min(ae,fct_plus(nz-1,n))
	   ae=min(ae,fct_minus(nz,n))
	end if   
      fct_aec_ver(nz,n)=ae*fct_aec_ver(nz,n)
   END DO
   ! the bottom flux is always zero 
END DO
!Horizontal
DO edge=1, myDim_edge2D
   enodes(1:2)=edges(:,edge)
   el=edge_tri(:,edge)
   nl1=nlevels(el(1))-1
   
   nl2=0
   if(el(2)>0) then
      nl2=nlevels(el(2))-1
   end if  
   DO nz=1, max(nl1,nl2)
      ae=1.0
      flux=fct_aec(nz,edge)
      if(flux>=0.) then
          ae=min(ae,fct_plus(nz,enodes(1)))
          ae=min(ae,fct_minus(nz,enodes(2)))
      else
          ae=min(ae,fct_minus(nz,enodes(1)))
          ae=min(ae,fct_plus(nz,enodes(2)))
      endif
      fct_aec(nz,edge)=ae*fct_aec(nz,edge)
   END DO
END DO
!==========================
! c. Update the solution
!==========================

   ! Vertical part
DO n=1, myDim_nod2d
   DO nz=1,nlevels_nod2D(n)-1  
   dttf(nz,n)=dttf(nz,n)-ttf(nz,n)+fct_LO(nz,n)+ dt* &
             (fct_aec_ver(nz,n)-fct_aec_ver(nz+1,n))/(zbar(nz)-zbar(nz+1))/area(nz,n)    
   END DO
END DO
   ! Horizontal part
DO edge=1, myDim_edge2D
   enodes(1:2)=edges(:,edge)
   el=edge_tri(:,edge)
   nl1=nlevels(el(1))-1
   nl2=0
   if(el(2)>0) nl2=nlevels(el(2))-1
  DO nz=1, max(nl1,nl2)
      dttf(nz,enodes(1))=dttf(nz,enodes(1))+fct_aec(nz,edge)*dt/area(nz,enodes(1))
      dttf(nz,enodes(2))=dttf(nz,enodes(2))-fct_aec(nz,edge)*dt/area(nz,enodes(2))
   END DO
END DO  
END SUBROUTINE fct
!===========================================================================
SUBROUTINE fct_muscl_LH(ttf, ttfold, num_ord)
USE o_MESH
USE o_ARRAYS
USE o_PARAM
USE g_PARSUP
IMPLICIT NONE
 integer       :: el(2), enodes(2), n, nz, edge
 integer       :: nl1, nl2, n2
 real(kind=WP) :: c1, deltaX1, deltaY1, deltaX2, deltaY2 
 real(kind=WP) :: tvert(nl), a, qc, qu, qd, un1, un2, Tupw1
 real(kind=WP) :: Tmean1, Tmean2, Tmean, num_ord
 real(kind=WP) :: ttf(nl-1,myDim_nod2D+eDim_nod2D), ttfold(nl-1,myDim_nod2D+eDim_nod2D)
 ! ----------------------------------------
 ! It is assumed that velocity is at n+1/2, hence only tracer field 
 ! is AB2 interpolated to n+1/2. 
 ! ttf contains tracer at step n
 ! ttfold is AB interpolated (n+1/2)
 ! num_ord is the fraction of fourth-order contribution in the HO solution
 ! The result is the low-order solution
 ! and vertical and horizontal antidiffusive fluxes
 ! They are put into fct_LO
 !                   fct_aec
 !                   fct_aec_ver
 ! ----------------------------------------
 ! =================
 ! Clean the low-order solution
 ! =================          
 fct_LO=0d0 
 ! =================
 ! Horizontal advection
 ! =================
  DO edge=1, myDim_edge2D
   enodes=edges(:,edge)   
   el=edge_tri(:,edge)
   nl1=nlevels(el(1))-1
   deltaX1=edge_cross_dxdy(1,edge)
   deltaY1=edge_cross_dxdy(2,edge)
   a=r_earth*elem_cos(el(1))
   nl2=0
   if(el(2)>0) then
     deltaX2=edge_cross_dxdy(3,edge)
     deltaY2=edge_cross_dxdy(4,edge)
     nl2=nlevels(el(2))-1
     a=0.5_WP*(a+r_earth*elem_cos(el(2)))
   end if     
   n2=min(nl1,nl2)
   DO nz=1, n2
   ! ============
   ! MUSCL-type recoonstruction (high order) and upwind (low order)
   ! ============
      Tmean2=ttfold(nz, enodes(2))- &
      (2d0*(ttfold(nz, enodes(2))-ttfold(nz,enodes(1)))+ &
      edge_dxdy(1,edge)*a*edge_up_dn_grad(2,nz,edge)+ &
      edge_dxdy(2,edge)*r_earth*edge_up_dn_grad(4,nz,edge))/6d0  
   
      Tmean1=ttfold(nz, enodes(1))+ &
      (2.0_8*(ttfold(nz, enodes(2))-ttfold(nz,enodes(1)))+ &
      edge_dxdy(1,edge)*a*edge_up_dn_grad(1,nz,edge)+ &
      edge_dxdy(2,edge)*r_earth*edge_up_dn_grad(3,nz,edge))/6.0_8    
   
   un1=  -V_n(nz,el(1))*deltaX1+ U_n(nz,el(1))*deltaY1   ! normal outer velocity
   un2=   V_n(nz,el(2))*deltaX2- U_n(nz,el(2))*deltaY2   ! from 1 to 2
   un1=un1+un2
   Tupw1=0.5_8*(ttf(nz, enodes(1))*(un1+abs(un1))+ttf(nz, enodes(2))*(un1-abs(un1)))
   c1=-Tupw1 
   
   fct_LO(nz,enodes(1))=fct_LO(nz,enodes(1))+c1            ! Low-order upwind 
   fct_LO(nz,enodes(2))=fct_LO(nz,enodes(2))-c1            !
   c1=(un1+abs(un1))*Tmean1+(un1-abs(un1))*Tmean2
   c1=-0.5_8*((1.0_8-num_ord)*c1+un1*num_ord*(Tmean1+Tmean2))
   fct_aec(nz,edge)=c1+Tupw1                               ! Antidiffusive edge flux
   END DO
   if(nl1>nl2) then
   DO nz=n2+1, nl1
       Tmean2=ttfold(nz, enodes(2))- &
      (2.0_8*(ttfold(nz, enodes(2))-ttfold(nz,enodes(1)))+ &
      edge_dxdy(1,edge)*a*edge_up_dn_grad(2,nz,edge)+ &
      edge_dxdy(2,edge)*r_earth*edge_up_dn_grad(4,nz,edge))/6.0_8    
   
      Tmean1=ttfold(nz, enodes(1))+ &
      (2.0_8*(ttfold(nz, enodes(2))-ttfold(nz,enodes(1)))+ &
      edge_dxdy(1,edge)*a*edge_up_dn_grad(1,nz,edge)+ &
      edge_dxdy(2,edge)*r_earth*edge_up_dn_grad(3,nz,edge))/6.0_8    
   
      un1=  -V_n(nz,el(1))*deltaX1+ U_n(nz,el(1))*deltaY1   ! normal outer velocity
                                                              ! from 1 to 2
      Tupw1=0.5_8*(ttf(nz, enodes(1))*(un1+abs(un1))+ttf(nz, enodes(2))*(un1-abs(un1)))
      c1=-Tupw1
      fct_LO(nz,enodes(1))=fct_LO(nz,enodes(1))+c1 
      fct_LO(nz,enodes(2))=fct_LO(nz,enodes(2))-c1
      c1=(un1+abs(un1))*Tmean1+(un1-abs(un1))*Tmean2
      c1=-0.5_8*((1.0_8-num_ord)*c1+un1*num_ord*(Tmean1+Tmean2))
      fct_aec(nz,edge)=c1+Tupw1                               ! Antidiffusive edge flux 
   END DO
   else
   DO nz=n2+1, nl2
      Tmean2=ttfold(nz, enodes(2))- &
      (2.0_8*(ttfold(nz, enodes(2))-ttfold(nz,enodes(1)))+ &
      edge_dxdy(1,edge)*a*edge_up_dn_grad(2,nz,edge)+ &
      edge_dxdy(2,edge)*r_earth*edge_up_dn_grad(4,nz,edge))/6.0_8    
   
      Tmean1=ttfold(nz, enodes(1))+ &
      (2.0_8*(ttfold(nz, enodes(2))-ttfold(nz,enodes(1)))+ &
      edge_dxdy(1,edge)*a*edge_up_dn_grad(1,nz,edge)+ &
      edge_dxdy(2,edge)*r_earth*edge_up_dn_grad(3,nz,edge))/6.0_8    
                                                              ! normal outer velocity
      un1=   V_n(nz,el(2))*deltaX2- U_n(nz,el(2))*deltaY2   ! from 1 to 2
      Tupw1=0.5_8*(ttf(nz, enodes(1))*(un1+abs(un1))+ttf(nz, enodes(2))*(un1-abs(un1)))
      c1=-Tupw1
      fct_LO(nz,enodes(1))=fct_LO(nz,enodes(1))+c1 
      fct_LO(nz,enodes(2))=fct_LO(nz,enodes(2))-c1
      c1=(un1+abs(un1))*Tmean1+(un1-abs(un1))*Tmean2
      c1=-0.5_8*((1.0_8-num_ord)*c1+un1*num_ord*(Tmean1+Tmean2))
      fct_aec(nz,edge)=c1+Tupw1                               ! Antidiffusive edge flux 
   END DO
   end if
 END DO

! ===================
! Vertical advection
! ===================
DO n=1, myDim_nod2D
    nl1=nlevels_nod2D(n)
    ! -----
    ! upper level --- linear free surface
    ! -----
    tvert(1)= -Wvel(1,n)*ttf(1,n)*area(1,n)
    ! -----
    ! No flux through the bottom
    ! -----
    tvert(nl1)=0.0_8;
    ! ------
    ! Centered differences for the second and 
    ! last but one
    ! ------ 
    nz=2  
    Tupw1=0.5_8*(ttf(nz,n)*(Wvel(nz,n)+abs(Wvel(nz,n)))+ &
                 ttf(nz-1,n)*(Wvel(nz,n)-abs(Wvel(nz,n))))
    Tmean=0.5_8*(ttfold(nz-1,n)+ttfold(nz,n))*Wvel(nz,n) 
    tvert(nz)= -Tupw1*area(nz,n)
    fct_aec_ver(nz,n)=(-Tmean+Tupw1)*area(nz,n)
    nz=nl1-1
    Tupw1=0.5_8*(ttf(nz,n)*(Wvel(nz,n)+abs(Wvel(nz,n)))+ &
                 ttf(nz-1,n)*(Wvel(nz,n)-abs(Wvel(nz,n))))
    Tmean=0.5_8*(ttfold(nz-1,n)+ttfold(nz,n))*Wvel(nz,n) 
    tvert(nz)= -Tupw1*area(nz,n)
    fct_aec_ver(nz,n)=(-Tmean+Tupw1)*area(nz,n)
    ! the remaining levels    
    DO nz=3,nl1-2
    ! ------
    ! First-order upwind estimate
    ! ------
    Tupw1=0.5*(ttf(nz,n)*(Wvel(nz,n)+abs(Wvel(nz,n)))+ &
              ttf(nz-1,n)*(Wvel(nz,n)-abs(Wvel(nz,n))))
    ! -----
    ! centered (4th order)
    ! -----
    qc=(ttfold(nz-1,n)-ttfold(nz,n))/(Z(nz-1)-Z(nz))
    qu=(ttfold(nz,n)-ttfold(nz+1,n))/(Z(nz)-Z(nz+1))    
    qd=(ttfold(nz-2,n)-ttfold(nz-1,n))/(Z(nz-2)-Z(nz-1))
    Tmean1=ttfold(nz,n)+(2*qc+qu)*(zbar(nz)-Z(nz))/3.0_8
    Tmean2=ttfold(nz-1,n)+(2*qc+qd)*(zbar(nz)-Z(nz-1))/3.0_8
    Tmean=(Wvel(nz,n)+abs(Wvel(nz,n)))*Tmean1+(Wvel(nz,n)-abs(Wvel(nz,n)))*Tmean2
    tvert(nz)=-Tupw1*area(nz,n)
    fct_aec_ver(nz,n)=(-0.5_8*(num_ord*(Tmean1+Tmean2)*Wvel(nz,n)+(1.0_8-num_ord)*Tmean)+Tupw1)*area(nz,n)
    end do
    fct_aec_ver(1,n)=0.0_8
    fct_aec_ver(nl1,n)=0.0_8
   DO nz=1,nlevels_nod2D(n)-1
      fct_LO(nz,n)=ttf(nz,n)+dt*(fct_LO(nz,n)+ &
                      (tvert(nz)-tvert(nz+1))/(zbar(nz)-zbar(nz+1)))/area(nz,n)
   END DO
 ENDDO
 call exchange_nod3D(fct_LO) 
 ! Summary:   
 ! fct_LO contains full low-order solution
 ! fct_aec contains antidiffusive component of horizontal flux 
 ! fct_aec_ver contains antidiffusive component of vertical fluxes
end subroutine fct_muscl_LH
!===========================================================================.


