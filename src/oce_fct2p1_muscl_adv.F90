!Contains:
!	adv_tracer_fct34
!	fct_init
!	fct_hor
!	fct_ver
!	fct_muscl34_solve_LH
! ==========================================================================
SUBROUTINE adv_tracer_fct2p1_34(ttf,dttf,ttfold,tr_num,tt_xy,tt_xynodes)
USE o_ARRAYS
USE o_MESH
USE g_PARSUP
IMPLICIT NONE
real(kind=WP) :: ttf(nl-1, myDim_nod2D+eDim_nod2D),dttf(nl-1, myDim_nod2D+eDim_nod2D)
real(kind=WP) :: ttfold(nl-1, myDim_nod2D+eDim_nod2D)
real*8      :: tt_xy(2,nl-1,myDim_elem2D+eDim_elem2D+eXDim_elem2D)
real*8      :: tt_xynodes(2,nl-1,myDim_nod2D+eDim_nod2D)
integer tr_num
call fct_muscl_LH(ttf, ttfold, tr_num, 0.5)
call fct_hor(ttf,dttf)
call fct_ver(ttf, dttf,ttfold, tr_num)

END SUBROUTINE adv_tracer_fct2p1_34
!===========================================================================
SUBROUTINE fct_hor(ttf,dttf)

! Update of low-order (Upwind) to high-order (MUSCL) with account for horizontal 
! fluxes. 

USE o_MESH
USE o_ARRAYS
USE o_PARAM
USE g_PARSUP
USE g_CONFIG
USE g_comm_auto
IMPLICIT NONE
integer       ::  n, nz, k, elem, enodes(3), nl1, nl2, edge, el(2)
real(kind=WP) :: flux, ae 
real(kind=WP) :: ttf(nl-1, myDim_nod2D+eDim_nod2D),dttf(nl-1, myDim_nod2D+eDim_nod2D)
! ==========
! a1. max, min between old solution and updated low-order solution 
!    (we do not do this screening in FEOM's FCT)
! ==========

DO n=1,myDim_nod2D+eDim_nod2D
                                 !! n=myList_nod2D(m)
   DO nz=1, nlevels_nod2D(n)-1 
      fct_ttf_max(nz,n)=max(fct_LO(nz,n), ttf(nz,n))
      fct_ttf_min(nz,n)=min(fct_LO(nz,n), ttf(nz,n))
   END DO
END DO       

! ==========
! a2. Bounds on elements
! (use U_rhs and V-rhs as placeholder for temporary storage)
! ==========

DO elem=1, myDim_elem2D
                                 !!elem=myList_elem2D(m)
   enodes=elem2D_nodes(:,elem)
   DO nz=1, nlevels(elem)-1
      UV_rhs(1,nz,elem)=maxval(fct_ttf_max(nz,enodes))
      UV_rhs(2,nz,elem)=minval(fct_ttf_min(nz,enodes))
   END DO
if (nlevels(elem)<=nl-1) then
   DO nz=nlevels(elem),nl-1
      UV_rhs(1,nz,elem)=-1000.0
      UV_rhs(2,nz,elem)=1000.0
   END DO
endif
END DO

! ==========
! a3. Bounds on clusters and admissible increments
! ==========
DO n=1, myDim_nod2D
                                   !! n=myList_nod2D(m)     
   DO nz=1,nlevels_nod2D(n)-1
      fct_ttf_max(nz,n)= &
          maxval(UV_rhs(1,nz,nod_in_elem2D(1:nod_in_elem2D_num(n),n)))-fct_LO(nz,n)
      fct_ttf_min(nz,n)= &
          minval(UV_rhs(2,nz,nod_in_elem2D(1:nod_in_elem2D_num(n),n)))-fct_LO(nz,n)
   END DO
END DO

! ==========
! b1. Split positive and negative contributions from edges
! ==========

DO n=1, myDim_nod2D
                                   !! n=myList_nod2D(m)     
   DO nz=1,nlevels_nod2D(n)-1
      fct_plus(nz,n)=0.
      fct_minus(nz,n)=0.
   END DO
END DO

DO edge=1, myDim_edge2D                   !! P (c) edge=1, edge2d
                                   !! edge=myList_edge2D(m)
   enodes(1:2)=edges(:,edge)   
   el=edge_tri(:,edge)
   nl1=nlevels(el(1))-1
   
   nl2=0
   if(el(2)>0) then
   nl2=nlevels(el(2))-1
   end if   
   DO nz=1, max(nl1,nl2)
      fct_plus(nz,enodes(1))=fct_plus(nz,enodes(1))+max(0.0_8,fct_aec(nz,edge))
      fct_minus(nz,enodes(1))=fct_minus(nz,enodes(1))+min(0.0_8,fct_aec(nz,edge))   
      fct_plus(nz,enodes(2))=fct_plus(nz,enodes(2))+max(0.0_8,-fct_aec(nz,edge))
      fct_minus(nz,enodes(2))=fct_minus(nz,enodes(2))+min(0.0_8,-fct_aec(nz,edge))   
      ! They have to be multiplied with dt and divided by area(nz,n) to get
      ! positive and negative contributions. We do it on the next step.
   END DO   
END DO 

! ==========
! b2. Limiting factors
! ==========
DO n=1,myDim_nod2D
                                               !! n=myList_nod2D(m)
   DO nz=1,nlevels_nod2D(n)-1
     flux=fct_plus(nz,n)*dt/area(nz,n)
     if (abs(flux)>0.) then    
        fct_plus(nz,n)=min(1.0_8,fct_ttf_max(nz,n)/flux)
     else
        fct_plus(nz,n)=0.
     end if
     flux=fct_minus(nz,n)*dt/area(nz,n)
     if (abs(flux)>0.) then     
        fct_minus(nz,n)=min(1.0_8,fct_ttf_min(nz,n)/flux)
     else
        fct_minus(nz,n)=0.
     end if
  END DO
END DO  
  ! fct_minus and fct_plus must be known to neighbouring PE
  call exchange_nod(fct_plus)
  call exchange_nod(fct_minus)	   
 
!========================	 
! b3. Limiting
!========================	 
DO edge=1, myDim_edge2D
                                      !! edge=myList_edge2D(m)
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
	end if
      fct_aec(nz,edge)=ae*fct_aec(nz,edge)
   END DO
END DO


!==========================
! c. Update the solution (see fct_ver)
!==========================
!write(*,*) 'FCT_HOR', maxval(fct_aec), minval(fct_aec)
!=============
! Summary:
! fct_aec contains now limited horizontal antidiffusive fluxes associated with edges.
! Tracer field is not updated here, but will be at the end of vertical part.
!=============
END SUBROUTINE fct_hor
!===========================================================================
SUBROUTINE fct_ver(ttf, dttf, ttfold, tr_num)
! FCT in vertical
! Add vertical antidiffusive fluxes where possible
! 
USE o_MESH
USE o_ARRAYS
USE o_PARAM
USE g_PARSUP
USE g_CONFIG
use g_comm_auto
IMPLICIT NONE
integer       :: n, nz, k, elem, enodes(3), num, el(2), nl1, nl2, edge, tr_num
!real(kind=WP) :: flux, ae,tvert(nl), tvert_(nl) 
real(kind=WP) :: flux, ae,tvert(nl-1) 
real(kind=WP) :: ttf(nl-1, myDim_nod2D+eDim_nod2D),dttf(nl-1, myDim_nod2D+eDim_nod2D)
real(kind=WP) :: ttfold(nl-1, myDim_nod2D+eDim_nod2D)
real*8 :: flux_eps=1d-16

! ==========
! a1. max, min between old solution and updated low-order solution 
!    (we do not do  this screening in FEOM's FCT)
! ==========
DO n=1,myDim_nod2D
                             !! n=myList_nod2D(m)
   DO nz=1, nlevels_nod2D(n)-1 
      fct_ttf_max(nz,n)=max(fct_LO(nz,n), ttf(nz,n))
      fct_ttf_min(nz,n)=min(fct_LO(nz,n), ttf(nz,n))
   END DO
END DO       

! ==========
! a2. Admissible increments on elements
!     (only layers below the first and above the last layer)
! ==========
DO n=1,myDim_nod2D
                             !! n=myList_nod2D(m)
   tvert=fct_ttf_max(:,n)
   DO nz=2, nlevels_nod2D(n)-2 
      fct_ttf_max(nz,n)=maxval(tvert(nz-1:nz+1))-fct_LO(nz,n)
   END DO
   tvert=fct_ttf_min(:,n)
   DO nz=2, nlevels_nod2D(n)-2 
      fct_ttf_min(nz,n)=minval(tvert(nz-1:nz+1))-fct_LO(nz,n)
   END DO
END DO       
! ==========
! b1. Split positive and negative antidiffusive contributions
! ==========
DO n=1, myDim_nod2D
   DO nz=1,nlevels_nod2D(n)-1
      fct_plus(nz,n)=0.
      fct_minus(nz,n)=0.
   END DO
END DO

DO n=1, myDim_nod2D
   DO nz=2,nlevels_nod2D(n)-2
      fct_plus(nz,n)=fct_plus(nz,n)  + &
                     dt*(max(0.0_8,fct_aec_ver(nz,n))+max(0.0_8,-fct_aec_ver(nz+1,n))) &
                     /(zbar(nz)-zbar(nz+1))/area(nz,n)
      fct_minus(nz,n)=fct_minus(nz,n)+ &
                     dt*(min(0.0_8,fct_aec_ver(nz,n))+min(0.0_8,-fct_aec_ver(nz+1,n))) &
                     /(zbar(nz)-zbar(nz+1))/area(nz,n)
   END DO
END DO

! ==========
! b2. Limiting factors
! ==========
DO n=1,myDim_nod2D
   DO nz=2,nlevels_nod2D(n)-2
     flux=fct_plus(nz,n)
     if (abs(flux)>flux_eps) then    
        fct_plus(nz,n)=min(1.0_8,fct_ttf_max(nz,n)/flux)
     else
        fct_plus(nz,n)=0.
     end if
     flux=fct_minus(nz,n)
     if (abs(flux)>flux_eps) then     
        fct_minus(nz,n)=min(1.0_8,fct_ttf_min(nz,n)/flux)
     else
        fct_minus(nz,n)=0.
     end if
  END DO
END DO  
  
!========================	 
! b3. Limiting   
!========================	 
DO n=1, myDim_nod2D
   if(nlevels_nod2D(n)<=3) cycle    ! there are only 1 or 2 layers, no limiting
   nz=2
   ae=1.0
   flux=fct_aec_ver(nz,n)
   if(flux>=0.0) then 
      ae=min(ae,fct_plus(nz,n))
   else
      ae=min(ae,fct_minus(nz,n))
   end if
   fct_aec_ver(nz,n)=ae*fct_aec_ver(nz,n)    
   DO nz=3,nlevels_nod2D(n)-2
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
   nz=nlevels_nod2D(n)-1
   ae=1.0
   flux=fct_aec_ver(nz,n)
   if(flux>=0.0) then 
      ae=min(ae,fct_minus(nz-1,n))
   else
      ae=min(ae,fct_plus(nz-1,n))
   end if
   fct_aec_ver(nz,n)=ae*fct_aec_ver(nz,n)    
END DO
!------------------------------------------------------------------   
DO n=1, myDim_nod2d
!   nz=1
!   dttf(nz,n)=dttf(nz,n)-ttf(nz,n)+fct_LO(nz,n)+ dt* &
!             (fct_aec_ver(nz,n)-fct_aec_ver(nz+1,n))/(zbar(nz)+eta_n(n)-zbar(nz+1))/area(nz,n)    
   DO nz=1,nlevels_nod2D(n)-1  
   dttf(nz,n)=dttf(nz,n)-ttf(nz,n)+fct_LO(nz,n)+ dt* &
             (fct_aec_ver(nz,n)-fct_aec_ver(nz+1,n))/(zbar(nz)-zbar(nz+1))/area(nz,n)    
   END DO
END DO
   ! Horizontal part
DO edge=1, myDim_edge2D
                                      !! edge=myList_edge2D(m)
   enodes(1:2)=edges(:,edge)
   el=edge_tri(:,edge)
   nl1=nlevels(el(1))-1
   
   nl2=0
   if(el(2)>0) then
   nl2=nlevels(el(2))-1
   end if  
   DO nz=1, max(nl1,nl2)
      dttf(nz,enodes(1))=dttf(nz,enodes(1))+fct_aec(nz,edge)*dt/area(nz,enodes(1))
      dttf(nz,enodes(2))=dttf(nz,enodes(2))-fct_aec(nz,edge)*dt/area(nz,enodes(2))
   END DO
END DO  
!=======================================================================
END SUBROUTINE fct_ver
