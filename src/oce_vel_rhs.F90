SUBROUTINE compute_vel_rhs
USE o_MESH
USE o_ARRAYS
USE o_PARAM
USE g_PARSUP
USE g_CONFIG
use g_comm_auto
IMPLICIT NONE
integer          :: elem, elnodes(3), nz 
real(kind=WP)    :: eta(3), ff, gg, mm 
real(kind=WP)    :: Fx, Fy, pre(3)
logical, save    :: lfirst=.true.
real(kind=WP)    :: t1, t2, t3, t4

t1=MPI_Wtime()
! =================
! Take care of the AB part
! =================
do elem=1, myDim_elem2D
   do nz=1,nl-1 
      UV_rhs(1,nz,elem)=-(0.5_WP+epsilon)*UV_rhsAB(1,nz,elem)   
      UV_rhs(2,nz,elem)=-(0.5_WP+epsilon)*UV_rhsAB(2,nz,elem)
    end do
end do
! ====================
! Sea level and pressure contribution   -\nabla(\eta +hpressure/rho_0)
! and the Coriolis force + metric terms
! ====================
DO elem=1, myDim_elem2D
   elnodes=elem2D_nodes(:,elem)
   !eta=g*eta_n(elnodes)*(1-theta)        !! this place needs update (1-theta)!!!
   eta=g*eta_n(elnodes)
   gg=elem_area(elem)
   ff=coriolis(elem)*gg
   !mm=metric_factor(elem)*gg
   DO nz=1,nlevels(elem)-1
      pre=-(eta+hpressure(nz,elnodes)/density_0)
      Fx=sum(gradient_sca(1:3,elem)*pre)
      Fy=sum(gradient_sca(4:6,elem)*pre)
      UV_rhs(1,nz,elem)=UV_rhs(1,nz,elem)+Fx*gg 
      UV_rhs(2,nz,elem)=UV_rhs(2,nz,elem)+Fy*gg 
      UV_rhsAB(1,nz,elem)= UV(2,nz,elem)*ff! + mm*UV(1,nz,elem)*UV(2,nz,elem)
      UV_rhsAB(2,nz,elem)=-UV(1,nz,elem)*ff! - mm*UV(1,nz,elem)*UV(2,nz,elem)
   END DO
END DO
t2=MPI_Wtime() 
! ====================
! Horizontal advection
! ====================
if (mom_adv==1) then
   call momentum_adv_p1
   call momentum_adv_vert
elseif (mom_adv==2) then
   call momentum_adv_scalar
end if
t3=MPI_Wtime() 
! =======================
! Update the rhs   
! =======================
ff=(1.5_WP+epsilon)
if (lfirst.and.(.not.r_restart)) then
   ff=1.0_WP
   lfirst=.false.
end if

DO elem=1, myDim_elem2D
   DO nz=1,nlevels(elem)-1   
      UV_rhs(1,nz,elem)=dt*(UV_rhs(1,nz,elem)+UV_rhsAB(1,nz,elem)*ff)/elem_area(elem)
      UV_rhs(2,nz,elem)=dt*(UV_rhs(2,nz,elem)+UV_rhsAB(2,nz,elem)*ff)/elem_area(elem)
   END DO 
END DO
! =======================  
! U_rhs contains all contributions to velocity from old time steps   
! =======================
t4=MPI_Wtime() 
if (mype==0) then 
   write(*,*) 'Momentum:   ', t4-t1
   write(*,*) 'pres., Cor: ', t2-t1
   write(*,*) 'h adv       ', t3-t2
   write(*,*) 'vert. part  ', t4-t3
end if     
END SUBROUTINE compute_vel_rhs
! ===========================================================================
SUBROUTINE momentum_adv_p1
USE o_MESH
USE o_ARRAYS
USE o_PARAM
USE g_PARSUP
use g_comm_auto
IMPLICIT NONE
real(kind=WP)  :: tvol, tx, ty, uu, vv
integer        :: n, elem, nz, k, elnodes(3), m
integer        :: nl1, nl2, nodes(2), el(2), ed
real(kind=WP)  :: un, u1, v1, xe, ye

! ======  
! Set boundary velocity to zero
! ======  
DO ed=1, myDim_edge2D+eDim_edge2D
   if(myList_edge2D(ed)>edge2D_in) cycle
   nodes=edges(:,ed)   
   el=edge_tri(:,ed)
   nl1=nlevels(el(1))-1
   nl2=nlevels(el(2))-1
   
   xe=edge_dxdy(1,ed)*r_earth*0.5_WP*(elem_cos(el(1))+elem_cos(el(2)))
   ye=edge_dxdy(2,ed)*r_earth
                           ! since velocities are obtained by averaging, 
			   ! one cannot introduce clean metrics here!    
 
   DO nz=1, min(nl1,nl2)   ! Only faces that do not belong to
                           ! vertical walls can contribute to
			   ! the momentum advection 
      !====== 
      ! The piece below gives second order spatial accuracy for
      ! the momentum fluxes. 
      !======
   
      u1=0.5_WP*(Unode(1,nz,nodes(1))+Unode(1,nz,nodes(2)))
      v1=0.5_WP*(Unode(2,nz,nodes(1))+Unode(2,nz,nodes(2)))
   
      !======
      ! Normal velocity at edge ed directed to el(2)
      ! (outer to el(1)) multiplied with the length of the edge
      !======
      un=u1*ye-v1*xe
      UV_rhsAB(1,nz,el(1))=UV_rhsAB(1,nz,el(1))-u1*un
      UV_rhsAB(2,nz,el(1))=UV_rhsAB(2,nz,el(1))-v1*un
      UV_rhsAB(1,nz,el(2))=UV_rhsAB(1,nz,el(2))+u1*un
      UV_rhsAB(2,nz,el(2))=UV_rhsAB(2,nz,el(2))+v1*un
   END DO
END DO  
END subroutine momentum_adv_p1
! ===================================================================
SUBROUTINE momentum_adv_vert
! 
! Vertical momentum advection and viscosity
! For advection, quadratic upwind reconstruction is used.
! 
USE o_PARAM
USE o_MESH
USE o_ARRAYS
USE g_PARSUP
USE g_CONFIG
IMPLICIT NONE

integer          :: nz, elem, elnodes(3)
real(kind=WP)    :: w, uvertAB(2, nl)
real(kind=WP)    :: umean, vmean, a, b, c, d, dg, da, db

! =======================
! Vertical momentum advection 
! (a complementary subroutine for momentum_adv_p1)
! =======================
uvertAB=0.0_WP

DO elem=1, myDim_elem2D
   elnodes=elem2D_nodes(:,elem)
   DO nz=2, nlevels(elem)-1
      w=sum(Wvel(nz,elnodes))*elem_area(elem)/3.0_WP
      umean=0.5_WP*(UV(1,nz-1,elem)+UV(1,nz,elem))
      vmean=0.5_WP*(UV(2,nz-1,elem)+UV(2,nz,elem))
      if (w>0) then
         if (nz==nlevels(elem)-1) then
            umean=0.5_WP*(UV(1,nz-1,elem)+UV(1,nz,elem))
  	    vmean=0.5_WP*(UV(2,nz-1,elem)+UV(2,nz,elem))
	                               ! or replace this with first
	                               ! order upwind  
         else
	    a=Z(nz-1)-zbar(nz)
            b=zbar(nz)-Z(nz)
            c=zbar(nz)-Z(nz+1)
            d=(c+a)*(b**2-a**2)-(c**2-a**2)*(b+a)
            dg=a*b*(a+b)/d
            db=-a*c*(a+c)/d
            da=1.0_WP-dg-db
            umean=UV(1,nz-1,elem)*da+UV(1,nz,elem)*db+UV(1,nz+1,elem)*dg
            vmean=UV(2,nz-1,elem)*da+UV(2,nz,elem)*db+UV(2,nz+1,elem)*dg
         end if
      end if

      if (w<0) then
         if (nz==2) then
            umean=0.5_WP*(UV(1,nz-1,elem)+UV(1,nz,elem))
            vmean=0.5_WP*(UV(2,nz-1,elem)+UV(2,nz,elem))
 	 else  
            a=zbar(nz)-Z(nz)
            b=Z(nz-1)-zbar(nz)
            c=Z(nz-2)-zbar(nz)
            d=(c+a)*(b**2-a**2)-(c**2-a**2)*(b+a)
            dg=a*b*(a+b)/d
            db=-a*c*(a+c)/d
            da=1.0_WP-dg-db
            umean=UV(1,nz,elem)*da+UV(1,nz-1,elem)*db+UV(1,nz-2,elem)*dg
            vmean=UV(2,nz,elem)*da+UV(2,nz-1,elem)*db+UV(2,nz-2,elem)*dg
         end if
      end if
      uvertAB(1,nz)= -umean*w 
      uvertAB(2,nz)= -vmean*w 
   END DO
   w=sum(Wvel(1,elnodes))*elem_area(elem)/3.0_WP 
   uvertAB(1,1)= -w*UV(1,1,elem)
   uvertAB(2,1)= -w*UV(2,1,elem) 
   ! + sign here because it is subtracted!
   uvertAB(1,nlevels(elem))=0.0_WP 
   uvertAB(2,nlevels(elem))=0.0_WP 
   DO nz=1,nlevels(elem)-1
      UV_rhsAB(1,nz,elem)=UV_rhsAB(1,nz,elem)+(uvertAB(1,nz)-uvertAB(1,nz+1))/(zbar(nz)-zbar(nz+1))
      UV_rhsAB(2,nz,elem)=UV_rhsAB(2,nz,elem)+(uvertAB(2,nz)-uvertAB(2,nz+1))/(zbar(nz)-zbar(nz+1))
   END DO
END DO
END SUBROUTINE momentum_adv_vert
! ===================================================================================
SUBROUTINE momentum_adv_scalar
USE o_MESH
USE o_ARRAYS
USE o_PARAM
USE g_PARSUP
use g_comm_auto
IMPLICIT NONE

real(kind=WP)               :: tvol, tx, ty, uu, vv
integer                     :: n, elem, nz, k, elnodes(3), m
integer                     :: nl1, nl2, nodes(2), el(2), ed
real(kind=WP)               :: un, x1, y1, x2, y2, wu(1:nl),wv(1:nl),wuw(1:nl),wvw(1:nl) 
real(kind=WP), allocatable  ::  Unode_rhs(:,:,:)
allocate(Unode_rhs(2,nl-1,myDim_nod2d+eDim_nod2D))
Unode_rhs=0.0_WP
! ======
! Momentum advection on scalar control volumes
! ======
! ============
! The horizontal part
! ============

DO ed=1, myDim_edge2D
   nodes=edges(:,ed)   
   el=edge_tri(:,ed)
   nl1=nlevels(el(1))-1
   x1=edge_cross_dxdy(1,ed)
   y1=edge_cross_dxdy(2,ed)
   if (el(2)>0) then
      nl2=nlevels(el(2))-1
      x2=edge_cross_dxdy(3,ed)
      y2=edge_cross_dxdy(4,ed)
   else
      nl2=0
      x2=0.0_WP
      y2=0.0_WP
   end if
   DO nz=1, nl1
      un=(UV(2,nz,el(1))*x1- UV(1,nz,el(1))*y1)
      Unode_rhs(1,nz,nodes(1))=Unode_rhs(1,nz,nodes(1))+un*UV(1,nz,el(1))/area(nz,nodes(1))
      Unode_rhs(1,nz,nodes(2))=Unode_rhs(1,nz,nodes(2))-un*UV(1,nz,el(1))/area(nz,nodes(2))
      Unode_rhs(2,nz,nodes(1))=Unode_rhs(2,nz,nodes(1))+un*UV(2,nz,el(1))/area(nz,nodes(1))
      Unode_rhs(2,nz,nodes(2))=Unode_rhs(2,nz,nodes(2))-un*UV(2,nz,el(1))/area(nz,nodes(2))
   END DO
   if (el(2)>0)then
      DO nz=1, nl2
         un=-(UV(2,nz,el(2))*x2- UV(1,nz,el(2))*y2)
         Unode_rhs(1,nz,nodes(1))=Unode_rhs(1,nz,nodes(1))+un*UV(1,nz,el(2))/area(nz,nodes(1))
         Unode_rhs(1,nz,nodes(2))=Unode_rhs(1,nz,nodes(2))-un*UV(1,nz,el(2))/area(nz,nodes(2))
         Unode_rhs(2,nz,nodes(1))=Unode_rhs(2,nz,nodes(1))+un*UV(2,nz,el(2))/area(nz,nodes(1))
         Unode_rhs(2,nz,nodes(2))=Unode_rhs(2,nz,nodes(2))-un*UV(2,nz,el(2))/area(nz,nodes(2))
      END DO
   END IF
END DO
! =============
! The vertical part
! =============
DO ed=1, myDim_edge2D
   nodes=edges(:,ed)   
   el=edge_tri(:,ed)
   nl1=nlevels(el(1))-1
   if (el(2)>0) then
      nl2=nlevels(el(2))-1
   end if
   ! Here 1/6 because the 1/6 of area is related to the edge
   wu(1)=UV(1,1,el(1))*elem_area(el(1))/6.0_WP
   wu(2:nl1)=(UV(1,2:nl1,el(1))+UV(1,1:nl1-1,el(1)))*0.5_WP*elem_area(el(1))/6.0_WP
   wu(nl1+1)=0.0_WP
   
   wv(1)=UV(2,1,el(1))*elem_area(el(1))/6.0_WP
   wv(2:nl1)=(UV(2,2:nl1,el(1))+UV(2,1:nl1-1,el(1)))*0.5_WP*elem_area(el(1))/6.0_WP
   wv(nl1+1)=0.0_WP
   
   wuw(1:nl1+1)=wu*Wvel_e(1:nl1+1,nodes(1))
   wvw(1:nl1+1)=wv*Wvel_e(1:nl1+1,nodes(1))
   DO nz=1,nl1 
      Unode_rhs(1,nz,nodes(1))=Unode_rhs(1,nz,nodes(1))-(wuw(nz)-wuw(nz+1))/ &
           (zbar(nz)-zbar(nz+1))/area(nz,nodes(1))
      Unode_rhs(2,nz,nodes(1))=Unode_rhs(2,nz,nodes(1))-(wvw(nz)-wvw(nz+1))/ &
           (zbar(nz)-zbar(nz+1))/area(nz,nodes(1))
   END DO 
   
   wuw(1:nl1+1)=wu*Wvel_e(1:nl1+1,nodes(2))
   wvw(1:nl1+1)=wv*Wvel_e(1:nl1+1,nodes(2))
   DO nz=1,nl1 
      Unode_rhs(1,nz,nodes(2))=Unode_rhs(1,nz,nodes(2))-(wuw(nz)-wuw(nz+1))/ &
           (zbar(nz)-zbar(nz+1))/area(nz,nodes(2))
      Unode_rhs(2,nz,nodes(2))=Unode_rhs(2,nz,nodes(2))-(wvw(nz)-wvw(nz+1))/ &
           (zbar(nz)-zbar(nz+1))/area(nz,nodes(2))
   END DO 
   
   if (el(2)>0) then
      wu(1)=UV(1,1,el(2))*elem_area(el(2))/6.0_WP
      wu(2:nl2)=(UV(1,2:nl2,el(2))+UV(1,1:nl2-1,el(2)))*0.5_WP*elem_area(el(2))/6.0_WP
      wu(nl2+1)=0.0_WP
   
      wv(1)=UV(2,1,el(2))*elem_area(el(2))/6.0_WP
      wv(2:nl2)=(UV(2,2:nl2,el(2))+UV(2,1:nl2-1,el(2)))*0.5_WP*elem_area(el(2))/6.0_WP
      wv(nl2+1)=0.0_WP
   
      wuw(1:nl2+1)=wu*Wvel_e(1:nl2+1,nodes(1))
      wvw(1:nl2+1)=wv*Wvel_e(1:nl2+1,nodes(1))
      DO nz=1,nl2 
         Unode_rhs(1,nz,nodes(1))=Unode_rhs(1,nz,nodes(1))-(wuw(nz)-wuw(nz+1))/ &
               (zbar(nz)-zbar(nz+1))/area(nz,nodes(1))
         Unode_rhs(2,nz,nodes(1))=Unode_rhs(2,nz,nodes(1))-(wvw(nz)-wvw(nz+1))/ &
               (zbar(nz)-zbar(nz+1))/area(nz,nodes(1))
      END DO 
   
      wuw(1:nl2+1)=wu*Wvel_e(1:nl2+1,nodes(2))
      wvw(1:nl2+1)=wv*Wvel_e(1:nl2+1,nodes(2))
      DO nz=1,nl2 
         Unode_rhs(1,nz,nodes(2))=Unode_rhs(1,nz,nodes(2))-(wuw(nz)-wuw(nz+1))/ &
                (zbar(nz)-zbar(nz+1))/area(nz,nodes(2))
         Unode_rhs(2,nz,nodes(2))=Unode_rhs(2,nz,nodes(2))-(wvw(nz)-wvw(nz+1))/ &
                (zbar(nz)-zbar(nz+1))/area(nz,nodes(2))
      END DO 
   end if
END DO

call exchange_nod(Unode_rhs)
DO elem=1, myDim_elem2D
   elnodes=elem2D_nodes(:,elem)
   DO nz=1,nlevels(elem)-1
      UV_rhsAB(1,nz,elem)=UV_rhsAB(1,nz,elem)+elem_area(elem)*sum(Unode_rhs(1,nz,elnodes))/3.0_WP
      UV_rhsAB(2,nz,elem)=UV_rhsAB(2,nz,elem)+elem_area(elem)*sum(Unode_rhs(2,nz,elnodes))/3.0_WP
   END DO
END DO
deallocate(Unode_rhs)
END subroutine momentum_adv_scalar
! ===================================================================
