SUBROUTINE h_viscosity
!
! Coefficient of horizontal viscosity is a combination of 
! Smagorinsky contribution (with Smag2)
! Modified Leith (with Div_c)
! and background (with A_hor)
! Background is scaled as sqrt(A/A_0), others are scaled
! in natural way by construction
!
USE o_MESH
USE o_ARRAYS
USE o_PARAM
USE g_PARSUP
use g_comm_auto
IMPLICIT NONE
real(kind=WP)  :: Smag2, div_elem(3), d1, d2, xe, ye, vi, vimax
integer        :: elem, nz, elnod(3)
real(kind=WP)  :: a_horizontal, dz2_inv(nl-1)

  Smag2 = Smag_c*Smag_c 
   
!NR precompute auxilliary array: depending only on nz, not on elem.
!NR Compare with original computation in comments.
  do nz=1,nl-1
     dz2_inv(nz) = Div_c /((zbar(nz)-zbar(nz+1))**2)
  enddo

  ! Fill in viscosity:
  DO  elem=1, myDim_elem2D
      vi=sqrt(elem_area(elem)/scale_area)
      vimax=vi*A_hor_max
      vi=vi*A_hor 
      elnod=elem2D_nodes(:,elem)

      do nz=1,nlevels(elem)-1

	div_elem= Wvel(nz,elnod) - Wvel(nz+1,elnod)
        xe      = sum(gradient_sca(1:3,elem)*div_elem)
	ye      = sum(gradient_sca(4:6,elem)*div_elem)

        d1      = vel_grad(nz,1,elem) - vel_grad(nz,4,elem)
        d2      = vel_grad(nz,2,elem) + vel_grad(nz,3,elem)

        Visc(nz,elem) = min(elem_area(elem)* &
	    sqrt(Smag2*(d1*d1+d2*d2) + elem_area(elem)*(xe*xe+ye*ye)*dz2_inv(nz)) &
	     + vi, vimax)
      end do

      do nz=nlevels(elem), nl-1
        Visc(nz, elem)=0.0_WP
      end do		    
  END DO 
  call exchange_elem(Visc)
END subroutine h_viscosity  
! =======================================================================
SUBROUTINE h_viscosity2
!
! Similar to h_viscosity, but with an account for the 
! Leith viscosity. It needs vorticity, which is only computed for 
! the vector invariant form of momentum advection 
!
! Coefficient of horizontal viscosity is a combination of 
! Smagorinsky contribution (with Smag2)
! Modified Leith (with Div_c)
! and background (with A_hor)
! Background is scaled as sqrt(A/A_0), others are scaled
! in natural way by construction
!
USE o_MESH
USE o_ARRAYS
USE o_PARAM
USE g_PARSUP
use g_comm_auto
IMPLICIT NONE
real(kind=WP)  :: Smag2, div_elem(3), rot_elem(3), d1, d2, xe, ye, vi, vimax
integer        :: elem, nz, elnod(3)
real(kind=WP)  :: leith_x, leith_y
real(kind=WP)  :: a_horizontal, dz2_inv(nl-1)

  if(mom_adv.ne.4) then 
  call h_viscosity
  return
  end if         ! to use the Leith viscosity we need to have vorticity array 
                 ! allocated and filled. 

  Smag2 = Smag_c*Smag_c 
   
!NR precompute auxilliary array: depending only on nz, not on elem.
!NR Compare with original computation in comments.
  do nz=1,nl-1
     dz2_inv(nz) = Div_c /((zbar(nz)-zbar(nz+1))**2)
  enddo

  ! Fill in viscosity:
  DO  elem=1, myDim_elem2D
      vi=sqrt(elem_area(elem)/scale_area)
      vimax=vi*A_hor_max
      vi=vi*A_hor 
      elnod=elem2D_nodes(:,elem)

      do nz=1,nlevels(elem)-1

	div_elem= Wvel(nz,elnod) - Wvel(nz+1,elnod)
        xe      = sum(gradient_sca(1:3,elem)*div_elem)
	ye      = sum(gradient_sca(4:6,elem)*div_elem)

        d1      = vel_grad(nz,1,elem) - vel_grad(nz,4,elem)
        d2      = vel_grad(nz,2,elem) + vel_grad(nz,3,elem)

        rot_elem=vorticity(nz,elnod)
        leith_x=sum(gradient_sca(1:3,elem)*rot_elem)
	leith_y=sum(gradient_sca(4:6,elem)*rot_elem)

        Visc(nz,elem) = min(elem_area(elem)* &
	    sqrt(Smag2*(d1*d1+d2*d2) + elem_area(elem)*(xe*xe+ye*ye)*dz2_inv(nz)+&
		Leith_c*(leith_x*leith_x+leith_y*leith_y))+&
		+ vi, vimax)
      end do

      do nz=nlevels(elem), nl-1
        Visc(nz, elem)=0.0_WP
      end do		    
  END DO 
  call exchange_elem(Visc)
END subroutine h_viscosity2
! =======================================================================
SUBROUTINE viscosity_filtx
USE o_MESH
USE o_ARRAYS
USE o_PARAM
USE g_PARSUP
USE g_config
IMPLICIT NONE

real(kind=8)  :: u1, v1, tau_inv
integer       :: nz, ed, el(2)

   
 ! subroutine applies filter to velocity, which on uniform 
 ! equilateral mesh will be equivalent to the Laplacian operator
 ! It adds to the rhs(0) tau_filt*(u1+u2+u3-3*u0)
 ! on triangles and tau_filt*(u1+u2+u3+u4-4*u0) on quads. 
 ! The contribution from boundary edges is neglected, it can 
 ! be done in function of boundary conditions. 
 ! The result should be equivalent to appr. L/3 on triangles and 
 ! to L on quads. 
 
 
 tau_inv=dt*tau_c/3600.0/24.0     ! SET IT experimentally 
 
 DO ed=1, myDim_edge2D+eDim_edge2D
    if(myList_edge2D(ed)>edge2D_in) cycle
    el=edge_tri(:,ed)
    DO  nz=1,minval(nlevels(el))-1 
     u1=tau_inv*(UV(1,nz,el(1))-UV(1,nz,el(2)))
     v1=tau_inv*(UV(2,nz,el(1))-UV(2,nz,el(2)))
     UV_rhs(1,nz,el(1))=UV_rhs(1,nz,el(1))-u1
     UV_rhs(1,nz,el(2))=UV_rhs(1,nz,el(2))+u1
     UV_rhs(2,nz,el(1))=UV_rhs(2,nz,el(1))-v1
     UV_rhs(2,nz,el(2))=UV_rhs(2,nz,el(2))+v1
    END DO 
 END DO
    
end subroutine viscosity_filtx
! ===================================================================
SUBROUTINE viscosity_filt2x
USE o_MESH
USE o_ARRAYS
USE o_PARAM
USE g_PARSUP
USE g_config
USE g_comm_auto
IMPLICIT NONE

real(kind=8)  :: u1, v1, tau_inv, s
integer       :: ed, el(2), nz
real(kind=8), allocatable  :: UV_c(:,:,:), UV_f(:,:,:)
 
 ! Filter is applied twice. It should be approximately 
 ! equivalent to biharmonic operator with the coefficient
 ! (tau_c/day)a^3/9. Scaling inside is found to help 
 ! with smoothness in places of mesh transition. *(it makes a^3 from a^4) 
ed=myDim_elem2D+eDim_elem2D
allocate(UV_c(2,nl-1,ed)) ! to store the filtered velocity
allocate(UV_f(2,nl-1,ed)) ! to store the contributions into the RHS


 UV_c=0.0_8
 UV_f=0.0_8
 tau_inv=dt*tau_c/3600.0/24.0     ! SET IT experimentally 
  
 DO ed=1, myDim_edge2D+eDim_edge2D
    if(myList_edge2D(ed)>edge2D_in) cycle
    el=edge_tri(:,ed)
    DO  nz=1,minval(nlevels(el))-1
        u1=(UV(1,nz,el(1))-UV(1,nz,el(2)))
        v1=(UV(2,nz,el(1))-UV(2,nz,el(2)))

        UV_c(1,nz,el(1))=UV_c(1,nz,el(1))-u1
        UV_c(1,nz,el(2))=UV_c(1,nz,el(2))+u1
        UV_c(2,nz,el(1))=UV_c(2,nz,el(1))-v1
        UV_c(2,nz,el(2))=UV_c(2,nz,el(2))+v1
    END DO 
 END DO
 
 ! ============ 
 ! Contribution from boundary edges (Dirichlet boundary conditions)
 ! ============
! DO ed=1, myDim_edge2D+eDim_edge2D
!    if(myList_edge2D(ed)<=edge2D_in) cycle
!    el=edge_tri(:, ed)
!    DO  nz=1, nlevels(el(1))-1
!        UV_c(1,nz,el(1))=UV_c(1,nz,el(1))-2.0_WP*UV(1,nz,el(1))
!        UV_c(2,nz,el(1))=UV_c(2,nz,el(1))-2.0_WP*UV(2,nz,el(1))
!    END DO
! END DO

 Do ed=1,myDim_elem2D
    Do nz=1,nlevels(ed)-1
     UV_c(1,nz,ed)=-UV_c(1,nz,ed)*tau_inv*sqrt(scale_area/elem_area(ed))
     UV_c(2,nz,ed)=-UV_c(2,nz,ed)*tau_inv*sqrt(scale_area/elem_area(ed))
    END DO
 end do

 call exchange_elem(UV_c(1,:,:))
 call exchange_elem(UV_c(2,:,:))

 DO ed=1, myDim_edge2D+eDim_edge2D
    if(myList_edge2D(ed)>edge2D_in) cycle
    el=edge_tri(:,ed)
    DO  nz=1,minval(nlevels(el))-1 
        u1=(UV_c(1,nz,el(1))-UV_c(1,nz,el(2)))
        v1=(UV_c(2,nz,el(1))-UV_c(2,nz,el(2)))

        UV_f(1,nz,el(1))=UV_f(1,nz,el(1))-u1
        UV_f(1,nz,el(2))=UV_f(1,nz,el(2))+u1
        UV_f(2,nz,el(1))=UV_f(2,nz,el(1))-v1
        UV_f(2,nz,el(2))=UV_f(2,nz,el(2))+v1
    END DO 
 END DO
 
 DO ed=1, myDim_elem2D
    DO  nz=1, nlevels(ed)-1 
        u1=sqrt(UV_f(1,nz,ed)**2+UV_f(2,nz,ed)**2)+1.e-5
        v1=sqrt(UV(1,nz,ed)**2+UV(2,nz,ed)**2)
        ! we limit the maximum contribution from the filter such, that the update is less than the N (N=2 currently) times velocity
        ! this is done to force the CFL, which is otherwise exceeded in some points
	! some other criteria is welcome (i.e. like computing the eigenvalues from filtering)
        UV_rhs(1,nz,ed)=UV_rhs(1,nz,ed)+UV_f(1,nz,ed)*min(1.0_WP, 2.0_WP*v1/u1)
        UV_rhs(2,nz,ed)=UV_rhs(2,nz,ed)+UV_f(2,nz,ed)*min(1.0_WP, 2.0_WP*v1/u1)
    END DO 
 END DO

 deallocate(UV_f, UV_c)    
end subroutine viscosity_filt2x
! =======================================================================
SUBROUTINE laplacian_viscosity
USE o_MESH
USE o_ARRAYS
USE o_PARAM
USE g_PARSUP
IMPLICIT NONE

real(kind=WP)  :: x1, x2, y1, y2, xe, ye, u1, v1, xed
real(kind=WP)  :: ah, eT1, eS1, eT2, eS2, mcos
integer              :: nl1, nl2, nz, ed, el(2)
integer              :: elem, elnodes(3)
   
 DO  ed=1,myDim_edge2D+eDim_edge2D   !! m=1,myDim_edge2D+eDim_edge2D   
                                     !! ed=myList_edge2D(m)
   el=edge_tri(:,ed)
   xe=edge_dxdy(1,ed)*r_earth
   ye=edge_dxdy(2,ed)*r_earth
   
 if(myList_edge2D(ed)<=edge2D_in) then    ! internal edge, elem on both sides
   nl1=nlevels(el(1))-1
   nl2=nlevels(el(2))-1
   DO nz=1, min(nl1,nl2)   ! Only faces that do not belong to vert. walls
   eT1=vel_grad(nz,1, el(1))-vel_grad(nz,4, el(1))
   eT2=vel_grad(nz,1, el(2))-vel_grad(nz,4, el(2))
   eS1=vel_grad(nz,2, el(1))+vel_grad(nz,3, el(1))
   eS2=vel_grad(nz,2, el(2))+vel_grad(nz,3, el(2))
   
   u1=(eT1*ye-eS1*xe*elem_cos(el(1)))*Visc(nz,el(1))
   u1=0.5_WP*(u1+(eT2*ye-eS2*xe*elem_cos(el(2)))*Visc(nz,el(2)))
   v1=(eS1*ye+eT1*xe*elem_cos(el(1)))*Visc(nz,el(1))
   v1=0.5_WP*(v1+(eS2*ye+eT2*xe*elem_cos(el(2)))*Visc(nz,el(2)))
   
   UV_rhs(1,nz,el(1))=UV_rhs(1,nz,el(1))+u1
   UV_rhs(2,nz,el(1))=UV_rhs(2,nz,el(1))+v1
   UV_rhs(1,nz,el(2))=UV_rhs(1,nz,el(2))-u1
   UV_rhs(2,nz,el(2))=UV_rhs(2,nz,el(2))-v1
   END DO
  
   !======
   ! Horizontal viscosity, the contribution 
   ! from vertical walls
   !======
   DO nz=nl2+1,nl1
   Ah=Visc(nz,el(1))
   eT1=vel_grad(nz,1, el(1))-vel_grad(nz,4, el(1))
   eS1=vel_grad(nz,2, el(1))+vel_grad(nz,3, el(1))
   xed=xe*elem_cos(el(1))
   u1=Ah*(eT1*ye-eS1*xed)
   v1=Ah*(eS1*ye+eT1*xed)
   
   UV_rhs(1,nz,el(1))=UV_rhs(1,nz,el(1))+u1
   UV_rhs(2,nz,el(1))=UV_rhs(2,nz,el(1))+v1
   END DO
   
   
   DO nz=nl1+1,nl2
   Ah=Visc(nz,el(2))
   eT2=vel_grad(nz,1, el(2))-vel_grad(nz,4, el(2))
   eS2=vel_grad(nz,2, el(2))+vel_grad(nz,3, el(2))
   xed=xe*elem_cos(el(2))
   u1=Ah*(eT2*ye-eS2*xed)
   v1=Ah*(eS2*ye+eT2*xed)
   
   UV_rhs(1,nz,el(2))=UV_rhs(1,nz,el(2))-u1
   UV_rhs(2,nz,el(2))=UV_rhs(2,nz,el(2))-v1
   END DO
   
 else
   ! ============ 
   ! Viscosity contribution from boundary edges
   ! ============
   DO nz=1, nlevels(el(1))-1
   Ah=Visc(nz,el(1))
   eT1=vel_grad(nz,1, el(1))-vel_grad(nz,4, el(1))
   eS1=vel_grad(nz,2, el(1))+vel_grad(nz,3, el(1))
   xed=xe*elem_cos(el(1))
   u1=Ah*(eT1*ye-eS1*xed)
   v1=Ah*(eS1*ye+eT1*xed)
   if (free_slip) then
   ! remove tangent component
   mcos=(u1*ye-v1*xe)/(xe*xe+ye*ye) 
   ! Projection; division because xe, ye are not normalized
   u1=mcos*ye
   v1=-mcos*xe
   end if
   
   UV_rhs(1,nz,el(1))=UV_rhs(1,nz,el(1))+u1
   UV_rhs(2,nz,el(1))=UV_rhs(2,nz,el(1))+v1
   END DO
   end if
 END DO 
end subroutine laplacian_viscosity
! ===========================================================================
SUBROUTINE gradients_p1
USE o_MESH
USE o_ARRAYS
USE o_PARAM
USE g_PARSUP
  use g_comm_auto
IMPLICIT NONE
real(kind=WP)  :: tvol, tx, ty
integer       :: n, elem, nz, k, elnodes(3), m

! ======  
! Set boundary velocity to zero
! ======

! ======
! Compute gradients
! ======
  
 DO elem=1, myDim_elem2D   !! m=1, myDim_elem2D 
                           !! elem=myList_elem2D(m) 
   elnodes=elem2D_nodes(:,elem)
   DO nz=1, nlevels(elem)-1   
      vel_grad(nz, 1, elem)=sum(gradient_sca(1:3,elem)*Unode(1,nz,elnodes))
      vel_grad(nz, 2, elem)=sum(gradient_sca(4:6,elem)*Unode(1,nz,elnodes))
      vel_grad(nz, 3, elem)=sum(gradient_sca(1:3,elem)*Unode(2,nz,elnodes))
      vel_grad(nz, 4, elem)=sum(gradient_sca(4:6,elem)*Unode(2,nz,elnodes))
   END DO
 END DO

 do m=1,4
   call exchange_elem(vel_grad(:,m,:))   
 enddo
END subroutine gradients_p1
! ===================================================================
SUBROUTINE momentum_adv_p1
USE o_MESH
USE o_ARRAYS
USE o_PARAM
USE g_PARSUP
  use g_comm_auto
IMPLICIT NONE
real(kind=WP)  :: tvol, tx, ty, uu, vv
integer       :: n, elem, nz, k, elnodes(3), m
integer       :: nl1, nl2, nodes(2), el(2), ed
real(kind=WP)  :: un, u1, v1, xe, ye
! ======  
! Set boundary velocity to zero
! ======  
 
 DO ed=1, myDim_edge2D+eDim_edge2D   !! m=1, myDim_edge2D+eDim_edge2D 
                                     !! ed=myList_edge2D(m)
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
SUBROUTINE momentum_adv_scalar
USE o_MESH
USE o_ARRAYS
USE o_PARAM
USE g_PARSUP
  use g_comm_auto
IMPLICIT NONE
real(kind=WP)  :: tvol, tx, ty, uu, vv
integer       :: n, elem, nz, k, elnodes(3), m
integer       :: nl1, nl2, nodes(2), el(2), ed
real(kind=WP)  :: un, x1, y1, x2, y2, wu(1:nl),wv(1:nl),wuw(1:nl),wvw(1:nl) 
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
   if(el(2)>0) then
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
   if(el(2)>0)then
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
   if(el(2)>0) then
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
   
   if(el(2)>0) then
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

 call exchange_nod(Unode_rhs(:	,:,:))
! call exchange_nod(Unode_rhs(1,:,:))
! call exchange_nod(Unode_rhs(2,:,:))     
 DO elem=1,  myDim_elem2D     !! m=1,  myDim_elem2D   
                              !! elem=myList_elem2D(m)
   elnodes=elem2D_nodes(:,elem)
   DO nz=1,nlevels(elem)-1
   UV_rhsAB(1,nz,elem)=UV_rhsAB(1,nz,elem)+elem_area(elem)*sum(Unode_rhs(1,nz,elnodes))/3.0_WP
   UV_rhsAB(2,nz,elem)=UV_rhsAB(2,nz,elem)+elem_area(elem)*sum(Unode_rhs(2,nz,elnodes))/3.0_WP
   END DO
 END DO
deallocate(Unode_rhs)
END subroutine momentum_adv_scalar
 ! ===================================================================
