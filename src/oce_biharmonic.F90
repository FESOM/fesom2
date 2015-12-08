! Biharmonic viscosity routines.
! Contains: biharmonic_viscosity (only biharmonic operator)
!           biPlusharmonic_viscosity (harmonic+biharmonic)     
!           vel_lapl_gradiens ==vel_gradients, but simplified (no account for free slip). 
!       
! True viscosity operator is implemented in laplacian_viscosity (div sigma_ij, 
! with sigma_ij=2*A_hor(epsilon_ij-delta_ij*div u/2). If we apply it twice, the resultant
! biharmonic operator does not kill noise if eddies are strong. We use operator
! Lu=nabla nabla u, corrected as suggested by Blazek, 2001, (the correction concerns the
! way nabla u is computed on cell-centered discretization) and apply it twice:
! bh_visc=-L Abh Lu.
!
! If Laplacian and biharmonic are used simultaneously, we assemble 
! lapl_visc=nabla A_hor nabla u and bh_visc in the same subroutine to 
! increase numerical efficiency. Note, however, that operators introduced in that way are
! not zero under solid body rotations, while the operator in laplacian_viscosity is.
! 
! Note also that we ignore metric differentiation terms because our viscosity
! is just a parameterization. There is no problem with including them.
!
! sergey.danilov@awi.de  2012
! ==============================================================
SUBROUTINE vel_lapl_gradients(UV_c)
  ! Similar to vel_gradients, the difference is  
  ! the argument and the implementation of free slip option
  USE o_MESH
  USE o_ARRAYS
  USE o_PARAM
  USE g_PARSUP
  use g_comm_auto
  IMPLICIT NONE
  real(kind=WP)    :: x(3), y(3), u, v, r1, r2
  integer         :: elem, el, j, nz
  real(kind=WP)    :: UV_c(nl-1,2,myDim_elem2D+eDim_elem2D)

  DO elem=1,myDim_elem2D  !! m=1,myDim_elem2D

     vel_grad(1:nlevels(elem)-1,:,elem) = 0._WP
     !vel_grad(nz,:,n)=(ux,uy,vx,vy)
     DO j=1,3
        el=elem_neighbors(j,elem)

        if (el>0) then
           ! ======================
           ! fill in virtual values u, v
           ! ======================
           ! ======================
           ! Filling in velocity gradients:
           ! ======================
           if(nlevels(el) < nlevels(elem)) then
              DO nz=1, nlevels(el)-1
                 u = UV_c(nz,1,el)-UV_c(nz,1,elem)
                 v = UV_c(nz,2,el)-UV_c(nz,2,elem)
                vel_grad(nz,1,elem) = vel_grad(nz,1,elem)+gradient_vec(j,elem)*u
                vel_grad(nz,2,elem) = vel_grad(nz,2,elem)+gradient_vec(j+3,elem)*u
                vel_grad(nz,3,elem) = vel_grad(nz,3,elem)+gradient_vec(j,elem)*v
                vel_grad(nz,4,elem) = vel_grad(nz,4,elem)+gradient_vec(j+3,elem)*v
              END DO

              DO nz=nlevels(el),nlevels(elem)-1
                 u = -2.*UV_c(nz,1,elem)
                 v = -2.*UV_c(nz,2,elem)
                vel_grad(nz,1,elem) = vel_grad(nz,1,elem)+gradient_vec(j,elem)*u
                vel_grad(nz,2,elem) = vel_grad(nz,2,elem)+gradient_vec(j+3,elem)*u
                vel_grad(nz,3,elem) = vel_grad(nz,3,elem)+gradient_vec(j,elem)*v
                vel_grad(nz,4,elem) = vel_grad(nz,4,elem)+gradient_vec(j+3,elem)*v  
              END DO

           else
              DO nz=1, nlevels(elem)-1
                 u = UV_c(nz,1,el)-UV_c(nz,1,elem)
                 v = UV_c(nz,2,el)-UV_c(nz,2,elem)
                vel_grad(nz,1,elem) = vel_grad(nz,1,elem)+gradient_vec(j,elem)*u
                vel_grad(nz,2,elem) = vel_grad(nz,2,elem)+gradient_vec(j+3,elem)*u
                vel_grad(nz,3,elem) = vel_grad(nz,3,elem)+gradient_vec(j,elem)*v
                vel_grad(nz,4,elem) = vel_grad(nz,4,elem)+gradient_vec(j+3,elem)*v
              END DO
           end if

        else

           ! ===============
           ! Boundary element
           ! ===============
           ! ======================
           ! Filling in velocity gradients:
           ! ======================
           DO nz=1, nlevels(elem)-1
                 u = -2.*UV_c(nz,1,elem)
                 v = -2.*UV_c(nz,2,elem)
             vel_grad(nz,1,elem) = vel_grad(nz,1,elem) + gradient_vec(j,elem)*u
             vel_grad(nz,2,elem) = vel_grad(nz,2,elem) + gradient_vec(j+3,elem)*u
             vel_grad(nz,3,elem) = vel_grad(nz,3,elem) + gradient_vec(j,elem)*v
             vel_grad(nz,4,elem) = vel_grad(nz,4,elem) + gradient_vec(j+3,elem)*v
           END DO
        end if
     end do   ! cycle over neighbor elements 	 
  END DO


  call exchange_elem(vel_grad)

END SUBROUTINE vel_lapl_gradients
!
!===========================================================================
! With Laplacian regularization --- does not help in reality.
SUBROUTINE biharmonic_viscosity
  USE o_MESH
  USE o_ARRAYS
  USE o_PARAM
  USE g_PARSUP
  use g_comm_auto
  IMPLICIT NONE
  real(kind=WP)          :: xe, ye, u1, v1, Ah  
  real(kind=WP)          :: x1, x2, y1, y2, tx, ty, tl, tt, g1, g2
  integer               :: ed, nodes(2), el(2), nl1, nl2, nz, elem, m
  real(kind=WP), allocatable  :: UV_c(:,:,:)

  ed=myDim_elem2D+eDIm_elem2D
  allocate(UV_c(nl-1,2,ed))

!  DO elem=1, myDim_elem2D  
!     DO nz=1, nlevels(elem)-1
!        UV_c(nz,:,elem)=0.0_WP
!     END DO
!  END DO

  UV_c=0.0_WP

  DO ed=1, myDim_edge2D+eDim_edge2D 
     if(myList_edge2D(ed)>edge2D_in) cycle
     nodes=edges(:,ed)   
     el=edge_tri(:,ed)
     nl1=nlevels(el(1))-1
     nl2=nlevels(el(2))-1
     xe=edge_dxdy(1,ed)*r_earth*0.5_WP*(elem_cos(el(1))+elem_cos(el(2)))
     ye=edge_dxdy(2,ed)*r_earth

     x1=-edge_cross_dxdy(1,ed)
     y1=-edge_cross_dxdy(2,ed)
     x2=-edge_cross_dxdy(3,ed)
     y2=-edge_cross_dxdy(4,ed)
     tx=x1-x2
     ty=y1-y2
     tl=1./sqrt(tx**2+ty**2)
     tx=tx*tl
     ty=ty*tl

     DO nz=1, min(nl1,nl2)   
        g1=0.5_WP*(vel_grad(nz,1,el(1))+vel_grad(nz,1,el(2)))
        g2=0.5_WP*(vel_grad(nz,2,el(1))+vel_grad(nz,2,el(2)))
        tt=g1*tx+g2*ty-(UV(1,nz,el(2))-UV(1,nz,el(1)))*tl
        g1=g1-tx*tt
        g2=g2-ty*tt
        u1=(g1*ye-g2*xe)
        g1=0.5_WP*(vel_grad(nz,3,el(1))+vel_grad(nz,3,el(2)))
        g2=0.5_WP*(vel_grad(nz,4,el(1))+vel_grad(nz,4,el(2)))
        tt=g1*tx+g2*ty-(UV(2,nz,el(2))-UV(2,nz,el(1)))*tl
        g1=g1-tx*tt
        g2=g2-ty*tt
        v1=(g1*ye-g2*xe)

        UV_c(nz,1,el(1))=UV_c(nz,1,el(1))+u1
        UV_c(nz,2,el(1))=UV_c(nz,2,el(1))+v1
        UV_c(nz,1,el(2))=UV_c(nz,1,el(2))-u1
        UV_c(nz,2,el(2))=UV_c(nz,2,el(2))-v1

     END DO

     if(nl1>nl2) then 
        DO nz=nl2+1,nl1
           u1=vel_grad(nz,1,el(1))*ye-vel_grad(nz,2,el(1))*xe
           v1=vel_grad(nz,3,el(1))*ye-vel_grad(nz,4,el(1))*xe
           UV_c(nz,1,el(1))=UV_c(nz,1,el(1))+u1
           UV_c(nz,2,el(1))=UV_c(nz,2,el(1))+v1
        END DO
     end if
     if(nl1<nl2) then 
        DO nz=nl1+1,nl2
           u1=vel_grad(nz,1,el(2))*ye-vel_grad(nz,2,el(2))*xe
           v1=vel_grad(nz,3,el(2))*ye-vel_grad(nz,4,el(2))*xe
           UV_c(nz,1,el(2))=UV_c(nz,1,el(2))-u1
           UV_c(nz,2,el(2))=UV_c(nz,2,el(2))-v1
        END DO
     END IF

  END DO

  ! ============ 
  ! Contribution from boundary edges
  ! ============
  DO ed=1, myDim_edge2D+eDim_edge2D  !! m=1, myDim_edge2D+eDim_edge2D  
     !! ed=myList_edge2D(m)
     if(myList_edge2D(ed)<=edge2D_in) cycle
     el=edge_tri(:,ed)
     xe=edge_dxdy(1,ed)*r_earth*elem_cos(el(1)) 
     ye=edge_dxdy(2,ed)*r_earth

     DO nz=1, nlevels(el(1))-1
        u1=(vel_grad(nz,1,el(1)))*ye-(vel_grad(nz,2,el(1)))*xe
        v1=(vel_grad(nz,3,el(1)))*ye-(vel_grad(nz,4,el(1)))*xe
        UV_c(nz,1,el(1))=UV_c(nz,1,el(1))+u1
        UV_c(nz,2,el(1))=UV_c(nz,2,el(1))+v1
     END DO
  END DO

  ! =============  
  ! Multiply with biharmonic viscosity
  ! =============
  DO elem=1, myDim_elem2D     !! m=1, myDim_elem2D       
     !! elem=myList_elem2D(m)
     !Ah=-Abh0*elem_area(elem)/scale_area**2   ! we have to divide on elem_area
     ! - takes into account the fact that 
     ! biharmonic goes with different sign.   
     Ah=-Abh0*sqrt(elem_area(elem)/scale_area)/scale_area
     DO nz=1, nlevels(elem)-1
        UV_c(nz,1,elem)=UV_c(nz,1,elem)*(Ah-Visc(nz,elem))
        UV_c(nz,2,elem)=UV_c(nz,2,elem)*(Ah-Visc(nz,elem))
     END DO
  END DO

 ! Visc should be multiplied with elem_area to become biharmonic viscosity.
 ! However, UV_c contains contributions over the area, so in the result 
 ! no multiplication is required. For the same reason Ah does not contain
 ! this multiplication too.  


  call exchange_elem(UV_c)
!  call exchange_elem(UV_c(:,1,:))
!  call exchange_elem(UV_c(:,2,:))
  call vel_lapl_gradients(UV_c)

  ! =============
  ! Apply Laplace operator once more
  ! =============
  DO ed=1, myDim_edge2D+eDim_edge2D !! m=1, myDim_edge2D+eDim_edge2D
     !! ed=myList_edge2D(m)
     if(myList_edge2D(ed)>edge2D_in) cycle
     nodes=edges(:,ed)   
     el=edge_tri(:,ed)
     nl1=nlevels(el(1))-1
     nl2=nlevels(el(2))-1
     xe=edge_dxdy(1,ed)*r_earth*0.5_WP*(elem_cos(el(1))+elem_cos(el(2)))
     ye=edge_dxdy(2,ed)*r_earth
     x1=-edge_cross_dxdy(1,ed)
     y1=-edge_cross_dxdy(2,ed)
     x2=-edge_cross_dxdy(3,ed)
     y2=-edge_cross_dxdy(4,ed)

     tx=x1-x2
     ty=y1-y2
     tl=1./sqrt(tx**2+ty**2)
     tx=tx*tl
     ty=ty*tl

     DO nz=1, min(nl1,nl2)   
        g1=0.5_WP*(vel_grad(nz,1,el(1))+vel_grad(nz,1,el(2)))
        g2=0.5_WP*(vel_grad(nz,2,el(1))+vel_grad(nz,2,el(2)))
        tt=g1*tx+g2*ty-(UV_c(nz,1,el(2))-UV_c(nz,1,el(1)))*tl
        g1=g1-tx*tt
        g2=g2-ty*tt
        u1=(g1*ye-g2*xe)
        g1=0.5_WP*(vel_grad(nz,3,el(1))+vel_grad(nz,3,el(2)))
        g2=0.5_WP*(vel_grad(nz,4,el(1))+vel_grad(nz,4,el(2)))
        tt=g1*tx+g2*ty-(UV_c(nz,2,el(2))-UV_c(nz,2,el(1)))*tl
        g1=g1-tx*tt
        g2=g2-ty*tt
        v1=(g1*ye-g2*xe)

        UV_rhs(1,nz,el(1))=UV_rhs(1,nz,el(1))+u1
        UV_rhs(2,nz,el(1))=UV_rhs(2,nz,el(1))+v1
        UV_rhs(1,nz,el(2))=UV_rhs(1,nz,el(2))-u1
        UV_rhs(2,nz,el(2))=UV_rhs(2,nz,el(2))-v1

     END DO
     if(nl1>nl2) then 
        DO nz=nl2+1,nl1
           u1=vel_grad(nz,1,el(1))*ye-vel_grad(nz,2,el(1))*xe
           v1=vel_grad(nz,3,el(1))*ye-vel_grad(nz,4,el(1))*xe
           UV_rhs(1,nz,el(1))=UV_rhs(1,nz,el(1))+u1
           UV_rhs(2,nz,el(1))=UV_rhs(2,nz,el(1))+v1
        END DO
     end if
     if(nl1<nl2) then 
        DO nz=nl1+1,nl2
           u1=vel_grad(nz,1,el(2))*ye-vel_grad(nz,2,el(2))*xe
           v1=vel_grad(nz,3,el(2))*ye-vel_grad(nz,4,el(2))*xe
           UV_rhs(1,nz,el(2))=UV_rhs(1,nz,el(2))-u1
           UV_rhs(2,nz,el(2))=UV_rhs(2,nz,el(2))-v1
        END DO
     END IF

  END DO
  ! ============ 
  ! Contribution from boundary edges
  ! ============
  DO ed=1, myDim_edge2D+eDim_edge2D !! m=1, myDim_edge2D+eDim_edge2D 
     !! ed=myList_edge2D(m)
     if(myList_edge2D(ed)<=edge2D_in) cycle
     nodes=edges(:,ed)   
     el=edge_tri(:,ed)
     xe=edge_dxdy(1,ed)*r_earth*elem_cos(el(1)) 
     ye=edge_dxdy(2,ed)*r_earth

     DO nz=1, nlevels(el(1))-1
        u1=(vel_grad(nz,1,el(1)))*ye-(vel_grad(nz,2,el(1)))*xe
        v1=(vel_grad(nz,3,el(1)))*ye-(vel_grad(nz,4,el(1)))*xe
        UV_rhs(1,nz,el(1))=UV_rhs(1,nz,el(1))+u1
        UV_rhs(2,nz,el(1))=UV_rhs(2,nz,el(1))+v1
     END DO
  END DO

  deallocate(UV_c)
END subroutine biharmonic_viscosity
!
! =========================================================================
!
SUBROUTINE biPlusharmonic_viscosity
  !
  ! Compute harmonic and biharmonic in a single sweep. See comments in
  ! biharmonic_viscosity
  !
  USE o_MESH
  USE o_ARRAYS
  USE o_PARAM
  USE g_PARSUP
  use g_comm_auto
  IMPLICIT NONE
  real(kind=WP)    :: xe, ye, u1, v1, Ah, uc1, vc1, xx  
  real(kind=WP)    :: x1, x2, y1, y2,tl, tt, tx, ty, g1, g2
  integer          :: ed, nl1, nl2, nz, elem, m, el1, el2
  real(kind=WP)    :: UV_c(nl-1, 2, myDim_elem2D+eDIm_elem2D)

  UV_c=0.0_WP
  

  call h_viscosity

  DO ed=1, myDim_edge2D+eDim_edge2D !! m=1, myDim_edge2D+eDim_edge2D
     !! ed=myList_edge2D(m)

     if(myList_edge2D(ed) <= edge2D_in) then

     el1 = edge_tri(1,ed)
     el2 = edge_tri(2,ed)
     nl1= nlevels(el1)-1
     nl2= nlevels(el2)-1
     xe = edge_dxdy(1,ed)*r_earth*0.5_WP*(elem_cos(el1)+elem_cos(el2))
     ye = edge_dxdy(2,ed)*r_earth

     tx = edge_cross_dxdy(3,ed) - edge_cross_dxdy(1,ed)
     ty = edge_cross_dxdy(4,ed) - edge_cross_dxdy(2,ed)
     tl = 1./sqrt(tx**2+ty**2)
     tx = tx*tl
     ty = ty*tl

     DO nz=1, min(nl1,nl2)   
        g1 = 0.5_WP*(vel_grad(nz,1,el1)+vel_grad(nz,1,el2))
        g2 = 0.5_WP*(vel_grad(nz,2,el1)+vel_grad(nz,2,el2))
        tt = g1*tx+g2*ty-(UV(1,nz,el2)-UV(1,nz,el1))*tl
        u1 = (g1-tx*tt)*ye - (g2-ty*tt)*xe

        g1 = 0.5_WP*(vel_grad(nz,3,el1)+vel_grad(nz,3,el2))
        g2 = 0.5_WP*(vel_grad(nz,4,el1)+vel_grad(nz,4,el2))
        tt = g1*tx + g2*ty - (UV(2,nz,el2)-UV(2,nz,el1))*tl
        v1 = (g1-tx*tt)*ye - (g2-ty*tt)*xe

        Ah = (Visc(nz,el1)+Visc(nz,el2))*0.5

        UV_c(nz,1,el1)=UV_c(nz,1,el1)+u1
        UV_c(nz,2,el1)=UV_c(nz,2,el1)+v1
        UV_c(nz,1,el2)=UV_c(nz,1,el2)-u1
        UV_c(nz,2,el2)=UV_c(nz,2,el2)-v1
        UV_rhs(1,nz,el1)=UV_rhs(1,nz,el1)+u1*Ah
        UV_rhs(2,nz,el1)=UV_rhs(2,nz,el1)+v1*Ah
        UV_rhs(1,nz,el2)=UV_rhs(1,nz,el2)-u1*Ah
        UV_rhs(2,nz,el2)=UV_rhs(2,nz,el2)-v1*Ah

     END DO

     if(nl1>nl2) then 
        DO nz=nl2+1,nl1
           u1=vel_grad(nz,1,el1)*ye-vel_grad(nz,2,el1)*xe
           v1=vel_grad(nz,3,el1)*ye-vel_grad(nz,4,el1)*xe
           UV_c(nz,1,el1)=UV_c(nz,1,el1)+u1
           UV_c(nz,2,el1)=UV_c(nz,2,el1)+v1
           UV_rhs(1,nz,el1)=UV_rhs(1,nz,el1)+u1*Visc(nz, el1)
           UV_rhs(2,nz,el1)=UV_rhs(2,nz,el1)+v1*Visc(nz, el1)
        END DO
     elseif(nl1<nl2) then 
        DO nz=nl1+1,nl2
           u1=vel_grad(nz,1,el2)*ye-vel_grad(nz,2,el2)*xe
           v1=vel_grad(nz,3,el2)*ye-vel_grad(nz,4,el2)*xe
           UV_c(nz,1,el2)=UV_c(nz,1,el2)-u1
           UV_c(nz,2,el2)=UV_c(nz,2,el2)-v1
           UV_rhs(1,nz,el2)=UV_rhs(1,nz,el2)-u1*Visc(nz, el2)
           UV_rhs(2,nz,el2)=UV_rhs(2,nz,el2)-v1*Visc(nz, el2)

        END DO
     endif
     else
  ! ============ 
  ! Contribution from boundary edges
  ! ============
  
     el1=edge_tri(1,ed)
     xe=edge_dxdy(1,ed)*r_earth*elem_cos(el1) 
     ye=edge_dxdy(2,ed)*r_earth

     DO nz=1, nlevels(el1)-1
        u1=(vel_grad(nz,1,el1))*ye-(vel_grad(nz,2,el1))*xe
        v1=(vel_grad(nz,3,el1))*ye-(vel_grad(nz,4,el1))*xe
        UV_c(nz,1,el1)=UV_c(nz,1,el1)+u1
        UV_c(nz,2,el1)=UV_c(nz,2,el1)+v1
        UV_rhs(1,nz,el1)=UV_rhs(1,nz,el1)+u1*Visc(nz, el1)
        UV_rhs(2,nz,el1)=UV_rhs(2,nz,el1)+v1*Visc(nz, el1)
     END DO
  endif
  END DO

  ! =============  
  ! Multiply with biharmonic viscosity
  ! =============
  DO elem=1, myDim_elem2D                    ! U_c, V_c are area-weighted 
     ! we have to divide on elem_area
     ! - takes into account the fact that 
     ! biharmonic goes with different sign.   
     Ah=-Abh0*sqrt(elem_area(elem)/scale_area)/scale_area
     DO nz=1, nlevels(elem)-1
        UV_c(nz,1,elem) = UV_c(nz,1,elem)*(Ah-Visc(nz,elem))
        UV_c(nz,2,elem) = UV_c(nz,2,elem)*(Ah-Visc(nz,elem))
     END DO
  END DO

 ! Visc should be multiplied with elem_area to become biharmonic viscosity.
 ! However, UV_c contains contributions over the area, so in the result 
 ! no multiplication is required. For the same reason Ah does not contain
 ! this multiplication too.  


  call exchange_elem(UV_c(:,1,:))
  call exchange_elem(UV_c(:,2,:))
  call vel_lapl_gradients(UV_c)

  ! =============
  ! Apply Laplace operator once again
  ! =============
  DO ed=1, myDim_edge2D+eDim_edge2D

     if(myList_edge2D(ed) <= edge2D_in) then

     el1=edge_tri(1,ed)
     el2=edge_tri(2,ed)
     nl1=nlevels(el1)-1
     nl2=nlevels(el2)-1
     xe=edge_dxdy(1,ed)*r_earth*0.5_WP*(elem_cos(el1)+elem_cos(el2))
     ye=edge_dxdy(2,ed)*r_earth

     tx = edge_cross_dxdy(3,ed) - edge_cross_dxdy(1,ed)
     ty = edge_cross_dxdy(4,ed) - edge_cross_dxdy(2,ed)
     tl = 1./sqrt(tx**2+ty**2)
     tx = tx*tl
     ty = ty*tl

     DO nz=1, min(nl1,nl2)   
        g1 = 0.5_WP*(vel_grad(nz,1,el1)+vel_grad(nz,1,el2))
        g2 = 0.5_WP*(vel_grad(nz,2,el1)+vel_grad(nz,2,el2))
        tt = g1*tx+g2*ty-(UV_c(nz,1,el2)-UV_c(nz,1,el1))*tl
        u1 = (g1-tx*tt)*ye - (g2-ty*tt)*xe

        g1 = 0.5_WP*(vel_grad(nz,3,el1)+vel_grad(nz,3,el2))
        g2 = 0.5_WP*(vel_grad(nz,4,el1)+vel_grad(nz,4,el2))
        tt = g1*tx+g2*ty-(UV_c(nz,2,el2)-UV_c(nz,2,el1))*tl
        v1 = (g1-tx*tt)*ye - (g2-ty*tt)*xe

        UV_rhs(1,nz,el1) = UV_rhs(1,nz,el1)+u1
        UV_rhs(2,nz,el1) = UV_rhs(2,nz,el1)+v1
        UV_rhs(1,nz,el2) = UV_rhs(1,nz,el2)-u1
        UV_rhs(2,nz,el2) = UV_rhs(2,nz,el2)-v1

     END DO
     if(nl1>nl2) then 
        DO nz=nl2+1,nl1
           u1 = vel_grad(nz,1,el1)*ye-vel_grad(nz,2,el1)*xe
           v1 = vel_grad(nz,3,el1)*ye-vel_grad(nz,4,el1)*xe
           UV_rhs(1,nz,el1) = UV_rhs(1,nz,el1)+u1
           UV_rhs(2,nz,el1) = UV_rhs(2,nz,el1)+v1
        END DO
     elseif (nl1<nl2) then 
        DO nz=nl1+1,nl2
           u1 = vel_grad(nz,1,el2)*ye-vel_grad(nz,2,el2)*xe
           v1 = vel_grad(nz,3,el2)*ye-vel_grad(nz,4,el2)*xe
           UV_rhs(1,nz,el2) = UV_rhs(1,nz,el2)-u1
           UV_rhs(2,nz,el2) = UV_rhs(2,nz,el2)-v1
        END DO
     endif

  else
  ! ============ 
  ! Contribution from boundary edges
  ! ============

     el1=edge_tri(1,ed)
     xe=edge_dxdy(1,ed)*r_earth*elem_cos(el1) 
     ye=edge_dxdy(2,ed)*r_earth

     DO nz=1, nlevels(el1)-1
        u1=(vel_grad(nz,1,el1))*ye-(vel_grad(nz,2,el1))*xe
        v1=(vel_grad(nz,3,el1))*ye-(vel_grad(nz,4,el1))*xe
        UV_rhs(1,nz,el1)=UV_rhs(1,nz,el1)+u1
        UV_rhs(2,nz,el1)=UV_rhs(2,nz,el1)+v1
     END DO
  end if
  END DO

END subroutine biPlusharmonic_viscosity
! ==============================================================
