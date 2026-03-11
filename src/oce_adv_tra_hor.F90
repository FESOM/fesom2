!===============================================================================================================================
!**************** routines for horizontal tracer advection ***********************
module oce_adv_tra_hor_interfaces
  interface
! (low order upwind)
! returns flux given at edges which contributes with
! plus sign into 1st. node and with the minus sign into the 2nd node
! IF init_zero=.TRUE.  : flux will be set to zero before computation
! IF init_zero=.FALSE. : flux=flux-input flux
! flux is not multiplied with dt
    subroutine adv_tra_hor_upw1(vel, ttf, partit, mesh, flux, o_init_zero)
      use MOD_MESH
      use MOD_TRACER
      USE MOD_PARTIT
      USE MOD_PARSUP
      type(t_partit),intent(in), target :: partit
      type(t_mesh),  intent(in), target :: mesh
      real(kind=WP), intent(in)         :: ttf(   mesh%nl-1, partit%myDim_nod2D+partit%eDim_nod2D)
      real(kind=WP), intent(in)         :: vel(2, mesh%nl-1, partit%myDim_elem2D+partit%eDim_elem2D)
      real(kind=WP), intent(inout)      :: flux(  mesh%nl-1, partit%myDim_edge2D)
      logical, optional                 :: o_init_zero
    end subroutine adv_tra_hor_upw1
!===============================================================================
! MUSCL
! returns flux given at edges which contributes with
! plus sign into 1st. node and with the minus sign into the 2nd node
! IF init_zero=.TRUE.  : flux will be set to zero before computation
! IF init_zero=.FALSE. : flux=flux-input flux
! flux is not multiplied with dt
    subroutine adv_tra_hor_muscl(vel, ttf, partit, mesh, num_ord, flux, edge_up_dn_grad, nboundary_lay, o_init_zero)
      use MOD_MESH
      USE MOD_PARTIT
      USE MOD_PARSUP
      type(t_partit),intent(in), target :: partit
      type(t_mesh),  intent(in), target :: mesh
      real(kind=WP), intent(in)         :: num_ord    ! num_ord is the fraction of fourth-order contribution in the solution
      real(kind=WP), intent(in)         :: ttf(   mesh%nl-1, partit%myDim_nod2D+partit%eDim_nod2D)
      real(kind=WP), intent(in)         :: vel(2, mesh%nl-1, partit%myDim_elem2D+partit%eDim_elem2D)
      real(kind=WP), intent(inout)      :: flux(  mesh%nl-1, partit%myDim_edge2D)
      integer,       intent(in)         :: nboundary_lay(partit%myDim_nod2D+partit%eDim_nod2D)
      real(kind=WP), intent(in)         :: edge_up_dn_grad(4, mesh%nl-1, partit%myDim_edge2D)
      logical, optional                 :: o_init_zero
    end subroutine adv_tra_hor_muscl
! a not stable version of MUSCL (reconstruction in the vicinity of bottom topography is not upwind)
! it runs with FCT option only
    subroutine adv_tra_hor_mfct(vel, ttf, partit, mesh, num_ord, flux, edge_up_dn_grad,                 o_init_zero)
      use MOD_MESH
      USE MOD_PARTIT
      USE MOD_PARSUP
      type(t_partit),intent(inout), target :: partit
      type(t_mesh),  intent(in), target :: mesh
      real(kind=WP), intent(in)         :: num_ord    ! num_ord is the fraction of fourth-order contribution in the solution
      real(kind=WP), intent(in)         :: ttf(   mesh%nl-1, partit%myDim_nod2D+partit%eDim_nod2D)
      real(kind=WP), intent(in)         :: vel(2, mesh%nl-1, partit%myDim_elem2D+partit%eDim_elem2D)
      real(kind=WP), intent(inout)      :: flux(  mesh%nl-1, partit%myDim_edge2D)
      real(kind=WP), intent(in)         :: edge_up_dn_grad(4, mesh%nl-1, partit%myDim_edge2D)
      logical, optional                 :: o_init_zero
    end subroutine adv_tra_hor_mfct
    
    ! superbee advection num_ord=0: 2nd order in space num_ord=1: 2nd order in space
    ! and time through  Direct space-time scheme
    subroutine adv_tra_hor_spbee(vel, ttf, partit, mesh, num_ord, flux, edge_up_dn_grad, flag_posdef, o_init_zero)
      use MOD_MESH
      USE MOD_PARTIT
      USE MOD_PARSUP
      type(t_partit), intent(inout), target :: partit
      type(t_mesh)  , intent(in)   , target :: mesh
      real(kind=WP) , intent(in)            :: num_ord    ! num_ord is the fraction of fourth-order contribution in the solution
      real(kind=WP) , intent(in)            :: ttf(   mesh%nl-1, partit%myDim_nod2D+partit%eDim_nod2D)
      real(kind=WP) , intent(in)            :: vel(2, mesh%nl-1, partit%myDim_elem2D+partit%eDim_elem2D)
      real(kind=WP) , intent(inout)         :: flux(  mesh%nl-1, partit%myDim_edge2D)
      real(kind=WP) , intent(in)            :: edge_up_dn_grad(4, mesh%nl-1, partit%myDim_edge2D)
      logical       , intent(in)            :: flag_posdef
      logical       , optional              :: o_init_zero
    end subroutine adv_tra_hor_spbee
  end interface
end module oce_adv_tra_hor_interfaces
!
!
!===============================================================================
subroutine adv_tra_hor_upw1(vel, ttf, partit, mesh, flux, o_init_zero)
    use MOD_MESH
    USE MOD_PARTIT
    USE MOD_PARSUP
    use g_comm_auto
    implicit none
    type(t_partit),intent(in), target :: partit
    type(t_mesh),  intent(in), target :: mesh
    real(kind=WP), intent(in)         :: ttf(   mesh%nl-1, partit%myDim_nod2D+partit%eDim_nod2D)
    real(kind=WP), intent(in)         :: vel(2, mesh%nl-1, partit%myDim_elem2D+partit%eDim_elem2D)
    real(kind=WP), intent(inout)      :: flux(  mesh%nl-1, partit%myDim_edge2D)
    logical, optional                 :: o_init_zero
    logical                           :: l_init_zero
    real(kind=WP)                     :: deltaX1, deltaY1, deltaX2, deltaY2
    real(kind=WP)                     :: a, vflux
    integer                           :: el(2), enodes(2), nz, edge
    integer                           :: nu12, nl12, nl1, nl2, nu1, nu2

#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"

    l_init_zero=.true.
    if (present(o_init_zero)) then
       l_init_zero=o_init_zero
    end if
    if (l_init_zero) then
#ifndef ENABLE_OPENACC
!$OMP PARALLEL DO
#else
       !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) DEFAULT(PRESENT) VECTOR_LENGTH(acc_vl)
#endif
       do edge=1, myDim_edge2D
          do nz=1, mesh%nl-1
             flux(nz,edge)=0.0_WP
          end do
       end do
#ifndef ENABLE_OPENACC
!$OMP END PARALLEL DO
#else
       !$ACC END PARALLEL LOOP
#endif
    end if

    ! The result is the low-order solution horizontal fluxes
    ! They are put into flux
    !___________________________________________________________________________
#ifndef ENABLE_OPENACC
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(edge, deltaX1, deltaY1, deltaX2, deltaY2, &
!$OMP                       a, vflux, el, enodes, nz, nu12, nl12, nl1, nl2, nu1, nu2)
!$OMP DO
#else
    !$ACC PARALLEL LOOP GANG PRIVATE(enodes, el) DEFAULT(PRESENT) VECTOR_LENGTH(acc_vl)
#endif
    do edge=1, myDim_edge2D
        ! local indice of nodes that span up edge ed
        enodes=edges(:,edge)

        ! local index of element that contribute to edge
        el=edge_tri(:,edge)

        ! number of layers -1 at elem el(1)
        nl1=nlevels(el(1))-1

        ! index off surface layer in case of cavity !=1
        nu1=ulevels(el(1))

        ! edge_cross_dxdy(1:2,ed)... dx,dy distance from element centroid el(1) to
        ! center of edge --> needed to calc flux perpedicular to edge from elem el(1)
        deltaX1=edge_cross_dxdy(1,edge)
        deltaY1=edge_cross_dxdy(2,edge)
        a=r_earth*elem_cos(el(1))

        !_______________________________________________________________________
        ! same parameter but for other element el(2) that contributes to edge ed
        ! if el(2)==0 than edge is boundary edge
        nl2=0
        nu2=0
        if(el(2)>0) then
            deltaX2=edge_cross_dxdy(3,edge)
            deltaY2=edge_cross_dxdy(4,edge)
            ! number of layers -1 at elem el(2)
            nl2=nlevels(el(2))-1
            nu2=ulevels(el(2))
            a=0.5_WP*(a+r_earth*elem_cos(el(2)))
        end if

        !_______________________________________________________________________
        ! nl12 ... minimum number of layers -1 between element el(1) & el(2) that
        ! contribute to edge ed
        ! nu12 ... upper index of layers between element el(1) & el(2) that
        ! contribute to edge ed
        ! be carefull !!! --> if ed is a boundary edge than el(1)~=0 and el(2)==0
        !                     that means nl1>0, nl2==0, n2=min(nl1,nl2)=0 !!!
        nl12=min(nl1,nl2)
        nu12=max(nu1,nu2)

        !_______________________________________________________________________
        ! (A) goes only into this loop when the edge has only facing element
        ! el(1) --> so the edge is a boundary edge --> this is for ocean
        ! surface in case of cavity
        !$ACC LOOP VECTOR
        do nz=nu1, nu12-1
           !____________________________________________________________________
           ! volume flux across the segments
           vflux=(-VEL(2,nz,el(1))*deltaX1 + VEL(1,nz,el(1))*deltaY1)*helem(nz,el(1))

           !____________________________________________________________________
           ! 1st. low order upwind solution
           flux(nz, edge)=-0.5_WP*(                                                     &
                         ttf(nz, enodes(1))*(vflux+abs(vflux))+ &
                         ttf(nz, enodes(2))*(vflux-abs(vflux))  &
                         )-flux(nz, edge)
        end do
        !$ACC END LOOP

        !_______________________________________________________________________
        ! (B) goes only into this loop when the edge has only facing elemenmt
        ! el(2) --> so the edge is a boundary edge --> this is for ocean
        ! surface in case of cavity
        if (nu2 > 0) then
            !$ACC LOOP VECTOR
            do nz=nu2, nu12-1
                !___________________________________________________________
                ! volume flux across the segments
                vflux=(VEL(2,nz,el(2))*deltaX2 - VEL(1,nz,el(2))*deltaY2)*helem(nz,el(2))

                !___________________________________________________________
                ! 1st. low order upwind solution
                flux(nz, edge)=-0.5_WP*(                                           &
                            ttf(nz, enodes(1))*(vflux+abs(vflux))+ &
                            ttf(nz, enodes(2))*(vflux-abs(vflux)))-flux(nz, edge)
            end do
            !$ACC END LOOP
        end if

        !_______________________________________________________________________
        ! (C) Both segments
        ! loop over depth layers from top (nu12) to nl12
        ! be carefull !!! --> if ed is a boundary edge, el(2)==0 than nl12=0 so
        !                     you wont enter in this loop
        !$ACC LOOP VECTOR
        do nz=nu12, nl12
            !___________________________________________________________________
            ! 1st. low order upwind solution
            ! here already assumed that ed is NOT! a boundary edge so el(2) should exist
            vflux=(-VEL(2,nz,el(1))*deltaX1 + VEL(1,nz,el(1))*deltaY1)*helem(nz,el(1)) &
                  +(VEL(2,nz,el(2))*deltaX2 - VEL(1,nz,el(2))*deltaY2)*helem(nz,el(2))

            flux(nz, edge)=-0.5_WP*(                                                     &
                         ttf(nz, enodes(1))*(vflux+abs(vflux))+ &
                         ttf(nz, enodes(2))*(vflux-abs(vflux)))-flux(nz, edge)
        end do
        !$ACC END LOOP

        !_______________________________________________________________________
        ! (D) remaining segments on the left or on the right
        !$ACC LOOP VECTOR
        do nz=nl12+1, nl1
           !____________________________________________________________________
           ! volume flux across the segments
           vflux=(-VEL(2,nz,el(1))*deltaX1 + VEL(1,nz,el(1))*deltaY1)*helem(nz,el(1))
           !____________________________________________________________________
           ! 1st. low order upwind solution
           flux(nz, edge)=-0.5_WP*(                                                     &
                         ttf(nz, enodes(1))*(vflux+abs(vflux))+ &
                         ttf(nz, enodes(2))*(vflux-abs(vflux))  &
                         )-flux(nz, edge)
        end do
        !$ACC END LOOP

        !_______________________________________________________________________
        ! (E) remaining segments on the left or on the right
        !$ACC LOOP VECTOR
        do nz=nl12+1, nl2
                !_______________________________________________________________
                ! volume flux across the segments
                vflux=(VEL(2,nz,el(2))*deltaX2 - VEL(1,nz,el(2))*deltaY2)*helem(nz,el(2))
                !_______________________________________________________________
                ! 1st. low order upwind solution
                flux(nz, edge)=-0.5_WP*(                                           &
                             ttf(nz, enodes(1))*(vflux+abs(vflux))+ &
                             ttf(nz, enodes(2))*(vflux-abs(vflux)))-flux(nz, edge)
        end do
        !$ACC END LOOP
    end do
#ifndef ENABLE_OPENACC
!$OMP END DO
!$OMP END PARALLEL
#else
    !$ACC END PARALLEL LOOP
#endif
end subroutine adv_tra_hor_upw1
!
!
!===============================================================================
subroutine adv_tra_hor_muscl(vel, ttf, partit, mesh, num_ord, flux, edge_up_dn_grad, nboundary_lay, o_init_zero)
    use MOD_MESH
    use MOD_TRACER
    USE MOD_PARTIT
    USE MOD_PARSUP
    use g_comm_auto
    implicit none
    type(t_partit),intent(in), target :: partit
    type(t_mesh),  intent(in), target :: mesh
    real(kind=WP), intent(in)         :: num_ord    ! num_ord is the fraction of fourth-order contribution in the solution
    real(kind=WP), intent(in)         :: ttf(   mesh%nl-1, partit%myDim_nod2D+partit%eDim_nod2D)
    real(kind=WP), intent(in)         :: vel(2, mesh%nl-1, partit%myDim_elem2D+partit%eDim_elem2D)
    real(kind=WP), intent(inout)      :: flux(  mesh%nl-1, partit%myDim_edge2D)
    integer,       intent(in)         :: nboundary_lay(partit%myDim_nod2D+partit%eDim_nod2D)
    real(kind=WP), intent(in)         :: edge_up_dn_grad(4, mesh%nl-1, partit%myDim_edge2D)
    logical, optional                 :: o_init_zero
    logical                           :: l_init_zero
    real(kind=WP)                     :: deltaX1, deltaY1, deltaX2, deltaY2
    real(kind=WP)                     :: Tmean1, Tmean2, cHO
    real(kind=WP)                     :: c_lo(2)
    real(kind=WP)                     :: a, vflux
    integer                           :: el(2), enodes(2), nz, edge
    integer                           :: nu12, nl12, nl1, nl2, nu1, nu2

#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"

    l_init_zero=.true.
    if (present(o_init_zero)) then
       l_init_zero=o_init_zero
    end if
    if (l_init_zero) then
!$OMP PARALLEL DO
       do edge=1, myDim_edge2D
          flux(:,edge)=0.0_WP
       end do
!$OMP END PARALLEL DO
    end if

    ! The result is the low-order solution horizontal fluxes
    ! They are put into flux
    !___________________________________________________________________________
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(edge, deltaX1, deltaY1, deltaX2, deltaY2, Tmean1, Tmean2, cHO, &
!$OMP                                      c_lo, a, vflux, el, enodes, nz, nu12, nl12, nl1, nl2, nu1, nu2)
!$OMP DO
    do edge=1, myDim_edge2D
        ! local indice of nodes that span up edge ed
        enodes=edges(:,edge)

        ! local index of element that contribute to edge
        el=edge_tri(:,edge)

        ! number of layers -1 at elem el(1)
        nl1=nlevels(el(1))-1

        ! index off surface layer in case of cavity !=1
        nu1=ulevels(el(1))

        ! edge_cross_dxdy(1:2,ed)... dx,dy distance from element centroid el(1) to
        ! center of edge --> needed to calc flux perpedicular to edge from elem el(1)
        deltaX1=edge_cross_dxdy(1,edge)
        deltaY1=edge_cross_dxdy(2,edge)
        a=r_earth*elem_cos(el(1))

        !_______________________________________________________________________
        ! same parameter but for other element el(2) that contributes to edge ed
        ! if el(2)==0 than edge is boundary edge
        nl2=0
        nu2=0
        if(el(2)>0) then
            deltaX2=edge_cross_dxdy(3,edge)
            deltaY2=edge_cross_dxdy(4,edge)
            ! number of layers -1 at elem el(2)
            nl2=nlevels(el(2))-1
            nu2=ulevels(el(2))
            a=0.5_WP*(a+r_earth*elem_cos(el(2)))
        end if

        !_______________________________________________________________________
        ! n2 ... minimum number of layers -1 between element el(1) & el(2) that
        ! contribute to edge ed
        ! nu12 ... upper index of layers between element el(1) & el(2) that
        ! contribute to edge ed
        ! be carefull !!! --> if ed is a boundary edge than el(1)~=0 and el(2)==0
        !                     that means nl1>0, nl2==0, n2=min(nl1,nl2)=0 !!!
        nl12=min(nl1,nl2)
        nu12=max(nu1,nu2)

        !_______________________________________________________________________
        ! (A) goes only into this loop when the edge has only facing element
        ! el(1) --> so the edge is a boundary edge --> this is for ocean
        ! surface in case of cavity
        do nz=nu1, nu12-1
           c_lo(1)=real(max(sign(1, nboundary_lay(enodes(1))-nz), 0),WP)
           c_lo(2)=real(max(sign(1, nboundary_lay(enodes(2))-nz), 0),WP)

           !____________________________________________________________________
           Tmean2=ttf(nz, enodes(2))- &
                  (2.0_WP*(ttf(nz, enodes(2))-ttf(nz,enodes(1)))+ &
                  edge_dxdy(1,edge)*a*edge_up_dn_grad(2,nz,edge)+ &
                  edge_dxdy(2,edge)*r_earth*edge_up_dn_grad(4,nz,edge))/6.0_WP*c_lo(2)

           Tmean1=ttf(nz, enodes(1))+ &
                  (2.0_WP*(ttf(nz, enodes(2))-ttf(nz,enodes(1)))+ &
                  edge_dxdy(1,edge)*a*edge_up_dn_grad(1,nz,edge)+ &
                  edge_dxdy(2,edge)*r_earth*edge_up_dn_grad(3,nz,edge))/6.0_WP*c_lo(1)

           !____________________________________________________________________
           ! volume flux across the segments
           vflux=(-VEL(2,nz,el(1))*deltaX1 + VEL(1,nz,el(1))*deltaY1)*helem(nz,el(1))
           cHO=(vflux+abs(vflux))*Tmean1 + (vflux-abs(vflux))*Tmean2
           flux(nz,edge)=-0.5_WP*(1.0_WP-num_ord)*cHO - vflux*num_ord*0.5_WP*(Tmean1+Tmean2)-flux(nz,edge)
        end do

        !_______________________________________________________________________
        ! (B) goes only into this loop when the edge has only facing elemenmt
        ! el(2) --> so the edge is a boundary edge --> this is for ocean
        ! surface in case of cavity
        if (nu2 > 0) then
            do nz=nu2, nu12-1
                c_lo(1)=real(max(sign(1, nboundary_lay(enodes(1))-nz), 0),WP)
                c_lo(2)=real(max(sign(1, nboundary_lay(enodes(2))-nz), 0),WP)

                !_______________________________________________________________
                Tmean2=ttf(nz, enodes(2))- &
                        (2.0_WP*(ttf(nz, enodes(2))-ttf(nz,enodes(1)))+ &
                        edge_dxdy(1,edge)*a*edge_up_dn_grad(2,nz,edge)+ &
                        edge_dxdy(2,edge)*r_earth*edge_up_dn_grad(4,nz,edge))/6.0_WP*c_lo(2)

                Tmean1=ttf(nz, enodes(1))+ &
                        (2.0_WP*(ttf(nz, enodes(2))-ttf(nz,enodes(1)))+ &
                        edge_dxdy(1,edge)*a*edge_up_dn_grad(1,nz,edge)+ &
                        edge_dxdy(2,edge)*r_earth*edge_up_dn_grad(3,nz,edge))/6.0_WP*c_lo(1)

                !_______________________________________________________________
                ! volume flux across the segments
                vflux=(VEL(2,nz,el(2))*deltaX2 - VEL(1,nz,el(2))*deltaY2)*helem(nz,el(2))
                cHO=(vflux+abs(vflux))*Tmean1 + (vflux-abs(vflux))*Tmean2
                flux(nz,edge)=-0.5_WP*(1.0_WP-num_ord)*cHO - vflux*num_ord*0.5_WP*(Tmean1+Tmean2)-flux(nz,edge)
            end do
        end if

        !_______________________________________________________________________
        ! (C) Both segments
        ! loop over depth layers from top to n2
        ! be carefull !!! --> if ed is a boundary edge, el(2)==0 than n2=0 so
        !                     you wont enter in this loop
        do nz=nu12, nl12
            c_lo(1)=real(max(sign(1, nboundary_lay(enodes(1))-nz), 0),WP)
            c_lo(2)=real(max(sign(1, nboundary_lay(enodes(2))-nz), 0),WP)

           !___________________________________________________________________
           ! MUSCL-type reconstruction
           ! check if upwind or downwind triagle is necessary
           !
           ! cross product between velocity vector and cross vector edge-elem-center
           ! cross product > 0 --> angle vec_v and (dx,dy) --> [0   180] --> upwind triangle
           ! cross product < 0 --> angle vec_v and (dx,dy) --> [180 360] --> downwind triangle
           !
           !       o                  o      !     o                  o
           !      / \                / \     !    / \                / \
           !     /   \    \ vec_v   /   \    !   /   \        /     /   \
           !    /  up \    \       / dn  \   !  /  up \      /     / dn  \
           !   o-------o----+---->o-------o  ! o-------o----+---->o-------o
           !           1   /      2          !         1     \vec_v
           !              /vec_v             !                \
           !   --> downwind triangle         ! --> upwind triangle
           !
           !  edge_up_dn_grad(1,nz,edge) ... gradTR_x upwind
           !  edge_up_dn_grad(2,nz,edge) ... gradTR_x downwind
           !  edge_up_dn_grad(3,nz,edge) ... gradTR_y upwind
           !  edge_up_dn_grad(4,nz,edge) ... gradTR_y downwind

           !___________________________________________________________________
           ! use downwind triangle to interpolate Tracer to edge center with
           ! fancy scheme --> Linear upwind reconstruction
           ! T_n+0.5 = T_n+1 - 1/2*deltax*GRADIENT
           ! --> GRADIENT = 2/3 GRAD_edgecenter + 1/3 GRAD_downwindtri
           ! T_n+0.5 = T_n+1 - 2/6*(T_n+1-T_n) + 1/6*gradT_down
           ! --> edge_up_dn_grad ... contains already elemental tracer gradient
           !     of up and dn wind triangle
           ! --> Tmean2 ... edge center interpolated Tracer using tracer
           !     gradient info from upwind triangle
           Tmean2=ttf(nz, enodes(2))- &
                  (2.0_WP*(ttf(nz, enodes(2))-ttf(nz,enodes(1)))+ &
                  edge_dxdy(1,edge)*a*edge_up_dn_grad(2,nz,edge)+ &
                  edge_dxdy(2,edge)*r_earth*edge_up_dn_grad(4,nz,edge))/6.0_WP*c_lo(2)

            ! use upwind triangle to interpolate Tracer to edge center with
            ! fancy scheme --> Linear upwind reconstruction
            ! T_n+0.5 = T_n + 1/2*deltax*GRADIENT
            ! --> GRADIENT = 2/3 GRAD_edgecenter + 1/3 GRAD_downwindtri
            ! T_n+0.5 = T_n + 2/6*(T_n+1-T_n) + 1/6*gradT_down
            ! --> Tmean1 ... edge center interpolated Tracer using tracer
            !     gradient info from downwind triangle
            Tmean1=ttf(nz, enodes(1))+ &
                   (2.0_WP*(ttf(nz, enodes(2))-ttf(nz,enodes(1)))+ &
                   edge_dxdy(1,edge)*a*edge_up_dn_grad(1,nz,edge)+ &
                   edge_dxdy(2,edge)*r_earth*edge_up_dn_grad(3,nz,edge))/6.0_WP*c_lo(1)

            !___________________________________________________________________
            ! volume flux along the edge segment ed
            ! netto volume flux along segment that comes from edge node 1 and 2
            !
            !
            !                         C1 (centroid el(1)) --> (u1,v1)
            !                         x
            !                         ^
            !               (dx1,dy1) |
            !                         |---> vec_n1 (dy1,-dx1)--> project vec_u1 onto vec_n1 --> -v1*dx1+u1*dy1 -->
            !                         |                                                                          |
            !    enodes(1) o----------O---------o enodes(2)                                                      |-> calculate volume flux out of/in
            !          vflux_________/|                                                                          |   the volume of enode1(enode2) through
            !                         |---> vec_n2 (dy2,-dx2)--> project vec_u2 onto vec_n2 --> -v2*dx2+u2*dy2 -->   sections of dx1,dy1 and dx2,dy2
            !               (dx2,dy2) |                                                                              --> vflux
            !                         v
            !                         x
            !                         C2 (centroid el(2)) --> (u2,v2)

            ! here already assumed that ed is NOT! a boundary edge so el(2) should exist
            vflux=(-VEL(2,nz,el(1))*deltaX1 + VEL(1,nz,el(1))*deltaY1)*helem(nz,el(1)) &
                  +(VEL(2,nz,el(2))*deltaX2 - VEL(1,nz,el(2))*deltaY2)*helem(nz,el(2))

            !___________________________________________________________________
            ! (1-num_ord) is done with 3rd order upwind
            cHO=(vflux+abs(vflux))*Tmean1 + (vflux-abs(vflux))*Tmean2
            flux(nz,edge)=-0.5_WP*(1.0_WP-num_ord)*cHO - vflux*num_ord*0.5_WP*(Tmean1+Tmean2)-flux(nz,edge)
        end do

        !_______________________________________________________________________
        ! (D) remaining segments on the left or on the right
        do nz=nl12+1, nl1
           c_lo(1)=real(max(sign(1, nboundary_lay(enodes(1))-nz), 0),WP)
           c_lo(2)=real(max(sign(1, nboundary_lay(enodes(2))-nz), 0),WP)

           !____________________________________________________________________
           Tmean2=ttf(nz, enodes(2))- &
                  (2.0_WP*(ttf(nz, enodes(2))-ttf(nz,enodes(1)))+ &
                  edge_dxdy(1,edge)*a*edge_up_dn_grad(2,nz,edge)+ &
                  edge_dxdy(2,edge)*r_earth*edge_up_dn_grad(4,nz,edge))/6.0_WP*c_lo(2)

           Tmean1=ttf(nz, enodes(1))+ &
                  (2.0_WP*(ttf(nz, enodes(2))-ttf(nz,enodes(1)))+ &
                  edge_dxdy(1,edge)*a*edge_up_dn_grad(1,nz,edge)+ &
                  edge_dxdy(2,edge)*r_earth*edge_up_dn_grad(3,nz,edge))/6.0_WP*c_lo(1)

           !____________________________________________________________________
           ! volume flux across the segments
           vflux=(-VEL(2,nz,el(1))*deltaX1 + VEL(1,nz,el(1))*deltaY1)*helem(nz,el(1))
           cHO=(vflux+abs(vflux))*Tmean1 + (vflux-abs(vflux))*Tmean2
           flux(nz,edge)=-0.5_WP*(1.0_WP-num_ord)*cHO - vflux*num_ord*0.5_WP*(Tmean1+Tmean2)-flux(nz,edge)
        end do

        !_______________________________________________________________________
        ! (E) remaining segments on the left or on the right
        do nz=nl12+1, nl2
           c_lo(1)=real(max(sign(1, nboundary_lay(enodes(1))-nz), 0),WP)
           c_lo(2)=real(max(sign(1, nboundary_lay(enodes(2))-nz), 0),WP)

           !____________________________________________________________________
           Tmean2=ttf(nz, enodes(2))- &
                  (2.0_WP*(ttf(nz, enodes(2))-ttf(nz,enodes(1)))+ &
                  edge_dxdy(1,edge)*a*edge_up_dn_grad(2,nz,edge)+ &
                  edge_dxdy(2,edge)*r_earth*edge_up_dn_grad(4,nz,edge))/6.0_WP*c_lo(2)

           Tmean1=ttf(nz, enodes(1))+ &
                  (2.0_WP*(ttf(nz, enodes(2))-ttf(nz,enodes(1)))+ &
                  edge_dxdy(1,edge)*a*edge_up_dn_grad(1,nz,edge)+ &
                  edge_dxdy(2,edge)*r_earth*edge_up_dn_grad(3,nz,edge))/6.0_WP*c_lo(1)

           !____________________________________________________________________
           ! volume flux across the segments
           vflux=(VEL(2,nz,el(2))*deltaX2 - VEL(1,nz,el(2))*deltaY2)*helem(nz,el(2))
           cHO=(vflux+abs(vflux))*Tmean1 + (vflux-abs(vflux))*Tmean2
           flux(nz,edge)=-0.5_WP*(1.0_WP-num_ord)*cHO - vflux*num_ord*0.5_WP*(Tmean1+Tmean2)-flux(nz,edge)
        end do
    end do
!$OMP END DO
!$OMP END PARALLEL
end subroutine adv_tra_hor_muscl
!
!
!===============================================================================
    subroutine adv_tra_hor_mfct(vel, ttf, partit, mesh, num_ord, flux, edge_up_dn_grad,                 o_init_zero)
    use MOD_MESH
    use MOD_TRACER
    USE MOD_PARTIT
    USE MOD_PARSUP
    use g_comm_auto
    implicit none
    type(t_partit),intent(inout), target :: partit
    type(t_mesh),  intent(in), target :: mesh
    real(kind=WP), intent(in)         :: num_ord    ! num_ord is the fraction of fourth-order contribution in the solution
    real(kind=WP), intent(in)         :: ttf(   mesh%nl-1, partit%myDim_nod2D+partit%eDim_nod2D)
    real(kind=WP), intent(in)         :: vel(2, mesh%nl-1, partit%myDim_elem2D+partit%eDim_elem2D)
    real(kind=WP), intent(inout)      :: flux(  mesh%nl-1, partit%myDim_edge2D)
    real(kind=WP), intent(in)         :: edge_up_dn_grad(4, mesh%nl-1, partit%myDim_edge2D)
    logical, optional                 :: o_init_zero
    logical                           :: l_init_zero
    real(kind=WP)                     :: deltaX1, deltaY1, deltaX2, deltaY2
    real(kind=WP)                     :: Tmean1, Tmean2, cHO
    real(kind=WP)                     :: a, vflux
    integer                           :: el(2), enodes(2), nz, edge
    integer                           :: nu12, nl12, nl1, nl2, nu1, nu2

#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"

    l_init_zero=.true.
    if (present(o_init_zero)) then
       l_init_zero=o_init_zero
    end if
    if (l_init_zero) then
#ifndef ENABLE_OPENACC
!$OMP PARALLEL DO
#else
       !$ACC PARALLEL LOOP GANG DEFAULT(PRESENT) VECTOR_LENGTH(acc_vl)
#endif
       do edge=1, myDim_edge2D
          flux(:,edge)=0.0_WP
       end do
#ifndef ENABLE_OPENACC
!$OMP END PARALLEL DO
#else
       !$ACC END PARALLEL LOOP
#endif
    end if

    ! The result is the low-order solution horizontal fluxes
    ! They are put into flux
    !___________________________________________________________________________
#ifndef ENABLE_OPENACC
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(edge, deltaX1, deltaY1, deltaX2, deltaY2, Tmean1, Tmean2, cHO, &
!$OMP                                     a, vflux, el, enodes, nz, nu12, nl12, nl1, nl2, nu1, nu2)
!$OMP DO
#else
    !$ACC PARALLEL LOOP GANG PRIVATE(enodes, el) DEFAULT(PRESENT) VECTOR_LENGTH(acc_vl)
#endif
    do edge=1, myDim_edge2D
        ! local indice of nodes that span up edge ed
        enodes=edges(:,edge)

        ! local index of element that contribute to edge
        el=edge_tri(:,edge)

        ! number of layers -1 at elem el(1)
        nl1=nlevels(el(1))-1

        ! index off surface layer in case of cavity !=1
        nu1=ulevels(el(1))

        ! edge_cross_dxdy(1:2,ed)... dx,dy distance from element centroid el(1) to
        ! center of edge --> needed to calc flux perpedicular to edge from elem el(1)
        deltaX1=edge_cross_dxdy(1,edge)
        deltaY1=edge_cross_dxdy(2,edge)
        a=r_earth*elem_cos(el(1))

        !_______________________________________________________________________
        ! same parameter but for other element el(2) that contributes to edge ed
        ! if el(2)==0 than edge is boundary edge
        nl2=0
        nu2=0
        if(el(2)>0) then
            deltaX2=edge_cross_dxdy(3,edge)
            deltaY2=edge_cross_dxdy(4,edge)
            ! number of layers -1 at elem el(2)
            nl2=nlevels(el(2))-1
            nu2=ulevels(el(2))
            a=0.5_WP*(a+r_earth*elem_cos(el(2)))
        end if

        !_______________________________________________________________________
        ! n2 ... minimum number of layers -1 between element el(1) & el(2) that
        ! contribute to edge ed
        ! nu12 ... upper index of layers between element el(1) & el(2) that
        ! contribute to edge ed
        ! be carefull !!! --> if ed is a boundary edge than el(1)~=0 and el(2)==0
        !                     that means nl1>0, nl2==0, n2=min(nl1,nl2)=0 !!!
        nl12=min(nl1,nl2)
        nu12=max(nu1,nu2)

        !_______________________________________________________________________
        ! (A) goes only into this loop when the edge has only facing element
        ! el(1) --> so the edge is a boundary edge --> this is for ocean
        ! surface in case of cavity
        !$ACC LOOP VECTOR
        do nz=nu1, nu12-1
           !____________________________________________________________________
           Tmean2=ttf(nz, enodes(2))- &
                  (2.0_WP*(ttf(nz, enodes(2))-ttf(nz,enodes(1)))+ &
                  edge_dxdy(1,edge)*a*edge_up_dn_grad(2,nz,edge)+ &
                  edge_dxdy(2,edge)*r_earth*edge_up_dn_grad(4,nz,edge))/6.0_WP

           Tmean1=ttf(nz, enodes(1))+ &
                  (2.0_WP*(ttf(nz, enodes(2))-ttf(nz,enodes(1)))+ &
                  edge_dxdy(1,edge)*a*edge_up_dn_grad(1,nz,edge)+ &
                  edge_dxdy(2,edge)*r_earth*edge_up_dn_grad(3,nz,edge))/6.0_WP

           !____________________________________________________________________
           ! volume flux across the segments
           vflux=(-VEL(2,nz,el(1))*deltaX1 + VEL(1,nz,el(1))*deltaY1)*helem(nz,el(1))
           cHO=(vflux+abs(vflux))*Tmean1 + (vflux-abs(vflux))*Tmean2
           flux(nz,edge)=-0.5_WP*(1.0_WP-num_ord)*cHO - vflux*num_ord*0.5_WP*(Tmean1+Tmean2)-flux(nz,edge)
        end do
        !$ACC END LOOP

        !_______________________________________________________________________
        ! (B) goes only into this loop when the edge has only facing elemenmt
        ! el(2) --> so the edge is a boundary edge --> this is for ocean
        ! surface in case of cavity
        if (nu2 > 0) then
            !$ACC LOOP VECTOR
            do nz=nu2,nu12-1
            !___________________________________________________________________
            Tmean2=ttf(nz, enodes(2))- &
                    (2.0_WP*(ttf(nz, enodes(2))-ttf(nz,enodes(1)))+ &
                    edge_dxdy(1,edge)*a*edge_up_dn_grad(2,nz,edge)+ &
                    edge_dxdy(2,edge)*r_earth*edge_up_dn_grad(4,nz,edge))/6.0_WP

            Tmean1=ttf(nz, enodes(1))+ &
                    (2.0_WP*(ttf(nz, enodes(2))-ttf(nz,enodes(1)))+ &
                    edge_dxdy(1,edge)*a*edge_up_dn_grad(1,nz,edge)+ &
                    edge_dxdy(2,edge)*r_earth*edge_up_dn_grad(3,nz,edge))/6.0_WP
            !___________________________________________________________________
            ! volume flux across the segments
            vflux=(VEL(2,nz,el(2))*deltaX2 - VEL(1,nz,el(2))*deltaY2)*helem(nz,el(2))
            cHO=(vflux+abs(vflux))*Tmean1 + (vflux-abs(vflux))*Tmean2
            flux(nz,edge)=-0.5_WP*(1.0_WP-num_ord)*cHO - vflux*num_ord*0.5_WP*(Tmean1+Tmean2)-flux(nz,edge)
            end do
            !$ACC END LOOP
        end if

        !_______________________________________________________________________
        ! (C) Both segments
        ! loop over depth layers from top to n2
        ! be carefull !!! --> if ed is a boundary edge, el(2)==0 than n2=0 so
        !                     you wont enter in this loop
        !$ACC LOOP VECTOR
        do nz=nu12, nl12
           !___________________________________________________________________
           ! MUSCL-type reconstruction
           ! check if upwind or downwind triagle is necessary
           !
           ! cross product between velocity vector and cross vector edge-elem-center
           ! cross product > 0 --> angle vec_v and (dx,dy) --> [0   180] --> upwind triangle
           ! cross product < 0 --> angle vec_v and (dx,dy) --> [180 360] --> downwind triangle
           !
           !       o                  o      !     o                  o
           !      / \                / \     !    / \                / \
           !     /   \    \ vec_v   /   \    !   /   \        /     /   \
           !    /  up \    \       / dn  \   !  /  up \      /     / dn  \
           !   o-------o----+---->o-------o  ! o-------o----+---->o-------o
           !           1   /      2          !         1     \vec_v
           !              /vec_v             !                \
           !   --> downwind triangle         ! --> upwind triangle
           !
           !  edge_up_dn_grad(1,nz,edge) ... gradTR_x upwind
           !  edge_up_dn_grad(2,nz,edge) ... gradTR_x downwind
           !  edge_up_dn_grad(3,nz,edge) ... gradTR_y upwind
           !  edge_up_dn_grad(4,nz,edge) ... gradTR_y downwind

           !___________________________________________________________________
           ! use downwind triangle to interpolate Tracer to edge center with
           ! fancy scheme --> Linear upwind reconstruction
           ! T_n+0.5 = T_n+1 - 1/2*deltax*GRADIENT
           ! --> GRADIENT = 2/3 GRAD_edgecenter + 1/3 GRAD_downwindtri
           ! T_n+0.5 = T_n+1 - 2/6*(T_n+1-T_n) + 1/6*gradT_down
           ! --> edge_up_dn_grad ... contains already elemental tracer gradient
           !     of up and dn wind triangle
           ! --> Tmean2 ... edge center interpolated Tracer using tracer
           !     gradient info from upwind triangle
           Tmean2=ttf(nz, enodes(2))- &
                  (2.0_WP*(ttf(nz, enodes(2))-ttf(nz,enodes(1)))+ &
                  edge_dxdy(1,edge)*a*edge_up_dn_grad(2,nz,edge)+ &
                  edge_dxdy(2,edge)*r_earth*edge_up_dn_grad(4,nz,edge))/6.0_WP

            ! use upwind triangle to interpolate Tracer to edge center with
            ! fancy scheme --> Linear upwind reconstruction
            ! T_n+0.5 = T_n + 1/2*deltax*GRADIENT
            ! --> GRADIENT = 2/3 GRAD_edgecenter + 1/3 GRAD_downwindtri
            ! T_n+0.5 = T_n + 2/6*(T_n+1-T_n) + 1/6*gradT_down
            ! --> Tmean1 ... edge center interpolated Tracer using tracer
            !     gradient info from downwind triangle
            Tmean1=ttf(nz, enodes(1))+ &
                   (2.0_WP*(ttf(nz, enodes(2))-ttf(nz,enodes(1)))+ &
                   edge_dxdy(1,edge)*a*edge_up_dn_grad(1,nz,edge)+ &
                   edge_dxdy(2,edge)*r_earth*edge_up_dn_grad(3,nz,edge))/6.0_WP
            !___________________________________________________________________
            ! volume flux along the edge segment ed
            ! netto volume flux along segment that comes from edge node 1 and 2
            !
            !
            !                         C1 (centroid el(1)) --> (u1,v1)
            !                         x
            !                         ^
            !               (dx1,dy1) |
            !                         |---> vec_n1 (dy1,-dx1)--> project vec_u1 onto vec_n1 --> -v1*dx1+u1*dy1 -->
            !                         |                                                                          |
            !    enodes(1) o----------O---------o enodes(2)                                                      |-> calculate volume flux out of/in
            !          vflux_________/|                                                                          |   the volume of enode1(enode2) through
            !                         |---> vec_n2 (dy2,-dx2)--> project vec_u2 onto vec_n2 --> -v2*dx2+u2*dy2 -->   sections of dx1,dy1 and dx2,dy2
            !               (dx2,dy2) |                                                                              --> vflux
            !                         v
            !                         x
            !                         C2 (centroid el(2)) --> (u2,v2)

            ! here already assumed that ed is NOT! a boundary edge so el(2) should exist
            vflux=(-VEL(2,nz,el(1))*deltaX1 + VEL(1,nz,el(1))*deltaY1)*helem(nz,el(1)) &
                  +(VEL(2,nz,el(2))*deltaX2 - VEL(1,nz,el(2))*deltaY2)*helem(nz,el(2))

            !___________________________________________________________________
            ! (1-num_ord) is done with 3rd order upwind
            cHO=(vflux+abs(vflux))*Tmean1 + (vflux-abs(vflux))*Tmean2
            flux(nz,edge)=-0.5_WP*(1.0_WP-num_ord)*cHO - vflux*num_ord*0.5_WP*(Tmean1+Tmean2)-flux(nz,edge)
        end do
        !$ACC END LOOP

        !_______________________________________________________________________
        ! (D) remaining segments on the left or on the right
        !$ACC LOOP VECTOR
        do nz=nl12+1, nl1
           !____________________________________________________________________
           Tmean2=ttf(nz, enodes(2))- &
                  (2.0_WP*(ttf(nz, enodes(2))-ttf(nz,enodes(1)))+ &
                  edge_dxdy(1,edge)*a*edge_up_dn_grad(2,nz,edge)+ &
                  edge_dxdy(2,edge)*r_earth*edge_up_dn_grad(4,nz,edge))/6.0_WP

           Tmean1=ttf(nz, enodes(1))+ &
                  (2.0_WP*(ttf(nz, enodes(2))-ttf(nz,enodes(1)))+ &
                  edge_dxdy(1,edge)*a*edge_up_dn_grad(1,nz,edge)+ &
                  edge_dxdy(2,edge)*r_earth*edge_up_dn_grad(3,nz,edge))/6.0_WP

           !____________________________________________________________________
           ! volume flux across the segments
           vflux=(-VEL(2,nz,el(1))*deltaX1 + VEL(1,nz,el(1))*deltaY1)*helem(nz,el(1))
           cHO=(vflux+abs(vflux))*Tmean1 + (vflux-abs(vflux))*Tmean2
           flux(nz,edge)=-0.5_WP*(1.0_WP-num_ord)*cHO - vflux*num_ord*0.5_WP*(Tmean1+Tmean2)-flux(nz,edge)
        end do
        !$ACC END LOOP

        !_______________________________________________________________________
        ! (E) remaining segments on the left or on the right
        !$ACC LOOP VECTOR
        do nz=nl12+1, nl2
           !____________________________________________________________________
           Tmean2=ttf(nz, enodes(2))- &
                  (2.0_WP*(ttf(nz, enodes(2))-ttf(nz,enodes(1)))+ &
                  edge_dxdy(1,edge)*a*edge_up_dn_grad(2,nz,edge)+ &
                  edge_dxdy(2,edge)*r_earth*edge_up_dn_grad(4,nz,edge))/6.0_WP

           Tmean1=ttf(nz, enodes(1))+ &
                  (2.0_WP*(ttf(nz, enodes(2))-ttf(nz,enodes(1)))+ &
                  edge_dxdy(1,edge)*a*edge_up_dn_grad(1,nz,edge)+ &
                  edge_dxdy(2,edge)*r_earth*edge_up_dn_grad(3,nz,edge))/6.0_WP

           !____________________________________________________________________
           ! volume flux across the segments
           vflux=(VEL(2,nz,el(2))*deltaX2 - VEL(1,nz,el(2))*deltaY2)*helem(nz,el(2))
           cHO=(vflux+abs(vflux))*Tmean1 + (vflux-abs(vflux))*Tmean2
           flux(nz,edge)=-0.5_WP*(1.0_WP-num_ord)*cHO - vflux*num_ord*0.5_WP*(Tmean1+Tmean2)-flux(nz,edge)
        end do
        !$ACC END LOOP
    end do
#ifndef ENABLE_OPENACC
!$OMP END DO
!$OMP END PARALLEL
#else
    !$ACC END PARALLEL LOOP
#endif 
end subroutine adv_tra_hor_mfct



!
!
!_______________________________________________________________________________
! horizontal advection second order in space (2nd oprder in time) using superbee
! slope limiter
subroutine adv_tra_hor_spbee(             &
            vel                         , & 
            ttf                         , &
            partit                      , &
            mesh                        , &
            num_ord                     , &
            flux                        , &
            edge_up_dn_grad             , &
            flag_posdef                 , &
            o_init_zero                   &
            )
    use MOD_MESH
    use MOD_TRACER
    use MOD_PARTIT
    use MOD_PARSUP
    use g_config, only: dt
    use g_comm_auto
    implicit none
    !___INPUT/OUTPUT VARIABLES__________________________________________________
    type(t_partit)  , intent(inout), target :: partit
    type(t_mesh)    , intent(in)   , target :: mesh
    real(kind=WP)   , intent(in)            :: ttf(               mesh%nl-1, partit%myDim_nod2D+partit%eDim_nod2D)
    real(kind=WP)   , intent(in)            :: vel(2,             mesh%nl-1, partit%myDim_elem2D+partit%eDim_elem2D)
    real(kind=WP)   , intent(inout)         :: flux(              mesh%nl-1, partit%myDim_edge2D)
    real(kind=WP)   , intent(in)            :: edge_up_dn_grad(4, mesh%nl-1, partit%myDim_edge2D)
    real(kind=WP)   , intent(in)            :: num_ord
    logical         , intent(in)            :: flag_posdef
    logical         , optional              :: o_init_zero
    !___LOCAL VARIABLES_________________________________________________________
    logical                                 :: l_init_zero, flag_2ndord_time
    real(kind=WP)                           :: dx1, dy1, dx2, dy2, dx0, dy0, dxdy12(2), dt_over_edlen
    real(kind=WP)                           :: n_x, n_y, nlen, inv_nlen
    real(kind=WP)                           :: u1, u2, v1, v2, Ue
    real(kind=WP)                           :: vflux, cfl, T12vflux
    integer                                 :: edel(2), ednodes(2), edge, edel0
    integer                                 :: nz, nl1, nl2, nl12, nu1, nu2, nu12, nzs, nze
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"
    
    !___________________________________________________________________________
    if (num_ord == 1.0_WP) then
        flag_2ndord_time = .True.
    else
        flag_2ndord_time = .False.
    end if 
    
    !___________________________________________________________________________
    l_init_zero=.true.
    if (present(o_init_zero)) then
       l_init_zero=o_init_zero
    end if
    if (l_init_zero) then
#ifndef ENABLE_OPENACC
        !$OMP PARALLEL DO
#else
        !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) DEFAULT(PRESENT) VECTOR_LENGTH(acc_vl)
#endif
       do edge=1, myDim_edge2D
          do nz=1, mesh%nl-1
             flux(nz,edge)=0.0_WP
          end do
       end do
#ifndef ENABLE_OPENACC
        !$OMP END PARALLEL DO
#else
        !$ACC END PARALLEL LOOP
#endif
    end if
    
    !___________________________________________________________________________
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(edge, ednodes, edel, edel0, &
!$OMP                                  nz, nl1, nl2, nl12, nu1, nu2, nu12, nzs, nze, &
!$OMP                                  dx1, dy1, dx2, dy2, dxdy12, dx0, dy0, &
!$OMP                                  vflux, cfl, dt_over_edlen, T12vflux,  &
!$OMP                                  u1, u2, v1, v2, n_x, n_y, Ue)
!$OMP DO        
    do edge=1, myDim_edge2D
        !_______________________________________________________________________
        ! local indice of nodes that span up edge ed
        ednodes= edges(:,edge)
            
        ! local index of element that contribute to edge
        edel     = edge_tri(:,edge)
        
        nl1      = nlevels(edel(1))-1
        ! nu1 ... upper index of ocean default = 1 but can consider cavity !=1
        nu1      = ulevels(edel(1))
        
        ! edge_cross_dxdy(1:2,ed)... dx,dy distance from center of edge to
        ! element centroid el(1) --> needed to calc flux perpedicular to edge from elem el(1)
        dx1    = edge_cross_dxdy(1,edge)
        dy1    = edge_cross_dxdy(2,edge)
            
        ! length of edge dx dy
        dxdy12 = edge_dxdy(:,edge)*r_earth  
            
        !_______________________________________________________________________
        ! same parameter but for other element el(2) that contributes to edge ed
        ! if el(2)==0 than edge is boundary edge
        nl2    = 0
        nu2    = 0
        dx2    = 0.0_WP
        dy2    = 0.0_WP
        if(edel(2)>0) then
            nl2      = nlevels(edel(2))-1
            nu2      = ulevels(edel(2))            
            dx2      = edge_cross_dxdy(3,edge)
            dy2      = edge_cross_dxdy(4,edge)
            dxdy12(1)= dxdy12(1) * (elem_cos(edel(1)) + elem_cos(edel(2))) * 0.5_WP
                
            ! compute mean edge-segment normal vector
            !                         (dx1,dy1)
            !                             ^ 
            !                             │
            !                             ├──> (dy1,-dx1)  
            !                             │
            !                 ●━━━━━━━━━━━┿━━━━━━━━━━►● 
            !                             │
            !               (dy2,-dx2) <──┼──> (-dy2,dx2)  
            !                             │ 
            !                             v  
            !                         (dx2,dy2)
            ! 
            n_x      = (  dy1  + (-dy2))/2
            n_y      = ((-dx1) +   dx2 )/2
                
        else
            dxdy12(1)= dxdy12(1) * elem_cos(edel(1))   
            ! compute edge-segment normal vector
            n_x      =  dy1
            n_y      = -dx1
        end if
        
        ! normalize n_vec=(n_x, n_y)
        nlen     = sqrt(n_x**2 + n_y**2)
        inv_nlen = 1/max(nlen, tiny(1.0_WP))
        n_x      = n_x * inv_nlen
        n_y      = n_y * inv_nlen
        
        dt_over_edlen = dt/sqrt(dxdy12(1)**2 + dxdy12(2)**2)
            
        !_______________________________________________________________________
        ! compute volume fluxes
        nu12=max(nu1,nu2)
        nl12=min(nl1,nl2)
        
        !_______________________________________________________________________
        ! loop over upper uncommon edge levels between edel(1) and edel(2)
        ! --> edge is boundary edge, edel(2) does not exist
        !     In this case the entire advection of this edge is only handled by the first
        !     loop the other two loop will be skiped
        if (edel(2) <= 0) then
            nzs   = nu1
            nze   = nl1
            edel0 = edel(1)
            dx0   = -dx1
            dy0   =  dy1
        ! --> edge is interior edge
        elseif (nu1 < nu2) then
            nzs   = nu1
            nze   = nu12-1
            edel0 = edel(1)
            dx0   = -dx1
            dy0   =  dy1
        ! --> edge is interior edge
        else ! (nu1 > nu2)
            nzs   = nu2
            nze   = nu12-1
            edel0 = edel(2)
            dx0   =  dx2
            dy0   = -dy2
        end if
        do nz=nzs, nze
            ! compute volume flux across the segments from el(1)
            u1        = vel(1, nz, edel0)
            v1        = vel(2, nz, edel0)
            vflux     = (v1*dx0 + u1*dy0)*helem(nz,edel0)
            
            ! compute approximated edge centered and along edge projected 
            ! mean velocity --> need to add component to make second order
            ! in time
            Ue        = abs( (u1 * n_x) + (v1 * n_y) )
            cfl       = Ue*dt_over_edlen
            
            ! compute upwind, downwind edge centered tracers volume flux, based on 
            ! gradient projected ttfm1 and ttfp2 values and applied superbee limiter
            T12vflux  = spbee_compute_T12vflux(ttf(nz, ednodes(1)), ttf(nz, ednodes(2)),  &
                                                edge_up_dn_grad(:, nz, edge), & 
                                                vflux, cfl, dxdy12, flag_2ndord_time, flag_posdef)
                                                
            flux(nz, edge)= -T12vflux - flux(nz, edge)
        end do ! --> do nz=nu0, nu12-1
        
        !_______________________________________________________________________
        ! if edge is boundary edge skip this code block everything is handled by the 
        ! first loop
        if (edel(2)>0) then
            !_______________________________________________________________________
            ! loop over common edge depth levels between edel(1) and edel(2)
            do nz = nu12, nl12
                ! compute volume flux across the segments from el(1)
                u1        = vel(1, nz, edel(1))
                v1        = vel(2, nz, edel(1))
                u2        = vel(1, nz, edel(2))
                v2        = vel(2, nz, edel(2))  
                vflux     =  (-v1*dx1 + u1*dy1)*helem(nz,edel(1)) + ( v2*dx2 - u2*dy2)*helem(nz,edel(2))
                
                ! compute approximated edge centered and along edge projected 
                ! mean velocity --> need to add component to make second order
                ! in time
                Ue        = abs( 0.5_WP*( u1 + u2 )*n_x + 0.5_WP*( v1 + v2 )*n_y )
                cfl       = Ue*dt_over_edlen
                
                ! compute upwind, downwind edge centered tracers volume flux, based on 
                ! gradient projected ttfm1 and ttfp2 values and applied superbee limiter
                T12vflux  = spbee_compute_T12vflux(ttf(nz, ednodes(1)), ttf(nz, ednodes(2)), &
                                                edge_up_dn_grad(:, nz, edge), & 
                                                vflux, cfl, dxdy12, flag_2ndord_time, flag_posdef)
                                                
                flux(nz, edge)= -T12vflux - flux(nz, edge)
            end do !--> do nz=nu12, nl12
            
            !_______________________________________________________________________
            ! loop over lower uncommon edge levels between edel(1) and edel(2)
            if (nl1 > nl2) then
                nzs   = nl12+1 
                nze   = nl1
                edel0 = edel(1)
                dx0   = -dx1
                dy0   =  dy1
            else ! (nl1 < nl2)
                nzs   = nl12+1 
                nze   = nl2
                edel0 = edel(2)
                dx0   =  dx2
                dy0   = -dy2
            end if
            do nz = nzs, nze
                ! compute volume flux across the segments from el(1)
                u1        = vel(1, nz, edel0)
                v1        = vel(2, nz, edel0)
                vflux     = (v1*dx0 + u1*dy0)*helem(nz,edel0)
                
                ! compute approximated edge centered and along edge projected 
                ! mean velocity --> need to add component to make second order
                ! in time
                Ue        = abs( u1 * n_x + v1 * n_y)
                cfl       = Ue*dt_over_edlen
                
                ! compute upwind, downwind edge centered tracers volume flux, based on 
                ! gradient projected ttfm1 and ttfp2 values and applied superbee limiter
                T12vflux  = spbee_compute_T12vflux(ttf(nz, ednodes(1)), ttf(nz, ednodes(2)), &
                                                    edge_up_dn_grad(:, nz, edge), & 
                                                    vflux, cfl, dxdy12, flag_2ndord_time, flag_posdef)
                                                
                flux(nz, edge)= -T12vflux - flux(nz, edge)
            end do ! --> do nz=nl12+1,nl1
        end if ! --> if (edel(2)>0) then
    end do !--> do edge=1, myDim_edge2D
!$OMP END DO
!$OMP END PARALLEL

    contains 
    
    !
    !
    !___________________________________________________________________________
    ! superbee slope limiter
    pure elemental real(kind=WP) function spbee_limiter(R) result(Cr)
        real(kind=WP), intent(in) :: R
        Cr = max(0._WP, min( min(2._WP*R, 0.5_WP+R/2._WP), 2._WP ))
    end function spbee_limiter
    
    !
    !
    !___________________________________________________________________________
    ! compute edge centered superbee limited tracer volume flux
    pure real(kind=WP) function spbee_compute_T12vflux(&
                                ttf0                 , &
                                ttfp1                , &
                                ttf_grad_ed          , &
                                vflux                , &
                                cfl                  , &
                                dxdy12               , & 
                                flag_2ndord_time     , &
                                flag_posdef            & 
                                ) result(T12vflux)
        real(kind=WP), intent(in) :: ttf0, ttfp1, ttf_grad_ed(4),  &
                                     vflux, cfl, dxdy12(2)
        logical      , intent(in) :: flag_2ndord_time, flag_posdef             
        real(kind=WP)             :: dttf0p1, ttfp2, ttfm1, R, Cr, vfabs, Tmean1, Tmean2
        
        !_______________________________________________________________________
        ! compute tracer difference allong edge
        dttf0p1   = ttfp1-ttf0
            
        !_______________________________________________________________________
        ! tracer Slope Ratio Calculation for upwind augmented point
        !
        !                             o
        !            ○._            .´ `.            .-○
        !            |  `._       .´  ^ dx1,dy1   .-´  |
        !            |     `._  .´    ├─>  `.  .-´     |
        !            |   □    `●━━━━━━┿━━━━━►●´----□---|-->●ttf_p2
        !            |      .-´│`.  <─┤    .´│`._      |
        !            |   .-´  dx2,dy2 v  .´  │   `._   |
        !            ○.-´      │    `. .´    │      `._○ 
        !                      │      o      │
        !                      ├------------>┤
        !                      │    dxdy12   │
        !                      │   dttf0p1   │
        !                      v             v
        !                    ttf_0      ttf_grad_ed(3:4)
        !                                   ttf_p1
        ttfp2   =   ttf0                                & 
                    + 2.0_WP * dxdy12(1)*ttf_grad_ed(3) &
                    + 2.0_WP * dxdy12(2)*ttf_grad_ed(4)  
        
        ! considering here we want to advect energy we can ensure the 
        ! variable to be positiv definit
        ttfp2   = merge(max(0.0_WP,ttfp2), ttfp2, flag_posdef)
        
        ! compute tracer slope 
        R       = (ttfp1-ttfp2)/(-dttf0p1+small)
        
        ! apply superbee limiter
        Cr      = spbee_limiter(R)     
        
        ! construct edge centered tracer value
        ! T_(i+0.5) = T_i+1 - [ (T_i+1 - T_i)/dx ]_Limited * dx/2
        !  --> this is seconds order in space, but first order in time 
        !  
        ! T_(i+0.5) = T_i+1 - [ (T_i+1 - T_i)/dx ]_Limited * (dx/2 - u_ed*n_ed*dt/2)
        ! T_(i+0.5) = T_i+1 - [ (T_i+1 - T_i)/dx ]_Limited * dx/2 *(1 - CFL_h)
        !  --> becomes second order in space and time  (Direct space-time scheme)
        !  --> u_ed  = (vel(el(1)) + vel(el(2)))/2
        !  --> n_ed  = (n_1 + n_2)/2
        !  --> CFL_h = u_ed*n_ed*dt/dx
        Tmean2  = ttfp1 + 0.5_WP * Cr * (-dttf0p1) * (1.0_WP-merge(cfl, 0.0_WP, flag_2ndord_time)) 
        
        !_______________________________________________________________________
        ! tracer Slope Ratio Calculation for downwind augmented point
        !
        !                             o
        !            ○._            .´ `.            .-○
        !            |  `._       .´  ^ dx1,dy1   .-´  |
        !            |     `._  .´    ├─>  `.  .-´     |
        !  ttf_m1●<--|---□----`●━━━━━━┿━━━━━►●´    □   |
        !            |      .-´│`.  <─┤    .´│`._      |
        !            |   .-´  dx2,dy2 v  .´  │   `._   |
        !            ○.-´      │    `. .´    │      `._○ 
        !                      │      o      │
        !                      ├------------>┤
        !                      │    dxdy12   │
        !                      │   dttf0p1   │
        !                      v             v
        !                 ttf_grad_ed(1:2) ttf_p1 
        !                    ttf0    
        ttfm1   =   ttfp1                               & 
                    - 2.0_WP * dxdy12(1)*ttf_grad_ed(1) &
                    - 2.0_WP * dxdy12(2)*ttf_grad_ed(2)  
        
        ! considering here we want to advect energy we can ensure the 
        ! variable to be positiv definit
        ttfm1   = merge(max(0.0_WP,ttfm1), ttfm1, flag_posdef)
        
        ! compute tracer slope 
        R       = (ttf0-ttfm1)/(dttf0p1+small)
        
        ! apply superbee limiter
        Cr      = spbee_limiter(R)     
        
        ! construct edge centered tracer value
        ! T_(i+0.5) = T_i + [ (T_i+1 - T_i)/dx ]_Limited * dx/2
        !  --> this is seconds order in space, but first order in time 
        !  
        ! T_(i+0.5) = T_i + [ (T_i+1 - T_i)/dx ]_Limited * (dx/2 - u_ed*n_ed*dt/2)
        ! T_(i+0.5) = T_i + [ (T_i+1 - T_i)/dx ]_Limited * dx/2 *(1 - CFL_h)
        !  --> becomes second order in space and time  (Direct space-time scheme)
        !  --> u_ed  = (vel(el(1)) + vel(el(2)))/2
        !  --> n_ed  = (n_1 + n_2)/2
        !  --> CFL_h = u_ed*n_ed*dt/dx
        Tmean1  = ttf0 + 0.5_WP * Cr * (dttf0p1) * (1.0_WP-merge(cfl, 0.0_WP, flag_2ndord_time)) 
        
        ! compute upwind/downwind related tracer volume flux 
        vfabs   = abs(vflux)
        T12vflux= 0.5_WP*( (vflux+vfabs)*Tmean1+ (vflux-vfabs)*Tmean2)
    end function spbee_compute_T12vflux
    
end subroutine adv_tra_hor_spbee
    
    
