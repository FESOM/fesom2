! Reference kernels: verbatim copies of FESOM's adv_tra_hor_upw1 and adv_tra_hor_mfct
! (src/oce_adv_tra_hor.F90 at branch dwarf-tradv), and the CPU loops of the low-order tendency gather
! (do_oce_adv_tra) and the horizontal part of oce_tra_adv_flux2dtracer (src/oce_adv_tra_driver.F90).
module kernels_ref
    use o_PARAM
    use MOD_MESH
    use MOD_PARTIT
    implicit none
contains
subroutine adv_tra_hor_upw1(vel, ttf, partit, mesh, flux, o_init_zero)
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
       do edge=1, partit%myDim_edge2D
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
    do edge=1, partit%myDim_edge2D
        ! local indice of nodes that span up edge ed
        enodes=mesh%edges(:,edge)

        ! local index of element that contribute to edge
        el=mesh%edge_tri(:,edge)

        ! number of layers -1 at elem el(1)
        nl1=mesh%nlevels(el(1))-1

        ! index off surface layer in case of cavity !=1
        nu1=mesh%ulevels(el(1))

        ! edge_cross_dxdy(1:2,ed)... dx,dy distance from element centroid el(1) to
        ! center of edge --> needed to calc flux perpedicular to edge from elem el(1)
        deltaX1=mesh%edge_cross_dxdy(1,edge)
        deltaY1=mesh%edge_cross_dxdy(2,edge)
        a=r_earth*mesh%elem_cos(el(1))

        !_______________________________________________________________________
        ! same parameter but for other element el(2) that contributes to edge ed
        ! if el(2)==0 than edge is boundary edge
        nl2=0
        nu2=0
        if(el(2)>0) then
            deltaX2=mesh%edge_cross_dxdy(3,edge)
            deltaY2=mesh%edge_cross_dxdy(4,edge)
            ! number of layers -1 at elem el(2)
            nl2=mesh%nlevels(el(2))-1
            nu2=mesh%ulevels(el(2))
            a=0.5_WP*(a+r_earth*mesh%elem_cos(el(2)))
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
           vflux=(-VEL(2,nz,el(1))*deltaX1 + VEL(1,nz,el(1))*deltaY1)*mesh%helem(nz,el(1))

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
                vflux=(VEL(2,nz,el(2))*deltaX2 - VEL(1,nz,el(2))*deltaY2)*mesh%helem(nz,el(2))

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
            vflux=(-VEL(2,nz,el(1))*deltaX1 + VEL(1,nz,el(1))*deltaY1)*mesh%helem(nz,el(1)) &
                  +(VEL(2,nz,el(2))*deltaX2 - VEL(1,nz,el(2))*deltaY2)*mesh%helem(nz,el(2))

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
           vflux=(-VEL(2,nz,el(1))*deltaX1 + VEL(1,nz,el(1))*deltaY1)*mesh%helem(nz,el(1))
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
                vflux=(VEL(2,nz,el(2))*deltaX2 - VEL(1,nz,el(2))*deltaY2)*mesh%helem(nz,el(2))
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

    subroutine adv_tra_hor_mfct(vel, ttf, partit, mesh, num_ord, flux, edge_up_dn_grad,                 o_init_zero)
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
       do edge=1, partit%myDim_edge2D
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
    do edge=1, partit%myDim_edge2D
        ! local indice of nodes that span up edge ed
        enodes=mesh%edges(:,edge)

        ! local index of element that contribute to edge
        el=mesh%edge_tri(:,edge)

        ! number of layers -1 at elem el(1)
        nl1=mesh%nlevels(el(1))-1

        ! index off surface layer in case of cavity !=1
        nu1=mesh%ulevels(el(1))

        ! edge_cross_dxdy(1:2,ed)... dx,dy distance from element centroid el(1) to
        ! center of edge --> needed to calc flux perpedicular to edge from elem el(1)
        deltaX1=mesh%edge_cross_dxdy(1,edge)
        deltaY1=mesh%edge_cross_dxdy(2,edge)
        a=r_earth*mesh%elem_cos(el(1))

        !_______________________________________________________________________
        ! same parameter but for other element el(2) that contributes to edge ed
        ! if el(2)==0 than edge is boundary edge
        nl2=0
        nu2=0
        if(el(2)>0) then
            deltaX2=mesh%edge_cross_dxdy(3,edge)
            deltaY2=mesh%edge_cross_dxdy(4,edge)
            ! number of layers -1 at elem el(2)
            nl2=mesh%nlevels(el(2))-1
            nu2=mesh%ulevels(el(2))
            a=0.5_WP*(a+r_earth*mesh%elem_cos(el(2)))
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
                  mesh%edge_dxdy(1,edge)*a*edge_up_dn_grad(2,nz,edge)+ &
                  mesh%edge_dxdy(2,edge)*r_earth*edge_up_dn_grad(4,nz,edge))/6.0_WP

           Tmean1=ttf(nz, enodes(1))+ &
                  (2.0_WP*(ttf(nz, enodes(2))-ttf(nz,enodes(1)))+ &
                  mesh%edge_dxdy(1,edge)*a*edge_up_dn_grad(1,nz,edge)+ &
                  mesh%edge_dxdy(2,edge)*r_earth*edge_up_dn_grad(3,nz,edge))/6.0_WP

           !____________________________________________________________________
           ! volume flux across the segments
           vflux=(-VEL(2,nz,el(1))*deltaX1 + VEL(1,nz,el(1))*deltaY1)*mesh%helem(nz,el(1))
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
                    mesh%edge_dxdy(1,edge)*a*edge_up_dn_grad(2,nz,edge)+ &
                    mesh%edge_dxdy(2,edge)*r_earth*edge_up_dn_grad(4,nz,edge))/6.0_WP

            Tmean1=ttf(nz, enodes(1))+ &
                    (2.0_WP*(ttf(nz, enodes(2))-ttf(nz,enodes(1)))+ &
                    mesh%edge_dxdy(1,edge)*a*edge_up_dn_grad(1,nz,edge)+ &
                    mesh%edge_dxdy(2,edge)*r_earth*edge_up_dn_grad(3,nz,edge))/6.0_WP
            !___________________________________________________________________
            ! volume flux across the segments
            vflux=(VEL(2,nz,el(2))*deltaX2 - VEL(1,nz,el(2))*deltaY2)*mesh%helem(nz,el(2))
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
                  mesh%edge_dxdy(1,edge)*a*edge_up_dn_grad(2,nz,edge)+ &
                  mesh%edge_dxdy(2,edge)*r_earth*edge_up_dn_grad(4,nz,edge))/6.0_WP

            ! use upwind triangle to interpolate Tracer to edge center with
            ! fancy scheme --> Linear upwind reconstruction
            ! T_n+0.5 = T_n + 1/2*deltax*GRADIENT
            ! --> GRADIENT = 2/3 GRAD_edgecenter + 1/3 GRAD_downwindtri
            ! T_n+0.5 = T_n + 2/6*(T_n+1-T_n) + 1/6*gradT_down
            ! --> Tmean1 ... edge center interpolated Tracer using tracer
            !     gradient info from downwind triangle
            Tmean1=ttf(nz, enodes(1))+ &
                   (2.0_WP*(ttf(nz, enodes(2))-ttf(nz,enodes(1)))+ &
                   mesh%edge_dxdy(1,edge)*a*edge_up_dn_grad(1,nz,edge)+ &
                   mesh%edge_dxdy(2,edge)*r_earth*edge_up_dn_grad(3,nz,edge))/6.0_WP
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
            vflux=(-VEL(2,nz,el(1))*deltaX1 + VEL(1,nz,el(1))*deltaY1)*mesh%helem(nz,el(1)) &
                  +(VEL(2,nz,el(2))*deltaX2 - VEL(1,nz,el(2))*deltaY2)*mesh%helem(nz,el(2))

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
                  mesh%edge_dxdy(1,edge)*a*edge_up_dn_grad(2,nz,edge)+ &
                  mesh%edge_dxdy(2,edge)*r_earth*edge_up_dn_grad(4,nz,edge))/6.0_WP

           Tmean1=ttf(nz, enodes(1))+ &
                  (2.0_WP*(ttf(nz, enodes(2))-ttf(nz,enodes(1)))+ &
                  mesh%edge_dxdy(1,edge)*a*edge_up_dn_grad(1,nz,edge)+ &
                  mesh%edge_dxdy(2,edge)*r_earth*edge_up_dn_grad(3,nz,edge))/6.0_WP

           !____________________________________________________________________
           ! volume flux across the segments
           vflux=(-VEL(2,nz,el(1))*deltaX1 + VEL(1,nz,el(1))*deltaY1)*mesh%helem(nz,el(1))
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
                  mesh%edge_dxdy(1,edge)*a*edge_up_dn_grad(2,nz,edge)+ &
                  mesh%edge_dxdy(2,edge)*r_earth*edge_up_dn_grad(4,nz,edge))/6.0_WP

           Tmean1=ttf(nz, enodes(1))+ &
                  (2.0_WP*(ttf(nz, enodes(2))-ttf(nz,enodes(1)))+ &
                  mesh%edge_dxdy(1,edge)*a*edge_up_dn_grad(1,nz,edge)+ &
                  mesh%edge_dxdy(2,edge)*r_earth*edge_up_dn_grad(3,nz,edge))/6.0_WP

           !____________________________________________________________________
           ! volume flux across the segments
           vflux=(VEL(2,nz,el(2))*deltaX2 - VEL(1,nz,el(2))*deltaY2)*mesh%helem(nz,el(2))
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

! low-order tendency: each node gathers the fluxes of its incident edges (do_oce_adv_tra, CPU path)
subroutine lo_gather(adv_flux_hor, fct_LO_tend, partit, mesh)
    type(t_partit), intent(in), target :: partit
    type(t_mesh),   intent(in), target :: mesh
    real(kind=WP),      intent(in)    :: adv_flux_hor(mesh%nl-1, partit%myDim_edge2D)
    real(kind=WP_full), intent(inout) :: fct_LO_tend(mesh%nl-1, partit%myDim_nod2D+partit%eDim_nod2D)
    integer :: n, nz, k, edge, el(2), nl1, nu1, nl2, nu2, nl12, nu12
!$OMP PARALLEL DO
        do n=1, partit%myDim_nod2D+partit%eDim_nod2D
           do nz=1, mesh%nl - 1
              fct_LO_tend(nz,n) = 0.0_WP_full
           end do
        end do
!$OMP END PARALLEL DO
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(n, k, edge, el, nl1, nu1, nl2, nu2, nu12, nl12, nz)
        do n=1, partit%myDim_nod2D+partit%eDim_nod2D
            do k=1, mesh%nod_in_edge2D_num(n)
                edge=mesh%nod_in_edge2D(k,n)
                el=mesh%edge_tri(:,edge)
                nl1=mesh%nlevels(el(1))-1
                nu1=mesh%ulevels(el(1))
                nl2=0
                nu2=0
                if(el(2)>0) then
                    nl2=mesh%nlevels(el(2))-1
                    nu2=mesh%ulevels(el(2))
                end if
                nl12 = max(nl1,nl2)
                nu12 = nu1
                if (nu2>0) nu12 = min(nu1,nu2)
                do nz=nu12, nl12
                    fct_LO_tend(nz, n)=fct_LO_tend(nz, n)+mesh%nod_in_edge2D_sgn(k,n)*real(adv_flux_hor(nz, edge), WP_full)
                end do
            end do
        end do
!$OMP END PARALLEL DO
end subroutine lo_gather

! horizontal part of oce_tra_adv_flux2dtracer (CPU path)
subroutine hor_gather(dt, dttf_h, flux_h, partit, mesh)
    real(kind=WP), intent(in) :: dt
    type(t_partit), intent(in), target :: partit
    type(t_mesh),   intent(in), target :: mesh
    real(kind=WP_full), intent(inout) :: dttf_h(mesh%nl-1, partit%myDim_nod2D+partit%eDim_nod2D)
    real(kind=WP),      intent(in)    :: flux_h(mesh%nl-1, partit%myDim_edge2D)
    integer :: n, nz, k, edge, el(2), nl1, nu1, nl2, nu2, nl12, nu12
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(n, k, edge, el, nl1, nu1, nl2, nu2, nu12, nl12, nz)
    do n=1, partit%myDim_nod2D+partit%eDim_nod2D
        do k=1, mesh%nod_in_edge2D_num(n)
            edge=mesh%nod_in_edge2D(k,n)
            el=mesh%edge_tri(:,edge)
            nl1=mesh%nlevels(el(1))-1
            nu1=mesh%ulevels(el(1))
            nl2=0
            nu2=0
            if(el(2)>0) then
                nl2=mesh%nlevels(el(2))-1
                nu2=mesh%ulevels(el(2))
            end if
            nl12 = max(nl1,nl2)
            nu12 = nu1
            if (nu2>0) nu12 = min(nu1,nu2)
            do nz=nu12, nl12
                dttf_h(nz,n)=dttf_h(nz,n)+mesh%nod_in_edge2D_sgn(k,n)*real(flux_h(nz,edge), WP_full)*dt/mesh%areasvol(nz,n)
            end do
        end do
    end do
!$OMP END PARALLEL DO
end subroutine hor_gather
end module kernels_ref
