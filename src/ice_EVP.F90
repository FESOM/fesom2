module ice_EVP_interfaces
    interface
        subroutine stress_tensor(ice, partit, mesh)
        USE MOD_ICE
        USE MOD_PARTIT
        USE MOD_PARSUP
        USE MOD_MESH
        type(t_ice)   , intent(inout), target :: ice
        type(t_partit), intent(inout), target :: partit
        type(t_mesh)  , intent(in)   , target :: mesh
        end subroutine

        subroutine stress2rhs(ice, partit, mesh)
        USE MOD_ICE
        USE MOD_PARTIT
        USE MOD_PARSUP
        USE MOD_MESH
        type(t_ice)   , intent(inout), target :: ice
        type(t_partit), intent(inout), target :: partit
        type(t_mesh)  , intent(in)   , target :: mesh
        end subroutine
    end interface
end module

module ice_EVPdynamics_interface
    interface
        subroutine EVPdynamics(ice, partit, mesh)
        USE MOD_ICE
        USE MOD_PARTIT
        USE MOD_PARSUP
        USE MOD_MESH
        type(t_ice)   , intent(inout), target :: ice
        type(t_partit), intent(inout), target :: partit
        type(t_mesh)  , intent(in)   , target :: mesh
        end subroutine
    end interface
end module

!
! Contains routines of EVP dynamics
!
!_______________________________________________________________________________
! EVP rheology. The routine computes stress tensor components based on ice
! velocity field. They are stored as elemental arrays (sigma11, sigma22 and
! sigma12). The ocean velocity is at nodal locations.
subroutine stress_tensor(ice, partit, mesh)
    USE MOD_ICE
    USE MOD_PARTIT
    USE MOD_PARSUP
    USE MOD_MESH
    use o_param
    use g_CONFIG
#if defined (__icepack)
    use icedrv_main,   only: rdg_conv_elem, rdg_shear_elem, strength
#endif
    implicit none
    type(t_partit), intent(inout), target :: partit
    type(t_ice)   , intent(inout), target :: ice
    type(t_mesh)  , intent(in)   , target :: mesh
    !___________________________________________________________________________
    integer         :: el
    real(kind=WP)   :: det1, det2, dte, vale, r1, r2, r3, si1, si2
    real(kind=WP)   :: zeta, delta, delta_inv, d1, d2
    !___________________________________________________________________________
    ! pointer on necessary derived types
    real(kind=WP), dimension(:), pointer  :: u_ice, v_ice
    real(kind=WP), dimension(:), pointer  :: eps11, eps12, eps22
    real(kind=WP), dimension(:), pointer  :: sigma11, sigma12, sigma22
    real(kind=WP), dimension(:), pointer  :: ice_strength
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"
    u_ice       => ice%uice(:)
    v_ice       => ice%vice(:)
    eps11       => ice%work%eps11(:)
    eps12       => ice%work%eps12(:)
    eps22       => ice%work%eps22(:)
    sigma11     => ice%work%sigma11(:)
    sigma12     => ice%work%sigma12(:)
    sigma22     => ice%work%sigma22(:)
    ice_strength=> ice%work%ice_strength(:)
    !___________________________________________________________________________
    vale = 1.0_WP/(ice%ellipse**2)
    dte  = ice%ice_dt/(1.0_WP*ice%evp_rheol_steps)
    det1 = 1.0_WP/(1.0_WP + 0.5_WP*ice%Tevp_inv*dte)
    det2 = 1.0_WP/(1.0_WP + 0.5_WP*ice%Tevp_inv*dte) !*ellipse**2

#ifndef ENABLE_OPENACC
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(el, r1, r2, r3, si1, si2, zeta, delta, delta_inv, d1, d2)
#else
!$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT)
#endif
    do el=1,myDim_elem2D
        !_______________________________________________________________________
        ! if element contains cavity node skip it
        if (ulevels(el) > 1) cycle

        ! ===== Check if there is ice on elem
        ! There is no ice in elem
        ! if (any(m_ice(elnodes)<= 0.) .or. any(a_ice(elnodes) <=0.)) CYCLE
        if (ice_strength(el) > 0.) then
            ! =====
            ! ===== Deformation rate tensor on element elem:
                !du/dx
            eps11(el) = sum(gradient_sca(1:3,el)*U_ice(elem2D_nodes(1:3,el))) &
                - metric_factor(el) * sum(V_ice(elem2D_nodes(1:3,el)))/3.0_WP

            eps22(el) = sum(gradient_sca(4:6, el)*V_ice(elem2D_nodes(1:3,el)))

            eps12(el) = 0.5_WP*(sum(gradient_sca(4:6,el)*U_ice(elem2D_nodes(1:3,el))) &
                        + sum(gradient_sca(1:3,el)*V_ice(elem2D_nodes(1:3,el))) &
                        + metric_factor(el) * sum(U_ice(elem2D_nodes(1:3,el)))/3.0_WP)
            ! ===== moduli:
            delta = sqrt((eps11(el)*eps11(el) + eps22(el)*eps22(el))*(1.0_WP+vale) + 4.0_WP*vale*eps12(el)*eps12(el) + &
                                2.0_WP*eps11(el)*eps22(el)*(1.0_WP-vale))

            ! =======================================
            ! ===== Here the EVP rheology piece starts
            ! =======================================

            ! ===== viscosity zeta should exceed zeta_min
            ! (done via limiting delta from above)

            !if(delta>pressure/ice%zeta_min) delta=pressure/ice%zeta_min
                !It does not work properly by
            !creating response where ice_strength is small
                ! Uncomment and test if necessary

            ! ===== if delta is too small or zero, viscosity will too large (unlimited)
            ! (limit delta_inv)
            delta_inv = 1.0_WP/max(delta,ice%delta_min)
            zeta = ice_strength(el)*delta_inv
            ! ===== Limiting pressure/Delta  (zeta): it may still happen that pressure/Delta
            ! is too large in some regions and CFL criterion is violated.
            ! The regularization below was introduced by Hunke,
            ! but seemingly is not used in the current CICE.
            ! Without it divergence and zeta can be noisy (but code
            ! remains stable), using it reduces viscosities too strongly.
            ! It is therefore commented

            !if (zeta>ice%clim_evp*voltriangle(el)) then
            !zeta=ice%clim_evp*voltriangle(el)
            !end if

            zeta = zeta*ice%Tevp_inv

            r1  = zeta*(eps11(el)+eps22(el)) - ice_strength(el)*ice%Tevp_inv
            r2  = zeta*(eps11(el)-eps22(el))*vale
            r3  = zeta*eps12(el)*vale

            si1 = det1*(sigma11(el) + sigma22(el) + dte*r1)
            si2 = det2*(sigma11(el) - sigma22(el) + dte*r2)

            sigma12(el) = det2*(sigma12(el)+dte*r3)
            sigma11(el) = 0.5_WP*(si1+si2)
            sigma22(el) = 0.5_WP*(si1-si2)

#if defined (__icepack)
            rdg_conv_elem(el)  = -min((eps11(el)+eps22(el)),0.0_WP)
            rdg_shear_elem(el) = 0.5_WP*(delta - abs(eps11(el)+eps22(el)))
#endif
        endif
    end do
#ifndef ENABLE_OPENACC
!$OMP END PARALLEL DO
#else
!$ACC END PARALLEL LOOP
#endif
end subroutine stress_tensor
!
!
!_______________________________________________________________________________
! EVP implementation:
! Computes the divergence of stress tensor and puts the result into the
! rhs vectors
subroutine stress2rhs(ice, partit, mesh)
    USE MOD_ICE
    USE MOD_PARTIT
    USE MOD_PARSUP
    USE MOD_MESH
    USE o_PARAM
    IMPLICIT NONE
    type(t_ice)   , intent(inout), target :: ice
    type(t_partit), intent(inout), target :: partit
    type(t_mesh)  , intent(in)   , target :: mesh
    !___________________________________________________________________________
    INTEGER                   :: n, el,  k
    REAL(kind=WP)             :: val3
    !___________________________________________________________________________
    ! pointer on necessary derived types
    real(kind=WP), dimension(:), pointer  :: sigma11, sigma12, sigma22
    real(kind=WP), dimension(:), pointer  :: u_rhs_ice, v_rhs_ice, rhs_a, rhs_m
    real(kind=WP), dimension(:), pointer  :: inv_areamass, ice_strength
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"
    sigma11      => ice%work%sigma11(:)
    sigma12      => ice%work%sigma12(:)
    sigma22      => ice%work%sigma22(:)
    u_rhs_ice    => ice%uice_rhs(:)
    v_rhs_ice    => ice%vice_rhs(:)
    rhs_a        => ice%data(1)%values_rhs(:)
    rhs_m        => ice%data(2)%values_rhs(:)
    inv_areamass => ice%work%inv_areamass(:)
    ice_strength => ice%work%ice_strength(:)

    !___________________________________________________________________________
    val3=1/3.0_WP

#ifndef ENABLE_OPENACC
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(n, el, k)
!$OMP DO
#else
    !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT)
#endif
    DO  n=1, myDim_nod2D
        U_rhs_ice(n)=0.0_WP
        V_rhs_ice(n)=0.0_WP
    END DO

#ifndef ENABLE_OPENACC
!$OMP END DO
#else
    !$ACC END PARALLEL LOOP
#endif

#ifndef ENABLE_OPENACC
!$OMP DO
#else
    !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT)
    #if !defined(DISABLE_OPENACC_ATOMICS)
       !$ACC ATOMIC UPDAATE
    #else
    !$ACC UPDATE SELF(u_rhs_ice, v_rhs_ice, sigma11, sigma12, sigma22)
    #endif
#endif
    do el=1,myDim_elem2D
        ! ===== Skip if ice is absent
        !   if (any(m_ice(elnodes)<= 0.) .or. any(a_ice(elnodes) <=0.)) CYCLE
        !_______________________________________________________________________
        ! if element contains cavity node skip it
        if (ulevels(el) > 1) cycle

        !_______________________________________________________________________
        if (ice_strength(el) > 0._WP) then
            DO k=1,3

#ifndef ENABLE_OPENACC
#if defined(_OPENMP) && !defined(__openmp_reproducible)
        call omp_set_lock  (partit%plock(elem2D_nodes(k,el)))
#else
!$OMP ORDERED
#endif
#endif
#ifdef ENABLE_OPENACC
#if !defined(DISABLE_OPENACC_ATOMICS)
                !$ACC ATOMIC UPDATE
#endif
#endif
                U_rhs_ice(elem2D_nodes(k,el)) = U_rhs_ice(elem2D_nodes(k,el)) &
                - elem_area(el) * &
                    (sigma11(el)*gradient_sca(k,el) + sigma12(el)*gradient_sca(k+3,el) &
                    +sigma12(el)*val3*metric_factor(el))            !metrics

#ifdef ENABLE_OPENACC
#if !defined(DISABLE_OPENACC_ATOMICS)
                !$ACC ATOMIC UPDATE
#endif
#endif
                V_rhs_ice(elem2D_nodes(k,el)) = V_rhs_ice(elem2D_nodes(k,el)) &
                    - elem_area(el) * &
                    (sigma12(el)*gradient_sca(k,el) + sigma22(el)*gradient_sca(k+3,el) &
                    -sigma11(el)*val3*metric_factor(el))

#ifndef ENABLE_OPENACC
#if defined(_OPENMP) && !defined(__openmp_reproducible)
        call omp_unset_lock(partit%plock(elem2D_nodes(k,el)))
#else
!$OMP END ORDERED
#endif
#endif
            END DO
        endif
    end do
#ifdef ENABLE_OPENACC
   #if !defined(DISABLE_OPENACC_ATOMICS)
    !$ACC UPDATE DEVICE(u_rhs_ice, v_rhs_ice)
    #endif
#endif

#ifndef ENABLE_OPENACC
!$OMP END DO
#else
    !$ACC END PARALLEL LOOP
#endif

#ifndef ENABLE_OPENACC
!$OMP DO
#else
    !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT)
#endif
    DO n=1, myDim_nod2D
        !_______________________________________________________________________
        ! if cavity node skip it
        if (ulevels_nod2d(n)>1) cycle

        !_______________________________________________________________________
        if (inv_areamass(n) > 0._WP) then
            U_rhs_ice(n) = U_rhs_ice(n)*inv_areamass(n) + rhs_a(n)
            V_rhs_ice(n) = V_rhs_ice(n)*inv_areamass(n) + rhs_m(n)
        else
            U_rhs_ice(n) = 0._WP
            V_rhs_ice(n) = 0._WP
        endif
    END DO
#ifndef ENABLE_OPENACC
!$OMP END DO
!$OMP END PARALLEL
#else
    !$ACC END PARALLEL LOOP
#endif
end subroutine stress2rhs
!
!
!_______________________________________________________________________________
! EVP implementation. Does subcycling and boundary conditions.
! Velocities at nodes
subroutine EVPdynamics(ice, partit, mesh)
    USE MOD_ICE
    USE MOD_PARTIT
    USE MOD_PARSUP
    USE MOD_MESH
    USE o_PARAM
    USE o_ARRAYS
    USE g_CONFIG
    USE g_comm_auto
    use ice_EVP_interfaces
#if defined (__icepack)
    use icedrv_main,   only: rdg_conv_elem, rdg_shear_elem, strength
    use icedrv_main,   only: icepack_to_fesom
#endif
    IMPLICIT NONE
    type(t_ice)   , intent(inout), target :: ice
    type(t_partit), intent(inout), target :: partit
    type(t_mesh)  , intent(in)   , target :: mesh
    !___________________________________________________________________________
    integer         :: steps, shortstep
    real(kind=WP)   :: rdt, asum, msum, r_a, r_b
    real(kind=WP)   :: drag, det, umod, rhsu, rhsv
    integer         :: n, ed, ednodes(2), el,  elnodes(3)
    real(kind=WP)   :: ax, ay, aa, elevation_dx, elevation_dy

    real(kind=WP)   :: elevation_elem(3), p_ice(3)
    integer         :: use_pice

    real(kind=WP)   :: eta, delta
    integer         :: k
    real(kind=WP)   :: vale, dx(3), dy(3), val3
    real(kind=WP)   :: det1, det2, r1, r2, r3, si1, si2, dte
    real(kind=WP)   :: zeta, delta_inv, d1, d2
    INTEGER         :: elem
    !_______________________________________________________________________________
    ! pointer on necessary derived types
    real(kind=WP), dimension(:), pointer  :: u_ice, v_ice
    real(kind=WP), dimension(:), pointer  :: a_ice, m_ice, m_snow
    real(kind=WP), dimension(:), pointer  :: u_ice_old, v_ice_old
    real(kind=WP), dimension(:), pointer  :: u_rhs_ice, v_rhs_ice, rhs_a, rhs_m
    real(kind=WP), dimension(:), pointer  :: u_w, v_w, elevation
    real(kind=WP), dimension(:), pointer  :: stress_atmice_x, stress_atmice_y
    real(kind=WP), dimension(:), pointer  :: inv_areamass, inv_mass, ice_strength
#if defined (__icepack)
    real(kind=WP), dimension(:), pointer  :: a_ice_old, m_ice_old, m_snow_old
#endif
    real(kind=WP)              , pointer  :: inv_rhowat, rhosno, rhoice
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"
    u_ice           => ice%uice(:)
    v_ice           => ice%vice(:)
    a_ice           => ice%data(1)%values(:)
    m_ice           => ice%data(2)%values(:)
    m_snow          => ice%data(3)%values(:)
    u_ice_old       => ice%uice_old(:)
    v_ice_old       => ice%vice_old(:)
    u_rhs_ice       => ice%uice_rhs(:)
    v_rhs_ice       => ice%vice_rhs(:)
    rhs_a           => ice%data(1)%values_rhs(:)
    rhs_m           => ice%data(2)%values_rhs(:)
    u_w             => ice%srfoce_u(:)
    v_w             => ice%srfoce_v(:)
    elevation       => ice%srfoce_ssh(:)
    stress_atmice_x => ice%stress_atmice_x(:)
    stress_atmice_y => ice%stress_atmice_y(:)
#if defined (__icepack)
    a_ice_old       => ice%data(1)%values_old(:)
    m_ice_old       => ice%data(2)%values_old(:)
    m_snow_old      => ice%data(3)%values_old(:)
#endif
    rhosno          => ice%thermo%rhosno
    rhoice          => ice%thermo%rhoice
    inv_rhowat      => ice%thermo%inv_rhowat

    inv_areamass    => ice%work%inv_areamass(:)
    inv_mass        => ice%work%inv_mass(:)
    ice_strength    => ice%work%ice_strength(:)

    !___________________________________________________________________________
    ! If Icepack is used, always update the tracers
#if defined (__icepack)
    a_ice_old(:)  = a_ice(:)
    m_ice_old(:)  = a_ice(:)
    m_snow_old(:) = m_snow(:)
    call icepack_to_fesom (nx_in=(myDim_nod2D+eDim_nod2D), &
                            aice_out=a_ice,                 &
                            vice_out=m_ice,                 &
                            vsno_out=m_snow)
#endif

    !___________________________________________________________________________
    rdt=ice%ice_dt/(1.0*ice%evp_rheol_steps)
    ax=cos(ice%theta_io)
    ay=sin(ice%theta_io)

    !___________________________________________________________________________
    ! Precompute values that are never changed during the iteration
#ifndef ENABLE_OPENACC
!$OMP PARALLEL DO
#else
    !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT)
#endif
    do n=1, myDim_nod2D+eDim_nod2D
       inv_areamass(n) =0.0_WP
       inv_mass(n)     =0.0_WP
       rhs_a(n)        =0.0_WP
       rhs_m(n)        =0.0_WP
    end do
#ifndef ENABLE_OPENACC
!$OMP END PARALLEL DO
#else
    !$ACC END PARALLEL LOOP
#endif

#ifndef ENABLE_OPENACC
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(n)
#else
    !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT)
#endif
    do n=1,myDim_nod2D
        !_______________________________________________________________________
        ! if cavity node skip it
        if (ulevels_nod2d(n)>1) cycle

        !_______________________________________________________________________
        if ((rhoice*m_ice(n)+rhosno*m_snow(n)) > 1.e-3_WP) then
            inv_areamass(n) = 1._WP/(area(1,n)*(rhoice*m_ice(n)+rhosno*m_snow(n)))
        else
            inv_areamass(n) = 0._WP
        endif

        if (a_ice(n) < 0.01_WP) then
            ! Skip if ice is absent
            inv_mass(n) = 0._WP
        else
            inv_mass(n) = (rhoice*m_ice(n)+rhosno*m_snow(n))/a_ice(n)
            inv_mass(n) = 1.0_WP/max(inv_mass(n), 9.0_WP)        ! Limit the mass
                                        ! if it is too small
        endif
        rhs_a(n)=0.0_WP       ! these are used as temporal storage here
        rhs_m(n)=0.0_WP       ! for the contribution due to ssh
    enddo
#ifndef ENABLE_OPENACC
!$OMP END PARALLEL DO
#else
    !$ACC END PARALLEL LOOP
#endif
    !___________________________________________________________________________
    use_pice=0
    if (use_floatice .and.  .not. trim(which_ale)=='linfs') use_pice=1
    if ( .not. trim(which_ALE)=='linfs') then
        ! for full free surface include pressure from ice mass

#ifndef ENABLE_OPENACC
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(el, elnodes, msum, asum, aa, p_ice, elevation_elem, elevation_dx, elevation_dy)
!$OMP DO
#else

#if !defined(DISABLE_OPENACC_ATOMICS)
        !$ACC PARALLEL LOOP GANG VECTOR PRIVATE(elnodes, elevation_elem, p_ice) DEFAULT(PRESENT)
#else
        !$ACC UPDATE SELF(rhs_a, rhs_m, m_ice, a_ice)
#endif
#endif
        do el = 1,myDim_elem2D
            elnodes = elem2D_nodes(:,el)
            ice_strength(el)=0.0_WP
            !___________________________________________________________________
            ! if element has any cavity node skip it
            if (ulevels(el) > 1) cycle

            !___________________________________________________________________
            if (any(m_ice(elnodes)<=0._WP) .or. &
                any(a_ice(elnodes)<=0._WP)) then

                ! There is no ice in elem
                ice_strength(el) = 0._WP

            !___________________________________________________________________
            else
                msum = sum(m_ice(elnodes))/3.0_WP
                asum = sum(a_ice(elnodes))/3.0_WP

                !_______________________________________________________________
                ! Hunke and Dukowicz c*h*p*
#if defined (__icepack)
                ice_strength(el) = ice%pstar*msum*exp(-ice%c_pressure*(1.0_WP-asum))
#else
                ice_strength(el) = ice%pstar*msum*exp(-ice%c_pressure*(1.0_WP-asum))
#endif
                ice_strength(el) = 0.5_WP*ice_strength(el)

                !_______________________________________________________________
                ! use rhs_m and rhs_a for storing the contribution from elevation:
                aa = 9.81_WP*elem_area(el)/3.0_WP

                !_______________________________________________________________
                ! add and limit pressure from ice weight in case of floating ice
                ! like in FESOM 1.4
                p_ice=(rhoice*m_ice(elnodes)+rhosno*m_snow(elnodes))*inv_rhowat
                do n=1,3
                    p_ice(n)=min(p_ice(n),max_ice_loading)
                end do
                !!PS p_ice= 0.0_WP

                !_______________________________________________________________
                elevation_elem = elevation(elnodes)
                elevation_dx   = sum(gradient_sca(1:3,el)*(elevation_elem+p_ice*use_pice))
                elevation_dy   = sum(gradient_sca(4:6,el)*(elevation_elem+p_ice*use_pice))

                !_______________________________________________________________
                do k = 1, 3
#if !defined(DISABLE_OPENACC_ATOMICS)
                    !$ACC ATOMIC UPDATE
#endif
                    rhs_a(elnodes(k)) = rhs_a(elnodes(k))-aa*elevation_dx
#if !defined(DISABLE_OPENACC_ATOMICS)
                    !$ACC ATOMIC UPDATE
#endif
                    rhs_m(elnodes(k)) = rhs_m(elnodes(k))-aa*elevation_dy
                end do
            end if
        enddo
#ifndef ENABLE_OPENACC
!$OMP END DO
!$OMP END PARALLEL
#else
#if !defined(DISABLE_OPENACC_ATOMICS)
        !$ACC END PARALLEL LOOP
#else
        !$ACC UPDATE DEVICE(rhs_a, rhs_m, ice_strength)
#endif
#endif
    else
        ! for linear free surface
#ifndef ENABLE_OPENACC
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(el, elnodes, msum, asum, aa, elevation_elem, elevation_dx, elevation_dy)
#else
#if !defined(DISABLE_OPENACC_ATOMICS)
            !$ACC PARALLEL LOOP GANG VECTOR PRIVATE(elnodes) DEFAULT(PRESENT)
#else
            !$ACC UPDATE SELF(rhs_a, rhs_m, ice_strength, m_ice, a_ice)
#endif
#endif
        do el = 1,myDim_elem2D
            ice_strength(el)=0.0_WP
            elnodes = elem2D_nodes(:,el)
            !___________________________________________________________________
            ! if element has any cavity node skip it
            if (ulevels(el) > 1) cycle

            !___________________________________________________________________
            if (any(m_ice(elnodes) <= 0._WP) .or. &
                any(a_ice(elnodes) <=0._WP)) then

                ! There is no ice in elem
                ice_strength(el) = 0._WP
            else
                msum = sum(m_ice(elnodes))/3.0_WP
                asum = sum(a_ice(elnodes))/3.0_WP

                ! ===== Hunke and Dukowicz c*h*p*
#if defined (__icepack)
                ice_strength(el) = ice%pstar*msum*exp(-ice%c_pressure*(1.0_WP-asum))
#else
                ice_strength(el) = ice%pstar*msum*exp(-ice%c_pressure*(1.0_WP-asum))
#endif
                ice_strength(el) = 0.5_WP*ice_strength(el)

                ! use rhs_m and rhs_a for storing the contribution from elevation:
                aa = 9.81_WP*elem_area(el)/3.0_WP

                elevation_dx = sum(gradient_sca(1:3,el)*elevation(elnodes))
                elevation_dy = sum(gradient_sca(4:6,el)*elevation(elnodes))

                do k = 1, 3
#if !defined(DISABLE_OPENACC_ATOMICS)
                    !$ACC ATOMIC UPDATE
#endif
                    rhs_a(elnodes(k)) = rhs_a(elnodes(k))-aa*elevation_dx
#if !defined(DISABLE_OPENACC_ATOMICS)
                    !$ACC ATOMIC UPDATE
#endif
                    rhs_m(elnodes(k)) = rhs_m(elnodes(k))-aa*elevation_dy
                end do
            end if
        enddo
#ifndef ENABLE_OPENACC
!$OMP END PARALLEL DO
#else
#if !defined(DISABLE_OPENACC_ATOMICS)
            !$ACC END PARALLEL LOOP
#else
            !$ACC UPDATE DEVICE(rhs_a, rhs_m, ice_strength)
#endif
#endif
    endif ! --> if ( .not. trim(which_ALE)=='linfs') then
#ifndef ENABLE_OPENACC
!$OMP PARALLEL DO
#else

    !___________________________________________________________________________
    !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT)
#endif
    do n=1,myDim_nod2D
        if (ulevels_nod2d(n)>1) cycle
        rhs_a(n) = rhs_a(n)/area(1,n)
        rhs_m(n) = rhs_m(n)/area(1,n)
    enddo
#ifndef ENABLE_OPENACC
!$OMP END PARALLEL DO
#else
    !$ACC END PARALLEL LOOP
#endif
    !___________________________________________________________________________
    ! End of Precomputing --> And the ice stepping starts
#if defined (__icepack)
    rdg_conv_elem(:)  = 0.0_WP
    rdg_shear_elem(:) = 0.0_WP
#endif
    do shortstep=1, ice%evp_rheol_steps
        !_______________________________________________________________________
        !TODO: temporary workaround for cray16.0.1.1 bug
#if defined(_CRAYFTN)
	!dir$ noinline
#endif
        call stress_tensor(ice, partit, mesh)
#if defined(_CRAYFTN)
	!dir$ noinline
#endif
        call stress2rhs(ice, partit, mesh)

        !_______________________________________________________________________
#ifndef ENABLE_OPENACC
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(n, ed, umod, drag, rhsu, rhsv, r_a, r_b, det)
!$OMP DO
#else
        !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT)
#endif
        do n=1,myDim_nod2D+eDim_nod2D
           U_ice_old(n) = U_ice(n) !PS
           V_ice_old(n) = V_ice(n) !PS
        end do
#ifndef ENABLE_OPENACC
!$OMP END DO
#else
        !$ACC END PARALLEL LOOP
#endif
#ifndef ENABLE_OPENACC
!$OMP DO
#else
        !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT)
#endif
        do n=1,myDim_nod2D
            !___________________________________________________________________
            ! if cavity node skip it
            if ( ulevels_nod2d(n)>1 ) cycle

            !___________________________________________________________________
            if (a_ice(n) >= 0.01_WP) then               ! Skip if ice is absent
                umod = sqrt((U_ice(n)-U_w(n))**2+(V_ice(n)-V_w(n))**2)
                drag = ice%cd_oce_ice*umod*density_0*inv_mass(n)

                rhsu = U_ice(n) +rdt*(drag*(ax*U_w(n) - ay*V_w(n))+ &
                        inv_mass(n)*stress_atmice_x(n) + U_rhs_ice(n))
                rhsv = V_ice(n) +rdt*(drag*(ax*V_w(n) + ay*U_w(n))+ &
                        inv_mass(n)*stress_atmice_y(n) + V_rhs_ice(n))

                r_a      = 1._WP + ax*drag*rdt
                r_b      = rdt*(mesh%coriolis_node(n) + ay*drag)
                det      = 1.0_WP/(r_a*r_a + r_b*r_b)
                U_ice(n) = det*(r_a*rhsu +r_b*rhsv)
                V_ice(n) = det*(r_a*rhsv -r_b*rhsu)
            else  ! Set velocities to 0 if ice is absent
                U_ice(n) = 0.0_WP
                V_ice(n) = 0.0_WP
            end if
        end do
#ifndef ENABLE_OPENACC
!$OMP END DO
#else
        !$ACC END PARALLEL LOOP
#endif
        !_______________________________________________________________________
        ! apply sea ice velocity boundary condition

#ifndef ENABLE_OPENACC
!$OMP DO
#else
        ! With the binary data of np2 goes only inside the first if
        !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT)
#endif
        DO  ed=1,myDim_edge2D
            !___________________________________________________________________
            ! apply coastal sea ice velocity boundary conditions
            if(myList_edge2D(ed) > edge2D_in) then
                U_ice(edges(1:2,ed))=0.0_WP
                V_ice(edges(1:2,ed))=0.0_WP
            endif

            !___________________________________________________________________
            ! apply sea ice velocity boundary conditions at cavity-ocean edge
            if (use_cavity) then
!                 if ( (ulevels(edge_tri(1,ed))>1) .or. &
!                     ( edge_tri(2,ed)>0 .and. ulevels(edge_tri(2,ed))>1) ) then
! #if defined(_OPENMP)  && !defined(__openmp_reproducible)
!                     call omp_set_lock  (partit%plock(edges(1,ed)))
! #else
! !$OMP ORDERED
! #endif
!                     U_ice(edges(1,ed))=0.0_WP
!                     V_ice(edges(1,ed))=0.0_WP
!
! #if defined(_OPENMP) && !defined(__openmp_reproducible)
!                     call omp_unset_lock(partit%plock(edges(1,ed)))
!                     call omp_set_lock  (partit%plock(edges(2,ed)))
! #endif
!                     U_ice(edges(2,ed))=0.0_WP
!                     V_ice(edges(2,ed))=0.0_WP
!
! #if defined(_OPENMP)  && !defined(__openmp_reproducible)
!                     call omp_unset_lock(partit%plock(edges(2,ed)))
! #else
! !$OMP END ORDERED
! #endif
!                 end if
                if (ulevels(edge_tri(1,ed))>1) then
#ifndef ENABLE_OPENACC
#if defined(_OPENMP)  && !defined(__openmp_reproducible)
                    call omp_set_lock  (partit%plock(edges(1,ed)))
#else
!$OMP ORDERED
#endif
#endif
                    U_ice(edges(1,ed))=0.0_WP
                    V_ice(edges(1,ed))=0.0_WP

#ifndef ENABLE_OPENACC
#if defined(_OPENMP) && !defined(__openmp_reproducible)
                    call omp_unset_lock(partit%plock(edges(1,ed)))
                    call omp_set_lock  (partit%plock(edges(2,ed)))
#endif
#endif
                    U_ice(edges(2,ed))=0.0_WP
                    V_ice(edges(2,ed))=0.0_WP

#ifndef ENABLE_OPENACC
#if defined(_OPENMP)  && !defined(__openmp_reproducible)
                    call omp_unset_lock(partit%plock(edges(2,ed)))
#else
!$OMP END ORDERED
#endif
#endif 
                elseif ( edge_tri(2,ed)>0) then
                    if (ulevels(edge_tri(2,ed))>1) then
#ifndef ENABLE_OPENACC
#if defined(_OPENMP)  && !defined(__openmp_reproducible)
                    call omp_set_lock  (partit%plock(edges(1,ed)))
#else
!$OMP ORDERED
#endif
#endif
                    U_ice(edges(1,ed))=0.0_WP
                    V_ice(edges(1,ed))=0.0_WP

#ifndef ENABLE_OPENACC
#if defined(_OPENMP) && !defined(__openmp_reproducible)
                    call omp_unset_lock(partit%plock(edges(1,ed)))
                    call omp_set_lock  (partit%plock(edges(2,ed)))
#endif
#endif
                    U_ice(edges(2,ed))=0.0_WP
                    V_ice(edges(2,ed))=0.0_WP

#ifndef ENABLE_OPENACC
#if defined(_OPENMP)  && !defined(__openmp_reproducible)
                    call omp_unset_lock(partit%plock(edges(2,ed)))
#else
!$OMP END ORDERED
#endif
#endif

                    end if
                end if
            end if
        end do
#ifndef ENABLE_OPENACC
!$OMP END DO
!$OMP END PARALLEL
#else
        !$ACC END PARALLEL LOOP
#endif

!write(*,*) partit%mype, shortstep, 'CP4'
        !_______________________________________________________________________
        call exchange_nod(U_ice,V_ice,partit, luse_g2g = .true.)

!ifndef ENABLE_OPENACC
!$OMP BARRIER
!endif
    END DO !--> do shortstep=1, ice%evp_rheol_steps

end subroutine EVPdynamics