module ice_fct_interfaces
    interface
        subroutine ice_mass_matrix_fill(ice, partit, mesh)
        USE MOD_ICE
        USE MOD_PARTIT
        USE MOD_PARSUP
        USE MOD_MESH
        type(t_ice),    intent(inout), target :: ice
        type(t_partit), intent(inout), target :: partit
        type(t_mesh),   intent(in),    target :: mesh
        end subroutine

        subroutine ice_solve_high_order(ice, partit, mesh)
        USE MOD_ICE
        USE MOD_PARTIT
        USE MOD_PARSUP
        USE MOD_MESH
        type(t_ice),    intent(inout), target :: ice
        type(t_partit), intent(inout), target :: partit
        type(t_mesh),   intent(in),    target :: mesh
        end subroutine

        subroutine ice_solve_low_order(ice, partit, mesh)
        USE MOD_ICE
        USE MOD_PARTIT
        USE MOD_PARSUP
        USE MOD_MESH
        type(t_ice),    intent(inout), target :: ice
        type(t_partit), intent(inout), target :: partit
        type(t_mesh),   intent(in),    target :: mesh
        end subroutine

        subroutine ice_fem_fct(tr_array_id, ice, partit, mesh)
        USE MOD_ICE
        USE MOD_PARTIT
        USE MOD_PARSUP
        USE MOD_MESH
        integer   :: tr_array_id
        type(t_ice),    intent(inout), target :: ice
        type(t_partit), intent(inout), target :: partit
        type(t_mesh),   intent(in),    target :: mesh
        end subroutine

        subroutine ice_TG_rhs_div(ice, partit, mesh)
        USE MOD_ICE
        USE MOD_PARTIT
        USE MOD_PARSUP
        USE MOD_MESH
        type(t_ice),    intent(inout), target :: ice
        type(t_partit), intent(inout), target :: partit
        type(t_mesh),   intent(in),    target :: mesh
        end subroutine

        subroutine ice_TG_rhs(ice, partit, mesh)
        USE MOD_ICE
        USE MOD_PARTIT
        USE MOD_PARSUP
        USE MOD_MESH
        type(t_ice),    intent(inout), target :: ice
        type(t_partit), intent(inout), target :: partit
        type(t_mesh),   intent(in),    target :: mesh
        end subroutine

        subroutine ice_update_for_div(ice, partit, mesh)
        USE MOD_ICE
        USE MOD_PARTIT
        USE MOD_PARSUP
        USE MOD_MESH
        type(t_ice),    intent(inout), target :: ice
        type(t_partit), intent(inout), target :: partit
        type(t_mesh),   intent(in),    target :: mesh
        end subroutine
    end interface
end module
!
!
!_______________________________________________________________________________
! This file collect subroutines implementing FE-FCT
! advection scheme by Loehner et al.
! There is a tunable parameter ice_gamma_fct.
! Increasing it leads to positivity preserving solution.
!
! Driving routine is fct_ice_solve. It calles other routines
! that do low-order and figh order solutions and then combine them in a flux
! corrected way. Taylor-Galerkin scheme is used as a high-order one.
!
! The code is adapted from  FESOM
!
!
!_______________________________________________________________________________
subroutine ice_TG_rhs(ice, partit, mesh)
    use MOD_MESH
    USE MOD_PARTIT
    USE MOD_PARSUP
    USE MOD_ICE
    use o_PARAM
    USE g_CONFIG
    implicit none
    type(t_ice),    intent(inout), target :: ice
    type(t_partit), intent(inout), target :: partit
    type(t_mesh),   intent(in),    target :: mesh
    !___________________________________________________________________________
    real(kind=WP)   :: diff, entries(3),  um, vm, vol, dx(3), dy(3)
    integer         :: n, q, row, elem, elnodes(3)
    !___________________________________________________________________________

!! Juha: Use associate blocks instead of pointers to work around gpu memcopy issues with Cray
#include "associate_part_ass_combined.h"
#include "associate_mesh_ass_combined.h"

    associate (    &
#if defined (__oifs) || defined (__ifsinterface)
        ice_temp => ice%data(4)%values(:), &
        rhs_temp => ice%data(4)%values_rhs(:), &
#endif
        u_ice    => ice%uice(:), &
        v_ice    => ice%vice(:), &
        a_ice    => ice%data(1)%values(:), &
        m_ice    => ice%data(2)%values(:), &
        m_snow   => ice%data(3)%values(:), &
        rhs_a    => ice%data(1)%values_rhs(:), &
        rhs_m    => ice%data(2)%values_rhs(:), &
        rhs_ms   => ice%data(3)%values_rhs(:) )
    !___________________________________________________________________________
    ! Taylor-Galerkin (Lax-Wendroff) rhs
#ifndef ENABLE_OPENACC
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(n, q, row, elem, elnodes, diff, entries,  um, vm, vol, dx, dy)
!$OMP DO
#else
    !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT)
#endif
    DO row=1, myDim_nod2D
        rhs_m(row)=0._WP
        rhs_a(row)=0._WP
        rhs_ms(row)=0._WP
#if defined (__oifs) || defined (__ifsinterface)
        rhs_temp(row)=0._WP
#endif
    END DO

#ifndef ENABLE_OPENACC
!$OMP END DO
#else
    !$ACC END PARALLEL LOOP
#endif
    ! Velocities at nodes


#ifndef ENABLE_OPENACC
!$OMP DO
#else
    !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) private(n, q, row, elem, elnodes, diff, entries,  um, vm, vol, dx, dy)
#endif
    do elem=1,myDim_elem2D          !assembling rhs over elements
        elnodes=elem2D_nodes(:,elem)
        !_______________________________________________________________________
        ! if cavity element skip it
        if (ulevels(elem)>1) cycle

        !derivatives
        dx=gradient_sca(1:3,elem)
        dy=gradient_sca(4:6,elem)
        vol=elem_area(elem)
        !um=sum(U_ice(elnodes))/3.0_WP
        !vm=sum(V_ice(elnodes))/3.0_WP
        um=sum(U_ice(elnodes))
        vm=sum(V_ice(elnodes))

        !diffusivity
        diff=ice%ice_diff*sqrt(elem_area(elem)/scale_area)
        !$ACC LOOP SEQ
        DO n=1,3
            row=elnodes(n)
	    !$ACC LOOP SEQ
            DO q = 1,3
                !entries(q)= vol*dt*((dx(n)*um+dy(n)*vm)/3.0_WP - &
                !            diff*(dx(n)*dx(q)+ dy(n)*dy(q))- &
                !	       0.5*dt*(um*dx(n)+vm*dy(n))*(um*dx(q)+vm*dy(q)))
                entries(q)= vol*ice%ice_dt*((dx(n)*(um+u_ice(elnodes(q)))+ &
                            dy(n)*(vm+v_ice(elnodes(q))))/12.0_WP - &
                            diff*(dx(n)*dx(q)+ dy(n)*dy(q))- &
                            0.5_WP*ice%ice_dt*(um*dx(n)+vm*dy(n))*(um*dx(q)+vm*dy(q))/9.0_WP)
            END DO
	    !$ACC END LOOP
            rhs_m(row)=rhs_m(row)+sum(entries*m_ice(elnodes))
            rhs_a(row)=rhs_a(row)+sum(entries*a_ice(elnodes))
            rhs_ms(row)=rhs_ms(row)+sum(entries*m_snow(elnodes))
#if defined (__oifs) || defined (__ifsinterface)
            rhs_temp(row)=rhs_temp(row)+sum(entries*ice_temp(elnodes))
#endif
        END DO
	!$ACC END LOOP
    end do

#ifndef ENABLE_OPENACC
!$OMP END DO
!$OMP END PARALLEL
#else
    !$ACC END PARALLEL LOOP
#endif

    ! Juha: close associate blocks
    end associate
    end associate
    end associate

end subroutine ice_TG_rhs
!
!
!_______________________________________________________________________________
subroutine ice_fct_solve(ice, partit, mesh)
  USE MOD_ICE
  USE MOD_PARTIT
  USE MOD_PARSUP
  USE MOD_MESH
  use ice_fct_interfaces
  implicit none
  type(t_ice),    intent(inout), target :: ice
  type(t_partit), intent(inout), target :: partit
  type(t_mesh),   intent(in),    target :: mesh
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"
  !_____________________________________________________________________________
  ! Driving routine
  call ice_solve_high_order(ice, partit, mesh)   ! uses arrays of low-order solutions as temp
                                    ! storage. It should preceed the call of low
                                    ! order solution.
  call ice_solve_low_order(ice, partit, mesh)

  call ice_fem_fct(1, ice, partit, mesh)    ! m_ice
  call ice_fem_fct(2, ice, partit, mesh)    ! a_ice
  call ice_fem_fct(3, ice, partit, mesh)    ! m_snow

#if defined (__oifs) || defined (__ifsinterface)
  call ice_fem_fct(4, ice, partit, mesh)    ! ice_temp
#endif

end subroutine ice_fct_solve
!
!
!_______________________________________________________________________________
subroutine ice_solve_low_order(ice, partit, mesh)

    !============================
    ! Low-order solution
    !============================
    !
    ! It is assumed that m_ice, a_ice and m_snow from the previous time step
    ! are known at 1:myDim_nod2D+eDim_nod2D.
    ! We add diffusive contribution to the rhs. The diffusion operator
    ! is implemented as the difference between the consistent and lumped mass
    ! matrices acting on the field from the previous time step. The consistent
    ! mass matrix on the lhs is replaced with the lumped one.
    USE MOD_ICE
    USE MOD_PARTIT
    USE MOD_PARSUP
    USE MOD_MESH
    use g_comm_auto
    implicit none
    type(t_ice),    intent(inout), target :: ice
    type(t_partit), intent(inout), target :: partit
    type(t_mesh),   intent(in),    target :: mesh
    !___________________________________________________________________________
    integer       :: row, clo, clo2, cn, location(100)
    real(kind=WP) :: gamma
    !___________________________________________________________________________

    ! Juha: Use associate blocks instead of pointers to work around gpu memcopy issues with Cray
#include "associate_part_ass_combined.h"
#include "associate_mesh_ass_combined.h"

    associate( &
#if defined (__oifs) || defined (__ifsinterface)
    ice_temp     => ice%data(4)%values(:), &
    rhs_temp     => ice%data(4)%values_rhs(:), &
    m_templ      => ice%data(4)%valuesl(:), &
#endif
    a_ice        => ice%data(1)%values(:), &
    m_ice        => ice%data(2)%values(:), &
    m_snow       => ice%data(3)%values(:), &
    rhs_a        => ice%data(1)%values_rhs(:), &
    rhs_m        => ice%data(2)%values_rhs(:), &
    rhs_ms       => ice%data(3)%values_rhs(:), &
    a_icel       => ice%data(1)%valuesl(:), &
    m_icel       => ice%data(2)%valuesl(:), &
    m_snowl      => ice%data(3)%valuesl(:), &
    mass_matrix  => ice%work%fct_massmatrix(:) &
    )


    !___________________________________________________________________________
    gamma=ice%ice_gamma_fct         ! Added diffusivity parameter
                                ! Adjust it to ensure posivity of solution

#ifndef ENABLE_OPENACC
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(row, clo, clo2, cn, location)
#else
    !$ACC PARALLEL LOOP GANG VECTOR PRESENT(ssh_stiff, ssh_stiff%rowptr) PRIVATE(location) DEFAULT(PRESENT)
#endif
    do row=1,myDim_nod2D
        !_______________________________________________________________________
        ! if there is cavity no ice fxt low order
        if (ulevels_nod2D(row)>1) cycle

        !_______________________________________________________________________
        clo=ssh_stiff%rowptr(row)-ssh_stiff%rowptr(1)+1
        clo2=ssh_stiff%rowptr(row+1)-ssh_stiff%rowptr(1)
        cn=clo2-clo+1
        location(1:cn)=nn_pos(1:cn, row)
        m_icel(row)=(rhs_m(row)+gamma*sum(mass_matrix(clo:clo2)* &
                    m_ice(location(1:cn))))/area(1,row) + &
                    (1.0_WP-gamma)*m_ice(row)
        a_icel(row)=(rhs_a(row)+gamma*sum(mass_matrix(clo:clo2)* &
                    a_ice(location(1:cn))))/area(1,row) + &
                    (1.0_WP-gamma)*a_ice(row)
        m_snowl(row)=(rhs_ms(row)+gamma*sum(mass_matrix(clo:clo2)* &
                    m_snow(location(1:cn))))/area(1,row) + &
                    (1.0_WP-gamma)*m_snow(row)
#if defined (__oifs) || defined (__ifsinterface)
        m_templ(row)=(rhs_temp(row)+gamma*sum(mass_matrix(clo:clo2)* &
                  ice_temp(location(1:cn))))/area(1,row) + &
                  (1.0_WP-gamma)*ice_temp(row)
#endif
    end do
#ifndef ENABLE_OPENACC
!$OMP END PARALLEL DO
#else 
   !$ACC END PARALLEL LOOP
#endif
    ! Low-order solution must be known to neighbours
    call exchange_nod(m_icel,a_icel,m_snowl, partit, luse_g2g = .true.)
#if defined (__oifs) || defined (__ifsinterface)
    call exchange_nod(m_templ, partit, luse_g2g = .true.)
#endif

#ifndef ENABLE_OPENACC
!$OMP BARRIER
#endif

    ! Juha: close associate blocks
    end associate
    end associate
    end associate

end subroutine ice_solve_low_order
!
!
!_______________________________________________________________________________
subroutine ice_solve_high_order(ice, partit, mesh)
    USE MOD_ICE
    USE MOD_TRACER
    USE MOD_PARTIT
    USE MOD_PARSUP
    USE MOD_MESH
    use o_PARAM
    use g_comm_auto
    implicit none
    type(t_ice)   , intent(inout), target :: ice
    type(t_partit), intent(inout), target :: partit
    type(t_mesh)  , intent(in)   , target :: mesh
    !___________________________________________________________________________
    integer                               :: n,clo,clo2,cn,location(100),row
    real(kind=WP)                         :: rhs_new
    integer                               :: num_iter_solve=3
    !___________________________________________________________________________

    ! Juha: Use associate blocks instead of pointers to work around gpu memcopy issues with Cray
#include "associate_part_ass_combined.h"
#include "associate_mesh_ass_combined.h"

    associate( &
#if defined (__oifs) || defined (__ifsinterface)
    rhs_temp     => ice%data(4)%values_rhs(:), &
    m_templ      => ice%data(4)%valuesl(:), &
    dm_temp      => ice%data(4)%dvalues(:), &
#endif
    rhs_a        => ice%data(1)%values_rhs(:), &
    rhs_m        => ice%data(2)%values_rhs(:), &
    rhs_ms       => ice%data(3)%values_rhs(:), &
    a_icel       => ice%data(1)%valuesl(:), &
    m_icel       => ice%data(2)%valuesl(:), &
    m_snowl      => ice%data(3)%valuesl(:), &
    da_ice       => ice%data(1)%dvalues(:), &
    dm_ice       => ice%data(2)%dvalues(:), &
    dm_snow      => ice%data(3)%dvalues(:), &
    mass_matrix  => ice%work%fct_massmatrix(:)  &
    )

    !___________________________________________________________________________
    ! Does Taylor-Galerkin solution
    !
    !the first approximation
#ifndef ENABLE_OPENACC
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(row)
#else
    !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT)
#endif
    do row=1,myDim_nod2D
        ! if cavity node skip it
        if (ulevels_nod2d(row)>1) cycle

        dm_ice(row)=rhs_m(row)/area(1,row)
        da_ice(row)=rhs_a(row)/area(1,row)
        dm_snow(row)=rhs_ms(row)/area(1,row)
#if defined (__oifs) || defined (__ifsinterface)
        dm_temp(row)=rhs_temp(row)/area(1,row)
#endif
    end do
#ifndef ENABLE_OPENACC
!$OMP END PARALLEL DO
#else
    !$ACC END PARALLEL LOOP
#endif
    call exchange_nod(dm_ice, da_ice, dm_snow, partit, luse_g2g = .true.)
#if defined (__oifs) || defined (__ifsinterface)
    call exchange_nod(dm_temp, partit, luse_g2g = .true.)
#endif /* (__oifs) */
#ifndef ENABLE_OPENACC
!$OMP BARRIER
#endif
    !___________________________________________________________________________
    !iterate

    do n=1,num_iter_solve-1
#ifndef ENABLE_OPENACC
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(n, clo, clo2, cn, location, row, rhs_new)
!$OMP DO
#else
        !$ACC PARALLEL LOOP GANG VECTOR PRESENT(ssh_stiff, ssh_stiff%rowptr) PRIVATE(location) DEFAULT(PRESENT)
#endif
        do row=1,myDim_nod2D
            ! if cavity node skip it
            if (ulevels_nod2d(row)>1) cycle
            !___________________________________________________________________
            clo  = ssh_stiff%rowptr(row)-ssh_stiff%rowptr(1)+1
            clo2 = ssh_stiff%rowptr(row+1)-ssh_stiff%rowptr(1)
            cn   = clo2-clo+1
            location(1:cn)=nn_pos(1:cn,row)
            !___________________________________________________________________
            rhs_new     = rhs_m(row) - sum(mass_matrix(clo:clo2)*dm_ice(location(1:cn)))
            m_icel(row) = dm_ice(row)+rhs_new/area(1,row)
            rhs_new     = rhs_a(row) - sum(mass_matrix(clo:clo2)*da_ice(location(1:cn)))
            a_icel(row) = da_ice(row)+rhs_new/area(1,row)
            rhs_new     = rhs_ms(row) - sum(mass_matrix(clo:clo2)*dm_snow(location(1:cn)))
            m_snowl(row)= dm_snow(row)+rhs_new/area(1,row)
#if defined (__oifs) || defined (__ifsinterface)
            rhs_new     = rhs_temp(row) - sum(mass_matrix(clo:clo2)*dm_temp(location(1:cn)))
            m_templ(row)= dm_temp(row)+rhs_new/area(1,row)
#endif
        end do
#ifndef ENABLE_OPENACC
!$OMP END DO
#else
        !$ACC END PARALLEL LOOP
#endif
        !_______________________________________________________________________
#ifndef ENABLE_OPENACC
!$OMP DO
#else
        !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT)
#endif
        do row=1,myDim_nod2D
            ! if cavity node skip it
            if (ulevels_nod2d(row)>1) cycle
            dm_ice(row)=m_icel(row)
            da_ice(row)=a_icel(row)
            dm_snow(row)=m_snowl(row)
#if defined (__oifs) || defined (__ifsinterface)
            dm_temp(row)=m_templ(row)
#endif
        end do
#ifndef ENABLE_OPENACC
!$OMP END DO
!$OMP END PARALLEL
#else
        !$ACC END PARALLEL LOOP
#endif
        !_______________________________________________________________________
        call exchange_nod(dm_ice, da_ice, dm_snow, partit, luse_g2g = .true.)
#if defined (__oifs) || defined (__ifsinterface)
        call exchange_nod(dm_temp, partit, luse_g2g = .true.)
#endif /* (__oifs) */
#ifndef ENABLE_OPENACC
!$OMP BARRIER
#endif
    end do

    ! Juha: close associate blocks
    end associate
    end associate
    end associate

end subroutine ice_solve_high_order
!
!
!_______________________________________________________________________________
! Flux corrected transport algorithm for tracer advection
! It is based on Loehner et al. (Finite-element flux-corrected
! transport (FEM-FCT) for the Euler and Navier-Stokes equation,
! Int. J. Numer. Meth. Fluids, 7 (1987), 1093--1109) as described by Kuzmin and
! Turek. (kuzmin@math.uni-dortmund.de)
subroutine ice_fem_fct(tr_array_id, ice, partit, mesh)
    USE MOD_ICE
    USE MOD_PARTIT
    USE MOD_PARSUP
    USE MOD_MESH
    use o_PARAM
    use g_comm_auto
    implicit none
    type(t_ice)   , intent(inout), target :: ice
    type(t_partit), intent(inout), target :: partit
    type(t_mesh)  , intent(in)   , target :: mesh
    !___________________________________________________________________________
    integer       :: tr_array_id
    integer       :: icoef(3,3), n, q, elem, elnodes(3), row
    real(kind=WP) :: vol, flux, ae, gamma
    !___________________________________________________________________________

    !! Juha: Use associate blocks instead of pointers to work around gpu memcopy issues with Cray
#include "associate_part_ass_combined.h"
#include "associate_mesh_ass_combined.h"

    associate( &
#if defined (__oifs) || defined (__ifsinterface)
    ice_temp  => ice%data(4)%values(:), &
    m_templ   => ice%data(4)%valuesl(:), &
    dm_temp   => ice%data(4)%dvalues(:), &
#endif
    a_ice     => ice%data(1)%values(:), &
    m_ice     => ice%data(2)%values(:), &
    m_snow    => ice%data(3)%values(:), &
    a_icel    => ice%data(1)%valuesl(:), &
    m_icel    => ice%data(2)%valuesl(:), &
    m_snowl   => ice%data(3)%valuesl(:), &
    da_ice    => ice%data(1)%dvalues(:), &
    dm_ice    => ice%data(2)%dvalues(:), &
    dm_snow   => ice%data(3)%dvalues(:), &
    icefluxes => ice%work%fct_fluxes(:,:), &
    icepplus  => ice%work%fct_plus(:), &
    icepminus => ice%work%fct_minus(:), &
    tmax      => ice%work%fct_tmax(:), &
    tmin      => ice%work%fct_tmin(:) &
    )


    !___________________________________________________________________________
    ! It should coinside with gamma in ts_solve_low_order
    gamma=ice%ice_gamma_fct

    !___________________________________________________________________________
    ! Compute elemental antidiffusive fluxes to nodes
    ! This is the most unpleasant part ---
    ! it takes memory and time. For every element
    ! we need its antidiffusive contribution to
    ! each of its 3 nodes
#ifndef ENABLE_OPENACC
!$OMP PARALLEL DO
#else
    !$ACC DATA CREATE(icoef, elnodes)
    !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT)
#endif
    do n = 1, myDim_nod2D + eDim_nod2D
        tmax(n) = 0.0_WP
        tmin(n) = 0.0_WP
    end do
#ifndef ENABLE_OPENACC
!$OMP END PARALLEL  DO
#else
    !$ACC END PARALLEL LOOP
#endif
    ! Auxiliary elemental operator (mass matrix- lumped mass matrix)

    !$ACC KERNELS
    icoef = 1
    !$ACC END KERNELS
    !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT)
    do n=1,3   ! three upper nodes
        ! Cycle over rows  row=elnodes(n)
        icoef(n,n)=-2
    end do
    !$ACC END PARALLEL LOOP


#ifndef ENABLE_OPENACC
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(n, q, elem, elnodes, row, vol, flux, ae)
!$OMP DO
#else
    !$ACC PARALLEL LOOP GANG VECTOR PRIVATE(elnodes) DEFAULT(PRESENT)
#endif
    do elem=1, myDim_elem2D
        !_______________________________________________________________________
        elnodes=elem2D_nodes(:,elem)

        !_______________________________________________________________________
        ! if cavity cycle over
        if(ulevels(elem)>1) cycle !LK89140

        !_______________________________________________________________________
        vol=elem_area(elem)
        if (tr_array_id==1) then
            do q=1,3
                icefluxes(elem,q)=-sum(icoef(:,q)*(gamma*m_ice(elnodes) + &
                            dm_ice(elnodes)))*(vol/area(1,elnodes(q)))/12.0_WP
            end do
        end if

        if (tr_array_id==2) then
            do q=1,3
                icefluxes(elem,q)=-sum(icoef(:,q)*(gamma*a_ice(elnodes) + &
                            da_ice(elnodes)))*(vol/area(1,elnodes(q)))/12.0_WP
            end do
        end if

        if (tr_array_id==3) then
            do q=1,3
                icefluxes(elem,q)=-sum(icoef(:,q)*(gamma*m_snow(elnodes) + &
                            dm_snow(elnodes)))*(vol/area(1,elnodes(q)))/12.0_WP
            end do
        end if

#if defined (__oifs) || defined (__ifsinterface)
        if (tr_array_id==4) then
            do q=1,3
                icefluxes(elem,q)=-sum(icoef(:,q)*(gamma*ice_temp(elnodes) + &
                            dm_temp(elnodes)))*(vol/area(1,elnodes(q)))/12.0_WP
            end do
        end if
#endif
    end do
#ifndef ENABLE_OPENACC
!$OMP END DO
#else
    !$ACC END PARALLEL LOOP
#endif
    !___________________________________________________________________________
    ! Screening the low-order solution
    ! TO BE ADDED IF FOUND NECESSARY
    ! Screening means comparing low-order solutions with the
    ! solution on the previous time step and using whichever
    ! is greater/smaller in computations of max/min below

    !___________________________________________________________________________
    ! Cluster min/max
    if (tr_array_id==1) then
#ifndef ENABLE_OPENACC
!$OMP DO
#else
        !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT)
#endif
        do row=1, myDim_nod2D
            if (ulevels_nod2d(row)>1) cycle
            n=nn_num(row)
            tmax(row)=maxval(m_icel(nn_pos(1:n,row)))
            tmin(row)=minval(m_icel(nn_pos(1:n,row)))
            ! Admissible increments
            tmax(row)=tmax(row)-m_icel(row)
            tmin(row)=tmin(row)-m_icel(row)
        end do
#ifndef ENABLE_OPENACC
!$OMP END DO
#else
        !$ACC END PARALLEL LOOP
#endif
    end if

    if (tr_array_id==2) then
#ifndef ENABLE_OPENACC
!$OMP DO
#else
        !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT)
#endif
        do row=1, myDim_nod2D
            if (ulevels_nod2d(row)>1) cycle
            n=nn_num(row)
            tmax(row)=maxval(a_icel(nn_pos(1:n,row)))
            tmin(row)=minval(a_icel(nn_pos(1:n,row)))
            ! Admissible increments
            tmax(row)=tmax(row)-a_icel(row)
            tmin(row)=tmin(row)-a_icel(row)
        end do
#ifndef ENABLE_OPENACC
!$OMP END DO
#else
        !$ACC END PARALLEL LOOP
#endif
    end if

    if (tr_array_id==3) then
#ifndef ENABLE_OPENACC
!$OMP DO
#else
        !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT)
#endif
        do row=1, myDim_nod2D
            if (ulevels_nod2d(row)>1) cycle
            n=nn_num(row)
            tmax(row)=maxval(m_snowl(nn_pos(1:n,row)))
            tmin(row)=minval(m_snowl(nn_pos(1:n,row)))
            ! Admissible increments
            tmax(row)=tmax(row)-m_snowl(row)
            tmin(row)=tmin(row)-m_snowl(row)
        end do
#ifndef ENABLE_OPENACC
!$OMP END DO
#else
        !$ACC END PARALLEL LOOP
#endif
    end if

#if defined (__oifs) || defined (__ifsinterface)
    if (tr_array_id==4) then
#ifndef ENABLE_OPENACC
!$OMP DO
#else
        !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT)
#endif
        do row=1, myDim_nod2D
            if (ulevels_nod2d(row)>1) cycle
            n=nn_num(row)
            tmax(row)=maxval(m_templ(nn_pos(1:n,row)))
            tmin(row)=minval(m_templ(nn_pos(1:n,row)))
            ! Admissible increments
            tmax(row)=tmax(row)-m_templ(row)
            tmin(row)=tmin(row)-m_templ(row)
        end do
#ifndef ENABLE_OPENACC
!$OMP END DO
#else
        !$ACC END PARALLEL LOOP
#endif
    end if
#endif

    !___________________________________________________________________________
    ! Sums of positive/negative fluxes to node row
#ifndef ENABLE_OPENACC
!$OMP DO
#else
    !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT)
#endif
    do n=1, myDim_nod2D+eDim_nod2D
       icepplus (n)=0._WP
       icepminus(n)=0._WP
    end do
#ifndef ENABLE_OPENACC
!$OMP END DO
#else
    !$ACC END PARALLEL LOOP
#endif

#ifndef ENABLE_OPENACC
!$OMP DO
#else
#if !defined(DISABLE_OPENACC_ATOMICS)
    !$ACC PARALLEL LOOP GANG VECTOR PRIVATE(elnodes) DEFAULT(PRESENT)
#else
    !$ACC UPDATE SELF(icefluxes, icepplus, icepminus)
#endif
#endif
    do elem=1, myDim_elem2D
        ! if cavity cycle over
        if(ulevels(elem)>1) cycle !LK89140

        !_______________________________________________________________________
        elnodes=elem2D_nodes(:,elem)
        do q=1,3
            n=elnodes(q)
            flux=icefluxes(elem,q)
#ifndef ENABLE_OPENACC
#if defined(_OPENMP)  && !defined(__openmp_reproducible)
        call omp_set_lock  (partit%plock(n))
#else
!$OMP ORDERED
#endif
#endif
            if (flux>0) then
#if !defined(DISABLE_OPENACC_ATOMICS)
                !$ACC ATOMIC UPDATE
#endif
                icepplus(n)=icepplus(n)+flux
            else
#if !defined(DISABLE_OPENACC_ATOMICS)
                !$ACC ATOMIC UPDATE
#endif
                icepminus(n)=icepminus(n)+flux
            end if
#ifndef ENABLE_OPENACC
#if defined(_OPENMP)  && !defined(__openmp_reproducible)
        call omp_unset_lock(partit%plock(n))
#else
!$OMP END ORDERED
#endif
#endif
        end do
    end do
#ifndef ENABLE_OPENACC
!$OMP END DO
#else
#if !defined(DISABLE_OPENACC_ATOMICS)
    !$ACC END PARALLEL LOOP
#else
    !$ACC UPDATE DEVICE(icepplus, icepminus)
#endif
#endif


    !___________________________________________________________________________
    ! The least upper bound for the correction factors
#ifndef ENABLE_OPENACC
!$OMP DO
#else
    !$ACC PARALLEL LOOP GANG VECTOR PRESENT(icepplus, icepminus) DEFAULT(PRESENT)
#endif
    do n=1,myDim_nod2D
        ! if cavity cycle over
        if(ulevels_nod2D(n)>1) cycle !LK89140

        flux=icepplus(n)
        if (abs(flux)>0) then
            icepplus(n)=min(1.0_WP,tmax(n)/max(flux,1.e-12))
        else
            icepplus(n)=0._WP
        end if

        flux=icepminus(n)
        if (abs(flux)>0) then
            icepminus(n)=min(1.0_WP,tmin(n)/min(flux,-1.e-12))
        else
            icepminus(n)=0._WP
        end if
    end do
#ifndef ENABLE_OPENACC
!$OMP END DO
#else
    !$ACC END PARALLEL LOOP
#endif 
   ! pminus and pplus are to be known to neighbouting PE
!$ACC wait

#if defined(_OPENMP)
!$OMP MASTER
#endif
    call exchange_nod(icepminus, icepplus, partit, luse_g2g = .true.)
#if defined(_OPENMP)
!$OMP END MASTER
!$OMP BARRIER
#endif
    !___________________________________________________________________________
    ! Limiting
#ifndef ENABLE_OPENACC
!$OMP DO
#else
    !$ACC PARALLEL LOOP GANG VECTOR PRESENT(icepplus, icepminus) PRIVATE(elnodes) DEFAULT(PRESENT)
#endif
    do elem=1, myDim_elem2D
        ! if cavity cycle over
        if(ulevels(elem)>1) cycle !LK89140

        !_______________________________________________________________________
        elnodes=elem2D_nodes(:,elem)
        ae=1.0_WP
        do q=1,3
            n=elnodes(q)
            flux=icefluxes(elem,q)
            if(flux>=0._WP) ae=min(ae,icepplus(n))
            if(flux<0._WP) ae=min(ae,icepminus(n))
        end do
        icefluxes(elem,:)=ae*icefluxes(elem,:)
    end do
#ifndef ENABLE_OPENACC
!$OMP END DO
#else
    !$ACC END PARALLEL LOOP
#endif
    !___________________________________________________________________________
    ! Update the solution
    if(tr_array_id==1) then
#ifndef ENABLE_OPENACC
!$OMP DO
#else
        !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT)
#endif
        do n=1,myDim_nod2D
            if(ulevels_nod2D(n)>1) cycle !LK89140
            m_ice(n)=m_icel(n)
        end do
#ifndef ENABLE_OPENACC
!$OMP END DO
#else
        !$ACC END PARALLEL LOOP
#endif

#ifndef ENABLE_OPENACC
!$OMP DO
#else
#if !defined(DISABLE_OPENACC_ATOMICS)
        !$ACC PARALLEL LOOP GANG VECTOR PRIVATE(elnodes) DEFAULT(PRESENT)
#else
        !$ACC UPDATE SELF(m_ice, icefluxes)
#endif
#endif
        do elem=1, myDim_elem2D
            ! if cavity cycle over
            if(ulevels(elem)>1) cycle !LK89140

            elnodes=elem2D_nodes(:,elem)
            do q=1,3
                n=elnodes(q)
#ifndef ENABLE_OPENACC
#if defined(_OPENMP)  && !defined(__openmp_reproducible)
                call omp_set_lock  (partit%plock(n))
#else
!$OMP ORDERED
#endif
#endif
#if !defined(DISABLE_OPENACC_ATOMICS)
                !$ACC ATOMIC UPDATE
#endif
                m_ice(n)=m_ice(n)+icefluxes(elem,q)
#ifndef ENABLE_OPENACC
#if defined(_OPENMP)  && !defined(__openmp_reproducible)
                call omp_unset_lock(partit%plock(n))
#else
!$OMP END ORDERED
#endif
#endif
            end do
        end do
#ifndef ENABLE_OPENACC
!$OMP END DO
#else
#if !defined(DISABLE_OPENACC_ATOMICS)
        !$ACC END PARALLEL LOOP
#else
        !$ACC UPDATE DEVICE(m_ice)
#endif
#endif
    end if

    if(tr_array_id==2) then
#ifndef ENABLE_OPENACC
!$OMP DO
#else
        !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT)
#endif
        do n=1,myDim_nod2D
            if(ulevels_nod2D(n)>1) cycle !LK89140
            a_ice(n)=a_icel(n)
        end do
#ifndef ENABLE_OPENACC
!$OMP END DO
!$OMP DO
#else
        !$ACC END PARALLEL LOOP
#if !defined(DISABLE_OPENACC_ATOMICS)
        !$ACC PARALLEL LOOP GANG VECTOR PRIVATE(elnodes) DEFAULT(PRESENT)
#else
        !$ACC UPDATE SELF(a_ice, icefluxes)
#endif
#endif
        do elem=1, myDim_elem2D
            ! if cavity cycle over
            if(ulevels(elem)>1) cycle !LK89140

            elnodes=elem2D_nodes(:,elem)
            do q=1,3
                n=elnodes(q)
#ifndef ENABLE_OPENACC
#if defined(_OPENMP)  && !defined(__openmp_reproducible)
                call omp_set_lock  (partit%plock(n))
#else
!$OMP ORDERED
#endif
#endif
#if !defined(DISABLE_OPENACC_ATOMICS)
                !$ACC ATOMIC UPDATE
#endif
                a_ice(n)=a_ice(n)+icefluxes(elem,q)
#ifndef ENABLE_OPENACC
#if defined(_OPENMP) && !defined(__openmp_reproducible)
                call omp_unset_lock(partit%plock(n))
#else
!$OMP END ORDERED
#endif
#endif
            end do
        end do
#ifndef ENABLE_OPENACC
!$OMP END DO
#else
#if !defined(DISABLE_OPENACC_ATOMICS)
        !$ACC END PARALLEL LOOP
#else
        !$ACC UPDATE DEVICE(a_ice)
#endif
#endif
    end if

    if(tr_array_id==3) then
#ifndef ENABLE_OPENACC
!$OMP DO
#else
        !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT)
#endif
        do n=1,myDim_nod2D
            if(ulevels_nod2D(n)>1) cycle !LK89140
            m_snow(n)=m_snowl(n)
        end do
#ifndef ENABLE_OPENACC
!$OMP END DO
!$OMP DO
#else
        !$ACC END PARALLEL LOOP
#if !defined(DISABLE_OPENACC_ATOMICS)
        !$ACC PARALLEL LOOP GANG VECTOR PRIVATE(elnodes) DEFAULT(PRESENT)
#else
        !$ACC UPDATE SELF(m_snow, icefluxes)
#endif
#endif
        do elem=1, myDim_elem2D
            ! if cavity cycle over
            if(ulevels(elem)>1) cycle !LK89140

            elnodes=elem2D_nodes(:,elem)
            do q=1,3
                n=elnodes(q)
#ifndef ENABLE_OPENACC
#if defined(_OPENMP)  && !defined(__openmp_reproducible)
                call omp_set_lock  (partit%plock(n))
#else
!$OMP ORDERED
#endif
#endif
#if !defined(DISABLE_OPENACC_ATOMICS)
                !$ACC ATOMIC UPDATE
#endif
                m_snow(n)=m_snow(n)+icefluxes(elem,q)
#ifndef ENABLE_OPENACC
#if defined(_OPENMP) && !defined(__openmp_reproducible)
                call omp_unset_lock(partit%plock(n))
#else
!$OMP END ORDERED
#endif
#endif
            end do
        end do
#ifndef ENABLE_OPENACC
!$OMP END DO
#else
#if !defined(DISABLE_OPENACC_ATOMICS)
        !$ACC END PARALLEL LOOP
#else
        !$ACC UPDATE DEVICE(m_snow)
#endif
#endif
    end if

#if defined (__oifs) || defined (__ifsinterface)
    if(tr_array_id==4) then
#ifndef ENABLE_OPENACC
!$OMP DO
#else
        !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT)
#endif
        do n=1,myDim_nod2D
            if(ulevels_nod2D(n)>1) cycle !LK89140
            ice_temp(n)=m_templ(n)
        end do
#ifndef ENABLE_OPENACC
!$OMP END DO
!$OMP DO
#else
        !$ACC END PARALLEL LOOP
#endif
#if !defined(DISABLE_OPENACC_ATOMICS)
        !$ACC PARALLEL LOOP GANG VECTOR PRIVATE(elnodes) DEFAULT(PRESENT)
#else
        !$ACC UPDATE SELF(ice_temp, icefluxes)
#endif
        do elem=1, myDim_elem2D
            ! if cavity cycle over
            if(ulevels(elem)>1) cycle !LK89140

            elnodes=elem2D_nodes(:,elem)
            do q=1,3
                n=elnodes(q)
#ifndef ENABLE_OPENACC
#if defined(_OPENMP)  && !defined(__openmp_reproducible)
                call omp_set_lock  (partit%plock(n))
#else
!$OMP ORDERED
#endif
#endif
#if !defined(DISABLE_OPENACC_ATOMICS)
                !$ACC ATOMIC UPDATE
#endif
                ice_temp(n)=ice_temp(n)+icefluxes(elem,q)
#ifndef ENABLE_OPENACC
#if defined(_OPENMP)  && !defined(__openmp_reproducible)
                call omp_unset_lock(partit%plock(n))
#else
!$OMP END ORDERED
#endif
#endif
            end do
        end do
#if !defined(DISABLE_OPENACC_ATOMICS)
        !$ACC END PARALLEL LOOP
#else
        !$ACC UPDATE DEVICE(ice_temp)
#endif
#ifndef ENABLE_OPENACC
!$OMP END DO
#endif
    end if
#endif
#ifndef ENABLE_OPENACC
!$OMP END PARALLEL
#endif
    call exchange_nod(m_ice, a_ice, m_snow, partit, luse_g2g = .true.)
#if defined (__oifs) || defined (__ifsinterface)
    call exchange_nod(ice_temp, partit, luse_g2g = .true.)
#endif

!$ACC END DATA

!$OMP BARRIER

    !Juha: for associate blocks
    end associate
    end associate
    end associate

end subroutine ice_fem_fct
!
!
!_______________________________________________________________________________
! Used in ice_fct inherited from FESOM
SUBROUTINE ice_mass_matrix_fill(ice, partit, mesh)
    USE MOD_ICE
    USE MOD_PARTIT
    USE MOD_PARSUP
    USE MOD_MESH
    implicit none
    type(t_ice)   , intent(inout), target :: ice
    type(t_partit), intent(inout), target :: partit
    type(t_mesh)  , intent(in)   , target :: mesh
    !___________________________________________________________________________
    integer                             :: n, k, row
    integer                             :: elem, elnodes(3), q, offset, ipos
    real(kind=WP)                       :: aa
    integer                             :: flag=0, iflag=0
    !___________________________________________________________________________
    ! pointer on necessary derived types
    real(kind=WP), dimension(:), pointer  :: mass_matrix
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"
    mass_matrix => ice%work%fct_massmatrix(:)
    !
    ! a)
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(n, k, row, elem, elnodes, q, offset, ipos, aa)
!$OMP DO
    DO elem=1,myDim_elem2D
        elnodes=elem2D_nodes(:,elem)

        !_______________________________________________________________________
        do n=1,3
            row=elnodes(n)
            if(row>myDim_nod2D) cycle
            !___________________________________________________________________
            ! Global-to-local neighbourhood correspondence
            ! we have to modify col_pos construction for OMP compatibility. The MPI version might become a bit slower :(
            ! loop over number of neghbouring nodes of node-row
            offset=ssh_stiff%rowptr(row)-ssh_stiff%rowptr(1)
            do q=1, 3
               !_______________________________________________________________
               ! if element is cavity cycle over
               if(ulevels(elem)>1) cycle
               do k=1, nn_num(row)
                  if (nn_pos(k,row)==elnodes(q)) then
                     ipos=offset+k
                     exit
                  end if
                  if (k==nn_num(row)) write(*,*) 'FATAL ERROR'
               end do
#if defined(_OPENMP)  && !defined(__openmp_reproducible)
               call omp_set_lock  (partit%plock(row)) ! it shall be sufficient to block writing into the same row of SSH_stiff
#else
!$OMP ORDERED
#endif
               mass_matrix(ipos)=mass_matrix(ipos)+elem_area(elem)/12.0_WP
               if(q==n) then
                   mass_matrix(ipos)=mass_matrix(ipos)+elem_area(elem)/12.0_WP
               end if
#if defined(_OPENMP)  && !defined(__openmp_reproducible)
               call omp_unset_lock(partit%plock(row))
#else
!$OMP END ORDERED
#endif
           END DO
        end do
    END DO
!$OMP END DO
    ! TEST: area==sum of row entries in mass_matrix:
!$OMP DO
    DO q=1,myDim_nod2D
        ! if cavity cycle over
        if(ulevels_nod2d(q)>1) cycle

        !_______________________________________________________________________
        offset=ssh_stiff%rowptr(q)-ssh_stiff%rowptr(1)+1
        n=ssh_stiff%rowptr(q+1)-ssh_stiff%rowptr(1)
        aa=sum(mass_matrix(offset:n))
        !!PS if(abs(area(1,q)-aa)>.1_WP) then
        if(abs(area(ulevels_nod2d(q),q)-aa)>.1_WP) then
!$OMP CRITICAL
            iflag=q
            flag=1
!$OMP END CRITICAL
        endif
    END DO
!$OMP END DO
!$OMP END PARALLEL
    if(flag>0) then
        offset=ssh_stiff%rowptr(iflag)-ssh_stiff%rowptr(1)+1
        n=ssh_stiff%rowptr(iflag+1)-ssh_stiff%rowptr(1)
#if !defined(__openmp_reproducible)
        aa=0
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(row) REDUCTION(+:aa)
!$OMP DO
        do row=offset, n
           aa=aa+mass_matrix(row)
        end do
!$OMP END DO
!$OMP END PARALLEL
#else
     aa = sum(mass_matrix(offset:n))
#endif
        write(*,*) '#### MASS MATRIX PROBLEM', mype, iflag, aa, area(1,iflag), ulevels_nod2D(iflag)
    endif
END SUBROUTINE ice_mass_matrix_fill
!
!
!_______________________________________________________________________________
subroutine ice_TG_rhs_div(ice, partit, mesh)
    USE MOD_ICE
    USE MOD_PARTIT
    USE MOD_PARSUP
    USE MOD_MESH
    use o_PARAM
    USE g_CONFIG
    implicit none
    type(t_ice)   , intent(inout), target :: ice
    type(t_partit), intent(inout), target :: partit
    type(t_mesh)  , intent(in)   , target :: mesh
    !___________________________________________________________________________
    real(kind=WP)            :: diff, entries(3),  um, vm, vol, dx(3), dy(3), tmp_sum
    integer                  :: n, q, row, elem, elnodes(3)
    real(kind=WP)            :: c1, c2, c3, c4, cx1, cx2, cx3, cx4, entries2(3)
    !___________________________________________________________________________
    real(kind=WP), dimension(:), pointer  :: rhs_a, rhs_m, rhs_ms
!! Juha: Use associate blocks instead of pointers to work around gpu memcopy issues with Cray
#include "associate_part_ass_combined.h"
#include "associate_mesh_ass_combined.h"

    rhs_a       => ice%data(1)%values_rhs(:)  !! Juha: leavign these here because for some reason acc atomic updates wouldnt work, look into later
    rhs_m       => ice%data(2)%values_rhs(:)
    rhs_ms      => ice%data(3)%values_rhs(:)
    
    associate(  &
#if defined (__oifs) || defined (__ifsinterface)
    ice_temp    => ice%data(4)%values(:), &
    rhs_temp    => ice%data(4)%values_rhs(:), &
    rhs_tempdiv => ice%data(4)%values_div_rhs(:), & 
#endif
    u_ice       => ice%uice(:), &
    v_ice       => ice%vice(:), &
    a_ice       => ice%data(1)%values(:), &
    m_ice       => ice%data(2)%values(:), &
    m_snow      => ice%data(3)%values(:), &
    !rhs_a       => ice%data(1)%values_rhs(:), &
    !rhs_m       => ice%data(2)%values_rhs(:), &
    !rhs_ms      => ice%data(3)%values_rhs(:), &
    rhs_adiv    => ice%data(1)%values_div_rhs(:), &
    rhs_mdiv    => ice%data(2)%values_div_rhs(:), &
    rhs_msdiv   => ice%data(3)%values_div_rhs(:) & 
    )
    


    !___________________________________________________________________________
    ! Computes the rhs in a Taylor-Galerkin way (with upwind type of
    ! correction for the advection operator)
    ! In this version I tr to split divergent term off, so that FCT works without it.
#ifndef ENABLE_OPENACC
!$OMP PARALLEL DO
#else
    !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT)
#endif
    do row=1, myDim_nod2D
                    !! row=myList_nod2D(m)
        rhs_m(row)=0.0_WP
        rhs_a(row)=0.0_WP
        rhs_ms(row)=0.0_WP
#if defined (__oifs) || defined (__ifsinterface)
        rhs_temp(row)=0.0_WP
#endif
        rhs_mdiv(row)=0.0_WP
        rhs_adiv(row)=0.0_WP
        rhs_msdiv(row)=0.0_WP
#if defined (__oifs) || defined (__ifsinterface)
        rhs_tempdiv(row)=0.0_WP
#endif
    end do
#ifndef ENABLE_OPENACC
!$OMP END PARALLEL DO
#else
    !$ACC END PARALLEL LOOP
#endif

#ifndef ENABLE_OPENACC
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(diff, entries, um, vm, vol, dx, dy, n, q, row, elem, elnodes, c1, c2, c3, c4, cx1, cx2, cx3, cx4, entries2)
!$OMP DO
#else
#if !defined(DISABLE_OPENACC_ATOMICS)
    !$ACC PARALLEL LOOP GANG VECTOR PRIVATE(elnodes, dx, dy, entries, entries2) DEFAULT(PRESENT)
#else
    !$ACC UPDATE SELF(rhs_a, rhs_m, rhs_ms, rhs_adiv, rhs_mdiv, rhs_msdiv, u_ice, v_ice, m_ice, a_ice, m_snow)
#endif
#endif
    do elem=1,myDim_elem2D          !assembling rhs over elements
        elnodes=elem2D_nodes(:,elem)

        ! if cavity element skip it
        if (ulevels(elem)>1) cycle

        !derivatives
        dx=gradient_sca(1:3,elem)
        dy=gradient_sca(4:6,elem)
        vol=elem_area(elem)
        um=sum(u_ice(elnodes))
        vm=sum(v_ice(elnodes))
        ! this is exact computation (no assumption of u=const on elements used
        ! in the standard version)
        c1=(um*um+sum(u_ice(elnodes)*u_ice(elnodes)))/12.0_WP
        c2=(vm*vm+sum(v_ice(elnodes)*v_ice(elnodes)))/12.0_WP
        c3=(um*vm+sum(v_ice(elnodes)*u_ice(elnodes)))/12.0_WP
        c4=sum(dx*u_ice(elnodes)+dy*v_ice(elnodes))
        do n=1,3
            row=elnodes(n)
                !!PS         if(ulevels_nod2D(row)>1) cycle !LK89140
            do q = 1,3
                entries(q)= vol*ice%ice_dt*((1.0_WP-0.5_WP*ice%ice_dt*c4)*(dx(n)*(um+u_ice(elnodes(q)))+ &
                            dy(n)*(vm+v_ice(elnodes(q))))/12.0_WP - &
                            0.5_WP*ice%ice_dt*(c1*dx(n)*dx(q)+c2*dy(n)*dy(q)+c3*(dx(n)*dy(q)+dx(q)*dy(n))))
                        !um*dx(n)+vm*dy(n))*(um*dx(q)+vm*dy(q))/9.0)
                entries2(q)=0.5_WP*ice%ice_dt*(dx(n)*(um+u_ice(elnodes(q)))+ &
                            dy(n)*(vm+v_ice(elnodes(q)))-dx(q)*(um+u_ice(row))- &
                            dy(q)*(vm+v_ice(row)))
            end do

            !___________________________________________________________________
            cx1=vol*ice%ice_dt*c4*(sum(m_ice(elnodes))+m_ice(elnodes(n))+sum(entries2*m_ice(elnodes)))/12.0_WP
            cx2=vol*ice%ice_dt*c4*(sum(a_ice(elnodes))+a_ice(elnodes(n))+sum(entries2*a_ice(elnodes)))/12.0_WP
            cx3=vol*ice%ice_dt*c4*(sum(m_snow(elnodes))+m_snow(elnodes(n))+sum(entries2*m_snow(elnodes)))/12.0_WP
#if defined (__oifs) || defined (__ifsinterface)
            cx4=vol*ice%ice_dt*c4*(sum(ice_temp(elnodes))+ice_temp(elnodes(n))+sum(entries2*ice_temp(elnodes)))/12.0_WP
#endif

            !___________________________________________________________________
#ifndef ENABLE_OPENACC
#if defined(_OPENMP)  && !defined(__openmp_reproducible)
                call omp_set_lock  (partit%plock(row))
#else
!$OMP ORDERED
#endif
#endif
            tmp_sum = sum(entries*m_ice(elnodes))
#if !defined(DISABLE_OPENACC_ATOMICS)
            !$ACC ATOMIC WRITE
#endif
            rhs_m(row)=rhs_m(row)+tmp_sum+cx1

            tmp_sum = sum(entries*a_ice(elnodes))
#if !defined(DISABLE_OPENACC_ATOMICS)
            !$ACC ATOMIC WRITE
#endif
            rhs_a(row)=rhs_a(row)+tmp_sum+cx2

            tmp_sum = sum(entries*m_snow(elnodes))
#if !defined(DISABLE_OPENACC_ATOMICS)
            !$ACC ATOMIC WRITE
#endif
            rhs_ms(row)=rhs_ms(row)+tmp_sum+cx3

#if defined (__oifs) || defined (__ifsinterface)
            tmp_sum = sum(entries*ice_temp(elnodes))
#if !defined(DISABLE_OPENACC_ATOMICS)
            !$ACC ATOMIC UPDATE
#endif
            rhs_temp(row)=rhs_temp(row)+tmp_sum+cx4
#endif

            !___________________________________________________________________
#if !defined(DISABLE_OPENACC_ATOMICS)
            !$ACC ATOMIC UPDATE
#endif
            rhs_mdiv(row)=rhs_mdiv(row)-cx1
#if !defined(DISABLE_OPENACC_ATOMICS)
            !$ACC ATOMIC UPDATE
#endif
            rhs_adiv(row)=rhs_adiv(row)-cx2
#if !defined(DISABLE_OPENACC_ATOMICS)
            !$ACC ATOMIC UPDATE
#endif
            rhs_msdiv(row)=rhs_msdiv(row)-cx3
#if defined (__oifs) || defined (__ifsinterface)
#if !defined(DISABLE_OPENACC_ATOMICS)
            !$ACC ATOMIC UPDATE
#endif
            rhs_tempdiv(row)=rhs_tempdiv(row)-cx4
#endif /* (__oifs) */
#ifndef ENABLE_OPENACC
#if defined(_OPENMP)  && !defined(__openmp_reproducible)
                call omp_unset_lock(partit%plock(row))
#else
!$OMP END ORDERED
#endif
#endif
        end do
    end do
#ifndef ENABLE_OPENACC
!$OMP END DO
!$OMP END PARALLEL
#else
#if !defined(DISABLE_OPENACC_ATOMICS)
    !$ACC END PARALLEL LOOP
#else
    !$ACC UPDATE DEVICE(rhs_a, rhs_m, rhs_ms, rhs_adiv, rhs_mdiv, rhs_msdiv)
#endif
#endif

    !! Juha: close associate blocks
    end associate
    end associate
    end associate

end subroutine ice_TG_rhs_div
!
!
!_______________________________________________________________________________
subroutine ice_update_for_div(ice, partit, mesh)
    use MOD_ICE
    USE MOD_PARTIT
    USE MOD_PARSUP
    USE MOD_MESH
    use o_PARAM
    use g_CONFIG
    use g_comm_auto
    implicit none
    type(t_ice)   , intent(inout), target   :: ice
    type(t_partit), intent(inout), target   :: partit
    type(t_mesh)  , intent(in)   , target   :: mesh
    !___________________________________________________________________________
    integer                                 :: n,clo,clo2,cn,location(100),row, tmp_it
    real(kind=WP)                           :: rhs_new
    integer                                 :: num_iter_solve=3
    !___________________________________________________________________________

! Juha : Use associate blocks instead of pointers to work around gpu memcopy issues with Cray
#include "associate_part_ass_combined.h"
#include "associate_mesh_ass_combined.h"

    associate( &
#if defined (__oifs) || defined (__ifsinterface)
    ice_temp     => ice%data(4)%values(:), &
    m_templ      => ice%data(4)%valuesl(:), &
    dm_temp      => ice%data(4)%dvalues(:), &
    rhs_tempdiv  => ice%data(4)%values_div_rhs(:), &
#endif
    a_ice        => ice%data(1)%values(:), &
    m_ice        => ice%data(2)%values(:), &
    m_snow       => ice%data(3)%values(:), &
    rhs_adiv     => ice%data(1)%values_div_rhs(:), &
    rhs_mdiv     => ice%data(2)%values_div_rhs(:), &
    rhs_msdiv    => ice%data(3)%values_div_rhs(:), &
    a_icel       => ice%data(1)%valuesl(:), &
    m_icel       => ice%data(2)%valuesl(:), &
    m_snowl      => ice%data(3)%valuesl(:), &
    da_ice       => ice%data(1)%dvalues(:), &
    dm_ice       => ice%data(2)%dvalues(:), &
    dm_snow      => ice%data(3)%dvalues(:), &
    mass_matrix  => ice%work%fct_massmatrix(:) &
    )


    !___________________________________________________________________________
    ! Does Taylor-Galerkin solution
    ! the first approximation
#ifndef ENABLE_OPENACC
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(row)
#else
    !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT)
#endif
    do row=1,myDim_nod2D
        !! row=myList_nod2D(m)
        ! if cavity node skip it
        if (ulevels_nod2d(row)>1) cycle

        dm_ice(row) =rhs_mdiv(row) /area(1,row)
        da_ice(row) =rhs_adiv(row) /area(1,row)
        dm_snow(row)=rhs_msdiv(row)/area(1,row)
#if defined (__oifs) || defined (__ifsinterface)
        dm_temp(row)=rhs_tempdiv(row)/area(1,row)
#endif
    end do
#ifndef ENABLE_OPENACC
!$OMP END PARALLEL DO
#else
    !$ACC END PARALLEL LOOP
#endif
    call exchange_nod(dm_ice, partit, luse_g2g = .true.)
    call exchange_nod(da_ice, partit, luse_g2g = .true.)
    call exchange_nod(dm_snow, partit, luse_g2g = .true.)
#if defined (__oifs) || defined (__ifsinterface)
    call exchange_nod(dm_temp, partit, luse_g2g = .true.)
#endif /* (__oifs) */
#ifndef ENABLE_OPENACC
!$OMP BARRIER
#endif
    !___________________________________________________________________________
    !iterate
    do n=1,num_iter_solve-1
#ifndef ENABLE_OPENACC
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(row, n, clo, clo2, cn, location, rhs_new)
!$OMP DO
#else
        !$ACC PARALLEL LOOP GANG VECTOR PRIVATE(location) DEFAULT(PRESENT)
#endif
        do row=1,myDim_nod2D
            ! if cavity node skip it
            if (ulevels_nod2d(row)>1) cycle
                  !! row=myList_nod2D(m)
            !___________________________________________________________________
            clo=ssh_stiff%rowptr(row)-ssh_stiff%rowptr(1)+1
            clo2=ssh_stiff%rowptr(row+1)-ssh_stiff%rowptr(1)
            cn=clo2-clo+1
            do tmp_it = 1, cn
                location(tmp_it)=nn_pos(tmp_it, row)
            end do

            !___________________________________________________________________
            rhs_new     = rhs_mdiv(row) - sum(mass_matrix(clo:clo2)*dm_ice(location(1:cn)))
            m_icel(row) = dm_ice(row)+rhs_new/area(1,row)
            rhs_new     = rhs_adiv(row) - sum(mass_matrix(clo:clo2)*da_ice(location(1:cn)))
            a_icel(row) = da_ice(row)+rhs_new/area(1,row)
            rhs_new     = rhs_msdiv(row) - sum(mass_matrix(clo:clo2)*dm_snow(location(1:cn)))
            m_snowl(row)= dm_snow(row)+rhs_new/area(1,row)
#if defined (__oifs) || defined (__ifsinterface)
            rhs_new     = rhs_tempdiv(row) - sum(mass_matrix(clo:clo2)*dm_temp(location(1:cn)))
            m_templ(row)= dm_temp(row)+rhs_new/area(1,row)
#endif
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
        do row=1,myDim_nod2D
            ! if cavity node skip it
            if (ulevels_nod2d(row)>1) cycle
            dm_ice(row)  = m_icel(row)
            da_ice(row)  = a_icel(row)
            dm_snow(row) = m_snowl(row)
#if defined (__oifs) || defined (__ifsinterface)
            dm_temp(row) = m_templ(row)
#endif
        end do
#ifndef ENABLE_OPENACC
!$OMP END DO
!$OMP END PARALLEL
#else
        !$ACC END PARALLEL LOOP
#endif
        call exchange_nod(dm_ice, partit, luse_g2g = .true.)
        call exchange_nod(da_ice, partit, luse_g2g = .true.)
        call exchange_nod(dm_snow, partit, luse_g2g = .true.)
#if defined (__oifs) || defined (__ifsinterface)
        call exchange_nod(dm_temp, partit, luse_g2g = .true.)
#endif /* (__oifs) */
#ifndef ENABLE_OPENACC
!$OMP BARRIER
#endif
    end do

#ifndef ENABLE_OPENACC
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(row)
#else
    !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT)
#endif
    do row=1, myDim_nod2D+eDim_nod2D
       m_ice(row)   = m_ice (row)+dm_ice (row)
       a_ice(row)   = a_ice (row)+da_ice (row)
       m_snow(row)  = m_snow(row)+dm_snow(row)
#if defined (__oifs) || defined (__ifsinterface)
       ice_temp(row)= ice_temp(row)+dm_temp(row)
#endif
    end do
#ifndef ENABLE_OPENACC
!$OMP END PARALLEL DO
#else
    !$ACC END PARALLEL LOOP
#endif

    ! Juha: close associate blocks
    end associate
    end associate
    end associate

end subroutine ice_update_for_div
! =============================================================
