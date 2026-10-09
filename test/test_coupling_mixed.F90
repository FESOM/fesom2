program test_coupling_mixed
    ! Forcing and ice-ocean coupling on mixed triangle/quad meshes (plan
    ! 2026-10-08-hybrid-meshes, Task 13): the node->element surface stress of oce_fluxes_mom
    ! (mod_ice_oce_coupling) and forcing_stress_surf (mod_forcing_bulk) is the kite-weighted
    ! mean sum(elem_w*stress_node_surf) over the 4 vertex slots, and the oce_fluxes balancing
    ! uses the mixed-mesh node areas. Mesh from FESOM3_MESH_DIR (default the NEMO-mixed
    ! cut-out, unrotated; FESOM3_TEST_ROTATE=1 for pi/CORE2) at 1 or N ranks and, at 1 rank,
    ! the analytic mixed patch. Node fields are hashes of the global node id (identical on
    ! owned and halo copies), so the identities below are exact to roundoff:
    !   C1 a uniform node stress gives the same stress on every owned element (sum elem_w
    !      = 1; the former /3 over the 4 slots of a quad gives 4/3)
    !   C2 momentum input is the same on nodes and elements: sum_e elem_area*stress_surf ==
    !      sum_n area*stress_node_surf (global; each element counted once, on the owner of
    !      its first vertex), with ice (a_ice in [0,1], ice-ocean drag) — holds because
    !      area(n) = sum_e elem_area*elem_w(n)
    !   C3 triangles: the kite mean equals FESOM2's sum/3 to roundoff (CORE2/pi regression)
    !   C4 (1 rank) forcing_stress_surf == oce_fluxes_mom at a_ice = 0, bitwise
    !   C5 oce_fluxes balancing: integral of virtual_salt and of relax_salt == 0, and the
    !      freshwater balance integral(water_flux) == integral(-flx_fw + flux) with
    !      ocean_area == sum of the owned areasvol
    ! A mutation check prints the error of the former /3 element mean on the quads (C1 must
    ! fail with it on a mesh that has quads).
    use mpi
    use test_util,         only: setup_test_mesh, check, finish, max_all, sum_all, hval
    use mod_precision,        only: WP
    use mod_mesh,             only: t_mesh
    use mod_partit,           only: t_partit
    use mod_partitioning,     only: par_init
    use mod_mesh_analytic,    only: generate_mixed_patch
    use mod_part_bounds,      only: owned_bounds
    use mod_ice,              only: t_ice
    use mod_ice_thermo,       only: t_atmflux
    use mod_tracer,           only: t_tracer
    use mod_ice_oce_coupling, only: oce_fluxes_mom, oce_fluxes
    use mod_forcing_bulk,     only: forcing_stress_surf
    implicit none

    type(t_partit) :: partit, apartit
    type(t_mesh)   :: mesh, amesh

    call par_init(partit)
    call setup_test_mesh(mesh, partit, '/work/ab0995/a270029/meshes/nemo025mix_cut', cartesian=.false., partition_np1=.false.)
    call check_coupling(mesh, partit, 'mesh')

    if (partit%npes == 1) then
        call generate_mixed_patch(amesh, apartit, nx=13, ny=11, nl=6, &
                                  Lx=12000.0_WP, Ly=10000.0_WP, max_depth=100.0_WP)
        call check_coupling(amesh, apartit, 'patch')
    end if

    call finish(partit, 'test_coupling_mixed')

contains

    subroutine check_coupling(m, p, tag)
        type(t_mesh),     intent(inout) :: m
        type(t_partit),   intent(inout) :: p
        character(len=*), intent(in)    :: tag
        type(t_ice)     :: ice
        type(t_atmflux) :: atm
        type(t_tracer)  :: tracers
        integer :: nNodO, nNodL, nEdgeO, nElemO, n, e, gid, ntri, nquad, k
        real(kind=WP), allocatable :: ss(:,:), ss2(:,:), fin(:)
        real(kind=WP) :: err(6), mut, iel(2), inod(2), sc(2), v(8), w(8), smax
        real(kind=WP), parameter :: SX = 0.1_WP, SY = -0.05_WP

        call owned_bounds(m, nNodO, nNodL, nEdgeO, nElemO, p)
        ntri  = count(m%elem2D_nnodes(1:nElemO) == 3)
        nquad = count(m%elem2D_nnodes(1:nElemO) == 4)
        call alloc_fields(ice, atm, tracers, nNodL, m%nl)
        allocate(ss(2, m%elem2D), ss2(2, m%elem2D), fin(nNodL))

        !--- C1: uniform atm-ocean stress, no ice ---
        ice%data(1)%values = 0.0_WP
        atm%stress_atmoce_x = SX; atm%stress_atmoce_y = SY
        ss = 0.0_WP
        call oce_fluxes_mom(ice, atm, ss, m, p)
        err(1) = 0.0_WP; mut = 0.0_WP
        do e = 1, nElemO
            err(1) = max(err(1), abs(ss(1,e) - SX)/abs(SX), abs(ss(2,e) - SY)/abs(SY))
            if (m%elem2D_nnodes(e) == 4) &        ! the former element mean, for the record
                mut = max(mut, abs(sum(atm%stress_node_surf(1, m%elem2D_nodes(1:4,e)))/3.0_WP - SX)/abs(SX))
        end do
        call max_all(err(1), p); call max_all(mut, p)

        !--- C2/C3: non-uniform stress, partial ice cover with ice-ocean drag ---
        do n = 1, nNodL
            gid = n
            if (p%npes > 1) gid = p%myList_nod2D(n)
            ice%data(1)%values(n) = 0.5_WP + 0.5_WP*hval(gid, 1)
            ice%uice(n)     = 0.2_WP*hval(gid, 2); ice%vice(n)     = 0.2_WP*hval(gid, 3)
            ice%srfoce_u(n) = 0.1_WP*hval(gid, 4); ice%srfoce_v(n) = 0.1_WP*hval(gid, 5)
            atm%stress_atmoce_x(n) = 0.1_WP*hval(gid, 6)
            atm%stress_atmoce_y(n) = 0.1_WP*hval(gid, 7)
        end do
        ss = 0.0_WP
        call oce_fluxes_mom(ice, atm, ss, m, p)
        iel = 0.0_WP; inod = 0.0_WP; sc = 0.0_WP
        do e = 1, nElemO
            ! owned elements are shared (every element with an owned vertex): count each
            ! once, on the rank that owns its first vertex
            if (m%elem2D_nodes(1,e) > nNodO) cycle
            iel = iel + real(m%elem_area(e), WP)*ss(:,e)
        end do
        do n = 1, nNodO
            inod = inod + real(m%area(n), WP)*atm%stress_node_surf(:,n)
            sc = sc + real(m%area(n), WP)*abs(atm%stress_node_surf(:,n))
        end do
        v(1:2) = iel; v(3:4) = inod; v(5:6) = sc
        call sum_all(v(1:6), p)
        err(2) = max(abs(v(1) - v(3))/v(5), abs(v(2) - v(4))/v(6))
        err(3) = 0.0_WP; smax = maxval(abs(atm%stress_node_surf(:, 1:nNodL)))
        do e = 1, nElemO
            if (m%elem2D_nnodes(e) /= 3) cycle
            do k = 1, 2
                err(3) = max(err(3), abs(ss(k,e) - sum(atm%stress_node_surf(k, m%elem2D_nodes(1:3,e)))/3.0_WP))
            end do
        end do
        call max_all(err(3), p); call max_all(smax, p)
        err(3) = err(3)/smax

        !--- C4: forcing_stress_surf (1-rank driver path) == oce_fluxes_mom at a_ice = 0 ---
        err(4) = 0.0_WP
        if (p%npes == 1) then
            ice%data(1)%values = 0.0_WP
            ss = 0.0_WP; ss2 = 0.0_WP
            call oce_fluxes_mom(ice, atm, ss, m, p)
            call forcing_stress_surf(atm%stress_atmoce_x, atm%stress_atmoce_y, ss2, m)
            err(4) = maxval(abs(ss(:, 1:m%elem2D) - ss2(:, 1:m%elem2D)))
        end if

        !--- C5: oce_fluxes balancing on the mixed-mesh node areas ---
        do n = 1, nNodL
            gid = n
            if (p%npes > 1) gid = p%myList_nod2D(n)
            ice%data(1)%values_old(n) = 0.5_WP + 0.5_WP*hval(gid, 8)
            ice%flx_fw(n)  = 1.0e-7_WP*hval(gid, 9)
            ice%flx_h(n)   = 100.0_WP*hval(gid, 10)
            tracers%data(2)%values(:, n) = 34.0_WP + hval(gid, 11)
            atm%Ssurf(n)          = 34.5_WP + 0.5_WP*hval(gid, 12)
            atm%evaporation(n)    = -2.0e-8_WP*(1.0_WP + hval(gid, 13))
            atm%ice_sublimation(n) = 1.0e-9_WP*hval(gid, 14)
            atm%prec_rain(n)      = 3.0e-8_WP*(1.0_WP + hval(gid, 15))
            atm%prec_snow(n)      = 1.0e-8_WP*(1.0_WP + hval(gid, 16))
            atm%runoff(n)         = 5.0e-9_WP*(1.0_WP + hval(gid, 17))
            fin(n) = -ice%flx_fw(n) + atm%evaporation(n) - atm%ice_sublimation(n) + atm%prec_rain(n) &
                     + atm%prec_snow(n)*(1.0_WP - ice%data(1)%values_old(n)) + atm%runoff(n)
        end do
        atm%use_virt_salt = .true.; atm%ref_sss_local = .true.
        call oce_fluxes(ice, tracers, atm, m, p)
        w = 0.0_WP
        do n = 1, nNodO
            w(1) = w(1) + real(m%areasvol(n), WP)*atm%virtual_salt(n)
            w(2) = w(2) + real(m%areasvol(n), WP)*abs(atm%virtual_salt(n))
            w(3) = w(3) + real(m%areasvol(n), WP)*atm%relax_salt(n)
            w(4) = w(4) + real(m%areasvol(n), WP)*abs(atm%relax_salt(n))
            w(5) = w(5) + real(m%areasvol(n), WP)*(atm%water_flux(n) - fin(n))
            w(6) = w(6) + real(m%areasvol(n), WP)*abs(fin(n))
            w(7) = w(7) + real(m%areasvol(n), WP)
        end do
        call sum_all(w(1:7), p)
        err(5) = max(abs(w(1))/w(2), abs(w(3))/w(4), abs(w(5))/w(6))
        err(6) = abs(w(7) - real(m%ocean_area, WP))/w(7)

        if (p%mype == 0) write(*,'(a,a,i0,a,i0,a,6es10.2,a,es10.2)') tag, ': owned triangles=', ntri, &
            ' quads=', nquad, '  C1 C2 C3 C4 C5 ocean_area rel. errors:', err, &
            '  (former /3 mean on quads: C1 error', mut, ')'
        call check(err(1) < 1.0e-14_WP, tag//': C1 uniform node stress -> same element stress', p%mype)
        call check(err(2) < 1.0e-13_WP, tag//': C2 momentum input sum_e A*stress_surf == sum_n area*stress', p%mype)
        call check(err(3) < 1.0e-15_WP, tag//': C3 triangle kite mean == sum/3 (roundoff)', p%mype)
        call check(err(4) == 0.0_WP,    tag//': C4 forcing_stress_surf == oce_fluxes_mom at a_ice=0', p%mype)
        call check(err(5) < 1.0e-13_WP, tag//': C5 virtual_salt / relax_salt / freshwater balanced', p%mype)
        call check(err(6) < 1.0e-13_WP, tag//': C5 ocean_area == sum of the owned areasvol', p%mype)
        if (nquad > 0) call check(mut > 0.3_WP, tag//': mutation: the former /3 mean fails C1 on quads', p%mype)
        deallocate(ss, ss2, fin)
    end subroutine check_coupling

    subroutine alloc_fields(ice, atm, tracers, nn, nl)
        type(t_ice),     intent(inout) :: ice
        type(t_atmflux), intent(inout) :: atm
        type(t_tracer),  intent(inout) :: tracers
        integer,         intent(in)    :: nn, nl
        allocate(ice%uice(nn), ice%vice(nn), ice%srfoce_u(nn), ice%srfoce_v(nn))
        allocate(ice%stress_iceoce_x(nn), ice%stress_iceoce_y(nn), ice%flx_fw(nn), ice%flx_h(nn))
        allocate(ice%data(ice%num_itracers))
        allocate(ice%data(1)%values(nn), ice%data(1)%values_old(nn))
        allocate(ice%thermo%thdgr(nn), ice%thermo%thdgrsn(nn))
        allocate(atm%stress_atmoce_x(nn), atm%stress_atmoce_y(nn), atm%stress_node_surf(2, nn))
        allocate(atm%heat_flux(nn), atm%water_flux(nn), atm%heat_flux_in(nn))
        allocate(atm%virtual_salt(nn), atm%relax_salt(nn), atm%Ssurf(nn))
        allocate(atm%evaporation(nn), atm%ice_sublimation(nn))
        allocate(atm%prec_rain(nn), atm%prec_snow(nn), atm%runoff(nn))
        allocate(tracers%data(2), tracers%data(2)%values(nl-1, nn))
        ice%uice = 0.0_WP; ice%vice = 0.0_WP; ice%srfoce_u = 0.0_WP; ice%srfoce_v = 0.0_WP
        ice%stress_iceoce_x = 0.0_WP; ice%stress_iceoce_y = 0.0_WP
        ice%flx_fw = 0.0_WP; ice%flx_h = 0.0_WP
        ice%data(1)%values = 0.0_WP; ice%data(1)%values_old = 0.0_WP
        ice%thermo%thdgr = 0.0_WP; ice%thermo%thdgrsn = 0.0_WP
        atm%stress_node_surf = 0.0_WP; atm%heat_flux = 0.0_WP; atm%water_flux = 0.0_WP
        atm%heat_flux_in = 0.0_WP; atm%virtual_salt = 0.0_WP; atm%relax_salt = 0.0_WP
        tracers%data(2)%values = 0.0_WP
    end subroutine alloc_fields

end program test_coupling_mixed
