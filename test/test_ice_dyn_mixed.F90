program test_ice_dyn_mixed
    ! Sea-ice EVP / mEVP dynamics on mixed triangle/quad meshes (plan
    ! 2026-10-08-hybrid-meshes, Task 15): strain rates over the 4 slots of gradient_sca with
    ! kite-weighted metric means, the stress divergence (stress2rhs and the fused mEVP loop)
    ! with node shares elem_area*elem_w(k), the SSH-tilt term with the same shares, msum/asum
    ! via elem_w and the kite-share-weighted ocean-velocity node gather (ocean2ice). Mesh from
    ! FESOM3_MESH_DIR (default the NEMO-mixed cut-out, unrotated; FESOM3_TEST_ROTATE=1 for
    ! pi/CORE2) at 1 or N ranks and, at 1 rank, the analytic triangle, mixed and
    ! parallelogram patches (planar, metric 0). Node/element fields are hashes of the global
    ! id (identical on owned and halo copies):
    !   D1 strain rates (stress_tensor): planar patches -- uniform drift and rigid rotation
    !      give zero strain, a linear deformation its exact strain; sphere -- rigid rotation
    !      about the mesh pole (u = Omega*r*cos(lat), v = 0) gives zero strain up to the
    !      discretisation error (the metric terms cancel du/dy: < 0.1 of the no-metric value,
    !      both reported; measured 3.6e-4/0.10 NEMO-mixed cut, 3.3e-2/0.50 pi near the
    !      rotated pole); triangles == the former P1 formulas (1:3, /3)
    !   D2 stress2rhs is the exact negative adjoint of the strain rates (global, any mesh):
    !      sum_n (u*Ru + v*Rv) == -sum_e elem_area*(s11*e11 + s22*e22 + 2*s12*e12); the
    !      former /3 metric share breaks it on quads. Triangles: stress2rhs == the former
    !      formula at all-triangle nodes
    !   D3 SSH tilt (EVP and mEVP): sum_n area*rhs_a == -g*sum_e elem_area*d(eta)/dx (global,
    !      any mesh); planar patches with a linear eta: rhs_a == -g*a at every node
    !   D4 uniform strain is force-free (planar patches, EVP and mEVP, one subcycle): uniform
    !      ice with a linear velocity drifting with the ocean has a uniform stress, which the
    !      4-slot scatter must cancel around every interior node -> Coriolis-only update there
    !   D5 ocean2ice: a uniform surface velocity gives the same node velocity (owned + halo);
    !      triangles: the kite-share gather == the former area-weighted mean at all-triangle
    !      nodes
    !   D6 (np > 1) shared elements (docs/plans/completed/20261012-shared-element-sync.md): a border
    !      element is owned by every rank that owns one of its vertices and only the ocean uv
    !      is synced between the copies; the sigma copies must stay bitwise equal by
    !      themselves after EVP and after mEVP (3 subcycles, hashed ice velocity, a, m and
    !      SSH), and sigma must be nonzero on some shared element (not vacuous)
    use mpi
    use test_util,         only: setup_test_mesh, check, finish, max_all, sum_all, hval
    use mod_precision,     only: WP
    use mod_constants,     only: r_earth, g
    use mod_mesh,          only: t_mesh, MAX_NV
    use mod_partit,        only: t_partit
    use mod_partitioning,  only: par_init
    use mod_mesh_analytic, only: generate_analytic_mesh, generate_mixed_patch, generate_quad_patch
    use mod_part_bounds,   only: owned_bounds
    use mod_ice,           only: t_ice
    use mod_dyn,           only: t_dyn
    use mod_tracer,        only: t_tracer
    use mod_ice_dyn,       only: ocean2ice, EVPdynamics_solve, stress_tensor, stress2rhs
    use mod_ice_setup,     only: build_bc_index_nod2D
    use mod_halo,          only: shared_copy_max_elem
    implicit none

    type(t_partit) :: partit, tpartit, apartit, qpartit
    type(t_mesh)   :: mesh, tmesh, amesh, qmesh

    call par_init(partit)
    call setup_test_mesh(mesh, partit, '/work/ab0995/a270029/meshes/nemo025mix_cut', cartesian=.false., partition_np1=.false.)
    call check_ice(mesh, partit, 'mesh', .false.)

    if (partit%npes == 1) then
        call generate_analytic_mesh(tmesh, tpartit, nx=13, ny=11, nl=4, &
                                    Lx=12000.0_WP, Ly=10000.0_WP, max_depth=100.0_WP)
        call check_ice(tmesh, tpartit, 'patch_tri', .true.)
        call generate_mixed_patch(amesh, apartit, nx=13, ny=11, nl=4, &
                                  Lx=12000.0_WP, Ly=10000.0_WP, max_depth=100.0_WP)
        call check_ice(amesh, apartit, 'patch_mixed', .true.)
        call generate_quad_patch(qmesh, qpartit, nx=13, ny=11, nl=4, &
                                 Lx=12000.0_WP, Ly=10000.0_WP, max_depth=100.0_WP, skew=0.3_WP)
        call check_ice(qmesh, qpartit, 'patch_quad', .true.)
    end if

    call finish(partit, 'test_ice_dyn_mixed')

contains

    subroutine check_ice(m, p, tag, planar)
        type(t_mesh),     intent(inout) :: m
        type(t_partit),   intent(inout) :: p
        character(len=*), intent(in)    :: tag
        logical,          intent(in)    :: planar
        type(t_ice)    :: ice
        type(t_dyn)    :: dyn
        type(t_tracer) :: tracers
        integer :: nNodO, nNodL, nEdgeO, nElemO, ne, n, e, k, nd, ntri, nquad, iv
        integer, allocatable :: ngid(:), egid(:)
        logical, allocatable :: alltri(:)
        real(kind=WP), allocatable :: X(:), Y(:), ro(:), rabs(:)
        real(kind=WP) :: d1(4), d1s(2), d1t, d2, d2t, d3(2), d4(2), d4s(2), d5(2), v(4), sc, sc2
        real(kind=WP) :: d6(2), d6s(2)
        real(kind=WP) :: e11, e22, e12, ex, wa, uo, vo, rb, cb, ue, ve
        real(kind=WP), allocatable :: u0(:), v0(:)
        real(kind=WP), parameter :: OM = 1.0e-6_WP, LA = 1.0e-6_WP, LB = 2.0e-6_WP, &
                                    LC = -0.5e-6_WP, LD = 0.7e-6_WP, EA = 1.0e-6_WP, EB = -2.0e-6_WP

        call owned_bounds(m, nNodO, nNodL, nEdgeO, nElemO, p)
        ne = size(m%elem2D_nodes, 2)
        ntri  = count(m%elem2D_nnodes(1:nElemO) == 3)
        nquad = count(m%elem2D_nnodes(1:nElemO) == 4)
        allocate(ngid(nNodL), egid(ne), X(nNodL), Y(nNodL), ro(nNodL), rabs(nNodL), alltri(nNodO))
        do n = 1, nNodL
            ngid(n) = n
            if (p%npes > 1) ngid(n) = p%myList_nod2D(n)
            X(n) = real(m%coord_nod2D(1, n), WP)*r_earth   ! planar patches: metres
            Y(n) = real(m%coord_nod2D(2, n), WP)*r_earth
        end do
        do e = 1, ne
            egid(e) = e
            if (p%npes > 1) egid(e) = p%myList_elem2D(e)
        end do
        do n = 1, nNodO
            alltri(n) = all(m%elem2D_nnodes(m%nod_in_elem2D(1:m%nod_in_elem2D_num(n), n)) == 3)
        end do
        call alloc_ice(ice, m, p, nNodL, ne)

        !--- D1: strain rates ---
        ice%work%ice_strength = 0.0_WP
        ice%work%ice_strength(1:nElemO) = 1.0_WP
        d1 = 0.0_WP; d1s = 0.0_WP; d1t = 0.0_WP
        if (planar) then
            sc = 0.0_WP
            do e = 1, nElemO
                sc = max(sc, sum(abs(real(m%gradient_sca(:, e), WP))))
            end do
            ! uniform drift, relative to the size of the single terms
            ice%uice = 0.3_WP; ice%vice = -0.2_WP
            call strain(ice, m, p)
            d1(1) = max_eps(ice, nElemO)/(0.3_WP*sc)
            ! rigid rotation about the patch centre
            do n = 1, nNodL
                ice%uice(n) = -OM*(Y(n) - 5000.0_WP)
                ice%vice(n) =  OM*(X(n) - 6000.0_WP)
            end do
            call strain(ice, m, p)
            d1(2) = max_eps(ice, nElemO)/OM
            ! linear deformation: exact strain
            do n = 1, nNodL
                ice%uice(n) = LA*X(n) + LB*Y(n)
                ice%vice(n) = LC*X(n) + LD*Y(n)
            end do
            call strain(ice, m, p)
            do e = 1, nElemO
                d1(3) = max(d1(3), abs(ice%work%eps11(e) - LA), abs(ice%work%eps22(e) - LD), &
                        abs(ice%work%eps12(e) - 0.5_WP*(LB + LC)))
            end do
            d1(3) = d1(3)/LB
        else
            ! rigid rotation about the mesh pole: zero strain up to discretisation error
            do n = 1, nNodL
                ice%uice(n) = OM*r_earth*cos(real(m%coord_nod2D(2, n), WP))
                ice%vice(n) = 0.0_WP
            end do
            call strain(ice, m, p)
            d1s(1) = max_eps(ice, nElemO)/OM
            do e = 1, nElemO               ! the same without the metric terms, for scale
                d1s(2) = max(d1s(2), 0.5_WP*abs(sum(real(m%gradient_sca(5:8, e), WP) &
                                                    *ice%uice(m%elem2D_nodes(:, e)))))
            end do
            d1s(2) = d1s(2)/OM
            call max_all(d1s(1), p); call max_all(d1s(2), p)
        end if
        ! triangles == the former P1 formulas (hashed velocities)
        call hashed_velocity(ice, ngid, nNodL)
        call strain(ice, m, p)
        sc = 0.0_WP
        do e = 1, nElemO
            if (m%elem2D_nnodes(e) /= 3) cycle
            call former_strain(m, ice, e, e11, e22, e12)
            d1t = max(d1t, abs(ice%work%eps11(e) - e11), abs(ice%work%eps22(e) - e22), &
                      abs(ice%work%eps12(e) - e12))
            sc = max(sc, abs(e11), abs(e22), abs(e12))
        end do
        call max_all(d1t, p); call max_all(sc, p)
        if (sc > 0.0_WP) d1t = d1t/sc
        call max_all(d1(1), p); call max_all(d1(2), p); call max_all(d1(3), p)

        !--- D2: stress2rhs = - adjoint of the strain rates ---
        ! eps of the hashed velocities is in place; hashed element stresses; raw rhs
        do e = 1, ne
            ice%work%sigma11(e) = 1.0e4_WP*hval(egid(e), 11)
            ice%work%sigma22(e) = 1.0e4_WP*hval(egid(e), 12)
            ice%work%sigma12(e) = 1.0e4_WP*hval(egid(e), 13)
        end do
        ice%work%inv_areamass = 1.0_WP
        ice%data(1)%values_rhs = 0.0_WP; ice%data(2)%values_rhs = 0.0_WP
        call stress2rhs(ice, m, p)
        v = 0.0_WP
        do n = 1, nNodO
            v(1) = v(1) + ice%uice(n)*ice%uice_rhs(n) + ice%vice(n)*ice%vice_rhs(n)
            v(2) = v(2) + abs(ice%uice(n)*ice%uice_rhs(n)) + abs(ice%vice(n)*ice%vice_rhs(n))
        end do
        do e = 1, nElemO
            if (m%elem2D_nodes(1, e) > nNodO) cycle   ! each element once: owner of its 1st vertex
            wa = real(m%elem_area(e), WP)
            v(3) = v(3) + wa*(ice%work%sigma11(e)*ice%work%eps11(e) + ice%work%sigma22(e)*ice%work%eps22(e) &
                              + 2.0_WP*ice%work%sigma12(e)*ice%work%eps12(e))
            v(4) = v(4) + wa*(abs(ice%work%sigma11(e)*ice%work%eps11(e)) + abs(ice%work%sigma22(e)*ice%work%eps22(e)) &
                              + 2.0_WP*abs(ice%work%sigma12(e)*ice%work%eps12(e)))
        end do
        call sum_all(v, p)
        d2 = abs(v(1) + v(3))/max(v(2), v(4))
        ! triangles: the former stress2rhs (k = 1..3, metric share elem_area/3) at all-triangle nodes
        ro = 0.0_WP; rabs = 0.0_WP; d2t = 0.0_WP; sc2 = 0.0_WP
        do e = 1, nElemO
            if (m%elem2D_nnodes(e) /= 3) cycle
            do k = 1, 3
                nd = m%elem2D_nodes(k, e)
                ro(nd) = ro(nd) - real(m%elem_area(e), WP)*(ice%work%sigma11(e)*real(m%gradient_sca(k, e), WP) &
                         + ice%work%sigma12(e)*real(m%gradient_sca(k+4, e), WP) &
                         + ice%work%sigma12(e)*(1.0_WP/3.0_WP)*real(m%metric_factor(e), WP))
                rabs(nd) = rabs(nd) + real(m%elem_area(e), WP)*(abs(ice%work%sigma11(e)*real(m%gradient_sca(k, e), WP)) &
                         + abs(ice%work%sigma12(e)*real(m%gradient_sca(k+4, e), WP)) &
                         + abs(ice%work%sigma12(e)*real(m%metric_factor(e), WP))/3.0_WP)
            end do
        end do
        do n = 1, nNodO
            if (.not. alltri(n)) cycle
            d2t = max(d2t, abs(ice%uice_rhs(n) - ro(n)))
            sc2 = max(sc2, rabs(n))
        end do
        call max_all(d2t, p); call max_all(sc2, p)
        if (sc2 > 0.0_WP) d2t = d2t/sc2

        !--- D3: SSH tilt, EVP (whichEVP 0) and mEVP (1), 2 subcycles ---
        d3 = 0.0_WP
        do iv = 0, 1
            call ice_state_uniform(ice, nNodL)
            do n = 1, nNodL
                if (planar) then
                    ice%srfoce_ssh(n) = EA*X(n) + EB*Y(n)
                else
                    ice%srfoce_ssh(n) = 0.1_WP*hval(ngid(n), 21)
                end if
            end do
            ice%whichEVP = iv; ice%evp_rheol_steps = 2
            call EVPdynamics_solve(ice, m, p)
            v = 0.0_WP
            do n = 1, nNodO
                v(1) = v(1) + real(m%area(n), WP)*ice%data(1)%values_rhs(n)
                v(2) = v(2) + real(m%area(n), WP)*abs(ice%data(1)%values_rhs(n))
            end do
            do e = 1, nElemO
                if (m%elem2D_nodes(1, e) > nNodO) cycle
                ex = sum(real(m%gradient_sca(1:4, e), WP)*ice%srfoce_ssh(m%elem2D_nodes(:, e)))
                v(3) = v(3) - g*real(m%elem_area(e), WP)*ex
            end do
            call sum_all(v, p)
            d3(iv+1) = abs(v(1) - v(3))/v(2)
            if (planar) then
                do n = 1, nNodO
                    d3(iv+1) = max(d3(iv+1), abs(ice%data(1)%values_rhs(n) + g*EA)/abs(g*EA), &
                                             abs(ice%data(2)%values_rhs(n) + g*EB)/abs(g*EB))
                end do
            end if
            call max_all(d3(iv+1), p)
        end do

        !--- D4: uniform strain is force-free at interior nodes (planar patches) ---
        ! uniform ice with a linear velocity field (uniform strain -> uniform stress, also in
        ! mEVP, whose stress vanishes at rest), drifting with the ocean (no drag), flat sea
        ! surface, no wind: after ONE subcycle every interior node has the Coriolis-only
        ! update of its velocity. (A rest state is no test: mEVP has zero stress at rest,
        ! and the delta_min-regularised standard EVP amplifies roundoff at rest -- O(1e-2)
        ! m/s after 120 subcycles on the plain triangle patch, independent of this change.)
        d4 = 0.0_WP; d4s = 0.0_WP
        if (planar) then
            allocate(u0(nNodL), v0(nNodL))
            do iv = 0, 1
                call ice_state_uniform(ice, nNodL)
                ice%srfoce_ssh = 0.0_WP
                do n = 1, nNodL
                    u0(n) = 1.0e-5_WP*(X(n) - 6000.0_WP) + 2.0e-5_WP*(Y(n) - 5000.0_WP)
                    v0(n) = -0.5e-5_WP*(X(n) - 6000.0_WP) + 0.7e-5_WP*(Y(n) - 5000.0_WP)
                end do
                ice%uice = u0; ice%vice = v0; ice%srfoce_u = u0; ice%srfoce_v = v0
                ice%whichEVP = iv; ice%evp_rheol_steps = 1
                call EVPdynamics_solve(ice, m, p)
                d4s(iv+1) = maxval(abs(ice%work%sigma11(1:nElemO)))
                do n = 1, nNodO
                    if (ice%bc_index_nod2D(n) < 0.5_WP) cycle      ! boundary node (velocity 0)
                    rb = ice%ice_dt*real(m%coriolis_node(n), WP)
                    if (iv == 0) then
                        cb = 1.0_WP
                    else
                        cb = 1.0_WP + ice%beta_evp
                    end if
                    ue = (cb*cb*u0(n) + rb*cb*v0(n))/(cb*cb + rb*rb)
                    ve = (cb*cb*v0(n) - rb*cb*u0(n))/(cb*cb + rb*rb)
                    d4(iv+1) = max(d4(iv+1), abs(ice%uice(n) - ue), abs(ice%vice(n) - ve))
                end do
                d4(iv+1) = d4(iv+1)/maxval(abs(u0))
            end do
            deallocate(u0, v0)
        end if

        !--- D5: ocean2ice node gather ---
        allocate(dyn%uv(2, m%nl-1, ne))
        allocate(tracers%data(2))
        allocate(tracers%data(1)%values(m%nl-1, nNodL), tracers%data(2)%values(m%nl-1, nNodL))
        tracers%data(1)%values = 1.0_WP; tracers%data(2)%values = 35.0_WP
        if (.not. allocated(m%hbar)) then
            allocate(m%hbar(nNodL)); m%hbar = 0.0_WP
        end if
        dyn%uv(1, :, :) = 0.1_WP; dyn%uv(2, :, :) = -0.05_WP
        call ocean2ice(ice, dyn, tracers, m, p)
        d5 = 0.0_WP
        do n = 1, nNodL
            d5(1) = max(d5(1), abs(ice%srfoce_u(n) - 0.1_WP)/0.1_WP, abs(ice%srfoce_v(n) + 0.05_WP)/0.05_WP)
        end do
        do e = 1, ne
            dyn%uv(1, 1, e) = 0.1_WP*hval(egid(e), 31)
            dyn%uv(2, 1, e) = 0.1_WP*hval(egid(e), 32)
        end do
        call ocean2ice(ice, dyn, tracers, m, p)
        do n = 1, nNodO
            if (.not. alltri(n)) cycle
            uo = 0.0_WP; vo = 0.0_WP; wa = 0.0_WP
            do k = 1, m%nod_in_elem2D_num(n)
                e = m%nod_in_elem2D(k, n)
                uo = uo + dyn%uv(1, 1, e)*real(m%elem_area(e), WP)
                vo = vo + dyn%uv(2, 1, e)*real(m%elem_area(e), WP)
                wa = wa + real(m%elem_area(e), WP)
            end do
            d5(2) = max(d5(2), abs(ice%srfoce_u(n) - uo/wa), abs(ice%srfoce_v(n) - vo/wa))
        end do
        d5(2) = d5(2)/0.1_WP
        call max_all(d5(1), p); call max_all(d5(2), p)

        !--- D6 (np > 1): the sigma copies of the shared elements stay bitwise equal ---
        d6 = 0.0_WP; d6s = 0.0_WP
        if (p%npes > 1) then
            do iv = 0, 1
                call ice_state_uniform(ice, nNodL)
                call hashed_velocity(ice, ngid, nNodL)        ! mEVP sigma stays 0 at rest
                do n = 1, nNodL
                    ice%data(1)%values(n) = 0.5_WP + 0.4_WP*hval(ngid(n), 41)
                    ice%data(2)%values(n) = 2.0_WP + 1.0_WP*hval(ngid(n), 42)
                    ice%srfoce_ssh(n)     = 0.1_WP*hval(ngid(n), 21)
                end do
                ice%whichEVP = iv; ice%evp_rheol_steps = 3
                call EVPdynamics_solve(ice, m, p)
                d6(iv+1) = max(shared_copy_max_elem(ice%work%sigma11, p), &
                               shared_copy_max_elem(ice%work%sigma12, p), &
                               shared_copy_max_elem(ice%work%sigma22, p))
                call max_all(d6(iv+1), p)
                do k = 1, size(p%com_elem2D_shared%rlist)
                    e = p%com_elem2D_shared%rlist(k)
                    d6s(iv+1) = max(d6s(iv+1), abs(ice%work%sigma11(e)), abs(ice%work%sigma12(e)), &
                                    abs(ice%work%sigma22(e)))
                end do
                call max_all(d6s(iv+1), p)
            end do
        end if

        if (p%mype == 0) then
            write(*,'(a,a,i0,a,i0)') tag, ': owned triangles=', ntri, ' quads=', nquad
            if (planar) then
                write(*,'(a,a,3es10.2)') tag, ': D1 strain of uniform drift / rigid rotation / linear deformation:', d1(1:3)
                write(*,'(a,a,2es10.2,a,2es10.2,a)') tag, ': D4 uniform strain, interior update EVP / mEVP:', d4, &
                    '  (max|sigma11|', d4s, ' Pa)'
            else
                write(*,'(a,a,2es10.2)') tag, ': D1 rigid rotation about the pole max|eps|/Omega (without metric):', d1s
            end if
            write(*,'(a,a,es10.2)') tag, ': D1 triangles vs former strain formulas:', d1t
            write(*,'(a,a,2es10.2)') tag, ': D2 adjoint identity / triangles vs former stress2rhs:', d2, d2t
            write(*,'(a,a,2es10.2)') tag, ': D3 SSH tilt EVP / mEVP:', d3
            write(*,'(a,a,2es10.2)') tag, ': D5 ocean2ice uniform / triangles vs former gather:', d5
            if (p%npes > 1) write(*,'(a,a,2es10.2,a,2es10.2,a)') tag, &
                ': D6 shared-copy spread of sigma EVP / mEVP:', d6, '  (max|sigma| on shared', d6s, ' Pa)'
        end if
        if (planar) then
            call check(d1(1) < 1.0e-13_WP, tag//': D1 uniform drift has zero strain rate', p%mype)
            call check(d1(2) < 1.0e-10_WP, tag//': D1 rigid rotation has zero strain rate', p%mype)
            call check(d1(3) < 1.0e-10_WP, tag//': D1 linear deformation gives its exact strain rate', p%mype)
            call check(d4(1) < 1.0e-10_WP .and. d4s(1) > 1.0_WP, &
                       tag//': D4 EVP: uniform strain is force-free at interior nodes', p%mype)
            call check(d4(2) < 1.0e-10_WP .and. d4s(2) > 1.0_WP, &
                       tag//': D4 mEVP: uniform strain is force-free at interior nodes', p%mype)
        else
            call check(d1s(1) < 0.1_WP*d1s(2), &
                       tag//': D1 rigid rotation about the pole: metric terms cancel du/dy', p%mype)
        end if
        call check(d1t < 1.0e-13_WP, tag//': D1 triangles: strain rates == former P1 formulas', p%mype)
        call check(d2  < 1.0e-13_WP, tag//': D2 stress2rhs == -adjoint of the strain rates', p%mype)
        call check(d2t < 1.0e-13_WP, tag//': D2 triangles: stress2rhs == former formula', p%mype)
        call check(d3(1) < 1.0e-12_WP, tag//': D3 EVP SSH tilt node shares', p%mype)
        call check(d3(2) < 1.0e-12_WP, tag//': D3 mEVP SSH tilt node shares', p%mype)
        call check(d5(1) < 1.0e-14_WP, tag//': D5 ocean2ice uniform velocity gathered exactly', p%mype)
        call check(d5(2) < 1.0e-14_WP, tag//': D5 triangles: ocean2ice == former area-weighted mean', p%mype)
        if (p%npes > 1) then
            call check(d6(1) == 0.0_WP .and. d6s(1) > 0.0_WP, &
                       tag//': D6 EVP: sigma copies of shared elements bitwise equal (and nonzero)', p%mype)
            call check(d6(2) == 0.0_WP .and. d6s(2) > 0.0_WP, &
                       tag//': D6 mEVP: sigma copies of shared elements bitwise equal (and nonzero)', p%mype)
        end if
        deallocate(ngid, egid, X, Y, ro, rabs, alltri)
    end subroutine check_ice

    subroutine strain(ice, m, p)
        type(t_ice),    intent(inout) :: ice
        type(t_mesh),   intent(in)    :: m
        type(t_partit), intent(in)    :: p
        ice%work%eps11 = 0.0_WP; ice%work%eps22 = 0.0_WP; ice%work%eps12 = 0.0_WP
        call stress_tensor(ice, m, p)
    end subroutine strain

    real(kind=WP) function max_eps(ice, nElemO)
        type(t_ice), intent(in) :: ice
        integer,     intent(in) :: nElemO
        max_eps = max(maxval(abs(ice%work%eps11(1:nElemO))), maxval(abs(ice%work%eps22(1:nElemO))), &
                      maxval(abs(ice%work%eps12(1:nElemO))))
    end function max_eps

    subroutine former_strain(m, ice, e, e11, e22, e12)
        ! FESOM2 stress_tensor strain rates of a triangle (1:3 slices, metric mean /3)
        type(t_mesh),  intent(in)  :: m
        type(t_ice),   intent(in)  :: ice
        integer,       intent(in)  :: e
        real(kind=WP), intent(out) :: e11, e22, e12
        integer :: nn(3)
        nn = m%elem2D_nodes(1:3, e)
        e11 = sum(m%gradient_sca(1:3, e)*ice%uice(nn)) - m%metric_factor(e)*sum(ice%vice(nn))/3.0_WP
        e22 = sum(m%gradient_sca(5:7, e)*ice%vice(nn))
        e12 = 0.5_WP*(sum(m%gradient_sca(5:7, e)*ice%uice(nn)) + sum(m%gradient_sca(1:3, e)*ice%vice(nn)) &
                      + m%metric_factor(e)*sum(ice%uice(nn))/3.0_WP)
    end subroutine former_strain

    subroutine hashed_velocity(ice, ngid, nn)
        type(t_ice), intent(inout) :: ice
        integer,     intent(in)    :: ngid(:), nn
        integer :: n
        do n = 1, nn
            ice%uice(n) = 0.3_WP*hval(ngid(n), 1)
            ice%vice(n) = 0.3_WP*hval(ngid(n), 2)
        end do
    end subroutine hashed_velocity

    subroutine ice_state_uniform(ice, nn)
        ! uniform ice (a = 1, m = 2 m, snow 0.3 m) at rest, no wind, no ocean current
        type(t_ice), intent(inout) :: ice
        integer,     intent(in)    :: nn
        ice%data(1)%values(1:nn) = 1.0_WP
        ice%data(2)%values(1:nn) = 2.0_WP
        ice%data(3)%values(1:nn) = 0.3_WP
        ice%uice = 0.0_WP; ice%vice = 0.0_WP; ice%uice_aux = 0.0_WP; ice%vice_aux = 0.0_WP
        ice%srfoce_u = 0.0_WP; ice%srfoce_v = 0.0_WP
        ice%stress_atmice_x = 0.0_WP; ice%stress_atmice_y = 0.0_WP
        ice%work%sigma11 = 0.0_WP; ice%work%sigma12 = 0.0_WP; ice%work%sigma22 = 0.0_WP
    end subroutine ice_state_uniform

    subroutine alloc_ice(ice, m, p, nn, ne)
        ! the EVP/mEVP subset of ice_allocate + ice_setup (ice_dt = the CORE2 1800 s,
        ! namelist doubles for cd_oce_ice / delta_min) and the mEVP boundary mask
        ! (mod_ice_setup build_bc_index_nod2D)
        type(t_ice),    intent(inout) :: ice
        type(t_mesh),   intent(in)    :: m
        type(t_partit), intent(in)    :: p
        integer,        intent(in)    :: nn, ne
        integer :: k
        allocate(ice%uice(nn), ice%vice(nn), ice%uice_old(nn), ice%vice_old(nn))
        allocate(ice%uice_rhs(nn), ice%vice_rhs(nn), ice%uice_aux(nn), ice%vice_aux(nn))
        allocate(ice%stress_atmice_x(nn), ice%stress_atmice_y(nn), ice%bc_index_nod2D(nn))
        allocate(ice%srfoce_u(nn), ice%srfoce_v(nn), ice%srfoce_temp(nn), ice%srfoce_salt(nn), ice%srfoce_ssh(nn))
        allocate(ice%data(ice%num_itracers))
        do k = 1, ice%num_itracers
            allocate(ice%data(k)%values(nn), ice%data(k)%values_rhs(nn))
            ice%data(k)%values = 0.0_WP; ice%data(k)%values_rhs = 0.0_WP
        end do
        allocate(ice%work%sigma11(ne), ice%work%sigma12(ne), ice%work%sigma22(ne))
        allocate(ice%work%eps11(ne), ice%work%eps12(ne), ice%work%eps22(ne))
        allocate(ice%work%ice_strength(ne), ice%work%inv_areamass(nn), ice%work%inv_mass(nn))
        ice%uice = 0.0_WP; ice%vice = 0.0_WP; ice%uice_old = 0.0_WP; ice%vice_old = 0.0_WP
        ice%uice_rhs = 0.0_WP; ice%vice_rhs = 0.0_WP; ice%uice_aux = 0.0_WP; ice%vice_aux = 0.0_WP
        ice%stress_atmice_x = 0.0_WP; ice%stress_atmice_y = 0.0_WP
        ice%srfoce_u = 0.0_WP; ice%srfoce_v = 0.0_WP; ice%srfoce_temp = 0.0_WP
        ice%srfoce_salt = 0.0_WP; ice%srfoce_ssh = 0.0_WP
        ice%work%sigma11 = 0.0_WP; ice%work%sigma12 = 0.0_WP; ice%work%sigma22 = 0.0_WP
        ice%work%eps11 = 0.0_WP; ice%work%eps12 = 0.0_WP; ice%work%eps22 = 0.0_WP
        ice%work%ice_strength = 0.0_WP; ice%work%inv_areamass = 0.0_WP; ice%work%inv_mass = 0.0_WP
        ice%ice_dt     = 1800.0_WP
        ice%Tevp_inv   = 3.0_WP/ice%ice_dt
        ice%cd_oce_ice = 0.0055_WP
        ice%delta_min  = 1.0e-11_WP
        call build_bc_index_nod2D(ice, m, p)
    end subroutine alloc_ice

end program test_ice_dyn_mixed
