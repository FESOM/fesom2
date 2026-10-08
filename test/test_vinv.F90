program test_vinv
    ! Vector-invariant momentum advection (momadv_opt==1, src/oce/oce_dyn_vinv.F90).
    !
    ! WRITTEN BEFORE THE OPERATORS (TDD), because there is no oracle: the FESOM2 v2.7.3
    ! momadv_opt==1 branch aborts ("not adapted mom_adv advection typ for ALE"), so these
    ! assertions ARE the specification. See
    ! docs/plans/2026-10-01-momadv-vector-invariant.md.
    !
    ! WHY cartesian=.true.
    ! --------------------
    ! On the SPHERICAL metric the median-dual contour does not close: edge_cross_dxdy(1:2)
    ! carries elem_cos(el1) and (3:4) carries elem_cos(el2) at the SAME edge midpoint
    ! (src/mesh/mod_mesh_areas.F90:336-347), so a uniform field leaves a residual of
    ! relative size (L/R)*tan(lat) ~ 3% on pi -- about 1e9 x round-off. No exactness
    ! assertion is possible there. compute_geometry(..., cartesian=.true.) sets
    ! elem_cos = 1 and metric_factor = 0 (mod_mesh_areas.F90:239-242), the contour closes,
    ! and the discrete circulation integral becomes EXACT for a velocity field linear in
    ! (x,y): for straight segments the quadrature of the closed line integral is exact, so
    ! zeta = circulation/area reproduces curl(u) to round-off. That is what makes V2 below
    ! a real test rather than a tolerance guess.
    !
    ! WHAT IS ASSERTED
    ! ----------------
    !   V1 uniform flow     -> zeta == 0 to round-off. Catches a sign error or a missing
    !                          one-sided level range: a constant field must cancel edge by
    !                          edge.
    !   V2 linear shear     -> u = alpha*y, v = 0  =>  zeta == -alpha to ROUND-OFF.
    !                          THE load-bearing test. Getting -alpha exactly requires
    !                          area(n) == sum over wet adjacent elements of (1/3)*elem_area,
    !                          so this is what pins the area(nz,n) -> area(n) deviation.
    !                          The field depends on LATITUDE ONLY, so it is continuous
    !                          across pi's cyclic seam (a field linear in longitude would
    !                          not be, and would give garbage there).
    !   V3 linearity        -> zeta(-u) == -zeta(u) exactly.
    !   V4 non-degeneracy   -> zeta /= 0 somewhere. Catches a stub or a no-op (the whole
    !                          point of writing this test before the implementation).
    !   V5 finiteness       -> no NaN/Inf anywhere, including at bathymetry steps.
    !
    ! RESTRICTIONS on every assertion: INTERIOR nodes only (the dual contour is open at a
    ! boundary node, so zeta there is not curl(u) -- pi has ~455 boundary edges), and only
    ! levels nz < nlevels_nod2D_min(n) where ALL adjacent elements are wet (at a bathymetry
    ! step the one-sided ranges deliberately drop the dry element's contribution while
    ! area(n) keeps its full weight, so zeta is not curl(u) there either -- that is the
    ! documented zero-with-full-weight rule, not an error).
    use mpi
    use mod_precision,     only: WP, MP
    use mod_constants,     only: r_earth
    use mod_mesh,          only: t_mesh
    use mod_dyn,           only: t_dyn
    use mod_partit,        only: t_partit
    use mod_partitioning,  only: par_init, par_ex, set_partition
    use mod_param_phys,    only: rvo_upwind
    use mod_mesh_rotate,   only: trim_cyclic
    use mod_mesh_read,     only: read_mesh
    use mod_mesh_areas,    only: compute_geometry
    use mod_part_bounds,   only: owned_bounds
    use oce_dyn_vinv,      only: relative_vorticity, momentum_adv_vinv
    use mod_halo,          only: exchange_elem
    implicit none

    character(len=512) :: mesh_dir
    type(t_partit) :: partit
    type(t_mesh)   :: mesh
    type(t_dyn)    :: dyn
    integer :: nfail, nsw, ierr
    integer :: nNodO, nNodL, nEdgeO, nElemO, nElemF, nl, n, nz
    real(kind=WP) :: zmax
    real(kind=WP), allocatable :: zsave(:,:)
    logical,       allocatable :: interior(:), elem_int(:), elem_int2(:)
    real(kind=WP) :: scal
    integer       :: nzall
    real(kind=WP), allocatable :: rhsA(:,:,:,:)
    real(kind=WP), parameter   :: ALPHA = 1.0e-5_WP   ! du/dy [1/s]

    nfail = 0
    call get_environment_variable('FESOM3_MESH_DIR', mesh_dir)
    if (len_trim(mesh_dir) == 0) &
        mesh_dir = '/home/a/a270088/port2/fesom2/tests/data/MESHES/pi'
    call par_init(partit)

    if (partit%npes > 1) call set_partition(partit, trim(mesh_dir))
    call read_mesh(mesh, partit, trim(mesh_dir), 50.0_WP, 15.0_WP, -90.0_WP, 360.0_WP, &
                   force_rotation=.true., n_cw_swaps=nsw)
    ! cartesian: makes the discrete curl exact for a linear field (see header)
    call compute_geometry(mesh, partit, cartesian=.true.)
    call owned_bounds(mesh, nNodO, nNodL, nEdgeO, nElemO, partit)
    if (partit%npes == 1) then
        nElemF = mesh%elem2D
    else
        nElemF = partit%myDim_elem2D + partit%eDim_elem2D + partit%eXDim_elem2D
    end if
    nl = mesh%nl

    allocate(dyn%uv(2, nl-1, nElemF))
    allocate(dyn%work%vorticity(nl-1, nNodL))
    allocate(zsave(nl-1, nNodL))
    allocate(interior(nNodL))
    dyn%uv = 0.0_WP; dyn%work%vorticity = 0.0_WP

    ! interior = touched by no boundary edge (edge_tri(2,ed) == edge_tri(1,ed) is the
    ! partition-robust boundary test;the ed > edge2D_in form needs myList_edge2D at npes>1)
    interior = .true.
    do n = 1, nEdgeO                 ! LOCAL owned edges; mesh%edge2D is the GLOBAL count
        if (mesh%edge_tri(2,n) == mesh%edge_tri(1,n)) then
            if (mesh%edges(1,n) <= nNodL) interior(mesh%edges(1,n)) = .false.
            if (mesh%edges(2,n) <= nNodL) interior(mesh%edges(2,n)) = .false.
        end if
    end do
    write(*,'(a,i7,a,i7)') '  interior nodes: ', count(interior(1:nNodO)), ' of ', nNodO

    !=========================================================================
    ! A - elem_neighbors adjacency (build_elem_adjacency, Task 1 of rvo_upwind)
    !=========================================================================
    call check_true('A1 elem_neighbors allocated+resolved', &
                    allocated(mesh%elem_neighbors) .and. all(mesh%elem_neighbors(:,1:nElemO) >= 1))
    block
        integer :: e, k, nb, nshared, i, j, nbnd_slots, nbnd_edges
        logical :: mutual
        nbnd_slots = 0
        do e = 1, nElemO
            do k = 1, 3
                nb = mesh%elem_neighbors(k, e)
                if (nb == e) then
                    nbnd_slots = nbnd_slots + 1
                    cycle
                end if
                if (nb > nElemO) cycle                     ! halo neighbour: no local data
                ! A2 mutuality
                mutual = any(mesh%elem_neighbors(:, nb) == e)
                if (.not. mutual) call check_true('A2 neighbour relation mutual', .false.)
                ! A4 slot alignment: e and nb share EXACTLY the slot-k vertex pair
                nshared = 0
                do i = 1, 3
                    do j = 1, 3
                        if (mesh%elem2D_nodes(i, e) == mesh%elem2D_nodes(j, nb)) nshared = nshared + 1
                    end do
                end do
                if (nshared /= 2) call check_true('A4 neighbours share exactly 2 vertices', .false.)
                if (.not. (any(mesh%elem2D_nodes(:, nb) == mesh%elem2D_nodes(k, e)) .and. &
                           any(mesh%elem2D_nodes(:, nb) == mesh%elem2D_nodes(mod(k,3)+1, e)))) &
                    call check_true('A4 slot k matches the shared vertex pair', .false.)
            end do
        end do
        ! A3 boundary accounting, PER RANK: self-slots of owned elements == boundary edges
        ! whose el1 is owned. Do NOT expect the cross-rank sum to equal the global boundary
        ! count: FESOM2 element ownership OVERLAPS at partition seams (myDim elements =
        ! elements touching an owned node), so a boundary-strip element -- and hence its
        ! boundary edge -- is legitimately counted by several ranks (pi dist_2: 6 such).
        nbnd_edges = 0
        do i = 1, size(mesh%edge_tri, 2)
            if (mesh%edge_tri(2, i) == mesh%edge_tri(1, i) .and. mesh%edge_tri(1, i) >= 1 .and. &
                mesh%edge_tri(1, i) <= nElemO) nbnd_edges = nbnd_edges + 1
        end do
        write(*,'(a,i6,a,i6)') '  A3 boundary slots = ', nbnd_slots, '  boundary edges = ', nbnd_edges
        call check_true('A3 boundary slots == owned boundary edges', nbnd_slots == nbnd_edges)
    end block

    !-------------------------------------------------------------------------
    ! A5 - the halo assumption rvo_upwind stands on (revalidation of 6e38cf6).
    ! momentum_adv_vinv fills omega_e over OWNED elements, calls the STANDARD
    ! exchange_elem (com_elem2D = the eDim ring) and then reads omega_e at every
    ! elem_neighbors id. That is right only if every edge-neighbour of an owned element
    ! lies in the ring exchange_elem fills -- never in eXDim, which only
    ! exchange_elem_full touches and which would otherwise be read as the INITIAL value
    ! (silently, no bounds error: the array spans nElemF). By construction it holds
    ! (FESOM2 gen_comm.F90: com_elem2D = "only those sharing an edge", com_elem2D_full =
    ! "all neighbors"); A5 makes it a tested invariant using the production pattern
    ! itself: owned entries <- own global id, halo entries <- a sentinel, ONE
    ! exchange_elem, and every neighbour slot must then hold that neighbour's global id.
    ! A surviving sentinel = an unfilled ring or a receive list pointing elsewhere; a
    ! wrong id = a sender/receiver mismatch. The halo-slot count printed is the number
    ! of slots that pass ONLY because the exchange filled them (they start as the
    ! sentinel), i.e. the test's teeth; at npes > 1 it must be > 0 on every rank.
    block
        real(kind=WP), allocatable :: f(:,:)
        real(kind=WP), parameter   :: SENTINEL = -1.0_WP
        integer, allocatable :: gid(:)
        integer :: e, k, nb, nhalo, nbad, nring
        ! global element id of each LOCAL id; the identity at npes == 1, where this test
        ! never calls set_partition and partit%myList_elem2D stays unallocated
        allocate(gid(nElemF))
        if (partit%npes > 1) then
            gid = partit%myList_elem2D(1:nElemF)
        else
            gid = [(e, e = 1, nElemF)]
        end if
        allocate(f(nl-1, nElemF))
        f = SENTINEL
        do e = 1, nElemO
            f(:, e) = real(gid(e), WP)
        end do
        if (partit%npes > 1) call exchange_elem(f, partit)
        nhalo = 0; nbad = 0
        do e = 1, nElemO
            do k = 1, 3
                nb = mesh%elem_neighbors(k, e)
                if (nb == e) cycle
                if (nb > nElemO) nhalo = nhalo + 1
                if (f(1, nb)    /= real(gid(nb), WP) .or. &
                    f(nl-1, nb) /= real(gid(nb), WP)) nbad = nbad + 1
            end do
        end do
        if (partit%npes > 1) then
            nring = partit%myDim_elem2D + partit%eDim_elem2D
        else
            nring = nElemF                 ! identity partition: no halo at all
        end if
        write(*,'(a,i6,a,i6,a,i7)') '  A5 halo neighbour slots = ', nhalo, &
            '  unreadable/mismatched = ', nbad, '  eDim ring end = ', nring
        call check_true('A5 every neighbour readable after the STANDARD exchange_elem', nbad == 0)
        call check_true('A5 every neighbour id inside the eDim ring (never eXDim)', &
                        all(mesh%elem_neighbors(:, 1:nElemO) <= nring))
        if (partit%npes > 1) &
            call check_true('A5 not vacuous: this rank has halo neighbour slots', nhalo > 0)
        deallocate(f, gid)
    end block

    !=========================================================================
    ! V1 - uniform flow must give zeta == 0 EXACTLY
    !=========================================================================
    dyn%uv(1,:,:) =  0.17_WP
    dyn%uv(2,:,:) = -0.43_WP
    call relative_vorticity(dyn, mesh, partit)
    zmax = 0.0_WP
    do n = 1, nNodO
        if (.not. interior(n)) cycle
        do nz = mesh%ulevels_nod2D(n), mesh%nlevels_nod2D_min(n)-1
            zmax = max(zmax, abs(dyn%work%vorticity(nz,n)))
        end do
    end do
    call max_all(zmax)
    write(*,'(a,es12.4)') '  V1 uniform flow   : max|zeta| = ', zmax
    call check_true('V1 uniform flow -> zeta == 0 (interior, all-wet)', zmax < 1.0e-16_WP)

    !=========================================================================
    ! V2 / V4 / V5 - linear shear u = alpha*y : zeta must be EXACTLY -alpha
    !=========================================================================
    call set_shear(1.0_WP)
    call relative_vorticity(dyn, mesh, partit)
    zsave = dyn%work%vorticity

    zmax = owned_absmax(dyn%work%vorticity)
    write(*,'(a,es12.4)') '  V4 non-degeneracy : max|zeta| = ', zmax
    call check_true('V4 zeta /= 0 somewhere (operator is not a no-op)', zmax > 1.0e-12_WP)

    zmax = 0.0_WP
    do n = 1, nNodO
        if (.not. interior(n)) cycle
        do nz = mesh%ulevels_nod2D(n), mesh%nlevels_nod2D_min(n)-1
            zmax = max(zmax, abs(dyn%work%vorticity(nz,n) + ALPHA))
        end do
    end do
    call max_all(zmax)
    write(*,'(a,es12.4,a,es12.4)') '  V2 linear shear   : max|zeta-(-alpha)| = ', zmax, &
                                   '   relative = ', zmax/ALPHA
    ! DIAGNOSTIC: is the error systematic (a factor) or localised (a few bad nodes)?
    block
        integer :: nok, nbad, nworst
        real(kind=WP) :: r, racc, rworst, latw
        nok = 0; nbad = 0; racc = 0.0_WP; rworst = 0.0_WP; nworst = 0
        do n = 1, nNodO
            if (.not. interior(n)) cycle
            do nz = mesh%ulevels_nod2D(n), mesh%nlevels_nod2D_min(n)-1
                r = dyn%work%vorticity(nz,n)/(-ALPHA)
                racc = racc + r; nok = nok + 1
                if (abs(r - 1.0_WP) > 0.01_WP) nbad = nbad + 1
                if (abs(r - 1.0_WP) > rworst) then
                    rworst = abs(r - 1.0_WP); nworst = n
                end if
            end do
        end do
        if (nok > 0) then
            latw = real(mesh%coord_nod2D(2, max(nworst,1)), WP)*57.29578_WP
            write(*,'(a,i7,a,f9.4,a,i7,a,f7.2,a,f8.2)') &
                '     levels=', nok, '  mean zeta/(-alpha)=', racc/real(nok,WP), &
                '  off-by->1%=', nbad, '  worst ratio off=', rworst, ' at lat=', latw
        end if
    end block
    call check_true('V2 u=alpha*y -> zeta == -alpha to round-off (interior, all-wet)', &
                    zmax/ALPHA < 1.0e-11_WP)

    call check_true('V5 finite', all(dyn%work%vorticity == dyn%work%vorticity) .and. &
                                 all(abs(dyn%work%vorticity) < huge(1.0_WP)))

    !=========================================================================
    ! V3 - linearity: reversing the flow must flip zeta exactly
    !=========================================================================
    call set_shear(-1.0_WP)
    call relative_vorticity(dyn, mesh, partit)
    zmax = owned_absmax(dyn%work%vorticity + zsave)
    write(*,'(a,es12.4)') '  V3 linearity      : max|zeta(-u)+zeta(u)| = ', zmax
    call check_true('V3 zeta(-u) == -zeta(u)', zmax < 1.0e-16_WP)

    !=========================================================================
    ! V6 - full operator: uniform u with NONZERO w must give ZERO tendency.
    ! This is the Block C net. For a uniform u: zeta == 0 (V1), grad(KE) == 0 at interior
    ! elements, and the vertical flux must telescope EXACTLY -- including across the bottom
    ! face, which is where qq's uvert(nl1+1)=0 breaks under bottom-at-vertices. A `w == 0`
    ! test cannot see that, since every vertical term is proportional to w.
    !=========================================================================
    allocate(mesh%helem(nl-1, nElemF))
    allocate(dyn%uv_rhsAB(1, 2, nl-1, nElemF))
    allocate(dyn%w_e(nl, nNodL))
    allocate(elem_int(nElemO))
    allocate(rhsA(1, 2, nl-1, nElemF))
    do nz = 1, nl-1
        mesh%helem(nz,:) = 10.0_MP + 2.0_MP*real(nz, MP)   ! any positive ALE thickness
    end do
    do n = 1, nNodL
        do nz = 1, nl
            dyn%w_e(nz,n) = 1.0e-5_WP*real(nl-nz, WP)      ! nonzero AND depth-varying
        end do
    end do
    do n = 1, nElemO
        elem_int(n) = all(interior(mesh%elem2D_nodes(1:3,n)))
    end do
    ! one ring deeper: the upwind blend reads the NEIGHBOURS' vertex vorticity, so a
    ! constant-field inertness assertion must exclude elements whose neighbours touch the
    ! boundary too (and, at np>1, elements with a halo neighbour whose vertices cannot be
    ! inspected locally -- elem2D_nodes is owned-only).
    allocate(elem_int2(nElemO))
    do n = 1, nElemO
        elem_int2(n) = elem_int(n)
        do nz = 1, 3
            block
                integer :: nbb
                nbb = mesh%elem_neighbors(nz, n)
                if (nbb == n .or. nbb > nElemO) then
                    elem_int2(n) = .false.
                else if (.not. elem_int(nbb)) then
                    elem_int2(n) = .false.
                end if
            end block
        end do
    end do

    dyn%uv(1,:,:) =  0.17_WP
    dyn%uv(2,:,:) = -0.43_WP
    dyn%uv_rhsAB  = 0.0_WP
    call momentum_adv_vinv(dyn, mesh, partit)

    ! Restricted to levels where ALL elements adjacent to ALL THREE nodes are wet. Below
    ! that, KE(nz,n) sums only over the elements wet at nz while dividing by the FULL
    ! area(n) -- the zero-with-full-weight rule -- so KE is deliberately not uniform there
    ! and grad(KE) /= 0 even for a uniform u. That is intended, not an error.
    zmax = 0.0_WP; scal = 0.0_WP
    do n = 1, nElemO
        if (.not. elem_int(n)) cycle
        nzall = minval(mesh%nlevels_nod2D_min(mesh%elem2D_nodes(1:3,n))) - 1
        do nz = mesh%ulevels(n), min(mesh%nlevels(n)-1, nzall)
            zmax = max(zmax, abs(dyn%uv_rhsAB(1,1,nz,n)), abs(dyn%uv_rhsAB(1,2,nz,n)))
            ! scale of the individual vertical terms, so the tolerance is relative
            scal = max(scal, abs(dyn%w_e(nz,mesh%elem2D_nodes(1,n))*0.43_WP) &
                             *real(mesh%elem_area(n),WP)/real(mesh%helem(nz,n),WP))
        end do
    end do
    call max_all(zmax); call max_all(scal)
    write(*,'(a,es12.4,a,es12.4)') '  V6 uniform u, w/=0: max|UV_rhsAB| = ', zmax, &
                                   '   relative = ', zmax/max(scal, tiny(1.0_WP))
    call check_true('V6 uniform u + nonzero w -> zero tendency (interior elems, incl bottom face)', &
                    zmax/max(scal, tiny(1.0_WP)) < 1.0e-12_WP)

    !=========================================================================
    ! V7 - the full operator is not a no-op
    !=========================================================================
    call set_shear(1.0_WP)
    dyn%uv_rhsAB = 0.0_WP
    call momentum_adv_vinv(dyn, mesh, partit)
    zmax = 0.0_WP
    do n = 1, nElemO
        do nz = mesh%ulevels(n), mesh%nlevels(n)-1
            zmax = max(zmax, abs(dyn%uv_rhsAB(1,1,nz,n)), abs(dyn%uv_rhsAB(1,2,nz,n)))
        end do
    end do
    call max_all(zmax)
    write(*,'(a,es12.4)') '  V7 full operator  : max|UV_rhsAB| = ', zmax
    call check_true('V7 full operator is not a no-op', zmax > 0.0_WP)

    !=========================================================================
    ! V8-V11 - upwind-blended face vorticity (rvo_upwind), written BEFORE the blend
    ! exists (TDD): V11's single-neighbour equality is the failing specification.
    !=========================================================================
    ! V8: rvo_upwind = 0 must be exactly the reference behaviour
    rvo_upwind = 0.0_WP
    call set_shear(1.0_WP)
    dyn%uv_rhsAB = 0.0_WP
    call momentum_adv_vinv(dyn, mesh, partit)
    rhsA = dyn%uv_rhsAB
    dyn%uv_rhsAB = 0.0_WP
    call momentum_adv_vinv(dyn, mesh, partit)
    call check_true('V8 rvo_upwind=0 bit-identical', all(dyn%uv_rhsAB == rhsA))

    ! V9a: constant zeta (linear shear) -> blend inert at rvo_upwind=1 (weighted average
    ! of equal values re-forms the same value up to round-off)
    rvo_upwind = 1.0_WP
    dyn%uv_rhsAB = 0.0_WP
    call momentum_adv_vinv(dyn, mesh, partit)
    zmax = 0.0_WP; scal = 0.0_WP
    do n = 1, nElemO
        if (.not. elem_int2(n)) cycle
        nzall = minval(mesh%nlevels_nod2D_min(mesh%elem2D_nodes(1:3,n))) - 1
        do nz = 1, 3
            nzall = min(nzall, minval(mesh%nlevels_nod2D_min( &
                        mesh%elem2D_nodes(1:3, mesh%elem_neighbors(nz,n)))) - 1)
        end do
        do nz = mesh%ulevels(n), min(mesh%nlevels(n)-1, nzall)
            zmax = max(zmax, abs(dyn%uv_rhsAB(1,1,nz,n)-rhsA(1,1,nz,n)), &
                             abs(dyn%uv_rhsAB(1,2,nz,n)-rhsA(1,2,nz,n)))
            scal = max(scal, abs(rhsA(1,1,nz,n)), abs(rhsA(1,2,nz,n)))
        end do
    end do
    call max_all(zmax); call max_all(scal)
    write(*,'(a,es12.4)') '  V9a const-zeta blend inertness: rel = ', zmax/max(scal,tiny(1.0_WP))
    call check_true('V9a constant zeta -> blend inert at rvo_upwind=1', &
                    zmax/max(scal,tiny(1.0_WP)) < 1.0e-12_WP)

    ! V9b: uniform u + nonzero w, rvo_upwind=1 -> the V6 zero-tendency bound must hold
    dyn%uv(1,:,:) =  0.17_WP
    dyn%uv(2,:,:) = -0.43_WP
    dyn%uv_rhsAB  = 0.0_WP
    call momentum_adv_vinv(dyn, mesh, partit)
    zmax = 0.0_WP; scal = 0.0_WP
    do n = 1, nElemO
        if (.not. elem_int2(n)) cycle
        nzall = minval(mesh%nlevels_nod2D_min(mesh%elem2D_nodes(1:3,n))) - 1
        do nz = 1, 3
            nzall = min(nzall, minval(mesh%nlevels_nod2D_min( &
                        mesh%elem2D_nodes(1:3, mesh%elem_neighbors(nz,n)))) - 1)
        end do
        do nz = mesh%ulevels(n), min(mesh%nlevels(n)-1, nzall)
            zmax = max(zmax, abs(dyn%uv_rhsAB(1,1,nz,n)), abs(dyn%uv_rhsAB(1,2,nz,n)))
            scal = max(scal, abs(dyn%w_e(nz,mesh%elem2D_nodes(1,n))*0.43_WP) &
                             *real(mesh%elem_area(n),WP)/real(mesh%helem(nz,n),WP))
        end do
    end do
    call max_all(zmax); call max_all(scal)
    write(*,'(a,es12.4)') '  V9b uniform-u zero tendency at rvo=1: rel = ', zmax/max(scal,tiny(1.0_WP))
    call check_true('V9b uniform u + w/=0 -> zero tendency at rvo_upwind=1', &
                    zmax/max(scal,tiny(1.0_WP)) < 1.0e-12_WP)

    ! V10 + V11 share the test-side reconstruction of normals/weights
    call upwind_checks()
    rvo_upwind = 0.0_WP

    if (nfail == 0) then
        write(*,'(a)') 'test_vinv: OK'
    else
        write(*,'(a,i0,a)') 'test_vinv: ', nfail, ' FAILURE(S)'
    end if
    call par_ex(partit%MPI_COMM_FESOM, partit%mype)   ! NO third arg: abort is presence-based
    if (nfail /= 0) error stop 1

contains

    subroutine upwind_checks()
        ! V10: the inward-normal/weight convention, proven against geometry: every edge
        ! the TEST's reconstruction marks as inflow (w>0) must have its neighbour's
        ! centroid UPSTREAM of the element's: (c_nb - c_e).u < 0. Pure geometry; the
        ! model couples to it through V11, which uses the SAME convention.
        ! V11: end-to-end weighting through the public output. u = C*lat^2 (zeta ~ lat,
        ! VARYING -- a constant-zeta field is blend-inert and cannot test selection).
        ! delta_zb is recovered from UV_rhsAB(rvo=1)-UV_rhsAB(rvo=0):
        !   drhs_x = dzb*V*area, drhs_y = -dzb*U*area  =>  dzb = (dx*V-dy*U)/((U2+V2)*area)
        ! and compared with the test-side blend of omega_e over the wet owned inflow
        ! neighbours: exactly one -> equality to round-off; several -> max principle.
        real(kind=WP), parameter :: CU = 0.3_WP, U0x = 0.31_WP, U0y = 0.17_WP
        ! (rhsB allocated on first assignment via sourced allocation below)
        real(kind=WP) :: nxk(3), nyk(3), dk(3), wk(3), cex, cey, cnx, cny
        real(kind=WP) :: latm, u_e, v_e, zb0, dzb, ea, sumw, zup, zlo, zhi, err, errmax
        real(kind=WP), allocatable :: rhsB(:,:,:,:)
        integer :: e, k, nb, nwet, nupw(1), nz1b, n1, n2, n3, nbad10, none_cnt
        ! --- fields: u = CU*lat^2, v = 0 (rotated frame) ---
        do e = 1, nElemO
            latm = (real(mesh%coord_nod2D(2,mesh%elem2D_nodes(1,e)),WP) + &
                    real(mesh%coord_nod2D(2,mesh%elem2D_nodes(2,e)),WP) + &
                    real(mesh%coord_nod2D(2,mesh%elem2D_nodes(3,e)),WP))/3.0_WP
            dyn%uv(1,:,e) = CU*latm*latm
            dyn%uv(2,:,e) = 0.0_WP
        end do
        rvo_upwind = 0.0_WP
        dyn%uv_rhsAB = 0.0_WP
        call momentum_adv_vinv(dyn, mesh, partit)
        rhsA = dyn%uv_rhsAB
        rvo_upwind = 1.0_WP
        dyn%uv_rhsAB = 0.0_WP
        call momentum_adv_vinv(dyn, mesh, partit)
        rhsB = dyn%uv_rhsAB

        nbad10 = 0; errmax = 0.0_WP; none_cnt = 0
        do e = 1, nElemO
            n1 = mesh%elem2D_nodes(1,e); n2 = mesh%elem2D_nodes(2,e); n3 = mesh%elem2D_nodes(3,e)
            call elem_normals(e, nxk, nyk, cex, cey)
            u_e = dyn%uv(1,1,e); v_e = dyn%uv(2,1,e)
            ! ---- V10, velocity-free: the slot-k neighbour must lie on the OUTWARD side
            ! of the slot-k inward normal, (c_nb - c_e).n_in < 0, in the SAME metre
            ! metric as the normals. (The first draft asserted (c_nb-c_e).u < 0 for
            ! inflow edges -- too strong: u.n_in > 0 does not bound the TANGENTIAL part
            ! of the centroid offset, so ~20% of irregular-triangle edges legitimately
            ! violated it. The velocity-free form is the actual orientation property;
            ! d = u.n_in is then correct by construction.) ----
            do k = 1, 3
                nb = mesh%elem_neighbors(k, e)
                if (nb == e .or. nb > nElemO) cycle
                call elem_centroid(nb, cex, cey, cnx, cny)
                if ( ((cnx-cex)*real(mesh%elem_cos(e),WP)*r_earth)*nxk(k) &
                   + ((cny-cey)*r_earth)*nyk(k) >= 0.0_WP ) nbad10 = nbad10 + 1
            end do
            ! ---- V11 at the surface level of all-wet elements ----
            nz1b = mesh%ulevels(e)
            if (nz1b > minval(mesh%nlevels_nod2D_min(mesh%elem2D_nodes(1:3,e))) - 1) cycle
            if (u_e*u_e + v_e*v_e < 1.0e-12_WP) cycle
            zb0 = (dyn%work%vorticity(nz1b,n1)+dyn%work%vorticity(nz1b,n2) &
                  +dyn%work%vorticity(nz1b,n3))/3.0_WP
            sumw = 0.0_WP; zup = 0.0_WP; nwet = 0
            zlo = zb0; zhi = zb0
            do k = 1, 3
                nb = mesh%elem_neighbors(k, e)
                if (nb == e .or. nb > nElemO) then
                    if (nb /= e) nwet = -999       ! halo neighbour: skip this element
                    cycle
                end if
                if (nz1b < mesh%ulevels(nb) .or. nz1b > mesh%nlevels(nb)-1) cycle
                dk(k) = u_e*nxk(k) + v_e*nyk(k)
                wk(k) = dk(k) + abs(dk(k))
                if (wk(k) > 0.0_WP) then
                    nwet = nwet + 1; nupw(1) = nb
                    sumw = sumw + wk(k)
                    zup  = zup + wk(k)*omega_of(nb, nz1b)
                    zlo  = min(zlo, omega_of(nb, nz1b)); zhi = max(zhi, omega_of(nb, nz1b))
                end if
            end do
            if (nwet <= 0) then
                none_cnt = none_cnt + 1; cycle
            end if
            ea  = real(mesh%elem_area(e), WP)
            dzb = ((dyn%uv_rhsAB(1,1,nz1b,e)-rhsA(1,1,nz1b,e))*v_e  &
                 - (dyn%uv_rhsAB(1,2,nz1b,e)-rhsA(1,2,nz1b,e))*u_e) &
                 /((u_e*u_e+v_e*v_e)*ea)
            if (nwet == 1) then
                err = abs(dzb - (omega_of(nupw(1),nz1b) - zb0)) &
                      /max(abs(omega_of(nupw(1),nz1b))+abs(zb0), 1.0e-12_WP)
                errmax = max(errmax, err)
            else
                if (zb0+dzb < zlo-1.0e-10_WP*abs(zlo) .or. zb0+dzb > zhi+1.0e-10_WP*abs(zhi)) &
                    call check_true('V11 max principle (2-neighbour blend in range)', .false.)
            end if
        end do
        ! V12: the blend is LINEAR in rvo_upwind pointwise:
        !   rhs(rvo) - rhs(0) = rvo * [ rhs(1) - rhs(0) ]   (exactly, up to one multiply)
        rvo_upwind = 0.5_WP
        dyn%uv_rhsAB = 0.0_WP
        call momentum_adv_vinv(dyn, mesh, partit)
        err = 0.0_WP; zhi = 0.0_WP
        do e = 1, nElemO
            do k = mesh%ulevels(e), mesh%nlevels(e)-1
                err = max(err, &
                    abs((dyn%uv_rhsAB(1,1,k,e)-rhsA(1,1,k,e)) - 0.5_WP*(rhsB(1,1,k,e)-rhsA(1,1,k,e))), &
                    abs((dyn%uv_rhsAB(1,2,k,e)-rhsA(1,2,k,e)) - 0.5_WP*(rhsB(1,2,k,e)-rhsA(1,2,k,e))))
                zhi = max(zhi, abs(rhsB(1,1,k,e)-rhsA(1,1,k,e)), abs(rhsB(1,2,k,e)-rhsA(1,2,k,e)))
            end do
        end do
        call max_all(err); call max_all(zhi)
        write(*,'(a,es12.4)') '  V12 linear-in-coefficient: rel = ', err/max(zhi,tiny(1.0_WP))
        call check_true('V12 blend linear in rvo_upwind', err/max(zhi,tiny(1.0_WP)) < 1.0e-12_WP)
        write(*,'(a,i6,a,es12.4,a,i6)') '  V10 bad orientations = ', nbad10, &
            '   V11 1-neighbour errmax = ', errmax, '   no-inflow elems = ', none_cnt
        call check_true('V10 inflow neighbours are upstream', nbad10 == 0)
        call check_true('V11 single-neighbour blend equals upwind omega', errmax < 1.0e-9_WP)
    end subroutine upwind_checks

    real(kind=WP) function omega_of(e, nz)
        integer, intent(in) :: e, nz
        omega_of = (dyn%work%vorticity(nz, mesh%elem2D_nodes(1,e)) &
                  + dyn%work%vorticity(nz, mesh%elem2D_nodes(2,e)) &
                  + dyn%work%vorticity(nz, mesh%elem2D_nodes(3,e)))/3.0_WP
    end function omega_of

    subroutine elem_centroid(e, refx, refy, cx, cy)
        ! centroid in the metre-scaled plane about (refx,refy) to stay cyclic-safe
        integer,       intent(in)  :: e
        real(kind=WP), intent(in)  :: refx, refy
        real(kind=WP), intent(out) :: cx, cy
        real(kind=WP) :: dxs, dys
        integer :: k
        cx = 0.0_WP; cy = 0.0_WP
        do k = 1, 3
            dxs = real(mesh%coord_nod2D(1,mesh%elem2D_nodes(k,e)),WP) - refx
            call trim_cyclic(dxs)
            dys = real(mesh%coord_nod2D(2,mesh%elem2D_nodes(k,e)),WP) - refy
            cx = cx + dxs; cy = cy + dys
        end do
        cx = refx + cx/3.0_WP; cy = refy + cy/3.0_WP
    end subroutine elem_centroid

    subroutine elem_normals(e, nxk, nyk, cex, cey)
        ! the three INWARD edge normals of element e, metre-scaled with elem_cos(e);
        ! edge k connects elnodes(k), elnodes(k+1); inward = toward elnodes(k+2).
        ! Mirrors the model's convention; V10 then checks the convention against
        ! geometry, so it is proven rather than shared-by-copy.
        integer,       intent(in)  :: e
        real(kind=WP), intent(out) :: nxk(3), nyk(3), cex, cey
        real(kind=WP) :: x(3), y(3), tx, ty, px, py, dref
        integer :: k, k2, k3
        do k = 1, 3
            x(k) = real(mesh%coord_nod2D(1, mesh%elem2D_nodes(k,e)), WP)
            y(k) = real(mesh%coord_nod2D(2, mesh%elem2D_nodes(k,e)), WP)
        end do
        ! cyclic-safe local frame about vertex 1
        do k = 2, 3
            tx = x(k) - x(1); call trim_cyclic(tx); x(k) = x(1) + tx
        end do
        cex = sum(x)/3.0_WP; cey = sum(y)/3.0_WP
        do k = 1, 3
            k2 = mod(k,3)+1; k3 = mod(k+1,3)+1
            tx = (x(k2)-x(k))*real(mesh%elem_cos(e),WP)*r_earth
            ty = (y(k2)-y(k))*r_earth
            px =  ty; py = -tx
            dref = px*(x(k3)-x(k))*real(mesh%elem_cos(e),WP)*r_earth + py*(y(k3)-y(k))*r_earth
            if (dref < 0.0_WP) then
                px = -px; py = -py
            end if
            nxk(k) = px; nyk(k) = py
        end do
    end subroutine elem_normals

    subroutine set_shear(sgn)
        ! u = sgn*ALPHA*y, v = 0 with y = r_earth*lat : LINEAR in the cartesian metric, so
        ! the discrete curl is exact and zeta == -sgn*ALPHA. Depends on latitude only, so
        ! it is continuous across pi's cyclic seam. Element value from the mean of its 3
        ! corner latitudes (exact for a linear field).
        !
        ! coord_nod2D, NOT geo_coord_nod2D: the mesh is read with force_rotation=.true., so
        ! the geometry (elem_area, edge_cross_dxdy, element centres) lives in the ROTATED
        ! computational frame. A field linear in GEOGRAPHIC latitude is not linear in that
        ! frame, and the test then fails worst near the rotated pole -- which is exactly how
        ! this was found (mean zeta/(-alpha) = 0.836, worst at rotated lat 84 deg).
        real(kind=WP), intent(in) :: sgn
        integer :: el, k
        real(kind=WP) :: latm
        ! 1..nElemO, NOT nElemF: elem2D_nodes is owned-only, and both triangles of an
        ! owned edge are owned (LESSONS L31), so the edge loop never reads a halo element.
        do el = 1, nElemO
            latm = 0.0_WP
            do k = 1, 3
                latm = latm + real(mesh%coord_nod2D(2, mesh%elem2D_nodes(k,el)), WP)
            end do
            latm = latm/3.0_WP
            dyn%uv(1,:,el) = sgn*ALPHA*r_earth*latm
            dyn%uv(2,:,el) = 0.0_WP
        end do
    end subroutine set_shear

    real(kind=WP) function owned_absmax(f)
        real(kind=WP), intent(in) :: f(:,:)
        integer :: nn, kk
        owned_absmax = 0.0_WP
        do nn = 1, nNodO
            do kk = mesh%ulevels_nod2D(nn), mesh%nlevels_nod2D(nn)-1
                owned_absmax = max(owned_absmax, abs(f(kk,nn)))
            end do
        end do
        call max_all(owned_absmax)
    end function owned_absmax

    subroutine sum_all(x)
        real(kind=WP), intent(inout) :: x
        real(kind=WP) :: t
        if (partit%npes > 1) then
            call MPI_Allreduce(x, t, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
                               partit%MPI_COMM_FESOM, ierr)
            x = t
        end if
    end subroutine sum_all

    subroutine max_all(x)
        real(kind=WP), intent(inout) :: x
        real(kind=WP) :: t
        if (partit%npes > 1) then
            call MPI_Allreduce(x, t, 1, MPI_DOUBLE_PRECISION, MPI_MAX, &
                               partit%MPI_COMM_FESOM, ierr)
            x = t
        end if
    end subroutine max_all

    subroutine sum_int(k)
        integer, intent(inout) :: k
        integer :: t
        if (partit%npes > 1) then
            call MPI_Allreduce(k, t, 1, MPI_INTEGER, MPI_SUM, &
                               partit%MPI_COMM_FESOM, ierr)
            k = t
        end if
    end subroutine sum_int

    subroutine check_true(name, cond)
        character(len=*), intent(in) :: name
        logical,          intent(in) :: cond
        if (.not. cond) then
            nfail = nfail + 1
            if (partit%mype == 0) write(*,'(a)') '  FAIL: '//name
        end if
    end subroutine check_true

end program test_vinv
