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
    !                          area(n) == sum over wet adjacent elements of the kite shares
    !                          elem_w*elem_area (1/3 on a triangle),
    !                          so this is what pins the area(nz,n) -> area(n) deviation.
    !                          The field depends on LATITUDE ONLY, so it is continuous
    !                          across pi's cyclic seam (a field linear in longitude would
    !                          not be, and would give garbage there).
    !   V3 linearity        -> zeta(-u) == -zeta(u) exactly.
    !   V4 non-degeneracy   -> zeta /= 0 somewhere. Catches a stub or a no-op (the whole
    !                          point of writing this test before the implementation).
    !   V5 finiteness       -> no NaN/Inf anywhere, including at bathymetry steps.
    !   V13 (np > 1)        -> the uv_rhsAB copies of the shared (border) elements stay
    !                          bitwise equal with the upwind blend on (shared_copy_guard).
    !
    ! MESHES (hybrid meshes, Task 10): FESOM3_MESH_DIR (default pi, rotated 50/15/-90;
    ! FESOM3_TEST_ROTATE=0 reads it unrotated, e.g. the NEMO-mixed cut-out), or, with
    ! FESOM3_TEST_PATCH=mixed at np 1, the analytic mixed triangle/quad patch
    ! (generate_mixed_patch). All element loops run over the 4 vertex slots (a triangle's
    ! slot 4 repeats slot 1 with kite weight 0); element values of the analytic fields
    ! are taken at the element CENTROID (area-weighted, as mod_mesh_areas), which is
    ! where the element velocity of a linear field must sit for V2 to be exact.
    !
    ! EXACTNESS ON QUADS. The median-dual circulation with an element-constant velocity
    ! u_e = u(centroid) is exact for linear fields on triangles and PARALLELOGRAMS only:
    ! per element the residual against the exact contour integral telescopes around the
    ! node iff the centroid is the midpoint of the diagonal opposite the node (always
    ! true on a triangle, true for both diagonals only on a parallelogram). So V2 and V9a
    ! assert round-off only at "exact" nodes (every adjacent element a triangle or a
    ! parallelogram: pi, the NEMO cut-outs -- lon/lat rectangles -- and the parallelogram
    ! patch); at the other interior nodes (the jittered mixed patch) V2 asserts first-order
    ! consistency, |zeta/(-alpha) - 1| < 0.1 (measured 6e-2 with 15% jitter).
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
    use mod_mesh,          only: t_mesh, MAX_NV
    use mod_dyn,           only: t_dyn
    use mod_partit,        only: t_partit
    use mod_partitioning,  only: par_init, par_ex, set_partition
    use mod_param_phys,    only: rvo_upwind
    use mod_mesh_rotate,   only: trim_cyclic
    use mod_mesh_read,     only: read_mesh
    use mod_mesh_areas,    only: compute_geometry
    use mod_mesh_analytic, only: generate_mixed_patch, generate_quad_patch
    use mod_part_bounds,   only: owned_bounds
    use oce_dyn_vinv,      only: relative_vorticity, momentum_adv_vinv
    use mod_halo,          only: exchange_elem, shared_copy_max_elem
    implicit none

    character(len=512) :: mesh_dir, cbuf
    type(t_partit) :: partit
    type(t_mesh)   :: mesh
    type(t_dyn)    :: dyn
    integer :: nfail, nsw, ierr, nex, ngen, nv9
    logical :: jittered
    integer :: nNodO, nNodL, nEdgeO, nElemO, nElemF, nl, n, nz
    real(kind=WP) :: zmax, zmax2
    real(kind=WP), allocatable :: zsave(:,:)
    logical,       allocatable :: interior(:), elem_int(:), elem_int2(:), exact(:)
    real(kind=WP) :: scal
    integer       :: nzall
    real(kind=WP), allocatable :: rhsA(:,:,:,:)
    real(kind=WP), parameter   :: ALPHA = 1.0e-5_WP   ! du/dy [1/s]

    nfail = 0
    call get_environment_variable('FESOM3_MESH_DIR', mesh_dir)
    if (len_trim(mesh_dir) == 0) &
        mesh_dir = '/home/a/a270088/port2/fesom2/tests/data/MESHES/pi'
    call par_init(partit)

    cbuf = ''
    call get_environment_variable('FESOM3_TEST_PATCH', cbuf)
    jittered = trim(cbuf) == 'mixed'      ! no exact interior node by design (see header)
    if (trim(cbuf) == 'mixed') then
        if (partit%npes /= 1) error stop 'test_vinv: FESOM3_TEST_PATCH needs np 1'
        ! analytic mixed patch (cartesian geometry built by the generator)
        call generate_mixed_patch(mesh, partit, nx=13, ny=11, nl=6, &
                                  Lx=12000.0_WP, Ly=10000.0_WP, max_depth=100.0_WP)
        write(*,'(a,i0,a,i0)') '  analytic mixed patch: triangles=', count(mesh%elem2D_nnodes == 3), &
                               ' quads=', count(mesh%elem2D_nnodes == 4)
    else if (trim(cbuf) == 'quad') then
        if (partit%npes /= 1) error stop 'test_vinv: FESOM3_TEST_PATCH needs np 1'
        ! analytic all-quad patch of (skewed) parallelograms
        call generate_quad_patch(mesh, partit, nx=11, ny=9, nl=6, Lx=10000.0_WP, &
                                 Ly=8000.0_WP, max_depth=100.0_WP, skew=0.35_WP)
        write(*,'(a,i0)') '  analytic parallelogram patch: quads=', count(mesh%elem2D_nnodes == 4)
    else
        if (partit%npes > 1) call set_partition(partit, trim(mesh_dir))
        cbuf = ''
        call get_environment_variable('FESOM3_TEST_ROTATE', cbuf)
        if (trim(cbuf) == '0') then
            call read_mesh(mesh, partit, trim(mesh_dir), 0.0_WP, 0.0_WP, 0.0_WP, 360.0_WP, &
                           force_rotation=.false., n_cw_swaps=nsw)
        else
            call read_mesh(mesh, partit, trim(mesh_dir), 50.0_WP, 15.0_WP, -90.0_WP, 360.0_WP, &
                           force_rotation=.true., n_cw_swaps=nsw)
        end if
        ! cartesian: makes the discrete curl exact for a linear field (see header)
        call compute_geometry(mesh, partit, cartesian=.true.)
    end if
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

    ! exact = every adjacent element a triangle or a parallelogram (see header); the
    ! elements around an owned node are owned, so their vertex lists are local
    allocate(exact(nNodL))
    exact = .false.
    do n = 1, nNodO
        exact(n) = .true.
        do nz = 1, mesh%nod_in_elem2D_num(n)
            if (.not. is_parallelogram(mesh%nod_in_elem2D(nz, n))) exact(n) = .false.
        end do
    end do
    write(*,'(a,i7,a,i7)') '  exact interior nodes (triangles/parallelograms only): ', &
        count(interior(1:nNodO) .and. exact(1:nNodO)), ' of ', count(interior(1:nNodO))

    !=========================================================================
    ! A - elem_neighbors adjacency (build_elem_adjacency, Task 1 of rvo_upwind)
    !=========================================================================
    call check_true('A1 elem_neighbors allocated+resolved', &
                    allocated(mesh%elem_neighbors) .and. all(mesh%elem_neighbors(:,1:nElemO) >= 1))
    block
        integer :: e, k, nb, nshared, i, j, nbnd_slots, nbnd_edges, nv
        logical :: mutual
        nbnd_slots = 0
        do e = 1, nElemO
            nv = mesh%elem2D_nnodes(e)
            ! a triangle's slot 4 repeats slot 1 (checked here, then skipped)
            if (nv == 3 .and. mesh%elem_neighbors(4, e) /= mesh%elem_neighbors(1, e)) &
                call check_true('A4 triangle neighbour slot 4 repeats slot 1', .false.)
            do k = 1, nv
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
                do i = 1, nv
                    do j = 1, mesh%elem2D_nnodes(nb)
                        if (mesh%elem2D_nodes(i, e) == mesh%elem2D_nodes(j, nb)) nshared = nshared + 1
                    end do
                end do
                if (nshared /= 2) call check_true('A4 neighbours share exactly 2 vertices', .false.)
                if (.not. (any(mesh%elem2D_nodes(:, nb) == mesh%elem2D_nodes(k, e)) .and. &
                           any(mesh%elem2D_nodes(:, nb) == mesh%elem2D_nodes(mod(k,nv)+1, e)))) &
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
            do k = 1, mesh%elem2D_nnodes(e)
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

    !-------------------------------------------------------------------------
    ! A6 - the node halo assumption of every element->node and node->element operator
    ! here: all vertices of an owned element (4 slots on a quad) lie in owned+halo
    ! (1..nNodL), never in the double halo, so vort/KE/w_e are readable there after the
    ! standard exchange_nod.
    write(*,'(a,i7,a,i7)') '  A6 max vertex id of owned elements = ', &
        maxval(mesh%elem2D_nodes(:, 1:nElemO)), '  nNodL = ', nNodL
    call check_true('A6 owned-element vertices inside owned+halo nodes', &
                    maxval(mesh%elem2D_nodes(:, 1:nElemO)) <= nNodL)

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

    zmax = 0.0_WP; zmax2 = 0.0_WP; nex = 0; ngen = 0
    do n = 1, nNodO
        if (.not. interior(n)) cycle
        do nz = mesh%ulevels_nod2D(n), mesh%nlevels_nod2D_min(n)-1
            if (exact(n)) then
                zmax = max(zmax, abs(dyn%work%vorticity(nz,n) + ALPHA)); nex = nex + 1
            else
                zmax2 = max(zmax2, abs(dyn%work%vorticity(nz,n) + ALPHA)); ngen = ngen + 1
            end if
        end do
    end do
    call max_all(zmax); call max_all(zmax2); call sum_int(nex); call sum_int(ngen)
    write(*,'(a,i0,a,i0)') '  V2 levels checked: exact nodes ', nex, ', general quads ', ngen
    ! non-vacuity: the round-off set on every mesh but the jittered patch, the first-order
    ! set there
    if (jittered) then
        call check_true('V2 general-quad levels checked (non-vacuity)', ngen > 0)
    else
        call check_true('V2 exact-node levels checked (non-vacuity)', nex > 0)
    end if
    write(*,'(a,es12.4,a,es12.4)') '  V2 linear shear   : max|zeta-(-alpha)| = ', zmax, &
                                   '   relative = ', zmax/ALPHA
    write(*,'(a,es12.4)') '  V2 (general quads): relative = ', zmax2/ALPHA
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
    call check_true('V2 u=alpha*y -> zeta == -alpha to round-off (interior, all-wet, exact nodes)', &
                    zmax/ALPHA < 1.0e-11_WP)
    call check_true('V2 u=alpha*y -> zeta == -alpha to first order at general quads', &
                    zmax2/ALPHA < 0.1_WP)

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
        elem_int(n) = all(interior(mesh%elem2D_nodes(:,n)))
    end do
    ! one ring deeper: the upwind blend reads the NEIGHBOURS' vertex vorticity, so a
    ! constant-field inertness assertion must exclude elements whose neighbours touch the
    ! boundary too (and, at np>1, elements with a halo neighbour whose vertices cannot be
    ! inspected locally -- elem2D_nodes is owned-only).
    allocate(elem_int2(nElemO))
    do n = 1, nElemO
        elem_int2(n) = elem_int(n)
        do nz = 1, MAX_NV
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

    ! Restricted to levels where ALL elements adjacent to ALL the element's nodes are wet. Below
    ! that, KE(nz,n) sums only over the elements wet at nz while dividing by the FULL
    ! area(n) -- the zero-with-full-weight rule -- so KE is deliberately not uniform there
    ! and grad(KE) /= 0 even for a uniform u. That is intended, not an error.
    zmax = 0.0_WP; scal = 0.0_WP
    do n = 1, nElemO
        if (.not. elem_int(n)) cycle
        nzall = minval(mesh%nlevels_nod2D_min(mesh%elem2D_nodes(:,n))) - 1
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
    ! of equal values re-forms the same value up to round-off). Elements whose own and
    ! neighbours' vertices are exact nodes only (zeta is constant only there), and
    ! normalised by the size of the vorticity term |u|*alpha*area: the full tendency of a
    ! pure shear is itself zero (zeta x u cancels grad KE) on a planar mesh.
    rvo_upwind = 1.0_WP
    dyn%uv_rhsAB = 0.0_WP
    call momentum_adv_vinv(dyn, mesh, partit)
    zmax = 0.0_WP; scal = 0.0_WP; nv9 = 0
    do n = 1, nElemO
        if (.not. elem_int2(n)) cycle
        if (.not. all(exact(mesh%elem2D_nodes(:,n)))) cycle
        if (.not. all([(all(exact(mesh%elem2D_nodes(:, mesh%elem_neighbors(nz,n)))), nz = 1, MAX_NV)])) cycle
        nv9 = nv9 + 1
        nzall = minval(mesh%nlevels_nod2D_min(mesh%elem2D_nodes(:,n))) - 1
        do nz = 1, MAX_NV
            nzall = min(nzall, minval(mesh%nlevels_nod2D_min( &
                        mesh%elem2D_nodes(:, mesh%elem_neighbors(nz,n)))) - 1)
        end do
        do nz = mesh%ulevels(n), min(mesh%nlevels(n)-1, nzall)
            zmax = max(zmax, abs(dyn%uv_rhsAB(1,1,nz,n)-rhsA(1,1,nz,n)), &
                             abs(dyn%uv_rhsAB(1,2,nz,n)-rhsA(1,2,nz,n)))
            scal = max(scal, sqrt(dyn%uv(1,nz,n)**2 + dyn%uv(2,nz,n)**2)*ALPHA &
                             *real(mesh%elem_area(n),WP))
        end do
    end do
    call max_all(zmax); call max_all(scal); call sum_int(nv9)
    write(*,'(a,es12.4,a,i0)') '  V9a const-zeta blend inertness: rel = ', zmax/max(scal,tiny(1.0_WP)), &
                               '  elements ', nv9
    if (.not. jittered) call check_true('V9a elements checked (non-vacuity)', nv9 > 0)
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
        nzall = minval(mesh%nlevels_nod2D_min(mesh%elem2D_nodes(:,n))) - 1
        do nz = 1, MAX_NV
            nzall = min(nzall, minval(mesh%nlevels_nod2D_min( &
                        mesh%elem2D_nodes(:, mesh%elem_neighbors(nz,n)))) - 1)
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
    ! V13 (np > 1): shared-element copies of uv_rhsAB stay equal with the blend on
    if (partit%npes > 1) call shared_copy_guard()
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
        real(kind=WP) :: nxk(MAX_NV), nyk(MAX_NV), dk(MAX_NV), wk(MAX_NV), cex, cey, cnx, cny
        real(kind=WP) :: latm, u_e, v_e, zb0, dzb, ea, sumw, zup, zlo, zhi, err, errmax
        real(kind=WP), allocatable :: rhsB(:,:,:,:)
        integer :: e, k, nb, nwet, nupw(1), nz1b, nbad10, none_cnt
        ! --- fields: u = CU*lat^2, v = 0 (rotated frame) ---
        do e = 1, nElemO
            call elem_centroid(e, real(mesh%coord_nod2D(1,mesh%elem2D_nodes(1,e)),WP), &
                               0.0_WP, cnx, latm)
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
            call elem_normals(e, nxk, nyk, cex, cey)
            u_e = dyn%uv(1,1,e); v_e = dyn%uv(2,1,e)
            ! ---- V10, velocity-free: the slot-k neighbour must lie on the OUTWARD side
            ! of the slot-k inward normal, (c_nb - c_e).n_in < 0, in the SAME metre
            ! metric as the normals. (The first draft asserted (c_nb-c_e).u < 0 for
            ! inflow edges -- too strong: u.n_in > 0 does not bound the TANGENTIAL part
            ! of the centroid offset, so ~20% of irregular-triangle edges legitimately
            ! violated it. The velocity-free form is the actual orientation property;
            ! d = u.n_in is then correct by construction.) ----
            do k = 1, mesh%elem2D_nnodes(e)
                nb = mesh%elem_neighbors(k, e)
                if (nb == e .or. nb > nElemO) cycle
                call elem_centroid(nb, cex, cey, cnx, cny)
                if ( ((cnx-cex)*real(mesh%elem_cos(e),WP)*r_earth)*nxk(k) &
                   + ((cny-cey)*r_earth)*nyk(k) >= 0.0_WP ) nbad10 = nbad10 + 1
            end do
            ! ---- V11 at the surface level of all-wet elements ----
            nz1b = mesh%ulevels(e)
            if (nz1b > minval(mesh%nlevels_nod2D_min(mesh%elem2D_nodes(:,e))) - 1) cycle
            if (u_e*u_e + v_e*v_e < 1.0e-12_WP) cycle
            zb0 = omega_of(e, nz1b)
            sumw = 0.0_WP; zup = 0.0_WP; nwet = 0
            zlo = zb0; zhi = zb0
            do k = 1, MAX_NV
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

    subroutine shared_copy_guard()
        ! V13: a shared (border) element is owned by every rank that owns one of its
        ! vertices, and only uv is synced to the designated owner's copy (update_vel,
        ! docs/plans/completed/20261012-shared-element-sync.md). With equal uv copies the uv_rhsAB
        ! copies must stay BITWISE equal by themselves: every input of momentum_adv_vinv is
        ! per element from the same global geometry, or a nodal field (vorticity, KE, w_e)
        ! that its owner computes and exchanges. V11's field (u = CU*lat^2 on the owned
        ! elements, identical on every copy) has varying vorticity, so the upwind blend
        ! (rvo_upwind = 1) selects neighbours -- the branchy path that turned the roundoff
        ! spread of the copies into O(1) in the CORE2 leak. Not vacuous: the blend must
        ! change uv_rhsAB on some shared element by more than 1e-6 relative.
        real(kind=WP), allocatable :: rhs0(:,:,:,:)
        real(kind=WP) :: spread, dmax, rmax
        integer :: p, k, e
        rvo_upwind = 0.0_WP
        dyn%uv_rhsAB = 0.0_WP
        call momentum_adv_vinv(dyn, mesh, partit)
        rhs0 = dyn%uv_rhsAB
        rvo_upwind = 1.0_WP
        dyn%uv_rhsAB = 0.0_WP
        call momentum_adv_vinv(dyn, mesh, partit)
        spread = shared_copy_max_elem(dyn%uv_rhsAB(1,:,:,:), partit)
        call max_all(spread)
        ! the shared elements of this rank: its receive slots and the elements it sends
        dmax = 0.0_WP; rmax = 0.0_WP
        do p = 1, 2
            do k = 1, merge(size(partit%com_elem2D_shared%rlist), size(partit%com_elem2D_shared%slist), p == 1)
                if (p == 1) then
                    e = partit%com_elem2D_shared%rlist(k)
                else
                    e = partit%com_elem2D_shared%slist(k)
                end if
                dmax = max(dmax, maxval(abs(dyn%uv_rhsAB(1,:,:,e) - rhs0(1,:,:,e))))
                rmax = max(rmax, maxval(abs(rhs0(1,:,:,e))))
            end do
        end do
        call max_all(dmax); call max_all(rmax)
        if (partit%mype == 0) write(*,'(a,es12.4,a,es12.4)') '  V13 shared-copy spread of uv_rhsAB = ', &
            spread, '   blend change on shared elements (rel) = ', dmax/max(rmax, tiny(1.0_WP))
        call check_true('V13 uv_rhsAB copies of shared elements bitwise equal', spread == 0.0_WP)
        call check_true('V13 blend active on shared elements (non-vacuity)', &
                        dmax > 1.0e-6_WP*rmax .and. rmax > 0.0_WP)
    end subroutine shared_copy_guard

    logical function is_parallelogram(e)
        ! a triangle, or a quad whose diagonals bisect each other (to 1e-9 of its size)
        integer, intent(in) :: e
        real(kind=WP) :: x(MAX_NV), y(MAX_NV)
        integer :: k
        is_parallelogram = .true.
        if (mesh%elem2D_nnodes(e) == 3) return
        do k = 1, MAX_NV
            x(k) = real(mesh%coord_nod2D(1,mesh%elem2D_nodes(k,e)),WP) - &
                   real(mesh%coord_nod2D(1,mesh%elem2D_nodes(1,e)),WP)
            call trim_cyclic(x(k))
            y(k) = real(mesh%coord_nod2D(2,mesh%elem2D_nodes(k,e)),WP) - &
                   real(mesh%coord_nod2D(2,mesh%elem2D_nodes(1,e)),WP)
        end do
        is_parallelogram = abs(x(1)+x(3)-x(2)-x(4)) + abs(y(1)+y(3)-y(2)-y(4)) &
                           <= 1.0e-9_WP*(maxval(abs(x)) + maxval(abs(y)))
    end function is_parallelogram

    real(kind=WP) function omega_of(e, nz)
        ! the kite-weighted vertex mean, as the model's face vorticity
        integer, intent(in) :: e, nz
        integer :: k
        omega_of = 0.0_WP
        do k = 1, MAX_NV
            omega_of = omega_of + real(mesh%elem_w(k,e),WP)*dyn%work%vorticity(nz, mesh%elem2D_nodes(k,e))
        end do
    end function omega_of

    subroutine elem_centroid(e, refx, refy, cx, cy)
        ! area-weighted centroid of the sub-triangles (1,3,4) and (1,2,3) (a triangle:
        ! the vertex mean), built from the vertex offsets to (refx,refy) -- the lon
        ! offsets trimmed, so it is cyclic-safe -- and returned as refx/refy + offset
        integer,       intent(in)  :: e
        real(kind=WP), intent(in)  :: refx, refy
        real(kind=WP), intent(out) :: cx, cy
        real(kind=WP) :: x(MAX_NV), y(MAX_NV), s1, s2
        integer :: k
        do k = 1, MAX_NV
            x(k) = real(mesh%coord_nod2D(1,mesh%elem2D_nodes(k,e)),WP) - refx
            call trim_cyclic(x(k))
            y(k) = real(mesh%coord_nod2D(2,mesh%elem2D_nodes(k,e)),WP) - refy
        end do
        s1 = abs((x(3)-x(4))*(y(1)-y(4)) - (y(3)-y(4))*(x(1)-x(4)))
        s2 = abs((x(3)-x(2))*(y(1)-y(2)) - (y(3)-y(2))*(x(1)-x(2)))
        cx = refx + ((x(1)+x(3)+x(4))*s1 + (x(1)+x(2)+x(3))*s2)/(3.0_WP*(s1+s2))
        cy = refy + ((y(1)+y(3)+y(4))*s1 + (y(1)+y(2)+y(3))*s2)/(3.0_WP*(s1+s2))
    end subroutine elem_centroid

    subroutine elem_normals(e, nxk, nyk, cex, cey)
        ! the INWARD side normals of element e (4 slots), metre-scaled with elem_cos(e);
        ! side k connects elnodes(k), elnodes(k+1) (slot 4 -> 1; a triangle's side 4 is
        ! degenerate, normal 0); inward = toward the centroid (cex, cey).
        ! Mirrors the model's convention; V10 then checks the convention against
        ! geometry, so it is proven rather than shared-by-copy.
        integer,       intent(in)  :: e
        real(kind=WP), intent(out) :: nxk(MAX_NV), nyk(MAX_NV), cex, cey
        real(kind=WP) :: x(MAX_NV), y(MAX_NV), tx, ty, px, py, dref
        integer :: k, k2
        do k = 1, MAX_NV
            x(k) = real(mesh%coord_nod2D(1, mesh%elem2D_nodes(k,e)), WP)
            y(k) = real(mesh%coord_nod2D(2, mesh%elem2D_nodes(k,e)), WP)
        end do
        ! cyclic-safe local frame about vertex 1
        do k = 2, MAX_NV
            tx = x(k) - x(1); call trim_cyclic(tx); x(k) = x(1) + tx
        end do
        call elem_centroid(e, x(1), 0.0_WP, cex, cey)
        do k = 1, MAX_NV
            k2 = mod(k,MAX_NV)+1
            tx = (x(k2)-x(k))*real(mesh%elem_cos(e),WP)*r_earth
            ty = (y(k2)-y(k))*r_earth
            px =  ty; py = -tx
            dref = px*(cex-x(k))*real(mesh%elem_cos(e),WP)*r_earth + py*(cey-y(k))*r_earth
            if (dref < 0.0_WP) then
                px = -px; py = -py
            end if
            nxk(k) = px; nyk(k) = py
        end do
    end subroutine elem_normals

    subroutine set_shear(sgn)
        ! u = sgn*ALPHA*y, v = 0 with y = r_earth*lat : LINEAR in the cartesian metric, so
        ! the discrete curl is exact and zeta == -sgn*ALPHA. Depends on latitude only, so
        ! it is continuous across pi's cyclic seam. Element value at the element centroid
        ! (a triangle: the mean of its 3 corner latitudes).
        !
        ! coord_nod2D, NOT geo_coord_nod2D: the mesh is read with force_rotation=.true., so
        ! the geometry (elem_area, edge_cross_dxdy, element centres) lives in the ROTATED
        ! computational frame. A field linear in GEOGRAPHIC latitude is not linear in that
        ! frame, and the test then fails worst near the rotated pole -- which is exactly how
        ! this was found (mean zeta/(-alpha) = 0.836, worst at rotated lat 84 deg).
        real(kind=WP), intent(in) :: sgn
        integer :: el
        real(kind=WP) :: latm, lonm
        ! 1..nElemO, NOT nElemF: elem2D_nodes is owned-only, and both triangles of an
        ! owned edge are owned (LESSONS L31), so the edge loop never reads a halo element.
        do el = 1, nElemO
            call elem_centroid(el, real(mesh%coord_nod2D(1,mesh%elem2D_nodes(1,el)),WP), &
                               0.0_WP, lonm, latm)
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
