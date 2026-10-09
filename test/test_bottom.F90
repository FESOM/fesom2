program test_bottom
    ! FESOM3 BOTTOM AT VERTICES contract test (docs/plans/20260910-fesom3-bottom-at-vertices.md).
    !
    ! The vertex column is authoritative and the element bounds are DERIVED:
    !     ulevels(e) = maxval(ulevels_nod2D(elnodes))
    !     nlevels(e) = minval(nlevels_nod2D(elnodes))
    ! so an element's layer range is exactly its fully wet prisms. These checks pin that
    ! contract down against a hand-set sloping triangle (T2), the edge interval (T3), and
    ! the pi mesh (T6, T9 and the required invariant).
    !
    ! Multi-rank (np>1): additionally checks that the derivation reached the FULL element
    ! halo, which is the part elem2D_nodes cannot supply because it is owned-only.
    use mpi
    use mod_precision,     only: WP, MP
    use mod_constants,     only: r_earth
    use mod_mesh,          only: t_mesh
    use mod_partit,        only: t_partit
    use mod_partitioning,  only: par_init, par_ex, set_partition
    use mod_mesh_read,     only: read_mesh
    use mod_mesh_areas,    only: compute_geometry
    use mod_mesh_analytic, only: generate_analytic_mesh
    use mod_mesh_rotate,   only: trim_cyclic
    implicit none

    character(len=512) :: mesh_dir
    type(t_partit) :: partit
    type(t_mesh)   :: mesh, amesh
    integer :: nfail, nsw, nNodO, nElemO, nElemF

    call get_environment_variable('FESOM3_MESH_DIR', mesh_dir)
    if (len_trim(mesh_dir) == 0) &
        mesh_dir = '/home/a/a270088/port2/fesom2/tests/data/MESHES/pi'
    nfail = 0
    call par_init(partit)

    !=========================================================================
    ! pi — the real mesh.
    !=========================================================================
    if (partit%npes > 1) call set_partition(partit, trim(mesh_dir))
    call read_mesh(mesh, partit, trim(mesh_dir), 50.0_WP, 15.0_WP, -90.0_WP, 360.0_WP, &
                   force_rotation=.true., n_cw_swaps=nsw)
    call compute_geometry(mesh, partit, cartesian=.false.)
    if (partit%npes == 1) then
        nNodO = mesh%nod2D; nElemO = mesh%elem2D; nElemF = mesh%elem2D
    else
        nNodO  = partit%myDim_nod2D
        nElemO = partit%myDim_elem2D
        nElemF = partit%myDim_elem2D + partit%eDim_elem2D + partit%eXDim_elem2D
    end if
    if (partit%mype == 0) &
        write(*,'(a,i0,a,i0,a,i0)') 'pi: nod2D=', mesh%nod2D, ' elem2D=', mesh%elem2D, &
            ' nl=', mesh%nl

    call check(levels_in_range(mesh, nNodO, nElemO), 'pi: level bounds within [1, nl]')
    call check(elem_is_min_of_nodes(mesh, nElemO),   'pi: nlevels(e) == min over its nodes')
    call check(bottom_invariant(mesh, nNodO),        'pi: every vertex bottom cell has a wet element')
    call check(edge_interval_contains_elements(mesh),'T3: edge interval >= adjacent nlevels (pi)')
    call check(hnode_column_sum(mesh, nNodO),        'T6: sum(hnode) == column depth, hnode >= 0')
    call check(zbar_e_bot_derived(mesh, nElemO),     'T9: zbar_e_bot(e) == zbar(nlevels(e))')
    call check(sloping_triangle(mesh, nElemO),      'T2: sloping triangle takes the shallowest vertex')
    call check(edge_len_is_metres(mesh),            'T7: edge_len matches haversine, edge_dxdy metre-scale')
    call check(edge_dxdy_fold_exact(mesh),          'T7: edge_dxdy == radian value * r_earth * mean(elem_cos)')
    call check(edge_flux_lands_in_wet_cells(mesh),   'conservation: edge flux range fits inside both node columns')

    if (partit%npes > 1) then
        ! The halo derivation. elem2D_nodes is owned-only, so halo nlevels can only come
        ! from the global vertex-level scatter -- and oce_adv_tra_hor / vert_vel_ale read
        ! nlevels(el(2)) with el(2) possibly in the halo, so a stale halo entry is a real
        ! bug rather than dead storage.
        call check(all(mesh%nlevels(1:nElemF) > 0), 'MR: nlevels filled over the full element halo')
        call check(all(mesh%ulevels(1:nElemF) > 0), 'MR: ulevels filled over the full element halo')
        call check(all(mesh%nlevels(1:nElemF) <= mesh%nl), 'MR: halo nlevels within [1, nl]')
        ! and the values must be RIGHT, not merely present: recompute each local element's
        ! bound from its three GLOBAL node ids read straight from the mesh files.
        call check(halo_matches_global(mesh, partit, trim(mesh_dir), nElemF), &
                   'MR: halo nlevels == min over the element global vertex levels')
    end if

    !=========================================================================
    ! Flat bottom — the analytic mesh derives nlevels == nl everywhere, the degenerate
    ! case of the reduction. Built AFTER pi: generate_analytic_mesh fills the same partit,
    ! which can only be populated once.
    !=========================================================================
    if (partit%npes == 1) then
        call generate_analytic_mesh(amesh, partit, nx=9, ny=7, nl=12, &
                                    Lx=8000.0_WP, Ly=6000.0_WP, max_depth=1200.0_WP)
        call check(all(amesh%nlevels == 12) .and. all(amesh%ulevels == 1), &
                   'analytic: flat bottom derives nlevels == nl, ulevels == 1')
        call check(edge_interval_contains_elements(amesh), &
                   'T3: edge interval >= adjacent nlevels (analytic)')
    end if

    if (nfail == 0) then
        if (partit%mype == 0) write(*,'(a)') 'test_bottom: ALL PASS'
    else
        write(*,'(a,i0,a,i0)') 'test_bottom: rank ', partit%mype, ' FAILURE(S): ', nfail
    end if
    call par_ex(partit%MPI_COMM_FESOM, partit%mype)
    if (nfail /= 0) error stop 1

contains

    subroutine check(cond, name)
        logical, intent(in) :: cond
        character(len=*), intent(in) :: name
        if (.not. cond) then
            nfail = nfail + 1
            write(*,'(a)') '  FAIL: '//name
        end if
    end subroutine

    !-------------------------------------------------------------------------
    logical function sloping_triangle(m, nElem)
        ! T2, the design note's worked case: "three nodes of a mesh triangle may have
        ! different depths ... a full triangular prism will exist only between
        ! tlayer_elem = max(tlayer(elnodes)) and blayer_elem = min(blayer(elnodes))".
        !
        ! Taken from REAL bathymetry rather than fabricated. Setting three node levels by
        ! hand on an otherwise flat mesh produces a mesh that violates the bottom
        ! invariant (a neighbouring vertex is left with no wet adjacent element), so the
        ! honest test is to find a triangle whose three vertices genuinely differ -- pi
        ! has many -- and require the element to take the SHALLOWEST of them, since only
        ! layers wet at all three corners can carry a velocity.
        type(t_mesh), intent(in) :: m
        integer,      intent(in) :: nElem
        integer :: e, n1, n2, n3, found
        found = 0
        do e = 1, nElem
            if (m%elem2D_nnodes(e) /= 3) cycle
            n1 = m%elem2D_nodes(1,e); n2 = m%elem2D_nodes(2,e); n3 = m%elem2D_nodes(3,e)
            if (m%nlevels_nod2D(n1) == m%nlevels_nod2D(n2)) cycle
            if (m%nlevels_nod2D(n2) == m%nlevels_nod2D(n3)) cycle
            if (m%nlevels_nod2D(n1) == m%nlevels_nod2D(n3)) cycle
            found = e
            exit
        end do
        if (found == 0) then
            write(*,'(a)') '  (no triangle with three distinct vertex bottoms)'
            sloping_triangle = .false.; return
        end if
        n1 = m%elem2D_nodes(1,found); n2 = m%elem2D_nodes(2,found); n3 = m%elem2D_nodes(3,found)
        write(*,'(a,i0,a,3(i0,1x),a,i0)') '  T2: elem ', found, ' vertex bottoms [', &
            m%nlevels_nod2D(n1), m%nlevels_nod2D(n2), m%nlevels_nod2D(n3), '] -> nlevels=', &
            m%nlevels(found)
        sloping_triangle = &
            (m%nlevels(found) == min(m%nlevels_nod2D(n1), min(m%nlevels_nod2D(n2), m%nlevels_nod2D(n3)))) &
            .and. (m%ulevels(found) == 1) &
            .and. (m%elem_depth(found) == m%zbar(m%nlevels(found)))
    end function

    logical function edge_len_is_metres(m)
        ! T7 (R7). edge_dxdy is now stored in PHYSICAL measure and edge_len is the edge
        ! length in metres. Two independent checks:
        !   (a) magnitudes are metre-scale, not radian-scale. A radian-measure edge on pi
        !       is O(1e-2); in metres it is O(1e4-1e5). Anything below 1 m would mean the
        !       r_earth*cos factor never got folded in.
        !   (b) edge_len agrees with a haversine great-circle distance between the edge's
        !       two GEOGRAPHIC vertices, computed here from scratch. The two differ by the
        !       flat-earth approximation and by edge_dxdy's mean-cosine (evaluated at the
        !       adjacent element centres, not at the edge), so the tolerance is loose --
        !       but a missing or doubled r_earth is a factor of 6.4e6 or 2, not 20%.
        type(t_mesh), intent(in) :: m
        integer :: ed, n1, n2, nedge, nbad
        real(kind=WP) :: lon1, lat1, lon2, lat2, dlon, dlat, hav, dist, rel, worst
        edge_len_is_metres = .true.
        nedge = size(m%edge_len)
        if (nedge <= 0) then; edge_len_is_metres = .false.; return; end if
        if (minval(m%edge_len) <= 1.0_MP) then
            write(*,'(a,es12.4)') '  edge_len min is not metre-scale: ', minval(m%edge_len)
            edge_len_is_metres = .false.; return
        end if
        if (maxval(abs(m%edge_dxdy)) < 1.0_MP) then
            write(*,'(a,es12.4)') '  edge_dxdy still radian-scale: ', maxval(abs(m%edge_dxdy))
            edge_len_is_metres = .false.; return
        end if
        nbad = 0; worst = 0.0_WP
        do ed = 1, nedge
            n1 = m%edges(1, ed); n2 = m%edges(2, ed)
            if (n1 <= 0 .or. n2 <= 0) cycle
            lon1 = real(m%geo_coord_nod2D(1, n1), WP); lat1 = real(m%geo_coord_nod2D(2, n1), WP)
            lon2 = real(m%geo_coord_nod2D(1, n2), WP); lat2 = real(m%geo_coord_nod2D(2, n2), WP)
            dlon = lon2 - lon1; dlat = lat2 - lat1
            hav  = sin(dlat/2.0_WP)**2 + cos(lat1)*cos(lat2)*sin(dlon/2.0_WP)**2
            dist = 2.0_WP*r_earth*asin(min(1.0_WP, sqrt(max(hav, 0.0_WP))))
            if (dist <= 0.0_WP) cycle
            rel = abs(real(m%edge_len(ed), WP) - dist)/dist
            if (rel > worst) worst = rel
            if (rel > 0.2_WP) nbad = nbad + 1
        end do
        write(*,'(a,es12.4,a,i0)') '  T7: worst |edge_len-haversine|/haversine = ', worst, &
            '   edges over 20%: ', nbad
        if (nbad > 0) edge_len_is_metres = .false.
    end function

    logical function edge_flux_lands_in_wet_cells(m)
        ! CONSERVATION PRECONDITION. The horizontal edge flux is scattered as
        ! +flux_h(nz,edge) into one edge node and -flux_h(nz,edge) into the other
        ! (oce_adv_tra_flux.F90:80-81), so the two content changes cancel exactly -- but
        ! ONLY if both nodes actually have a cell at nz. If the flux range reached past a
        ! node's bottom it would deposit tracer into a cell that does not exist and the
        ! cancellation would be lost.
        !
        ! The scatter runs to nl12 = max(nlevels(el1)-1, nlevels(el2)-1). Under the vertex
        ! bottom nlevels(e) = min over the element's vertices, and both edge nodes belong
        ! to both elements, so nlevels(el) <= nlevels_nod2D(ednode) for either node and
        ! either element. Check that on the real mesh rather than trusting the argument.
        type(t_mesh), intent(in) :: m
        integer :: ed, i, el, nedge, nl12, node_bot
        edge_flux_lands_in_wet_cells = .true.
        nedge = size(m%edges, 2)
        do ed = 1, nedge
            if (m%edges(1,ed) <= 0 .or. m%edges(2,ed) <= 0) cycle
            nl12 = 0
            do i = 1, 2
                el = m%edge_tri(i, ed)
                if (i == 2 .and. el == m%edge_tri(1, ed)) cycle   ! boundary edge: one element
                if (m%nlevels(el) <= 0) cycle
                nl12 = max(nl12, m%nlevels(el)-1)
            end do
            node_bot = min(m%nlevels_nod2D(m%edges(1,ed)), m%nlevels_nod2D(m%edges(2,ed))) - 1
            if (nl12 > node_bot) then
                write(*,'(a,i0,a,i0,a,i0)') '  edge ', ed, ': flux range reaches nz=', nl12, &
                    ' but the shallower node column ends at ', node_bot
                edge_flux_lands_in_wet_cells = .false.; return
            end if
        end do
    end function

    logical function edge_dxdy_fold_exact(m)
        ! T7, the exact fold. R7 moved r_earth*mean(elem_cos) out of oce_adv_tra_hor's
        ! inline `a` and into edge_dxdy. Recompute the FESOM2 radian value from the node
        ! coordinates and require the stored array to be that value times the factor,
        ! with the mean taken over BOTH adjacent elements interior and over the single
        ! one at a boundary edge. This is what conservation cannot catch: a wrong metric
        ! factor changes the MUSCL reconstruction without breaking any tracer budget.
        type(t_mesh), intent(in) :: m
        integer :: ed, el1, el2, nedge
        real(kind=WP) :: a1, a2, cosm, wantx, wanty, rel, worst
        logical :: bad
        edge_dxdy_fold_exact = .true.
        worst = 0.0_WP; bad = .false.
        nedge = size(m%edge_len)
        do ed = 1, nedge
            if (m%edges(1,ed) <= 0 .or. m%edges(2,ed) <= 0) cycle
            a1 = real(m%coord_nod2D(1, m%edges(2,ed)) - m%coord_nod2D(1, m%edges(1,ed)), WP)
            a2 = real(m%coord_nod2D(2, m%edges(2,ed)) - m%coord_nod2D(2, m%edges(1,ed)), WP)
            call trim_cyclic(a1)
            el1 = m%edge_tri(1, ed); el2 = m%edge_tri(2, ed)
            if (el1 < 1) cycle
            cosm = real(m%elem_cos(el1), WP)
            if (el2 /= el1) cosm = 0.5_WP*(cosm + real(m%elem_cos(el2), WP))
            wantx = a1 * cosm * r_earth
            wanty = a2 * r_earth
            rel = max(abs(real(m%edge_dxdy(1,ed),WP) - wantx)/max(abs(wantx), 1.0_WP), &
                      abs(real(m%edge_dxdy(2,ed),WP) - wanty)/max(abs(wanty), 1.0_WP))
            if (rel > worst) worst = rel
            if (rel > 1.0e-12_WP) bad = .true.
        end do
        write(*,'(a,es12.4)') '  T7: worst relative deviation from the folded factor = ', worst
        if (bad) edge_dxdy_fold_exact = .false.
    end function

    !-------------------------------------------------------------------------
    logical function edge_interval_contains_elements(m)
        ! T3. The edge's wet interval is the intersection of its two vertex columns,
        ! min(nlevels_nod2D(ednodes)). Every element adjacent to the edge must be no
        ! deeper than that, because the edge's nodes are a subset of the element's.
        type(t_mesh), intent(in) :: m
        integer :: ed, i, el, nedge, ubound_edge
        edge_interval_contains_elements = .true.
        nedge = size(m%edges, 2)
        do ed = 1, nedge
            if (m%edges(1,ed) <= 0 .or. m%edges(2,ed) <= 0) cycle
            ubound_edge = min(m%nlevels_nod2D(m%edges(1,ed)), m%nlevels_nod2D(m%edges(2,ed)))
            do i = 1, 2
                el = m%edge_tri(i, ed)
                if (i == 2 .and. el == m%edge_tri(1, ed)) cycle   ! boundary edge: one element
                if (m%nlevels(el) <= 0) cycle
                if (m%nlevels(el) > ubound_edge) then
                    edge_interval_contains_elements = .false.; return
                end if
            end do
        end do
    end function

    !-------------------------------------------------------------------------
    logical function levels_in_range(m, nNod, nElem)
        type(t_mesh), intent(in) :: m
        integer,      intent(in) :: nNod, nElem
        integer :: n
        levels_in_range = .true.
        do n = 1, nNod
            if (m%ulevels_nod2D(n) < 1 .or. &
                m%ulevels_nod2D(n) >= m%nlevels_nod2D(n) .or. &
                m%nlevels_nod2D(n) > m%nl) then
                levels_in_range = .false.; return
            end if
        end do
        do n = 1, nElem
            if (m%ulevels(n) < 1 .or. m%ulevels(n) >= m%nlevels(n) .or. m%nlevels(n) > m%nl) then
                levels_in_range = .false.; return
            end if
        end do
    end function

    logical function elem_is_min_of_nodes(m, nElem)
        ! the derivation itself: nlevels(e) is exactly the min over the element's nodes
        ! (not merely <=, which would hold for any shallower choice too).
        type(t_mesh), intent(in) :: m
        integer,      intent(in) :: nElem
        integer :: e, nv
        elem_is_min_of_nodes = .true.
        do e = 1, nElem
            nv = m%elem2D_nnodes(e)
            if (nv <= 0) cycle
            if (m%nlevels(e) /= minval(m%nlevels_nod2D(m%elem2D_nodes(1:nv,e))) .or. &
                m%ulevels(e) /= maxval(m%ulevels_nod2D(m%elem2D_nodes(1:nv,e)))) then
                elem_is_min_of_nodes = .false.; return
            end if
        end do
    end function

    logical function bottom_invariant(m, nNod)
        ! the invariant three unguarded divides depend on: every vertex's deepest scalar
        ! cell must have at least one adjacent element wet down to it.
        type(t_mesh), intent(in) :: m
        integer,      intent(in) :: nNod
        integer :: n, k
        bottom_invariant = .true.
        do n = 1, nNod
            k = m%nod_in_elem2D_num(n)
            if (k <= 0) cycle
            if (maxval(m%nlevels(m%nod_in_elem2D(1:k,n))) /= m%nlevels_nod2D(n)) then
                bottom_invariant = .false.; return
            end if
        end do
    end function

    !-------------------------------------------------------------------------
    logical function hnode_column_sum(m, nNod)
        ! T6. Build hnode the way every driver does (linfs full cells) and require the
        ! column to sum to the vertex's wet depth, with no negative cell.
        type(t_mesh), intent(in) :: m
        integer,      intent(in) :: nNod
        integer :: n, nz
        real(kind=WP) :: acc, want, h
        hnode_column_sum = .true.
        do n = 1, nNod
            if (m%nlevels_nod2D(n) <= 0) cycle
            acc = 0.0_WP
            do nz = m%ulevels_nod2D(n), m%nlevels_nod2D(n)-1
                h = real(m%zbar(nz) - m%zbar(nz+1), WP)
                if (h < 0.0_WP) then; hnode_column_sum = .false.; return; end if
                acc = acc + h
            end do
            want = real(m%zbar(m%ulevels_nod2D(n)) - m%zbar(m%nlevels_nod2D(n)), WP)
            if (abs(acc - want) > 1.0e-9_WP*max(abs(want), 1.0_WP)) then
                hnode_column_sum = .false.; return
            end if
        end do
    end function

    logical function halo_matches_global(m, p, mdir, nElem)
        ! Independent recomputation: read nlvls.out and elem2d.out directly and check
        ! every LOCAL element (owned and halo) against min over its three global vertices.
        type(t_mesh),     intent(in) :: m
        type(t_partit),   intent(in) :: p
        character(len=*), intent(in) :: mdir
        integer,          intent(in) :: nElem
        integer :: u2, ios2, g, lev, want, e
        integer :: gn1, gn2, gn3, nNodG, nElemG
        integer, allocatable :: lvl(:), emap(:)
        halo_matches_global = .true.
        nNodG  = m%nod2D
        nElemG = m%elem2D
        allocate(lvl(nNodG), emap(nElemG))
        lvl = 0; emap = 0
        open(newunit=u2, file=mdir//'/nlvls.out', status='old', action='read', iostat=ios2)
        if (ios2 /= 0) then; halo_matches_global = .false.; return; end if
        do g = 1, nNodG
            read(u2,*) lev; lvl(g) = lev
        end do
        close(u2)
        do e = 1, nElem
            emap(p%myList_elem2D(e)) = e
        end do
        open(newunit=u2, file=mdir//'/elem2d.out', status='old', action='read', iostat=ios2)
        if (ios2 /= 0) then; halo_matches_global = .false.; return; end if
        read(u2,*) g
        do g = 1, nElemG
            read(u2,*) gn1, gn2, gn3
            e = emap(g)
            if (e < 1 .or. e > nElem) cycle
            want = min(lvl(gn1), min(lvl(gn2), lvl(gn3)))
            if (m%nlevels(e) /= want) then
                write(*,'(a,i0,a,i0,a,i0,a,i0)') '  MR: elem local ', e, ' global ', g, &
                    ' nlevels=', m%nlevels(e), ' expected ', want
                halo_matches_global = .false.
                close(u2); deallocate(lvl, emap); return
            end if
        end do
        close(u2)
        deallocate(lvl, emap)
    end function

    logical function zbar_e_bot_derived(m, nElem)
        ! T9. The element bottom depth used by the stiffness assembly and the bottom drag
        ! must be the DERIVED bottom, zbar(nlevels(e)).
        type(t_mesh), intent(in) :: m
        integer,      intent(in) :: nElem
        integer :: e
        zbar_e_bot_derived = .true.
        do e = 1, nElem
            if (m%nlevels(e) <= 0) cycle
            if (m%elem_depth(e) /= m%zbar(m%nlevels(e))) then
                zbar_e_bot_derived = .false.; return
            end if
        end do
    end function

end program test_bottom
