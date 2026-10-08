module oce_muscl_adv
    ! MUSCL-type advection setup, transcribed from FESOM2 v2.7.3 oce_muscl_adv.F90.
    !   muscl_adv_init            -> calls find_up_downwind_elements
    !   find_up_downwind_elements -> edge_up_dn_tri (upwind/downwind element per edge)
    !
    ! Reference: Abalakin, Dervieux, Kozubskaya (2002), INRIA RR-4459; the concept
    ! of upwind/downwind triangles to a given edge (sergey.danilov@awi.de 2012).
    ! Hybrid meshes (triangles + quads): the search runs over elements as in the qq
    ! prototype (qq/oce_muscl_adv.F90 find_up_downwind_elements): at the anchor vertex the
    ! candidate sector is the angle spanned by the element's two sides to the NEXT and the
    ! PREVIOUS vertex (cyclic over the nv unique vertices; nv = 3 when slot 4 repeats
    ! slot 1). The edge-array name edge_up_dn_tri is kept (t_tracer_work, binary layout).
    !
    ! SCOPE / clean-architecture deviations (D7, plan "USE-globals -> explicit args"):
    !  * 1-rank: coord_elem(:,n,el)==coord_nod2D(:,elem2D_nodes(n,el)) and the FESOM2
    !    global id myList_nod2D(node)==node, so the e_nodes(n,el)==myList_nod2D(...) tests
    !    become elem2D_nodes(n,el)==node.
    !  * muscl_adv_init's nn_num/nn_pos block (FESOM2 oce_muscl_adv.F90:59-90,106-127)
    !    is OMITTED: those node-neighbour arrays size off SSH_stiff (not built until
    !    M2.6) and are used only by quadratic reconstruction, NOT by horizontal
    !    advection. FESOM2's nboundary_lay loop is also dropped: it fed only the clamp of
    !    FESOM2's adv_tra_hor_muscl, which FESOM3 does not port (single unclamped kernel).
    !  * FESOM2's fill_up_dn_grad (per-edge stored up/dn gradient) is not ported: the MUSCL
    !    kernel takes the up/downwind element gradient on the fly (oce_adv_tra_hor).
    use mod_precision,   only: WP
    use mod_mesh,        only: t_mesh, MAX_NV
    use mod_partit,      only: t_partit
    use mod_tracer,      only: t_tracer_work
    use mod_mesh_rotate, only: get_cyclic_length
    use mod_part_bounds, only: owned_bounds, is_multirank
    use mod_halo,        only: exchange_elem_full
    implicit none
    private
    public :: muscl_adv_init, find_up_downwind_elements

    ! Angle tolerance [rad] of the "x along a side" test. FESOM2 compared exactly
    ! (ab == ax, ax == 0); on structured quad grids the edge direction is often exactly
    ! collinear with a side of the neighbouring elements, where an exact compare can miss
    ! both candidates by roundoff. qq uses 0.01; a roundoff-sized value keeps the triangle
    ! choice unchanged away from exact collinearity.
    real(kind=WP), parameter :: ang_tol = 1.0e-10_WP

contains

    !---------------------------------------------------------------------------
    subroutine muscl_adv_init(twork, mesh, partit)
        ! oce_muscl_adv.F90:32-158 (up/downwind elements only — see module header).
        type(t_mesh),        intent(in)    :: mesh
        type(t_tracer_work), intent(inout) :: twork
        type(t_partit),      intent(in), optional :: partit
        ! find upwind and downwind element for each local edge
        call find_up_downwind_elements(twork, mesh, partit)
    end subroutine muscl_adv_init


    !---------------------------------------------------------------------------
    subroutine find_up_downwind_elements(twork, mesh, partit)
        ! oce_muscl_adv.F90:162-352 (qq oce_muscl_adv.F90:56-255). For each edge, find the
        ! element around node 1 whose vertex angle contains the direction -edge (upwind),
        ! and around node 2 the element containing +edge (downwind). Decompose b (the side
        ! to the next vertex) and the search direction x along c (the side to the previous
        ! vertex) and the 90deg-CCW normal, then compare the atan2 angles.
        !
        ! M2.12b (optional partit): a halo node's element list (nod_in_elem2D) reaches
        ! eXDim elements whose elem2D_nodes are NOT stored locally and whose nodes may
        ! be beyond the local node list. So at npes>1 we build, exactly as FESOM2
        ! (oce_muscl_adv.F90:192-235), coord_elem(2,4,*) and e_nodes(4,*) over the FULL
        ! element halo via exchange_elem_full, and pick the element sides from those +
        ! the GLOBAL anchor id (myList_nod2D). At 1-rank (partit absent) the sides come
        ! from elem2D_nodes/coord_nod2D.
        type(t_mesh),        intent(in)    :: mesh
        type(t_tracer_work), intent(inout) :: twork
        type(t_partit),      intent(in), optional :: partit
        integer       :: n, k, ednodes(2), elem
        real(kind=WP) :: x(2), b(2), c(2), cr, bx, by, xx, xy, ab, ax, cl
        integer       :: nNodO, nNodL, nEdgeO, nElemO, nElemF, kk, nn
        logical       :: lmr
        real(kind=WP), allocatable :: coord_elem(:,:,:), tmp(:)
        integer,       allocatable :: e_nodes(:,:), tmp_i(:)

        lmr = is_multirank(partit)
        call owned_bounds(mesh, nNodO, nNodL, nEdgeO, nElemO, partit)

        if (.not. allocated(twork%edge_up_dn_tri))  allocate(twork%edge_up_dn_tri(2, nEdgeO))
        twork%edge_up_dn_tri = 0
        cl = get_cyclic_length()

        if (lmr) then
            ! element vertex coords + GLOBAL node ids over the full element halo
            nElemF = partit%myDim_elem2D + partit%eDim_elem2D + partit%eXDim_elem2D
            allocate(coord_elem(2, MAX_NV, nElemF), e_nodes(MAX_NV, nElemF))
            allocate(tmp(nElemF), tmp_i(nElemF))
            do nn = 1, MAX_NV
                do kk = 1, 2
                    tmp = 0.0_WP
                    do elem = 1, nElemO
                        tmp(elem) = mesh%coord_nod2D(kk, mesh%elem2D_nodes(nn, elem))
                    end do
                    call exchange_elem_full(tmp, partit)
                    coord_elem(kk, nn, :) = tmp
                end do
                tmp_i = 0
                do elem = 1, nElemO
                    tmp_i(elem) = partit%myList_nod2D(mesh%elem2D_nodes(nn, elem))
                end do
                call exchange_elem_full(tmp_i, partit)
                e_nodes(nn, :) = tmp_i
            end do
            deallocate(tmp, tmp_i)
        end if

        do n = 1, nEdgeO
            ednodes = mesh%edges(:, n)
            x = mesh%coord_nod2D(:, ednodes(2)) - mesh%coord_nod2D(:, ednodes(1))
            if (x(1) >  cl/2.0_WP) x(1) = x(1) - cl
            if (x(1) < -cl/2.0_WP) x(1) = x(1) + cl

            ! Find upwind (in the sense of x) element, i.e. which contains -x:
            x = -x
            do k = 1, mesh%nod_in_elem2D_num(ednodes(1))
                elem = mesh%nod_in_elem2D(k, ednodes(1))
                call get_bc(elem, ednodes(1), b, c)
                if (b(1) >  cl/2.0_WP) b(1) = b(1) - cl
                if (b(1) < -cl/2.0_WP) b(1) = b(1) + cl
                if (c(1) >  cl/2.0_WP) c(1) = c(1) - cl
                if (c(1) < -cl/2.0_WP) c(1) = c(1) + cl
                ! Decompose b and x into parts along c and along (-cy,cx) (90deg CCW)
                cr = sum(c*c)
                bx = sum(b*c)/cr
                by = (-b(1)*c(2)+b(2)*c(1))/cr
                xx = sum(x*c)/cr
                xy = (-x(1)*c(2)+x(2)*c(1))/cr
                ab = atan2(by, bx)
                ax = atan2(xy, xx)
                ! Since b,c are sides of a convex element, |ab|<pi, so atan2 is what is needed
                if ((ab > 0.0_WP) .and. (ax > 0.0_WP) .and. (ax < ab)) then
                    twork%edge_up_dn_tri(1, n) = elem; cycle
                end if
                if ((ab < 0.0_WP) .and. (ax < 0.0_WP) .and. (ax > ab)) then
                    twork%edge_up_dn_tri(1, n) = elem; cycle
                end if
                if ((abs(ab-ax) < ang_tol) .or. (abs(ax) < ang_tol)) then
                    twork%edge_up_dn_tri(1, n) = elem; cycle
                end if
            end do

            ! Find downwind element (direction +x)
            x = -x
            do k = 1, mesh%nod_in_elem2D_num(ednodes(2))
                elem = mesh%nod_in_elem2D(k, ednodes(2))
                call get_bc(elem, ednodes(2), b, c)
                if (b(1) >  cl/2.0_WP) b(1) = b(1) - cl
                if (b(1) < -cl/2.0_WP) b(1) = b(1) + cl
                if (c(1) >  cl/2.0_WP) c(1) = c(1) - cl
                if (c(1) < -cl/2.0_WP) c(1) = c(1) + cl
                cr = sum(c*c)
                bx = sum(b*c)/cr
                by = (-b(1)*c(2)+b(2)*c(1))/cr
                xx = sum(x*c)/cr
                xy = (-x(1)*c(2)+x(2)*c(1))/cr
                ab = atan2(by, bx)
                ax = atan2(xy, xx)
                if ((ab > 0.0_WP) .and. (ax > 0.0_WP) .and. (ax < ab)) then
                    twork%edge_up_dn_tri(2, n) = elem; cycle
                end if
                if ((ab < 0.0_WP) .and. (ax < 0.0_WP) .and. (ax > ab)) then
                    twork%edge_up_dn_tri(2, n) = elem; cycle
                end if
                if ((abs(ab-ax) < ang_tol) .or. (abs(ax) < ang_tol)) then
                    twork%edge_up_dn_tri(2, n) = elem; cycle
                end if
            end do
        end do

        if (lmr) deallocate(coord_elem, e_nodes)

        ! For edges touching the boundary, up/downwind elements may be absent; the
        ! MUSCL kernels then use a zero reconstruction increment on that side.

    contains
        subroutine get_bc(el, anchor_loc, bb, cc)
            ! Element sides from vertex `anchor_loc`. 1-rank: via elem2D_nodes/
            ! coord_nod2D. Multi-rank: via coord_elem + the GLOBAL anchor id matched
            ! against e_nodes (FESOM2 oce_muscl_adv.F90:251-303, qq :110-136).
            integer,       intent(in)  :: el, anchor_loc
            real(kind=WP), intent(out) :: bb(2), cc(2)
            real(kind=WP) :: xy4(2, MAX_NV)
            integer :: j
            if (.not. lmr) then
                do j = 1, MAX_NV
                    xy4(:, j) = mesh%coord_nod2D(:, mesh%elem2D_nodes(j, el))
                end do
                call elem_sides(mesh%elem2D_nodes(:, el), xy4, anchor_loc, bb, cc)
            else
                call elem_sides(e_nodes(:, el), coord_elem(:, :, el), partit%myList_nod2D(anchor_loc), bb, cc)
            end if
        end subroutine get_bc
    end subroutine find_up_downwind_elements

    pure subroutine elem_sides(ids, xy4, anchor, b, c)
        ! The two element sides emanating from the vertex `anchor`: b to the NEXT vertex,
        ! c to the PREVIOUS one, cyclic over the nv unique vertices (nv = 3 when slot 4
        ! repeats slot 1, the triangle convention). ids are the 4 vertex ids (local or
        ! global, matching `anchor`), xy4 their coordinates. For a triangle this is FESOM2's
        ! b/c at vertices 1 and 3 and b/c swapped at vertex 2 (the same sector).
        integer,       intent(in)  :: ids(MAX_NV), anchor
        real(kind=WP), intent(in)  :: xy4(2, MAX_NV)
        real(kind=WP), intent(out) :: b(2), c(2)
        integer :: nv, j, ja, jn, jp
        nv = MAX_NV
        if (ids(MAX_NV) == ids(1)) nv = 3
        ja = nv                               ! as FESOM2: not the first ones -> the last
        do j = 1, nv - 1
            if (ids(j) == anchor) then
                ja = j
                exit
            end if
        end do
        jn = mod(ja, nv) + 1
        jp = mod(ja + nv - 2, nv) + 1
        b = xy4(:, jn) - xy4(:, ja)
        c = xy4(:, jp) - xy4(:, ja)
    end subroutine elem_sides

end module oce_muscl_adv
