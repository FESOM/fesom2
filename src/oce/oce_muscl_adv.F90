module oce_muscl_adv
    ! MUSCL-type advection setup, transcribed from FESOM2 v2.7.3 oce_muscl_adv.F90.
    !   muscl_adv_init            -> calls find_up_downwind_triangles
    !   find_up_downwind_triangles-> edge_up_dn_tri (upwind/downwind tri per edge)
    !   fill_up_dn_grad           -> per-edge up/dn gradient (ORACLE only, see the routine)
    !
    ! Reference: Abalakin, Dervieux, Kozubskaya (2002), INRIA RR-4459; the concept
    ! of upwind/downwind triangles to a given edge (sergey.danilov@awi.de 2012).
    !
    ! SCOPE / clean-architecture deviations (D7, plan "USE-globals -> explicit args"):
    !  * 1-rank only (myDim_* == global; eDim==eXDim==0). The FESOM2 halo exchanges
    !    (exchange_elem of coord_elem/e_nodes in find_up_downwind_triangles) are
    !    no-ops at 1-rank and are dropped here; multi-rank is deferred to M1.5/M2.12.
    !    coord_elem(:,n,el)==coord_nod2D(:,elem2D_nodes(n,el)) and the FESOM2 global
    !    id myList_nod2D(node)==node, so the e_nodes(n,el)==myList_nod2D(...) tests
    !    become elem2D_nodes(n,el)==node.
    !  * muscl_adv_init's nn_num/nn_pos block (FESOM2 oce_muscl_adv.F90:59-90,106-127)
    !    is OMITTED: those node-neighbour arrays size off SSH_stiff (not built until
    !    M2.6) and are used only by quadratic reconstruction, NOT by horizontal
    !    advection. FESOM2's nboundary_lay loop is also dropped: it fed only the clamp of
    !    FESOM2's adv_tra_hor_muscl, which FESOM3 does not port (single unclamped kernel).
    !  * tr_xy is an explicit argument to fill_up_dn_grad (FESOM2 reads global o_ARRAYS).
    use mod_precision,   only: WP
    use mod_mesh,        only: t_mesh
    use mod_partit,      only: t_partit
    use mod_tracer,      only: t_tracer_work
    use mod_mesh_rotate, only: get_cyclic_length
    use mod_part_bounds, only: owned_bounds, is_multirank
    use mod_halo,        only: exchange_elem_full
    implicit none
    private
    public :: muscl_adv_init, find_up_downwind_triangles, fill_up_dn_grad

contains

    !---------------------------------------------------------------------------
    subroutine muscl_adv_init(twork, mesh, partit)
        ! oce_muscl_adv.F90:32-158 (up/downwind triangles only — see module header).
        type(t_mesh),        intent(in)    :: mesh
        type(t_tracer_work), intent(inout) :: twork
        type(t_partit),      intent(in), optional :: partit
        ! find upwind and downwind triangle for each local edge
        call find_up_downwind_triangles(twork, mesh, partit)
    end subroutine muscl_adv_init


    !---------------------------------------------------------------------------
    subroutine find_up_downwind_triangles(twork, mesh, partit)
        ! oce_muscl_adv.F90:162-352. For each edge, find the triangle around node 1
        ! that contains the direction -edge (upwind), and around node 2 the triangle
        ! containing +edge (downwind). Decompose b (a triangle side) and the search
        ! direction x along c (the other side) and the 90deg-CCW normal, then compare
        ! the atan2 angles.
        !
        ! M2.12b (optional partit): a halo node's element list (nod_in_elem2D) reaches
        ! eXDim elements whose elem2D_nodes are NOT stored locally and whose nodes may
        ! be beyond the local node list. So at npes>1 we build, exactly as FESOM2
        ! (oce_muscl_adv.F90:192-235), coord_elem(2,3,*) and e_nodes(3,*) over the FULL
        ! element halo via exchange_elem_full, and pick the triangle sides from those +
        ! the GLOBAL anchor id (myList_nod2D). At 1-rank (partit absent) the proven path
        ! through tri_sides(elem2D_nodes/coord_nod2D) is used unchanged.
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
            allocate(coord_elem(2, 3, nElemF), e_nodes(3, nElemF))
            allocate(tmp(nElemF), tmp_i(nElemF))
            do nn = 1, 3
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

            ! Find upwind (in the sense of x) triangle, i.e. which contains -x:
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
                ! Since b,c are triangle sides, |ab|<pi, so atan2 is what is needed
                if ((ab > 0.0_WP) .and. (ax > 0.0_WP) .and. (ax < ab)) then
                    twork%edge_up_dn_tri(1, n) = elem; cycle
                end if
                if ((ab < 0.0_WP) .and. (ax < 0.0_WP) .and. (ax > ab)) then
                    twork%edge_up_dn_tri(1, n) = elem; cycle
                end if
                if ((ab == ax) .or. (ax == 0.0_WP)) then
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
                if ((ab == ax) .or. (ax == 0.0_WP)) then
                    twork%edge_up_dn_tri(2, n) = elem; cycle
                end if
            end do
        end do

        if (lmr) deallocate(coord_elem, e_nodes)

        ! For edges touching the boundary, up/downwind elements may be absent; the
        ! MUSCL kernels then use a zero reconstruction increment on that side.

    contains
        subroutine get_bc(el, anchor_loc, bb, cc)
            ! Triangle sides from vertex `anchor_loc`. 1-rank: via elem2D_nodes/
            ! coord_nod2D (tri_sides). Multi-rank: via coord_elem + the GLOBAL anchor
            ! id matched against e_nodes (FESOM2 oce_muscl_adv.F90:251-303).
            integer,       intent(in)  :: el, anchor_loc
            real(kind=WP), intent(out) :: bb(2), cc(2)
            integer :: anchor_g
            if (.not. lmr) then
                call tri_sides(mesh, el, anchor_loc, bb, cc)
                return
            end if
            anchor_g = partit%myList_nod2D(anchor_loc)
            if (e_nodes(1, el) == anchor_g) then
                bb = coord_elem(:, 2, el) - coord_elem(:, 1, el)
                cc = coord_elem(:, 3, el) - coord_elem(:, 1, el)
            elseif (e_nodes(2, el) == anchor_g) then
                bb = coord_elem(:, 1, el) - coord_elem(:, 2, el)
                cc = coord_elem(:, 3, el) - coord_elem(:, 2, el)
            else
                bb = coord_elem(:, 1, el) - coord_elem(:, 3, el)
                cc = coord_elem(:, 2, el) - coord_elem(:, 3, el)
            end if
        end subroutine get_bc
    end subroutine find_up_downwind_triangles

    pure subroutine tri_sides(mesh, elem, anchor, b, c)
        ! The two triangle sides emanating from the vertex 'anchor' of 'elem', in
        ! FESOM2's b/c convention (oce_muscl_adv.F90:251-260): if anchor is vertex 1
        ! then b=v2-v1, c=v3-v1; if vertex 2 then b=v1-v2, c=v3-v2; else b=v1-v3,
        ! c=v2-v3. (1-rank: coord_elem == coord_nod2D(:,elem2D_nodes(:,elem)).)
        type(t_mesh),  intent(in)  :: mesh
        integer,       intent(in)  :: elem, anchor
        real(kind=WP), intent(out) :: b(2), c(2)
        integer :: e1, e2, e3
        e1 = mesh%elem2D_nodes(1, elem)
        e2 = mesh%elem2D_nodes(2, elem)
        e3 = mesh%elem2D_nodes(3, elem)
        if (e1 == anchor) then
            b = mesh%coord_nod2D(:, e2) - mesh%coord_nod2D(:, e1)
            c = mesh%coord_nod2D(:, e3) - mesh%coord_nod2D(:, e1)
        elseif (e2 == anchor) then
            b = mesh%coord_nod2D(:, e1) - mesh%coord_nod2D(:, e2)
            c = mesh%coord_nod2D(:, e3) - mesh%coord_nod2D(:, e2)
        else
            b = mesh%coord_nod2D(:, e1) - mesh%coord_nod2D(:, e3)
            c = mesh%coord_nod2D(:, e2) - mesh%coord_nod2D(:, e3)
        end if
    end subroutine tri_sides

    !---------------------------------------------------------------------------
    subroutine fill_up_dn_grad(eudg, twork, tr_xy, mesh, partit)
        ! oce_muscl_adv.F90:356-525. Per edge, build the up/downwind elemental tracer
        ! gradient eudg(1:4,nz,edge): (1,3)=upwind (x,y), (2,4)=downwind.
        ! NOT USED BY THE MODEL: FESOM2 stores this array in twork%edge_up_dn_grad; FESOM3
        ! forms the reconstruction on the fly from tr_xy of the up/downwind triangle (the
        ! reference RK3 rule, no Miura node average). Kept as FESOM2's reference for
        ! test_muscl_onthefly (part R) and for the legacy FESOM2 dump drivers
        ! (fesom_advhordump*). eudg is zeroed here, so levels the fill does not write are
        ! 0, as in FESOM2's once-zeroed array.
        ! On shared levels take the gradient straight from edge_up_dn_tri; on
        ! not-shared levels (and on boundary edges) area-weighted-average tr_xy over
        ! the triangles around each edge node (standard Miura).
        ! M2.12b: optional partit -> loop OWNED edges. The body reads nod_in_elem2D
        ! (completed for halo nodes by the find_neighbors dance), tr_xy and elem_area
        ! at halo elements (both halo-exchanged by the caller / compute_geometry).
        type(t_mesh),        intent(in)    :: mesh
        real(kind=WP),       intent(out)   :: eudg(4, mesh%nl-1, *)   ! (4, nl-1, nEdgeO)
        type(t_tracer_work), intent(in)    :: twork
        real(kind=WP),       intent(in)    :: tr_xy(2, mesh%nl-1, mesh%elem2D)
        type(t_partit),      intent(in), optional :: partit
        integer       :: edge, nz, elem, k, ednodes(2), nzmin, nzmax
        real(kind=WP) :: tvol, tx, ty
        integer       :: nNodO, nNodL, nEdgeO, nElemO

        call owned_bounds(mesh, nNodO, nNodL, nEdgeO, nElemO, partit)
        eudg(:, :, 1:nEdgeO) = 0.0_WP
        do edge = 1, nEdgeO
            ednodes = mesh%edges(:, edge)
            !___ edge has both upwind and downwind triangle on the surface __________
            if ((twork%edge_up_dn_tri(1, edge) /= 0) .and. (twork%edge_up_dn_tri(2, edge) /= 0)) then
                nzmin = maxval(mesh%ulevels_nod2D_max(ednodes))
                nzmax = minval(mesh%nlevels_nod2D_min(ednodes))

                ! not-shared upper levels of edge node 1
                do nz = mesh%ulevels_nod2D(ednodes(1)), nzmin-1
                    tvol = 0.0_WP; tx = 0.0_WP; ty = 0.0_WP
                    do k = 1, mesh%nod_in_elem2D_num(ednodes(1))
                        elem = mesh%nod_in_elem2D(k, ednodes(1))
                        if (mesh%nlevels(elem)-1 < nz .or. nz < mesh%ulevels(elem)) cycle
                        tvol = tvol + mesh%elem_area(elem)
                        tx = tx + tr_xy(1, nz, elem)*mesh%elem_area(elem)
                        ty = ty + tr_xy(2, nz, elem)*mesh%elem_area(elem)
                    end do
                    eudg(1, nz, edge) = tx/tvol
                    eudg(3, nz, edge) = ty/tvol
                end do
                ! not-shared upper levels of edge node 2
                do nz = mesh%ulevels_nod2D(ednodes(2)), nzmin-1
                    tvol = 0.0_WP; tx = 0.0_WP; ty = 0.0_WP
                    do k = 1, mesh%nod_in_elem2D_num(ednodes(2))
                        elem = mesh%nod_in_elem2D(k, ednodes(2))
                        if (mesh%nlevels(elem)-1 < nz .or. nz < mesh%ulevels(elem)) cycle
                        tvol = tvol + mesh%elem_area(elem)
                        tx = tx + tr_xy(1, nz, elem)*mesh%elem_area(elem)
                        ty = ty + tr_xy(2, nz, elem)*mesh%elem_area(elem)
                    end do
                    eudg(2, nz, edge) = tx/tvol
                    eudg(4, nz, edge) = ty/tvol
                end do
                ! shared levels: take gradient straight from up/downwind triangle
                do nz = nzmin, nzmax-1
                    eudg(1:2, nz, edge) = tr_xy(1, nz, twork%edge_up_dn_tri(:, edge))
                    eudg(3:4, nz, edge) = tr_xy(2, nz, twork%edge_up_dn_tri(:, edge))
                end do
                ! not-shared lower levels of edge node 1
                do nz = nzmax, mesh%nlevels_nod2D(ednodes(1))-1
                    tvol = 0.0_WP; tx = 0.0_WP; ty = 0.0_WP
                    do k = 1, mesh%nod_in_elem2D_num(ednodes(1))
                        elem = mesh%nod_in_elem2D(k, ednodes(1))
                        if (mesh%nlevels(elem)-1 < nz .or. nz < mesh%ulevels(elem)) cycle
                        tvol = tvol + mesh%elem_area(elem)
                        tx = tx + tr_xy(1, nz, elem)*mesh%elem_area(elem)
                        ty = ty + tr_xy(2, nz, elem)*mesh%elem_area(elem)
                    end do
                    eudg(1, nz, edge) = tx/tvol
                    eudg(3, nz, edge) = ty/tvol
                end do
                ! not-shared lower levels of edge node 2
                do nz = nzmax, mesh%nlevels_nod2D(ednodes(2))-1
                    tvol = 0.0_WP; tx = 0.0_WP; ty = 0.0_WP
                    do k = 1, mesh%nod_in_elem2D_num(ednodes(2))
                        elem = mesh%nod_in_elem2D(k, ednodes(2))
                        if (mesh%nlevels(elem)-1 < nz .or. nz < mesh%ulevels(elem)) cycle
                        tvol = tvol + mesh%elem_area(elem)
                        tx = tx + tr_xy(1, nz, elem)*mesh%elem_area(elem)
                        ty = ty + tr_xy(2, nz, elem)*mesh%elem_area(elem)
                    end do
                    eudg(2, nz, edge) = tx/tvol
                    eudg(4, nz, edge) = ty/tvol
                end do
            !___ edge has only one triangle on the surface (boundary edge) __________
            else
                ! Only linear reconstruction part (standard Miura at both nodes)
                nzmin = mesh%ulevels_nod2D(ednodes(1))
                nzmax = mesh%nlevels_nod2D(ednodes(1))
                do nz = nzmin, nzmax-1
                    tvol = 0.0_WP; tx = 0.0_WP; ty = 0.0_WP
                    do k = 1, mesh%nod_in_elem2D_num(ednodes(1))
                        elem = mesh%nod_in_elem2D(k, ednodes(1))
                        if (mesh%nlevels(elem)-1 < nz .or. nz < mesh%ulevels(elem)) cycle
                        tvol = tvol + mesh%elem_area(elem)
                        tx = tx + tr_xy(1, nz, elem)*mesh%elem_area(elem)
                        ty = ty + tr_xy(2, nz, elem)*mesh%elem_area(elem)
                    end do
                    eudg(1, nz, edge) = tx/tvol
                    eudg(3, nz, edge) = ty/tvol
                end do
                nzmin = mesh%ulevels_nod2D(ednodes(2))
                nzmax = mesh%nlevels_nod2D(ednodes(2))
                do nz = nzmin, nzmax-1
                    tvol = 0.0_WP; tx = 0.0_WP; ty = 0.0_WP
                    do k = 1, mesh%nod_in_elem2D_num(ednodes(2))
                        elem = mesh%nod_in_elem2D(k, ednodes(2))
                        if (mesh%nlevels(elem)-1 < nz .or. nz < mesh%ulevels(elem)) cycle
                        tvol = tvol + mesh%elem_area(elem)
                        tx = tx + tr_xy(1, nz, elem)*mesh%elem_area(elem)
                        ty = ty + tr_xy(2, nz, elem)*mesh%elem_area(elem)
                    end do
                    eudg(2, nz, edge) = tx/tvol
                    eudg(4, nz, edge) = ty/tvol
                end do
            end if
        end do
    end subroutine fill_up_dn_grad


end module oce_muscl_adv
