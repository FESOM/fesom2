module muscl_oracle
    ! THE STORED-ARRAY KERNEL: FESOM3's MUSCL horizontal-flux kernel (FESOM2's MFCT, no
    ! boundary clamp) as it was before the on-the-fly refactor (verbatim, commit 48cc272: it
    ! reads a per-edge up/downwind gradient array edge_up_dn_grad). Fed with the array built
    ! by the reference RK3 rule it gives that rule's fluxes: test_muscl_onthefly part R, the
    ! specification of the production kernel.
    use mod_precision,   only: WP
    use mod_mesh,        only: t_mesh
    use mod_partit,      only: t_partit
    use mod_part_bounds, only: owned_bounds
    implicit none
    private
    public :: oracle_muscl
contains

    subroutine oracle_muscl(vel, ttf, mesh, num_ord, flux, edge_up_dn_grad, o_init_zero, partit)
        ! MUSCL (FESOM2 adv_tra_hor_mfct, oce_adv_tra_hor.F90:546-834), no boundary clamp.
        ! M2.12b: optional partit -> loop over OWNED edges (myDim_edge2D).
        type(t_mesh),  intent(in)    :: mesh
        real(kind=WP), intent(in)    :: num_ord
        real(kind=WP), intent(in)    :: ttf(mesh%nl-1, mesh%nod2D)
        real(kind=WP), intent(in)    :: vel(2, mesh%nl-1, mesh%elem2D)
        real(kind=WP), intent(inout) :: flux(mesh%nl-1, mesh%edge2D)
        real(kind=WP), intent(in)    :: edge_up_dn_grad(4, mesh%nl-1, mesh%edge2D)
        logical, optional, intent(in) :: o_init_zero
        type(t_partit), intent(in), optional :: partit
        logical       :: l_init_zero
        real(kind=WP) :: deltaX1, deltaY1, deltaX2, deltaY2, vflux
        integer       :: el(2), enodes(2), nz, edge, nu12, nl12, nl1, nl2, nu1, nu2
        integer       :: nNodO, nNodL, nEdgeO, nElemO

        call owned_bounds(mesh, nNodO, nNodL, nEdgeO, nElemO, partit)
        l_init_zero = .true.
        if (present(o_init_zero)) l_init_zero = o_init_zero
        if (l_init_zero) then
            do edge = 1, nEdgeO
                flux(:, edge) = 0.0_WP
            end do
        end if

        do edge = 1, nEdgeO
            enodes = mesh%edges(:, edge)
            el     = mesh%edge_tri(:, edge)
            nl1    = mesh%nlevels(el(1)) - 1
            nu1    = mesh%ulevels(el(1))
            deltaX1 = mesh%edge_cross_dxdy(1, edge)
            deltaY1 = mesh%edge_cross_dxdy(2, edge)
            nl2 = 0; nu2 = 0
            if (el(2) /= el(1)) then
                deltaX2 = mesh%edge_cross_dxdy(3, edge)
                deltaY2 = mesh%edge_cross_dxdy(4, edge)
                nl2 = mesh%nlevels(el(2)) - 1
                nu2 = mesh%ulevels(el(2))
            end if
            nl12 = min(nl1, nl2)
            nu12 = max(nu1, nu2)
            ! (A)
            do nz = nu1, nu12-1
                vflux = (-vel(2,nz,el(1))*deltaX1 + vel(1,nz,el(1))*deltaY1) * mesh%helem(nz,el(1))
                call flux_ho(nz, vflux)
            end do
            ! (B)
            if (nu2 > 0) then
                do nz = nu2, nu12-1
                    vflux = (vel(2,nz,el(2))*deltaX2 - vel(1,nz,el(2))*deltaY2) * mesh%helem(nz,el(2))
                    call flux_ho(nz, vflux)
                end do
            end if
            ! (C)
            do nz = nu12, nl12
                vflux = (-vel(2,nz,el(1))*deltaX1 + vel(1,nz,el(1))*deltaY1) * mesh%helem(nz,el(1)) &
                      + ( vel(2,nz,el(2))*deltaX2 - vel(1,nz,el(2))*deltaY2) * mesh%helem(nz,el(2))
                call flux_ho(nz, vflux)
            end do
            ! (D)
            do nz = nl12+1, nl1
                vflux = (-vel(2,nz,el(1))*deltaX1 + vel(1,nz,el(1))*deltaY1) * mesh%helem(nz,el(1))
                call flux_ho(nz, vflux)
            end do
            ! (E)
            do nz = nl12+1, nl2
                vflux = (vel(2,nz,el(2))*deltaX2 - vel(1,nz,el(2))*deltaY2) * mesh%helem(nz,el(2))
                call flux_ho(nz, vflux)
            end do
        end do

    contains
        subroutine flux_ho(nz, vflux)
            integer,       intent(in) :: nz
            real(kind=WP), intent(in) :: vflux
            real(kind=WP) :: Tmean1, Tmean2, cHO
            Tmean2 = ttf(nz, enodes(2)) - &
                     (2.0_WP*(ttf(nz, enodes(2))-ttf(nz, enodes(1))) + &
                      mesh%edge_dxdy(1,edge)*edge_up_dn_grad(2,nz,edge) + &
                      mesh%edge_dxdy(2,edge)*edge_up_dn_grad(4,nz,edge))/6.0_WP
            Tmean1 = ttf(nz, enodes(1)) + &
                     (2.0_WP*(ttf(nz, enodes(2))-ttf(nz, enodes(1))) + &
                      mesh%edge_dxdy(1,edge)*edge_up_dn_grad(1,nz,edge) + &
                      mesh%edge_dxdy(2,edge)*edge_up_dn_grad(3,nz,edge))/6.0_WP
            cHO = (vflux+abs(vflux))*Tmean1 + (vflux-abs(vflux))*Tmean2
            flux(nz,edge) = -0.5_WP*(1.0_WP-num_ord)*cHO - vflux*num_ord*0.5_WP*(Tmean1+Tmean2) - flux(nz,edge)
        end subroutine flux_ho
    end subroutine oracle_muscl

end module muscl_oracle

program test_muscl_onthefly
    ! MUSCL horizontal advection with the up/downwind reconstruction formed ON THE FLY
    ! (docs/plans/completed/2026-10-06-muscl-onthefly.md).
    !
    ! THE RULE (the reference RK3 code, qq/oce_stepRK3.F90 t_hor_adv_muscl_RK3): for each
    ! edge, side k = 1 takes the elemental gradient tr_xy of the UPWIND triangle
    ! edge_up_dn_tri(1, edge), side k = 2 that of the DOWNWIND triangle, on the levels where
    ! that triangle is wet (ulevels..nlevels-1); a missing or dry triangle gives a ZERO
    ! increment. No stored per-edge array, no Miura node average (FESOM2's fill_up_dn_grad,
    ! removed from FESOM3; the comparison with it is recorded in docs/HANDOFF.md).
    !
    ! PART R - the specification: adv_tra_hor_muscl == the stored-array
    !   kernel (module muscl_oracle) fed with the array built by the rule above, every owned
    !   edge and level, both velocity signs, BITWISE; on pi as read (R1) and with synthesised
    !   cavity columns (R2: ulevels_nod2D = 3 on every 10th global node with >= 14 levels,
    !   element ulevels and ulevels_nod2D(_max) rebuilt consistently and exchanged).
    !   Hybrid meshes (Task 11): FESOM3_MESH_DIR selects the mesh (default pi, rotated);
    !   FESOM3_TEST_ROTATE=0 reads it unrotated (the NEMO quad / mixed cut-outs), where the
    !   up/downwind elements are triangles and quads (find_up_downwind_elements).
    use mpi
    use mod_precision,     only: WP, MP
    use mod_mesh,          only: t_mesh
    use mod_partit,        only: t_partit
    use mod_tracer,        only: t_tracer_work
    use mod_partitioning,  only: par_init, par_ex, set_partition
    use mod_mesh_read,     only: read_mesh
    use mod_mesh_areas,    only: compute_geometry
    use mod_part_bounds,   only: owned_bounds, is_multirank
    use mod_halo,          only: exchange_elem_full, exchange_nod
    use oce_tracer_grad,   only: tracer_gradient_elements
    use oce_muscl_adv,     only: muscl_adv_init
    use oce_adv_tra_hor,   only: adv_tra_hor_muscl
    use muscl_oracle,      only: oracle_muscl
    implicit none

    character(len=512) :: mesh_dir, cbuf
    type(t_partit)      :: partit
    type(t_mesh)        :: mesh
    type(t_tracer_work) :: twork
    integer :: nfail, nsw
    integer :: nNodO, nNodL, nEdgeO, nElemO, nElemF, nl
    real(kind=WP), allocatable :: ttf(:,:), tr_xy(:,:,:)
    real(kind=WP), allocatable :: vel(:,:,:), eref(:,:,:), f_ref(:,:), f_new(:,:)
    real(kind=WP), parameter   :: num_ord = 0.25_WP      ! any mix of the 3rd/4th-order parts
    ! cavity synthesis (the constants of test_wsplit / test_wimpl_tra: same column set)
    integer, parameter :: cav_every    = 10      ! every 10th global node ...
    integer, parameter :: cav_ulev     = 3       ! ... gets its surface at level 3 ...
    integer, parameter :: cav_nlev_min = 14      ! ... if it has >= 14 levels (>= 11 layers below)

    nfail = 0
    call get_environment_variable('FESOM3_MESH_DIR', mesh_dir)
    if (len_trim(mesh_dir) == 0) &
        mesh_dir = '/home/a/a270088/port2/fesom2/tests/data/MESHES/pi'
    call par_init(partit)
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
    call compute_geometry(mesh, partit, cartesian=.false.)
    call owned_bounds(mesh, nNodO, nNodL, nEdgeO, nElemO, partit)
    if (partit%mype == 0) write(*,'(a,i0,a,i0)') '  owned elements on rank 0: triangles=', &
        count(mesh%elem2D_nnodes(1:nElemO) == 3), ' quads=', count(mesh%elem2D_nnodes(1:nElemO) == 4)
    nElemF = mesh%elem2D
    if (partit%npes > 1) nElemF = partit%myDim_elem2D + partit%eDim_elem2D + partit%eXDim_elem2D
    nl = mesh%nl
    allocate(ttf(nl-1, nNodL), tr_xy(2, nl-1, nElemF))
    allocate(vel(2, nl-1, nElemF), eref(4, nl-1, nEdgeO))
    allocate(f_ref(nl-1, nEdgeO), f_new(nl-1, nEdgeO))
    allocate(mesh%helem(nl-1, nElemF))
    call muscl_adv_init(twork, mesh, partit)      ! edge_up_dn_tri
    call build_flow()

    call build_gradient()
    call part_r('R1 mesh as read      ')
    call synth_cavity()
    call build_gradient()
    call part_r('R2 mesh with cavities')

    if (partit%mype == 0) then
        if (nfail == 0) then
            write(*,'(a)') 'test_muscl_onthefly: OK'
        else
            write(*,'(a,i0,a)') 'test_muscl_onthefly: ', nfail, ' FAILURE(S)'
        end if
    end if
    call par_ex(partit%MPI_COMM_FESOM, partit%mype)   ! NO third arg: abort is presence-based
    if (nfail /= 0) error stop 1

contains
    subroutine check_true(name, cond)
        character(len=*), intent(in) :: name
        logical,          intent(in) :: cond
        if (.not. cond) then
            nfail = nfail + 1
            if (partit%mype == 0) write(*,'(a)') '  FAIL: '//name
        end if
    end subroutine check_true

    function gsum(i) result(s)
        integer, intent(in) :: i
        integer :: s, ierr
        s = i
        if (partit%npes > 1) call MPI_Allreduce(i, s, 1, MPI_INTEGER, MPI_SUM, &
                                                partit%MPI_COMM_FESOM, ierr)
    end function gsum

    !=========================================================================
    ! tracer field and its elemental gradient (full halo), as init_tracers_AB builds it
    !=========================================================================
    subroutine build_gradient()
        ! smooth large-scale structure + a level- and node-dependent wiggle, all functions
        ! of the node COORDINATES (identical on owner and halo copies)
        integer :: n, nz
        real(kind=WP) :: x, y
        do n = 1, nNodL
            x = mesh%coord_nod2D(1, n); y = mesh%coord_nod2D(2, n)
            do nz = 1, nl-1
                ttf(nz, n) = 10.0_WP*cos(y)*sin(x + 0.1_WP*real(nz, WP)) &
                           + 0.37_WP*sin(17.0_WP*x*y + 1.3_WP*real(nz, WP)) - 0.05_WP*real(nz, WP)
            end do
        end do
        tr_xy = 0.0_WP
        call tracer_gradient_elements(ttf, tr_xy, mesh, partit)
        if (is_multirank(partit)) call exchange_elem_full(tr_xy, partit)
    end subroutine build_gradient

    subroutine build_flow()
        ! element velocity with structure in both components and every level; helem any
        ! positive thickness (the kernels only multiply the edge velocity by it)
        integer :: e, nz
        real(kind=WP) :: a
        do e = 1, nElemF
            a = 0.731_WP*real(mod(e, 97), WP)
            do nz = 1, nl-1
                vel(1, nz, e) = 0.3_WP*sin(a + 0.2_WP*real(nz, WP)) + 0.05_WP
                vel(2, nz, e) = 0.2_WP*cos(1.7_WP*a - 0.1_WP*real(nz, WP))
                mesh%helem(nz, e) = 10.0_MP + 2.0_MP*real(nz, MP)
            end do
        end do
    end subroutine build_flow

    !=========================================================================
    ! the reference rule as a per-edge array (side k: x in component k, y in k+2)
    !=========================================================================
    subroutine build_ref_array()
        integer :: edge, k, tri, nz
        eref = 0.0_WP
        do edge = 1, nEdgeO
            do k = 1, 2
                tri = twork%edge_up_dn_tri(k, edge)
                if (tri <= 0) cycle
                do nz = mesh%ulevels(tri), mesh%nlevels(tri)-1
                    eref(k,   nz, edge) = tr_xy(1, nz, tri)
                    eref(k+2, nz, edge) = tr_xy(2, nz, tri)
                end do
            end do
        end do
    end subroutine build_ref_array

    !=========================================================================
    ! PART R: production kernels == the reference rule, bitwise
    !=========================================================================
    subroutine part_r(label)
        character(len=*), intent(in) :: label
        integer :: isign, n_bad_m, n_flux
        call build_ref_array()
        n_bad_m = 0; n_flux = 0
        do isign = 1, 2
            if (isign == 2) vel = -vel
            call oracle_muscl(vel, ttf, mesh, num_ord, f_ref, eref, o_init_zero=.true., partit=partit)
            call adv_tra_hor_muscl(vel, ttf, mesh, num_ord, f_new, tr_xy, twork%edge_up_dn_tri, &
                                   o_init_zero=.true., partit=partit)
            n_bad_m = n_bad_m + count(transfer(f_ref, 0_8, size(f_ref)) /= transfer(f_new, 0_8, size(f_new)))
            n_flux  = n_flux  + count(f_ref /= 0.0_WP)
        end do
        vel = -vel                                   ! restore
        n_bad_m = gsum(n_bad_m); n_flux = gsum(n_flux)
        if (partit%mype == 0) write(*,'(2a,i0,a,i0)') '  ', label//': nonzero fluxes ', n_flux, &
            '  bitwise mismatches vs reference rule: ', n_bad_m
        call check_true(trim(label)//': adv_tra_hor_muscl == reference rule bitwise', n_bad_m == 0)
        call check_true(trim(label)//': fluxes nonzero (not vacuous)', n_flux > 0)
    end subroutine part_r

    !=========================================================================
    ! Cavity columns, element-consistent
    !=========================================================================
    subroutine synth_cavity()
        ! node ulevels by the GLOBAL index (owner and halo copies agree); element ulevels
        ! = maxval over its vertices (owned elements, then the full element halo);
        ! ulevels_nod2D = minval and ulevels_nod2D_max = maxval over the elements around
        ! the node (owned, then exchange_nod) -- FESOM2's cavity convention.
        integer :: n, g, e, k, ncav
        ncav = 0
        do n = 1, nNodL
            g = n
            if (partit%npes > 1) g = partit%myList_nod2D(n)
            if (mod(g, cav_every) == 0 .and. mesh%nlevels_nod2D(n) >= cav_nlev_min) then
                mesh%ulevels_nod2D(n) = cav_ulev
                if (n <= nNodO) ncav = ncav + 1
            end if
        end do
        do e = 1, nElemO
            mesh%ulevels(e) = maxval(mesh%ulevels_nod2D(mesh%elem2D_nodes(:, e)))   ! 4 slots
        end do
        if (is_multirank(partit)) call exchange_elem_full(mesh%ulevels, partit)
        ! FESOM2's cavity convention: a node's surface is the SHALLOWEST surface of its
        ! elements (ulevels_nod2D = minval(ulevels around)), so a node whose elements all
        ! touch a cavity vertex moves down too; without this its top levels would have no
        ! wet element (FESOM2's fill would then divide 0/0)
        do n = 1, nNodO
            k = mesh%nod_in_elem2D_num(n)
            mesh%ulevels_nod2D(n)     = minval(mesh%ulevels(mesh%nod_in_elem2D(1:k, n)))
            mesh%ulevels_nod2D_max(n) = maxval(mesh%ulevels(mesh%nod_in_elem2D(1:k, n)))
        end do
        if (is_multirank(partit)) then
            call exchange_nod(mesh%ulevels_nod2D, partit)
            call exchange_nod(mesh%ulevels_nod2D_max, partit)
        end if
        ncav = gsum(ncav)
        if (partit%mype == 0) write(*,'(a,i0)') '  synthesised owned cavity columns (ulevels_nod2D = 3): ', ncav
        call check_true('cavity columns synthesised', ncav > 0)
        do e = 1, nElemO
            call check_e(mesh%nlevels(e)-1 >= mesh%ulevels(e))
        end do
    end subroutine synth_cavity

    subroutine check_e(ok)
        logical, intent(in) :: ok
        if (.not. ok) call check_true('synthesised element keeps a wet level', .false.)
    end subroutine check_e
end program test_muscl_onthefly
