module mod_ice_setup
    ! M3a: sea-ice foundation — allocation + mass matrix + cold-start initial state.
    ! Faithful transcription of FESOM2 v2.7.3 ice_init (MOD_ICE.F90:572) +
    ! ice_mass_matrix_fill (ice_fct.F90:1145) + ice_initial_state (ice_setup_step.F90:358),
    ! standard EVP path (whichEVP=0; NO icepack/meltponds/cavity/oasis-yac/oifs).
    !
    ! M2.12 optional-partit pattern: partit absent OR npes==1 -> the proven 1-rank path
    ! VERBATIM (owned_bounds returns the global mesh counts). The ice mass matrix rides
    ! the M2.6 ssh_stiff CSR (rowptr_loc/colind_loc) exactly like FESOM2 rides ssh_stiff
    ! (rowptr/nn_num/nn_pos) — same sparsity, same nza, so fct_massmatrix is byte-identical.
    !
    ! Hybrid meshes (plan docs/plans/completed/2026-10-08-hybrid-meshes.md, Task 16): the mass matrix
    ! uses the user's element contributions (mod_ice_fct ice_mcoef: triangle S/6, S/12;
    ! quad (S/36)(4,2,1,2) rotated per row) on the SSH sparsity incl. quad diagonals, and
    ! its row sums are stored once as the FCT lumped mass (fct_lumpedmass).
    use mpi
    use mod_precision,    only: WP
    use mod_mesh,         only: t_mesh, MAX_NV
    use mod_partit,       only: t_partit
    use mod_tracer,       only: t_tracer
    use mod_ice,          only: t_ice
    use mod_part_bounds,  only: owned_bounds, local_dims, is_multirank
    use mod_halo,         only: exchange_nod
    use mod_ice_fct,      only: ice_mcoef
    implicit none
    private
    public :: ice_setup, ice_allocate, ice_mass_matrix_fill, ice_initial_state
    public :: build_bc_index_nod2D   ! also used by test_ice_dyn_mixed

contains

    !___________________________________________________________________________
    ! ice_setup (FESOM2 ice_setup_step.F90:51): allocate -> ice_dt/Tevp_inv/Clim_evp
    ! -> mass matrix -> cold-start IC. dt is the ocean timestep (ice_dt = ice_ave_steps*dt).
    subroutine ice_setup(ice, tracers, mesh, dt, partit)
        type(t_ice),    intent(inout)         :: ice
        type(t_tracer), intent(in)            :: tracers
        type(t_mesh),   intent(in)            :: mesh
        real(kind=WP),  intent(in)            :: dt
        type(t_partit), intent(in), optional  :: partit

        call ice_allocate(ice, mesh, partit)

        ! DO not change (FESOM2 ice_setup:74-78). evp_rheol_steps is integer -> the
        ! (evp_rheol_steps/ice_dt) divide is int/real promotion, kept verbatim.
        ice%ice_dt   = real(ice%ice_ave_steps, WP) * dt
        ice%Tevp_inv = 3.0_WP / ice%ice_dt
        ice%Clim_evp = ice%Clim_evp * (ice%evp_rheol_steps/ice%ice_dt)**2 / ice%Tevp_inv

        call build_bc_index_nod2D(ice, mesh, partit)
        call ice_mass_matrix_fill(ice, mesh, partit)
        call ice_initial_state(ice, tracers, mesh, partit)
    end subroutine ice_setup

    !___________________________________________________________________________
    ! build_bc_index_nod2D (FESOM2 ice_init MOD_ICE.F90:888-895). The node boundary
    ! mask used by the mEVP velocity solve (whichEVP/=0): 1 in the interior, 0 at the
    ! two nodes of every boundary edge (global id > edge2D_in). Built unconditionally
    ! (FESOM2 "also for whichEVP==0"). owned edges only; an owned node's incident edges
    ! are all owned (invariant, M2.12c-2) so owned-node bc_index is complete (no exchange,
    ! exactly like FESOM2). FESOM2 keeps this in mesh%bc_index_nod2D; here in ice (the
    ! FESOM3 t_mesh type-bound I/O trips an ifort cascade when a component is added).
    subroutine build_bc_index_nod2D(ice, mesh, partit)
        type(t_ice),    intent(inout)         :: ice
        type(t_mesh),   intent(in)            :: mesh
        type(t_partit), intent(in), optional  :: partit
        integer :: n, nNodO, nNodL, nEdgeO, nElemO
        logical :: lmr

        call owned_bounds(mesh, nNodO, nNodL, nEdgeO, nElemO, partit)
        lmr = is_multirank(partit)

        ice%bc_index_nod2D = 1.0_WP
        do n = 1, nEdgeO
            if (lmr) then
                if (partit%myList_edge2D(n) <= mesh%edge2D_in) cycle
            else
                if (n <= mesh%edge2D_in) cycle
            end if
            ice%bc_index_nod2D(mesh%edges(1:2,n)) = 0.0_WP
        end do
    end subroutine build_bc_index_nod2D

    !___________________________________________________________________________
    ! ice_allocate (FESOM2 ice_init MOD_ICE.F90:711-905, standard-EVP subset).
    ! node arrays sized nNodL (owned+halo); element arrays nElemL (owned+eDim);
    ! fct_massmatrix sized ssh_stiff%nza (owned CSR rows) — REQUIRES ssh_stiff built
    ! (init_stiff_mat_ale) before this call.
    subroutine ice_allocate(ice, mesh, partit)
        type(t_ice),    intent(inout)         :: ice
        type(t_mesh),   intent(in)            :: mesh
        type(t_partit), intent(in), optional  :: partit
        integer :: node_size, elem_size, n
        integer :: nNodO, nNodL, nEdgeO, nEdgeL, nElemO, nElemL, nElemF

        call local_dims(mesh, partit, nNodO, nNodL, nEdgeO, nEdgeL, nElemO, nElemL, nElemF)
        node_size = nNodL
        elem_size = nElemL

        ! velocity + stress (nodes). uice_aux/vice_aux are the mEVP (whichEVP/=0)
        ! auxiliary solver velocities (FESOM2 ice_init allocates them when whichEVP/=0);
        ! allocated unconditionally here (cheap, lets the EVP dispatcher pick at runtime).
        allocate(ice%uice(node_size), ice%uice_rhs(node_size), ice%uice_old(node_size))
        allocate(ice%vice(node_size), ice%vice_rhs(node_size), ice%vice_old(node_size))
        allocate(ice%uice_aux(node_size), ice%vice_aux(node_size))
        allocate(ice%stress_atmice_x(node_size), ice%stress_iceoce_x(node_size))
        allocate(ice%stress_atmice_y(node_size), ice%stress_iceoce_y(node_size))
        allocate(ice%h_ice(node_size), ice%h_snow(node_size))
        allocate(ice%bc_index_nod2D(node_size))
        ice%uice = 0.0_WP; ice%uice_rhs = 0.0_WP; ice%uice_old = 0.0_WP
        ice%vice = 0.0_WP; ice%vice_rhs = 0.0_WP; ice%vice_old = 0.0_WP
        ice%uice_aux = 0.0_WP; ice%vice_aux = 0.0_WP
        ice%bc_index_nod2D = 1.0_WP
        ice%stress_atmice_x = 0.0_WP; ice%stress_iceoce_x = 0.0_WP
        ice%stress_atmice_y = 0.0_WP; ice%stress_iceoce_y = 0.0_WP
        ice%h_ice = 0.0_WP; ice%h_snow = 0.0_WP

        ! surface ocean arrays (nodes)
        allocate(ice%srfoce_u(node_size), ice%srfoce_v(node_size))
        allocate(ice%srfoce_temp(node_size), ice%srfoce_salt(node_size), ice%srfoce_ssh(node_size))
        ice%srfoce_u = 0.0_WP; ice%srfoce_v = 0.0_WP
        ice%srfoce_temp = 0.0_WP; ice%srfoce_salt = 0.0_WP; ice%srfoce_ssh = 0.0_WP

        ! freshwater & heat flux (nodes)
        allocate(ice%flx_fw(node_size), ice%flx_h(node_size))
        ice%flx_fw = 0.0_WP; ice%flx_h = 0.0_WP

        ! ice tracers data(1:3) = a_ice/m_ice/m_snow
        allocate(ice%data(ice%num_itracers))
        do n = 1, ice%num_itracers
            allocate(ice%data(n)%values(node_size), ice%data(n)%values_old(node_size))
            allocate(ice%data(n)%values_rhs(node_size), ice%data(n)%values_div_rhs(node_size))
            allocate(ice%data(n)%dvalues(node_size), ice%data(n)%valuesl(node_size))
            ice%data(n)%ID             = n
            ice%data(n)%values         = 0.0_WP
            ice%data(n)%values_old     = 0.0_WP
            ice%data(n)%values_rhs     = 0.0_WP
            ice%data(n)%values_div_rhs = 0.0_WP
            ice%data(n)%dvalues        = 0.0_WP
            ice%data(n)%valuesl        = 0.0_WP
        end do

        ! work arrays: FCT (nodes) + mass matrix (CSR) + stress/strain (elements)
        allocate(ice%work%fct_tmax(node_size), ice%work%fct_tmin(node_size))
        allocate(ice%work%fct_plus(node_size), ice%work%fct_minus(node_size))
        allocate(ice%work%fct_fluxes(elem_size, MAX_NV))
        ice%work%fct_tmax = 0.0_WP; ice%work%fct_tmin = 0.0_WP
        ice%work%fct_plus = 0.0_WP; ice%work%fct_minus = 0.0_WP
        ice%work%fct_fluxes = 0.0_WP

        allocate(ice%work%fct_massmatrix(mesh%ssh_stiff%nza), ice%work%fct_lumpedmass(node_size))
        ice%work%fct_massmatrix = 0.0_WP; ice%work%fct_lumpedmass = 0.0_WP

        allocate(ice%work%sigma11(elem_size), ice%work%sigma12(elem_size), ice%work%sigma22(elem_size))
        allocate(ice%work%eps11(elem_size), ice%work%eps12(elem_size), ice%work%eps22(elem_size))
        ice%work%sigma11 = 0.0_WP; ice%work%sigma12 = 0.0_WP; ice%work%sigma22 = 0.0_WP
        ice%work%eps11 = 0.0_WP; ice%work%eps12 = 0.0_WP; ice%work%eps22 = 0.0_WP

        allocate(ice%work%ice_strength(elem_size), ice%work%inv_areamass(node_size), ice%work%inv_mass(node_size))
        ice%work%ice_strength = 0.0_WP; ice%work%inv_areamass = 0.0_WP; ice%work%inv_mass = 0.0_WP

        ! thermo work arrays (nodes)
        allocate(ice%thermo%ustar(node_size), ice%thermo%t_skin(node_size))
        allocate(ice%thermo%thdgr(node_size), ice%thermo%thdgrsn(node_size))
        allocate(ice%thermo%thdgra(node_size), ice%thermo%thdgr_old(node_size))
        allocate(ice%thermo%dyngr(node_size), ice%thermo%dyngrsn(node_size), ice%thermo%dyngra(node_size))
        ice%thermo%ustar = 0.0_WP; ice%thermo%t_skin = 0.0_WP
        ice%thermo%thdgr = 0.0_WP; ice%thermo%thdgrsn = 0.0_WP
        ice%thermo%thdgra = 0.0_WP; ice%thermo%thdgr_old = 0.0_WP
        ice%thermo%dyngr = 0.0_WP; ice%thermo%dyngrsn = 0.0_WP; ice%thermo%dyngra = 0.0_WP

        ! NOTE: bc_index_nod2D (FESOM2 ice_init:889) is deferred to M3b (EVP) — the only
        ! consumer is the EVP boundary handling, not the IC/mass-matrix gated here.
    end subroutine ice_allocate

    !___________________________________________________________________________
    ! ice_mass_matrix_fill (FESOM2 ice_fct.F90:1145). Builds the consistent ice FCT mass
    ! matrix on the ssh_stiff CSR sparsity and its row sums, the lumped mass. FESOM2 locates
    ! the CSR slot for the (row, elnodes(q)) pair via nn_num/nn_pos; here we scan the
    ! equivalent ssh_stiff%colind_loc(rowptr_loc(row):rowptr_loc(row+1)-1) — the SAME
    ! ordered neighbour list (init_stiff_mat_ale fills colind_loc = n_pos(1:n_num,row)).
    ! A column that is not in the row stops the run (the sparsity does not match the mesh).
    ! Element contributions elem_area*ice_mcoef(:,:,nv) over the nv unique vertices
    ! (FESOM2: S/12 off-diagonal, S/6 diagonal = the triangle case). The lumped mass
    ! fct_lumpedmass = row sums (S/3 per triangle vertex, S/4 per quad vertex; FESOM2 used
    ! the node area instead and only warned where they differ) is exchanged to the halo
    ! (the antidiffusive fluxes of owned elements divide by it at halo vertices). A node
    ! with no wet element (all its elements cavity) has a zero row; it keeps the node area
    ! there (FESOM2 behaviour) so that the FCT divides stay finite.
    subroutine ice_mass_matrix_fill(ice, mesh, partit)
        type(t_ice),    intent(inout)         :: ice
        type(t_mesh),   intent(in)            :: mesh
        type(t_partit), intent(in), optional  :: partit
        integer :: n, row, elem, elnodes(MAX_NV), nv, q, ipos, nini, nend
        integer :: nNodO, nNodL, nEdgeO, nElemO

        call owned_bounds(mesh, nNodO, nNodL, nEdgeO, nElemO, partit)

        associate(ssh_stiff => mesh%ssh_stiff)
        do elem = 1, nElemO
            if (mesh%ulevels(elem) > 1) cycle   ! cavity element -> no mass
            elnodes = mesh%elem2D_nodes(:, elem)
            nv      = mesh%elem2D_nnodes(elem)
            do n = 1, nv
                row = elnodes(n)
                if (row > nNodO) cycle
                nini = ssh_stiff%rowptr_loc(row)
                nend = ssh_stiff%rowptr_loc(row+1) - 1
                do q = 1, nv
                    do ipos = nini, nend
                        if (ssh_stiff%colind_loc(ipos) == elnodes(q)) exit
                    end do
                    if (ipos > nend) call column_not_found(row, elnodes(q), elem)
                    ice%work%fct_massmatrix(ipos) = ice%work%fct_massmatrix(ipos) + &
                                                    mesh%elem_area(elem)*ice_mcoef(n,q,nv)
                end do
            end do
        end do

        ! lumped mass = row sums
        do row = 1, nNodO
            nini = ssh_stiff%rowptr_loc(row)
            nend = ssh_stiff%rowptr_loc(row+1) - 1
            ice%work%fct_lumpedmass(row) = sum(ice%work%fct_massmatrix(nini:nend))
            if (ice%work%fct_lumpedmass(row) <= 0.0_WP) ice%work%fct_lumpedmass(row) = mesh%area(row)
        end do
        end associate
        if (is_multirank(partit)) call exchange_nod(ice%work%fct_lumpedmass, partit)

    contains

        subroutine column_not_found(row, col, elem)
            integer, intent(in) :: row, col, elem
            integer :: grow, gcol, gelem, comm, ierr
            grow = row; gcol = col; gelem = elem; comm = MPI_COMM_WORLD
            if (present(partit)) then
                comm = partit%MPI_COMM_FESOM
                if (is_multirank(partit)) then
                    grow  = partit%myList_nod2D(row)
                    gcol  = partit%myList_nod2D(col)
                    gelem = partit%myList_elem2D(elem)
                end if
            end if
            write(*,'(a,i0,a,i0,a,i0,a)') 'ice_mass_matrix_fill: column of node ', gcol, &
                ' (element ', gelem, ') not found in the ssh_stiff CSR row of node ', grow, &
                ': the SSH sparsity does not match the mesh'
            call MPI_Abort(comm, 1, ierr)
        end subroutine column_not_found
    end subroutine ice_mass_matrix_fill

    !___________________________________________________________________________
    ! ice_initial_state (FESOM2 ice_setup_step.F90:358), cold start (ini_ice_from_file
    ! = .false.). a_ice/m_ice/m_snow/uice/vice = 0 everywhere, then where surface ocean
    ! temperature < 0: NH m_ice=1/m_snow=0.1, SH m_ice=2/m_snow=0.5, a_ice=0.9. Cavity
    ! nodes (ulevels_nod2D>1) skipped. SST<0 and lat>0 are SIGN tests -> robust to ULP.
    subroutine ice_initial_state(ice, tracers, mesh, partit)
        type(t_ice),    intent(inout)         :: ice
        type(t_tracer), intent(in)            :: tracers
        type(t_mesh),   intent(in)            :: mesh
        type(t_partit), intent(in), optional  :: partit
        integer :: i, nNodO, nNodL, nEdgeO, nElemO

        call owned_bounds(mesh, nNodO, nNodL, nEdgeO, nElemO, partit)

        ice%data(1)%values = 0.0_WP   ! a_ice
        ice%data(2)%values = 0.0_WP   ! m_ice
        ice%data(3)%values = 0.0_WP   ! m_snow
        ice%uice = 0.0_WP
        ice%vice = 0.0_WP

        do i = 1, nNodL
            ! if cavity, no sea ice, no initial state
            if (mesh%ulevels_nod2D(i) > 1) cycle
            if (tracers%data(1)%values(1, i) < 0.0_WP) then
                if (mesh%geo_coord_nod2D(2, i) > 0._WP) then
                    ice%data(2)%values(i) = 1.0_WP
                    ice%data(3)%values(i) = 0.1_WP
                else
                    ice%data(2)%values(i) = 2.0_WP
                    ice%data(3)%values(i) = 0.5_WP
                end if
                ice%data(1)%values(i) = 0.9_WP
                ice%uice(i) = 0.0_WP
                ice%vice(i) = 0.0_WP
            end if
        end do
    end subroutine ice_initial_state

end module mod_ice_setup
