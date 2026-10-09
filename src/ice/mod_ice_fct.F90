module mod_ice_fct
    ! M3c: sea-ice FCT advection. Faithful transcription of FESOM2 v2.7.3 ice_fct.F90 —
    ! ice_TG_rhs (Taylor-Galerkin rhs + ice_diff stabilising diffusion) + ice_fct_solve
    ! (ice_solve_high_order -> ice_solve_low_order -> ice_fem_fct Zalesak limiter x3
    ! tracers). Advects the 3 ice tracers (1=a_ice, 2=m_ice, 3=m_snow) with the EVP ice
    ! velocity uice/vice. Standard path (whichEVP=0; NO icepack/meltponds/cavity/oifs).
    !
    ! CSR BASIS (the M3a key finding): FESOM2 locates the (row, neighbour) mass-matrix slot
    ! via mesh%nn_num/nn_pos; FESOM3 uses the ssh_stiff CSR (rowptr_loc/colind_loc) — the
    ! SAME ordered neighbour list ice_mass_matrix_fill scanned to build fct_massmatrix. M3a
    ! PROVED fct_massmatrix byte-matches FESOM2's mass_matrix position-by-position, which
    ! means colind_loc(rowptr_loc(row):rowptr_loc(row+1)-1) IS FESOM2's nn_pos(1:nn_num,row)
    ! as an ORDERED list (incl. self). So pairing fct_massmatrix(clo:clo2) with
    ! colind_loc(clo:clo2) reproduces FESOM2's mass_matrix(clo:clo2)*field(nn_pos) — same
    ! operands, same summation order => byte-identical. (FESOM3's own mesh%nn_pos is built by
    ! muscl_adv_init and may order neighbours differently, so it is NOT used here.)
    !
    ! M2.12 optional-partit pattern (= every M2/M3 kernel): partit absent OR npes==1 -> the
    ! proven 1-rank path VERBATIM (owned_bounds returns global counts, is_multirank=.false.
    ! so every exchange is skipped). present+npes>1 -> owned/halo loop bounds + the FESOM2
    ! exchanges (exchange_nod on valuesl / dvalues / icepplus-minus / final values). The
    ! per-element / per-node arithmetic is UNCHANGED so codegen — and the byte-match —
    ! is preserved; only loop bounds change + guarded exchanges are added. (Multi-rank ice
    ! advection is gated at M3f; M3c is the CORE2 1-rank gate.)
    !
    ! BIT-IDENTITY NOTES (the L9 transitive-gate pattern):
    !  - gradient_sca/elem_area/area are geometry-gate proven (max|delta|=0, M2.11a);
    !    fct_massmatrix is M3a-proven; uice/vice are M3b-proven; a_ice/m_ice/m_snow come from
    !    the cold-start IC (M3a). Every operand is already byte-pinned.
    !  - the per-node rhs accumulation (ice_TG_rhs) and the per-node flux scatter (ice_fem_fct)
    !    over OWNED elements are FP-order-sensitive; the 1-rank global element order matches
    !    FESOM2 exactly (M3a pinned elem2D_nodes order; M3b's stress2rhs proved the same
    !    accumulation pattern byte-matches; the oracle build is ENABLE_OPENMP=OFF so the !$OMP
    !    ORDERED/atomic scatters run serially in element order = FESOM3's serial path).
    !  - the runtime divides (rhs/lumped, vol/lumped, tmax/max(flux,1e-12)) are scalar with byte-
    !    identical operands (-no-prec-div reciprocals agree, L7/L29). Pointers mirror the
    !    oracle's access pattern (mass_matrix => ice%work%fct_massmatrix) to keep codegen.
    !  - ice_diff=0.0 / ice_gamma_fct=0.5 are set by the driver to the CORE2 namelist values
    !    (NOT the t_ice single->WP defaults 10.0/0.25); the FCT reads ice%ice_diff /
    !    ice%ice_gamma_fct so the override propagates. scale_area only multiplies the (zero)
    !    diffusion, so it cannot move results while ice_diff=0 -- but it is a &oce_dyn
    !    NAMELIST parameter, and the CORE2 reference sets it to 5.8e9, NOT the o_PARAM
    !    default 2.0e8 that an earlier version of this comment called "correct".
    !    read_param_phys() is never called by the lifecycle drivers, so the default stood;
    !    at ice_diff/=0 that would make this diffusion sqrt(5.8e9/2.0e8) = 5.4x the
    !    reference's. The native drivers now pin scale_area = 5.8e9 with the rest of the
    !    CORE2 namelist values.
    !  - hybrid meshes (plan docs/plans/completed/2026-10-08-hybrid-meshes.md, Task 16, the user's
    !    specification): unified 4-slot formulas, so triangles agree with FESOM2 to roundoff
    !    only. The consistent mass matrix is assembled from the element contributions
    !    elem_area*ice_mcoef(:,:,nv) (triangle S/6 diagonal, S/12 off-diagonal; quad, regular-
    !    square simplification, row j (S/36)*(4,2,1,2) rotated); its row sums (S/3 per
    !    triangle vertex, S/4 per quad vertex) are the lumped mass ice%work%fct_lumpedmass,
    !    which replaces the node area in the low-order step, the high-order iteration and
    !    the antidiffusive fluxes (ice mass is conserved as sum(lumped*m)). Antidiffusive
    !    element operator = consistent - lumped (ice_acoef: triangle (-2,1,1)*S/12, quad
    !    (-5,2,1,2)*S/36). TG rhs: advection dx(n)*(M_e u)_q, diffusion over the 4 gradient
    !    slots, TG term with the kite-weighted element-mean velocity sum(elem_w*u). A
    !    triangle's slot 4 (repeated 1st vertex) has coefficient 0 in every table and in
    !    gradient_sca; the limiter loops run over the elem2D_nnodes unique vertices.
    use mod_precision,    only: WP
    use mod_mesh,         only: t_mesh, MAX_NV
    use mod_ice,          only: t_ice
    use mod_partit,       only: t_partit
    use mod_part_bounds,  only: owned_bounds, is_multirank
    use mod_halo,         only: exchange_nod
    use mod_param_phys,   only: scale_area
    implicit none
    private
    public :: ice_TG_rhs, ice_fct_solve
    public :: ice_mcoef, ice_acoef

    ! Element mass contributions per unit element area, ice_mcoef(j,k,nv) (symmetric):
    ! triangle (nv=3) 1/6 diagonal, 1/12 off-diagonal, slot 4 = 0; quad (nv=4) row j
    ! (4,2,1,2)/36 rotated so that the 4 sits on the diagonal. Row sums 1/3 and 1/4.
    real(kind=WP), parameter :: ice_mcoef(MAX_NV, MAX_NV, 3:4) = reshape([ &
        2.0_WP/12.0_WP, 1.0_WP/12.0_WP, 1.0_WP/12.0_WP, 0.0_WP,         &
        1.0_WP/12.0_WP, 2.0_WP/12.0_WP, 1.0_WP/12.0_WP, 0.0_WP,         &
        1.0_WP/12.0_WP, 1.0_WP/12.0_WP, 2.0_WP/12.0_WP, 0.0_WP,         &
        0.0_WP,         0.0_WP,         0.0_WP,         0.0_WP,         &
        4.0_WP/36.0_WP, 2.0_WP/36.0_WP, 1.0_WP/36.0_WP, 2.0_WP/36.0_WP, &
        2.0_WP/36.0_WP, 4.0_WP/36.0_WP, 2.0_WP/36.0_WP, 1.0_WP/36.0_WP, &
        1.0_WP/36.0_WP, 2.0_WP/36.0_WP, 4.0_WP/36.0_WP, 2.0_WP/36.0_WP, &
        2.0_WP/36.0_WP, 1.0_WP/36.0_WP, 2.0_WP/36.0_WP, 4.0_WP/36.0_WP], &
        [MAX_NV, MAX_NV, 2])
    ! Antidiffusive element operator = consistent - lumped (ice_mcoef minus its row sums on
    ! the diagonal): triangle (-2,1,1)/12, quad (-5,2,1,2)/36, rotated per row.
    real(kind=WP), parameter :: ice_acoef(MAX_NV, MAX_NV, 3:4) = reshape([ &
        -2.0_WP/12.0_WP,  1.0_WP/12.0_WP,  1.0_WP/12.0_WP,  0.0_WP,         &
         1.0_WP/12.0_WP, -2.0_WP/12.0_WP,  1.0_WP/12.0_WP,  0.0_WP,         &
         1.0_WP/12.0_WP,  1.0_WP/12.0_WP, -2.0_WP/12.0_WP,  0.0_WP,         &
         0.0_WP,          0.0_WP,          0.0_WP,          0.0_WP,         &
        -5.0_WP/36.0_WP,  2.0_WP/36.0_WP,  1.0_WP/36.0_WP,  2.0_WP/36.0_WP, &
         2.0_WP/36.0_WP, -5.0_WP/36.0_WP,  2.0_WP/36.0_WP,  1.0_WP/36.0_WP, &
         1.0_WP/36.0_WP,  2.0_WP/36.0_WP, -5.0_WP/36.0_WP,  2.0_WP/36.0_WP, &
         2.0_WP/36.0_WP,  1.0_WP/36.0_WP,  2.0_WP/36.0_WP, -5.0_WP/36.0_WP], &
        [MAX_NV, MAX_NV, 2])

contains

    !___________________________________________________________________________
    ! ice_TG_rhs (FESOM2 ice_fct.F90:91). Taylor-Galerkin (Lax-Wendroff) rhs assembled
    ! over elements into the nodal rhs (values_rhs of each ice tracer), with an ice_diff
    ! stabilising diffusion. The advection velocity is uice/vice. Skip cavity elements
    ! (ulevels>1). No exchange in FESOM2 (the rhs at the halo is each rank's own; the
    ! downstream solves exchange). Per element (4 slots): advection dx(n)*(M_e u)_q with
    ! the element mass contributions (FESOM2 S*dx(n)*(sum(u)+u_q)/12 = the triangle case),
    ! diffusion S*diff*grad(n).grad(q), TG term with the kite-weighted mean velocity
    ! (FESOM2 sum(u)/3).
    subroutine ice_TG_rhs(ice, mesh, partit)
        type(t_ice),    intent(inout), target :: ice
        type(t_mesh),   intent(in),    target :: mesh
        type(t_partit), intent(in), optional  :: partit
        real(kind=WP) :: diff, entries(MAX_NV), um, vm, vol, dx(MAX_NV), dy(MAX_NV)
        real(kind=WP) :: mu(MAX_NV), mv(MAX_NV)
        integer       :: n, q, nv, row, elem, elnodes(MAX_NV)
        integer       :: nNodO, nNodL, nEdgeO, nElemO
        real(kind=WP), dimension(:), pointer :: u_ice, v_ice
        real(kind=WP), dimension(:), pointer :: a_ice, m_ice, m_snow
        real(kind=WP), dimension(:), pointer :: rhs_a, rhs_m, rhs_ms

        u_ice  => ice%uice(:)
        v_ice  => ice%vice(:)
        a_ice  => ice%data(1)%values(:)
        m_ice  => ice%data(2)%values(:)
        m_snow => ice%data(3)%values(:)
        rhs_a  => ice%data(1)%values_rhs(:)
        rhs_m  => ice%data(2)%values_rhs(:)
        rhs_ms => ice%data(3)%values_rhs(:)
        call owned_bounds(mesh, nNodO, nNodL, nEdgeO, nElemO, partit)

        do row = 1, nNodO
            rhs_m(row)  = 0._WP
            rhs_a(row)  = 0._WP
            rhs_ms(row) = 0._WP
        end do

        do elem = 1, nElemO          ! assembling rhs over elements
            elnodes = mesh%elem2D_nodes(:, elem)
            ! if cavity element skip it
            if (mesh%ulevels(elem) > 1) cycle
            nv  = mesh%elem2D_nnodes(elem)

            ! derivatives (4 slots; a triangle's slot 4 is 0)
            dx  = mesh%gradient_sca(1:4, elem)
            dy  = mesh%gradient_sca(5:8, elem)
            vol = mesh%elem_area(elem)
            ! element mass contributions applied to the velocity, (M_e u)_q (ice_mcoef is
            ! symmetric), and the kite-weighted element-mean velocity
            do q = 1, MAX_NV
                mu(q) = vol*sum(ice_mcoef(:,q,nv)*u_ice(elnodes))
                mv(q) = vol*sum(ice_mcoef(:,q,nv)*v_ice(elnodes))
            end do
            um  = sum(mesh%elem_w(:,elem)*u_ice(elnodes))
            vm  = sum(mesh%elem_w(:,elem)*v_ice(elnodes))

            ! diffusivity
            diff = ice%ice_diff*sqrt(mesh%elem_area(elem)/scale_area)
            do n = 1, MAX_NV
                row = elnodes(n)
                do q = 1, MAX_NV
                    entries(q) = ice%ice_dt*(dx(n)*mu(q) + dy(n)*mv(q) - &
                                 vol*diff*(dx(n)*dx(q) + dy(n)*dy(q)) - &
                                 0.5_WP*ice%ice_dt*vol*(um*dx(n)+vm*dy(n))*(um*dx(q)+vm*dy(q)))
                end do
                rhs_m(row)  = rhs_m(row)  + sum(entries*m_ice(elnodes))
                rhs_a(row)  = rhs_a(row)  + sum(entries*a_ice(elnodes))
                rhs_ms(row) = rhs_ms(row) + sum(entries*m_snow(elnodes))
            end do
        end do
    end subroutine ice_TG_rhs

    !___________________________________________________________________________
    ! ice_fct_solve (FESOM2 ice_fct.F90:210). Driving routine: high-order (Taylor-
    ! Galerkin) solution -> low-order solution -> Zalesak FCT limiter per tracer.
    ! ice_solve_high_order MUST precede ice_solve_low_order (it uses valuesl as temp).
    subroutine ice_fct_solve(ice, mesh, partit)
        type(t_ice),    intent(inout), target :: ice
        type(t_mesh),   intent(in),    target :: mesh
        type(t_partit), intent(in), optional  :: partit

        call ice_solve_high_order(ice, mesh, partit)   ! uses valuesl as temp storage
        call ice_solve_low_order (ice, mesh, partit)

        call ice_fem_fct(1, ice, mesh, partit)         ! m_ice
        call ice_fem_fct(2, ice, mesh, partit)         ! a_ice
        call ice_fem_fct(3, ice, mesh, partit)         ! m_snow
    end subroutine ice_fct_solve

    !___________________________________________________________________________
    ! ice_solve_low_order (FESOM2 ice_fct.F90:243). Low-order solution: add the diffusive
    ! contribution (consistent - lumped mass matrix acting on the previous field) to the
    ! rhs, with the lumped mass matrix (fct_lumpedmass; FESOM2: the node area) on the lhs.
    ! gamma=ice_gamma_fct. Result in valuesl (a_icel/m_icel/m_snowl), exchanged to the halo.
    subroutine ice_solve_low_order(ice, mesh, partit)
        type(t_ice),    intent(inout), target :: ice
        type(t_mesh),   intent(in),    target :: mesh
        type(t_partit), intent(in), optional  :: partit
        integer       :: row, clo, clo2, cn, location(100)
        integer       :: nNodO, nNodL, nEdgeO, nElemO
        logical       :: lmr
        real(kind=WP) :: gamma
        real(kind=WP), dimension(:), pointer :: a_ice, m_ice, m_snow
        real(kind=WP), dimension(:), pointer :: rhs_a, rhs_m, rhs_ms
        real(kind=WP), dimension(:), pointer :: a_icel, m_icel, m_snowl
        real(kind=WP), dimension(:), pointer :: mass_matrix, lumped

        a_ice       => ice%data(1)%values(:)
        m_ice       => ice%data(2)%values(:)
        m_snow      => ice%data(3)%values(:)
        rhs_a       => ice%data(1)%values_rhs(:)
        rhs_m       => ice%data(2)%values_rhs(:)
        rhs_ms      => ice%data(3)%values_rhs(:)
        a_icel      => ice%data(1)%valuesl(:)
        m_icel      => ice%data(2)%valuesl(:)
        m_snowl     => ice%data(3)%valuesl(:)
        mass_matrix => ice%work%fct_massmatrix(:)
        lumped      => ice%work%fct_lumpedmass(:)
        call owned_bounds(mesh, nNodO, nNodL, nEdgeO, nElemO, partit)
        lmr = is_multirank(partit)

        gamma = ice%ice_gamma_fct        ! Added diffusivity parameter

        associate(ssh_stiff => mesh%ssh_stiff)
        do row = 1, nNodO
            ! if there is cavity no ice fct low order
            if (mesh%ulevels_nod2D(row) > 1) cycle

            clo  = ssh_stiff%rowptr_loc(row)   - ssh_stiff%rowptr_loc(1) + 1
            clo2 = ssh_stiff%rowptr_loc(row+1) - ssh_stiff%rowptr_loc(1)
            cn   = clo2 - clo + 1
            location(1:cn) = ssh_stiff%colind_loc(clo:clo2)
            m_icel(row)  = (rhs_m(row) +gamma*sum(mass_matrix(clo:clo2)* &
                            m_ice(location(1:cn))))/lumped(row) + &
                           (1.0_WP-gamma)*m_ice(row)
            a_icel(row)  = (rhs_a(row) +gamma*sum(mass_matrix(clo:clo2)* &
                            a_ice(location(1:cn))))/lumped(row) + &
                           (1.0_WP-gamma)*a_ice(row)
            m_snowl(row) = (rhs_ms(row)+gamma*sum(mass_matrix(clo:clo2)* &
                            m_snow(location(1:cn))))/lumped(row) + &
                           (1.0_WP-gamma)*m_snow(row)
        end do
        end associate

        ! Low-order solution must be known to neighbours
        if (lmr) then
            call exchange_nod(m_icel,  partit)
            call exchange_nod(a_icel,  partit)
            call exchange_nod(m_snowl, partit)
        end if
    end subroutine ice_solve_low_order

    !___________________________________________________________________________
    ! ice_solve_high_order (FESOM2 ice_fct.F90:347). Taylor-Galerkin high-order solution
    ! by num_iter_solve=3 Jacobi-style iterations of the consistent mass matrix (the
    ! difference correction) preconditioned by the lumped mass (fct_lumpedmass). Result in
    ! dvalues (da_ice/dm_ice/dm_snow), with valuesl used as the per-iteration temp.
    subroutine ice_solve_high_order(ice, mesh, partit)
        type(t_ice),    intent(inout), target :: ice
        type(t_mesh),   intent(in),    target :: mesh
        type(t_partit), intent(in), optional  :: partit
        integer       :: n, clo, clo2, cn, location(100), row
        integer       :: nNodO, nNodL, nEdgeO, nElemO
        integer       :: num_iter_solve = 3
        logical       :: lmr
        real(kind=WP) :: rhs_new
        real(kind=WP), dimension(:), pointer :: rhs_a, rhs_m, rhs_ms
        real(kind=WP), dimension(:), pointer :: a_icel, m_icel, m_snowl
        real(kind=WP), dimension(:), pointer :: da_ice, dm_ice, dm_snow
        real(kind=WP), dimension(:), pointer :: mass_matrix, lumped

        rhs_a       => ice%data(1)%values_rhs(:)
        rhs_m       => ice%data(2)%values_rhs(:)
        rhs_ms      => ice%data(3)%values_rhs(:)
        a_icel      => ice%data(1)%valuesl(:)
        m_icel      => ice%data(2)%valuesl(:)
        m_snowl     => ice%data(3)%valuesl(:)
        da_ice      => ice%data(1)%dvalues(:)
        dm_ice      => ice%data(2)%dvalues(:)
        dm_snow     => ice%data(3)%dvalues(:)
        mass_matrix => ice%work%fct_massmatrix(:)
        lumped      => ice%work%fct_lumpedmass(:)
        call owned_bounds(mesh, nNodO, nNodL, nEdgeO, nElemO, partit)
        lmr = is_multirank(partit)

        ! the first approximation
        do row = 1, nNodO
            if (mesh%ulevels_nod2D(row) > 1) cycle
            dm_ice(row)  = rhs_m(row) /lumped(row)
            da_ice(row)  = rhs_a(row) /lumped(row)
            dm_snow(row) = rhs_ms(row)/lumped(row)
        end do
        if (lmr) then
            call exchange_nod(dm_ice,  partit)
            call exchange_nod(da_ice,  partit)
            call exchange_nod(dm_snow, partit)
        end if

        ! iterate
        associate(ssh_stiff => mesh%ssh_stiff)
        do n = 1, num_iter_solve-1
            do row = 1, nNodO
                if (mesh%ulevels_nod2D(row) > 1) cycle
                clo  = ssh_stiff%rowptr_loc(row)   - ssh_stiff%rowptr_loc(1) + 1
                clo2 = ssh_stiff%rowptr_loc(row+1) - ssh_stiff%rowptr_loc(1)
                cn   = clo2 - clo + 1
                location(1:cn) = ssh_stiff%colind_loc(clo:clo2)
                rhs_new      = rhs_m(row)  - sum(mass_matrix(clo:clo2)*dm_ice(location(1:cn)))
                m_icel(row)  = dm_ice(row) + rhs_new/lumped(row)
                rhs_new      = rhs_a(row)  - sum(mass_matrix(clo:clo2)*da_ice(location(1:cn)))
                a_icel(row)  = da_ice(row) + rhs_new/lumped(row)
                rhs_new      = rhs_ms(row) - sum(mass_matrix(clo:clo2)*dm_snow(location(1:cn)))
                m_snowl(row) = dm_snow(row)+ rhs_new/lumped(row)
            end do

            do row = 1, nNodO
                if (mesh%ulevels_nod2D(row) > 1) cycle
                dm_ice(row)  = m_icel(row)
                da_ice(row)  = a_icel(row)
                dm_snow(row) = m_snowl(row)
            end do

            if (lmr) then
                call exchange_nod(dm_ice,  partit)
                call exchange_nod(da_ice,  partit)
                call exchange_nod(dm_snow, partit)
            end if
        end do
        end associate
    end subroutine ice_solve_high_order

    !___________________________________________________________________________
    ! ice_fem_fct (FESOM2 ice_fct.F90:498). Zalesak flux-corrected transport limiter for
    ! tracer tr_array_id (1=m_ice, 2=a_ice, 3=m_snow): elemental antidiffusive fluxes ->
    ! cluster admissible increments (tmax/tmin) -> sum of positive/negative fluxes per node
    ! -> correction factors (icepplus/icepminus) -> per-element limiting -> update the
    ! solution (valuesl + the limited flux scatter). gamma=ice_gamma_fct.
    subroutine ice_fem_fct(tr_array_id, ice, mesh, partit)
        integer,        intent(in)            :: tr_array_id
        type(t_ice),    intent(inout), target :: ice
        type(t_mesh),   intent(in),    target :: mesh
        type(t_partit), intent(in), optional  :: partit
        integer       :: n, q, nv, elem, elnodes(MAX_NV), row, clo, clo2
        integer       :: nNodO, nNodL, nEdgeO, nElemO
        logical       :: lmr
        real(kind=WP) :: vol, flux, ae, gamma
        real(kind=WP), dimension(:),   pointer :: a_ice, m_ice, m_snow
        real(kind=WP), dimension(:),   pointer :: a_icel, m_icel, m_snowl
        real(kind=WP), dimension(:),   pointer :: da_ice, dm_ice, dm_snow
        real(kind=WP), dimension(:),   pointer :: icepplus, icepminus, tmax, tmin, lumped
        real(kind=WP), dimension(:,:), pointer :: icefluxes

        a_ice     => ice%data(1)%values(:)
        m_ice     => ice%data(2)%values(:)
        m_snow    => ice%data(3)%values(:)
        a_icel    => ice%data(1)%valuesl(:)
        m_icel    => ice%data(2)%valuesl(:)
        m_snowl   => ice%data(3)%valuesl(:)
        da_ice    => ice%data(1)%dvalues(:)
        dm_ice    => ice%data(2)%dvalues(:)
        dm_snow   => ice%data(3)%dvalues(:)
        icefluxes => ice%work%fct_fluxes(:,:)
        icepplus  => ice%work%fct_plus(:)
        icepminus => ice%work%fct_minus(:)
        tmax      => ice%work%fct_tmax(:)
        tmin      => ice%work%fct_tmin(:)
        lumped    => ice%work%fct_lumpedmass(:)
        call owned_bounds(mesh, nNodO, nNodL, nEdgeO, nElemO, partit)
        lmr = is_multirank(partit)

        ! It should coincide with gamma in ice_solve_low_order
        gamma = ice%ice_gamma_fct

        !_______________________________________________________________________
        ! Compute elemental antidiffusive fluxes to nodes
        do n = 1, nNodL
            tmax(n) = 0.0_WP
            tmin(n) = 0.0_WP
        end do
        ! Auxiliary elemental operator (mass matrix - lumped mass matrix) = ice_acoef,
        ! applied at the nv unique vertices (ice_acoef is symmetric)
        do elem = 1, nElemO
            elnodes = mesh%elem2D_nodes(:, elem)
            ! if cavity cycle over
            if (mesh%ulevels(elem) > 1) cycle

            nv  = mesh%elem2D_nnodes(elem)
            vol = mesh%elem_area(elem)
            if (tr_array_id == 1) then
                do q = 1, nv
                    icefluxes(elem,q) = -sum(ice_acoef(:,q,nv)*(gamma*m_ice(elnodes) + &
                                dm_ice(elnodes)))*(vol/lumped(elnodes(q)))
                end do
            end if
            if (tr_array_id == 2) then
                do q = 1, nv
                    icefluxes(elem,q) = -sum(ice_acoef(:,q,nv)*(gamma*a_ice(elnodes) + &
                                da_ice(elnodes)))*(vol/lumped(elnodes(q)))
                end do
            end if
            if (tr_array_id == 3) then
                do q = 1, nv
                    icefluxes(elem,q) = -sum(ice_acoef(:,q,nv)*(gamma*m_snow(elnodes) + &
                                dm_snow(elnodes)))*(vol/lumped(elnodes(q)))
                end do
            end if
        end do

        !_______________________________________________________________________
        ! Cluster min/max
        associate(ssh_stiff => mesh%ssh_stiff)
        if (tr_array_id == 1) then
            do row = 1, nNodO
                if (mesh%ulevels_nod2D(row) > 1) cycle
                clo  = ssh_stiff%rowptr_loc(row)   - ssh_stiff%rowptr_loc(1) + 1
                clo2 = ssh_stiff%rowptr_loc(row+1) - ssh_stiff%rowptr_loc(1)
                tmax(row) = max(maxval(m_icel(ssh_stiff%colind_loc(clo:clo2))), &
                                maxval(m_ice (ssh_stiff%colind_loc(clo:clo2))))
                tmin(row) = min(minval(m_icel(ssh_stiff%colind_loc(clo:clo2))), &
                                minval(m_ice (ssh_stiff%colind_loc(clo:clo2))))
                ! Admissible increments
                tmax(row) = tmax(row) - m_icel(row)
                tmin(row) = tmin(row) - m_icel(row)
            end do
        end if
        if (tr_array_id == 2) then
            do row = 1, nNodO
                if (mesh%ulevels_nod2D(row) > 1) cycle
                clo  = ssh_stiff%rowptr_loc(row)   - ssh_stiff%rowptr_loc(1) + 1
                clo2 = ssh_stiff%rowptr_loc(row+1) - ssh_stiff%rowptr_loc(1)
                tmax(row) = max(maxval(a_icel(ssh_stiff%colind_loc(clo:clo2))), &
                                maxval(a_ice (ssh_stiff%colind_loc(clo:clo2))))
                tmin(row) = min(minval(a_icel(ssh_stiff%colind_loc(clo:clo2))), &
                                minval(a_ice (ssh_stiff%colind_loc(clo:clo2))))
                tmax(row) = tmax(row) - a_icel(row)
                tmin(row) = tmin(row) - a_icel(row)
            end do
        end if
        if (tr_array_id == 3) then
            do row = 1, nNodO
                if (mesh%ulevels_nod2D(row) > 1) cycle
                clo  = ssh_stiff%rowptr_loc(row)   - ssh_stiff%rowptr_loc(1) + 1
                clo2 = ssh_stiff%rowptr_loc(row+1) - ssh_stiff%rowptr_loc(1)
                tmax(row) = max(maxval(m_snowl(ssh_stiff%colind_loc(clo:clo2))), &
                                maxval(m_snow (ssh_stiff%colind_loc(clo:clo2))))
                tmin(row) = min(minval(m_snowl(ssh_stiff%colind_loc(clo:clo2))), &
                                minval(m_snow (ssh_stiff%colind_loc(clo:clo2))))
                tmax(row) = tmax(row) - m_snowl(row)
                tmin(row) = tmin(row) - m_snowl(row)
            end do
        end if
        end associate

        !_______________________________________________________________________
        ! Sums of positive/negative fluxes to node row
        do n = 1, nNodL
            icepplus (n) = 0._WP
            icepminus(n) = 0._WP
        end do

        do elem = 1, nElemO
            if (mesh%ulevels(elem) > 1) cycle
            elnodes = mesh%elem2D_nodes(:, elem)
            do q = 1, mesh%elem2D_nnodes(elem)
                n    = elnodes(q)
                flux = icefluxes(elem,q)
                if (flux > 0) then
                    icepplus(n) = icepplus(n) + flux
                else
                    icepminus(n) = icepminus(n) + flux
                end if
            end do
        end do

        !_______________________________________________________________________
        ! The least upper bound for the correction factors
        do n = 1, nNodO
            if (mesh%ulevels_nod2D(n) > 1) cycle
            flux = icepplus(n)
            if (abs(flux) > 0) then
                icepplus(n) = min(1.0_WP, tmax(n)/max(flux, 1.e-12))
            else
                icepplus(n) = 0._WP
            end if

            flux = icepminus(n)
            if (abs(flux) > 0) then
                icepminus(n) = min(1.0_WP, tmin(n)/min(flux, -1.e-12))
            else
                icepminus(n) = 0._WP
            end if
        end do
        ! pminus and pplus are to be known to neighbouring PE
        if (lmr) then
            call exchange_nod(icepminus, partit)
            call exchange_nod(icepplus,  partit)
        end if

        !_______________________________________________________________________
        ! Limiting
        do elem = 1, nElemO
            if (mesh%ulevels(elem) > 1) cycle
            elnodes = mesh%elem2D_nodes(:, elem)
            ae = 1.0_WP
            do q = 1, mesh%elem2D_nnodes(elem)
                n    = elnodes(q)
                flux = icefluxes(elem,q)
                if (flux >= 0._WP) ae = min(ae, icepplus(n))
                if (flux <  0._WP) ae = min(ae, icepminus(n))
            end do
            icefluxes(elem,:) = ae*icefluxes(elem,:)
        end do

        !_______________________________________________________________________
        ! Update the solution
        if (tr_array_id == 1) then
            do n = 1, nNodO
                if (mesh%ulevels_nod2D(n) > 1) cycle
                m_ice(n) = m_icel(n)
            end do
            do elem = 1, nElemO
                if (mesh%ulevels(elem) > 1) cycle
                elnodes = mesh%elem2D_nodes(:, elem)
                do q = 1, mesh%elem2D_nnodes(elem)
                    n = elnodes(q)
                    m_ice(n) = m_ice(n) + icefluxes(elem,q)
                end do
            end do
        end if
        if (tr_array_id == 2) then
            do n = 1, nNodO
                if (mesh%ulevels_nod2D(n) > 1) cycle
                a_ice(n) = a_icel(n)
            end do
            do elem = 1, nElemO
                if (mesh%ulevels(elem) > 1) cycle
                elnodes = mesh%elem2D_nodes(:, elem)
                do q = 1, mesh%elem2D_nnodes(elem)
                    n = elnodes(q)
                    a_ice(n) = a_ice(n) + icefluxes(elem,q)
                end do
            end do
        end if
        if (tr_array_id == 3) then
            do n = 1, nNodO
                if (mesh%ulevels_nod2D(n) > 1) cycle
                m_snow(n) = m_snowl(n)
            end do
            do elem = 1, nElemO
                if (mesh%ulevels(elem) > 1) cycle
                elnodes = mesh%elem2D_nodes(:, elem)
                do q = 1, mesh%elem2D_nnodes(elem)
                    n = elnodes(q)
                    m_snow(n) = m_snow(n) + icefluxes(elem,q)
                end do
            end do
        end if

        if (lmr) then
            call exchange_nod(m_ice,  partit)
            call exchange_nod(a_ice,  partit)
            call exchange_nod(m_snow, partit)
        end if
    end subroutine ice_fem_fct

end module mod_ice_fct
