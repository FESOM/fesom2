module mod_ice_dyn
    ! M3b: sea-ice EVP dynamics + ocean->ice coupling. Faithful transcription of
    ! FESOM2 v2.7.3 ocean2ice (ice_oce_coupling.F90:194) + EVPdynamics / stress_tensor /
    ! stress2rhs (ice_EVP.F90:414/130/261), standard EVP path (whichEVP=0; NO
    ! icepack / meltponds / cavity / oasis-yac / oifs).
    !
    ! Pinned to the reduced-M2/M3 gated path (= the M2 dynamics kernels):
    !   - which_ALE='linfs'  => use_pice=0  -> the EVP elevation term drops the floating-
    !     ice loading p_ice (FESOM2's `if (.not. linfs)` branch is the unported non-linfs
    !     path; the `else` linfs branch is transcribed here).
    !   - whichEVP=0 (standard EVP); use_cavity=.false. (no cavity-edge velocity BC);
    !     ice%ice_update=.true. (ocean2ice copies the instantaneous surface state).
    !
    ! M2.12 optional-partit pattern (= every M2 dynamics kernel): partit absent OR
    ! npes==1 -> the proven 1-rank path VERBATIM (owned_bounds returns global counts,
    ! is_multirank=.false. so every exchange is skipped). present+npes>1 -> owned/halo
    ! loop bounds + the FESOM2 exchanges (exchange_nod u_w/v_w in ocean2ice, exchange_nod
    ! u_ice/v_ice per EVP subcycle). The per-element / per-node arithmetic is UNCHANGED so
    ! codegen — and the byte-match — is preserved; only loop bounds change + guarded
    ! exchanges are added. (Multi-rank ice is gated at M3f; M3b is the CORE2 1-rank gate.)
    !
    ! BIT-IDENTITY NOTES (the L9 transitive-gate pattern — every operand already byte-
    ! pinned, so faithful transcription byte-matches like the M2 dynamics kernels):
    !  - gradient_sca / metric_factor / elem_area / area / coriolis_node are geometry-gate
    !    proven (max|delta|=0, M2.11a); a_ice/m_ice/m_snow come from the cold-start IC
    !    (M3a, gated); the prescribed u_w/elevation/stress_atmice are byte-identical
    !    analytic functions of coord_nod2D (= the driver / oracle shim).
    !  - the per-node stress-divergence accumulation (stress2rhs) over OWNED elements and
    !    the per-node ocean-velocity area-average (ocean2ice) over nod_in_elem2D are FP-
    !    order-sensitive; the 1-rank global element order matches FESOM2's exactly (the
    !    area gate pinned nod_in_elem2D order; M2.11a pinned elem2D_nodes order).
    !  - the runtime divides (1/(area*mass), 1/max(.,9), rhs/area, 1/(r_a^2+r_b^2),
    !    1/max(delta,delta_min)) are byte-identical on both sides (operands geom/IC-proven),
    !    so the -no-prec-div reciprocals agree (L7/L10/L14).
    !  - cd_oce_ice / delta_min are set by the driver to the namelist DOUBLES (0.0055_WP /
    !    1.0e-11_WP) NOT the t_ice single->WP defaults (the M3a precision note); the EVP
    !    arithmetic reads ice%cd_oce_ice / ice%delta_min so the override propagates.
    !  - hybrid meshes (plan docs/plans/completed/2026-10-08-hybrid-meshes.md, Task 15): unified 4-slot
    !    formulas, so triangles agree with FESOM2 to roundoff only. Strain rates and the
    !    stress divergence use gradient_sca over the 4 vertex slots (x 1:4, y 5:8; a
    !    triangle's slot 4 has coefficient 0); every element mean (metric terms, msum/asum)
    !    is the kite-weighted sum(elem_w*f) and every element->node share (metric terms of
    !    stress2rhs, SSH tilt) is elem_area*elem_w(k) (FESOM2: elem_area/3). The shares sum
    !    to the node area (area = sum elem_area*elem_w), and stress2rhs stays the exact
    !    negative adjoint of the strain rates (sum_n u.rhs = -sum_e elem_area*sigma:eps).
    !    The ocean-velocity node gather is kite-share weighted (sum A*w*UV / sum A*w).
    use mod_precision,   only: WP
    use mod_mesh,        only: t_mesh, MAX_NV
    use mod_dyn,         only: t_dyn
    use mod_tracer,      only: t_tracer
    use mod_ice,         only: t_ice
    use mod_partit,      only: t_partit
    use mod_part_bounds, only: owned_bounds, is_multirank
    use mod_halo,        only: exchange_nod
    use mod_constants,   only: density_0, g
    implicit none
    private
    public :: ocean2ice, EVPdynamics, EVPdynamics_m, EVPdynamics_solve
    public :: stress_tensor, stress2rhs          ! for test_ice_dyn_mixed

contains

    !___________________________________________________________________________
    ! ocean2ice (FESOM2 ice_oce_coupling.F90:194). Transmits the ocean surface state
    ! to the ice model: srfoce_temp/salt/ssh = surface tracer / hbar (owned+halo, no
    ! exchange — already valid at the halo), srfoce_u/v = kite-share-weighted node average
    ! (elem_area*nod_in_elem2D_w) of the surface element velocity UV(:,1,:) over the wet
    ! elements of nod_in_elem2D (owned, then exchange_nod).
    ! ice_update=.true. branch (ice_steps_since_upd averaging is the ice_ave_steps>1 path).
    subroutine ocean2ice(ice, dynamics, tracers, mesh, partit)
        type(t_ice),    intent(inout), target :: ice
        type(t_dyn),    intent(in),    target :: dynamics
        type(t_tracer), intent(in),    target :: tracers
        type(t_mesh),   intent(in),    target :: mesh
        type(t_partit), intent(in), optional  :: partit
        integer :: n, elem, k
        integer :: nNodO, nNodL, nEdgeO, nElemO
        real(kind=WP) :: uw, vw, vol, wa
        real(kind=WP), dimension(:,:),   pointer :: temp, salt
        real(kind=WP), dimension(:,:,:), pointer :: UV
        real(kind=WP), dimension(:),     pointer :: u_w, v_w, T_oc_array, S_oc_array, elevation

        temp       => tracers%data(1)%values(:,:)
        salt       => tracers%data(2)%values(:,:)
        UV         => dynamics%uv(:,:,:)
        u_w        => ice%srfoce_u(:)
        v_w        => ice%srfoce_v(:)
        T_oc_array => ice%srfoce_temp(:)
        S_oc_array => ice%srfoce_salt(:)
        elevation  => ice%srfoce_ssh(:)
        call owned_bounds(mesh, nNodO, nNodL, nEdgeO, nElemO, partit)

        ! the arrays in the ice model are renamed (ice_update=.true.): owned+halo, no
        ! exchange (temp/salt/hbar already valid at the halo).
        do n = 1, nNodL
            if (mesh%ulevels_nod2D(n) > 1) cycle
            T_oc_array(n) = temp(1,n)
            S_oc_array(n) = salt(1,n)
            elevation(n)  = mesh%hbar(n)
        end do

        do n = 1, nNodL
            u_w(n) = 0.0_WP
            v_w(n) = 0.0_WP
        end do

        do n = 1, nNodO
            if (mesh%ulevels_nod2D(n) > 1) cycle
            uw  = 0.0_WP
            vw  = 0.0_WP
            vol = 0.0_WP
            do k = 1, mesh%nod_in_elem2D_num(n)
                elem = mesh%nod_in_elem2D(k, n)
                if (mesh%ulevels(elem) > 1) cycle
                wa  = mesh%elem_area(elem)*mesh%nod_in_elem2D_w(k, n)
                vol = vol + wa
                uw  = uw  + UV(1,1,elem)*wa
                vw  = vw  + UV(2,1,elem)*wa
            end do
            uw = uw/vol
            vw = vw/vol
            u_w(n) = uw          ! ice_update=.true.
            v_w(n) = vw
        end do

        if (is_multirank(partit)) then
            call exchange_nod(u_w, partit)
            call exchange_nod(v_w, partit)
        end if
    end subroutine ocean2ice

    !___________________________________________________________________________
    ! stress_tensor (FESOM2 ice_EVP.F90:130). EVP rheology: computes the elemental
    ! stress tensor (sigma11/22/12) from the nodal ice velocity via the deformation-rate
    ! tensor (eps11/22/12) and the elastic-viscous-plastic update. Elements only; skip
    ! cavity (ulevels>1) and ice-free (ice_strength<=0) elements.
    subroutine stress_tensor(ice, mesh, partit)
        type(t_ice),    intent(inout), target :: ice
        type(t_mesh),   intent(in),    target :: mesh
        type(t_partit), intent(in), optional  :: partit
        integer :: el, nNodO, nNodL, nEdgeO, nElemO
        real(kind=WP) :: det1, det2, dte, vale, r1, r2, r3, si1, si2
        real(kind=WP) :: delta, delta_inv, zeta
        real(kind=WP), dimension(:), pointer :: u_ice, v_ice
        real(kind=WP), dimension(:), pointer :: eps11, eps12, eps22
        real(kind=WP), dimension(:), pointer :: sigma11, sigma12, sigma22
        real(kind=WP), dimension(:), pointer :: ice_strength

        u_ice        => ice%uice(:)
        v_ice        => ice%vice(:)
        eps11        => ice%work%eps11(:)
        eps12        => ice%work%eps12(:)
        eps22        => ice%work%eps22(:)
        sigma11      => ice%work%sigma11(:)
        sigma12      => ice%work%sigma12(:)
        sigma22      => ice%work%sigma22(:)
        ice_strength => ice%work%ice_strength(:)
        call owned_bounds(mesh, nNodO, nNodL, nEdgeO, nElemO, partit)

        vale = 1.0_WP/(ice%ellipse**2)
        dte  = ice%ice_dt/(1.0_WP*ice%evp_rheol_steps)
        det1 = 1.0_WP/(1.0_WP + 0.5_WP*ice%Tevp_inv*dte)
        det2 = 1.0_WP/(1.0_WP + 0.5_WP*ice%Tevp_inv*dte) !*ellipse**2

        do el = 1, nElemO
            ! if element contains cavity node skip it
            if (mesh%ulevels(el) > 1) cycle

            if (ice_strength(el) > 0.0_WP) then
                ! deformation rate tensor on element el (4 vertex slots; metric terms with
                ! the kite-weighted element mean)
                eps11(el) = sum(mesh%gradient_sca(1:4,el)*u_ice(mesh%elem2D_nodes(:,el))) &
                    - mesh%metric_factor(el) * sum(mesh%elem_w(:,el)*v_ice(mesh%elem2D_nodes(:,el)))

                eps22(el) = sum(mesh%gradient_sca(5:8,el)*v_ice(mesh%elem2D_nodes(:,el)))

                eps12(el) = 0.5_WP*(sum(mesh%gradient_sca(5:8,el)*u_ice(mesh%elem2D_nodes(:,el))) &
                            + sum(mesh%gradient_sca(1:4,el)*v_ice(mesh%elem2D_nodes(:,el))) &
                            + mesh%metric_factor(el) * sum(mesh%elem_w(:,el)*u_ice(mesh%elem2D_nodes(:,el))))

                ! moduli
                delta = sqrt((eps11(el)*eps11(el) + eps22(el)*eps22(el))*(1.0_WP+vale) &
                             + 4.0_WP*vale*eps12(el)*eps12(el) &
                             + 2.0_WP*eps11(el)*eps22(el)*(1.0_WP-vale))

                ! limit delta from below -> viscosity zeta bounded above
                delta_inv = 1.0_WP/max(delta, ice%delta_min)
                zeta = ice_strength(el)*delta_inv
                zeta = zeta*ice%Tevp_inv

                r1 = zeta*(eps11(el)+eps22(el)) - ice_strength(el)*ice%Tevp_inv
                r2 = zeta*(eps11(el)-eps22(el))*vale
                r3 = zeta*eps12(el)*vale

                si1 = det1*(sigma11(el) + sigma22(el) + dte*r1)
                si2 = det2*(sigma11(el) - sigma22(el) + dte*r2)

                sigma12(el) = det2*(sigma12(el)+dte*r3)
                sigma11(el) = 0.5_WP*(si1+si2)
                sigma22(el) = 0.5_WP*(si1-si2)
            end if
        end do
    end subroutine stress_tensor

    !___________________________________________________________________________
    ! stress2rhs (FESOM2 ice_EVP.F90:261). Divergence of the stress tensor -> the
    ! momentum rhs (u_rhs_ice/v_rhs_ice). Accumulates the elemental stress divergence
    ! into the element vertices (owned elements; 4 slots, metric term with the node share
    ! elem_area*elem_w(k)), then scales by inv_areamass and adds the elevation
    ! contribution (rhs_a/rhs_m). Nodes with no ice mass -> rhs 0.
    subroutine stress2rhs(ice, mesh, partit)
        type(t_ice),    intent(inout), target :: ice
        type(t_mesh),   intent(in),    target :: mesh
        type(t_partit), intent(in), optional  :: partit
        integer :: n, el, k, nNodO, nNodL, nEdgeO, nElemO
        real(kind=WP), dimension(:), pointer :: sigma11, sigma12, sigma22
        real(kind=WP), dimension(:), pointer :: u_rhs_ice, v_rhs_ice, rhs_a, rhs_m
        real(kind=WP), dimension(:), pointer :: inv_areamass, ice_strength

        sigma11      => ice%work%sigma11(:)
        sigma12      => ice%work%sigma12(:)
        sigma22      => ice%work%sigma22(:)
        u_rhs_ice    => ice%uice_rhs(:)
        v_rhs_ice    => ice%vice_rhs(:)
        rhs_a        => ice%data(1)%values_rhs(:)
        rhs_m        => ice%data(2)%values_rhs(:)
        inv_areamass => ice%work%inv_areamass(:)
        ice_strength => ice%work%ice_strength(:)
        call owned_bounds(mesh, nNodO, nNodL, nEdgeO, nElemO, partit)

        do n = 1, nNodO
            u_rhs_ice(n) = 0.0_WP
            v_rhs_ice(n) = 0.0_WP
        end do

        do el = 1, nElemO
            ! if element contains cavity node skip it
            if (mesh%ulevels(el) > 1) cycle

            if (ice_strength(el) > 0.0_WP) then
                ! 4 vertex slots: a triangle's slot 4 (= its 1st vertex) has gradient and
                ! kite weight 0, so it adds exactly 0
                do k = 1, MAX_NV
                    u_rhs_ice(mesh%elem2D_nodes(k,el)) = u_rhs_ice(mesh%elem2D_nodes(k,el)) &
                        - mesh%elem_area(el) * &
                          (sigma11(el)*mesh%gradient_sca(k,el) + sigma12(el)*mesh%gradient_sca(k+4,el) &
                           + sigma12(el)*mesh%elem_w(k,el)*mesh%metric_factor(el))   ! metrics

                    v_rhs_ice(mesh%elem2D_nodes(k,el)) = v_rhs_ice(mesh%elem2D_nodes(k,el)) &
                        - mesh%elem_area(el) * &
                          (sigma12(el)*mesh%gradient_sca(k,el) + sigma22(el)*mesh%gradient_sca(k+4,el) &
                           - sigma11(el)*mesh%elem_w(k,el)*mesh%metric_factor(el))
                end do
            end if
        end do

        do n = 1, nNodO
            ! if cavity node skip it
            if (mesh%ulevels_nod2D(n) > 1) cycle

            if (inv_areamass(n) > 0.0_WP) then
                u_rhs_ice(n) = u_rhs_ice(n)*inv_areamass(n) + rhs_a(n)
                v_rhs_ice(n) = v_rhs_ice(n)*inv_areamass(n) + rhs_m(n)
            else
                u_rhs_ice(n) = 0.0_WP
                v_rhs_ice(n) = 0.0_WP
            end if
        end do
    end subroutine stress2rhs

    !___________________________________________________________________________
    ! EVPdynamics (FESOM2 ice_EVP.F90:414). The EVP subcycling driver: precompute the
    ! per-node inverse area-mass / inverse mass + the elemental ice_strength and the ssh-
    ! gradient contribution to the velocity rhs, then iterate evp_rheol_steps subcycles of
    ! (stress_tensor -> stress2rhs -> implicit velocity update [ocean drag + Coriolis +
    ! stress rhs + wind-on-ice stress] -> coastal BC -> exchange_nod(u_ice/v_ice)).
    subroutine EVPdynamics(ice, mesh, partit)
        type(t_ice),    intent(inout), target :: ice
        type(t_mesh),   intent(in),    target :: mesh
        type(t_partit), intent(in), optional  :: partit
        integer :: shortstep, n, ed, el, k, elnodes(MAX_NV)
        integer :: nNodO, nNodL, nEdgeO, nElemO
        logical :: lmr, is_bnd
        real(kind=WP) :: rdt, asum, msum, r_a, r_b, drag, det, umod, rhsu, rhsv
        real(kind=WP) :: ax, ay, aa, elevation_dx, elevation_dy
        real(kind=WP), dimension(:), pointer :: u_ice, v_ice, a_ice, m_ice, m_snow
        real(kind=WP), dimension(:), pointer :: u_ice_old, v_ice_old
        real(kind=WP), dimension(:), pointer :: u_rhs_ice, v_rhs_ice, rhs_a, rhs_m
        real(kind=WP), dimension(:), pointer :: u_w, v_w, elevation
        real(kind=WP), dimension(:), pointer :: stress_atmice_x, stress_atmice_y
        real(kind=WP), dimension(:), pointer :: inv_areamass, inv_mass, ice_strength
        real(kind=WP) :: rhosno, rhoice

        u_ice           => ice%uice(:)
        v_ice           => ice%vice(:)
        a_ice           => ice%data(1)%values(:)
        m_ice           => ice%data(2)%values(:)
        m_snow          => ice%data(3)%values(:)
        u_ice_old       => ice%uice_old(:)
        v_ice_old       => ice%vice_old(:)
        u_rhs_ice       => ice%uice_rhs(:)
        v_rhs_ice       => ice%vice_rhs(:)
        rhs_a           => ice%data(1)%values_rhs(:)
        rhs_m           => ice%data(2)%values_rhs(:)
        u_w             => ice%srfoce_u(:)
        v_w             => ice%srfoce_v(:)
        elevation       => ice%srfoce_ssh(:)
        stress_atmice_x => ice%stress_atmice_x(:)
        stress_atmice_y => ice%stress_atmice_y(:)
        inv_areamass    => ice%work%inv_areamass(:)
        inv_mass        => ice%work%inv_mass(:)
        ice_strength    => ice%work%ice_strength(:)
        rhosno = ice%thermo%rhosno
        rhoice = ice%thermo%rhoice
        call owned_bounds(mesh, nNodO, nNodL, nEdgeO, nElemO, partit)
        lmr = is_multirank(partit)

        rdt = ice%ice_dt/(1.0_WP*ice%evp_rheol_steps)
        ax  = cos(ice%theta_io)
        ay  = sin(ice%theta_io)

        ! precompute values that are never changed during the iteration
        do n = 1, nNodL
            inv_areamass(n) = 0.0_WP
            inv_mass(n)     = 0.0_WP
            rhs_a(n)        = 0.0_WP
            rhs_m(n)        = 0.0_WP
        end do

        do n = 1, nNodO
            ! if cavity node skip it
            if (mesh%ulevels_nod2D(n) > 1) cycle

            if ((rhoice*m_ice(n)+rhosno*m_snow(n)) > 1.e-3_WP) then
                inv_areamass(n) = 1.0_WP/(mesh%area(n)*(rhoice*m_ice(n)+rhosno*m_snow(n)))
            else
                inv_areamass(n) = 0.0_WP
            end if

            if (a_ice(n) < 0.01_WP) then
                inv_mass(n) = 0.0_WP                     ! Skip if ice is absent
            else
                inv_mass(n) = (rhoice*m_ice(n)+rhosno*m_snow(n))/a_ice(n)
                inv_mass(n) = 1.0_WP/max(inv_mass(n), 9.0_WP)     ! Limit the mass
            end if
            rhs_a(n) = 0.0_WP        ! temporal storage for the ssh contribution
            rhs_m(n) = 0.0_WP
        end do

        ! ice_strength + ssh-gradient rhs. which_ALE='linfs' => use_pice=0, so the
        ! floating-ice loading p_ice term is dropped (the FESOM2 `else` linfs branch).
        do el = 1, nElemO
            ice_strength(el) = 0.0_WP
            elnodes = mesh%elem2D_nodes(:,el)
            ! if element has any cavity node skip it
            if (mesh%ulevels(el) > 1) cycle

            ! (a triangle's repeated 4th slot leaves any() unchanged)
            if (any(m_ice(elnodes) <= 0.0_WP) .or. &
                any(a_ice(elnodes) <= 0.0_WP)) then
                ! There is no ice in elem
                ice_strength(el) = 0.0_WP
            else
                msum = sum(mesh%elem_w(:,el)*m_ice(elnodes))
                asum = sum(mesh%elem_w(:,el)*a_ice(elnodes))

                ! Hunke and Dukowicz c*h*p*
                ice_strength(el) = ice%pstar*msum*exp(-ice%c_pressure*(1.0_WP-asum))
                ice_strength(el) = 0.5_WP*ice_strength(el)

                ! use rhs_m and rhs_a for storing the contribution from elevation; node
                ! share elem_area*elem_w(k) (a triangle's slot 4 has weight 0)
                aa = 9.81_WP*mesh%elem_area(el)
                elevation_dx = sum(mesh%gradient_sca(1:4,el)*elevation(elnodes))
                elevation_dy = sum(mesh%gradient_sca(5:8,el)*elevation(elnodes))

                do k = 1, MAX_NV
                    rhs_a(elnodes(k)) = rhs_a(elnodes(k)) - aa*mesh%elem_w(k,el)*elevation_dx
                    rhs_m(elnodes(k)) = rhs_m(elnodes(k)) - aa*mesh%elem_w(k,el)*elevation_dy
                end do
            end if
        end do

        do n = 1, nNodO
            if (mesh%ulevels_nod2D(n) > 1) cycle
            rhs_a(n) = rhs_a(n)/mesh%area(n)
            rhs_m(n) = rhs_m(n)/mesh%area(n)
        end do

        !_______________________________________________________________________
        ! End of precomputing -> the ice stepping starts
        do shortstep = 1, ice%evp_rheol_steps
            call stress_tensor(ice, mesh, partit)
            call stress2rhs(ice, mesh, partit)

            do n = 1, nNodL
                u_ice_old(n) = u_ice(n)
                v_ice_old(n) = v_ice(n)
            end do

            do n = 1, nNodO
                ! if cavity node skip it
                if (mesh%ulevels_nod2D(n) > 1) cycle

                if (a_ice(n) >= 0.01_WP) then            ! Skip if ice is absent
                    umod = sqrt((u_ice(n)-u_w(n))**2 + (v_ice(n)-v_w(n))**2)
                    drag = ice%cd_oce_ice*umod*density_0*inv_mass(n)

                    rhsu = u_ice(n) + rdt*(drag*(ax*u_w(n) - ay*v_w(n)) + &
                            inv_mass(n)*stress_atmice_x(n) + u_rhs_ice(n))
                    rhsv = v_ice(n) + rdt*(drag*(ax*v_w(n) + ay*u_w(n)) + &
                            inv_mass(n)*stress_atmice_y(n) + v_rhs_ice(n))

                    r_a      = 1.0_WP + ax*drag*rdt
                    r_b      = rdt*(mesh%coriolis_node(n) + ay*drag)
                    det      = 1.0_WP/(r_a*r_a + r_b*r_b)
                    u_ice(n) = det*(r_a*rhsu + r_b*rhsv)
                    v_ice(n) = det*(r_a*rhsv - r_b*rhsu)
                else                                      ! Set velocities to 0 if ice is absent
                    u_ice(n) = 0.0_WP
                    v_ice(n) = 0.0_WP
                end if
            end do

            ! apply coastal sea-ice velocity boundary condition (zero on boundary edges).
            ! use_cavity=.false. -> no cavity-ocean edge BC.
            do ed = 1, nEdgeO
                if (lmr) then
                    is_bnd = (partit%myList_edge2D(ed) > mesh%edge2D_in)
                else
                    is_bnd = (ed > mesh%edge2D_in)
                end if
                if (is_bnd) then
                    u_ice(mesh%edges(1:2,ed)) = 0.0_WP
                    v_ice(mesh%edges(1:2,ed)) = 0.0_WP
                end if
            end do

            if (lmr) then
                call exchange_nod(u_ice, partit)
                call exchange_nod(v_ice, partit)
            end if
        end do
    end subroutine EVPdynamics

    !___________________________________________________________________________
    ! EVPdynamics_solve: dispatch on whichEVP (mirrors FESOM2 ice_timestep
    ! ice_setup_step.F90:209 SELECT CASE). 0=standard EVP, 1=modified EVP (mEVP).
    ! aEVP (2) needs alpha_evp_array/beta_evp_array + find_alpha_field_a — deferred.
    subroutine EVPdynamics_solve(ice, mesh, partit)
        type(t_ice),    intent(inout)         :: ice
        type(t_mesh),   intent(in),    target :: mesh
        type(t_partit), intent(in), optional  :: partit
        select case (ice%whichEVP)
        case (0)
            call EVPdynamics(ice, mesh, partit)
        case (1)
            call EVPdynamics_m(ice, mesh, partit)
        case default
            error stop 'EVPdynamics_solve: only whichEVP=0 (EVP) and 1 (mEVP) ported (aEVP=2 deferred)'
        end select
    end subroutine EVPdynamics_solve

    !___________________________________________________________________________
    ! EVPdynamics_m (FESOM2 ice_maEVP.F90:429). Modified EVP (mEVP, whichEVP=1), the
    ! Bouillon et al. (Ocean Modelling 2013) pseudo-time iteration: the elastic
    ! relaxation uses the stability constants alpha_evp (stress) / beta_evp (velocity)
    ! instead of the standard Tevp_inv, iterating an AUXILIARY velocity u_ice_aux/v_ice_aux
    ! toward the implicit solution; u_ice/v_ice are frozen during the subcycling and
    ! copied from the aux at the end. This transcribes the INLINED form FESOM2 actually
    ! executes (stress_tensor_m + stress2rhs_m + ssh2rhs fused into the loop), NOT the
    ! standalone stress_tensor_m/stress2rhs_m (which are unused in this path).
    !
    ! Pinned to the reduced-M2/M3 path: which_ALE='linfs' => the ssh2rhs floating-ice
    ! p_ice term drops; use_cavity=.false. => no cavity-edge BC. Same optional-partit
    ! pattern. The mEVP velocity BC is the node mask bc_index_nod2D (det*bc_index zeroes
    ! boundary nodes) PLUS the edge zeroing (redundant for no-cavity, kept faithfully).
    subroutine EVPdynamics_m(ice, mesh, partit)
        type(t_ice),    intent(inout), target :: ice
        type(t_mesh),   intent(in),    target :: mesh
        type(t_partit), intent(in), optional  :: partit
        integer :: shortstep, steps, n, ed, el, i, k, row, elnodes(MAX_NV)
        integer :: nNodO, nNodL, nEdgeO, nElemO
        logical :: lmr, is_bnd
        real(kind=WP) :: vale, det1, det2, rdt
        real(kind=WP) :: dx(MAX_NV), dy(MAX_NV), w(MAX_NV), msum, asum, metric
        real(kind=WP) :: eps1, eps2, delta, pressure, vol, aa, bb
        real(kind=WP) :: umod, drag, det, rhsu, rhsv
        real(kind=WP), allocatable :: inv_thickness(:), mass(:), pressure_fac(:)
        logical,       allocatable :: ice_nod(:), ice_el(:)
        real(kind=WP), dimension(:), pointer :: u_ice, v_ice, a_ice, m_ice, m_snow
        real(kind=WP), dimension(:), pointer :: eps11, eps12, eps22, sigma11, sigma12, sigma22
        real(kind=WP), dimension(:), pointer :: u_rhs_ice, v_rhs_ice, rhs_a, rhs_m
        real(kind=WP), dimension(:), pointer :: u_w, v_w, elevation
        real(kind=WP), dimension(:), pointer :: stress_atmice_x, stress_atmice_y
        real(kind=WP), dimension(:), pointer :: u_ice_aux, v_ice_aux
        real(kind=WP) :: rhoice, rhosno

        u_ice           => ice%uice(:)
        v_ice           => ice%vice(:)
        a_ice           => ice%data(1)%values(:)
        m_ice           => ice%data(2)%values(:)
        m_snow          => ice%data(3)%values(:)
        eps11           => ice%work%eps11(:)
        eps12           => ice%work%eps12(:)
        eps22           => ice%work%eps22(:)
        sigma11         => ice%work%sigma11(:)
        sigma12         => ice%work%sigma12(:)
        sigma22         => ice%work%sigma22(:)
        u_rhs_ice       => ice%uice_rhs(:)
        v_rhs_ice       => ice%vice_rhs(:)
        rhs_a           => ice%data(1)%values_rhs(:)
        rhs_m           => ice%data(2)%values_rhs(:)
        u_w             => ice%srfoce_u(:)
        v_w             => ice%srfoce_v(:)
        elevation       => ice%srfoce_ssh(:)
        stress_atmice_x => ice%stress_atmice_x(:)
        stress_atmice_y => ice%stress_atmice_y(:)
        u_ice_aux       => ice%uice_aux(:)
        v_ice_aux       => ice%vice_aux(:)
        rhoice = ice%thermo%rhoice
        rhosno = ice%thermo%rhosno
        call owned_bounds(mesh, nNodO, nNodL, nEdgeO, nElemO, partit)
        lmr = is_multirank(partit)
        allocate(inv_thickness(nNodO), mass(nNodO), ice_nod(nNodO))
        allocate(pressure_fac(nElemO), ice_el(nElemO))

        vale = 1.0_WP/(ice%ellipse**2)
        det2 = 1.0_WP/(1.0_WP+ice%alpha_evp)
        det1 = ice%alpha_evp*det2
        rdt  = ice%ice_dt
        steps = ice%evp_rheol_steps

        u_ice_aux = u_ice        ! Initialize solver variables
        v_ice_aux = v_ice

        ! ssh2rhs (inlined, linfs branch): elevation contribution to rhs_a/rhs_m, node
        ! share elem_area*elem_w(n) (a triangle's slot 4 has weight 0).
        do row = 1, nNodO
            rhs_a(row) = 0.0_WP
            rhs_m(row) = 0.0_WP
        end do
        do el = 1, nElemO
            elnodes = mesh%elem2D_nodes(:,el)
            ! if element has any cavity node skip it
            if (mesh%ulevels(el) > 1) cycle
            vol = mesh%elem_area(el)
            dx  = mesh%gradient_sca(1:4,el)
            dy  = mesh%gradient_sca(5:8,el)
            bb  = g*vol
            aa  = bb*sum(dx*elevation(elnodes))
            bb  = bb*sum(dy*elevation(elnodes))
            do n = 1, MAX_NV
                rhs_a(elnodes(n)) = rhs_a(elnodes(n)) - aa*mesh%elem_w(n,el)
                rhs_m(elnodes(n)) = rhs_m(elnodes(n)) - bb*mesh%elem_w(n,el)
            end do
        end do

        ! precompute thickness (the inverse) + mass (scaled by area) + scale rhs.
        do i = 1, nNodO
            inv_thickness(i) = 0.0_WP
            mass(i)          = 0.0_WP
            ice_nod(i)       = .false.
            ! if cavity node skip it
            if (mesh%ulevels_nod2D(i) > 1) cycle

            if (a_ice(i) >= 0.01_WP) then
                inv_thickness(i) = (rhoice*m_ice(i)+rhosno*m_snow(i))/a_ice(i)
                inv_thickness(i) = 1.0_WP/max(inv_thickness(i), 9.0_WP)   ! Limit the mass

                mass(i) = (m_ice(i)*rhoice+m_snow(i)*rhosno)
                mass(i) = mass(i)/((1.0_WP+mass(i)*mass(i))*mesh%area(i))

                ! scale rhs_a, rhs_m, too.
                rhs_a(i) = rhs_a(i)/mesh%area(i)
                rhs_m(i) = rhs_m(i)/mesh%area(i)

                ice_nod(i) = .true.
            end if
        end do

        ! precompute pressure factor (det2 folded in) + ice_el flag.
        do el = 1, nElemO
            elnodes = mesh%elem2D_nodes(:,el)
            pressure_fac(el) = 0.0_WP
            ice_el(el)       = .false.
            ! if element has any cavity node skip it
            if (mesh%ulevels(el) > 1) cycle

            msum = sum(mesh%elem_w(:,el)*m_ice(elnodes))
            if (msum > 0.01_WP) then
                ice_el(el) = .true.
                asum = sum(mesh%elem_w(:,el)*a_ice(elnodes))
                pressure_fac(el) = det2*ice%pstar*msum*exp(-ice%c_pressure*(1.0_WP-asum))
            end if
        end do

        do row = 1, nNodO
            u_rhs_ice(row) = 0.0_WP
            v_rhs_ice(row) = 0.0_WP
        end do

        !_______________________________________________________________________
        ! mEVP subcycling main loop
        do shortstep = 1, steps
            ! stress_tensor_m + stress2rhs_m fused (per owned element).
            do el = 1, nElemO
                if (mesh%ulevels(el) > 1) cycle
                if (ice_el(el)) then
                    elnodes = mesh%elem2D_nodes(:,el)
                    dx = mesh%gradient_sca(1:4,el)
                    dy = mesh%gradient_sca(5:8,el)
                    w  = mesh%elem_w(:,el)
                    metric = mesh%metric_factor(el)              ! metrics

                    ! deformation rate tensor on element el (metric terms with the
                    ! kite-weighted element mean)
                    eps11(el) = sum(dx*u_ice_aux(elnodes)) - sum(w*v_ice_aux(elnodes))*metric
                    eps22(el) = sum(dy*v_ice_aux(elnodes))
                    eps12(el) = 0.5_WP*(sum(dy*u_ice_aux(elnodes) + dx*v_ice_aux(elnodes)) &
                                    + sum(w*u_ice_aux(elnodes))*metric)

                    eps1 = eps11(el) + eps22(el)
                    eps2 = eps11(el) - eps22(el)
                    delta = sqrt(eps1**2 + vale*(eps2**2 + 4.0_WP*eps12(el)**2))
                    pressure = pressure_fac(el)/(delta + ice%delta_min)

                    sigma12(el) = det1*sigma12(el) +        pressure*eps12(el)*vale
                    sigma11(el) = det1*sigma11(el) + 0.5_WP*pressure*(eps1 - delta + eps2*vale)
                    sigma22(el) = det1*sigma22(el) + 0.5_WP*pressure*(eps1 - delta - eps2*vale)

                    ! scatter the stress divergence into the rhs (owned nodes only; metric
                    ! term with the node share elem_area*elem_w(k); a triangle's slot 4 adds 0)
                    do k = 1, MAX_NV
                        if (elnodes(k) <= nNodO) then
                            u_rhs_ice(elnodes(k)) = u_rhs_ice(elnodes(k)) - mesh%elem_area(el)* &
                                    (sigma11(el)*dx(k)+sigma12(el)*dy(k) + sigma12(el)*w(k)*metric)
                            v_rhs_ice(elnodes(k)) = v_rhs_ice(elnodes(k)) - mesh%elem_area(el)* &
                                    (sigma12(el)*dx(k)+sigma22(el)*dy(k) - sigma11(el)*w(k)*metric)
                        end if
                    end do
                end if
            end do

            ! velocity solve (Coriolis + ocean drag implicit; mEVP beta relaxation).
            do i = 1, nNodO
                if (mesh%ulevels_nod2D(i) > 1) cycle
                if (ice_nod(i)) then
                    u_rhs_ice(i) = u_rhs_ice(i)*mass(i) + rhs_a(i)
                    v_rhs_ice(i) = v_rhs_ice(i)*mass(i) + rhs_m(i)

                    umod = sqrt((u_ice_aux(i)-u_w(i))**2 + (v_ice_aux(i)-v_w(i))**2)
                    drag = rdt*ice%cd_oce_ice*umod*density_0*inv_thickness(i)

                    rhsu = u_ice(i)+drag*u_w(i)+rdt*(inv_thickness(i)*stress_atmice_x(i)+u_rhs_ice(i)) + ice%beta_evp*u_ice_aux(i)
                    rhsv = v_ice(i)+drag*v_w(i)+rdt*(inv_thickness(i)*stress_atmice_y(i)+v_rhs_ice(i)) + ice%beta_evp*v_ice_aux(i)

                    det = ice%bc_index_nod2D(i) / ((1.0_WP+ice%beta_evp+drag)**2 + (rdt*mesh%coriolis_node(i))**2)
                    u_ice_aux(i) = det*((1.0_WP+ice%beta_evp+drag)*rhsu + rdt*mesh%coriolis_node(i)*rhsv)
                    v_ice_aux(i) = det*((1.0_WP+ice%beta_evp+drag)*rhsv - rdt*mesh%coriolis_node(i)*rhsu)
                end if
            end do

            ! coastal sea-ice velocity boundary condition (use_cavity=.false.).
            do ed = 1, nEdgeO
                if (lmr) then
                    is_bnd = (partit%myList_edge2D(ed) > mesh%edge2D_in)
                else
                    is_bnd = (ed > mesh%edge2D_in)
                end if
                if (is_bnd) then
                    u_ice_aux(mesh%edges(1:2,ed)) = 0.0_WP
                    v_ice_aux(mesh%edges(1:2,ed)) = 0.0_WP
                end if
            end do

            if (lmr) then
                call exchange_nod(u_ice_aux, partit)
                call exchange_nod(v_ice_aux, partit)
            end if

            do row = 1, nNodO
                u_rhs_ice(row) = 0.0_WP
                v_rhs_ice(row) = 0.0_WP
            end do
        end do

        ! copy the converged auxiliary velocity into the ice velocity (owned+halo).
        do row = 1, nNodL
            u_ice(row) = u_ice_aux(row)
            v_ice(row) = v_ice_aux(row)
        end do

        deallocate(inv_thickness, mass, ice_nod, pressure_fac, ice_el)
    end subroutine EVPdynamics_m

end module mod_ice_dyn
