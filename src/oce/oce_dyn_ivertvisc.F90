module oce_dyn_ivertvisc
    ! Implicit vertical viscosity for momentum (the vertical-implicit diffusion solve
    ! run AFTER viscosity_filter in the timestep, FESOM2 oce_ale.F90:3874-3876 ->
    ! impl_vert_visc_ale). Transcribed from FESOM2 v2.7.3 oce_ale.F90:3315-3526 (the
    ! non-subcycl, non-vtransp branch — use_ssh_se_subcycl=.false.).
    !
    ! Per element column it assembles and solves a tridiagonal system
    !       (I + dt*[vertical viscosity + vertical advection]) u^{n+1} = rhs
    ! by the Thomas algorithm (forward elimination + back substitution). The system
    ! is written in INCREMENT form: the routine takes the explicit UV_rhs (the post-
    ! viscosity_filter right-hand side) and the old velocity UV, subtracts the implicit
    ! operator applied to UV from the rhs ("model solves for the difference to timestep
    ! N"), solves, and OVERWRITES UV_rhs with the solution. UV itself is NOT modified
    ! here (the velocity update UV += UV_rhs happens later in the timestep). So the
    ! kernel's single output is UV_rhs.
    !
    ! Tridiagonal coefficients per layer nz (interior):
    !   a(nz) = -Av(nz)  /(Z(nz-1)-Z(nz))  *zinv      (sub-diagonal, viscous up-coupling)
    !   c(nz) = -Av(nz+1)/(Z(nz)-Z(nz+1))  *zinv      (super-diagonal, viscous down-coupling)
    !   b(nz) = -a(nz)-c(nz)+1                         (diagonal)
    !   zinv  = dt/(zbar(nz)-zbar(nz+1))               (1/layer-thickness * dt)
    ! plus the vertical-advection upwind update from the implicit vertical velocity
    ! w_i averaged to the prism faces (wu at top face nz, wd at bottom face nz+1):
    !   a += min(0,wu)*zinv ; b += max(0,wu)*zinv ; b -= min(0,wd)*zinv ; c -= max(0,wd)*zinv
    ! That is the upwind FLUX form d(w_i u)/dz: its row sums are (wu-wd)*zinv, so a
    ! uniform u gets the tendency -u*dw_i/dz. For the scalar scheme (momadv_opt==2) that
    ! is consistent -- its explicit part is flux form too, and the u*dw/dz piece is
    ! balanced by the horizontal u*div(u) through continuity. The vector-invariant
    ! scheme (momadv_opt==1, oce_dyn_vinv) is ADVECTIVE form in both directions
    ! (zeta x u + grad KE; Block C does w_e du/dz as d(w_e u)/dz - u dw_e/dz), so the
    ! flux-form operator would leave a spurious u*dw_i/dz wherever the split is active
    ! (CFL_z > wsplit_mincfl, oce_wsplit). For momadv_opt==1 the diagonal is therefore
    ! corrected by
    !   b -= (wu-wd)*zinv            (bottom row: b -= wu*zinv, its wd is dropped)
    ! which turns the operator into the upwind ADVECTIVE form w_i du/dz (inflow-face
    ! velocity), zero row sums, constants preserved implicitly, diagonal dominance kept.
    ! Structurally gated: the ==2 path is untouched, and with use_wsplit=.false. w_i==0.
    ! Net: test/test_ivertvisc.F90 (closed forms for both options); end to end in the
    ! fesom_conserve wsplit / momadv-vinv-wsplit gate configs (FESOM3_WSPLIT).
    ! Boundary rows: the surface row carries the wind-stress flux
    !   ur(top) += zinv*stress_surf(1,elem)/density_0 ;  vr(top) += zinv*stress_surf(2,elem)/density_0
    ! and the bottom row the quadratic bottom drag
    !   friction = -C_d*sqrt(UV(1,bot)^2+UV(2,bot)^2) ;  ur(bot) += zinv*friction*UV(1,bot) ...
    !
    ! INPUTS not yet produced by the dynamics core, passed as explicit dummy arguments
    ! (honest about provenance; deferred to their own milestones):
    !   Av(nl,elem2D)         vertical (eddy) viscosity on elements — from PP/other
    !                         vertical mixing (M2.8). Prescribed analytically for the
    !                         M2.5 byte-gate so the TDMA is gated independently of mixing.
    !   stress_surf(2,elem2D) surface wind stress on elements — from the atmospheric
    !                         forcing (M2.10). Prescribed analytically for the gate.
    ! When M2.8/M2.10 land, the timestep caller sources Av from the mixing work field
    ! and stress_surf from the forcing module; this kernel is unchanged.
    !
    ! Pinned to the pi / reduced-M2 gated path (v1-out / UNGATED branches dropped, to be
    ! re-added with their own gates later):
    !   - use_ssh_se_subcycl=.false.  -> this routine (not impl_vert_visc_ale_vtransp).
    !   - toy_ocean=.false. (pi)      -> the standard quadratic bottom drag friction
    !     (the dbgyre/neverworld2/soufflet constant-C_d toy branches are dropped).
    !   - ldiag_ke=.false. (gate)     -> the ke_wind/ke_drag kinetic-energy diagnostics
    !     are dropped (FESOM3 t_dyn carries no ke_* arrays — M2 diagnostics scope; cf.
    !     compute_vel_rhs which drops ldiag_ke the same way).
    !
    ! BIT-IDENTITY NOTES (the L9 transitive-gate pattern — every operand is already
    ! byte-pinned, so faithful transcription byte-matches like M1.1-M2.4):
    !  - The Thomas solve is a strictly SEQUENTIAL recurrence (forward then backward);
    !    there is NO summation/scatter order ambiguity (unlike the edge-based operators).
    !    Per-element columns are independent. So the only requirement is that every
    !    operand entering each column is byte-identical to FESOM2 — which holds: UV,
    !    UV_rhs (post-viscosity, M2.4-visc-gated), w_i/Av/stress_surf (prescribed, byte-
    !    identical), helem/zbar_e_bot/ulevels/nlevels (geometry-proven), C_d/density_0/dt.
    !  - The many runtime divisors (zinv, the Z/zbar differences, m, b(nz), the Thomas
    !    1/b and 1/m) are byte-identical on both sides (same operands), so the Intel
    !    -no-prec-div reciprocal matches (L7/L10/L14).
    !  - Hybrid meshes: wu/wd = the kite-weighted element mean sum(elem_w*w_i) over the 4
    !    vertex slots (FESOM2: sum over the 3 nodes / 3; triangles equal to roundoff) — a
    !    fixed 4-term sum in index order, so deterministic and order-stable.
    !  - SINGLE-LAYER COLUMN CAVEAT (deferred to M2.11 CORE2): for nlevels(elem)==2 the
    !    "last row" block reads Z_n(nzmax-2)=Z_n(0) and UV(:,nzmax-2)=UV(:,0) — benign
    !    out-of-bounds reads (the values are overwritten by the "first row" block / are
    !    multiplied by a(top)=0). FESOM2's Release oracle tolerates them; -check all would
    !    trap. pi has NO such columns (min nlevels = 5, i.e. >= 4 layers), so the exact
    !    transcription is both byte-identical AND -check all clean here. A mesh with
    !    nlevels==2 shelf columns (CORE2) will need a guard; matched to FESOM2's behaviour
    !    then. nzmin=ulevels=1 on pi (no cavity); cavity (nzmin>1) is transcribed but
    !    ungated until M2.11.
    !  - Local intermediates are real(WP). FESOM2 hard-codes them real(kind=WP); at the DP
    !    anchor WP==MP so reading the MP-typed mesh helem/zbar_e_bot into WP arithmetic is
    !    byte-identical (as in all prior kernels).
    use mod_precision,  only: WP
    use mod_constants,  only: density_0
    use mod_param_phys, only: C_d
    use mod_mesh,       only: t_mesh, MAX_NV
    use mod_dyn,        only: t_dyn
    use mod_partit,     only: t_partit
    use mod_part_bounds, only: owned_bounds
    implicit none
    private
    public :: impl_vert_visc_ale

contains

    subroutine impl_vert_visc_ale(dynamics, mesh, dt, Av, stress_surf, partit)
        ! M2.12c: optional partit -> owned element loop (FESOM2 oce_ale.F90:3350
        ! do elem=1,myDim_elem2D). The per-element Thomas solve is column-local with NO
        ! exchange; UV_rhs is overwritten at owned elements (the halo UV_rhs keeps its
        ! post-viscosity value, exactly as FESOM2).
        type(t_dyn),   intent(inout), target :: dynamics
        type(t_mesh),  intent(in),    target :: mesh
        real(kind=WP), intent(in)            :: dt
        real(kind=WP), intent(in)            :: Av(:,:)          ! (nl,  elem2D) vertical viscosity (M2.8 mixing)
        real(kind=WP), intent(in)            :: stress_surf(:,:) ! (2,   elem2D) surface wind stress (M2.10 forcing)
        type(t_partit), intent(in), optional :: partit
        !______________________________________________________________________
        real(kind=WP) :: a(mesh%nl-1), b(mesh%nl-1), c(mesh%nl-1), ur(mesh%nl-1), vr(mesh%nl-1)
        real(kind=WP) :: cp(mesh%nl-1), up(mesh%nl-1), vp(mesh%nl-1)
        integer       :: nz, elem, nzmin, nzmax, elnodes(MAX_NV)
        integer       :: nNodO, nNodL, nEdgeO, nElemO
        real(kind=WP) :: zinv, m, friction, wu, wd
        real(kind=WP) :: zbar_n(mesh%nl), Z_n(mesh%nl-1)
        real(kind=WP), dimension(:,:,:), pointer :: UV, UV_rhs
        real(kind=WP), dimension(:,:),   pointer :: Wvel_i

        UV     => dynamics%uv
        UV_rhs => dynamics%uv_rhs
        Wvel_i => dynamics%w_i
        call owned_bounds(mesh, nNodO, nNodL, nEdgeO, nElemO, partit)

        !______________________________________________________________________
        do elem = 1, nElemO
            elnodes = mesh%elem2D_nodes(:, elem)
            nzmin   = mesh%ulevels(elem)
            nzmax   = mesh%nlevels(elem)

            !__________________________________________________________________
            ! per-element column interface (zbar_n) and mid (Z_n) depths from the
            ! element thickness helem + partial-cell bottom zbar_e_bot (built bottom-up).
            zbar_n = 0.0_WP
            Z_n    = 0.0_WP
            zbar_n(nzmax)   = mesh%zbar_e_bot(elem)
            Z_n(nzmax-1)    = zbar_n(nzmax) + mesh%helem(nzmax-1,elem)/2.0_WP
            do nz = nzmax-1, nzmin+1, -1
                zbar_n(nz)  = zbar_n(nz+1) + mesh%helem(nz,elem)
                Z_n(nz-1)   = zbar_n(nz)   + mesh%helem(nz-1,elem)/2.0_WP
            end do
            zbar_n(nzmin)   = zbar_n(nzmin+1) + mesh%helem(nzmin,elem)

            !__________________________________________________________________
            ! Operator — regular (interior) rows
            do nz = nzmin+1, nzmax-2
                zinv  = 1.0_WP*dt/(zbar_n(nz)-zbar_n(nz+1))
                a(nz) = -Av(nz,  elem)/(Z_n(nz-1)-Z_n(nz))  *zinv
                c(nz) = -Av(nz+1,elem)/(Z_n(nz)  -Z_n(nz+1))*zinv
                b(nz) = -a(nz)-c(nz)+1.0_WP
                ! update from the vertical advection
                wu = sum(mesh%elem_w(:,elem)*Wvel_i(nz,   elnodes))
                wd = sum(mesh%elem_w(:,elem)*Wvel_i(nz+1, elnodes))
                a(nz) = a(nz)+min(0._WP, wu)*zinv
                b(nz) = b(nz)+max(0._WP, wu)*zinv
                b(nz) = b(nz)-min(0._WP, wd)*zinv
                c(nz) = c(nz)-max(0._WP, wd)*zinv
                ! advective form for the vector-invariant scheme (see header)
                if (dynamics%momadv_opt == 1) b(nz) = b(nz)-(wu-wd)*zinv
            end do
            ! The last row
            zinv       = 1.0_WP*dt/(zbar_n(nzmax-1)-zbar_n(nzmax))
            a(nzmax-1) = -Av(nzmax-1,elem)/(Z_n(nzmax-2)-Z_n(nzmax-1))*zinv
            b(nzmax-1) = -a(nzmax-1)+1.0_WP
            c(nzmax-1) = 0.0_WP
            wu         = sum(mesh%elem_w(:,elem)*Wvel_i(nzmax-1, elnodes))
            a(nzmax-1) = a(nzmax-1)+min(0._WP, wu)*zinv
            b(nzmax-1) = b(nzmax-1)+max(0._WP, wu)*zinv
            if (dynamics%momadv_opt == 1) b(nzmax-1) = b(nzmax-1)-wu*zinv
            ! The first row
            zinv       = 1.0_WP*dt/(zbar_n(nzmin)-zbar_n(nzmin+1))
            c(nzmin)   = -Av(nzmin+1,elem)/(Z_n(nzmin)-Z_n(nzmin+1))*zinv
            a(nzmin)   = 0.0_WP
            b(nzmin)   = -c(nzmin)+1.0_WP
            wu         = sum(mesh%elem_w(:,elem)*Wvel_i(nzmin,   elnodes))
            wd         = sum(mesh%elem_w(:,elem)*Wvel_i(nzmin+1, elnodes))
            b(nzmin)   = b(nzmin)+wu*zinv
            b(nzmin)   = b(nzmin)-min(0._WP, wd)*zinv
            c(nzmin)   = c(nzmin)-max(0._WP, wd)*zinv
            if (dynamics%momadv_opt == 1) b(nzmin)   = b(nzmin)-(wu-wd)*zinv

            !__________________________________________________________________
            ! The rhs
            ur(nzmin:nzmax-1) = UV_rhs(1,nzmin:nzmax-1,elem)
            vr(nzmin:nzmax-1) = UV_rhs(2,nzmin:nzmax-1,elem)
            ! first row: surface wind-stress forcing (zinv is the first-row zinv above)
            ur(nzmin) = ur(nzmin)+zinv*stress_surf(1,elem)/density_0
            vr(nzmin) = vr(nzmin)+zinv*stress_surf(2,elem)/density_0
            ! (ldiag_ke ke_wind diagnostic dropped — see header)
            ! last row: quadratic bottom drag
            zinv     = 1.0_WP*dt/(zbar_n(nzmax-1)-zbar_n(nzmax))
            friction = -C_d*sqrt(UV(1,nzmax-1,elem)**2 + UV(2,nzmax-1,elem)**2)
            ur(nzmax-1) = ur(nzmax-1)+zinv*friction*UV(1,nzmax-1,elem)
            vr(nzmax-1) = vr(nzmax-1)+zinv*friction*UV(2,nzmax-1,elem)
            ! (ldiag_ke ke_drag diagnostic dropped — see header)

            ! solve for the difference to timestep N: subtract the operator applied to UV
            do nz = nzmin+1, nzmax-2
                ur(nz) = ur(nz)-a(nz)*UV(1,nz-1,elem)-(b(nz)-1.0_WP)*UV(1,nz,elem)-c(nz)*UV(1,nz+1,elem)
                vr(nz) = vr(nz)-a(nz)*UV(2,nz-1,elem)-(b(nz)-1.0_WP)*UV(2,nz,elem)-c(nz)*UV(2,nz+1,elem)
            end do
            ur(nzmin)   = ur(nzmin)-(b(nzmin)-1.0_WP)*UV(1,nzmin,elem)-c(nzmin)*UV(1,nzmin+1,elem)
            vr(nzmin)   = vr(nzmin)-(b(nzmin)-1.0_WP)*UV(2,nzmin,elem)-c(nzmin)*UV(2,nzmin+1,elem)
            ur(nzmax-1) = ur(nzmax-1)-a(nzmax-1)*UV(1,nzmax-2,elem)-(b(nzmax-1)-1.0_WP)*UV(1,nzmax-1,elem)
            vr(nzmax-1) = vr(nzmax-1)-a(nzmax-1)*UV(2,nzmax-2,elem)-(b(nzmax-1)-1.0_WP)*UV(2,nzmax-1,elem)

            !__________________________________________________________________
            ! Thomas sweep — forward elimination
            cp(nzmin) = c(nzmin)/b(nzmin)
            up(nzmin) = ur(nzmin)/b(nzmin)
            vp(nzmin) = vr(nzmin)/b(nzmin)
            do nz = nzmin+1, nzmax-1
                m      = b(nz)-cp(nz-1)*a(nz)
                cp(nz) = c(nz)/m
                up(nz) = (ur(nz)-up(nz-1)*a(nz))/m
                vp(nz) = (vr(nz)-vp(nz-1)*a(nz))/m
            end do
            ! back substitution
            ur(nzmax-1) = up(nzmax-1)
            vr(nzmax-1) = vp(nzmax-1)
            do nz = nzmax-2, nzmin, -1
                ur(nz) = up(nz)-cp(nz)*ur(nz+1)
                vr(nz) = vp(nz)-cp(nz)*vr(nz+1)
            end do

            !__________________________________________________________________
            ! write the solution back into UV_rhs (the kernel output)
            do nz = nzmin, nzmax-1
                UV_rhs(1,nz,elem) = ur(nz)
                UV_rhs(2,nz,elem) = vr(nz)
            end do
        end do
    end subroutine impl_vert_visc_ale

end module oce_dyn_ivertvisc
