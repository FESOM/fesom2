module oce_mo_conv
    ! Convective adjustment + (deferred) Monin-Obukhov / wind-mixing enhancements.
    ! Transcribed from FESOM2 v2.7.3 oce_mo_conv.F90:4-121 (subroutine mo_convect), with
    ! the convective term REPLACED by the ROMS ramp (see DEVIATION below).
    !
    ! M2.8b SCOPE — the convective (static-instability) adjustment, called AFTER the
    ! mixing scheme in the step (FESOM2 oce_ale.F90:3729, FESOM3 mod_step_oce.F90:190/
    ! 195/198 — it runs after KPP, TKE and PP alike, it is NOT selected by mix_scheme).
    ! Where N^2 < 0 the column is statically unstable, so the vertical diffusivity Kv
    ! (nodes) and viscosity Av (elements) are enhanced to mimic convective overturning.
    !
    ! DEVIATION FROM FESOM2 (deliberate). FESOM2 applies a hard step:
    !     Kv(nz,node) = max(Kv(nz,node), instabmix_kv)  where bvfreq(nz,node)    < 0
    !     Av(nz,elem) = max(Av(nz,elem), instabmix_kv)  where any(bvfreq(elnodes) < 0)
    ! — a discontinuous coefficient keyed on the SIGN of N^2, so N^2 = -1e-12 got exactly
    ! the treatment of N^2 = -1e-4. This version uses the ROMS LMD_CONVEC ramp instead
    ! (ROMS/Nonlinear/lmd_vmix.F):
    !     cff    = max(N2, n2ref);  cff = min(1, (n2ref - cff)/n2ref)
    !     nu_sxc = (1 - cff^2)^3                    ! 0 at N2 = 0, 1 at N2 <= n2ref
    !     Kv     = Kv + instabmix_kv*nu_sxc         ! ADDITIVE, as in ROMS
    !     Av     = Av + instabmix_kv*maxval(nu_sxc over the element's vertices)
    ! Continuous in N^2, so no discontinuous coefficient enters the implicit vertical
    ! solve; additive rather than max() so there is no second kink where the incoming
    ! scheme value crosses the floor; maxval over the corner nodes is the continuous
    ! generalisation of FESOM2's any(... < 0). A stable cell gets nu_sxc = 0 EXACTLY, so
    ! stable columns are bit-unchanged. In the saturated limit under PP the net effect is
    ! Kv = 0.010 + 0.1 instead of max(0.010, 0.1) — a 10% increase, not a regime change.
    ! Ramp width instabmix_n2ref (mod_param_phys); see the note there on why the default
    ! is narrower than the ROMS value.
    !
    ! reads:  dyn%work%bvfreq_raw (nl,nod2D) -- N^2 as pressure_bv computed it, BEFORE the
    !         horizontal smooth_nod pass; falls back to dyn%work%bvfreq if a driver did not
    !         allocate the raw copy. This is a SECOND deliberate deviation from FESOM2,
    !         which tests the smoothed field. smooth_nod averages 1/3 own + 2/3 neighbour
    !         patch, so it can flip the SIGN of N^2: a convecting node surrounded by
    !         stratified water comes out stable (convection suppressed where the column
    !         really IS unstable), and a stratified node beside a convecting patch comes
    !         out unstable (convection triggered where it is not) -- the instability region
    !         is displaced outward by one element patch, ~100 km on core2. A static-
    !         instability test must see the true local stratification. ROMS never smooths
    !         bvf anywhere (rho_eos.F computes it pointwise; lmd_vmix.F takes it intent(in)).
    !         GM/Redi and PP/KPP/TKE still read the SMOOTHED dyn%work%bvfreq -- only this
    !         test was repointed.
    !         Also reads mesh level arrays. in/out: dyn%work%Kv (nl,nod2D) +
    !         dyn%work%Av (nl,elem2D), modified IN PLACE (the scheme output is incoming).
    ! params: use_instabmix, instabmix_kv, instabmix_n2ref, use_windmix, windmix_kv,
    !         windmix_nl.
    !
    ! DEFERRED to M2.10 (forcing): the use_momix Monin-Obukhov / TB04 block (FESOM2
    ! oce_mo_conv.F90:34-74 + the Av momix term :114) reads forcing/ice fields not ported
    ! yet (water_flux / heat_flux / stress_node_surf / u_ice/v_ice/a_ice / the mo /
    ! mixlength arrays / mo_length / pmlktmo) — OMITTED here (NOT just guarded; the fields
    ! do not exist). pi runs use_momix=.true. in production, but the M2.8b gate FORCES
    ! use_momix=.false. (pre-forcing). The use_windmix block (uses only params + Kv/Av) IS
    ! transcribed, guarded off (use_windmix=.false. — the FESOM2 default; pi does not set it).
    !
    ! 1-rank: node loops 1..nod2D (FESOM2 myDim+eDim, eDim=0); elem loop 1..elem2D (FESOM2
    ! myDim). No halo exchange (purely local). NOT byte-identical to FESOM2 any more: the
    ! ramp replaces the step by design, so a FESOM2 oracle gate covering mo_convect
    ! (run_step_gate.sh, run_lifecycle_{tke,kpp}*_gate_core2.sh) differs on purpose. Those
    ! comparisons are already invalid on this branch — bottom-at-vertices moved the element
    ! bottom, so no FESOM2 run shares this mesh discretisation. The live net is
    ! tools/run_conserve_pi.sh, which is unaffected: vertical diffusion is conservative for
    ! ANY Kv, so changing the coefficient cannot move the heat/salt/volume budgets.
    ! The maxval runs over all 4 vertex slots (a triangle's repeated 4th slot leaves the
    ! max unchanged). cavity nzmin>1 path transcribed, ungated on pi.
    use mod_precision,  only: WP
    use mod_mesh,       only: t_mesh, MAX_NV
    use mod_dyn,        only: t_dyn
    use mod_param_phys, only: use_instabmix, instabmix_kv, instabmix_n2ref, &
                              use_windmix, windmix_kv, windmix_nl
    use mod_partit,     only: t_partit
    use mod_part_bounds, only: owned_bounds
    implicit none
    private
    public :: mo_convect

contains

    elemental function convec_ramp(n2, n2ref, n2ref_inv) result(f)
        ! ROMS LMD_CONVEC shape factor: 0 for N^2 >= 0, 1 for N^2 <= n2ref (< 0), cubic
        ! in between. n2ref_inv = 1/n2ref is passed in so the inner loops carry no divide
        ! (the variable-divisor SIMD trap the pressure/mixing kernels also avoid).
        real(kind=WP), intent(in) :: n2, n2ref, n2ref_inv
        real(kind=WP)             :: f
        f = max(n2, n2ref)
        f = min(1.0_WP, (n2ref - f)*n2ref_inv)
        f = 1.0_WP - f*f
        f = f*f*f
    end function convec_ramp

    subroutine mo_convect(dyn, mesh, partit)
        ! M2.12c: optional partit -> the FESOM2 bounds (oce_mo_conv.F90): the node pass
        ! runs 1..myDim_nod2D+eDim_nod2D (the elem pass reads bvfreq at the element's
        ! corner nodes, halo-valid from pressure_bv, so no exchange) and the elem pass
        ! runs 1..myDim_elem2D (Av consumed at owned elements).
        type(t_mesh), intent(in)            :: mesh
        type(t_dyn),  intent(inout), target :: dyn
        type(t_partit), intent(in), optional :: partit
        integer       :: node, elem, nz, nzmin, nzmax, elnodes(MAX_NV)
        integer       :: nNodO, nNodL, nEdgeO, nElemO
        real(kind=WP) :: n2ref, n2ref_inv
        real(kind=WP), dimension(:,:), pointer :: Kv, Av, bvf

        Kv     => dyn%work%Kv
        Av     => dyn%work%Av
        ! unsmoothed N^2 where the driver provides it (see header); the fallback keeps
        ! pressure_bv callers that pass no raw copy working unchanged.
        if (allocated(dyn%work%bvfreq_raw)) then
            bvf => dyn%work%bvfreq_raw
        else
            bvf => dyn%work%bvfreq
        end if
        call owned_bounds(mesh, nNodO, nNodL, nEdgeO, nElemO, partit)

        ! the ramp width must be strictly negative; clamp once, never in the inner loop
        n2ref     = min(instabmix_n2ref, -tiny(1.0_WP))
        n2ref_inv = 1.0_WP/n2ref

        !_______________________________________________________________________
        ! enhance vertical DIFFUSIVITY (nodes): convective adjustment + wind mixing.
        do node = 1, nNodL
            nzmax = mesh%nlevels_nod2D(node)
            nzmin = mesh%ulevels_nod2D(node)
            do nz = nzmin+1, nzmax-1
                ! convection: ramp in on static instability (0 exactly where N^2 >= 0)
                if (use_instabmix) Kv(nz,node) = Kv(nz,node) &
                    + instabmix_kv*convec_ramp(bvf(nz,node), n2ref, n2ref_inv)
                ! cavity: no wind mixing
                if (nzmin > 1) cycle
                ! enhanced near-surface wind mixing (FESOM1.4 style)
                if (use_windmix .and. nz <= windmix_nl+1) Kv(nz,node) = max(Kv(nz,node), windmix_kv)
            end do
        end do

        !_______________________________________________________________________
        ! enhance vertical VISCOSITY (elements): convective adjustment + wind mixing.
        ! (the use_momix Monin-Obukhov viscosity term is DEFERRED to M2.10 — see header.)
        do elem = 1, nElemO
            elnodes = mesh%elem2D_nodes(:,elem)
            nzmax = mesh%nlevels(elem)
            nzmin = mesh%ulevels(elem)
            do nz = nzmin+1, nzmax-1
                ! convection: strongest of the element's corner ramps (was any(N^2<0))
                if (use_instabmix) Av(nz,elem) = Av(nz,elem) &
                    + instabmix_kv*maxval(convec_ramp(bvf(nz,elnodes), n2ref, n2ref_inv))
                ! cavity: no Monin-Obukhov / wind mixing
                if (nzmin > 1) cycle
                ! enhanced near-surface wind mixing
                if (use_windmix .and. nz <= windmix_nl+1) Av(nz,elem) = max(Av(nz,elem), windmix_kv)
            end do
        end do
    end subroutine mo_convect

end module oce_mo_conv
