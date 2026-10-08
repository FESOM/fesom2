module mod_ice_oce_coupling
    ! M3e: sea-ice -> ocean coupling-out (THE PAYOFF). Faithful transcription of FESOM2
    ! v2.7.3 ice_oce_coupling.F90 — oce_fluxes_mom (ice-ocean + atm-ocean momentum stress
    ! -> stress_surf, the dynamics surface BC) + oce_fluxes (the air-sea heat/freshwater/salt
    ! budget -> heat_flux/water_flux/virtual_salt/relax_salt, the tracer/SSH surface BCs).
    ! These REPLACE the M2.11c-2 PRESCRIBED flux dump: until now the ocean step read the
    ! surface fluxes from a FESOM2 oracle dump; M3e produces them natively from the M3d ice
    ! thermo (ice%flx_h/flx_fw + the t_atmflux evaporation/prec/runoff) + the EVP ice velocity.
    !
    ! REDUCED CONFIG (matches the reduced-M2 oracle, like every M2/M3 gate): NO __icepack /
    ! use_cavity / use_icebergs / lwiso / use_landice_water / use_age_tracer / __oasis. The
    ! gated path is use_virt_salt=.true. (which_ALE='linfs'), ref_sss_local=.true., l_snow=
    ! .true., open_water_albedo=0. The dens_flux MOC diagnostic (oce_fluxes:691-699) is
    ! DEFERRED — it needs sw_alpha/sw_beta (EOS) and vcpw and feeds only the MOC diagnostic,
    ! neither of which is in M3 scope; it does not affect the gated surface fluxes.
    !
    ! BIT-IDENTITY (the load-bearing pieces):
    !  - integrate_nod_2D is the EXACT FESOM2 gen_support.F90:318 reduction: an explicit
    !    SEQUENTIAL loop over OWNED nodes accumulating data(row)*areasvol(ulevels,row), then
    !    (multi-rank) MPI_AllREDUCE. A sequential do-loop — NOT sum() — so the order matches
    !    FESOM2 bit-for-bit under -fp-model precise; at 1-rank the allreduce is identity.
    !  - ocean_area (the balancing divisor net/ocean_area) is mesh%ocean_area, now built by
    !    the same sequential areasvol(1,n) loop FESOM2 uses (mod_mesh_areas.F90).
    !  - the per-node arithmetic + the left-to-right additive order of the freshwater flux are
    !    transcribed VERBATIM; the /ocean_area divide matches (shared runtime operand, L10).
    !  - hybrid meshes: the node->element stress_surf mean is the kite-weighted one over the
    !    4 vertex slots (sum(elem_w*stress_node_surf); FESOM2: plain /3 over 3 nodes), so
    !    triangles agree with FESOM2 to roundoff only. Area-consistent: sum_e elem_area*
    !    stress_surf == sum_n area*stress_node_surf (area = sum elem_area*elem_w).
    !
    ! M2.12 optional-partit pattern: partit absent OR npes==1 -> the proven 1-rank path
    ! VERBATIM (owned_bounds returns global counts, is_multirank=.false.); present+npes>1 ->
    ! owned/halo loop bounds + the cross-rank allreduce_sum in integrate_nod_2D. The momentum-
    ! stress / surface-flux halo exchanges (a_ice/uice/vice/srfoce/stress at the halo, the
    ! exchange_nod(u_w,v_w) ocean2ice already does upstream) are the M3f multi-rank concern;
    ! at 1-rank nNodL==nod2D so every read is local. The per-node arithmetic is UNCHANGED so
    ! codegen — and the byte-match — is preserved.
    use mod_precision,   only: WP
    use mod_constants,   only: density_0
    use mod_mesh,        only: t_mesh, MAX_NV
    use mod_ice,         only: t_ice
    use mod_ice_thermo,  only: t_atmflux
    use mod_tracer,      only: t_tracer
    use mod_partit,      only: t_partit
    use mod_part_bounds, only: owned_bounds, is_multirank
    use mod_halo,        only: allreduce_sum
    implicit none
    private
    public :: oce_fluxes_mom, oce_fluxes

contains

    !___________________________________________________________________________
    ! oce_fluxes_mom (FESOM2 ice_oce_coupling.F90:53). Total surface momentum stress
    ! (ice-ocean drag blended with atm-ocean stress by ice concentration) on nodes
    ! (stress_node_surf) -> averaged to elements (stress_surf, the dynamics surface BC).
    subroutine oce_fluxes_mom(ice, atm, stress_surf, mesh, partit)
        ! mesh declared BEFORE stress_surf so its bound mesh%elem2D resolves (Intel #6415, L27).
        type(t_mesh),    intent(in),    target :: mesh
        type(t_ice),     intent(inout), target :: ice
        type(t_atmflux), intent(inout), target :: atm
        real(kind=WP),   intent(inout)         :: stress_surf(2, mesh%elem2D)
        type(t_partit),  intent(in),    optional :: partit
        integer :: n, elem, elnodes(MAX_NV), nNodO, nNodL, nEdgeO, nElemO
        real(kind=WP) :: aux
        real(kind=WP), dimension(:),   pointer :: u_ice, v_ice, a_ice, u_w, v_w
        real(kind=WP), dimension(:),   pointer :: stress_iceoce_x, stress_iceoce_y
        real(kind=WP), dimension(:),   pointer :: stress_atmoce_x, stress_atmoce_y
        real(kind=WP), dimension(:,:), pointer :: stress_node_surf

        call owned_bounds(mesh, nNodO, nNodL, nEdgeO, nElemO, partit)
        u_ice           => ice%uice(:)
        v_ice           => ice%vice(:)
        a_ice           => ice%data(1)%values(:)
        u_w             => ice%srfoce_u(:)
        v_w             => ice%srfoce_v(:)
        stress_iceoce_x => ice%stress_iceoce_x(:)
        stress_iceoce_y => ice%stress_iceoce_y(:)
        stress_atmoce_x => atm%stress_atmoce_x(:)
        stress_atmoce_y => atm%stress_atmoce_y(:)
        stress_node_surf => atm%stress_node_surf(:,:)

        !_______________________________________________________________________
        ! total surface stress (iceoce+atmoce) on nodes (owned+halo)
        do n = 1, nNodL
            ! if cavity node skip it
            if (mesh%ulevels_nod2D(n) > 1) cycle

            if (a_ice(n) > 0.001_WP) then
                aux = sqrt((u_ice(n)-u_w(n))**2 + (v_ice(n)-v_w(n))**2)*density_0*ice%cd_oce_ice
                stress_iceoce_x(n) = aux * (u_ice(n)-u_w(n))
                stress_iceoce_y(n) = aux * (v_ice(n)-v_w(n))
            else
                stress_iceoce_x(n) = 0.0_WP
                stress_iceoce_y(n) = 0.0_WP
            end if

            stress_node_surf(1,n) = stress_iceoce_x(n)*a_ice(n) + stress_atmoce_x(n)*(1.0_WP-a_ice(n))
            stress_node_surf(2,n) = stress_iceoce_y(n)*a_ice(n) + stress_atmoce_y(n)*(1.0_WP-a_ice(n))
        end do

        !_______________________________________________________________________
        ! total surface stress (iceoce+atmoce) on elements (owned), kite-weighted mean
        do elem = 1, nElemO
            ! if cavity element skip it
            if (mesh%ulevels(elem) > 1) cycle

            elnodes = mesh%elem2D_nodes(:,elem)
            stress_surf(1,elem) = sum(mesh%elem_w(:,elem)*stress_node_surf(1,elnodes))
            stress_surf(2,elem) = sum(mesh%elem_w(:,elem)*stress_node_surf(2,elnodes))
        end do
    end subroutine oce_fluxes_mom

    !___________________________________________________________________________
    ! oce_fluxes (FESOM2 ice_oce_coupling.F90:289). Air-sea heat/freshwater/salt budget.
    ! Heat/water flux from the ice thermo (sign-flipped), virtual salt + SSS relaxation
    ! (both globally balanced to zero net), and the freshwater flux balanced into water_flux.
    subroutine oce_fluxes(ice, tracers, atm, mesh, partit)
        type(t_ice),     intent(inout), target :: ice
        type(t_tracer),  intent(in),    target :: tracers
        type(t_atmflux), intent(inout), target :: atm
        type(t_mesh),    intent(in),    target :: mesh
        type(t_partit),  intent(in),    optional :: partit
        integer :: n, nNodO, nNodL, nEdgeO, nElemO
        real(kind=WP) :: rsss, net
        real(kind=WP), allocatable :: flux(:)
        real(kind=WP), dimension(:,:), pointer :: salt
        real(kind=WP), dimension(:),   pointer :: a_ice_old
        real(kind=WP), dimension(:),   pointer :: fresh_wa_flux, net_heat_flux
        real(kind=WP), dimension(:),   pointer :: heat_flux, water_flux, heat_flux_in
        real(kind=WP), dimension(:),   pointer :: virtual_salt, relax_salt, Ssurf
        real(kind=WP), dimension(:),   pointer :: evaporation, ice_sublimation
        real(kind=WP), dimension(:),   pointer :: prec_rain, prec_snow, runoff
        real(kind=WP),                 pointer :: rhoice, rhosno, inv_rhowat
        real(kind=WP), dimension(:),   pointer :: thdgr, thdgrsn

        call owned_bounds(mesh, nNodO, nNodL, nEdgeO, nElemO, partit)

        salt          => tracers%data(2)%values(:,:)
        a_ice_old     => ice%data(1)%values_old(:)
        fresh_wa_flux => ice%flx_fw(:)
        net_heat_flux => ice%flx_h(:)
        rhoice        => ice%thermo%rhoice
        rhosno        => ice%thermo%rhosno
        inv_rhowat    => ice%thermo%inv_rhowat
        thdgr         => ice%thermo%thdgr(:)
        thdgrsn       => ice%thermo%thdgrsn(:)
        heat_flux       => atm%heat_flux(:)
        water_flux      => atm%water_flux(:)
        heat_flux_in    => atm%heat_flux_in(:)
        virtual_salt    => atm%virtual_salt(:)
        relax_salt      => atm%relax_salt(:)
        Ssurf           => atm%Ssurf(:)
        evaporation     => atm%evaporation(:)
        ice_sublimation => atm%ice_sublimation(:)
        prec_rain       => atm%prec_rain(:)
        prec_snow       => atm%prec_snow(:)
        runoff          => atm%runoff(:)

        allocate(flux(nNodL))
        do n = 1, nNodL
            flux(n) = 0.0_WP
        end do

        !_______________________________________________________________________
        ! heat and freshwater flux (standard, no __icepack)
        do n = 1, nNodL
            heat_flux(n)  = -net_heat_flux(n)
            water_flux(n) = -fresh_wa_flux(n)
        end do

        !_______________________________________________________________________
        ! save total heat flux (heat_flux_in) since heat_flux will be altered by sw_pene
        do n = 1, nNodL
            heat_flux_in(n) = heat_flux(n)
        end do

        !_______________________________________________________________________
        ! balance virtual salt flux (use_virt_salt -> linfs)
        if (atm%use_virt_salt) then
            rsss = atm%ref_sss
            do n = 1, nNodL
                if (atm%ref_sss_local) rsss = salt(mesh%ulevels_nod2D(n), n)
                virtual_salt(n) = rsss*water_flux(n)
            end do

            call integrate_nod_2D(virtual_salt, net, mesh, partit)
            net = net/mesh%ocean_area
            do n = 1, nNodL
                if (mesh%ulevels_nod2D(n) > 1) cycle  ! cavity node
                virtual_salt(n) = virtual_salt(n) - net
            end do
        end if

        !_______________________________________________________________________
        ! balance SSS restoring to climatology (no cavity)
        do n = 1, nNodL
            relax_salt(n) = atm%surf_relax_S*(Ssurf(n) - salt(mesh%ulevels_nod2D(n), n))
        end do
        call integrate_nod_2D(relax_salt, net, mesh, partit)
        net = net/mesh%ocean_area
        do n = 1, nNodL
            if (mesh%ulevels_nod2D(n) > 1) cycle  ! cavity node
            relax_salt(n) = relax_salt(n) - net
        end do

        !_______________________________________________________________________
        ! enforce the total freshwater flux be zero. Standard path (no __icepack):
        ! snow scaled by (1-a_ice_old) (previous-step ice concentration). The additive
        ! order is transcribed VERBATIM (-fp-model precise preserves it).
        do n = 1, nNodL
            flux(n) = evaporation(n)                      &
                      -ice_sublimation(n)                 &
                      +prec_rain(n)                       &
                      +prec_snow(n)*(1.0_WP-a_ice_old(n)) &
                      +runoff(n)
        end do

        ! levitating sea ice (zlevel/zstar, .not. use_virt_salt): add the thermodynamic
        ! growth rates. SKIPPED here (use_virt_salt -> linfs balances it via virtual_salt).
        if (.not. atm%use_virt_salt) then
            do n = 1, nNodL
                flux(n) = flux(n) - thdgr(n)*rhoice*inv_rhowat - thdgrsn(n)*rhosno*inv_rhowat
            end do
        end if

        !_______________________________________________________________________
        ! compute total global net freshwater flux, balance it into water_flux (no cavity).
        ! '+' because water_flux = -fresh_wa_flux flipped the sign but evap/prec/runoff keep theirs.
        call integrate_nod_2D(flux, net, mesh, partit)
        net = net/mesh%ocean_area
        do n = 1, nNodL
            water_flux(n) = water_flux(n) + net
        end do

        !_______________________________________________________________________
        ! dens_flux (MOC diagnostic, oce_fluxes:691-699) DEFERRED — needs sw_alpha/sw_beta
        ! (EOS) + vcpw and feeds only the MOC diagnostic; out of M3 scope, not a surface BC.

        deallocate(flux)
    end subroutine oce_fluxes

    !___________________________________________________________________________
    ! integrate_nod_2D (FESOM2 gen_support.F90:318). Global node integral with the
    ! surface scalar-cell area weight: sum over OWNED nodes of data*areasvol(ulevels,n),
    ! then cross-rank MPI_SUM. Explicit sequential loop (not sum()) so the accumulation
    ! order is byte-identical to FESOM2; at 1-rank the reduction collapses to the local sum.
    subroutine integrate_nod_2D(data, int2D, mesh, partit)
        real(kind=WP),  intent(in)            :: data(:)
        real(kind=WP),  intent(out)           :: int2D
        type(t_mesh),   intent(in)            :: mesh
        type(t_partit), intent(in), optional  :: partit
        integer :: row, nNodO, nNodL, nEdgeO, nElemO
        real(kind=WP) :: lval

        call owned_bounds(mesh, nNodO, nNodL, nEdgeO, nElemO, partit)
        lval = 0.0_WP
        do row = 1, nNodO
            lval = lval + data(row)*mesh%areasvol(row)
        end do
        int2D = lval
        if (is_multirank(partit)) call allreduce_sum(int2D, partit)
    end subroutine integrate_nod_2D

end module mod_ice_oce_coupling
