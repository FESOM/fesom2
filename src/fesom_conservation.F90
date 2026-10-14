!===============================================================================
! Budgets of the ocean's volume, heat and salt against the fluxes the model
! applies at its surface, on the model grid. Switched on with conservation_freq=N
! in &run_config; every N steps rank 0 prints
!
!   CONS <step> <name> <integral> <change since start> <flux since start> <residual> <residual/|change|>
!
! with
!   volume  V = sum areasvol * hnode                          [m3]
!   heat    H = vcpw * sum T * hnode * areasvol               [J]
!   salt    S = sum (S + S_ref_anomaly) * hnode * areasvol    [psu m3]
!
! The fluxes are integrated step by step with the expressions the model uses:
! compute_ssh_rhs_ale for the volume, bc_surface and the shortwave penetration
! in diff_ver_part_impl_ale for heat and salt, relax_to_clim for the nudging. The
! residual is therefore what the discretisation does not account for.
!
! Probes inside the tracer step (conservation_probe) give the thickness-weighted
! integral of temperature and salinity after each stage; their differences are
! printed as CONSP lines and attribute a residual to a stage.
!
! All sums are WP_full, in a fixed order over the owned nodes (no threads), with
! one MPI reduction, as in the monitor record.
!
! Not accounted for (they appear in the residual): icebergs, hosing, the KPP
! non-local salt term, transient tracers, and anything a coupled model adds
! outside the surface-flux arrays.
!===============================================================================
module fesom_conservation_module
    use o_PARAM,          only: WP, WP_full, vcpw, S_ref_anomaly, clim_relax
    use o_ARRAYS,         only: heat_flux, water_flux, virtual_salt, relax_salt, &
                                Tclim, Sclim, relax2clim, is_nonlinfs
    use g_forcing_arrays, only: sw_3d, real_salt_flux
    use g_config,         only: dt, use_sw_pene, toy_ocean, conservation_freq, use_cavity_fw2press
    use MOD_MESH
    use MOD_PARTIT
    use MOD_TRACER

    implicit none
    private
    public :: conservation_before_step, conservation_after_step, conservation_probe

    integer, parameter :: NQ = 3
    character(len=6), parameter :: qname(NQ) = ['volume', 'heat  ', 'salt  ']

    ! Integrals at the start of the run, and the local parts of the surface fluxes
    ! and of the nudging accumulated since then.
    logical,            save :: started = .false.
    real(kind=WP_full), save :: i0(NQ)    = 0.0_WP_full
    real(kind=WP_full), save :: flux(NQ)  = 0.0_WP_full
    real(kind=WP_full), save :: relax(NQ) = 0.0_WP_full

    ! Probes: global integral of tracers 1 and 2 at each named point of the last
    ! tracer step, in the order of the first step.
    integer, parameter :: NPROBE = 12
    character(len=12),  save :: pname(NPROBE) = ''
    real(kind=WP_full), save :: pval(NPROBE, 2) = 0.0_WP_full
    integer,            save :: np = 0

contains
!
!
!===============================================================================
! Before the ocean step, with this step's surface fluxes in place: record the
! start integrals on the first call, then add this step's fluxes.
subroutine conservation_before_step(partit, tracers, mesh)
    type(t_partit), intent(inout), target :: partit
    type(t_tracer), intent(in),    target :: tracers
    type(t_mesh),   intent(in),    target :: mesh
    real(kind=WP_full) :: a, fv, fh, fs, ts, sw, dtf
    integer            :: n, nz, nzmin, nzmax

    if (.not. started) then
        call integrals(tracers, partit, mesh, i0)
        started = .true.
    end if
    dtf = real(dt, WP_full)
    fv = 0.0_WP_full; fh = 0.0_WP_full; fs = 0.0_WP_full
    do n = 1, partit%myDim_nod2D
        nzmin = mesh%ulevels_nod2D(n)
        nzmax = mesh%nlevels_nod2D(n) - 1
        a  = real(mesh%areasvol(nzmin, n), WP_full)
        ts = real(tracers%data(1)%values(nzmin, n), WP_full)
        ! volume: under an ice-shelf cavity melt water enters only with use_cavity_fw2press
        if (nzmin == 1 .or. use_cavity_fw2press) fv = fv - dtf*real(water_flux(n), WP_full)*a
        ! heat, in K m3 until the end of the loop
        fh = fh - dtf*(real(heat_flux(n), WP_full)/real(vcpw, WP_full) &
                  + ts*real(water_flux(n), WP_full)*real(is_nonlinfs, WP_full))*a
        if (use_sw_pene .and. .not. toy_ocean) then
            sw = 0.0_WP_full
            do nz = nzmin, nzmax
                sw = sw + real(sw_3d(nz, n), WP_full)*real(mesh%areasvol(nz, n), WP_full) &
                        - real(sw_3d(nz+1, n), WP_full)*real(mesh%area(nz+1, n), WP_full)
            end do
            fh = fh + dtf*sw
        end if
        fs = fs + dtf*(real(virtual_salt(n), WP_full) + real(relax_salt(n), WP_full) &
                  + (real(real_salt_flux(n), WP_full) + real(S_ref_anomaly, WP_full)*real(water_flux(n), WP_full)) &
                    *real(is_nonlinfs, WP_full))*a
    end do
    flux(1) = flux(1) + fv
    flux(2) = flux(2) + fh*real(vcpw, WP_full)
    flux(3) = flux(3) + fs
end subroutine conservation_before_step
!
!
!===============================================================================
! After the ocean step: add this step's nudging, reconstructed from the state it
! produced, and every conservation_freq steps print the budgets.
subroutine conservation_after_step(istep, partit, tracers, mesh)
    integer,        intent(in)            :: istep
    type(t_partit), intent(inout), target :: partit
    type(t_tracer), intent(in),    target :: tracers
    type(t_mesh),   intent(in),    target :: mesh
    real(kind=WP_full) :: inow(NQ), l(2*NQ), g(2*NQ), change, fl, res, rdt, w
    integer            :: n, nz, k, ierr

    if (.not. started) return
    ! relax_to_clim sets T = T_b + r dt (Tclim - T_b), so T_b = (T - r dt Tclim)/(1 - r dt)
    if (clim_relax > 1.0e-8_WP) then
        do n = 1, partit%myDim_nod2D
            rdt = real(relax2clim(n), WP_full)*real(dt, WP_full)
            if (rdt == 0.0_WP_full) cycle
            do nz = mesh%ulevels_nod2D(n), mesh%nlevels_nod2D(n) - 1
                w = real(mesh%areasvol(nz, n), WP_full)*real(mesh%hnode(nz, n), WP_full)
                relax(2) = relax(2) + w*real(vcpw, WP_full)*rdt*(real(Tclim(nz, n), WP_full) &
                         - real(tracers%data(1)%values(nz, n), WP_full))/(1.0_WP_full - rdt)
                relax(3) = relax(3) + w*rdt*(real(Sclim(nz, n), WP_full) &
                         - real(tracers%data(2)%values(nz, n), WP_full))/(1.0_WP_full - rdt)
            end do
        end do
    end if
    if (mod(istep, conservation_freq) /= 0) return

    call integrals(tracers, partit, mesh, inow)
    l(1:NQ) = flux
    l(NQ+1:2*NQ) = relax
    call MPI_Allreduce(l, g, 2*NQ, MPI_WP_FULL, MPI_SUM, partit%MPI_COMM_FESOM, ierr)
    if (partit%mype /= 0) return
    ! stage contributions of the last tracer step: temperature in J, salinity in psu m3
    do k = 2, np
        write(*, '(" CONSP ",I8,1X,A12,2(1X,ES23.15E3))') istep, pname(k), &
            (pval(k, 1) - pval(k-1, 1))*real(vcpw, WP_full), pval(k, 2) - pval(k-1, 2)
    end do
    do k = 1, NQ
        change = inow(k) - i0(k)
        fl     = g(k) + g(NQ+k)
        res    = change - fl
        write(*, '(" CONS ",I8,1X,A6,5(1X,ES23.15E3))') istep, qname(k), inow(k), change, fl, res, &
            res/max(abs(change), tiny(1.0_WP_full))
    end do
end subroutine conservation_after_step
!
!
!===============================================================================
! Global integral of tracer tr_num (1 or 2) at a named point of the tracer step.
! Phase 'a', before the division by the new thickness: sum areasvol*(hnode*T + del_ttf);
! phase 'b', after it: sum areasvol*hnode_new*T.
subroutine conservation_probe(name, phase, tr_num, tracers, partit, mesh)
    character(len=*), intent(in)            :: name
    character(len=1), intent(in)            :: phase
    integer,          intent(in)            :: tr_num
    type(t_tracer),   intent(in),    target :: tracers
    type(t_partit),   intent(inout), target :: partit
    type(t_mesh),     intent(in),    target :: mesh
    real(kind=WP_full) :: l, g, a
    integer            :: n, nz, k, ierr

    if (conservation_freq <= 0 .or. tr_num > 2) return
    k = 0
    do n = 1, np
        if (pname(n) == name) k = n
    end do
    if (k == 0) then
        if (np >= NPROBE) return
        np = np + 1
        k  = np
        pname(k) = name
    end if
    l = 0.0_WP_full
    do n = 1, partit%myDim_nod2D
        do nz = mesh%ulevels_nod2D(n), mesh%nlevels_nod2D(n) - 1
            a = real(mesh%areasvol(nz, n), WP_full)
            if (phase == 'a') then
                l = l + a*(real(mesh%hnode(nz, n), WP_full)*real(tracers%data(tr_num)%values(nz, n), WP_full) &
                         + real(tracers%work%del_ttf(nz, n), WP_full))
            else
                l = l + a*real(mesh%hnode_new(nz, n), WP_full)*real(tracers%data(tr_num)%values(nz, n), WP_full)
            end if
        end do
    end do
    call MPI_Allreduce(l, g, 1, MPI_WP_FULL, MPI_SUM, partit%MPI_COMM_FESOM, ierr)
    pval(k, tr_num) = g
end subroutine conservation_probe
!
!
!===============================================================================
! Volume, heat and salt of the ocean over all ranks.
subroutine integrals(tracers, partit, mesh, g)
    type(t_tracer),     intent(in),    target :: tracers
    type(t_partit),     intent(inout), target :: partit
    type(t_mesh),       intent(in),    target :: mesh
    real(kind=WP_full), intent(out)           :: g(NQ)
    real(kind=WP_full) :: l(NQ), w
    integer            :: n, nz, ierr

    l = 0.0_WP_full
    do n = 1, partit%myDim_nod2D
        do nz = mesh%ulevels_nod2D(n), mesh%nlevels_nod2D(n) - 1
            w = real(mesh%areasvol(nz, n), WP_full)*real(mesh%hnode(nz, n), WP_full)
            l(1) = l(1) + w
            l(2) = l(2) + w*real(tracers%data(1)%values(nz, n), WP_full)
            l(3) = l(3) + w*(real(tracers%data(2)%values(nz, n), WP_full) + real(S_ref_anomaly, WP_full))
        end do
    end do
    l(2) = l(2)*real(vcpw, WP_full)
    call MPI_Allreduce(l, g, NQ, MPI_WP_FULL, MPI_SUM, partit%MPI_COMM_FESOM, ierr)
end subroutine integrals

end module fesom_conservation_module
