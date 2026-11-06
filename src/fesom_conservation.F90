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
! In OpenACC builds the host copies are not current between the stages, so the probes are off
! there (no CONSP lines); the budgets use the end-of-step state, which is on the host.
!
! All sums are WP_full, in a fixed order over the owned nodes (no threads), with
! one MPI reduction, as in the monitor record.
!
! Not accounted for (they appear in the residual): icebergs, hosing, the KPP
! non-local salt term, transient tracers, and anything a coupled model adds
! outside the surface-flux arrays.
!
! With REcoM (__recom) also DIC, alkalinity and total carbon (all REcoM carbon
! pools, organic and calcite, plus the benthic carbon and calcite), in mmol:
!
!   CONSB  <step> <name> <inventory> <change> <terms> <residual> <residual/inventory> <floor> <ulp>
!   CONSBT <step> <name> <term> <integral of the term since start>
!   CONSPB <step> <stage> <DIC> <Alk>   stage contributions in the last tracer step
!
! inventory = sum areasvol*(values*hnode + carry) over the owned nodes and wet
! levels. It is summed as a double word (two real64), so that the change since the
! start is not limited by the spacing of one real64 number at the inventory (ulp,
! 1.1e-16 of it). What is left is the rounding of the cell contents themselves:
! floor = eps*sqrt(sum cell content**2), the noise of one evaluation. The model's
! own real64 arithmetic rounds at that level every step, so a residual of a few
! floors that does not grow faster than sqrt(steps) is closure in real64.
! The terms are taken where the model applies them:
!   recom           change of the water pools in one call of REcoM, real64 at the
!                   interface (sources minus sinks and REcoM's lower bounds; they
!                   cannot be told apart without instrumenting REcoM itself)
!   benthos_recom   change of the benthic pools in that call (decay)
!   surface_bc      surface term as applied in diff_ver_part_impl_ale
!   benthic_return  dtr_bf of diff_ver_recom_expl
!   sinking         vert_sink of ver_sinking_recom
!   seafloor        str_bf of ver_sinking_recom_benthos
!   deposit         change of the benthic pools in the ocean step
! Not in the sum: the components of surface_bc from REcoM's real64 fluxes (air-sea,
! virtual, river, erosion, restoring, loopback), surface_bc minus them, and the
! benthic mismatches benthic_return+benthos_recom and seafloor+deposit (carbon).
! Not covered: explicit vertical diffusion (i_vert_diff false) and more than one
! FESOM group (__usetp with num_fesom_groups > 1).
!===============================================================================
module fesom_conservation_module
    use o_PARAM,          only: WP, WP_full, vcpw, S_ref_anomaly, clim_relax, use_tracer_carry
    use o_ARRAYS,         only: heat_flux, water_flux, virtual_salt, relax_salt, &
                                Tclim, Sclim, relax2clim, is_nonlinfs
    use g_forcing_arrays, only: sw_3d, real_salt_flux
    use g_config,         only: dt, use_sw_pene, toy_ocean, conservation_freq, use_cavity_fw2press
    use MOD_MESH
    use MOD_PARTIT
    use MOD_TRACER
#if defined(__recom)
    use recom_declarations, only: tracer_ids, is_riverinput, is_erosioninput
    use recom_config,       only: use_MEDUSA, add_loopback, benthos_num
    use recom_glovar,       only: tracers_info_type, GloCO2flux_seaicemask, virtual_dic, virtual_alk, &
                                  relax_alk, RiverDIC2D, RiverAlk2D, RiverDOC2D, ErosionTOC2D, &
                                  lb_flux, Benthos
#endif

    implicit none
    private
    public :: conservation_before_step, conservation_after_step, conservation_probe
#if defined(__recom)
    public :: conservation_recom, conservation_bgc_column, conservation_bgc_bc
#endif

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

#if defined(__recom)
    ! DIC, alkalinity, total carbon
    integer, parameter :: NB = 3
    character(len=6), parameter :: bname(NB) = ['DIC   ', 'Alk   ', 'carbon']
    integer, parameter :: T_RECOM = 1, T_BENREC = 2, T_BC = 3, T_RETURN = 4, T_SINK = 5, &
                          T_FLOOR = 6, T_DEPOSIT = 7, NTB = 7
    character(len=16), parameter :: tname(NTB) = ['recom           ', 'benthos_recom   ', &
        'surface_bc      ', 'benthic_return  ', 'sinking         ', 'seafloor        ', 'deposit         ']
    integer, parameter :: C_AIRSEA = 1, C_VIRT = 2, C_RIVER = 3, C_EROS = 4, C_RESTORE = 5, &
                          C_LOOP = 6, NCB = 6
    character(len=16), parameter :: cname(NCB) = ['  airsea        ', '  virtual       ', &
        '  river         ', '  erosion       ', '  restoring     ', '  loopback      ']
    ! weight of each tracer in each budget, and the tracers with a weight
    real(kind=WP_full), allocatable, save :: wb(:,:)
    integer,            allocatable, save :: btr(:)
    logical,                         save :: benthic_c = .false., has_doc = .false.
    ! the state REcoM reads, for the change in one call; the benthic carbon then
    ! and at the start of the ocean step
    real(kind=WP_full), allocatable, save :: before(:,:,:), ben_r(:), ben_s(:)
    ! terms and surface components since the start, local
    real(kind=WP_full), save :: bterm(NTB, NB) = 0.0_WP_full
    real(kind=WP_full), save :: bcomp(NCB, NB) = 0.0_WP_full
    ! start inventories as double words (rank 0)
    real(kind=WP_full), save :: b0h(NB) = 0.0_WP_full, b0l(NB) = 0.0_WP_full
    ! probes of DIC and alkalinity as double words (rank 0)
    real(kind=WP_full), save :: pbh(NPROBE, 2) = 0.0_WP_full, pbl(NPROBE, 2) = 0.0_WP_full
#endif

contains
!
!
!===============================================================================
! The start integrals. With REcoM they are taken before its first call
! (conservation_recom), so that REcoM's change in the first step of a run is a
! term like in every other step: the start inventory of a run is then the end
! inventory of the run it restarts from.
subroutine conservation_start(partit, tracers, mesh)
    type(t_partit), intent(inout), target :: partit
    type(t_tracer), intent(in),    target :: tracers
    type(t_mesh),   intent(in),    target :: mesh
#if defined(__recom)
    real(kind=WP_full) :: g(3*NB)
#endif

    call integrals(tracers, partit, mesh, i0)
#if defined(__recom)
    call bgc_setup(tracers, partit, mesh)
    call bgc_inventory(tracers, partit, mesh, g)
    b0h = g(1:NB)
    b0l = g(NB+1:2*NB)
#endif
    started = .true.
end subroutine conservation_start
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
#if defined(__recom)
    real(kind=WP_full) :: c(NCB, NB), x
#endif

    if (.not. started) call conservation_start(partit, tracers, mesh)
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

#if defined(__recom)
    ! the components of bc_surface for DIC (1002) and alkalinity (1003) from
    ! REcoM's real64 fluxes; carbon has DIC's and the river DOC (1012)
    c = 0.0_WP_full
    x = 0.0_WP_full
    do n = 1, partit%myDim_nod2D
        nzmin = mesh%ulevels_nod2D(n)
        a = dtf*real(mesh%areasvol(nzmin, n), WP_full)
        c(C_AIRSEA, 1) = c(C_AIRSEA, 1) + GloCO2flux_seaicemask(n)*a
        c(C_VIRT, 1)   = c(C_VIRT, 1)   + virtual_dic(n)*a
        c(C_VIRT, 2)   = c(C_VIRT, 2)   + virtual_alk(n)*a
        c(C_RIVER, 1)  = c(C_RIVER, 1)  + RiverDIC2D(n)*is_riverinput*a
        c(C_RIVER, 2)  = c(C_RIVER, 2)  + RiverAlk2D(n)*is_riverinput*a
        x              = x              + RiverDOC2D(n)*is_riverinput*a
        c(C_EROS, 1)   = c(C_EROS, 1)   + ErosionTOC2D(n)*is_erosioninput*a
        c(C_RESTORE, 2) = c(C_RESTORE, 2) + relax_alk(n)*a
        if (use_MEDUSA .and. add_loopback) then
            c(C_LOOP, 1) = c(C_LOOP, 1) + (lb_flux(n, 2) + lb_flux(n, 5))*a
            c(C_LOOP, 2) = c(C_LOOP, 2) + (lb_flux(n, 3) + 2.0_WP_full*lb_flux(n, 5))*a
        end if
    end do
    c(:, 3) = c(:, 1)
    if (has_doc) c(C_RIVER, 3) = c(C_RIVER, 3) + x
    bcomp = bcomp + c
    if (benthic_c) then
        do n = 1, partit%myDim_nod2D
            ben_s(n) = benthic_carbon(n, mesh)
        end do
    end if
#endif
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
#if defined(__recom)
    if (benthic_c) then
        w = 0.0_WP_full
        do n = 1, partit%myDim_nod2D
            w = w + (benthic_carbon(n, mesh) - ben_s(n))
        end do
        bterm(T_DEPOSIT, 3) = bterm(T_DEPOSIT, 3) + w
    end if
#endif
    if (mod(istep, conservation_freq) /= 0) return

    call integrals(tracers, partit, mesh, inow)
    l(1:NQ) = flux
    l(NQ+1:2*NQ) = relax
    call MPI_Allreduce(l, g, 2*NQ, MPI_WP_FULL, MPI_SUM, partit%MPI_COMM_FESOM, ierr)
    if (partit%mype == 0) then
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
    end if
#if defined(__recom)
    call bgc_report(istep, tracers, partit, mesh)
    if (partit%mype == 0) then
        ! DIC and alkalinity [mmol]
        do k = 2, np
            write(*, '(" CONSPB ",I8,1X,A12,2(1X,ES23.15E3))') istep, pname(k), &
                (pbh(k, 1) - pbh(k-1, 1)) + (pbl(k, 1) - pbl(k-1, 1)), &
                (pbh(k, 2) - pbh(k-1, 2)) + (pbl(k, 2) - pbl(k-1, 2))
        end do
    end if
#endif
end subroutine conservation_after_step
!
!
!===============================================================================
! Global integral of tracer tr_num (1 or 2) at a named point of the tracer step.
! Phase 'a', before the division by the new thickness: sum areasvol*(hnode*T + del_ttf);
! phase 'b', after it: sum areasvol*hnode_new*T. Both include the carry (use_tracer_carry).
! With REcoM the same for DIC and alkalinity, as double words.
subroutine conservation_probe(name, phase, tr_num, tracers, partit, mesh)
    character(len=*), intent(in)            :: name
    character(len=1), intent(in)            :: phase
    integer,          intent(in)            :: tr_num
    type(t_tracer),   intent(in),    target :: tracers
    type(t_partit),   intent(inout), target :: partit
    type(t_mesh),     intent(in),    target :: mesh
    real(kind=WP_full) :: l, g, a, lo, gg(2)
    integer            :: n, nz, k, ierr, slot

    if (conservation_freq <= 0) return
#ifdef ENABLE_OPENACC
    return
#endif
    slot = 0
    if (tr_num <= 2) slot = tr_num
#if defined(__recom)
    if (tracers%data(tr_num)%ID == 1002) slot = 3
    if (tracers%data(tr_num)%ID == 1003) slot = 4
#endif
    if (slot == 0) return
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
#if defined(__recom)
    if (slot > 2) then
        ! DIC, alkalinity: cell contents with the carry, summed as a double word
        l = 0.0_WP_full
        lo = 0.0_WP_full
        do n = 1, partit%myDim_nod2D
            do nz = mesh%ulevels_nod2D(n), mesh%nlevels_nod2D(n) - 1
                if (phase == 'a') then
                    a = mesh%hnode(nz, n)*real(tracers%data(tr_num)%values(nz, n), WP_full) &
                      + tracers%work%del_ttf(nz, n)
                else
                    a = mesh%hnode_new(nz, n)*real(tracers%data(tr_num)%values(nz, n), WP_full)
                end if
                if (use_tracer_carry) a = a + real(tracers%work%carry(nz, n, tr_num), WP_full)
                call two_sum(l, lo, a*real(mesh%areasvol(nz, n), WP_full))
            end do
        end do
        call dd_gather([l], [lo], gg(1:1), gg(2:2), 1, partit)
        pbh(k, slot-2) = gg(1)
        pbl(k, slot-2) = gg(2)
        return
    end if
#endif
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
    if (use_tracer_carry) then
        do n = 1, partit%myDim_nod2D
            do nz = mesh%ulevels_nod2D(n), mesh%nlevels_nod2D(n) - 1
                l = l + real(mesh%areasvol(nz, n), WP_full)*real(tracers%work%carry(nz, n, tr_num), WP_full)
            end do
        end do
    end if
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
    if (use_tracer_carry) then
        do n = 1, partit%myDim_nod2D
            do nz = mesh%ulevels_nod2D(n), mesh%nlevels_nod2D(n) - 1
                w = real(mesh%areasvol(nz, n), WP_full)
                l(2) = l(2) + w*real(tracers%work%carry(nz, n, 1), WP_full)
                l(3) = l(3) + w*real(tracers%work%carry(nz, n, 2), WP_full)
            end do
        end do
    end if
    l(2) = l(2)*real(vcpw, WP_full)
    call MPI_Allreduce(l, g, NQ, MPI_WP_FULL, MPI_SUM, partit%MPI_COMM_FESOM, ierr)
end subroutine integrals
!
!
!===============================================================================
! s + e += x with the rounding error of s + x kept in e (Knuth's two-sum)
subroutine two_sum(s, e, x)
    real(kind=WP_full), intent(inout) :: s, e
    real(kind=WP_full), intent(in)    :: x
    real(kind=WP_full) :: t, z
    t = s + x
    z = t - s
    e = e + ((s - (t - z)) + (x - z))
    s = t
end subroutine two_sum
#if defined(__recom)
!
!
!===============================================================================
! Global double-word sums of m local double words (hi, lo), on rank 0: gathered
! and added in rank order.
subroutine dd_gather(hi, lo, gh, gl, m, partit)
    integer,            intent(in)            :: m
    real(kind=WP_full), intent(in)            :: hi(m), lo(m)
    real(kind=WP_full), intent(out)           :: gh(m), gl(m)
    type(t_partit),     intent(inout), target :: partit
    real(kind=WP_full), allocatable :: r(:,:)
    real(kind=WP_full) :: s(2*m)
    integer :: p, k, ierr

    s(1:m) = hi
    s(m+1:2*m) = lo
    if (partit%mype == 0) then
        allocate(r(2*m, partit%npes))
    else
        allocate(r(1, 1))
    end if
    call MPI_Gather(s, 2*m, MPI_WP_FULL, r, 2*m, MPI_WP_FULL, 0, partit%MPI_COMM_FESOM, ierr)
    gh = 0.0_WP_full
    gl = 0.0_WP_full
    if (partit%mype == 0) then
        do p = 1, partit%npes
            do k = 1, m
                call two_sum(gh(k), gl(k), r(k, p))
                gl(k) = gl(k) + r(m+k, p)
            end do
        end do
    end if
    deallocate(r)
end subroutine dd_gather
!
!
!===============================================================================
! Which tracers count in which budget: DIC, alkalinity, and for total carbon all
! carbon pools of REcoM (inactive species have the ID -1); benthic carbon and
! calcite are pools 2 and 4.
subroutine bgc_setup(tracers, partit, mesh)
    type(t_tracer), intent(in),    target :: tracers
    type(t_partit), intent(inout), target :: partit
    type(t_mesh),   intent(in),    target :: mesh
    integer :: tr, id, k, nsel
    integer :: cid(14)

    if (allocated(wb)) return
    cid = [tracer_ids%dissolved_inorganic_carbon, tracer_ids%dissolved_organic_carbon,     &
           tracer_ids%phytoplankton_carbon, tracer_ids%diatom_carbon,                       &
           tracer_ids%heterotroph_carbon, tracer_ids%detrital_carbon,                       &
           tracer_ids%phytoplankton_calcite, tracer_ids%detrital_calcite,                   &
           tracer_ids%macrozooplankton_carbon, tracer_ids%macrozooplankton_detrital_carbon, &
           tracer_ids%macrozooplankton_detrital_calcite, tracer_ids%microzooplankton_carbon, &
           tracer_ids%coccolithophore_carbon, tracer_ids%phaeocystis_carbon]
    allocate(wb(NB, tracers%num_tracers))
    wb = 0.0_WP_full
    nsel = 0
    do tr = 1, tracers%num_tracers
        id = tracers%data(tr)%ID
        if (id == tracer_ids%dissolved_inorganic_carbon) wb(1, tr) = 1.0_WP_full
        if (id == tracer_ids%alkalinity)                 wb(2, tr) = 1.0_WP_full
        if (any(cid == id) .and. id > 0)                 wb(3, tr) = 1.0_WP_full
        if (id == tracer_ids%dissolved_organic_carbon)   has_doc = .true.
        if (any(wb(:, tr) /= 0.0_WP_full)) nsel = nsel + 1
    end do
    allocate(btr(nsel))
    k = 0
    do tr = 1, tracers%num_tracers
        if (any(wb(:, tr) /= 0.0_WP_full)) then
            k = k + 1
            btr(k) = tr
        end if
    end do
    allocate(before(mesh%nl-1, partit%myDim_nod2D, nsel))
    before = 0.0_WP_full
    benthic_c = benthos_num >= 4
    allocate(ben_r(partit%myDim_nod2D), ben_s(partit%myDim_nod2D))
    ben_r = 0.0_WP_full
    ben_s = 0.0_WP_full
    if (partit%mype == 0) then
        write(*, '(" CONSB budget tracers (ID: DIC Alk carbon):")')
        do k = 1, nsel
            write(*, '(" CONSB   ",I6,3(1X,F3.0))') tracers%data(btr(k))%ID, wb(:, btr(k))
        end do
        if (benthic_c) write(*, '(" CONSB   benthos pools 2 (C) and 4 (CaCO3) in carbon")')
    end if
end subroutine bgc_setup
!
!
!===============================================================================
! benthic carbon and calcite of node n [mmol]: REcoM's pools are per area(ulevels)
function benthic_carbon(n, mesh) result(c)
    integer,      intent(in)         :: n
    type(t_mesh), intent(in), target :: mesh
    real(kind=WP_full) :: c
    c = (Benthos(n, 2) + Benthos(n, 4))*real(mesh%area(mesh%ulevels_nod2D(n), n), WP_full)
end function benthic_carbon
!
!
!===============================================================================
! Inventories of the budgets as global double words (rank 0): g(1:NB) high
! words, g(NB+1:2NB) low words, g(2NB+1:3NB) sum of the squared cell contents.
subroutine bgc_inventory(tracers, partit, mesh, g)
    type(t_tracer),     intent(in),    target :: tracers
    type(t_partit),     intent(inout), target :: partit
    type(t_mesh),       intent(in),    target :: mesh
    real(kind=WP_full), intent(out)           :: g(3*NB)
    real(kind=WP_full) :: hi(2*NB), lo(2*NB), gh(2*NB), gl(2*NB), c, w
    integer :: k, tr, n, nz, q

    hi = 0.0_WP_full
    lo = 0.0_WP_full
    do k = 1, size(btr)
        tr = btr(k)
        do n = 1, partit%myDim_nod2D
            do nz = mesh%ulevels_nod2D(n), mesh%nlevels_nod2D(n) - 1
                c = real(tracers%data(tr)%values(nz, n), WP_full)*mesh%hnode(nz, n)
                if (use_tracer_carry) c = c + real(tracers%work%carry(nz, n, tr), WP_full)
                c = c*real(mesh%areasvol(nz, n), WP_full)
                do q = 1, NB
                    w = wb(q, tr)
                    if (w == 0.0_WP_full) cycle
                    call two_sum(hi(q), lo(q), w*c)
                    hi(NB+q) = hi(NB+q) + (w*c)**2
                end do
            end do
        end do
    end do
    if (benthic_c) then
        do n = 1, partit%myDim_nod2D
            c = benthic_carbon(n, mesh)
            call two_sum(hi(3), lo(3), c)
            hi(NB+3) = hi(NB+3) + c**2
        end do
    end if
    call dd_gather(hi, lo, gh, gl, 2*NB, partit)
    g(1:NB)        = gh(1:NB)
    g(NB+1:2*NB)   = gl(1:NB)
    g(2*NB+1:3*NB) = gh(NB+1:2*NB) + gl(NB+1:2*NB)
end subroutine bgc_inventory
!
!
!===============================================================================
! Around the call of REcoM (fesom_module): 'before' keeps the state REcoM reads
! (info%data_pointers, real64 in either precision), 'after' adds REcoM's change
! of it, cell by cell, before recom_tracers_back rounds it into the WP tracers.
subroutine conservation_recom(stage, info, tracers, partit, mesh)
    character(len=*),        intent(in)            :: stage
    type(tracers_info_type), intent(in)            :: info
    type(t_tracer),          intent(in),    target :: tracers
    type(t_partit),          intent(inout), target :: partit
    type(t_mesh),            intent(in),    target :: mesh
    real(kind=WP_full), pointer :: p(:,:)
    real(kind=WP_full) :: s, l(NB)
    integer :: k, tr, n, nz

    if (.not. started) then
        if (stage /= 'before') return
        call conservation_start(partit, tracers, mesh)
    end if
    if (stage == 'before') then
        do k = 1, size(btr)
            p => info%data_pointers(btr(k))%tracer_data
            do n = 1, partit%myDim_nod2D
                do nz = mesh%ulevels_nod2D(n), mesh%nlevels_nod2D(n) - 1
                    before(nz, n, k) = p(nz, n)
                end do
            end do
        end do
        if (benthic_c) then
            do n = 1, partit%myDim_nod2D
                ben_r(n) = benthic_carbon(n, mesh)
            end do
        end if
        return
    end if
    l = 0.0_WP_full
    do k = 1, size(btr)
        tr = btr(k)
        p => info%data_pointers(tr)%tracer_data
        s = 0.0_WP_full
        do n = 1, partit%myDim_nod2D
            do nz = mesh%ulevels_nod2D(n), mesh%nlevels_nod2D(n) - 1
                s = s + (p(nz, n) - before(nz, n, k))*mesh%hnode(nz, n)*real(mesh%areasvol(nz, n), WP_full)
            end do
        end do
        l = l + wb(:, tr)*s
    end do
    bterm(T_RECOM, :) = bterm(T_RECOM, :) + l
    if (benthic_c) then
        s = 0.0_WP_full
        do n = 1, partit%myDim_nod2D
            s = s + (benthic_carbon(n, mesh) - ben_r(n))
        end do
        bterm(T_BENREC, 3) = bterm(T_BENREC, 3) + s
    end if
end subroutine conservation_recom
!
!
!===============================================================================
! REcoM's column terms of tracer tr_num as diff_tracers_ale applies them, as
! field*hnode content: 'return' (dtr_bf), 'sinking' (vert_sink), 'seafloor' (str_bf)
subroutine conservation_bgc_column(term, tr_num, field, partit, mesh)
    character(len=*),   intent(in)            :: term
    integer,            intent(in)            :: tr_num
    real(kind=WP_full), intent(in)            :: field(:,:)
    type(t_partit),     intent(inout), target :: partit
    type(t_mesh),       intent(in),    target :: mesh
    real(kind=WP_full) :: s
    integer :: n, nz, t

    if (.not. started) return
    if (all(wb(:, tr_num) == 0.0_WP_full)) return
    select case (term)
    case ('return')
        t = T_RETURN
    case ('sinking')
        t = T_SINK
    case default
        t = T_FLOOR
    end select
    s = 0.0_WP_full
    do n = 1, partit%myDim_nod2D
        do nz = mesh%ulevels_nod2D(n), mesh%nlevels_nod2D(n) - 1
            s = s + field(nz, n)*mesh%hnode(nz, n)*real(mesh%areasvol(nz, n), WP_full)
        end do
    end do
    bterm(t, :) = bterm(t, :) + wb(:, tr_num)*s
end subroutine conservation_bgc_column
!
!
!===============================================================================
! The surface term of tracer tr_num as diff_ver_part_impl_ale adds it, bc*areasvol(nzmin)
subroutine conservation_bgc_bc(tr_num, bc, partit, mesh)
    integer,        intent(in)            :: tr_num
    real(kind=WP),  intent(in)            :: bc(:)
    type(t_partit), intent(inout), target :: partit
    type(t_mesh),   intent(in),    target :: mesh
    real(kind=WP_full) :: s
    integer :: n

    if (.not. started) return
    if (all(wb(:, tr_num) == 0.0_WP_full)) return
    s = 0.0_WP_full
    do n = 1, partit%myDim_nod2D
        s = s + real(bc(n), WP_full)*real(mesh%areasvol(mesh%ulevels_nod2D(n), n), WP_full)
    end do
    bterm(T_BC, :) = bterm(T_BC, :) + wb(:, tr_num)*s
end subroutine conservation_bgc_bc
!
!
!===============================================================================
! CONSB and CONSBT lines (rank 0; collective)
subroutine bgc_report(istep, tracers, partit, mesh)
    integer,        intent(in)            :: istep
    type(t_tracer), intent(in),    target :: tracers
    type(t_partit), intent(inout), target :: partit
    type(t_mesh),   intent(in),    target :: mesh
    real(kind=WP_full) :: g(3*NB), lt(NTB*NB + NCB*NB), gt(NTB*NB + NCB*NB)
    real(kind=WP_full) :: term(NTB, NB), comp(NCB, NB), inv, change, fl, res, floor
    integer :: q, t, ierr
    character(len=*), parameter :: fb = '(" CONSB ",I8,1X,A6,7(1X,ES23.15E3))'
    character(len=*), parameter :: ft = '(" CONSBT ",I8,1X,A6,1X,A16,1X,ES23.15E3)'

    call bgc_inventory(tracers, partit, mesh, g)
    lt(1:NTB*NB) = reshape(bterm, [NTB*NB])
    lt(NTB*NB+1:) = reshape(bcomp, [NCB*NB])
    call MPI_Reduce(lt, gt, size(lt), MPI_WP_FULL, MPI_SUM, 0, partit%MPI_COMM_FESOM, ierr)
    if (partit%mype /= 0) return
    term = reshape(gt(1:NTB*NB), [NTB, NB])
    comp = reshape(gt(NTB*NB+1:), [NCB, NB])
    do q = 1, NB
        inv    = g(q) + g(NB+q)
        change = (g(q) - b0h(q)) + (g(NB+q) - b0l(q))
        fl     = sum(term(:, q))
        res    = change - fl
        floor  = spacing(1.0_WP_full)*sqrt(g(2*NB+q))
        write(*, fb) istep, bname(q), inv, change, fl, res, res/max(abs(inv), tiny(1.0_WP_full)), &
                     floor, spacing(inv)
        do t = 1, NTB
            if (term(t, q) /= 0.0_WP_full) write(*, ft) istep, bname(q), tname(t), term(t, q)
        end do
        do t = 1, NCB
            if (comp(t, q) /= 0.0_WP_full) write(*, ft) istep, bname(q), cname(t), comp(t, q)
        end do
        write(*, ft) istep, bname(q), '  bc-components ', term(T_BC, q) - sum(comp(:, q))
        if (q == 3 .and. benthic_c) then
            write(*, ft) istep, bname(q), 'return+decay    ', term(T_RETURN, q) + term(T_BENREC, q)
            write(*, ft) istep, bname(q), 'seafloor+deposit', term(T_FLOOR, q) + term(T_DEPOSIT, q)
        end if
    end do
end subroutine bgc_report
#endif

end module fesom_conservation_module
