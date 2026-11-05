!> FESOM_F adapter for RECOM4 (RECOM4 M8, #46, #47, #48): the only FESOM_F code
!! that talks to RECOM4. RECOM4 is on only when FESOM_F is built with USE_RECOM4
!! (which defines __RECOM4) and FESOM3_RECOM_NML is set at run time.
!! Otherwise every routine returns at once and recom4_adapter_ntr() = 0, so
!! FESOM_F runs exactly as without RECOM4.
!!
!! Run-time variables:
!!   FESOM3_RECOM_NML          path of namelist.recom; RECOM4 is on when it is set
!!   FESOM3_RECOM_IC_DIR       folder with REcoM's input files (FESOM2's ClimateDataPath):
!!                             initial conditions and the dust climatology
!!   FESOM3_RECOM_DUMP_DIR     column dumps for gate item 1 (#48) into this folder
!!   FESOM3_RECOM_DUMP_POINTS  name:lon:lat[,name:lon:lat ...], degrees; the nearest
!!                             open-ocean node to each point is dumped
!!   FESOM3_RECOM_DUMP_STEPS   number of steps dumped (default 11)
!!   FESOM3_RECOM_CHECK        gate items 3 and 4 (#48): total N and the DICremin reset
module mod_recom4_adapter
    use mod_precision,     only: WP
    use mod_mesh,          only: t_mesh
    use mod_partit,        only: t_partit
    use mod_tracer,        only: t_tracer
    use mod_partitioning,  only: par_ex
    use oce_initial_state, only: t_ic3d_config
    use mod_forcing_read,  only: t_atm_forcing
#ifdef __RECOM4
    use recom4_api,     only: recom4_state, recom4_init, recom4_ntr, recom4_ncol,    &
                              recom4_ndiag2d, recom4_ndiag3d, recom4_ndiag3dw, recom4_tracer_id, &
                              recom4_initial_values, recom4_step, recom4_transport,  &
                              recom4_host_settings, recom4_get_host_settings,        &
                              recom4_col_profile, recom4_col_benthos, recom4_col_decay, &
                              RECOM4_D3_PAR, RECOM4_D3_GPP_PHY, RECOM4_D3_NPP_PHY,   &
                              RECOM4_D3_GPP_DIA, RECOM4_D3_NPP_DIA,                  &
                              RECOM4_D3_VERT_SINK, RECOM4_D3_DTR_BF, RECOM4_D3_STR_BF, &
                              RECOM4_D2_CO2FLUX, RECOM4_D2_CO2FLUX_DAY, RECOM4_D2_O2FLUX, &
                              RECOM4_D2_O2FLUX_DAY, RECOM4_D2_PCO2SURF, RECOM4_D2_DPCO2SURF, &
                              RECOM4_D2_PH_SURF, RECOM4_D2_K0, RECOM4_D2_KW660,       &
                              RECOM4_D2_DECAY_BEN_N, RECOM4_D3W_WSINK_DET1, RECOM4_D3W_WSINK_DET2
    use recom4_text_io, only: recom4_read_text
    use recom4_dump_io, only: dump_open_write, dump_write, dump_close
    use iso_fortran_env, only: real64
    use netcdf
    use mod_halo,       only: exchange_nod
    use mod_clock,      only: daynew, ndpyr
    use mpi
#endif
    implicit none
    private
    public :: recom4_adapter_init, recom4_adapter_on, recom4_adapter_ntr
    public :: recom4_adapter_register, recom4_adapter_ic, recom4_adapter_initial_values
    public :: recom4_adapter_forcing, recom4_adapter_monthly
    public :: recom4_adapter_step, recom4_adapter_bc_surface, recom4_adapter_transport

    logical :: on   = .false.
    integer :: ntr  = 0
    integer :: islp = 0      ! field number of sea-level pressure in FESOM_F's forcing
#ifdef __RECOM4
    type(recom4_state) :: state
    ! Per-node arrays of the host (IMPLEMENTATION.md section 1), all local nodes
    real(WP), allocatable :: col(:, :)          ! (ncol, nNodL)  per-column memory
    real(WP), allocatable :: sflux(:, :)        ! (ntr, nNodL)   surface flux of the step
    real(WP), allocatable :: diag2d(:, :)       ! (ndiag2d, nNodL)
    real(WP), allocatable :: diag3d(:, :, :)    ! (nl-1, ndiag3d, nNodL)
    real(WP), allocatable :: diag3dw(:, :, :)   ! (nl, ndiag3dw, nNodL), absolute levels:
                                                ! wsink_det1/2 from their transport calls (M9)
    real(WP), allocatable :: dust_fe(:)         ! (nNodL)  dust iron [umol/m2/s], monthly
    type(recom4_host_settings) :: hs
    integer  :: ialk = 0                        ! RECOM4 index of alkalinity (ID 1003)
    integer  :: itr_of_id(1000:1099) = 0        ! RECOM4 index of each REcoM ID
    real(WP) :: ocean_area = 0.0_WP             ! open-ocean area, cavities excluded [m2]
    real(WP), allocatable :: alk_surf(:)        ! (nNodL) initial surface alkalinity
    real(WP), allocatable :: relax_alk(:)       ! (nNodL) alkalinity restoring of the step

    ! Column dumps for gate item 1 (#48): off unless FESOM3_RECOM_DUMP_DIR is set
    integer, parameter :: NDUMP_MAX = 4
    integer :: ndump = 0, dump_nsteps = 0, cur_step = 0
    integer :: dump_node(NDUMP_MAX) = 0
    character(len=32) :: dump_name(NDUMP_MAX) = ''
    integer :: u_step(NDUMP_MAX), u_surf(NDUMP_MAX), u_surfo(NDUMP_MAX), u_tran(NDUMP_MAX)
    ! the carbonate profiles of col, in RECOM4's order (recom4_col_profile)
    character(len=8), parameter :: carb_names(8) = [character(len=8) ::                &
        'co2', 'ph', 'pco2', 'hco3', 'co3', 'omegac', 'kspc', 'rhosw']

    ! Gate items 3 and 4 (#48): off unless FESOM3_RECOM_CHECK is set
    logical  :: check_on = .false.
    real(WP) :: total_n0 = 0.0_WP              ! total N at the first check [mmol]
    integer, parameter :: NIDS_N = 9          ! the N tracers (2p3z2d; missing ones skipped)
    integer, parameter :: ids_n(NIDS_N) = [1001, 1004, 1007, 1009, 1011, 1013, 1023, 1025, 1029]
#endif

contains

    logical function recom4_adapter_on()
        recom4_adapter_on = on
    end function recom4_adapter_on

    integer function recom4_adapter_ntr()
        recom4_adapter_ntr = ntr
    end function recom4_adapter_ntr

    ! Before the tracers are allocated. nl: FESOM_F's number of levels;
    ! nNodL: the local nodes, halo included.
    subroutine recom4_adapter_init(nl, nNodL, partit)
        integer,        intent(in) :: nl, nNodL
        type(t_partit), intent(in) :: partit
#ifdef __RECOM4
        character(len=1024) :: path
        character(len=:), allocatable :: text(:)
        character(len=512) :: errmsg
        integer :: len_env, ios, ierr, i

        call get_environment_variable('FESOM3_RECOM_NML', path, length=len_env, status=ios)
        if (ios /= 0 .or. len_env == 0) return

        call recom4_read_text(trim(path), text, ierr, errmsg)
        if (ierr /= 0) call fail('reading ' // trim(path) // ': ' // trim(errmsg), partit)
        call recom4_init(text, nl - 1, state, ierr, errmsg)
        if (ierr /= 0) call fail('recom4_init: ' // trim(errmsg), partit)

        call recom4_get_host_settings(state, hs, ierr, errmsg)
        if (ierr /= 0) call fail('recom4_get_host_settings: ' // trim(errmsg), partit)
        ! Switches the adapter does not implement yet: stop rather than ignore them
        if (hs%use_atbox)          call fail('use_atbox is not supported by the adapter yet', partit)
        if (.not. hs%constant_CO2) call fail('constant_CO2 = .false. is not supported by the adapter yet', partit)
        if (hs%useAeolianN)        call fail('useAeolianN is not supported by the adapter yet', partit)
        if (hs%useRivers)          call fail('useRivers is not supported by the adapter yet', partit)
        if (hs%useErosion)         call fail('useErosion is not supported by the adapter yet', partit)
        if (hs%use_virt_tracers)   call fail('use_virt_tracers is not supported by the adapter yet', partit)

        on  = .true.
        ntr = recom4_ntr(state)
        do i = 1, ntr
            itr_of_id(recom4_tracer_id(state, i)) = i
            if (recom4_tracer_id(state, i) == 1003) ialk = i
        end do
        allocate(alk_surf(nNodL), relax_alk(nNodL))
        alk_surf  = 0.0_WP
        relax_alk = 0.0_WP
        allocate(col(recom4_ncol(state), nNodL), sflux(ntr, nNodL))
        allocate(diag2d(recom4_ndiag2d(state), nNodL))
        allocate(diag3d(nl - 1, recom4_ndiag3d(state), nNodL))
        allocate(diag3dw(nl, recom4_ndiag3dw(state), nNodL))
        diag3dw = 0.0_WP
        col    = 0.0_WP
        sflux  = 0.0_WP
        diag2d = 0.0_WP
        diag3d = 0.0_WP
        allocate(dust_fe(nNodL))
        dust_fe = 0.0_WP
        if (partit%mype == 0) write(*,'(a,i0,2a)') 'recom4_adapter: RECOM4 on, ', ntr, &
            ' tracers, ', trim(path)

        call get_environment_variable('FESOM3_RECOM_CHECK', path, length=len_env, status=ios)
        check_on = (ios == 0 .and. len_env > 0)
        if (check_on .and. partit%mype == 0) write(*,'(a)') &
            'recom4_adapter: checks of total N and the DICremin reset on'
#endif
    end subroutine recom4_adapter_init

    ! After tracers%data(:) is allocated with 2 + ntr entries: the IDs of data(3:).
    subroutine recom4_adapter_register(tracers)
        type(t_tracer), intent(inout) :: tracers
#ifdef __RECOM4
        integer :: i

        if (.not. on) return
        do i = 1, ntr
            tracers%data(2 + i)%ID = recom4_tracer_id(state, i)
        end do
#endif
    end subroutine recom4_adapter_register

    ! Before do_ic3d: adds the 7 REcoM fields of the reference &tracer_init3d
    ! (FESOM2-REcoM test3 namelist.tra), so that one do_ic3d call reads all.
    subroutine recom4_adapter_ic(ic, partit)
        type(t_ic3d_config), intent(inout) :: ic
        type(t_partit),      intent(in)    :: partit
#ifdef __RECOM4
        integer, parameter :: nic = 7
        integer, parameter :: ids(nic) = [1019, 1022, 1018, 1003, 1002, 1001, 1037]
        character(len=64), parameter :: files(nic) = [character(len=64) ::             &
            'fe_pisces_opa_eq_init_3D_changed_name.nc',                               &
            'woa18_all_o00_01_mmol_fesom2.nc',                                        &
            'woa13_all_i00_01_fesom2.nc',                                             &
            'GLODAPv2.2016b.TAlk_fesom2_mmol_fix_z_Fillvalue.nc',                     &
            'GLODAPv2.2016b.PI_TCO2_fesom2_mmol_fix_z_Fillvalue.nc',                  &
            'woa13_all_n00_01_fesom2.nc',                                             &
            'dicremin_fromAnnualAOU_WOA2023_noneg.nc']
        character(len=16), parameter :: vars(nic) = [character(len=16) ::              &
            'Fe', 'oxygen_mmol', 'i_an', 'TAlk_mmol', 'PI_TCO2_mmol', 'n_an', 'DICremin']
        character(len=1024) :: dir
        integer :: len_env, ios, i

        if (.not. on) return
        call get_environment_variable('FESOM3_RECOM_IC_DIR', dir, length=len_env, status=ios)
        if (ios /= 0 .or. len_env == 0) call fail('FESOM3_RECOM_IC_DIR is not set', partit)
        if (ic%n_ic3d + nic > size(ic%idlist)) call fail('too many initial-condition fields', partit)
        do i = 1, nic
            ic%n_ic3d = ic%n_ic3d + 1
            ic%idlist(ic%n_ic3d)   = ids(i)
            ic%filelist(ic%n_ic3d) = trim(dir) // '/' // trim(files(i))
            ic%varlist(ic%n_ic3d)  = vars(i)
        end do
#endif
    end subroutine recom4_adapter_ic

    ! After do_ic3d: every other tracer, as REcoM's recom_init. All local nodes,
    ! halo included (REcoM sets every node), so no halo exchange is needed.
    ! The vent test uses the static mid-depths Z: at a cold start the surface
    ! is flat, so Z equals REcoM's Z_3d_n.
    subroutine recom4_adapter_initial_values(tracers, mesh, nNodL, partit)
        type(t_tracer), intent(inout) :: tracers
        type(t_mesh),   intent(in)    :: mesh
        integer,        intent(in)    :: nNodL
        type(t_partit), intent(in)    :: partit
#ifdef __RECOM4
        real(WP), allocatable :: tr(:, :)
        real(WP), allocatable :: vmin(:), vmax(:)
        character(len=512) :: errmsg
        integer :: n, i, nzmin, nzmax, nz, ierr
        real(WP) :: lval

        if (.not. on) return
        allocate(tr(mesh%nl - 1, ntr))
        allocate(vmin(ntr), vmax(ntr))
        vmin = huge(1.0_WP)
        vmax = -huge(1.0_WP)
        do n = 1, nNodL
            nzmin = mesh%ulevels_nod2D(n)
            nzmax = mesh%nlevels_nod2D(n) - 1
            nz    = nzmax - nzmin + 1
            do i = 1, ntr
                tr(1:nz, i) = tracers%data(2 + i)%values(nzmin:nzmax, n)
            end do
            call recom4_initial_values(state, nz, mesh%Z(nzmin:nzmax),               &
                                       mesh%geo_coord_nod2D(1, n), mesh%geo_coord_nod2D(2, n), &
                                       tr(1:nz, :), ierr, errmsg)
            if (ierr /= 0) call fail('recom4_initial_values: ' // trim(errmsg), partit)
            do i = 1, ntr
                tracers%data(2 + i)%values(nzmin:nzmax, n) = tr(1:nz, i)
                vmin(i) = min(vmin(i), minval(tr(1:nz, i)))
                vmax(i) = max(vmax(i), maxval(tr(1:nz, i)))
            end do
        end do

        ! Alkalinity restoring (REcoM bio_fluxes): FESOM2 keeps the initial surface
        ! alkalinity, level 1, of every local node (oce_setup_step.F90), and the open-
        ! ocean area. At the first step Alk_surf - Alk = 0, so the restoring starts at 0.
        do n = 1, nNodL
            alk_surf(n) = tracers%data(2 + ialk)%values(1, n)
        end do
        lval = 0.0_WP
        do n = 1, partit%myDim_nod2D
            if (mesh%ulevels_nod2D(n) == 1) lval = lval + mesh%area(n)
        end do
        call MPI_Allreduce(lval, ocean_area, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
                           partit%MPI_COMM_FESOM, ierr)

        if (partit%mype == 0) then
            do i = 1, ntr
                write(*,'(a,i0,2es12.4)') 'recom4_adapter: initial (rank 0, wet levels) ID ', &
                    recom4_tracer_id(state, i), vmin(i), vmax(i)
            end do
        end if

        call setup_dumps(mesh, partit)
#endif
    end subroutine recom4_adapter_initial_values

    ! After FESOM_F has filled its forcing table: sea-level pressure as one more
    ! field, only with RECOM4 on. FESOM_F's bulk formulae do not use it; REcoM's
    ! gas exchange does. JRA55-do psl is in Pa, REcoM's unit.
    subroutine recom4_adapter_forcing(frc, forcing_dir, forc_set, partit)
        type(t_atm_forcing), intent(inout) :: frc
        character(len=*),    intent(in)    :: forcing_dir, forc_set
        type(t_partit),      intent(in)    :: partit
#ifdef __RECOM4
        if (.not. on) return
        if (trim(forc_set) /= 'JRA55') &
            call fail('RECOM4 needs FESOM3_FORCING=JRA55 (sea-level pressure psl)', partit)
        frc%nfld = frc%nfld + 1
        islp     = frc%nfld
        frc%f(islp)%file_base = trim(forcing_dir) // '/psl.'
        frc%f(islp)%varname   = 'psl'
#endif
    end subroutine recom4_adapter_forcing

    ! From the driver's roll_monthly_clim, with its month index: FESOM2's
    ! sbc_do_recom updates the dust at the same steps (the first step, then
    ! the last step of each month) with the same month (the current one at the
    ! first step, else the next). Albani monthly climatology, given on the
    ! CORE2 nodes [umol/m2/s], as fe_data_source = 'Albani' in the reference.
    subroutine recom4_adapter_monthly(imonth, nNodL, partit)
        integer,        intent(in) :: imonth, nNodL
        type(t_partit), intent(in) :: partit
#ifdef __RECOM4
        character(len=*), parameter :: fname = 'DustClimMonthlyAlbani.nc'
        character(len=*), parameter :: vname = 'DustClim'
        character(len=1024) :: dir
        real(WP), allocatable :: slice(:)
        integer :: ncid, varid, dimids(2), nnod_file, len_env, ios, n, gid

        if (.not. on) return
        call get_environment_variable('FESOM3_RECOM_IC_DIR', dir, length=len_env, status=ios)
        if (ios /= 0 .or. len_env == 0) call fail('FESOM3_RECOM_IC_DIR is not set', partit)

        call nc(nf90_open(trim(dir) // '/' // fname, nf90_nowrite, ncid), 'open ' // fname, partit)
        call nc(nf90_inq_varid(ncid, vname, varid), 'variable ' // vname, partit)
        call nc(nf90_inquire_variable(ncid, varid, dimids=dimids), 'dimensions of ' // vname, partit)
        call nc(nf90_inquire_dimension(ncid, dimids(1), len=nnod_file), 'nodes of ' // vname, partit)
        allocate(slice(nnod_file))
        call nc(nf90_get_var(ncid, varid, slice, start=[1, imonth], count=[nnod_file, 1]), &
                'read ' // vname, partit)
        call nc(nf90_close(ncid), 'close ' // fname, partit)

        do n = 1, nNodL
            gid = n
            if (partit%npes > 1) gid = partit%myList_nod2D(n)
            if (gid > nnod_file) call fail(fname // ' has fewer nodes than the mesh', partit)
            dust_fe(n) = slice(gid)
        end do
        if (partit%mype == 0) write(*,'(a,i0,2a)') 'recom4_adapter: dust iron of month ', &
            imonth, ' from ', fname
#endif
    end subroutine recom4_adapter_monthly

    ! Hook 1, after oce_fluxes and the shortwave penetration, before step_oce:
    ! recom4_step for every owned column, as REcoM's recom() loop
    ! (recom_main.F90 at 6e0b28aa), then the halo of the REcoM tracers.
    subroutine recom4_adapter_step(mstep, dt, tracers, mesh, partit, ice_frac, shortwave, &
                                   u_wind, v_wind, frc)
        integer,             intent(in)    :: mstep
        real(WP),            intent(in)    :: dt
        type(t_tracer),      intent(inout) :: tracers
        type(t_mesh),        intent(in)    :: mesh
        type(t_partit),      intent(in)    :: partit
        real(WP),            intent(in)    :: ice_frac(:), shortwave(:), u_wind(:), v_wind(:)
        type(t_atm_forcing), intent(in)    :: frc
#ifdef __RECOM4
        real(WP), allocatable :: tr(:, :), zero(:), rflux(:)
        character(len=512) :: errmsg
        real(WP) :: wind10
        integer  :: n, i, nzmin, nzmax, nz, ierr, kd

        if (.not. on) return
        cur_step = mstep
        if (check_on) call check_total_n(mstep, tracers, mesh, partit)
        if (hs%restore_alkalinity) call alkalinity_restoring(tracers, mesh, partit)

        allocate(tr(mesh%nl - 1, ntr), zero(ntr), rflux(ntr))
        zero = 0.0_WP
        do n = 1, partit%myDim_nod2D
            nzmin = mesh%ulevels_nod2D(n)
            nzmax = mesh%nlevels_nod2D(n) - 1
            nz    = nzmax - nzmin + 1
            do i = 1, ntr
                tr(1:nz, i) = tracers%data(2 + i)%values(nzmin:nzmax, n)
            end do
            rflux = 0.0_WP
            if (hs%restore_alkalinity) rflux(ialk) = relax_alk(n)
            wind10 = sqrt(u_wind(n)**2 + v_wind(n)**2)

            kd = dump_index(n)
            if (kd > 0) call dump_step_in(kd, n, nzmin, nzmax, nz, mstep, dt, mesh, tracers, &
                                          shortwave(n), ice_frac(n), frc%atmdata(islp, n), &
                                          wind10, tr(1:nz, :))

            call recom4_step(state, nz, mesh%hnode(nzmin:nzmax, n),                       &
                             mesh%zbar_3d_n(nzmin:nzmax + 1, n), mesh%Z_3d_n(nzmin:nzmax, n), &
                             tracers%data(1)%values(nzmin:nzmax, n),                    &
                             tracers%data(2)%values(nzmin:nzmax, n),                    &
                             mesh%geo_coord_nod2D(1, n), mesh%geo_coord_nod2D(2, n),    &
                             daynew, ndpyr, mstep, dt,                                  &
                             shortwave(n), ice_frac(n), frc%atmdata(islp, n), wind10,   &
                             hs%CO2_for_spinup, dust_fe(n), 0.0_WP, nzmin > 1,          &
                             zero, rflux, zero, zero, 0.0_WP,                           &
                             tr(1:nz, :), col(:, n), sflux(:, n), diag2d(:, n),         &
                             diag3d(1:nz, :, n), ierr, errmsg)
            if (ierr /= 0) call fail('recom4_step: ' // trim(errmsg), partit)

            if (kd > 0) call dump_step_out(kd, n, nz, mstep, tr(1:nz, :), rflux)

            do i = 1, ntr
                tracers%data(2 + i)%values(nzmin:nzmax, n) = tr(1:nz, i)
            end do
        end do

        if (check_on) call check_total_n(-mstep, tracers, mesh, partit)
        if (check_on) call check_dicremin_top(mstep, tracers, mesh, partit)

        ! REcoM's loop covers owned nodes only and "needs exchange_nod in the end"
        if (partit%npes > 1) then
            do i = 1, ntr
                call exchange_nod(tracers%data(2 + i)%values, partit)
            end do
        end if
#endif
    end subroutine recom4_adapter_step

    ! Hook 3: what the host adds to the top wet level of REcoM tracer id at
    ! owned node n in one step, dt * surface_flux (IMPLEMENTATION.md section 2).
    ! 0 when RECOM4 is off.
    real(WP) function recom4_adapter_bc_surface(id, n, dt)
        integer,  intent(in) :: id, n
        real(WP), intent(in) :: dt

        recom4_adapter_bc_surface = 0.0_WP
#ifdef __RECOM4
        if (.not. on) return
        if (id < lbound(itr_of_id, 1) .or. id > ubound(itr_of_id, 1)) return
        if (itr_of_id(id) == 0) return
        recom4_adapter_bc_surface = dt * sflux(itr_of_id(id), n)
#endif
    end function recom4_adapter_bc_surface

    ! Hook 2, in diff_tracers_ale before the ALE reconstruct, for tracer tr_num:
    ! recom4_transport for every owned column (benthic return, sinking,
    ! deposition), as REcoM's block in FESOM2's diff_tracers_ale. tr holds the
    ! current values of all REcoM tracers (those earlier in this tracer loop at
    ! their end-of-step state); only tr(:, itr) changes. FESOM_F's tracer cells
    ! are straight prisms: area and areasvol are the same at every level, the
    ! interface below the deepest layer has area 0 (as in FESOM2), and the whole
    ! seafloor lies under the bottom layer. A no-op for T and S.
    subroutine recom4_adapter_transport(tr_num, dt, nNodO, tracers, mesh, partit)
        integer,        intent(in)           :: tr_num, nNodO
        real(WP),       intent(in)           :: dt
        type(t_tracer), intent(inout)        :: tracers
        type(t_mesh),   intent(in)           :: mesh
        type(t_partit), intent(in), optional :: partit
#ifdef __RECOM4
        real(WP), allocatable :: tr(:, :), area_i(:), avol(:), seafl(:), d3t(:, :), d3wt(:, :)
        character(len=512) :: errmsg
        integer :: n, i, itr, nzmin, nzmax, nz, ierr, kd

        if (.not. on) return
        if (tr_num <= 2) return
        itr = tr_num - 2
        allocate(tr(mesh%nl - 1, ntr), area_i(mesh%nl), avol(mesh%nl - 1), seafl(mesh%nl - 1))
        ! recom4_transport sets its whole diag3d / diag3dw to 0 at every call;
        ! scratch arrays keep the step's diagnostics in diag3d (M9)
        allocate(d3t(mesh%nl - 1, recom4_ndiag3d(state)), d3wt(mesh%nl, recom4_ndiag3dw(state)))
        do n = 1, nNodO
            nzmin = mesh%ulevels_nod2D(n)
            nzmax = mesh%nlevels_nod2D(n) - 1
            nz    = nzmax - nzmin + 1
            do i = 1, ntr
                tr(1:nz, i) = tracers%data(2 + i)%values(nzmin:nzmax, n)
            end do
            area_i(1:nz + 1) = mesh%area(n)
            area_i(nz + 1)   = 0.0_WP      ! no water below the deepest layer (FESOM2: in_area ends in 0)
            avol(1:nz)       = mesh%areasvol(n)
            seafl(1:nz)      = 0.0_WP
            seafl(nz)        = mesh%area(n)

            kd = dump_index(n)
            if (kd > 0) call dump_transport_in(kd, n, nzmin, nzmax, nz, itr, dt, mesh, tracers, &
                                               area_i(1:nz + 1), avol(1:nz), seafl(1:nz), tr(1:nz, :))

            call recom4_transport(state, itr, nz, mesh%hnode(nzmin:nzmax, n),               &
                                  mesh%hnode_new(nzmin:nzmax, n),                           &
                                  mesh%zbar_3d_n(nzmin:nzmax + 1, n), mesh%Z_3d_n(nzmin:nzmax, n), &
                                  area_i(1:nz + 1), avol(1:nz), seafl(1:nz),                &
                                  tracers%data(1)%values(nzmin:nzmax, n),                   &
                                  tracers%data(2)%values(nzmin:nzmax, n),                   &
                                  mesh%geo_coord_nod2D(1, n), mesh%geo_coord_nod2D(2, n),   &
                                  dt, tr(1:nz, :), col(:, n), d3t(1:nz, :), d3wt(1:nz + 1, :), ierr, errmsg)
            if (ierr /= 0) then
                if (present(partit)) call fail('recom4_transport: ' // trim(errmsg), partit)
                error stop 'recom4_adapter: recom4_transport failed'
            end if

            ! the transport's own columns (6-8, column dumps), and the sinking
            ! speeds of the two detritus classes, only from their own calls
            diag3d(1:nz, RECOM4_D3_VERT_SINK:RECOM4_D3_STR_BF, n) = &
                d3t(1:nz, RECOM4_D3_VERT_SINK:RECOM4_D3_STR_BF)
            if (itr == itr_of_id(1021)) diag3dw(nzmin:nzmax + 1, RECOM4_D3W_WSINK_DET1, n) = &
                d3wt(1:nz + 1, RECOM4_D3W_WSINK_DET1)
            if (itr == itr_of_id(1028)) diag3dw(nzmin:nzmax + 1, RECOM4_D3W_WSINK_DET2, n) = &
                d3wt(1:nz + 1, RECOM4_D3W_WSINK_DET2)

            if (kd > 0) call dump_transport_out(kd, n, nz, itr, tr(1:nz, :))

            tracers%data(tr_num)%values(nzmin:nzmax, n) = tr(1:nz, itr)
        end do
#endif
    end subroutine recom4_adapter_transport

#ifdef __RECOM4
    ! REcoM's bio_fluxes (recom_main.F90 at 6e0b28aa): relaxation of the surface
    ! alkalinity (level 1) toward its initial value at open-ocean nodes, minus its
    ! area-weighted global mean. FESOM_F's area is the same at every level, so
    ! REcoM's areasvol(ulevels(n), n) is area(n) here.
    subroutine alkalinity_restoring(tracers, mesh, partit)
        type(t_tracer), intent(in) :: tracers
        type(t_mesh),   intent(in) :: mesh
        type(t_partit), intent(in) :: partit
        real(WP) :: lval, net
        integer  :: n, ierr

        do n = 1, size(relax_alk)
            if (mesh%ulevels_nod2D(n) == 1) then
                relax_alk(n) = hs%surf_relax_Alk * (alk_surf(n) - tracers%data(2 + ialk)%values(1, n))
            end if
        end do
        lval = 0.0_WP
        do n = 1, partit%myDim_nod2D
            lval = lval + relax_alk(n) * mesh%area(n)
        end do
        net = 0.0_WP
        call MPI_Allreduce(lval, net, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
                           partit%MPI_COMM_FESOM, ierr)
        do n = 1, size(relax_alk)
            if (mesh%ulevels_nod2D(n) == 1) relax_alk(n) = relax_alk(n) - net / ocean_area
        end do
    end subroutine alkalinity_restoring

    ! Gate item 3: the ocean's total nitrogen at the start of hook 1 [mmol], the N
    ! tracers of every owned wet cell (concentration * hnode * area) plus the
    ! benthic N pool (mmol/m2 * area), and its relative change since the first call.
    ! Printed on rank 0 in FESOM_F's format for global sums (es).
    subroutine check_total_n(mstep, tracers, mesh, partit)
        integer,        intent(in) :: mstep
        type(t_tracer), intent(in) :: tracers
        type(t_mesh),   intent(in) :: mesh
        type(t_partit), intent(in) :: partit
        real(WP) :: lval(2), tot(2)
        integer  :: n, j, i, nzmin, nzmax, ierr

        lval = 0.0_WP
        do n = 1, partit%myDim_nod2D
            nzmin = mesh%ulevels_nod2D(n)
            nzmax = mesh%nlevels_nod2D(n) - 1
            do j = 1, NIDS_N
                i = itr_of_id(ids_n(j))
                if (i == 0) cycle
                lval(1) = lval(1) + sum(tracers%data(2 + i)%values(nzmin:nzmax, n) &
                                        * mesh%hnode(nzmin:nzmax, n)) * mesh%area(n)
            end do
            lval(2) = lval(2) + col(recom4_col_benthos(state, 1), n) * mesh%area(n)
        end do
        call MPI_Allreduce(lval, tot, 2, MPI_DOUBLE_PRECISION, MPI_SUM, partit%MPI_COMM_FESOM, ierr)
        if (total_n0 == 0.0_WP) total_n0 = tot(1) + tot(2)
        if (partit%mype == 0) write(*,'(a,i0,a,es23.16,a,es23.16,a,es23.16,a,es10.3)')       &
            'recom4_adapter: step ', mstep, '  water N ', tot(1), '  benthic N ', tot(2),    &
            '  total N ', tot(1) + tot(2), '  rel. change ', (tot(1) + tot(2) - total_n0) / total_n0
    end subroutine check_total_n
    ! Gate item 4: after hook 1, the largest DICremin at the top wet level, which
    ! REcoM resets to 0 before its sources minus sinks, next to the largest at any level.
    subroutine check_dicremin_top(mstep, tracers, mesh, partit)
        integer,        intent(in) :: mstep
        type(t_tracer), intent(in) :: tracers
        type(t_mesh),   intent(in) :: mesh
        type(t_partit), intent(in) :: partit
        real(WP) :: lmax(2), gmax(2)
        integer  :: n, i, nzmin, nzmax, ierr

        i = itr_of_id(1037)
        if (i == 0) return
        lmax = 0.0_WP
        do n = 1, partit%myDim_nod2D
            nzmin = mesh%ulevels_nod2D(n)
            nzmax = mesh%nlevels_nod2D(n) - 1
            lmax(1) = max(lmax(1), abs(tracers%data(2 + i)%values(nzmin, n)))
            lmax(2) = max(lmax(2), maxval(abs(tracers%data(2 + i)%values(nzmin:nzmax, n))))
        end do
        call MPI_Allreduce(lmax, gmax, 2, MPI_DOUBLE_PRECISION, MPI_MAX, partit%MPI_COMM_FESOM, ierr)
        if (partit%mype == 0) write(*,'(a,i0,a,es10.3,a,es10.3)') 'recom4_adapter: step ', &
            mstep, '  DICremin max, top wet level ', gmax(1), ', any level ', gmax(2)
    end subroutine check_dicremin_top

    ! ---------------------------------------------------------------- column dumps
    ! The nearest open-ocean node to each point of FESOM3_RECOM_DUMP_POINTS, over
    ! all ranks (as the FESOM2-REcoM reference dumps chose their columns, see
    ! tools/ref/SOURCE); the owning rank opens the four dump files of the column.
    subroutine setup_dumps(mesh, partit)
        type(t_mesh),   intent(in) :: mesh
        type(t_partit), intent(in) :: partit
        character(len=1024) :: dir, pts, env
        character(len=256)  :: item
        real(WP) :: plon, plat, d, deg, best(2), gbest(2)
        integer  :: len_env, ios, i0, i1, c1, c2, n, nbest, ierr

        call get_environment_variable('FESOM3_RECOM_DUMP_DIR', dir, length=len_env, status=ios)
        if (ios /= 0 .or. len_env == 0) return
        call get_environment_variable('FESOM3_RECOM_DUMP_POINTS', pts, length=len_env, status=ios)
        if (ios /= 0 .or. len_env == 0) &
            call fail('FESOM3_RECOM_DUMP_DIR is set, but FESOM3_RECOM_DUMP_POINTS is not', partit)
        dump_nsteps = 11
        call get_environment_variable('FESOM3_RECOM_DUMP_STEPS', env, length=len_env, status=ios)
        if (ios == 0 .and. len_env > 0) read(env, *) dump_nsteps

        deg = acos(-1.0_WP) / 180.0_WP
        i0 = 1
        do while (i0 <= len_trim(pts))
            i1 = index(pts(i0:), ',')
            if (i1 == 0) then
                item = pts(i0:len_trim(pts))
                i0   = len_trim(pts) + 1
            else
                item = pts(i0:i0 + i1 - 2)
                i0   = i0 + i1
            end if
            c1 = index(item, ':')
            c2 = index(item, ':', back=.true.)
            if (c1 == 0 .or. c2 == c1) &
                call fail('FESOM3_RECOM_DUMP_POINTS: expected name:lon:lat, got ' // trim(item), partit)
            read(item(c1 + 1:c2 - 1), *) plon
            read(item(c2 + 1:), *) plat

            best  = [huge(1.0_WP), real(partit%mype, WP)]
            nbest = 0
            do n = 1, partit%myDim_nod2D
                if (mesh%ulevels_nod2D(n) /= 1) cycle
                d = gc_angle(plon * deg, plat * deg, mesh%geo_coord_nod2D(1, n), mesh%geo_coord_nod2D(2, n))
                if (d < best(1)) then
                    best(1) = d
                    nbest   = n
                end if
            end do
            call MPI_Allreduce(best, gbest, 1, MPI_2DOUBLE_PRECISION, MPI_MINLOC, &
                               partit%MPI_COMM_FESOM, ierr)
            if (nint(gbest(2)) /= partit%mype) cycle

            if (ndump == NDUMP_MAX) call fail('more than 4 dump columns on one rank', partit)
            ndump = ndump + 1
            dump_node(ndump) = nbest
            dump_name(ndump) = item(1:c1 - 1)
            call open_dump(trim(dir) // '/' // trim(dump_name(ndump)) // '_step.bin',        u_step(ndump),  partit)
            call open_dump(trim(dir) // '/' // trim(dump_name(ndump)) // '_surface.bin',     u_surf(ndump),  partit)
            call open_dump(trim(dir) // '/' // trim(dump_name(ndump)) // '_surface_out.bin', u_surfo(ndump), partit)
            call open_dump(trim(dir) // '/' // trim(dump_name(ndump)) // '_transport.bin',   u_tran(ndump),  partit)
            write(*,'(3a,i0,a,i0,a,f10.4,a,f10.4,a,i0,a)') 'recom4_adapter: dump column ',     &
                trim(dump_name(ndump)), ': rank ', partit%mype, ', local node ', nbest,        &
                ', lon ', mesh%geo_coord_nod2D(1, nbest) / deg, ', lat ',                      &
                mesh%geo_coord_nod2D(2, nbest) / deg, ', ',                                    &
                mesh%nlevels_nod2D(nbest) - mesh%ulevels_nod2D(nbest), ' levels'
        end do
    end subroutine setup_dumps

    subroutine open_dump(path, u, partit)
        character(len=*), intent(in)  :: path
        integer,          intent(out) :: u
        type(t_partit),   intent(in)  :: partit
        character(len=512) :: errmsg
        integer :: ierr

        call dump_open_write(path, u, ierr, errmsg)
        if (ierr /= 0) call fail('opening ' // path // ': ' // trim(errmsg), partit)
    end subroutine open_dump

    ! Index of node n among this rank's dump columns, 0 if none or past the dump steps
    integer function dump_index(n)
        integer, intent(in) :: n
        integer :: k

        dump_index = 0
        if (cur_step < 1 .or. cur_step > dump_nsteps) return
        do k = 1, ndump
            if (dump_node(k) == n) dump_index = k
        end do
    end function dump_index

    ! Hook 1, before recom4_step: what the column driver reads (drivers/recom4_column.F90)
    subroutine dump_step_in(k, n, nzmin, nzmax, nz, m, dt, mesh, tracers, sw, ice, slp, wind10, tr)
        integer,        intent(in) :: k, n, nzmin, nzmax, nz, m
        real(WP),       intent(in) :: dt, sw, ice, slp, wind10, tr(:, :)
        type(t_mesh),   intent(in) :: mesh
        type(t_tracer), intent(in) :: tracers
        integer :: u, j, i

        u = u_step(k)
        call w0(u, 'in_nzmin',  m, real(nzmin, WP))
        call w0(u, 'in_nzmax',  m, real(nzmax, WP))
        call w0(u, 'in_mstep',  m, real(m, WP))
        call w0(u, 'in_daynew', m, real(daynew, WP))
        call w0(u, 'in_ndpyr',  m, real(ndpyr, WP))
        call w0(u, 'in_dt',     m, dt)
        call w0(u, 'in_lon',    m, mesh%geo_coord_nod2D(1, n))
        call w0(u, 'in_lat',    m, mesh%geo_coord_nod2D(2, n))
        call w1(u, 'in_hnode',  m, mesh%hnode(nzmin:nzmax, n))
        call w1(u, 'in_zbar',   m, mesh%zbar_3d_n(nzmin:nzmax + 1, n))
        call w1(u, 'in_zmid',   m, mesh%Z_3d_n(nzmin:nzmax, n))
        call w1(u, 'in_temp',   m, tracers%data(1)%values(nzmin:nzmax, n))
        call w1(u, 'in_salt',   m, tracers%data(2)%values(nzmin:nzmax, n))
        call w0(u, 'in_shortwave_raw', m, sw)
        call w0(u, 'in_ice_frac',      m, ice)
        call w0(u, 'in_press_air',     m, slp)
        call w0(u, 'in_wind10',        m, wind10)
        call w0(u, 'in_atm_co2',       m, hs%CO2_for_spinup)
        call w0(u, 'in_dust_fe_raw',   m, dust_fe(n))
        call w0(u, 'in_dust_n_raw',    m, 0.0_WP)
        do i = 1, size(carb_names)
            j = recom4_col_profile(state, i)
            call w1(u, 'in_' // trim(carb_names(i)), m, col(j:j + nz - 1, n))
        end do
        call w1(u, 'in_benthos', m, col(recom4_col_benthos(state, 1):recom4_col_benthos(state, 4), n))
        call w2(u, 'in_tr', m, tr)
    end subroutine dump_step_in

    ! Hook 1, after recom4_step: its results, and the surface terms
    subroutine dump_step_out(k, n, nz, m, tr, rflux)
        integer,  intent(in) :: k, n, nz, m
        real(WP), intent(in) :: tr(:, :), rflux(:)
        character(len=16) :: id
        integer :: u, i, j

        u = u_step(k)
        call w2(u, 'out_tr',      m, tr)
        call w1(u, 'out_par',     m, diag3d(1:nz, RECOM4_D3_PAR, n))
        call w1(u, 'out_npp_phy', m, diag3d(1:nz, RECOM4_D3_NPP_PHY, n))
        call w1(u, 'out_gpp_phy', m, diag3d(1:nz, RECOM4_D3_GPP_PHY, n))
        call w1(u, 'out_npp_dia', m, diag3d(1:nz, RECOM4_D3_NPP_DIA, n))
        call w1(u, 'out_gpp_dia', m, diag3d(1:nz, RECOM4_D3_GPP_DIA, n))
        ! the air-sea gas exchange, one value per column (as drivers/recom4_column.F90)
        call w0(u, 'out_co2flux',     m, diag2d(RECOM4_D2_CO2FLUX, n))
        call w0(u, 'out_co2flux_day', m, diag2d(RECOM4_D2_CO2FLUX_DAY, n))
        call w0(u, 'out_o2flux',      m, diag2d(RECOM4_D2_O2FLUX, n))
        call w0(u, 'out_o2flux_day',  m, diag2d(RECOM4_D2_O2FLUX_DAY, n))
        call w0(u, 'out_pco2surf',    m, diag2d(RECOM4_D2_PCO2SURF, n))
        call w0(u, 'out_dpco2surf',   m, diag2d(RECOM4_D2_DPCO2SURF, n))
        call w0(u, 'out_ph_surf',     m, diag2d(RECOM4_D2_PH_SURF, n))
        call w0(u, 'out_k0',          m, diag2d(RECOM4_D2_K0, n))
        call w0(u, 'out_kw660',       m, diag2d(RECOM4_D2_KW660, n))
        ! the benthic pools after the step, and their decay per day from diag2d
        call w1(u, 'out_benthos',       m, col(recom4_col_benthos(state, 1):recom4_col_benthos(state, 4), n))
        call w1(u, 'out_decay_benthos', m, diag2d(RECOM4_D2_DECAY_BEN_N:RECOM4_D2_DECAY_BEN_N + 3, n))
        do i = 1, size(carb_names)
            j = recom4_col_profile(state, i)
            call w1(u, 'out_' // trim(carb_names(i)), m, col(j:j + nz - 1, n))
        end do
        ! surf_<ID>: the replay reads terms 2-5 (virtual, restoring, river, erosion);
        ! RECOM4 does not return its own term separately, so term 1 is 0 here and
        ! the total, for the comparison, goes to _surface_out.bin
        do i = 1, ntr
            write(id, '(i0)') recom4_tracer_id(state, i)
            call w1(u_surf(k), 'surf_' // trim(id), m, &
                    [0.0_WP, 0.0_WP, rflux(i), 0.0_WP, 0.0_WP, sflux(i, n)])
            call w0(u_surfo(k), 'out_surf_' // trim(id), m, sflux(i, n))
        end do
    end subroutine dump_step_out

    ! Hook 2, before recom4_transport; the column's inputs once per step, at the
    ! first REcoM tracer, then the tracer itself
    subroutine dump_transport_in(k, n, nzmin, nzmax, nz, itr, dt, mesh, tracers, area_i, avol, seafl, tr)
        integer,        intent(in) :: k, n, nzmin, nzmax, nz, itr
        real(WP),       intent(in) :: dt, area_i(:), avol(:), seafl(:), tr(:, :)
        type(t_mesh),   intent(in) :: mesh
        type(t_tracer), intent(in) :: tracers
        character(len=16) :: id
        integer :: u, m

        u = u_tran(k)
        m = cur_step
        if (itr == 1) then
            call w0(u, 'in_nzmin',     m, real(nzmin, WP))
            call w0(u, 'in_nzmax',     m, real(nzmax, WP))
            call w0(u, 'in_mstep',     m, real(m, WP))
            call w0(u, 'in_dt',        m, dt)
            call w0(u, 'in_lat',       m, mesh%geo_coord_nod2D(2, n))
            call w1(u, 'in_hnode',     m, mesh%hnode(nzmin:nzmax, n))
            call w1(u, 'in_hnode_new', m, mesh%hnode_new(nzmin:nzmax, n))
            call w1(u, 'in_zbar',      m, mesh%zbar_3d_n(nzmin:nzmax + 1, n))
            call w1(u, 'in_zmid',      m, mesh%Z_3d_n(nzmin:nzmax, n))
            call w1(u, 'in_area',      m, area_i)
            call w1(u, 'in_areasvol',  m, avol)
            call w1(u, 'in_seafloor',  m, seafl)
            call w1(u, 'in_temp',      m, tracers%data(1)%values(nzmin:nzmax, n))
            call w1(u, 'in_salt',      m, tracers%data(2)%values(nzmin:nzmax, n))
            call w1(u, 'in_benthos',       m, col(recom4_col_benthos(state, 1):recom4_col_benthos(state, 4), n))
            call w1(u, 'in_decay_benthos', m, col(recom4_col_decay(state, 1):recom4_col_decay(state, 4), n))
        end if
        write(id, '(i0)') recom4_tracer_id(state, itr)
        call w1(u, 'in_tr_' // trim(id), m, tr(1:nz, itr))
    end subroutine dump_transport_in

    ! Hook 2, after recom4_transport; closes the files after the last dumped step
    subroutine dump_transport_out(k, n, nz, itr, tr)
        integer,  intent(in) :: k, n, nz, itr
        real(WP), intent(in) :: tr(:, :)
        character(len=16) :: id
        integer :: u, m

        u = u_tran(k)
        m = cur_step
        write(id, '(i0)') recom4_tracer_id(state, itr)
        call w1(u, 'out_vert_sink_' // trim(id), m, diag3d(1:nz, RECOM4_D3_VERT_SINK, n))
        call w1(u, 'out_dtr_bf_'    // trim(id), m, diag3d(1:nz, RECOM4_D3_DTR_BF, n))
        call w1(u, 'out_str_bf_'    // trim(id), m, diag3d(1:nz, RECOM4_D3_STR_BF, n))
        call w1(u, 'out_tr_sink_'   // trim(id), m, tr(1:nz, itr))
        if (itr == ntr .and. m == dump_nsteps) then
            call dump_close(u_step(k))
            call dump_close(u_surf(k))
            call dump_close(u_surfo(k))
            call dump_close(u_tran(k))
        end if
    end subroutine dump_transport_out

    subroutine w0(u, name, m, x)
        integer,          intent(in) :: u, m
        character(len=*), intent(in) :: name
        real(WP),         intent(in) :: x
        character(len=256) :: errmsg
        integer :: ierr

        call dump_write(u, name, m, real(x, real64), ierr, errmsg)
        if (ierr /= 0) error stop 'recom4_adapter: writing a column dump failed'
    end subroutine w0

    subroutine w1(u, name, m, v)
        integer,          intent(in) :: u, m
        character(len=*), intent(in) :: name
        real(WP),         intent(in) :: v(:)
        character(len=256) :: errmsg
        integer :: ierr

        call dump_write(u, name, m, real(v, real64), ierr, errmsg)
        if (ierr /= 0) error stop 'recom4_adapter: writing a column dump failed'
    end subroutine w1

    subroutine w2(u, name, m, a)
        integer,          intent(in) :: u, m
        character(len=*), intent(in) :: name
        real(WP),         intent(in) :: a(:, :)
        character(len=256) :: errmsg
        integer :: ierr

        call dump_write(u, name, m, real(a, real64), ierr, errmsg)
        if (ierr /= 0) error stop 'recom4_adapter: writing a column dump failed'
    end subroutine w2

    ! Great-circle angle between two points (radians), haversine form
    pure real(WP) function gc_angle(lon1, lat1, lon2, lat2)
        real(WP), intent(in) :: lon1, lat1, lon2, lat2
        real(WP) :: s1, s2

        s1 = sin(0.5_WP * (lat2 - lat1))
        s2 = sin(0.5_WP * (lon2 - lon1))
        gc_angle = 2.0_WP * asin(min(1.0_WP, sqrt(s1 * s1 + cos(lat1) * cos(lat2) * s2 * s2)))
    end function gc_angle

    ! Stops with a message if a netCDF call failed.
    subroutine nc(status, what, partit)
        integer,          intent(in) :: status
        character(len=*), intent(in) :: what
        type(t_partit),   intent(in) :: partit

        if (status /= nf90_noerr) call fail('netCDF: ' // what // ': ' // trim(nf90_strerror(status)), partit)
    end subroutine nc
#endif

    ! Prints msg and aborts every rank (a failing rank must not leave the others
    ! waiting in a collective), as FESOM_F's own die().
    subroutine fail(msg, partit)
        character(len=*), intent(in) :: msg
        type(t_partit),   intent(in) :: partit

        write(*,'(a,i0,2a)') 'recom4_adapter (rank ', partit%mype, '): ', msg
        call par_ex(partit%MPI_COMM_FESOM, partit%mype, 1)
    end subroutine fail

end module mod_recom4_adapter
