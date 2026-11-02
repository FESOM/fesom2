!> FESOM_F adapter for RECOM4 (RECOM4 M8, #46): the only FESOM_F code that
!! talks to RECOM4. RECOM4 is on only when FESOM_F is built with USE_RECOM4
!! (which defines __RECOM4) and FESOM3_RECOM_NML is set at run time.
!! Otherwise every routine returns at once and recom4_adapter_ntr() = 0, so
!! FESOM_F runs exactly as without RECOM4.
!!
!! Run-time variables:
!!   FESOM3_RECOM_NML     path of namelist.recom; RECOM4 is on when it is set
!!   FESOM3_RECOM_IC_DIR  folder with the REcoM initial-condition files
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
                              recom4_ndiag2d, recom4_ndiag3d, recom4_tracer_id,      &
                              recom4_initial_values
    use recom4_text_io, only: recom4_read_text
#endif
    implicit none
    private
    public :: recom4_adapter_init, recom4_adapter_on, recom4_adapter_ntr
    public :: recom4_adapter_register, recom4_adapter_ic, recom4_adapter_initial_values
    public :: recom4_adapter_forcing

    logical :: on  = .false.
    integer :: ntr = 0
    integer :: islp = 0      ! field number of sea-level pressure in FESOM_F's forcing
#ifdef __RECOM4
    type(recom4_state) :: state
    ! Per-node arrays of the host (IMPLEMENTATION.md section 1), all local nodes
    real(WP), allocatable :: col(:, :)          ! (ncol, nNodL)  per-column memory
    real(WP), allocatable :: sflux(:, :)        ! (ntr, nNodL)   surface flux of the step
    real(WP), allocatable :: diag2d(:, :)       ! (ndiag2d, nNodL)
    real(WP), allocatable :: diag3d(:, :, :)    ! (nl-1, ndiag3d, nNodL)
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
        integer :: len_env, ios, ierr

        call get_environment_variable('FESOM3_RECOM_NML', path, length=len_env, status=ios)
        if (ios /= 0 .or. len_env == 0) return

        call recom4_read_text(trim(path), text, ierr, errmsg)
        if (ierr /= 0) call fail('reading ' // trim(path) // ': ' // trim(errmsg), partit)
        call recom4_init(text, nl - 1, state, ierr, errmsg)
        if (ierr /= 0) call fail('recom4_init: ' // trim(errmsg), partit)

        on  = .true.
        ntr = recom4_ntr(state)
        allocate(col(recom4_ncol(state), nNodL), sflux(ntr, nNodL))
        allocate(diag2d(recom4_ndiag2d(state), nNodL))
        allocate(diag3d(nl - 1, recom4_ndiag3d(state), nNodL))
        col    = 0.0_WP
        sflux  = 0.0_WP
        diag2d = 0.0_WP
        diag3d = 0.0_WP
        if (partit%mype == 0) write(*,'(a,i0,2a)') 'recom4_adapter: RECOM4 on, ', ntr, &
            ' tracers, ', trim(path)
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
        character(len=512) :: errmsg
        integer :: n, i, nzmin, nzmax, nz, ierr
        real(WP), allocatable :: vmin(:), vmax(:)

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
        if (partit%mype == 0) then
            do i = 1, ntr
                write(*,'(a,i0,2es12.4)') 'recom4_adapter: initial (rank 0, wet levels) ID ', &
                    recom4_tracer_id(state, i), vmin(i), vmax(i)
            end do
        end if
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

    ! Prints msg and aborts every rank (a failing rank must not leave the others
    ! waiting in a collective), as FESOM_F's own die().
    subroutine fail(msg, partit)
        character(len=*), intent(in) :: msg
        type(t_partit),   intent(in) :: partit

        write(*,'(a,i0,2a)') 'recom4_adapter (rank ', partit%mype, '): ', msg
        call par_ex(partit%MPI_COMM_FESOM, partit%mype, 1)
    end subroutine fail

end module mod_recom4_adapter
