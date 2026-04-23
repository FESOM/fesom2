!==============================================================================
! fesom_dump_shim.F90
!
! Reference-dump instrumentation for the FESOM2 -> C port.
!
! When the environment variable FESOM_DUMP_FILE is set to a path prefix, this
! module opens a per-rank binary file at <prefix>.<mype:05d> and writes a
! stream of records describing the value of named fields at a small set of
! probe nodes (hardcoded global IDs) after each substep of oce_timestep_ale.
!
! The C port reads the same file and asserts per-substep parity, turning
! later debugging sessions from days into hours. Per-rank files avoid an
! MPI gather; on serial Phase-1 runs, only one file is produced.
!
! Activation:
!   FESOM_DUMP_FILE       (required, no default)  -- enables shim
!   FESOM_DUMP_MAXSTEPS   (optional, default 10)  -- stop recording after step N
!
! Binary record layout (little-endian, stream access, no header):
!
!   int32   step                ( >= 1; 0 reserved for pre-loop init dumps )
!   int32   substep_id          ( DUMP_SUBSTEP_* below )
!   int32   probe_global_id     ( Fortran 1-based node ID )
!   int32   nlevels
!   char[24] field_name         ( blank-padded ASCII, e.g. "density" )
!   real64  values[nlevels]     ( column from surface to nlevels-th interface )
!
! Total per record: 40 + 8*nlevels bytes. The reader streams records to EOF
! and groups them by (step, substep_id, probe_global_id, field_name).
!==============================================================================
module fesom_dump_shim
  use, intrinsic :: iso_fortran_env, only: int32, real64
  implicit none
  private

  ! ---- Substep IDs (mirror FRESH_START.md section 5) -----------------------
  integer, parameter, public :: DUMP_SUBSTEP_INIT          = 0
  integer, parameter, public :: DUMP_SUBSTEP_PRESSURE_BV   = 1
  integer, parameter, public :: DUMP_SUBSTEP_SW_AB         = 2
  integer, parameter, public :: DUMP_SUBSTEP_PGF           = 3
  integer, parameter, public :: DUMP_SUBSTEP_MIXING        = 4
  integer, parameter, public :: DUMP_SUBSTEP_VEL_RHS       = 5
  integer, parameter, public :: DUMP_SUBSTEP_VISC_FILTER   = 6
  integer, parameter, public :: DUMP_SUBSTEP_IMPL_VISC     = 7
  integer, parameter, public :: DUMP_SUBSTEP_SSH_RHS       = 8
  integer, parameter, public :: DUMP_SUBSTEP_SSH_SOLVE     = 9
  integer, parameter, public :: DUMP_SUBSTEP_UPDATE_VEL    = 10
  integer, parameter, public :: DUMP_SUBSTEP_HBAR          = 11
  integer, parameter, public :: DUMP_SUBSTEP_ETA_N         = 12
  integer, parameter, public :: DUMP_SUBSTEP_ALE           = 13
  integer, parameter, public :: DUMP_SUBSTEP_GM_BOLUS      = 14
  integer, parameter, public :: DUMP_SUBSTEP_TRACERS       = 15
  integer, parameter, public :: DUMP_SUBSTEP_THICKNESS     = 16

  ! ---- Probe nodes ---------------------------------------------------------
  ! Hardcoded global Fortran 1-based node IDs. Valid on both pi (3140 nodes)
  ! and CORE2 (126858 nodes). A rank silently skips IDs not in its myList.
  ! Expand or move to env-var-driven lat/lon resolution in a later iteration.
  integer, parameter, public :: DUMP_NPROBES = 5
  integer, dimension(DUMP_NPROBES), parameter, public :: DUMP_PROBE_GIDS = &
       [ 1001, 1500, 2000, 2500, 3000 ]

  ! ---- Internal state ------------------------------------------------------
  logical, save :: shim_inited = .false.
  logical, save :: shim_active = .false.
  integer, save :: dump_unit   = -1
  integer, save :: max_steps   = 10
  integer, save :: probe_locals(DUMP_NPROBES) = -1   ! local indices, -1 = absent

  public :: dump_shim_init, dump_shim_record_node, dump_shim_record_node_2d, &
            dump_shim_finalize

contains

  !---------------------------------------------------------------------------
  ! Initialise. Safe to call repeatedly; only the first call has effect.
  ! Reads env vars, resolves probe global IDs to local indices on this rank,
  ! and opens the per-rank dump file. If FESOM_DUMP_FILE is unset, the module
  ! stays a no-op for the lifetime of the run.
  !---------------------------------------------------------------------------
  subroutine dump_shim_init(mype, myDim_nod2D, myList_nod2D)
    integer, intent(in) :: mype, myDim_nod2D
    integer, intent(in) :: myList_nod2D(:)

    character(len=512) :: env_path, env_max, fname
    integer            :: i, nn, ios, env_len

    if (shim_inited) return
    shim_inited = .true.

    call get_environment_variable('FESOM_DUMP_FILE', env_path, length=env_len, status=ios)
    if (ios /= 0 .or. env_len == 0) return

    call get_environment_variable('FESOM_DUMP_MAXSTEPS', env_max, status=ios)
    if (ios == 0 .and. len_trim(env_max) > 0) then
      read(env_max, *, iostat=ios) max_steps
      if (ios /= 0) max_steps = 10
    end if

    do i = 1, DUMP_NPROBES
      probe_locals(i) = -1
      do nn = 1, myDim_nod2D
        if (myList_nod2D(nn) == DUMP_PROBE_GIDS(i)) then
          probe_locals(i) = nn
          exit
        end if
      end do
    end do

    write(fname, '(A,".",I5.5)') trim(env_path), mype
    open(newunit=dump_unit, file=trim(fname), status='replace', &
         form='unformatted', access='stream', action='write', iostat=ios)
    if (ios /= 0) then
      print *, '[fesom_dump_shim] FAILED to open ', trim(fname), ' iostat=', ios
      return
    end if

    shim_active = .true.
    if (mype == 0) then
      print *, '[fesom_dump_shim] active: prefix=', trim(env_path), &
               ' max_steps=', max_steps, ' nprobes=', DUMP_NPROBES
    end if
  end subroutine dump_shim_init

  !---------------------------------------------------------------------------
  ! Record one column-valued node field at all probe nodes resident on this
  ! rank. field_3d is (nl, nod_total); nlevels_nod2D(:) gives the column
  ! length per node. Truncates the column to nlevels_nod2D(local_id) entries.
  !---------------------------------------------------------------------------
  subroutine dump_shim_record_node(substep_id, n, field_name, field_3d, nlevels_nod2D)
    integer,          intent(in) :: substep_id, n
    character(len=*), intent(in) :: field_name
    real(real64),     intent(in) :: field_3d(:,:)     ! (nl, nod_total)
    integer,          intent(in) :: nlevels_nod2D(:)

    integer            :: i, nlev, lid
    character(len=24)  :: name24

    if (.not. shim_active) return
    if (n > max_steps)     return

    name24 = field_name      ! Fortran auto-blank-pads on assignment

    do i = 1, DUMP_NPROBES
      lid = probe_locals(i)
      if (lid <= 0) cycle
      nlev = nlevels_nod2D(lid)
      write(dump_unit) int(n,             int32), &
                       int(substep_id,    int32), &
                       int(DUMP_PROBE_GIDS(i), int32), &
                       int(nlev,          int32), &
                       name24, &
                       field_3d(1:nlev, lid)
    end do
  end subroutine dump_shim_record_node

  !---------------------------------------------------------------------------
  ! Record one scalar-per-node field at all probe nodes resident on this rank.
  ! field_2d is (nod_total). Encoded as nlevels=1 records so the C reader
  ! handles it identically to column fields.
  !---------------------------------------------------------------------------
  subroutine dump_shim_record_node_2d(substep_id, n, field_name, field_2d)
    integer,          intent(in) :: substep_id, n
    character(len=*), intent(in) :: field_name
    real(real64),     intent(in) :: field_2d(:)

    integer            :: i, lid
    character(len=24)  :: name24
    real(real64)       :: scratch(1)

    if (.not. shim_active) return
    if (n > max_steps)     return

    name24 = field_name

    do i = 1, DUMP_NPROBES
      lid = probe_locals(i)
      if (lid <= 0) cycle
      scratch(1) = field_2d(lid)
      write(dump_unit) int(n,             int32), &
                       int(substep_id,    int32), &
                       int(DUMP_PROBE_GIDS(i), int32), &
                       int(1,             int32), &
                       name24, &
                       scratch
    end do
  end subroutine dump_shim_record_node_2d

  !---------------------------------------------------------------------------
  ! Close the dump file. Safe to call when shim is inactive.
  !---------------------------------------------------------------------------
  subroutine dump_shim_finalize()
    if (.not. shim_active) return
    close(dump_unit)
    shim_active = .false.
  end subroutine dump_shim_finalize

end module fesom_dump_shim
