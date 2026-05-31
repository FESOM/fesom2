!> Ship-track / mooring-array curtain output for FESOM-2.
!>
!> See ../runtime/fesom-2.7/meltpond_diag/ship_track_curtain_plan.md for
!> the full design. The model walks each user-supplied polyline CSV,
!> densifies the geodesic between consecutive vertices at
!> track_resolution_km, snaps every dense sample to its nearest
!> strictly-owned mesh node via MPI_MINLOC, collapses consecutive
!> duplicate nodes, and writes the resulting sequence of unique mesh
!> nodes at all depths as a NetCDF "curtain" via XIOS.
!>
!> Multi-track: track_files / track_vars / track_names are
!> semicolon-separated strings. One file_group at a single shared
!> output cadence (track_output_freq), one <file> per track inside it,
!> one runtime XIOS domain per track.
!>
!> Activation: ltracks=.true. plus a non-empty track_files. XML override
!> in io_xios_init wins over the (skipped in XIOS mode) namelist.
module io_tracks_module
  ! Configuration variables are always declared so namelist.io can parse
  ! the &nml_tracks block in any build. The XIOS-using subroutines below
  ! are compiled in only when __XIOS is defined.
#if defined(__XIOS)
  use xios
  use ixml_tree,   only: xios_add_axistogrid, xios_add_domaintogrid
  use mpi
  use mod_mesh,    only: t_mesh
  use mod_partit,  only: t_partit
  use mod_tracer,  only: t_tracer
  use o_param,     only: WP, rad
#endif
  implicit none
  private

  logical,             save, public :: ltracks               = .false.
  character(len=4096), save, public :: track_files           = ''
  character(len=512),  save, public :: track_vars            = 'temp'
  character(len=512),  save, public :: track_names           = ''
  real(kind=8),        save, public :: track_resolution_km   = 5.0_8
  character(len=16),   save, public :: track_output_freq     = '1h'

  public :: io_tracks_is_on
#if defined(__XIOS)
  public :: io_tracks_register_xios
  public :: io_tracks_send

  ! Per-track state. One element per user-supplied CSV.
  type :: track_t
    character(len=64)    :: name       = ''
    character(len=512)   :: csv_path   = ''
    character(len=32)    :: var_name   = 'temp'
    integer              :: tracer_idx = 1            ! 1=temp, 2=salt
    integer              :: N          = 0            ! unique mesh nodes after dedupe
    integer              :: nL         = 0            ! per-rank owned count
    character(len=64)    :: field_id   = ''
    integer, allocatable :: track_pos(:)              ! length nL
    integer, allocatable :: track_lid(:)              ! length nL
  end type track_t

  type(track_t), allocatable, save :: tracks(:)
  integer,                    save :: n_tracks = 0
  integer,                    save :: nz_track = 0   ! shared vertical levels
#endif

contains

  !> True iff tracks output is enabled this run.
  logical function io_tracks_is_on() result(on)
    on = ltracks
  end function io_tracks_is_on


#if defined(__XIOS)


  !> Parse semicolon-separated lists, then for each track: read CSV,
  !> densify the polyline, MPI_MINLOC ownership, dedupe, build per-rank
  !> sparse lists, register XIOS domain/grid/field. One shared file_group
  !> is created once and each track's <file> is added inside it.
  subroutine io_tracks_register_xios(mesh, partit, nz_cell)
    type(t_mesh),   intent(in) :: mesh
    type(t_partit), intent(in) :: partit
    integer,        intent(in) :: nz_cell
    logical, save :: already = .false.
    integer :: i
    integer :: ierr

    if (.not. io_tracks_is_on()) return
    if (already) return
    already = .true.

    nz_track = nz_cell

    call parse_track_config(partit)
    if (n_tracks == 0) then
       if (partit%mype == 0) write(*,'(a)') &
            '[TRACKS] ERROR: ltracks=.true. but track_files is empty.'
       call mpi_abort(partit%MPI_COMM_FESOM, 1, ierr)
    end if

    call register_track_filegroup(partit)

    do i = 1, n_tracks
       call process_one_track(mesh, partit, nz_cell, tracks(i))
    end do
  end subroutine io_tracks_register_xios


  !> Loop tracks, pack the per-track variable from the tracer array at
  !> owned mesh nodes, send via XIOS. Call once per timestep alongside
  !> the existing 0D send block in fesom_module.F90.
  subroutine io_tracks_send(tracers, partit)
    type(t_tracer), intent(in) :: tracers
    type(t_partit), intent(in) :: partit
    integer :: i

    if (.not. io_tracks_is_on()) return
    if (n_tracks == 0) return

    do i = 1, n_tracks
       call send_one_track(tracers, partit, tracks(i))
    end do
  end subroutine io_tracks_send


  ! ====================================================================
  ! Per-track pipeline
  ! ====================================================================

  subroutine process_one_track(mesh, partit, nz_cell, t)
    type(t_mesh),   intent(in)    :: mesh
    type(t_partit), intent(in)    :: partit
    integer,        intent(in)    :: nz_cell
    type(track_t),  intent(inout) :: t

    real(kind=8), allocatable :: lon_csv(:), lat_csv(:)
    real(kind=8), allocatable :: lon_dense(:), lat_dense(:)
    integer,      allocatable :: gid_track(:), owner_track(:), lid_owner(:)
    integer,      allocatable :: i_index(:)
    real(kind=8), allocatable :: lon_out(:), lat_out(:)
    character(len=64) :: domain_id, grid_id, file_id, file_name
    integer :: N_csv, N_dense, N_uniq, ierr

    ! --- 0. tracer slot ---------------------------------------------------
    select case (trim(t%var_name))
    case ('temp')
       t%tracer_idx = 1
    case ('salt')
       t%tracer_idx = 2
    case default
       if (partit%mype == 0) then
          write(*,'(a,a,a,a,a)') '[TRACKS] ERROR: track "', trim(t%name), &
               '" unsupported track_var="', trim(t%var_name),             &
               '". Supported: temp, salt.'
       end if
       call mpi_abort(partit%MPI_COMM_FESOM, 1, ierr)
    end select
    t%field_id = trim(t%var_name) // '_track_' // trim(t%name)

    ! --- 1. parse + densify on rank 0, bcast dense (lon, lat) -------------
    if (partit%mype == 0) then
       call parse_csv(trim(t%csv_path), N_csv, lon_csv, lat_csv, partit)
       write(*,'(a,a,a,a,a,i0,a)') '[TRACKS] ', trim(t%name), ': parsed ', &
            trim(t%csv_path), ' -> ', N_csv, ' polyline vertices'
       call densify_polyline(N_csv, lon_csv, lat_csv, track_resolution_km, &
                             N_dense, lon_dense, lat_dense)
       write(*,'(a,a,a,f6.2,a,i0,a)') '[TRACKS] ', trim(t%name),           &
            ': densified at ', track_resolution_km, ' km -> ',             &
            N_dense, ' query points'
       deallocate(lon_csv, lat_csv)
    end if
    call mpi_bcast(N_dense, 1, MPI_INTEGER, 0, partit%MPI_COMM_FESOM, ierr)
    if (partit%mype /= 0) allocate(lon_dense(N_dense), lat_dense(N_dense))
    call mpi_bcast(lon_dense, N_dense, MPI_DOUBLE_PRECISION, 0, partit%MPI_COMM_FESOM, ierr)
    call mpi_bcast(lat_dense, N_dense, MPI_DOUBLE_PRECISION, 0, partit%MPI_COMM_FESOM, ierr)

    ! --- 2. nearest-neighbour resolution + dedupe -------------------------
    allocate(gid_track(N_dense), owner_track(N_dense), lid_owner(N_dense))
    call resolve_owners(mesh, partit, N_dense, lon_dense, lat_dense,       &
                        gid_track, owner_track, lid_owner)
    deallocate(lon_dense, lat_dense)
    call dedupe_adjacent(N_dense, gid_track, owner_track, lid_owner, N_uniq)
    t%N = N_uniq
    if (partit%mype == 0) then
       write(*,'(a,a,a,i0,a)') '[TRACKS] ', trim(t%name),                  &
            ': deduplicated to ', t%N, ' unique mesh nodes'
    end if

    ! --- 3. per-rank sparse lists, monotonic by track_pos -----------------
    call build_local_lists(partit, t%N, owner_track, lid_owner,            &
                           t%nL, t%track_pos, t%track_lid)

    ! --- 4. XIOS attr arrays for this rank --------------------------------
    allocate(i_index(t%nL), lon_out(t%nL), lat_out(t%nL))
    call build_xios_attrs(mesh, t, lon_out, lat_out, i_index)

    ! --- 5. domain + grid + field + add file in the shared file_group -----
    domain_id = 'track_' // trim(t%name)
    grid_id   = 'grid_3d_track_' // trim(t%name)
    file_id   = 'f_' // trim(t%field_id)
    file_name = trim(t%field_id) // '.fesom'

    call register_xios_objects_for_track(partit, t%N, t%nL,                &
                                         trim(domain_id), trim(grid_id),   &
                                         trim(t%field_id), trim(file_id),  &
                                         trim(file_name),                  &
                                         i_index, lon_out, lat_out)

    if (partit%mype == 0) then
       write(*,'(a,a,a,a,a,i0,a)') '[TRACKS] ', trim(t%name),              &
            ': registered runtime XML (domain=', trim(domain_id),          &
            ', N=', t%N, ' cells)'
       write(*,'(a,a,a,a)')        '[TRACKS]   field=', trim(t%field_id),  &
                                   '  file=', trim(file_name)
    end if

    deallocate(gid_track, owner_track, lid_owner)
    deallocate(i_index, lon_out, lat_out)
  end subroutine process_one_track


  subroutine send_one_track(tracers, partit, t)
    type(t_tracer), intent(in) :: tracers
    type(t_partit), intent(in) :: partit
    type(track_t),  intent(in) :: t
    real(kind=8), allocatable  :: buf(:,:)
    integer :: j, kk

    if (len_trim(t%field_id) == 0)                            return
    if (.not. xios_is_valid_field(trim(t%field_id)))          return
    if (.not. xios_field_is_active(trim(t%field_id), .true.)) return

    allocate(buf(nz_track, t%nL))
    do j = 1, t%nL
       do kk = 1, nz_track
          buf(kk, j) = real(tracers%data(t%tracer_idx)%values(kk, t%track_lid(j)), kind=8)
       end do
    end do
    call xios_send_field(trim(t%field_id), buf)
    deallocate(buf)
  end subroutine send_one_track


  ! ====================================================================
  ! Config parsing
  ! ====================================================================

  !> Parse the comma- (or semicolon-) separated config strings into the
  !> tracks(:) array. All three lists must be the same length; track_names
  !> defaults to "track<i>" if left empty.
  subroutine parse_track_config(partit)
    type(t_partit), intent(in) :: partit
    character(len=512), allocatable :: files(:), vars(:), names(:)
    integer :: nf, nv, nn, i, ierr
    character(len=4) :: tag

    call split_semi_or_comma(track_files, nf, files)
    call split_semi_or_comma(track_vars,  nv, vars)
    call split_semi_or_comma(track_names, nn, names)

    n_tracks = nf
    if (n_tracks == 0) return

    if (nv == 1 .and. n_tracks > 1) then
       ! Single var supplied -> broadcast to all tracks
       deallocate(vars)
       allocate(vars(n_tracks))
       vars(:) = trim(track_vars)
       nv = n_tracks
    end if
    if (nv /= n_tracks .and. partit%mype == 0) then
       write(*,'(a,i0,a,i0)') '[TRACKS] ERROR: track_vars has ', nv,       &
            ' entries but track_files has ', n_tracks
    end if

    if (nn == 0) then
       ! Default names: "track1", "track2", ...
       allocate(names(n_tracks))
       do i = 1, n_tracks
          write(tag, '(i0)') i
          names(i) = 'track' // trim(tag)
       end do
    else if (nn /= n_tracks .and. partit%mype == 0) then
       write(*,'(a,i0,a,i0)') '[TRACKS] ERROR: track_names has ', nn,      &
            ' entries but track_files has ', n_tracks
    end if

    if (nv /= n_tracks .or. (nn /= 0 .and. nn /= n_tracks)) then
       call mpi_abort(partit%MPI_COMM_FESOM, 1, ierr)
    end if

    if (allocated(tracks)) deallocate(tracks)
    allocate(tracks(n_tracks))
    do i = 1, n_tracks
       tracks(i)%name     = trim(adjustl(names(i)))
       tracks(i)%csv_path = trim(adjustl(files(i)))
       tracks(i)%var_name = trim(adjustl(vars(i)))
    end do
    if (allocated(files)) deallocate(files)
    if (allocated(vars))  deallocate(vars)
    if (allocated(names)) deallocate(names)

    if (partit%mype == 0) then
       write(*,'(a,i0,a)') '[TRACKS] parsed config: ', n_tracks, ' track(s):'
       do i = 1, n_tracks
          write(*,'(a,i0,a,a,a,a,a,a,a)') '[TRACKS]   ', i, '. name="',    &
               trim(tracks(i)%name), '" var="', trim(tracks(i)%var_name),  &
               '" csv="', trim(tracks(i)%csv_path), '"'
       end do
       write(*,'(a,a,a,f6.2,a)') '[TRACKS] shared cadence=',               &
            trim(track_output_freq), ', resolution=',                       &
            track_resolution_km, ' km'
    end if
  end subroutine parse_track_config


  !> Split a string on ';' or ',' (whichever appears first); trim tokens.
  !> Returns n=0 and unallocated tokens if input is empty.
  subroutine split_semi_or_comma(s_in, n, tokens)
    character(len=*),               intent(in)  :: s_in
    integer,                        intent(out) :: n
    character(len=512), allocatable, intent(out) :: tokens(:)
    character(len=len(s_in)) :: s
    character :: sep
    integer :: i, j, k, lenT, ntok

    s = adjustl(s_in)
    lenT = len_trim(s)
    if (lenT == 0) then
       n = 0
       return
    end if
    if (index(s(1:lenT), ';') > 0) then
       sep = ';'
    else
       sep = ','
    end if

    ntok = 1
    do i = 1, lenT
       if (s(i:i) == sep) ntok = ntok + 1
    end do
    n = ntok
    allocate(tokens(n))

    j = 1
    k = 1
    do i = 1, lenT + 1
       if (i > lenT .or. s(i:i) == sep) then
          tokens(k) = adjustl(s(j:i-1))
          j = i + 1
          k = k + 1
       end if
    end do
  end subroutine split_semi_or_comma


  !> Parse a cadence string like '1h', '6h', '1d', '1mo', '1y' into the
  !> matching xios_duration. Default fallback: 1 hour.
  function parse_xios_freq(s) result(d)
    character(len=*), intent(in) :: s
    type(xios_duration) :: d
    character(len=16) :: t
    integer :: n, idx, ios
    t = adjustl(s)
    idx = index(t, 'mo')
    if (idx > 0) then
       read(t(1:idx-1), *, iostat=ios) n
       if (ios == 0) then; d = xios_duration(month=real(n, 8)); return; end if
    end if
    idx = index(t, 'y')
    if (idx > 0) then
       read(t(1:idx-1), *, iostat=ios) n
       if (ios == 0) then; d = xios_duration(year=real(n, 8)); return; end if
    end if
    idx = index(t, 'd')
    if (idx > 0) then
       read(t(1:idx-1), *, iostat=ios) n
       if (ios == 0) then; d = xios_duration(day=real(n, 8)); return; end if
    end if
    idx = index(t, 'h')
    if (idx > 0) then
       read(t(1:idx-1), *, iostat=ios) n
       if (ios == 0) then; d = xios_duration(hour=real(n, 8)); return; end if
    end if
    d = xios_duration(hour=1.0_8)
  end function parse_xios_freq


  ! ====================================================================
  ! XIOS runtime registration
  ! ====================================================================

  !> Create the shared file_group "fesom_tracks" once with the global
  !> output_freq + split_freq + enabled=true. Each per-track <file> is
  !> later added inside this group.
  subroutine register_track_filegroup(partit)
    type(t_partit), intent(in) :: partit
    type(xios_filegroup) :: filgrp_root, my_filgrp

    call xios_get_handle("file_definition", filgrp_root)
    call xios_add_child(filgrp_root, my_filgrp, "fesom_tracks")
    call xios_set_attr(my_filgrp,                                          &
         output_freq = parse_xios_freq(trim(track_output_freq)),            &
         split_freq  = xios_duration(year = 1.0_8),                         &
         enabled     = .true.)
    if (partit%mype == 0) then
       write(*,'(a,a,a)') '[TRACKS] file_group "fesom_tracks" created with output_freq=', &
                          trim(track_output_freq), ' split_freq=1y'
    end if
  end subroutine register_track_filegroup


  !> Per-track: register domain, grid, field, plus a <file> inside the
  !> existing "fesom_tracks" file_group with field_ref + freq_op.
  subroutine register_xios_objects_for_track(partit, N_glo, nL_loc,        &
                                             domain_id, grid_id,            &
                                             field_id_in, file_id, file_name, &
                                             i_index, lon_out, lat_out)
    type(t_partit),   intent(in) :: partit
    integer,          intent(in) :: N_glo, nL_loc
    character(len=*), intent(in) :: domain_id, grid_id
    character(len=*), intent(in) :: field_id_in, file_id, file_name
    integer,          intent(in) :: i_index(nL_loc)
    real(kind=8),     intent(in) :: lon_out(nL_loc), lat_out(nL_loc)

    type(xios_domaingroup) :: dgrp
    type(xios_domain)      :: dom, dom_in_grid
    type(xios_gridgroup)   :: ggrp
    type(xios_grid)        :: grd
    type(xios_axis)        :: ax_in_grid
    type(xios_fieldgroup)  :: fgrp
    type(xios_field)       :: fld
    type(xios_filegroup)   :: my_filgrp
    type(xios_file)        :: fil
    type(xios_field)       :: fld_in_file
    type(xios_duration)    :: freq

    freq = parse_xios_freq(trim(track_output_freq))

    ! domain
    call xios_get_handle("domain_definition", dgrp)
    call xios_add_child(dgrp, dom, domain_id)
    call xios_set_attr(dom,                              &
         type        = "unstructured",                    &
         ni_glo      = N_glo,                              &
         ni          = nL_loc,                             &
         i_index     = i_index,                            &
         lonvalue_1d = lon_out,                            &
         latvalue_1d = lat_out,                            &
         data_dim    = 1,                                  &
         data_ni     = nL_loc)

    ! grid (nz x domain)
    call xios_get_handle("grid_definition", ggrp)
    call xios_add_child(ggrp, grd, grid_id)
    call xios_add_axistogrid(grd, ax_in_grid)
    call xios_set_attr(ax_in_grid, axis_ref = "nz")
    call xios_add_domaintogrid(grd, dom_in_grid)
    call xios_set_attr(dom_in_grid, domain_ref = domain_id)

    ! field
    call xios_get_handle("field_definition", fgrp)
    call xios_add_child(fgrp, fld, field_id_in)
    call xios_set_attr(fld,                              &
         grid_ref  = grid_id,                             &
         operation = "average",                           &
         prec      = 8,                                   &
         long_name = "ship-track/mooring curtain output")

    ! file inside the shared fesom_tracks group
    call xios_get_handle("fesom_tracks", my_filgrp)
    call xios_add_child(my_filgrp, fil, file_id)
    call xios_set_attr(fil, name = file_name)
    call xios_add_child(fil, fld_in_file)
    call xios_set_attr(fld_in_file,                                        &
         field_ref = field_id_in,                                           &
         operation = "instant",                                             &
         freq_op   = freq)

    call xios_solve_inheritance()
  end subroutine register_xios_objects_for_track


  ! ====================================================================
  ! Private helpers (CSV, geometry, MPI, per-rank packing)
  ! ====================================================================

  subroutine parse_csv(path, N, lon_out, lat_out, partit)
    character(len=*),              intent(in)  :: path
    integer,                       intent(out) :: N
    real(kind=8), allocatable,     intent(out) :: lon_out(:), lat_out(:)
    type(t_partit),                intent(in)  :: partit

    integer, parameter :: MAX_WPTS = 100000
    real(kind=8) :: lon_buf(MAX_WPTS), lat_buf(MAX_WPTS)
    character(len=512) :: line
    integer :: u, ios, line_no, k
    integer :: comma1, ierr
    character(len=128) :: ftok, stok
    real(kind=8) :: lon, lat

    open(newunit=u, file=trim(path), status='old', action='read', iostat=ios)
    if (ios /= 0) then
       write(*,'(a,a,a,i0)') '[TRACKS] ERROR: cannot open track file "',   &
                              trim(path), '" iostat=', ios
       call mpi_abort(partit%MPI_COMM_FESOM, 1, ierr)
    end if

    N = 0
    line_no = 0
    do
       read(u, '(a)', iostat=ios) line
       if (ios /= 0) exit
       line_no = line_no + 1

       if (line_no == 1 .and. len(line) >= 3) then
          if (iachar(line(1:1)) == 239 .and. &
              iachar(line(2:2)) == 187 .and. &
              iachar(line(3:3)) == 191) then
             line = line(4:)
          end if
       end if

       line = adjustl(line)
       if (len_trim(line) == 0) cycle
       if (line(1:1) == '#') cycle

       comma1 = index(line, ',')
       if (comma1 < 2) then
          call csv_fatal(path, line_no, 'expected "," separator')
          call mpi_abort(partit%MPI_COMM_FESOM, 1, ierr)
       end if
       ftok = adjustl(line(1:comma1-1))
       stok = adjustl(line(comma1+1:))

       if (index(ftok, ',') > 0) then
          call csv_fatal(path, line_no, 'decimal point required (no decimal comma)')
          call mpi_abort(partit%MPI_COMM_FESOM, 1, ierr)
       end if

       read(ftok, *, iostat=ios) lon
       if (ios /= 0) then
          call csv_fatal(path, line_no, 'cannot parse lon')
          call mpi_abort(partit%MPI_COMM_FESOM, 1, ierr)
       end if
       read(stok, *, iostat=ios) lat
       if (ios /= 0) then
          call csv_fatal(path, line_no, 'cannot parse lat')
          call mpi_abort(partit%MPI_COMM_FESOM, 1, ierr)
       end if

       if (N >= MAX_WPTS) then
          write(*,'(a,a,a,i0)') '[TRACKS] ERROR: ', trim(path),            &
                                 ' exceeds MAX_WPTS=', MAX_WPTS
          call mpi_abort(partit%MPI_COMM_FESOM, 1, ierr)
       end if

       do k = 1, N
          if (abs(lon_buf(k) - lon) < 1.0e-12_8 .and. &
              abs(lat_buf(k) - lat) < 1.0e-12_8) then
             write(*,'(a,a,a,i0,a,i0)') '[TRACKS] ERROR: ', trim(path),    &
                  ' line ', line_no, ' duplicates waypoint ', k
             call mpi_abort(partit%MPI_COMM_FESOM, 1, ierr)
          end if
       end do

       N = N + 1
       lon_buf(N) = lon
       lat_buf(N) = lat
    end do
    close(u)

    if (N == 0) then
       write(*,'(a,a,a)') '[TRACKS] ERROR: ', trim(path), ' has no waypoints'
       call mpi_abort(partit%MPI_COMM_FESOM, 1, ierr)
    end if

    allocate(lon_out(N), lat_out(N))
    lon_out(1:N) = lon_buf(1:N)
    lat_out(1:N) = lat_buf(1:N)
  end subroutine parse_csv


  subroutine csv_fatal(path, line_no, msg)
    character(len=*), intent(in) :: path, msg
    integer,          intent(in) :: line_no
    write(*,'(a,a,a,i0,a,a)') '[TRACKS] ERROR: ', trim(path), &
                              ':', line_no, ': ', trim(msg)
  end subroutine csv_fatal


  subroutine resolve_owners(mesh, partit, N, lon_in, lat_in,               &
                            gid_track, owner_track, lid_owner)
    type(t_mesh),   intent(in)  :: mesh
    type(t_partit), intent(in)  :: partit
    integer,        intent(in)  :: N
    real(kind=8),   intent(in)  :: lon_in(N), lat_in(N)
    integer,        intent(out) :: gid_track(N), owner_track(N), lid_owner(N)
    integer :: k, i, ierr
    real(kind=8) :: local_best_dist, d, lon_k, lat_k
    integer :: local_best_lid
    real(kind=8) :: pair_in(2), pair_out(2)
    integer :: bcast_buf(2)

    lid_owner(:) = 0
    do k = 1, N
       lon_k = lon_in(k)
       lat_k = lat_in(k)
       local_best_dist = huge(0.0_8)
       local_best_lid  = 0
       do i = 1, partit%myDim_nod2D
          d = gcdist_km(lon_k, lat_k,                                       &
                        real(mesh%geo_coord_nod2D(1, i), kind=8) / rad,    &
                        real(mesh%geo_coord_nod2D(2, i), kind=8) / rad)
          if (d < local_best_dist) then
             local_best_dist = d
             local_best_lid  = i
          end if
       end do
       pair_in(1) = local_best_dist
       pair_in(2) = real(partit%mype, kind=8)
       call mpi_allreduce(pair_in, pair_out, 1, MPI_2DOUBLE_PRECISION,     &
                          MPI_MINLOC, partit%MPI_COMM_FESOM, ierr)
       owner_track(k) = nint(pair_out(2))
       if (partit%mype == owner_track(k)) then
          bcast_buf(1) = partit%myList_nod2D(local_best_lid)
          bcast_buf(2) = local_best_lid
       else
          bcast_buf(:) = 0
       end if
       call mpi_bcast(bcast_buf, 2, MPI_INTEGER, owner_track(k),           &
                      partit%MPI_COMM_FESOM, ierr)
       gid_track(k) = bcast_buf(1)
       if (partit%mype == owner_track(k)) lid_owner(k) = bcast_buf(2)
    end do
  end subroutine resolve_owners


  function gcdist_km(lon1, lat1, lon2, lat2) result(d)
    real(kind=8), intent(in) :: lon1, lat1, lon2, lat2
    real(kind=8) :: d
    real(kind=8), parameter :: deg2rad = 3.141592653589793_8 / 180.0_8
    real(kind=8) :: phi1, phi2, dl, c
    phi1 = lat1 * deg2rad
    phi2 = lat2 * deg2rad
    dl   = (lon2 - lon1) * deg2rad
    c    = sin(phi1)*sin(phi2) + cos(phi1)*cos(phi2)*cos(dl)
    c    = max(-1.0_8, min(1.0_8, c))
    d    = 6371.0_8 * acos(c)
  end function gcdist_km


  subroutine densify_polyline(N_in, lon_in, lat_in, res_km,                &
                              N_out, lon_out, lat_out)
    integer,                   intent(in)  :: N_in
    real(kind=8),              intent(in)  :: lon_in(N_in), lat_in(N_in)
    real(kind=8),              intent(in)  :: res_km
    integer,                   intent(out) :: N_out
    real(kind=8), allocatable, intent(out) :: lon_out(:), lat_out(:)

    real(kind=8), parameter :: deg2rad = 3.141592653589793_8 / 180.0_8
    real(kind=8) :: lon1, lat1, lon2, lat2, d_km
    real(kind=8) :: p1(3), p2(3), pv(3), angle, sinang, t
    integer :: seg, ns, j, k, cap
    real(kind=8), allocatable :: buf_lon(:), buf_lat(:)

    if (N_in <= 0) then
       N_out = 0
       allocate(lon_out(0), lat_out(0))
       return
    end if
    if (N_in == 1) then
       N_out = 1
       allocate(lon_out(1), lat_out(1))
       lon_out(1) = lon_in(1)
       lat_out(1) = lat_in(1)
       return
    end if

    cap = N_in + max(1, int(40000.0_8 / max(res_km, 0.1_8))) * (N_in - 1)
    allocate(buf_lon(cap), buf_lat(cap))

    k = 1
    buf_lon(k) = lon_in(1)
    buf_lat(k) = lat_in(1)

    do seg = 2, N_in
       lon1 = lon_in(seg-1); lat1 = lat_in(seg-1)
       lon2 = lon_in(seg);   lat2 = lat_in(seg)
       d_km = gcdist_km(lon1, lat1, lon2, lat2)
       ns = max(1, int(ceiling(d_km / max(res_km, 1.0e-3_8))))

       p1(1) = cos(lat1*deg2rad) * cos(lon1*deg2rad)
       p1(2) = cos(lat1*deg2rad) * sin(lon1*deg2rad)
       p1(3) = sin(lat1*deg2rad)
       p2(1) = cos(lat2*deg2rad) * cos(lon2*deg2rad)
       p2(2) = cos(lat2*deg2rad) * sin(lon2*deg2rad)
       p2(3) = sin(lat2*deg2rad)
       angle  = acos(max(-1.0_8, min(1.0_8, p1(1)*p2(1)+p1(2)*p2(2)+p1(3)*p2(3))))
       sinang = sin(angle)

       do j = 1, ns
          t = real(j, kind=8) / real(ns, kind=8)
          if (sinang < 1.0e-12_8) then
             pv = p2
          else
             pv = (sin((1.0_8 - t) * angle) * p1                           &
                 + sin(t          * angle) * p2) / sinang
          end if
          k = k + 1
          if (k > cap) then
             write(*,'(a)') '[TRACKS] BUG: densify buffer overflow'
             exit
          end if
          buf_lat(k) = asin(max(-1.0_8, min(1.0_8, pv(3)))) / deg2rad
          buf_lon(k) = atan2(pv(2), pv(1)) / deg2rad
       end do
    end do

    N_out = k
    allocate(lon_out(N_out), lat_out(N_out))
    lon_out = buf_lon(1:N_out)
    lat_out = buf_lat(1:N_out)
    deallocate(buf_lon, buf_lat)
  end subroutine densify_polyline


  subroutine dedupe_adjacent(N_in, gid, owner, lid, N_out)
    integer, intent(in)    :: N_in
    integer, intent(inout) :: gid(N_in), owner(N_in), lid(N_in)
    integer, intent(out)   :: N_out
    integer :: i
    if (N_in == 0) then
       N_out = 0
       return
    end if
    N_out = 1
    do i = 2, N_in
       if (gid(i) /= gid(N_out)) then
          N_out = N_out + 1
          gid(N_out)   = gid(i)
          owner(N_out) = owner(i)
          lid(N_out)   = lid(i)
       end if
    end do
  end subroutine dedupe_adjacent


  subroutine build_local_lists(partit, N, owner_track, lid_owner,          &
                               nL_out, track_pos_out, track_lid_out)
    type(t_partit), intent(in)  :: partit
    integer,        intent(in)  :: N, owner_track(N), lid_owner(N)
    integer,        intent(out) :: nL_out
    integer,    allocatable, intent(out) :: track_pos_out(:), track_lid_out(:)
    integer :: k, j, ierr

    nL_out = count(owner_track == partit%mype)
    allocate(track_pos_out(nL_out), track_lid_out(nL_out))

    j = 0
    do k = 1, N
       if (owner_track(k) /= partit%mype) cycle
       j = j + 1
       track_pos_out(j) = k
       track_lid_out(j) = lid_owner(k)
    end do
    do j = 2, nL_out
       if (track_pos_out(j) < track_pos_out(j-1)) then
          write(*,'(a)') '[TRACKS] BUG: track_pos out of order; sort needed.'
          call mpi_abort(partit%MPI_COMM_FESOM, 1, ierr)
       end if
    end do
  end subroutine build_local_lists


  subroutine build_xios_attrs(mesh, t, lon_out, lat_out, i_index)
    type(t_mesh),  intent(in)  :: mesh
    type(track_t), intent(in)  :: t
    real(kind=8),  intent(out) :: lon_out(t%nL), lat_out(t%nL)
    integer,       intent(out) :: i_index(t%nL)
    integer :: j, lid
    do j = 1, t%nL
       lid = t%track_lid(j)
       lon_out(j) = real(mesh%geo_coord_nod2D(1, lid), kind=8) / rad
       lat_out(j) = real(mesh%geo_coord_nod2D(2, lid), kind=8) / rad
       i_index(j) = t%track_pos(j) - 1
    end do
  end subroutine build_xios_attrs

#endif
end module io_tracks_module
