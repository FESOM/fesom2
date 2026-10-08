module mod_io_means
    ! FESOM3 field output as per-variable-per-period Zarr v2 stores (M9 Stage 2, Tasks 2.1-2.6).
    !
    ! Each registered field gets its OWN store, split by year (or month): <out_dir>/<name>.fesom.<P>.zarr
    ! (<P> = <YYYY> or <YYYY>_<MM>), holding the data var (time[,nz],nod2) + a CF `time` coord + embedded
    ! `lon`/`lat` node coords (+ `nz`/`nz1` for 3-D). Mirrors FESOM2's per-variable-per-year mean files
    ! (<name>.fesom.<year>.nc), but as Zarr.
    !
    ! Reuses the proven Stage-1 stack verbatim:
    !   - mod_io_decomp : canonical-order MPI_Alltoallv redistribution -> writer subset (no rank-0
    !                     gather), so the store is partition-INDEPENDENT (dist_2 == dist_8).
    !   - mod_io_zarr   : the Zarr v2 writer; store-create ordering = rank 0 creates store + defines
    !                     all arrays -> barrier -> writers write their chunks.
    !
    ! CONFIGURABLE (Task 2.6) — knobs from `namelist.io` (&nml_general + &nml_list), FESOM3_* env overrides:
    !   GLOBAL  : n_writers, chunk_time, chunk_vert, chunk_horiz, compressor (none|lz4), filesplit_freq
    !             (y|m), vec_frame (geographic|native).
    !   PER FIELD: precision (4|8 -> <f4|<f8), mean|snap, and the output cadence freq + unit (y/m/d/h/s),
    !             evaluated each step by the FESOM2-ported events (annual/monthly/daily/hourly/step).
    !
    ! GROWING time dim: each field keeps its OWN record counter `t` (per period) and store handles. With
    ! chunk_time=1 a record == one fresh data chunk [t,..,c] + time chunk [t] + a .zarray shape[0] bump
    ! (zarr_rewrite_zarray). With chunk_time=N>1 the record lands in time-chunk t/N slot mod(t,N): the
    ! partial chunk already on disk is READ (zarr_read_chunk), the slot overwritten, the whole chunk
    ! re-written (read-modify-write). chunk_vert splits a 3-D field's nz dim into ceil(nz/cv) vert chunks.
    !
    ! MEAN vs SNAPSHOT (Task 2.2) — transcribed from io_meandata.F90:update_means/compute_means:
    !   - every step `means_accumulate` does, IN THE OUTPUT PRECISION, mean: acc += value, count++ ;
    !     snapshot: acc = value, count = 1  (io_meandata.F90:2107/2127/2142).
    !   - at output `means_output` writes acc / count (divide in the output precision: real64 for <f8,
    !     real32 for <f4 — io_meandata.F90:2335/2353), then zeroes acc + count for the next interval.
    !   So float32 means accumulate AND divide in float32, byte-matching FESOM2's r4 mean stream.
    !
    ! Field types: node SCALARS 2-D (ssh/sst/...) + 3-D (T/S/w), and node VECTOR pairs (unod/vnod, the
    ! velocity at nodes — FESOM2's default velocity output, dynamics%uvnode) which are r2g-rotated to
    ! geographic at write time (Task 2.5, means_define_vector3d). Element-based u/v (dynamics%uv) is the
    ! alternative FESOM2 offers but never enables by default — left as a future elem-decomp addition.
    !
    ! Caller protocol (all ranks, same order):
    !   means_read_namelist(path, cfg, list, nlist, ok)  optional: parse namelist.io -> global cfg + rows
    !   means_init(io, dir, mesh, partit, ...)            once (cfg knobs as args; FESOM3_* env overrides)
    !   means_define_node2d(io, name, long, units, ...)   per field, once (mean / precision / freq / unit)
    !   ... each step ...
    !   means_accumulate(io, 'ssh', eta_n(1:myDim)) ...   EVERY step (mean sums; snapshot overwrites)
    !   means_output(io, istep, clk)                      EVERY step: events decide which fields write now
    use mpi
    use, intrinsic :: iso_fortran_env, only: real64, real32
    use mod_precision,   only: WP
    use mod_constants,   only: rad
    use mod_mesh,        only: t_mesh
    use mod_partit,      only: t_partit
    use mod_part_bounds, only: is_multirank, local_dims
    use mod_mesh_rotate, only: vector_r2g
    use mod_io_zarr
    use mod_io_decomp
    use mod_io_coords    ! shared lon/lat + _ARRAY_DIMENSIONS + UGRID-attr embedding (Task 2.1)
    implicit none
    private
    public :: t_io_means, t_io_config, t_io_entry, t_means_clock
    public :: means_init, means_read_namelist
    public :: means_define_node2d, means_define_node3d, means_define_vector3d
    public :: means_define_elem2d, means_define_elem3d, means_define_vector3d_elem
    public :: means_accumulate, means_has, means_output, means_finalize
    public :: MEANS_MAXF

    interface means_accumulate
        module procedure means_accumulate_2d, means_accumulate_3d
    end interface means_accumulate

    integer, parameter :: MEANS_MAXF = 128
    ! netCDF standard fill values (== FESOM2 io_meandata.F90 NC_FILL_*); below-bottom levels get this,
    ! xarray masks _FillValue -> NaN (CF). Plan Task 2.4: nlevels-based masking (cleaner than FESOM2's
    ! value-based abs(acc)<1e-30 quirk; valid-level VALUES still byte-match FESOM2).
    real(real64), parameter :: NC_FILL_DOUBLE = 9.9692099683868690e+36_real64
    real(real32), parameter :: NC_FILL_FLOAT  = 9.9692099683868690e+36_real32

    ! ---- namelist.io schema (Task 2.6) -----------------------------------------------------------
    ! &nml_general global knobs (FESOM3 analog of FESOM2 io_meandata.F90:nml_general).
    type :: t_io_config
        integer           :: n_writers      = 0
        integer           :: chunk_time     = 1
        integer           :: chunk_vert     = 0          ! 0 => full depth single chunk
        integer           :: chunk_horiz    = 500000
        character(len=16) :: compressor     = 'none'     ! 'none' | 'lz4'
        character(len=1)  :: filesplit_freq = 'y'        ! 'y' | 'm'
        character(len=16) :: vec_frame      = 'geographic'
    end type t_io_config

    ! &nml_list row: '<id>', <freq>, '<unit>', <precision 4|8>, '<mean|snap>'. Component ORDER matches
    ! the row order so the Fortran namelist read fills one entry per 5 values (FESOM2's io_entry trick,
    ! extended with the 5th 'op' field). Default init => unread tail rows stay 'unknown' (the sentinel).
    type :: t_io_entry
        character(len=15) :: id        = 'unknown'
        integer           :: freq      = 1
        character(len=1)  :: unit      = 's'
        integer           :: precision = 4
        character(len=4)  :: op        = 'snap'
    end type t_io_entry

    ! Clock snapshot the per-field output events read (FESOM3-explicit instead of g_clock globals, so
    ! mod_io_means stays clock-independent and the smoke driver can synthesize it). ndim_month =
    ! num_day_in_month(fleapyear, month) (the caller resolves the leap lookup).
    type :: t_means_clock
        integer      :: year = 0, yearstart = 0, daynew = 1, ndpyr = 365
        integer      :: month = 1, day_in_month = 1, ndim_month = 31
        real(real64) :: timenew = 0.0_real64
    end type t_means_clock

    ! one registered output field (its per-period store + array handles + accumulator + cadence)
    type :: t_mean_field
        character(len=64)  :: name = '', long_name = '', units = '', std = ''
        character(len=8)   :: dtype = '<f4'        ! float32 default (Task 2.6 exposes the knob)
        logical            :: is_mean = .false.     ! mean stream (Task 2.2) vs snapshot (Task 2.1)
        integer            :: ndim = 2              ! 2 = node scalar (time,nod2); 3 = (time,nz,nod2)
        integer            :: nlev = 1              ! vertical size (3-D): nl-1 (layers) or nl (levels)
        integer            :: voff = 0              ! valid levels per entity = nlevels[_nod2D] - voff
        character(len=8)   :: vdim = ''             ! 'nz1' (layers) / 'nz' (levels)
        logical            :: is_elem = .false.     ! .true. => element-based (Task 2.7); else node-based
        character(len=8)   :: hdim = 'nod2'         ! horizontal dim name: 'nod2' (node) / 'elem' (element)
        integer            :: freq = 1              ! output cadence (Task 2.6): freq + unit -> event
        character(len=1)   :: unit = 's'            ! 'y'|'m'|'d'|'h'|'s'
        type(t_zarr_store) :: store                 ! current-period store
        type(t_zarr_array) :: a_data, a_time, a_lon, a_lat, a_vert
        logical            :: open  = .false.       ! store created for the current period
        integer            :: cur_period = -2000000000 ! sentinel "no period opened yet"
        integer            :: t = 0                  ! record index within the current period (0-based)
        ! owned accumulator, in the OUTPUT precision (FESOM2-faithful) — one set per field
        real(real64), allocatable :: acc8(:),    acc8_3d(:,:)
        real(real32), allocatable :: acc4(:),    acc4_3d(:,:)
        integer                   :: count = 0
        ! vector pairing (Task 2.5): the (x,y) components of a velocity/wind pair are r2g-rotated
        ! TOGETHER at write time (io_meandata.F90:io_r2g). is_vec_x marks the x-component (drives the
        ! paired write); is_vec_y the y-component (its data is written by its partner). vec_partner =
        ! the other component's field index.
        logical :: is_vec_x = .false., is_vec_y = .false.
        integer :: vec_partner = 0
    end type t_mean_field

    type :: t_io_means
        character(len=512) :: out_dir  = '.'
        character(len=32)  :: calendar = 'standard'
        integer            :: nf = 0
        type(t_mean_field) :: f(MEANS_MAXF)
        type(t_io_decomp)  :: Dn                     ! node decomp (node fields)
        type(t_io_decomp)  :: De                     ! element decomp (Task 2.7, element fields)
        integer            :: nNodO = 0, nElemO = 0
        integer            :: mype = 0, comm = MPI_COMM_SELF, npes = 1
        logical            :: mr = .false.
        integer            :: nl = 0                 ! mesh%nl (level count)
        ! global writer knobs (Task 2.6) — set in means_init from cfg args; FESOM3_* env overrides
        integer            :: chunk_time = 1         ! time-chunk length (>1 => RMW append)
        integer            :: chunk_vert = 0         ! vertical chunk (0 => full depth single chunk)
        character(len=16)  :: compressor = 'none'    ! data-array codec: 'none' | 'lz4'
        character(len=1)   :: filesplit  = 'y'       ! 'y' per-year / 'm' per-month stores
        real(WP), allocatable :: lon_owned(:), lat_owned(:)  ! cached owned node coords (deg, geographic)
        real(WP), allocatable :: rlon_owned(:), rlat_owned(:)! cached owned ROTATED node coords (rad) for r2g
        integer,  allocatable :: nlev_owned(:)               ! cached owned nlevels_nod2D (3-D node mask)
        ! Task 2.7 element coords: geographic centroid (deg, embed) + ROTATED centroid (rad, r2g) +
        ! element nlevels. The centroid = mod_mesh_areas elem_center of the ROTATED vertices (area-
        ! weighted, cyclic-wrapped; FESOM2 io_r2g:3004 took the plain 3-vertex mean, see mod_io_coords).
        real(WP), allocatable :: elon_owned(:), elat_owned(:)
        real(WP), allocatable :: relon_owned(:), relat_owned(:)
        integer,  allocatable :: enlev_owned(:)              ! cached owned nlevels (3-D element mask)
        real(WP), allocatable :: depth_nz(:), depth_nz1(:)   ! cached vertical coords (-zbar / -Z)
        ! vector output frame (Task 2.5): .true. => r2g-rotate vector pairs to GEOGRAPHIC before write
        ! (= FESOM2 vec_autorotate=.true., the production default); .false. => write NATIVE rotated-mesh
        ! components. FESOM3 default = geographic (the scientifically-useful frame).
        logical               :: autorotate = .true.
    end type t_io_means

contains

    ! ----------------------------------------------------------------- namelist.io (Task 2.6)

    ! Read namelist.io: &nml_general -> cfg (global knobs), &nml_list io_list -> list(1:nlist) rows.
    ! Mirrors FESOM2 io_meandata.F90:213-241 (open + read nml_general + read nml_list + 'unknown'
    ! sentinel scan), with the FESOM3 schema (chunk/compressor/filesplit/vec_frame + the 5th op field).
    ! ok=.false. when the file is absent (caller falls back to defaults/env). All ranks call it.
    subroutine means_read_namelist(path, cfg, list, nlist, ok)
        character(len=*), intent(in)  :: path
        type(t_io_config), intent(out):: cfg
        type(t_io_entry), intent(out) :: list(:)
        integer,          intent(out) :: nlist
        logical,          intent(out) :: ok
        integer            :: n_writers, chunk_time, chunk_vert, chunk_horiz, io_listsize
        character(len=16)  :: compressor, filesplit_freq, vec_frame
        type(t_io_entry)   :: io_list(MEANS_MAXF)
        integer            :: u, ios, i
        namelist /nml_general/ n_writers, chunk_time, chunk_vert, chunk_horiz, &
                               compressor, filesplit_freq, vec_frame, io_listsize
        namelist /nml_list/ io_list
        ! defaults (match t_io_config)
        n_writers = 0; chunk_time = 1; chunk_vert = 0; chunk_horiz = 500000
        compressor = 'none'; filesplit_freq = 'y'; vec_frame = 'geographic'; io_listsize = 0
        ok = .false.; nlist = 0
        open(newunit=u, file=trim(path), status='old', action='read', &
             form='formatted', access='sequential', iostat=ios)
        if (ios /= 0) return
        read(u, nml=nml_general, iostat=ios)   ! tolerate a missing/empty &nml_general (keep defaults)
        if (ios /= 0) rewind(u)
        read(u, nml=nml_list, iostat=ios)
        call zarr_check(ios == 0, 'means_read_namelist: cannot parse &nml_list in '//trim(path))
        close(u)
        cfg%n_writers = n_writers; cfg%chunk_time = max(1, chunk_time)
        cfg%chunk_vert = chunk_vert; cfg%chunk_horiz = chunk_horiz
        cfg%compressor = compressor; cfg%filesplit_freq = filesplit_freq(1:1); cfg%vec_frame = vec_frame
        do i = 1, min(size(io_list), size(list))
            if (trim(io_list(i)%id) == 'unknown') exit
            list(i) = io_list(i); nlist = i
        end do
        ok = .true.
    end subroutine means_read_namelist

    ! ----------------------------------------------------------------- setup / registration

    subroutine means_init(io, out_dir, mesh, partit, chunk_horiz, n_writers, calendar, vec_frame, &
                          chunk_time, chunk_vert, compressor, filesplit_freq)
        type(t_io_means), intent(out)          :: io
        character(len=*), intent(in)           :: out_dir
        type(t_mesh),     intent(in)           :: mesh
        type(t_partit),   intent(in), optional :: partit
        integer,          intent(in), optional :: chunk_horiz, n_writers, chunk_time, chunk_vert
        character(len=*), intent(in), optional :: calendar, vec_frame, compressor, filesplit_freq
        integer :: C, nw, i, nNodL, nEdgeO, nEdgeL, nElemO, nElemL, nElemF
        character(len=16) :: frame, cbuf
        io%out_dir  = out_dir
        io%nf       = 0
        io%calendar = 'standard'
        if (present(calendar)) io%calendar = calendar
        ! vector output frame: default geographic (autorotate); arg overrides default; env overrides arg.
        io%autorotate = .true.
        if (present(vec_frame)) io%autorotate = (trim(vec_frame) /= 'native')
        frame = ''
        call get_environment_variable('FESOM3_VEC_FRAME', frame)
        if (len_trim(frame) > 0) io%autorotate = (trim(frame) /= 'native')
        io%mr   = is_multirank(partit)
        io%mype = 0; io%comm = MPI_COMM_SELF; io%npes = 1
        if (present(partit)) then
            io%mype = partit%mype; io%comm = partit%MPI_COMM_FESOM; io%npes = partit%npes
        end if
        ! global writer knobs: arg (from cfg) overrides default; FESOM3_* env overrides arg.
        C = 500000; if (present(chunk_horiz)) C = chunk_horiz
        nw = 0;     if (present(n_writers))   nw = n_writers
        io%chunk_time = 1; if (present(chunk_time)) io%chunk_time = max(1, chunk_time)
        io%chunk_vert = 0; if (present(chunk_vert)) io%chunk_vert = chunk_vert
        io%compressor = 'none'; if (present(compressor)) io%compressor = compressor
        io%filesplit  = 'y'
        if (present(filesplit_freq)) then
            if (len_trim(filesplit_freq) > 0) io%filesplit = filesplit_freq(1:1)
        end if
        call read_env_int('FESOM3_CHUNK_HORIZ', C)
        call read_env_int('FESOM3_N_WRITERS',  nw)
        call read_env_int('FESOM3_CHUNK_TIME', io%chunk_time);  io%chunk_time = max(1, io%chunk_time)
        call read_env_int('FESOM3_CHUNK_VERT', io%chunk_vert)
        cbuf = ''; call get_environment_variable('FESOM3_COMPRESSOR', cbuf)
        if (len_trim(cbuf) > 0) io%compressor = cbuf
        cbuf = ''; call get_environment_variable('FESOM3_FILESPLIT', cbuf)
        if (len_trim(cbuf) > 0) io%filesplit = cbuf(1:1)
        call local_dims(mesh, partit, io%nNodO, nNodL, nEdgeO, nEdgeL, io%nElemO, nElemL, nElemF)
        call decomp_init_entity(io%Dn, C, nw, DECOMP_NODE, mesh, partit)
        call decomp_init_entity(io%De, C, nw, DECOMP_ELEM, mesh, partit)   ! Task 2.7 element decomp
        io%nl = mesh%nl
        ! cache owned node lon/lat (deg) + nlevels (3-D mask) — embedded into each per-period store.
        ! Also cache the ROTATED node coords (radians) that vector_r2g needs at write time (flag_coord=0,
        ! exactly as io_meandata.F90:io_r2g passes mesh%coord_nod2D for node-based vectors).
        allocate(io%lon_owned(max(1,io%nNodO)),  io%lat_owned(max(1,io%nNodO)), &
                 io%rlon_owned(max(1,io%nNodO)), io%rlat_owned(max(1,io%nNodO)), &
                 io%nlev_owned(max(1,io%nNodO)))
        call io_coords_compute(DECOMP_NODE, mesh, io%nNodO, io%lon_owned, io%lat_owned, &
                               io%rlon_owned, io%rlat_owned, io%nlev_owned)
        ! Task 2.7: owned ELEMENT-centroid coords. ROTATED centroid = elem_center (area-weighted
        ! centroid of the sub-triangles (1,3,4)+(1,2,3), cyclic-wrapped; hybrid-mesh plan Task 14);
        ! embedded GEOGRAPHIC centroid = its r2g image (deg). Element nlevels drives the 3-D mask.
        allocate(io%elon_owned(max(1,io%nElemO)),  io%elat_owned(max(1,io%nElemO)), &
                 io%relon_owned(max(1,io%nElemO)), io%relat_owned(max(1,io%nElemO)), &
                 io%enlev_owned(max(1,io%nElemO)))
        call io_coords_compute(DECOMP_ELEM, mesh, io%nElemO, io%elon_owned, io%elat_owned, &
                               io%relon_owned, io%relat_owned, io%enlev_owned)
        ! vertical coords (CF positive-down depths): nz = -zbar(1:nl), nz1 = -Z(1:nl-1)
        allocate(io%depth_nz(io%nl), io%depth_nz1(io%nl-1))
        do i = 1, io%nl;   io%depth_nz(i)  = real(-mesh%zbar(i), WP); end do
        do i = 1, io%nl-1; io%depth_nz1(i) = real(-mesh%Z(i),    WP); end do
    end subroutine means_init

    ! Shared field registration core (node OR element, 2-D OR 3-D). on_full_levels (3-D only):
    ! .true. => nl levels (vdim 'nz', valid 1..nlevels); .false. => nl-1 layers (vdim 'nz1', valid
    ! 1..nlevels-1). Below-bottom => NC_FILL (CF mask). Accumulator is sized to the OWNED entity count.
    subroutine add_field(io, name, long_name, units, ndim, is_elem, on_full_levels, std, precision, &
                         mean, freq, unit)
        type(t_io_means), intent(inout)        :: io
        character(len=*), intent(in)           :: name, long_name, units
        integer,          intent(in)           :: ndim
        logical,          intent(in)           :: is_elem, on_full_levels
        character(len=*), intent(in), optional :: std, precision, unit
        logical,          intent(in), optional :: mean
        integer,          intent(in), optional :: freq
        integer :: k, n
        call zarr_check(io%nf < MEANS_MAXF, 'means: too many fields ('//trim(name)//')')
        io%nf = io%nf + 1; k = io%nf
        io%f(k)%name = name; io%f(k)%long_name = long_name; io%f(k)%units = units
        io%f(k)%std = ''; if (present(std)) io%f(k)%std = std
        io%f(k)%dtype = field_dtype(precision)
        io%f(k)%is_mean = .false.; if (present(mean)) io%f(k)%is_mean = mean
        call set_cadence(io%f(k), freq, unit)
        io%f(k)%open = .false.; io%f(k)%count = 0
        io%f(k)%ndim = ndim
        io%f(k)%is_elem = is_elem
        io%f(k)%hdim = 'nod2'; if (is_elem) io%f(k)%hdim = 'elem'
        n = max(1, entity_owned(io, is_elem))
        if (ndim == 3) then
            if (on_full_levels) then
                io%f(k)%nlev = io%nl;   io%f(k)%voff = 0; io%f(k)%vdim = 'nz'
            else
                io%f(k)%nlev = io%nl-1; io%f(k)%voff = 1; io%f(k)%vdim = 'nz1'
            end if
            if (io%f(k)%dtype == '<f8') then
                allocate(io%f(k)%acc8_3d(io%f(k)%nlev, n)); io%f(k)%acc8_3d = 0.0_real64
            else
                allocate(io%f(k)%acc4_3d(io%f(k)%nlev, n)); io%f(k)%acc4_3d = 0.0_real32
            end if
        else
            if (io%f(k)%dtype == '<f8') then
                allocate(io%f(k)%acc8(n)); io%f(k)%acc8 = 0.0_real64
            else
                allocate(io%f(k)%acc4(n)); io%f(k)%acc4 = 0.0_real32
            end if
        end if
    end subroutine add_field

    subroutine means_define_node2d(io, name, long_name, units, std, precision, mean, freq, unit)
        type(t_io_means), intent(inout)        :: io
        character(len=*), intent(in)           :: name, long_name, units
        character(len=*), intent(in), optional :: std, precision, unit
        logical,          intent(in), optional :: mean
        integer,          intent(in), optional :: freq
        call add_field(io, name, long_name, units, 2, .false., .false., std, precision, mean, freq, unit)
    end subroutine means_define_node2d

    ! Register a 3-D node field (T/S on layers nz1, w on levels nz). See add_field for on_full_levels.
    subroutine means_define_node3d(io, name, long_name, units, on_full_levels, std, precision, mean, &
                                   freq, unit)
        type(t_io_means), intent(inout)        :: io
        character(len=*), intent(in)           :: name, long_name, units
        logical,          intent(in)           :: on_full_levels
        character(len=*), intent(in), optional :: std, precision, unit
        logical,          intent(in), optional :: mean
        integer,          intent(in), optional :: freq
        call add_field(io, name, long_name, units, 3, .false., on_full_levels, std, precision, mean, freq, unit)
    end subroutine means_define_node3d

    ! Task 2.7: register a 2-D ELEMENT field (data on owned elements; embed elem-centroid lon/lat).
    subroutine means_define_elem2d(io, name, long_name, units, std, precision, mean, freq, unit)
        type(t_io_means), intent(inout)        :: io
        character(len=*), intent(in)           :: name, long_name, units
        character(len=*), intent(in), optional :: std, precision, unit
        logical,          intent(in), optional :: mean
        integer,          intent(in), optional :: freq
        call add_field(io, name, long_name, units, 2, .true., .false., std, precision, mean, freq, unit)
    end subroutine means_define_elem2d

    ! Task 2.7: register a 3-D ELEMENT field (e.g. Av on levels nz, u/v on layers nz1).
    subroutine means_define_elem3d(io, name, long_name, units, on_full_levels, std, precision, mean, &
                                   freq, unit)
        type(t_io_means), intent(inout)        :: io
        character(len=*), intent(in)           :: name, long_name, units
        logical,          intent(in)           :: on_full_levels
        character(len=*), intent(in), optional :: std, precision, unit
        logical,          intent(in), optional :: mean
        integer,          intent(in), optional :: freq
        call add_field(io, name, long_name, units, 3, .true., on_full_levels, std, precision, mean, freq, unit)
    end subroutine means_define_elem3d

    ! Register a 3-D node VECTOR pair (Task 2.5): two node3d fields (x=zonal, y=meridional) linked so
    ! means_output r2g-rotates them TOGETHER to geographic (when autorotate) — the io_meandata io_r2g
    ! analog for unod/vnod. Each component is accumulated independently (means_accumulate as usual); the
    ! pairing matters only at write time. on_full_levels follows means_define_node3d (false => nl-1/nz1).
    ! Both components share freq/unit/precision (written together).
    subroutine means_define_vector3d(io, name_x, name_y, long_x, long_y, units, on_full_levels, &
                                     std_x, std_y, precision, mean, freq, unit)
        type(t_io_means), intent(inout)        :: io
        character(len=*), intent(in)           :: name_x, name_y, long_x, long_y, units
        logical,          intent(in)           :: on_full_levels
        character(len=*), intent(in), optional :: std_x, std_y, precision, unit
        logical,          intent(in), optional :: mean
        integer,          intent(in), optional :: freq
        integer :: kx, ky
        call means_define_node3d(io, name_x, long_x, units, on_full_levels, &
                                 std=std_x, precision=precision, mean=mean, freq=freq, unit=unit)
        kx = io%nf
        call means_define_node3d(io, name_y, long_y, units, on_full_levels, &
                                 std=std_y, precision=precision, mean=mean, freq=freq, unit=unit)
        ky = io%nf
        call link_vec_pair(io, kx, ky)
    end subroutine means_define_vector3d

    ! Task 2.7: a 3-D ELEMENT vector pair (e.g. u/v, bolus_u/bolus_v), r2g-rotated at the elem centroid.
    subroutine means_define_vector3d_elem(io, name_x, name_y, long_x, long_y, units, on_full_levels, &
                                          std_x, std_y, precision, mean, freq, unit)
        type(t_io_means), intent(inout)        :: io
        character(len=*), intent(in)           :: name_x, name_y, long_x, long_y, units
        logical,          intent(in)           :: on_full_levels
        character(len=*), intent(in), optional :: std_x, std_y, precision, unit
        logical,          intent(in), optional :: mean
        integer,          intent(in), optional :: freq
        integer :: kx, ky
        call means_define_elem3d(io, name_x, long_x, units, on_full_levels, &
                                 std=std_x, precision=precision, mean=mean, freq=freq, unit=unit)
        kx = io%nf
        call means_define_elem3d(io, name_y, long_y, units, on_full_levels, &
                                 std=std_y, precision=precision, mean=mean, freq=freq, unit=unit)
        ky = io%nf
        call link_vec_pair(io, kx, ky)
    end subroutine means_define_vector3d_elem

    subroutine link_vec_pair(io, kx, ky)
        type(t_io_means), intent(inout) :: io
        integer,          intent(in)    :: kx, ky
        io%f(kx)%is_vec_x = .true.; io%f(kx)%vec_partner = ky
        io%f(ky)%is_vec_y = .true.; io%f(ky)%vec_partner = kx
    end subroutine link_vec_pair

    ! True if a field named `name` is registered (lets the caller accumulate only listed fields).
    logical function means_has(io, name)
        type(t_io_means), intent(in) :: io
        character(len=*), intent(in) :: name
        means_has = (find_field(io, name) > 0)
    end function means_has

    ! ----------------------------------------------------------------- accumulate

    ! Accumulate a 2-D field's owned value for THIS step (in the output precision). mean => running sum
    ! + count; snapshot => overwrite (count=1). Call every step (FESOM2 update_means). nO follows the
    ! field's entity (node nNodO / element nElemO).
    subroutine means_accumulate_2d(io, name, owned)
        type(t_io_means), intent(inout) :: io
        character(len=*), intent(in)    :: name
        real(WP),         intent(in)    :: owned(:)
        integer :: k, i, nO
        k = find_field(io, name)
        call zarr_check(k > 0, 'means_accumulate: unknown field '//trim(name))
        call zarr_check(io%f(k)%ndim == 2, 'means_accumulate: 2-D call on 3-D field '//trim(name))
        nO = entity_owned(io, io%f(k)%is_elem)
        if (io%f(k)%is_mean) then
            if (io%f(k)%dtype == '<f8') then
                do i = 1, nO; io%f(k)%acc8(i) = io%f(k)%acc8(i) + real(owned(i), real64); end do
            else
                do i = 1, nO; io%f(k)%acc4(i) = io%f(k)%acc4(i) + real(owned(i), real32); end do
            end if
            io%f(k)%count = io%f(k)%count + 1
        else
            if (io%f(k)%dtype == '<f8') then
                do i = 1, nO; io%f(k)%acc8(i) = real(owned(i), real64); end do
            else
                do i = 1, nO; io%f(k)%acc4(i) = real(owned(i), real32); end do
            end if
            io%f(k)%count = 1
        end if
    end subroutine means_accumulate_2d

    ! Accumulate a 3-D field owned(1:nlev, 1:nO) for THIS step (output precision; nO = node/element).
    subroutine means_accumulate_3d(io, name, owned)
        type(t_io_means), intent(inout) :: io
        character(len=*), intent(in)    :: name
        real(WP),         intent(in)    :: owned(:,:)
        integer :: k, i, L, nlev, nO
        k = find_field(io, name)
        call zarr_check(k > 0, 'means_accumulate: unknown field '//trim(name))
        call zarr_check(io%f(k)%ndim == 3, 'means_accumulate: 3-D call on 2-D field '//trim(name))
        nlev = io%f(k)%nlev; nO = entity_owned(io, io%f(k)%is_elem)
        if (io%f(k)%is_mean) then
            if (io%f(k)%dtype == '<f8') then
                do i = 1, nO; do L = 1, nlev
                    io%f(k)%acc8_3d(L,i) = io%f(k)%acc8_3d(L,i) + real(owned(L,i), real64); end do; end do
            else
                do i = 1, nO; do L = 1, nlev
                    io%f(k)%acc4_3d(L,i) = io%f(k)%acc4_3d(L,i) + real(owned(L,i), real32); end do; end do
            end if
            io%f(k)%count = io%f(k)%count + 1
        else
            if (io%f(k)%dtype == '<f8') then
                do i = 1, nO; do L = 1, nlev
                    io%f(k)%acc8_3d(L,i) = real(owned(L,i), real64); end do; end do
            else
                do i = 1, nO; do L = 1, nlev
                    io%f(k)%acc4_3d(L,i) = real(owned(L,i), real32); end do; end do
            end if
            io%f(k)%count = 1
        end if
    end subroutine means_accumulate_3d

    ! ----------------------------------------------------------------- output (events + write)

    ! Per-step output driver: for each registered field, decide via its freq/unit event whether it is
    ! due THIS step; if so, (re)open its per-period store on a period change and write the in-flight
    ! record (divide by count, redistribute to canonical order, writers write chunks, rank 0 appends
    ! the time coord + bumps shape), then reset the accumulator. Collective: event_due + period are
    ! deterministic from clk (identical on every rank) so all ranks write the same fields in lockstep.
    subroutine means_output(io, istep, clk)
        type(t_io_means),   intent(inout) :: io
        integer,            intent(in)    :: istep
        type(t_means_clock), intent(in)   :: clk
        integer :: k, ierr, period
        real(real64) :: time_sec
        logical :: any_written
        time_sec = real(clk%daynew - 1, real64)*86400.0_real64 + clk%timenew
        period   = period_key(io, clk)
        any_written = .false.
        do k = 1, io%nf
            if (io%f(k)%is_vec_y) cycle                 ! written by its vec_x partner
            if (.not. event_due(io%f(k)%unit, io%f(k)%freq, istep, clk)) cycle
            call ensure_period_open(io, k, period, clk)
            call write_one_field(io, k, time_sec)
            any_written = .true.
        end do
        if (io%mr .and. any_written) call MPI_Barrier(io%comm, ierr)
    end subroutine means_output

    ! End of run. v1 keeps each <var>/.zarray current every record, so stores open with
    ! consolidated=False already; consolidated .zmetadata for field stores is a Task 2.6 polish.
    subroutine means_finalize(io)
        type(t_io_means), intent(inout) :: io
        integer :: ierr
        if (io%mr) call MPI_Barrier(io%comm, ierr)
    end subroutine means_finalize

    ! ----------------------------------------------------------------- write internals

    ! Dispatch one field's record to the right writer (scalar 2-D / scalar 3-D / vector pair), after
    ! resolving the entity context: D = node/element decomp, nO = owned count, nlevown = the per-entity
    ! valid-level count (3-D mask), rlon/rlat = ROTATED entity coords for the vector r2g.
    subroutine write_one_field(io, k, time_sec)
        type(t_io_means), intent(inout), target :: io
        integer,          intent(in)            :: k
        real(real64),     intent(in)            :: time_sec
        type(t_io_decomp), pointer :: D
        integer                    :: nO
        integer,  pointer          :: nlevown(:)
        real(WP), pointer          :: rlon(:), rlat(:)
        if (io%f(k)%is_elem) then
            D => io%De; nO = io%nElemO; nlevown => io%enlev_owned
            rlon => io%relon_owned; rlat => io%relat_owned
        else
            D => io%Dn; nO = io%nNodO; nlevown => io%nlev_owned
            rlon => io%rlon_owned; rlat => io%rlat_owned
        end if
        if (io%f(k)%is_vec_x) then
            call write_vec_pair(io, k, io%f(k)%vec_partner, time_sec, D, nO, nlevown, rlon, rlat)
        else if (io%f(k)%ndim == 2) then
            call write_scalar_2d(io, k, time_sec, D, nO)
        else
            call write_scalar_3d(io, k, time_sec, D, nO, nlevown)
        end if
    end subroutine write_one_field

    subroutine write_scalar_2d(io, k, time_sec, D, nO)
        type(t_io_means),   intent(inout) :: io
        integer,            intent(in)    :: k, nO
        real(real64),       intent(in)    :: time_sec
        type(t_io_decomp),  intent(in)    :: D
        real(WP), allocatable :: owned_mean(:), buf(:)
        integer :: i, cnt
        allocate(owned_mean(max(1,nO)), buf(max(1, D%w_nbuf)))
        cnt = max(1, io%f(k)%count)
        if (io%f(k)%dtype == '<f8') then
            do i = 1, nO; owned_mean(i) = io%f(k)%acc8(i) / real(cnt, real64); end do
        else
            do i = 1, nO; owned_mean(i) = real(io%f(k)%acc4(i) / real(cnt, real32), WP); end do
        end if
        call decomp_redistribute(D, owned_mean(1:nO), buf, 0.0_WP)   ! collective
        call emit_chunks_2d(io, D, io%f(k)%store, io%f(k)%a_data, buf, io%f(k)%t)
        if (io%mype == 0) call finish_record(io, k, time_sec)
        call reset_field(io%f(k))
        io%f(k)%t = io%f(k)%t + 1
    end subroutine write_scalar_2d

    ! 3-D: divide (output precision), mask below-bottom (L > nlevels[_nod2D] - voff) to NC_FILL,
    ! redistribute level-by-level, writers write the (chunk_time, chunk_vert) chunks.
    subroutine write_scalar_3d(io, k, time_sec, D, nO, nlevown)
        type(t_io_means),   intent(inout) :: io
        integer,            intent(in)    :: k, nO, nlevown(:)
        real(real64),       intent(in)    :: time_sec
        type(t_io_decomp),  intent(in)    :: D
        real(WP), allocatable :: mean3(:,:), buf3(:,:)
        integer :: i, L, cnt, nlev, nvalid
        nlev = io%f(k)%nlev
        cnt  = max(1, io%f(k)%count)
        allocate(mean3(nlev, max(1,nO)), buf3(nlev, max(1, D%w_nbuf)))
        do i = 1, nO
            nvalid = nlevown(i) - io%f(k)%voff
            do L = 1, nlev
                if (L <= nvalid) then
                    if (io%f(k)%dtype == '<f8') then
                        mean3(L,i) = io%f(k)%acc8_3d(L,i) / real(cnt, real64)
                    else
                        mean3(L,i) = real(io%f(k)%acc4_3d(L,i) / real(cnt, real32), WP)
                    end if
                else
                    mean3(L,i) = real(NC_FILL_DOUBLE, WP)
                end if
            end do
        end do
        call decomp_redistribute(D, mean3(1:nlev,1:nO), buf3, real(NC_FILL_DOUBLE, WP))
        call emit_chunks_3d(io, D, io%f(k)%store, io%f(k)%a_data, buf3, nlev, io%f(k)%t)
        if (io%mype == 0) call finish_record(io, k, time_sec)
        call reset_field(io%f(k))
        io%f(k)%t = io%f(k)%t + 1
    end subroutine write_scalar_3d

    ! 3-D vector pair (kx=zonal, ky=meridional), node OR element. FESOM2 ORDER (io_meandata.F90): io_r2g
    ! rotates the ACCUMULATED SUM in place (L2265) BEFORE compute_means divides by addcounter (L2335) —
    ! and rotate-then-divide /= divide-then-rotate in floating point, so we match it: rotate the SUM,
    ! THEN divide. <f8: rotate the r8 sum, divide in r8. <f4: promote the r4 sum to WP(=r8), rotate,
    ! demote to r4 (io_r2g r4 branch L3028-3033), divide in r4. Below-bottom (nlevels mask) -> NC_FILL,
    ! not rotated. native frame (autorotate=.false.) skips rotation. rlon/rlat = the ROTATED entity
    ! coords (node coord_nod2D, or the element centroid elem_center — io_r2g flag_coord=0).
    subroutine write_vec_pair(io, kx, ky, time_sec, D, nO, nlevown, rlon, rlat)
        type(t_io_means),   intent(inout) :: io
        integer,            intent(in)    :: kx, ky, nO, nlevown(:)
        real(real64),       intent(in)    :: time_sec
        type(t_io_decomp),  intent(in)    :: D
        real(WP),           intent(in)    :: rlon(:), rlat(:)
        real(WP), allocatable :: mx(:,:), my(:,:), bx(:,:), by(:,:)
        integer :: i, L, cnt, nlev, nvalid
        real(real64) :: sx8, sy8
        real(real32) :: sx4, sy4
        real(WP)     :: tx, ty
        nlev = io%f(kx)%nlev
        cnt  = max(1, io%f(kx)%count)
        allocate(mx(nlev, max(1,nO)), my(nlev, max(1,nO)), &
                 bx(nlev, max(1, D%w_nbuf)), by(nlev, max(1, D%w_nbuf)))
        do i = 1, nO
            nvalid = nlevown(i) - io%f(kx)%voff
            do L = 1, nlev
                if (L > nvalid) then
                    mx(L,i) = real(NC_FILL_DOUBLE, WP); my(L,i) = real(NC_FILL_DOUBLE, WP)
                    cycle
                end if
                if (io%f(kx)%dtype == '<f8') then
                    sx8 = io%f(kx)%acc8_3d(L,i); sy8 = io%f(ky)%acc8_3d(L,i)
                    if (io%autorotate) then
                        tx = real(sx8, WP); ty = real(sy8, WP)
                        call vector_r2g(tx, ty, rlon(i), rlat(i), 0)
                        sx8 = real(tx, real64); sy8 = real(ty, real64)
                    end if
                    mx(L,i) = real(sx8 / real(cnt, real64), WP)
                    my(L,i) = real(sy8 / real(cnt, real64), WP)
                else
                    sx4 = io%f(kx)%acc4_3d(L,i); sy4 = io%f(ky)%acc4_3d(L,i)
                    if (io%autorotate) then
                        tx = real(sx4, WP); ty = real(sy4, WP)
                        call vector_r2g(tx, ty, rlon(i), rlat(i), 0)
                        sx4 = real(tx, real32); sy4 = real(ty, real32)
                    end if
                    mx(L,i) = real(sx4 / real(cnt, real32), WP)
                    my(L,i) = real(sy4 / real(cnt, real32), WP)
                end if
            end do
        end do
        call decomp_redistribute(D, mx(1:nlev,1:nO), bx, real(NC_FILL_DOUBLE, WP))
        call decomp_redistribute(D, my(1:nlev,1:nO), by, real(NC_FILL_DOUBLE, WP))
        call emit_chunks_3d(io, D, io%f(kx)%store, io%f(kx)%a_data, bx, nlev, io%f(kx)%t)
        call emit_chunks_3d(io, D, io%f(ky)%store, io%f(ky)%a_data, by, nlev, io%f(ky)%t)
        if (io%mype == 0) then
            call finish_record(io, kx, time_sec)
            call finish_record(io, ky, time_sec)
        end if
        call reset_field(io%f(kx)); call reset_field(io%f(ky))
        io%f(kx)%t = io%f(kx)%t + 1; io%f(ky)%t = io%f(ky)%t + 1
    end subroutine write_vec_pair

    ! rank-0: append the time-coord value for record t (chunk_time-aware RMW) and bump the data/time
    ! .zarray shape[0] to t+1.
    subroutine finish_record(io, k, time_sec)
        type(t_io_means), intent(inout) :: io
        integer,          intent(in)    :: k
        real(real64),     intent(in)    :: time_sec
        real(WP) :: tbuf(io%chunk_time)
        integer  :: ct, tc, slot
        ct = io%chunk_time; tc = io%f(k)%t / ct; slot = mod(io%f(k)%t, ct)
        if (ct > 1 .and. slot > 0) then
            call zarr_read_chunk(io%f(k)%store, io%f(k)%a_time, [tc], tbuf)
        else
            tbuf = 0.0_WP
        end if
        tbuf(slot+1) = real(time_sec, WP)
        call zarr_write_chunk(io%f(k)%store, io%f(k)%a_time, [tc], tbuf)
        io%f(k)%a_data%dims(1) = io%f(k)%t + 1
        io%f(k)%a_time%dims(1) = io%f(k)%t + 1
        call zarr_rewrite_zarray(io%f(k)%store, io%f(k)%a_data)
        call zarr_rewrite_zarray(io%f(k)%store, io%f(k)%a_time)
    end subroutine finish_record

    ! Writers write their 2-D (time,<hdim>) data chunks for record t. chunk_time=1 => fresh chunk
    ! [t,c]; chunk_time>1 => read-modify-write the time-chunk [t/ct, c] at slot mod(t,ct). D = the
    ! field's entity decomp (node Dn / element De).
    subroutine emit_chunks_2d(io, D, store, arr, buf, t)
        type(t_io_means),   intent(in) :: io
        type(t_io_decomp),  intent(in) :: D
        type(t_zarr_store), intent(in) :: store
        type(t_zarr_array), intent(in) :: arr
        real(WP),           intent(in) :: buf(:)
        integer,            intent(in) :: t
        real(WP), allocatable :: chunk2(:,:)
        integer :: c, lo, ct, tc, slot, C0
        C0 = D%C; ct = io%chunk_time; tc = t/ct; slot = mod(t, ct)
        allocate(chunk2(ct, C0))
        do c = D%w_first_chunk, D%w_last_chunk
            lo = (c - D%w_first_chunk)*C0 + 1
            if (ct > 1 .and. slot > 0) then
                call zarr_read_chunk(store, arr, [tc, c], chunk2)
            else
                chunk2 = real(arr%fill, WP)
            end if
            chunk2(slot+1, 1:C0) = buf(lo:lo + C0 - 1)
            call zarr_write_chunk(store, arr, [tc, c], chunk2)
        end do
    end subroutine emit_chunks_2d

    ! Writers write their 3-D (time,nz,<hdim>) data chunks for record t, splitting nz into ceil(nlev/cv)
    ! vertical chunks (cv = effective chunk_vert). chunk_time RMW as in emit_chunks_2d.
    subroutine emit_chunks_3d(io, D, store, arr, buf3, nlev, t)
        type(t_io_means),   intent(in) :: io
        type(t_io_decomp),  intent(in) :: D
        type(t_zarr_store), intent(in) :: store
        type(t_zarr_array), intent(in) :: arr
        real(WP),           intent(in) :: buf3(:,:)
        integer,            intent(in) :: nlev, t
        real(WP), allocatable :: chunk3(:,:,:)
        integer :: c, lo, ct, tc, slot, C0, cv, nvc, vc, L0, cvn
        C0 = D%C; ct = io%chunk_time; tc = t/ct; slot = mod(t, ct)
        cv = vchunk_eff(io%chunk_vert, nlev); nvc = (nlev + cv - 1)/cv
        allocate(chunk3(ct, cv, C0))
        do c = D%w_first_chunk, D%w_last_chunk
            lo = (c - D%w_first_chunk)*C0 + 1
            do vc = 0, nvc - 1
                L0 = vc*cv; cvn = min(cv, nlev - L0)
                if (ct > 1 .and. slot > 0) then
                    call zarr_read_chunk(store, arr, [tc, vc, c], chunk3)
                else
                    chunk3 = real(arr%fill, WP)
                end if
                chunk3(slot+1, 1:cvn, 1:C0) = buf3(L0+1:L0+cvn, lo:lo + C0 - 1)
                if (cvn < cv) chunk3(slot+1, cvn+1:cv, 1:C0) = real(arr%fill, WP)
                call zarr_write_chunk(store, arr, [tc, vc, c], chunk3)
            end do
        end do
    end subroutine emit_chunks_3d

    subroutine reset_field(f)
        type(t_mean_field), intent(inout) :: f
        if (allocated(f%acc8))    f%acc8    = 0.0_real64
        if (allocated(f%acc4))    f%acc4    = 0.0_real32
        if (allocated(f%acc8_3d)) f%acc8_3d = 0.0_real64
        if (allocated(f%acc4_3d)) f%acc4_3d = 0.0_real32
        f%count = 0
    end subroutine reset_field

    ! ----------------------------------------------------------------- period stores / events

    ! Lazily (re)create field k's per-period store on a period change (collective: same clk on every
    ! rank). A vec_x also opens its partner so the pair writes into matching stores.
    subroutine ensure_period_open(io, k, period, clk)
        type(t_io_means),   intent(inout) :: io
        integer,            intent(in)    :: k, period
        type(t_means_clock), intent(in)   :: clk
        integer :: kp
        if (io%f(k)%cur_period == period) return
        call open_field_store(io, k, clk)
        io%f(k)%cur_period = period; io%f(k)%t = 0
        if (io%f(k)%is_vec_x) then
            kp = io%f(k)%vec_partner
            call open_field_store(io, kp, clk)
            io%f(kp)%cur_period = period; io%f(kp)%t = 0
        end if
    end subroutine ensure_period_open

    ! Create + define field k's store for the current period, then writers embed the static lon/lat
    ! (+ rank 0 the vertical coord). Store-create ordering: rank 0 defines -> barrier -> writers embed.
    ! D + own_lon/own_lat = the field's entity decomp + centroid coords (node or element).
    subroutine open_field_store(io, k, clk)
        type(t_io_means),   intent(inout), target :: io
        integer,            intent(in)    :: k
        type(t_means_clock), intent(in)   :: clk
        integer :: ierr
        character(len=4)              :: cy
        character(len=2)              :: cm
        character(len=:), allocatable :: tunits, suffix
        type(t_io_decomp), pointer    :: D
        real(WP),          pointer    :: own_lon(:), own_lat(:)
        if (io%f(k)%is_elem) then
            D => io%De; own_lon => io%elon_owned; own_lat => io%elat_owned
        else
            D => io%Dn; own_lon => io%lon_owned; own_lat => io%lat_owned
        end if
        write(cy, '(i4.4)') clk%year
        tunits = 'seconds since '//cy//'-01-01 00:00:00'     ! year-start ref (FESOM2 keeps yearnew)
        if (io%filesplit == 'm') then
            write(cm, '(i2.2)') clk%month
            suffix = cy//'_'//cm
        else
            suffix = cy
        end if
        io%f(k)%store%path = trim(io%out_dir)//'/'//trim(io%f(k)%name)//'.fesom.'//suffix//'.zarr'
        call def_field_store(io%f(k), D%N, D%C, io%mype, tunits, trim(io%calendar), &
                             io%chunk_time, vchunk_eff(io%chunk_vert, io%f(k)%nlev), trim(io%compressor))
        if (io%mr) call MPI_Barrier(io%comm, ierr)
        call io_coords_put(D, io%f(k)%store, io%f(k)%a_lon, own_lon)
        call io_coords_put(D, io%f(k)%store, io%f(k)%a_lat, own_lat)
        ! vertical coord is global/replicated -> rank 0 writes it whole (like mesh.diag nz/nz1)
        if (io%mype == 0 .and. io%f(k)%ndim == 3) then
            if (trim(io%f(k)%vdim) == 'nz') then
                call zarr_write_whole(io%f(k)%store, io%f(k)%a_vert, io%depth_nz(1:io%nl))
            else
                call zarr_write_whole(io%f(k)%store, io%f(k)%a_vert, io%depth_nz1(1:io%nl-1))
            end if
        end if
        io%f(k)%open = .true.
        if (io%mr) call MPI_Barrier(io%comm, ierr)
    end subroutine open_field_store

    ! the per-period store key: year ('y' filesplit) or year*100+month ('m').
    integer function period_key(io, clk)
        type(t_io_means),   intent(in) :: io
        type(t_means_clock), intent(in) :: clk
        if (io%filesplit == 'm') then
            period_key = clk%year*100 + clk%month
        else
            period_key = clk%year
        end if
    end function period_key

    ! Output event (FESOM2 gen_events.F90 ported verbatim): is a field with cadence (unit,freq) due
    ! at this clock instant? 'y' annual / 'm' monthly / 'd' daily / 'h' hourly / 's' per-step.
    logical function event_due(unit, freq, istep, clk)
        character(len=*),   intent(in) :: unit
        integer,            intent(in) :: freq, istep
        type(t_means_clock), intent(in) :: clk
        integer :: fr
        fr = max(1, freq)
        select case (unit(1:1))
        case ('y')
            event_due = (mod(clk%year - clk%yearstart + 1, fr) == 0 .and. &
                         clk%daynew == clk%ndpyr .and. clk%timenew == 86400.0_real64)
        case ('m')
            event_due = (mod(clk%month, fr) == 0 .and. &
                         clk%day_in_month == clk%ndim_month .and. clk%timenew == 86400.0_real64)
        case ('d')
            event_due = (mod(clk%daynew, fr) == 0 .and. clk%timenew == 86400.0_real64)
        case ('h')
            event_due = (mod(clk%timenew, 3600.0_real64*real(fr, real64)) == 0.0_real64)
        case ('s')
            event_due = (mod(istep, fr) == 0)
        case default
            call zarr_check(.false., 'means: unknown output freq unit '//unit)
            event_due = .false.
        end select
    end function event_due

    ! ----------------------------------------------------------------- store / coord definition

    ! Define one field store: data (time[,nz],nod2) growing from shape[0]=0 (time chunk = chunk_time,
    ! vert chunk = chunk_vert_eff, codec = compressor); time (time,); lon/lat (nod2,); nz/nz1 (3-D).
    ! Array handles init on all ranks; dirs/.zarray/.zattrs written on rank 0 (store-create ordering).
    subroutine def_field_store(f, N, C, mype, tunits, calendar, chunk_time, chunk_vert_eff, compressor)
        type(t_mean_field), intent(inout) :: f
        integer,            intent(in)    :: N, C, mype, chunk_time, chunk_vert_eff
        character(len=*),   intent(in)    :: tunits, calendar, compressor
        type(t_zarr_attrs) :: gat, at
        ! handles (all ranks). 2-D: fill_value:null (a valid 0 isn't NaN-masked). 3-D: NC_FILL so
        ! xarray masks below-bottom -> NaN (CF). time/lon/lat/vert coords: null fill, codec none.
        if (f%ndim == 2) then
            call zarr_array_init(f%a_data, trim(f%name), [0, N], [chunk_time, C], f%dtype, &
                                 has_fill=.false., codec=compressor)
        else
            call zarr_array_init(f%a_data, trim(f%name), [0, f%nlev, N], &
                                 [chunk_time, chunk_vert_eff, C], f%dtype, &
                                 fill=NC_FILL_DOUBLE, has_fill=.true., codec=compressor)
            call zarr_array_init(f%a_vert, trim(f%vdim), [f%nlev], [f%nlev], '<f8', has_fill=.false.)
        end if
        call zarr_array_init(f%a_time, 'time', [0], [chunk_time], '<f8', has_fill=.false.)
        call io_coords_init_lonlat(f%a_lon, f%a_lat, N, C)
        if (mype /= 0) return
        call zarr_attrs_init(gat)
        call zattr_str(gat, 'Conventions', 'CF-1.8')
        call zattr_str(gat, 'description', 'FESOM3 field output (Zarr v2, M9)')
        call zarr_create_store(f%store, trim(f%store%path), gat)
        ! data variable (horizontal dim = f%hdim: 'nod2' node / 'elem' element)
        call zarr_attrs_init(at)
        if (f%ndim == 2) then
            call zattr_str_arr(at, '_ARRAY_DIMENSIONS', [character(len=8) :: 'time', f%hdim])
        else
            call zattr_str_arr(at, '_ARRAY_DIMENSIONS', [character(len=8) :: 'time', f%vdim, f%hdim])
        end if
        call zattr_str(at, 'long_name', trim(f%long_name))
        if (len_trim(f%units) > 0) call zattr_str(at, 'units', trim(f%units))
        if (len_trim(f%std)   > 0) call zattr_str(at, 'standard_name', trim(f%std))
        if (f%is_mean) then
            call zattr_str(at, 'cell_methods', 'time: mean')
        else
            call zattr_str(at, 'cell_methods', 'time: point')
        end if
        call zattr_str(at, 'coordinates', 'lon lat')
        call zarr_define_array(f%store, f%a_data, at)
        ! vertical coordinate (3-D only)
        if (f%ndim == 3) then
            call zarr_attrs_init(at)
            call zattr_str_arr(at, '_ARRAY_DIMENSIONS', [character(len=8) :: f%vdim])
            call zattr_str(at, 'long_name', 'depth'); call zattr_str(at, 'units', 'meters')
            call zattr_str(at, 'positive', 'down')
            call zarr_define_array(f%store, f%a_vert, at)
        end if
        ! time coordinate (CF)
        call zarr_attrs_init(at)
        call zattr_str_arr(at, '_ARRAY_DIMENSIONS', [character(len=8) :: 'time'])
        call zattr_str(at, 'long_name', 'time')
        call zattr_str(at, 'units', tunits)
        call zattr_str(at, 'calendar', calendar)
        call zarr_define_array(f%store, f%a_time, at)
        ! lon / lat entity coords (node or element-centroid, dim = f%hdim) — shared embedding
        call io_coords_define_lonlat(f%store, f%a_lon, f%a_lat, f%hdim)
    end subroutine def_field_store

    ! ----------------------------------------------------------------- small helpers

    integer function find_field(io, name) result(k)
        type(t_io_means), intent(in) :: io
        character(len=*), intent(in) :: name
        integer :: i
        k = -1
        do i = 1, io%nf
            if (trim(io%f(i)%name) == trim(name)) then; k = i; return; end if
        end do
    end function find_field

    ! map a precision arg ('double'|'<f8'|'8' => <f8; else <f4) to a Zarr dtype. float32 default.
    function field_dtype(precision) result(dt)
        character(len=*), intent(in), optional :: precision
        character(len=8) :: dt
        dt = '<f4'
        if (present(precision)) then
            if (trim(precision) == 'double' .or. trim(precision) == '<f8' .or. trim(precision) == '8') &
                dt = '<f8'
        end if
    end function field_dtype

    subroutine set_cadence(f, freq, unit)
        type(t_mean_field), intent(inout)      :: f
        integer,            intent(in), optional :: freq
        character(len=*),   intent(in), optional :: unit
        f%freq = 1; f%unit = 's'
        if (present(freq)) f%freq = max(1, freq)
        if (present(unit)) then
            if (len_trim(unit) > 0) f%unit = unit(1:1)
        end if
    end subroutine set_cadence

    ! owned-entity count for a field: element nElemO (Task 2.7) or node nNodO.
    integer function entity_owned(io, is_elem)
        type(t_io_means), intent(in) :: io
        logical,          intent(in) :: is_elem
        if (is_elem) then; entity_owned = io%nElemO; else; entity_owned = io%nNodO; end if
    end function entity_owned

    ! effective vertical chunk: full depth when chunk_vert is 0 or >= nlev, else chunk_vert.
    integer function vchunk_eff(cv, nlev)
        integer, intent(in) :: cv, nlev
        if (cv <= 0 .or. cv >= nlev) then
            vchunk_eff = max(1, nlev)
        else
            vchunk_eff = cv
        end if
    end function vchunk_eff

    subroutine read_env_int(name, val)
        character(len=*), intent(in)    :: name
        integer,          intent(inout) :: val
        character(len=64) :: buf
        integer :: ios, tmp
        call get_environment_variable(name, buf, status=ios)
        if (ios == 0 .and. len_trim(buf) > 0) then
            read(buf, *, iostat=ios) tmp
            if (ios == 0) val = tmp
        end if
    end subroutine read_env_int

end module mod_io_means
