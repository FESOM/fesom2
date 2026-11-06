module mod_io_restart
    ! FESOM3 restart / checkpoint WRITE path (restart Stage 3, Task 3.2).
    !
    ! A checkpoint is an immutable folder  <RestartOutPath>/fesom.<YYYY>.<DDD>.<SSSSS>/  holding ONE
    ! single-variable, single-entity SNAPSHOT Zarr store per registered prognostic field, plus a tiny
    ! checkpoint.json (rank-0 provenance). Each <name>.zarr is an M9 field store in every respect
    ! (same lon/lat + _ARRAY_DIMENSIONS + UGRID embedding via mod_io_coords) EXCEPT it has NO time
    ! dimension and NO mean accumulation: the data var is (entity) for a 2-D field or (nlev, entity)
    ! for a 3-D field. So ushow / xarray render a checkpoint store verbatim, and history is just
    ! separate per-checkpoint folders (immutable, atomic, prunable) instead of records in a file.
    !
    ! Reuses the proven M9 stack verbatim — NOT mod_io_means (whose growing time-dim / per-period /
    ! cadence / accumulation are the wrong tool for an immutable snapshot):
    !   - mod_io_decomp : canonical-order MPI_Alltoallv redistribution -> writer subset (no rank-0
    !                     gather), so the store is partition-INDEPENDENT (np=2 write reads at np=8).
    !   - mod_io_coords : the shared lon/lat + _ARRAY_DIMENSIONS + CF/UGRID-attr embedding (Task 2.1).
    !   - mod_io_zarr   : the Zarr v2 chunk writer; store-create ordering = rank 0 creates store +
    !                     defines all arrays -> barrier -> writers write their chunks.
    !
    ! REGISTRATION, not binary dumps: a field descriptor holds a Fortran POINTER to the LIVE model
    ! array (eta_n, uv components, tracers, ice arrays, sigma), exactly like M9's register_output_var.
    ! The reader (Task 4.1) reads back INTO the same pointers. Pointer dummies (not copies) so a strided
    ! model section such as dyn%uv(1,:,:) associates without a copy and stays live for read-back.
    !
    ! Caller protocol (all ranks, same order):
    !   restart_init(R, mesh, partit, ...)                          once (decomps + cached coords)
    !   restart_register_field(R, name, units, entity, p2d=.. )     per field, once (associates a pointer)
    !   ... at a checkpoint instant ...
    !   restart_write(R, dir, year, day, time_sec, globalstep)      writes the folder + per-field stores
    !   restart_finalize(R)                                         end-of-run barrier
    !
    ! ATOMICITY (Task 3.3): restart_write stages a whole checkpoint into fesom.<tag>.tmp/ and only then
    ! atomically renames it into place and flips a one-line restart.latest pointer (write .tmp + rename),
    ! so a crash leaves either the OLD valid checkpoint or the NEW one — never a half-written folder and
    ! never a partial pointer. A reader (restart_resolve_latest) ALWAYS follows restart.latest and never
    ! scans the directory, so a stray crashed-write tmp folder or an unpointed finalized fesom folder is ignored.
    ! A keep-N prune (restart_keep) trims the oldest immutable checkpoints after each finalize.
    !
    ! SCOPE: Task 3.2 = the per-field writer + folder + checkpoint.json; Task 3.3 = atomic finalize +
    ! restart.latest + keep-N prune + restart_resolve_latest; Task 3.4 (this file now) =
    ! restart_register_state, which maps the FULL oce+ice prognostic state (incl. EVP sigma and the
    ! real(MP) mesh%hbar/hnode via a lossless WP staging copy) to the registry. The READ path (which
    ! calls restart_resolve_latest and copies WP->MP back for mp_src fields) is Stage 4 (Task 4.1).
    use mpi
    use, intrinsic :: iso_fortran_env, only: real64
    use mod_precision,   only: WP, MP
    use mod_mesh,        only: t_mesh
    use mod_partit,      only: t_partit
    use mod_dyn,         only: t_dyn        ! full-state registration (Task 3.4)
    use mod_tracer,      only: t_tracer
    use mod_ice,         only: t_ice
    use mod_part_bounds, only: is_multirank, local_dims
    use mod_io_zarr
    use mod_io_decomp
    use mod_io_coords    ! shared lon/lat + _ARRAY_DIMENSIONS + UGRID-attr embedding (Task 2.1)
    use mod_io_posix     ! POSIX fs shims: rename/fsync/rmtree/listdir for atomic finalize + prune (Task 3.1)
    use mod_halo,        only: exchange_nod, exchange_elem, exchange_elem_full   ! read-back halo recon (Task 4.1)
    implicit none
    private

    ! Re-export the entity selectors so a caller needs only `use mod_io_restart`.
    public :: DECOMP_NODE, DECOMP_ELEM
    public :: t_restart, t_restart_field, RESTART_MAXF
    public :: restart_init, restart_register_field, restart_register_field_mp, restart_register_state
    public :: restart_write_field, restart_write, restart_finalize
    public :: restart_resolve_latest
    ! READ path (Stage 4, Task 4.1) + the per-field halo-exchange variant selectors
    public :: restart_read, restart_read_field, restart_halo_exchange_all
    public :: restart_skip_field_on_read, restart_skip_tracers_on_read
    public :: RESTART_HALO_NONE, RESTART_HALO_NODE, RESTART_HALO_ELEM, RESTART_HALO_ELEM_FULL

    integer, parameter :: RESTART_MAXF = 64
    integer, parameter :: RESTART_FORMAT_VERSION = 1

    ! Halo-exchange variant a field's read-back uses to reconstruct its halo from owner values, matching
    ! the variant the field's in-step consumer uses (a wrong variant leaves stale halo cells):
    !   NONE       no exchange  — consumer reads owned-only (uv_rhsAB, EVP sigma); halo stays fresh-allocate
    !   NODE       exchange_nod — every node field (node halo owner-consistent at end of step)
    !   ELEM       exchange_elem      (eDim com_elem2D) — dyn%uv u/v (production update_vel; full CORRUPTS)
    !   ELEM_FULL  exchange_elem_full (eDim+eXDim)      — element fields whose consumer needs the full halo
    integer, parameter :: RESTART_HALO_NONE = 0, RESTART_HALO_NODE = 1, &
                          RESTART_HALO_ELEM = 2, RESTART_HALO_ELEM_FULL = 3

    ! One registered prognostic field: its identity + a live POINTER to the model array it serializes.
    ! ndim=2 uses p2d (entity,); ndim=3 uses p3d (nlev, entity). on_full_levels picks nl levels ('nz')
    ! vs nl-1 layers ('nz1'), the M9 distinction — a wrong level count breaks the shape/gate.
    type :: t_restart_field
        character(len=64) :: name   = ''
        character(len=32) :: units  = ''
        integer           :: entity = DECOMP_NODE      ! DECOMP_NODE | DECOMP_ELEM
        integer           :: ndim   = 2                ! 2 or 3
        integer           :: halo   = RESTART_HALO_NONE ! read-back halo-exchange variant (Task 4.1)
        logical           :: on_full_levels = .false.  ! 3-D: nl levels (true) vs nl-1 layers (false)
        character(len=8)  :: dtype  = '<f8'            ! restart default = full precision
        integer           :: nlev   = 1                ! vertical size (3-D): nl or nl-1
        character(len=8)  :: hdim   = 'nod2'           ! 'nod2' (node) / 'elem' (element)
        character(len=8)  :: vdim   = ''               ! 'nz' (levels) / 'nz1' (layers)
        real(WP), pointer :: p2d(:)   => null()        ! live array (2-D field, WP source)
        real(WP), pointer :: p3d(:,:) => null()        ! live array (3-D field, WP source)
        ! MP-source path (Task 3.4 precision gotcha): mesh%hbar / mesh%hnode are real(MP), which may
        ! differ from WP (the decomp/Zarr path is WP). We never pass an MP array to a WP dummy; instead
        ! we hold the live MP pointer + a flag, copy MP->WP into a local staging buffer at write time,
        ! and (Task 4.1) copy WP->MP back at read time. dtype is forced <f8 so WP(real64) >= MP is
        ! lossless. mp_src=.false. fields use p2d/p3d directly (no copy).
        logical           :: mp_src   = .false.
        real(MP), pointer :: pmp2d(:)   => null()      ! live array (2-D field, MP source)
        real(MP), pointer :: pmp3d(:,:) => null()      ! live array (3-D field, MP source)
        ! RECOM4 M10 (#57): registered and written, but not read by restart_read
        logical           :: skip_read = .false.
    end type t_restart_field

    type :: t_restart
        integer               :: nf = 0
        type(t_restart_field) :: f(RESTART_MAXF)
        type(t_io_decomp)     :: Dn                    ! node decomp (node fields)
        type(t_io_decomp)     :: De                    ! element decomp (element fields)
        integer               :: nNodO = 0, nElemO = 0
        integer               :: nl = 0
        integer               :: mype = 0, comm = MPI_COMM_SELF, npes = 1
        logical               :: mr = .false.
        ! global writer knobs (mirror mod_io_means)
        integer               :: chunk_vert = 0        ! vertical chunk (0 => full depth single chunk)
        character(len=16)     :: compressor = 'none'   ! data-array codec: 'none' | 'lz4'
        integer               :: restart_keep = 0      ! keep-N prune: keep newest N checkpoints (0 = keep all)
        ! cached owned coords for the embed (node + element-centroid), like means_init
        real(WP), allocatable :: lon_n(:), lat_n(:), rlon_n(:), rlat_n(:)
        integer,  allocatable :: nlev_n(:)
        real(WP), allocatable :: lon_e(:), lat_e(:), rlon_e(:), rlat_e(:)
        integer,  allocatable :: nlev_e(:)
        real(WP), allocatable :: depth_nz(:), depth_nz1(:)  ! -zbar / -Z (CF positive-down depths)
    end type t_restart

contains

    ! ----------------------------------------------------------------- setup / registration

    subroutine restart_init(R, mesh, partit, chunk_horiz, n_writers, chunk_vert, compressor, restart_keep)
        type(t_restart),  intent(out)          :: R
        type(t_mesh),     intent(in)           :: mesh
        type(t_partit),   intent(in), optional :: partit
        integer,          intent(in), optional :: chunk_horiz, n_writers, chunk_vert
        character(len=*), intent(in), optional :: compressor
        integer,          intent(in), optional :: restart_keep
        integer :: C, nw, i, nNodL, nEdgeO, nEdgeL, nElemL, nElemF
        character(len=16) :: cbuf
        R%nf = 0
        R%mr   = is_multirank(partit)
        R%mype = 0; R%comm = MPI_COMM_SELF; R%npes = 1
        if (present(partit)) then
            R%mype = partit%mype; R%comm = partit%MPI_COMM_FESOM; R%npes = partit%npes
        end if
        ! writer knobs: arg overrides default; FESOM3_* env overrides arg (mirror means_init).
        C  = 500000; if (present(chunk_horiz)) C  = chunk_horiz
        nw = 0;      if (present(n_writers))   nw = n_writers
        R%chunk_vert = 0; if (present(chunk_vert)) R%chunk_vert = chunk_vert
        R%compressor = 'none'; if (present(compressor)) R%compressor = compressor
        R%restart_keep = 0; if (present(restart_keep)) R%restart_keep = restart_keep
        call read_env_int('FESOM3_CHUNK_HORIZ', C)
        call read_env_int('FESOM3_N_WRITERS',  nw)
        call read_env_int('FESOM3_CHUNK_VERT', R%chunk_vert)
        call read_env_int('FESOM3_RESTART_KEEP', R%restart_keep)
        cbuf = ''; call get_environment_variable('FESOM3_COMPRESSOR', cbuf)
        if (len_trim(cbuf) > 0) R%compressor = cbuf
        call local_dims(mesh, partit, R%nNodO, nNodL, nEdgeO, nEdgeL, R%nElemO, nElemL, nElemF)
        call decomp_init_entity(R%Dn, C, nw, DECOMP_NODE, mesh, partit)
        call decomp_init_entity(R%De, C, nw, DECOMP_ELEM, mesh, partit)
        R%nl = mesh%nl
        ! cached owned coords (geographic deg lon/lat for the embed + rotated rad + nlevels mask).
        allocate(R%lon_n(max(1,R%nNodO)),  R%lat_n(max(1,R%nNodO)), &
                 R%rlon_n(max(1,R%nNodO)), R%rlat_n(max(1,R%nNodO)), R%nlev_n(max(1,R%nNodO)))
        call io_coords_compute(DECOMP_NODE, mesh, R%nNodO, R%lon_n, R%lat_n, R%rlon_n, R%rlat_n, R%nlev_n)
        allocate(R%lon_e(max(1,R%nElemO)),  R%lat_e(max(1,R%nElemO)), &
                 R%rlon_e(max(1,R%nElemO)), R%rlat_e(max(1,R%nElemO)), R%nlev_e(max(1,R%nElemO)))
        call io_coords_compute(DECOMP_ELEM, mesh, R%nElemO, R%lon_e, R%lat_e, R%rlon_e, R%rlat_e, R%nlev_e)
        ! vertical coords (CF positive-down depths): nz = -zbar(1:nl), nz1 = -Z(1:nl-1)
        allocate(R%depth_nz(R%nl), R%depth_nz1(R%nl-1))
        do i = 1, R%nl;   R%depth_nz(i)  = real(-mesh%zbar(i), WP); end do
        do i = 1, R%nl-1; R%depth_nz1(i) = real(-mesh%Z(i),    WP); end do
    end subroutine restart_init

    ! Register one prognostic field by associating a LIVE pointer to its model array. Pass p2d for a
    ! 2-D field, p3d for a 3-D field (on_full_levels selects nl vs nl-1). entity = DECOMP_NODE/ELEM.
    ! The pointer is stored, not copied, so the writer reads the live values and the reader (Task 4.1)
    ! can write back into the same array. Caller passes a Fortran pointer (e.g. peta => dyn%eta_n).
    subroutine restart_register_field(R, name, units, entity, p2d, p3d, on_full_levels, precision, halo)
        type(t_restart),  intent(inout)        :: R
        character(len=*), intent(in)           :: name, units
        integer,          intent(in)           :: entity
        real(WP), pointer, intent(in), optional :: p2d(:)
        real(WP), pointer, intent(in), optional :: p3d(:,:)
        logical,          intent(in), optional :: on_full_levels
        character(len=*), intent(in), optional :: precision
        integer,          intent(in), optional :: halo     ! override the entity-default read-back variant
        integer :: k
        k = restart_new_field(R, name, units, entity, precision)
        if (present(halo)) R%f(k)%halo = halo
        call restart_set_levels(R, k, present(p3d), on_full_levels)
        if (present(p3d)) then
            R%f(k)%p3d => p3d
        else if (present(p2d)) then
            R%f(k)%p2d => p2d
        end if
    end subroutine restart_register_field

    ! Register an MP-source field (mesh%hbar / mesh%hnode). Identical to restart_register_field but the
    ! live pointer is real(MP): we never bind it to a WP dummy. dtype is forced <f8 so the WP staging
    ! copy at write time (and the WP->MP copy at read time, Task 4.1) round-trips MP exactly when
    ! WP=real64 >= MP. Pass pmp2d for a 2-D field, pmp3d for a 3-D field (on_full_levels picks nl vs nl-1).
    subroutine restart_register_field_mp(R, name, units, entity, pmp2d, pmp3d, on_full_levels)
        type(t_restart),   intent(inout)        :: R
        character(len=*),  intent(in)           :: name, units
        integer,           intent(in)           :: entity
        real(MP), pointer, intent(in), optional :: pmp2d(:)
        real(MP), pointer, intent(in), optional :: pmp3d(:,:)
        logical,           intent(in), optional :: on_full_levels
        integer :: k
        k = restart_new_field(R, name, units, entity, '8')      ! MP -> always full precision
        call restart_set_levels(R, k, present(pmp3d), on_full_levels)
        R%f(k)%mp_src = .true.
        if (present(pmp3d)) then
            R%f(k)%pmp3d => pmp3d
        else if (present(pmp2d)) then
            R%f(k)%pmp2d => pmp2d
        end if
    end subroutine restart_register_field_mp

    ! Bump the registry and set the identity fields (name/units/entity/hdim/dtype) common to the WP and
    ! MP registration entry points; returns the new field index k.
    integer function restart_new_field(R, name, units, entity, precision) result(k)
        type(t_restart),  intent(inout)        :: R
        character(len=*), intent(in)           :: name, units
        integer,          intent(in)           :: entity
        character(len=*), intent(in), optional :: precision
        call zarr_check(R%nf < RESTART_MAXF, 'restart: too many fields ('//trim(name)//')')
        R%nf = R%nf + 1; k = R%nf
        R%f(k)%name   = name
        R%f(k)%units  = units
        R%f(k)%entity = entity
        R%f(k)%hdim   = 'nod2'; if (entity == DECOMP_ELEM) R%f(k)%hdim = 'elem'
        R%f(k)%dtype  = restart_dtype(precision)
        ! default read-back halo: node -> exchange_nod; element -> none (overridable, e.g. u/v -> ELEM)
        R%f(k)%halo   = RESTART_HALO_NONE
        if (entity == DECOMP_NODE) R%f(k)%halo = RESTART_HALO_NODE
    end function restart_new_field

    ! Resolve the level kind (2-D vs 3-D, nl vs nl-1) for field k from the pointer-presence + on_full_levels.
    subroutine restart_set_levels(R, k, is3d, on_full_levels)
        type(t_restart), intent(inout)        :: R
        integer,         intent(in)           :: k
        logical,         intent(in)           :: is3d
        logical,         intent(in), optional :: on_full_levels
        logical :: ofl
        if (is3d) then
            R%f(k)%ndim = 3
            ofl = .false.; if (present(on_full_levels)) ofl = on_full_levels
            R%f(k)%on_full_levels = ofl
            if (ofl) then
                R%f(k)%nlev = R%nl;   R%f(k)%vdim = 'nz'
            else
                R%f(k)%nlev = R%nl-1; R%f(k)%vdim = 'nz1'
            end if
        else
            R%f(k)%ndim = 2
            R%f(k)%nlev = 1
        end if
    end subroutine restart_set_levels

    ! ----------------------------------------------------------------- full prognostic state (Task 3.4)

    ! Register EVERY prognostic field of the full oce+ice state by associating a live pointer to each
    ! model array, mirroring FESOM2's ini_ocean_io / ini_ice_io field lists (the verified safe superset,
    ! F-F) and adding the EVP sigma stress (F-C, which FESOM2 omits). This is the one call Stage 5 makes
    ! with the real structs; the reader (Task 4.1) iterates the same registry and writes back into these
    ! pointers (mp_src fields via a WP->MP copy). Conditionals (AB_order, num_tracers, mix_scheme) are
    ! read from the structs / args so a different config registers a different — but self-consistent — set.
    !
    !   store names (== FESOM2 oracle, except eta_n which FESOM3 names per the plan field-set table):
    !     node 2-D : eta_n d_eta hbar(MP) ssh_rhs_old | area hice hsnow uice vice t_skin
    !     node 3-D : hnode(MP, nl-1) | <tr> <tr>_AB <tr>_M1 [<tr>_M2 if tracer AB_order==3]
    !                w w_expl w_impl (FULL nl) | tke (FULL nl, if mix_scheme==5)
    !     elem 3-D : u v urhs_AB vrhs_AB [urhs_AB3 vrhs_AB3 if dyn%AB_order==3]   (nl-1)
    !     elem 2-D : sigma11 sigma12 sigma22
    subroutine restart_register_state(R, dyn, tracers, ice, mesh, mix_scheme)
        type(t_restart), intent(inout)      :: R
        type(t_dyn),     intent(in), target :: dyn
        type(t_tracer),  intent(in), target :: tracers
        type(t_ice),     intent(in), target :: ice
        type(t_mesh),    intent(in), target :: mesh
        integer,         intent(in), optional :: mix_scheme
        real(WP), pointer :: p2(:), p3(:,:)
        real(MP), pointer :: q2(:), q3(:,:)
        integer           :: j, ms
        character(len=32) :: tn
        ms = 0; if (present(mix_scheme)) ms = mix_scheme

        ! ---- ocean NODE 2-D ----
        p2 => dyn%eta_n;       call restart_register_field(R, 'eta_n',       'm', DECOMP_NODE, p2d=p2)
        ! d_eta: SSH increment is the CG solve's INITIAL GUESS (oce_ssh_solve.F90:87 x0=dynamics%d_eta).
        ! The CG converges only to a relative tolerance (soltol=1e-5), NOT machine precision, so the
        ! converged d_eta depends on the carried x0. Without it the resumed SSH solve (substep 9) drifts at
        ! the ~1e-5 tolerance floor -> d_eta/eta_n diverge. Carried prognostic state. NODE halo variant.
        p2 => dyn%d_eta;       call restart_register_field(R, 'd_eta',       'm', DECOMP_NODE, p2d=p2)
        q2 => mesh%hbar;       call restart_register_field_mp(R, 'hbar',      'm', DECOMP_NODE, pmp2d=q2)
        p2 => dyn%ssh_rhs_old; call restart_register_field(R, 'ssh_rhs_old', '',  DECOMP_NODE, p2d=p2)

        ! ---- ocean NODE 3-D layers (nl-1): ALE layer thickness (MP) ----
        q3 => mesh%hnode;      call restart_register_field_mp(R, 'hnode',     'm', DECOMP_NODE, pmp3d=q3)

        ! ---- ocean ELEMENT 3-D layers (nl-1): velocity + Adams-Bashforth history ----
        ! u/v reconstruct over the eDim com_elem2D halo (production update_vel; full halo CORRUPTS the
        ! trajectory — oce_ale.F90:140); uv_rhsAB is owned-only (compute_vel_rhs) => RESTART_HALO_NONE.
        p3 => dyn%uv(1,:,:);         call restart_register_field(R, 'u',       'm/s', DECOMP_ELEM, p3d=p3, halo=RESTART_HALO_ELEM)
        p3 => dyn%uv(2,:,:);         call restart_register_field(R, 'v',       'm/s', DECOMP_ELEM, p3d=p3, halo=RESTART_HALO_ELEM)
        p3 => dyn%uv_rhsAB(1,1,:,:); call restart_register_field(R, 'urhs_AB', 'm/s', DECOMP_ELEM, p3d=p3)
        p3 => dyn%uv_rhsAB(1,2,:,:); call restart_register_field(R, 'vrhs_AB', 'm/s', DECOMP_ELEM, p3d=p3)
        if (dyn%AB_order == 3) then
            p3 => dyn%uv_rhsAB(2,1,:,:); call restart_register_field(R, 'urhs_AB3', 'm/s', DECOMP_ELEM, p3d=p3)
            p3 => dyn%uv_rhsAB(2,2,:,:); call restart_register_field(R, 'vrhs_AB3', 'm/s', DECOMP_ELEM, p3d=p3)
        end if

        ! ---- tracers (NODE 3-D layers nl-1): values + AB interp + valuesold history (M1 mandatory) ----
        do j = 1, tracers%num_tracers
            tn = tracer_name(tracers%data(j)%ID, j)
            p3 => tracers%data(j)%values;          call restart_register_field(R, trim(tn),         '', DECOMP_NODE, p3d=p3)
            p3 => tracers%data(j)%valuesAB;         call restart_register_field(R, trim(tn)//'_AB',  '', DECOMP_NODE, p3d=p3)
            p3 => tracers%data(j)%valuesold(1,:,:); call restart_register_field(R, trim(tn)//'_M1',  '', DECOMP_NODE, p3d=p3)
            if (tracers%data(j)%AB_order == 3) then
                p3 => tracers%data(j)%valuesold(2,:,:); call restart_register_field(R, trim(tn)//'_M2', '', DECOMP_NODE, p3d=p3)
            end if
        end do

        ! ---- vertical velocities (NODE 3-D, FULL levels nl) ----
        p3 => dyn%w;   call restart_register_field(R, 'w',      'm/s', DECOMP_NODE, p3d=p3, on_full_levels=.true.)
        p3 => dyn%w_e; call restart_register_field(R, 'w_expl', 'm/s', DECOMP_NODE, p3d=p3, on_full_levels=.true.)
        p3 => dyn%w_i; call restart_register_field(R, 'w_impl', 'm/s', DECOMP_NODE, p3d=p3, on_full_levels=.true.)

        ! ---- optional TKE (NODE 3-D, FULL levels nl) — only when TKE mixing is active (mix_scheme==5) ----
        if (ms == 5 .and. allocated(dyn%work%tke)) then
            p3 => dyn%work%tke; call restart_register_field(R, 'tke', 'm2/s2', DECOMP_NODE, p3d=p3, on_full_levels=.true.)
        end if

        ! ---- ice NODE 2-D: tracers (a_ice/m_ice/m_snow) + velocity ----
        p2 => ice%data(1)%values; call restart_register_field(R, 'area',  '',    DECOMP_NODE, p2d=p2)
        p2 => ice%data(2)%values; call restart_register_field(R, 'hice',  'm',   DECOMP_NODE, p2d=p2)
        p2 => ice%data(3)%values; call restart_register_field(R, 'hsnow', 'm',   DECOMP_NODE, p2d=p2)
        p2 => ice%uice;           call restart_register_field(R, 'uice',  'm/s', DECOMP_NODE, p2d=p2)
        p2 => ice%vice;           call restart_register_field(R, 'vice',  'm/s', DECOMP_NODE, p2d=p2)
        ! t_skin: thermo skin temperature is genuine carried prognostic state — the Newton ice-surface
        ! solver seeds from it each step (mod_ice_thermo.F90 read@222 t=t_skin(i), write@260 t_skin(i)=t).
        ! Without it the resumed flx_fw -> water_flux -> ssh_rhs diverges (substep 8). NODE halo variant.
        p2 => ice%thermo%t_skin;  call restart_register_field(R, 't_skin', 'degC', DECOMP_NODE, p2d=p2)

        ! ---- ice ELEMENT 2-D: EVP stress tensor (F-C: serialized so the gate is max|Δ|=0 on ice too) ----
        p2 => ice%work%sigma11; call restart_register_field(R, 'sigma11', '', DECOMP_ELEM, p2d=p2)
        p2 => ice%work%sigma12; call restart_register_field(R, 'sigma12', '', DECOMP_ELEM, p2d=p2)
        p2 => ice%work%sigma22; call restart_register_field(R, 'sigma22', '', DECOMP_ELEM, p2d=p2)
    end subroutine restart_register_state

    ! Tracer store name == FESOM2 ini_ocean_io CASE(id) (io_restart.F90:167-212): id keys the canonical
    ! name (1=temp, 2=salt, passive species by tracer ID), default 'tra'//j for an unlisted passive tracer.
    ! RECOM4 M10 (#57): mark a registered field as written but not read. Used for a
    ! biology cold start on a physics checkpoint (REcoM_restart = .false.): the field
    ! stays registered, so the run's own checkpoints carry it, but restart_read leaves
    ! its live array as the cold init set it. Stops if no registered field has the name.
    subroutine restart_skip_field_on_read(R, name)
        type(t_restart),  intent(inout) :: R
        character(len=*), intent(in)    :: name
        integer :: k
        logical :: found
        found = .false.
        do k = 1, R%nf
            if (trim(R%f(k)%name) == trim(name)) then
                R%f(k)%skip_read = .true.
                found = .true.
            end if
        end do
        call zarr_check(found, 'restart_skip_field_on_read: no registered field '//trim(name))
    end subroutine restart_skip_field_on_read

    ! The same for tracers jfirst .. num_tracers: the value and history stores
    ! (<tr>, <tr>_AB, <tr>_M1 and, for AB_order 3, <tr>_M2), named as
    ! restart_register_state names them.
    subroutine restart_skip_tracers_on_read(R, tracers, jfirst)
        type(t_restart), intent(inout) :: R
        type(t_tracer),  intent(in)    :: tracers
        integer,         intent(in)    :: jfirst
        integer           :: j
        character(len=32) :: tn
        do j = jfirst, tracers%num_tracers
            tn = tracer_name(tracers%data(j)%ID, j)
            call restart_skip_field_on_read(R, trim(tn))
            call restart_skip_field_on_read(R, trim(tn)//'_AB')
            call restart_skip_field_on_read(R, trim(tn)//'_M1')
            if (tracers%data(j)%AB_order == 3) call restart_skip_field_on_read(R, trim(tn)//'_M2')
        end do
    end subroutine restart_skip_tracers_on_read

    function tracer_name(id, j) result(nm)
        integer, intent(in) :: id, j
        character(len=32)   :: nm
        select case (id)
        case (1);   nm = 'temp'
        case (2);   nm = 'salt'
        case (6);   nm = 'sf6'
        case (11);  nm = 'cfc11'
        case (12);  nm = 'cfc12'
        case (14);  nm = 'r14c'
        case (39);  nm = 'r39ar'
        case (101); nm = 'h2o18'
        case (102); nm = 'hDo16'
        case (103); nm = 'h2o16'
        case default
            write(nm, '(a3,i4.4)') 'tra', j     ! FESOM2 oracle: write(trname,'(A3,i4.4)') 'tra_', j
        end select
    end function tracer_name

    ! ----------------------------------------------------------------- checkpoint write

    ! Write one full checkpoint ATOMICALLY (Task 3.3). The on-disk publish sequence guarantees that a
    ! crash at any point leaves either the previous valid checkpoint or this one — never a torn folder
    ! and never a partial restart.latest:
    !
    !   1. rank 0 clears any stale fesom.<tag>.tmp/ (from a crashed earlier write at THIS tag) and
    !      mkdir's a clean staging dir.                                   <- barrier ->
    !   2. EVERY writer writes its field stores into fesom.<tag>.tmp/, then rank 0 writes
    !      checkpoint.json there.                                         <- barrier ->
    !   3. rank 0 fsyncs the staging dir, then rename(2)s fesom.<tag>.tmp/ -> fesom.<tag>/ (atomic
    !      within one filesystem: the immutable checkpoint appears all-at-once), then fsyncs the parent.
    !   4. rank 0 atomically flips restart.latest to the new folder NAME (write restart.latest.tmp +
    !      rename), then runs the keep-N prune (best-effort; warns, never aborts).
    !
    ! folder/tag = fesom.<YYYY>.<DDD>.<SSSSS> (zero-padded; SSSSS = int(time_sec) = sec-of-day / timenew;
    ! lexical order == chronological). time_sec is the new-clock sec-of-day; Stage 5 passes timenew.
    subroutine restart_write(R, checkpoint_dir, year, day, time_sec, globalstep)
        type(t_restart),  intent(in), target :: R
        character(len=*), intent(in)         :: checkpoint_dir
        integer,          intent(in)         :: year, day
        real(real64),     intent(in)         :: time_sec
        integer,          intent(in), optional :: globalstep
        integer :: k, ierr, gstep, irc
        character(len=4) :: cy
        character(len=3) :: cd
        character(len=5) :: cs
        character(len=:), allocatable :: tag, fname, final_folder, tmp_folder
        gstep = 0; if (present(globalstep)) gstep = globalstep
        write(cy, '(i4.4)') year
        write(cd, '(i3.3)') day
        write(cs, '(i5.5)') int(time_sec)
        tag          = cy//'.'//cd//'.'//cs                ! YYYY.DDD.SSSSS (lexical == chronological)
        fname        = 'fesom.'//tag                       ! the FINAL folder NAME (what restart.latest holds)
        final_folder = trim(checkpoint_dir)//'/'//fname
        tmp_folder   = trim(checkpoint_dir)//'/'//fname//'.tmp'
        ! 1) rank 0: clear any stale tmp from a crashed write at this tag, then (re)create clean tmp.
        if (R%mype == 0) then
            irc = posix_rmtree(tmp_folder)                 ! best-effort (ENOENT if none) -> ignore status
            call zarr_mkdir(tmp_folder)                    ! recursive mkdir -p (also creates the parent)
        end if
        if (R%mr) call MPI_Barrier(R%comm, ierr)           ! tmp folder must exist before any store write
        ! 2) every writer stages its field stores INTO the tmp folder; rank 0 stages checkpoint.json.
        do k = 1, R%nf
            call restart_write_field(R, R%f(k), trim(tmp_folder)//'/'//trim(R%f(k)%name)//'.zarr')
        end do
        if (R%mype == 0) call write_checkpoint_json(R, trim(tmp_folder)//'/checkpoint.json', year, day, &
                                                    time_sec, gstep)
        if (R%mr) call MPI_Barrier(R%comm, ierr)           ! all stores + json fully staged before publish
        ! 3) + 4) rank 0: durably publish (fsync+rename), flip restart.latest atomically, then prune.
        if (R%mype == 0) then
            irc = posix_fsync_dir(tmp_folder)              ! durability of the staged dir before rename
            irc = posix_rename(tmp_folder, final_folder)   ! ATOMIC publish: fesom.<tag>.tmp -> fesom.<tag>
            call zarr_check(irc == 0, 'restart: rename tmp->final failed: '//trim(final_folder))
            irc = posix_fsync_dir(checkpoint_dir)          ! durably record the new dir entry
            call update_restart_latest(checkpoint_dir, fname)
            call restart_prune(R, checkpoint_dir, fname)
        end if
        if (R%mr) call MPI_Barrier(R%comm, ierr)
    end subroutine restart_write

    ! Atomically point restart.latest at <folder_name> (the bare folder NAME, e.g. fesom.2000.001.03600).
    ! Write restart.latest.tmp then rename(2) it onto restart.latest: rename is atomic within one
    ! filesystem, so a concurrent or crashing reader sees the OLD or NEW pointer — never a partial line.
    ! rank-0 only.
    subroutine update_restart_latest(restart_dir, folder_name)
        character(len=*), intent(in) :: restart_dir, folder_name
        character(len=:), allocatable :: latest, latest_tmp
        integer :: irc
        latest     = trim(restart_dir)//'/restart.latest'
        latest_tmp = trim(restart_dir)//'/restart.latest.tmp'
        call write_text_file(latest_tmp, trim(folder_name))
        irc = posix_rename(latest_tmp, latest)
        call zarr_check(irc == 0, 'restart: rename restart.latest.tmp -> restart.latest failed')
    end subroutine update_restart_latest

    ! keep-N prune (rank-0, best-effort). After a checkpoint is finalized, enumerate the immutable
    ! checkpoint folders under restart_dir, sort lexically (== chronologically — the tag is fixed-width
    ! zero-padded), and posix_rmtree the oldest beyond restart_keep. restart_keep<=0 keeps all. The
    ! CURRENT pointer target (keep_name) is ALWAYS protected, even if a stray later-named folder would
    ! rank ahead of it, so the prune can never delete the checkpoint restart.latest just committed.
    ! Stray *.tmp/ dirs do NOT match the strict checkpoint pattern (is_checkpoint_name), so they are
    ! never pruned and never resolved — they are left as-is (documented: the reader only trusts the
    ! pointer). ANY failure warns and continues — a prune problem must never abort a good checkpoint.
    subroutine restart_prune(R, restart_dir, keep_name)
        type(t_restart),  intent(in) :: R
        character(len=*), intent(in) :: restart_dir, keep_name
        integer, parameter :: MAXLIST = 4096
        character(len=64), allocatable :: names(:)
        character(len=64) :: swap
        integer :: n, i, j, irc
        logical :: ok
        character(len=:), allocatable :: folder
        if (R%restart_keep <= 0) return                    ! 0 => keep all
        allocate(names(MAXLIST))
        call posix_listdir(restart_dir, names, n, ok)
        if (.not. ok) then
            write(*,'(a)') '[restart] WARN: cannot list '//trim(restart_dir)//' for keep-N prune (skipped)'
            return
        end if
        ! keep only strict checkpoint folder names (excludes fesom.clock, *.zarr, restart.latest, *.tmp)
        j = 0
        do i = 1, n
            if (is_checkpoint_name(names(i))) then
                j = j + 1; names(j) = names(i)
            end if
        end do
        n = j
        if (n <= R%restart_keep) return                    ! nothing to prune
        ! selection sort ascending (lexical == chronological); n is tiny, clarity over speed.
        do i = 1, n - 1
            do j = i + 1, n
                if (names(j) < names(i)) then
                    swap = names(i); names(i) = names(j); names(j) = swap
                end if
            end do
        end do
        ! rmtree the oldest (n - restart_keep), but NEVER the active pointer target (keep_name).
        do i = 1, n - R%restart_keep
            if (trim(names(i)) == trim(keep_name)) cycle   ! protect the just-committed checkpoint
            folder = trim(restart_dir)//'/'//trim(names(i))
            irc = posix_rmtree(folder)
            if (irc /= 0) then
                write(*,'(a)') '[restart] WARN: keep-N prune could not fully remove '//trim(folder)
            else
                write(*,'(a)') '[restart] keep-N prune removed old checkpoint '//trim(names(i))
            end if
        end do
    end subroutine restart_prune

    ! True iff `name` is EXACTLY a checkpoint folder name: 'fesom.' + 4 digits + '.' + 3 digits + '.' +
    ! 5 digits (20 chars). This strict shape excludes fesom.clock, fesom.mesh.diag.zarr, restart.latest,
    ! and any fesom.<tag>.tmp staging dir — so the prune only ever touches real finalized checkpoints.
    logical function is_checkpoint_name(name)
        character(len=*), intent(in) :: name
        character(len=20) :: s
        integer :: i
        is_checkpoint_name = .false.
        if (len_trim(name) /= 20) return
        s = name(1:20)
        if (s(1:6) /= 'fesom.')               return
        if (s(11:11) /= '.' .or. s(15:15) /= '.') return
        do i = 7, 10
            if (s(i:i) < '0' .or. s(i:i) > '9') return
        end do
        do i = 12, 14
            if (s(i:i) < '0' .or. s(i:i) > '9') return
        end do
        do i = 16, 20
            if (s(i:i) < '0' .or. s(i:i) > '9') return
        end do
        is_checkpoint_name = .true.
    end function is_checkpoint_name

    ! Resolve the newest finalized checkpoint by FOLLOWING restart.latest — never by scanning the dir —
    ! so a stray crashed-write fesom.<tag>.tmp/ or an unpointed finalized fesom.<tag>/ is ignored by
    ! construction. Returns folder_out = <restart_dir>/<name-in-restart.latest> and ok=.true. iff
    ! restart.latest exists, is readable, names a non-empty folder, and that folder carries a
    ! checkpoint.json (the manifest that confirms a complete checkpoint). This is the entry point the
    ! Task 4.1 read path calls; any rank may call it (a pure file read).
    subroutine restart_resolve_latest(restart_dir, folder_out, ok)
        character(len=*),              intent(in)  :: restart_dir
        character(len=:), allocatable, intent(out) :: folder_out
        logical,                       intent(out) :: ok
        character(len=512) :: name
        character(len=:), allocatable :: latest
        integer :: u, ios
        logical :: ex
        ok = .false.
        folder_out = ''
        latest = trim(restart_dir)//'/restart.latest'
        inquire(file=latest, exist=ex)
        if (.not. ex) return
        open(newunit=u, file=latest, status='old', action='read', form='formatted', iostat=ios)
        if (ios /= 0) return
        read(u, '(a)', iostat=ios) name
        close(u)
        if (ios /= 0) return
        name = adjustl(name)
        if (len_trim(name) == 0) return
        folder_out = trim(restart_dir)//'/'//trim(name)
        ! trust the pointer, but confirm the target is a complete checkpoint (has its manifest).
        inquire(file=folder_out//'/checkpoint.json', exist=ex)
        ok = ex
    end subroutine restart_resolve_latest

    ! Write ONE field's snapshot store (general: node OR element, 2-D OR 3-D). Resolves the field's
    ! entity decomp + cached coords from R, creates a single-variable single-entity store (data var
    ! (entity) or (nlev, entity), <f8 by default, codec from the namelist), embeds lon/lat +
    ! _ARRAY_DIMENSIONS + UGRID attrs via the Stage 2 helper, then decomp_redistributes the live array
    ! to canonical writer buffers and writes the data chunks. Store-create ordering: rank 0 defines ->
    ! barrier -> every writer writes its chunks.
    subroutine restart_write_field(R, f, store_path)
        type(t_restart),       intent(in), target :: R
        type(t_restart_field), intent(in)         :: f
        character(len=*),      intent(in)         :: store_path
        type(t_io_decomp), pointer :: D
        real(WP),          pointer :: own_lon(:), own_lat(:)
        integer                    :: nO
        type(t_zarr_store) :: store
        type(t_zarr_array) :: a_data, a_vert, a_lon, a_lat
        type(t_zarr_attrs) :: gat, at
        real(WP), allocatable :: buf(:), buf3(:,:), stg2(:), stg3(:,:)
        integer :: c, lo, ierr, cv, nvc, vc, L0, cvn
        ! resolve the entity context (node vs element-centroid)
        if (f%entity == DECOMP_ELEM) then
            D => R%De; nO = R%nElemO; own_lon => R%lon_e; own_lat => R%lat_e
        else
            D => R%Dn; nO = R%nNodO; own_lon => R%lon_n; own_lat => R%lat_n
        end if
        cv = vchunk_eff(R%chunk_vert, max(1, f%nlev))
        ! ---- array handles (all ranks) ----
        if (f%ndim == 2) then
            call zarr_array_init(a_data, trim(f%name), [D%N], [D%C], trim(f%dtype), &
                                 has_fill=.false., codec=trim(R%compressor))
        else
            call zarr_array_init(a_data, trim(f%name), [f%nlev, D%N], [cv, D%C], trim(f%dtype), &
                                 has_fill=.false., codec=trim(R%compressor))
            call zarr_array_init(a_vert, trim(f%vdim), [f%nlev], [f%nlev], '<f8', has_fill=.false.)
        end if
        call io_coords_init_lonlat(a_lon, a_lat, D%N, D%C)
        ! the chunk writers (every writer) reference store%path -> set on ALL ranks
        store%path = trim(store_path)
        ! ---- rank 0 defines the store + arrays + attrs (store-create ordering) ----
        if (R%mype == 0) then
            call zarr_attrs_init(gat)
            call zattr_str(gat, 'Conventions', 'CF-1.8')
            call zattr_str(gat, 'description', 'FESOM3 restart checkpoint field (Zarr v2 snapshot)')
            call zarr_create_store(store, trim(store_path), gat)
            call zarr_attrs_init(at)
            if (f%ndim == 2) then
                call zattr_str_arr(at, '_ARRAY_DIMENSIONS', [character(len=8) :: f%hdim])
            else
                call zattr_str_arr(at, '_ARRAY_DIMENSIONS', [character(len=8) :: f%vdim, f%hdim])
            end if
            if (len_trim(f%units) > 0) call zattr_str(at, 'units', trim(f%units))
            call zattr_str(at, 'coordinates', 'lon lat')
            call zarr_define_array(store, a_data, at)
            if (f%ndim == 3) then
                call zarr_attrs_init(at)
                call zattr_str_arr(at, '_ARRAY_DIMENSIONS', [character(len=8) :: f%vdim])
                call zattr_str(at, 'long_name', 'depth'); call zattr_str(at, 'units', 'meters')
                call zattr_str(at, 'positive', 'down')
                call zarr_define_array(store, a_vert, at)
            end if
            call io_coords_define_lonlat(store, a_lon, a_lat, f%hdim)
        end if
        if (R%mr) call MPI_Barrier(R%comm, ierr)
        ! ---- every writer redistributes the live field + writes its data chunks ----
        if (f%ndim == 2) then
            allocate(buf(max(1, D%w_nbuf)))
            if (f%mp_src) then                                              ! MP->WP staged (lossless)
                allocate(stg2(max(1, nO))); stg2(1:nO) = real(f%pmp2d(1:nO), WP)
                call decomp_redistribute(D, stg2(1:nO), buf, 0.0_WP)        ! collective
            else
                call decomp_redistribute(D, f%p2d(1:nO), buf, 0.0_WP)       ! collective
            end if
            do c = D%w_first_chunk, D%w_last_chunk
                lo = (c - D%w_first_chunk)*D%C + 1
                call zarr_write_chunk(store, a_data, [c], buf(lo:lo + D%C - 1))
            end do
        else
            allocate(buf3(f%nlev, max(1, D%w_nbuf)))
            if (f%mp_src) then                                              ! MP->WP staged (lossless)
                allocate(stg3(f%nlev, max(1, nO)))
                stg3(1:f%nlev, 1:nO) = real(f%pmp3d(1:f%nlev, 1:nO), WP)
                call decomp_redistribute(D, stg3(1:f%nlev, 1:nO), buf3, 0.0_WP)  ! collective
            else
                call decomp_redistribute(D, f%p3d(1:f%nlev,1:nO), buf3, 0.0_WP)  ! collective
            end if
            nvc = (f%nlev + cv - 1)/cv
            do c = D%w_first_chunk, D%w_last_chunk
                lo = (c - D%w_first_chunk)*D%C + 1
                do vc = 0, nvc - 1
                    L0 = vc*cv; cvn = min(cv, f%nlev - L0)
                    call zarr_write_chunk(store, a_data, [vc, c], buf3(L0+1:L0+cvn, lo:lo + D%C - 1))
                end do
            end do
        end if
        ! embedded lon/lat coord data (shared helper redistributes + writes per-chunk on each writer)
        call io_coords_put(D, store, a_lon, own_lon)
        call io_coords_put(D, store, a_lat, own_lat)
        ! vertical coord is global/replicated -> rank 0 writes it whole
        if (R%mype == 0 .and. f%ndim == 3) then
            if (trim(f%vdim) == 'nz') then
                call zarr_write_whole(store, a_vert, R%depth_nz(1:R%nl))
            else
                call zarr_write_whole(store, a_vert, R%depth_nz1(1:R%nl-1))
            end if
        end if
        if (R%mr) call MPI_Barrier(R%comm, ierr)
    end subroutine restart_write_field

    subroutine restart_finalize(R)
        type(t_restart), intent(in) :: R
        integer :: ierr
        if (R%mr) call MPI_Barrier(R%comm, ierr)
    end subroutine restart_finalize

    ! ----------------------------------------------------------------- checkpoint READ (Stage 4, Task 4.1)
    !
    ! The exact INVERSE of the write path. For every registered field: writer ranks read their owned
    ! canonical chunks (mirroring the write chunk loop), decomp_gather (the literal transpose of the
    ! write's decomp_redistribute) pulls them back to the compute partition's owned entities, those are
    ! written INTO the live model array's owned slots, and the halo is reconstructed by the SAME exchange
    ! the field's in-step consumer uses. Owned values round-trip max|Δ|=0 (the gather is exact); halos
    ! become owner-derived copies => the restored full array is bit-identical to a straight-through run
    ! at the same instant. The store holds only canonical OWNED values, so partition-independence falls
    ! out: a different np builds a different decomp plan over the SAME canonical chunks.

    ! Restore the FULL prognostic state from the newest finalized checkpoint under restart_dir. ABORTS
    ! (clear message) only if restart.latest is missing/unreadable, names an incomplete checkpoint, or a
    ! registered field's store is missing/corrupt — exactly as clock_init aborts on a missing .clock. On
    ! a checkpoint.json time that differs from the current clock it WARNs and CONTINUES (a legitimate
    ! dt-change restart trips it; oracle io_restart.F90:904-914 warns, never aborts). mesh is accepted
    ! for signature symmetry with the lifecycle (halo reconstruction needs only partit's com structures);
    ! the optional clock args drive the time-vs-clock safety warning.
    subroutine restart_read(R, restart_dir, mesh, partit, clock_year, clock_day, clock_time_sec)
        type(t_restart),  intent(in), target  :: R
        character(len=*), intent(in)           :: restart_dir
        type(t_mesh),     intent(in)           :: mesh
        type(t_partit),   intent(in)           :: partit
        integer,          intent(in), optional :: clock_year, clock_day
        real(real64),     intent(in), optional :: clock_time_sec
        character(len=:), allocatable :: folder
        logical :: ok, have_json
        integer :: k, jyear, jday, ierr
        real(real64) :: jtime
        call restart_resolve_latest(restart_dir, folder, ok)
        call zarr_check(ok, 'restart_read: restart.latest missing/unreadable or checkpoint incomplete under '// &
                            trim(restart_dir))
        if (R%mype == 0) write(*,'(a)') '[restart] reading checkpoint '//folder
        ! time-vs-clock safety check (WARN, never abort)
        if (present(clock_time_sec)) then
            call read_checkpoint_json(folder//'/checkpoint.json', jyear, jday, jtime, have_json)
            if (have_json .and. abs(jtime - clock_time_sec) > 0.0_real64 .and. R%mype == 0) then
                write(*,'(a,es24.16,a,es24.16,a)') '[restart] WARN: checkpoint time_sec=', jtime, &
                    ' differs from clock time_sec=', clock_time_sec, ' (continuing)'
            end if
        end if
        ! restore + halo-reconstruct every field
        do k = 1, R%nf
            if (R%f(k)%skip_read) cycle      ! RECOM4 M10: same flag on every rank
            call restart_read_field(R, R%f(k), folder//'/'//trim(R%f(k)%name)//'.zarr', partit)
        end do
        if (R%mr) call MPI_Barrier(R%comm, ierr)
    end subroutine restart_read

    ! Read ONE field's snapshot store back into its live model array, then reconstruct its halo. WRITER
    ! ranks read their owned canonical chunks ([w_first_chunk..w_last_chunk]; 3-D loops the vertical
    ! chunks exactly like the writer); ALL ranks then decomp_gather the canonical buffer to owned values,
    ! write them into the live array's owned slots (MP-source fields via a lossless WP->MP copy), and
    ! finally halo-exchange (restart_halo_exchange_field) so the restored array matches a straight-through
    ! run over the extent that variant refreshes.
    subroutine restart_read_field(R, f, store_path, partit)
        type(t_restart),       intent(in), target :: R
        type(t_restart_field), intent(in)         :: f
        character(len=*),      intent(in)         :: store_path
        type(t_partit),        intent(in)         :: partit
        type(t_io_decomp), pointer :: D
        integer                    :: nO
        type(t_zarr_store) :: store
        type(t_zarr_array) :: a_data
        real(WP), allocatable :: buf(:), buf3(:,:), bo(:), bo3(:,:), rdc(:,:)
        integer :: c, lo, cv, nvc, vc, L0, cvn
        logical :: ex
        ! entity context (node vs element decomp + owned count)
        if (f%entity == DECOMP_ELEM) then; D => R%De; nO = R%nElemO
        else;                              D => R%Dn; nO = R%nNodO; end if
        ! abort on a missing/corrupt store (mirrors clock_init's missing-.clock abort)
        inquire(file=trim(store_path)//'/.zgroup', exist=ex)
        call zarr_check(ex, 'restart_read: missing/corrupt store '//trim(store_path))
        store%path = trim(store_path)
        if (f%ndim == 2) then
            ! data array (entity,) chunked (C); build the handle IDENTICALLY to the writer so chunk
            ! paths + chunk byte sizes match exactly.
            call zarr_array_init(a_data, trim(f%name), [D%N], [D%C], trim(f%dtype), &
                                 has_fill=.false., codec=trim(R%compressor))
            allocate(buf(max(1, D%w_nbuf)), bo(max(1, nO)))
            do c = D%w_first_chunk, D%w_last_chunk                ! writers only (else empty range)
                lo = (c - D%w_first_chunk)*D%C + 1
                call zarr_read_chunk(store, a_data, [c], buf(lo:lo + D%C - 1))
            end do
            call decomp_gather(D, buf, bo)                        ! collective: canonical -> owned
            if (f%mp_src) then
                f%pmp2d(1:nO) = real(bo(1:nO), MP)               ! WP -> MP (lossless at WP>=MP)
            else
                f%p2d(1:nO) = bo(1:nO)
            end if
        else
            ! data array (nlev, entity) chunked (cv, C)
            cv = vchunk_eff(R%chunk_vert, max(1, f%nlev))
            call zarr_array_init(a_data, trim(f%name), [f%nlev, D%N], [cv, D%C], trim(f%dtype), &
                                 has_fill=.false., codec=trim(R%compressor))
            allocate(buf3(f%nlev, max(1, D%w_nbuf)), bo3(f%nlev, max(1, nO)), rdc(cv, D%C))
            nvc = (f%nlev + cv - 1)/cv
            do c = D%w_first_chunk, D%w_last_chunk
                lo = (c - D%w_first_chunk)*D%C + 1
                do vc = 0, nvc - 1
                    L0 = vc*cv; cvn = min(cv, f%nlev - L0)
                    call zarr_read_chunk(store, a_data, [vc, c], rdc) ! reads the FULL (cv,C) chunk
                    buf3(L0+1:L0+cvn, lo:lo + D%C - 1) = rdc(1:cvn, 1:D%C)  ! keep the valid rows
                end do
            end do
            call decomp_gather(D, buf3, bo3)
            if (f%mp_src) then
                f%pmp3d(1:f%nlev, 1:nO) = real(bo3(1:f%nlev, 1:nO), MP)
            else
                f%p3d(1:f%nlev, 1:nO) = bo3(1:f%nlev, 1:nO)
            end if
        end if
        ! reconstruct the halo with the field's chosen variant (no-op at np=1: no halo cells exist)
        call restart_halo_exchange_field(f, partit)
    end subroutine restart_read_field

    ! Reconstruct field f's halo IN-PLACE with the exchange variant its in-step consumer uses, so the
    ! restored array is bit-identical to a straight-through run over the extent that variant refreshes.
    ! MP-source fields (mesh%hbar/hnode) exchange via a WP staging copy (exchange_* is WP-typed): copy
    ! MP->WP, exchange, copy WP->MP (lossless at WP=real64>=MP). f is intent(in) but its POINTER targets
    ! (the live arrays) are modified — allowed. At single rank there are no halo cells (return early).
    subroutine restart_halo_exchange_field(f, partit)
        type(t_restart_field), intent(in) :: f
        type(t_partit),        intent(in) :: partit
        real(WP), allocatable :: s2(:), s3(:,:)
        integer :: n1, n2
        if (f%halo == RESTART_HALO_NONE) return                  ! owned-only consumer: halo stays as-is
        if (.not. is_multirank(partit)) return                   ! np=1: array is wholly owned
        select case (f%halo)
        case (RESTART_HALO_NODE)
            if (f%mp_src) then
                if (f%ndim == 2) then
                    n1 = size(f%pmp2d); allocate(s2(n1)); s2 = real(f%pmp2d, WP)
                    call exchange_nod(s2, partit);                     f%pmp2d = real(s2, MP)
                else
                    n1 = size(f%pmp3d,1); n2 = size(f%pmp3d,2); allocate(s3(n1,n2)); s3 = real(f%pmp3d, WP)
                    call exchange_nod(s3, partit);                     f%pmp3d = real(s3, MP)
                end if
            else
                if (f%ndim == 2) then; call exchange_nod(f%p2d, partit)
                else;                  call exchange_nod(f%p3d, partit); end if
            end if
        case (RESTART_HALO_ELEM)        ! eDim com_elem2D (dyn%uv u/v); broadcast-only => per-component == blk
            if (f%ndim == 2) then; call exchange_elem(f%p2d, partit)
            else;                  call exchange_elem(f%p3d, partit); end if
        case (RESTART_HALO_ELEM_FULL)   ! eDim+eXDim com_elem2D_full (rank-1 element fields only —
            !                             exchange_elem_full has no (nlev,elem) rank-2 variant)
            if (f%ndim == 2) then; call exchange_elem_full(f%p2d, partit)
            else; call zarr_check(.false., 'restart: ELEM_FULL halo on a (nlev,elem) field is '// &
                                           'unsupported (no rank-2 exchange_elem_full): '//trim(f%name)); end if
        end select
    end subroutine restart_halo_exchange_field

    ! Reconstruct EVERY registered field's halo. Used by the read path's per-field loop indirectly, and
    ! directly by the Task 4.1 round-trip gate to make its REFERENCE halos owner-derived with the exact
    ! same logic the read path uses (so reference and read-back halos are built identically).
    subroutine restart_halo_exchange_all(R, partit)
        type(t_restart), intent(in) :: R
        type(t_partit),  intent(in) :: partit
        integer :: k
        do k = 1, R%nf
            call restart_halo_exchange_field(R%f(k), partit)
        end do
    end subroutine restart_halo_exchange_all

    ! Minimal line-based reader for the manifest keys the safety check needs (year/day/time_sec) from the
    ! plain document write_checkpoint_json emits (one "key": value per line). ok=.false. if unreadable.
    ! Not a general JSON parser — only the keys we wrote.
    subroutine read_checkpoint_json(path, year, day, time_sec, ok)
        character(len=*), intent(in)  :: path
        integer,          intent(out) :: year, day
        real(real64),     intent(out) :: time_sec
        logical,          intent(out) :: ok
        integer :: u, ios
        character(len=256) :: line
        year = 0; day = 0; time_sec = 0.0_real64; ok = .false.
        open(newunit=u, file=trim(path), status='old', action='read', form='formatted', iostat=ios)
        if (ios /= 0) return
        do
            read(u, '(a)', iostat=ios) line
            if (ios /= 0) exit
            if (index(line, '"year"')     > 0) call json_value_int(line, year)
            if (index(line, '"day"')      > 0) call json_value_int(line, day)
            if (index(line, '"time_sec"') > 0) call json_value_real(line, time_sec)
        end do
        close(u)
        ok = .true.
    end subroutine read_checkpoint_json

    subroutine json_value_int(line, v)
        character(len=*), intent(in)  :: line
        integer,          intent(out) :: v
        character(len=64) :: s
        integer :: p, ios
        v = 0
        p = index(line, ':'); if (p == 0) return
        s = adjustl(line(p+1:))
        p = index(s, ','); if (p > 0) s(p:) = ' '
        read(s, *, iostat=ios) v
    end subroutine json_value_int

    subroutine json_value_real(line, v)
        character(len=*), intent(in)  :: line
        real(real64),     intent(out) :: v
        character(len=64) :: s
        integer :: p, ios
        v = 0.0_real64
        p = index(line, ':'); if (p == 0) return
        s = adjustl(line(p+1:))
        p = index(s, ','); if (p > 0) s(p:) = ' '
        read(s, *, iostat=ios) v
    end subroutine json_value_real

    ! ----------------------------------------------------------------- checkpoint.json (rank 0)

    ! Plain Fortran string writes (no JSON library): the small provenance manifest the reader uses for
    ! the time-vs-clock safety check + npes provenance. format_version/year/day/globalstep/npes_wrote
    ! are integers; time_sec a real (the sec-of-day passed in); fesom_git from $FESOM3_GIT else 'unknown'.
    subroutine write_checkpoint_json(R, path, year, day, time_sec, gstep)
        type(t_restart),  intent(in) :: R
        character(len=*), intent(in) :: path
        integer,          intent(in) :: year, day, gstep
        real(real64),     intent(in) :: time_sec
        character(len=64) :: git
        character(len=32) :: rbuf
        character(len=:), allocatable :: doc, rstr
        write(rbuf, '(es24.16)') time_sec
        rstr = trim(adjustl(rbuf))
        git = ''; call get_environment_variable('FESOM3_GIT', git)
        if (len_trim(git) == 0) git = 'unknown'
        doc = '{'//char(10)// &
              '  "format_version": '//itoa(RESTART_FORMAT_VERSION)//','//char(10)// &
              '  "year": '//itoa(year)//','//char(10)// &
              '  "day": '//itoa(day)//','//char(10)// &
              '  "time_sec": '//rstr//','//char(10)// &
              '  "globalstep": '//itoa(gstep)//','//char(10)// &
              '  "fesom_git": "'//trim(git)//'",'//char(10)// &
              '  "npes_wrote": '//itoa(R%npes)//char(10)// &
              '}'
        call write_text_file(path, doc)
    end subroutine write_checkpoint_json

    ! ----------------------------------------------------------------- small helpers

    ! map a precision arg ('double'|'<f8'|'8' => <f8; '<f4'|'4'|'single' => <f4) to a Zarr dtype.
    ! restart DEFAULT = <f8 (full precision; lossy restart is not a v1 option).
    function restart_dtype(precision) result(dt)
        character(len=*), intent(in), optional :: precision
        character(len=8) :: dt
        dt = '<f8'
        if (present(precision)) then
            if (trim(precision) == '<f4' .or. trim(precision) == '4' .or. trim(precision) == 'single') &
                dt = '<f4'
        end if
    end function restart_dtype

    ! effective vertical chunk: full depth when chunk_vert is 0 or >= nlev, else chunk_vert.
    integer function vchunk_eff(cv, nlev)
        integer, intent(in) :: cv, nlev
        if (cv <= 0 .or. cv >= nlev) then
            vchunk_eff = max(1, nlev)
        else
            vchunk_eff = cv
        end if
    end function vchunk_eff

    function itoa(i) result(s)
        integer, intent(in)           :: i
        character(len=:), allocatable :: s
        character(len=32) :: tmp
        write(tmp, '(i0)') i
        s = trim(tmp)
    end function itoa

    ! Write a text blob verbatim + a trailing newline via stream (any length is safe). Mirrors
    ! mod_io_zarr's private write_text_file (kept local so mod_io_restart owns its manifest I/O).
    subroutine write_text_file(path, text)
        character(len=*), intent(in) :: path, text
        integer :: u, ios
        open(newunit=u, file=trim(path), status='replace', action='write', &
             form='unformatted', access='stream', iostat=ios)
        call zarr_check(ios == 0, 'restart: open(write) '//trim(path))
        write(u) text//char(10)
        close(u)
    end subroutine write_text_file

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

end module mod_io_restart
