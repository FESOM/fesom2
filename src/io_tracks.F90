!> Ship-track / mooring-array curtain output for FESOM-2, v2 (geometry-aware).
!>
!> Each user-supplied polyline CSV is converted to a list of "crossed
!> edges": the actual edges of the FESOM mesh whose interior the polyline
!> passes through. For each crossed edge a sample value is built by
!> linearly interpolating the two endpoint values
!>     V_at_crossing = V(n1) * (1 - t) + V(n2) * t,
!> with t in [0,1] from the line-edge intersection. This replaces the
!> earlier snap-to-nearest-node approach which produced mesh-scale
!> zig-zag artefacts and could not compute transport.
!>
!> Geometry pipeline (lives in mod_tracks_geometry, pure Fortran):
!>   * bbox + polar widening
!>   * signed-distance line-edge intersection
!>   * alternating triangle path with up/down section dx/dy
!>   * cumulative great-circle distance
!>
!> MPI layout:
!>   * Rank 0 gathers the global topology (edges, edge_tri, face_nodes,
!>     edge_cross_dxdy, lon/lat) using FESOM's existing gather_* primitives
!>   * Rank 0 runs analyse_transect once per track
!>   * The resulting transect_t (M <= ~few thousand) is broadcast to all
!>     ranks
!>   * Each rank precomputes which edge endpoints it strict-owns
!>   * At each output step: per-rank contributions are summed via
!>     MPI_Allreduce; slots where no rank contributed (boundary edge or
!>     both endpoints dry at this level) become NaN
!>
!> Activation: ltracks=.true. plus a non-empty track_files. XML override
!> in io_xios_init wins over the (skipped in XIOS mode) namelist.
module io_tracks_module
  ! Configuration variables are always declared so namelist.io can parse
  ! the &nml_tracks block in any build. The XIOS-using subroutines below
  ! are compiled in only when __XIOS is defined.
#if defined(__XIOS)
  use xios
  use ixml_tree,           only: xios_add_axistogrid, xios_add_domaintogrid
  use mpi
  use mod_mesh,            only: t_mesh
  use mod_partit,          only: t_partit
  use mod_tracer,          only: t_tracer
  use o_param,             only: WP, rad
  use g_comm_auto,         only: gather_nod, gather_elem, gather_edge
  use mod_tracks_geometry, only: transect_t, analyse_transect, free_transect
  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
#endif
  implicit none
  private

  logical,             save, public :: ltracks               = .false.
  character(len=4096), save, public :: track_files           = ''
  character(len=512),  save, public :: track_vars            = 'temp'
  character(len=512),  save, public :: track_names           = ''
  character(len=16),   save, public :: track_output_freq     = '1h'

  public :: io_tracks_is_on
#if defined(__XIOS)
  public :: io_tracks_register_xios
  public :: io_tracks_send

  ! Per-track state.
  type :: track_t
    character(len=64)    :: name       = ''
    character(len=512)   :: csv_path   = ''
    character(len=32)    :: var_name   = 'temp'
    integer              :: tracer_idx = 1            ! 1=temp, 2=salt
    character(len=64)    :: field_id   = ''

    ! Geometry — broadcast from rank 0 to all ranks after analyse_transect.
    integer              :: M          = 0            ! # crossed edges (global)
    integer, allocatable :: edge_cut_ni(:, :)         ! (2, M) global node ids
    real(WP), allocatable :: edge_cut_lint(:)         ! (M)
    real(WP), allocatable :: lon_mid(:), lat_mid(:)   ! (M) midpoint coords (deg)
    real(WP), allocatable :: dist_km(:)               ! (M)

    ! Per-rank sampling lookups. lid_*(m) > 0 iff this rank strict-owns
    ! endpoint *; the lid value indexes nodal arrays (1..myDim_nod2D).
    integer, allocatable :: lid_n1(:)                 ! (M) 0 if not owned
    integer, allocatable :: lid_n2(:)                 ! (M) 0 if not owned

    ! XIOS-side: each rank "owns" a stride of m's for write purposes.
    ! nL = # m's this rank pushes to XIOS; my_m(k) = global m index (1-based)
    ! and i_index(k) = m - 1 (0-based for XIOS).
    integer              :: nL         = 0
    integer, allocatable :: my_m(:)                   ! (nL)
    integer, allocatable :: i_index(:)                ! (nL)
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

  !> Parse semicolon-separated lists, gather the global mesh on rank 0,
  !> run analyse_transect per track, broadcast geometry, register the
  !> per-track XIOS domain/grid/field.
  subroutine io_tracks_register_xios(mesh, partit, nz_cell)
    type(t_mesh),   intent(in)    :: mesh
    type(t_partit), intent(inout) :: partit
    integer,        intent(in)    :: nz_cell
    logical, save :: already = .false.
    integer :: i, ierr

    ! Global mesh assembled on rank 0
    real(WP), allocatable :: lon_glo(:), lat_glo(:)
    integer,  allocatable :: edges_glo(:, :), edge_tri_glo(:, :), elem_nodes_glo(:, :)
    real(WP), allocatable :: edge_cross_dxdy_glo(:, :)

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

    call gather_global_mesh(mesh, partit, &
                            lon_glo, lat_glo, edges_glo, edge_tri_glo, &
                            elem_nodes_glo, edge_cross_dxdy_glo)

    call register_track_filegroup(partit)

    do i = 1, n_tracks
       call process_one_track(mesh, partit, nz_cell, tracks(i), &
                              lon_glo, lat_glo, edges_glo, edge_tri_glo, &
                              elem_nodes_glo, edge_cross_dxdy_glo)
    end do

    ! Drop the global temps once all tracks have their geometry.
    if (allocated(lon_glo))            deallocate(lon_glo)
    if (allocated(lat_glo))            deallocate(lat_glo)
    if (allocated(edges_glo))          deallocate(edges_glo)
    if (allocated(edge_tri_glo))       deallocate(edge_tri_glo)
    if (allocated(elem_nodes_glo))     deallocate(elem_nodes_glo)
    if (allocated(edge_cross_dxdy_glo)) deallocate(edge_cross_dxdy_glo)
  end subroutine io_tracks_register_xios


  !> Loop tracks, sample, MPI_Allreduce, send via XIOS.
  subroutine io_tracks_send(tracers, mesh, partit)
    type(t_tracer), intent(in) :: tracers
    type(t_mesh),   intent(in) :: mesh
    type(t_partit), intent(in) :: partit
    integer :: i

    if (.not. io_tracks_is_on()) return
    if (n_tracks == 0) return

    do i = 1, n_tracks
       call send_one_track(tracers, mesh, partit, tracks(i))
    end do
  end subroutine io_tracks_send


  ! ====================================================================
  ! Rank-0 mesh gather
  ! ====================================================================

  !> Use FESOM's gather_edge / gather_elem / gather_nod primitives to
  !> assemble the global topology on rank 0. Connectivity arrays carry
  !> LOCAL ids per rank; translate to GLOBAL ids before gather (the
  !> pattern at io_mesh_info.F90:441-526).
  !>
  !> Non-rank-0 receive dummy size-1 arrays — the gather routines write
  !> to arr_global only on rank 0.
  subroutine gather_global_mesh(mesh, partit, &
                                lon_glo, lat_glo, edges_glo, edge_tri_glo, &
                                elem_nodes_glo, edge_cross_dxdy_glo)
    type(t_mesh),   intent(in)    :: mesh
    type(t_partit), intent(inout) :: partit
    real(WP), allocatable, intent(out) :: lon_glo(:), lat_glo(:)
    integer,  allocatable, intent(out) :: edges_glo(:, :), edge_tri_glo(:, :)
    integer,  allocatable, intent(out) :: elem_nodes_glo(:, :)
    real(WP), allocatable, intent(out) :: edge_cross_dxdy_glo(:, :)

    integer :: nod2D_g, edge2D_g, elem2D_g, k, j
    integer :: myDim_nod2D, myDim_edge2D, myDim_elem2D
    integer, allocatable :: lbuf_i(:), tmp_col(:)
    real(WP), allocatable :: lbuf_r(:), tmp_col_r(:)

    nod2D_g  = mesh%nod2D
    edge2D_g = mesh%edge2D
    elem2D_g = mesh%elem2D

    myDim_nod2D   = partit%myDim_nod2D
    myDim_edge2D  = partit%myDim_edge2D
    myDim_elem2D  = partit%myDim_elem2D

    if (partit%mype == 0) then
       allocate(lon_glo(nod2D_g), lat_glo(nod2D_g))
       allocate(edges_glo(2, edge2D_g))
       allocate(edge_tri_glo(2, edge2D_g))
       allocate(elem_nodes_glo(3, elem2D_g))
       allocate(edge_cross_dxdy_glo(4, edge2D_g))
    else
       allocate(lon_glo(1), lat_glo(1))
       allocate(edges_glo(1, 1))
       allocate(edge_tri_glo(1, 1))
       allocate(elem_nodes_glo(1, 1))
       allocate(edge_cross_dxdy_glo(1, 1))
    end if

    ! ---- geo_coord_nod2D (radians on the model side; convert to degrees later) ----
    allocate(lbuf_r(myDim_nod2D), tmp_col_r(nod2D_g))
    do k = 1, myDim_nod2D
       lbuf_r(k) = mesh%geo_coord_nod2D(1, k)
    end do
    if (partit%mype == 0) then
       call gather_nod(lbuf_r, tmp_col_r, partit)
       lon_glo = tmp_col_r / rad     ! rad-per-degree -> degrees
    else
       call gather_nod(lbuf_r, tmp_col_r, partit)
    end if
    do k = 1, myDim_nod2D
       lbuf_r(k) = mesh%geo_coord_nod2D(2, k)
    end do
    if (partit%mype == 0) then
       call gather_nod(lbuf_r, tmp_col_r, partit)
       lat_glo = tmp_col_r / rad
    else
       call gather_nod(lbuf_r, tmp_col_r, partit)
    end if
    deallocate(lbuf_r, tmp_col_r)

    ! ---- edges (LOCAL node ids -> GLOBAL via myList_nod2D) ----
    allocate(lbuf_i(myDim_edge2D), tmp_col(edge2D_g))
    do j = 1, 2
       do k = 1, myDim_edge2D
          lbuf_i(k) = partit%myList_nod2D(mesh%edges(j, k))
       end do
       call gather_edge(lbuf_i, tmp_col, partit)
       if (partit%mype == 0) edges_glo(j, :) = tmp_col
    end do

    ! ---- edge_tri (LOCAL elem ids -> GLOBAL via myList_elem2D; -1 = boundary) ----
    do j = 1, 2
       do k = 1, myDim_edge2D
          if (mesh%edge_tri(j, k) > 0) then
             lbuf_i(k) = partit%myList_elem2D(mesh%edge_tri(j, k))
          else
             lbuf_i(k) = -1
          end if
       end do
       call gather_edge(lbuf_i, tmp_col, partit)
       if (partit%mype == 0) edge_tri_glo(j, :) = tmp_col
    end do
    deallocate(lbuf_i, tmp_col)

    ! ---- elem_nodes (LOCAL node ids -> GLOBAL) ----
    allocate(lbuf_i(myDim_elem2D), tmp_col(elem2D_g))
    do j = 1, 3
       do k = 1, myDim_elem2D
          lbuf_i(k) = partit%myList_nod2D(mesh%elem2D_nodes(j, k))
       end do
       call gather_elem(lbuf_i, tmp_col, partit)
       if (partit%mype == 0) elem_nodes_glo(j, :) = tmp_col
    end do
    deallocate(lbuf_i, tmp_col)

    ! ---- edge_cross_dxdy (no id translation; values only) ----
    allocate(lbuf_r(myDim_edge2D), tmp_col_r(edge2D_g))
    do j = 1, 4
       do k = 1, myDim_edge2D
          lbuf_r(k) = mesh%edge_cross_dxdy(j, k)
       end do
       call gather_edge(lbuf_r, tmp_col_r, partit)
       if (partit%mype == 0) edge_cross_dxdy_glo(j, :) = tmp_col_r
    end do
    deallocate(lbuf_r, tmp_col_r)
  end subroutine gather_global_mesh


  ! ====================================================================
  ! Per-track pipeline
  ! ====================================================================

  subroutine process_one_track(mesh, partit, nz_cell, t, &
                               lon_glo, lat_glo, edges_glo, edge_tri_glo, &
                               elem_nodes_glo, edge_cross_dxdy_glo)
    type(t_mesh),   intent(in)    :: mesh
    type(t_partit), intent(inout) :: partit
    integer,        intent(in)    :: nz_cell
    type(track_t),  intent(inout) :: t
    real(WP),       intent(in)    :: lon_glo(:), lat_glo(:)
    integer,        intent(in)    :: edges_glo(:, :), edge_tri_glo(:, :)
    integer,        intent(in)    :: elem_nodes_glo(:, :)
    real(WP),       intent(in)    :: edge_cross_dxdy_glo(:, :)

    real(WP), allocatable :: lon_csv(:), lat_csv(:)
    type(transect_t)      :: geom
    integer :: N_csv, ierr

    ! --- tracer slot ---------------------------------------------------
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

    ! --- parse CSV on rank 0; analyse_transect there; broadcast --------
    if (partit%mype == 0) then
       call parse_csv(trim(t%csv_path), N_csv, lon_csv, lat_csv, partit)
       write(*,'(a,a,a,a,a,i0,a)') '[TRACKS] ', trim(t%name), ': parsed ', &
            trim(t%csv_path), ' -> ', N_csv, ' polyline vertices'

       call analyse_transect(size(lon_glo),     lon_glo, lat_glo,             &
                             size(edges_glo, 2), edges_glo, edge_tri_glo,     &
                             edge_cross_dxdy_glo,                              &
                             size(elem_nodes_glo, 2), elem_nodes_glo,         &
                             N_csv, lon_csv, lat_csv, trim(t%name), geom)

       deallocate(lon_csv, lat_csv)
       write(*,'(a,a,a,i0,a)') '[TRACKS] ', trim(t%name),                    &
            ': geometry M=', geom%M, ' crossed edges'
    end if

    call broadcast_transect_min(geom, partit)
    if (geom%M == 0) then
       if (partit%mype == 0) write(*,'(a,a,a)') '[TRACKS] WARNING: track "', &
            trim(t%name), '" has no mesh crossings — skipping.'
       call free_transect(geom)
       return
    end if

    ! --- store the small geometry pieces we need at sample time ---------
    t%M = geom%M
    allocate(t%edge_cut_ni(2, t%M))
    allocate(t%edge_cut_lint(t%M))
    allocate(t%lon_mid(t%M), t%lat_mid(t%M))
    allocate(t%dist_km(t%M))
    t%edge_cut_ni   = geom%edge_cut_ni
    t%edge_cut_lint = geom%edge_cut_lint
    t%lon_mid       = geom%edge_cut_midP(1, :)
    t%lat_mid       = geom%edge_cut_midP(2, :)
    t%dist_km       = geom%edge_cut_dist

    call free_transect(geom)

    ! --- per-rank: for each m, find local strict-owned lid for n1, n2 ---
    call build_sampling_lookups(mesh, partit, t)

    ! --- per-rank XIOS assignment: m mod npes -> my_m / i_index ---------
    call assign_xios_cells(partit, t)

    ! --- register XIOS objects ------------------------------------------
    call register_xios_objects_for_track(partit, t)

    if (partit%mype == 0) then
       write(*,'(a,a,a,i0,a)') '[TRACKS] ', trim(t%name),                    &
            ': registered runtime XML, M=', t%M, ' cells'
       write(*,'(a,a,a,a)')    '[TRACKS]   field=', trim(t%field_id),         &
                               '  file=', trim(t%field_id) // '.fesom'
    end if
  end subroutine process_one_track


  !> Broadcast only the fields we sample at runtime (edge_cut_ni,
  !> edge_cut_lint, edge_cut_midP, edge_cut_dist) plus M. The richer
  !> transect_t produced by analyse_transect carries path-side arrays
  !> that Step 1 doesn't need; we drop them after geometry runs.
  subroutine broadcast_transect_min(geom, partit)
    type(transect_t), intent(inout) :: geom
    type(t_partit),   intent(in)    :: partit
    integer :: ierr, M_b

    M_b = geom%M
    call mpi_bcast(M_b, 1, MPI_INTEGER, 0, partit%MPI_COMM_FESOM, ierr)
    if (M_b == 0) then
       geom%M = 0
       return
    end if

    if (partit%mype /= 0) then
       if (allocated(geom%edge_cut_ni))   deallocate(geom%edge_cut_ni)
       if (allocated(geom%edge_cut_lint)) deallocate(geom%edge_cut_lint)
       if (allocated(geom%edge_cut_midP)) deallocate(geom%edge_cut_midP)
       if (allocated(geom%edge_cut_dist)) deallocate(geom%edge_cut_dist)
       allocate(geom%edge_cut_ni(2, M_b))
       allocate(geom%edge_cut_lint(M_b))
       allocate(geom%edge_cut_midP(2, M_b))
       allocate(geom%edge_cut_dist(M_b))
       geom%M = M_b
    end if
    call mpi_bcast(geom%edge_cut_ni,   2*M_b, MPI_INTEGER,          0, partit%MPI_COMM_FESOM, ierr)
    call mpi_bcast(geom%edge_cut_lint, M_b,   MPI_DOUBLE_PRECISION, 0, partit%MPI_COMM_FESOM, ierr)
    call mpi_bcast(geom%edge_cut_midP, 2*M_b, MPI_DOUBLE_PRECISION, 0, partit%MPI_COMM_FESOM, ierr)
    call mpi_bcast(geom%edge_cut_dist, M_b,   MPI_DOUBLE_PRECISION, 0, partit%MPI_COMM_FESOM, ierr)
  end subroutine broadcast_transect_min


  !> For each crossed edge m, record this rank's lid for n1 and n2 if it
  !> strict-owns them, else 0. Strict-own = present in
  !> myList_nod2D(1:myDim_nod2D) only — never the halo range.
  !> Implementation: build a global->local lookup once (size nod2D), use
  !> for all m's, drop. Memory: 4*nod2D bytes per rank — cheap.
  subroutine build_sampling_lookups(mesh, partit, t)
    type(t_mesh),   intent(in)    :: mesh
    type(t_partit), intent(in)    :: partit
    type(track_t),  intent(inout) :: t
    integer, allocatable :: g2l(:)
    integer :: m, g, k

    allocate(g2l(mesh%nod2D))
    g2l = 0
    do k = 1, partit%myDim_nod2D
       g2l(partit%myList_nod2D(k)) = k
    end do

    allocate(t%lid_n1(t%M), t%lid_n2(t%M))
    do m = 1, t%M
       g = t%edge_cut_ni(1, m)
       if (g >= 1 .and. g <= mesh%nod2D) then
          t%lid_n1(m) = g2l(g)
       else
          t%lid_n1(m) = 0
       end if
       g = t%edge_cut_ni(2, m)
       if (g >= 1 .and. g <= mesh%nod2D) then
          t%lid_n2(m) = g2l(g)
       else
          t%lid_n2(m) = 0
       end if
    end do
    deallocate(g2l)
  end subroutine build_sampling_lookups


  !> Round-robin assignment of global edge indices m=1..M to ranks for
  !> XIOS write purposes. Each rank's nL = ceil((M - mype) / npes).
  subroutine assign_xios_cells(partit, t)
    type(t_partit), intent(in)    :: partit
    type(track_t),  intent(inout) :: t
    integer :: m, k, nL

    nL = 0
    do m = 1, t%M
       if (mod(m - 1, partit%npes) == partit%mype) nL = nL + 1
    end do
    t%nL = nL
    allocate(t%my_m(nL), t%i_index(nL))
    k = 0
    do m = 1, t%M
       if (mod(m - 1, partit%npes) == partit%mype) then
          k = k + 1
          t%my_m(k)   = m
          t%i_index(k) = m - 1     ! XIOS expects 0-based
       end if
    end do
  end subroutine assign_xios_cells


  !> Sample one track: per-rank contributions to (nz_track x M), MPI_Allreduce,
  !> mask boundary/dry-bottom slots as NaN, push this rank's stride to XIOS.
  subroutine send_one_track(tracers, mesh, partit, t)
    type(t_tracer), intent(in) :: tracers
    type(t_mesh),   intent(in) :: mesh
    type(t_partit), intent(in) :: partit
    type(track_t),  intent(in) :: t
    real(WP), allocatable :: contrib(:, :), wsum(:, :), sampled(:, :), wsum_tot(:, :), buf(:, :)
    integer  :: m, kk, ierr, nz_wet, k, lid
    real(WP) :: tt, nan, w

    if (len_trim(t%field_id) == 0)                            return
    if (.not. xios_is_valid_field(trim(t%field_id)))          return
    if (.not. xios_field_is_active(trim(t%field_id), .true.)) return
    if (t%M == 0)                                             return

    nan = ieee_value(0.0_WP, ieee_quiet_nan)

    ! Sum both the weighted contributions AND the contributed weights so
    ! we can rescale at the end. This makes wet/dry edges fall back to
    ! the surviving endpoint's value rather than the biased partial sum.
    !   both wet:  sampled = V1*(1-t)+V2*t   (weight sum = 1)
    !   only n1:   sampled = V1               (weight sum = 1-t, value = V1*(1-t)/(1-t))
    !   only n2:   sampled = V2
    !   both dry:  wsum=0 -> NaN
    allocate(contrib(nz_track, t%M), wsum(nz_track, t%M))
    contrib = 0.0_WP
    wsum    = 0.0_WP

    do m = 1, t%M
       tt = t%edge_cut_lint(m)

       lid = t%lid_n1(m)
       if (lid > 0) then
          if (lid > partit%myDim_nod2D) call mpi_abort(partit%MPI_COMM_FESOM, 1, ierr)
          nz_wet = max(0, mesh%nlevels_nod2D(lid) - 1)
          w = 1.0_WP - tt
          do kk = 1, min(nz_track, nz_wet)
             contrib(kk, m) = contrib(kk, m) + &
                  real(tracers%data(t%tracer_idx)%values(kk, lid), WP) * w
             wsum   (kk, m) = wsum   (kk, m) + w
          end do
       end if

       lid = t%lid_n2(m)
       if (lid > 0) then
          if (lid > partit%myDim_nod2D) call mpi_abort(partit%MPI_COMM_FESOM, 1, ierr)
          nz_wet = max(0, mesh%nlevels_nod2D(lid) - 1)
          w = tt
          do kk = 1, min(nz_track, nz_wet)
             contrib(kk, m) = contrib(kk, m) + &
                  real(tracers%data(t%tracer_idx)%values(kk, lid), WP) * w
             wsum   (kk, m) = wsum   (kk, m) + w
          end do
       end if
    end do

    allocate(sampled(nz_track, t%M), wsum_tot(nz_track, t%M))
    call MPI_Allreduce(contrib, sampled,  nz_track*t%M, MPI_DOUBLE_PRECISION, &
                       MPI_SUM, partit%MPI_COMM_FESOM, ierr)
    call MPI_Allreduce(wsum,    wsum_tot, nz_track*t%M, MPI_DOUBLE_PRECISION, &
                       MPI_SUM, partit%MPI_COMM_FESOM, ierr)
    deallocate(contrib, wsum)

    ! Rescale by the actual contributed weight, NaN where no endpoint contributed.
    do m = 1, t%M
       do kk = 1, nz_track
          if (wsum_tot(kk, m) > 0.0_WP) then
             sampled(kk, m) = sampled(kk, m) / wsum_tot(kk, m)
          else
             sampled(kk, m) = nan
          end if
       end do
    end do
    deallocate(wsum_tot)

    ! Push this rank's stride of m's.
    allocate(buf(nz_track, t%nL))
    do k = 1, t%nL
       do kk = 1, nz_track
          buf(kk, k) = sampled(kk, t%my_m(k))
       end do
    end do
    call xios_send_field(trim(t%field_id), buf)
    deallocate(buf, sampled)
  end subroutine send_one_track


  ! ====================================================================
  ! Config parsing (unchanged from v1 modulo track_resolution_km drop)
  ! ====================================================================

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
       deallocate(vars)
       allocate(vars(n_tracks))
       vars(:) = trim(track_vars)
       nv = n_tracks
    end if
    if (nv /= n_tracks .and. partit%mype == 0) then
       write(*,'(a,i0,a,i0)') '[TRACKS] ERROR: track_vars has ', nv,         &
            ' entries but track_files has ', n_tracks
    end if

    if (nn == 0) then
       allocate(names(n_tracks))
       do i = 1, n_tracks
          write(tag, '(i0)') i
          names(i) = 'track' // trim(tag)
       end do
    else if (nn /= n_tracks .and. partit%mype == 0) then
       write(*,'(a,i0,a,i0)') '[TRACKS] ERROR: track_names has ', nn,        &
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
          write(*,'(a,i0,a,a,a,a,a,a,a)') '[TRACKS]   ', i, '. name="',      &
               trim(tracks(i)%name), '" var="', trim(tracks(i)%var_name),    &
               '" csv="', trim(tracks(i)%csv_path), '"'
       end do
       write(*,'(a,a)') '[TRACKS] shared cadence=', trim(track_output_freq)
    end if
  end subroutine parse_track_config


  subroutine split_semi_or_comma(s_in, n, tokens)
    character(len=*),                intent(in)  :: s_in
    integer,                         intent(out) :: n
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


  function parse_xios_freq(s) result(d)
    character(len=*), intent(in) :: s
    type(xios_duration) :: d
    character(len=16) :: tok
    integer :: n, idx, ios
    tok = adjustl(s)
    idx = index(tok, 'mo')
    if (idx > 0) then
       read(tok(1:idx-1), *, iostat=ios) n
       if (ios == 0) then; d = xios_duration(month=real(n, 8)); return; end if
    end if
    idx = index(tok, 'y')
    if (idx > 0) then
       read(tok(1:idx-1), *, iostat=ios) n
       if (ios == 0) then; d = xios_duration(year=real(n, 8)); return; end if
    end if
    idx = index(tok, 'd')
    if (idx > 0) then
       read(tok(1:idx-1), *, iostat=ios) n
       if (ios == 0) then; d = xios_duration(day=real(n, 8)); return; end if
    end if
    idx = index(tok, 'h')
    if (idx > 0) then
       read(tok(1:idx-1), *, iostat=ios) n
       if (ios == 0) then; d = xios_duration(hour=real(n, 8)); return; end if
    end if
    d = xios_duration(hour=1.0_8)
  end function parse_xios_freq


  ! ====================================================================
  ! XIOS runtime registration
  ! ====================================================================

  subroutine register_track_filegroup(partit)
    type(t_partit), intent(in) :: partit
    type(xios_filegroup) :: filgrp_root, my_filgrp

    call xios_get_handle("file_definition", filgrp_root)
    call xios_add_child(filgrp_root, my_filgrp, "fesom_tracks")
    call xios_set_attr(my_filgrp,                                              &
         output_freq = parse_xios_freq(trim(track_output_freq)),                &
         split_freq  = xios_duration(year = 1.0_8),                             &
         enabled     = .true.)
    if (partit%mype == 0) then
       write(*,'(a,a,a)') '[TRACKS] file_group "fesom_tracks" created with output_freq=', &
                          trim(track_output_freq), ' split_freq=1y'
    end if
  end subroutine register_track_filegroup


  subroutine register_xios_objects_for_track(partit, t)
    type(t_partit), intent(in) :: partit
    type(track_t),  intent(in) :: t

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
    character(len=64)      :: domain_id, grid_id, file_id, file_name
    real(WP), allocatable  :: lon_loc(:), lat_loc(:)
    integer :: k

    freq = parse_xios_freq(trim(track_output_freq))

    domain_id = 'track_' // trim(t%name)
    grid_id   = 'grid_3d_track_' // trim(t%name)
    file_id   = 'f_' // trim(t%field_id)
    file_name = trim(t%field_id) // '.fesom'

    allocate(lon_loc(t%nL), lat_loc(t%nL))
    do k = 1, t%nL
       lon_loc(k) = t%lon_mid(t%my_m(k))
       lat_loc(k) = t%lat_mid(t%my_m(k))
    end do

    call xios_get_handle("domain_definition", dgrp)
    call xios_add_child(dgrp, dom, trim(domain_id))
    call xios_set_attr(dom,                              &
         type        = "unstructured",                    &
         ni_glo      = t%M,                                &
         ni          = t%nL,                               &
         i_index     = t%i_index,                          &
         lonvalue_1d = lon_loc,                            &
         latvalue_1d = lat_loc,                            &
         data_dim    = 1,                                  &
         data_ni     = t%nL)

    call xios_get_handle("grid_definition", ggrp)
    call xios_add_child(ggrp, grd, trim(grid_id))
    call xios_add_axistogrid(grd, ax_in_grid)
    call xios_set_attr(ax_in_grid, axis_ref = "nz")
    call xios_add_domaintogrid(grd, dom_in_grid)
    call xios_set_attr(dom_in_grid, domain_ref = trim(domain_id))

    call xios_get_handle("field_definition", fgrp)
    call xios_add_child(fgrp, fld, trim(t%field_id))
    call xios_set_attr(fld,                              &
         grid_ref  = trim(grid_id),                       &
         operation = "average",                           &
         prec      = 8,                                   &
         long_name = "ship-track/mooring curtain (geometry-aware)")

    call xios_get_handle("fesom_tracks", my_filgrp)
    call xios_add_child(my_filgrp, fil, trim(file_id))
    call xios_set_attr(fil, name = trim(file_name))
    call xios_add_child(fil, fld_in_file)
    call xios_set_attr(fld_in_file,                                            &
         field_ref = trim(t%field_id),                                          &
         operation = "instant",                                                 &
         freq_op   = freq)

    call xios_solve_inheritance()
    deallocate(lon_loc, lat_loc)
  end subroutine register_xios_objects_for_track


  ! ====================================================================
  ! CSV parser (unchanged)
  ! ====================================================================

  subroutine parse_csv(path, N, lon_out, lat_out, partit)
    character(len=*),              intent(in)  :: path
    integer,                       intent(out) :: N
    real(WP), allocatable,         intent(out) :: lon_out(:), lat_out(:)
    type(t_partit),                intent(in)  :: partit

    integer, parameter :: MAX_WPTS = 100000
    real(WP) :: lon_buf(MAX_WPTS), lat_buf(MAX_WPTS)
    character(len=512) :: line
    integer :: u, ios, line_no, k
    integer :: comma1, ierr
    character(len=128) :: ftok, stok
    real(WP) :: lon, lat

    open(newunit=u, file=trim(path), status='old', action='read', iostat=ios)
    if (ios /= 0) then
       write(*,'(a,a,a,i0)') '[TRACKS] ERROR: cannot open track file "',     &
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
          write(*,'(a,a,a,i0)') '[TRACKS] ERROR: ', trim(path),               &
                                 ' exceeds MAX_WPTS=', MAX_WPTS
          call mpi_abort(partit%MPI_COMM_FESOM, 1, ierr)
       end if

       do k = 1, N
          if (abs(lon_buf(k) - lon) < 1.0e-12_WP .and. &
              abs(lat_buf(k) - lat) < 1.0e-12_WP) then
             write(*,'(a,a,a,i0,a,i0)') '[TRACKS] ERROR: ', trim(path),       &
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

#endif

end module io_tracks_module
