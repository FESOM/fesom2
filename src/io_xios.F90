!> XIOS output integration for FESOM-2 in coupled AWI-ESM3 (OASIS + dedicated
!> xios_server.exe ranks). Follows the NEMO/OIFS pattern: OASIS has already
!> split MPI_COMM_WORLD into per-component local comms; FESOM hands its OASIS
!> local comm (partit%MPI_COMM_FESOM) to XIOS via local_comm=, and the separate
!> xios_server.exe binary registers with OASIS independently.
!>
!> Lifecycle (cf. arpifs/module/suxios.F90 in OpenIFS 48r1, NEMO 4.2 nemogcm):
!>   xios_initialize("fesom", local_comm=parent_comm)    ! no return_comm split
!>   xios_context_initialize("fesom", parent_comm)
!>   xios_set_domain_attr("nodes",    ...)               ! fills XML blanks
!>   xios_set_domain_attr("elements", ...)
!>   xios_set_axis_attr  ("nz",       ...)
!>   xios_close_context_definition()
!>
!> Ordering: OASIS3-MCT is initialised in fesom_module.F90 before par_init /
!> mesh_setup, therefore before this routine is invoked. Matches NEMO's safe
!> EC-Earth ordering (OASIS first, then XIOS on the OASIS local comm).
module io_xios_module
#if defined(__XIOS)
  use xios
  use mod_mesh,    only: T_MESH
  use mod_partit,  only: T_PARTIT
  use o_param,     only: WP, rad
  use g_config,    only: dt, use_cavity
  use g_clock,     only: yearnew, daynew, timenew
  use diagnostics,          only: ldiag_solver, lcurt_stress_surf,       &
                                  ldiag_curl_vel3, ldiag_Ri,             &
                                  ldiag_TurbFlux, ldiag_salt3D,          &
                                  ldiag_dMOC, ldiag_DVD, ldiag_forc,     &
                                  ldiag_extflds, ldiag_destine,          &
                                  ldiag_ice, ldiag_trflx,                &
                                  ldiag_uvw_sqr, ldiag_trgrd_xyz
  use cmor_variables_diag,  only: ldiag_cmor
  implicit none
  private

  public :: io_xios_init, io_xios_close
  public :: io_xios_update_calendar
  public :: io_xios_send_2d_r8, io_xios_send_3d_r8
  public :: io_xios_send_2d_r4, io_xios_send_3d_r4
  public :: io_xios_is_on
  public :: io_xios_set_ice_conc, io_xios_is_ice_field
  public :: io_xios_apply_ice_mask_2d_r4, io_xios_apply_ice_mask_2d_r8
  public :: io_xios_apply_ice_mask_2d_elem_r4, io_xios_apply_ice_mask_2d_elem_r8

  ! NetCDF default _FillValue constants (used for sender-side ice masking;
  ! XIOS must also have detect_missing_value=true on the corresponding field
  ! group so these slots are excluded from the running mean).
  real(kind=8), parameter :: NC_FILL_DOUBLE = 9.9692099683868690e+36_8
  real(kind=4), parameter :: NC_FILL_FLOAT  = 9.9692099683868690e+36_4
  real(kind=8), parameter :: ICE_CONC_EPS   = 1.0e-6_8

  ! Pointer to ice-area-fraction array (ice%data(1)%values). Registered by
  ! ini_mean_io once MOD_ICE is in scope.
  real(kind=WP), pointer, save :: p_ice_conc(:) => null()

  ! Pointer to elem2D_nodes(1:3, 1:elem2D_local) — used to map node-based
  ! a_ice onto element fields (element is ice-covered iff any of its 3
  ! vertex nodes has a_ice >= eps).
  integer, pointer, save :: p_elem2D_nodes(:,:) => null()

  type(xios_context), save :: ctx_hdl
  logical,            save :: xios_on = .false.

  ! Indices (1-based, into the local myDim_elem2D list) of elements owned by
  ! this rank — i.e. those sent to XIOS. Populated in io_xios_init.
  integer,     allocatable, save, target :: owned_elem_local(:)
  integer,                  save :: n_owned_elem = 0

  public :: io_xios_owned_elem_local, io_xios_n_owned_elem

contains

  function io_xios_owned_elem_local() result(p)
    integer, pointer :: p(:)
    p => null()
    if (allocated(owned_elem_local)) p => owned_elem_local
  end function

  integer function io_xios_n_owned_elem() result(n)
    n = n_owned_elem
  end function

  logical function io_xios_is_on() result(r)
    r = xios_on
  end function


  !> Initialise the XIOS client side.
  !> parent_comm:  the FESOM OASIS local communicator (partit%MPI_COMM_FESOM).
  !> client_comm:  returned equal to parent_comm (kept for API compat; OASIS
  !>               has already done the world split, so no further split here).
  subroutine io_xios_init(mesh, partit, parent_comm, client_comm)
    type(T_MESH),    intent(in), target :: mesh
    type(T_PARTIT),  intent(in), target :: partit
    integer,         intent(in)         :: parent_comm
    integer,         intent(out)        :: client_comm

    integer                       :: i, e, n1, n2, n3, nn, ne, ne_owned, nz_cell
    real(kind=8), allocatable     :: lon_n(:),  lat_n(:)
    real(kind=8), allocatable     :: lon_e(:),  lat_e(:)
    integer,      allocatable     :: i_index_n(:), i_index_e(:)
    real(kind=8), allocatable     :: z_mid(:)
    real(kind=8)                  :: xlon, ylat

    nn = partit%myDim_nod2D
    ne = partit%myDim_elem2D

    ! --- 1. initialise XIOS client on the FESOM OASIS local comm ------------
    call xios_initialize("fesom", local_comm=parent_comm)
    client_comm = parent_comm    ! OASIS already split; keep API compat

    ! --- 2. open the context on the FESOM communicator -----------------------
    call xios_context_initialize("fesom", parent_comm)
    call xios_get_handle("fesom", ctx_hdl)
    call xios_set_current_context(ctx_hdl)

    ! --- 2b. optional XML override of FESOM diagnostic flags (OIFS pattern).
    ! Must run before init_cmor_diag / compute_diag(mode=0) allocates arrays.
    ! Each flag retains its namelist value if absent from the XML.
    block
      logical :: ov
      ov = xios_getvar("ldiag_solver",      ldiag_solver)
      if (ov .and. partit%mype==0) write(*,*) '[XIOS] ldiag_solver=',      ldiag_solver
      ov = xios_getvar("lcurt_stress_surf", lcurt_stress_surf)
      if (ov .and. partit%mype==0) write(*,*) '[XIOS] lcurt_stress_surf=', lcurt_stress_surf
      ov = xios_getvar("ldiag_curl_vel3",   ldiag_curl_vel3)
      if (ov .and. partit%mype==0) write(*,*) '[XIOS] ldiag_curl_vel3=',   ldiag_curl_vel3
      ov = xios_getvar("ldiag_Ri",          ldiag_Ri)
      if (ov .and. partit%mype==0) write(*,*) '[XIOS] ldiag_Ri=',          ldiag_Ri
      ov = xios_getvar("ldiag_TurbFlux",    ldiag_TurbFlux)
      if (ov .and. partit%mype==0) write(*,*) '[XIOS] ldiag_TurbFlux=',    ldiag_TurbFlux
      ov = xios_getvar("ldiag_salt3D",      ldiag_salt3D)
      if (ov .and. partit%mype==0) write(*,*) '[XIOS] ldiag_salt3D=',      ldiag_salt3D
      ov = xios_getvar("ldiag_dMOC",        ldiag_dMOC)
      if (ov .and. partit%mype==0) write(*,*) '[XIOS] ldiag_dMOC=',        ldiag_dMOC
      ov = xios_getvar("ldiag_DVD",         ldiag_DVD)
      if (ov .and. partit%mype==0) write(*,*) '[XIOS] ldiag_DVD=',         ldiag_DVD
      ov = xios_getvar("ldiag_forc",        ldiag_forc)
      if (ov .and. partit%mype==0) write(*,*) '[XIOS] ldiag_forc=',        ldiag_forc
      ov = xios_getvar("ldiag_extflds",     ldiag_extflds)
      if (ov .and. partit%mype==0) write(*,*) '[XIOS] ldiag_extflds=',     ldiag_extflds
      ov = xios_getvar("ldiag_destine",     ldiag_destine)
      if (ov .and. partit%mype==0) write(*,*) '[XIOS] ldiag_destine=',     ldiag_destine
      ov = xios_getvar("ldiag_ice",         ldiag_ice)
      if (ov .and. partit%mype==0) write(*,*) '[XIOS] ldiag_ice=',         ldiag_ice
      ov = xios_getvar("ldiag_trflx",       ldiag_trflx)
      if (ov .and. partit%mype==0) write(*,*) '[XIOS] ldiag_trflx=',       ldiag_trflx
      ov = xios_getvar("ldiag_uvw_sqr",     ldiag_uvw_sqr)
      if (ov .and. partit%mype==0) write(*,*) '[XIOS] ldiag_uvw_sqr=',     ldiag_uvw_sqr
      ov = xios_getvar("ldiag_trgrd_xyz",   ldiag_trgrd_xyz)
      if (ov .and. partit%mype==0) write(*,*) '[XIOS] ldiag_trgrd_xyz=',   ldiag_trgrd_xyz
      ov = xios_getvar("ldiag_cmor",        ldiag_cmor)
      if (ov .and. partit%mype==0) write(*,*) '[XIOS] ldiag_cmor=',        ldiag_cmor
    end block

    ! --- 3. calendar: left to XML (calendar_type on the context element) -----

    ! --- 4. unstructured node-based domain "nodes" ---------------------------
    allocate(lon_n(nn), lat_n(nn), i_index_n(nn))
    do i = 1, nn
       lon_n(i)    = mesh%geo_coord_nod2D(1, i) / rad   ! rad -> deg
       lat_n(i)    = mesh%geo_coord_nod2D(2, i) / rad
       i_index_n(i)= partit%myList_nod2D(i) - 1         ! 1-based -> 0-based
    end do

    call xios_set_domain_attr("nodes",              &
         type        = "unstructured",              &
         ni_glo      = mesh%nod2D,                  &
         ni          = nn,                          &
         ibegin      = 0,                           &
         i_index     = i_index_n,                   &
         lonvalue_1d = lon_n,                       &
         latvalue_1d = lat_n,                       &
         data_dim    = 1,                           &
         data_ni     = nn)

    deallocate(lon_n, lat_n, i_index_n)

    ! --- 5. unstructured element (triangle-centre) domain "elements" --------
    ! FESOM's myList_elem2D overlaps across ranks (any element touching a node
    ! owned by this rank is in the list), which makes XIOS complain about
    ! duplicated i_index. Restrict to elements STRICTLY owned by this rank,
    ! where ownership = min(partit%part) over the element's 3 vertex nodes.
    ! All ranks agree on the same partit%part array, so each global element
    ! ends up claimed by exactly one rank.
    allocate(lon_e(ne), lat_e(ne), i_index_e(ne), owned_elem_local(ne))
    ne_owned = 0
    do e = 1, ne
       n1 = mesh%elem2D_nodes(1, e)
       n2 = mesh%elem2D_nodes(2, e)
       n3 = mesh%elem2D_nodes(3, e)
       if (min(partit%part(n1), partit%part(n2), partit%part(n3)) &
                                                    /= partit%mype) cycle
       ne_owned = ne_owned + 1
       owned_elem_local(ne_owned) = e
       xlon = ( mesh%geo_coord_nod2D(1, n1) &
              + mesh%geo_coord_nod2D(1, n2) &
              + mesh%geo_coord_nod2D(1, n3) ) / 3.0_WP
       ylat = ( mesh%geo_coord_nod2D(2, n1) &
              + mesh%geo_coord_nod2D(2, n2) &
              + mesh%geo_coord_nod2D(2, n3) ) / 3.0_WP
       lon_e(ne_owned)     = xlon / rad
       lat_e(ne_owned)     = ylat / rad
       i_index_e(ne_owned) = partit%myList_elem2D(e) - 1
    end do
    n_owned_elem = ne_owned

    ! Sort the strictly-owned elements by global index. XIOS server2 maps
    ! i_index -> server rank assuming roughly monotonic indices per client;
    ! unsorted scatter causes uneven server chunks (e.g. some servers ending
    ! up with a tiny handful of elements while the domain global array is
    ! redistributed as a whole -> writeDomain_ 'intern vs input' mismatch).
    call sort_owned_by_gid(i_index_e, lon_e, lat_e, owned_elem_local, ne_owned)

    ! Omit ibegin: per XIOS 2.5 domain.cpp:479-488, supplying i_index causes
    ! XIOS to set ibegin = i_index(0). Explicitly passing ibegin=0 alongside
    ! i_index confuses the server-side band distribution attribute check.
    call xios_set_domain_attr("elements",                   &
         type        = "unstructured",                      &
         ni_glo      = mesh%elem2D,                         &
         ni          = ne_owned,                            &
         i_index     = i_index_e(1:ne_owned),               &
         lonvalue_1d = lon_e(1:ne_owned),                   &
         latvalue_1d = lat_e(1:ne_owned),                   &
         data_dim    = 1,                                   &
         data_ni     = ne_owned)

    deallocate(lon_e, lat_e, i_index_e)

    ! --- 6. vertical axes ---------------------------------------------------
    !   nz  : cell-centred,  nl-1 layers (for T, S, u, v, unod, vnod, ...)
    !   nz1 : interface/w-grid, nl levels (for w, bolus_w, Kv, N2, ...)
    nz_cell = mesh%nl - 1
    allocate(z_mid(nz_cell))
    ! mesh%zbar is NEGATIVE (z up, surface=0). Legacy/CF convention expects
    ! depths positive-down, matching axis attribute positive="down".
    do i = 1, nz_cell
       z_mid(i) = -0.5_WP * (mesh%zbar(i) + mesh%zbar(i+1))
    end do
    call xios_set_axis_attr("nz",  n_glo = nz_cell, value = z_mid)
    deallocate(z_mid)
    call xios_set_axis_attr("nz1", n_glo = mesh%nl, value = -mesh%zbar(1:mesh%nl))

    ! --- 7. timestep (required by XIOS before close_context_definition) -----
    call xios_set_timestep(timestep = xios_duration(second = dt))

    ! --- 7b. calendar origin / start date (drives filename date suffix) -----
    ! Without these the split suffix comes out as "0000-0000". Matches the
    ! OIFS suxios.F90 pattern. time_origin = Jan 1 of the run's start year,
    ! start_date = time_origin + (daynew-1) days + timenew seconds.
    call xios_set_time_origin(time_origin = &
         xios_date(yearnew, 1, 1, 0, 0, 0))
    call xios_set_start_date(start_date =   &
         xios_date(yearnew, 1, 1, 0, 0, 0) + &
         xios_duration(day = real(daynew - 1, kind=8), &
                       second = real(timenew, kind=8)))

    ! --- 7c. masks for land / below-bottom / cavity --------------------------
    ! XIOS averages uninitialised buffer slots as if they were valid, which
    ! pollutes land / below-bottom / ice-shelf points with zeros or garbage.
    ! Declare the valid-point masks here so XIOS writes _FillValue to the
    ! invalid slots instead.
    !
    ! Non-cavity run:
    !   - mesh has no land nodes/elements, so 2D masks collapse to all-true.
    !   - 3D masks select k in [1 .. nlevels-1] (or 1 .. nlevels for interfaces).
    ! Cavity run (use_cavity=.true.):
    !   - ulevels_* can be > 1 → some 2D nodes/elements are DRY at surface.
    !   - 3D mask lower bound is ulevels_*, not 1.
    call io_xios_set_masks(mesh, partit, ne_owned)

    ! --- 8. close context definition ----------------------------------------
    call xios_close_context_definition()

    xios_on = .true.
  end subroutine io_xios_init


  !> Called every FESOM timestep from the main loop.
  subroutine io_xios_update_calendar(step)
    integer, intent(in) :: step
    if (.not. xios_on) return
    call xios_set_current_context(ctx_hdl)
    call xios_update_calendar(step)
  end subroutine


  !> Send a 2D (node-only) field. Shape: (ni=myDim_nod2D).
  !> Caller (io_meandata.F90) must already have sliced to owned nodes.
  subroutine io_xios_send_2d_r8(name, buf)
    character(len=*), intent(in) :: name
    real(kind=8),     intent(in) :: buf(:)
    logical :: is_valid
    if (.not. xios_on) return
    is_valid = xios_is_valid_field(trim(name))
    if (.not. is_valid) return
    call xios_send_field(trim(name), buf)
  end subroutine

  subroutine io_xios_send_2d_r4(name, buf)
    character(len=*), intent(in) :: name
    real(kind=4),     intent(in) :: buf(:)
    logical :: is_valid
    if (.not. xios_on) return
    is_valid = xios_is_valid_field(trim(name))
    if (.not. is_valid) return
    call xios_send_field(trim(name), buf)
  end subroutine


  !> Send a 3D (node x level) field. Shape: (nz_cell, myDim_nod2D).
  subroutine io_xios_send_3d_r8(name, buf)
    character(len=*), intent(in) :: name
    real(kind=8),     intent(in) :: buf(:,:)
    logical :: is_valid
    if (.not. xios_on) return
    is_valid = xios_is_valid_field(trim(name))
    if (.not. is_valid) return
    call xios_send_field(trim(name), buf)
  end subroutine

  subroutine io_xios_send_3d_r4(name, buf)
    character(len=*), intent(in) :: name
    real(kind=4),     intent(in) :: buf(:,:)
    logical :: is_valid
    if (.not. xios_on) return
    is_valid = xios_is_valid_field(trim(name))
    if (.not. is_valid) return
    call xios_send_field(trim(name), buf)
  end subroutine


  !> Build and push wet/bottom masks to the FESOM domains and grids.
  !> Must be called between xios_set_domain_attr("nodes"/"elements") and
  !> xios_close_context_definition().
  subroutine io_xios_set_masks(mesh, partit, ne_owned)
    type(T_MESH),   intent(in), target :: mesh
    type(T_PARTIT), intent(in), target :: partit
    integer,        intent(in)         :: ne_owned

    integer              :: nn, nz_cell, nl, n, e, k, ee
    logical, allocatable :: m1d_n(:), m1d_e(:)
    logical, allocatable :: m3d_n_nz(:,:),   m3d_n_nz1(:,:)
    logical, allocatable :: m3d_e_nz(:,:),   m3d_e_nz1(:,:)
    integer              :: u_n, b_n, u_e, b_e

    nn      = partit%myDim_nod2D
    nz_cell = mesh%nl - 1
    nl      = mesh%nl

    ! -- 2D node mask --
    allocate(m1d_n(nn))
    if (use_cavity) then
       do n = 1, nn
          m1d_n(n) = (mesh%ulevels_nod2D(n) == 1)
       end do
    else
       m1d_n(:) = .true.
    end if
    call xios_set_domain_attr("nodes", mask_1d = m1d_n)

    ! -- 2D element mask (strictly-owned subset) --
    allocate(m1d_e(ne_owned))
    if (use_cavity) then
       do ee = 1, ne_owned
          e = owned_elem_local(ee)
          m1d_e(ee) = (mesh%ulevels(e) == 1)
       end do
    else
       m1d_e(:) = .true.
    end if
    call xios_set_domain_attr("elements", mask_1d = m1d_e)

    ! -- 3D node masks, axis-first (nz, nn) to match send buffer layout --
    allocate(m3d_n_nz (nz_cell, nn))
    allocate(m3d_n_nz1(nl,      nn))
    m3d_n_nz (:,:) = .false.
    m3d_n_nz1(:,:) = .false.
    do n = 1, nn
       u_n = mesh%ulevels_nod2D(n)
       b_n = mesh%nlevels_nod2D(n)
       do k = u_n, min(b_n - 1, nz_cell)
          m3d_n_nz(k, n) = .true.
       end do
       do k = u_n, min(b_n, nl)
          m3d_n_nz1(k, n) = .true.
       end do
    end do
    call xios_set_grid_attr("grid_3d_nod",     mask_2d = m3d_n_nz)
    call xios_set_grid_attr("grid_3d_nod_nz1", mask_2d = m3d_n_nz1)

    ! -- 3D element masks, axis-first (nz, ne_owned) --
    allocate(m3d_e_nz (nz_cell, ne_owned))
    allocate(m3d_e_nz1(nl,      ne_owned))
    m3d_e_nz (:,:) = .false.
    m3d_e_nz1(:,:) = .false.
    do ee = 1, ne_owned
       e   = owned_elem_local(ee)
       u_e = mesh%ulevels(e)
       b_e = mesh%nlevels(e)
       do k = u_e, min(b_e - 1, nz_cell)
          m3d_e_nz(k, ee) = .true.
       end do
       do k = u_e, min(b_e, nl)
          m3d_e_nz1(k, ee) = .true.
       end do
    end do
    call xios_set_grid_attr("grid_3d_elem",     mask_2d = m3d_e_nz)
    call xios_set_grid_attr("grid_3d_elem_nz1", mask_2d = m3d_e_nz1)

    deallocate(m1d_n, m1d_e, m3d_n_nz, m3d_n_nz1, m3d_e_nz, m3d_e_nz1)
  end subroutine io_xios_set_masks


  !> Insertion sort owned elements by global index (i_index), carrying
  !> lon/lat/local-id in lock-step. ne_owned is typically ~1000-3000 per rank,
  !> so O(n^2) is fine here.
  subroutine sort_owned_by_gid(idx, lon, lat, loc, n)
    integer, intent(inout) :: idx(:), loc(:)
    real(kind=8), intent(inout) :: lon(:), lat(:)
    integer, intent(in) :: n
    integer :: i, j, k_idx, k_loc
    real(kind=8) :: k_lon, k_lat
    do i = 2, n
       k_idx = idx(i); k_lon = lon(i); k_lat = lat(i); k_loc = loc(i)
       j = i - 1
       do while (j >= 1)
          if (idx(j) <= k_idx) exit
          idx(j+1) = idx(j); lon(j+1) = lon(j); lat(j+1) = lat(j); loc(j+1) = loc(j)
          j = j - 1
       end do
       idx(j+1) = k_idx; lon(j+1) = k_lon; lat(j+1) = k_lat; loc(j+1) = k_loc
    end do
  end subroutine


  !> Called once from ini_mean_io after MOD_ICE is in scope.
  subroutine io_xios_set_ice_conc(p, elem_nodes)
    real(kind=WP), target, intent(in) :: p(:)
    integer,       target, intent(in) :: elem_nodes(:,:)
    p_ice_conc     => p
    p_elem2D_nodes => elem_nodes
  end subroutine

  !> Ice-tagged field names: these get sender-side masking where a_ice is
  !> below ICE_CONC_EPS. Keep in sync with field_def_fesom.xml ice subgroups.
  logical function io_xios_is_ice_field(name) result(r)
    character(len=*), intent(in) :: name
    select case (trim(name))
    case ('a_ice', 'm_ice', 'm_snow', 'h_ice', 'h_snow', 'ist', &
          'uice', 'vice', 'apnd', 'hpnd', 'ipnd', &
          'thdgrice', 'thdgrsnw', 'thdgrarea', 'dyngrice', 'dyngrarea', &
          'fw_ice', 'fw_snw', 'atmice_x', 'atmice_y', &
          'iceoce_x', 'iceoce_y', 'strength_ice', 'sgm11', 'sgm12', 'sgm22', &
          'qcon')
       r = .true.
    case default
       r = .false.
    end select
  end function

  subroutine io_xios_apply_ice_mask_2d_r8(buf)
    real(kind=8), intent(inout) :: buf(:)
    integer :: i, n
    if (.not. associated(p_ice_conc)) return
    n = min(size(buf), size(p_ice_conc))
    do i = 1, n
       if (real(p_ice_conc(i), kind=8) < ICE_CONC_EPS) buf(i) = NC_FILL_DOUBLE
    end do
  end subroutine

  subroutine io_xios_apply_ice_mask_2d_r4(buf)
    real(kind=4), intent(inout) :: buf(:)
    integer :: i, n
    if (.not. associated(p_ice_conc)) return
    n = min(size(buf), size(p_ice_conc))
    do i = 1, n
       if (real(p_ice_conc(i), kind=4) < real(ICE_CONC_EPS, kind=4)) buf(i) = NC_FILL_FLOAT
    end do
  end subroutine

  !> Element-based ice mask: element is ice-covered iff any of its 3 vertex
  !> nodes has a_ice >= eps. buf is indexed by owned-element order; local
  !> element indices come from owned_elem_local.
  subroutine io_xios_apply_ice_mask_2d_elem_r8(buf)
    real(kind=8), intent(inout) :: buf(:)
    integer :: ee, e, n1, n2, n3
    real(kind=8) :: amax
    if (.not. associated(p_ice_conc) .or. .not. associated(p_elem2D_nodes)) return
    if (.not. allocated(owned_elem_local)) return
    do ee = 1, min(size(buf), n_owned_elem)
       e  = owned_elem_local(ee)
       n1 = p_elem2D_nodes(1, e)
       n2 = p_elem2D_nodes(2, e)
       n3 = p_elem2D_nodes(3, e)
       amax = max(real(p_ice_conc(n1), kind=8), &
                  real(p_ice_conc(n2), kind=8), &
                  real(p_ice_conc(n3), kind=8))
       if (amax < ICE_CONC_EPS) buf(ee) = NC_FILL_DOUBLE
    end do
  end subroutine

  subroutine io_xios_apply_ice_mask_2d_elem_r4(buf)
    real(kind=4), intent(inout) :: buf(:)
    integer :: ee, e, n1, n2, n3
    real(kind=4) :: amax
    if (.not. associated(p_ice_conc) .or. .not. associated(p_elem2D_nodes)) return
    if (.not. allocated(owned_elem_local)) return
    do ee = 1, min(size(buf), n_owned_elem)
       e  = owned_elem_local(ee)
       n1 = p_elem2D_nodes(1, e)
       n2 = p_elem2D_nodes(2, e)
       n3 = p_elem2D_nodes(3, e)
       amax = max(real(p_ice_conc(n1), kind=4), &
                  real(p_ice_conc(n2), kind=4), &
                  real(p_ice_conc(n3), kind=4))
       if (amax < real(ICE_CONC_EPS, kind=4)) buf(ee) = NC_FILL_FLOAT
    end do
  end subroutine


  subroutine io_xios_close()
    if (.not. xios_on) return
    call xios_context_finalize()
    call xios_finalize()
    xios_on = .false.
  end subroutine

#else
  ! Stub module when built without XIOS; keeps `use io_xios_module` legal.
  ! All symbols imported by io_meandata.F90 must exist here as no-ops so
  ! the standalone (non-XIOS) build links without modification.
  implicit none
  private
  public :: io_xios_is_on
  public :: io_xios_send_2d_r8, io_xios_send_3d_r8
  public :: io_xios_send_2d_r4, io_xios_send_3d_r4
  public :: io_xios_owned_elem_local, io_xios_n_owned_elem
  public :: io_xios_set_ice_conc, io_xios_is_ice_field
  public :: io_xios_apply_ice_mask_2d_r4, io_xios_apply_ice_mask_2d_r8
  public :: io_xios_apply_ice_mask_2d_elem_r4, io_xios_apply_ice_mask_2d_elem_r8
contains
  logical function io_xios_is_on() result(r)
    r = .false.
  end function

  subroutine io_xios_send_2d_r8(name, buf)
    character(len=*), intent(in) :: name
    real(kind=8),     intent(in) :: buf(:)
  end subroutine

  subroutine io_xios_send_2d_r4(name, buf)
    character(len=*), intent(in) :: name
    real(kind=4),     intent(in) :: buf(:)
  end subroutine

  subroutine io_xios_send_3d_r8(name, buf)
    character(len=*), intent(in) :: name
    real(kind=8),     intent(in) :: buf(:,:)
  end subroutine

  subroutine io_xios_send_3d_r4(name, buf)
    character(len=*), intent(in) :: name
    real(kind=4),     intent(in) :: buf(:,:)
  end subroutine

  function io_xios_owned_elem_local() result(p)
    integer, pointer :: p(:)
    p => null()
  end function

  integer function io_xios_n_owned_elem() result(n)
    n = 0
  end function

  subroutine io_xios_set_ice_conc(p, elem_nodes)
    real(kind=8), target, intent(in) :: p(:)
    integer,      target, intent(in) :: elem_nodes(:,:)
  end subroutine

  logical function io_xios_is_ice_field(name) result(r)
    character(len=*), intent(in) :: name
    r = .false.
  end function

  subroutine io_xios_apply_ice_mask_2d_r8(buf)
    real(kind=8), intent(inout) :: buf(:)
  end subroutine

  subroutine io_xios_apply_ice_mask_2d_r4(buf)
    real(kind=4), intent(inout) :: buf(:)
  end subroutine

  subroutine io_xios_apply_ice_mask_2d_elem_r8(buf)
    real(kind=8), intent(inout) :: buf(:)
  end subroutine

  subroutine io_xios_apply_ice_mask_2d_elem_r4(buf)
    real(kind=4), intent(inout) :: buf(:)
  end subroutine
#endif
end module io_xios_module
