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
  use mpi
  use mod_mesh,    only: T_MESH
  use mod_partit,  only: T_PARTIT
  use o_param,     only: WP, rad
  use g_config,    only: dt
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
  public :: io_xios_set_wet_ptrs
  public :: io_xios_apply_wet_2d_r4,      io_xios_apply_wet_2d_r8
  public :: io_xios_apply_wet_2d_elem_r4, io_xios_apply_wet_2d_elem_r8
  public :: io_xios_apply_wet_3d_r4,      io_xios_apply_wet_3d_r8
  public :: io_xios_apply_wet_3d_elem_r4, io_xios_apply_wet_3d_elem_r8

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

  ! Mesh pointers registered by ini_mean_io for sender-side wet/bottom
  ! masking. Written to NC_FILL before xios_send_field so XIOS averages
  ! with detect_missing_value="true" exclude dry slots. Node arrays are
  ! sized myDim_nod2D; element arrays are sized myDim_elem2D (the caller
  ! packs via owned_elem_local).
  integer, pointer, save :: p_ulevels_nod(:) => null()
  integer, pointer, save :: p_nlevels_nod(:) => null()
  integer, pointer, save :: p_ulevels_elem(:) => null()
  integer, pointer, save :: p_nlevels_elem(:) => null()
  integer,           save :: p_nl = 0

  type(xios_context), save :: ctx_hdl
  logical,            save :: xios_on = .false.

  ! Indices (1-based, into the local myDim_elem2D list) of elements owned by
  ! this rank — i.e. those sent to XIOS. Populated in io_xios_init.
  integer,     allocatable, save, target :: owned_elem_local(:)
  integer,                  save :: n_owned_elem = 0
  integer,                  save :: init_call_count = 0

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
    block
      integer :: wrank, wsize, prank, psize, ierr
      call MPI_Comm_rank(MPI_COMM_WORLD, wrank, ierr)
      call MPI_Comm_size(MPI_COMM_WORLD, wsize, ierr)
      call MPI_Comm_rank(parent_comm,     prank, ierr)
      call MPI_Comm_size(parent_comm,     psize, ierr)
      write(*,'("IO_XIOS_COMM mype=",I5," wrank=",I5,"/",I5," prank=",I5,"/",I5," calls=",I0)') &
            partit%mype, wrank, wsize, prank, psize, init_call_count + 1
      flush(6)
    end block
    init_call_count = init_call_count + 1

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
    ! Use myDim_elem2D_shrinked for strictly-unique element ownership across
    ! ranks (FESOM rule: element owned iff elem2D_nodes(1,e) <= myDim_nod2D;
    ! sum over ranks = elem2D_glo). myInd_elem2D_shrinked maps shrinked idx
    ! -> local element index in 1..myDim_elem2D.
    ne_owned = partit%myDim_elem2D_shrinked
    allocate(lon_e(ne_owned), lat_e(ne_owned), i_index_e(ne_owned), owned_elem_local(ne_owned))
    do i = 1, ne_owned
       e = partit%myInd_elem2D_shrinked(i)
       owned_elem_local(i) = e
       xlon = ( mesh%geo_coord_nod2D(1, mesh%elem2D_nodes(1, e)) &
              + mesh%geo_coord_nod2D(1, mesh%elem2D_nodes(2, e)) &
              + mesh%geo_coord_nod2D(1, mesh%elem2D_nodes(3, e)) ) / 3.0_WP
       ylat = ( mesh%geo_coord_nod2D(2, mesh%elem2D_nodes(1, e)) &
              + mesh%geo_coord_nod2D(2, mesh%elem2D_nodes(2, e)) &
              + mesh%geo_coord_nod2D(2, mesh%elem2D_nodes(3, e)) ) / 3.0_WP
       lon_e(i)     = xlon / rad
       lat_e(i)     = ylat / rad
       i_index_e(i) = partit%myList_elem2D(e) - 1
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


  !> Register mesh pointers for sender-side wet/bottom masking. Called
  !> once from ini_mean_io (MOD_MESH is in scope there).
  subroutine io_xios_set_wet_ptrs(ulevels_nod, nlevels_nod, &
                                  ulevels_elem, nlevels_elem, nl)
    integer, target, intent(in) :: ulevels_nod(:),  nlevels_nod(:)
    integer, target, intent(in) :: ulevels_elem(:), nlevels_elem(:)
    integer,         intent(in) :: nl
    p_ulevels_nod  => ulevels_nod
    p_nlevels_nod  => nlevels_nod
    p_ulevels_elem => ulevels_elem
    p_nlevels_elem => nlevels_elem
    p_nl            = nl
  end subroutine

  !> 2D node field: NC_FILL where ulevels_nod2D(n) > 1 (cavity only; no-op
  !> for non-cavity runs where all nodes are wet at surface).
  subroutine io_xios_apply_wet_2d_r8(buf)
    real(kind=8), intent(inout) :: buf(:)
    integer :: n
    if (.not. associated(p_ulevels_nod)) return
    do n = 1, min(size(buf), size(p_ulevels_nod))
       if (p_ulevels_nod(n) > 1) buf(n) = NC_FILL_DOUBLE
    end do
  end subroutine

  subroutine io_xios_apply_wet_2d_r4(buf)
    real(kind=4), intent(inout) :: buf(:)
    integer :: n
    if (.not. associated(p_ulevels_nod)) return
    do n = 1, min(size(buf), size(p_ulevels_nod))
       if (p_ulevels_nod(n) > 1) buf(n) = NC_FILL_FLOAT
    end do
  end subroutine

  !> 2D element field (strictly-owned subset): NC_FILL where any of the
  !> 3 vertex nodes has ulevels > 1 (cavity element).
  subroutine io_xios_apply_wet_2d_elem_r8(buf)
    real(kind=8), intent(inout) :: buf(:)
    integer :: ee, e
    if (.not. associated(p_ulevels_elem)) return
    if (.not. allocated(owned_elem_local)) return
    do ee = 1, min(size(buf), n_owned_elem)
       e = owned_elem_local(ee)
       if (p_ulevels_elem(e) > 1) buf(ee) = NC_FILL_DOUBLE
    end do
  end subroutine

  subroutine io_xios_apply_wet_2d_elem_r4(buf)
    real(kind=4), intent(inout) :: buf(:)
    integer :: ee, e
    if (.not. associated(p_ulevels_elem)) return
    if (.not. allocated(owned_elem_local)) return
    do ee = 1, min(size(buf), n_owned_elem)
       e = owned_elem_local(ee)
       if (p_ulevels_elem(e) > 1) buf(ee) = NC_FILL_FLOAT
    end do
  end subroutine

  !> 3D node field, axis-first (nz, nn). Valid range is
  !> ulevels_nod2D(n) <= k <= bound(n), where bound = nlevels-1 for mid-
  !> layer axis (nz), nlevels for interface axis (nz1). Caller shape tells
  !> us which: size(buf,1)==mesh%nl => interface, else => mid-layer.
  subroutine io_xios_apply_wet_3d_r8(buf)
    real(kind=8), intent(inout) :: buf(:,:)
    integer :: n, k, nz, nn, ub, un, bn
    logical :: is_interface
    if (.not. associated(p_ulevels_nod) .or. .not. associated(p_nlevels_nod)) return
    nz = size(buf, 1); nn = size(buf, 2)
    is_interface = (p_nl > 0 .and. nz == p_nl)
    do n = 1, min(nn, size(p_ulevels_nod))
       un = p_ulevels_nod(n)
       bn = p_nlevels_nod(n)
       ub = bn; if (.not. is_interface) ub = bn - 1
       do k = 1, nz
          if (k < un .or. k > ub) buf(k, n) = NC_FILL_DOUBLE
       end do
    end do
  end subroutine

  subroutine io_xios_apply_wet_3d_r4(buf)
    real(kind=4), intent(inout) :: buf(:,:)
    integer :: n, k, nz, nn, ub, un, bn
    logical :: is_interface
    if (.not. associated(p_ulevels_nod) .or. .not. associated(p_nlevels_nod)) return
    nz = size(buf, 1); nn = size(buf, 2)
    is_interface = (p_nl > 0 .and. nz == p_nl)
    do n = 1, min(nn, size(p_ulevels_nod))
       un = p_ulevels_nod(n)
       bn = p_nlevels_nod(n)
       ub = bn; if (.not. is_interface) ub = bn - 1
       do k = 1, nz
          if (k < un .or. k > ub) buf(k, n) = NC_FILL_FLOAT
       end do
    end do
  end subroutine

  !> 3D element field, axis-first (nz, ne_owned). Same bound logic as
  !> node variant, indexed via owned_elem_local.
  subroutine io_xios_apply_wet_3d_elem_r8(buf)
    real(kind=8), intent(inout) :: buf(:,:)
    integer :: ee, e, k, nz, ne, ub, ue, be
    logical :: is_interface
    if (.not. associated(p_ulevels_elem) .or. .not. associated(p_nlevels_elem)) return
    if (.not. allocated(owned_elem_local)) return
    nz = size(buf, 1); ne = size(buf, 2)
    is_interface = (p_nl > 0 .and. nz == p_nl)
    do ee = 1, min(ne, n_owned_elem)
       e = owned_elem_local(ee)
       ue = p_ulevels_elem(e)
       be = p_nlevels_elem(e)
       ub = be; if (.not. is_interface) ub = be - 1
       do k = 1, nz
          if (k < ue .or. k > ub) buf(k, ee) = NC_FILL_DOUBLE
       end do
    end do
  end subroutine

  subroutine io_xios_apply_wet_3d_elem_r4(buf)
    real(kind=4), intent(inout) :: buf(:,:)
    integer :: ee, e, k, nz, ne, ub, ue, be
    logical :: is_interface
    if (.not. associated(p_ulevels_elem) .or. .not. associated(p_nlevels_elem)) return
    if (.not. allocated(owned_elem_local)) return
    nz = size(buf, 1); ne = size(buf, 2)
    is_interface = (p_nl > 0 .and. nz == p_nl)
    do ee = 1, min(ne, n_owned_elem)
       e = owned_elem_local(ee)
       ue = p_ulevels_elem(e)
       be = p_nlevels_elem(e)
       ub = be; if (.not. is_interface) ub = be - 1
       do k = 1, nz
          if (k < ue .or. k > ub) buf(k, ee) = NC_FILL_FLOAT
       end do
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
  public :: io_xios_set_wet_ptrs
  public :: io_xios_apply_wet_2d_r4,      io_xios_apply_wet_2d_r8
  public :: io_xios_apply_wet_2d_elem_r4, io_xios_apply_wet_2d_elem_r8
  public :: io_xios_apply_wet_3d_r4,      io_xios_apply_wet_3d_r8
  public :: io_xios_apply_wet_3d_elem_r4, io_xios_apply_wet_3d_elem_r8
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

  subroutine io_xios_set_wet_ptrs(ulevels_nod, nlevels_nod, &
                                  ulevels_elem, nlevels_elem, nl)
    integer, target, intent(in) :: ulevels_nod(:),  nlevels_nod(:)
    integer, target, intent(in) :: ulevels_elem(:), nlevels_elem(:)
    integer,         intent(in) :: nl
  end subroutine

  subroutine io_xios_apply_wet_2d_r8(buf)
    real(kind=8), intent(inout) :: buf(:)
  end subroutine
  subroutine io_xios_apply_wet_2d_r4(buf)
    real(kind=4), intent(inout) :: buf(:)
  end subroutine
  subroutine io_xios_apply_wet_2d_elem_r8(buf)
    real(kind=8), intent(inout) :: buf(:)
  end subroutine
  subroutine io_xios_apply_wet_2d_elem_r4(buf)
    real(kind=4), intent(inout) :: buf(:)
  end subroutine
  subroutine io_xios_apply_wet_3d_r8(buf)
    real(kind=8), intent(inout) :: buf(:,:)
  end subroutine
  subroutine io_xios_apply_wet_3d_r4(buf)
    real(kind=4), intent(inout) :: buf(:,:)
  end subroutine
  subroutine io_xios_apply_wet_3d_elem_r8(buf)
    real(kind=8), intent(inout) :: buf(:,:)
  end subroutine
  subroutine io_xios_apply_wet_3d_elem_r4(buf)
    real(kind=4), intent(inout) :: buf(:,:)
  end subroutine
#endif
end module io_xios_module
