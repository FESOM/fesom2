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
    do i = 1, nz_cell
       z_mid(i) = 0.5_WP * (mesh%zbar(i) + mesh%zbar(i+1))
    end do
    call xios_set_axis_attr("nz",  n_glo = nz_cell, value = z_mid)
    deallocate(z_mid)
    call xios_set_axis_attr("nz1", n_glo = mesh%nl, value = mesh%zbar(1:mesh%nl))

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


  subroutine io_xios_close()
    if (.not. xios_on) return
    call xios_context_finalize()
    call xios_finalize()
    xios_on = .false.
  end subroutine

#else
  ! Stub module when built without XIOS; keeps `use io_xios_module` legal.
  implicit none
  private
  public :: io_xios_is_on
contains
  logical function io_xios_is_on() result(r)
    r = .false.
  end function
#endif
end module io_xios_module
