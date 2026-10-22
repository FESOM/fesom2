module cpl_fields
  !============================================================================
  ! The fields FESOM exchanges with its external components, by name.
  !
  ! define_cpl_fields registers the field set of the configured partner once,
  ! in the order the coupler announces and exchanges them. The send and
  ! receive loops dispatch on the name instead of the position, so a field
  ! keeps its meaning whatever optional fields are registered before it, and
  ! nsend/nrecv follow from what was registered.
  !
  ! Modelled on CPLNG2_ADD_FLD in OpenIFS. The names are the ones the
  ! namcouple (OASIS) or coupling.yaml (YAC) refers to.
  !============================================================================
  implicit none
  save
  private

  public :: define_cpl_fields
  public :: nsend, nrecv, cpl_send, cpl_recv
  public :: cpl_send_collection_size, cpl_recv_collection_size

  integer                                      :: nsend = 0
  integer                                      :: nrecv = 0
  character(len=32), allocatable, dimension(:) :: cpl_send
  character(len=32), allocatable, dimension(:) :: cpl_recv
  ! Number of 2D fields bundled under one name. Always 1 for OASIS; YAC
  ! exchanges some fields as collections.
  integer,           allocatable, dimension(:) :: cpl_send_collection_size
  integer,           allocatable, dimension(:) :: cpl_recv_collection_size

contains

  !____________________________________________________________________________
  !> Register the field set of the configured partner. Called once, after the
  !> namelists are read and before the coupler definition phase.
  subroutine define_cpl_fields()
    use cpl_config, only: is_coupled_to_echam, is_coupled_to_oifs, &
                          is_coupled_to_icon_a
    use g_config,   only: lwiso, use_icebergs, use_atm_ice_tskin, use_atm_mslp

    nsend = 0
    nrecv = 0
    allocate(cpl_send(0), cpl_recv(0))
    allocate(cpl_send_collection_size(0), cpl_recv_collection_size(0))

    if (is_coupled_to_oifs) then
       ! OpenIFS (AWI-CM3, AWI-ESM3) over OASIS3-MCT 5
       call add_send('sst_feom')   ! sea surface temperature [K]
       call add_send('sie_feom')   ! sea ice concentration [0-1]
       call add_send('snt_feom')   ! snow thickness [m]
       call add_send('ist_feom')   ! ice surface temperature * concentration [K]
       call add_send('sia_feom')   ! sea ice albedo [0-1]
       call add_send('u_feom')     ! eastward surface velocity [m/s]
       call add_send('v_feom')     ! northward surface velocity [m/s]
       call add_send('sit_feom')   ! effective sea ice thickness [m]
#if defined (__recom)
       call add_send('FCO2_feom')  ! CO2 flux [kgCO2 m-2 s-1]
#endif
       call add_common_atm_recv()
       call add_recv('calv_oce')   ! solid runoff (calving)
       call add_recv('u10w_oce')   ! zonal 10 m wind
       call add_recv('v10w_oce')   ! meridional 10 m wind
#if defined (__recom)
       call add_recv('XCO2_oce')   ! atmospheric CO2 mass mixing ratio
#endif
       if (use_atm_ice_tskin) call add_recv('tsk_ico')  ! ice-tile skin temperature
       if (use_atm_mslp)      call add_recv('mslp_oce') ! mean sea-level pressure

    else if (is_coupled_to_echam) then
       ! ECHAM6 (AWI-CM2, AWI-ESM2) over OASIS3-MCT 2.8
       call add_send('sst_feom')   ! sea surface temperature [degC]
       call add_send('sit_feom')   ! effective sea ice thickness [m]
       call add_send('sie_feom')   ! sea ice concentration [0-1]
       call add_send('snt_feom')   ! snow thickness [m]
       if (lwiso) then
          call add_send('o18w_oce') ! delta 18O of surface water
          call add_send('hdow_oce') ! delta D of surface water
          call add_send('o16w_oce') ! delta H216O of surface water
          call add_send('o18i_oce') ! delta 18O of sea ice
          call add_send('hdoi_oce') ! delta D of sea ice
          call add_send('o16i_oce') ! delta H216O of sea ice
       end if
       call add_common_atm_recv()
       if (lwiso) then
          call add_recv('w1_oce')  ! H218O precipitation over water
          call add_recv('w2_oce')  ! HDO precipitation over water
          call add_recv('w3_oce')  ! H216O precipitation over water
          call add_recv('i1_oce')  ! H218O snowfall over sea ice
          call add_recv('i2_oce')  ! HDO snowfall over sea ice
          call add_recv('i3_oce')  ! H216O snowfall over sea ice
       end if
       if (use_icebergs) then
          call add_recv('u10w_oce') ! zonal 10 m wind
          call add_recv('v10w_oce') ! meridional 10 m wind
       end if

    else if (is_coupled_to_icon_a) then
       ! ICON-A over YAC
       call add_send('sst_feom')                     ! sea surface temperature [K]
       call add_send('ocean_sea_ice_bundle', 3)      ! ice thickness, snow thickness, concentration
       call add_recv('taux', 2)                      ! zonal stress over water, ice
       call add_recv('tauy', 2)                      ! meridional stress over water, ice
       call add_recv('surface_fresh_water_flux', 3)  ! rain, snow, evaporation
       call add_recv('total_heat_flux', 4)           ! shortwave, longwave, sensible, latent
       call add_recv('atmosphere_sea_ice_bundle', 2) ! heat flux into the ice
       call add_recv('river_runoff')                 ! river runoff
    end if
  end subroutine define_cpl_fields

  !____________________________________________________________________________
  !> The atmosphere-to-ocean fields ECHAM6 and OpenIFS both send, in the order
  !> both namcouples use.
  subroutine add_common_atm_recv()
    call add_recv('taux_oce')   ! zonal stress over water
    call add_recv('tauy_oce')   ! meridional stress over water
    call add_recv('taux_ico')   ! zonal stress over ice
    call add_recv('tauy_ico')   ! meridional stress over ice
    call add_recv('prec_oce')   ! total precipitation
    call add_recv('snow_oce')   ! snowfall
    call add_recv('evap_oce')   ! evaporation over water
    call add_recv('subl_oce')   ! sublimation over ice
    call add_recv('heat_oce')   ! net heat flux over water
    call add_recv('heat_ico')   ! net heat flux over ice
    call add_recv('heat_swo')   ! shortwave over water
    call add_recv('hydr_oce')   ! runoff
  end subroutine add_common_atm_recv

  !____________________________________________________________________________
  subroutine add_send(name, collection_size)
    character(len=*), intent(in)           :: name
    integer,          intent(in), optional :: collection_size

    call append(cpl_send, cpl_send_collection_size, nsend, name, collection_size)
  end subroutine add_send

  !____________________________________________________________________________
  subroutine add_recv(name, collection_size)
    character(len=*), intent(in)           :: name
    integer,          intent(in), optional :: collection_size

    call append(cpl_recv, cpl_recv_collection_size, nrecv, name, collection_size)
  end subroutine add_recv

  !____________________________________________________________________________
  subroutine append(names, sizes, n, name, collection_size)
    character(len=32), allocatable, intent(inout) :: names(:)
    integer,           allocatable, intent(inout) :: sizes(:)
    integer,                        intent(inout) :: n
    character(len=*),               intent(in)    :: name
    integer,          optional,     intent(in)    :: collection_size

    integer :: size_of_this

    size_of_this = 1
    if (present(collection_size)) size_of_this = collection_size

    names = [character(len=32) :: names, name]
    sizes = [sizes, size_of_this]
    n     = n + 1
  end subroutine append

end module cpl_fields
