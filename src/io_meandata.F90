module io_MEANDATA

  use o_PARAM, only : WP
  use, intrinsic :: iso_fortran_env, only: real64, real32
  use io_data_strategy_module

  implicit none
#include "netcdf.inc"
  private
  public output
!
!--------------------------------------------------------------------------------------------
!
  integer, parameter  :: i_real8=8, i_real4=4

  integer, parameter :: root_rank = 0


  type Meandata
    private
    integer                                            :: ndim
    integer                                            :: glsize(2)
    integer                                            :: accuracy
    real(real64), allocatable, dimension(:,:) :: local_values_r8
    real(real32), allocatable, dimension(:,:) :: local_values_r4
    real(real64), allocatable :: aux_r8(:)
    integer                                            :: addcounter=0
    real(kind=WP), pointer                             :: ptr3(:,:) ! todo: use netcdf types, not WP
    character(500)                                     :: filename
    character(100)                                     :: name
    character(500)                                     :: description
    character(100)                                     :: units
    character(100)                                     :: dimname(2)
    integer                                            :: ncid
    integer                                            :: rec_count=0
    integer                                            :: recID, tID
    integer                                            :: dimID(2), varID
    integer                                            :: freq=1
    character                                          :: freq_unit='m'
    logical                                            :: is_in_use=.false.
    logical :: is_elem_based = .false.
    class(data_strategy_type), allocatable :: data_strategy
    integer :: callback_level = 0
  contains
    final destructor
  end type  
!
!--------------------------------------------------------------------------------------------
!
  type(Meandata), save, allocatable, target :: io_stream(:)
  integer, save                             :: io_NSTREAMS=0
  real(kind=WP)                             :: ctime !current time in seconds from the beginning of the year
!
!--------------------------------------------------------------------------------------------
!
  integer, save                  :: io_listsize=0
  type io_entry
        CHARACTER(len=15)        :: id        ='unknown   '
        INTEGER                  :: freq      =0
        CHARACTER                :: unit      =''
        INTEGER                  :: precision =0
  end type

  type(io_entry), save, allocatable, target   :: io_list(:)
!
!--------------------------------------------------------------------------------------------
! generic interface was required to associate variables of unknown rank with the pointers of the same rank
! this allows for automatic streaming of associated variables into the netcdf file
  INTERFACE def_stream
            MODULE PROCEDURE def_stream2D, def_stream3D
  END INTERFACE
!
!--------------------------------------------------------------------------------------------
!
  contains
!
!--------------------------------------------------------------------------------------------
!

  subroutine destructor(this)
    type(Meandata), intent(inout) :: this
    ! EO args
    call assert_nf(nf_close(this%ncid), __LINE__)
  end subroutine


subroutine ini_mean_io(mesh)
  use g_cvmix_tke
  use g_cvmix_idemix
  use g_cvmix_kpp
  use g_cvmix_tidal
  use g_PARSUP
  use diagnostics
  use i_PARAM, only: whichEVP
  implicit none
  integer                   :: i, j
  integer, save             :: nm_io_unit  = 103       ! unit to open namelist file, skip 100-102 for cray
  integer                   :: iost
  integer,dimension(12)     :: sel_forcvar=0
  character(len=10)         :: id_string

  type(t_mesh), intent(in) , target :: mesh
  namelist /nml_listsize/ io_listsize
  namelist /nml_list    / io_list

#include  "associate_mesh.h"

  ! OPEN and read namelist for I/O
  open( unit=nm_io_unit, file='namelist.io', form='formatted', access='sequential', status='old', iostat=iost )
  if (iost == 0) then
  if (mype==0) WRITE(*,*) '     file   : ', 'namelist.io',' open ok'
     else
  if (mype==0) WRITE(*,*) 'ERROR: --> bad opening file   : ', 'namelist.io',' ; iostat=',iost
     call par_ex
     stop
  endif
  READ(nm_io_unit, nml=nml_listsize, iostat=iost )
  allocate(io_list(io_listsize))
  READ(nm_io_unit, nml=nml_list,     iostat=iost )
  close(nm_io_unit )


  do i=1, io_listsize
     if (trim(io_list(i)%id)=='unknown   ') then
        if (mype==0) write(*,*) 'io_listsize will be changed from ', io_listsize, ' to ', i-1, '!'
        io_listsize=i-1
        EXIT
     end if
  end do
  
DO i=1, io_listsize
SELECT CASE (trim(io_list(i)%id))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!2D streams!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CASE ('sst       ')
    call def_stream(nod2D, myDim_nod2D, 'sst',      'sea surface temperature',        'C',      tr_arr(1,1:myDim_nod2D,1), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, mesh)
CASE ('sss       ')
    call def_stream(nod2D, myDim_nod2D, 'sss',      'sea surface salinity',           'psu',    tr_arr(1,1:myDim_nod2D,2), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, mesh)
CASE ('ssh       ')
    call def_stream(nod2D, myDim_nod2D, 'ssh',      'sea surface elevation',          'm',      eta_n,                     io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, mesh)
CASE ('vve_5     ')
    call def_stream(nod2D, myDim_nod2D, 'vve_5',    'vertical velocity at 5th level', 'm/s',    Wvel(5,:),                 io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, mesh)
    
!___________________________________________________________________________________________________________________________________
! output sea ice 
CASE ('uice      ')
    if (use_ice) then
    call def_stream(nod2D, myDim_nod2D, 'uice',     'ice velocity x',                 'm/s',    u_ice,                     io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, mesh)
    end if
CASE ('vice      ')
    if (use_ice) then
    call def_stream(nod2D, myDim_nod2D, 'vice',     'ice velocity y',                 'm/s',    v_ice,                     io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, mesh)
    end if
CASE ('a_ice     ')
    if (use_ice) then
    call def_stream(nod2D, myDim_nod2D, 'a_ice',    'ice concentration',              '%',      a_ice(1:myDim_nod2D),      io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, mesh)
    end if
CASE ('m_ice     ')
    if (use_ice) then
    call def_stream(nod2D, myDim_nod2D, 'm_ice',    'ice height',                     'm',      m_ice(1:myDim_nod2D),      io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, mesh)
    end if
CASE ('thdgr     ')
    if (use_ice) then
    call def_stream(nod2D, myDim_nod2D, 'thdgr',    'growth rate ice',                 'm/s',    thdgr(1:myDim_nod2D),      io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, mesh)
    end if
CASE ('thdgrsn   ')
    if (use_ice) then
    call def_stream(nod2D, myDim_nod2D, 'thdgrsn',  'growth rate ice',                 'm/s',    thdgrsn(1:myDim_nod2D),    io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, mesh)
    end if
CASE ('m_snow    ')
    if (use_ice) then
    call def_stream(nod2D, myDim_nod2D, 'm_snow',   'snow height',                     'm',      m_snow(1:myDim_nod2D),     io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, mesh)
    end if
    
!___________________________________________________________________________________________________________________________________
! output mixed layer depth    
CASE ('MLD1      ')
    call def_stream(nod2D, myDim_nod2D, 'MLD1',     'Mixed Layer Depth',               'm',      MLD1(1:myDim_nod2D),       io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, mesh)
CASE ('MLD2      ')
    call def_stream(nod2D, myDim_nod2D, 'MLD2',     'Mixed Layer Depth',               'm',      MLD2(1:myDim_nod2D),       io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, mesh)
    
!___________________________________________________________________________________________________________________________________
! output surface forcing
CASE ('fh        ')
    call def_stream(nod2D, myDim_nod2D, 'fh',       'heat flux',                       'W',      heat_flux(:),              io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, mesh)
CASE ('fw        ')
    call def_stream(nod2D, myDim_nod2D, 'fw',       'fresh water flux',                'm/s',    water_flux(:),             io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, mesh)
CASE ('atmice_x  ')
    call def_stream(nod2D, myDim_nod2D, 'atmice_x', 'stress atmice x',                 'N/m2',   stress_atmice_x(:),        io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, mesh)
CASE ('atmice_y  ')
    call def_stream(nod2D, myDim_nod2D, 'atmice_y', 'stress atmice y',                 'N/m2',   stress_atmice_y(:),        io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, mesh)
CASE ('atmoce_x  ')
    call def_stream(nod2D, myDim_nod2D, 'atmoce_x', 'stress atmoce x',                 'N/m2',   stress_atmoce_x(:),        io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, mesh)
CASE ('atmoce_y  ')
    call def_stream(nod2D, myDim_nod2D, 'atmoce_y', 'stress atmoce y',                 'N/m2',   stress_atmoce_y(:),        io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, mesh)
CASE ('alpha     ')
    call def_stream(nod2D, myDim_nod2D, 'alpha',    'thermal expansion',               'none',   sw_alpha(1,:),             io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, mesh)
CASE ('beta      ')
    call def_stream(nod2D, myDim_nod2D, 'beta',     'saline contraction',              'none',   sw_beta (1,:),             io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, mesh)
CASE ('runoff    ')
    sel_forcvar(10)= 1
    call def_stream(nod2D, myDim_nod2D, 'runoff',   'river runoff',                    'none',   runoff(:),                 io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, mesh)
CASE ('evap      ')
    sel_forcvar(7) = 1
    call def_stream(nod2D, myDim_nod2D, 'evap',     'evaporation',                     'm/s',    evaporation(:),            io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, mesh)
CASE ('prec      ')
    sel_forcvar(5) = 1
    call def_stream(nod2D, myDim_nod2D, 'prec',     'precicipation rain',              'm/s',    prec_rain(:),              io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, mesh)
CASE ('snow      ')
    sel_forcvar(6) = 1
    call def_stream(nod2D, myDim_nod2D, 'snow',     'precicipation snow',              'm/s',    prec_snow(:),              io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, mesh)
CASE ('tair      ')
    sel_forcvar(3) = 1
    call def_stream(nod2D, myDim_nod2D, 'tair',     'surface air temperature',         '°C',     Tair(:),                   io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, mesh)
CASE ('shum      ')
    sel_forcvar(4) = 1
    call def_stream(nod2D, myDim_nod2D, 'shum',     'specific humidity',               '',       shum(:),                   io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, mesh)
CASE ('swr ')
    sel_forcvar(8) = 1
    call def_stream(nod2D, myDim_nod2D, 'swr',      'short wave radiation',            'W/m^2',  shortwave(:),              io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, mesh)
CASE ('lwr ')
    sel_forcvar(9) = 1
    call def_stream(nod2D, myDim_nod2D, 'lwr',      'long wave radiation',             'W/m^2',  longwave(:),               io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, mesh)
CASE ('uwind ')
    sel_forcvar(1) = 1
    call def_stream(nod2D, myDim_nod2D, 'uwind',    '10m zonal surface wind velocity', 'm/s',    u_wind(:),               io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, mesh)
CASE ('vwind ')
    sel_forcvar(2) = 1
    call def_stream(nod2D, myDim_nod2D, 'vwind',    '10m merid. surface wind velocity','m/s',    v_wind(:),               io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, mesh)

    
!___________________________________________________________________________________________________________________________________
! output KPP vertical mixing schemes
CASE ('kpp_obldepth   ')
    if     (mix_scheme_nmb==1 .or. mix_scheme_nmb==17) then! fesom KPP
        call def_stream(nod2D, myDim_nod2D,    'kpp_obldepth',    'KPP ocean boundary layer depth', 'm',   hbl(:),          io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, mesh)
    elseif (mix_scheme_nmb==3 .or. mix_scheme_nmb==37) then ! cvmix KPP
        call def_stream(nod2D, myDim_nod2D,    'kpp_obldepth',    'KPP ocean boundary layer depth', 'm',   kpp_obldepth(:), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, mesh)
    end if
CASE ('kpp_sbuoyflx')
    if     (mix_scheme_nmb==1 .or. mix_scheme_nmb==17) then ! fesom KPP
        call def_stream(nod2D, myDim_nod2D,    'kpp_sbuoyflx',    'surface buoyancy flux',   'm2/s3',  Bo(:),             io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, mesh)
    elseif (mix_scheme_nmb==3 .or. mix_scheme_nmb==37) then ! cvmix KPP
        call def_stream(nod2D, myDim_nod2D,    'kpp_sbuoyflx',    'surface buoyancy flux',   'm2/s3',  kpp_sbuoyflx(:), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, mesh)
    end if
CASE ('tx_sur    ')
    sel_forcvar(11) = 1
    call def_stream(elem2D, myDim_elem2D,  'tx_sur',    'zonal wind str. to ocean',       'm/s2',   stress_surf(1, :),         io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, mesh)
CASE ('ty_sur    ')
    sel_forcvar(12) = 1
    call def_stream(elem2D, myDim_elem2D,  'ty_sur',    'meridional wind str. to ocean',  'm/s2',   stress_surf(2, :),         io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, mesh)
CASE ('curl_surf ')
    if (lcurt_stress_surf) then
    call def_stream(nod2D, myDim_nod2D,    'curl_surf', 'vorticity of the surface stress','none',   curl_stress_surf(:),       io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, mesh)
    end if
    
!___________________________________________________________________________________________________________________________________
! output Ferrari/GM parameterisation 2D  
CASE ('fer_C     ')
    if (Fer_GM) then
    call def_stream(nod2D,  myDim_nod2D,   'fer_C',     'GM,   depth independent speed',  'm/s' ,   fer_c(:),                  io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, mesh)
    end if
    
!___________________________________________________________________________________________________________________________________    
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>   3D streams   <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!___________________________________________________________________________________________________________________________________
CASE ('temp      ')
    call def_stream((/nl-1, nod2D/),  (/nl-1, myDim_nod2D/),  'temp',      'temperature', 'C',      tr_arr(:,:,1),             io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, mesh)
CASE ('salt      ')
    call def_stream((/nl-1, nod2D/),  (/nl-1, myDim_nod2D/),  'salt',      'salinity',    'psu',    tr_arr(:,:,2),             io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, mesh)
CASE ('otracers  ')
    do j=3, num_tracers
    write (id_string, "(I3.3)") tracer_id(j)
    call def_stream((/nl-1, nod2D/),  (/nl-1, myDim_nod2D/),  'tra_'//id_string, 'pasive tracer ID='//id_string, 'n/a', tr_arr(:,:,j), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, mesh)
    end do
CASE ('slope_x   ')
    call def_stream((/nl-1,  nod2D/), (/nl-1, myDim_nod2D/),  'slope_x',   'neutral slope X',    'none', slope_tapered(1,:,:), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, mesh)
CASE ('slope_y   ')
    call def_stream((/nl-1,  nod2D/), (/nl-1, myDim_nod2D/),  'slope_y',   'neutral slope Y',    'none', slope_tapered(2,:,:), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, mesh)
CASE ('slope_z   ')
    call def_stream((/nl-1,  nod2D/), (/nl-1, myDim_nod2D/),  'slope_z',   'neutral slope Z',    'none', slope_tapered(3,:,:), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, mesh)
CASE ('N2        ')
    call def_stream((/nl,    nod2D/), (/nl,   myDim_nod2D/),  'N2',        'brunt väisälä',      '1/s2', bvfreq(:,:),          io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, mesh)
CASE ('Kv        ')
    call def_stream((/nl,    nod2D/), (/nl,   myDim_nod2D/),  'Kv',        'vertical diffusivity Kv',  'm2/s', Kv(:,:),              io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, mesh)
CASE ('u         ')
    call def_stream((/nl-1, elem2D/), (/nl-1, myDim_elem2D/), 'u',         'horizontal velocity','m/s',  uv(1,:,:),            io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, mesh)
CASE ('v         ')
    call def_stream((/nl-1, elem2D/), (/nl-1, myDim_elem2D/), 'v',         'meridional velocity','m/s',  uv(2,:,:),            io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, mesh)
CASE ('w         ')
    call def_stream((/nl,    nod2D/), (/nl,   myDim_nod2D/),  'w',         'vertical velocity',  'm/s',  Wvel(:,:),            io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, mesh)
CASE ('Av        ')
    call def_stream((/nl,   elem2D/), (/nl,   myDim_elem2D/), 'Av',        'vertical viscosity Av',  'm2/s', Av(:,:),              io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, mesh)
    
!___________________________________________________________________________________________________________________________________
! output Ferrari/GM parameterisation
CASE ('bolus_u   ')
    if (Fer_GM) then
    call def_stream((/nl-1, elem2D/), (/nl-1, myDim_elem2D/), 'bolus_u',   'GM bolus velocity U','m/s',  fer_uv(1,:,:),        io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, mesh)
    end if
CASE ('bolus_v   ')
    if (Fer_GM) then  
    call def_stream((/nl-1, elem2D/), (/nl-1, myDim_elem2D/), 'bolus_v',   'GM bolus velocity V','m/s',  fer_uv(2,:,:),        io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, mesh)
    end if
CASE ('bolus_w   ')
    if (Fer_GM) then
    call def_stream((/nl  , nod2D /), (/nl,   myDim_nod2D /), 'bolus_w',   'GM bolus velocity W','m/s',  fer_Wvel(:,:),        io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, mesh)
    end if
CASE ('fer_K     ')
    if (Fer_GM) then
    call def_stream((/nl  , nod2D /), (/nl,   myDim_nod2D /), 'fer_K',     'GM, stirring diff.','m2/s',  fer_k(:,:),           io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, mesh)
    end if
CASE ('fer_scal  ')
    if (Fer_GM) then
    call def_stream(        nod2D   ,         myDim_nod2D   , 'fer_scal',  'GM surface scaling','',  fer_scal(:),           io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, mesh)
    end if
CASE ('dMOC      ')
    if (ldiag_dMOC) then
       call def_stream((/std_dens_N, elem2D/),  (/std_dens_N, myDim_elem2D/), 'U_rho_x_DZ',     'fluxes for density MOC', 'fluxes', std_dens_UVDZ(1,:,:),   1, 'm', i_real4, mesh)
       call def_stream((/std_dens_N, elem2D/),  (/std_dens_N, myDim_elem2D/), 'V_rho_x_DZ',     'fluxes for density MOC', 'fluxes', std_dens_UVDZ(2,:,:),   1, 'm', i_real4, mesh)
       call def_stream((/std_dens_N, elem2D/),  (/std_dens_N, myDim_elem2D/), 'RHO_Z',          'drho/dz',                'kg/m4' , std_dens_RHOZ(:,:),     1, 'm', i_real4, mesh)
       call def_stream((/nl-1,       nod2D /),  (/nl-1,       myDim_nod2D /), 'density_dMOC',   'density'               , 'm',      density_dmoc(:,:),      1, 'm', i_real4, mesh)
    end if
!___________________________________________________________________________________________________________________________________
CASE ('pgf_x     ')    
    call def_stream((/nl-1, elem2D/), (/nl-1, myDim_elem2D/), 'pgf_x', 'zonal pressure gradient force'     , 'm/s^2', pgf_x(:,:), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, mesh)
CASE ('pgf_y     ')    
    call def_stream((/nl-1, elem2D/), (/nl-1, myDim_elem2D/), 'pgf_y', 'meridional pressure gradient force', 'm/s^2', pgf_y(:,:), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, mesh)
!___________________________________________________________________________________________________________________________________    
CASE DEFAULT
    if (mype==0) write(*,*) 'stream ', io_list(i)%id, ' is not defined !'
END SELECT
END DO


!3D
  if (ldiag_energy) then
     call def_stream((/nl,   nod2D/),  (/nl,     myDim_nod2D/), 'rhof',     'in-situ density at faces',    'kg/m3',     rhof(:,:),  1, 'm', i_real8, mesh)
     call def_stream((/nl,   nod2D/),  (/nl,     myDim_nod2D/), 'wrhof',    'vertical velocity x density', 'kg/(s*m2)', wrhof(:,:), 1, 'm', i_real8, mesh)
     call def_stream((/nl-1, nod2D/),  (/nl-1,   myDim_nod2D/), 'uu',   'u times u', 'm2/s2', u_x_u(:,:), 1, 'm', i_real8, mesh)
     call def_stream((/nl-1, nod2D/),  (/nl-1,   myDim_nod2D/), 'uv',   'u times v', 'm2/s2', u_x_v(:,:), 1, 'm', i_real8, mesh)
     call def_stream((/nl-1, nod2D/),  (/nl-1,   myDim_nod2D/), 'vv',   'v times v', 'm2/s2', v_x_v(:,:), 1, 'm', i_real8, mesh)
     call def_stream((/nl,  elem2D/),  (/nl-1,   myDim_elem2D/),'uw',   'u times w', 'm2/s2', u_x_w(:,:), 1, 'm', i_real8, mesh)
     call def_stream((/nl,  elem2D/),  (/nl-1,   myDim_elem2D/),'vw',   'v times w', 'm2/s2', v_x_w(:,:), 1, 'm', i_real8, mesh)
     call def_stream((/nl-1, nod2D/),  (/nl-1,   myDim_nod2D/), 'dudx', 'du/dx',     '1/s',   dudx(:,:),  1, 'm', i_real8, mesh)
     call def_stream((/nl-1, nod2D/),  (/nl-1,   myDim_nod2D/), 'dudy', 'du/dy',     '1/s',   dudy(:,:),  1, 'm', i_real8, mesh)
     call def_stream((/nl-1, nod2D/),  (/nl-1,   myDim_nod2D/), 'dvdx', 'dv/dx',     '1/s',   dvdx(:,:),  1, 'm', i_real8, mesh)
     call def_stream((/nl-1, nod2D/),  (/nl-1,   myDim_nod2D/), 'dvdy', 'dv/dy',     '1/s',   dvdy(:,:),  1, 'm', i_real8, mesh)
     call def_stream((/nl,  elem2D/),  (/nl,     myDim_elem2D/), 'dudz', 'du/dz',     '1/s',   dudz(:,:),  1, 'm', i_real8, mesh)
     call def_stream((/nl,  elem2D/),  (/nl,     myDim_elem2D/), 'dvdz', 'dv/dz',     '1/s',   dvdz(:,:),  1, 'm', i_real8, mesh)
     call def_stream((/nl,  elem2D/),  (/nl,     myDim_elem2D/), 'av_dudz', 'int(Av * du/dz)',        'm3/s2',   av_dudz(:,:),    1, 'm', i_real4, mesh)
     call def_stream((/nl,  elem2D/),  (/nl,     myDim_elem2D/), 'av_dvdz', 'int(Av * dv/dz)',        'm3/s2',   av_dvdz(:,:),    1, 'm', i_real4, mesh)
     call def_stream((/nl,  elem2D/),  (/nl,     myDim_elem2D/), 'av_dudz_sq',  'Av * (du/dz)^2',     'm^2/s^3', av_dudz_sq(:,:), 1, 'm', i_real4, mesh)
     call def_stream((/nl,   elem2D/), (/nl,     myDim_elem2D/), 'Av',    'Vertical mixing A',         'm2/s',            Av(:,:), 1, 'm', i_real4, mesh)
     call def_stream((/nl-1, nod2D/),  (/nl-1,   myDim_nod2D/),  'unod',  'horizontal velocity at nodes', 'm/s', Unode(1,:,:), 1, 'm', i_real8, mesh)
     call def_stream((/nl-1, nod2D/),  (/nl-1,   myDim_nod2D/),  'vnod',  'meridional velocity at nodes', 'm/s', Unode(2,:,:), 1, 'm', i_real8, mesh)
    
     call def_stream((/nl-1, elem2D/), (/nl-1,   myDim_elem2D/), 'um',  'horizontal velocity', 'm/s', uv(1,:,:),     1, 'm', i_real4, mesh)
     call def_stream((/nl-1, elem2D/), (/nl-1,   myDim_elem2D/), 'vm',  'meridional velocity', 'm/s', uv(2,:,:),     1, 'm', i_real4, mesh)
     call def_stream((/nl, nod2D/),    (/nl,     myDim_nod2D/),  'wm',  'vertical velocity',   'm/s', Wvel(:,:),     1, 'm', i_real8, mesh)

     call def_stream(elem2D, myDim_elem2D,   'utau_surf',  '(u, tau) at the surface',      'N/(m s)', utau_surf(1:myDim_elem2D),     1, 'm', i_real4, mesh)
     call def_stream(elem2D, myDim_elem2D,   'utau_bott',  '(u, tau) at the bottom',       'N/(m s)', utau_bott(1:myDim_elem2D),     1, 'm', i_real4, mesh)
     call def_stream(elem2D, myDim_elem2D,   'u_bott',     'bottom velocity',                  'm/s', u_bott(1:myDim_elem2D),        1, 'm', i_real4, mesh)
     call def_stream(elem2D, myDim_elem2D,   'v_bott',     'bottom velocity',                  'm/s', v_bott(1:myDim_elem2D),        1, 'm', i_real4, mesh)
     call def_stream(elem2D, myDim_elem2D,   'u_surf',     'surface velocity',                 'm/s', u_surf(1:myDim_elem2D),        1, 'm', i_real4, mesh)
     call def_stream(elem2D, myDim_elem2D,   'v_surf',     'surface velocity',                 'm/s', u_surf(1:myDim_elem2D),        1, 'm', i_real4, mesh)
     call def_stream(elem2D, myDim_elem2D,   'tx_bot',     'bottom stress x',                 'N/m2', stress_bott(1, 1:myDim_elem2D),1, 'm', i_real4, mesh)
     call def_stream(elem2D, myDim_elem2D,   'ty_bot',     'bottom stress y',                 'N/m2', stress_bott(2, 1:myDim_elem2D),1, 'm', i_real4, mesh)
     if (sel_forcvar(11)==0) call def_stream(elem2D, myDim_elem2D,   'tx_sur',     'zonal wind stress to ocean',     'm/s2', stress_surf(1, 1:myDim_elem2D),1, 'm', i_real4, mesh) ; sel_forcvar(11)=1
     if (sel_forcvar(12)==0) call def_stream(elem2D, myDim_elem2D,   'ty_sur',     'meridional wind stress to ocean','m/s2', stress_surf(2, 1:myDim_elem2D),1, 'm', i_real4, mesh) ; sel_forcvar(12)=1
  end if

#if defined (__oifs)
  call def_stream(nod2D, myDim_nod2D, 'alb',    'ice albedo',             'none',   ice_alb(:),                    1, 'm', i_real4, mesh)
#endif
    
    if (mix_scheme_nmb==5 .or. mix_scheme_nmb==56) then
        ! TKE diagnostic 
        call def_stream((/nl,nod2D/), (/nl,myDim_nod2D/), 'tke'     , 'turbulent kinetic energy'                    , 'm^2/s^2', tke(:,:)     , 1, 'y', i_real4, mesh)
        call def_stream((/nl,nod2D/), (/nl,myDim_nod2D/), 'tke_Ttot', 'total production of turbulent kinetic energy', 'm^2/s^3', tke_Ttot(:,:), 1, 'y', i_real4, mesh)
        call def_stream((/nl,nod2D/), (/nl,myDim_nod2D/), 'tke_Tbpr', 'TKE production by buoyancy'                  , 'm^2/s^3', tke_Tbpr(:,:), 1, 'y', i_real4, mesh)
        call def_stream((/nl,nod2D/), (/nl,myDim_nod2D/), 'tke_Tspr', 'TKE production by shear'                     , 'm^2/s^3', tke_Tspr(:,:), 1, 'y', i_real4, mesh)
        call def_stream((/nl,nod2D/), (/nl,myDim_nod2D/), 'tke_Tdif', 'TKE production by vertical diffusion'        , 'm^2/s^3', tke_Tdif(:,:), 1, 'y', i_real4, mesh)
        call def_stream((/nl,nod2D/), (/nl,myDim_nod2D/), 'tke_Tdis', 'TKE production by dissipation'               , 'm^2/s^3', tke_Tdis(:,:), 1, 'y', i_real4, mesh)
        call def_stream((/nl,nod2D/), (/nl,myDim_nod2D/), 'tke_Twin', 'TKE production by wind'                      , 'm^2/s^3', tke_Twin(:,:), 1, 'y', i_real4, mesh)
        call def_stream((/nl,nod2D/), (/nl,myDim_nod2D/), 'tke_Tbck', 'background forcing for TKE'                  , 'm^2/s^3', tke_Tbck(:,:), 1, 'y', i_real4, mesh)
        call def_stream((/nl,nod2D/), (/nl,myDim_nod2D/), 'tke_Lmix', 'mixing length scale of TKE'                  , 'm'      , tke_Lmix(:,:), 1, 'y', i_real4, mesh)
        call def_stream((/nl,nod2D/), (/nl,myDim_nod2D/), 'tke_Pr'  , 'Prantl number'                               , ''       , tke_Pr(:,:)  , 1, 'y', i_real4, mesh)
        if (mix_scheme_nmb==56) then
            ! TKE-IDEMIX diagnostic 
            call def_stream((/nl,nod2D/), (/nl,myDim_nod2D/), 'tke_Tiwf', 'TKE production by internal waves (IDEMIX)', 'm^2/s^3', tke_Tiwf(:,:), 1, 'y', i_real4, mesh)
        end if 
    end if 
    
    if (mod(mix_scheme_nmb,10)==6) then
        ! IDEMIX Internal-Wave-Energy diagnostics
        call def_stream((/nl,nod2D/), (/nl,myDim_nod2D/), 'iwe'     , 'internal wave energy'                    , 'm^2/s^2', iwe(:,:)     , 1, 'y', i_real4, mesh)
        call def_stream((/nl,nod2D/), (/nl,myDim_nod2D/), 'iwe_Ttot', 'total production of internal wave energy', 'm^2/s^2', iwe_Ttot(:,:), 1, 'y', i_real4, mesh)
        call def_stream((/nl,nod2D/), (/nl,myDim_nod2D/), 'iwe_Tdif', 'IWE production by vertical diffusion'    , 'm^2/s^3', iwe_Tdif(:,:), 1, 'y', i_real4, mesh)
        call def_stream((/nl,nod2D/), (/nl,myDim_nod2D/), 'iwe_Tdis', 'IWE production by dissipation'           , 'm^2/s^3', iwe_Tdis(:,:), 1, 'y', i_real4, mesh)
        call def_stream((/nl,nod2D/), (/nl,myDim_nod2D/), 'iwe_Tsur', 'IWE production from surface forcing'     , 'm^2/s^2', iwe_Tsur(:,:), 1, 'y', i_real4, mesh)
        call def_stream((/nl,nod2D/), (/nl,myDim_nod2D/), 'iwe_Tbot', 'IWE production from bottom forcing'      , 'm^2/s^2', iwe_Tbot(:,:), 1, 'y', i_real4, mesh)
        call def_stream((/nl,nod2D/), (/nl,myDim_nod2D/), 'iwe_c0'  , 'IWE vertical group velocity'             , 'm/s'    , iwe_c0(:,:)  , 1, 'y', i_real4, mesh)
        call def_stream((/nl,nod2D/), (/nl,myDim_nod2D/), 'iwe_v0'  , 'IWE horizontal group velocity'           , 'm/s'    , iwe_c0(:,:)  , 1, 'y', i_real4, mesh)
    end if
    
    if (mod(mix_scheme_nmb,10)==7) then
        ! cvmix_TIDAL diagnostics
        call def_stream((/nl,nod2D/), (/nl,myDim_nod2D/), 'tidal_Kv'  , 'tidal diffusivity'           , 'm^2/s'    , tidal_Kv(:,:)  , 1, 'y', i_real4, mesh)
        call def_stream((/nl,nod2D/), (/nl,myDim_nod2D/), 'tidal_Av'  , 'tidal viscosity'             , 'm^2/s'    , tidal_Av(:,:)  , 1, 'y', i_real4, mesh)
        call def_stream(     nod2D  ,      myDim_nod2D  , 'tidal_forcbot', 'near tidal bottom forcing', 'W/m^2'    , tidal_forc_bottom_2D  , 100, 'y', i_real4, mesh)
    end if
    
  !___________________________________________________________________________________________________________________________________
  ! output Redi parameterisation
  if (Redi) then
     call def_stream((/nl-1  , nod2D /), (/nl-1,   myDim_nod2D /), 'Redi_K',   'Redi diffusion coefficient', 'm2/s', Ki(:,:),    1, 'y', i_real4, mesh)
  end if

  !___________________________________________________________________________________________________________________________________
  ! output Monin-Obukov (TB04) mixing length
  if (use_momix) then
     call def_stream(nod2D, myDim_nod2D, 'momix_length',   'Monin-Obukov mixing length', 'm', mixlength(:),    1, 'm', i_real4, mesh)
  end if
  
    !___________________________________________________________________________________________________________________________________
    if (ldiag_curl_vel3) then
        call def_stream((/nl-1, nod2D/),  (/nl-1, myDim_nod2D/),  'curl_u',     'relative vorticity',          '1/s',   vorticity,                   1, 'm', i_real4, mesh)
    end if

    !___________________________________________________________________________________________________________________________________
    if (whichEVP==1) then
    end if
    
    if (whichEVP==2) then
        call def_stream(elem2D, myDim_elem2D, 'alpha_EVP', 'alpha in EVP', 'n/a', alpha_evp_array,  1, 'd', i_real4, mesh)
        call def_stream(nod2D,  myDim_nod2D,  'beta_EVP',  'beta in EVP',  'n/a', beta_evp_array,   1, 'd', i_real4, mesh)
    end if
  
    !___________________________________________________________________________
    if (ldiag_dvd) then
        call def_stream((/nl-1, nod2D/), (/nl-1, myDim_nod2D/), 'dvd_temp_h', 'horiz. dvd of temperature', '°C/s' , tr_dvd_horiz(:,:,1), 1, 'm', i_real4, mesh)
        call def_stream((/nl-1, nod2D/), (/nl-1, myDim_nod2D/), 'dvd_temp_v', 'vert. dvd of temperature' , '°C/s' , tr_dvd_vert(:,:,1) , 1, 'm', i_real4, mesh)
        call def_stream((/nl-1, nod2D/), (/nl-1, myDim_nod2D/), 'dvd_salt_h', 'horiz. dvd of salinity'   , 'psu/s', tr_dvd_horiz(:,:,2), 1, 'm', i_real4, mesh)
        call def_stream((/nl-1, nod2D/), (/nl-1, myDim_nod2D/), 'dvd_salt_v', 'vert. dvd of salinity'    , 'psu/s', tr_dvd_vert(:,:,2) , 1, 'm', i_real4, mesh)
    end if 
    
    !___________________________________________________________________________
    if (ldiag_forc) then
        if (sel_forcvar( 1)==0) call def_stream(nod2D , myDim_nod2D , 'uwind' , '10m zonal surface wind velocity', 'm/s'  , u_wind(:)        , 1, 'm', i_real4, mesh)
        if (sel_forcvar( 2)==0) call def_stream(nod2D , myDim_nod2D , 'vwind' , '10m merid surface wind velocity', 'm/s'  , v_wind(:)        , 1, 'm', i_real4, mesh)
        if (sel_forcvar( 3)==0) call def_stream(nod2D , myDim_nod2D , 'tair'  , 'surface air temperature'        , '°C'   , Tair(:)          , 1, 'm', i_real4, mesh)
        if (sel_forcvar( 4)==0) call def_stream(nod2D , myDim_nod2D , 'shum'  , 'specific humidity'              , ''     , shum(:)          , 1, 'm', i_real4, mesh)
        if (sel_forcvar( 5)==0) call def_stream(nod2D , myDim_nod2D , 'prec'  , 'precicipation rain'             , 'm/s'  , prec_rain(:)     , 1, 'm', i_real4, mesh)
        if (sel_forcvar( 6)==0) call def_stream(nod2D , myDim_nod2D , 'snow'  , 'precicipation snow'             , 'm/s'  , prec_snow(:)     , 1, 'm', i_real4, mesh)
        if (sel_forcvar( 7)==0) call def_stream(nod2D , myDim_nod2D , 'evap'  , 'evaporation'                    , 'm/s'  , evaporation(:)   , 1, 'm', i_real4, mesh)
        if (sel_forcvar( 8)==0) call def_stream(nod2D , myDim_nod2D , 'swr'   , 'short wave radiation'           , 'W/m^2', shortwave(:)     , 1, 'm', i_real4, mesh)
        if (sel_forcvar( 9)==0) call def_stream(nod2D , myDim_nod2D , 'lwr'   , 'long wave radiation'            , 'W/m^2', longwave(:)      , 1, 'm', i_real4, mesh)
        if (sel_forcvar(10)==0) call def_stream(nod2D , myDim_nod2D , 'runoff', 'river runoff'                   , 'none' , runoff(:)        , 1, 'm', i_real4, mesh)
        if (sel_forcvar(11)==0) call def_stream(elem2D, myDim_elem2D, 'tx_sur', 'zonal wind str. to ocean'       , 'm/s^2', stress_surf(1, :), 1, 'm', i_real4, mesh)
        if (sel_forcvar(12)==0) call def_stream(elem2D, myDim_elem2D, 'ty_sur', 'meridional wind str. to ocean'  , 'm/s^2', stress_surf(2, :), 1, 'm', i_real4, mesh)
        call def_stream(nod2D , myDim_nod2D , 'cd','wind drag coef. '             , '', cd_atm_oce_arr(:), 1, 'm', i_real4, mesh)
        call def_stream(nod2D , myDim_nod2D , 'ch','transfer coeff. sensible heat', '', ch_atm_oce_arr(:), 1, 'm', i_real4, mesh)
        call def_stream(nod2D , myDim_nod2D , 'ce','transfer coeff. evaporation ' , '', ce_atm_oce_arr(:), 1, 'm', i_real4, mesh)
    end if
    
    
end subroutine
!
!--------------------------------------------------------------------------------------------
!
function mesh_dimname_from_dimsize(size, mesh) result(name)
  use mod_mesh
  use g_PARSUP
  use diagnostics
  implicit none
  integer       :: size
  type(t_mesh) mesh
  character(50) :: name

  if (size==mesh%nod2D) then
    name='nod2'
  elseif (size==mesh%elem2D) then
    name='elem'
  elseif (size==mesh%nl) then
    name='nz'
  elseif (size==mesh%nl-1) then
    name='nz1'
  elseif (size==std_dens_N) then
    name='ndens'
  else
    name='unknown'
    if (mype==0) write(*,*) 'WARNING: unknown dimension in mean I/O with size of ', size
  end if
end function
!
!--------------------------------------------------------------------------------------------
!
subroutine create_new_file(entry)
  use g_clock
  use g_PARSUP
  implicit none
  character(2000)               :: att_text

  type(Meandata), intent(inout) :: entry
  ! Serial output implemented so far
  if (mype/=root_rank) return
  ! create an ocean output file
  write(*,*) 'initializing I/O file for ', trim(entry%name)

  call assert_nf( nf_create(entry%filename, IOR(NF_NOCLOBBER,IOR(NF_NETCDF4,NF_CLASSIC_MODEL)), entry%ncid), __LINE__)

!___Create mesh related dimensions__________________________________________
  if (entry%ndim==1) then
     call assert_nf( nf_def_dim(entry%ncid, entry%dimname(1), entry%glsize(2), entry%dimID(1)), __LINE__)
  else if (entry%ndim==2) then
     call assert_nf( nf_def_dim(entry%ncid, entry%dimname(1), entry%glsize(1), entry%dimID(1)), __LINE__)
     call assert_nf( nf_def_dim(entry%ncid, entry%dimname(2), entry%glsize(2), entry%dimID(2)), __LINE__)
  end if
!___Create time related dimensions__________________________________________
  call assert_nf( nf_def_dim(entry%ncid, 'time', NF_UNLIMITED, entry%recID), __LINE__)
!___Define the time and iteration variables_________________________________
  call assert_nf( nf_def_var(entry%ncid, 'time', NF_DOUBLE, 1, entry%recID, entry%tID), __LINE__)

  att_text='time'
  call assert_nf( nf_put_att_text(entry%ncid, entry%tID, 'long_name', len_trim(att_text), trim(att_text)), __LINE__)
  write(att_text, '(a14,I4.4,a1,I2.2,a1,I2.2,a6)') 'seconds since ', yearold, '-', 1, '-', 1, ' 0:0:0'
  call assert_nf( nf_put_att_text(entry%ncid, entry%tID, 'units', len_trim(att_text), trim(att_text)), __LINE__)

  call assert_nf( nf_def_var(entry%ncid, trim(entry%name), entry%data_strategy%netcdf_type(), entry%ndim+1, &
                                    (/entry%dimid(1:entry%ndim), entry%recID/), entry%varID), __LINE__)
  call assert_nf( nf_put_att_text(entry%ncid, entry%varID, 'description', len_trim(entry%description), entry%description), __LINE__)
  call assert_nf( nf_put_att_text(entry%ncid, entry%varID, 'units',       len_trim(entry%units),       entry%units), __LINE__)
  call assert_nf( nf_close(entry%ncid), __LINE__)
end subroutine
!
!--------------------------------------------------------------------------------------------
!
subroutine assoc_ids(entry)
  use g_PARSUP
  implicit none

  type(Meandata), intent(inout) :: entry
  integer                       :: j

  write(*,*) 'associating mean I/O file ', trim(entry%filename)

  do j=1, entry%ndim
     call assert_nf( nf_inq_dimid(entry%ncid, entry%dimname(j), entry%dimID(j)), __LINE__)
  end do
!___Associate time related dimensions_______________________________________
  call assert_nf( nf_inq_dimid (entry%ncid, 'time', entry%recID), __LINE__)
  call assert_nf( nf_inq_dimlen(entry%ncid, entry%recID, entry%rec_count), __LINE__)
!___Associate the time and iteration variables______________________________
  call assert_nf( nf_inq_varid(entry%ncid, 'time', entry%tID), __LINE__)
!___Associate physical variables____________________________________________
  call assert_nf( nf_inq_varid(entry%ncid, entry%name, entry%varID), __LINE__)
end subroutine
!
!--------------------------------------------------------------------------------------------
!
subroutine write_mean(entry, entry_index)
  use mod_mesh
  use g_PARSUP
  use io_gather_module
  implicit none
  type(Meandata), intent(inout) :: entry
  real(real32),   allocatable   :: aux_r4(:)
  integer                       :: size1, size2
  integer                       :: lev
  integer, intent(in) :: entry_index


  ! Serial output implemented so far
  if (mype==root_rank) then
     write(*,*) 'writing mean record for ', trim(entry%name), '; rec. count = ', entry%rec_count
     call assert_nf( nf_put_vara_double(entry%ncid, entry%Tid, entry%rec_count, 1, ctime, 1), __LINE__)
  end if
! !_______writing 2D and 3D fields________________________________________________
  size1=entry%glsize(1)
  size2=entry%glsize(2)
!___________writing 8 byte real_________________________________________ 
  if (entry%accuracy == i_real8) then
     if (mype==root_rank) allocate(entry%aux_r8(size2))
     do lev=1, size1
       if(.not. entry%is_elem_based) then
         call gather_nod2D (entry%local_values_r8(lev,1:size(entry%local_values_r8,dim=2)), entry%aux_r8, root_rank)
       else
         call gather_elem2D(entry%local_values_r8(lev,1:size(entry%local_values_r8,dim=2)), entry%aux_r8, root_rank)
       end if
        if (mype==root_rank) then
          entry%callback_level = lev
          call write_netcdf_callback(entry_index)
        end if
     end do
     if (mype==root_rank) deallocate(entry%aux_r8)
!___________writing real 4 byte real _________________________________________ 
  elseif (entry%accuracy == i_real4) then
     if (mype==root_rank) allocate(aux_r4(size2))
     do lev=1, size1
       if(.not. entry%is_elem_based) then
         call gather_real4_nod2D(entry%local_values_r4(lev,1:size(entry%local_values_r4,dim=2)), aux_r4, root_rank)
       else
         call gather_real4_elem2D(entry%local_values_r4(lev,1:size(entry%local_values_r4,dim=2)), aux_r4, root_rank)
       end if
        if (mype==root_rank) then
           if (entry%ndim==1) then
             call assert_nf( nf_put_vara_real(entry%ncid, entry%varID, (/1, entry%rec_count/), (/size2, 1/), aux_r4, 1), __LINE__)
           elseif (entry%ndim==2) then
             call assert_nf( nf_put_vara_real(entry%ncid, entry%varID, (/lev, 1, entry%rec_count/), (/1, size2, 1/), aux_r4, 1), __LINE__)
           end if
        end if
     end do
     if (mype==root_rank) deallocate(aux_r4)
  endif

end subroutine


subroutine write_netcdf_callback(entry_index)
  integer, intent(in) :: entry_index
  ! EO args
  type(Meandata), pointer :: entry
  integer size2
  integer lev
  
  entry=>io_stream(entry_index)
  size2 = entry%glsize(2)
  lev = entry%callback_level

  if (entry%ndim==1) then
    call assert_nf( nf_put_vara_double(entry%ncid, entry%varID, (/1, entry%rec_count/), (/size2, 1/), entry%aux_r8, 1), __LINE__)
  elseif (entry%ndim==2) then
    call assert_nf( nf_put_vara_double(entry%ncid, entry%varID, (/lev, 1, entry%rec_count/), (/1, size2, 1/), entry%aux_r8, 1), __LINE__)
  end if
end subroutine
!
!--------------------------------------------------------------------------------------------
!
subroutine update_means
  use g_PARSUP
  implicit none
  type(Meandata), pointer :: entry
  integer                 :: n

  do n=1, io_NSTREAMS
     entry=>io_stream(n)
!_____________ compute in 8 byte accuracy _________________________
    if (entry%accuracy == i_real8) then
      entry%local_values_r8 = entry%local_values_r8 + entry%ptr3(1:size(entry%local_values_r8,dim=1),1:size(entry%local_values_r8,dim=2))

!_____________ compute in 4 byte accuracy _________________________
    elseif (entry%accuracy == i_real4) then
      entry%local_values_r4 = entry%local_values_r4 + real(entry%ptr3(1:size(entry%local_values_r4,dim=1),1:size(entry%local_values_r4,dim=2)),real32)
     endif

     entry%addcounter=entry%addcounter+1
  end do
end subroutine
!
!--------------------------------------------------------------------------------------------
!
subroutine output(istep, mesh)
  use g_clock
  use mod_mesh
  use g_PARSUP
  implicit none

  integer       :: istep
  logical, save :: lfirst=.true.
  integer       :: n, k
  logical       :: do_output
  type(Meandata), pointer :: entry
  type(t_mesh), intent(in) , target :: mesh
  character(:), allocatable :: filepath
  real(real64)                  :: rtime !timestamp of the record

  ctime=timeold+(dayold-1.)*86400
  if (lfirst) call ini_mean_io(mesh)

  call update_means

  do n=1, io_NSTREAMS
     entry=>io_stream(n)
     !check whether restart will be written
     do_output=.false.

     if (entry%freq_unit.eq.'y') then
        call annual_event(do_output)

     else if (entry%freq_unit == 'm') then 
        call monthly_event(do_output) 

     else if (entry%freq_unit == 'd') then
        call daily_event(do_output, entry%freq)

     else if (entry%freq_unit == 'h') then
        call hourly_event(do_output, entry%freq)

     else if (entry%freq_unit == 's') then
        call step_event(do_output, istep, entry%freq)

     else
        write(*,*) 'You did not specify a supported outputflag.'
        write(*,*) 'The program will stop to give you opportunity to do it.'
        call par_ex(1)
        stop
     endif

     if (do_output) then
        filepath = trim(ResultPath)//trim(entry%name)//'.'//trim(runid)//'.'//cyearnew//'.nc'
        if(mype == root_rank) then
          if(filepath /= trim(entry%filename)) then
            if("" /= trim(entry%filename)) call assert_nf(nf_close(entry%ncid), __LINE__)   
            entry%filename = filepath
            call create_new_file(entry)
            call assert_nf( nf_open(entry%filename, nf_write, entry%ncid), __LINE__)
            call assoc_ids(entry)
          end if

          !___if the time rtime at the rec_count is larger than ctime we look for the closest record with the timestamp less than ctime
          do k=entry%rec_count, 1, -1
             call assert_nf( nf_get_vara_double(entry%ncid, entry%tID, k, 1, rtime, 1), __LINE__)
             if (ctime > rtime) then
                entry%rec_count=k+1
                exit ! a proper rec_count detected, exit the loop
             end if
             if (k==1) then
                write(*,*) 'I/O '//trim(entry%name)//' WARNING: the existing output file will be overwritten'//'; ', entry%rec_count, ' records in the file;'
                entry%rec_count=1
                exit ! no appropriate rec_count detected
             end if
          end do
          entry%rec_count=max(entry%rec_count, 1)
          write(*,*) trim(entry%name)//': current mean I/O counter = ', entry%rec_count
        end if

        call do_output_callback(n)
     endif
  end do
  lfirst=.false.
end subroutine


subroutine do_output_callback(entry_index)
  integer, intent(in) :: entry_index
  ! EO args
  type(Meandata), pointer :: entry

  entry=>io_stream(entry_index)

  if (entry%accuracy == i_real8) then
     entry%local_values_r8 = entry%local_values_r8 /real(entry%addcounter,real64)  ! compute_means
     call write_mean(entry, entry_index)
     entry%local_values_r8 = 0._real64 ! clean_meanarrays

  elseif (entry%accuracy == i_real4) then
     entry%local_values_r4 = entry%local_values_r4 /real(entry%addcounter,real32) ! compute_means
     call write_mean(entry, entry_index)
     entry%local_values_r4 = 0._real32 ! clean_meanarrays

  endif  ! accuracy

  entry%addcounter   = 0  ! clean_meanarrays
end subroutine

!
!--------------------------------------------------------------------------------------------
!
subroutine def_stream3D(glsize, lcsize, name, description, units, data, freq, freq_unit, accuracy, mesh)
  use mod_mesh
  use g_PARSUP
  implicit none
  integer,               intent(in)    :: glsize(2), lcsize(2)
  character(len=*),      intent(in)    :: name, description, units
  real(kind=WP), target, intent(inout) :: data(:,:)
  integer,               intent(in)    :: freq
  character,             intent(in)    :: freq_unit
  integer,               intent(in)    :: accuracy
  type(Meandata),        allocatable   :: tmparr(:)
  type(Meandata),        pointer       :: entry
  type(t_mesh), intent(in), target     :: mesh
  integer i
  
  do i = 1, rank(data)
    if ((ubound(data, dim = i)<=0)) then
      if (mype==0) then
        write(*,*) 'WARNING: adding I/O stream for ', trim(name), ' failed (contains 0 dimension)'
        write(*,*) 'upper bound is: ', ubound(data, dim = i)
      end if
      return
    end if    
  end do

  if (mype==0) then
     write(*,*) 'adding I/O stream for ', trim(name)
  end if

   ! add this instance to io_stream array
  if ( .not. allocated(io_stream)) then
    allocate(io_stream(1))
  else
    allocate(tmparr(size(io_stream)+1))
    tmparr(1:size(io_stream)) = io_stream
    deallocate(io_stream)
    call move_alloc(tmparr, io_stream)
  end if
  entry=>io_stream(size(io_stream))
  entry%ptr3 => data                      !2D! entry%ptr3(1:1,1:size(data)) => data
  entry%accuracy = accuracy

  if (accuracy == i_real8) then
    allocate(data_strategy_nf_double_type :: entry%data_strategy)
     allocate(entry%local_values_r8(lcsize(1), lcsize(2)))          !2D! allocate(entry%local_values_r8(1, lcsize))
     entry%local_values_r8 = 0._real64
  elseif (accuracy == i_real4) then
    allocate(data_strategy_nf_float_type :: entry%data_strategy)
     allocate(entry%local_values_r4(lcsize(1), lcsize(2)))          !2D! allocate(entry%local_values_r4(1, lcsize))
     entry%local_values_r4 = 0._real32
  else
     if (mype==0) write(*,*) 'not supported output accuracy:',accuracy,'for',trim(name)
     call par_ex
     stop
  endif ! accuracy

  entry%ndim=2                            !2D! entry%ndim=1
  entry%glsize=glsize                     !2D! entry%glsize=(/1, glsize/)
  entry%name = name
  entry%description = description
  entry%units = units
  entry%filename = ""

  entry%dimname(1)=mesh_dimname_from_dimsize(glsize(1), mesh)     !2D! mesh_dimname_from_dimsize(glsize, mesh)
  entry%dimname(2)=mesh_dimname_from_dimsize(glsize(2), mesh)     !2D! entry%dimname(2)='unknown'
  entry%freq=freq
  entry%freq_unit=freq_unit
  entry%addcounter   = 0
  entry%is_in_use=.true.
  io_NSTREAMS=io_NSTREAMS+1

  if(entry%glsize(1)==mesh%nod2D  .or. entry%glsize(2)==mesh%nod2D) then
    entry%is_elem_based = .false.
  else if(entry%glsize(1)==mesh%elem2D .or. entry%glsize(2)==mesh%elem2D) then
    entry%is_elem_based = .true.
  else
    if(mype == 0) print *,"can not determine if ",trim(name)," is node or elem based"
    stop
  end if
end subroutine
!
!--------------------------------------------------------------------------------------------
!
subroutine def_stream2D(glsize, lcsize, name, description, units, data, freq, freq_unit, accuracy, mesh)
  use mod_mesh
  use g_PARSUP
  implicit none
  integer,               intent(in)    :: glsize, lcsize
  character(len=*),      intent(in)    :: name, description, units
  real(kind=WP), target, intent(inout) :: data(:)
  integer,               intent(in)    :: freq
  character,             intent(in)    :: freq_unit
  integer,               intent(in)    :: accuracy
  type(Meandata),        allocatable   :: tmparr(:)
  type(Meandata),        pointer       :: entry
  type(t_mesh), intent(in), target     :: mesh
  integer i
  
  do i = 1, rank(data)
    if ((ubound(data, dim = i)<=0)) then
      if (mype==0) then
        write(*,*) 'WARNING: adding I/O stream for ', trim(name), ' failed (contains 0 dimension)'
        write(*,*) 'upper bound is: ', ubound(data, dim = i)
      end if
      return
    end if    
  end do

  if (mype==0) then
     write(*,*) 'adding I/O stream for ', trim(name)
  end if

   ! add this instance to io_stream array
  if ( .not. allocated(io_stream)) then
    allocate(io_stream(1))
  else
    allocate(tmparr(size(io_stream)+1))
    tmparr(1:size(io_stream)) = io_stream
    deallocate(io_stream)
    call move_alloc(tmparr, io_stream)
  end if
  entry=>io_stream(size(io_stream))
  entry%ptr3(1:1,1:size(data)) => data
  entry%accuracy = accuracy
  
  if (accuracy == i_real8) then
    allocate(data_strategy_nf_double_type :: entry%data_strategy)
     allocate(entry%local_values_r8(1, lcsize))
     entry%local_values_r8 = 0._real64 
  elseif (accuracy == i_real4) then
    allocate(data_strategy_nf_float_type :: entry%data_strategy)
     allocate(entry%local_values_r4(1, lcsize))
     entry%local_values_r4 = 0._real32
  else
     if (mype==0) write(*,*) 'not supported output accuracy:',accuracy,'for',trim(name)
     call par_ex
     stop
  endif ! accuracy

  entry%ndim=1
  entry%glsize=(/1, glsize/)
  entry%name = name
  entry%description = description
  entry%units = units
  entry%filename = ""

  entry%dimname(1)=mesh_dimname_from_dimsize(glsize, mesh)
  entry%dimname(2)='unknown'
  entry%freq=freq
  entry%freq_unit=freq_unit
  entry%addcounter   = 0
  entry%is_in_use=.true.
  io_NSTREAMS=io_NSTREAMS+1

  if(entry%glsize(1)==mesh%nod2D  .or. entry%glsize(2)==mesh%nod2D) then
    entry%is_elem_based = .false.
  else if(entry%glsize(1)==mesh%elem2D .or. entry%glsize(2)==mesh%elem2D) then
    entry%is_elem_based = .true.
  else
    if(mype == 0) print *,"can not determine if ",trim(name)," is node or elem based"
    stop
  end if
end subroutine


  subroutine assert_nf(status, line)
    integer, intent(in) :: status
    integer, intent(in) :: line
    ! EO args
    include "netcdf.inc" ! old netcdf fortran interface required?
    if(status /= NF_NOERR) then
      print *, "error in line ",line, __FILE__, ' ', trim(nf_strerror(status))
      stop 1
    endif   
  end subroutine

end module

