module io_MEANDATA
  USE MOD_PARTIT
  USE MOD_PARSUP
  use o_PARAM, only : WP
  use, intrinsic :: iso_fortran_env, only: real64, real32
  use io_data_strategy_module
  use async_threads_module

  implicit none
#include "netcdf.inc"
  private
  public :: def_stream2D, def_stream3D, output, finalize_output
!
!--------------------------------------------------------------------------------------------
!
  integer, parameter  :: i_real8=8, i_real4=4

  type Meandata
    private
    type(t_partit), pointer                            :: p_partit
    integer                                            :: ndim
    integer                                            :: glsize(2)
    integer                                            :: accuracy
    real(real64), allocatable, dimension(:,:) :: local_values_r8
    real(real32), allocatable, dimension(:,:) :: local_values_r4
    real(real64), allocatable :: aux_r8(:)
    real(real32), allocatable :: aux_r4(:)
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
    integer                                            :: dimID(2), dimvarID(2), varID
    integer                                            :: freq=1
    character                                          :: freq_unit='m'
    logical                                            :: is_in_use=.false.
    logical :: is_elem_based = .false.
    logical :: flip
    class(data_strategy_type), allocatable :: data_strategy
    integer :: comm
    type(thread_type) thread
    integer :: root_rank = 0
    logical :: thread_running = .false.
    real(real64), allocatable, dimension(:,:) :: local_values_r8_copy
    real(real32), allocatable, dimension(:,:) :: local_values_r4_copy
    real(kind=WP) :: ctime_copy
    integer :: mype_workaround
  contains
    final destructor
  end type  
!
!--------------------------------------------------------------------------------------------
!
  type(Meandata), save, target :: io_stream(150) ! todo: find a way to increase the array withhout move_alloc to keep the derived types in Meandata intact
  integer, save                             :: io_NSTREAMS=0
  real(kind=WP)                             :: ctime !current time in seconds from the beginning of the year
!
!--------------------------------------------------------------------------------------------
!
  integer, save                  :: io_listsize   =0
  logical, save                  :: vec_autorotate=.FALSE.
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


subroutine ini_mean_io(ice, dynamics, tracers, partit, mesh)
  use MOD_MESH
  use MOD_TRACER
  USE MOD_PARTIT
  USE MOD_PARSUP
  USE MOD_DYN
  USE MOD_ICE
  use g_cvmix_tke
  use g_cvmix_idemix
  use g_cvmix_kpp
  use g_cvmix_tidal
  use diagnostics
  implicit none
  integer                   :: i, j
  integer, save             :: nm_io_unit  = 103       ! unit to open namelist file, skip 100-102 for cray
  integer                   :: iost
  integer,dimension(15)     :: sel_forcvar=0
  character(len=10)         :: id_string

  type(t_mesh), intent(in) , target :: mesh
  type(t_partit), intent(inout), target :: partit
  type(t_tracer), intent(in)   , target :: tracers
  type(t_dyn)   , intent(in)   , target :: dynamics
  type(t_ice)   , intent(in)   , target :: ice
  namelist /nml_general / io_listsize, vec_autorotate
  namelist /nml_list    / io_list

#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"

  ! OPEN and read namelist for I/O
  open( unit=nm_io_unit, file='namelist.io', form='formatted', access='sequential', status='old', iostat=iost )
  if (iost == 0) then
  if (mype==0) WRITE(*,*) '     file   : ', 'namelist.io',' open ok'
     else
  if (mype==0) WRITE(*,*) 'ERROR: --> bad opening file   : ', 'namelist.io',' ; iostat=',iost
     call par_ex(partit%MPI_COMM_FESOM, partit%mype)
     stop
  endif
  READ(nm_io_unit, nml=nml_general,  iostat=iost )
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
    call def_stream(nod2D, myDim_nod2D, 'sst',      'sea surface temperature',        'C',      tracers%data(1)%values(1,1:myDim_nod2D), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
CASE ('sss       ')
    call def_stream(nod2D, myDim_nod2D, 'sss',      'sea surface salinity',           'psu',    tracers%data(2)%values(1,1:myDim_nod2D), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
CASE ('ssh       ')
    call def_stream(nod2D, myDim_nod2D, 'ssh',      'sea surface elevation',          'm',      dynamics%eta_n,                     io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
CASE ('vve_5     ')
    call def_stream(nod2D, myDim_nod2D, 'vve_5',    'vertical velocity at 5th level', 'm/s',    dynamics%w(5,:),                 io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)

CASE ('ssh_rhs       ')
    call def_stream(nod2D, myDim_nod2D, 'ssh_rhs',      'ssh rhs',          '?',      dynamics%ssh_rhs,                     io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
CASE ('ssh_rhs_old   ')
    call def_stream(nod2D, myDim_nod2D, 'ssh_rhs_old',      'ssh rhs',          '?',      dynamics%ssh_rhs_old,                     io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)

!___________________________________________________________________________________________________________________________________
! output sea ice 
CASE ('uice      ')
    if (use_ice) then
    call def_stream(nod2D, myDim_nod2D, 'uice',     'ice velocity x',                 'm/s',    ice%uice(:),                     io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    end if
CASE ('vice      ')
    if (use_ice) then
    call def_stream(nod2D, myDim_nod2D, 'vice',     'ice velocity y',                 'm/s',    ice%vice(:),                     io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    end if
CASE ('a_ice     ')
    if (use_ice) then
    call def_stream(nod2D, myDim_nod2D, 'a_ice',    'ice concentration',              '%',      ice%data(1)%values(1:myDim_nod2D),      io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    end if
CASE ('m_ice     ')
    if (use_ice) then
    call def_stream(nod2D, myDim_nod2D, 'm_ice',    'ice height',                     'm',      ice%data(2)%values(1:myDim_nod2D),      io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    end if
CASE ('thdgr     ')
    if (use_ice) then
    call def_stream(nod2D, myDim_nod2D, 'thdgr',    'thermodynamic growth rate ice',    'm/s',    ice%thermo%thdgr(1:myDim_nod2D),      io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    end if
CASE ('thdgrsn   ')
    if (use_ice) then
    call def_stream(nod2D, myDim_nod2D, 'thdgrsn',  'thermodynamic growth rate snow',   'm/s',    ice%thermo%thdgrsn(1:myDim_nod2D),    io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    end if
CASE ('flice     ')
    if (use_ice) then
    call def_stream(nod2D, myDim_nod2D,  'flice',    'flooding growth rate ice',       'm/s',    flice(:),                  io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    end if
CASE ('m_snow    ')
    if (use_ice) then
    call def_stream(nod2D, myDim_nod2D, 'm_snow',   'snow height',                     'm',      ice%data(3)%values(1:myDim_nod2D),     io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    end if
    
!___________________________________________________________________________________________________________________________________
! output mixed layer depth    
CASE ('MLD1      ')
    call def_stream(nod2D, myDim_nod2D, 'MLD1',     'Mixed Layer Depth',               'm',      MLD1(1:myDim_nod2D),       io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
CASE ('MLD2      ')
    call def_stream(nod2D, myDim_nod2D, 'MLD2',     'Mixed Layer Depth',               'm',      MLD2(1:myDim_nod2D),       io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    
!___________________________________________________________________________________________________________________________________
! output surface forcing
CASE ('fh        ')
    call def_stream(nod2D, myDim_nod2D, 'fh',       'heat flux',                       'W',      heat_flux_in(:),           io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
CASE ('fw        ')
    call def_stream(nod2D, myDim_nod2D, 'fw',       'fresh water flux',                'm/s',    water_flux(:),             io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
CASE ('atmice_x  ')
    call def_stream(nod2D, myDim_nod2D, 'atmice_x', 'stress atmice x',                 'N/m2',   ice%stress_atmice_x(:),        io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
CASE ('atmice_y  ')
    call def_stream(nod2D, myDim_nod2D, 'atmice_y', 'stress atmice y',                 'N/m2',   ice%stress_atmice_y(:),        io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
CASE ('atmoce_x  ')
    call def_stream(nod2D, myDim_nod2D, 'atmoce_x', 'stress atmoce x',                 'N/m2',   stress_atmoce_x(:),        io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
CASE ('atmoce_y  ')
    call def_stream(nod2D, myDim_nod2D, 'atmoce_y', 'stress atmoce y',                 'N/m2',   stress_atmoce_y(:),        io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
CASE ('iceoce_x  ')
    call def_stream(nod2D, myDim_nod2D, 'iceoce_x', 'stress iceoce x',                 'N/m2',   ice%stress_iceoce_x(:),        io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
CASE ('iceoce_y  ')
    call def_stream(nod2D, myDim_nod2D, 'iceoce_y', 'stress iceoce y',                 'N/m2',   ice%stress_iceoce_y(:),        io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
CASE ('alpha     ')
    call def_stream(nod2D, myDim_nod2D, 'alpha',    'thermal expansion',               'none',   sw_alpha(1,:),             io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
CASE ('beta      ')
    call def_stream(nod2D, myDim_nod2D, 'beta',     'saline contraction',              'none',   sw_beta (1,:),             io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
CASE ('dens_flux ')
    call def_stream(nod2D, myDim_nod2D , 'dflux',   'density flux',               'kg/(m3*s)',   dens_flux(:),              io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
CASE ('runoff    ')
    sel_forcvar(10)= 1
    call def_stream(nod2D, myDim_nod2D, 'runoff',   'river runoff',                    'none',   runoff(:),                 io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
CASE ('evap      ')
    sel_forcvar(7) = 1
    call def_stream(nod2D, myDim_nod2D, 'evap',     'evaporation',                     'm/s',    evaporation(:),            io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
CASE ('prec      ')
    sel_forcvar(5) = 1
    call def_stream(nod2D, myDim_nod2D, 'prec',     'precicipation rain',              'm/s',    prec_rain(:),              io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
CASE ('snow      ')
    sel_forcvar(6) = 1
    call def_stream(nod2D, myDim_nod2D, 'snow',     'precicipation snow',              'm/s',    prec_snow(:),              io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
CASE ('tair      ')
    sel_forcvar(3) = 1
    call def_stream(nod2D, myDim_nod2D, 'tair',     'surface air temperature',         '°C',     Tair(:),                   io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
CASE ('shum      ')
    sel_forcvar(4) = 1
    call def_stream(nod2D, myDim_nod2D, 'shum',     'specific humidity',               '',       shum(:),                   io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
CASE ('swr ')
    sel_forcvar(8) = 1
    call def_stream(nod2D, myDim_nod2D, 'swr',      'short wave radiation',            'W/m^2',  shortwave(:),              io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
CASE ('lwr ')
    sel_forcvar(9) = 1
    call def_stream(nod2D, myDim_nod2D, 'lwr',      'long wave radiation',             'W/m^2',  longwave(:),               io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
CASE ('uwind ')
    sel_forcvar(1) = 1
    call def_stream(nod2D, myDim_nod2D, 'uwind',    '10m zonal surface wind velocity', 'm/s',    u_wind(:),               io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
CASE ('vwind ')
    sel_forcvar(2) = 1
    call def_stream(nod2D, myDim_nod2D, 'vwind',    '10m merid. surface wind velocity','m/s',    v_wind(:),               io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)

    
!___________________________________________________________________________________________________________________________________
! output KPP vertical mixing schemes
CASE ('kpp_obldepth   ')
    if     (mix_scheme_nmb==1 .or. mix_scheme_nmb==17) then! fesom KPP
        call def_stream(nod2D, myDim_nod2D,    'kpp_obldepth',    'KPP ocean boundary layer depth', 'm',   hbl(:),          io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    elseif (mix_scheme_nmb==3 .or. mix_scheme_nmb==37) then ! cvmix KPP
        call def_stream(nod2D, myDim_nod2D,    'kpp_obldepth',    'KPP ocean boundary layer depth', 'm',   kpp_obldepth(:), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    end if
CASE ('kpp_sbuoyflx')
    if     (mix_scheme_nmb==1 .or. mix_scheme_nmb==17) then ! fesom KPP
        call def_stream(nod2D, myDim_nod2D,    'kpp_sbuoyflx',    'surface buoyancy flux',   'm2/s3',  Bo(:),             io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    elseif (mix_scheme_nmb==3 .or. mix_scheme_nmb==37) then ! cvmix KPP
        call def_stream(nod2D, myDim_nod2D,    'kpp_sbuoyflx',    'surface buoyancy flux',   'm2/s3',  kpp_sbuoyflx(:), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    end if
CASE ('tx_sur    ')
    sel_forcvar(11) = 1
    call def_stream(elem2D, myDim_elem2D,  'tx_sur',    'zonal wind str. to ocean',       'm/s2',   stress_surf(1, :),         io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
CASE ('ty_sur    ')
    sel_forcvar(12) = 1
    call def_stream(elem2D, myDim_elem2D,  'ty_sur',    'meridional wind str. to ocean',  'm/s2',   stress_surf(2, :),         io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
CASE ('curl_surf ')
    if (lcurt_stress_surf) then
    call def_stream(nod2D, myDim_nod2D,    'curl_surf', 'vorticity of the surface stress','none',   curl_stress_surf(:),       io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    
    end if
!___________________________________________________________________________________________________________________________________
! output Ferrari/GM parameterisation 2D  
CASE ('fer_C     ')
    if (Fer_GM) then
    call def_stream(nod2D,  myDim_nod2D,   'fer_C',     'GM,   depth independent speed',  'm/s' ,   fer_c(:),                  io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    end if
    
!___________________________________________________________________________________________________________________________________    
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>   3D streams   <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!___________________________________________________________________________________________________________________________________
CASE ('temp      ')
    call def_stream((/nl-1, nod2D/),  (/nl-1, myDim_nod2D/),  'temp',      'temperature', 'C',      tracers%data(1)%values(:,:),             io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
CASE ('salt      ')
    call def_stream((/nl-1, nod2D/),  (/nl-1, myDim_nod2D/),  'salt',      'salinity',    'psu',    tracers%data(2)%values(:,:),             io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
CASE ('otracers  ')
    do j=3, tracers%num_tracers
    write (id_string, "(I3.3)") tracers%data(j)%ID
    call def_stream((/nl-1, nod2D/),  (/nl-1, myDim_nod2D/),  'tra_'//id_string, 'pasive tracer ID='//id_string, 'n/a', tracers%data(j)%values(:,:), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    end do
CASE ('slope_x   ')
    call def_stream((/nl-1,  nod2D/), (/nl-1, myDim_nod2D/),  'slope_x',   'neutral slope X',    'none', slope_tapered(1,:,:), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
CASE ('slope_y   ')
    call def_stream((/nl-1,  nod2D/), (/nl-1, myDim_nod2D/),  'slope_y',   'neutral slope Y',    'none', slope_tapered(2,:,:), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
CASE ('slope_z   ')
    call def_stream((/nl-1,  nod2D/), (/nl-1, myDim_nod2D/),  'slope_z',   'neutral slope Z',    'none', slope_tapered(3,:,:), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
CASE ('N2        ')
    call def_stream((/nl,    nod2D/), (/nl,   myDim_nod2D/),  'N2',        'brunt väisälä',      '1/s2', bvfreq(:,:),          io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
CASE ('Kv        ')
    call def_stream((/nl,    nod2D/), (/nl,   myDim_nod2D/),  'Kv',        'vertical diffusivity Kv',  'm2/s', Kv(:,:),        io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
CASE ('u         ')
    call def_stream((/nl-1, elem2D/), (/nl-1, myDim_elem2D/), 'u',         'horizontal velocity','m/s',           dynamics%uv(1,:,:),     io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
CASE ('v         ')
    call def_stream((/nl-1, elem2D/), (/nl-1, myDim_elem2D/), 'v',         'meridional velocity','m/s',           dynamics%uv(2,:,:),     io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
CASE ('unod      ')
    call def_stream((/nl-1, nod2D/),  (/nl-1,   myDim_nod2D/),'unod',      'horizontal velocity at nodes', 'm/s', dynamics%uvnode(1,:,:), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
CASE ('vnod      ')
    call def_stream((/nl-1, nod2D/),  (/nl-1,   myDim_nod2D/),'vnod',      'meridional velocity at nodes', 'm/s', dynamics%uvnode(2,:,:), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
CASE ('w         ')
    call def_stream((/nl,    nod2D/), (/nl,   myDim_nod2D/),  'w',         'vertical velocity',  'm/s',           dynamics%w(:,:),        io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
CASE ('Av        ')
    call def_stream((/nl,   elem2D/), (/nl,   myDim_elem2D/), 'Av',        'vertical viscosity Av',  'm2/s',      Av(:,:),                io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
CASE ('u_dis_tend')
    if(dynamics%opt_visc==8) then
    call def_stream((/nl-1, elem2D/), (/nl-1, myDim_elem2D/), 'u_dis_tend',    'horizontal velocity viscosity tendency', 'm/s', UV_dis_tend(1,:,:), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    end if
CASE ('v_dis_tend')
    if(dynamics%opt_visc==8) then
    call def_stream((/nl-1, elem2D/), (/nl-1, myDim_elem2D/), 'v_dis_tend',    'meridional velocity viscosity tendency', 'm/s', UV_dis_tend(2,:,:), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh) 
    end if
CASE ('u_back_tend')
    if(dynamics%opt_visc==8) then    
    call def_stream((/nl-1, elem2D/), (/nl-1, myDim_elem2D/), 'u_back_tend',    'horizontal velocity backscatter tendency', 'm2/s2', UV_back_tend(1,:,:), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    end if
CASE ('v_back_tend') 
    if(dynamics%opt_visc==8) then
    call def_stream((/nl-1, elem2D/), (/nl-1, myDim_elem2D/), 'v_back_tend',    'meridional velocity backscatter tendency', 'm2/s2', UV_back_tend(2,:,:), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh) 
    end if
CASE ('u_total_tend')
    if(dynamics%opt_visc==8) then
    call def_stream((/nl-1, elem2D/), (/nl-1, myDim_elem2D/), 'u_total_tend',    'horizontal velocity total viscosity tendency', 'm/s', UV_total_tend(1,:,:), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    end if
CASE ('v_total_tend')
    if(dynamics%opt_visc==8) then
    call def_stream((/nl-1, elem2D/), (/nl-1, myDim_elem2D/), 'v_total_tend',    'meridional velocity total viscosity tendency', 'm/s', UV_total_tend(2,:,:), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh) 
    end if
!___________________________________________________________________________________________________________________________________
! output Ferrari/GM parameterisation
CASE ('bolus_u   ')
    if (Fer_GM) then
    call def_stream((/nl-1, elem2D/), (/nl-1, myDim_elem2D/), 'bolus_u',   'GM bolus velocity U','m/s',  dynamics%fer_uv(1,:,:),        io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    end if
CASE ('bolus_v   ')
    if (Fer_GM) then  
    call def_stream((/nl-1, elem2D/), (/nl-1, myDim_elem2D/), 'bolus_v',   'GM bolus velocity V','m/s',  dynamics%fer_uv(2,:,:),        io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    end if
CASE ('bolus_w   ')
    if (Fer_GM) then
    call def_stream((/nl  , nod2D /), (/nl,   myDim_nod2D /), 'bolus_w',   'GM bolus velocity W','m/s',  dynamics%fer_w(:,:),        io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    end if
CASE ('fer_K     ')
    if (Fer_GM) then
    call def_stream((/nl  , nod2D /), (/nl,   myDim_nod2D /), 'fer_K',     'GM, stirring diff.','m2/s',  fer_k(:,:),           io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    end if
CASE ('fer_scal  ')
    if (Fer_GM) then
    call def_stream(        nod2D   ,         myDim_nod2D   , 'fer_scal',  'GM surface scaling','',  fer_scal(:),           io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    end if
CASE ('dMOC      ')
    if (ldiag_dMOC) then
       call def_stream((/std_dens_N, elem2D/),  (/std_dens_N, myDim_elem2D/), 'U_rho_x_DZ',     'fluxes for density MOC', 'fluxes', std_dens_UVDZ(1,:,:),   1, 'y', i_real4, partit, mesh)
       call def_stream((/std_dens_N, elem2D/),  (/std_dens_N, myDim_elem2D/), 'V_rho_x_DZ',     'fluxes for density MOC', 'fluxes', std_dens_UVDZ(2,:,:),   1, 'y', i_real4, partit, mesh)
       call def_stream((/std_dens_N, elem2D/),  (/std_dens_N, myDim_elem2D/), 'std_heat_flux',  'HF bouyancy flux      ', 'kg*m/s' ,std_dens_flux(1,:,:),   1, 'y', i_real4, partit, mesh)
       call def_stream((/std_dens_N, elem2D/),  (/std_dens_N, myDim_elem2D/), 'std_rest_flux',  'RESTOR. bouyancy flux ', 'kg*m/s' ,std_dens_flux(2,:,:),   1, 'y', i_real4, partit, mesh)
       call def_stream((/std_dens_N, elem2D/),  (/std_dens_N, myDim_elem2D/), 'std_frwt_flux',  'FW bouyancy flux      ', 'kg*m/s' ,std_dens_flux(3,:,:),   1, 'y', i_real4, partit, mesh)
       call def_stream((/std_dens_N, elem2D/),  (/std_dens_N, myDim_elem2D/), 'std_dens_dVdT',  'dV/dT',                  'm3/s'   ,std_dens_dVdT(:,:),     1, 'y', i_real4, partit, mesh)
       call def_stream((/std_dens_N, nod2D /),  (/std_dens_N,  myDim_nod2D/), 'std_dens_DIV',   'm3/s',                   'm3/s'   ,std_dens_DIV(:,:),      1, 'y', i_real4, partit, mesh)
       call def_stream((/std_dens_N, elem2D/),  (/std_dens_N, myDim_elem2D/), 'std_dens_Z',     'm',                      'm'      ,std_dens_Z(:,:),        1, 'y', i_real4, partit, mesh)
       call def_stream((/nl-1,       nod2D /),  (/nl-1,       myDim_nod2D /), 'density_dMOC',   'density'               , 'm',      density_dmoc(:,:),      1, 'y', i_real4, partit, mesh)
       call def_stream(elem2D,                                myDim_elem2D  , 'density_flux_e', 'density flux at elems ', 'm',      dens_flux_e(:),         1, 'y', i_real4, partit, mesh)
    end if
!___________________________________________________________________________________________________________________________________
CASE ('pgf_x     ')    
    call def_stream((/nl-1, elem2D/), (/nl-1, myDim_elem2D/), 'pgf_x', 'zonal pressure gradient force'     , 'm/s^2', pgf_x(:,:), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
CASE ('pgf_y     ')    
    call def_stream((/nl-1, elem2D/), (/nl-1, myDim_elem2D/), 'pgf_y', 'meridional pressure gradient force', 'm/s^2', pgf_y(:,:), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
!___________________________________________________________________________________________________________________________________    

#if defined (__oifs)
CASE ('alb       ')
  call def_stream(nod2D, myDim_nod2D, 'alb',    'ice albedo',              'none',   ice%atmcoupl%ice_alb(:),                   io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
CASE ('ist       ')
  call def_stream(nod2D, myDim_nod2D, 'ist',    'ice surface temperature', 'K',      ice%data(4)%values(:),                  io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
CASE ('qsi       ')
  call def_stream(nod2D, myDim_nod2D, 'qsi',    'ice heat flux',           'W/m^2',  ice%atmcoupl%ice_flx_h(:),             io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
CASE ('qso       ')
  call def_stream(nod2D, myDim_nod2D, 'qso',    'oce heat flux',           'W/m^2',  ice%atmcoupl%oce_flx_h(:),             io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
#endif
!___________________________________________________________________________________________________________________________________


CASE DEFAULT
    if (mype==0) write(*,*) 'stream ', io_list(i)%id, ' is not defined !'
END SELECT
END DO


!3D
  if (ldiag_energy) then
     call def_stream((/nl,   nod2D/),  (/nl,     myDim_nod2D/), 'rhof',     'in-situ density at faces',    'kg/m3',     rhof(:,:),  1, 'm', i_real8, partit, mesh)
     call def_stream((/nl,   nod2D/),  (/nl,     myDim_nod2D/), 'wrhof',    'vertical velocity x density', 'kg/(s*m2)', wrhof(:,:), 1, 'm', i_real8, partit, mesh)
     call def_stream((/nl-1, nod2D/),  (/nl-1,   myDim_nod2D/), 'uu',   'u times u', 'm2/s2', u_x_u(:,:), 1, 'm', i_real8, partit, mesh)
     call def_stream((/nl-1, nod2D/),  (/nl-1,   myDim_nod2D/), 'uv',   'u times v', 'm2/s2', u_x_v(:,:), 1, 'm', i_real8, partit, mesh)
     call def_stream((/nl-1, nod2D/),  (/nl-1,   myDim_nod2D/), 'vv',   'v times v', 'm2/s2', v_x_v(:,:), 1, 'm', i_real8, partit, mesh)
     call def_stream((/nl,  elem2D/),  (/nl-1,   myDim_elem2D/),'uw',   'u times w', 'm2/s2', u_x_w(:,:), 1, 'm', i_real8, partit, mesh)
     call def_stream((/nl,  elem2D/),  (/nl-1,   myDim_elem2D/),'vw',   'v times w', 'm2/s2', v_x_w(:,:), 1, 'm', i_real8, partit, mesh)
     call def_stream((/nl-1, nod2D/),  (/nl-1,   myDim_nod2D/), 'dudx', 'du/dx',     '1/s',   dudx(:,:),  1, 'm', i_real8, partit, mesh)
     call def_stream((/nl-1, nod2D/),  (/nl-1,   myDim_nod2D/), 'dudy', 'du/dy',     '1/s',   dudy(:,:),  1, 'm', i_real8, partit, mesh)
     call def_stream((/nl-1, nod2D/),  (/nl-1,   myDim_nod2D/), 'dvdx', 'dv/dx',     '1/s',   dvdx(:,:),  1, 'm', i_real8, partit, mesh)
     call def_stream((/nl-1, nod2D/),  (/nl-1,   myDim_nod2D/), 'dvdy', 'dv/dy',     '1/s',   dvdy(:,:),  1, 'm', i_real8, partit, mesh)
     call def_stream((/nl,  elem2D/),  (/nl,     myDim_elem2D/), 'dudz', 'du/dz',     '1/s',   dudz(:,:),  1, 'm', i_real8, partit, mesh)
     call def_stream((/nl,  elem2D/),  (/nl,     myDim_elem2D/), 'dvdz', 'dv/dz',     '1/s',   dvdz(:,:),  1, 'm', i_real8, partit, mesh)
     call def_stream((/nl,  elem2D/),  (/nl,     myDim_elem2D/), 'av_dudz', 'int(Av * du/dz)',        'm3/s2',   av_dudz(:,:),    1, 'm', i_real4, partit, mesh)
     call def_stream((/nl,  elem2D/),  (/nl,     myDim_elem2D/), 'av_dvdz', 'int(Av * dv/dz)',        'm3/s2',   av_dvdz(:,:),    1, 'm', i_real4, partit, mesh)
     call def_stream((/nl,  elem2D/),  (/nl,     myDim_elem2D/), 'av_dudz_sq',  'Av * (du/dz)^2',     'm^2/s^3', av_dudz_sq(:,:), 1, 'm', i_real4, partit, mesh)
     call def_stream((/nl,   elem2D/), (/nl,     myDim_elem2D/), 'Av',    'Vertical mixing A',         'm2/s',            Av(:,:), 1, 'm', i_real4, partit, mesh)
   
     call def_stream((/nl-1, elem2D/), (/nl-1,   myDim_elem2D/), 'um',  'horizontal velocity', 'm/s', dynamics%uv(1,:,:),     1, 'm', i_real4, partit, mesh)
     call def_stream((/nl-1, elem2D/), (/nl-1,   myDim_elem2D/), 'vm',  'meridional velocity', 'm/s', dynamics%uv(2,:,:),     1, 'm', i_real4, partit, mesh)
     call def_stream((/nl, nod2D/),    (/nl,     myDim_nod2D/),  'wm',  'vertical velocity',   'm/s', dynamics%w(:,:),     1, 'm', i_real8, partit, mesh)

     call def_stream(elem2D, myDim_elem2D,   'utau_surf',  '(u, tau) at the surface',      'N/(m s)', utau_surf(1:myDim_elem2D),     1, 'm', i_real4, partit, mesh)
     call def_stream(elem2D, myDim_elem2D,   'utau_bott',  '(u, tau) at the bottom',       'N/(m s)', utau_bott(1:myDim_elem2D),     1, 'm', i_real4, partit, mesh)
     call def_stream(elem2D, myDim_elem2D,   'u_bott',     'bottom velocity',                  'm/s', u_bott(1:myDim_elem2D),        1, 'm', i_real4, partit, mesh)
     call def_stream(elem2D, myDim_elem2D,   'v_bott',     'bottom velocity',                  'm/s', v_bott(1:myDim_elem2D),        1, 'm', i_real4, partit, mesh)
     call def_stream(elem2D, myDim_elem2D,   'u_surf',     'surface velocity',                 'm/s', u_surf(1:myDim_elem2D),        1, 'm', i_real4, partit, mesh)
     call def_stream(elem2D, myDim_elem2D,   'v_surf',     'surface velocity',                 'm/s', u_surf(1:myDim_elem2D),        1, 'm', i_real4, partit, mesh)
     call def_stream(elem2D, myDim_elem2D,   'tx_bot',     'bottom stress x',                 'N/m2', stress_bott(1, 1:myDim_elem2D),1, 'm', i_real4, partit, mesh)
     call def_stream(elem2D, myDim_elem2D,   'ty_bot',     'bottom stress y',                 'N/m2', stress_bott(2, 1:myDim_elem2D),1, 'm', i_real4, partit, mesh)
     if (sel_forcvar(11)==0) call def_stream(elem2D, myDim_elem2D,   'tx_sur',     'zonal wind stress to ocean',     'm/s2', stress_surf(1, 1:myDim_elem2D),1, 'm', i_real4, partit, mesh) ; sel_forcvar(11)=1
     if (sel_forcvar(12)==0) call def_stream(elem2D, myDim_elem2D,   'ty_sur',     'meridional wind stress to ocean','m/s2', stress_surf(2, 1:myDim_elem2D),1, 'm', i_real4, partit, mesh) ; sel_forcvar(12)=1
  end if

    
    if (mix_scheme_nmb==5 .or. mix_scheme_nmb==56) then
        ! TKE diagnostic 
        call def_stream((/nl,nod2D/), (/nl,myDim_nod2D/), 'tke'     , 'turbulent kinetic energy'                    , 'm^2/s^2', tke(:,:)     , 1, 'y', i_real4, partit, mesh)
        call def_stream((/nl,nod2D/), (/nl,myDim_nod2D/), 'tke_Ttot', 'total production of turbulent kinetic energy', 'm^2/s^3', tke_Ttot(:,:), 1, 'y', i_real4, partit, mesh)
        call def_stream((/nl,nod2D/), (/nl,myDim_nod2D/), 'tke_Tbpr', 'TKE production by buoyancy'                  , 'm^2/s^3', tke_Tbpr(:,:), 1, 'y', i_real4, partit, mesh)
        call def_stream((/nl,nod2D/), (/nl,myDim_nod2D/), 'tke_Tspr', 'TKE production by shear'                     , 'm^2/s^3', tke_Tspr(:,:), 1, 'y', i_real4, partit, mesh)
        call def_stream((/nl,nod2D/), (/nl,myDim_nod2D/), 'tke_Tdif', 'TKE production by vertical diffusion'        , 'm^2/s^3', tke_Tdif(:,:), 1, 'y', i_real4, partit, mesh)
        call def_stream((/nl,nod2D/), (/nl,myDim_nod2D/), 'tke_Tdis', 'TKE production by dissipation'               , 'm^2/s^3', tke_Tdis(:,:), 1, 'y', i_real4, partit, mesh)
        call def_stream((/nl,nod2D/), (/nl,myDim_nod2D/), 'tke_Twin', 'TKE production by wind'                      , 'm^2/s^3', tke_Twin(:,:), 1, 'y', i_real4, partit, mesh)
        call def_stream((/nl,nod2D/), (/nl,myDim_nod2D/), 'tke_Tbck', 'background forcing for TKE'                  , 'm^2/s^3', tke_Tbck(:,:), 1, 'y', i_real4, partit, mesh)
        call def_stream((/nl,nod2D/), (/nl,myDim_nod2D/), 'tke_Lmix', 'mixing length scale of TKE'                  , 'm'      , tke_Lmix(:,:), 1, 'y', i_real4, partit, mesh)
        call def_stream((/nl,nod2D/), (/nl,myDim_nod2D/), 'tke_Pr'  , 'Prantl number'                               , ''       , tke_Pr(:,:)  , 1, 'y', i_real4, partit, mesh)
        if (mix_scheme_nmb==56) then
            ! TKE-IDEMIX diagnostic 
            call def_stream((/nl,nod2D/), (/nl,myDim_nod2D/), 'tke_Tiwf', 'TKE production by internal waves (IDEMIX)', 'm^2/s^3', tke_Tiwf(:,:), 1, 'y', i_real4, partit, mesh)
        end if 
    end if 
    
    if (mod(mix_scheme_nmb,10)==6) then
        ! IDEMIX Internal-Wave-Energy diagnostics
        call def_stream((/nl,nod2D/), (/nl,myDim_nod2D/), 'iwe'     , 'internal wave energy'                    , 'm^2/s^2', iwe(:,:)     , 1, 'y', i_real4, partit, mesh)
        call def_stream((/nl,nod2D/), (/nl,myDim_nod2D/), 'iwe_Ttot', 'total production of internal wave energy', 'm^2/s^2', iwe_Ttot(:,:), 1, 'y', i_real4, partit, mesh)
        call def_stream((/nl,nod2D/), (/nl,myDim_nod2D/), 'iwe_Tdif', 'IWE production by vertical diffusion'    , 'm^2/s^3', iwe_Tdif(:,:), 1, 'y', i_real4, partit, mesh)
        call def_stream((/nl,nod2D/), (/nl,myDim_nod2D/), 'iwe_Tdis', 'IWE production by dissipation'           , 'm^2/s^3', iwe_Tdis(:,:), 1, 'y', i_real4, partit, mesh)
        call def_stream((/nl,nod2D/), (/nl,myDim_nod2D/), 'iwe_Tsur', 'IWE production from surface forcing'     , 'm^2/s^2', iwe_Tsur(:,:), 1, 'y', i_real4, partit, mesh)
        call def_stream((/nl,nod2D/), (/nl,myDim_nod2D/), 'iwe_Tbot', 'IWE production from bottom forcing'      , 'm^2/s^2', iwe_Tbot(:,:), 1, 'y', i_real4, partit, mesh)
        call def_stream((/nl,nod2D/), (/nl,myDim_nod2D/), 'iwe_c0'  , 'IWE vertical group velocity'             , 'm/s'    , iwe_c0(:,:)  , 1, 'y', i_real4, partit, mesh)
        call def_stream((/nl,nod2D/), (/nl,myDim_nod2D/), 'iwe_v0'  , 'IWE horizontal group velocity'           , 'm/s'    , iwe_c0(:,:)  , 1, 'y', i_real4, partit, mesh)
    end if
    
    if (mod(mix_scheme_nmb,10)==7) then
        ! cvmix_TIDAL diagnostics
        call def_stream((/nl,nod2D/), (/nl,myDim_nod2D/), 'tidal_Kv'  , 'tidal diffusivity'           , 'm^2/s'    , tidal_Kv(:,:)  , 1, 'y', i_real4, partit, mesh)
        call def_stream((/nl,nod2D/), (/nl,myDim_nod2D/), 'tidal_Av'  , 'tidal viscosity'             , 'm^2/s'    , tidal_Av(:,:)  , 1, 'y', i_real4, partit, mesh)
        call def_stream(     nod2D  ,      myDim_nod2D  , 'tidal_forcbot', 'near tidal bottom forcing', 'W/m^2'    , tidal_forc_bottom_2D  , 100, 'y', i_real4, partit, mesh)
    end if
    
  !___________________________________________________________________________________________________________________________________
  ! output Redi parameterisation
  if (Redi) then
     call def_stream((/nl-1  , nod2D /), (/nl-1,   myDim_nod2D /), 'Redi_K',   'Redi diffusion coefficient', 'm2/s', Ki(:,:),    1, 'y', i_real4, partit, mesh)
  end if

  !___________________________________________________________________________________________________________________________________
  ! output Monin-Obukov (TB04) mixing length
  if (use_momix) then
     call def_stream(nod2D, myDim_nod2D, 'momix_length',   'Monin-Obukov mixing length', 'm', mixlength(:),    1, 'm', i_real4, partit, mesh)
  end if
  
    !___________________________________________________________________________________________________________________________________
    if (ldiag_curl_vel3) then
        call def_stream((/nl-1, nod2D/),  (/nl-1, myDim_nod2D/),  'curl_u',     'relative vorticity',          '1/s',   vorticity,                   1, 'm', i_real4, partit, mesh)
    end if

    !___________________________________________________________________________________________________________________________________
    if (ice%whichEVP==1) then
    end if
    
    if (ice%whichEVP==2) then
        call def_stream(elem2D, myDim_elem2D, 'alpha_EVP', 'alpha in EVP', 'n/a', ice%alpha_evp_array,  1, 'd', i_real4, partit, mesh)
        call def_stream(nod2D,  myDim_nod2D,  'beta_EVP',  'beta in EVP',  'n/a', ice%beta_evp_array,   1, 'd', i_real4, partit, mesh)
    end if
  
    !___________________________________________________________________________
    if (ldiag_dvd) then
        call def_stream((/nl-1, nod2D/), (/nl-1, myDim_nod2D/), 'dvd_temp_h', 'horiz. dvd of temperature', '°C/s' , tracers%work%tr_dvd_horiz(:,:,1), 1, 'm', i_real4, partit, mesh)
        call def_stream((/nl-1, nod2D/), (/nl-1, myDim_nod2D/), 'dvd_temp_v', 'vert. dvd of temperature' , '°C/s' , tracers%work%tr_dvd_vert(:,:,1) , 1, 'm', i_real4, partit, mesh)
        call def_stream((/nl-1, nod2D/), (/nl-1, myDim_nod2D/), 'dvd_salt_h', 'horiz. dvd of salinity'   , 'psu/s', tracers%work%tr_dvd_horiz(:,:,2), 1, 'm', i_real4, partit, mesh)
        call def_stream((/nl-1, nod2D/), (/nl-1, myDim_nod2D/), 'dvd_salt_v', 'vert. dvd of salinity'    , 'psu/s', tracers%work%tr_dvd_vert(:,:,2) , 1, 'm', i_real4, partit, mesh)
    end if 
    
    !___________________________________________________________________________
    if (ldiag_forc) then
        if (sel_forcvar( 1)==0) call def_stream(nod2D , myDim_nod2D , 'uwind' , '10m zonal surface wind velocity', 'm/s'  , u_wind(:)        , 1, 'm', i_real4, partit, mesh)
        if (sel_forcvar( 2)==0) call def_stream(nod2D , myDim_nod2D , 'vwind' , '10m merid surface wind velocity', 'm/s'  , v_wind(:)        , 1, 'm', i_real4, partit, mesh)
        if (sel_forcvar( 3)==0) call def_stream(nod2D , myDim_nod2D , 'tair'  , 'surface air temperature'        , '°C'   , Tair(:)          , 1, 'm', i_real4, partit, mesh)
        if (sel_forcvar( 4)==0) call def_stream(nod2D , myDim_nod2D , 'shum'  , 'specific humidity'              , ''     , shum(:)          , 1, 'm', i_real4, partit, mesh)
        if (sel_forcvar( 5)==0) call def_stream(nod2D , myDim_nod2D , 'prec'  , 'precicipation rain'             , 'm/s'  , prec_rain(:)     , 1, 'm', i_real4, partit, mesh)
        if (sel_forcvar( 6)==0) call def_stream(nod2D , myDim_nod2D , 'snow'  , 'precicipation snow'             , 'm/s'  , prec_snow(:)     , 1, 'm', i_real4, partit, mesh)
        if (sel_forcvar( 7)==0) call def_stream(nod2D , myDim_nod2D , 'evap'  , 'evaporation'                    , 'm/s'  , evaporation(:)   , 1, 'm', i_real4, partit, mesh)
        if (sel_forcvar( 8)==0) call def_stream(nod2D , myDim_nod2D , 'swr'   , 'short wave radiation'           , 'W/m^2', shortwave(:)     , 1, 'm', i_real4, partit, mesh)
        if (sel_forcvar( 9)==0) call def_stream(nod2D , myDim_nod2D , 'lwr'   , 'long wave radiation'            , 'W/m^2', longwave(:)      , 1, 'm', i_real4, partit, mesh)
        if (sel_forcvar(10)==0) call def_stream(nod2D , myDim_nod2D , 'runoff', 'river runoff'                   , 'none' , runoff(:)        , 1, 'm', i_real4, partit, mesh)
        if (sel_forcvar(11)==0) call def_stream(elem2D, myDim_elem2D, 'tx_sur', 'zonal wind str. to ocean'       , 'm/s^2', stress_surf(1, :), 1, 'm', i_real4, partit, mesh)
        if (sel_forcvar(12)==0) call def_stream(elem2D, myDim_elem2D, 'ty_sur', 'meridional wind str. to ocean'  , 'm/s^2', stress_surf(2, :), 1, 'm', i_real4, partit, mesh)
        call def_stream(nod2D , myDim_nod2D , 'cd',    'wind drag coef. '             , '',     cd_atm_oce_arr(:), 1, 'm', i_real4, partit, mesh)
        call def_stream(nod2D , myDim_nod2D , 'ch',    'transfer coeff. sensible heat', '',     ch_atm_oce_arr(:), 1, 'm', i_real4, partit, mesh)
        call def_stream(nod2D , myDim_nod2D , 'ce',    'transfer coeff. evaporation ' , '',     ce_atm_oce_arr(:), 1, 'm', i_real4, partit, mesh)
#if defined (__oasis)
        call def_stream(nod2D,  myDim_nod2D,  'subli', 'sublimation',                   'm/s',  sublimation(:),   1, 'm',  i_real4, partit, mesh)
#endif
    end if
    
    
end subroutine
!
!--------------------------------------------------------------------------------------------
!
function mesh_dimname_from_dimsize(size, partit, mesh) result(name)
  use mod_mesh
  USE MOD_PARTIT
  USE MOD_PARSUP
  use diagnostics
#if defined (__icepack)
  use icedrv_main,   only: ncat ! number of ice thickness cathegories
#endif
  implicit none
  integer       :: size
  type(t_mesh)  , intent(in) :: mesh
  type(t_partit), intent(in) :: partit
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
#if defined (__icepack)
  elseif (size==ncat) then
    name='ncat'
#endif
  else
    name='unknown'
    if (partit%mype==0) write(*,*) 'WARNING: unknown dimension in mean I/O with size of ', size
  end if
end function
!
!--------------------------------------------------------------------------------------------
!
subroutine create_new_file(entry, ice, dynamics, partit, mesh)
  use g_clock
  use mod_mesh
  USE MOD_PARTIT
  USE MOD_PARSUP
  USE MOD_DYN
  USE MOD_ICE
  use fesom_version_info_module
  use g_config
  use o_PARAM

  implicit none
  character(2000)               :: att_text
  type(t_mesh)  , intent(in) :: mesh
  type(t_partit), intent(in) :: partit
  type(t_dyn)   , intent(in) :: dynamics
  type(t_ice)   , intent(in) :: ice
 
  type(Meandata), intent(inout) :: entry
  character(len=*), parameter :: global_attributes_prefix = "FESOM_"
  ! Serial output implemented so far
  if (partit%mype/=entry%root_rank) return
  ! create an ocean output file
  write(*,*) 'initializing I/O file for ', trim(entry%name)

  call assert_nf( nf_create(entry%filename, IOR(NF_NOCLOBBER,IOR(NF_NETCDF4,NF_CLASSIC_MODEL)), entry%ncid), __LINE__)

!___Create mesh related dimensions__________________________________________
  if (entry%ndim==1) then
     call assert_nf( nf_def_dim(entry%ncid, entry%dimname(1), entry%glsize(2), entry%dimID(1)), __LINE__)
  else if (entry%ndim==2) then
     call assert_nf( nf_def_dim(entry%ncid,  entry%dimname(1), entry%glsize(1), entry%dimID(1)), __LINE__)
     call assert_nf( nf_def_var(entry%ncid,  entry%dimname(1), NF_DOUBLE,   1,  entry%dimID(1), entry%dimvarID(1)), __LINE__)
     if (entry%dimname(1)=='nz') then
       call assert_nf( nf_put_att_text(entry%ncid, entry%dimvarID(1), 'long_name', len_trim('depth at layer interface'),'depth at layer interface'), __LINE__)
     elseif (entry%dimname(1)=='nz1') then
       call assert_nf( nf_put_att_text(entry%ncid, entry%dimvarID(1), 'long_name', len_trim('depth at layer midpoint'),'depth at layer midpoint'), __LINE__)
     elseif (entry%dimname(1)=='ncat') then
       call assert_nf( nf_put_att_text(entry%ncid, entry%dimvarID(1), 'long_name', len_trim('sea-ice thickness class'),'sea-ice thickness class'), __LINE__)
     else
       if (partit%mype==0) write(*,*) 'WARNING: unknown first dimension in 2d mean I/O data'
     end if 
     call assert_nf( nf_put_att_text(entry%ncid, entry%dimvarID(1), 'units', len_trim('m'),'m'), __LINE__)
     call assert_nf( nf_put_att_text(entry%ncid, entry%dimvarID(1), 'positive', len_trim('down'),'down'), __LINE__)
     call assert_nf( nf_put_att_text(entry%ncid, entry%dimvarID(1), 'axis', len_trim('Z'),'Z'), __LINE__)
     
     call assert_nf( nf_def_dim(entry%ncid, entry%dimname(2), entry%glsize(2), entry%dimID(2)), __LINE__)
  end if
!___Create time related dimensions__________________________________________
  call assert_nf( nf_def_dim(entry%ncid, 'time', NF_UNLIMITED, entry%recID), __LINE__)
!___Define the time and iteration variables_________________________________
  call assert_nf( nf_def_var(entry%ncid, 'time', NF_DOUBLE, 1, entry%recID, entry%tID), __LINE__)
  att_text='time'
  call assert_nf( nf_put_att_text(entry%ncid, entry%tID, 'long_name', len_trim(att_text), trim(att_text)), __LINE__)
  call assert_nf( nf_put_att_text(entry%ncid, entry%tID, 'standard_name', len_trim(att_text), trim(att_text)), __LINE__)
  write(att_text, '(a14,I4.4,a1,I2.2,a1,I2.2,a6)') 'seconds since ', yearold, '-', 1, '-', 1, ' 0:0:0'
  call assert_nf( nf_put_att_text(entry%ncid, entry%tID, 'units', len_trim(att_text), trim(att_text)), __LINE__)
  call assert_nf( nf_put_att_text(entry%ncid, entry%tID, 'axis', len_trim('T'), trim('T')), __LINE__)
  call assert_nf( nf_put_att_text(entry%ncid, entry%tID, 'stored_direction', len_trim('increasing'), trim('increasing')), __LINE__)
  
  call assert_nf( nf_def_var(entry%ncid, trim(entry%name), entry%data_strategy%netcdf_type(), entry%ndim+1, (/entry%dimid(entry%ndim:1:-1), entry%recID/), entry%varID), __LINE__)


  call assert_nf( nf_put_att_text(entry%ncid, entry%varID, 'description', len_trim(entry%description), entry%description), __LINE__)
  call assert_nf( nf_put_att_text(entry%ncid, entry%varID, 'long_name', len_trim(entry%description), entry%description), __LINE__)
  call assert_nf( nf_put_att_text(entry%ncid, entry%varID, 'units',       len_trim(entry%units),       entry%units), __LINE__)
  
 
!___Global attributes________  
  call assert_nf( nf_put_att_text(entry%ncid, NF_GLOBAL, global_attributes_prefix//'model', len_trim('FESOM2'),'FESOM2'), __LINE__)
  call assert_nf( nf_put_att_text(entry%ncid, NF_GLOBAL, global_attributes_prefix//'website', len_trim('fesom.de'), trim('fesom.de')), __LINE__)
 
  call assert_nf( nf_put_att_text(entry%ncid, NF_GLOBAL, global_attributes_prefix//'git_SHA', len_trim(fesom_git_sha()), fesom_git_sha()), __LINE__)
  call assert_nf( nf_put_att_text(entry%ncid, NF_GLOBAL, global_attributes_prefix//'MeshPath', len_trim(MeshPath), trim(MeshPath)), __LINE__)
  call assert_nf( nf_put_att_text(entry%ncid, NF_GLOBAL, global_attributes_prefix//'mesh_representative_checksum', len(mesh%representative_checksum), mesh%representative_checksum), __LINE__)
  call assert_nf( nf_put_att_text(entry%ncid, NF_GLOBAL, global_attributes_prefix//'ClimateDataPath', len_trim(ClimateDataPath), trim(ClimateDataPath)), __LINE__)
  call assert_nf( nf_put_att_text(entry%ncid, NF_GLOBAL, global_attributes_prefix//'which_ALE', len_trim(which_ALE), trim(which_ALE)), __LINE__)
  call assert_nf( nf_put_att_text(entry%ncid, NF_GLOBAL, global_attributes_prefix//'mix_scheme', len_trim(mix_scheme), trim(mix_scheme)), __LINE__)
! call assert_nf( nf_put_att_text(entry%ncid, NF_GLOBAL, global_attributes_prefix//'tra_adv_hor', len_trim(tra_adv_hor), trim(tra_adv_hor)), __LINE__)
! call assert_nf( nf_put_att_text(entry%ncid, NF_GLOBAL, global_attributes_prefix//'tra_adv_ver', len_trim(tra_adv_ver), trim(tra_adv_ver)), __LINE__)
! call assert_nf( nf_put_att_text(entry%ncid, NF_GLOBAL, global_attributes_prefix//'tra_adv_lim', len_trim(tra_adv_lim), trim(tra_adv_lim)), __LINE__)
 
 
  call assert_nf( nf_put_att_int(entry%ncid, NF_GLOBAL, global_attributes_prefix//'use_partial_cell', NF_INT, 1,  use_partial_cell), __LINE__)
  call assert_nf( nf_put_att_int(entry%ncid, NF_GLOBAL, global_attributes_prefix//'force_rotation', NF_INT, 1,  force_rotation), __LINE__)
  call assert_nf( nf_put_att_int(entry%ncid, NF_GLOBAL, global_attributes_prefix//'include_fleapyear', NF_INT, 1,  include_fleapyear), __LINE__)
  call assert_nf( nf_put_att_int(entry%ncid, NF_GLOBAL, global_attributes_prefix//'use_floatice', NF_INT, 1,  use_floatice), __LINE__)
  call assert_nf( nf_put_att_int(entry%ncid, NF_GLOBAL, global_attributes_prefix//'whichEVP'         , NF_INT, 1,  ice%whichEVP), __LINE__)
  call assert_nf( nf_put_att_int(entry%ncid, NF_GLOBAL, global_attributes_prefix//'evp_rheol_steps'  , NF_INT, 1,  ice%evp_rheol_steps), __LINE__)
  call assert_nf( nf_put_att_int(entry%ncid, NF_GLOBAL, global_attributes_prefix//'opt_visc'         , NF_INT, 1,  dynamics%opt_visc), __LINE__)
  call assert_nf( nf_put_att_int(entry%ncid, NF_GLOBAL, global_attributes_prefix//'use_wsplit'       , NF_INT, 1,  dynamics%use_wsplit), __LINE__)
  call assert_nf( nf_put_att_int(entry%ncid, NF_GLOBAL, global_attributes_prefix//'use_partial_cell', NF_INT, 1,  use_partial_cell), __LINE__)
  call assert_nf( nf_put_att_int(entry%ncid, NF_GLOBAL, global_attributes_prefix//'autorotate_back_to_geo', NF_INT, 1,  vec_autorotate), __LINE__)
 
 
  
!___This ends definition part of the file, below filling in variables is possible
  call assert_nf( nf_enddef(entry%ncid), __LINE__)
  if (entry%dimname(1)=='nz') then
      call assert_nf( nf_put_var_double(entry%ncid, entry%dimvarID(1), abs(mesh%zbar)), __LINE__)
  elseif (entry%dimname(1)=='nz1') then
      call assert_nf( nf_put_var_double(entry%ncid, entry%dimvarID(1), abs(mesh%Z)), __LINE__)
  else
      if (partit%mype==0) write(*,*) 'WARNING: unknown first dimension in 2d mean I/O data'
  end if 
 
  call assert_nf( nf_close(entry%ncid), __LINE__)
end subroutine
!
!--------------------------------------------------------------------------------------------
!
subroutine assoc_ids(entry)
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
  use io_gather_module
  implicit none
  type(Meandata), intent(inout) :: entry
  integer, intent(in) :: entry_index
  integer tag
  integer                       :: i, size1, size2, size_gen, size_lev, order
  integer                       :: c, lev
  integer mpierr


  ! Serial output implemented so far
  if (entry%p_partit%mype==entry%root_rank) then
     write(*,*) 'writing mean record for ', trim(entry%name), '; rec. count = ', entry%rec_count
     call assert_nf( nf_put_vara_double(entry%ncid, entry%Tid, entry%rec_count, 1, entry%ctime_copy, 1), __LINE__)
  end if
! !_______writing 2D and 3D fields________________________________________________
  size1=entry%glsize(1)
  size2=entry%glsize(2)
  tag = 2 ! we can use a fixed tag here as we have an individual communicator for each output field
!___________writing 8 byte real_________________________________________ 
  if (entry%accuracy == i_real8) then
     if(entry%p_partit%mype==entry%root_rank) then
       if(.not. allocated(entry%aux_r8)) allocate(entry%aux_r8(size2))
     end if
     do lev=1, size1
#ifdef ENABLE_ALEPH_CRAYMPICH_WORKAROUNDS
        ! aleph cray-mpich workaround
        call MPI_Barrier(entry%comm, mpierr)
#endif
       if(.not. entry%is_elem_based) then
         call gather_nod2D (entry%local_values_r8_copy(lev,1:size(entry%local_values_r8_copy,dim=2)), entry%aux_r8, entry%root_rank, tag, entry%comm, entry%p_partit)
       else
         call gather_elem2D(entry%local_values_r8_copy(lev,1:size(entry%local_values_r8_copy,dim=2)), entry%aux_r8, entry%root_rank, tag, entry%comm, entry%p_partit)
       end if
        if (entry%p_partit%mype==entry%root_rank) then
          if (entry%ndim==1) then
            call assert_nf( nf_put_vara_double(entry%ncid, entry%varID, (/1, entry%rec_count/), (/size2, 1/), entry%aux_r8, 1), __LINE__)
          elseif (entry%ndim==2) then
            call assert_nf( nf_put_vara_double(entry%ncid, entry%varID, (/1, lev, entry%rec_count/), (/size2, 1, 1/), entry%aux_r8, 1), __LINE__)
          end if
        end if
     end do

!___________writing 4 byte real _________________________________________ 
  else if (entry%accuracy == i_real4) then
     if(entry%p_partit%mype==entry%root_rank) then
       if(.not. allocated(entry%aux_r4)) allocate(entry%aux_r4(size2))
     end if
     do lev=1, size1
#ifdef ENABLE_ALEPH_CRAYMPICH_WORKAROUNDS
        ! aleph cray-mpich workaround
        call MPI_Barrier(entry%comm, mpierr)
#endif
       if(.not. entry%is_elem_based) then
         call gather_real4_nod2D (entry%local_values_r4_copy(lev,1:size(entry%local_values_r4_copy,dim=2)), entry%aux_r4, entry%root_rank, tag, entry%comm, entry%p_partit)
       else
         call gather_real4_elem2D(entry%local_values_r4_copy(lev,1:size(entry%local_values_r4_copy,dim=2)), entry%aux_r4, entry%root_rank, tag, entry%comm, entry%p_partit)
       end if
        if (entry%p_partit%mype==entry%root_rank) then
           if (entry%ndim==1) then
             call assert_nf( nf_put_vara_real(entry%ncid, entry%varID, (/1, entry%rec_count/), (/size2, 1/), entry%aux_r4, 1), __LINE__)
           elseif (entry%ndim==2) then
             call assert_nf( nf_put_vara_real(entry%ncid, entry%varID, (/1, lev, entry%rec_count/), (/size2, 1, 1/), entry%aux_r4, 1), __LINE__)
           end if
        end if
     end do
  end if

end subroutine


subroutine update_means
  implicit none
  type(Meandata), pointer :: entry
  integer                 :: n
  integer                 :: I, J

  DO n=1, io_NSTREAMS
     entry=>io_stream(n)
!_____________ compute in 8 byte accuracy _________________________
     IF (entry%accuracy == i_real8) then
        IF (entry%flip) then
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I,J)
           DO J=1, size(entry%local_values_r8,dim=2)
              DO I=1, size(entry%local_values_r8,dim=1)
                 entry%local_values_r8(I,J)=entry%local_values_r8(I,J)+entry%ptr3(J,I)
              END DO
           END DO
!$OMP END PARALLEL DO
        ELSE
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I,J)
           DO J=1, size(entry%local_values_r8,dim=2)
              DO I=1, size(entry%local_values_r8,dim=1)
                 entry%local_values_r8(I,J)=entry%local_values_r8(I,J)+entry%ptr3(I,J)
              END DO
           END DO
!$OMP END PARALLEL DO
        END IF
!_____________ compute in 4 byte accuracy _________________________
     ELSE IF (entry%accuracy == i_real4) then
        IF (entry%flip) then
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I,J)
           DO J=1, size(entry%local_values_r4,dim=2)
              DO I=1, size(entry%local_values_r4,dim=1)
                 entry%local_values_r4(I,J)=entry%local_values_r4(I,J)+real(entry%ptr3(J,I), real32)
              END DO
           END DO
!$OMP END PARALLEL DO
        ELSE
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I, J)
           DO J=1, size(entry%local_values_r4,dim=2)
              DO I=1, size(entry%local_values_r4,dim=1)
                 entry%local_values_r4(I,J)=entry%local_values_r4(I,J)+real(entry%ptr3(I,J), real32)
              END DO
           END DO
!$OMP END PARALLEL DO
        END IF
     END IF
     entry%addcounter=entry%addcounter+1
  END DO
end subroutine
!
!--------------------------------------------------------------------------------------------
!
subroutine output(istep, ice, dynamics, tracers, partit, mesh)
  use g_clock
  use mod_mesh
  USE MOD_PARTIT
  USE MOD_PARSUP
  use MOD_DYN
  use MOD_ICE
  use mod_tracer
  use io_gather_module
#if defined (__icepack)
  use icedrv_main,    only: init_io_icepack
#endif
  implicit none
  integer       :: istep
  logical, save :: lfirst=.true.
  integer       :: n, k
  integer       :: i, j !for OMP loops
  logical       :: do_output
  type(Meandata), pointer :: entry
  type(t_mesh), intent(in) , target :: mesh
  type(t_partit), intent(inout), target :: partit
  type(t_tracer), intent(in)   , target :: tracers
  type(t_dyn)   , intent(in)   , target :: dynamics
  type(t_ice)   , intent(inout), target :: ice
  
  character(:), allocatable :: filepath
  real(real64)                  :: rtime !timestamp of the record

  ctime=timeold+(dayold-1.)*86400
  if (lfirst) then
     call ini_mean_io(ice, dynamics, tracers, partit, mesh)
     call init_io_gather(partit)
#if defined (__icepack)
     call init_io_icepack(mesh) !icapack has its copy of p_partit => partit
#endif
  end if

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
        call par_ex(partit%MPI_COMM_FESOM, partit%mype, 1)
        stop
     endif

     if (do_output) then
        if (vec_autorotate) call io_r2g(n, partit, mesh) ! automatically detect if a vector field and rotate if makes sense!
        if(entry%thread_running) call entry%thread%join()
        entry%thread_running = .false.

        filepath = trim(ResultPath)//trim(entry%name)//'.'//trim(runid)//'.'//cyearnew//'.nc'
        if(partit%mype == entry%root_rank) then
          if(filepath /= trim(entry%filename)) then
            if("" /= trim(entry%filename)) call assert_nf(nf_close(entry%ncid), __LINE__)   
            entry%filename = filepath
            ! use any existing file with this name or create a new one
            if( nf_open(entry%filename, nf_write, entry%ncid) /= nf_noerr ) then
              call create_new_file(entry, ice, dynamics, partit, mesh)
              call assert_nf( nf_open(entry%filename, nf_write, entry%ncid), __LINE__)
            end if
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

        if (entry%accuracy == i_real8) then
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I,J)
           DO J=1, size(entry%local_values_r8,dim=2)
              DO I=1, size(entry%local_values_r8,dim=1)
                 entry%local_values_r8_copy(I,J) = entry%local_values_r8(I,J) /real(entry%addcounter,real64)  ! compute_means
                 entry%local_values_r8(I,J) = 0._real64 ! clean_meanarrays
              END DO
           END DO
!$OMP END PARALLEL DO
        else if (entry%accuracy == i_real4) then
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I,J)
           DO J=1, size(entry%local_values_r4,dim=2)
              DO I=1, size(entry%local_values_r4,dim=1)
                 entry%local_values_r4_copy(I,J) = entry%local_values_r4(I,J) /real(entry%addcounter,real32)  ! compute_means
                 entry%local_values_r4(I,J) = 0._real32 ! clean_meanarrays
              END DO
           END DO
!$OMP END PARALLEL DO
        end if
        entry%addcounter   = 0  ! clean_meanarrays
        entry%ctime_copy = ctime
        call entry%thread%run()
        entry%thread_running = .true.
     endif
  end do
  lfirst=.false.
end subroutine


subroutine do_output_callback(entry_index)
use mod_mesh
USE MOD_PARTIT
USE MOD_PARSUP
  integer, intent(in) :: entry_index
  ! EO args
  type(Meandata), pointer :: entry


  entry=>io_stream(entry_index)
  entry%p_partit%mype=entry%mype_workaround ! for the thread callback, copy back the value of our mype as a workaround for errors with the cray envinronment (at least with ftn 2.5.9 and cray-mpich 7.5.3)

  call write_mean(entry, entry_index)
  if(entry%p_partit%mype == entry%root_rank) call assert_nf( nf_sync(entry%ncid), __LINE__ ) ! flush the file to disk after each write
end subroutine


subroutine finalize_output()
  integer i
  type(Meandata), pointer :: entry

  do i=1, io_NSTREAMS
    entry=>io_stream(i)
    if(entry%thread_running) call entry%thread%join()
    entry%thread_running = .false.    
  end do
end subroutine
!
!--------------------------------------------------------------------------------------------
!
subroutine def_stream3D(glsize, lcsize, name, description, units, data, freq, freq_unit, accuracy, partit, mesh, flip_array)
  use mod_mesh
  USE MOD_PARTIT
  USE MOD_PARSUP
  implicit none
  type(t_partit),        intent(inout), target :: partit
  integer,               intent(in)    :: glsize(2), lcsize(2)
  character(len=*),      intent(in)    :: name, description, units
  real(kind=WP), target, intent(in)    :: data(:,:)
  integer,               intent(in)    :: freq
  character,             intent(in)    :: freq_unit
  integer,               intent(in)    :: accuracy
  type(Meandata),        allocatable   :: tmparr(:)
  type(Meandata),        pointer       :: entry
  type(t_mesh), intent(in), target     :: mesh
  logical, optional, intent(in)        :: flip_array
  integer i
 
#if !defined(__PGI)  
  do i = 1, rank(data)
    if ((ubound(data, dim = i)<=0)) then
      if (partit%mype==0) then
        write(*,*) 'WARNING: adding I/O stream for ', trim(name), ' failed (contains 0 dimension)'
        write(*,*) 'upper bound is: ', ubound(data, dim = i)
      end if
      return
    end if    
  end do
#endif

  if (partit%mype==0) then
     write(*,*) 'adding I/O stream 3D for ', trim(name)
  end if

  call associate_new_stream(name, entry)

  ! 3d specific
  entry%ptr3 => data                      !2D! entry%ptr3(1:1,1:size(data)) => data

  if (present(flip_array)) then
      if (flip_array) then
          entry%flip = .true.
      else
          entry%flip = .false.
      end if
  else
      entry%flip = .false.
  end if

  entry%ndim=2
  entry%glsize=glsize                     !2D! entry%glsize=(/1, glsize/)

  if (accuracy == i_real8) then
    allocate(entry%local_values_r8(lcsize(1), lcsize(2)))
    entry%local_values_r8 = 0._real64
  elseif (accuracy == i_real4) then
    allocate(entry%local_values_r4(lcsize(1), lcsize(2)))
    entry%local_values_r4 = 0._real32
  end if

  entry%dimname(1)=mesh_dimname_from_dimsize(glsize(1), partit, mesh)     !2D! mesh_dimname_from_dimsize(glsize, mesh)
  entry%dimname(2)=mesh_dimname_from_dimsize(glsize(2), partit, mesh)     !2D! entry%dimname(2)='unknown'
  ! non dimension specific
  call def_stream_after_dimension_specific(entry, name, description, units, freq, freq_unit, accuracy, partit, mesh)
end subroutine
!
!--------------------------------------------------------------------------------------------
!
subroutine def_stream2D(glsize, lcsize, name, description, units, data, freq, freq_unit, accuracy, partit, mesh)
  use mod_mesh
  USE MOD_PARTIT
  USE MOD_PARSUP
  implicit none
  integer,               intent(in)    :: glsize, lcsize
  character(len=*),      intent(in)    :: name, description, units
  real(kind=WP), target, intent(in)    :: data(:)
  integer,               intent(in)    :: freq
  character,             intent(in)    :: freq_unit
  integer,               intent(in)    :: accuracy
  type(Meandata),        allocatable   :: tmparr(:)
  type(Meandata),        pointer       :: entry
  type(t_mesh),          intent(in)    :: mesh
  type(t_partit),        intent(inout) :: partit
  integer i

#if !defined(__PGI)   
  do i = 1, rank(data)
    if ((ubound(data, dim = i)<=0)) then
      if (partit%mype==0) then
        write(*,*) 'WARNING: adding I/O stream for ', trim(name), ' failed (contains 0 dimension)'
        write(*,*) 'upper bound is: ', ubound(data, dim = i)
      end if
      return
    end if    
  end do
#endif

  if (partit%mype==0) then
     write(*,*) 'adding I/O stream 2D for ', trim(name)
  end if

  call associate_new_stream(name, entry)
  
  ! 2d specific
  entry%ptr3(1:1,1:size(data)) => data(:)

  if (accuracy == i_real8) then
    allocate(entry%local_values_r8(1, lcsize))
    entry%local_values_r8 = 0._real64
  elseif (accuracy == i_real4) then
    allocate(entry%local_values_r4(1, lcsize))
    entry%local_values_r4 = 0._real32
  end if

  entry%ndim=1
  entry%glsize=(/1, glsize/)

  entry%dimname(1)=mesh_dimname_from_dimsize(glsize, partit, mesh)
  entry%dimname(2)='unknown'

  ! non dimension specific
  call def_stream_after_dimension_specific(entry, name, description, units, freq, freq_unit, accuracy, partit, mesh)
end subroutine


  subroutine associate_new_stream(name, entry)
    type(Meandata), pointer :: entry
    character(len=*), intent(in) :: name
    integer i

    entry => null()
    
    ! check if we already have this variable
    do i=1, io_NSTREAMS
      if(trim(io_stream(i)%name) .eq. name) then
          print *,"variable '"//name//"' already exists, &
              &check if you define it multiple times, for example in namelist.io, &
              &namelist.icepack, io_meandata.F90 or other place that add I/O stream."
          call assert(.false., __LINE__) 
      end if
    end do
        
    ! add this instance to io_stream array
    io_NSTREAMS = io_NSTREAMS +1
    call assert(size(io_stream) >= io_NSTREAMS, __LINE__)
    entry=>io_stream(io_NSTREAMS)
  end subroutine


  subroutine def_stream_after_dimension_specific(entry, name, description, units, freq, freq_unit, accuracy, partit, mesh)
    use mod_mesh
    USE MOD_PARTIT
    USE MOD_PARSUP
    use io_netcdf_workaround_module
    type(Meandata), intent(inout) :: entry
    character(len=*),      intent(in)    :: name, description, units
    integer,               intent(in)    :: freq
    character,             intent(in)    :: freq_unit
    integer,               intent(in)    :: accuracy
    type(t_mesh), intent(in), target     :: mesh
    type(t_partit), intent(inout), target :: partit
    ! EO args
    logical async_netcdf_allowed
    integer provided_mpi_thread_support_level
    integer entry_index
    integer err
    
    entry_index = io_NSTREAMS
    
    entry%accuracy = accuracy

    if (accuracy == i_real8) then
      allocate(data_strategy_nf_double_type :: entry%data_strategy)
    elseif (accuracy == i_real4) then
      allocate(data_strategy_nf_float_type :: entry%data_strategy)
    else
       if (partit%mype==0) write(*,*) 'not supported output accuracy:',accuracy,'for',trim(name)
       call par_ex(partit%MPI_COMM_FESOM, partit%mype)
       stop
    endif ! accuracy

    entry%name = name
    entry%description = description
    entry%units = units
    entry%filename = ""

    entry%freq=freq
    entry%freq_unit=freq_unit
    entry%addcounter   = 0
    entry%is_in_use=.true.

    if(entry%glsize(1)==mesh%nod2D  .or. entry%glsize(2)==mesh%nod2D) then
      entry%is_elem_based = .false.
    else if(entry%glsize(1)==mesh%elem2D .or. entry%glsize(2)==mesh%elem2D) then
      entry%is_elem_based = .true.
    else
      if(partit%mype == 0) print *,"can not determine if ",trim(name)," is node or elem based"
      stop
    end if

    if (accuracy == i_real8) then
      allocate(entry%local_values_r8_copy(size(entry%local_values_r8, dim=1), size(entry%local_values_r8, dim=2)))
    else if (accuracy == i_real4) then
      allocate(entry%local_values_r4_copy(size(entry%local_values_r4, dim=1), size(entry%local_values_r4, dim=2)))
    end if

    ! set up async output
    
    entry%root_rank = next_io_rank(partit%MPI_COMM_FESOM, async_netcdf_allowed, partit)

    call MPI_Comm_dup(partit%MPI_COMM_FESOM, entry%comm, err)

    call entry%thread%initialize(do_output_callback, entry_index)
    if(.not. async_netcdf_allowed) call entry%thread%disable_async()
  
    ! check if we have multi thread support available in the MPI library
    ! tough MPI_THREAD_FUNNELED should be enough here, at least on cray-mpich 7.5.3 async mpi calls fail if we do not have support level 'MPI_THREAD_MULTIPLE'
    ! on cray-mpich we only get level 'MPI_THREAD_MULTIPLE' if 'MPICH_MAX_THREAD_SAFETY=multiple' is set in the environment
    call MPI_Query_thread(provided_mpi_thread_support_level, err)
    if(provided_mpi_thread_support_level < MPI_THREAD_MULTIPLE) call entry%thread%disable_async()
    
    entry%mype_workaround = partit%mype ! make a copy of the mype variable as there is an error with the cray compiler or environment which voids the global mype for our threads
    entry%p_partit=>partit
  end subroutine


  subroutine assert_nf(status, line)
    integer, intent(in) :: status
    integer, intent(in) :: line
    ! EO args
    include "netcdf.inc" ! old netcdf fortran interface required?
    if(status /= NF_NOERR) then
      print *, "error in line ",line, __FILE__, ' ', trim(nf_strerror(status))
      stop 1
    end if   
  end subroutine


  subroutine assert(val, line)
    logical, intent(in) :: val
    integer, intent(in) :: line
    ! EO args
    if(.NOT. val) then
      print *, "error in line ",line, __FILE__
      stop 1
    end if
  end subroutine


  subroutine io_r2g(n, partit, mesh)
    USE MOD_MESH
    USE MOD_PARTIT
    USE g_rotate_grid
    implicit none
    integer,        intent(in)            :: n
    type(t_partit), intent(inout), target :: partit
    type(t_mesh),   intent(in),    target :: mesh
    integer                               :: I, J
    type(Meandata), pointer               :: entry_x, entry_y
    real(kind=WP)                         :: temp_x, temp_y
    real(kind=WP)                         :: xmean, ymean
    logical                               :: do_rotation

    if (n==io_NSTREAMS) RETURN
    entry_x=>io_stream(n)
    entry_y=>io_stream(n+1)
    IF (.NOT. (entry_x%freq_unit==entry_y%freq_unit) .and. ((entry_x%freq==entry_y%freq))) RETURN
    IF (entry_x%accuracy /= entry_y%accuracy) RETURN
    do_rotation=.FALSE.
! we need to improve the logistic here in order to use this routinely. a new argument in def_stream
! will be needed.
    IF ((trim(entry_x%name)=='u'       ) .AND. ((trim(entry_y%name)=='v'       ))) do_rotation=.TRUE.
    IF ((trim(entry_x%name)=='uice'    ) .AND. ((trim(entry_y%name)=='vice'    ))) do_rotation=.TRUE.
    IF ((trim(entry_x%name)=='unod'    ) .AND. ((trim(entry_y%name)=='vnod'    ))) do_rotation=.TRUE.
    IF ((trim(entry_x%name)=='tau_x'   ) .AND. ((trim(entry_y%name)=='tau_y   '))) do_rotation=.TRUE.
    IF ((trim(entry_x%name)=='atmice_x') .AND. ((trim(entry_y%name)=='atmice_y'))) do_rotation=.TRUE.
    IF ((trim(entry_x%name)=='atmoce_x') .AND. ((trim(entry_y%name)=='atmoce_y'))) do_rotation=.TRUE.    
    IF ((trim(entry_x%name)=='iceoce_x') .AND. ((trim(entry_y%name)=='iceoce_y'))) do_rotation=.TRUE.    

    IF (.NOT. (do_rotation)) RETURN
   
    IF (partit%mype==0) THEN
       write(*,*) trim(entry_x%name)//' and '//trim(entry_y%name)//' will be rotated before output!'
    END IF

    IF ((entry_x%accuracy == i_real8) .AND. (entry_y%accuracy == i_real8)) THEN
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I, J, xmean, ymean)
    DO J=1, size(entry_x%local_values_r8,dim=2)
       if (entry_x%is_elem_based) then
          xmean=sum(mesh%coord_nod2D(1, mesh%elem2D_nodes(:, J)))/3._WP
          ymean=sum(mesh%coord_nod2D(2, mesh%elem2D_nodes(:, J)))/3._WP
       else
          xmean=mesh%coord_nod2D(1, J)
          ymean=mesh%coord_nod2D(2, J)
       end if
       DO I=1, size(entry_x%local_values_r8,dim=1)
          call vector_r2g(entry_x%local_values_r8(I,J), entry_y%local_values_r8(I,J), xmean, ymean, 0)
       END DO
    END DO
!$OMP END PARALLEL DO
    END IF

    IF ((entry_x%accuracy == i_real4) .AND. (entry_y%accuracy == i_real4)) THEN
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I, J, temp_x, temp_y, xmean, ymean)
    DO J=1, size(entry_x%local_values_r4,dim=2)
       if (entry_x%is_elem_based) then
          xmean=sum(mesh%coord_nod2D(1, mesh%elem2D_nodes(:, J)))/3._WP
          ymean=sum(mesh%coord_nod2D(2, mesh%elem2D_nodes(:, J)))/3._WP
       else
          xmean=mesh%coord_nod2D(1, J)
          ymean=mesh%coord_nod2D(2, J)
       end if
       DO I=1, size(entry_x%local_values_r4,dim=1)
          temp_x=real(entry_x%local_values_r4(I,J), real64)
          temp_y=real(entry_y%local_values_r4(I,J), real64)
          call vector_r2g(temp_x, temp_y, xmean, ymean, 0)
          entry_x%local_values_r4(I,J)=real(temp_x, real32)
          entry_y%local_values_r4(I,J)=real(temp_y, real32)
       END DO
    END DO
!$OMP END PARALLEL DO
    END IF
  end subroutine
end module
