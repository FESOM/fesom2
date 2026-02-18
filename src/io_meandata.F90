module io_MEANDATA
  USE MOD_PARTIT
  USE MOD_PARSUP
#if defined(__recom)
  use recom_glovar
  use recom_config
  use recom_ciso
#endif
  USE g_clock
  use o_PARAM, only : WP
  use, intrinsic :: iso_fortran_env, only: real64, real32
  use io_data_strategy_module
  use async_threads_module
  use netcdf

  implicit none
  private
  public :: def_stream, def_stream2D, def_stream3D, output, finalize_output
!
!--------------------------------------------------------------------------------------------
!
  integer, parameter  :: i_real8=8, i_real4=4

  ! NetCDF standard fill values for missing/invalid data
  real(real32), parameter :: NC_FILL_FLOAT  = 9.9692099683868690e+36_real32
  real(real64), parameter :: NC_FILL_DOUBLE = 9.9692099683868690e+36_real64

  type Meandata
    private
    type(t_partit), pointer                            :: p_partit
    integer                                            :: ndim
    integer                                            :: glsize(2)
    integer                                            :: shrinked_size
    integer, allocatable, dimension(:)                 :: shrinked_indx
    integer                                            :: accuracy
    real(real64), allocatable, dimension(:,:) :: local_values_r8
    real(real32), allocatable, dimension(:,:) :: local_values_r4
    real(real64), allocatable :: aux_r8(:)
    real(real32), allocatable :: aux_r4(:)
    integer                                            :: addcounter =0
    integer                                            :: lastcounter=0 ! before addcounter is set to 0
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
    character(500) :: long_description="" 
    character(100) :: defined_on=""
    character(100) :: mesh="fesom_mesh"
    ! to be passed to MULTIO (time window for the accumulations)
    integer :: currentDate,  currentTime
    integer :: previousDate, previousTime
    integer :: startDate, startTime
  contains
    final destructor
  end type Meandata
!
!--------------------------------------------------------------------------------------------
!
  type(Meandata), save, target   :: io_stream(150) ! todo: find a way to increase the array withhout move_alloc to keep the derived types in Meandata intact
  integer, save                  :: io_NSTREAMS=0
  real(kind=WP)                  :: ctime !current time in seconds from the beginning of the year
!
!--------------------------------------------------------------------------------------------
!
  integer, save                  :: io_listsize   =0
  logical, save                  :: vec_autorotate=.FALSE.
  logical, save                  :: lnextGEMS     =.FALSE.
  integer, save                  :: nlev_upper=1
  character(len=1), save         :: filesplit_freq='y'
  integer, save                  :: compression_level=0
  type io_entry
        CHARACTER(len=15)        :: id        ='unknown   '
        INTEGER                  :: freq      =0
        CHARACTER                :: unit      =''
        INTEGER                  :: precision =0
  end type io_entry 

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
  REAL(real64), DIMENSION(:), ALLOCATABLE, TARGET :: multio_temporary_array

  contains
!
!--------------------------------------------------------------------------------------------
!

!
!
!_______________________________________________________________________________
! not sure why this is needed --> seems to become method of meandata stream object
! type Meandata
!   private
!   ...  
!   contains
!       final destructor
! end type  
subroutine destructor(this)
    type(Meandata), intent(inout) :: this
    ! EO args
    call assert_nf(nf90_close(this%ncid), __LINE__)
end subroutine destructor
!
!
!_______________________________________________________________________________
! define 2d/3d meandata stream parameter
subroutine ini_mean_io(ice, dynamics, tracers, partit, mesh)
    !------------------------------------------
    ! LA 2023-01-31 add iceberg params
    use iceberg_params
    !------------------------------------------
    use MOD_MESH
    use MOD_TRACER
    use MOD_PARTIT
    use MOD_PARSUP
    use MOD_DYN
    use MOD_ICE
    use o_ARRAYS       
    use o_mixing_KPP_mod
    use g_backscatter
    use diagnostics
    use g_forcing_arrays
#if defined (__cvmix)    
    use g_cvmix_tke
    use g_cvmix_idemix
    use g_cvmix_kpp
    use g_cvmix_tidal
#endif    
#if defined(__recom)
    use recom_glovar
    use recom_config
    use recom_ciso
#endif
    use g_forcing_param, only: use_virt_salt, use_landice_water, use_age_tracer !---fwf-code, age-code
    use g_config, only : use_cavity, lwiso !---wiso-code
    use mod_transit, only : index_transit_r14c, index_transit_r39ar, index_transit_f11, index_transit_f12, index_transit_sf6

    implicit none
    integer                   :: i, j
    integer, save             :: nm_io_unit  = 103       ! unit to open namelist file, skip 100-102 for cray
    integer                   :: iost
    integer,dimension(15)     :: sel_forcvar=0
    integer                   :: sel_dmoc=0, sel_trgrd_xyz=0, sel_dvd=0, sel_redi=0
    character(len=10)         :: id_string

    type(t_mesh), intent(in) , target :: mesh
    type(t_partit), intent(inout), target :: partit
    type(t_tracer), intent(in)   , target :: tracers
    type(t_dyn)   , intent(in)   , target :: dynamics
    type(t_ice)   , intent(in)   , target :: ice
    namelist /nml_general / io_listsize, vec_autorotate, lnextGEMS, nlev_upper, filesplit_freq, compression_level
    namelist /nml_list    / io_list

#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"

    !___________________________________________________________________________
    ! OPEN and read namelist for I/O
    open( unit=nm_io_unit, file='namelist.io', form='formatted', access='sequential', status='old', iostat=iost )
    if (iost == 0) then
      if (mype==0) WRITE(*,*) '     file   : ', 'namelist.io',' open ok'
    else
      if (mype==0) WRITE(*,*) 'ERROR: --> file not found   : ', 'namelist.io',' ; iostat=',iost
      call par_ex(partit%MPI_COMM_FESOM, partit%mype, 1)
      stop
    endif

    READ(nm_io_unit, nml=nml_general,  iostat=iost )
    if (iost/=0) then
       if (mype==0) WRITE(*,*) 'ERROR: in reading nml_general block in namelist.io, invalid formatting.'
       call par_ex(partit%MPI_COMM_FESOM, partit%mype, 1)
    endif

    allocate(io_list(io_listsize))
    READ(nm_io_unit, nml=nml_list,     iostat=iost )
    close(nm_io_unit )

    !___________________________________________________________________________
    ! TODO: unknown variable found then write clearly in log, saves lot of frustration.
    do i=1, io_listsize
        if (trim(io_list(i)%id)=='unknown   ') then
            if (mype==0) write(*,*) 'io_listsize will be changed from ', io_listsize, ' to ', i-1, '!'
            io_listsize=i-1
            EXIT
        end if
    end do

!_______________________________________________________________________________    
DO i=1, io_listsize
SELECT CASE (trim(io_list(i)%id))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!2D streams!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CASE ('sst       ')
    call def_stream(nod2D, myDim_nod2D, 'sst',      'sea surface temperature',        'C', tracers%data(1)%values(1,1:myDim_nod2D), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh, "Sea surface temperature")
CASE ('sss       ')
    call def_stream(nod2D, myDim_nod2D, 'sss',      'sea surface salinity',           'psu', tracers%data(2)%values(1,1:myDim_nod2D), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh, "Sea surface salinity")
CASE ('ssh       ')
    call def_stream(nod2D, myDim_nod2D, 'ssh',      'sea surface elevation',          'm',      dynamics%eta_n,                          io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
CASE ('vve_5     ')
    call def_stream(nod2D, myDim_nod2D, 'vve_5',    'vertical velocity at 5th level', 'm/s',    dynamics%w(5,:),                         io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
CASE ('t_star        ')
    call def_stream(nod2D, myDim_nod2D,'t_star'        , 'air temperature'      , 'C'  , t_star(:)       , io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
CASE ('qsr        ')
    call def_stream(nod2D, myDim_nod2D,'qsr'        , 'solar radiation'      , 'W/s^2'  , qsr_c(:)       , io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)


! 2d ssh diagnostic variables
CASE ('ssh_rhs    ')
    call def_stream(nod2D, myDim_nod2D,  'ssh_rhs'    , 'ssh rhs'           , 'm/s', dynamics%ssh_rhs    , io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
CASE ('ssh_rhs_old')
    call def_stream(nod2D, myDim_nod2D,  'ssh_rhs_old', 'ssh rhs'           , 'm/s', dynamics%ssh_rhs_old, io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
CASE ('d_eta      ')
    call def_stream(nod2D, myDim_nod2D,  'd_eta'      , 'dssh (from solver)', 'm'  , dynamics%d_eta      , io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
CASE ('hbar       ')
    call def_stream(nod2D, myDim_nod2D,  'hbar'       , 'ssh n+0.5 tstep'   , 'm'  , mesh%hbar           , io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
CASE ('hbar_old   ')
    call def_stream(nod2D, myDim_nod2D,  'hbar_old'   , 'ssh n-0.5 tstep'   , 'm'  , mesh%hbar_old       , io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
CASE ('dhe        ')
    call def_stream(elem2D, myDim_elem2D,'dhe'        , 'dhbar @ elem'      , 'm'  , mesh%dhe            , io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    
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
    call def_stream(nod2D, myDim_nod2D, 'a_ice',    'ice concentration',              'fraction',      ice%data(1)%values(1:myDim_nod2D),      io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    end if
CASE ('m_ice     ')
    if (use_ice) then
    call def_stream(nod2D, myDim_nod2D, 'm_ice',    'ice height per unit area',       'm',      ice%data(2)%values(1:myDim_nod2D),      io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    end if

! ice thermodynamic growth rate: ice, snow, area    
CASE ('thdgrice  ')
    if (use_ice) then
    call def_stream(nod2D, myDim_nod2D, 'thdgrice' , 'thermodynamic growth rate ice',               'm/s', ice%thermo%thdgr(1:myDim_nod2D),      io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    end if
CASE ('thdgrsnw  ')
    if (use_ice) then
    call def_stream(nod2D, myDim_nod2D, 'thdgrsnw' , 'thermodynamic growth rate snow',              'm/s', ice%thermo%thdgrsn(1:myDim_nod2D),    io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    end if
CASE ('thdgrarea ')
    if (use_ice) then
    call def_stream(nod2D, myDim_nod2D, 'thdgrarea', 'thermodynamic growth rate ice concentration', 'frac/s', ice%thermo%thdgra(1:myDim_nod2D),    io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    end if
    
! ice dynamic growth rate: ice, snow, area    
CASE ('dyngrice  ')
    if (use_ice) then
    call def_stream(nod2D, myDim_nod2D, 'dyngrice' , 'dynamic growth rate ice',               'm/s', ice%thermo%dyngr(1:myDim_nod2D),      io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    end if
CASE ('dyngrsnw  ')
    if (use_ice) then
    call def_stream(nod2D, myDim_nod2D, 'dyngrsnw' , 'dynamic growth rate snow',              'm/s', ice%thermo%dyngrsn(1:myDim_nod2D),    io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    end if
CASE ('dyngrarea ')
    if (use_ice) then
    call def_stream(nod2D, myDim_nod2D, 'dyngrarea', 'dynamic growth rate ice concentration', 'frac/s', ice%thermo%dyngra(1:myDim_nod2D),    io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    end if

CASE ('flice     ')
    if (use_ice) then
    call def_stream(nod2D, myDim_nod2D,  'flice',    'flooding growth rate ice',       'm/s',    flice(:),                  io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    end if
CASE ('m_snow    ')
    if (use_ice) then
    call def_stream(nod2D, myDim_nod2D, 'm_snow',   'snow height per unit area',      'm',      ice%data(3)%values(1:myDim_nod2D),     io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    end if
CASE ('h_ice ')
    if (use_ice) then
    call def_stream(nod2D, myDim_nod2D, 'h_ice',    'ice thickness over ice-covered fraction',   'm',     ice%h_ice(1:myDim_nod2D), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    end if
CASE ('h_snow    ')
    if (use_ice) then
    call def_stream(nod2D, myDim_nod2D, 'h_snow',   'snow thickness over ice-covered fraction',  'm',     ice%h_snow(1:myDim_nod2D), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    end if
! Melt pond variables
CASE ('apnd      ')
    if (use_ice .and. ice%thermo%use_meltponds) then
    call def_stream(nod2D, myDim_nod2D, 'apnd',     'melt pond area fraction',                  'frac',  ice%thermo%apnd(1:myDim_nod2D), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    end if
CASE ('hpnd      ')
    if (use_ice .and. ice%thermo%use_meltponds) then
    call def_stream(nod2D, myDim_nod2D, 'hpnd',     'melt pond depth',                          'm',     ice%thermo%hpnd(1:myDim_nod2D), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    end if
CASE ('ipnd      ')
    if (use_ice .and. ice%thermo%use_meltponds) then
    call def_stream(nod2D, myDim_nod2D, 'ipnd',     'melt pond ice lid thickness',              'm',     ice%thermo%ipnd(1:myDim_nod2D), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    end if    
! Debug ice variables    
CASE ('strength_ice')
    if (use_ice) then
    call def_stream(elem2D, myDim_elem2D, 'strength_ice', 'ice strength', '?', ice%work%ice_strength(1:myDim_elem2D), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    end if
CASE ('inv_areamass')
    if (use_ice) then
    call def_stream(nod2D,  myDim_nod2D,  'inv_areamass', 'inv_areamass', '?', ice%work%inv_areamass(1:myDim_nod2D) , io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    end if
CASE ('rhs_a')
    if (use_ice) then
    call def_stream(nod2D,  myDim_nod2D,  'rhs_a'    , 'rhs_a'    , '?', ice%data(1)%values_rhs(1:myDim_nod2D), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    end if
CASE ('rhs_m')
    if (use_ice) then
    call def_stream(nod2D,  myDim_nod2D,  'rhs_m'    , 'rhs_m'    , '?', ice%data(2)%values_rhs(1:myDim_nod2D), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    end if
CASE ('sgm11')
    if (use_ice) then
    call def_stream(elem2D, myDim_elem2D, 'sgm11'    , 'sgm11'    , '?', ice%work%sigma11(1:myDim_elem2D)     , io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    end if
CASE ('sgm12')
    if (use_ice) then
    call def_stream(elem2D, myDim_elem2D, 'sgm12'    , 'sgm12'    , '?', ice%work%sigma12(1:myDim_elem2D)     , io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    end if
CASE ('sgm22')
    if (use_ice) then
    call def_stream(elem2D, myDim_elem2D, 'sgm22'    , 'sgm22'    , '?', ice%work%sigma22(1:myDim_elem2D)     , io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    end if
CASE ('eps11')
    if (use_ice) then
    call def_stream(elem2D, myDim_elem2D, 'eps11'    , 'eps11'    , '?', ice%work%eps11(1:myDim_elem2D)       , io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    end if
CASE ('eps12')
    if (use_ice) then
    call def_stream(elem2D, myDim_elem2D, 'eps12'    , 'eps12'    , '?', ice%work%eps12(1:myDim_elem2D)       , io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    end if
CASE ('eps22')
    if (use_ice) then
    call def_stream(elem2D, myDim_elem2D, 'eps22'    , 'eps22'    , '?', ice%work%eps22(1:myDim_elem2D)       , io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    end if
CASE ('u_rhs_ice')
    if (use_ice) then
    call def_stream(nod2D,  myDim_nod2D , 'u_rhs_ice', 'u_rhs_ice', '?', ice%uice_rhs(1:myDim_nod2D)          , io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    end if
CASE ('v_rhs_ice')
    if (use_ice) then
    call def_stream(nod2D,  myDim_nod2D , 'v_rhs_ice', 'v_rhs_ice', '?', ice%vice_rhs(1:myDim_nod2D)          , io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    end if
CASE ('metric_fac')
    if (use_ice) then
    call def_stream(elem2D, myDim_elem2D, 'metric_fac', 'metric_fac',               '?',    mesh%metric_factor(1:myDim_elem2D),     io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    end if
CASE ('elevat_ice')
    if (use_ice) then
    call def_stream(nod2D,  myDim_nod2D , 'elevat_ice', 'elevat_ice',               '?',      ice%srfoce_ssh(1:myDim_nod2D),     io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    end if
CASE ('uwice')
    if (use_ice) then
    call def_stream(nod2D,  myDim_nod2D , 'uwice'     , 'uwice'     ,               '?',      ice%srfoce_u(1:myDim_nod2D),     io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    end if
CASE ('vwice')
    if (use_ice) then
    call def_stream(nod2D,  myDim_nod2D ,  'vwice',   'vwice',               '?',      ice%srfoce_v(1:myDim_nod2D),     io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    end if
CASE ('twice')
    if (use_ice) then
    call def_stream(nod2D,  myDim_nod2D ,  'twice',   'twice',               '?',      ice%srfoce_temp(1:myDim_nod2D),     io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    end if
CASE ('swice')
    if (use_ice) then
    call def_stream(nod2D,  myDim_nod2D ,  'swice',   'swice',               '?',      ice%srfoce_salt(1:myDim_nod2D),     io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    end if
    
!_______________________________________________________________________________
! output mixed layer depth    
CASE ('MLD1      ')
    call def_stream(nod2D, myDim_nod2D, 'MLD1',     'Mixed Layer Depth',               'm',      MLD1(1:myDim_nod2D),       io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
CASE ('MLD2      ')
    call def_stream(nod2D, myDim_nod2D, 'MLD2',     'Mixed Layer Depth',               'm',      MLD2(1:myDim_nod2D),       io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
CASE ('MLD3      ')
    call def_stream(nod2D, myDim_nod2D, 'MLD3',     'Mixed Layer Depth',               'm',      MLD3(1:myDim_nod2D),       io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)

!_______________________________________________________________________________
! output heat content (for destine)
CASE ('hc300m')
    if (ldiag_destine) then
        call def_stream(nod2D, myDim_nod2D, 'hc300m', 'Vertically integrated heat content upper 300m',   'J m**-2', heatcontent(1:myDim_nod2D,1), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    end if 
CASE ('hc700m')
    if (ldiag_destine) then
        call def_stream(nod2D, myDim_nod2D, 'hc700m', 'Vertically integrated heat content upper 700m',   'J m**-2', heatcontent(1:myDim_nod2D,2), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    end if 
CASE ('hc')
    if (ldiag_destine) then
        call def_stream(nod2D, myDim_nod2D, 'hc',     'Vertically integrated heat content total column', 'J m**-2', heatcontent(1:myDim_nod2D,3), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    end if 
!_______________________________________________________________________________
!---wiso-code
! output water isotopes in sea ice
CASE ('h2o18_ice ')
    if (lwiso) then
    call def_stream(nod2D,  myDim_nod2D,  'h2o18_ice',      'h2o18 concentration in sea ice',    'kmol/m**3',    tr_arr_ice(:,1),    io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    endif
CASE ('hDo16_ice ')
    if (lwiso) then
    call def_stream(nod2D,  myDim_nod2D,  'hDo16_ice',      'hDo16 concentration in sea ice',    'kmol/m**3',    tr_arr_ice(:,2),    io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    endif
CASE ('h2o16_ice ')
    if (lwiso) then
    call def_stream(nod2D,  myDim_nod2D,  'h2o16_ice',      'h2o16 concentration in sea ice',    'kmol/m**3',    tr_arr_ice(:,3),    io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    endif
!---wiso-code-end

!---fwf-code-begin
CASE ('landice   ')
    if (use_landice_water) then
    call def_stream(nod2D,  myDim_nod2D,  'landice',      'freshwater flux',    'm/s',    runoff_landice,    io_list(i)%freq,    io_list(i)%unit, io_list(i)%precision, partit, mesh)
    endif
!---fwf-code-end

!---age-code-begin
CASE ('age       ')
    if (use_age_tracer) then
    call def_stream((/nl-1, nod2D/),  (/nl-1, myDim_nod2D/),  'age',      'water age tracer',    'year',    tracers%data(index_age_tracer)%values(:,:),             io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    end if
!---age-code-end

!_______________________________________________________________________________
! output surface forcing
CASE ('fh        ')
    call def_stream(nod2D, myDim_nod2D, 'fh'       , 'heat flux',                             'W/m2', heat_flux_in(:),           io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
CASE ('fh_sen    ')
    call def_stream(nod2D, myDim_nod2D, 'fh_sen'   , 'sensible heat flux',                    'W/m2', heat_flux_in(:),           io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
CASE ('fh_lat    ')
    call def_stream(nod2D, myDim_nod2D, 'fh_lat'   , 'latent heat flux',                      'W/m2', heat_flux_in(:),           io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
CASE ('fh_radtot ')
    call def_stream(nod2D, myDim_nod2D, 'fh_radtot', 'total radiation heat flux',             'W/m2', heat_flux_in(:),           io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
CASE ('fh_swr    ')
    call def_stream(nod2D, myDim_nod2D, 'fh_swr'   , 'shortwave radiation heat flux',         'W/m2', heat_flux_in(:),           io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
CASE ('fh_lwr    ')
    call def_stream(nod2D, myDim_nod2D, 'fh_lwr'   , 'longwave radiation heat flux',          'W/m2', heat_flux_in(:),           io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
CASE ('fh_lwrout ')
    call def_stream(nod2D, myDim_nod2D, 'fh_lwrout', 'outgoing longwave radiation heat flux', 'W/m2', heat_flux_in(:),           io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)

CASE ('fw        ')
    call def_stream(nod2D, myDim_nod2D, 'fw'       , 'fresh water flux',                'm/s',    water_flux(:),             io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
CASE ('fw_ice    ')
    call def_stream(nod2D, myDim_nod2D, 'fw_ice'   , 'fresh water flux from ice',       'm/s',    fw_ice(:),             io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
CASE ('fw_snw    ')
    call def_stream(nod2D, myDim_nod2D, 'fw_snw'   , 'fresh water flux from snow',      'm/s',    fw_snw(:),             io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)

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
    call def_stream(nod2D, myDim_nod2D, 'dflux',    'density flux',               'kg/(m3*s)',   dens_flux(:),              io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
CASE ('runoff    ')
    sel_forcvar(10)= 1
    call def_stream(nod2D, myDim_nod2D, 'runoff',   'river runoff',                    'm/s',    runoff(:),                 io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
CASE ('evap      ')
    sel_forcvar(7) = 1
    call def_stream(nod2D, myDim_nod2D, 'evap',     'evaporation',                     'm/s',    evaporation(:),            io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
CASE ('prec      ')
    sel_forcvar(5) = 1
    call def_stream(nod2D, myDim_nod2D, 'prec',     'precipitation rain',              'm/s',    prec_rain(:),              io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
CASE ('snow      ')
    sel_forcvar(6) = 1
    call def_stream(nod2D, myDim_nod2D, 'snow',     'precipitation snow',              'm/s',    prec_snow(:),              io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
CASE ('tair      ')
    sel_forcvar(3) = 1
    call def_stream(nod2D, myDim_nod2D, 'tair',     'surface air temperature',         '°C',     Tair(:),                   io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
CASE ('shum      ')
    sel_forcvar(4) = 1
    call def_stream(nod2D, myDim_nod2D, 'shum',     'specific humidity',               '',       shum(:),                   io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
CASE ('swr       ')
    sel_forcvar(8) = 1
    call def_stream(nod2D, myDim_nod2D, 'swr',      'short wave radiation',            'W/m^2',  shortwave(:),              io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
CASE ('lwr       ')
    sel_forcvar(9) = 1
    call def_stream(nod2D, myDim_nod2D, 'lwr',      'long wave radiation',             'W/m^2',  longwave(:),               io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
CASE ('uwind     ')
    sel_forcvar(1) = 1
    call def_stream(nod2D, myDim_nod2D, 'uwind',    '10m zonal surface wind velocity', 'm/s',    u_wind(:),               io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
CASE ('vwind     ')
    sel_forcvar(2) = 1
    call def_stream(nod2D, myDim_nod2D, 'vwind',    '10m merid. surface wind velocity','m/s',    v_wind(:),               io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)

!___________________________________________________________________________________________________________________________________    
CASE ('virtsalt  ')
    sel_forcvar(13) = 1
    call def_stream(nod2D , myDim_nod2D , 'virtsalt' , 'virtual salt flux'          , 'm/s*psu', virtual_salt(:), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
CASE ('relaxsalt ')
    sel_forcvar(14) = 1
    call def_stream(nod2D , myDim_nod2D , 'relaxsalt', 'relaxation salt flux'       , 'm/s*psu', relax_salt(:) , io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
CASE ('realsalt  ')
    sel_forcvar(15) = 1
    call def_stream(nod2D , myDim_nod2D , 'realsalt' , 'real salt flux from sea ice', 'm/s*psu', real_salt_flux(:)  , io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)

CASE ('ice_rejectsalt')
    if (SPP) call def_stream(nod2D , myDim_nod2D , 'ice_rejectsalt' , 'salt flux from plum parameterisation ', 'm/s*psu', ice_rejected_salt(:), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    
!___________________________________________________________________________________________________________________________________
! output KPP vertical mixing schemes
CASE ('kpp_obldepth   ')
    if     (mix_scheme_nmb==1 .or. mix_scheme_nmb==17) then! fesom KPP
        call def_stream(nod2D, myDim_nod2D,    'kpp_obldepth',    'KPP ocean boundary layer depth', 'm',   hbl(:),          io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
#if defined (__cvmix)    
    elseif (mix_scheme_nmb==3 .or. mix_scheme_nmb==37) then ! cvmix KPP
        call def_stream(nod2D, myDim_nod2D,    'kpp_obldepth',    'KPP ocean boundary layer depth', 'm',   kpp_obldepth(:), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
#endif        
    end if
CASE ('kpp_sbuoyflx')
    if     (mix_scheme_nmb==1 .or. mix_scheme_nmb==17) then ! fesom KPP
        call def_stream(nod2D, myDim_nod2D,    'kpp_sbuoyflx',    'surface buoyancy flux',   'm2/s3',  Bo(:),             io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
#if defined (__cvmix)            
    elseif (mix_scheme_nmb==3 .or. mix_scheme_nmb==37) then ! cvmix KPP
        call def_stream(nod2D, myDim_nod2D,    'kpp_sbuoyflx',    'surface buoyancy flux',   'm2/s3',  kpp_sbuoyflx(:), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
#endif        
    end if
CASE ('tx_sur    ')
    sel_forcvar(11) = 1
    call def_stream(elem2D, myDim_elem2D,  'tx_sur',    'zonal wind str. to ocean',       'N/m2',   stress_surf(1, :),         io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
CASE ('ty_sur    ')
    sel_forcvar(12) = 1
    call def_stream(elem2D, myDim_elem2D,  'ty_sur',    'meridional wind str. to ocean',  'N/m2',   stress_surf(2, :),         io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
CASE ('curl_surf ')
    if (lcurt_stress_surf) then
    call def_stream(nod2D, myDim_nod2D,    'curl_surf', 'vorticity of the surface stress', 'none',   curl_stress_surf(:),       io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    
    end if
!___________________________________________________________________________________________________________________________________
! output RECOM 2D
#if defined(__recom)

!===============================================================================
! RECOM BIOGEOCHEMICAL MODEL - OUTPUT VARIABLE DEFINITIONS
!===============================================================================
! This module defines output streams for the REcoM (Regulated Ecosystem Model)
! biogeochemical tracer variables. Variables are organized by functional groups.
!===============================================================================

! =====================================================================
! CARBONATE CHEMISTRY & AIR-SEA GAS EXCHANGE
! =====================================================================
! These variables control and diagnose the exchange of CO2 and O2 between
! the ocean and atmosphere, as well as ocean acidification processes.
! =====================================================================

CASE ('alphaCO2  ')
    ! =====================================================================
    ! Variable: alphaCO2
    ! Description: CO2 solubility coefficient (Henry's Law constant)
    ! Function: Determines how much CO2 dissolves in seawater at equilibrium
    ! Dependencies: Temperature and salinity dependent
    ! Units: mol/kg/atm
    ! =====================================================================
    if (use_REcoM) then
        call def_stream(nod2D, myDim_nod2D, &
                       'alphaCO2', 'CO2 solubility', 'mol/kg/atm', alphaCO2(:), &
                       io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    end if

CASE ('Kw        ')
    ! =====================================================================
    ! Variable: Kw (also known as k660 or gas transfer velocity)
    ! Description: Air-sea piston velocity for gas exchange
    ! Function: Controls the rate of gas transfer across the air-sea interface
    ! Dependencies: Wind speed, sea state, and gas-specific Schmidt number
    ! Units: m/s
    ! =====================================================================
    if (use_REcoM) then
        call def_stream(nod2D, myDim_nod2D, &
                       'Kw', 'Air-sea piston velocity', 'm/s', PistonVelocity(:), &
                       io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    end if

CASE ('dpCO2s    ')
    ! =====================================================================
    ! Variable: dpCO2s (Delta pCO2)
    ! Description: Air-sea pCO2 gradient (ocean minus atmosphere)
    ! Function: Driving force for CO2 flux; positive = ocean supersaturated
    ! Sign convention: Positive when ocean releases CO2 to atmosphere
    ! Units: μatm (microatmospheres)
    ! =====================================================================
    if (use_REcoM) then
        call def_stream(nod2D, myDim_nod2D, &
                       'dpCO2s', 'Difference of oceanic pCO2 minus atmospheric pCO2', 'uatm', GlodPCO2surf(:), &
                       io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    end if

CASE ('pCO2s     ')
    ! =====================================================================
    ! Variable: pCO2s
    ! Description: Oceanic partial pressure of CO2 at the surface
    ! Function: Measures CO2 concentration in surface waters
    ! Context: Atmospheric pCO2 ≈ 420 μatm (as of 2024)
    ! Units: μatm (microatmospheres)
    ! =====================================================================
    if (use_REcoM) then
        call def_stream(nod2D, myDim_nod2D, &
                        'pCO2s', 'Partial pressure of oceanic CO2', 'uatm', GloPCO2surf(:), &
                        io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    end if

CASE ('CO2f      ')
    ! =====================================================================
    ! Variable: CO2f
    ! Description: Net CO2 flux across the air-sea interface
    ! Function: Quantifies ocean carbon uptake (or release)
    ! Sign convention: Positive = flux into ocean (ocean is a CO2 sink)
    ! Calculation: Flux = Kw × alphaCO2 × dpCO2s
    ! Units: mmolC/m²/d (millimoles carbon per square meter per day)
    ! =====================================================================
    if (use_REcoM) then
        call def_stream(nod2D, myDim_nod2D, &
                       'CO2f', 'CO2-flux into the surface water', 'mmolC/m2/d', GloCO2flux(:), &
                       io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    end if

CASE ('O2f       ')
    ! =====================================================================
    ! Variable: O2f
    ! Description: Net O2 flux across the air-sea interface
    ! Function: Balances photosynthetic O2 production and respiratory consumption
    ! Sign convention: Positive = flux into ocean
    ! Context: Opposite sign to CO2 flux in productive regions
    ! Units: mmolO/m²/d (millimoles oxygen per square meter per day)
    ! =====================================================================
    if (use_REcoM) then
        call def_stream(nod2D, myDim_nod2D, &
                       'O2f', 'O2-flux into the surface water', 'mmolO/m2/d', GloO2flux(:), &
                       io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    end if

CASE ('Hp        ')
    ! =====================================================================
    ! Variable: Hp (H⁺ concentration)
    ! Description: Hydrogen ion concentration in surface water
    ! Function: Direct measure of ocean acidity (pH = -log₁₀[H⁺])
    ! Context: pH ≈ 8.1 in modern ocean → [H⁺] ≈ 10⁻⁸·¹ mol/kg
    ! Relevance: Ocean acidification monitoring
    ! Units: mol/kg
    ! =====================================================================
    if (use_REcoM) then
        call def_stream(nod2D, myDim_nod2D, &
                       'Hp', 'Mean of H-plus ions in the surface water', 'mol/kg', GloHplus(:), &
                       io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    end if

! =====================================================================
! ATMOSPHERIC DEPOSITION
! =====================================================================
! External nutrient inputs to the ocean surface from atmospheric sources
! =====================================================================

CASE ('aFe       ')
    ! =====================================================================
    ! Variable: aFe
    ! Description: Atmospheric iron deposition flux
    ! Function: Primary iron source for open ocean phytoplankton
    ! Sources: Mineral dust (especially from deserts), combustion aerosols
    ! Importance: Iron limits primary production in HNLC regions
    ! Units: μmolFe/m²/s (micromoles iron per square meter per second)
    ! =====================================================================
    if (use_REcoM) then
        call def_stream(nod2D, myDim_nod2D, &
                       'aFe', 'Atmospheric iron input', 'umolFe/m2/s', AtmFeInput(:), &
                       io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    end if

CASE ('aN        ')
    ! =====================================================================
    ! Variable: aN
    ! Description: Atmospheric dissolved inorganic nitrogen (DIN) deposition
    ! Function: External nitrogen input to surface ocean
    ! Sources: NOₓ from combustion, ammonia from agriculture, natural sources
    ! Forms: Primarily nitrate (NO₃⁻) and ammonium (NH₄⁺)
    ! Units: mmolN/m²/s (millimoles nitrogen per square meter per second)
    ! =====================================================================
    if (use_REcoM) then
        call def_stream(nod2D, myDim_nod2D, &
                       'aN', 'Atmospheric DIN input', 'mmolN/m2/s', AtmNInput(:), &
                       io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    end if

CASE ('riverDIN    ')
    if (use_REcoM) then
    call def_stream(nod2D,  myDim_nod2D,   'riverDIN','riverine DIN input','mmolN/m2/s', RiverDIN2D(:), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    end if
    
CASE ('riverDIC    ')
    if (use_REcoM) then
    call def_stream(nod2D,  myDim_nod2D,   'riverDIC','riverine DIC input','mmolC/m2/s', RiverDIC2D(:), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    end if

CASE ('riverDOCsl    ')
    if (use_REcoM) then
    call def_stream(nod2D,  myDim_nod2D,   'riverDOCsl','riverine semi-labile DOC input','mmolC/m2/s', RiverDOCsl2D(:), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    end if

CASE ('riverDOCl    ')
    if (use_REcoM) then
    call def_stream(nod2D,  myDim_nod2D,   'riverDOCl','riverine labile DOC input','mmolC/m2/s', RiverDOCl2D(:), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    end if

CASE ('riverPOC    ')
    if (use_REcoM) then
    call def_stream(nod2D,  myDim_nod2D,   'riverPOC','riverine POC input','mmolC/m2/s', RiverPOC2D(:), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    end if

CASE ('riverDFe    ')
    if (use_REcoM) then
    call def_stream(nod2D,  myDim_nod2D,   'riverDFe','riverine DFe input','mmolC/m2/s', RiverFe(:), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    end if


! =====================================================================
! BENTHIC (SEAFLOOR) POOLS
! =====================================================================
! Accumulation of particulate organic and inorganic material in sediments
! These pools act as temporary storage before remineralization/burial
! =====================================================================

CASE ('benN      ')
    ! =====================================================================
    ! Variable: benN
    ! Description: Benthic nitrogen pool
    ! Function: Sedimentary organic nitrogen awaiting remineralization
    ! Process: Sinking particles → benthos → remineralization → DIN
    ! Role: Nutrient recycling in shallow waters; burial in deep ocean
    ! Units: mmol (millimoles per grid cell)
    ! =====================================================================
    if (use_REcoM) then
        call def_stream(nod2D, myDim_nod2D, &
                       'benN', 'Benthos Nitrogen', 'mmol', Benthos(:,1), &
                       io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    end if

CASE ('benC      ')
    ! =====================================================================
    ! Variable: benC
    ! Description: Benthic carbon pool
    ! Function: Sedimentary organic carbon storage
    ! Process: POC sinking → benthic accumulation → remineralization/burial
    ! Importance: Long-term carbon sequestration pathway
    ! Units: mmol (millimoles per grid cell)
    ! =====================================================================
    if (use_REcoM) then
        call def_stream(nod2D, myDim_nod2D, &
                       'benC', 'Benthos Carbon', 'mmol', Benthos(:,2), &
                        io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    end if

CASE ('benSi     ')
    ! =====================================================================
    ! Variable: benSi
    ! Description: Benthic silicon pool (biogenic silica)
    ! Function: Opal (SiO₂) from diatom frustules in sediments
    ! Process: Diatom mortality → opal dissolution in water column/benthos
    ! Fate: Slow dissolution returns Si to dissolved silicate pool
    ! Units: mmol (millimoles per grid cell)
    ! =====================================================================
    if (use_REcoM) then
        call def_stream(nod2D, myDim_nod2D, &
                       'benSi', 'Benthos silicon', 'mmol', Benthos(:,3), &
                       io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    end if

CASE ('benCalc   ')
    ! =====================================================================
    ! Variable: benCalc
    ! Description: Benthic calcite pool (calcium carbonate)
    ! Function: CaCO₃ from coccolithophore shells and other sources
    ! Process: Calcifier mortality → CaCO₃ sinking → accumulation/dissolution
    ! Importance: Carbonate chemistry buffering; geological carbon burial
    ! Units: mmol (millimoles per grid cell)
    ! =====================================================================
    if (use_REcoM) then
        call def_stream(nod2D,  myDim_nod2D, &
                       'benCalc','Benthos calcite','mmol', Benthos(:,4), &
                       io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    end if

CASE ('denb      ')
    ! Variable: denb
    ! Description: Benthic denitrification rate
    ! Function: Conversion of nitrate to N₂ gas in anoxic sediments
    ! Process: NO₃⁻ → NO₂⁻ → NO → N₂O → N₂ (anaerobic respiration pathway)
    ! Role: Permanent nitrogen loss from ocean; regulates ocean N inventory
    ! Units: mmol/m² (millimoles per square meter)
    if (use_REcoM) then
        call def_stream(nod2D, myDim_nod2D, &
                       'denb', 'Benthic denitrification rate', 'mmol/m2', DenitBen(:), &
                       io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    end if

! =====================================================================
! BURIAL FLUXES
! =====================================================================
CASE ('BurialN      ')
    ! Variable: BurialN
    ! Description: Nitrogen burial rate
    ! Function: Permanent removal of nitrogen from biogeochemical cycling
    ! Process: Benthic N → deep sediments → geological timescale storage
    ! Role: Long-term nitrogen sink; reduces ocean fixed-N inventory
    ! Units: mmolN/m² (millimoles nitrogen per square meter)
    if (use_REcoM) then
        call def_stream(nod2D, myDim_nod2D, &
                       'BurialN', 'Benthic Burial rate of Nitrogen', 'mmolN/m2', Burial(1,:), &
                       io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    end if

CASE ('BurialC      ')
    ! Variable: BurialC
    ! Description: Carbon burial rate
    ! Function: Permanent sequestration of organic carbon in sediments
    ! Process: Benthic C → deep sediments → long-term CO₂ drawdown
    ! Role: Major carbon sink; climate regulation on geological timescales
    ! Units: mmolC/m² (millimoles carbon per square meter)
    if (use_REcoM) then
        call def_stream(nod2D, myDim_nod2D, &
                       'BurialC', 'Benthic Burial rate of Carbon', 'mmolC/m2', Burial(2,:), &
                       io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    end if

CASE ('BurialSi     ')
    ! Variable: BurialSi
    ! Description: Biogenic silica burial rate
    ! Function: Permanent removal of opal from dissolution-precipitation cycle
    ! Process: Benthic opal → deep sediments → geological storage
    ! Role: Silicon sink; sedimentary opal record of diatom productivity
    ! Units: mmolSi/m² (millimoles silicon per square meter)
    if (use_REcoM) then
        call def_stream(nod2D, myDim_nod2D, &
                       'BurialSi', 'Benthic Burial rate', 'mmolSi/m2', Burial(3,:), &
                       io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    end if

CASE ('BurialCalc     ')
    ! Variable: BurialCalc
    ! Description: Calcite burial rate
    ! Function: Permanent sequestration of CaCO₃ in sediments
    ! Process: Benthic calcite → deep sediments → limestone formation
    ! Role: Long-term alkalinity/DIC sink; carbonate compensation depth control
    ! Units: mmolC/m² (millimoles carbon as CaCO₃ per square meter)
    if (use_REcoM) then
        call def_stream(nod2D, myDim_nod2D, &
                       'BurialCalc', 'Benthic Burial rate of Calcite', 'mmolC/m2', Burial(4,:), &
                       io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    end if

! =====================================================================
! OCEAN → SEDIMENT FLUXES
! =====================================================================
CASE ('wFNb     ')
    ! Variable: wFNb
    ! Description: Sinking particulate nitrogen flux to benthos
    ! Function: Downward transport of organic nitrogen to seafloor
    ! Process: Detritus/phytoplankton → sinking → benthic deposition
    ! Role: Primary nitrogen input to sediments; fuels benthic remineralization
    ! Units: mmolN/(m²·d) (millimoles nitrogen per square meter per day)
    if (use_REcoM) then
        call def_stream(nod2D, myDim_nod2D, &
                       'wFNb', 'Sinking flux into benthic N', 'mmolN/(m2*d)', Ocean_2_Sed_flux(:,1), &
                       io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    end if

CASE ('wFCb     ')
    ! Variable: wFCb
    ! Description: Sinking particulate organic carbon flux to benthos
    ! Function: Downward transport of organic carbon to seafloor
    ! Process: POC/phytoplankton → sinking → benthic deposition
    ! Role: Primary carbon input to sediments; energy source for benthos
    ! Units: mmolC/(m²·d) (millimoles carbon per square meter per day)
    if (use_REcoM) then
        call def_stream(nod2D, myDim_nod2D, &
                       'wFCb', 'Sinking flux into benthic C', 'mmolC/(m2*d)', Ocean_2_Sed_flux(:,2), &
                       io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    end if

CASE ('wFSib     ')
    ! Variable: wFSib
    ! Description: Sinking biogenic silica flux to benthos
    ! Function: Downward transport of opal (diatom frustules) to seafloor
    ! Process: Diatom detritus → sinking → benthic deposition
    ! Role: Primary silicon input to sediments; creates opal-rich sediments
    ! Units: mmolSi/(m²·d) (millimoles silicon per square meter per day)
    if (use_REcoM) then
        call def_stream(nod2D, myDim_nod2D, &
                       'wFSib', 'Sinking flux into benthic Si', 'mmolSi/(m2*d)', Ocean_2_Sed_flux(:,3), &
                       io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    end if

CASE ('wFCalcb     ')
    ! Variable: wFCalcb
    ! Description: Sinking calcite flux to benthos
    ! Function: Downward transport of CaCO₃ (coccoliths, foraminifera) to seafloor
    ! Process: Calcifying organisms → sinking → benthic deposition
    ! Role: Primary carbonate input; lysocline/CCD depth control
    ! Units: mmolC/(m²·d) (millimoles carbon as CaCO₃ per square meter per day)
    if (use_REcoM) then
        call def_stream(nod2D, myDim_nod2D, &
                       'wFCalcb', 'Sinking flux into benthic Carb', 'mmolC/(m2*d)', Ocean_2_Sed_flux(:,4), &
                       io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    end if

! =====================================================================
! SEDIMENT → OCEAN FLUXES (REMINERALIZATION/DISSOLUTION)
! =====================================================================
CASE ('ReNb     ')
    ! Variable: ReNb
    ! Description: Benthic nitrogen remineralization flux
    ! Function: Release of dissolved inorganic nitrogen from sediments
    ! Process: Benthic organic N → bacterial degradation → NH₄⁺/NO₃⁻ → ocean
    ! Role: Nutrient regeneration; supports primary production in shallow waters
    ! Units: mmolN/(m²·s) (millimoles nitrogen per square meter per second)
    if (use_REcoM) then
        call def_stream(nod2D, myDim_nod2D, &
                       'ReNb', 'benthic N release to the ocean (remineralization)', 'mmolN/(m2*s)', &
                       Sed_2_Ocean_Flux(:,1), &
                       io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    end if

CASE ('ReCb     ')
    ! Variable: ReCb
    ! Description: Benthic carbon remineralization + calcite dissolution flux
    ! Function: Release of dissolved inorganic carbon from sediments
    ! Process: Organic C → remineralization → DIC + CaCO₃ dissolution → ocean
    ! Role: CO₂ source; carbonate system buffering; alkalinity regulation
    ! Units: mmolC/(m²·s) (millimoles carbon per square meter per second)
    if (use_REcoM) then
        call def_stream(nod2D, myDim_nod2D, &
                       'ReCb', 'benthic C release to the ocean (remineralization + calcification)', &
                       'mmolC/(m2*s)', Sed_2_Ocean_Flux(:,2), &
                       io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    end if

CASE ('ReSib     ')
    ! Variable: ReSib
    ! Description: Benthic silica dissolution flux
    ! Function: Release of dissolved silicon from opal dissolution
    ! Process: Benthic opal → porewater dissolution → Si(OH)₄ → ocean
    ! Role: Silicon regeneration; supports diatom productivity
    ! Units: mmolSi/(m²·s) (millimoles silicon per square meter per second)
    if (use_REcoM) then
        call def_stream(nod2D, myDim_nod2D, &
                       'ReSib', 'benthic Si release to the ocean (dissolution)', 'mmolSi/(m2*s)', &
                       Sed_2_Ocean_Flux(:,4), &
                       io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    end if

CASE ('ReAlkb     ')
    ! Variable: ReAlkb
    ! Description: Benthic alkalinity flux
    ! Function: Release of alkalinity from organic matter and CaCO₃ processing
    ! Process: Denitrification + sulfate reduction + CaCO₃ dissolution → Alk
    ! Role: pH buffering; carbonate system regulation; CO₂ uptake capacity
    ! Units: mmolC-eq/(m²·s) (millimole carbon-equivalents per square meter per second)
    if (use_REcoM) then
        call def_stream(nod2D, myDim_nod2D, &
                       'ReAlkb', 'benthic Alk release to the ocean', 'mmolC-eq/(m2*s)', &
                       Sed_2_Ocean_Flux(:,3), &
                       io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    end if

CASE ('ReFeb     ')
    ! Variable: ReFeb
    ! Description: Benthic iron flux
    ! Function: Release of dissolved iron from reducing sediments
    ! Process: Benthic Fe-oxides → anoxic reduction → Fe²⁺ → ocean
    ! Role: Micronutrient supply; supports primary production in Fe-limited regions
    ! Units: mmolFe/(m²·s) (millimoles iron per square meter per second)
    if (use_REcoM) then
        call def_stream(nod2D, myDim_nod2D, &
                       'ReFeb', 'benthic Fe release to the ocean', 'mmolFe/(m2*s)', &
                       Sed_2_Ocean_Flux(:,5), &
                       io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    end if

CASE ('ReO2b     ')
    ! Variable: ReO2b
    ! Description: Benthic oxygen consumption flux
    ! Function: Oxygen demand from aerobic remineralization in sediments
    ! Process: Organic matter + O₂ → CO₂ + nutrients (negative flux = consumption)
    ! Role: Sediment oxygen demand; creates anoxic conditions; hypoxia indicator
    ! Units: mmolO₂/(m²·s) (millimoles oxygen per square meter per second)
    if (use_REcoM) then
        call def_stream(nod2D, myDim_nod2D, &
                       'ReO2b', 'benthic O2 release to the ocean', 'mmolO2/(m2*s)', &
                       Sed_2_Ocean_Flux(:,6), &
                       io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    end if

! =====================================================================
! PHYTOPLANKTON PRIMARY PRODUCTION
! =====================================================================
! Production and nutrient assimilation by different phytoplankton groups
! NPP = Net Primary Production (photosynthesis - respiration)
! GPP = Gross Primary Production (total photosynthesis)
! NNA = Net Nitrogen Assimilation (nitrogen uptake - excretion)
! ChlDeg = Chlorophyll degradation rate
! =====================================================================

! ---------------------------------------------------------------------
! SMALL PHYTOPLANKTON
! ---------------------------------------------------------------------
! Small phytoplankton (typically < 20 μm; e.g., small flagellates)
! Ecological role: Dominant in oligotrophic (nutrient-poor) regions
! ---------------------------------------------------------------------

CASE ('NPPn      ')
    ! =====================================================================
    ! Variable: NPPn
    ! Description: Net Primary Production of nanophytoplankton
    ! Function: Carbon fixation minus autotrophic respiration
    ! Integration: Vertically integrated over water column
    ! Units: mmolC/m²/d (millimoles carbon per square meter per day)
    ! =====================================================================
    if (use_REcoM) then
        call def_stream(nod2D, myDim_nod2D, &
                       'NPPn', 'Mean NPP nanophytoplankton', 'mmolC/m2/d', NPPn, &
                       io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    end if

CASE ('GPPn     ')
    ! =====================================================================
    ! Variable: GPPn
    ! Description: Gross Primary Production of nanophytoplankton
    ! Function: Total photosynthetic carbon fixation (before respiration)
    ! Relationship: GPPn = NPPn + autotrophic respiration
    ! Units: mmolC/m²/d
    ! =====================================================================
    if (use_REcoM) then
        call def_stream(nod2D, myDim_nod2D, &
                       'GPPn', 'Mean GPP nanophytoplankton', 'mmolC/m2/d', GPPn, &
                       io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    end if

CASE ('NNAn     ')
    ! =====================================================================
    ! Variable: NNAn
    ! Description: Net Nitrogen Assimilation by nanophytoplankton
    ! Function: Nitrogen uptake (NO₃⁻, NH₄⁺) minus excretion
    ! Stoichiometry: Related to NPPn via Redfield ratio (C:N ≈ 6.6:1)
    ! Units: mmolN/m²/d
    ! =====================================================================
    if (use_REcoM) then
        call def_stream(nod2D, myDim_nod2D, & 
                       'NNAn', 'Net N-assimilation nanophytoplankton', 'mmolN/m2/d', NNAn, &
                       io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    end if

CASE ('Chldegn  ')
    ! =====================================================================
    ! Variable: Chldegn
    ! Description: Chlorophyll degradation rate for nanophytoplankton
    ! Function: Loss rate of chlorophyll (photopigment bleaching)
    ! Causes: Photo-oxidation, senescence, nutrient stress
    ! Application: Chlorophyll evolution = synthesis - degradation
    ! Units: 1/d (per day; first-order rate constant)
    ! =====================================================================
    if (use_REcoM) then
        call def_stream(nod2D, myDim_nod2D, & 
                       'ChlDegn', 'Chlorophyll degradation nanophytoplankton', '1/d', Chldegn, &
                       io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    end if

! ---------------------------------------------------------------------
! DIATOMS
! ---------------------------------------------------------------------
! Silica-shelled phytoplankton (typically 20-200 μm)
! Ecological role: Bloom formers in nutrient-rich upwelling regions
! Requirement: Need silicate (Si(OH)₄) for frustule construction
! ---------------------------------------------------------------------

CASE ('NPPd     ')
    ! =====================================================================
    ! Variable: NPPd
    ! Description: Net Primary Production of diatoms
    ! Function: Diatom carbon fixation minus respiration
    ! Ecological context: Dominates spring blooms, coastal upwelling
    ! Units: mmolC/m²/d
    ! =====================================================================
    if (use_REcoM) then
        call def_stream(nod2D,  myDim_nod2D, &
                       'NPPd', 'Mean NPP diatoms', 'mmolC/m2/d', NPPd, &
                       io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    end if

CASE ('GPPd     ')
    ! =====================================================================
    ! Variable: GPPd
    ! Description: Gross Primary Production of diatoms
    ! Function: Total diatom photosynthesis
    ! Importance: Major contributor to biological carbon pump
    ! Units: mmolC/m²/d
    ! =====================================================================
    if (use_REcoM) then
        call def_stream(nod2D, myDim_nod2D, &
                       'GPPd', 'Mean GPP diatoms', 'mmolC/m2/d', GPPd, &
                       io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    end if

CASE ('NNAd     ')
    ! =====================================================================
    ! Variable: NNAd
    ! Description: Net Nitrogen Assimilation by diatoms
    ! Function: Diatom nitrogen uptake
    ! Preference: Often prefer NH₄⁺ but can use NO₃⁻
    ! Units: mmolN/m²/d
    ! =====================================================================
    if (use_REcoM) then
        call def_stream(nod2D, myDim_nod2D, &
                       'NNAd', 'Net N-assimilation diatoms', 'mmolN/m2/d', NNAd, &
                       io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    end if

CASE ('Chldegd  ')
    ! =====================================================================
    ! Variable: Chldegd
    ! Description: Chlorophyll degradation rate for diatoms
    ! Function: Diatom chlorophyll loss rate
    ! Context: Can be elevated during silicate stress
    ! Units: 1/d
    ! =====================================================================
    if (use_REcoM) then
        call def_stream(nod2D, myDim_nod2D, &
                       'ChlDegd', 'Chlorophyll degradation diatoms', '1/d', Chldegd, &
                       io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    end if



! ---------------------------------------------------------------------
! COCCOLITHOPHORES
! ---------------------------------------------------------------------
! Calcite-producing phytoplankton (e.g., Emiliania huxleyi)
! Ecological role: CaCO₃ production affects alkalinity and air-sea CO₂
! Distribution: Temperate and subtropical regions
! ---------------------------------------------------------------------

CASE ('NPPc     ')
    ! =====================================================================
    ! Variable: NPPc
    ! Description: Net Primary Production of coccolithophores
    ! Function: Coccolith carbon fixation
    ! Special feature: Both organic C (photosynthesis) and inorganic C (CaCO₃)
    ! Units: mmolC/(m²*d)
    ! =====================================================================
    if (use_REcoM) then
        call def_stream(nod2D, myDim_nod2D, &
                       'NPPc', 'Mean NPP coccolithophores', 'mmolC/(m2*d)', NPPc, &
                       io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    end if

CASE ('GPPc     ')
    ! =====================================================================
    ! Variable: GPPc
    ! Description: Gross Primary Production of coccolithophores
    ! Function: Total coccolith photosynthesis
    ! Climate impact: Calcification increases surface pCO₂ (CO₂ source)
    ! Units: mmolC/(m²*d)
    ! =====================================================================
    if (use_REcoM) then
        call def_stream(nod2D, myDim_nod2D, &
                       'GPPc', 'Mean GPP coccolithophores','mmolC/(m2*d)', GPPc, &
                       io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    end if

CASE ('NNAc     ')
    ! =====================================================================
    ! Variable: NNAc
    ! Description: Net Nitrogen Assimilation by coccolithophores
    ! Function: Coccolith nitrogen uptake
    ! Competition: Less competitive for nutrients than diatoms
    ! Units: mmolN/(m²*d)
    ! =====================================================================
    if (use_REcoM) then
        call def_stream(nod2D, myDim_nod2D, &
                       'NNAc', 'Net N-assimilation coccolithophores', 'mmolN/(m2*d)', NNAc, &
                       io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    endif

CASE ('Chldegc  ')
    ! =====================================================================
    ! Variable: Chldegc
    ! Description: Chlorophyll degradation rate for coccolithophores
    ! Function: Coccolith chlorophyll loss
    ! Units: 1/d
    ! =====================================================================
    if (use_REcoM) then
        call def_stream(nod2D, myDim_nod2D, &  
                       'ChlDegc', 'Chlorophyll degradation coccolithophores', '1/d', Chldegc, &
                       io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    endif

! ---------------------------------------------------------------------
! PHAEOCYSTIS
! ---------------------------------------------------------------------
! Colony-forming phytoplankton (Phaeocystis spp.)
! Ecological role: Forms massive blooms in polar/subpolar regions
! Unique feature: Can exist as single cells or large mucilaginous colonies
! ---------------------------------------------------------------------

CASE ('NPPp     ')
    ! =====================================================================
    ! Variable: NPPp
    ! Description: Net Primary Production of Phaeocystis
    ! Function: Phaeocystis carbon fixation
    ! Regional importance: Dominant in Southern Ocean, North Sea blooms
    ! Units: mmolC/(m²*d)
    ! =====================================================================
    if (use_REcoM) then
        call def_stream(nod2D, myDim_nod2D, &
                       'NPPp', 'Mean NPP phaeocystis', 'mmolC/(m2*d)', NPPp, &
                       io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    end if

CASE ('GPPp     ')
    ! =====================================================================
    ! Variable: GPPp
    ! Description: Gross Primary Production of Phaeocystis
    ! Function: Total Phaeocystis photosynthesis
    ! Biogeochemical role: Colonial form less grazed; affects food web
    ! Units: mmolC/(m²*d)
    ! =====================================================================
    if (use_REcoM) then
        call def_stream(nod2D, myDim_nod2D, &
                       'GPPp', 'Mean GPP phaeocystis', 'mmolC/(m2*d)', GPPp, &
                       io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    end if

CASE ('NNAp     ')
    ! =====================================================================
    ! Variable: NNAp
    ! Description: Net Nitrogen Assimilation by Phaeocystis
    ! Function: Phaeocystis nitrogen uptake
    ! Nutrient preference: Can efficiently use various N forms
    ! Units: mmolN/(m²*d)
    ! =====================================================================
    if (use_REcoM) then
        call def_stream(nod2D, myDim_nod2D, &
                       'NNAp', 'Net N-assimilation phaeocystis', 'mmolN/(m2*d)', NNAp, &
                       io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    endif

CASE ('Chldegp  ')
    ! =====================================================================
    ! Variable: Chldegp
    ! Description: Chlorophyll degradation rate for Phaeocystis
    ! Function: Phaeocystis chlorophyll loss
    ! Units: 1/d
    ! =====================================================================
    if (use_REcoM) then
        call def_stream(nod2D, myDim_nod2D, &
                       'ChlDegp', 'Chlorophyll degradation phaeocystis', '1/d', Chldegp, &
                       io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    endif

! =====================================================================
! ZOOPLANKTON GRAZING FLUXES
! =====================================================================
! Predation fluxes by three zooplankton functional groups
! 
! IMPORTANT DISTINCTION:
! - Fluxes "without grazing efficiency" = loss from prey's perspective
!   (all carbon removed from prey pool)
! - "Total grazing" = actual ingestion by predator
!   (accounts for grazing efficiency < 1; rest becomes detritus)
! 
! Model structure: Three zooplankton size classes
!   - Microzooplankton: Small (e.g., ciliates, heterotrophic flagellates)
!   - Mesozooplankton: Medium (e.g., copepods)
!   - Macrozooplankton: Large (e.g., euphausiids, large copepods)
! =====================================================================

! ---------------------------------------------------------------------
! MESOZOOPLANKTON GRAZING
! ---------------------------------------------------------------------
! Medium-sized zooplankton (typically copepods, 0.2-2 mm)
! Diet: Omnivorous - phytoplankton, microzooplankton, detritus
! Ecological role: Link between primary producers and fish
! ---------------------------------------------------------------------

CASE ('grazmeso_tot')
    ! =====================================================================
    ! Variable: grazmeso_tot
    ! Description: Total mesozooplankton grazing (carbon ingested)
    ! Function: Sum of all prey items × grazing efficiency
    ! Fate of ingested C: Growth + respiration + excretion + egestion
    ! Units: mmolC/(m²d)
    ! =====================================================================
    call def_stream(nod2D, myDim_nod2D, &
                   'grazmeso_tot', 'Total grazing flux of mesozooplankton, dependent on grazing efficiency', 'mmolC/(m2d)', grazmeso_tot, &
                   io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)


CASE ('grazmeso_n')
    ! =====================================================================
    ! Variable: grazmeso_n
    ! Description: Mesozooplankton grazing on nanophytoplankton
    ! Function: Carbon loss from nanophyto pool due to mesozoo predation
    ! Note: Includes both ingested and sloppy feeding (→ DOM/detritus)
    ! Units: mmolC/(m²d)
    ! =====================================================================
    call def_stream(nod2D, myDim_nod2D, &
                   'grazmeso_n', 'Grazing flux of mesozooplankton on small phytoplankton without grazing efficiency (i.e., = loss small phytoplankton)', 'mmolC/(m2d)', grazmeso_n, &
                   io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)

CASE ('grazmeso_d')
    ! =====================================================================
    ! Variable: grazmeso_d
    ! Description: Mesozooplankton grazing on diatoms
    ! Function: Carbon loss from diatom pool
    ! Ecological note: Copepods are major diatom grazers in blooms
    ! Units: mmolC/(m²d)
    ! =====================================================================
    call def_stream(nod2D, myDim_nod2D, &
                   'grazmeso_d', 'Grazing flux of mesozooplankton on diatoms without grazing efficiency (i.e., = loss diatoms)', 'mmolC/(m2d)', grazmeso_d, &
                   io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)

CASE ('grazmeso_c')
    ! =====================================================================
    ! Variable: grazmeso_c
    ! Description: Mesozooplankton grazing on coccolithophores
    ! Function: Carbon loss from coccolithophore pool
    ! Units: mmolC/(m²d)
    ! =====================================================================
    call def_stream(nod2D, myDim_nod2D, &
                   'grazmeso_c', 'Grazing flux of mesozooplankton on coccolithophores without grazing efficiency (i.e., = loss coccolithophores)', 'mmolC/(m2d)', grazmeso_c, &
                   io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)

CASE ('grazmeso_p')
    ! =====================================================================
    ! Variable: grazmeso_p
    ! Description: Mesozooplankton grazing on Phaeocystis
    ! Function: Carbon loss from Phaeocystis pool
    ! Ecological note: Colonial Phaeocystis often poorly grazed
    ! Units: mmolC/(m²*d)
    ! =====================================================================
    call def_stream(nod2D, myDim_nod2D, &
                   'grazmeso_p', 'Grazing flux of mesozooplankton on phaeocystis without grazing efficiency (i.e., = loss phaeocystis)', 'mmolC/(m2*d)', grazmeso_p, &
                   io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)

CASE ('grazmeso_det')
    ! =====================================================================
    ! Variable: grazmeso_det
    ! Description: Mesozooplankton grazing on first detritus group
    ! Function: Detritivory by mesozooplankton
    ! Ecological role: Repackaging of detritus into fecal pellets
    ! Units: mmolC/(m²d)
    ! =====================================================================
    call def_stream(nod2D, myDim_nod2D, &
                   'grazmeso_det','Grazing flux of mesozooplankton on first detritus group without grazing efficiency (i.e., = loss first detritus group)', 'mmolC/(m2d)', grazmeso_det, &
                   io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)

CASE ('grazmeso_mic')
    ! =====================================================================
    ! Variable: grazmeso_mic
    ! Description: Mesozooplankton predation on microzooplankton
    ! Function: Trophic transfer from micro- to mesozooplankton
    ! Ecological significance: Microbial loop → classical food chain link
    ! Units: mmolC/(m²d)
    ! =====================================================================
    call def_stream(nod2D, myDim_nod2D, &
                   'grazmeso_mic', 'Grazing flux of mesozooplankton on microzooplankton without grazing efficiency (i.e., = loss microzooplankton)', 'mmolC/(m2d)', grazmeso_mic, &
                   io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)

CASE ('grazmeso_det2')
    ! =====================================================================
    ! Variable: grazmeso_det2
    ! Description: Mesozooplankton grazing on second detritus group
    ! Function: Feeding on different detritus pool (possibly larger/slower sinking)
    ! Model structure: Multiple detritus pools represent different size classes
    ! Units: mmolC/(m²*d)
    ! =====================================================================
    call def_stream(nod2D,  myDim_nod2D, &
                   'grazmeso_det2', 'Grazing flux of mesozooplankton on first detritus without grazing efficiency (i.e., = loss second detritus group)', 'mmolC/(m2*d)', grazmeso_det2, &
                   io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)

! ---------------------------------------------------------------------
! MACROZOOPLANKTON GRAZING
! ---------------------------------------------------------------------
! Large zooplankton (euphausiids/krill, large copepods, >2 mm)
! Diet: Phytoplankton, smaller zooplankton, detritus
! Ecological role: Major prey for fish, seabirds, marine mammals
! ---------------------------------------------------------------------

CASE ('grazmacro_tot')
    ! =====================================================================
    ! Variable: grazmacro_tot
    ! Description: Total macrozooplankton grazing (carbon ingested)
    ! Function: Total carbon intake by macrozooplankton
    ! Ecological importance: Critical for upper trophic levels
    ! Units: mmolC/(m²d)
    ! =====================================================================
    call def_stream(nod2D, myDim_nod2D, &
                   'grazmacro_tot', 'Total grazing flux of macrozooplankton, dependent on grazing efficiency', 'mmolC/(m2d)', grazmacro_tot, &
                   io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)

CASE ('grazmacro_n')
    ! =====================================================================
    ! Variable: grazmacro_n
    ! Description: Macrozooplankton grazing on nanophytoplankton
    ! Function: Loss from nanophyto due to macrozooplankton
    ! Units: mmolC/(m²d)
    ! =====================================================================
    call def_stream(nod2D, myDim_nod2D, &
                   'grazmacro_n','Grazing flux of macrozooplankton on small phytoplankton without grazing efficiency (i.e., = loss small phytoplankton)', 'mmolC/(m2d)', grazmacro_n, &
                    io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)

CASE ('grazmacro_d')
    ! =====================================================================
    ! Variable: grazmacro_d
    ! Description: Macrozooplankton grazing on diatoms
    ! Function: Loss from diatoms due to macrozooplankton
    ! Example: Antarctic krill feeding on diatom blooms
    ! Units: mmolC/(m²d)
    ! =====================================================================
    call def_stream(nod2D, myDim_nod2D, &
                   'grazmacro_d', 'Grazing flux of macrozooplankton on diatoms without grazing efficiency (i.e., = loss diatoms)','mmolC/(m2d)', grazmacro_d, &
                   io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)

CASE ('grazmacro_c')
    ! =====================================================================
    ! Variable: grazmacro_c
    ! Description: Macrozooplankton grazing on coccolithophores
    ! Function: Loss from coccolithophores due to macrozooplankton
    ! Units: mmolC/(m²d)
    ! =====================================================================
    call def_stream(nod2D, myDim_nod2D, &
                   'grazmacro_c', 'Grazing flux of macrozooplankton on coccolithophores without grazing efficiency (i.e., = loss coccolithophores)','mmolC/(m2d)', grazmacro_c, &
                   io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)

CASE ('grazmacro_p')
    ! =====================================================================
    ! Variable: grazmacro_p
    ! Description: Macrozooplankton grazing on Phaeocystis
    ! Function: Loss from Phaeocystis due to macrozooplankton
    ! Units: mmolC/(m²*d)
    ! =====================================================================
    call def_stream(nod2D, myDim_nod2D, &
                   'grazmacro_p', 'Grazing flux of macrozooplankton on phaeocystis without grazing efficiency (i.e., = loss phaeocystis)', 'mmolC/(m2*d)', grazmacro_p, &
                   io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)

CASE ('grazmacro_mes')
    ! =====================================================================
    ! Variable: grazmacro_mes
    ! Description: Macrozooplankton predation on mesozooplankton
    ! Function: Carnivory - trophic transfer to higher level
    ! Ecological note: Intraguild predation within zooplankton
    ! Units: mmolC/(m²d)
    ! =====================================================================
    call def_stream(nod2D, myDim_nod2D, &
                   'grazmacro_mes', 'Grazing flux of mesozooplankton on macrozooplankton without grazing efficiency (i.e., = loss mesozooplankton)', 'mmolC/(m2d)', grazmacro_mes, &
                   io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)

CASE ('grazmacro_det')
    ! =====================================================================
    ! Variable: grazmacro_det
    ! Description: Macrozooplankton grazing on first detritus group
    ! Function: Detritivory by macrozooplankton
    ! Units: mmolC/(m²d)
    ! =====================================================================
    call def_stream(nod2D, myDim_nod2D, &
                   'grazmacro_det', 'Grazing flux of macrozooplankton on first detritus group without grazing efficiency (i.e., = loss first detritus group)', 'mmolC/(m2d)', grazmacro_det, &
                   io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)

CASE ('grazmacro_mic')
    ! =====================================================================
    ! Variable: grazmacro_mic
    ! Description: Macrozooplankton predation on microzooplankton
    ! Function: Loss from microzooplankton to macrozooplankton
    ! Units: mmolC/(m²d)
    ! =====================================================================
    call def_stream(nod2D, myDim_nod2D, &
                   'grazmacro_mic', 'Grazing flux of macrozooplankton on microzooplankton without grazing efficiency (i.e., = loss microzooplankton)', 'mmolC/(m2d)', grazmacro_mic, &
                   io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)

CASE ('grazmacro_det2')
    ! =====================================================================
    ! Variable: grazmacro_det2
    ! Description: Macrozooplankton grazing on second detritus group
    ! Function: Feeding on alternative detritus pool
    ! Units: mmolC/(m²d)
    ! =====================================================================
    call def_stream(nod2D, myDim_nod2D, &
                   'grazmacro_det2', 'Grazing flux of macrozooplankton on first detritus without grazing efficiency (i.e., = loss second detritus group)', 'mmolC/(m2d)', grazmacro_det2, &
                   io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)

! ---------------------------------------------------------------------
! MICROZOOPLANKTON GRAZING
! ---------------------------------------------------------------------
! Small zooplankton (ciliates, heterotrophic flagellates, <200 μm)
! Diet: Primarily nanophytoplankton and bacteria
! Ecological role: Microbial loop - recycles nutrients rapidly
! ---------------------------------------------------------------------

CASE ('grazmicro_tot')
    ! =====================================================================
    ! Variable: grazmicro_tot
    ! Description: Total microzooplankton grazing (carbon ingested)
    ! Function: Total carbon intake by microzooplankton
    ! Ecological importance: Dominant grazers in oligotrophic systems
    ! Units: mmolC/(m²d)
    ! =====================================================================
    call def_stream(nod2D, myDim_nod2D, &
                   'grazmicro_tot', 'Total grazing flux of microzooplankton, dependent on grazing efficiency', 'mmolC/(m2d)', grazmicro_tot, &
                   io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)

CASE ('grazmicro_n')
    ! =====================================================================
    ! Variable: grazmicro_n
    ! Description: Microzooplankton grazing on nanophytoplankton
    ! Function: Primary grazing pressure on small phytoplankton
    ! Ecological note: Can control nanophyto biomass in stratified waters
    ! Units: mmolC/(m²d)
    ! =====================================================================
    call def_stream(nod2D, myDim_nod2D, &
                   'grazmicro_n', 'Grazing flux of microzooplankton on small phytoplankton without grazing efficiency (i.e., = loss small phytoplankton)', 'mmolC/(m2d)', grazmicro_n, &
                   io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)

CASE ('grazmicro_d')
    ! =====================================================================
    ! Variable: grazmicro_d
    ! Description: Microzooplankton grazing on diatoms
    ! Function: Microzoo predation on diatoms (less common, size-limited)
    ! Note: Typically graze smaller diatom species/early stages
    ! Units: mmolC/(m²d)
    ! =====================================================================
    call def_stream(nod2D, myDim_nod2D, &
                   'grazmicro_d', 'Grazing flux of microzooplankton on diatoms without grazing efficiency (i.e., = loss diatoms)', 'mmolC/(m2d)', grazmicro_d, &
                   io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)

CASE ('grazmicro_c')
    ! =====================================================================
    ! Variable: grazmicro_c
    ! Description: Microzooplankton grazing on coccolithophores
    ! Function: Microzoo predation on coccolithophores
    ! Units: mmolC/(m²d)
    ! =====================================================================
    call def_stream(nod2D,  myDim_nod2D,   'grazmicro_c', &
                   'Grazing flux of microzooplankton on coccolithophores without grazing efficiency (i.e., = loss coccolithophores)', 'mmolC/(m2d)', grazmicro_c, &
                   io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)

CASE ('grazmicro_p')
    ! =====================================================================
    ! Variable: grazmicro_p
    ! Description: Microzooplankton grazing on Phaeocystis
    ! Function: Microzoo predation on Phaeocystis
    ! Note: Primarily graze solitary cells, not large colonies
    ! Units: mmolC/(m²*d)
    ! =====================================================================
    call def_stream(nod2D, myDim_nod2D, &
                   'grazmicro_p', 'Grazing flux of microzooplankton on phaeocystis without grazing efficiency (i.e., = loss phaeocystis)', 'mmolC/(m2*d)', grazmicro_p, &
                   io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)

CASE ('REMOC     ')
    if (use_REcoM) then
    call def_stream(nod2D,  myDim_nod2D,   'REMOC','Total remineralization of DOC','mmolN/(m2*d)', REMOC, io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    end if

CASE ('REMOCt     ')
    if (use_REcoM) then
    call def_stream(nod2D,  myDim_nod2D,   'REMOCt','Remineralization of terrigenous DOC','mmolC/(m2*d)', REMOCt, io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    end if
    
CASE ('REMON     ')
    if (use_REcoM) then
    call def_stream(nod2D,  myDim_nod2D,   'REMON','Total remineralization of DON','mmolC/(m2*d)', REMON, io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    end if

CASE ('DISSOSi     ')
    if (use_REcoM) then
    call def_stream(nod2D,  myDim_nod2D,   'DISSOSi','Dissolution of Si','mmolSi/(m2*d)', DISSOSi, io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    end if

CASE ('DISSOC     ')
    if (use_REcoM) then
    call def_stream(nod2D,  myDim_nod2D,   'DISSOC','Dissolution of POC','mmolC/(m2*d)', DISSOC, io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    end if

CASE ('DISSON     ')
    if (use_REcoM) then
    call def_stream(nod2D,  myDim_nod2D,   'DISSON','Dissolution of PON','mmolN/(m2*d)', DISSON, io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    end if

#endif

!___________________________________________________________________________________________________________________________________    
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>   3D streams   <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!___________________________________________________________________________________________________________________________________
CASE ('temp      ')
    call def_stream((/nl-1, nod2D/),  (/nl-1, myDim_nod2D/),  'temp',      'temperature', 'C',      tracers%data(1)%values(:,:),             io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
CASE ('salt      ')
    call def_stream((/nl-1, nod2D/),  (/nl-1, myDim_nod2D/),  'salt',      'salinity',    'psu',    tracers%data(2)%values(:,:),             io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)

#if defined(__recom)

! =====================================================================
! CARBONATE CHEMISTRY VARIABLES
! Description: Variables related to ocean carbon chemistry and pH
! =====================================================================

CASE ('PAR       ')
    ! =====================================================================
    ! Variable: PAR
    ! Description: Photosynthetically Active Radiation
    ! Function: Light energy available for photosynthesis (400-700nm wavelength)
    ! Role: Primary driver of phytoplankton growth
    ! Units: W/m²
    ! =====================================================================
    if (use_REcoM) then
        call def_stream((/nl-1, nod2D/), (/nl-1, myDim_nod2D/), &
                       'PAR', 'PAR', 'W/m2', PAR3D(:,:), & 
                       io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    end if

CASE ('CO2       ')
    ! =====================================================================
    ! Variable: CO2
    ! Description: Aqueous CO2 concentration
    ! Function: Dissolved CO2 in seawater, primary inorganic carbon form
    ! Role: Substrate for phytoplankton photosynthesis, pH regulator
    ! Units: mol/m³
    ! =====================================================================
    if (use_REcoM) then
        call def_stream((/nl-1, nod2D/), (/nl-1, myDim_nod2D/), &
                       'CO2', 'Aqueous CO2 concentration', 'mol/m3', CO23D(:,:), &
                       io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    end if

CASE ('pH        ')
    ! =====================================================================
    ! Variable: pH
    ! Description: Acidity/alkalinity of seawater (total scale)
    ! Function: Measure of hydrogen ion concentration
    ! Role: Affects carbonate chemistry, organism physiology, calcification
    ! Units: total scale (dimensionless, typically 7.5-8.5 in ocean)
    ! =====================================================================
    if (use_REcoM) then
        call def_stream((/nl-1, nod2D/), (/nl-1, myDim_nod2D/), &
                       'pH', 'pH', 'total scale', pH3D(:,:), &
                       io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    end if

CASE ('pCO2      ')
    ! =====================================================================
    ! Variable: pCO2
    ! Description: Partial pressure of CO2 in seawater
    ! Function: CO2 fugacity in equilibrium with atmosphere
    ! Role: Determines air-sea CO2 flux direction and magnitude
    ! Units: μatm (microatmospheres)
    ! =====================================================================
    if (use_REcoM) then
        call def_stream((/nl-1, nod2D/), (/nl-1, myDim_nod2D/), &
                       'pCO2', 'CO2 partial pressure', 'uatm', pCO23D(:,:), &
                       io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    end if

CASE ('HCO3      ')
    ! =====================================================================
    ! Variable: HCO3
    ! Description: Bicarbonate ion concentration
    ! Function: Major dissolved inorganic carbon species (~90% of DIC)
    ! Role: pH buffer, carbon source for some phytoplankton
    ! Units: mol/m³
    ! =====================================================================
    if (use_REcoM) then
        call def_stream((/nl-1, nod2D/), (/nl-1, myDim_nod2D/), &
                       'HCO3', 'Bicarbonate ion concentration', 'mol/m3', HCO33D(:,:),  &
                       io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    end if

CASE ('CO3       ')
    ! =====================================================================
    ! Variable: CO3
    ! Description: Carbonate ion concentration
    ! Function: Secondary dissolved inorganic carbon species
    ! Role: Building block for calcium carbonate (CaCO3) formation
    ! Units: mol/m³
    ! =====================================================================
    if (use_REcoM) then
        call def_stream((/nl-1, nod2D/), (/nl-1, myDim_nod2D/), &
                       'CO3', 'Carbonate ion concentration', 'mol/m3', CO33D(:,:), &
                       io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    end if

CASE ('OmegaC    ')
    ! =====================================================================
    ! Variable: OmegaC
    ! Description: Calcite saturation state (Ω)
    ! Function: Ratio of [Ca²⁺][CO3²⁻] to calcite solubility product
    ! Role: Determines calcification vs dissolution; Ω>1 favors precipitation
    ! Context: Ocean acidification reduces Ω, threatening calcifying organisms
    ! Units: NN (dimensionless)
    ! =====================================================================
    if (use_REcoM) then
        call def_stream((/nl-1, nod2D/), (/nl-1, myDim_nod2D/), &
                       'OmegaC','calcite saturation state', 'NN', OmegaC3D(:,:), &
                       io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    end if

CASE ('kspc      ')
    ! =====================================================================
    ! Variable: kspc
    ! Description: Calcite solubility product
    ! Function: Equilibrium constant for CaCO3 dissolution
    ! Role: Temperature/pressure dependent; controls calcite stability
    ! Units: mol²/kg²
    ! =====================================================================
    if (use_REcoM) then
        call def_stream((/nl-1, nod2D/), (/nl-1, myDim_nod2D/), &
                       'kspc', 'calcite solubility product', 'mol^2/kg^2', kspc3D(:,:), &
                       io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    end if

CASE ('rhoSW     ')
    ! =====================================================================
    ! Variable: rhoSW
    ! Description: In-situ seawater density
    ! Function: Density at ambient temperature, salinity, and pressure
    ! Role: Physical property affecting stratification and mixing
    ! Units: mol/m³ (CHECK: unusual units, typically kg/m³)
    ! =====================================================================
    if (use_REcoM) then
        call def_stream((/nl-1, nod2D/),  (/nl-1, myDim_nod2D/), &
                       'rhoSW', 'in-situ density of seawater', 'mol/m3', rhoSW3D(:,:), &
                       io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    end if

! =====================================================================
! PARTICLE DYNAMICS - DENSITY AND SINKING
! Description: Physical properties of detrital particles
! =====================================================================

CASE ('rho_det1       ')
    ! =====================================================================
    ! Variable: rho_det1
    ! Description: Density of detrital particles in size class 1
    ! Function: Mass per unit volume of small/slow-sinking particles
    ! Role: Determines sinking speed via Stokes' law
    ! Context: Lower density = slower sinking = longer remineralization time
    ! Units: kg/m³
    ! =====================================================================
    call def_stream((/nl-1, nod2D/), (/nl-1, myDim_nod2D/), &
                   'rho_det1', 'rho of particles in class 1', 'kg/m3', rho_particle1(:,:), & 
                   io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)

CASE ('rho_det2       ')
    ! =====================================================================
    ! Variable: rho_det2
    ! Description: Density of detrital particles in size class 2
    ! Function: Mass per unit volume of large/fast-sinking particles
    ! Role: Determines sinking speed via Stokes' law
    ! Context: Higher density = faster sinking = efficient carbon export
    ! Units: kg/m³
    ! =====================================================================
    call def_stream((/nl-1, nod2D/), (/nl-1, myDim_nod2D/), &
                   'rho_det2', 'rho of particles in class 2', 'kg/m3', rho_particle2(:,:), &
                   io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)

CASE ('scaling_visc   ')
    ! =====================================================================
    ! Variable: scaling_visc
    ! Description: Viscosity-based scaling factor for sinking speed
    ! Function: Adjusts particle sinking based on water viscosity
    ! Role: Temperature/salinity affect viscosity, thus sinking rates
    ! Context: Warmer water = lower viscosity = faster sinking
    ! Units: n.d. (non-dimensional scaling factor)
    ! =====================================================================
    call def_stream((/nl-1, nod2D/), (/nl-1, myDim_nod2D/), &
                   'scaling_visc', 'scaling factor of particle sinking speed', 'n.d.', scaling_visc_3D(:,:), &
                   io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)

CASE ('wsink_det1')
    ! =====================================================================
    ! Variable: wsink_det1
    ! Description: Sinking velocity of particle class 1
    ! Function: Vertical descent rate of small detrital particles
    ! Role: Controls residence time in euphotic zone, carbon export efficiency
    ! Context: Typically 1-50 m/day for small aggregates
    ! Units: m/s
    ! =====================================================================
   call def_stream((/nl, nod2D/), (/nl, myDim_nod2D/), &
                  'wsink_det1', 'sinking speed of particles in class 1', 'm s-1', Sinkingvel1(:,:), &
                  io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)

CASE ('wsink_det2')
    ! =====================================================================
    ! Variable: wsink_det2
    ! Description: Sinking velocity of particle class 2
    ! Function: Vertical descent rate of large detrital particles
    ! Role: Fast export pathway for carbon to deep ocean
    ! Context: Typically 50-200+ m/day for large aggregates/fecal pellets
    ! Units: m/s
    ! =====================================================================
   call def_stream((/nl, nod2D/), (/nl, myDim_nod2D/), &
                  'wsink_det2', 'sinking speed of particles in class 2', 'm s-1', Sinkingvel2(:,:), &
                  io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)

! =====================================================================
! BENTHIC NUTRIENT FLUXES
! Description: Sediment-water interface nutrient release
! =====================================================================

CASE ('DIC_bf')
    ! =====================================================================
    ! Variable: DIC_bf
    ! Description: Dissolved Inorganic Carbon bottom flux
    ! Function: Release of CO2/HCO3/CO3 from sediment remineralization
    ! Role: Recycles organic carbon back to water column
    ! Context: Important in shallow seas, negligible in deep ocean
    ! Units: mmolC/m³
    ! =====================================================================
    call def_stream((/nl, nod2D/), (/nl, myDim_nod2D/), &
                   'DIC_bf', 'DIC bottom flux', 'mmolC/m3', dtr_bflux_dic(:,:), &
                   io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)

CASE ('DIN_bf')
    ! =====================================================================
    ! Variable: DIN_bf
    ! Description: Dissolved Inorganic Nitrogen bottom flux
    ! Function: Release of NH4⁺/NO3⁻ from sediment remineralization
    ! Role: Nutrient regeneration supporting primary production
    ! Context: Major N source in coastal/shelf systems
    ! Units: mmolN/m³
    ! =====================================================================
    call def_stream((/nl, nod2D/), (/nl, myDim_nod2D/), &
                   'DIN_bf', 'DIN bottom flux', 'mmolC/m3', dtr_bflux_din(:,:), &
                   io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)

CASE ('Alk_bf')
    ! =====================================================================
    ! Variable: Alk_bf
    ! Description: Alkalinity bottom flux
    ! Function: Release of acid-neutralizing capacity from sediments
    ! Role: Affects pH and carbonate chemistry of bottom waters
    ! Context: Dominated by NH4⁺ production and CaCO3 dissolution
    ! Units: mmolC/m³ (CHECK: should be meq/m³)
    ! =====================================================================
    call def_stream((/nl, nod2D/), (/nl, myDim_nod2D/), &
                   'Alk_bf', 'Alk bottom flux', 'mmolC/m3', dtr_bflux_alk(:,:), &
                   io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)

CASE ('DSi_bf')
    ! =====================================================================
    ! Variable: DSi_bf
    ! Description: Dissolved Silicate bottom flux
    ! Function: Release of Si(OH)4 from opal dissolution in sediments
    ! Role: Critical nutrient for diatom growth
    ! Context: Can control diatom vs non-diatom dominance
    ! Units: mmolSi/m³
    ! =====================================================================
    call def_stream((/nl, nod2D/), (/nl, myDim_nod2D/), &
                   'DSi_bf', 'DSi bottom flux', 'mmolC/m3',  dtr_bflux_dsi(:,:), &
                   io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)

! =====================================================================
! ZOOPLANKTON RESPIRATION
! Description: Metabolic CO2 production by different zooplankton size classes
! =====================================================================

CASE ('respmeso       ')
    ! =====================================================================
    ! Variable: respmeso
    ! Description: Respiration rate of mesozooplankton
    ! Function: Metabolic CO2 release by medium-sized grazers (0.2-2 mm)
    ! Role: Returns carbon to DIC pool, drives oxygen consumption
    ! Taxa: Copepods, euphausiids, chaetognaths
    ! Units: mmolC/m²/d
    ! =====================================================================
    if (use_REcoM) then
        call def_stream((/nl-1, nod2D/), (/nl-1, myDim_nod2D/), &
                       'respmeso', 'Respiration rate of mesozooplankton', 'mmolC/m2/d', respmeso(:,:), &
                       io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    endif

CASE ('respmacro      ')
    ! =====================================================================
    ! Variable: respmacro
    ! Description: Respiration rate of macrozooplankton
    ! Function: Metabolic CO2 release by large grazers (>2 mm)
    ! Role: Produces fast-sinking fecal pellets, vertical carbon transport
    ! Taxa: Large euphausiids, amphipods, jellyfish
    ! Units: mmolC/m²/d
    ! =====================================================================
    if (use_REcoM) then
        call def_stream((/nl-1, nod2D/), (/nl-1, myDim_nod2D/), &
                       'respmacro', 'Respiration rate of macrozooplankton', 'mmolC/m2/d', respmacro(:,:), &
                       io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    endif

CASE ('respmicro      ')
    ! =====================================================================
    ! Variable: respmicro
    ! Description: Respiration rate of microzooplankton
    ! Function: Metabolic CO2 release by small grazers (<0.2 mm)
    ! Role: Rapid nutrient recycling in microbial loop
    ! Taxa: Ciliates, heterotrophic dinoflagellates
    ! Units: mmolC/m²/d
    ! =====================================================================
    if (use_REcoM) then
        call def_stream((/nl-1, nod2D/), (/nl-1, myDim_nod2D/), &
                       'respmicro', 'Respiration rate of microzooplankton', 'mmolC/m2/d', respmicro(:,:), &
                       io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    endif

! =====================================================================
! SECTION 5: CALCIUM CARBONATE DYNAMICS
! Description: Calcification and dissolution processes
! =====================================================================

CASE ('calcdiss       ')
    ! =====================================================================
    ! Variable: calcdiss
    ! Description: Calcite (CaCO3) dissolution rate
    ! Function: Conversion of solid CaCO3 back to dissolved Ca²⁺ and CO3²⁻
    ! Role: Increases alkalinity, releases CO2, occurs when Ω < 1
    ! Context: Accelerated by ocean acidification
    ! Units: mmolC/m²/d
    ! =====================================================================
    if (use_REcoM) then
        call def_stream((/nl-1, nod2D/), (/nl-1, myDim_nod2D/), &
                       'calcdiss', 'Calcite dissolution', 'mmolC/m2/d', calcdiss(:,:), &
                       io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    endif

CASE ('calcif         ')
    ! =====================================================================
    ! Variable: calcif
    ! Description: Calcification rate (CaCO3 formation)
    ! Function: Precipitation of calcium carbonate shells/tests
    ! Role: Decreases alkalinity, releases CO2, ballasts organic carbon
    ! Organisms: Coccolithophores, foraminifera, pteropods
    ! Units: mmolC/m²/d
    ! =====================================================================
    if (use_REcoM) then
        call def_stream((/nl-1, nod2D/), (/nl-1, myDim_nod2D/), &
                       'calcif', 'Calcification', 'mmolC/m2/d', calcif(:,:), &
                       io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    endif

! =====================================================================
! PHYTOPLANKTON AGGREGATION
! Description: Formation of larger particles from small cells
! =====================================================================

CASE ('aggn           ')
    ! =====================================================================
    ! Variable: aggn
    ! Description: Aggregation of small (nano) phytoplankton
    ! Function: Collision-driven clumping of 2-20 μm cells
    ! Role: Converts slow-sinking cells to faster-sinking aggregates
    ! Context: Enhanced by stickiness (TEP), turbulence, high cell density
    ! Units: mmolC/m²/d
    ! =====================================================================
    if (use_REcoM) then
        call def_stream((/nl-1, nod2D/), (/nl-1, myDim_nod2D/), &
                       'aggn', 'Aggregation of small phytoplankton', 'mmolC/m2/d', aggn(:,:), &
                       io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    endif

CASE ('aggd           ')
    ! =====================================================================
    ! Variable: aggd
    ! Description: Aggregation of diatoms
    ! Function: Formation of diatom-dominated marine snow
    ! Role: Major pathway for carbon export during bloom collapse
    ! Context: Si:C ballasting enhances sinking speed (100-200 m/day)
    ! Units: mmolC/m²/d
    ! =====================================================================
    if (use_REcoM) then
        call def_stream((/nl-1, nod2D/), (/nl-1, myDim_nod2D/), &
                       'aggd', 'Aggregation of diatoms', 'mmolC/m2/d', aggd(:,:), &
                       io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    endif

CASE ('aggc           ')
    ! =====================================================================
    ! Variable: aggc
    ! Description: Aggregation of coccolithophores
    ! Function: Formation of coccolith-bearing aggregates
    ! Role: CaCO3 ballasting enhances export efficiency
    ! Context: Creates "white waters" visible from space during blooms
    ! Units: mmolC/m²/d
    ! =====================================================================
    if (use_REcoM) then
        call def_stream((/nl-1, nod2D/), (/nl-1, myDim_nod2D/), &
                       'aggc', 'Aggregation of coccolithophores', 'mmolC/m2/d', aggc(:,:), &
                       io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    endif

CASE ('aggp           ')
    ! =====================================================================
    ! Variable: aggp
    ! Description: Aggregation of Phaeocystis
    ! Function: Colony formation and aggregation of colonial haptophyte
    ! Role: Forms mucilaginous colonies up to 3 cm, efficient export
    ! Context: Dominant in polar/subpolar regions, nuisance blooms
    ! Units: mmolC/(m²*d)
    ! =====================================================================
    if (use_REcoM) then
        call def_stream((/nl-1, nod2D/), (/nl-1, myDim_nod2D/), &
                       'aggp', 'Aggregation of phaeocystis', 'mmolC/(m2*d)', aggp(:,:), &
                       io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    endif

! =====================================================================
! DISSOLVED ORGANIC CARBON (DOC) EXCRETION
! Description: Release of dissolved organics by phytoplankton
! =====================================================================

CASE ('docexn         ')
    ! =====================================================================
    ! Variable: docexn
    ! Description: DOC excretion by nanophytoplankton
    ! Function: Passive/active release of dissolved organic compounds
    ! Role: Fuels microbial loop, can be 5-20% of gross production
    ! Context: Increases under nutrient stress or viral lysis
    ! Units: mmolC/m²/d
    ! =====================================================================
    if (use_REcoM) then
        call def_stream((/nl-1, nod2D/), (/nl-1, myDim_nod2D/), &
                       'docexn', 'DOC excretion by small phytoplankton', 'mmolC/m2/d', docexn(:,:), &
                       io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    endif

CASE ('docexd         ')
    ! =====================================================================
    ! Variable: docexd
    ! Description: DOC excretion by diatoms
    ! Function: Release of organic compounds from diatom metabolism
    ! Role: Can be substantial during stationary phase or Si limitation
    ! Context: Includes transparent exopolymer particles (TEP)
    ! Units: mmolC/m²/d
    ! =====================================================================
    if (use_REcoM) then
        call def_stream((/nl-1, nod2D/), (/nl-1, myDim_nod2D/), &
                       'docexd', 'DOC excretion by diatoms', 'mmolC/m2/d', docexd(:,:), &
                       io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    endif

CASE ('docexc         ')
    ! =====================================================================
    ! Variable: docexc
    ! Description: DOC excretion by coccolithophores
    ! Function: Release of organic compounds during calcification/growth
    ! Role: Links organic and inorganic carbon cycles
    ! Context: May increase under high pCO2 (ocean acidification)
    ! Units: mmolC/m²/d
    ! =====================================================================
    if (use_REcoM) then
        call def_stream((/nl-1, nod2D/), (/nl-1, myDim_nod2D/), &
                       'docexc', 'DOC excretion by coccolithophores', 'mmolC/m2/d', docexc(:,:), &
                       io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    endif

CASE ('docexp         ')
    ! =====================================================================
    ! Variable: docexp
    ! Description: DOC excretion by Phaeocystis
    ! Function: Release of polysaccharides forming colonial matrix
    ! Role: Major mucus producer, can create foam on beaches
    ! Context: High C:N ratio DOC, relatively refractory
    ! Units: mmolC/(m²*d)
    ! =====================================================================
   if (use_REcoM) then
       call def_stream((/nl-1, nod2D/), (/nl-1, myDim_nod2D/), &
                      'docexp', 'DOC excretion by phaeocystis', 'mmolC/(m2*d)', docexp(:,:), &
                      io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
   endif

! =====================================================================
! PHYTOPLANKTON RESPIRATION
! Description: Metabolic CO2 release by autotrophs
! =====================================================================

CASE ('respn          ')
    ! =====================================================================
    ! Variable: respn
    ! Description: Respiration by nanophytoplankton
    ! Function: Maintenance metabolism and biosynthesis costs
    ! Role: Determines net vs gross primary production (NPP vs GPP)
    ! Context: Typically 10-30% of GPP, increases with temperature
    ! Units: mmolC/m²/d
    ! =====================================================================
   if (use_REcoM) then
       call def_stream((/nl-1, nod2D/), (/nl-1, myDim_nod2D/), &
                      'respn', 'Respiration by small phytoplankton', 'mmolC/m2/d', respn(:,:), &
                      io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
   endif

CASE ('respd          ')
    ! =====================================================================
    ! Variable: respd
    ! Description: Respiration by diatoms
    ! Function: Metabolic costs of rapid growth and Si metabolism
    ! Role: Higher in fast-growing bloom-forming species
    ! Context: Can increase under nutrient limitation
    ! Units: mmolC/m²/d
    ! =====================================================================
    if (use_REcoM) then
         call def_stream((/nl-1, nod2D/), (/nl-1, myDim_nod2D/), &
                        'respd', 'Respiration by diatoms', 'mmolC/m2/d', respd(:,:), & 
                        io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    endif

CASE ('respc          ')
    ! =====================================================================
    ! Variable: respc
    ! Description: Respiration by coccolithophores
    ! Function: Metabolic costs including calcification energy
    ! Role: Calcification requires ~20% additional energy vs non-calcifiers
    ! Context: May increase under ocean acidification stress
    ! Units: mmolC/(m²*d)
    ! =====================================================================
    if (use_REcoM) then
        call def_stream((/nl-1, nod2D/), (/nl-1, myDim_nod2D/), &
                       'respc', 'Respiration by coccolithophores', 'mmolC/(m2*d)', respc(:,:), &
                       io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    endif

CASE ('respp          ')
    ! =====================================================================
    ! Variable: respp
    ! Description: Respiration by Phaeocystis
    ! Function: Metabolic costs of colonial vs solitary forms
    ! Role: Colonial forms may have lower specific respiration rates
    ! Context: Life cycle shifts affect carbon cycling
    ! Units: mmolC/(m²*d)
    ! =====================================================================
    if (use_REcoM) then
        call def_stream((/nl-1, nod2D/), (/nl-1, myDim_nod2D/), &
                       'respp', 'Respiration by phaeocystis', 'mmolC/(m2*d)', respp(:,:), &
                       io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    endif

! =====================================================================
! NET PRIMARY PRODUCTION (NPP) BY FUNCTIONAL GROUP
! Description: Net carbon fixation after respiration losses
! =====================================================================

CASE ('NPPn3D         ')
    ! =====================================================================
    ! Variable: NPPn3D
    ! Description: Net Primary Production of nanophytoplankton
    ! Function: Gross photosynthesis minus respiration and exudation
    ! Role: Sustained baseline production in oligotrophic waters
    ! Ecology: Dominant in stratified, nutrient-poor conditions
    ! Units: mmolC/m²/d
    ! =====================================================================
    if (use_REcoM) then
    call def_stream((/nl-1, nod2D/), (/nl-1, myDim_nod2D/), &
                   'NPPn3D', 'Net primary production of small phytoplankton', 'mmolC/m2/d', NPPn3D(:,:), &
                   io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    endif

CASE ('NPPd3D         ')
    ! =====================================================================
    ! Variable: NPPd3D
    ! Description: Net Primary Production of diatoms
    ! Function: Net carbon fixation by silicified phytoplankton
    ! Role: Bloom-forming, episodic production pulses
    ! Ecology: Dominant in upwelling, mixed, nutrient-rich waters
    ! Units: mmolC/m²/d
    ! =====================================================================
    if (use_REcoM) then
        call def_stream((/nl-1, nod2D/), (/nl-1, myDim_nod2D/), &
                       'NPPd3D', 'Net primary production of diatoms', 'mmolC/m2/d', NPPd3D(:,:), &
                       io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    endif

CASE ('NPPc3D         ')
    ! =====================================================================
    ! Variable: NPPc3D
    ! Description: Net Primary Production of coccolithophores
    ! Function: Net carbon fixation by calcifying haptophytes
    ! Role: Links organic and inorganic carbon pumps
    ! Ecology: Bloom in stratified, post-bloom nutrient-depleted waters
    ! Units: mmolC/m²/d
    ! =====================================================================
    if (use_REcoM) then
        call def_stream((/nl-1, nod2D/), (/nl-1, myDim_nod2D/), &
                       'NPPc3D', 'Net primary production of coccolithophores', 'mmolC/m2/d', NPPc3D(:,:), &
                       io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    endif

CASE ('NPPp3D         ')
    ! =====================================================================
    ! Variable: NPPp3D
    ! Description: Net Primary Production of Phaeocystis
    ! Function: Net carbon fixation by colonial haptophyte
    ! Role: High-latitude spring bloom specialist
    ! Ecology: Tolerates low light, cold temperatures, ice-edge blooms
    ! Units: mmolC/(m²*d)
    ! =====================================================================
   if (use_REcoM) then
       call def_stream((/nl-1, nod2D/), (/nl-1, myDim_nod2D/), &
                      'NPPp3D', 'Net primary production of phaeocystis', 'mmolC/(m2*d)', NPPp3D(:,:), &
                      io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
   endif

! =====================================================================
! TEMPERATURE EFFECTS ON PHOTOSYNTHESIS
! Description: Temperature-dependent growth rate modifiers
! =====================================================================

CASE ('TTemp_diatoms          ')
    ! =====================================================================
    ! Variable: TTemp_diatoms
    ! Description: Temperature dependence factor for diatom photosynthesis
    ! Function: Q10-based or Arrhenius growth rate modifier
    ! Role: Accelerates growth in warmer waters
    ! Context: Diatoms often have lower temperature optimum than flagellates
    ! Units: per day (dimensionless multiplier)
    ! =====================================================================
    if (use_REcoM) then
        call def_stream((/nl-1, nod2D/), (/nl-1, myDim_nod2D/), &
                       'TTemp_diatoms', 'Temperature dependence of diatoms PS', 'per day', TTemp_diatoms(:,:), &
                       io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    endif

CASE ('TTemp_phyto          ')
    ! =====================================================================
    ! Variable: TTemp_phyto
    ! Description: Temperature dependence factor for nanophytoplankton
    ! Function: Growth rate scaling with temperature
    ! Role: Often have broader temperature tolerance than specialists
    ! Context: May benefit from warming in temperate/tropical regions
    ! Units: per day (dimensionless multiplier)
    ! =====================================================================
    if (use_REcoM) then 
        call def_stream((/nl-1, nod2D/), (/nl-1, myDim_nod2D/), &
                       'TTemp_phyto', 'Temperature dependence of small phytoplankton PS', 'per day', TTemp_phyto(:,:), &
                       io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    endif

CASE ('TTemp_cocco          ')
    ! =====================================================================
    ! Variable: TTemp_cocco
    ! Description: Temperature dependence factor for coccolithophores
    ! Function: Growth rate and calcification temperature sensitivity
    ! Role: Often have higher temperature optima than diatoms
    ! Context: May expand poleward with warming
    ! Units: per day (dimensionless multiplier)
    ! =====================================================================
    if (use_REcoM) then
        call def_stream((/nl-1, nod2D/), (/nl-1, myDim_nod2D/), &
                       'TTemp_cocco', 'Temperature dependence of coccolithophores PS', 'per day', TTemp_cocco(:,:), &
                       io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    endif

CASE ('TTemp_phaeo          ')
    ! =====================================================================
    ! Variable: TTemp_phaeo
    ! Description: Temperature dependence factor for Phaeocystis
    ! Function: Growth rate scaling in cold polar waters
    ! Role: Adapted to low temperatures (0-10°C optimum)
    ! Context: Climate warming may reduce competitive advantage
    ! Units: per day (dimensionless multiplier)
    ! =====================================================================
    if (use_REcoM) then
        call def_stream((/nl-1, nod2D/), (/nl-1, myDim_nod2D/), & 
                       'TTemp_phaeo', 'Temperature dependence of phaeocystis PS', 'per day', TTemp_phaeo(:,:), &
                       io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    endif

! =====================================================================
! SECTION 11: CO2 EFFECTS ON PHOTOSYNTHESIS
! Description: pCO2-dependent growth modifiers (ocean acidification impacts)
! =====================================================================

CASE ('TPhyCO2          ')
    ! =====================================================================
    ! Variable: TPhyCO2
    ! Description: Effect of CO2 on nanophytoplankton growth
    ! Function: Carbon concentration mechanism (CCM) efficiency modifier
    ! Role: Some species benefit from higher CO2 (reduced CCM energy cost)
    ! Context: Generally small positive or neutral effect
    ! Units: per day (dimensionless multiplier)
    ! =====================================================================
    if (use_REcoM) then
        call def_stream((/nl-1, nod2D/), (/nl-1, myDim_nod2D/), &
                       'TPhyCO2', 'Effect of CO2 of phyto growth/PS', 'per day', TPhyCO2(:,:), &
                       io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    endif

CASE ('TDiaCO2          ')
    ! =====================================================================
    ! Variable: TDiaCO2
    ! Description: Effect of CO2 on diatom growth
    ! Function: CCM efficiency and photorespiration modifier
    ! Role: Mixed responses, some species C3-like, others have efficient CCMs
    ! Context: High CO2 may reduce competitive advantage vs other groups
    ! Units: per day (dimensionless multiplier)
    ! =====================================================================
    if (use_REcoM) then
        call def_stream((/nl-1, nod2D/), (/nl-1, myDim_nod2D/), &
                       'TDiaCO2', 'Effect of CO2 of phyto growth/PS', 'per day', TDiaCO2(:,:), &
                       io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    endif

CASE ('TCoccoCO2          ')
    ! =====================================================================
    ! Variable: TCoccoCO2
    ! Description: Effect of CO2 on coccolithophore growth
    ! Function: Photosynthesis vs calcification trade-off under high CO2
    ! Role: High pCO2 may benefit photosynthesis but impair calcification
    ! Context: Net effect species-specific, generally negative at >1000 μatm
    ! Units: per day (dimensionless multiplier)
    ! =====================================================================
    if (use_REcoM) then
        call def_stream((/nl-1, nod2D/), (/nl-1, myDim_nod2D/), &
                       'TCoccoCO2', 'Effect of CO2 of phyto growth/PS', 'per day', TCoccoCO2(:,:), &
                       io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    endif

CASE ('TPhaeoCO2          ')
    ! =====================================================================
    ! Variable: TPhaeoCO2
    ! Description: Effect of CO2 on Phaeocystis growth
    ! Function: CCM efficiency in high-latitude species
    ! Role: May benefit from CO2 enrichment in cold waters
    ! Context: Limited experimental data compared to other groups
    ! Units: per day (dimensionless multiplier)
    ! =====================================================================
    if (use_REcoM) then
        call def_stream((/nl-1, nod2D/), (/nl-1, myDim_nod2D/), &
                       'TPhaeoCO2', 'Effect of CO2 of phyto growth/PS', 'per day', TPhaeoCO2(:,:), &
                       io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    endif

! =====================================================================
! SECTION 12: NUTRIENT LIMITATION FACTORS
! Description: Multi-nutrient limitation of photosynthesis
! =====================================================================

CASE ('TqLF_phyto          ')
    ! =====================================================================
    ! Variable: TqLF_phyto
    ! Description: Nutrient limitation factor for nanophytoplankton
    ! Function: Minimum of N, P, Fe limitation (Liebig's Law)
    ! Role: Determines realized vs maximum potential growth rate
    ! Context: Often N-limited in stratified waters, Fe-limited in HNLC
    ! Units: per day (0-1 dimensionless multiplier)
    ! =====================================================================
    if (use_REcoM) then
        call def_stream((/nl-1, nod2D/), (/nl-1, myDim_nod2D/), &
                       'TqLF_phyto', 'Nutrient limitation effect of phyto PS', 'per day', TqlimitFac_phyto(:,:), &
                       io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    endif

CASE ('TqLF_diatoms          ')
    ! =====================================================================
    ! Variable: TqLF_diatoms
    ! Description: Nutrient limitation factor for diatoms
    ! Function: Minimum of N, P, Si, Fe limitation
    ! Role: Si limitation unique to diatoms, controls bloom termination
    ! Context: Si:N ratio critical for diatom vs flagellate competition
    ! Units: per day (0-1 dimensionless multiplier)
    ! =====================================================================
    if (use_REcoM) then
    call def_stream((/nl-1, nod2D/), (/nl-1, myDim_nod2D/), &
                   'TqLF_diatoms', 'Nutrient limitation effect of diatoms PS', 'per day', TqlimitFac_diatoms(:,:), &
                   io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    endif

CASE ('TqLF_cocco          ')
    ! =====================================================================
    ! Variable: TqLF_cocco
    ! Description: Nutrient limitation factor for coccolithophores
    ! Function: Minimum of N, P limitation (no Si requirement)
    ! Role: Can exploit low-nutrient conditions after diatom bloom collapse
    ! Context: High C:N, C:P ratios during severe limitation
    ! Units: per day (0-1 dimensionless multiplier)
    ! =====================================================================
    if (use_REcoM) then
    call def_stream((/nl-1, nod2D/),  (/nl-1, myDim_nod2D/),   'TqLF_cocco','Nutrient limitation effect of diatoms PS', 'per day',TqlimitFac_cocco(:,:), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    endif

CASE ('TqLF_phaeo          ')
    ! =====================================================================
    ! Variable: TqLF_phaeo
    ! Description: Nutrient limitation factor for Phaeocystis
    ! Function: Minimum of N, P, Fe limitation
    ! Role: Can form luxury nutrient reserves in colonial form
    ! Context: Often blooms early when nutrients still replete
    ! Units: per day (0-1 dimensionless multiplier)
    ! =====================================================================
    if (use_REcoM) then
        call def_stream((/nl-1, nod2D/), (/nl-1, myDim_nod2D/), &
                       'TqLF_phaeo', 'Nutrient limitation effect of diatoms PS', 'per day', TqlimitFac_phaeo(:,:), &
                       io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    endif

! =====================================================================
! LIGHT LIMITATION FACTORS
! Description: Photosynthesis-irradiance relationship modifiers
! =====================================================================

CASE ('TCphotLL_phyto         ')
    ! =====================================================================
    ! Variable: TCphotLL_phyto
    ! Description: Light limitation factor for nanophytoplankton
    ! Function: P-I curve response (typically hyperbolic tangent or exponential)
    ! Role: Reduces growth in low light, photoinhibition at high light
    ! Context: Small cells acclimate faster to changing light
    ! Units: per day (0-1 dimensionless multiplier)
    ! =====================================================================
    if (use_REcoM) then
        call def_stream((/nl-1, nod2D/), (/nl-1, myDim_nod2D/), &
                       'TCphotLL_phyto', 'Light limitation on phyto PS', 'per day', TCphotLigLim_phyto(:,:), &
                       io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    endif

CASE ('TCphotLL_dia         ')
    ! =====================================================================
    ! Variable: TCphotLL_dia
    ! Description: Light limitation factor for diatoms
    ! Function: P-I curve with potentially higher light saturation point
    ! Role: Often shade-adapted in coastal turbid waters
    ! Context: Larger cells = slower photoacclimation
    ! Units: per day (0-1 dimensionless multiplier)
    ! =====================================================================
    if (use_REcoM) then
        call def_stream((/nl-1, nod2D/),  (/nl-1, myDim_nod2D/), &
                       'TCphotLL_dia', 'Light limitation on diatoms PS', 'per day', TCphotLigLim_diatoms(:,:), &
                       io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    endif

CASE ('TCphotLL_cocco          ')
    ! =====================================================================
    ! Variable: TCphotLL_cocco
    ! Description: Light limitation factor for coccolithophores
    ! Function: P-I curve for high-light adapted species
    ! Role: Often prefer high irradiance, stratified waters
    ! Context: Bloom in clear waters with deep light penetration
    ! Units: per day (0-1 dimensionless multiplier)
    ! =====================================================================
    if (use_REcoM) then
        call def_stream((/nl-1, nod2D/), (/nl-1, myDim_nod2D/), &
                       'TCphotLL_cocco', 'Light limitation on phyto PS', 'per day', TCphotLigLim_cocco(:,:), &
                       io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    endif

CASE ('TCphotLL_phaeo          ')
    ! =====================================================================
    ! Variable: TCphotLL_phaeo
    ! Description: Light limitation factor for Phaeocystis
    ! Function: P-I curve for low-light adapted polar species
    ! Role: Low light saturation point, efficient in weak polar light
    ! Context: Can grow under sea ice with <1% surface irradiance
    ! Units: per day (0-1 dimensionless multiplier)
    ! =====================================================================
    if (use_REcoM) then
        call def_stream((/nl-1, nod2D/), (/nl-1, myDim_nod2D/), &
                       'TCphotLL_phaeo', 'Light limitation on phyto PS', 'per day', TCphotLigLim_phaeo(:,:), &
                       io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    endif


! =====================================================================
! CARBON-SPECIFIC PHOTOSYNTHESIS RATES
! Description: Actual photosynthetic carbon fixation rates
! =====================================================================

CASE ('TCphot_phyto          ')
    ! =====================================================================
    ! Variable: TCphot_phyto
    ! Description: Carbon-specific photosynthesis rate of nanophytoplankton
    ! Function: Gross carbon fixation per unit biomass
    ! Role: Product of all limiting factors × maximum rate
    ! Context: Typically 0.5-3.0 d⁻¹ in natural populations
    ! Units: per day (d⁻¹)
    ! =====================================================================
    if (use_REcoM) then
        call def_stream((/nl-1, nod2D/), (/nl-1, myDim_nod2D/), &
                       'TCphot_phyto', 'Light limitation on phyto PS', 'per day', TCphot_phyto(:,:), &
                       io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    endif

CASE ('TCphot_diatoms          ')
    ! =====================================================================
    ! Variable: TCphot_diatoms
    ! Description: Carbon-specific photosynthesis rate of diatoms
    ! Function: Gross carbon fixation per unit biomass
    ! Role: Can reach very high rates (>5 d⁻¹) during blooms
    ! Context: r-strategists with high maximum growth rates
    ! Units: per day (d⁻¹)
    ! =====================================================================
    if (use_REcoM) then
        call def_stream((/nl-1, nod2D/), (/nl-1, myDim_nod2D/), &
                       'TCphot_diatoms', 'Light limitation on phyto PS', 'per day', TCphot_diatoms(:,:), &
                       io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    endif   

CASE ('TCphot_cocco          ')
    ! =====================================================================
    ! Variable: TCphot_cocco
    ! Description: Carbon-specific photosynthesis rate of coccolithophores
    ! Function: Gross carbon fixation including calcification costs
    ! Role: Lower than diatoms but sustained in stratified conditions
    ! Context: Additional energy for CaCO3 formation reduces net efficiency
    ! Units: per day (d⁻¹)
    ! =====================================================================
    if (use_REcoM) then
        call def_stream((/nl-1, nod2D/), (/nl-1, myDim_nod2D/), &
                       'TCphot_cocco', 'Light limitation on phyto PS', 'per day', TCphot_cocco(:,:), &
                       io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    endif

CASE ('TCphot_phaeo          ')
    ! =====================================================================
    ! Variable: TCphot_phaeo
    ! Description: Carbon-specific photosynthesis rate of Phaeocystis
    ! Function: Gross carbon fixation in colonial or solitary forms
    ! Role: Efficient at low temperature and light
    ! Context: Colonial forms may have lower per-cell rates but higher bloom biomass
    ! Units: per day (d⁻¹)
    ! =====================================================================
    if (use_REcoM) then
        call def_stream((/nl-1, nod2D/), (/nl-1, myDim_nod2D/), &
                       'TCphot_phaeo', 'Light limitation on phyto PS', 'per day', TCphot_phaeo(:,:), &
                       io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    endif

! =====================================================================
! SILICATE ASSIMILATION
! Description: Diatom-specific nutrient uptake
! =====================================================================

CASE ('TSi_assimDia          ')
    ! =====================================================================
    ! Variable: TSi_assimDia
    ! Description: Silicate assimilation rate by diatoms
    ! Function: Uptake of Si(OH)4 for frustule (shell) construction
    ! Role: Couples diatom growth to silicon cycle
    ! Context: Si:C ratio varies 0.1-1.0 depending on species and conditions
    ! Ecology: Si depletion terminates diatom blooms, shifts to flagellates
    ! Units: per day (d⁻¹)
    ! =====================================================================
    if (use_REcoM) then
        call def_stream((/nl-1, nod2D/), (/nl-1, myDim_nod2D/), &
                       'TSi_assimDia', 'Silicate assimilation', 'per day', TSi_assimDia(:,:), &
                       io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    endif

#endif

CASE ('otracers  ')
    do j=3, tracers%num_tracers
    write (id_string, "(I4.4)") tracers%data(j)%ID
#if defined(__recom)
      ! =====================================================================
      ! Base tracers (always present in all configurations)
      ! Tracer IDs: 1001-1022 (indices 3-24)
      ! =====================================================================
      if (tracers%data(j)%ID==1001) then
         if (use_REcoM) then
         call def_stream((/nl-1, nod2D/),  (/nl-1, myDim_nod2D/),  'DIN', 'Dissolved Inorganic Nitrogen', '[mmol/m3]', tracers%data(j)%values(:,:), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
         endif



            if (tracers%data(j)%ltra_diag) then
               ! horizontal advection
               call def_stream((/nl-1, nod2D/),  (/nl-1, myDim_nod2D/),  'DIN_hor_adv', 'Horizontal advection part of dissolved Inorganic N', '[mmol/m3]', tracers%work%tra_advhoriz(:,:,j), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
                ! horizontal advection LO
               call def_stream((/nl-1, nod2D/),  (/nl-1, myDim_nod2D/),  'DIN_lo_hor_adv', 'LO Horizontal advection part of dissolved Inorganic N', '[mmol/m3]', tracers%work%tra_advhoriz_LO(:,:,j), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)

               ! vertical advection
               call def_stream((/nl-1, nod2D/),  (/nl-1, myDim_nod2D/),  'DIN_ver_adv', 'Vertical advection part of dissolved Inorganic N', '[mmol/m3]', tracers%work%tra_advvert(:,:,j), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)

               ! vertical advection LO
               call def_stream((/nl-1, nod2D/),  (/nl-1, myDim_nod2D/),  'DIN_lo_ver_adv', 'LO Vertical advection part of dissolved Inorganic N', '[mmol/m3]', tracers%work%tra_advvert_LO(:,:,j), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)

               ! horizontal diffusion
               call def_stream((/nl-1, nod2D/),  (/nl-1, myDim_nod2D/),  'DIN_tra_diff_part_hor_redi', 'Horizontal diffusion of dissolved Inorganic N (includes Redi diffusivity if Redi=.true.)', '[mmol/m3]', tracers%work%tra_diff_part_hor_redi(:,:,j), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)

               if (tracers%data(j)%i_vert_diff) then
               ! vertical diffusion (Explicit)
                   call def_stream((/nl-1, nod2D/),  (/nl-1, myDim_nod2D/),  'DIN_tra_diff_part_ver_expl', 'Vertical diffusion of dissolved Inorganic N (Explicit)', '[mmol/m3]', tracers%work%tra_diff_part_ver_expl(:,:,j), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
               end if

               ! projection of horizontal Redi diffussivity onto vertical
               call def_stream((/nl-1, nod2D/),  (/nl-1, myDim_nod2D/),  'DIN_tra_diff_part_ver_redi_expl', 'Projection of horizontal Redi diffussivity onto vertical for dissolved Inorganic N (Explicit)', '[mmol/m3]', tracers%work%tra_diff_part_ver_redi_expl(:,:,j), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)

               ! vertical diffusion (Implicit)
               call def_stream((/nl-1, nod2D/),  (/nl-1, myDim_nod2D/),  'DIN_tra_diff_part_ver_impl', 'Vertical diffusion of dissolved Inorganic N (Implicit)', '[mmol/m3]', tracers%work%tra_diff_part_ver_impl(:,:,j), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)

               ! recom_sms
               call def_stream((/nl-1, nod2D/),  (/nl-1, myDim_nod2D/),  'DIN_recom_sms', 'Recom SMS', '[mmol/m3]', tracers%work%tra_recom_sms(:,:,j), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
            end if








      else if (tracers%data(j)%ID==1002) then
         if (use_REcoM) then
         call def_stream((/nl-1, nod2D/),  (/nl-1, myDim_nod2D/),  'DIC', 'Dissolved Inorganic C', '[mmol/m3]', tracers%data(j)%values(:,:), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)

            if (tracers%data(j)%ltra_diag) then
               ! horizontal advection
               call def_stream((/nl-1, nod2D/),  (/nl-1, myDim_nod2D/),  'DIC_hor_adv', 'Horizontal advection part of dissolved Inorganic C', '[mmol/m3]', tracers%work%tra_advhoriz(:,:,j), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
                ! horizontal advection LO
               call def_stream((/nl-1, nod2D/),  (/nl-1, myDim_nod2D/),  'DIC_lo_hor_adv', 'LO Horizontal advection part of dissolved Inorganic C', '[mmol/m3]', tracers%work%tra_advhoriz_LO(:,:,j), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)

               ! vertical advection
               call def_stream((/nl-1, nod2D/),  (/nl-1, myDim_nod2D/),  'DIC_ver_adv', 'Vertical advection part of dissolved Inorganic C', '[mmol/m3]', tracers%work%tra_advvert(:,:,j), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)

               ! vertical advection LO
               call def_stream((/nl-1, nod2D/),  (/nl-1, myDim_nod2D/),  'DIC_lo_ver_adv', 'LO Vertical advection part of dissolved Inorganic C', '[mmol/m3]', tracers%work%tra_advvert_LO(:,:,j), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)

               ! horizontal diffusion
               call def_stream((/nl-1, nod2D/),  (/nl-1, myDim_nod2D/),  'DIC_tra_diff_part_hor_redi', 'Horizontal diffusion of dissolved Inorganic C (includes Redi diffusivity if Redi=.true.)', '[mmol/m3]', tracers%work%tra_diff_part_hor_redi(:,:,j), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)

               if (tracers%data(j)%i_vert_diff) then
               ! vertical diffusion (Explicit)
                   call def_stream((/nl-1, nod2D/),  (/nl-1, myDim_nod2D/),  'DIC_tra_diff_part_ver_expl', 'Vertical diffusion of dissolved Inorganic C (Explicit)', '[mmol/m3]', tracers%work%tra_diff_part_ver_expl(:,:,j), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
               end if

               ! projection of horizontal Redi diffussivity onto vertical
               call def_stream((/nl-1, nod2D/),  (/nl-1, myDim_nod2D/),  'DIC_tra_diff_part_ver_redi_expl', 'Projection of horizontal Redi diffussivity onto vertical for dissolved Inorganic C (Explicit)', '[mmol/m3]', tracers%work%tra_diff_part_ver_redi_expl(:,:,j), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)

               ! vertical diffusion (Implicit)
               call def_stream((/nl-1, nod2D/),  (/nl-1, myDim_nod2D/),  'DIC_tra_diff_part_ver_impl', 'Vertical diffusion of dissolved Inorganic C (Implicit)', '[mmol/m3]', tracers%work%tra_diff_part_ver_impl(:,:,j), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)

               ! recom_sms
               call def_stream((/nl-1, nod2D/),  (/nl-1, myDim_nod2D/),  'DIC_recom_sms', 'Recom SMS', '[mmol/m3]', tracers%work%tra_recom_sms(:,:,j), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
            end if

         endif

      else if (tracers%data(j)%ID==1003) then
         if (use_REcoM) then
         call def_stream((/nl-1, nod2D/),  (/nl-1, myDim_nod2D/),  'Alk', 'Total Alkalinity', '[mmol/m3]', tracers%data(j)%values(:,:), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
!!!!24.06.2025
            if (tracers%data(j)%ltra_diag) then ! OG - tra_diag
               ! horizontal advection
               call def_stream((/nl-1, nod2D/),  (/nl-1, myDim_nod2D/),  'Alk_hor_adv', 'Horizontal advection part of Total Alkalinity', '[mmol/m3]', tracers%work%tra_advhoriz(:,:,j), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
                ! horizontal advection LO
               call def_stream((/nl-1, nod2D/),  (/nl-1, myDim_nod2D/),  'Alk_lo_hor_adv', 'LO Horizontal advection part of Total Alkalinity', '[mmol/m3]', tracers%work%tra_advhoriz_LO(:,:,j), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
               ! vertical advection
               call def_stream((/nl-1, nod2D/),  (/nl-1, myDim_nod2D/),  'Alk_ver_adv', 'Vertical advection part of Total Alkalinity', '[mmol/m3]', tracers%work%tra_advvert(:,:,j), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)

               ! horizontal diffusion
               call def_stream((/nl-1, nod2D/),  (/nl-1, myDim_nod2D/),  'Alk_tra_diff_part_hor_redi', 'Horizontal diffusion of Total Alkalinity (includes Redi diffusivity if Redi=.true.)', '[mmol/m3]', tracers%work%tra_diff_part_hor_redi(:,:,j), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)

               if (.not. tracers%data(j)%i_vert_diff) then
               ! vertical diffusion (Explicit)
                   call def_stream((/nl-1, nod2D/),  (/nl-1, myDim_nod2D/),  'Alk_tra_diff_part_ver_expl', 'Vertical diffusion of Total Alkalinity (Explicit)', '[mmol/m3]', tracers%work%tra_diff_part_ver_expl(:,:,j), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
               end if

               ! projection of horizontal Redi diffussivity onto vertical
               call def_stream((/nl-1, nod2D/),  (/nl-1, myDim_nod2D/),  'Alk_tra_diff_part_ver_redi_expl', 'Projection of horizontal Redi diffussivity onto vertical for Total Alkalinity (Explicit)', '[mmol/m3]', tracers%work%tra_diff_part_ver_redi_expl(:,:,j), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)

               ! vertical diffusion (Implicit)
               call def_stream((/nl-1, nod2D/),  (/nl-1, myDim_nod2D/),  'Alk_tra_diff_part_ver_impl', 'Vertical diffusion of Total Alkalinity (Implicit)', '[mmol/m3]', tracers%work%tra_diff_part_ver_impl(:,:,j), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)

               ! recom_sms
               call def_stream((/nl-1, nod2D/),  (/nl-1, myDim_nod2D/),  'Alk_recom_sms', 'Recom SMS', '[mmol/m3]', tracers%work%tra_recom_sms(:,:,j), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)

            end if
         endif

      ! =====================================================================
      ! Small phytoplankton
      ! =====================================================================
      else if (tracers%data(j)%ID==1004) then
         if (use_REcoM) then
         call def_stream((/nl-1, nod2D/),  (/nl-1, myDim_nod2D/),  'PhyN', 'Intracellular conc of Nitrogen in small phytoplankton', '[mmol/m3]', tracers%data(j)%values(:,:), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
         endif

      else if (tracers%data(j)%ID==1005) then
         if (use_REcoM) then
         call def_stream((/nl-1, nod2D/),  (/nl-1, myDim_nod2D/),  'PhyC', 'Intracellular conc of Carbon in small phytoplankton', '[mmol/m3]', tracers%data(j)%values(:,:), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
         endif

      else if (tracers%data(j)%ID==1006) then
         if (use_REcoM) then
         call def_stream((/nl-1, nod2D/),  (/nl-1, myDim_nod2D/),  'PhyChl', 'Current intracellular ChlA conc.', '[mg/m3]', tracers%data(j)%values(:,:), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
         endif

      ! =====================================================================
      ! Detritus 1
      ! =====================================================================
      else if (tracers%data(j)%ID==1007) then
         if (use_REcoM) then
         call def_stream((/nl-1, nod2D/),  (/nl-1, myDim_nod2D/),  'DetN', 'Conc of N in Detritus', '[mmol/m3]', tracers%data(j)%values(:,:), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
         endif

      else if (tracers%data(j)%ID==1008) then
         if (use_REcoM) then
         call def_stream((/nl-1, nod2D/),  (/nl-1, myDim_nod2D/),  'DetC', 'Conc of C in Detritus', '[mmol/m3]', tracers%data(j)%values(:,:), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
         endif

      ! =====================================================================
      ! Mesozooplankton
      ! =====================================================================
      else if (tracers%data(j)%ID==1009) then
         if (use_REcoM) then
         call def_stream((/nl-1, nod2D/),  (/nl-1, myDim_nod2D/),  'HetN', 'Conc of N in heterotrophs', '[mmol/m3]', tracers%data(j)%values(:,:), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
         endif

      else if (tracers%data(j)%ID==1010) then
         if (use_REcoM) then
         call def_stream((/nl-1, nod2D/),  (/nl-1, myDim_nod2D/),  'HetC', 'Conc of C in heterotrophs', '[mmol/m3]', tracers%data(j)%values(:,:), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
         endif

      ! =====================================================================
      ! Dissolved organic matter
      ! =====================================================================
      else if (tracers%data(j)%ID==1011) then
         if (use_REcoM) then
         call def_stream((/nl-1, nod2D/),  (/nl-1, myDim_nod2D/),  'DON', 'Dissolved organic N in the water', '[mmol/m3]', tracers%data(j)%values(:,:), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
         endif

      else if (tracers%data(j)%ID==1012) then
         if (use_REcoM) then
         call def_stream((/nl-1, nod2D/),  (/nl-1, myDim_nod2D/),  'DOC', 'Dissolved Organic C in the water', '[mmol/m3]', tracers%data(j)%values(:,:), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
         endif

      ! =====================================================================
      ! Diatoms
      ! =====================================================================

      else if (tracers%data(j)%ID==1013) then
         if (use_REcoM) then
         call def_stream((/nl-1, nod2D/),  (/nl-1, myDim_nod2D/),  'DiaN', 'DiaN', '[mmol/m3]', tracers%data(j)%values(:,:), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
         endif

      else if (tracers%data(j)%ID==1014) then
         if (use_REcoM) then
         call def_stream((/nl-1, nod2D/),  (/nl-1, myDim_nod2D/),  'DiaC', 'DiaC', '[mmol/m3]', tracers%data(j)%values(:,:), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
         endif

      else if (tracers%data(j)%ID==1015) then
         if (use_REcoM) then
         call def_stream((/nl-1, nod2D/),  (/nl-1, myDim_nod2D/),  'DiaChl', 'DiaChl', '[mg/m3]', tracers%data(j)%values(:,:), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
         endif

      else if (tracers%data(j)%ID==1016) then
         if (use_REcoM) then
         call def_stream((/nl-1, nod2D/),  (/nl-1, myDim_nod2D/),  'DiaSi', 'DiaSi', '[mmol/m3]', tracers%data(j)%values(:,:), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
         endif

      ! =====================================================================
      ! Detrital silica and other base tracers
      ! =====================================================================

      else if (tracers%data(j)%ID==1017) then
         if (use_REcoM) then
         call def_stream((/nl-1, nod2D/),  (/nl-1, myDim_nod2D/),  'DetSi', 'DetSi', '[mmol/m3]', tracers%data(j)%values(:,:), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
         endif

      else if (tracers%data(j)%ID==1018) then
         if (use_REcoM) then
         call def_stream((/nl-1, nod2D/),  (/nl-1, myDim_nod2D/),  'DSi', 'DSi', '[mmol/m3]', tracers%data(j)%values(:,:), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
         endif


            if (tracers%data(j)%ltra_diag) then
               ! horizontal advection
               call def_stream((/nl-1, nod2D/),  (/nl-1, myDim_nod2D/),  'DSi_hor_adv', 'Horizontal advection part of dissolved Silicic Acid', '[mmol/m3]', tracers%work%tra_advhoriz(:,:,j), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
                ! horizontal advection LO
               call def_stream((/nl-1, nod2D/),  (/nl-1, myDim_nod2D/),  'DSi_lo_hor_adv', 'LO Horizontal advection part of dissolved Silicic Acid', '[mmol/m3]', tracers%work%tra_advhoriz_LO(:,:,j), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)

               ! vertical advection
               call def_stream((/nl-1, nod2D/),  (/nl-1, myDim_nod2D/),  'DSi_ver_adv', 'Vertical advection part of dissolved Silicic Acid', '[mmol/m3]', tracers%work%tra_advvert(:,:,j), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)

               ! vertical advection LO
               call def_stream((/nl-1, nod2D/),  (/nl-1, myDim_nod2D/),  'DSi_lo_ver_adv', 'LO Vertical advection part of dissolved Silicic Acid', '[mmol/m3]', tracers%work%tra_advvert_LO(:,:,j), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)

               ! horizontal diffusion
               call def_stream((/nl-1, nod2D/),  (/nl-1, myDim_nod2D/),  'DSi_tra_diff_part_hor_redi', 'Horizontal diffusion of dissolved Silicic Acid (includes Redi diffusivity if Redi=.true.)', '[mmol/m3]', tracers%work%tra_diff_part_hor_redi(:,:,j), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)

               if (tracers%data(j)%i_vert_diff) then
               ! vertical diffusion (Explicit)
                   call def_stream((/nl-1, nod2D/),  (/nl-1, myDim_nod2D/),  'DSi_tra_diff_part_ver_expl', 'Vertical diffusion of dissolved Silicic Acid (Explicit)', '[mmol/m3]', tracers%work%tra_diff_part_ver_expl(:,:,j), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
               end if

               ! projection of horizontal Redi diffussivity onto vertical
               call def_stream((/nl-1, nod2D/),  (/nl-1, myDim_nod2D/),  'DSi_tra_diff_part_ver_redi_expl', 'Projection of horizontal Redi diffussivity onto vertical for dissolved Silicic Acid (Explicit)', '[mmol/m3]', tracers%work%tra_diff_part_ver_redi_expl(:,:,j), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)

               ! vertical diffusion (Implicit)
               call def_stream((/nl-1, nod2D/),  (/nl-1, myDim_nod2D/),  'DSi_tra_diff_part_ver_impl', 'Vertical diffusion of dissolved Silicic Acid (Implicit)', '[mmol/m3]', tracers%work%tra_diff_part_ver_impl(:,:,j), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)

               ! recom_sms
               call def_stream((/nl-1, nod2D/),  (/nl-1, myDim_nod2D/),  'DSi_recom_sms', 'Recom SMS', '[mmol/m3]', tracers%work%tra_recom_sms(:,:,j), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
            end if







      else if (tracers%data(j)%ID==1019) then
         if (use_REcoM) then
         call def_stream((/nl-1, nod2D/),  (/nl-1, myDim_nod2D/),  'DFe', 'DFe', '[mmol/m3]', tracers%data(j)%values(:,:), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
         endif

      else if (tracers%data(j)%ID==1020) then
         if (use_REcoM) then
         call def_stream((/nl-1, nod2D/),  (/nl-1, myDim_nod2D/),  'PhyCalc', 'PhyCalc', '[mmol/m3]', tracers%data(j)%values(:,:), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
         endif

      else if (tracers%data(j)%ID==1021) then
         if (use_REcoM) then
         call def_stream((/nl-1, nod2D/),  (/nl-1, myDim_nod2D/),  'DetCalc', 'DetCalc', '[mmol/m3]', tracers%data(j)%values(:,:), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
         endif

      else if (tracers%data(j)%ID==1022) then
         if (use_REcoM) then
         call def_stream((/nl-1, nod2D/),  (/nl-1, myDim_nod2D/),  'O2', 'O2', '[mmol/m3]', tracers%data(j)%values(:,:), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
         endif

      ! =====================================================================
      ! CONFIGURATION-SPECIFIC TRACERS (IDs 1023+)
      ! The meaning of these IDs changes based on model configuration!
      ! =====================================================================

      ! =====================================================================
      ! Tracer ID 1023:
      ! - In COCCOS-ONLY config: CoccoN (Coccolithophore Nitrogen)
      ! - In 3ZOO2DET configs: Zoo2N (Macrozooplankton Nitrogen)
      ! - In BASE config with rivers: DOCt (Terrestrial DOC)
      ! =====================================================================
      else if (tracers%data(j)%ID==1023) then
         if (use_REcoM .and. enable_coccos .and. .not. enable_3zoo2det) then
         ! Coccos-only configuration: This is Coccolithophore Nitrogen
         call def_stream((/nl-1, nod2D/),  (/nl-1, myDim_nod2D/),  'CoccoN', 'Coccolithophore Nitrogen', '[mmol/m3]', tracers%data(j)%values(:,:), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
         else if (use_REcoM .and. enable_3zoo2det) then
         ! 3zoo2det configurations: This is Macrozooplankton Nitrogen
         call def_stream((/nl-1, nod2D/),  (/nl-1, myDim_nod2D/),  'Zoo2N', 'Macrozooplankton Nitrogen', '[mmol/m3]', tracers%data(j)%values(:,:), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
         else if (use_REcoM .and. useRivers) then
         ! Base configuration with rivers: This is Terrestrial DOC
         call def_stream((/nl-1, nod2D/),  (/nl-1, myDim_nod2D/),  'DOCt', 'Terrestrial Dissolved Organic Carbon', '[mmol/m3]', tracers%data(j)%values(:,:), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
         endif

      !else if (tracers%data(j)%ID==1023) then
         !if (use_REcoM .and. enable_3zoo2det) then
         !call def_stream((/nl-1, nod2D/),  (/nl-1, myDim_nod2D/),  'Zoo2N', 'Intracellular conc of Nitrogen in second zooplankton', '[mmol/m3]', tracers%data(j)%values(:,:), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
         !endif

      ! =====================================================================
      ! Tracer ID 1024:
      ! - In COCCOS-ONLY config: CoccoC (Coccolithophore Carbon)
      ! - In 3ZOO2DET configs: Zoo2C (Macrozooplankton Carbon)
      ! =====================================================================
      else if (tracers%data(j)%ID==1024) then
         if (use_REcoM .and. enable_coccos .and. .not. enable_3zoo2det) then
         ! Coccos-only configuration: This is Coccolithophore Carbon
         call def_stream((/nl-1, nod2D/),  (/nl-1, myDim_nod2D/),  'CoccoC', 'Coccolithophore Carbon', '[mmol/m3]', tracers%data(j)%values(:,:), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
         else if (use_REcoM .and. enable_3zoo2det) then
         ! 3zoo2det configurations: This is Macrozooplankton Carbon
         call def_stream((/nl-1, nod2D/),  (/nl-1, myDim_nod2D/),  'Zoo2C', 'Macrozooplankton Carbon', '[mmol/m3]', tracers%data(j)%values(:,:), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
         endif

      !else if (tracers%data(j)%ID==1024) then
         !if (use_REcoM .and. enable_3zoo2det) then
         !call def_stream((/nl-1, nod2D/),  (/nl-1, myDim_nod2D/),  'Zoo2C', 'Intracellular conc of Carbon in second zooplankton', '[mmol/m3]', tracers%data(j)%values(:,:), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
         !endif

      ! =====================================================================
      ! Tracer ID 1025:
      ! - In COCCOS-ONLY config: CoccoChl (Coccolithophore Chlorophyll)
      ! - In 3ZOO2DET configs: Det2N (Second detritus pool Nitrogen)
      ! =====================================================================

      else if (tracers%data(j)%ID==1025) then
         if (use_REcoM .and. enable_coccos .and. .not. enable_3zoo2det) then
         ! Coccos-only configuration: This is Coccolithophore Chlorophyll
         call def_stream((/nl-1, nod2D/),  (/nl-1, myDim_nod2D/),  'CoccoChl', 'Coccolithophore Chlorophyll', '[mg/m3]', tracers%data(j)%values(:,:), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
         else if (use_REcoM .and. enable_3zoo2det) then
         ! 3zoo2det configurations: This is second detritus pool Nitrogen
         call def_stream((/nl-1, nod2D/),  (/nl-1, myDim_nod2D/),  'Det2N', 'Conc of N in second detritus pool', '[mmol/m3]', tracers%data(j)%values(:,:), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
         endif

      !else if (tracers%data(j)%ID==1025) then
         !if (use_REcoM .and. enable_3zoo2det) then
         !call def_stream((/nl-1, nod2D/),  (/nl-1, myDim_nod2D/),  'idetz2n', 'idetz2n', '[mmol/m3]', tracers%data(j)%values(:,:), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
         !endif

      ! =====================================================================
      ! Tracer ID 1026:
      ! - In COCCOS-ONLY config: PhaeoN (Phaeocystis Nitrogen)
      ! - In 3ZOO2DET configs: Det2C (Second detritus pool Carbon)
      ! =====================================================================
      else if (tracers%data(j)%ID==1026) then
         if (use_REcoM .and. enable_coccos .and. .not. enable_3zoo2det) then
         ! Coccos-only configuration: This is Phaeocystis Nitrogen
         call def_stream((/nl-1, nod2D/),  (/nl-1, myDim_nod2D/),  'PhaeoN', 'Phaeocystis Nitrogen', '[mmol/m3]', tracers%data(j)%values(:,:), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
         else if (use_REcoM .and. enable_3zoo2det) then
         ! 3zoo2det configurations: This is second detritus pool Carbon
         call def_stream((/nl-1, nod2D/),  (/nl-1, myDim_nod2D/),  'Det2C', 'Conc of C in second detritus pool', '[mmol/m3]', tracers%data(j)%values(:,:), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
         endif

      !else if (tracers%data(j)%ID==1026) then
         !if (use_REcoM .and. enable_3zoo2det) then
         !call def_stream((/nl-1, nod2D/),  (/nl-1, myDim_nod2D/),  'idetz2c', 'idetz2c', '[mmol/m3]', tracers%data(j)%values(:,:), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
         !endif

      ! =====================================================================
      ! Tracer ID 1027:
      ! - In COCCOS-ONLY config: PhaeoC (Phaeocystis Carbon)
      ! - In 3ZOO2DET configs: Det2Si (Second detritus pool Silica)
      ! =====================================================================
      else if (tracers%data(j)%ID==1027) then
         if (use_REcoM .and. enable_coccos .and. .not. enable_3zoo2det) then
         ! Coccos-only configuration: This is Phaeocystis Carbon
         call def_stream((/nl-1, nod2D/),  (/nl-1, myDim_nod2D/),  'PhaeoC', 'Phaeocystis Carbon', '[mmol/m3]', tracers%data(j)%values(:,:), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
         else if (use_REcoM .and. enable_3zoo2det) then
         ! 3zoo2det configurations: This is second detritus pool Silica
         call def_stream((/nl-1, nod2D/),  (/nl-1, myDim_nod2D/),  'Det2Si', 'Conc of Si in second detritus pool', '[mmol/m3]', tracers%data(j)%values(:,:), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
         endif

      !else if (tracers%data(j)%ID==1027) then
         !if (use_REcoM .and. enable_3zoo2det) then
         !call def_stream((/nl-1, nod2D/),  (/nl-1, myDim_nod2D/),  'idetz2si', 'idetz2si', '[mmol/m3]', tracers%data(j)%values(:,:), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
         !endif

      ! =====================================================================
      ! Tracer ID 1028:
      ! - In COCCOS-ONLY config: PhaeoChl (Phaeocystis Chlorophyll)
      ! - In 3ZOO2DET configs: Det2Calc (Second detritus pool Calcium)
      ! =====================================================================
      else if (tracers%data(j)%ID==1028) then
         if (use_REcoM .and. enable_coccos .and. .not. enable_3zoo2det) then
         ! Coccos-only configuration: This is Phaeocystis Chlorophyll
         call def_stream((/nl-1, nod2D/),  (/nl-1, myDim_nod2D/),  'PhaeoChl', 'Phaeocystis Chlorophyll', '[mg/m3]', tracers%data(j)%values(:,:), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
         else if (use_REcoM .and. enable_3zoo2det) then
         ! 3zoo2det configurations: This is second detritus pool Calcium
         call def_stream((/nl-1, nod2D/),  (/nl-1, myDim_nod2D/),  'Det2Calc', 'Conc of Calc in second detritus pool', '[mmol/m3]', tracers%data(j)%values(:,:), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
         endif

      !else if (tracers%data(j)%ID==1028) then
         !if (use_REcoM .and. enable_3zoo2det) then
         !call def_stream((/nl-1, nod2D/),  (/nl-1, myDim_nod2D/),  'idetz2calc', 'idetz2calc', '[mmol/m3]', tracers%data(j)%values(:,:), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
         !endif

      ! =====================================================================
      ! Tracer ID 1029:
      ! - In COCCOS-ONLY config with rivers: DOCt (Terrestrial DOC)
      ! - In 3ZOO2DET-ONLY config: Zoo3N (Microzooplankton Nitrogen)
      ! - In FULL model config: CoccoN (Coccolithophore Nitrogen)
      ! =====================================================================
      else if (tracers%data(j)%ID==1029) then
         if (use_REcoM .and. enable_coccos .and. enable_3zoo2det) then
         ! Full model: This is Coccolithophore Nitrogen
         call def_stream((/nl-1, nod2D/),  (/nl-1, myDim_nod2D/),  'CoccoN', 'Coccolithophore Nitrogen', '[mmol/m3]', tracers%data(j)%values(:,:), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
         else if (use_REcoM .and. enable_3zoo2det .and. .not. enable_coccos) then
         ! 3zoo2det-only configuration: This is Microzooplankton Nitrogen
         call def_stream((/nl-1, nod2D/),  (/nl-1, myDim_nod2D/),  'Zoo3N', 'Microzooplankton Nitrogen', '[mmol/m3]', tracers%data(j)%values(:,:), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
         else if (use_REcoM .and. enable_coccos .and. .not. enable_3zoo2det .and. useRivers) then
         ! Coccos-only configuration with rivers: This is Terrestrial DOC
         call def_stream((/nl-1, nod2D/),  (/nl-1, myDim_nod2D/),  'DOCt', 'Terrestrial Dissolved Organic Carbon', '[mmol/m3]', tracers%data(j)%values(:,:), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
         endif

      !else if (tracers%data(j)%ID==1029) then
         !if (use_REcoM .and. enable_coccos) then
         !call def_stream((/nl-1, nod2D/),  (/nl-1, myDim_nod2D/),  'CoccoN', 'CoccoN', '[mmol/m3]', tracers%data(j)%values(:,:), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
         !endif

      ! =====================================================================
      ! Tracer ID 1030:
      ! - In 3ZOO2DET-ONLY config: Zoo3C (Microzooplankton Carbon)
      ! - In FULL model config: CoccoC (Coccolithophore Carbon)
      ! =====================================================================
      else if (tracers%data(j)%ID==1030) then
         if (use_REcoM .and. enable_coccos .and. enable_3zoo2det) then
         ! Full model: This is Coccolithophore Carbon
         call def_stream((/nl-1, nod2D/),  (/nl-1, myDim_nod2D/),  'CoccoC', 'Coccolithophore Carbon', '[mmol/m3]', tracers%data(j)%values(:,:), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
         else if (use_REcoM .and. enable_3zoo2det .and. .not. enable_coccos) then
         ! 3zoo2det-only configuration: This is Microzooplankton Carbon
         call def_stream((/nl-1, nod2D/),  (/nl-1, myDim_nod2D/),  'Zoo3C', 'Microzooplankton Carbon', '[mmol/m3]', tracers%data(j)%values(:,:), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
         endif

      !else if (tracers%data(j)%ID==1030) then
         !if (use_REcoM .and. enable_coccos) then
         !call def_stream((/nl-1, nod2D/),  (/nl-1, myDim_nod2D/),  'CoccoC', 'CoccoC', '[mmol/m3]', tracers%data(j)%values(:,:), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
         !endif

      ! =====================================================================
      ! Tracer ID 1031:
      ! - In 3ZOO2DET-ONLY config with rivers: DOCt (Terrestrial DOC)
      ! - In FULL model config: CoccoChl (Coccolithophore Chlorophyll)
      ! =====================================================================
      else if (tracers%data(j)%ID==1031) then
         if (use_REcoM .and. enable_coccos .and. enable_3zoo2det) then
         ! Full model: This is Coccolithophore Chlorophyll
         call def_stream((/nl-1, nod2D/),  (/nl-1, myDim_nod2D/),  'CoccoChl', 'Coccolithophore Chlorophyll', '[mg/m3]', tracers%data(j)%values(:,:), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
         else if (use_REcoM .and. enable_3zoo2det .and. .not. enable_coccos .and. useRivers) then
         ! 3zoo2det-only configuration with rivers: This is Terrestrial DOC
         call def_stream((/nl-1, nod2D/),  (/nl-1, myDim_nod2D/),  'DOCt', 'Terrestrial Dissolved Organic Carbon', '[mmol/m3]', tracers%data(j)%values(:,:), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
         endif

      !else if (tracers%data(j)%ID==1031) then
         !if (use_REcoM .and. enable_coccos) then
         !call def_stream((/nl-1, nod2D/),  (/nl-1, myDim_nod2D/),  'CoccoChl', 'CoccoChl', '[mg/m3]', tracers%data(j)%values(:,:), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
         !endif

      ! =====================================================================
      ! Tracer IDs 1032-1034: Only in FULL model configuration
      ! These are Phaeocystis tracers
      ! =====================================================================
      else if (tracers%data(j)%ID==1032) then
         if (use_REcoM .and. enable_coccos .and. enable_3zoo2det) then
         ! Full model: This is Phaeocystis Nitrogen
         call def_stream((/nl-1, nod2D/),  (/nl-1, myDim_nod2D/),  'PhaeoN', 'Phaeocystis Nitrogen', '[mmol/m3]', tracers%data(j)%values(:,:), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
         endif

      !else if (tracers%data(j)%ID==1032) then
         !if (use_REcoM) then
         !call def_stream((/nl-1, nod2D/),  (/nl-1, myDim_nod2D/),  'Zoo3N', 'Zoo3N', '[mmol/m3]', tracers%data(j)%values(:,:), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
         !endif

      else if (tracers%data(j)%ID==1033) then
         if (use_REcoM .and. enable_coccos .and. enable_3zoo2det) then
         ! Full model: This is Phaeocystis Carbon
         call def_stream((/nl-1, nod2D/),  (/nl-1, myDim_nod2D/),  'PhaeoC', 'Phaeocystis Carbon', '[mmol/m3]', tracers%data(j)%values(:,:), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
         endif

      !else if (tracers%data(j)%ID==1033) then
         !if (use_REcoM) then
         !call def_stream((/nl-1, nod2D/),  (/nl-1, myDim_nod2D/),  'Zoo3C', 'Zoo3C', '[mmol/m3]', tracers%data(j)%values(:,:), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
         !endif

      else if (tracers%data(j)%ID==1034) then
         if (use_REcoM .and. enable_coccos .and. enable_3zoo2det) then
         ! Full model: This is Phaeocystis Chlorophyll
         call def_stream((/nl-1, nod2D/),  (/nl-1, myDim_nod2D/),  'PhaeoChl', 'Phaeocystis Chlorophyll', '[mg/m3]', tracers%data(j)%values(:,:), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
         endif

      ! =====================================================================
      ! Tracer IDs 1035-1036: Only in FULL model configuration
      ! These are Microzooplankton tracers
      ! =====================================================================
      else if (tracers%data(j)%ID==1035) then
         if (use_REcoM .and. enable_coccos .and. enable_3zoo2det) then
         ! Full model: This is Microzooplankton Nitrogen
         call def_stream((/nl-1, nod2D/),  (/nl-1, myDim_nod2D/),  'Zoo3N', 'Microzooplankton Nitrogen', '[mmol/m3]', tracers%data(j)%values(:,:), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
         endif

      else if (tracers%data(j)%ID==1036) then
         if (use_REcoM .and. enable_coccos .and. enable_3zoo2det) then
         ! Full model: This is Microzooplankton Carbon
         call def_stream((/nl-1, nod2D/),  (/nl-1, myDim_nod2D/),  'Zoo3C', 'Microzooplankton Carbon', '[mmol/m3]', tracers%data(j)%values(:,:), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
         endif
      ! =====================================================================
      ! Tracer ID 1037: Only in FULL model configuration with rivers
      ! This is Terrestrial DOC
      ! =====================================================================
      else if (tracers%data(j)%ID==1037) then
         if (use_REcoM .and. enable_coccos .and. enable_3zoo2det .and. useRivers) then
         ! Full model with rivers: This is Terrestrial DOC
         call def_stream((/nl-1, nod2D/),  (/nl-1, myDim_nod2D/),  'DOCt', 'Terrestrial Dissolved Organic Carbon', '[mmol/m3]', tracers%data(j)%values(:,:), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
        endif
      else
#endif

        if (use_transit) then
!         Transient tracers
          select case (tracers%data(j)%ID)
            case(6)
              call def_stream((/nl-1, nod2D/),  (/nl-1, myDim_nod2D/),  'sf6', 'sulfur hexafluoride', 'mol / m**3', tracers%data(index_transit_sf6)%values(:,:),  io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
            case(11)
              call def_stream((/nl-1, nod2D/),  (/nl-1, myDim_nod2D/),  'cfc11', 'chlorofluorocarbon CFC-11', 'mol / m**3', tracers%data(index_transit_f11)%values(:,:),  io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
            case(12)
              call def_stream((/nl-1, nod2D/),  (/nl-1, myDim_nod2D/),  'cfc12', 'chlorofluorocarbon CFC-12', 'mol / m**3', tracers%data(index_transit_f12)%values(:,:),  io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
            case(14)
              call def_stream((/nl-1, nod2D/),  (/nl-1, myDim_nod2D/),  'r14c', '14C / C ratio of DIC', 'none', tracers%data(index_transit_r14c)%values(:,:),  io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
            case(39)
              call def_stream((/nl-1, nod2D/),  (/nl-1, myDim_nod2D/),  'r39ar', '39Ar / Ar ratio', 'none', tracers%data(index_transit_r39ar)%values(:,:),  io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
            case default
              call def_stream((/nl-1, nod2D/),  (/nl-1, myDim_nod2D/),  'tra_'//id_string, 'passive tracer ID='//id_string, 'n/a', tracers%data(j)%values(:,:), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)        
          end select
        else
!         Other passive but not transient tracers
          call def_stream((/nl-1, nod2D/),  (/nl-1, myDim_nod2D/),  'tra_'//id_string, 'passive tracer ID='//id_string, 'n/a', tracers%data(j)%values(:,:), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
        end if

#if defined(__recom)
      end if
#endif

    end do

CASE ('sigma0      ')
    call def_stream((/nl-1, nod2D/),  (/nl-1, myDim_nod2D/),  'sigma0',      'potential density',    'kg/m3',    density_m_rho0(:,:),             io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)

!---wiso-code
!___________________________________________________________________________________________________________________________________
! output water isotopes in ocean water
CASE ('h2o18     ')
    if (lwiso) then
    call def_stream((/nl-1, nod2D/),  (/nl-1, myDim_nod2D/),  'h2o18', 'h2o18 concentration',    'kmol/m**3',    tracers%data(index_wiso_tracers(1))%values(:,:), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    end if
CASE ('hDo16     ')
    if (lwiso) then
    call def_stream((/nl-1, nod2D/),  (/nl-1, myDim_nod2D/),  'hDo16', 'hDo16 concentration',    'kmol/m**3',    tracers%data(index_wiso_tracers(2))%values(:,:), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    end if
CASE ('h2o16     ')
    if (lwiso) then
    call def_stream((/nl-1, nod2D/),  (/nl-1, myDim_nod2D/),  'h2o16', 'h2o16 concentration',    'kmol/m**3',    tracers%data(index_wiso_tracers(3))%values(:,:), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    end if
!---wiso-code-end
    
CASE ('slopetap_x   ')
    call def_stream((/nl-1,  nod2D/), (/nl-1, myDim_nod2D/),  'slopetap_x',   'neutral slope tapered X',    'none', slope_tapered(1,:,:), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
CASE ('slopetap_y   ')
    call def_stream((/nl-1,  nod2D/), (/nl-1, myDim_nod2D/),  'slopetap_y',   'neutral slope tapered Y',    'none', slope_tapered(2,:,:), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
CASE ('slopetap_z   ')
    call def_stream((/nl-1,  nod2D/), (/nl-1, myDim_nod2D/),  'slopetap_z',   'neutral slope tapered Z',    'none', slope_tapered(3,:,:), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)

CASE ('slope_x   ')
    call def_stream((/nl-1,  nod2D/), (/nl-1, myDim_nod2D/),  'slope_x',   'neutral slope X',    'none', neutral_slope(1,:,:), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
CASE ('slope_y   ')
    call def_stream((/nl-1,  nod2D/), (/nl-1, myDim_nod2D/),  'slope_y',   'neutral slope Y',    'none', neutral_slope(2,:,:), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
CASE ('slope_z   ')
    call def_stream((/nl-1,  nod2D/), (/nl-1, myDim_nod2D/),  'slope_z',   'neutral slope Z',    'none', neutral_slope(3,:,:), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)

CASE ('N2        ')
    call def_stream((/nl,    nod2D/), (/nl,   myDim_nod2D/),  'N2',        'brunt väisälä',      '1/s2', bvfreq(:,:),          io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
CASE ('Kv        ')
    call def_stream((/nl,    nod2D/), (/nl,   myDim_nod2D/),  'Kv',        'vertical diffusivity Kv',  'm2/s', Kv(:,:),        io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
CASE ('u         ')
    call def_stream((/nl-1, elem2D/), (/nl-1, myDim_elem2D/), 'u',         'zonal velocity','m/s',           dynamics%uv(1,:,:),     io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
CASE ('v         ')
    call def_stream((/nl-1, elem2D/), (/nl-1, myDim_elem2D/), 'v',         'meridional velocity','m/s',           dynamics%uv(2,:,:),     io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
CASE ('unod      ')
    call def_stream((/nl-1, nod2D/),  (/nl-1,   myDim_nod2D/),'unod',      'zonal velocity at nodes', 'm/s', dynamics%uvnode(1,:,:), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
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
    
!_______________________________________________________________________________
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
!PS CASE ('fer_GINsea')
!PS     if (Fer_GM .and. scaling_GINsea) then
!PS     call def_stream(        nod2D   ,         myDim_nod2D   , 'fer_GINsea',  'GM surface muliplicator','',  fer_GINsea_mask(:),           io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
!PS     end if
CASE ('fer_C     ')
    if (Fer_GM) then
    call def_stream(        nod2D   ,         myDim_nod2D   , 'fer_C',     'GM,   depth independent speed',  'm/s' ,   fer_c(:),                  io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    end if
CASE ('fer_gammax')
    if (Fer_GM) then
    call def_stream((/nl  , nod2D /), (/nl,   myDim_nod2D /), 'fer_gammax','GM, transport gamma_x',  'm^2/s^2', fer_gamma(1,:,:),       io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    end if
CASE ('fer_gammay')
    if (Fer_GM) then
    call def_stream((/nl  , nod2D /), (/nl,   myDim_nod2D /), 'fer_gammay','GM, transport gamma_y',  'm^2/s^2', fer_gamma(2,:,:),       io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    end if
CASE ('sigma_x   ')
    if (Fer_GM) then
    call def_stream((/nl-1, nod2D /), (/nl-1, myDim_nod2D /), 'sigma_x',   'zonal density gradient ', 'kg/m^4' , sigma_xy(1,:,:),        io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    end if
CASE ('sigma_y   ')
    if (Fer_GM) then
    call def_stream((/nl-1, nod2D /), (/nl-1, myDim_nod2D /), 'sigma_y',   'meridional density gradient ',  'kg/m^4' , sigma_xy(2,:,:),        io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    end if     
CASE ('cfl_z         ')
    call def_stream((/nl,    nod2D/), (/nl,   myDim_nod2D/),  'cfl_z',         'vertical CFL criteria',  '?',           dynamics%cfl_z(:,:),        io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
CASE ('fer_tapfac')
    if (Fer_GM) then
    call def_stream((/nl-1  , nod2D /), (/nl-1,   myDim_nod2D /), 'fer_tapfac','tapering factor',  '', fer_tapfac(:,:),       io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    end if    
CASE ('redi_K        ')    
    sel_redi=1
    if (Redi) then
        call def_stream((/nl-1  , nod2D /), (/nl-1,   myDim_nod2D /), 'redi_K',   'Redi diffusion coefficient', 'm2/s', Ki(:,:),   io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    end if
!_______________________________________________________________________________
! Density MOC diagnostics
CASE ('dMOC      ')
    sel_dmoc = 1
    if (ldiag_dMOC) then
       call def_stream((/std_dens_N, elem2D/),  (/std_dens_N, myDim_elem2D/), 'U_rho_x_DZ',     'fluxes for density MOC', 'fluxes', std_dens_UVDZ(1,:,:),   io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
       call def_stream((/std_dens_N, elem2D/),  (/std_dens_N, myDim_elem2D/), 'V_rho_x_DZ',     'fluxes for density MOC', 'fluxes', std_dens_UVDZ(2,:,:),   io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
       call def_stream((/std_dens_N, elem2D/),  (/std_dens_N, myDim_elem2D/), 'std_heat_flux',  'HF bouyancy flux      ', 'kg*m/s' ,std_dens_flux(1,:,:),   io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
       call def_stream((/std_dens_N, elem2D/),  (/std_dens_N, myDim_elem2D/), 'std_rest_flux',  'RESTOR. bouyancy flux ', 'kg*m/s' ,std_dens_flux(2,:,:),   io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
       call def_stream((/std_dens_N, elem2D/),  (/std_dens_N, myDim_elem2D/), 'std_frwt_flux',  'FW bouyancy flux      ', 'kg*m/s' ,std_dens_flux(3,:,:),   io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
       call def_stream((/std_dens_N, elem2D/),  (/std_dens_N, myDim_elem2D/), 'std_dens_dVdT',  'dV/dT',                  'm3/s'   ,std_dens_dVdT(:,:),     io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
       call def_stream((/std_dens_N, nod2D /),  (/std_dens_N,  myDim_nod2D/), 'std_dens_DIV',   'm3/s',                   'm3/s'   ,std_dens_DIV(:,:),      io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
       if (Fer_GM) call def_stream((/std_dens_N, nod2D /),  (/std_dens_N,  myDim_nod2D/), 'std_dens_DIVbolus',   'm3/s',  'm3/s'   ,std_dens_DIV_fer(:,:),  io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
       call def_stream((/std_dens_N, elem2D/),  (/std_dens_N, myDim_elem2D/), 'std_dens_Z',     'm',                      'm'      ,std_dens_Z(:,:),        io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
       call def_stream((/std_dens_N, elem2D/),  (/std_dens_N, myDim_elem2D/), 'std_dens_H'    , 'density thickness'     , 'm'     , std_dens_H(:,:),        io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
       call def_stream((/nl-1,       nod2D /),  (/nl-1,       myDim_nod2D /), 'density_dMOC',   'density'               , 'm',      density_dmoc(:,:),      io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
       call def_stream(elem2D,                                myDim_elem2D  , 'density_flux_e', 'density flux at elems ', 'm',      dens_flux_e(:),         io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    end if
    
!_______________________________________________________________________________
! PGF (pressure gradient force) diagnostic
CASE ('pgf_x     ')    
    call def_stream((/nl-1, elem2D/), (/nl-1, myDim_elem2D/), 'pgf_x', 'zonal pressure gradient force'     , 'm/s^2', pgf_x(:,:), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
CASE ('pgf_y     ')    
    call def_stream((/nl-1, elem2D/), (/nl-1, myDim_elem2D/), 'pgf_y', 'meridional pressure gradient force', 'm/s^2', pgf_y(:,:), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)

!___________________________________________________________________________________________________________________________________
! ALE layer thickness 
CASE ('hnode     ')    
    call def_stream((/nl-1, nod2D/) , (/nl-1, myDim_nod2D/) , 'hnode'    , 'vertice layer thickness'  , 'm', hnode(:,:), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
CASE ('hnode_new  ')
    call def_stream((/nl-1, nod2D/),  (/nl-1, myDim_nod2D/) , 'hnode_new', 'hnode_new'                , 'm', mesh%hnode_new(:,:), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)    
CASE ('helem     ')    
    call def_stream((/nl-1, elem2D/), (/nl-1, myDim_elem2D/), 'helem'    , 'elemental layer thickness', 'm', helem(:,:), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    
!___________________________________________________________________________________________________________________________________    
#if defined (__oifs) || defined (__ifsinterface)
CASE ('alb       ')
  call def_stream(nod2D, myDim_nod2D, 'alb',    'ice albedo',              'none',   ice%atmcoupl%ice_alb(:),               io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
CASE ('ist       ')
  call def_stream(nod2D, myDim_nod2D, 'ist',    'ice surface temperature', 'K',      ice%data(4)%values(:),                 io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
CASE ('qsi       ')
  call def_stream(nod2D, myDim_nod2D, 'qsi',    'ice heat flux',           'W/m^2',  ice%atmcoupl%ice_flx_h(:),             io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
CASE ('qso       ')
  call def_stream(nod2D, myDim_nod2D, 'qso',    'oce heat flux',           'W/m^2',  ice%atmcoupl%oce_flx_h(:),             io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
CASE ('runoff_liquid  ')
  call def_stream(nod2D, myDim_nod2D, 'runoff_liquid',  'liquid water runoff',     'm/s',    ice%atmcoupl%runoff_liquid(:),        io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
CASE ('runoff_solid  ')
  call def_stream(nod2D, myDim_nod2D, 'runoff_solid',  'solid water runoff',     'm/s',    ice%atmcoupl%runoff_solid(:),        io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
CASE ('enthalpy  ')
  call def_stream(nod2D, myDim_nod2D, 'enth',  'enthalpy of fusion of solid water runoff',     'W/m^2',    ice%atmcoupl%enthalpyoffuse(:),        io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
CASE ('qcon      ')
     call def_stream(nod2D, myDim_nod2D, 'qcon',  'conductive heat flux',   'W/m^2',    ice%atmcoupl%flx_qcon(:),           io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)  
CASE ('qres      ')
     call def_stream(nod2D, myDim_nod2D, 'qres',  'residual heat flux',     'W/m^2',    ice%atmcoupl%flx_qres(:),           io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
#endif

!------------------------------------------
! LA 2023-01-31 adding iceberg outputs
CASE ('icb       ')
  if (use_icebergs) then
    call def_stream(nod2D, myDim_nod2D, 'ibfwb',   'basal iceberg melting',            'm/s',    ibfwb(:),         1, 'm', i_real4, partit, mesh)
    call def_stream(nod2D, myDim_nod2D, 'ibfwbv',  'basal iceberg melting',            'm/s',    ibfwbv(:),        1, 'm', i_real4, partit, mesh)
    call def_stream(nod2D, myDim_nod2D, 'ibfwl',   'lateral iceberg melting',          'm/s',    ibfwl(:),         1, 'm', i_real4, partit, mesh)
    call def_stream(nod2D, myDim_nod2D, 'ibfwe',   'iceberg erosion',                  'm/s',    ibfwe(:),         1, 'm', i_real4, partit, mesh)
    call def_stream((/nl,nod2D/), (/nl,myDim_nod2D/), 'ibhf',    'heat flux from iceberg melting',   'W/m2',    ibhf_n(:,:),      1, 'm', i_real4, partit, mesh)
  end if

#if defined (__cvmix)    
!_______________________________________________________________________________
! TKE mixing diagnostic 
CASE ('TKE       ')
    if (mix_scheme_nmb==5 .or. mix_scheme_nmb==56) then
        call def_stream((/nl,nod2D/), (/nl,myDim_nod2D/), 'tke'     , 'turbulent kinetic energy'                    , 'm^2/s^2', tke(:,:)     , io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
        call def_stream((/nl,nod2D/), (/nl,myDim_nod2D/), 'tke_Ttot', 'total production of turbulent kinetic energy', 'm^2/s^3', tke_Ttot(:,:), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
        call def_stream((/nl,nod2D/), (/nl,myDim_nod2D/), 'tke_Tbpr', 'TKE production by buoyancy'                  , 'm^2/s^3', tke_Tbpr(:,:), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
        call def_stream((/nl,nod2D/), (/nl,myDim_nod2D/), 'tke_Tspr', 'TKE production by shear'                     , 'm^2/s^3', tke_Tspr(:,:), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
        call def_stream((/nl,nod2D/), (/nl,myDim_nod2D/), 'tke_Tdif', 'TKE production by vertical diffusion'        , 'm^2/s^3', tke_Tdif(:,:), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
        call def_stream((/nl,nod2D/), (/nl,myDim_nod2D/), 'tke_Tdis', 'TKE production by dissipation'               , 'm^2/s^3', tke_Tdis(:,:), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
        call def_stream((/nl,nod2D/), (/nl,myDim_nod2D/), 'tke_Twin', 'TKE production by wind'                      , 'm^2/s^3', tke_Twin(:,:), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
        call def_stream((/nl,nod2D/), (/nl,myDim_nod2D/), 'tke_Tbck', 'background forcing for TKE'                  , 'm^2/s^3', tke_Tbck(:,:), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
        call def_stream((/nl,nod2D/), (/nl,myDim_nod2D/), 'tke_Lmix', 'mixing length scale of TKE'                  , 'm'      , tke_Lmix(:,:), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
        call def_stream((/nl,nod2D/), (/nl,myDim_nod2D/), 'tke_Pr'  , 'Prantl number'                               , ''       , tke_Pr(:,:)  , io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
        if (mix_scheme_nmb==56) then
            ! TKE-IDEMIX diagnostic 
            call def_stream((/nl,nod2D/), (/nl,myDim_nod2D/), 'tke_Tiwf', 'TKE production by internal waves (IDEMIX)','m^2/s^3', tke_Tiwf(:,:), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
        end if 
    end if
    
!_______________________________________________________________________________
! IDEMIX mixing Internal-Wave-Energy diagnostics
CASE ('IDEMIX    ')
    if (mod(mix_scheme_nmb,10)==6) then
        call def_stream((/nl,elem2D/), (/nl,myDim_elem2D/), 'iwe'     , 'internal wave energy'                      , 'm^2/s^2', iwe(:,:)     , io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
        call def_stream((/nl,elem2D/), (/nl,myDim_elem2D/), 'iwe_Ttot', 'total production of internal wave energy'  , 'm^2/s^2', iwe_Ttot(:,:), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
        call def_stream((/nl,elem2D/), (/nl,myDim_elem2D/), 'iwe_Tdif', 'IWE production by vertical diffusion'      , 'm^2/s^3', iwe_Tdif(:,:), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
        call def_stream((/nl,elem2D/), (/nl,myDim_elem2D/), 'iwe_Tdis', 'IWE production by dissipation'             , 'm^2/s^3', iwe_Tdis(:,:), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
        call def_stream((/nl,elem2D/), (/nl,myDim_elem2D/), 'iwe_Tsur', 'IWE production from surface forcing'       , 'm^2/s^2', iwe_Tsur(:,:), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
        call def_stream((/nl,elem2D/), (/nl,myDim_elem2D/), 'iwe_Tbot', 'IWE production from bottom forcing'        , 'm^2/s^2', iwe_Tbot(:,:), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
        call def_stream((/nl,elem2D/), (/nl,myDim_elem2D/), 'iwe_Thdi', 'IWE production from hori. diffusion'       , 'm^2/s^2', iwe_Thdi(:,:), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
        call def_stream((/nl,elem2D/), (/nl,myDim_elem2D/), 'iwe_c0'  , 'IWE vertical group velocity'               , 'm/s'    , iwe_c0(:,:)  , io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
        call def_stream((/nl,elem2D/), (/nl,myDim_elem2D/), 'iwe_v0'  , 'IWE horizontal group velocity'             , 'm/s'    , iwe_c0(:,:)  , io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
        call def_stream(elem2D       , myDim_elem2D       , 'iwe_fbot', 'IDEMIX bottom forcing'                     , 'm^3/s^3', iwe_fbot(:)  , io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
        call def_stream(elem2D       , myDim_elem2D       , 'iwe_fsrf', 'IDEMIX surface forcing'                    , 'm^3/s^3', iwe_fsrf(:)  , io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    end if    

!_______________________________________________________________________________
! TIDAL mixing diagnostics
CASE ('TIDAL     ')
    if (mod(mix_scheme_nmb,10)==7) then
        ! cvmix_TIDAL diagnostics
        call def_stream((/nl,elem2D/), (/nl,myDim_elem2D/), 'tidal_Kv'  , 'tidal diffusivity'                       , 'm^2/s'  , tidal_Kv(:,:), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
        call def_stream((/nl,elem2D/), (/nl,myDim_elem2D/), 'tidal_Av'  , 'tidal viscosity'                         , 'm^2/s'  , tidal_Av(:,:), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
        call def_stream(elem2D       , myDim_elem2D       , 'tidal_fbot', 'near tidal bottom forcing'               , 'W/m^2'  , tidal_fbot   , io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    end if
#endif     

!_______________________________________________________________________________
! TIDAL mixing diagnostics
CASE ('FORC      ')
    if (ldiag_forc) then
        if (sel_forcvar( 1)==0) call def_stream(nod2D , myDim_nod2D , 'uwind' , '10m zonal surface wind velocity', 'm/s'  , u_wind(:)        , io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
        if (sel_forcvar( 2)==0) call def_stream(nod2D , myDim_nod2D , 'vwind' , '10m merid surface wind velocity', 'm/s'  , v_wind(:)        , io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
        if (sel_forcvar( 3)==0) call def_stream(nod2D , myDim_nod2D , 'tair'  , 'surface air temperature'        , '°C'   , Tair(:)          , io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
        if (sel_forcvar( 4)==0) call def_stream(nod2D , myDim_nod2D , 'shum'  , 'specific humidity'              , ''     , shum(:)          , io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
        if (sel_forcvar( 5)==0) call def_stream(nod2D , myDim_nod2D , 'prec'  , 'precicipation rain'             , 'm/s'  , prec_rain(:)     , io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
        if (sel_forcvar( 6)==0) call def_stream(nod2D , myDim_nod2D , 'snow'  , 'precicipation snow'             , 'm/s'  , prec_snow(:)     , io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
        if (sel_forcvar( 7)==0) call def_stream(nod2D , myDim_nod2D , 'evap'  , 'total evaporation'              , 'm/s'  , evaporation(:)   , io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
        if (sel_forcvar( 8)==0) call def_stream(nod2D , myDim_nod2D , 'swr'   , 'short wave radiation'           , 'W/m^2', shortwave(:)     , io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
        if (sel_forcvar( 9)==0) call def_stream(nod2D , myDim_nod2D , 'lwr'   , 'long wave radiation'            , 'W/m^2', longwave(:)      , io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
        if (sel_forcvar(10)==0) call def_stream(nod2D , myDim_nod2D , 'runoff', 'river runoff'                   , 'none' , runoff(:)        , io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
        if (sel_forcvar(11)==0) call def_stream(elem2D, myDim_elem2D, 'tx_sur', 'zonal wind str. to ocean'       , 'N/m^2', stress_surf(1, :), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
        if (sel_forcvar(12)==0) call def_stream(elem2D, myDim_elem2D, 'ty_sur', 'meridional wind str. to ocean'  , 'N/m^2', stress_surf(2, :), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
        
        call def_stream(nod2D , myDim_nod2D , 'cd',    'wind drag coef. '             , '',     cd_atm_oce_arr(:), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
        call def_stream(nod2D , myDim_nod2D , 'ch',    'transfer coeff. sensible heat', '',     ch_atm_oce_arr(:), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
        call def_stream(nod2D , myDim_nod2D , 'ce',    'transfer coeff. evaporation ' , '',     ce_atm_oce_arr(:), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
#if defined (__oasis)
        call def_stream(nod2D,  myDim_nod2D,  'subli', 'sublimation',                   'm/s',  sublimation(:)   , io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
#endif
        if ((use_virt_salt) .or. ( (.not. use_virt_salt) .and. (use_cavity) )) then
            if (sel_forcvar(13)==0) call def_stream(nod2D , myDim_nod2D , 'virtsalt', 'virtual salt flux'        , 'm/s*psu', virtual_salt(:)  , io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
        end if     
#if !defined (__oasis)        
        if (sel_forcvar(14)==0) call def_stream(nod2D , myDim_nod2D , 'relaxsalt', 'relaxation salt flux'        , 'm/s*psu', relax_salt(:)    , io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
#endif  
        if (sel_forcvar(15)==0) call def_stream(nod2D , myDim_nod2D , 'realsalt' , 'real salt flux from sea ice' , 'm/s*psu', real_salt_flux(:), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    end if 
    
!_______________________________________________________________________________
! Discrete Variance Decay (DVD) diagnostic after Klingbeil etal 2014, and 
! Tridib et al 2023
CASE ('DVD       ')
    sel_dvd=1
    if (ldiag_DVD) then
        !___temperature DVD_____________________________________________________
        call def_stream((/nl-1, nod2D/), (/nl-1, myDim_nod2D/), 'dvd_temp_KK_tot' ,  'tot. temperature DVD \n (Klingbeil et al. 2014)', 'K^2/s' , dvd_KK_tot(        :,:,1), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
        call def_stream((/nl-1, nod2D/), (/nl-1, myDim_nod2D/), 'dvd_temp_SD_tot' ,  'tot. temperature DVD \n (Banerjee et al. 2023)' , 'K^2/s' , dvd_SD_tot(        :,:,1), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
        call def_stream((/nl-1, nod2D/), (/nl-1, myDim_nod2D/), 'dvd_temp_SD_advh',  'temperature DVD horiz. adv.'                    , 'K^2/s' , dvd_SD_chi_adv_h(  :,:,1), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
        call def_stream((/nl-1, nod2D/), (/nl-1, myDim_nod2D/), 'dvd_temp_SD_advv',  'temperature DVD vert. adv. '                    , 'K^2/s' , dvd_SD_chi_adv_v(  :,:,1), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
        if (Redi .eqv. .False. .and. K_hor /= 0.0_WP) then
            call def_stream((/nl-1, nod2D/), (/nl-1, myDim_nod2D/), 'dvd_temp_SD_difh',  'temperature DVD horiz. diff. expl.'             , 'K^2/s' , dvd_SD_chi_dif_he( :,:,1), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
        end if
        call def_stream((/nl-1, nod2D/), (/nl-1, myDim_nod2D/), 'dvd_temp_SD_difvi', 'temperature DVD vert. diff. impl.'              , 'K^2/s' , dvd_SD_chi_dif_vi( :,:,1), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
        call def_stream((/nl-1, nod2D/), (/nl-1, myDim_nod2D/), 'dvd_temp_SD_difsbc', 'temperature DVD sbc.'                          , 'K^2/s' , dvd_SD_chi_dif_sbc(:,:,1), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
        if (.not. tracers%data(1)%i_vert_diff) then
            call def_stream((/nl-1, nod2D/), (/nl-1, myDim_nod2D/), 'dvd_temp_SD_difve', 'temperature DVD vert. diff. expl.'          , 'K^2/s' , dvd_SD_chi_dif_ve( :,:,1), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
        end if
        if (tracers%data(1)%smooth_bh_tra) then 
            call def_stream((/nl-1, nod2D/), (/nl-1, myDim_nod2D/), 'dvd_temp_SD_difhbh','temperature DVD horiz. diff. biharm'        , 'K^2/s' , dvd_SD_chi_dif_hbh(:,:,1), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
        end if 
        if (Redi) then 
            call def_stream((/nl-1, nod2D/), (/nl-1, myDim_nod2D/), 'dvd_temp_SD_difheR', 'temperature DVD horiz. redi diff. expl.'   , 'K^2/s' , dvd_SD_chi_dif_heR(:,:,1), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
            
            call def_stream((/nl-1, nod2D/), (/nl-1, myDim_nod2D/), 'dvd_temp_SD_difveR', 'temperature DVD vert. redi diff. expl.'    , 'K^2/s' , dvd_SD_chi_dif_veR(:,:,1), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
            call def_stream((/nl-1, nod2D/), (/nl-1, myDim_nod2D/), 'dvd_temp_SD_difviR', 'temperature DVD vert. redi diff. impl.'    , 'K^2/s' , dvd_SD_chi_dif_viR(:,:,1), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
        end if     
        
        !___salinity DVD________________________________________________________
        call def_stream((/nl-1, nod2D/), (/nl-1, myDim_nod2D/), 'dvd_salt_KK_tot' ,  'tot. salinity DVD \n (Klingbeil et al. 2014)'   , 'K^2/s' , dvd_KK_tot(        :,:,2), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
        call def_stream((/nl-1, nod2D/), (/nl-1, myDim_nod2D/), 'dvd_salt_SD_tot' ,  'tot. salinity DVD \n (Banerjee et al. 2023)'    , 'K^2/s' , dvd_SD_tot(        :,:,2), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
        call def_stream((/nl-1, nod2D/), (/nl-1, myDim_nod2D/), 'dvd_salt_SD_advh',  'salinity DVD horiz. adv.'                       , 'K^2/s' , dvd_SD_chi_adv_h(  :,:,2), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
        call def_stream((/nl-1, nod2D/), (/nl-1, myDim_nod2D/), 'dvd_salt_SD_advv',  'salinity DVD vert. adv. '                       , 'K^2/s' , dvd_SD_chi_adv_v(  :,:,2), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
        if (Redi .eqv. .False. .and. K_hor /= 0.0_WP) then
            call def_stream((/nl-1, nod2D/), (/nl-1, myDim_nod2D/), 'dvd_salt_SD_difh',  'salinity DVD horiz. diff. expl.'                , 'K^2/s' , dvd_SD_chi_dif_he( :,:,2), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
        end if 
        call def_stream((/nl-1, nod2D/), (/nl-1, myDim_nod2D/), 'dvd_salt_SD_difvi', 'salinity DVD vert. diff. impl.'                 , 'K^2/s' , dvd_SD_chi_dif_vi( :,:,2), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
        call def_stream((/nl-1, nod2D/), (/nl-1, myDim_nod2D/), 'dvd_salt_SD_difsbc','salinity DVD sbc.'                              , 'K^2/s' , dvd_SD_chi_dif_sbc(:,:,2), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
        if (.not. tracers%data(1)%i_vert_diff) then
            call def_stream((/nl-1, nod2D/), (/nl-1, myDim_nod2D/), 'dvd_salt_SD_difve', 'salinity DVD vert. diff. expl.'             , 'K^2/s' , dvd_SD_chi_dif_ve( :,:,2), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
        end if
        if (tracers%data(1)%smooth_bh_tra) then 
            call def_stream((/nl-1, nod2D/), (/nl-1, myDim_nod2D/), 'dvd_salt_SD_difhbh','salinity DVD horiz. diff. biharm'           , 'K^2/s' , dvd_SD_chi_dif_hbh(:,:,2), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
        end if 
        if (Redi) then 
            call def_stream((/nl-1, nod2D/), (/nl-1, myDim_nod2D/), 'dvd_salt_SD_difheR', 'salinity DVD horiz. redi diff. expl.'      , 'K^2/s' , dvd_SD_chi_dif_heR(:,:,2), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
            
            call def_stream((/nl-1, nod2D/), (/nl-1, myDim_nod2D/), 'dvd_salt_SD_difveR', 'salinity DVD vert. redi diff. expl.'       , 'K^2/s' , dvd_SD_chi_dif_veR(:,:,2), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
            call def_stream((/nl-1, nod2D/), (/nl-1, myDim_nod2D/), 'dvd_salt_SD_difviR', 'salinity DVD vert. redi diff. impl.'       , 'K^2/s' , dvd_SD_chi_dif_viR(:,:,2), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
        end if     
        
    end if !--> if (ldiag_DVD) then

!_______________________________________________________________________________
! Split-Explicite subcycling varaibles    
CASE ('SPLIT-EXPL')
    if (dynamics%use_ssh_se_subcycl) then 
        call def_stream(elem2D, myDim_elem2D , 'ubt'      , 'zonal barotrop. transport'          , '?' , dynamics%se_uvBT(1,:)      , io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
        call def_stream(elem2D, myDim_elem2D , 'vbt'      , 'merid. barotrop. transport'         , '?' , dynamics%se_uvBT(2,:)      , io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
        call def_stream(elem2D, myDim_elem2D , 'ubt_rhs'  , 'zonal barotrop. transport RHS'      , '?' , dynamics%se_uvBT_rhs(1,:)  , io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
        call def_stream(elem2D, myDim_elem2D , 'vbt_rhs'  , 'merid. barotrop. transport RHS'     , '?' , dynamics%se_uvBT_rhs(2,:)  , io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
        call def_stream(elem2D, myDim_elem2D , 'ubt_12'   , 'zonal barotrop. transport 12'       , '?' , dynamics%se_uvBT_12(1,:)   , io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
        call def_stream(elem2D, myDim_elem2D , 'vbt_12'   , 'merid. barotrop. transport 12'      , '?' , dynamics%se_uvBT_12(2,:)   , io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
        call def_stream(elem2D, myDim_elem2D , 'ubt_theta', 'zonal barotrop. transport theta'    , '?' , dynamics%se_uvBT_theta(1,:), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
        call def_stream(elem2D, myDim_elem2D , 'vbt_theta', 'merid. barotrop. transport theta'   , '?' , dynamics%se_uvBT_theta(2,:), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
        call def_stream(elem2D, myDim_elem2D , 'ubt_mean' , 'zonal barotrop. transport mean'     , '?' , dynamics%se_uvBT_mean(1,:) , io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
        call def_stream(elem2D, myDim_elem2D , 'vbt_mean' , 'merid. barotrop. transport mean'    , '?' , dynamics%se_uvBT_mean(2,:) , io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
        call def_stream((/nl-1,elem2D/), (/nl-1,myDim_elem2D/) , 'u_rhs' , 'zonal transport rhs' , '?' , dynamics%uv_rhs(1,:,:)     , io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
        call def_stream((/nl-1,elem2D/), (/nl-1,myDim_elem2D/) , 'v_rhs' , 'merid. transport rhs', '?' , dynamics%uv_rhs(2,:,:)     , io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
        call def_stream((/nl-1,elem2D/), (/nl-1,myDim_elem2D/) , 'uh' , 'zonal transport'        , '?' , dynamics%se_uvh(1,:,:)     , io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
        call def_stream((/nl-1,elem2D/), (/nl-1,myDim_elem2D/) , 'vh' , 'merid. transport'       , '?' , dynamics%se_uvh(2,:,:)     , io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
        if (dynamics%se_visc) then
            call def_stream(elem2D, myDim_elem2D , 'ubt_hvisc'    , 'zonal  hvisc. stabil.'      , '?' , dynamics%se_uvBT_stab_hvisc(1,:), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
            call def_stream(elem2D, myDim_elem2D , 'vbt_hvisc'    , 'merid. hvisc. stabil.'      , '?' , dynamics%se_uvBT_stab_hvisc(2,:), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
        end if 
    end if

!_______________________________________________________________________________
! compute squared velocities of u, v, w
CASE ('UVW_SQR   ')
    if (ldiag_uvw_sqr) then
        !___temperature DVD_____________________________________________________
        call def_stream((/nl-1, elem2D/), (/nl-1, myDim_elem2D/), 'u2' , 'squared zonal velocity'     , 'm^2/s^2' , uv2(1,:,:), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
        call def_stream((/nl-1, elem2D/), (/nl-1, myDim_elem2D/), 'v2' , 'squared meridional velocity', 'm^2/s^2' , uv2(2,:,:), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
        call def_stream((/nl-1, nod2D/) , (/nl-1, myDim_nod2D/) , 'w2' , 'squared vertical velocity'  , 'm^2/s^2' , wvel2(:,:), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    end if !--> if (ldiag_DVD) then    
    
!_______________________________________________________________________________
! compute horizontal and vertical tracer gradients
CASE ('TRGRD_XYZ ')
    sel_trgrd_xyz=1
    if (ldiag_trgrd_xyz) then
        call def_stream((/nl-1, elem2D/), (/nl-1, myDim_elem2D/),  'temp_grdx',   'zonal temperature gradient',        'K/m', trgrd_x(1,:,:), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
        call def_stream((/nl-1, elem2D/), (/nl-1, myDim_elem2D/),  'temp_grdy',   'meridional temperature gradient',   'K/m', trgrd_y(1,:,:), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
        call def_stream((/nl-1,  nod2D/), (/nl-1, myDim_nod2D /),  'temp_grdz',   'vertical temperature gradient',     'K/m', trgrd_z(1,:,:), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
        
        call def_stream((/nl-1, elem2D/), (/nl-1, myDim_elem2D/),  'salt_grdx',   'zonal salinity gradient',         'psu/m', trgrd_x(2,:,:), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
        call def_stream((/nl-1, elem2D/), (/nl-1, myDim_elem2D/),  'salt_grdy',   'meridional salinity gradient',    'psu/m', trgrd_y(2,:,:), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
        call def_stream((/nl-1,  nod2D/), (/nl-1, myDim_nod2D /),  'salt_grdz',   'vertical salinity gradient',      'psu/m', trgrd_z(2,:,:), io_list(i)%freq, io_list(i)%unit, io_list(i)%precision, partit, mesh)
    end if !--> if (ldiag_DVD) then    
    
!_______________________________________________________________________________
CASE DEFAULT
    if (mype==0) write(*,*) 'stream ', io_list(i)%id, ' is not defined !'
END SELECT ! --> SELECT CASE (trim(io_list(i)%id))
END DO ! --> DO i=1, io_listsize

    if (lnextGEMS) then        
        call def_stream((/nlev_upper,   nod2D/),  (/nlev_upper,   myDim_nod2D/),  'temp_upper',  'temperature', 'C',         tracers%data(1)%values(:nlev_upper,:),  3, 'h', 4, partit, mesh)
        call def_stream((/nlev_upper,   nod2D/),  (/nlev_upper,   myDim_nod2D/),  'salt_upper',  'salinity',    'psu',        tracers%data(2)%values(:nlev_upper,:),  3, 'h', 8, partit, mesh)
        call def_stream((/nlev_upper,   elem2D/), (/nlev_upper,   myDim_elem2D/), 'u_upper',     'zonal velocity','m/s',      dynamics%uv(1,:nlev_upper,:),           3, 'h', 4, partit, mesh)
        call def_stream((/nlev_upper,   elem2D/), (/nlev_upper,   myDim_elem2D/), 'v_upper',     'meridional velocity','m/s', dynamics%uv(2,:nlev_upper,:),           3, 'h', 4, partit, mesh)
        call def_stream((/nlev_upper+1, nod2D/),  (/nlev_upper+1, myDim_nod2D/),  'w_upper',     'vertical velocity',  'm/s', dynamics%w(:nlev_upper+1,:),            3, 'h', 4, partit, mesh)
    end if
    if (ldiag_ice) then
        call def_stream(nod2D,  myDim_nod2D,  'vol_ice',   'ice volume',   'm',  vol_ice(:),   1, 'd', 8, partit, mesh)
        call def_stream(nod2D,  myDim_nod2D,  'vol_snow',  'snow volume',  'm',  vol_snow(:),  1, 'd', 8, partit, mesh)
    end if

    !___________________________________________________________________________
    ! Richardson number diagnostics  without predefined freq, freq_unit, prec, -->
    ! default monthly output
    if (ldiag_Ri) then
        call def_stream((/nl,  nod2D/), (/nl, myDim_nod2D/), 'shear',   'du2/dz2',           '1/s2',     shear(:,:), 1, 'm', i_real8, partit, mesh)
        call def_stream((/nl,  nod2D/), (/nl, myDim_nod2D/), 'Ri',      'Richardson number', 'no dim.',  Ri(:,:),    1, 'm', i_real8, partit, mesh)
    end if
    !___________________________________________________________________________
    ! Turbulent flux diagnostics  without predefined freq, freq_unit, prec, -->
    ! default monthly output
    if (ldiag_turbflux) then
        call def_stream((/nl,  nod2D/), (/nl, myDim_nod2D/), 'KvdTdz',   'KvdTdz',           'K m/s',     KvdTdz(:,:), 1, 'm', i_real8, partit, mesh)
        call def_stream((/nl,  nod2D/), (/nl, myDim_nod2D/), 'KvdSdz',   'KvdSdz',           'PSU m/s',   KvdSdz(:,:), 1, 'm', i_real8, partit, mesh)
    end if
    !___________________________________________________________________________
    ! Tracers flux diagnostics without predefined freq, freq_unit, prec, -->
    ! default monthly output
    if (ldiag_trflx .and. sel_trgrd_xyz==0) then
        call def_stream((/nl-1,  elem2D/), (/nl-1, myDim_elem2D/), 'utemp',   'u*temp',           'm/s*°C',     tuv(1,:,:), 1, 'm', i_real8, partit, mesh)
        call def_stream((/nl-1,  elem2D/), (/nl-1, myDim_elem2D/), 'vtemp',   'v*temp',           'm/s*°C',     tuv(2,:,:), 1, 'm', i_real8, partit, mesh)
        call def_stream((/nl-1,  elem2D/), (/nl-1, myDim_elem2D/), 'usalt',   'u*salt',           'm/s*psu',   suv(1,:,:), 1, 'm', i_real8, partit, mesh)
        call def_stream((/nl-1,  elem2D/), (/nl-1, myDim_elem2D/), 'vsalt',   'v*salt',           'm/s*psu',   suv(2,:,:), 1, 'm', i_real8, partit, mesh)
    end if
    
    !___________________________________________________________________________
    ! Density MOC diagnostic, without predefined freq, freq_unit, prec, -->
    ! default annual output
    if (ldiag_dMOC .and. sel_dmoc == 0) then
       call def_stream((/std_dens_N, elem2D/),  (/std_dens_N, myDim_elem2D/), 'U_rho_x_DZ',     'fluxes for density MOC', 'fluxes', std_dens_UVDZ(1,:,:),   1, 'y', i_real4, partit, mesh)
       call def_stream((/std_dens_N, elem2D/),  (/std_dens_N, myDim_elem2D/), 'V_rho_x_DZ',     'fluxes for density MOC', 'fluxes', std_dens_UVDZ(2,:,:),   1, 'y', i_real4, partit, mesh)
       call def_stream((/std_dens_N, elem2D/),  (/std_dens_N, myDim_elem2D/), 'std_heat_flux',  'HF bouyancy flux      ', 'kg*m/s' ,std_dens_flux(1,:,:),   1, 'y', i_real4, partit, mesh)
       call def_stream((/std_dens_N, elem2D/),  (/std_dens_N, myDim_elem2D/), 'std_rest_flux',  'RESTOR. bouyancy flux ', 'kg*m/s' ,std_dens_flux(2,:,:),   1, 'y', i_real4, partit, mesh)
       call def_stream((/std_dens_N, elem2D/),  (/std_dens_N, myDim_elem2D/), 'std_frwt_flux',  'FW bouyancy flux      ', 'kg*m/s' ,std_dens_flux(3,:,:),   1, 'y', i_real4, partit, mesh)
       call def_stream((/std_dens_N, elem2D/),  (/std_dens_N, myDim_elem2D/), 'std_dens_dVdT',  'dV/dT',                  'm3/s'   ,std_dens_dVdT(:,:),     1, 'y', i_real4, partit, mesh)
       call def_stream((/std_dens_N, nod2D /),  (/std_dens_N,  myDim_nod2D/), 'std_dens_DIV',   'm3/s',                   'm3/s'   ,std_dens_DIV(:,:),      1, 'y', i_real4, partit, mesh)
       if (Fer_GM) call def_stream((/std_dens_N, nod2D /),  (/std_dens_N,  myDim_nod2D/), 'std_dens_DIVbolus',   'm3/s',                   'm3/s'   ,std_dens_DIV_fer(:,:),  1, 'y', i_real4, partit, mesh)
       call def_stream((/std_dens_N, elem2D/),  (/std_dens_N, myDim_elem2D/), 'std_dens_Z',     'm',                      'm'      ,std_dens_Z(:,:),        1, 'y', i_real4, partit, mesh)
       call def_stream((/std_dens_N, elem2D/),  (/std_dens_N, myDim_elem2D/), 'std_dens_H'    , 'density thickness'     , 'm'     , std_dens_H(:,:),        1, 'y', i_real4, partit, mesh)
       call def_stream((/nl-1,       nod2D /),  (/nl-1,       myDim_nod2D /), 'density_dMOC',   'density'               , 'm',      density_dmoc(:,:),      1, 'y', i_real4, partit, mesh)
       call def_stream(elem2D,                                myDim_elem2D  , 'density_flux_e', 'density flux at elems ', 'm',      dens_flux_e(:),         1, 'y', i_real4, partit, mesh)
    end if
    
    !___________________________________________________________________________
    ! Discrete Variance Decay Dignostic, without predefined freq, freq_unit, 
    ! prec, --> default annual output
    if (ldiag_DVD .and. sel_dvd==0) then
        !___temperature DVD_____________________________________________________
        call def_stream((/nl-1, nod2D/), (/nl-1, myDim_nod2D/), 'dvd_temp_KK_tot' ,  'tot. temperature DVD \n (Klingbeil et al. 2014)', 'K^2/s' , dvd_KK_tot(        :,:,1), 1, 'y', i_real4, partit, mesh)
        call def_stream((/nl-1, nod2D/), (/nl-1, myDim_nod2D/), 'dvd_temp_SD_tot' ,  'tot. temperature DVD \n (Banerjee et al. 2023)' , 'K^2/s' , dvd_SD_tot(        :,:,1), 1, 'y', i_real4, partit, mesh)
        call def_stream((/nl-1, nod2D/), (/nl-1, myDim_nod2D/), 'dvd_temp_SD_advh',  'temperature DVD horiz. adv.'                    , 'K^2/s' , dvd_SD_chi_adv_h(  :,:,1), 1, 'y', i_real4, partit, mesh)
        call def_stream((/nl-1, nod2D/), (/nl-1, myDim_nod2D/), 'dvd_temp_SD_advv',  'temperature DVD vert. adv. '                    , 'K^2/s' , dvd_SD_chi_adv_v(  :,:,1), 1, 'y', i_real4, partit, mesh)
        if (Redi .eqv. .False. .and. K_hor /= 0.0_WP) then
            call def_stream((/nl-1, nod2D/), (/nl-1, myDim_nod2D/), 'dvd_temp_SD_difh',  'temperature DVD horiz. diff. expl.'             , 'K^2/s' , dvd_SD_chi_dif_he( :,:,1), 1, 'y', i_real4, partit, mesh)
        end if
        call def_stream((/nl-1, nod2D/), (/nl-1, myDim_nod2D/), 'dvd_temp_SD_difvi', 'temperature DVD vert. diff. impl.'              , 'K^2/s' , dvd_SD_chi_dif_vi( :,:,1), 1, 'y', i_real4, partit, mesh)
        call def_stream((/nl-1, nod2D/), (/nl-1, myDim_nod2D/), 'dvd_temp_SD_difsbc', 'temperature DVD sbc.'                          , 'K^2/s' , dvd_SD_chi_dif_sbc(:,:,1), 1, 'y', i_real4, partit, mesh)
        if (.not. tracers%data(1)%i_vert_diff) then
            call def_stream((/nl-1, nod2D/), (/nl-1, myDim_nod2D/), 'dvd_temp_SD_difve', 'temperature DVD vert. diff. expl.'          , 'K^2/s' , dvd_SD_chi_dif_ve( :,:,1), 1, 'y', i_real4, partit, mesh)
        end if
        if (tracers%data(1)%smooth_bh_tra) then 
            call def_stream((/nl-1, nod2D/), (/nl-1, myDim_nod2D/), 'dvd_temp_SD_difhbh','temperature DVD horiz. diff. biharm'        , 'K^2/s' , dvd_SD_chi_dif_hbh(:,:,1), 1, 'y', i_real4, partit, mesh)
        end if 
        if (Redi) then 
            call def_stream((/nl-1, nod2D/), (/nl-1, myDim_nod2D/), 'dvd_temp_SD_difheR', 'temperature DVD horiz. redi diff. expl.'   , 'K^2/s' , dvd_SD_chi_dif_heR(:,:,1), 1, 'y', i_real4, partit, mesh)
            
            call def_stream((/nl-1, nod2D/), (/nl-1, myDim_nod2D/), 'dvd_temp_SD_difveR', 'temperature DVD vert. redi diff. expl.'    , 'K^2/s' , dvd_SD_chi_dif_veR(:,:,1), 1, 'y', i_real4, partit, mesh)
            call def_stream((/nl-1, nod2D/), (/nl-1, myDim_nod2D/), 'dvd_temp_SD_difviR', 'temperature DVD vert. redi diff. impl.'    , 'K^2/s' , dvd_SD_chi_dif_viR(:,:,1), 1, 'y', i_real4, partit, mesh)
        end if     
        
        !___salinity DVD________________________________________________________
        call def_stream((/nl-1, nod2D/), (/nl-1, myDim_nod2D/), 'dvd_salt_KK_tot' ,  'tot. salinity DVD \n (Klingbeil et al. 2014)'   , 'K^2/s' , dvd_KK_tot(        :,:,2), 1, 'y', i_real4, partit, mesh)
        call def_stream((/nl-1, nod2D/), (/nl-1, myDim_nod2D/), 'dvd_salt_SD_tot' ,  'tot. salinity DVD \n (Banerjee et al. 2023)'    , 'K^2/s' , dvd_SD_tot(        :,:,2), 1, 'y', i_real4, partit, mesh)
        call def_stream((/nl-1, nod2D/), (/nl-1, myDim_nod2D/), 'dvd_salt_SD_advh',  'salinity DVD horiz. adv.'                       , 'K^2/s' , dvd_SD_chi_adv_h(  :,:,2), 1, 'y', i_real4, partit, mesh)
        call def_stream((/nl-1, nod2D/), (/nl-1, myDim_nod2D/), 'dvd_salt_SD_advv',  'salinity DVD vert. adv. '                       , 'K^2/s' , dvd_SD_chi_adv_v(  :,:,2), 1, 'y', i_real4, partit, mesh)
        if (Redi .eqv. .False. .and. K_hor /= 0.0_WP) then
            call def_stream((/nl-1, nod2D/), (/nl-1, myDim_nod2D/), 'dvd_salt_SD_difh',  'salinity DVD horiz. diff. expl.'                , 'K^2/s' , dvd_SD_chi_dif_he( :,:,2), 1, 'y', i_real4, partit, mesh)
        end if 
        call def_stream((/nl-1, nod2D/), (/nl-1, myDim_nod2D/), 'dvd_salt_SD_difvi', 'salinity DVD vert. diff. impl.'                 , 'K^2/s' , dvd_SD_chi_dif_vi( :,:,2), 1, 'y', i_real4, partit, mesh)
        call def_stream((/nl-1, nod2D/), (/nl-1, myDim_nod2D/), 'dvd_salt_SD_difsbc','salinity DVD sbc.'                              , 'K^2/s' , dvd_SD_chi_dif_sbc(:,:,2), 1, 'y', i_real4, partit, mesh)
        if (.not. tracers%data(1)%i_vert_diff) then
            call def_stream((/nl-1, nod2D/), (/nl-1, myDim_nod2D/), 'dvd_salt_SD_difve', 'salinity DVD vert. diff. expl.'             , 'K^2/s' , dvd_SD_chi_dif_ve( :,:,2), 1, 'y', i_real4, partit, mesh)
        end if
        if (tracers%data(1)%smooth_bh_tra) then 
            call def_stream((/nl-1, nod2D/), (/nl-1, myDim_nod2D/), 'dvd_salt_SD_difhbh','salinity DVD horiz. diff. biharm'           , 'K^2/s' , dvd_SD_chi_dif_hbh(:,:,2), 1, 'y', i_real4, partit, mesh)
        end if 
        if (Redi) then 
            call def_stream((/nl-1, nod2D/), (/nl-1, myDim_nod2D/), 'dvd_salt_SD_difheR', 'salinity DVD horiz. redi diff. expl.'      , 'K^2/s' , dvd_SD_chi_dif_heR(:,:,2), 1, 'y', i_real4, partit, mesh)
            
            call def_stream((/nl-1, nod2D/), (/nl-1, myDim_nod2D/), 'dvd_salt_SD_difveR', 'salinity DVD vert. redi diff. expl.'       , 'K^2/s' , dvd_SD_chi_dif_veR(:,:,2), 1, 'y', i_real4, partit, mesh)
            call def_stream((/nl-1, nod2D/), (/nl-1, myDim_nod2D/), 'dvd_salt_SD_difviR', 'salinity DVD vert. redi diff. impl.'       , 'K^2/s' , dvd_SD_chi_dif_viR(:,:,2), 1, 'y', i_real4, partit, mesh)
        end if 
    end if
    
    !___________________________________________________________________________
    ! output Redi parameterisation
#if !defined(__MULTIO)
    if (Redi) then
        if (sel_redi==0) then 
            call def_stream((/nl-1  , nod2D /), (/nl-1,   myDim_nod2D /), 'redi_K',   'Redi diffusion coefficient', 'm2/s', Ki(:,:),    1, 'y', i_real4, partit, mesh)
        endif     
    end if
#endif

    !___________________________________________________________________________
    ! output Monin-Obukov (TB04) mixing length
    !if (use_momix) then
    !    call def_stream(nod2D, myDim_nod2D, 'momix_length',   'Monin-Obukov mixing length', 'm', mixlength(:),    1, 'm', i_real4, partit, mesh)
    !end if
  
    !___________________________________________________________________________
    if (ldiag_curl_vel3) then
        call def_stream((/nl-1, nod2D/),  (/nl-1, myDim_nod2D/),  'curl_u',     'relative vorticity',          '1/s',   curl_vel3,                   1, 'm', i_real4, partit, mesh)
    end if

    !___________________________________________________________________________
    if (ice%whichEVP==1) then
    end if
    
    if (ice%whichEVP==2) then
        call def_stream(elem2D, myDim_elem2D, 'alpha_EVP', 'alpha in EVP', 'n/a', ice%alpha_evp_array,  1, 'd', i_real4, partit, mesh)
        call def_stream(nod2D,  myDim_nod2D,  'beta_EVP',  'beta in EVP',  'n/a', ice%beta_evp_array,   1, 'd', i_real4, partit, mesh)
    end if
    
    !___________________________________________________________________________
    if (dynamics%ldiag_ke) then
       call def_stream((/nl-1, elem2D/), (/nl-1, myDim_elem2D/), 'ke_adv_u_xVEL', 'work of advection   [u]',        'm2/s2', dynamics%ke_adv_xVEL(1,:,:),  io_list(i)%freq, 'm', 8, partit, mesh)
       call def_stream((/nl-1, elem2D/), (/nl-1, myDim_elem2D/), 'ke_adv_v_xVEL', 'work of advection   [v]',        'm2/s2', dynamics%ke_adv_xVEL(2,:,:),  io_list(i)%freq, 'm', 8, partit, mesh)
       !TODO: @sidorenko to clean up 
       !call def_stream((/nl-1, elem2D/), (/nl-1, myDim_elem2D/), 'ke_cor_u_xVEL', 'work Coriolis  :)  [u]',         'm2/s2', dynamics%ke_cor_xVEL(1,:,:),  io_list(i)%freq, 'm', 8, partit, mesh)
       !call def_stream((/nl-1, elem2D/), (/nl-1, myDim_elem2D/), 'ke_cor_v_xVEL', 'work Coriolis  :)  [v]',         'm2/s2', dynamics%ke_cor_xVEL(2,:,:),  io_list(i)%freq, 'm', 8, partit, mesh)

       !call def_stream((/nl-1, elem2D/), (/nl-1, myDim_elem2D/), 'ke_pre_u_xVEL', 'work of pressure gradient  [x]', 'm2/s2', dynamics%ke_pre_xVEL(1,:,:),  io_list(i)%freq, 'm', 8, partit, mesh)
       !call def_stream((/nl-1, elem2D/), (/nl-1, myDim_elem2D/), 'ke_pre_v_xVEL', 'work of pressure gradient  [y]', 'm2/s2', dynamics%ke_pre_xVEL(2,:,:),  io_list(i)%freq, 'm', 8, partit, mesh)

       !call def_stream((/nl-1, elem2D/), (/nl-1, myDim_elem2D/), 'ke_hvis_u_xVEL', 'work of hor. visc. [x]',        'm2/s2', dynamics%ke_hvis_xVEL(1,:,:), io_list(i)%freq, 'm', 8, partit, mesh)
       !call def_stream((/nl-1, elem2D/), (/nl-1, myDim_elem2D/), 'ke_hvis_v_xVEL', 'work of hor. visc. [y]',        'm2/s2', dynamics%ke_hvis_xVEL(2,:,:), io_list(i)%freq, 'm', 8, partit, mesh)

       !call def_stream((/nl-1, elem2D/), (/nl-1, myDim_elem2D/), 'ke_vvis_u_xVEL', 'work of ver. visc. [x]',        'm2/s2', dynamics%ke_vvis_xVEL(1,:,:), io_list(i)%freq, 'm', 8, partit, mesh)
       !call def_stream((/nl-1, elem2D/), (/nl-1, myDim_elem2D/), 'ke_vvis_v_xVEL', 'work of ver. visc. [y]',        'm2/s2', dynamics%ke_vvis_xVEL(2,:,:), io_list(i)%freq, 'm', 8, partit, mesh)

       !call def_stream((/nl-1, elem2D/), (/nl-1, myDim_elem2D/), 'ke_du2',    'KE change [0.5 du2]',           'm2/s2', dynamics%ke_du2(1,:,:),  io_list(i)%freq, 'm', 8, partit, mesh)
       !call def_stream((/nl-1, elem2D/), (/nl-1, myDim_elem2D/), 'ke_dv2',    'KE change [0.5 dv2]',           'm2/s2', dynamics%ke_du2(2,:,:),  io_list(i)%freq, 'm', 8, partit, mesh)
       
       !call def_stream((/nl-1, elem2D/), (/nl-1, myDim_elem2D/), 'ke_cor_u', 'Coriolis  *dT   [X]',       'm/s', dynamics%ke_cor(1,:,:),  io_list(i)%freq, 'm', 8, partit, mesh)
       !call def_stream((/nl-1, elem2D/), (/nl-1, myDim_elem2D/), 'ke_cor_v', 'Coriolis  *dT   [Y]',       'm/s', dynamics%ke_cor(2,:,:),  io_list(i)%freq, 'm', 8, partit, mesh)

       !call def_stream((/nl-1, elem2D/), (/nl-1, myDim_elem2D/), 'ke_pre_u', 'pressure gradient *dT [x]', 'm/s', dynamics%ke_pre(1,:,:),  io_list(i)%freq, 'm', 8, partit, mesh)
       !call def_stream((/nl-1, elem2D/), (/nl-1, myDim_elem2D/), 'ke_pre_v', 'pressure gradient *dT [y]', 'm/s', dynamics%ke_pre(2,:,:),  io_list(i)%freq, 'm', 8, partit, mesh)

       !call def_stream((/nl-1, elem2D/), (/nl-1, myDim_elem2D/), 'ke_hvis_u', 'hor. visc. [x] *dT',       'm/s', dynamics%ke_hvis(1,:,:), io_list(i)%freq, 'm', 8, partit, mesh)
       !call def_stream((/nl-1, elem2D/), (/nl-1, myDim_elem2D/), 'ke_hvis_v', 'hor. visc. [y] *dT',       'm/s', dynamics%ke_hvis(2,:,:), io_list(i)%freq, 'm', 8, partit, mesh)

       !call def_stream((/nl-1, elem2D/), (/nl-1, myDim_elem2D/), 'ke_vvis_u', 'ver. visc. [x] *dT',       'm/s', dynamics%ke_vvis(1,:,:), io_list(i)%freq, 'm', 8, partit, mesh)
       !call def_stream((/nl-1, elem2D/), (/nl-1, myDim_elem2D/), 'ke_vvis_v', 'ver. visc. [y] *dT',       'm/s', dynamics%ke_vvis(2,:,:), io_list(i)%freq, 'm', 8, partit, mesh)

       call def_stream((/nl-1, nod2D/),  (/nl-1, myDim_nod2D/),  'W_x_RHO',   'buoyancy work = ke_pre x VEL',  'm2/s2', dynamics%ke_wrho(:,:),   io_list(i)%freq, 'y', 8, partit, mesh)

       call def_stream(elem2D, myDim_elem2D,   'ke_wind_x_xVEL',  'work of wind [x]', 'm2/s2', dynamics%ke_wind_xVEL(1,:),  io_list(i)%freq, 'y', 8, partit, mesh)
       call def_stream(elem2D, myDim_elem2D,   'ke_wind_y_xVEL',  'work of wind [y]', 'm2/s2', dynamics%ke_wind_xVEL(2,:),  io_list(i)%freq, 'y', 8, partit, mesh)

       !call def_stream(elem2D, myDim_elem2D,   'ke_drag_x_xVEL',  'work of drag [x]', 'm2/s2', dynamics%ke_drag_xVEL(1,:),  io_list(i)%freq, 'm', 8, partit, mesh)
       !call def_stream(elem2D, myDim_elem2D,   'ke_drag_y_xVEL',  'work of drag [y]', 'm2/s2', dynamics%ke_drag_xVEL(2,:),  io_list(i)%freq, 'm', 8, partit, mesh)
!      same as above but without multiplying with UMEAN (for later computation of turbulence fluxes) 
       call def_stream((/nl-1, elem2D/), (/nl-1, myDim_elem2D/), 'ke_adv_u', 'advection *dT   [u]',       'm/s', dynamics%ke_adv(1,:,:),  io_list(i)%freq, 'y', 8, partit, mesh)
       call def_stream((/nl-1, elem2D/), (/nl-1, myDim_elem2D/), 'ke_adv_v', 'advection *dT   [v]',       'm/s', dynamics%ke_adv(2,:,:),  io_list(i)%freq, 'y', 8, partit, mesh)


       call def_stream((/nl-1, elem2D/), (/nl-1, myDim_elem2D/), 'ke_Umean',   'mean U',                  'm/s', dynamics%ke_umean(1,:,:),  io_list(i)%freq, 'y', 8, partit, mesh)
       call def_stream((/nl-1, elem2D/), (/nl-1, myDim_elem2D/), 'ke_Vmean',   'mean V',                  'm/s', dynamics%ke_umean(2,:,:),  io_list(i)%freq, 'y', 8, partit, mesh)

       call def_stream((/nl-1, elem2D/), (/nl-1, myDim_elem2D/), 'ke_U2mean',   'U2 mean',                'm/s', dynamics%ke_u2mean(1,:,:),  io_list(i)%freq, 'y', 8, partit, mesh)
       call def_stream((/nl-1, elem2D/), (/nl-1, myDim_elem2D/), 'ke_V2mean',   'V2 mean',                'm/s', dynamics%ke_u2mean(2,:,:),  io_list(i)%freq, 'y', 8, partit, mesh)

       call def_stream((/nl-1, nod2D/),  (/nl-1, myDim_nod2D/),  'ke_dW',      'd/dz (vertical VEL)',     'm/s', dynamics%ke_dW(:,:),       io_list(i)%freq, 'y', 8, partit, mesh)
       call def_stream((/nl-1, nod2D/),  (/nl-1, myDim_nod2D/),  'ke_PFULL',   'full Pressure',           'm/s', dynamics%ke_Pfull(:,:),    io_list(i)%freq, 'y', 8, partit, mesh)

       call def_stream(elem2D, myDim_elem2D,   'ke_wind_x',      'wind [x] *dT', 'm/s', dynamics%ke_wind(1,:),  io_list(i)%freq, 'y', 8, partit, mesh)
       call def_stream(elem2D, myDim_elem2D,   'ke_wind_y',      'wind [y] *dT', 'm/s', dynamics%ke_wind(2,:),  io_list(i)%freq, 'y', 8, partit, mesh)

       !call def_stream(elem2D, myDim_elem2D,   'ke_drag_x',      'drag [x] *dT', 'm/s', dynamics%ke_drag(1,:),  io_list(i)%freq, 'm', 8, partit, mesh)
       !call def_stream(elem2D, myDim_elem2D,   'ke_drag_y',      'drag [y] *dT', 'm/s', dynamics%ke_drag(2,:),  io_list(i)%freq, 'm', 8, partit, mesh)
       ! surface fields for APE input computations...
       call def_stream(nod2D , myDim_nod2D ,   'ke_J',           'surface temperature flux [Js]','°C/s',               dynamics%ke_J(:),   io_list(i)%freq, 'm', 8, partit, mesh)
       call def_stream(nod2D , myDim_nod2D ,   'ke_G',           'surface salinity    flux [Gs]','PSU/s',              dynamics%ke_G(:),   io_list(i)%freq, 'm', 8, partit, mesh)
       call def_stream(nod2D , myDim_nod2D ,   'ke_D',           'surface density',              'kg/m^3',             dynamics%ke_D(:),   io_list(i)%freq, 'm', 8, partit, mesh)
       call def_stream(nod2D , myDim_nod2D ,   'ke_D2',          'surface density squared',      'kg^2/m^6',           dynamics%ke_D2(:),  io_list(i)%freq, 'm', 8, partit, mesh)
       call def_stream(nod2D , myDim_nod2D ,   'ke_JD',          'surface temperature flux [Js] * RHO','°C*kg/s/m^3',  dynamics%ke_JD(:),  io_list(i)%freq, 'm', 8, partit, mesh)
       call def_stream(nod2D , myDim_nod2D ,   'ke_GD',          'surface salinity    flux [Gs] * RHO','PSU*kg/s/m^3', dynamics%ke_GD(:),  io_list(i)%freq, 'm', 8, partit, mesh)
       call def_stream(nod2D , myDim_nod2D ,   'ke_swA',         'Thermal expansion coeff  (alpha)',   '1/°C',         dynamics%ke_swA(:), io_list(i)%freq, 'm', 8, partit, mesh)
       call def_stream(nod2D , myDim_nod2D ,   'ke_swB',         'Taline contraction coeff (beta)',    '1/PSU',        dynamics%ke_swB(:), io_list(i)%freq, 'm', 8, partit, mesh)
       ! fields required to compute the energy conversions between Pm<->Pe, as well as Pm & Pe
       call def_stream((/nl-1, nod2D/),  (/nl-1, myDim_nod2D/) ,   'ke_n0',     'dRHO/dz',         'kg/m^4',   dynamics%ke_n0(:,:),     io_list(i)%freq, 'm', 8, partit, mesh)
       call def_stream((/nl-1, elem2D/), (/nl-1, myDim_elem2D/) ,  'ke_Dx',     'dRHO/dx',         'kg/m^4',   dynamics%ke_Dx(:,:),     io_list(i)%freq, 'm', 8, partit, mesh)
       call def_stream((/nl-1, elem2D/), (/nl-1, myDim_elem2D/) ,  'ke_Dy',     'dRHO/dy',         'kg/m^4',   dynamics%ke_Dy(:,:),     io_list(i)%freq, 'm', 8, partit, mesh)
       call def_stream((/nl-1, elem2D/), (/nl-1, myDim_elem2D/) ,  'ke_DU',     'RHO*U',           'kg*s/m^2', dynamics%ke_DU(:,:),     io_list(i)%freq, 'm', 8, partit, mesh)
       call def_stream((/nl-1, elem2D/), (/nl-1, myDim_elem2D/) ,  'ke_DV',     'RHO*V',           'kg*s/m^2', dynamics%ke_DV(:,:),     io_list(i)%freq, 'm', 8, partit, mesh)
       call def_stream((/nl-1, elem2D/), (/nl-1, myDim_elem2D/) ,  'ke_elemD',  'RHO*    on elem', 'kg/m^3',   dynamics%ke_elemD(:,:),  io_list(i)%freq, 'm', 8, partit, mesh)
       call def_stream((/nl-1, elem2D/), (/nl-1, myDim_elem2D/) ,  'ke_elemD2', 'RHO*^2  on elem', 'kg2/m^6',  dynamics%ke_elemD2(:,:), io_list(i)%freq, 'm', 8, partit, mesh)
    end if    
end subroutine
!
!
!_______________________________________________________________________________
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
!
!_______________________________________________________________________________
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
    use diagnostics, only: std_dens

    implicit none
    character(2000)               :: att_text
    type(t_mesh)  , intent(in) :: mesh
    type(t_partit), intent(in) :: partit
    type(t_dyn)   , intent(in) :: dynamics
    type(t_ice)   , intent(in) :: ice
    
    type(Meandata), intent(inout) :: entry
    character(len=*), parameter :: global_attributes_prefix = "FESOM_"    
#if defined(__icepack)
    integer, allocatable :: ncat_arr(:)
    integer              :: ii
#endif

    ! Serial output implemented so far
    if (partit%mype/=entry%root_rank) return
    ! create an ocean output file
    write(*,*) 'initializing I/O file for ', trim(entry%name)
    
    !___________________________________________________________________________
    ! Create file
    call assert_nf( nf90_create(entry%filename, IOR(nf90_noclobber,IOR(nf90_netcdf4,nf90_classic_model)), entry%ncid), __LINE__)

    !___________________________________________________________________________
    ! Create mesh related dimensions
    if (entry%ndim==1) then
        call assert_nf( nf90_def_dim(entry%ncid, entry%dimname(1), entry%glsize(2), entry%dimID(1)), __LINE__)
        
    else if (entry%ndim==2) then
        call assert_nf( nf90_def_dim(entry%ncid,  entry%dimname(1), entry%glsize(1), entry%dimID(1)), __LINE__)
        if (entry%dimname(1)=='nz') then
            call assert_nf( nf90_def_var(entry%ncid,  entry%dimname(1), nf90_double,   (/entry%dimID(1)/), entry%dimvarID(1)), __LINE__)
            call assert_nf( nf90_put_att(entry%ncid, entry%dimvarID(1), 'long_name', 'depth at layer interface'), __LINE__)
            
        elseif (entry%dimname(1)=='nz1') then
            call assert_nf( nf90_def_var(entry%ncid,  entry%dimname(1), nf90_double,   (/entry%dimID(1)/), entry%dimvarID(1)), __LINE__)
            call assert_nf( nf90_put_att(entry%ncid, entry%dimvarID(1), 'long_name', 'depth at layer midpoint'), __LINE__)
            
        elseif (entry%dimname(1)=='ncat') then
            call assert_nf( nf90_def_var(entry%ncid,  entry%dimname(1), nf90_int,   (/entry%dimID(1)/), entry%dimvarID(1)), __LINE__)
            call assert_nf( nf90_put_att(entry%ncid, entry%dimvarID(1), 'long_name', 'sea-ice thickness class'), __LINE__)
        
        elseif (entry%dimname(1)=='ndens') then
            call assert_nf( nf90_def_var(entry%ncid,  entry%dimname(1), nf90_int,   (/entry%dimID(1)/), entry%dimvarID(1)), __LINE__)
            call assert_nf( nf90_put_att(entry%ncid, entry%dimvarID(1), 'long_name', 'sigma2 density class'), __LINE__)
        
        else
            if (partit%mype==0) write(*,*) 'WARNING: unknown first dimension in 2d mean I/O data'
            
        end if 
        call assert_nf( nf90_put_att(entry%ncid, entry%dimvarID(1), 'units', 'm'), __LINE__)
        call assert_nf( nf90_put_att(entry%ncid, entry%dimvarID(1), 'positive', 'down'), __LINE__)
        call assert_nf( nf90_put_att(entry%ncid, entry%dimvarID(1), 'axis', 'Z'), __LINE__)
        
        call assert_nf( nf90_def_dim(entry%ncid, entry%dimname(2), entry%glsize(2), entry%dimID(2)), __LINE__)
    end if
  
    !___________________________________________________________________________
    ! Create time related dimensions
    call assert_nf( nf90_def_dim(entry%ncid, 'time', nf90_unlimited, entry%recID), __LINE__)
  
    !___________________________________________________________________________
    ! Define the time and iteration variables
    call assert_nf( nf90_def_var(entry%ncid, 'time', nf90_double, (/entry%recID/), entry%tID), __LINE__)
    att_text='time'
    call assert_nf( nf90_put_att(entry%ncid, entry%tID, 'long_name', trim(att_text)), __LINE__)
    call assert_nf( nf90_put_att(entry%ncid, entry%tID, 'standard_name', trim(att_text)), __LINE__)
    write(att_text, '(a14,I4.4,a1,I2.2,a1,I2.2,a6)') 'seconds since ', yearold, '-', 1, '-', 1, ' 0:0:0'
    call assert_nf( nf90_put_att(entry%ncid, entry%tID, 'units', trim(att_text)), __LINE__)
    call assert_nf( nf90_put_att(entry%ncid, entry%tID, 'axis', 'T'), __LINE__)
    call assert_nf( nf90_put_att(entry%ncid, entry%tID, 'stored_direction', 'increasing'), __LINE__)

    call assert_nf( nf90_def_var(entry%ncid, trim(entry%name), entry%data_strategy%netcdf_type(), (/entry%dimid(entry%ndim:1:-1), entry%recID/), entry%varID), __LINE__)

    call assert_nf( nf90_def_var_deflate(entry%ncid, entry%varID, 0, 1, compression_level), __LINE__)
    call assert_nf( nf90_put_att(entry%ncid, entry%varID, 'description', entry%description), __LINE__)
    call assert_nf( nf90_put_att(entry%ncid, entry%varID, 'long_name', entry%description), __LINE__)
    call assert_nf( nf90_put_att(entry%ncid, entry%varID, 'units', entry%units), __LINE__)
    call assert_nf( nf90_put_att(entry%ncid, entry%varID, 'location', entry%defined_on), __LINE__)
    call assert_nf( nf90_put_att(entry%ncid, entry%varID, 'mesh', entry%mesh), __LINE__)

    ! Add _FillValue attribute for missing/invalid data (CF-compliant)
    if (entry%accuracy == i_real8) then
        call assert_nf( nf90_put_att(entry%ncid, entry%varID, '_FillValue', NC_FILL_DOUBLE), __LINE__)
    elseif (entry%accuracy == i_real4) then
        call assert_nf( nf90_put_att(entry%ncid, entry%varID, '_FillValue', NC_FILL_FLOAT), __LINE__)
    end if


  !___Global attributes________  
    call assert_nf( nf90_put_att(entry%ncid, nf90_global, 'title', 'FESOM2 output'), __LINE__)
    call assert_nf( nf90_put_att(entry%ncid, nf90_global, global_attributes_prefix//'model', 'FESOM2'), __LINE__)
    call assert_nf( nf90_put_att(entry%ncid, nf90_global, global_attributes_prefix//'website', 'fesom.de'), __LINE__)

    call assert_nf( nf90_put_att(entry%ncid, nf90_global, global_attributes_prefix//'git_SHA', fesom_git_sha()), __LINE__)
    call assert_nf( nf90_put_att(entry%ncid, nf90_global, global_attributes_prefix//'MeshPath', trim(MeshPath)), __LINE__)
    call assert_nf( nf90_put_att(entry%ncid, nf90_global, global_attributes_prefix//'mesh_representative_checksum', mesh%representative_checksum), __LINE__)
    call assert_nf( nf90_put_att(entry%ncid, nf90_global, global_attributes_prefix//'ClimateDataPath', trim(ClimateDataPath)), __LINE__)
    call assert_nf( nf90_put_att(entry%ncid, nf90_global, global_attributes_prefix//'which_ALE', trim(which_ALE)), __LINE__)
    call assert_nf( nf90_put_att(entry%ncid, nf90_global, global_attributes_prefix//'mix_scheme', trim(mix_scheme)), __LINE__)
    ! Global attributes already set above

  ! call assert_nf( nf90_put_att(entry%ncid, nf90_global, global_attributes_prefix//'tra_adv_hor', trim(tra_adv_hor)), __LINE__)
  ! call assert_nf( nf90_put_att(entry%ncid, nf90_global, global_attributes_prefix//'tra_adv_ver', trim(tra_adv_ver)), __LINE__)
  ! call assert_nf( nf90_put_att(entry%ncid, nf90_global, global_attributes_prefix//'tra_adv_lim', trim(tra_adv_lim)), __LINE__)
 
    call assert_nf( nf90_put_att(entry%ncid, nf90_global, global_attributes_prefix//'use_partial_cell', merge(1, 0, use_partial_cell)), __LINE__)
    call assert_nf( nf90_put_att(entry%ncid, nf90_global, global_attributes_prefix//'force_rotation', merge(1, 0, force_rotation)), __LINE__)
    call assert_nf( nf90_put_att(entry%ncid, nf90_global, global_attributes_prefix//'include_fleapyear', merge(1, 0, include_fleapyear)), __LINE__)
    call assert_nf( nf90_put_att(entry%ncid, nf90_global, global_attributes_prefix//'use_floatice', merge(1, 0, use_floatice)), __LINE__)
    call assert_nf( nf90_put_att(entry%ncid, nf90_global, global_attributes_prefix//'whichEVP', ice%whichEVP), __LINE__)
    call assert_nf( nf90_put_att(entry%ncid, nf90_global, global_attributes_prefix//'evp_rheol_steps', ice%evp_rheol_steps), __LINE__)
    call assert_nf( nf90_put_att(entry%ncid, nf90_global, global_attributes_prefix//'opt_visc', dynamics%opt_visc), __LINE__)
    call assert_nf( nf90_put_att(entry%ncid, nf90_global, global_attributes_prefix//'use_wsplit', merge(1, 0, dynamics%use_wsplit)), __LINE__)
    ! use_partial_cell already set above
    call assert_nf( nf90_put_att(entry%ncid, nf90_global, global_attributes_prefix//'autorotate_back_to_geo', merge(1, 0, vec_autorotate)), __LINE__)
 
    !___________________________________________________________________________
    ! This ends definition part of the file, below filling in variables is possible
    call assert_nf( nf90_enddef(entry%ncid), __LINE__)
    if (entry%dimname(1)=='nz') then
        call assert_nf( nf90_put_var(entry%ncid, entry%dimvarID(1), abs(mesh%zbar)), __LINE__)
    elseif (entry%dimname(1)=='nz1') then
        call assert_nf( nf90_put_var(entry%ncid, entry%dimvarID(1), abs(mesh%Z)), __LINE__)
#if defined(__icepack)    
    elseif (entry%dimname(1)=='ncat') then
        allocate(ncat_arr(entry%glsize(1)))
        do ii= 1, entry%glsize(1)
            ncat_arr(ii)=ii
        end do        
        call assert_nf( nf90_put_var(entry%ncid, entry%dimvarID(1), ncat_arr), __LINE__)    
        deallocate(ncat_arr)
#endif
    elseif (entry%dimname(1)=='ndens') then
        call assert_nf( nf90_put_var(entry%ncid, entry%dimvarID(1), std_dens), __LINE__)
    else
        if (partit%mype==0) write(*,*) 'WARNING: unknown first dimension in 2d mean I/O data'
    end if 
    
    !___________________________________________________________________________
    call assert_nf( nf90_close(entry%ncid), __LINE__)
end subroutine
!
!
!_______________________________________________________________________________
subroutine assoc_ids(entry)
    implicit none
    type(Meandata), intent(inout) :: entry
    integer                       :: j

    write(*,*) 'associating mean I/O file ', trim(entry%filename)

    do j=1, entry%ndim
        call assert_nf( nf90_inq_dimid(entry%ncid, entry%dimname(j), entry%dimID(j)), __LINE__)
    end do
    !___Associate time related dimensions_______________________________________
    call assert_nf( nf90_inq_dimid(entry%ncid, 'time', entry%recID), __LINE__)
    call assert_nf( nf90_inquire_dimension(entry%ncid, entry%recID, len=entry%rec_count), __LINE__)
    !___Associate the time and iteration variables______________________________
    call assert_nf( nf90_inq_varid(entry%ncid, 'time', entry%tID), __LINE__)
    !___Associate physical variables____________________________________________
    call assert_nf( nf90_inq_varid(entry%ncid, entry%name, entry%varID), __LINE__)
end subroutine
!
!
!_______________________________________________________________________________
! collect local mean output data (entry%local_values_r8_copy) into global 2d 
! array (entry%aux_r8) and use the root_rank CPU/Task to write them into netcdf
! file. In case of 3data write horizontal slices level wise.
subroutine write_mean(entry, entry_index)
    use mod_mesh
    use io_gather_module
    implicit none
    type(Meandata), intent(inout) :: entry
    integer, intent(in) :: entry_index
    integer tag
    integer                       :: i, size1, size2, size_gen, size_lev, order
    integer                       :: c, lev
    real(kind=8)                  :: t0,t1
    integer mpierr

    ! Serial output implemented so far
    !___________________________________________________________________________
    ! write new time index ctime_copy to file --> expand time array in nc file
    if (entry%p_partit%mype==entry%root_rank) then
        write(*,*) 'writing mean record for ', trim(entry%name), '; rec. count = ', entry%rec_count
        call assert_nf( nf90_put_var(entry%ncid, entry%Tid, entry%ctime_copy, start = (/entry%rec_count/) ), __LINE__)
    end if
  
    !_______writing 2D and 3D fields____________________________________________
    size1=entry%glsize(1)
    size2=entry%glsize(2)
    tag = 2 ! we can use a fixed tag here as we have an individual communicator for each output field
    
    !___________writing 8 byte real_____________________________________________
    if (entry%accuracy == i_real8) then
        
        !_______________________________________________________________________
        ! allocate global 2d array in which local data are gathered
        if(entry%p_partit%mype==entry%root_rank) then
            if(.not. allocated(entry%aux_r8)) allocate(entry%aux_r8(size2))
        else
            if(.not. allocated(entry%aux_r8)) allocate(entry%aux_r8(1))
        end if
        
        !_______________________________________________________________________
        ! loop over vertical layers --> do gather 3d variables layerwise in 2d
        ! slices
        do lev=1, size1
            !___________________________________________________________________
            ! local output variables are gahtered in 2d shaped entry%aux_r8 
            ! either for vertices or elements
            if(.not. entry%is_elem_based) then
                call gather_nod2D (entry%local_values_r8_copy(lev,1:size(entry%local_values_r8_copy,dim=2)), entry%aux_r8, entry%root_rank, tag, entry%comm, entry%p_partit)
            else
                call gather_elem2D(entry%local_values_r8_copy(lev,1:size(entry%local_values_r8_copy,dim=2)), entry%aux_r8, entry%root_rank, tag, entry%comm, entry%p_partit)
            end if
            
            !___________________________________________________________________
            ! use root_rank CPU/Task to write 2d slice into netcdf file for 3d 
            ! variables into specific layer position lev
            if (entry%p_partit%mype==entry%root_rank) then
                if (entry%ndim==1) then
                    call assert_nf( nf90_put_var(entry%ncid, entry%varID, entry%aux_r8, start=(/1, entry%rec_count/), count=(/size2, 1/)), __LINE__)
                elseif (entry%ndim==2) then
                    call assert_nf( nf90_put_var(entry%ncid, entry%varID, entry%aux_r8, start=(/1, lev, entry%rec_count/), count=(/size2, 1, 1/)), __LINE__)
                end if
            end if
        end do ! --> do lev=1, size1

    !___________writing 4 byte real ____________________________________________ 
    else if (entry%accuracy == i_real4) then
    
        !_______________________________________________________________________
        ! allocate global 2d array in which local data are gathered
        if(entry%p_partit%mype==entry%root_rank) then
            if(.not. allocated(entry%aux_r4)) allocate(entry%aux_r4(size2))
        else
            if(.not. allocated(entry%aux_r4)) allocate(entry%aux_r4(1))
        end if
        
        !_______________________________________________________________________
        ! loop over vertical layers --> do gather 3d variables layerwise in 2d
        ! slices
        do lev=1, size1
            !PS if (entry%p_partit%mype==entry%root_rank) t0=MPI_Wtime()  
            !___________________________________________________________________
            ! local output variables are gahtered in 2d shaped entry%aux_r8 
            ! either for vertices or elements
            if(.not. entry%is_elem_based) then
                call gather_real4_nod2D (entry%local_values_r4_copy(lev,1:size(entry%local_values_r4_copy,dim=2)), entry%aux_r4, entry%root_rank, tag, entry%comm, entry%p_partit)
            else
                call gather_real4_elem2D(entry%local_values_r4_copy(lev,1:size(entry%local_values_r4_copy,dim=2)), entry%aux_r4, entry%root_rank, tag, entry%comm, entry%p_partit)
            end if
            
            !___________________________________________________________________
            ! use root_rank CPU/Task to write 2d slice into netcdf file for 3d 
            ! variables into specific layer position lev
            if (entry%p_partit%mype==entry%root_rank) then
                if (entry%ndim==1) then
                    call assert_nf( nf90_put_var(entry%ncid, entry%varID, entry%aux_r4, start=(/1, entry%rec_count/), count=(/size2, 1/)), __LINE__)
                    !PS t1=MPI_Wtime()  
                    !PS if (entry%p_partit%flag_debug)  print *, achar(27)//'[31m'//' -I/O-> after nf90_put_var'//achar(27)//'[0m', entry%p_partit%mype, t1-t0
                elseif (entry%ndim==2) then
                    call assert_nf( nf90_put_var(entry%ncid, entry%varID, entry%aux_r4, start=(/1, lev, entry%rec_count/), count=(/size2, 1, 1/)), __LINE__)
                    !PS t1=MPI_Wtime()  
                    !PS if (entry%p_partit%flag_debug)  print *, achar(27)//'[31m'//' -I/O-> after nf90_put_var'//achar(27)//'[0m', entry%p_partit%mype, lev, t1-t0
                end if
            end if
        end do ! --> do lev=1, size1
    end if ! --> if (entry%accuracy == i_real8) then
end subroutine
!
!
!_______________________________________________________________________________
subroutine update_means
    implicit none
    type(Meandata), pointer :: entry
    integer                 :: n
    integer                 :: I, J

    DO n=1, io_NSTREAMS
        entry=>io_stream(n)
        
        !_____________ compute in 8 byte accuracy ______________________________
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
            
        !_____________ compute in 4 byte accuracy ______________________________
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
        END IF ! --> IF (entry%accuracy == i_real8) then
        
        entry%addcounter=entry%addcounter+1
    END DO ! --> DO n=1, io_NSTREAMS
end subroutine
!
!
!_______________________________________________________________________________
! main output routine called at the end of each time step --> here is decided if 
! output event is triggered
subroutine output(istep, ice, dynamics, tracers, partit, mesh)
    use g_clock
    use mod_mesh
    USE MOD_PARTIT
    USE MOD_PARSUP
    use MOD_DYN
    use MOD_ICE
    use mod_tracer
    use io_gather_module
#if defined(__MULTIO)
    use iom
#endif
#if defined (__icepack)
    use icedrv_main,    only: ini_mean_icepack_io
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
    character(:), allocatable             :: filepath
    real(real64)                          :: rtime !timestamp of the record
#if defined(__MULTIO)
    logical       :: output_done
    logical       :: trigger_flush
#endif

ctime=timeold+(dayold-1.)*86400
    
    !___________________________________________________________________________
    if (lfirst) then
        ! define output streams-->dimension, variable, long_name, units, array, freq, unit, precision
        !PS if (partit%flag_debug .and. partit%mype==0)  print *, achar(27)//'[32m'//' -I/O-> call ini_mean_io'//achar(27)//'[0m'
        call ini_mean_io(ice, dynamics, tracers, partit, mesh)
        
#if defined (__icepack)
        call ini_mean_icepack_io(mesh) !icapack has its copy of p_partit => partit
#endif

        !PS if (partit%flag_debug .and. partit%mype==0)  print *, achar(27)//'[33m'//' -I/O-> call init_io_gather'//achar(27)//'[0m'
        call init_io_gather(partit)

    end if ! --> if (lfirst) then
    
    !___________________________________________________________________________
    !PS if (partit%flag_debug .and. partit%mype==0)  print *, achar(27)//'[33m'//' -I/O-> call update_means'//achar(27)//'[0m'  
    call update_means

#if defined(__MULTIO)
    output_done = .false.
#endif

    !___________________________________________________________________________
    ! loop over defined streams
    do n=1, io_NSTREAMS
        !_______________________________________________________________________
        ! make pointer for entry onto io_stream object
        entry=>io_stream(n)
!#if defined(__MULTIO)
!        call mio_write_nod(mio, entry)
!        lfirst=.false.
!        return
!#endif

        !_______________________________________________________________________
        !check whether output will be written based on event frequency
        do_output=.false.
        if (entry%freq_unit.eq.'y') then
            call annual_event(do_output, entry%freq)
        else if (entry%freq_unit == 'm') then 
            call monthly_event(do_output, entry%freq) 
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

#if defined(__MULTIO)
        output_done = output_done .or. do_output
#endif
        
        !_______________________________________________________________________
        ! if its time for output --> do_output==.true.
        if (do_output) then
            if (vec_autorotate) call io_r2g(n, partit, mesh) ! automatically detect if a vector field and rotate if makes sense!
#if !defined(__MULTIO)
            if(entry%thread_running) call entry%thread%join()
            entry%thread_running = .false.
            
            ! define filepath
            if (filesplit_freq=='m') then
                filepath = trim(ResultPath)//trim(entry%name)//'.'//trim(runid)//'.'//cyearnew//'_'//cmonth//'.nc'
            else
                filepath = trim(ResultPath)//trim(entry%name)//'.'//trim(runid)//'.'//cyearnew//'.nc'
            endif

            !___________________________________________________________________
            ! only root rank task does output
            if(partit%mype == entry%root_rank) then
                !_______________________________________________________________
                ! create new output file ?!
                if(filepath /= trim(entry%filename)) then
                    if("" /= trim(entry%filename)) call assert_nf(nf90_close(entry%ncid), __LINE__)   
                    entry%filename = filepath
                    !___________________________________________________________
                    ! use any existing file with this name or create a new one
                    if( nf90_open(entry%filename, nf90_write, entry%ncid) /= nf90_noerr ) then
                        !PS if (partit%flag_debug)  print *, achar(27)//'[33m'//' -I/O-> call create_new_file'//achar(27)//'[0m'  
                        call create_new_file(entry, ice, dynamics, partit, mesh)
                        
                        !PS if (partit%flag_debug)  print *, achar(27)//'[33m'//' -I/O-> call assert_nf A'//achar(27)//'[0m'//',  k=',k, ', rootpart=', entry%root_rank    
                        call assert_nf( nf90_open(entry%filename, nf90_write, entry%ncid), __LINE__)
                    end if
                    
                    !___________________________________________________________
                    ! setup all dimension definition and attributes of the netcdf
                    ! file 
                    !PS if (partit%flag_debug)  print *, achar(27)//'[33m'//' -I/O-> call assoc_ids'//achar(27)//'[0m'  
                    call assoc_ids(entry)
                    
                end if ! --> if(filepath /= trim(entry%filename)) then
                
                !_______________________________________________________________
                ! if the time rtime at the rec_count is larger than ctime we 
                ! look for the closest record with the timestamp less than ctime
                do k=entry%rec_count, 1, -1
                    !PS if (partit%flag_debug)  print *, achar(27)//'[33m'//' -I/O-> call assert_nf B'//achar(27)//'[0m'//',  k=',k, ', rootpart=', entry%root_rank  
                    ! determine rtime from exiting file
                    call assert_nf( nf90_get_var(entry%ncid, entry%tID, rtime), __LINE__)
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
            end if ! --> if(partit%mype == entry%root_rank) then
#endif
            !___________________________________________________________________
            ! write double precision output
            if (entry%accuracy == i_real8) then
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I,J)
                DO J=1, size(entry%local_values_r8,dim=2)
                    DO I=1, size(entry%local_values_r8,dim=1)
                        ! Check if point has valid data (non-zero accumulated value)
                        ! Use small epsilon to account for floating point precision
                        if (abs(entry%local_values_r8(I,J)) < 1.0e-30_real64) then
                            entry%local_values_r8_copy(I,J) = NC_FILL_DOUBLE  ! No data - set to fill value
                        else
                            entry%local_values_r8_copy(I,J) = entry%local_values_r8(I,J) /real(entry%addcounter,real64)  ! compute_means
                        end if
                        entry%local_values_r8(I,J) = 0._real64 ! clean_meanarrays - reset to 0 for next accumulation
                    END DO ! --> DO I=1, size(entry%local_values_r8,dim=1)
                END DO ! --> DO J=1, size(entry%local_values_r8,dim=2)
!$OMP END PARALLEL DO
                
            !___________________________________________________________________
            ! write single precision output
            else if (entry%accuracy == i_real4) then
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I,J)
                DO J=1, size(entry%local_values_r4,dim=2)
                    DO I=1, size(entry%local_values_r4,dim=1)
                        ! Check if point has valid data (non-zero accumulated value)
                        ! Use small epsilon to account for floating point precision
                        if (abs(entry%local_values_r4(I,J)) < 1.0e-30_real32) then
                            entry%local_values_r4_copy(I,J) = NC_FILL_FLOAT  ! No data - set to fill value
                        else
                            entry%local_values_r4_copy(I,J) = entry%local_values_r4(I,J) /real(entry%addcounter,real32)  ! compute_means
                        end if
                        entry%local_values_r4(I,J) = 0._real32 ! clean_meanarrays - reset to 0 for next accumulation
                    END DO ! --> DO I=1, size(entry%local_values_r4,dim=1)
                END DO ! --> DO J=1, size(entry%local_values_r4,dim=2)
!$OMP END PARALLEL DO
            end if ! --> if (entry%accuracy == i_real8) then
            !___________________________________________________________________
            entry%lastcounter=entry%addcounter
            entry%addcounter   = 0  ! clean_meanarrays
            entry%ctime_copy = ctime

#if defined(__MULTIO)
!            if (n==1) then
            entry%rec_count = istep
            call send_data_to_multio(entry)
!            end if            
#else
            !___________________________________________________________________
            ! this is where the magic happens --> here do_output_callback is
            ! triggered as a method of the io_stream object --> call write_mean(...)
            call entry%thread%run()
            entry%thread_running = .true.
#endif
        endif ! --> if (do_output) then
    end do ! --> do n=1, io_NSTREAMS
    lfirst=.false.

#if defined(__MULTIO)
    if (output_done) then
        call iom_flush('N grid', istep)
    end if
#endif

end subroutine
!
!
!_______________________________________________________________________________
! this becomes the output callback functions that becomes part of the IOstream(idx).thread 
! object
! --> this callback function becomes initialised in def_stream_after_dimension_specific(...)
!     by call IOstream(idx)%thread%initialize(do_output_callback, entry_index)
! --> the execution of this function can be triggered by IOstream(idx)%thread%run()
! --> do_output_callback becomes an executable method of the object 
!     IOstream(idx)%thread
subroutine do_output_callback(entry_index)
    use mod_mesh
    USE MOD_PARTIT
    USE MOD_PARSUP
    integer, intent(in) :: entry_index ! index of variable in def_stream array 
    ! EO args
    type(Meandata), pointer :: entry

    entry=>io_stream(entry_index)
    entry%p_partit%mype=entry%mype_workaround ! for the thread callback, copy back the value of our mype as a workaround for errors with the cray envinronment (at least with ftn 2.5.9 and cray-mpich 7.5.3)
    !___________________________________________________________________________
    ! collect local mean output data (entry%local_values_r8_copy) into global 2d 
    ! array (entry%aux_r8) and use the root_rank CPU/Task to write them into netcdf
    ! file. In case of 3data write horizontal slices level wise.
    call write_mean(entry, entry_index)
  
    !___________________________________________________________________________
    ! The function NF_SYNC offers a way to synchronize the disk copy of a netCDF 
    ! dataset with in-memory buffers. There are two reasons you might want to 
    ! synchronize after writes:
    ! To minimize data loss in case of abnormal termination, or To make data 
    ! available to other processes for reading immediately after it is written. 
    if(entry%p_partit%mype == entry%root_rank) then 
        !PS if (entry%p_partit%flag_debug)  print *, achar(27)//'[31m'//' -I/O-> call nf_sync'//achar(27)//'[0m', entry%p_partit%mype
        call assert_nf( nf90_sync(entry%ncid), __LINE__ ) ! flush the file to disk after each write
    end if   
    
end subroutine
!
!
!_______________________________________________________________________________
! Why this is done?! --> executed in fesom_main.F90 in call fesom_finalize() 
! after the run lop is finished --> stop/cleanup threads?
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
!

!_______________________________________________________________________________
! build 3d meandata streaming object
subroutine def_stream3D(glsize, lcsize, name, description, units, data, freq, freq_unit, accuracy, partit, mesh, flip_array, long_description)
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
  character(len=*), optional, intent(in) :: long_description
  integer i

 
    !___________________________________________________________________________
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

    !___________________________________________________________________________
    if (partit%mype==0) then
        write(*,*) 'adding I/O stream 3D for ', trim(name)
    end if
    
    !___________________________________________________________________________
    ! initialise meandata streaming object
    call associate_new_stream(name, entry)
    entry%previousDate=-1
    entry%previousTime=-1
    entry%currentDate=yearold * 10000 + month * 100 + day_in_month
    entry%currentTime=INT(INT(timeold / 3600) * 10000 + (INT(timeold / 60) - INT(timeold / 3600) * 60) * 100 + (timeold-INT(timeold / 60) * 60))
    entry%startDate=entry%currentDate
    entry%startTime=entry%currentTime
    !___________________________________________________________________________
    ! fill up 3d meandata streaming object
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
    call def_stream_after_dimension_specific(entry, name, description, units, freq, freq_unit, accuracy, partit, mesh, long_description)

end subroutine
!
!
!_______________________________________________________________________________
! build 2d meandata streaming object
subroutine def_stream2D(glsize, lcsize, name, description, units, data, freq, freq_unit, accuracy, partit, mesh, long_description)
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
  character(len=*), optional, intent(in) :: long_description
  integer i

    !___________________________________________________________________________
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

    !___________________________________________________________________________
    if (partit%mype==0) then
        write(*,*) 'adding I/O stream 2D for ', trim(name)
    end if

    !___________________________________________________________________________
    ! initialise meandata streaming object
    call associate_new_stream(name, entry)
    entry%previousDate=-1
    entry%previousTime=-1
    entry%currentDate=yearold * 10000 + month * 100 + day_in_month
    entry%currentTime=INT(INT(timeold / 3600) * 10000 + (INT(timeold / 60) - INT(timeold / 3600) * 60) * 100 + (timeold-INT(timeold / 60) * 60))
    entry%startDate=entry%currentDate
    entry%startTime=entry%currentTime
    !___________________________________________________________________________
    ! fill up 3d meandata streaming object
    ! 2d specific
    entry%ptr3(1:1,1:size(data)) => data(:)

    if (accuracy == i_real8) then
        allocate(entry%local_values_r8(1, lcsize))
        entry%local_values_r8 = 0._real64
    elseif (accuracy == i_real4) then
        allocate(entry%local_values_r4(1, lcsize))
        entry%local_values_r4 = 0._real32
    end if

    ! non dimension specific
    entry%ndim=1
    entry%glsize=(/1, glsize/)
    entry%dimname(1)=mesh_dimname_from_dimsize(glsize, partit, mesh)
    entry%dimname(2)='unknown'

    ! non dimension specific
    call def_stream_after_dimension_specific(entry, name, description, units, freq, freq_unit, accuracy, partit, mesh, long_description)
end subroutine
!
!
!_______________________________________________________________________________
! initialse new meandata streaming object
subroutine associate_new_stream(name, entry)
    type(Meandata), pointer      :: entry
    character(len=*), intent(in) :: name
    integer i

    entry => null()
    
    !___________________________________________________________________________
    ! check if we already have this variable
    do i=1, io_NSTREAMS
        if(trim(io_stream(i)%name) .eq. name) then
            print *,"variable '"//name//"' already exists, &
                &check if you define it multiple times, for example in namelist.io, &
                &namelist.icepack, io_meandata.F90 or other place that add I/O stream."
            call assert(.false., __LINE__) 
        end if
    end do
    
    !___________________________________________________________________________
    ! add this instance to io_stream array
    io_NSTREAMS = io_NSTREAMS +1
    call assert(size(io_stream) >= io_NSTREAMS, __LINE__)
    entry=>io_stream(io_NSTREAMS)
end subroutine
!
!
!_______________________________________________________________________________
! further fill up 2d/3d meandata streaming object --> link output callback routine
! as stream object method
subroutine def_stream_after_dimension_specific(entry, name, description, units, freq, freq_unit, accuracy, partit, mesh, long_description)
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
    character(len=*), intent(in) :: long_description
    ! EO args
    logical async_netcdf_allowed
    integer provided_mpi_thread_support_level
    integer entry_index
    integer err
    
    entry_index = io_NSTREAMS
    
    !___________________________________________________________________________
    entry%accuracy = accuracy

    if (accuracy == i_real8) then
      allocate(data_strategy_nf_double_type :: entry%data_strategy)
    elseif (accuracy == i_real4) then
      allocate(data_strategy_nf_float_type :: entry%data_strategy)
    else
       if (partit%mype==0) write(*,*) 'not supported output accuracy:',accuracy,'for',trim(name)
       call par_ex(partit%MPI_COMM_FESOM, partit%mype, 1)
       stop
    endif ! accuracy

    !___________________________________________________________________________
    entry%name = name
    entry%description = description
    entry%long_description = long_description
    entry%mesh = "fesom_mesh";
    entry%units = units
    entry%filename = ""

    entry%freq=freq
    entry%freq_unit=freq_unit
    entry%addcounter   = 0
    entry%is_in_use=.true.

    !___________________________________________________________________________
    if(entry%glsize(1)==mesh%nod2D  .or. entry%glsize(2)==mesh%nod2D) then
      entry%is_elem_based = .false.
      entry%defined_on = "node"
      entry%shrinked_size=partit%myDim_nod2D
    else if(entry%glsize(1)==mesh%elem2D .or. entry%glsize(2)==mesh%elem2D) then
      entry%is_elem_based = .true.
      entry%defined_on = "face"
      entry%shrinked_size=partit%myDim_elem2D_shrinked
      allocate(entry%shrinked_indx(entry%shrinked_size))
      entry%shrinked_indx=partit%myInd_elem2D_shrinked
!      write(*,*) partit%mype, partit%myDim_elem2D, partit%myDim_elem2D_shrinked, partit%myDim_elem2D-partit%myDim_elem2D_shrinked
!      entry_index=0
!      call MPI_AllREDUCE(partit%myDim_elem2D_shrinked, entry_index, 1, MPI_INTEGER, MPI_SUM, partit%MPI_COMM_FESOM, err)
!      write(*,*) 'total elem=', mesh%elem2D, entry_index
    else
      if(partit%mype == 0) print *,"can not determine if ",trim(name)," is node or elem based"
      call par_ex(partit%MPI_COMM_FESOM, partit%mype, 1)
      stop
    end if
    
    !___________________________________________________________________________
    if (accuracy == i_real8) then
      allocate(entry%local_values_r8_copy(size(entry%local_values_r8, dim=1), size(entry%local_values_r8, dim=2)))
    else if (accuracy == i_real4) then
      allocate(entry%local_values_r4_copy(size(entry%local_values_r4, dim=1), size(entry%local_values_r4, dim=2)))
    end if

    !___________________________________________________________________________
    ! set up async output
    entry%root_rank = next_io_rank(partit%MPI_COMM_FESOM, async_netcdf_allowed, partit)

    call MPI_Comm_dup(partit%MPI_COMM_FESOM, entry%comm, err)

    ! initialise output callback routine as meandata stream object method
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
!
!
!_______________________________________________________________________________
subroutine assert_nf(status, line)
    integer, intent(in) :: status
    integer, intent(in) :: line
    ! EO args
    if(status /= nf90_noerr) then
        print *, "error in line ",line, __FILE__, ' ', trim(nf90_strerror(status))
        stop 1
    end if   
end subroutine
!
!
!_______________________________________________________________________________
subroutine assert(val, line)
    logical, intent(in) :: val
    integer, intent(in) :: line
    ! EO args
    if(.NOT. val) then
        print *, "error in line ",line, __FILE__
        stop 1
    end if
end subroutine
!
!
!_______________________________________________________________________________
! do vector rotation on the fly 
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
    IF ((trim(entry_x%name)=='utemp'   ) .AND. ((trim(entry_y%name)=='vtemp'   ))) do_rotation=.TRUE.
    IF ((trim(entry_x%name)=='usalt'   ) .AND. ((trim(entry_y%name)=='vsalt'   ))) do_rotation=.TRUE.

    IF (.NOT. (do_rotation)) RETURN
   
    IF (partit%mype==0) THEN
       write(*,*) trim(entry_x%name)//' and '//trim(entry_y%name)//' will be rotated before output!'
    END IF

    !___________________________________________________________________________
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

    !___________________________________________________________________________
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

#if defined(__MULTIO)
SUBROUTINE send_data_to_multio(entry)
    USE iom
    IMPLICIT NONE

    TYPE(Meandata), TARGET, INTENT(INOUT)           :: entry
    TYPE(iom_field_request)                         :: request
    INTEGER                                         :: numLevels, globalSize, lev, i, n

    numLevels = entry%glsize(1)
    globalSize = entry%glsize(2)

    request%name = trim(entry%name)
    entry%previousDate=entry%currentDate
    entry%previousTime=entry%currentTime
    entry%currentDate=yearnew * 10000 + month * 100 + day_in_month
    entry%currentTime=INT(INT(timenew / 3600) * 10000 + (INT(timenew / 60) - INT(timenew / 3600) * 60) * 100 + (timenew-INT(timenew / 60) * 60))

    request%previousDate=entry%previousDate
    request%previousTime=entry%previousTime
    request%currentDate =entry%currentDate
    request%currentTime =entry%currentTime
    request%startDate   =entry%startDate
    request%startTime   =entry%startTime
    request%lastcounter =entry%lastcounter
    request%sampleInterval=INT(dt)

    IF (.NOT. entry%is_elem_based) THEN
        request%gridType = "N grid"
    ELSE
        request%gridType = "C grid"
    END IF
    request%globalSize = globalSize
    request%step = entry%rec_count
    if (numLevels==1) then
        request%category="ocean-2d"
    else
        request%category="ocean-3d"
    end if
    ! loop over vertical layers --> do gather 3d variables layerwise in 2d slices
    DO lev=1, numLevels
        request%level = lev

        IF (entry%is_elem_based) THEN
            n = SIZE(entry%shrinked_indx)
        ELSE
            n = entry%shrinked_size
        END IF

        IF (ALLOCATED(multio_temporary_array) .AND. (SIZE(multio_temporary_array) .LT. n)) THEN
            DEALLOCATE(multio_temporary_array)
        ENDIF

        IF (.NOT. ALLOCATED(multio_temporary_array)) THEN
            ALLOCATE(multio_temporary_array(n))
        ENDIF

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
        DO i = 1, n
            IF (entry%is_elem_based) THEN
                IF (entry%accuracy == i_real8) THEN
                    multio_temporary_array(i) = entry%local_values_r8_copy(lev, entry%shrinked_indx(i))
                ELSE
                    multio_temporary_array(i) = entry%local_values_r4_copy(lev, entry%shrinked_indx(i))
                ENDIF
            ELSE
                IF (entry%accuracy == i_real8) THEN
                    multio_temporary_array(i) = entry%local_values_r8_copy(lev, i)
                ELSE
                    multio_temporary_array(i) = entry%local_values_r4_copy(lev, i)
                ENDIF
            ENDIF
        ENDDO
!$OMP END PARALLEL DO

        request%values(1:n) => multio_temporary_array(1:n)

        CALL iom_send_fesom_data(request)
    END DO
END SUBROUTINE
#endif
end module
