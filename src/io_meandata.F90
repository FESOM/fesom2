module io_MEANDATA

  use g_config
  use o_mesh
  use g_parsup
  use g_clock
  use g_comm_auto
  use o_ARRAYS
  use g_forcing_arrays
  use i_ARRAYS
  use o_mixing_KPP_mod
  use diagnostics
  use i_PARAM, only: whichEVP

  use, intrinsic :: ISO_FORTRAN_ENV

  implicit none
#include "netcdf.inc"
  private
  public :: def_stream3D, output
!
!--------------------------------------------------------------------------------------------
!
  integer, parameter  :: i_real8=8, i_real4=4, i_int2=2


  type Meandata
    private
    integer                                            :: ndim
    integer                                            :: lcsize(2)
    integer                                            :: glsize(2)
    integer                                            :: accuracy
    real(real64),  public, allocatable, dimension(:,:) :: local_values_r8
    real(real32),  public, allocatable, dimension(:,:) :: local_values_r4
    integer(int16),  public, allocatable, dimension(:,:) :: local_values_i2
    integer                                            :: addcounter=0
    real(kind=WP), pointer                             :: ptr2(:), ptr3(:,:)
    real(real32)                                       :: min_value        ! lower and upper bound, used to compute scale_factor, add_offset 
    real(real32)                                       :: max_value        !      keep for check if real life values remain in this interval
    real(real32)                                       :: scale_factor=1.  ! for netcdf4 conversion real4 <-> int2:
    real(real32)                                       :: add_offset=0.    !      real4_value = int2_value * scale_factor + add_offset
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
    integer                                            :: error_status(1000), error_count
    logical                                            :: is_in_use=.false.
  end type
!
!--------------------------------------------------------------------------------------------
!
  type(Meandata), save, allocatable, target :: io_stream(:)
  integer, save                             :: io_NSTREAMS=0
  real(kind=WP)                             :: ctime !current time in seconds from the beginning of the year
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
subroutine ini_mean_io
  implicit none
  integer           :: i
  character(len=10) :: id_string

!3D
  if (ldiag_energy) then
     call def_stream((/nl,   nod2D/),  (/nl,     myDim_nod2D/), 'rhof',     'in-situ density at faces',    'kg/m3',     rhof(:,:),  1, 'm', i_real8)
     call def_stream((/nl,   nod2D/),  (/nl,     myDim_nod2D/), 'wrhof',    'vertical velocity x density', 'kg/(s*m2)', wrhof(:,:), 1, 'm', i_real8)
     call def_stream((/nl-1, nod2D/),  (/nl-1,   myDim_nod2D/), 'uu',   'u times u', 'm2/s2', u_x_u(:,:), 1, 'm', i_real8)
     call def_stream((/nl-1, nod2D/),  (/nl-1,   myDim_nod2D/), 'uv',   'u times v', 'm2/s2', u_x_v(:,:), 1, 'm', i_real8)
     call def_stream((/nl-1, nod2D/),  (/nl-1,   myDim_nod2D/), 'vv',   'v times v', 'm2/s2', v_x_v(:,:), 1, 'm', i_real8)
     call def_stream((/nl,  elem2D/),  (/nl-1,   myDim_elem2D/),'uw',   'u times w', 'm2/s2', u_x_w(:,:), 1, 'm', i_real8)
     call def_stream((/nl,  elem2D/),  (/nl-1,   myDim_elem2D/),'vw',   'v times w', 'm2/s2', v_x_w(:,:), 1, 'm', i_real8)
     call def_stream((/nl-1, nod2D/),  (/nl-1,   myDim_nod2D/), 'dudx', 'du/dx',     '1/s',   dudx(:,:),  1, 'm', i_real8)
     call def_stream((/nl-1, nod2D/),  (/nl-1,   myDim_nod2D/), 'dudy', 'du/dy',     '1/s',   dudy(:,:),  1, 'm', i_real8)
     call def_stream((/nl-1, nod2D/),  (/nl-1,   myDim_nod2D/), 'dvdx', 'dv/dx',     '1/s',   dvdx(:,:),  1, 'm', i_real8)
     call def_stream((/nl-1, nod2D/),  (/nl-1,   myDim_nod2D/), 'dvdy', 'dv/dy',     '1/s',   dvdy(:,:),  1, 'm', i_real8)
     call def_stream((/nl,  elem2D/),  (/nl,     myDim_elem2D/), 'dudz', 'du/dz',     '1/s',   dudz(:,:),  1, 'm', i_real8)
     call def_stream((/nl,  elem2D/),  (/nl,     myDim_elem2D/), 'dvdz', 'dv/dz',     '1/s',   dvdz(:,:),  1, 'm', i_real8)
     call def_stream((/nl,  elem2D/),  (/nl,     myDim_elem2D/), 'av_dudz', 'int(Av * du/dz)',        'm3/s2',   av_dudz(:,:),    1, 'm', i_real4)
     call def_stream((/nl,  elem2D/),  (/nl,     myDim_elem2D/), 'av_dvdz', 'int(Av * dv/dz)',        'm3/s2',   av_dvdz(:,:),    1, 'm', i_real4)
     call def_stream((/nl,  elem2D/),  (/nl,     myDim_elem2D/), 'av_dudz_sq',  'Av * (du/dz)^2',     'm^2/s^3', av_dudz_sq(:,:), 1, 'm', i_real4)
     call def_stream((/nl,   elem2D/), (/nl,     myDim_elem2D/), 'Av',    'Vertical mixing A',         'm2/s',            Av(:,:), 1, 'm', i_real4)
     call def_stream((/nl-1, nod2D/),  (/nl-1,   myDim_nod2D/),  'unod',  'horizontal velocity at nodes', 'm/s', Unode(1,:,:), 1, 'm', i_real8)
     call def_stream((/nl-1, nod2D/),  (/nl-1,   myDim_nod2D/),  'vnod',  'meridional velocity at nodes', 'm/s', Unode(2,:,:), 1, 'm', i_real8)
    
    call def_stream((/nl-1, elem2D/), (/nl-1,   myDim_elem2D/), 'u',  'horizontal velocity', 'm/s', uv(1,:,:),     1, 'm', i_real4)
    call def_stream((/nl-1, elem2D/), (/nl-1,   myDim_elem2D/), 'v',  'meridional velocity', 'm/s', uv(2,:,:),     1, 'm', i_real4)
    call def_stream((/nl, nod2D/),    (/nl,     myDim_nod2D/),  'w',  'vertical velocity',   'm/s', Wvel(:,:),     1, 'm', i_real8)

    call def_stream(elem2D, myDim_elem2D,   'utau_surf',  '(u, tau) at the surface',      'N/(m s)', utau_surf(1:myDim_elem2D),     1, 'm', i_real4)
    call def_stream(elem2D, myDim_elem2D,   'utau_bott',  '(u, tau) at the bottom',       'N/(m s)', utau_bott(1:myDim_elem2D),     1, 'm', i_real4)
    call def_stream(elem2D, myDim_elem2D,   'u_bott',     'bottom velocity',                  'm/s', u_bott(1:myDim_elem2D),        1, 'm', i_real4)
    call def_stream(elem2D, myDim_elem2D,   'v_bott',     'bottom velocity',                  'm/s', v_bott(1:myDim_elem2D),        1, 'm', i_real4)
    call def_stream(elem2D, myDim_elem2D,   'u_surf',     'surface velocity',                 'm/s', u_surf(1:myDim_elem2D),        1, 'm', i_real4)
    call def_stream(elem2D, myDim_elem2D,   'v_surf',     'surface velocity',                 'm/s', u_surf(1:myDim_elem2D),        1, 'm', i_real4)
    call def_stream(elem2D, myDim_elem2D,   'tx_bot',     'bottom stress x',                 'N/m2', stress_bott(1, 1:myDim_elem2D),1, 'm', i_real4)
    call def_stream(elem2D, myDim_elem2D,   'ty_bot',     'bottom stress y',                 'N/m2', stress_bott(2, 1:myDim_elem2D),1, 'm', i_real4)
    call def_stream(elem2D, myDim_elem2D,   'tx_sur',     'zonal wind stress to ocean',      'm/s2', stress_surf(1, 1:myDim_elem2D),1, 'm', i_real4)
    call def_stream(elem2D, myDim_elem2D,   'ty_sur',     'meridional wind stress to ocean', 'm/s2', stress_surf(2, 1:myDim_elem2D),1, 'm', i_real4)

  else
     call def_stream((/nl-1, elem2D/), (/nl-1, myDim_elem2D/), 'u',    'horizontal velocity', 'm/s', uv(1,:,:),     1, 'y', i_real4)
     call def_stream((/nl-1, elem2D/), (/nl-1, myDim_elem2D/), 'v',    'meridional velocity', 'm/s', uv(2,:,:),     1, 'y', i_real4)
     call def_stream((/nl, nod2D/),    (/nl,   myDim_nod2D/),  'w',    'vertical velocity',   'm/s', Wvel(:,:),     1, 'y', i_real4)
     call def_stream((/nl,   elem2D/), (/nl,   myDim_elem2D/), 'Av',   'Vertical mixing A',   'm2/s',Av(:,:),       1, 'y', i_real4)
     call def_stream(elem2D, myDim_elem2D,   'tx_sur',     'zonal wind stress to ocean',      'm/s2', stress_surf(1, 1:myDim_elem2D),1, 'm', i_real4)
     call def_stream(elem2D, myDim_elem2D,   'ty_sur',     'meridional wind stress to ocean', 'm/s2', stress_surf(2, 1:myDim_elem2D),1, 'm', i_real4)

  end if

  call def_stream((/nl-1, nod2D/),  (/nl-1, myDim_nod2D/),  'temp', 'temperature',         'C',   tr_arr(:,:,1), 1, 'm', i_real4)
  call def_stream((/nl-1, nod2D/),  (/nl-1, myDim_nod2D/),  'salt', 'salinity',            'psu', tr_arr(:,:,2), 1, 'm', i_real4)
  
  do i=3, num_tracers
     write (id_string, "(I3.3)") tracer_id(i)
     call def_stream((/nl-1, nod2D/),  (/nl-1, myDim_nod2D/),  'tra_'//id_string, 'pasive tracer ID='//id_string, 'n/a', tr_arr(:,:,i), 1, 'm', i_real4)
  end do

  call def_stream((/nl-1, nod2D/),  (/nl-1, myDim_nod2D/),  'slope_x','neutral slope X',   'none',slope_tapered(1,:,:), 1, 'y', i_real4)
  call def_stream((/nl-1, nod2D/),  (/nl-1, myDim_nod2D/),  'slope_y','neutral slope Y',   'none',slope_tapered(2,:,:), 1, 'y', i_real4)
  call def_stream((/nl-1, nod2D/),  (/nl-1, myDim_nod2D/),  'slope_z','neutral slope Z',   'none',slope_tapered(3,:,:), 1, 'y', i_real4)
  call def_stream((/nl,   nod2D/),  (/nl,   myDim_nod2D/),  'N2',   'brunt väisälä',       '1/s2',bvfreq(:,:),   1, 'y', i_real4)
  call def_stream((/nl,   nod2D/),  (/nl,   myDim_nod2D/),  'Kv',   'Vertical mixing K',   'm2/s',Kv(:,:),       1, 'y', i_real4)

!2D
  call def_stream(nod2D, myDim_nod2D, 'ssh',   'sea surface elevation',   'm',      eta_n,                         1, 'm', i_real4)
!DS for frontier run
!  call def_stream(nod2D,  myDim_nod2D,  't100',  'temperature at 100m',                'C',    tr_arr(12,1:myDim_nod2D,1),   1, 'm', i_real4)
!  call def_stream(elem2D, myDim_elem2D, 'u100',  'horizontal velocity at 100m',      'm/s',    uv(1,12,1:myDim_elem2D),      1, 'm', i_real4)
!  call def_stream(elem2D, myDim_elem2D, 'v100',  'meridional velocity at 100m',      'm/s',    uv(2,12,1:myDim_elem2D),      1, 'm', i_real4)
!  call def_stream(nod2D,  myDim_nod2D,   't30',  'temperature at 30m',                 'C',    tr_arr(5,1:myDim_nod2D,1),    1, 'm', i_real4)
!  call def_stream(elem2D, myDim_elem2D,  'u30',  'horizontal velocity at 30m',       'm/s',    uv(1,5,1:myDim_elem2D),       1, 'm', i_real4)
!  call def_stream(elem2D, myDim_elem2D,  'v30',  'meridional velocity at 30m',       'm/s',    uv(2,5,1:myDim_elem2D),       1, 'm', i_real4)
!DS

  call def_stream(nod2D, myDim_nod2D, 'sst',   'sea surface temperature', 'C',      tr_arr(1,1:myDim_nod2D,1),     1, 'm', i_real4)
  call def_stream(nod2D, myDim_nod2D, 'sss',   'sea surface salinity',    'psu',    tr_arr(1,1:myDim_nod2D,2),     1, 'm', i_real4)
  call def_stream(nod2D, myDim_nod2D, 'vve',   'vertical velocity',       'm/s',    Wvel(5,:),                     1, 'm', i_real4)
  call def_stream(nod2D, myDim_nod2D, 'uice',  'ice velocity x',          'm/s',    u_ice,                         1, 'm', i_real4)
  call def_stream(nod2D, myDim_nod2D, 'vice',  'ice velocity y',          'm/s',    v_ice,                         1, 'm', i_real4)
  call def_stream(nod2D, myDim_nod2D, 'a_ice', 'ice concentration',       '%',      a_ice(1:myDim_nod2D),          1, 'm', i_real4)
  call def_stream(nod2D, myDim_nod2D, 'm_ice', 'ice height',              'm',      m_ice(1:myDim_nod2D),          1, 'm', i_real4)
  call def_stream(nod2D, myDim_nod2D, 'm_snow','snow height',             'm',      m_snow(1:myDim_nod2D),         1, 'm', i_real4)
  call def_stream(nod2D, myDim_nod2D, 'MLD1',   'Mixed Layer Depth',      'm',      MLD1(1:myDim_nod2D),           1, 'm', i_real4)
  call def_stream(nod2D, myDim_nod2D, 'MLD2',   'Mixed Layer Depth',      'm',      MLD2(1:myDim_nod2D),           1, 'm', i_real4)
  call def_stream(nod2D, myDim_nod2D, 'fh',     'heat flux',              'W',      heat_flux(:),                  1, 'm', i_real4)
  call def_stream(nod2D, myDim_nod2D, 'fw',     'fresh water flux',       'm/s',    water_flux(:),                 1, 'm', i_real4)
  call def_stream(nod2D, myDim_nod2D, 'atmice_x', 'stress atmice x',      'N/m2',   stress_atmice_x(:),            1, 'm', i_real4)
  call def_stream(nod2D, myDim_nod2D, 'atmice_y', 'stress atmice y',      'N/m2',   stress_atmice_y(:),            1, 'm', i_real4)
  call def_stream(nod2D, myDim_nod2D, 'atmoce_x', 'stress atmoce x',      'N/m2',   stress_atmoce_x(:),            1, 'm', i_real4)
  call def_stream(nod2D, myDim_nod2D, 'atmoce_y', 'stress atmoce y',      'N/m2',   stress_atmoce_y(:),            1, 'm', i_real4)
  call def_stream(nod2D, myDim_nod2D, 'iceoce_x', 'stress iceoce x',      'N/m2',   stress_iceoce_x(:),            1, 'm', i_real4)
  call def_stream(nod2D, myDim_nod2D, 'iceoce_y', 'stress iceoce y',      'N/m2',   stress_iceoce_y(:),            1, 'm', i_real4)
  call def_stream(nod2D, myDim_nod2D, 'alpha',  'thermal expansion',      'none',   sw_alpha(1,:),                 1, 'm', i_real4)
  call def_stream(nod2D, myDim_nod2D, 'beta',   'saline contraction',     'none',   sw_beta (1,:),                 1, 'm', i_real4)
  call def_stream(nod2D, myDim_nod2D, 'runoff', 'river runoff',           'none',   runoff(:),                     1, 'y', i_real4)
  call def_stream(nod2D, myDim_nod2D, 'evap',   'evaporation',            'm/s',    evaporation(:),                1, 'm', i_real4)
  call def_stream(nod2D, myDim_nod2D, 'prec',   'precicipation rain',     'm/s',    prec_rain(:),                  1, 'm', i_real4)
#if defined (__oifs)
  call def_stream(nod2D, myDim_nod2D, 'alb',    'ice albedo',             'none',   ice_alb(:),                    1, 'm', i_real4)
#endif

  if (trim(mix_scheme)=='KPP') then
     call def_stream(nod2D, myDim_nod2D, 'hbl',    'HBL KPP',                'none',   hbl(:),                        1, 'm', i_real4)
     call def_stream(nod2D, myDim_nod2D, 'Bo',     'surface boyancy flux',   'm2/s3',  Bo(:),                         1, 'm', i_real4)
  end if

  if (Redi) then
     call def_stream((/nl-1  , nod2D /), (/nl-1,   myDim_nod2D /), 'Redi_K',   'Redi diffusion coefficient', 'm2/s', Ki(:,:),    1, 'm', i_real4)
  end if

  if (Fer_GM) then
     call def_stream((/nl-1, elem2D/), (/nl-1, myDim_elem2D/), 'bolus_u', 'GM bolus velocity U',      'm/s',  fer_uv(1,:,:), 1, 'y', i_real4)
     call def_stream((/nl-1, elem2D/), (/nl-1, myDim_elem2D/), 'bolus_v', 'GM bolus velocity V',      'm/s',  fer_uv(2,:,:), 1, 'y', i_real4)
     call def_stream((/nl  , nod2D /), (/nl,   myDim_nod2D /), 'bolus_w', 'GM bolus velocity W',      'm/s',  fer_Wvel(:,:), 1, 'y', i_real4)
     call def_stream((/nl  , nod2D /), (/nl,   myDim_nod2D /), 'fer_K',   'GM, stirring diffusivity', 'm2/s', fer_k(:,:),    1, 'm', i_real4)

     call def_stream(nod2D, myDim_nod2D, 'fer_C', 'GM,   depth independent speed',  'm/s' ,  fer_c(1:myDim_nod2D),           1, 'm', i_real4)
     call def_stream(nod2D, myDim_nod2D, 'reso',  'mesh resolution',                'm',     mesh_resolution(1:myDim_nod2D), 1, 'm', i_real4)
  end if
  if (ldiag_solver) then
     call def_stream(nod2D, myDim_nod2D, 'rhs_diag',  'SSH_STIFF*d_eta', 'none',      rhs_diag(1:myDim_nod2D),     10, 's', i_real4)
     call def_stream(nod2D, myDim_nod2D, 'ssh_rhs',   'ssh_rhs',         'none',      ssh_rhs (1:myDim_nod2D),     10, 's', i_real4)
  end if
  if (lcurt_stress_surf) then
     call def_stream(nod2D, myDim_nod2D, 'curl_stress_surf',  'vorticity of the surface stress', 'none', curl_stress_surf(1:myDim_nod2D),  10, 'm', i_real4)
  end if
  if (ldiag_curl_vel3) then
     call def_stream((/nl-1, nod2D/),  (/nl-1, myDim_nod2D/),  'curl_u',     'relative vorticity',          '1/s',   vorticity,                   1, 'm', i_real4)
 !    call def_stream(nod2D,  myDim_nod2D,                      'curl_u100',  'relative vorticity at 100m',  '1/s',   vorticity(12,1:myDim_nod2D), 1, 'd', i_real4)
 !    call def_stream(nod2D,  myDim_nod2D,                      'curl_u280',  'relative vorticity at 280m',  '1/s',   vorticity(16,1:myDim_nod2D), 1, 'd', i_real4)
  end if

  if (whichEVP==2) then
     call def_stream(elem2D, myDim_elem2D, 'alpha_EVP', 'alpha in EVP', 'n/a', alpha_evp_array,  1, 'd', i_real4)
     call def_stream(nod2D,  myDim_nod2D,  'beta_EVP',  'beta in EVP',  'n/a', beta_evp_array,   1, 'd', i_real4)
  end if
end subroutine ini_mean_io
!
!--------------------------------------------------------------------------------------------
!
function get_dimname(n) result(s)
  implicit none
  integer       :: n
  character(50) :: s
  if (n==nod2D) then
     s='nod2'
  elseif (n==elem2D) then
     s='elem'
  elseif (n==nl) then
     s='nz'
  elseif (n==nl-1) then
     s='nz1'
  else
     s='unknown'
     if (mype==0) write(*,*) 'WARNING: unknown dimension in mean I/O with zise of ', n
  end if
  end function
!
!--------------------------------------------------------------------------------------------
!
subroutine create_new_file(entry)
  implicit none
  integer                       :: c, j
  character(2000)               :: att_text

  type(Meandata), intent(inout) :: entry
  ! Serial output implemented so far
  if (mype/=0) return
  c=1
  entry%error_status=0
  ! create an ocean output file
  write(*,*) 'initializing I/O file for ', trim(entry%name)

  entry%error_status(c) = nf_create(entry%filename, IOR(NF_NOCLOBBER,IOR(NF_NETCDF4,NF_CLASSIC_MODEL)), entry%ncid); c=c+1

  do j=1, entry%ndim
!___Create mesh related dimentions__________________________________________
     entry%error_status(c) = nf_def_dim(entry%ncid, entry%dimname(j), entry%glsize(j), entry%dimID(j)); c=c+1
  end do
!___Create time related dimentions__________________________________________
  entry%error_status(c) = nf_def_dim(entry%ncid, 'time', NF_UNLIMITED, entry%recID);                     c=c+1
!___Define the time and iteration variables_________________________________
  entry%error_status(c) = nf_def_var(entry%ncid, 'time', NF_DOUBLE, 1, entry%recID, entry%tID); c=c+1

  att_text='time'
  entry%error_status(c) = nf_put_att_text(entry%ncid, entry%tID, 'long_name', len_trim(att_text), trim(att_text)); c=c+1
  write(att_text, '(a14,I4.4,a1,I2.2,a1,I2.2,a6)'), 'seconds since ', yearold, '-', 1, '-', 1, ' 0:0:0'
  entry%error_status(c) = nf_put_att_text(entry%ncid, entry%tID, 'units', len_trim(att_text), trim(att_text)); c=c+1

  if (entry%accuracy == i_real8) then
     entry%error_status(c) = nf_def_var(entry%ncid, trim(entry%name), NF_DOUBLE, entry%ndim+1, &
                                      (/entry%dimid(1:entry%ndim), entry%recID/), entry%varID); c=c+1
  elseif (entry%accuracy == i_real4) then
     entry%error_status(c) = nf_def_var(entry%ncid, trim(entry%name), NF_REAL, entry%ndim+1, &
                                      (/entry%dimid(1:entry%ndim), entry%recID/), entry%varID); c=c+1
  elseif (entry%accuracy == i_int2) then
     entry%error_status(c) = nf_def_var(entry%ncid, trim(entry%name), NF_SHORT, entry%ndim+1, &
                                      (/entry%dimid(1:entry%ndim), entry%recID/), entry%varID); c=c+1

     entry%error_status(c) = nf_put_att_real(entry%ncid, entry%varID, 'scale_factor', NF_REAL, 1, entry%scale_factor); c=c+1
     entry%error_status(c) = nf_put_att_real(entry%ncid, entry%varID, 'add_offset',   NF_REAL, 1, entry%add_offset);   c=c+1
  endif
  entry%error_status(c)=nf_put_att_text(entry%ncid, entry%varID, 'description', len_trim(entry%description), entry%description); c=c+1
  entry%error_status(c)=nf_put_att_text(entry%ncid, entry%varID, 'units',       len_trim(entry%units),       entry%units);       c=c+1
  entry%error_status(c)=nf_close(entry%ncid); c=c+1
  entry%error_count=c-1
end subroutine create_new_file
!
!--------------------------------------------------------------------------------------------
!
subroutine assoc_ids(entry)
  implicit none

  type(Meandata), intent(inout) :: entry
  integer                       :: c, j, k
  real(real64)                  :: rtime !timestamp of the record
  ! Serial output implemented so far
  if (mype/=0) return
  c=1
  entry%error_status=0
  ! open existing netcdf file
  write(*,*) 'associating mean I/O file ', trim(entry%filename)

  entry%error_status(c) = nf_open(entry%filename, nf_nowrite, entry%ncid);

  if (entry%error_status(c) .ne. nf_noerr) then
     call create_new_file(entry) ! error status counter will be reset
     c=entry%error_count+1
     entry%error_status(c) = nf_open(entry%filename, nf_nowrite, entry%ncid); c=c+1
  end if

  do j=1, entry%ndim
!___Create mesh related dimentions__________________________________________
     entry%error_status(c) = nf_inq_dimid(entry%ncid, entry%dimname(j), entry%dimID(j)); c=c+1
  end do
!___Associate time related dimentions_______________________________________
  entry%error_status(c) = nf_inq_dimid (entry%ncid, 'time', entry%recID);          c=c+1
  entry%error_status(c) = nf_inq_dimlen(entry%ncid, entry%recID, entry%rec_count); c=c+1
!___Associate the time and iteration variables______________________________
  entry%error_status(c) = nf_inq_varid(entry%ncid, 'time', entry%tID); c=c+1
!___if the time rtime at the rec_count is larger than ctime we look for the closest record with the 
! timestamp less than ctime
  do k=entry%rec_count, 1, -1
     entry%error_status(c)=nf_get_vara_double(entry%ncid, entry%tID, k, 1, rtime, 1);
     if (ctime > rtime) then
        entry%rec_count=k+1
!       write(*,*) 'I/O '//trim(entry%name)//' : current record = ', entry%rec_count, '; ', entry%rec_count, ' records in the file;'
        exit ! a proper rec_count detected, exit the loop
     end if
     if (k==1) then
        write(*,*) 'I/O '//trim(entry%name)//' WARNING: the existing output file will be overwritten'//'; ', entry%rec_count, ' records in the file;'
        entry%rec_count=1
        exit ! no appropriate rec_count detected
     end if
  end do
  c=c+1 ! check will be made only for the last nf_get_vara_double

  entry%rec_count=max(entry%rec_count, 1)
!___Associate physical variables____________________________________________
  entry%error_status(c) = nf_inq_varid(entry%ncid, entry%name, entry%varID); c=c+1
  entry%error_status(c)=nf_close(entry%ncid); c=c+1
  entry%error_count=c-1
  write(*,*) trim(entry%name)//': current mean I/O counter = ', entry%rec_count
end subroutine assoc_ids
!
!--------------------------------------------------------------------------------------------
!
subroutine write_mean(entry)
  implicit none
  type(Meandata), intent(inout) :: entry
  real(real64)  , allocatable   :: aux_r8(:)
  real(real32),   allocatable   :: aux_r4(:)
  integer(int16), allocatable   :: aux_i2(:)
  
  integer                       :: i, size1, size2
  integer                       :: c, lev
  ! Serial output implemented so far
  if (mype==0) then
     c=1
     write(*,*) 'writing mean record for ', trim(entry%name), '; rec. count = ', entry%rec_count
     entry%error_status(c)=nf_open(entry%filename, nf_write, entry%ncid); c=c+1
     entry%error_status(c)=nf_put_vara_double(entry%ncid, entry%Tid, entry%rec_count, 1, ctime, 1); c=c+1
  end if
!_______writing 2D fields________________________________________________
     if (entry%ndim==1) then
        size1=entry%glsize(1)
!___________writing 8 byte real_________________________________________ 
        if (entry%accuracy == i_real8) then
           if (mype==0) allocate(aux_r8(size1))
           if (size1==nod2D)  call gather_nod (entry%local_values_r8(1:entry%lcsize(1),1), aux_r8)
           if (size1==elem2D) call gather_elem(entry%local_values_r8(1:entry%lcsize(1),1), aux_r8)
           if (mype==0) then
              entry%error_status(c)=nf_put_vara_double(entry%ncid, entry%varID, (/1, entry%rec_count/), (/size1, 1/), aux_r8, 1); c=c+1
           end if
           if (mype==0) deallocate(aux_r8)
        
!___________writing real 4 byte real _________________________________________ 
        elseif (entry%accuracy == i_real4) then
           if (mype==0) allocate(aux_r4(size1))
           if (size1==nod2D)  call gather_nod (entry%local_values_r4(1:entry%lcsize(1),1), aux_r4)
           if (size1==elem2D) call gather_elem(entry%local_values_r4(1:entry%lcsize(1),1), aux_r4)
           if (mype==0) then
              entry%error_status(c)=nf_put_vara_real(entry%ncid, entry%varID, (/1, entry%rec_count/), (/size1, 1/), aux_r4, 1); c=c+1
           end if
           if (mype==0) deallocate(aux_r4)

!___________writing real as 2 byte integer _________________________________________ 
        elseif (entry%accuracy == i_int2) then
           if (mype==0) allocate(aux_i2(size1))
           if (size1==nod2D)  call gather_nod (entry%local_values_i2(1:entry%lcsize(1),1), aux_i2)
           if (size1==elem2D) call gather_elem(entry%local_values_i2(1:entry%lcsize(1),1), aux_i2)
           if (mype==0) then
              entry%error_status(c)=nf_put_vara_int2(entry%ncid, entry%varID, (/1, entry%rec_count/), (/size1, 1/), aux_i2, 1); c=c+1
           end if
           if (mype==0) deallocate(aux_i2)
        else
           if (mype==0) write(*,*) 'not supported output accuracy for mean I/O file.'
           call par_ex
           stop
           
        endif

!_______writing 3D fields________________________________________________
     elseif (entry%ndim==2) then
        size1=entry%glsize(1)
        size2=entry%glsize(2)
!___________writing 8 byte real_________________________________________ 
        if (entry%accuracy == i_real8) then
           if (mype==0) allocate(aux_r8(size2))
           do lev=1, size1
              if (size1==nod2D  .or. size2==nod2D)  call gather_nod (entry%local_values_r8(lev,1:entry%lcsize(2)),  aux_r8)
              if (size1==elem2D .or. size2==elem2D) call gather_elem(entry%local_values_r8(lev,1:entry%lcsize(2)),  aux_r8)
              if (mype==0) then
                 entry%error_status(c)=nf_put_vara_double(entry%ncid, entry%varID, (/lev, 1, entry%rec_count/), (/1, size2, 1/), aux_r8, 1); c=c+1
              end if
           end do
           if (mype==0) deallocate(aux_r8)
!___________writing real 4 byte real _________________________________________ 
        elseif (entry%accuracy == i_real4) then
           if (mype==0) allocate(aux_r4(size2))
           do lev=1, size1
              if (size1==nod2D  .or. size2==nod2D)  call gather_nod (entry%local_values_r4(lev,1:entry%lcsize(2)), aux_r4)
              if (size1==elem2D .or. size2==elem2D) call gather_elem(entry%local_values_r4(lev,1:entry%lcsize(2)), aux_r4)
              if (mype==0) then
                 entry%error_status(c)=nf_put_vara_real(entry%ncid, entry%varID, (/lev, 1, entry%rec_count/), (/1, size2, 1/), aux_r4, 1); c=c+1
              end if
           end do
           if (mype==0) deallocate(aux_r4)
!___________writing real as 2 byte integer _________________________________________ 
        elseif (entry%accuracy == i_int2) then
           if (mype==0) allocate(aux_i2(size2))
           do lev=1, size1           
              if (size1==nod2D  .or. size2==nod2D)  call gather_nod (entry%local_values_i2(lev,1:entry%lcsize(2)), aux_i2)
              if (size1==elem2D .or. size2==elem2D) call gather_elem(entry%local_values_i2(lev,1:entry%lcsize(2)), aux_i2)
              if (mype==0) then
                 entry%error_status(c)=nf_put_vara_int2(entry%ncid, entry%varID, (/lev, 1, entry%rec_count/), (/1, size2, 1/), aux_i2, 1); c=c+1
              end if
           end do
           if (mype==0) deallocate(aux_i2)
        else
           if (mype==0) write(*,*) 'not supported output accuracy for mean I/O file.'
           call par_ex
           stop
           
        endif
     else
        if (mype==0) write(*,*) 'not supported shape of array in mean I/O file'
           call par_ex
           stop
     end if

  if (mype==0) entry%error_count=c-1
  call was_error(entry)
  if (mype==0) entry%error_status(1)=nf_close(entry%ncid);
  entry%error_count=1
  call was_error(entry)
end subroutine write_mean
!
!--------------------------------------------------------------------------------------------
!
subroutine update_means
  implicit none
  type(Meandata), pointer :: entry
  integer                 :: n

  do n=1, io_NSTREAMS
     entry=>io_stream(n)
!_____________ compute in 8 byte accuracy _________________________
     if (entry%accuracy == i_real8) then
        if (entry%ndim==1) then 
           entry%local_values_r8(1:entry%lcsize(1),1) = &
           entry%local_values_r8(1:entry%lcsize(1),1) + entry%ptr2(1:entry%lcsize(1))
        elseif (entry%ndim==2) then 
           entry%local_values_r8(1:entry%lcsize(1),1:entry%lcsize(2)) = &
           entry%local_values_r8(1:entry%lcsize(1),1:entry%lcsize(2)) + entry%ptr3(1:entry%lcsize(1),1:entry%lcsize(2))
        else
           if (mype==0) write(*,*) 'not supported size in update_means'
           call par_ex
           stop
        end if

!_____________ compute in 4 byte accuracy _________________________
     elseif (entry%accuracy == i_real4 .or. entry%accuracy == i_int2) then
        if (entry%ndim==1) then 
           entry%local_values_r4(1:entry%lcsize(1),1) = &
           entry%local_values_r4(1:entry%lcsize(1),1) + real(entry%ptr2(1:entry%lcsize(1)),real32)
        elseif (entry%ndim==2) then 
           entry%local_values_r4(1:entry%lcsize(1),1:entry%lcsize(2)) = &
           entry%local_values_r4(1:entry%lcsize(1),1:entry%lcsize(2)) + &
                 real(entry%ptr3(1:entry%lcsize(1),1:entry%lcsize(2)),real32)
        else
           if (mype==0) write(*,*) 'not supported size in update_means'
           call par_ex
           stop
        end if

     else

           if (mype==0) write(*,*) 'not supported output accuracy in update_means:',entry%accuracy,'for',trim(entry%name)
           call par_ex
           stop
     endif

     entry%addcounter=entry%addcounter+1
  end do
end subroutine update_means
!
!--------------------------------------------------------------------------------------------
!
subroutine output(istep)
  implicit none

  integer       :: istep
  logical, save :: lfirst=.true.
  integer       :: mpierr
  integer       :: n
  logical       :: do_output
  type(Meandata), pointer :: entry
  real(real64)  :: inv_addcounter_r8
  real(real32)  :: inv_addcounter_r4, inv_scale_factor

  ctime=timeold+(dayold-1.)*86400
  if (lfirst) call ini_mean_io

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
        entry%filename=trim(ResultPath)//trim(entry%name)//'.'//trim(runid)//'.'//cyearnew//'.nc'
        call assoc_ids(entry)
        if (entry%accuracy == i_real8) then
           inv_addcounter_r8 = 1._WP/real(entry%addcounter,real64)
           entry%local_values_r8 = entry%local_values_r8 * inv_addcounter_r8  ! compute_means
           call write_mean(entry)
           entry%local_values_r8 = 0. ! clean_meanarrays

        elseif (entry%accuracy == i_real4) then
           inv_addcounter_r4 = 1._WP/real(entry%addcounter,real32)
           entry%local_values_r4 = entry%local_values_r4 * inv_addcounter_r4 ! compute_means
           call write_mean(entry)
           entry%local_values_r4 = 0. ! clean_meanarrays

        elseif (entry%accuracy == i_int2) then
           ! compute_means and compress to int2, use scale_factor and add_offset to maintain as much accuracy as possible
           inv_addcounter_r4 = 1./real(entry%addcounter,real32)
           inv_scale_factor  = 1./entry%scale_factor 
           entry%local_values_i2 = nint( (entry%local_values_r4 * inv_addcounter_r4 - entry%add_offset)* inv_scale_factor,int16)
  
           ! write a warning if the values exceed the interval on which scale_factor and add_offset are based
           if (minval(entry%local_values_r4) < entry%min_value) then
              print *,'! WARNING ! Check output of ',trim(entry%name),' (MPI-task',mype,'):'
              print *,'            The minimum',minval(entry%local_values_r4),'is smaller than'
              print *,'            the range [',entry%min_value,',',entry%min_value,'] converted to int2.'
              print *,'            Adjust the interval or choose real4 in ini_mean, io_meandata.F90'
           endif
           if (maxval(entry%local_values_r4) > entry%max_value) then
              print *,'! WARNING ! Check output of ',trim(entry%name),' (MPI-task',mype,'):'
              print *,'            The maximum',maxval(entry%local_values_r4),'is larger than'
              print *,'            the range [',entry%min_value,',',entry%min_value,'] converted to int2.'
              print *,'            Adjust the interval or choose real4 in ini_mean, io_meandata.F90'
           endif
           call write_mean(entry)
           entry%local_values_r4 = 0. ! clean_meanarrays
        endif  ! accuracy

        entry%addcounter   = 0  ! clean_meanarrays
     endif
  end do
  lfirst=.false.
end subroutine output
!
!--------------------------------------------------------------------------------------------
!
subroutine def_stream3D(glsize, lcsize, name, description, units, data, freq, freq_unit, accuracy, minvalue, maxvalue)
  implicit none
  integer                              :: glsize(2), lcsize(2)
  character(len=*),      intent(in)    :: name, description, units
  real(kind=WP), target, intent(inout) :: data(:,:)
  integer,               intent(in)    :: freq
  character,             intent(in)    :: freq_unit
  integer,               intent(in)    :: accuracy
  real(kind=WP), intent(in), optional  :: minvalue, maxvalue
  type(Meandata),        allocatable   :: tmparr(:)
  type(Meandata),        pointer       :: entry

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
  entry%ptr3 => data
  entry%accuracy = accuracy

  if (accuracy == i_real8) then
     allocate(entry%local_values_r8(lcsize(1), lcsize(2)))
     entry%local_values_r8 = 0. 
  elseif (accuracy == i_real4) then
     allocate(entry%local_values_r4(lcsize(1), lcsize(2)))
     entry%local_values_r4 = 0.  
  elseif (accuracy == i_int2) then
     allocate(entry%local_values_r4(lcsize(1), lcsize(2)))
     allocate(entry%local_values_i2(lcsize(1), lcsize(2)))
     entry%local_values_r4 = 0. 
     entry%local_values_i2 = 0 
     if (present(minvalue) .and. present(maxvalue)) then
        entry%scale_factor = (maxvalue - minvalue) / real(2**16-1,int32)
        entry%add_offset   = 0.5*(maxvalue + minvalue)
        entry%min_value    = minvalue
        entry%max_value    = maxvalue
     else
        entry%scale_factor = 1.
        entry%add_offset   = 0.
        entry%min_value    = 0.
        entry%max_value    = real(2**16-1,int32)
        if (mype==0) then
           print *,'Warning: netcdf-output of',name,'is set to short integer,'
           print *,'but no interval [minvalue,maxvalue] is given, thus scale_factor=1., add_offset=0.'
           print *,'This may result in a huge loss of accuracy!' 
        endif
     endif
  endif ! accuracy

  entry%ndim=2
  entry%lcsize=lcsize
  entry%glsize=glsize
  entry%name = name
  entry%description = description
  entry%units = units

  entry%dimname(1)=get_dimname(glsize(1))
  entry%dimname(2)=get_dimname(glsize(2))
  entry%freq=freq
  entry%freq_unit=freq_unit
  ! clean_meanarrays
  entry%addcounter   = 0
  entry%is_in_use=.true.
  io_NSTREAMS=io_NSTREAMS+1

end subroutine def_stream3D
!
!--------------------------------------------------------------------------------------------
!
subroutine def_stream2D(glsize, lcsize, name, description, units, data, freq, freq_unit, accuracy, minvalue, maxvalue)
  implicit none
  integer,               intent(in)    :: glsize, lcsize
  character(len=*),      intent(in)    :: name, description, units
  real(kind=WP), target, intent(inout) :: data(:)
  integer,               intent(in)    :: freq
  character,             intent(in)    :: freq_unit
  integer,               intent(in)    :: accuracy
  integer, intent(in), optional        :: minvalue, maxvalue
  type(Meandata),        allocatable   :: tmparr(:)
  type(Meandata),        pointer       :: entry

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
  entry%ptr2 => data
  entry%accuracy = accuracy
  if (accuracy == i_real8) then
     allocate(entry%local_values_r8(lcsize, 1))
     ! clean_meanarrays
     entry%local_values_r8 = 0. 

  elseif (accuracy == i_real4) then
     allocate(entry%local_values_r4(lcsize, 1))
     ! clean_meanarrays
     entry%local_values_r4 = 0.  

  elseif (accuracy == i_int2) then
     allocate(entry%local_values_r4(lcsize, 1))
     allocate(entry%local_values_i2(lcsize, 1))
     ! clean_meanarrays
     entry%local_values_r4 = 0.  
     entry%local_values_i2 = 0  

     if (present(minvalue) .and. present(maxvalue)) then
        entry%scale_factor = (maxvalue - minvalue) / real(2**16-1,real32)
        entry%add_offset   = minvalue
     else
        entry%scale_factor = 1.
        entry%add_offset   = 0.
        if (mype==0) then
           print *,'Warning: netcdf-output of',name,'is set to short integer,'
           print *,'but no interval [minvalue,maxvalue] is given, thus scale_factor=1., add_offset=0.'
           print *,'This may result in a huge loss of accuracy!' 
        endif
     endif
  endif

  entry%ndim=1
  entry%lcsize=(/lcsize, 1/)
  entry%glsize=(/glsize, 1/)
  entry%name = name
  entry%description = description
  entry%units = units

  entry%dimname(1)=get_dimname(glsize)
  entry%dimname(2)='unknown'
  entry%freq=freq
  entry%freq_unit=freq_unit
  entry%addcounter   = 0
  entry%is_in_use=.true.
  io_NSTREAMS=io_NSTREAMS+1

end subroutine def_stream2D
!
!--------------------------------------------------------------------------------------------
!
subroutine was_error(entry)
  implicit none
  type(Meandata), intent(inout) :: entry
  integer                       :: k, status, ierror

  call MPI_BCast(entry%error_count, 1,  MPI_INTEGER, 0, MPI_COMM_FESOM, ierror)
  call MPI_BCast(entry%error_status(1), entry%error_count, MPI_INTEGER, 0, MPI_COMM_FESOM, ierror)

  do k=1, entry%error_count
     status=entry%error_status(k)
     if (status .ne. nf_noerr) then
        if (mype==0) call handle_err(status)
        call par_ex
        stop
     end if
  end do
end subroutine was_error

end module io_MEANDATA

