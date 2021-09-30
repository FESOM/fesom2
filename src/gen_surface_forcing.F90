MODULE g_sbf
   !!===========================================================================
   !! Ocean forcing:
   !!===========================================================================
   !! History: 0.1 ! 07/2015 I. Kuznetsov
   !!
   !! WARNING: in getcoeffld;nc_readTimeGrid
   !!   module will flip data in infiles for lat from -90 to 90 (NCEP-DOE Reanalysis 2 standart)
   !!   if you are going to use other coordinates in input files, please rewrite getcoeffld and nc_readTimeGrid functions.
   !!! WARNING : for now reading only files were time counts in hours from 1800 year; nc_readTimeGrid ; done FOR NCEP data
   !!! WARNING : move time forward for dt/2 , last is a last+dt from last -1 ; nc_readTimeGrid ; done FOR NCEP data

   !! Description:
   !!   read and interpolate atmpospheric forcing on model grid,
   !!     or use constants from namelist each time step
   !!
   !!   first initialization before first time step
   !!     model will read namelist, made some checks and prepare first interpolation coeficients
   !!   during model time steping, each time step sb_do is calling
   !!     during run model check if there is a time to read new data and construct new interpolation coefients
   !!     and interpolate data for new time step
   !!
   !!   taux and tuay defined and allocated outside of this module, but will be changed in this module
   !!   qns - Downward Non Solar heat flux over the ocean defined here and changed here, to use it outside:
   !!   USE fv_sbc
   !!   emp - evaporation minus precipitation defined here and changed here, to use it outside:
   !!   USE fv_sbc
   !!
   !!   Only NetCDF format is supported!
   !!   we assume that all NetCDF files have identical grid and time variable
   !!
   !! public:
   !!   sbc_ini  -- inizialization atmpospheric forcing
   !!   sbc_do   -- provide a sbc (surface boundary conditions) each time step
   !!
   USE o_ARRAYS
   USE MOD_MESH
   USE o_PARAM
   USE g_PARSUP
   USE g_comm_auto
   USE g_support
   USE g_rotate_grid
   USE g_config, only: dummy, ClimateDataPath, dt
   USE g_clock,  only: timeold, timenew, dayold, daynew, yearold, yearnew, cyearnew
   USE g_forcing_arrays,    only: runoff
   USE g_read_other_NetCDF, only: read_other_NetCDF, read_2ddata_on_grid_netcdf
   IMPLICIT NONE

   include 'netcdf.inc'

   public  sbc_ini  ! routine called before 1st time step (open files, read namelist,...)
   public  sbc_do   ! routine called each time step to provide a sbc fileds (wind,...)
   public  sbc_end  ! routine called after last time step
   public  julday   ! get julian day from date
   public  atmdata
   public  i_totfl, i_xwind, i_ywind, i_humi, i_qsr, i_qlw, i_tair, i_prec, i_mslp, i_cloud, i_snow
   public  l_xwind, l_ywind, l_humi, l_qsr, l_qlw, l_tair, l_prec, l_mslp, l_cloud, l_snow
   private

   integer :: i_totfl = 10! total number of fluxes
   integer :: i_xwind = 1 ! index of 10m wind velocity (x-component) [m/s]
   integer :: i_ywind = 2 ! index of 10m wind velocity (y-component) [m/s]
   integer :: i_humi  = 3 ! index of specific humidity               [kg/kg]
   integer :: i_qsr   = 4 ! index of solar heat                      [W/m2]
   integer :: i_qlw   = 5 ! index of Long wave                       [W/m2]
   integer :: i_tair  = 6 ! index of 2m air temperature              [degK]
   integer :: i_prec  = 7 ! index of total precipitation (rain+snow) [Kg/m^2/s]
   integer :: i_mslp  = 8 ! index of mean sea level pressure         [Pascals]
   integer :: i_cloud = 9 ! index of mean sea level pressure         [0-1]
   integer :: i_snow  = 10! index of mean sea level pressure         [Kg/m^2/s]

   logical :: l_totfl = .false. 
   logical :: l_xwind = .false. 
   logical :: l_ywind = .false. 
   logical :: l_humi  = .false.
   logical :: l_qsr   = .false.
   logical :: l_qlw   = .false.
   logical :: l_tair  = .false.
   logical :: l_prec  = .false.
   logical :: l_mslp  = .false.
   logical :: l_cloud = .false.
   logical :: l_snow  = .false.

   character(10),      save   :: runoff_data_source='CORE2'
   character(len=MAX_PATH), save   :: nm_runoff_file    ='runoff.nc'

   character(10),      save   :: sss_data_source   ='CORE2'
   character(len=MAX_PATH), save   :: nm_sss_data_file  ='PHC2_salx.nc'

   logical :: runoff_climatology =.false.

   real(wp), allocatable, save, dimension(:), public     :: qns   ! downward non solar heat over the ocean [W/m2]
   real(wp), allocatable, save, dimension(:), public     :: qsr   ! downward solar heat over the ocean [W/m2]
   real(wp), allocatable, save, dimension(:), public     :: emp   ! evaporation minus precipitation        [kg/m2/s]
   real(4), allocatable, dimension(:,:),     private, save, target     :: sbcdata1_,sbcdata2_
!   real(wp), allocatable, save, dimension(:), public     :: qns_2   ! downward non solar heat over the ocean [W/m2]
!   real(wp), allocatable, save, dimension(:), public     :: qsr_2   ! downward solar heat over the ocean [W/m2]
!   real(wp), allocatable, save, dimension(:), public     :: emp_2   ! evaporation minus precipitation        [kg/m2/s]

!   real(wp), allocatable, save, dimension(:), public     :: taux_node_2 ! wind at 10m        [m/s]
!   real(wp), allocatable, save, dimension(:), public     :: tauy_node_2 ! wind at 10m        [m/s]

!   windx and windy are defined in fv_var and allocated in main
!   real(wp), allocatable, save, dimension(:), public     :: windx ! wind at 10m        [m/s]
!   real(wp), allocatable, save, dimension(:), public     :: windy ! wind at 10m        [m/s]
!   mslp now defined in fv_var and allocated in main module
!   real(wp), allocatable, save, dimension(:), public     :: mslp  ! mean sea level pressure [Pascals]
!


   ! namelists
   integer, save  :: nm_sbc_unit     = 103       ! unit to open namelist file, skip 100-102 for cray fortran
   logical        :: ic_cyclic=.true.
  !============== namelistatmdata variables ================
   integer, save  :: nm_sbc       = 1        ! data  1= constant, 2=from file
   character(len=256), save   :: nm_xwind_file = 'xwind.dat' ! name of file with winds/stress, if netcdf file then provide only name from "nameyyyy.nc" yyyy.nc will be added by model
   character(len=256), save   :: nm_ywind_file = 'ywind.dat' ! name of file with winds/stress, if netcdf file then provide only name from "nameyyyy.nc" yyyy.nc will be added by model
   character(len=256), save   :: nm_humi_file  = 'humidity.dat' ! name of file with humidity,  if netcdf file then provide only name from "nameyyyy.nc" yyyy.nc will be added by model
   character(len=256), save   :: nm_qsr_file   = 'qsr.dat'   ! name of file with solar heat,   if netcdf file then provide only name from "nameyyyy.nc" yyyy.nc will be added by model
   character(len=256), save   :: nm_qlw_file   = 'qlw.dat'   ! name of file with Long wave,    if netcdf file then provide only name from "nameyyyy.nc" yyyy.nc will be added by model
   character(len=256), save   :: nm_tair_file  = 'tair.dat'  ! name of file with 2m air temperature, if netcdf file then provide only name from "nameyyyy.nc" yyyy.nc will be added by model
   character(len=256), save   :: nm_prec_file  = 'prec.dat'  ! name of file with total precipitation, if netcdf file then provide only name from "nameyyyy.nc" yyyy.nc will be added by model
   character(len=256), save   :: nm_snow_file  = 'snow.dat'  ! name of file with snow  precipitation, if netcdf file then provide only name from "nameyyyy.nc" yyyy.nc will be added by model
   character(len=256), save   :: nm_mslp_file  = 'mslp.dat'  ! name of file with mean sea level pressure, if netcdf file then provide only name from "nameyyyy.nc" yyyy.nc will be added by model
   character(len=256), save   :: nm_cloud_file = 'cloud.dat'  ! name of file with clouds, if netcdf file then provide only name from "nameyyyy.nc" yyyy.nc will be added by model

   character(len=34), save   :: nm_xwind_var = 'uwnd' ! name of variable in file with wind
   character(len=34), save   :: nm_ywind_var = 'vwnd' ! name of variable in file with wind
   character(len=34), save   :: nm_humi_var  = 'shum' ! name of variable in file with humidity
   character(len=34), save   :: nm_qsr_var   = 'dswrf'! name of variable in file with solar heat
   character(len=34), save   :: nm_qlw_var   = 'dlwrf'! name of variable in file with Long wave
   character(len=34), save   :: nm_tair_var  = 'air'  ! name of variable in file with 2m air temperature
   character(len=34), save   :: nm_prec_var  = 'prate'! name of variable in file with total precipitation
   character(len=34), save   :: nm_snow_var  = 'snow' ! name of variable in file with snow  precipitation
   character(len=34), save   :: nm_mslp_var  = 'mslp' ! name of variable in file with mean sea level pressure
   character(len=34), save   :: nm_cloud_var = 'cloud'! name of variable in file with clouds

   ! ========== netCDF time param
   integer, save :: nm_nc_iyear = 1948    ! initial year of time axis in netCDF (1948 like CoastDat,1800 NCEP)
   integer, save :: nm_nc_imm   = 1       ! initial month of time axis in netCDF
   integer, save :: nm_nc_idd   = 1       ! initial day of time axis in netCDF
   real,    save :: nm_nc_freq  = 86400.0 ! time units coef (86400 CoastDat, 24 NCEP)
   integer, save :: nm_nc_tmid  = 1       ! 1 if the time stamps are given at the mid points of the netcdf file, 0 otherwise!
   logical, save :: y_perpetual=.false.

   integer,save            :: warn       ! warning switch node/element coordinate out of forcing bounds

   real(wp), allocatable, save, dimension(:,:)   :: coef_b ! time inerp coef. b (x=a*t+b)
   real(wp), allocatable, save, dimension(:,:)   :: coef_a ! time inerp coef. a (x=a*t+b)

   real(wp), allocatable, save, dimension(:,:)   :: atmdata ! atmosperic data for current time step

   type, public ::   flfi_type    !flux file informations
      character(len = MAX_PATH)                 :: file_name ! file name
      character(len = 34)                  :: var_name  ! variable name in the NetCDF file
      character(len = 34)                  :: calendar  ! variable name in the NetCDF file
      integer                              :: nc_Nlon
      integer                              :: nc_Nlat
      integer                              :: nc_Ntime
      real(wp), allocatable, dimension(:)  :: nc_lon, nc_lat, nc_time
      ! time index for NC time array
      integer                              :: t_indx    ! now time index in nc_time array
      integer                              :: t_indx_p1 ! now time index +1 in nc_time array
      integer read_forcing_rootrank
      logical async_netcdf_allowed
      real(4), allocatable, dimension(:,:)     :: sbcdata_a
      integer sbcdata_a_t_index
      real(4), allocatable, dimension(:,:)     :: sbcdata_b
      integer sbcdata_b_t_index
      ! ========== interpolation coeficients
   end type flfi_type
   type(flfi_type), allocatable, save, target :: sbc_flfi(:)  !array for information about flux files
   integer,  allocatable, dimension(:,:)     :: bilin_indx_i ! indexs i for interpolation
   integer,  allocatable, dimension(:,:)     :: bilin_indx_j ! indexs j for interpolation
   !flip latitude from infiles (for example  NCEP-DOE Reanalysis 2 standart)
   integer, save              :: flip_lat ! 1 if we need to flip
!============== NETCDF ==========================================


CONTAINS
   SUBROUTINE nc_readTimeGrid(flf)
   ! Read time array and grid from nc file
      IMPLICIT NONE
 
     type(flfi_type),intent(inout) :: flf
      integer                      :: iost !I/O status     
      integer                      :: ncid      ! netcdf file id
      integer                      :: i
      ! ID dimensions and variables:
      integer                      :: id_lon
      integer                      :: id_lat
      integer                      :: id_lond
      integer                      :: id_latd
      integer                      :: id_time
      integer                      :: id_timed      
      integer                      :: nf_start(4)
      integer                      :: nf_edges(4)         
      integer                      :: ierror              ! return error code
      character(len=20)            :: aux_calendar
      integer                      :: aux_len

      !open file
      if (mype==0) then
         iost = nf_open(trim(flf%file_name),NF_NOWRITE,ncid)
      end if

      call MPI_BCast(iost, 1, MPI_INTEGER, 0, MPI_COMM_FESOM, ierror)
      call check_nferr(iost,flf%file_name)

      ! get dimensions
      if (mype==0) then
         iost = nf_inq_dimid(ncid,    "LAT",      id_latd)
         if (iost .ne. NF_NOERR) then
            iost = nf_inq_dimid(ncid, "lat",      id_latd)
         end if
         if (iost .ne. NF_NOERR) then
            iost = nf_inq_dimid(ncid, "latitude", id_latd)
         end if
         if (iost .ne. NF_NOERR) then
            iost = nf_inq_dimid(ncid, "LAT1",     id_latd)
         end if
      end if
      call MPI_BCast(iost, 1, MPI_INTEGER, 0, MPI_COMM_FESOM, ierror)
      call check_nferr(iost,flf%file_name)  

      if (mype==0) then 
         iost = nf_inq_dimid(ncid,    "LON",       id_lond)
         if (iost .ne. NF_NOERR) then
            iost = nf_inq_dimid(ncid, "lon",       id_lond)
         end if
         if (iost .ne. NF_NOERR) then
            iost = nf_inq_dimid(ncid, "longitude", id_lond)
         end if
         if (iost .ne. NF_NOERR) then
            iost = nf_inq_dimid(ncid, "LON1",      id_lond)
         end if
      end if
      call MPI_BCast(iost, 1, MPI_INTEGER, 0, MPI_COMM_FESOM, ierror)
      call check_nferr(iost,flf%file_name) 

      if (mype==0) then   
         iost = nf_inq_dimid(ncid, "TIME", id_timed)
         if      (iost .ne. NF_NOERR) then
                 iost = nf_inq_dimid(ncid, "time",  id_timed)
         end if
         if      (iost .ne. NF_NOERR) then
                 iost = nf_inq_dimid(ncid, "TIME1", id_timed)
         end if
      end if
      call MPI_BCast(iost, 1, MPI_INTEGER, 0, MPI_COMM_FESOM, ierror)
      call check_nferr(iost,flf%file_name)  

      ! get variable id
      if (mype==0) then
         iost = nf_inq_varid(ncid,    "LAT",      id_lat)
         if (iost .ne. NF_NOERR) then
            iost = nf_inq_varid(ncid, "lat",      id_lat)
         end if
         if (iost .ne. NF_NOERR) then
            iost = nf_inq_varid(ncid, "latitude", id_lat)
         end if
         if (iost .ne. NF_NOERR) then
            iost = nf_inq_varid(ncid, "LAT1",     id_lat)
         end if
      end if
      call MPI_BCast(iost, 1, MPI_INTEGER, 0, MPI_COMM_FESOM, ierror)
      call check_nferr(iost,flf%file_name)
      if (mype==0) then
         iost = nf_inq_varid(ncid,    "LON",       id_lon)
         if      (iost .ne. NF_NOERR) then
            iost = nf_inq_varid(ncid, "longitude", id_lon)
         end if
         if (iost .ne. NF_NOERR) then
            iost = nf_inq_varid(ncid, "lon",       id_lon)
         end if
         if (iost .ne. NF_NOERR) then
            iost = nf_inq_varid(ncid, "LON1",      id_lon)
         end if
      end if
      call MPI_BCast(iost, 1, MPI_INTEGER, 0, MPI_COMM_FESOM, ierror)
      call check_nferr(iost,flf%file_name)

      if (mype==0) then
         iost = nf_inq_varid(ncid, "TIME", id_time)
         if      (iost .ne. NF_NOERR) then
                 iost = nf_inq_varid(ncid, "time", id_time)
         end if
         if      (iost .ne. NF_NOERR) then
                 iost = nf_inq_varid(ncid, "TIME1",id_time)
         end if
      end if
      call MPI_BCast(iost, 1, MPI_INTEGER, 0, MPI_COMM_FESOM, ierror)
      call check_nferr(iost,flf%file_name)   
      ! get dimensions size
      if (mype==0) then
         iost = nf_inq_dimlen(ncid, id_latd, flf%nc_Nlat)
      end if
      call MPI_BCast(iost, 1, MPI_INTEGER, 0, MPI_COMM_FESOM, ierror)
      call check_nferr(iost,flf%file_name)
      if (mype==0) then      
         iost = nf_inq_dimlen(ncid, id_lond, flf%nc_Nlon)
      end if
      call MPI_BCast(iost, 1, MPI_INTEGER, 0, MPI_COMM_FESOM, ierror)
      call check_nferr(iost,flf%file_name)   
      if (mype==0) then      
         iost = nf_inq_dimlen(ncid, id_timed,flf%nc_Ntime)
      end if
      call MPI_BCast(iost, 1, MPI_INTEGER, 0, MPI_COMM_FESOM, ierror)
      call check_nferr(iost,flf%file_name) 
      flf%nc_Nlon=flf%nc_Nlon+2 !for the halo in case of periodic boundary
      call MPI_BCast(flf%nc_Nlon,   1, MPI_INTEGER, 0, MPI_COMM_FESOM, ierror)
      call MPI_BCast(flf%nc_Nlat,   1, MPI_INTEGER, 0, MPI_COMM_FESOM, ierror)
      call MPI_BCast(flf%nc_Ntime,  1, MPI_INTEGER, 0, MPI_COMM_FESOM, ierror)
         
      if (.not. allocated(flf%nc_time)) then
         allocate( flf%nc_lon(flf%nc_Nlon), flf%nc_lat(flf%nc_Nlat),&
                &       flf%nc_time(flf%nc_Ntime))
      else
      ! only the temporal axis is allowed to vary between the files
         deallocate(flf%nc_time)
           allocate(flf%nc_time(flf%nc_Ntime))
      end if  
    !____________________________________________________________________________   
    !read variables from file
    ! read lat
      if (mype==0) then
         nf_start(1)=1
         nf_edges(1)=flf%nc_Nlat
         iost = nf_get_vara_double(ncid, id_lat, nf_start, nf_edges, flf%nc_lat)
      end if
      call MPI_BCast(iost, 1, MPI_INTEGER, 0, MPI_COMM_FESOM, ierror)
      call check_nferr(iost,flf%file_name)
      
    ! read lon  
      if (mype==0) then
         nf_start(1)=1
         nf_edges(1)=flf%nc_Nlon-2
         iost = nf_get_vara_double(ncid, id_lon, nf_start, nf_edges, flf%nc_lon(2:flf%nc_Nlon-1))
         flf%nc_lon(1)        =flf%nc_lon(flf%nc_Nlon-1)
         flf%nc_lon(flf%nc_Nlon)  =flf%nc_lon(2)
      end if
      call MPI_BCast(iost, 1, MPI_INTEGER, 0, MPI_COMM_FESOM, ierror)      
      call check_nferr(iost,flf%file_name)
    !____________________________________________________________________________
    ! read time axis from file
      if (mype==0) then
         nf_start(1)=1
         nf_edges(1)=flf%nc_Ntime
         iost = nf_get_vara_double(ncid, id_time, nf_start, nf_edges, flf%nc_time)
         ! digg for calendar attribute in time axis variable         
      end if
      call MPI_BCast(flf%nc_time, flf%nc_Ntime,   MPI_DOUBLE_PRECISION, 0, MPI_COMM_FESOM, ierror)
      call MPI_BCast(iost, 1, MPI_INTEGER, 0, MPI_COMM_FESOM, ierror)
      call check_nferr(iost,flf%file_name)
      
      ! digg for calendar attribute in time axis variable
      if (mype==0 .and. use_flpyrcheck) then
         iost = nf_inq_attlen(ncid, id_time,'calendar',aux_len)
         iost = nf_get_att(ncid, id_time,'calendar',aux_calendar)
         aux_calendar = aux_calendar(1:aux_len)
         
         if (iost .ne. NF_NOERR) then
            flf%calendar='none'
            write(*,*) ' --> could not find/read calendar attribute in the time axis'
            write(*,*) '     of the forcing file (Is this right?). I assume there is'
            write(*,*) '     none and proceed in CORE2 style without leap years!'
         else
            flf%calendar=lowercase(aux_calendar)
            write(*,*) ' --> found calendar attr. in time axis: |',trim(flf%calendar),'|' 
         end if 
         
         ! check for calendar and include_fleapyear consistency
         if ((trim(flf%calendar).eq.'none')   .or. &
             (trim(flf%calendar).eq.'noleap') .or. &
             (trim(flf%calendar).eq.'365_days')) then
            if (include_fleapyear .eqv. .true.) then
                print *, achar(27)//'[33m'
                write(*,*) '____________________________________________________________'
                write(*,*) ' WARNING: It looks like you want to use CORE forcing, Right?'
                write(*,*) '          but setted include_fleapyear=.true.. CORE forcing '
                write(*,*) '          does not contain any leap years or particular '
                write(*,*) '          calender option (julian, gregorian). So if im right,'
                write(*,*) '          please go to namelist.config and set '
                write(*,*) '          include_fleapyear=.false. otherwise comment this '
                write(*,*) '          message block in gen_surface_forcing.F90.'
                write(*,*) '____________________________________________________________'
                print *, achar(27)//'[0m'
                call par_ex(0)
            end if
         elseif ((trim(flf%calendar).eq.'julian')    .or. &
                 (trim(flf%calendar).eq.'gregorian') .or. &
                 (trim(flf%calendar).eq.'proleptic_gregorian') .or. &
                 (trim(flf%calendar).eq.'standard')) then
            if (include_fleapyear .eqv. .false.) then
                print *, achar(27)//'[33m'
                write(*,*) '____________________________________________________________'
                write(*,*) ' WARNING: It looks like you want to use either JRA55, ERA,'
                write(*,*) '          NCEP or a similar forcing, Right?, but setted '
                write(*,*) '          include_fleapyear=.false. JRA55, ERA or NCEP contain'
                write(*,*) '          all fleapyears and use a specific calendar option '
                write(*,*) '          (julian, gregorian). So that the calendars in FESOM2.0'
                write(*,*) '          work properly, when using these forcings '
                write(*,*) '          include_fleapyear must be true. So if im right, please go'
                write(*,*) '          to namelist.config and set include_fleapyear=.true. '
                write(*,*) '          otherwise comment this message block in'
                write(*,*) '          gen_surface_forcing.F90'
                write(*,*) '____________________________________________________________'
                print *, achar(27)//'[0m'
                call par_ex(0)
            end if 
         else
            print *, achar(27)//'[31m'
            write(*,*) '____________________________________________________________'
            write(*,*) ' ERROR: I am not familiar with the found calendar option,'
            write(*,*) '        dont know what to do. Either talk to the FESOM2 developers'
            write(*,*) '        or add the calendar option by your self in ...'
            write(*,*) '        gen_surface_forcing.F90, line:364-367'
            write(*,*) '                                                            '
            write(*,*) '        elseif ((trim(flf%calendar).eq."julian")      .or. &'
            write(*,*) '                (trim(flf%calendar).eq."gregorian")   .or. &'
            write(*,*) '                (trim(flf%calendar).eq."NEW_CALENDAR").or. &'
            write(*,*) '                (trim(flf%calendar).eq."standard")) then    '
            write(*,*) '                                                            '
            write(*,*) '        The time axis calendar attribute can be checked for '
            write(*,*) '        example with ncdump -h forcing_file.nc '
            write(*,*) '____________________________________________________________'
            print *, achar(27)//'[0m'
            call par_ex(0)
         end if 
      end if
      
    ! transform time axis accorcing to calendar and include_fleapyear=.true./.false. flag  
      flf%nc_time = flf%nc_time / nm_nc_freq + julday(nm_nc_iyear,nm_nc_imm,nm_nc_idd)
      if (nm_nc_tmid/=1) then
         if (flf%nc_Ntime > 1) then
            do i = 1, flf%nc_Ntime-1
               flf%nc_time(i) = (flf%nc_time(i+1) + flf%nc_time(i))/2.0_WP
            end do
           flf%nc_time(flf%nc_Ntime) = flf%nc_time(flf%nc_Ntime) + (flf%nc_time(flf%nc_Ntime) - flf%nc_time(flf%nc_Ntime-1))/2.0
         end if
      end if
      call MPI_BCast(flf%nc_lon,   flf%nc_Nlon,   MPI_DOUBLE_PRECISION, 0, MPI_COMM_FESOM, ierror)
      call MPI_BCast(flf%nc_lat,   flf%nc_Nlat,   MPI_DOUBLE_PRECISION, 0, MPI_COMM_FESOM, ierror)
    
    !___________________________________________________________________________
    !flip lat and data in case of lat from -90 to 90
    !!!! WARNING this is temporal solution, needs some more checks
      flip_lat = 0
      if ( flf%nc_Nlat > 1 ) then
         if ( flf%nc_lat(1) > flf%nc_lat(flf%nc_Nlat) ) then
            flip_lat = 1
            flf%nc_lat=flf%nc_lat(flf%nc_Nlat:1:-1)
            if (mype==0) write(*,*) "fv_sbc: nc_readTimeGrid: FLIP lat and data while lat from -90 to 90"
         endif
      endif

      if (mype==0) then
         iost = nf_close(ncid)
      end if
      call MPI_BCast(iost, 1, MPI_INTEGER, 0, MPI_COMM_FESOM, ierror)
      call check_nferr(iost,flf%file_name)

      if (ic_cyclic) then
         flf%nc_lon(1)      =flf%nc_lon(1)-360._WP
         flf%nc_lon(flf%nc_Nlon)=flf%nc_lon(flf%nc_Nlon)+360._WP
      end if 
   END SUBROUTINE nc_readTimeGrid

   SUBROUTINE nc_sbc_ini_fillnames(yyyy)
      integer, intent(in)         :: yyyy
      character(len=4)            :: yyear

      write(yyear,"(I4)") yyyy
      if (y_perpetual)    yyear = ''

      !! ** Purpose : Fill names of sbc_flfi array (file names and variable names)

      !prepare proper nc file (add year and .nc to the end of the file name from namelist
      if (l_xwind) write(sbc_flfi(i_xwind)%file_name, *) trim(nm_xwind_file),trim(yyear),'.nc'
      if (l_ywind) write(sbc_flfi(i_ywind)%file_name, *) trim(nm_ywind_file),trim(yyear),'.nc'
      if (l_humi)  write(sbc_flfi(i_humi)%file_name,  *) trim(nm_humi_file), trim(yyear),'.nc'
      if (l_qsr)   write(sbc_flfi(i_qsr)%file_name,   *) trim(nm_qsr_file),  trim(yyear),'.nc'
      if (l_qlw)   write(sbc_flfi(i_qlw)%file_name,   *) trim(nm_qlw_file),  trim(yyear),'.nc'
      if (l_tair)  write(sbc_flfi(i_tair)%file_name,  *) trim(nm_tair_file), trim(yyear),'.nc'
      if (l_prec)  write(sbc_flfi(i_prec)%file_name,  *) trim(nm_prec_file), trim(yyear),'.nc'
      if (l_snow)  write(sbc_flfi(i_snow)%file_name,  *) trim(nm_snow_file), trim(yyear),'.nc'
      if (l_mslp)  write(sbc_flfi(i_mslp)%file_name,  *) trim(nm_mslp_file), trim(yyear),'.nc'
      if (l_cloud) write(sbc_flfi(i_cloud)%file_name, *) trim(nm_cloud_file),trim(yyear),'.nc'

      if (l_xwind) sbc_flfi(i_xwind)%file_name=ADJUSTL(trim(sbc_flfi(i_xwind)%file_name))
      if (l_ywind) sbc_flfi(i_ywind)%file_name=ADJUSTL(trim(sbc_flfi(i_ywind)%file_name))
      if (l_humi)  sbc_flfi(i_humi)%file_name=ADJUSTL(trim(sbc_flfi(i_humi)%file_name))
      if (l_qsr)   sbc_flfi(i_qsr)%file_name=ADJUSTL(trim(sbc_flfi(i_qsr)%file_name))
      if (l_qlw)   sbc_flfi(i_qlw)%file_name=ADJUSTL(trim(sbc_flfi(i_qlw)%file_name))
      if (l_tair)  sbc_flfi(i_tair)%file_name=ADJUSTL(trim(sbc_flfi(i_tair)%file_name))
      if (l_prec)  sbc_flfi(i_prec)%file_name=ADJUSTL(trim(sbc_flfi(i_prec)%file_name))
      if (l_snow)  sbc_flfi(i_snow)%file_name=ADJUSTL(trim(sbc_flfi(i_snow)%file_name))
      if (l_mslp)  sbc_flfi(i_mslp)%file_name=ADJUSTL(trim(sbc_flfi(i_mslp)%file_name))
      if (l_cloud) sbc_flfi(i_cloud)%file_name=ADJUSTL(trim(sbc_flfi(i_cloud)%file_name))

      if (l_xwind) sbc_flfi(i_xwind)%var_name=ADJUSTL(trim(nm_xwind_var))
      if (l_ywind) sbc_flfi(i_ywind)%var_name=ADJUSTL(trim(nm_ywind_var))
      if (l_humi)  sbc_flfi(i_humi)%var_name=ADJUSTL(trim(nm_humi_var))
      if (l_qsr)   sbc_flfi(i_qsr)%var_name=ADJUSTL(trim(nm_qsr_var))
      if (l_qlw)   sbc_flfi(i_qlw)%var_name=ADJUSTL(trim(nm_qlw_var))
      if (l_tair)  sbc_flfi(i_tair)%var_name=ADJUSTL(trim(nm_tair_var))
      if (l_prec)  sbc_flfi(i_prec)%var_name=ADJUSTL(trim(nm_prec_var))
      if (l_snow)  sbc_flfi(i_snow)%var_name=ADJUSTL(trim(nm_snow_var))
      if (l_mslp)  sbc_flfi(i_mslp)%var_name=ADJUSTL(trim(nm_mslp_var))
      if (l_cloud) sbc_flfi(i_cloud)%var_name=ADJUSTL(trim(nm_cloud_var))
   END SUBROUTINE nc_sbc_ini_fillnames

   SUBROUTINE nc_sbc_ini(rdate, mesh)
      !!---------------------------------------------------------------------
      !! ** Purpose : initialization of ocean forcing from NETCDF file
      !!----------------------------------------------------------------------

      IMPLICIT NONE
      real(wp),intent(in) :: rdate ! initialization date
      integer             :: idate
      integer             :: yyyy,mm,dd

      integer                  :: i
      integer                  :: sbc_alloc
      logical, save            :: lfirst=.true.
      integer                  :: elnodes(4) !4 nodes from one element
      integer                  :: numnodes   ! nu,ber of nodes in elem (3 for triangle, 4 for ... )
      real(wp)                 :: x, y       ! coordinates of elements
      integer                  :: fld_idx
      type(flfi_type), pointer :: flf
      type(t_mesh), intent(in)              , target :: mesh

#include "associate_mesh.h"

! used for interpolate on elements
!      ALLOCATE( bilin_indx_i(elem2D),bilin_indx_j(elem2D), &
!              & qns(elem2D), emp(elem2D), qsr(elem2D),     &
!                   &      STAT=sbc_alloc )
! used to inerpolate on nodes
      warn = 0


      ! get ini year; Fill names of sbc_flfi
      idate=int(rdate)
      call calendar_date(idate,yyyy,mm,dd)
      call nc_sbc_ini_fillnames(yyyy)
      ! we assume that all NetCDF files have identical grid and time variable
      do fld_idx = 1, i_totfl
         call nc_readTimeGrid(sbc_flfi(fld_idx))
      end do
      if (lfirst) then
      do fld_idx = 1, i_totfl
         flf=>sbc_flfi(fld_idx)
         ! prepare nearest coordinates in INfile , save to bilin_indx_i/j
         do i = 1, myDim_nod2D+eDim_nod2D
            x  = geo_coord_nod2D(1,i)/rad
            if (x < 0) x=x+360._WP
            y  = geo_coord_nod2D(2,i)/rad

            ! find nearest
            if ( x < flf%nc_lon(flf%nc_Nlon) .and. x >= flf%nc_lon(1) ) then
               call binarysearch(flf%nc_Nlon, flf%nc_lon, x, bilin_indx_i(fld_idx, i))
            else ! NO extrapolation in space
               if ( x < flf%nc_lon(1) ) then
                  bilin_indx_i(fld_idx, i)=-1
               else
                  bilin_indx_i(fld_idx, i)=0
               end if
            end if
            if ( y < flf%nc_lat(flf%nc_Nlat) .and. y >= flf%nc_lat(1) ) then
               call binarysearch(flf%nc_Nlat, flf%nc_lat, y, bilin_indx_j(fld_idx, i))
            else ! NO extrapolation in space
               if ( y < flf%nc_lat(1) ) then
                  bilin_indx_j(fld_idx, i)=-1
               else
                  bilin_indx_j(fld_idx, i)=0
               end if
            end if
            if (warn == 0) then
               if (bilin_indx_i(fld_idx, i) < 1 .or. bilin_indx_j(fld_idx, i) < 1) then
!                 WRITE(*,*) '     WARNING:  node/element coordinate out of forcing bounds,'
!                 WRITE(*,*) '        nearest value will be used as a constant field'
                  warn = 1
               end if
            end if
         end do
      end do
      lfirst=.false.
      end if
      do fld_idx = 1, i_totfl
         ! get first coefficients for time interpolation on model grid for all data
         call getcoeffld(fld_idx, rdate, mesh)
      end do
         ! interpolate in time
      call data_timeinterp(rdate)
   END SUBROUTINE nc_sbc_ini

   SUBROUTINE getcoeffld(fld_idx, rdate, mesh)
      use forcing_provider_async_module
      use io_netcdf_workaround_module
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE getcoeffld ***
      !!
      !! ** Purpose : read fields from files, interpolate on model mesh and prepare interpolation coefficients
      !! ** Method  :
      !! ** Action  :
      !!----------------------------------------------------------------------
      IMPLICIT NONE
      integer, intent(in)  :: fld_idx
      real(wp),intent(in)  :: rdate ! initialization date
      integer              :: iost  !I/O status
      integer              :: ncid      ! netcdf file id
      ! ID dimensions and variables:
      integer              :: id_data
      integer              :: nf_start(4)
      integer              :: nf_edges(4)
!      integer              :: zero_year,yyyy,mm,dd
!      character(len = 256) :: att_string ! attribute
      integer              :: i,j,ii, ip1, jp1, extrp
      integer              :: sbc_alloc, itot

      real(wp)             :: denom, x1, x2, y1, y2, x, y
      real(wp)             :: now_date

!     real(wp), allocatable, dimension(:,:)  :: sbcdata1,sbcdata2
      real(wp)             :: data1,data2
      real(wp)             :: delta_t   ! time(t_indx) - time(t_indx+1)

      integer              :: elnodes(4) !4 nodes from one element
      integer              :: numnodes   ! nu,ber of nodes in elem (3 for triangle, 4 for ... )
      integer              :: yyyy,mm,dd
      integer              :: ierror              ! return error code
      integer,   pointer   :: nc_Ntime, nc_Nlon, nc_Nlat, t_indx, t_indx_p1
      character(len=MAX_PATH), pointer   :: file_name
      character(len=34) , pointer   :: var_name
      real(wp),  pointer   :: nc_time(:), nc_lon(:), nc_lat(:)
      type(t_mesh), intent(in) , target :: mesh
      real(4), dimension(:,:), pointer :: sbcdata1, sbcdata2
      logical sbcdata1_from_cache, sbcdata2_from_cache
      integer rootrank

#include  "associate_mesh.h"

      ! fld_idx determines which ouf our forcing fields we use here
      nc_Ntime =>sbc_flfi(fld_idx)%nc_Ntime
      nc_Nlon  =>sbc_flfi(fld_idx)%nc_Nlon
      nc_Nlat  =>sbc_flfi(fld_idx)%nc_Nlat
      t_indx   =>sbc_flfi(fld_idx)%t_indx
      t_indx_p1=>sbc_flfi(fld_idx)%t_indx_p1
      file_name=>sbc_flfi(fld_idx)%file_name
      var_name =>sbc_flfi(fld_idx)%var_name
      nc_time  =>sbc_flfi(fld_idx)%nc_time
      nc_lon   =>sbc_flfi(fld_idx)%nc_lon
      nc_lat   =>sbc_flfi(fld_idx)%nc_lat

      if(.not. allocated(sbc_flfi(fld_idx)%sbcdata_a)) then
        allocate(sbc_flfi(fld_idx)%sbcdata_a(nc_Nlon,nc_Nlat))
        sbc_flfi(fld_idx)%sbcdata_a_t_index = -1
        allocate(sbc_flfi(fld_idx)%sbcdata_b(nc_Nlon,nc_Nlat))
        sbc_flfi(fld_idx)%sbcdata_b_t_index = -1
        sbc_flfi(fld_idx)%read_forcing_rootrank = next_io_rank(MPI_COMM_FESOM, sbc_flfi(fld_idx)%async_netcdf_allowed)
      end if
      rootrank = sbc_flfi(fld_idx)%read_forcing_rootrank

      ! no initialization of sbcdata required, the whole array will be overwritten anyway

      ! find time index in files
      now_date = rdate
      call binarysearch(nc_Ntime,nc_time,now_date,t_indx)
      if ( (t_indx < nc_Ntime) .and. (t_indx > 0) ) then
         t_indx_p1 = t_indx + 1
         delta_t   = nc_time(t_indx_p1) - nc_time(t_indx)
      elseif (t_indx > 0) then ! NO extrapolation to future
         t_indx    = nc_Ntime
         t_indx_p1 = t_indx
         delta_t = 1.0_wp
         if (mype==0) then
            write(*,*) 'WARNING: no temporal extrapolation into future (nearest neighbour is used): ', trim(var_name), ' !'
            write(*,*) trim(file_name)
            write(*,*) nc_time(1), nc_time(nc_Ntime), now_date
         end if
      elseif (t_indx < 1) then ! NO extrapolation back in time
         t_indx = 1
         t_indx_p1 = t_indx
         delta_t = 1.0_wp
         if (mype==0) then 
            write(*,*) 'WARNING: no temporal extrapolation back in time (nearest neighbour is used): ', trim(var_name), ' !'
            write(*,*) trim(file_name)
            write(*,*) nc_time(1), nc_time(nc_Ntime), now_date
         end if
      end if


      ! determine if we can use the broadcast cache
      if(yearold == yearnew) then ! todo: simplify if clause
        if(sbc_flfi(fld_idx)%sbcdata_a_t_index == t_indx) then
          sbcdata1_from_cache = .true.
          sbcdata1 => sbc_flfi(fld_idx)%sbcdata_a
          sbcdata2 => sbc_flfi(fld_idx)%sbcdata_b
          sbc_flfi(fld_idx)%sbcdata_b_t_index = t_indx_p1
        else if(sbc_flfi(fld_idx)%sbcdata_b_t_index == t_indx) then
          sbcdata1_from_cache = .true.
          sbcdata1 => sbc_flfi(fld_idx)%sbcdata_b
        
          sbcdata2 => sbc_flfi(fld_idx)%sbcdata_a ! 
          sbc_flfi(fld_idx)%sbcdata_a_t_index = t_indx_p1
        else
          sbcdata1_from_cache = .false.
          sbcdata1 => sbc_flfi(fld_idx)%sbcdata_a
          sbc_flfi(fld_idx)%sbcdata_a_t_index = t_indx

          sbcdata2_from_cache = .false.
          sbcdata2 => sbc_flfi(fld_idx)%sbcdata_b
          sbc_flfi(fld_idx)%sbcdata_b_t_index = t_indx_p1
        end if
      else
          sbcdata1_from_cache = .false.
          sbcdata1 => sbc_flfi(fld_idx)%sbcdata_a
          sbc_flfi(fld_idx)%sbcdata_a_t_index = t_indx

          sbcdata2_from_cache = .false.
          sbcdata2 => sbc_flfi(fld_idx)%sbcdata_b
          sbc_flfi(fld_idx)%sbcdata_b_t_index = t_indx_p1
      end if
      sbcdata2_from_cache = .false.         

      iost = 0
      !open file sbc_flfi
      if (mype==0) then
         !write(*,*) 'check: ', trim(file_name)
      end if

      !read data from file
      if (mype==rootrank) then
         nf_start(1)=1
         nf_edges(1)=nc_Nlon-2
         nf_start(2)=1
         nf_edges(2)=nc_Nlat
         nf_start(3)=t_indx
         nf_edges(3)=1
         if(.not. sbcdata1_from_cache) then
           call forcing_provider%get_forcingdata(i_totfl, fld_idx, sbc_flfi(fld_idx)%async_netcdf_allowed, trim(file_name), yearnew, trim(var_name), t_indx, sbcdata1(2:nc_Nlon-1,1:nc_Nlat))
         end if
         iost = 0
      end if

      if(mype == rootrank) then
        if(.not. sbcdata1_from_cache) then
           sbcdata1(1, 1:nc_Nlat)       = sbcdata1(nc_Nlon-1, 1:nc_Nlat)
           sbcdata1(nc_Nlon,1:nc_Nlat) = sbcdata1(2, 1:nc_Nlat)
        end if
      end if
      if(sbcdata1_from_cache) then ! use the cache instead of bcast
      else
        call MPI_BCast(sbcdata1(1:nc_Nlon,1:nc_Nlat), nc_Nlon*nc_Nlat, MPI_REAL, rootrank, MPI_COMM_FESOM, ierror)
     end if

       ! read next time step in file (check for +1 done before)
      if (mype==rootrank) then
        nf_start(3)=t_indx_p1
        nf_edges(3)=1
        if(.not. sbcdata2_from_cache) then
        call forcing_provider%get_forcingdata(i_totfl, fld_idx, sbc_flfi(fld_idx)%async_netcdf_allowed, trim(file_name), yearnew, trim(var_name), t_indx_p1, sbcdata2(2:nc_Nlon-1,1:nc_Nlat))
      end if
      iost = 0
      end if

      if (mype==rootrank) then
        if(.not. sbcdata2_from_cache) then
         sbcdata2(1, 1:nc_Nlat)       = sbcdata2(nc_Nlon-1, 1:nc_Nlat)
         sbcdata2(nc_Nlon, 1:nc_Nlat) = sbcdata2(2, 1:nc_Nlat)
        end if
      end if
     if(sbcdata2_from_cache) then ! use the cache instead of bcast
      !
     else
       call MPI_BCast(sbcdata2(1:nc_Nlon,1:nc_Nlat), nc_Nlon*nc_Nlat, MPI_REAL, rootrank, MPI_COMM_FESOM, ierror)
     end if

!      !flip data in case of lat from -90 to 90
!      !!!! WARNING
!      if ( flip_lat == 1 ) then
!         sbcdata1=sbcdata1(1:nc_Nlon,nc_Nlat:1:-1)
!         sbcdata2=sbcdata2(1:nc_Nlon,nc_Nlat:1:-1)
!      end if
      ! bilinear space interpolation, and time interpolation ,
      ! data is assumed to be sampled on a regular grid
!!$OMP PARALLEL
!!$OMP DO
      do ii = 1, myDim_nod2D+eDim_nod2D
         i = bilin_indx_i(fld_idx, ii)
         j = bilin_indx_j(fld_idx, ii)
         ip1 = i + 1
         jp1 = j + 1
         x  = geo_coord_nod2D(1,ii)/rad
         if (x < 0.0_WP) x=x+360._WP
         y  = geo_coord_nod2D(2,ii)/rad
         extrp = 0
         if ( i == 0 ) then
            i   = nc_Nlon
            ip1 = i
            extrp = extrp + 1
         end if
         if ( i == -1 ) then
            i   = 1
            ip1 = i
            extrp = extrp + 1
         end if
         if ( j == 0 ) then
            j   = nc_Nlat
            jp1 = j
            extrp = extrp + 2
         end if
         if ( j == -1 ) then
            j   = 1
            jp1 = j
            extrp = extrp + 2
         end if
         x1 = nc_lon(i)
         x2 = nc_lon(ip1)
         y1 = nc_lat(j)
         y2 = nc_lat(jp1)
         if ( extrp == 0 ) then
         ! if point inside forcing domain
            denom = (x2 - x1)*(y2 - y1)
            data1 = ( sbcdata1(i,j)   * (x2-x)*(y2-y)   + sbcdata1(ip1,j)    * (x-x1)*(y2-y) + &
                    sbcdata1(i,jp1) * (x2-x)*(y-y1)   + sbcdata1(ip1, jp1) * (x-x1)*(y-y1)     ) / denom
            data2 = ( sbcdata2(i,j)   * (x2-x)*(y2-y)   + sbcdata2(ip1,j)    * (x-x1)*(y2-y) + &
                    sbcdata2(i,jp1) * (x2-x)*(y-y1)   + sbcdata2(ip1, jp1) * (x-x1)*(y-y1)     ) / denom
         else if ( extrp == 1 ) then !  "extrapolation" in x direction
            denom = (y2 - y1)
            data1 = ( sbcdata1(i,j)   * (y2-y)   + sbcdata1(ip1, jp1) * (y-y1) ) / denom
            data2 = ( sbcdata2(i,j)   * (y2-y)   + sbcdata2(ip1, jp1) * (y-y1) ) / denom
         else if ( extrp == 2 ) then !  "extrapolation" in y direction
            denom = (x2 - x1)
            data1 = ( sbcdata1(i,j)   * (x2-x)   + sbcdata1(ip1, jp1) * (x-x1) ) / denom
            data2 = ( sbcdata2(i,j)   * (x2-x)   + sbcdata2(ip1, jp1) * (x-x1) ) / denom
         else if ( extrp == 3 ) then !  "extrapolation" in x and y direction
            data1 = sbcdata1(i,j)
            data2 = sbcdata2(i,j)
         end if
         ! calculate new coefficients for interpolations
         coef_a(fld_idx, ii) = ( data2 - data1 ) / delta_t !( nc_time(t_indx+1) - nc_time(t_indx) )
         coef_b(fld_idx, ii) = data1 - coef_a(fld_idx, ii) * nc_time(t_indx)

      end do !ii
!!$OMP END DO
!!$OMP END PARALLEL
   END SUBROUTINE getcoeffld

   SUBROUTINE data_timeinterp(rdate)
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE data_timeinterp ***
      !!
      !! ** Purpose : interpolation of fields(interpolated on model grid) from IN files in time
      !! ** Method  :
      !! ** Action  :
      !!----------------------------------------------------------------------
      IMPLICIT NONE
      real(wp),intent(in)    :: rdate  ! seconds

     ! assign data from interpolation to taux and tauy
      integer            :: fld_idx, i,j,ii

!!$OMP PARALLEL
!!$OMP DO
      do fld_idx = 1, i_totfl
         do i = 1, myDim_nod2D+eDim_nod2D
            ! store processed forcing data for fesom computation
            atmdata(fld_idx,i) = rdate * coef_a(fld_idx,i) + coef_b(fld_idx,i)
         end do !nod2D
      end do !fld_idx
!!$OMP END DO
!!$OMP END PARALLEL
   END SUBROUTINE data_timeinterp

   SUBROUTINE sbc_ini(mesh)
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE sbc_ini ***
      !!
      !! ** Purpose : inizialization of ocean forcing
      !! ** Method  :
      !! ** Action  :
      !!----------------------------------------------------------------------
      IMPLICIT NONE

      integer            :: idate ! initialization date
      real(wp)           :: rdate ! initialization date
      integer            :: iost  ! I/O status
      integer            :: sbc_alloc                   !: allocation status

      real(wp)           :: tx, ty
      type(t_mesh), intent(in)   , target :: mesh

      namelist /nam_sbc/ nm_xwind_file, nm_ywind_file, nm_humi_file, nm_qsr_file, &
                        nm_qlw_file, nm_tair_file, nm_prec_file, nm_snow_file, &
                        nm_mslp_file, nm_xwind_var, nm_ywind_var, nm_humi_var, &
                        nm_qsr_var, nm_qlw_var, nm_tair_var, nm_prec_var, nm_snow_var, &
                        nm_mslp_var, nm_cloud_var, nm_cloud_file, nm_nc_iyear, nm_nc_imm, nm_nc_idd, nm_nc_freq, nm_nc_tmid, y_perpetual, &
                        l_xwind, l_ywind, l_humi, l_qsr, l_qlw, l_tair, l_prec, l_mslp, l_cloud, l_snow, &
                        nm_runoff_file, runoff_data_source, runoff_climatology, nm_sss_data_file, sss_data_source
      ! OPEN and read namelist for SBC
      open( unit=nm_sbc_unit, file='namelist.forcing', form='formatted', access='sequential', status='old', iostat=iost )
      if (iost == 0) then
         if (mype==0) WRITE(*,*) '     file   : ', 'namelist_bc.nml',' open ok'
      else
         if (mype==0) WRITE(*,*) 'ERROR: --> bad opening file   : ', 'namelist_bc.nml',' ; iostat=',iost
         call par_ex
         stop
      endif
      READ( nm_sbc_unit, nml=nam_sbc, iostat=iost )
      close( nm_sbc_unit )
      
      if (mype==0) write(*,*) "Start: Ocean forcing inizialization."
      rdate = real(julday(yearnew,1,1))
      rdate = rdate+real(daynew-1,WP)+timenew/86400._WP
      idate = int(rdate)

      if (mype==0) then
         write(*,*) "Start: Ocean forcing inizialization."
         write(*,*) "Surface boundary conditions parameters:"
      end if

      i_totfl=0
      if (l_xwind) then
         if (mype==0) then
            write(*,*) "      nm_xwind_file = ", trim(nm_xwind_file) ," ! name of file with X winds"
            write(*,*) "      nm_xwind_var  = ", trim(nm_xwind_var)  ," ! name of variable in file with wind "
         end if
         i_totfl=i_totfl+1
         i_xwind=i_totfl
      end if

      if (l_ywind) then
         if (mype==0) then
            write(*,*) "      nm_ywind_file = ", trim(nm_ywind_file) ," ! name of file with Y winds"
            write(*,*) "      nm_ywind_var  = ", trim(nm_ywind_var)  ," ! name of variable in file with wind "
         end if
         i_totfl=i_totfl+1
         i_ywind=i_totfl
      end if

      if (l_humi) then
         if (mype==0) then
            write(*,*) "      nm_humi_file = ", trim(nm_humi_file) ," ! name of file with humidity"
            write(*,*) "      nm_humi_var   = ", trim(nm_humi_var) ," ! name of variable in file with humidity  "
         end if
         i_totfl=i_totfl+1
         i_humi =i_totfl
      end if

      if (l_qsr) then
         if (mype==0) then
            write(*,*) "      nm_qsr_file = ", trim(nm_qsr_file) ," ! name of file with solar heat "
            write(*,*) "      nm_qsr_var  = ", trim(nm_qsr_var)  ," ! name of variable in file with solar heat "
         end if
         i_totfl=i_totfl+1
         i_qsr  =i_totfl
      end if

      if (l_qlw) then
         if (mype==0) then
            write(*,*) "      nm_qlw_file   = ", trim(nm_qlw_file) ," ! name of file with Long wave "
            write(*,*) "      nm_qlw_var    = ", trim(nm_qlw_var)  ," ! name of variable in file with Long wave "
         end if
         i_totfl=i_totfl+1
         i_qlw =i_totfl     
      end if

      if (l_tair) then
         if (mype==0) then
            write(*,*) "      nm_tair_file  = ", trim(nm_tair_file) ," ! name of file with 2m air temperature "
            write(*,*) "      nm_tair_var   = ", trim(nm_tair_var)  ," ! name of variable in file with 2m air temperature "
         end if
         i_totfl=i_totfl+1
         i_tair  =i_totfl
      end if

      if (l_prec) then
         if (mype==0) then
            write(*,*) "      nm_prec_file  = ", trim(nm_prec_file) ," ! name of file with total precipitation "
            write(*,*) "      nm_prec_var   = ", trim(nm_prec_var)  ," ! name of variable in file with total precipitation  "
         end if
         i_totfl =i_totfl+1
         i_prec  =i_totfl
      end if

      if (l_snow) then
         if (mype==0) then
            write(*,*) "      nm_snow_file  = ", trim(nm_snow_file) ," ! name of file with snow precipitation "
            write(*,*) "      nm_snow_var   = ", trim(nm_snow_var) , " !name of variable in file with air_pressure_at_sea_level "
         end if
         i_totfl =i_totfl+1
         i_snow  =i_totfl
      end if

      if (l_mslp) then
         if (mype==0) then
            write(*,*) "      nm_mslp_file  = ", trim(nm_mslp_file)," !air_pressure_at_sea_level "
            write(*,*) "      nm_mslp_var   = ", trim(nm_mslp_var) ," !name of variable in file with air_pressure_at_sea_level "
         end if
         i_totfl =i_totfl+1
         i_mslp  =i_totfl
      end if

      if (mype==0) then
         write(*,*) 'total fluxes to read: ', i_totfl
      end if

      ALLOCATE( coef_a(i_totfl,myDim_nod2D+eDim_nod2D), coef_b(i_totfl,myDim_nod2D+eDim_nod2D), &
              & atmdata(i_totfl,myDim_nod2D+eDim_nod2D), &
                   &      STAT=sbc_alloc )
      coef_a       = 0.0_WP             
      coef_b       = 0.0_WP
      atmdata      = 0.0_WP

      ALLOCATE( bilin_indx_i(i_totfl, myDim_nod2D+eDim_nod2D), bilin_indx_j(i_totfl, myDim_nod2D+eDim_nod2D), &
              & qns(myDim_nod2D+eDim_nod2D), emp(myDim_nod2D+eDim_nod2D), qsr(myDim_nod2D+eDim_nod2D),  &
                   &      STAT=sbc_alloc )
      bilin_indx_i = 0.0_WP
      bilin_indx_j = 0.0_WP
      qns          = 0.0_WP
      emp          = 0.0_WP
      qsr          = 0.0_WP
      ALLOCATE(sbc_flfi(i_totfl))
      call nc_sbc_ini(rdate, mesh)
      !==========================================================================
      ! runoff    
      if (runoff_data_source=='CORE1' .or. runoff_data_source=='CORE2' ) then
         ! runoff in CORE is constant in time
         ! Warning: For a global mesh, conservative scheme is to be updated!!
         call read_other_NetCDF(nm_runoff_file, 'Foxx_o_roff', 1, runoff, .false., mesh) 
         runoff=runoff/1000.0_WP  ! Kg/s/m2 --> m/s
      end if

      if (mype==0) write(*,*) "DONE:  Ocean forcing inizialization."
      if (mype==0) write(*,*) 'Parts of forcing data (only constant in time fields) are read'
   END SUBROUTINE sbc_ini

   SUBROUTINE sbc_do(mesh)
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE sbc_do ***
      !!
      !! ** Purpose : provide at each time-step: wind stress, ...
      !! ** Method  :
      !! ** Action  :
      !!----------------------------------------------------------------------
      use g_clock
      IMPLICIT NONE

      real(wp)     :: rdate ! date
      integer      :: fld_idx, i
      logical      :: do_rotation, force_newcoeff, update_monthly_flag
      integer      :: yyyy, dd, mm
      integer,   pointer   :: nc_Ntime, t_indx, t_indx_p1
      real(wp),  pointer   :: nc_time(:)
      character(len=MAX_PATH)   :: filename
      type(t_mesh), intent(in) , target :: mesh
      
#include  "associate_mesh.h"

      force_newcoeff=.false.
      if (yearnew/=yearold) then
         rdate = real(julday(yearnew,1,1),WP)
         call calendar_date(int(rdate),yyyy,dd,mm)
         ! use next set of forcing files
         call nc_sbc_ini_fillnames(yyyy)
         ! we assume that all NetCDF files have identical grid and time variable
         do fld_idx = 1, i_totfl
            call nc_readTimeGrid(sbc_flfi(fld_idx))            
         end do
         force_newcoeff=.true.
      end if
      

      rdate = real(julday(yearnew,1,1),WP)
      rdate = rdate+real(daynew-1,WP)+timenew/86400._WP-dt/86400._WP/2._WP
      do_rotation=.false.

      do fld_idx = 1, i_totfl
         nc_time  =>sbc_flfi(fld_idx)%nc_time
         t_indx_p1=>sbc_flfi(fld_idx)%t_indx_p1
         t_indx   =>sbc_flfi(fld_idx)%t_indx
         nc_Ntime =>sbc_flfi(fld_idx)%nc_Ntime
         if ( ((rdate > nc_time(t_indx_p1)) .and. (nc_time(t_indx) < nc_time(nc_Ntime))) .or. force_newcoeff) then
            ! get new coefficients for time interpolation on model grid for all data
            call getcoeffld(fld_idx, rdate, mesh)
            if (fld_idx==i_xwind) do_rotation=.true.
         endif
      end do

      if (do_rotation) then
         do i=1, myDim_nod2D+eDim_nod2D
            call vector_g2r(coef_a(i_xwind,i), coef_a(i_ywind,i), coord_nod2D(1,i), coord_nod2D(2,i), 0)
            call vector_g2r(coef_b(i_xwind,i), coef_b(i_ywind,i), coord_nod2D(1,i), coord_nod2D(2,i), 0)
         end do
      end if
      
      !==========================================================================

      ! prepare a flag which checks whether to update monthly data (SSS, river runoff)
      update_monthly_flag=( (day_in_month==num_day_in_month(fleapyear,month) .AND. timenew==86400._WP) .OR. mstep==1  )

      ! read in SSS for applying SSS restoring
      if (surf_relax_S > 0._WP) then
         if (sss_data_source=='CORE1' .or. sss_data_source=='CORE2') then
            if (update_monthly_flag) then
               i=month
               if (mstep > 1) i=i+1 
               if (i > 12) i=1
               if (mype==0) write(*,*) 'Updating SSS restoring data for month ', i 
               call read_other_NetCDF(nm_sss_data_file, 'SALT', i, Ssurf, .true., mesh) 
            end if
         end if
      end if

     ! runoff  
     if(runoff_data_source=='Dai09' .or. runoff_data_source=='JRA55') then
       
       if(update_monthly_flag) then
         if(runoff_climatology) then
           !climatology monthly mean
           i=month
           if (mstep > 1) i=i+1 
           if (i > 12) i=1
           if (mype==0) write(*,*) 'Updating monthly climatology runoff for month ', i 
           filename=trim(nm_runoff_file)
           call read_2ddata_on_grid_NetCDF(filename,'runoff', i, runoff, mesh)

           !kg/m2/s -> m/s
           runoff=runoff/1000.0_WP

         else
           !monthly data
           i=month
           if (mstep > 1) i=i+1 
           if (i > 12) i=1
           if (mype==0) write(*,*) 'Updating monthly runoff for month ', i 
           filename=trim(nm_runoff_file)//cyearnew//'.nc' 
           call read_2ddata_on_grid_NetCDF(filename,'runoff', i, runoff, mesh)

           !kg/m2/s -> m/s
           runoff=runoff/1000.0_WP

         end if
       end if

     end if


      ! interpolate in time
      call data_timeinterp(rdate)
   END SUBROUTINE sbc_do


   SUBROUTINE err_call(iost,fname)
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE  err_call ***
      !!
      !! ** Purpose : call Error
      !! ** Method  :
      !! ** Action  :
      !!----------------------------------------------------------------------
   IMPLICIT NONE
      integer, intent(in)            :: iost
      character(len=MAX_PATH), intent(in) :: fname
      write(*,*) 'ERROR: I/O status=',iost,' file= ',fname
      STOP 'ERROR:  stop'


   END SUBROUTINE err_call

   FUNCTION julday(yyyy,mm,dd)

   IMPLICIT NONE
      integer, INTENT(IN) :: mm, dd, yyyy
      integer             :: julday
      ! In this routine julday returns the Julian Day Number that begins at noon of the calendar     
      !    date specified by month mm, day dd , and year yyyy, all integer variables. Positive year
      !    signifies A.D.; negative, B.C. Remember that the year after 1 B.C. was 1 A.D. (from Num. Rec.)
      integer, PARAMETER  :: IGREG=15+31*(10+12*1582)
      ! Gregorian Calendar adopted Oct. 15, 1582.
      integer             :: ja,jm,jy
      if (y_perpetual) then !to work with COREI forcing
         julday=0
         return
      end if
      if (include_fleapyear) then
         jy = yyyy
         if (jy == 0) STOP 'julday: there is no year zero'
         if (jy < 0) jy=jy+1
         if (mm > 2) then
            jm=mm+1
         else
            jy=jy-1
            jm=mm+13
         endif
         julday=int(365.25_WP*jy)+int(30.6001_WP*jm)+dd+1720995
         !Test whether to change to Gregorian Calendar.
         if (dd+31*(mm+12*yyyy) >= IGREG) then
            ja=int(0.01_WP*jy)
            julday=julday+2-ja+int(0.25_wp*ja)
         end if
       else
         julday=365*yyyy
       end if
   END FUNCTION julday


   SUBROUTINE calendar_date(julian,yyyy,mm,dd)

!  Converts a Julian day to a calendar date (year, month and day). Numerical Recipes
   IMPLICIT NONE
!
      integer,intent(in)  :: julian
      integer,intent(out) :: yyyy,mm,dd

      integer, parameter :: IGREG=2299161
      integer            :: ja,jb,jc,jd,je
      real(wp)           :: x
      !
      !-----------------------------------------------------------------------
      if (include_fleapyear) then
         if (julian >= IGREG ) then
            x = ((julian-1867216)-0.25_WP)/36524.25_WP
            ja = julian+1+int(x)-int(0.25*x)
         else
            ja = julian
         end if

         jb = ja+1524
         jc = int(6680 + ((jb-2439870)-122.1_WP)/365.25_WP)
         jd = int(365*jc+(0.25_WP*jc))
         je = int((jb-jd)/30.6001_WP)

         dd = jb-jd-int(30.6001_WP*je)
         mm = je-1
         if (mm > 12) mm = mm-12
         yyyy = jc - 4715
         if (mm > 2) yyyy = yyyy-1
         if (yyyy <= 0) yyyy = yyyy-1
      else
         yyyy=int((real(julian)+1.e-12_WP)/365._WP)
         mm=-1 !not supported (no need so far)
         dd=-1 !not supported (no need so far)
      end if
      return
   END SUBROUTINE calendar_date

   SUBROUTINE sbc_end
      IMPLICIT NONE
      integer      :: fld_idx      
      do fld_idx = 1, i_totfl     
         DEALLOCATE( sbc_flfi(fld_idx)%nc_lon, sbc_flfi(fld_idx)%nc_lat, sbc_flfi(fld_idx)%nc_time)
      end do
      DEALLOCATE( sbc_flfi )
      DEALLOCATE( coef_a, coef_b, atmdata, &
                  &  bilin_indx_i, bilin_indx_j,  &
                  &  qns, emp, qsr)
   END SUBROUTINE sbc_end

   SUBROUTINE check_nferr(iost,fname)
   IMPLICIT NONE
      character(len=MAX_PATH), intent(in) :: fname
      integer, intent(in) :: iost

      if (iost .ne. NF_NOERR) then
         write(*,*) 'ERROR: I/O status= "',trim(nf_strerror(iost)),'";',iost,' file= ',fname
         call par_ex
         stop
      endif
   END SUBROUTINE check_nferr

   SUBROUTINE binarysearch(length, array, value, ind)!, delta)
      ! Given an array and a value, returns the index of the element that
      ! is closest to, but less than, the given value.
      ! Uses a binary search algorithm.
      ! "delta" is the tolerance used to determine if two values are equal
      ! if ( abs(x1 - x2) <= delta) then
      !    assume x1 = x2
      ! endif
      !org. source from: https://github.com/cfinch/Shocksolution_Examples/blob/master/FORTRAN/BilinearInterpolation/interpolation.f90

      IMPLICIT NONE
      integer,  intent(in) :: length
      real(wp), dimension(length), intent(in) :: array
      real(wp), intent(in) :: value
   !   real, intent(in), optional :: delta

   !   integer :: binarysearch
      integer, intent(out) :: ind

      integer :: left, middle, right
      real(wp):: d

   !   if (present(delta) .eqv. .true.) then
   !      d = delta
   !   else
      d = 1e-9_WP
   !   endif
      left = 1
      right = length
      do
         if (left > right) then
            exit
         endif
         middle = nint((left+right) / 2.0_WP)
         if ( abs(array(middle) - value) <= d) then
            ind = middle
            return
         else if (array(middle) > value) then
            right = middle - 1
         else
            left = middle + 1
         end if
      end do
      ind = right
   END SUBROUTINE binarysearch

!-----------------------------------------------------------------------
! This subroutine taken from GOTM src, used to compare with our old fluxes,
!   original code from (see description of code)
!BOP
!
! !ROUTINE: Heat and momentum fluxes according to Fairall et al.
!
! !INTERFACE:
   subroutine fairall(sst,airt,u10,v10,precip,qs,qa,rhoa,evap,taux,tauy,qe,qh)
!
! !DESCRIPTION:
!  The surface momentum flux vector, $(\tau_x^s,\tau_y^s)$,
!  in [N\,m$^{-2}$],
!  the latent heat flux, $Q_e$,
!  and the sensible heat flux, $Q_h$, both in [W\,m$^{-2}$]
!  are calculated here according to the \cite{Fairalletal96a} bulk
!  formulae, which are build on the Liu-Katsaros-Businger
!  (\cite{Liuetal79}) method.
!  Cool skin and warm layer effects are considered according to the
!  suggestions of \cite{Fairalletal96b}.
!
!  The air temperature {\tt airt} and the sea surface temperature
!  {\tt sst} may be given in Kelvin or Celsius:
!  if they are $>$ 100 - Kelvin is assumed.
!
!  This piece of code has been adapted from the COARE code originally
!  written by David Rutgers and Frank Bradley - see
!  http://www.coaps.fsu.edu/COARE/flux\_algor/flux.html.
!
! !USES:
!   use airsea_variables, only: const06,rgas,rho_0,g,rho_0,kappa
!   use airsea_variables, only: qs,qa,rhoa
!   use airsea_variables, only: cpa,cpw
!   use airsea, only: rain_impact,calc_evaporation
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   real(wp), intent(in)                :: sst,airt,u10,v10,precip
   real(wp), intent(in)                :: qa ! specific humidity (kg/kg)
   real(wp), intent(in)                :: qs ! saturation specific humidity
   real(wp), intent(in)                :: rhoa ! air density
!
! !INPUT/OUTPUT PARAMETERS:
   real(wp), intent(out)               :: evap
!
! !OUTPUT PARAMETERS:
   real(wp), intent(out)               :: taux,tauy,qe,qh
!
! !REVISION HISTORY:
!  Original author(s): Adolf Stips
!
! !DEFINED PARAMETERS:
!  Fairall LKB roughness Reynolds number to Von Karman
   real(wp),parameter        :: fdg = 1.0          ! non-dimensional

!  Beta parameter evaluated from Fairall low windspeed turbulence data.
   real(wp),parameter        :: beta = 1.2         ! non-dimensional

!  Zabl      Height (m) of atmospheric boundary layer.
   real(wp),parameter        :: Zabl = 600.0       ! in [m]

   real(wp), parameter       :: r3 = 1.0/3.0
!
!  Liu et al. (1979) look-up table coefficients to compute roughness
!  Reynolds number for temperature (rt) and moisture (rq) as a
!  function of wind Reynolds number (rr):
!
!       rt = Liu_a(:,1) * Rr   ** Liu_b(:,1)    temperature
!       rq = Liu_a(:,2) * Rr   ** Liu_b(:,2)    moisture
!
   real(wp),parameter, dimension(8,2) :: Liu_a = reshape ( &
                 (/ 0.177,  1.376,    1.026,      1.625,   &
                    4.661, 34.904, 1667.190, 588000.0,     &
                    0.292,  1.808,    1.393,      1.956,   &
                    4.994, 30.709, 1448.680, 298000.0 /),  &
                 (/ 8, 2 /) )

   real(wp),parameter, dimension(8,2) :: Liu_b = reshape ( &
                 (/  0.0,    0.929, -0.599, -1.018,        &
                    -1.475, -2.067, -2.907, -3.935,        &
                     0.0,    0.826, -0.528, -0.870,        &
                    -1.297, -1.845, -2.682, -3.616 /),     &
                 (/ 8, 2 /) )

   real(wp),parameter, dimension(9) :: Liu_Rr =            &
                 (/    0.0,  0.11,   0.825,   3.0,         &
                      10.0, 30.0,  100.0,   300.0,         &
                    1000.0 /)
!
!  Height (m) of surface air temperature measurement.
   real(wp), parameter       ::  zt= 2.0
!  Height (m) of surface air humidity measurement
   real(wp), parameter       ::  zq= 2.0
!  Height (m) of surface winds measurement
   real(wp), parameter       ::  zw= 10.0
   integer,  parameter       :: itermax = 20

   real(wp), parameter       :: wgust=0.0


!   real(wp)                    :: psi

   real(wp), parameter         :: cpa=1008.
   real(wp), parameter         :: cpw=3985.
   real(wp), parameter         :: emiss=0.97
   real(wp), parameter         :: bolz=5.67e-8
   real(wp), parameter         :: kelvin=273.16
   real(wp), parameter         :: const06=0.62198
   real(wp), parameter         :: rgas = 287.1    !
   real(wp), parameter         :: g = 9.81        ! [m/s2]
   real(wp), parameter         :: rho_0 = 1025.   ! [kg/m3]
   real(wp), parameter         :: kappa = 0.41    ! von Karman
   logical                     :: rain_impact = .true.
   logical                     :: calc_evaporation =.true.

! !LOCAL VARIABLES:
   real(wp)                  :: tmp,cff,wgus
   real(wp)                  :: L
   real(wp)                  :: Cd
   real(wp)                  :: ta,ta_k,tw,tw_k
   integer                   :: ier,iter,k
   real(wp)                  :: vis_air
   real(wp)                  :: tpsi,qpsi,wpsi,ZWoL,oL,ZToL,ZQoL,ZoW,ZoT, ZoQ
   real(wp)                  :: Wstar,Tstar, Qstar, delQ, delT, rr,rt,rq
   real(wp)                  :: TVstar,Bf, upvel,delw,Wspeed, w
   real(wp)                  :: ri,cd_rain
   real(wp)                  :: x1,x2,x3
   real(wp)                  :: x
   real(wp)                  :: rainfall
   real(wp), parameter       :: eps=1.0e-12
!EOP
!-----------------------------------------------------------------------
!BOC
   evap = 0.0_wp
   w = sqrt(u10*u10+v10*v10)

   if (sst .lt. 100._WP) then
      tw  = sst
      tw_k= sst+kelvin
   else
      tw  = sst-kelvin
      tw_k= sst
   end if

   if (airt .lt. 100._WP) then
      ta_k  = airt + kelvin
      ta = airt
   else
      ta  = airt - kelvin
      ta_k = airt
   end if

!
!  Initialize.
!
   qe   = 0.0_wp
   qh   = 0.0_wp
   taux = 0.0_wp
   tauy = 0.0_wp
   delw=sqrt(w*w+wgust*wgust)
   if (delw .ne. 0.0_WP) then
!-----------------------------------------------------------------------
!     Compute Monin-Obukhov similarity parameters for wind (Wstar),
!     heat (Tstar), and moisture (Qstar), Liu et al. (1979).
!-----------------------------------------------------------------------

!     Kinematic viscosity of dry air (m2/s), Andreas (1989).
      vis_air=1.326e-5_WP*(1.0_WP+ta*(6.542e-3_WP+ta*(8.301e-6_WP-4.84e-9_WP*ta)))

!     Compute latent heat of vaporization (J/kg) at sea surface
      L = (2.501_WP-0.00237_WP*tw)*1.e6_WP
!
!     Assume that wind is measured relative to sea surface and include
!     gustiness.
!     Initialize.
      ier = 0
      delq=qa-qs
      delt=ta-tw

!     Initial guesses for Monin-Obukhov similarity scales.
      ZWoL=0.0_WP
      ZoW=0.0005_WP
      Wstar=0.04_WP*delw
      Tstar=0.04_WP*delt
      Qstar=0.04_WP*delq
      TVstar=Tstar*(1.0_WP+0.61_WP*qa)+0.61_WP*ta_k*Qstar

!     Compute Richardson number.
      ri=g*zw*(delt+0.61_WP*ta_k*delq)/(ta_k*delw*delw)

!     Fairall computes turbulent fluxes only if Ri< 0.25
      if ( ri .le. 0.25_WP) then
!        Iterate until convergence or when IER is negative.  It usually
!        converges within four iterations.
         do iter=1,itermax
            if ( ier .ge. 0 ) then
!              Compute Monin-Obukhov stability parameter, Z/L.
               oL=g*kappa*TVstar/(ta_k*(1.0_WP+0.61_WP*qa)*Wstar*Wstar)
               ZWoL=zw*oL
               ZToL=zt*oL
               ZQoL=zq*oL

!              Evaluate stability functions at Z/L.
               wpsi=psi(1,ZWoL)
               tpsi=psi(2,ZToL)
               qpsi=psi(2,ZQoL)

!              Compute wind scaling parameters, Wstar.
               ZoW=0.011_WP*Wstar*Wstar/g+0.11_WP*vis_air/Wstar
               Wstar=delw*kappa/(log(zw/ZoW)-wpsi)

!              Computes roughness Reynolds number for wind (Rr), heat (Rt),
!              and moisture (Rq). Use Liu et al. (1976) look-up table to
!              compute "Rt" and "Rq" as function of "Rr".
               rr=ZoW*Wstar/vis_air
               if ((rr .ge. 0.0_WP).and.(rr .lt. 1000.0_WP)) then
                  do k=1,8
                     if ((liu_rr(k).le.rr).and.(rr .lt. liu_rr(k+1))) then
                        rt=liu_a(k,1)*rr**liu_b(k,1)
                        rq=liu_a(k,2)*rr**liu_b(k,2)
                     end if
                  end do

!                Compute heat and moisture scaling parameters,
!                Tstar and Qstar.
                  cff=vis_air/Wstar
                  ZoT=rt*cff
                  ZoQ=rq*cff
                  cff=kappa*fdg
                  Tstar=(delt)*cff/(log(zt/ZoT)-tpsi)
                  Qstar=(delq)*cff/(log(zq/ZoQ)-qpsi)

!                 Compute gustiness in wind speed.
                  TVstar=Tstar*(1.0_WP+0.61_WP*qa)+0.61_WP*ta_k*Qstar
                  bf=-g/ta_k*Wstar*TVstar
                  if (bf .gt. 0) then
                     wgus=beta*(bf*Zabl)**r3
                  else
                     wgus=0.0_wp
                  end if
                  delw=sqrt(w*w+wgus*wgus)
               else
                  ier = -2
               end if
            end if
         end do

!        Compute transfer coefficients for momentun (Cd), heat (Ch),
!        and moisture (Ce).
         if (ier .ge. 0) then
            Wspeed=sqrt(w*w+wgus*wgus)
            Cd=Wstar*Wstar/(Wspeed*Wspeed)

!           Compute turbulent sensible heat flux (W/m2), qe.
!           out of ocean is negative
            qe=cpa*rhoa*Wstar*Tstar

!           compute sensible heatflux due to rain fall
            if (rain_impact) then
!              units of qs and qa - should be kg/kg
               rainfall=precip * 1000._WP ! (convert from m/s to kg/m2/s)
               x1 = 2.11e-5_WP*(ta_k/kelvin)**1.94_WP
               x2 = 0.02411_WP*(1.0_WP+ta*(3.309e-3_WP-1.44e-6_WP*ta))/(rhoa*cpa)
               x3 = qa * L /(rgas * ta_K * ta_K)
               cd_rain = 1.0_WP/(1.0_WP+const06*(x3*L*x1)/(cpa*x2))
               cd_rain = cd_rain*cpw*((tw-ta) + (qs-qa)*L/cpa)
               qe = qe - rainfall * cd_rain
            end if

!           Compute turbulent latent heat flux (W/m2), qh.
            qh=L*rhoa*Wstar*Qstar

!           Compute Webb correction (Webb effect) to latent heat flux
            upvel=-1.61_WP*Wstar*Qstar-(1.0_WP+1.61_WP*qa)*Wstar*Tstar/ta_k
            qh=qh-rhoa*L*upvel*qa

!           calculation of evaporation/condensation in m/s
            if (rain_impact .and. calc_evaporation) then
               evap = rhoa/rho_0*Wstar*Qstar
            else
               evap = 0.0_WP
            end if

!           Compute wind stress components (N/m2), Tau.
            cff=rhoa*Cd*Wspeed
            taux=(cff*u10)
            tauy=(cff*v10)

!           Compute momentum flux (N/m2) due to rainfall (kg/m2/s).
!           according to Caldwell and Elliott (1971, JPO)
            if ( rain_impact ) then
               tmp  = 0.85_WP * rainfall
               taux  = taux + tmp * u10
               tauy  = tauy + tmp * v10
            end if

         end if ! ier >0
      end if ! Ri < 0.25
   end if  !delw != 0.0

   return
   end subroutine fairall
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: Calculate the humidity \label{sec:humidity}
!
! !INTERFACE:
   subroutine humidity(hum_method,hum,airp,tw,ta,qa,qs,rhoa,ea,es)
!
! !DESCRIPTION:
!
! This routine calculated the saturation vapour pressure at SST and at
! air temperature, as well as the saturation specific humidty and the
! specific humidity. For the latter, four methods are implemented,
! and the method has to be chosen in the namelist file {\tt airsea.nml}
! as parameter {\tt hum\_method}, see \sect{sec:init-air-sea} for details.
!
!
! !USES:
!   use airsea_variables, only: kelvin,const06,rgas
!   use airsea_variables, only: es,ea,qs,qa,rhoa
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: hum_method
   real(wp), intent(in)                :: hum,airp,tw,ta

   real(wp), intent(out)               :: qa ! specific humidity (kg/kg)
   real(wp), intent(out)               :: ea ! actual water vapor pressure in Pascal
   real(wp), intent(out)               :: es ! saturation vapor pressure
   real(wp), intent(out)               :: qs ! saturation specific humidity
   real(wp), intent(out)               :: rhoa ! air density
!
! !OUTPUT PARAMETERS:
!
! !REVISION HISTORY:
!  Original author(s): Adolf Stips, Hans Burchard & Karsten Bolding
!
! !DEFINED PARAMETERS:
   real(wp), parameter         :: kelvin=273.16
   real(wp), parameter         :: const06=0.62198
   real(wp), parameter         :: rgas = 287.1
!  Note shift of indices for coefficients compared to Lowe (1977, J. Appl. Met.)
   real(wp), parameter       :: a1=6.107799961
   real(wp), parameter       :: a2=4.436518521e-1
   real(wp), parameter       :: a3=1.428945805e-2
   real(wp), parameter       :: a4=2.650648471e-4
   real(wp), parameter       :: a5=3.031240396e-6
   real(wp), parameter       :: a6=2.034080948e-8
   real(wp), parameter       :: a7=6.136820929e-11
!
! !LOCAL VARIABLES:
   real(wp)        :: rh,twet,twet_k,dew,dew_k
!EOP
!-----------------------------------------------------------------------
!BOC
!  saturation vapor pressure - using SST
   es = a1 +tw*(a2+tw*(a3+tw*(a4+tw*(a5+tw*(a6+tw*a7)))))
   es = es * 100.0_WP ! Conversion millibar --> Pascal

!  correction for seawater, following Kraus 1972
!  correcting for salt water assuming 98% RH
   es=0.98_WP * es
!  saturation specific humidity
   qs = const06*es/(airp-0.377_WP*es)

!  must be also calcuated for airtemperature, depending on humidity input
!  see ../ncdf/ncdf_meteo.F90 for defined constants
   select case (hum_method)
      case (1) ! relative humidity in % given
         rh = 0.01_WP * hum
!        saturation vapor pressure at that air temperature
         ea = a1 +ta*(a2+ta*(a3+ta*(a4+ta*(a5+ta*(a6+ta*a7)))))
         ea = ea * 100.0_WP ! Conversion millibar --> Pascal
!        get actual vapor pressure
         ea = rh * ea
!        convert to specific humidity
         qa = const06*ea/(airp-0.377_WP*ea)
      case (2)  ! Specific humidity from wet bulb temperature
!        calculate the SVP at wet bulb temp then
!        use the psychrometer formula to get vapour pressure
!        See Smithsonian Met tables 6th Edition pg 366 eqn 3
!        Make sure this is in degC
         if (hum .lt. 100.0_WP ) then
            twet_k=hum + kelvin
            twet=hum
         else
            twet=hum - kelvin
            twet_k=hum
         end if
!        saturation vapor pressure at wet bulb temperature
         ea = a1 +twet*(a2+twet*(a3+twet*(a4+twet*(a5+twet*(a6+twet*a7)))))
         ea = ea * 100.0_WP ! Conversion millibar --> Pascal
!        actual vapor pressure
         ea = ea - 6.6e-4_WP*(1+1.15e-3_WP*twet)*airp*(ta-twet)
!        specific humidity in kg/kg
         qa = const06*ea/(airp-0.377_WP*ea)
      case (3)  ! Specific humidity from dew point temperature
         if (hum .lt. 100._WP) then
            dew = hum
            dew_k = hum + kelvin
         else
            dew = hum - kelvin
            dew_k = hum
         end if
         ea = a1 +dew*(a2+dew*(a3+dew*(a4+dew*(a5+dew*(a6+dew*a7)))))
         ea = ea * 100.0_WP ! Conversion millibar --> Pascal
         qa = const06*ea/(airp-0.377_WP*ea)
      case (4)
!        specific humidity in kg/kg is given
         qa = hum
!        actual water vapor pressure in Pascal
         ea = qa *airp/(const06+0.378_WP*qa)
      case default
         stop 'not a valid hum_method bulk_fluxes()'
   end select

   rhoa = airp/(rgas*(ta+kelvin)*(1.0_WP+const06*qa))

   return
   end subroutine humidity
!EOC

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------


   function psi(iflag, ZoL)
!=======================================================================
!                                                                      !
!  This function evaluates the stability function, PSI, for wind       !
!  speed (iflag=1) or for air temperature and moisture (iflag=2)       !
!  profiles as function of the stability parameter, ZoL (z/L).         !
!                                                                      !
!  Reference:                                                          !
!                                                                      !
!    Liu, W.T., K.B. Katsaros, and J.A. Businger, 1979:  Bulk          !
!        parameterization of the air-sea exchange of heat and          !
!        water vapor including the molecular constraints at            !
!        the interface, J. Atmos. Sci, 36, 1722-1735.                  !
!                                                                      !
!=======================================================================
!
      real(wp)  :: psi
!
!  Imported variable declarations.
!
   integer,  intent(in) :: iflag
   real(wp), intent(in) :: ZoL
!
!  Local variable declarations.
!
   real(wp), parameter :: r3 = 1.0/3.0
   real(wp), parameter :: sqr3 = 1.7320508
   real(wp), parameter :: pi=3.141592653589
   real(wp)            :: Fw, chic, chik, psic, psik

!  Initialize for the zero "ZoL" case.
!
   psi=0.0_WP
!
!  Unstable conditions.
!
   if (ZoL .lt. 0.0_WP) then
      chik=(1.0_WP-16.0_WP*ZoL)**0.25_WP
      if (iflag .eq. 1) then
         psik=2.0_WP*LOG(0.5_WP*(1.0_WP+chik))+LOG(0.5_WP*(1.0_WP+chik*chik))-   &
              2.0_WP*ATAN(chik)+ 0.5_WP*pi
      else if (iflag .eq. 2) then
            psik=2.0_WP*LOG(0.5_WP*(1.0_WP+chik*chik))
      end if
!
!  For very unstable conditions, use free-convection (Fairall).
!
      chic=(1.0_WP-12.87_WP*ZoL)**r3
      psic=1.5_WP*LOG(r3*(1.0_WP+chic+chic*chic))-                    &
         sqr3*ATAN((1.0_WP+2.0_WP*chic)/sqr3)+ pi/sqr3
!
!  Match Kansas and free-convection forms with weighting Fw.
!
      Fw=1.0_WP/(1.0_WP+ZoL*ZoL)
      psi=Fw*psik+(1.0_WP-Fw)*psic
!
!  Stable conditions.
!
   else if (ZoL .gt. 0.0_WP) then
      psi=-4.7_WP*ZoL
   end if

   return
   end function psi

!EOC
!-----------------------------------------------------------------------
!Copyright (C) 2007 - Adolf Stips
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: Calculate the long-wave back radiation \label{sec:back-rad}
!
! !INTERFACE:
   subroutine back_radiation(method,dlat,tw,ta,cloud,ea,qa,qb)
!
! !DESCRIPTION:
!
! Here, the long-wave back radiation is calculated by means of one out
! of four methods, which depend on the value given to the parameter
! {\tt method}:
! {\tt method}=1: \cite{Clarketal74},
! {\tt method}=2: \cite{HastenrathLamb78},
! {\tt method}=3: \cite{Bignamietal95},
! {\tt method}=4: \cite{BerliandBerliand52}.
! It should ne noted that the latitude must here be given in degrees.
!
! !USES:
!   use airsea_variables, only: emiss,bolz
!   use airsea_variables, only: ea,qa
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: method
   real(wp), intent(in)                :: dlat,tw,ta,cloud
   real(wp), intent(in)                :: qa ! specific humidity (kg/kg)
   real(wp), intent(in)                :: ea ! actual water vapor pressure in Pascal
!
! !OUTPUT PARAMETERS:
   real(wp), intent(out)               :: qb
!
! !REVISION HISTORY:
!  Original author(s): Adols Stips, Hans Burchard & Karsten Bolding
!   real(wp), public, parameter         :: cpw=3985.
   real(wp), parameter         :: emiss=0.97
   real(wp), parameter         :: bolz=5.67e-8

! !LOCAL VARIABLES:

   integer, parameter   :: clark=1      ! Clark et al, 1974
   integer, parameter   :: hastenrath=2 ! Hastenrath and Lamb, 1978
   integer, parameter   :: bignami=3    ! Bignami et al., 1995 - Medsea
   integer, parameter   :: berliand=4   ! Berliand and Berliand, 1952 - ROMS


   real(wp), parameter, dimension(91)  :: cloud_correction_factor = (/ &
     0.497202,     0.501885,     0.506568,     0.511250,     0.515933, &
     0.520616,     0.525299,     0.529982,     0.534665,     0.539348, &
     0.544031,     0.548714,     0.553397,     0.558080,     0.562763, &
     0.567446,     0.572129,     0.576812,     0.581495,     0.586178, &
     0.590861,     0.595544,     0.600227,     0.604910,     0.609593, &
     0.614276,     0.618959,     0.623641,     0.628324,     0.633007, &
     0.637690,     0.642373,     0.647056,     0.651739,     0.656422, &
     0.661105,     0.665788,     0.670471,     0.675154,     0.679837, &
     0.684520,     0.689203,     0.693886,     0.698569,     0.703252, &
     0.707935,     0.712618,     0.717301,     0.721984,     0.726667, &
     0.731350,     0.736032,     0.740715,     0.745398,     0.750081, &
     0.754764,     0.759447,     0.764130,     0.768813,     0.773496, &
     0.778179,     0.782862,     0.787545,     0.792228,     0.796911, &
     0.801594,     0.806277,     0.810960,     0.815643,     0.820326, &
     0.825009,     0.829692,     0.834375,     0.839058,     0.843741, &
     0.848423,     0.853106,     0.857789,     0.862472,     0.867155, &
     0.871838,     0.876521,     0.881204,     0.885887,     0.890570, &
     0.895253,     0.899936,     0.904619,     0.909302,     0.913985, &
     0.918668 /)

   real(wp)                  :: ccf
   real(wp)                  :: x1,x2,x3
!
!EOP
!-----------------------------------------------------------------------
!BOC
!  calculate cloud correction factor,fortran counts from 1 !
   ccf= cloud_correction_factor(nint(abs(dlat))+1)

   select case(method)
      case(clark)
!        Clark et al. (1974) formula.
!        unit of ea is Pascal, must hPa
!        Black body defect term, clouds, water vapor correction
         x1=(1.0_WP-ccf*cloud*cloud)*(tw**4)
         x2=(0.39_WP-0.05_WP*sqrt(ea*0.01_WP))
!        temperature jump term
         x3=4.0_WP*(tw**3)*(tw-ta)
         qb=-emiss*bolz*(x1*x2+x3)
      case(hastenrath) ! qa in g(water)/kg(wet air)
!        Hastenrath and Lamb (1978) formula.
         x1=(1.0_WP-ccf*cloud*cloud)*(tw**4)
         x2=(0.39_WP-0.056_WP*sqrt(1000.0_WP*qa))
         x3=4.0_WP*(tw**3)*(tw-ta)
         qb=-emiss*bolz*(x1*x2+x3)
      case(bignami)
!        Bignami et al. (1995) formula (Med Sea).
!        unit of ea is Pascal, must hPa
         ccf=0.1762_WP
         x1=(1.0_WP+ccf*cloud*cloud)*ta**4
         x2=(0.653_WP+0.00535_WP*(ea*0.01_WP))
         x3= emiss*(tw**4)
         qb=-bolz*(-x1*x2+x3)
      case(berliand)
!        Berliand & Berliand (1952) formula (ROMS).
         x1=(1.0_WP-0.6823_WP*cloud*cloud)*ta**4
         x2=(0.39_WP-0.05_WP*sqrt(0.01_WP*ea))
         x3=4.0_WP*ta**3*(tw-ta)
         qb=-emiss*bolz*(x1*x2+x3)
      case default
   end select

   return
   end subroutine back_radiation
!EOC

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: Calculate the solar zenith angle \label{sec:swr}
!
! !INTERFACE:
   function solar_zenith_angle(yday,hh,dlon,dlat)
!
! !DESCRIPTION:
!  This subroutine calculates the solar zenith angle as being used both
!  in albedo_water() and short_wave_radiation(). The result is in degrees.
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: yday
   real(wp), intent(in)                :: hh
   real(wp), intent(in)                :: dlon,dlat
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding
!
! !LOCAL VARIABLES:
   real(wp), parameter       :: pi=3.14159265358979323846
   real(wp), parameter       :: deg2rad=pi/180.
   real(wp), parameter       :: rad2deg=180./pi

   real(wp)                  :: rlon,rlat
   real(wp)                  :: yrdays
   real(wp)                  :: th0,th02,th03,sundec
   real(wp)                  :: thsun,coszen

   real(wp)                  :: solar_zenith_angle
!EOP
!-----------------------------------------------------------------------
!BOC
!  from now on everything in radians
   rlon = deg2rad*dlon
   rlat = deg2rad*dlat

   yrdays=365.25_WP

   th0 = 2._WP*pi*yday/yrdays
   th02 = 2._WP*th0
   th03 = 3._WP*th0
!  sun declination :
   sundec = 0.006918_WP - 0.399912_WP*cos(th0) + 0.070257_WP*sin(th0)         &
           - 0.006758_WP*cos(th02) + 0.000907_WP*sin(th02)                 &
           - 0.002697_WP*cos(th03) + 0.001480_WP*sin(th03)
!  sun hour angle :
   thsun = (hh-12._WP)*15._WP*deg2rad + rlon

!  cosine of the solar zenith angle :
   coszen =sin(rlat)*sin(sundec)+cos(rlat)*cos(sundec)*cos(thsun)
   if (coszen .lt. 0.0_wp) coszen = 0.0_wp

   solar_zenith_angle = rad2deg*acos(coszen)

   return
   end function solar_zenith_angle
!EOC

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: Calculate the short--wave radiation \label{sec:swr}
!
! !INTERFACE:
   function short_wave_radiation(zenith_angle,yday,dlon,dlat,cloud)
!
! !DESCRIPTION:
!  This subroutine calculates the short--wave net radiation based on
!  solar zenith angle, year day, longitude, latitude, and fractional cloud cover.
!  No corrections for albedo - must be done by calls to albedo\_water() and
!  if ice is included albedo\_ice().
!  The basic formula for the short-wave radiation at the surface, $Q_s$,
!  has been taken from \cite{RosatiMiyacoda88}, who adapted the work
!  of \cite{Reed77} and \cite{SimpsonPaulson99}:
!
!  \begin{equation}
!  Q_s=Q_{tot} (1-0.62 C + 0.0019 \beta) (1-\alpha),
!  \end{equation}
!
!  with the total radiation reaching the surface under clear skies,
!  $Q_{tot}$, the fractional cloud cover, $C$, the solar noon altitude,
!  $\beta$, and the albedo, $\alpha$.
!  This piece of code has been taken the MOM-I (Modular Ocean Model)
!  version at the INGV (Istituto Nazionale di Geofisica e Vulcanologia,
!  see {\tt http://www.bo.ingv.it/}).
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   real(wp), intent(in)                :: zenith_angle
   integer, intent(in)                 :: yday
   real(wp), intent(in)                :: dlon,dlat
   real(wp), intent(in)                :: cloud
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding & Hans Burchard
!
! !LOCAL VARIABLES:
   real(wp), parameter       :: pi=3.14159265358979323846
   real(wp), parameter       :: deg2rad=pi/180.
   real(wp), parameter       :: rad2deg=180./pi

   real(wp), parameter       :: solar=1350.
   real(wp), parameter       :: eclips=23.439*deg2rad
   real(wp), parameter       :: tau=0.7
   real(wp), parameter       :: aozone=0.09

   real(wp)                  :: coszen,sunbet
   real(wp)                  :: qatten,qzer,qdir,qdiff,qtot,qshort
   real(wp)                  :: rlon,rlat,eqnx
   real(wp)                  :: yrdays

   real(wp)                  :: short_wave_radiation
!EOP
!-----------------------------------------------------------------------
!BOC
   coszen = cos(deg2rad*zenith_angle)
   if (coszen .le. 0.0_WP) then
      coszen = 0.0_WP
      qatten = 0.0_WP
   else
      qatten = tau**(1.0_wp/coszen)
   end if

   qzer  = coszen * solar
   qdir  = qzer * qatten
   qdiff = ((1.0_wp-aozone)*qzer - qdir) * 0.5_WP
   qtot  =  qdir + qdiff

!  from now on everything in radians
   rlon = deg2rad*dlon
   rlat = deg2rad*dlat

   yrdays=365._WP
   eqnx = (yday-81._WP)/yrdays*2._WP*pi
!  sin of the solar noon altitude in radians :
   sunbet=sin(rlat)*sin(eclips*sin(eqnx))+cos(rlat)*cos(eclips*sin(eqnx))
!  solar noon altitude in degrees :
   sunbet = asin(sunbet)*rad2deg

!  radiation as from Reed(1977), Simpson and Paulson(1979)
!  calculates SHORT WAVE FLUX ( watt/m*m )
!  Rosati,Miyakoda 1988 ; eq. 3.8
!  clouds from COADS perpetual data set
!#if 1
   qshort  = qtot*(1.0_WP-0.62_WP*cloud + .0019_WP*sunbet)
   if(qshort .gt. qtot ) then
      qshort  = qtot
   end if
!#else
!  original implementation
!   if(cloud .lt. 0.3) then
!      qshort  = qtot
!   else
!      qshort  = qtot*(1-0.62*cloud + 0.0019*sunbet)
!   endif
!#endif
   short_wave_radiation = qshort

   return
   end function short_wave_radiation
!EOC

    !___________________________________________________________________________
    ! make inserted string all in lower case and kick out weired mystery characters
    ! --> replaces 'space', '-' character with '_'
    function lowercase(string)
        implicit none
        character(len=:),allocatable :: lowercase, aux_string
        character(len=*)             :: string 
        character(len=48)            :: aux_string_end=''
        character(len=26)            :: cap  ='ABCDEFGHIJKLMNOPQRSTUVWXYZ'
        character(len=38)            :: small='abcdefghijklmnopqrstuvwxyz1234567890_'
        character(len=2)             :: replace='- '
        integer                      :: i, i1 ,pos_c, pos_s, pos_r
        i1 = 0
        aux_string = trim(string)
        do i=1,len_trim(aux_string)
            pos_c = index(cap,string(i:i))
            pos_s = index(small,string(i:i))
            pos_r = index(replace,string(i:i))
            ! there is problem in the JRA55 calendar attribut string, at the end of
            ! that string there is a character which is not seeable, which is no letter 
            ! and also no whitespace and which can not be removed with trim() --> 
            ! to get rid of that lowercase will use only character that are 
            ! found in either cap or small otherwise the sring comparison fails 
            if (pos_r .ne. 0) then
                ! replaces 'space', '-' character with '_'
                i1=i1+1
                aux_string_end(i1:i1)='_'
            elseif (pos_c .ne. 0 .and. pos_s .eq. 0) then
                i1=i1+1
                aux_string_end(i1:i1)=small(pos_c:pos_c)
            elseif (pos_c .eq. 0 .and. pos_s .ne. 0) then
                i1=i1+1
                aux_string_end(i1:i1)=aux_string(i:i)
            end if
        end do
        lowercase=trim(aux_string_end)
        return
    end function lowercase 

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------

END MODULE g_sbf
