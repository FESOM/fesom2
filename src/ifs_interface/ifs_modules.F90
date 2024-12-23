#define __MYFILE__ 'ifs_modules.F90'
#define key_mpp_mpi
! Set of modules needed by the interface to IFS.
!
! -Original code by Kristian Mogensen, ECMWF.

MODULE par_kind
  IMPLICIT NONE
  INTEGER, PUBLIC, PARAMETER ::          &  !: Floating point section
       sp = SELECTED_REAL_KIND( 6, 37),  &  !: single precision (real 4)
       dp = SELECTED_REAL_KIND(12,307),  &  !: double precision (real 8)
       wpIFS = SELECTED_REAL_KIND(12,307),  &  !: double precision (real 8)
       ik = SELECTED_INT_KIND(6)  , &          !: integer precision 
   jptime = SELECTED_INT_KIND(18)  ! used for time, can handle numbers up to 10^18
END MODULE par_kind

MODULE nctools

   ! Utility subroutines for netCDF access
   ! Modified    : MAB (nf90, handle_error, LINE&FILE)
   ! Modifled    : KSM (new shorter name)

   USE netcdf

   PUBLIC ldebug_netcdf, nchdlerr
   LOGICAL :: ldebug_netcdf = .FALSE.  ! Debug switch for netcdf

CONTAINS

   SUBROUTINE nchdlerr(status,lineno,filename)
      
      ! Error handler for netCDF access
      IMPLICIT NONE


      INTEGER :: status ! netCDF return status
      INTEGER :: lineno ! Line number (usually obtained from 
                        !  preprocessing __LINE__,__MYFILE__)
      CHARACTER(len=*),OPTIONAL :: filename

      IF (status/=nf90_noerr) THEN
         WRITE(*,*)'Netcdf error, code ',status
         IF (PRESENT(filename)) THEN
            WRITE(*,*)'In file ',filename,' in line ',lineno
         ELSE
            WRITE(*,*)'In line ',lineno
         END IF
         WRITE(*,'(2A)')' Error message : ',nf90_strerror(status)
         CALL abort
      ENDIF

   END SUBROUTINE nchdlerr

!----------------------------------------------------------------------
END MODULE nctools

MODULE scrippar
   INTEGER, PARAMETER :: scripdp = SELECTED_REAL_KIND(12,307)
   INTEGER, PARAMETER :: scriplen = 80
END MODULE scrippar

MODULE scripgrid

   USE nctools
   USE scrippar

   IMPLICIT NONE

   TYPE scripgridtype
      INTEGER :: grid_size
      INTEGER :: grid_corners
      INTEGER :: grid_rank
      INTEGER, ALLOCATABLE, DIMENSION(:) :: grid_dims
      REAL(scripdp), ALLOCATABLE, DIMENSION(:) :: grid_center_lat
      REAL(scripdp), ALLOCATABLE, DIMENSION(:) :: grid_center_lon
      INTEGER, ALLOCATABLE, DIMENSION(:) :: grid_imask
      REAL(scripdp), ALLOCATABLE, DIMENSION(:,:) :: grid_corner_lat
      REAL(scripdp), ALLOCATABLE, DIMENSION(:,:) :: grid_corner_lon
      CHARACTER(len=scriplen) :: grid_center_lat_units
      CHARACTER(len=scriplen) :: grid_center_lon_units
      CHARACTER(len=scriplen) :: grid_imask_units
      CHARACTER(len=scriplen) :: grid_corner_lat_units
      CHARACTER(len=scriplen) :: grid_corner_lon_units
      CHARACTER(len=scriplen) :: title
   END TYPE scripgridtype
   
CONTAINS

   SUBROUTINE scripgrid_read( cdfilename, grid )

      CHARACTER(len=*) :: cdfilename
      TYPE(scripgridtype) :: grid

      INTEGER :: ncid, dimid, varid

      CALL scripgrid_init(grid)

      CALL nchdlerr(nf90_open(TRIM(cdfilename),nf90_nowrite,ncid),&
         &          __LINE__,__MYFILE__)
      
      CALL nchdlerr(nf90_inq_dimid(ncid,'grid_size',dimid),&
         &          __LINE__,__MYFILE__)
      CALL nchdlerr(nf90_inquire_dimension(ncid,dimid,&
         &                                 len=grid%grid_size),&
         &          __LINE__,__MYFILE__)

      CALL nchdlerr(nf90_inq_dimid(ncid,'grid_corners',dimid),&
         &          __LINE__,__MYFILE__)
      CALL nchdlerr(nf90_inquire_dimension(ncid,dimid,&
         &                                 len=grid%grid_corners),&
         &          __LINE__,__MYFILE__)

      CALL nchdlerr(nf90_inq_dimid(ncid,'grid_rank',dimid),&
         &          __LINE__,__MYFILE__)
      CALL nchdlerr(nf90_inquire_dimension(ncid,dimid,&
         &                                 len=grid%grid_rank),&
         &          __LINE__,__MYFILE__)

      CALL scripgrid_alloc(grid)

      CALL nchdlerr(nf90_inq_varid(ncid,'grid_dims',varid), &
         &         __LINE__,__MYFILE__)
      CALL nchdlerr(nf90_get_var(ncid,varid,grid%grid_dims),&
         &          __LINE__,__MYFILE__)

      CALL nchdlerr(nf90_inq_varid(ncid,'grid_center_lat',varid), &
         &         __LINE__,__MYFILE__)
      CALL nchdlerr(nf90_get_att(ncid,varid,'units',grid%grid_center_lat_units),&
         &          __LINE__,__MYFILE__)
      CALL nchdlerr(nf90_get_var(ncid,varid,grid%grid_center_lat),&
         &          __LINE__,__MYFILE__)

      CALL nchdlerr(nf90_inq_varid(ncid,'grid_center_lon',varid), &
         &         __LINE__,__MYFILE__)
      CALL nchdlerr(nf90_get_att(ncid,varid,'units',grid%grid_center_lon_units),&
         &          __LINE__,__MYFILE__)
      CALL nchdlerr(nf90_get_var(ncid,varid,grid%grid_center_lon),&
         &          __LINE__,__MYFILE__)

      CALL nchdlerr(nf90_inq_varid(ncid,'grid_corner_lat',varid), &
         &         __LINE__,__MYFILE__)
      CALL nchdlerr(nf90_get_att(ncid,varid,'units',grid%grid_corner_lat_units),&
         &          __LINE__,__MYFILE__)
      CALL nchdlerr(nf90_get_var(ncid,varid,grid%grid_corner_lat),&
         &          __LINE__,__MYFILE__)

      CALL nchdlerr(nf90_inq_varid(ncid,'grid_corner_lon',varid), &
         &         __LINE__,__MYFILE__)
      CALL nchdlerr(nf90_get_att(ncid,varid,'units',grid%grid_corner_lon_units),&
         &          __LINE__,__MYFILE__)
      CALL nchdlerr(nf90_get_var(ncid,varid,grid%grid_corner_lon),&
         &          __LINE__,__MYFILE__)

      CALL nchdlerr(nf90_inq_varid(ncid,'grid_imask',varid), &
         &         __LINE__,__MYFILE__)
      CALL nchdlerr(nf90_get_att(ncid,varid,'units',grid%grid_imask_units),&
         &          __LINE__,__MYFILE__)
      CALL nchdlerr(nf90_get_var(ncid,varid,grid%grid_imask),&
         &          __LINE__,__MYFILE__)

      CALL nchdlerr(nf90_get_att(ncid,nf90_global,'title',grid%title),&
         &          __LINE__,__MYFILE__)

      CALL nchdlerr(nf90_close(ncid),__LINE__,__MYFILE__)
       
   END SUBROUTINE scripgrid_read

   SUBROUTINE scripgrid_write( cdgridfile, grid )

      CHARACTER(len=*) :: cdgridfile
      TYPE(scripgridtype) :: grid

      INTEGER :: ncid
      INTEGER :: ioldfill
      INTEGER :: idimsize,idimxsize,idimysize,idimcorners,idimrank
      INTEGER :: idims1rank(1),idims1size(1),idims2(2)
      INTEGER :: iddims,idcentlat,idcentlon,idimask,idcornlat,idcornlon
      INTEGER :: igriddims(2)

      ! Setup netcdf file

      CALL nchdlerr(nf90_create(TRIM(cdgridfile),nf90_clobber,ncid),&
         &          __LINE__,__MYFILE__)

      ! Define dimensions

      CALL nchdlerr(nf90_def_dim(ncid,'grid_size',&
         &                       grid%grid_size,idimsize),&
         &          __LINE__,__MYFILE__)
      CALL nchdlerr(nf90_def_dim(ncid,'grid_corners',&
         &                       grid%grid_corners,idimcorners),&
         &          __LINE__,__MYFILE__)
      CALL nchdlerr(nf90_def_dim(ncid,'grid_rank',&
         &                       grid%grid_rank,idimrank),&
         &          __LINE__,__MYFILE__)

      idims1rank(1) = idimrank

      idims1size(1) = idimsize

      idims2(1) = idimcorners
      idims2(2) = idimsize

      ! Define variables

      CALL nchdlerr(nf90_def_var(ncid,'grid_dims',&
         &                       nf90_int,idims1rank,iddims),&
         &                       __LINE__,__MYFILE__)

      CALL nchdlerr(nf90_def_var(ncid,'grid_center_lat',&
         &                       nf90_double,idims1size,idcentlat),&
         &                       __LINE__,__MYFILE__)
      CALL nchdlerr(nf90_put_att(ncid,idcentlat,'units',&
         &                       grid%grid_center_lat_units),&
         &          __LINE__,__MYFILE__)

      CALL nchdlerr(nf90_def_var(ncid,'grid_center_lon',&
         &                       nf90_double,idims1size,idcentlon),&
         &                       __LINE__,__MYFILE__)
      CALL nchdlerr(nf90_put_att(ncid,idcentlon,'units',&
         &                       grid%grid_center_lon_units),&
         &          __LINE__,__MYFILE__)

      CALL nchdlerr(nf90_def_var(ncid,'grid_imask',&
         &                       nf90_int,idims1size,idimask),&
         &                       __LINE__,__MYFILE__)
      CALL nchdlerr(nf90_put_att(ncid,idimask,'units',&
         &                       grid%grid_imask_units),&
         &          __LINE__,__MYFILE__)

      CALL nchdlerr(nf90_def_var(ncid,'grid_corner_lat',&
         &                       nf90_double,idims2,idcornlat),&
         &                       __LINE__,__MYFILE__)
      CALL nchdlerr(nf90_put_att(ncid,idcornlat,'units',&
         &                       grid%grid_corner_lat_units),&
         &          __LINE__,__MYFILE__)

      CALL nchdlerr(nf90_def_var(ncid,'grid_corner_lon',&
         &                       nf90_double,idims2,idcornlon),&
         &                       __LINE__,__MYFILE__)
      CALL nchdlerr(nf90_put_att(ncid,idcornlon,'units',&
         &                       grid%grid_corner_lon_units),&
         &          __LINE__,__MYFILE__)

      CALL nchdlerr(nf90_put_att(ncid,nf90_global,'title',&
         &                       TRIM(grid%title)),&
         &          __LINE__,__MYFILE__)

      ! End of netCDF definition phase
                                 
      CALL nchdlerr(nf90_enddef(ncid),__LINE__,__MYFILE__)

      ! Write variables

      
      CALL nchdlerr(nf90_put_var(ncid,iddims,grid%grid_dims),&
         &          __LINE__,__MYFILE__)
      
      CALL nchdlerr(nf90_put_var(ncid,idcentlat,&
         &                       grid%grid_center_lat),&
         &          __LINE__,__MYFILE__)
      
      CALL nchdlerr(nf90_put_var(ncid,idcentlon,&
         &                       grid%grid_center_lon),&
         &          __LINE__,__MYFILE__)
      
      CALL nchdlerr(nf90_put_var(ncid,idimask,&
         &                       grid%grid_imask), &
         &          __LINE__,__MYFILE__)
      
      CALL nchdlerr(nf90_put_var(ncid,idcornlat,&
         &                       grid%grid_corner_lat),&
         &          __LINE__,__MYFILE__)

      CALL nchdlerr(nf90_put_var(ncid,idcornlon,&
         &                       grid%grid_corner_lon),&
         &          __LINE__,__MYFILE__)
      
      ! Close file
      
      CALL nchdlerr(nf90_close(ncid),__LINE__,__MYFILE__)
      
   END SUBROUTINE scripgrid_write

   SUBROUTINE scripgrid_init( grid )

      TYPE(scripgridtype) :: grid

      grid%grid_size=0
      grid%grid_corners=0
      grid%grid_rank=0
      grid%grid_center_lat_units=''
      grid%grid_center_lon_units=''
      grid%grid_imask_units=''
      grid%grid_corner_lat_units=''
      grid%grid_corner_lon_units=''
      grid%title=''

   END SUBROUTINE scripgrid_init

   SUBROUTINE scripgrid_alloc( grid )

      TYPE(scripgridtype) :: grid 

      IF ( (grid%grid_size == 0) .OR. &
         & (grid%grid_corners == 0) .OR. &
         & (grid%grid_rank == 0) ) THEN
         WRITE(*,*)'scripgridtype not initialized'
         CALL abort
      ENDIF
     
      ALLOCATE( &
         & grid%grid_dims(grid%grid_rank), &
         & grid%grid_center_lat(grid%grid_size), &
         & grid%grid_center_lon(grid%grid_size), &         
         & grid%grid_corner_lat(grid%grid_corners, grid%grid_size), &
         & grid%grid_corner_lon(grid%grid_corners, grid%grid_size), &
         & grid%grid_imask(grid%grid_size) &
         & )
      
   END SUBROUTINE scripgrid_alloc
   
   SUBROUTINE scripgrid_dealloc( grid )

      TYPE(scripgridtype) :: grid 

      DEALLOCATE( &
         & grid%grid_dims, &
         & grid%grid_center_lat, &
         & grid%grid_center_lon, &
         & grid%grid_corner_lat, &
         & grid%grid_corner_lon, &
         & grid%grid_imask &
         & )
      
   END SUBROUTINE scripgrid_dealloc

END MODULE scripgrid

MODULE scripremap

#if defined key_mpp_mpi
   USE mpi
#endif
   USE nctools
   USE scrippar
   USE scripgrid

   IMPLICIT NONE

   TYPE scripremaptype
      INTEGER :: num_links
      INTEGER :: num_wgts
      TYPE(scripgridtype) :: src
      TYPE(scripgridtype) :: dst
      REAL(scripdp), ALLOCATABLE, DIMENSION(:) :: src_grid_area
      REAL(scripdp), ALLOCATABLE, DIMENSION(:) :: dst_grid_area
      REAL(scripdp), ALLOCATABLE, DIMENSION(:) :: src_grid_frac
      REAL(scripdp), ALLOCATABLE, DIMENSION(:) :: dst_grid_frac
      INTEGER, ALLOCATABLE, DIMENSION(:) :: src_address
      INTEGER, ALLOCATABLE, DIMENSION(:) :: dst_address
      REAL(scripdp), ALLOCATABLE, DIMENSION(:,:) :: remap_matrix
      CHARACTER(len=scriplen) :: src_grid_area_units
      CHARACTER(len=scriplen) :: dst_grid_area_units
      CHARACTER(len=scriplen) :: src_grid_frac_units
      CHARACTER(len=scriplen) :: dst_grid_frac_units
      CHARACTER(len=scriplen) :: title
      CHARACTER(len=scriplen) :: normalization
      CHARACTER(len=scriplen) :: map_method
      CHARACTER(len=scriplen) :: history
      CHARACTER(len=scriplen) :: conventions
   END TYPE scripremaptype

CONTAINS

   SUBROUTINE scripremap_read_work(cdfilename,remap)

      CHARACTER(len=*) :: cdfilename 
      TYPE(scripremaptype) :: remap

      INTEGER :: ncid, dimid, varid
      LOGICAL :: lcorners

      lcorners=.TRUE.
      
      CALL scripremap_init(remap)

      CALL nchdlerr(nf90_open(TRIM(cdfilename),nf90_nowrite,ncid),&
         &          __LINE__,__MYFILE__)
      
      CALL nchdlerr(nf90_inq_dimid(ncid,'src_grid_size',dimid),&
         &          __LINE__,__MYFILE__)
      CALL nchdlerr(nf90_inquire_dimension(ncid,dimid,&
         &                                 len=remap%src%grid_size),&
         &          __LINE__,__MYFILE__)

      CALL nchdlerr(nf90_inq_dimid(ncid,'dst_grid_size',dimid),&
         &          __LINE__,__MYFILE__)
      CALL nchdlerr(nf90_inquire_dimension(ncid,dimid,&
         &                                 len=remap%dst%grid_size),&
         &          __LINE__,__MYFILE__)

      
      IF (nf90_inq_dimid(ncid,'src_grid_corners',dimid)==nf90_noerr) THEN
         CALL nchdlerr(nf90_inquire_dimension(ncid,dimid,&
            &                                 len=remap%src%grid_corners),&
            &          __LINE__,__MYFILE__)
      ELSE
         lcorners=.FALSE.
         remap%src%grid_corners=1
      ENDIF

      IF (lcorners) THEN
         CALL nchdlerr(nf90_inq_dimid(ncid,'dst_grid_corners',dimid),&
            &          __LINE__,__MYFILE__)
         CALL nchdlerr(nf90_inquire_dimension(ncid,dimid,&
            &                                 len=remap%dst%grid_corners),&
            &          __LINE__,__MYFILE__)
      ELSE
         remap%dst%grid_corners=1
      ENDIF

      CALL nchdlerr(nf90_inq_dimid(ncid,'src_grid_rank',dimid),&
         &          __LINE__,__MYFILE__)
      CALL nchdlerr(nf90_inquire_dimension(ncid,dimid,&
         &                                 len=remap%src%grid_rank),&
         &          __LINE__,__MYFILE__)

      CALL nchdlerr(nf90_inq_dimid(ncid,'dst_grid_rank',dimid),&
         &          __LINE__,__MYFILE__)
      CALL nchdlerr(nf90_inquire_dimension(ncid,dimid,&
         &                                 len=remap%dst%grid_rank),&
         &          __LINE__,__MYFILE__)

      CALL nchdlerr(nf90_inq_dimid(ncid,'num_links',dimid),&
         &          __LINE__,__MYFILE__)
      CALL nchdlerr(nf90_inquire_dimension(ncid,dimid,&
         &                                 len=remap%num_links),&
         &          __LINE__,__MYFILE__)

      CALL nchdlerr(nf90_inq_dimid(ncid,'num_wgts',dimid),&
         &          __LINE__,__MYFILE__)
      CALL nchdlerr(nf90_inquire_dimension(ncid,dimid,&
         &                                 len=remap%num_wgts),&
         &          __LINE__,__MYFILE__)

      CALL scripremap_alloc(remap)
      
      CALL nchdlerr(nf90_inq_varid(ncid,'src_grid_dims',varid), &
         &         __LINE__,__MYFILE__)
      CALL nchdlerr(nf90_get_var(ncid,varid,remap%src%grid_dims),&
         &          __LINE__,__MYFILE__)

      CALL nchdlerr(nf90_inq_varid(ncid,'dst_grid_dims',varid), &
         &         __LINE__,__MYFILE__)
      CALL nchdlerr(nf90_get_var(ncid,varid,remap%dst%grid_dims),&
         &          __LINE__,__MYFILE__)

      CALL nchdlerr(nf90_inq_varid(ncid,'src_grid_center_lat',varid), &
         &         __LINE__,__MYFILE__)
      CALL nchdlerr(nf90_get_att(ncid,varid,'units',remap%src%grid_center_lat_units),&
         &          __LINE__,__MYFILE__)
      CALL nchdlerr(nf90_get_var(ncid,varid,remap%src%grid_center_lat),&
         &          __LINE__,__MYFILE__)

      CALL nchdlerr(nf90_inq_varid(ncid,'dst_grid_center_lat',varid), &
         &         __LINE__,__MYFILE__)
      CALL nchdlerr(nf90_get_att(ncid,varid,'units',remap%dst%grid_center_lat_units),&
         &          __LINE__,__MYFILE__)
      CALL nchdlerr(nf90_get_var(ncid,varid,remap%dst%grid_center_lat),&
         &          __LINE__,__MYFILE__)

      CALL nchdlerr(nf90_inq_varid(ncid,'src_grid_center_lon',varid), &
         &         __LINE__,__MYFILE__)
      CALL nchdlerr(nf90_get_att(ncid,varid,'units',remap%src%grid_center_lon_units),&
         &          __LINE__,__MYFILE__)
      CALL nchdlerr(nf90_get_var(ncid,varid,remap%src%grid_center_lon),&
         &          __LINE__,__MYFILE__)

      CALL nchdlerr(nf90_inq_varid(ncid,'dst_grid_center_lon',varid), &
         &         __LINE__,__MYFILE__)
      CALL nchdlerr(nf90_get_att(ncid,varid,'units',remap%dst%grid_center_lon_units),&
         &          __LINE__,__MYFILE__)
      CALL nchdlerr(nf90_get_var(ncid,varid,remap%dst%grid_center_lon),&
         &          __LINE__,__MYFILE__)

      IF (lcorners) THEN

         CALL nchdlerr(nf90_inq_varid(ncid,'src_grid_corner_lat',varid), &
            &         __LINE__,__MYFILE__)
         CALL nchdlerr(nf90_get_att(ncid,varid,'units',remap%src%grid_corner_lat_units),&
            &          __LINE__,__MYFILE__)
         CALL nchdlerr(nf90_get_var(ncid,varid,remap%src%grid_corner_lat),&
            &          __LINE__,__MYFILE__)

         CALL nchdlerr(nf90_inq_varid(ncid,'src_grid_corner_lon',varid), &
            &         __LINE__,__MYFILE__)
         CALL nchdlerr(nf90_get_att(ncid,varid,'units',remap%src%grid_corner_lon_units),&
            &          __LINE__,__MYFILE__)
         CALL nchdlerr(nf90_get_var(ncid,varid,remap%src%grid_corner_lon),&
            &          __LINE__,__MYFILE__)

         CALL nchdlerr(nf90_inq_varid(ncid,'dst_grid_corner_lat',varid), &
            &         __LINE__,__MYFILE__)
         CALL nchdlerr(nf90_get_att(ncid,varid,'units',remap%dst%grid_corner_lat_units),&
            &          __LINE__,__MYFILE__)
         CALL nchdlerr(nf90_get_var(ncid,varid,remap%dst%grid_corner_lat),&
            &          __LINE__,__MYFILE__)
         
         CALL nchdlerr(nf90_inq_varid(ncid,'dst_grid_corner_lon',varid), &
            &         __LINE__,__MYFILE__)
         CALL nchdlerr(nf90_get_att(ncid,varid,'units',remap%dst%grid_corner_lon_units),&
            &          __LINE__,__MYFILE__)
         CALL nchdlerr(nf90_get_var(ncid,varid,remap%dst%grid_corner_lon),&
            &          __LINE__,__MYFILE__)

      ELSE
         
         remap%src%grid_corner_lat(:,:) = 0.0
         remap%src%grid_corner_lon(:,:) = 0.0
         remap%dst%grid_corner_lat(:,:) = 0.0
         remap%dst%grid_corner_lon(:,:) = 0.0
         remap%src%grid_corner_lat_units = ''
         remap%src%grid_corner_lon_units = ''
         remap%dst%grid_corner_lat_units = ''
         remap%dst%grid_corner_lon_units = ''

      ENDIF

      CALL nchdlerr(nf90_inq_varid(ncid,'src_grid_imask',varid), &
         &         __LINE__,__MYFILE__)
      CALL nchdlerr(nf90_get_att(ncid,varid,'units',remap%src%grid_imask_units),&
         &          __LINE__,__MYFILE__)
      CALL nchdlerr(nf90_get_var(ncid,varid,remap%src%grid_imask),&
         &          __LINE__,__MYFILE__)

      CALL nchdlerr(nf90_inq_varid(ncid,'dst_grid_imask',varid), &
         &         __LINE__,__MYFILE__)
      CALL nchdlerr(nf90_get_att(ncid,varid,'units',remap%dst%grid_imask_units),&
         &          __LINE__,__MYFILE__)
      CALL nchdlerr(nf90_get_var(ncid,varid,remap%dst%grid_imask),&
         &          __LINE__,__MYFILE__)

      CALL nchdlerr(nf90_inq_varid(ncid,'src_grid_area',varid), &
         &         __LINE__,__MYFILE__)
      CALL nchdlerr(nf90_get_att(ncid,varid,'units',remap%src_grid_area_units),&
         &          __LINE__,__MYFILE__)
      CALL nchdlerr(nf90_get_var(ncid,varid,remap%src_grid_area),&
         &          __LINE__,__MYFILE__)

      CALL nchdlerr(nf90_inq_varid(ncid,'dst_grid_area',varid), &
         &         __LINE__,__MYFILE__)
      CALL nchdlerr(nf90_get_att(ncid,varid,'units',remap%dst_grid_area_units),&
         &          __LINE__,__MYFILE__)
      CALL nchdlerr(nf90_get_var(ncid,varid,remap%dst_grid_area),&
         &          __LINE__,__MYFILE__)

      CALL nchdlerr(nf90_inq_varid(ncid,'src_grid_frac',varid), &
         &         __LINE__,__MYFILE__)
      CALL nchdlerr(nf90_get_att(ncid,varid,'units',remap%src_grid_frac_units),&
         &          __LINE__,__MYFILE__)
      CALL nchdlerr(nf90_get_var(ncid,varid,remap%src_grid_frac),&
         &          __LINE__,__MYFILE__)

      CALL nchdlerr(nf90_inq_varid(ncid,'dst_grid_frac',varid), &
         &         __LINE__,__MYFILE__)
      CALL nchdlerr(nf90_get_att(ncid,varid,'units',remap%dst_grid_frac_units),&
         &          __LINE__,__MYFILE__)
      CALL nchdlerr(nf90_get_var(ncid,varid,remap%dst_grid_frac),&
         &          __LINE__,__MYFILE__)

      CALL nchdlerr(nf90_inq_varid(ncid,'src_address',varid), &
         &         __LINE__,__MYFILE__)
      CALL nchdlerr(nf90_get_var(ncid,varid,remap%src_address),&
         &          __LINE__,__MYFILE__)

      CALL nchdlerr(nf90_inq_varid(ncid,'dst_address',varid), &
         &         __LINE__,__MYFILE__)
      CALL nchdlerr(nf90_get_var(ncid,varid,remap%dst_address),&
         &          __LINE__,__MYFILE__)

      CALL nchdlerr(nf90_inq_varid(ncid,'remap_matrix',varid), &
         &         __LINE__,__MYFILE__)
      CALL nchdlerr(nf90_get_var(ncid,varid,remap%remap_matrix),&
         &          __LINE__,__MYFILE__)

      CALL nchdlerr(nf90_get_att(ncid,nf90_global,'title',remap%title),&
         &          __LINE__,__MYFILE__)
      CALL nchdlerr(nf90_get_att(ncid,nf90_global,'normalization',remap%normalization),&
         &          __LINE__,__MYFILE__)
      CALL nchdlerr(nf90_get_att(ncid,nf90_global,'map_method',remap%map_method),&
         &          __LINE__,__MYFILE__)
      CALL nchdlerr(nf90_get_att(ncid,nf90_global,'history',remap%history),&
         &          __LINE__,__MYFILE__)
      CALL nchdlerr(nf90_get_att(ncid,nf90_global,'conventions',remap%conventions),&
         &          __LINE__,__MYFILE__)
      CALL nchdlerr(nf90_get_att(ncid,nf90_global,'dest_grid',remap%dst%title),&
         &          __LINE__,__MYFILE__)
      CALL nchdlerr(nf90_get_att(ncid,nf90_global,'source_grid',remap%src%title),&
         &          __LINE__,__MYFILE__)

      CALL nchdlerr(nf90_close(ncid),__LINE__,__MYFILE__)

   END SUBROUTINE scripremap_read_work

   SUBROUTINE scripremap_read(cdfilename,remap)

      CHARACTER(len=*) :: cdfilename 
      TYPE(scripremaptype) :: remap

      CALL scripremap_read_work(cdfilename,remap)

   END SUBROUTINE scripremap_read


   SUBROUTINE scripremap_read_sgl(cdfilename,remap,&
      &                           mype,nproc,mycomm,linteronly)

      CHARACTER(len=*) :: cdfilename 
      TYPE(scripremaptype) :: remap
      INTEGER :: mype,nproc,mycomm
      LOGICAL :: linteronly

      INTEGER, DIMENSION(8) :: isizes
      INTEGER :: ierr, ip

      IF (mype==0) THEN
         CALL scripremap_read_work(cdfilename,remap)
#if defined key_mpp_mpi
         isizes(1)=remap%src%grid_size
         isizes(2)=remap%dst%grid_size
         isizes(3)=remap%src%grid_corners
         isizes(4)=remap%dst%grid_corners
         isizes(5)=remap%src%grid_rank
         isizes(6)=remap%dst%grid_rank
         isizes(7)=remap%num_links
         isizes(8)=remap%num_wgts
         CALL mpi_bcast( isizes, 8, mpi_integer, 0, mycomm, ierr)
      ELSE
         CALL mpi_bcast( isizes, 8, mpi_integer, 0, mycomm, ierr)
         CALL scripremap_init(remap)
         remap%src%grid_size=isizes(1)
         remap%dst%grid_size=isizes(2)
         remap%src%grid_corners=isizes(3)
         remap%dst%grid_corners=isizes(4)
         remap%src%grid_rank=isizes(5)
         remap%dst%grid_rank=isizes(6)
         remap%num_links=isizes(7)
         remap%num_wgts=isizes(8)
         CALL scripremap_alloc(remap)
#endif
      ENDIF

#if defined key_mpp_mpi
      
      IF (.NOT.linteronly) THEN

         CALL mpi_bcast( remap%src%grid_dims, remap%src%grid_rank, &
            & mpi_integer, 0, mycomm, ierr )
         CALL mpi_bcast( remap%src%grid_center_lat, remap%src%grid_size, &
            & mpi_double_precision, 0, mycomm, ierr )
         CALL mpi_bcast( remap%src%grid_center_lon, remap%src%grid_size, &
            & mpi_double_precision, 0, mycomm, ierr )
         CALL mpi_bcast( remap%src%grid_corner_lat, remap%src%grid_corners*remap%src%grid_size, &
            & mpi_double_precision, 0, mycomm, ierr )
         CALL mpi_bcast( remap%src%grid_corner_lon, remap%src%grid_corners*remap%src%grid_size, &
            & mpi_double_precision, 0, mycomm, ierr )

         CALL mpi_bcast( remap%dst%grid_dims, remap%dst%grid_rank, &
            & mpi_integer, 0, mycomm, ierr )
         CALL mpi_bcast( remap%dst%grid_center_lat, remap%dst%grid_size, &
            & mpi_double_precision, 0, mycomm, ierr )
         CALL mpi_bcast( remap%dst%grid_center_lon, remap%dst%grid_size, &
            & mpi_double_precision, 0, mycomm, ierr )
         CALL mpi_bcast( remap%dst%grid_corner_lat, remap%dst%grid_corners*remap%dst%grid_size, &
            & mpi_double_precision, 0, mycomm, ierr )
         CALL mpi_bcast( remap%dst%grid_corner_lon, remap%dst%grid_corners*remap%dst%grid_size, &
            & mpi_double_precision, 0, mycomm, ierr )

         CALL mpi_bcast( remap%src_grid_area, remap%src%grid_size, &
         & mpi_double_precision, 0, mycomm, ierr )
         CALL mpi_bcast( remap%dst_grid_area, remap%dst%grid_size, &
            & mpi_double_precision, 0, mycomm, ierr )
         CALL mpi_bcast( remap%src_grid_frac, remap%src%grid_size, &
            & mpi_double_precision, 0, mycomm, ierr )
         CALL mpi_bcast( remap%dst_grid_frac, remap%dst%grid_size, &
            & mpi_double_precision, 0, mycomm, ierr )

         CALL mpi_bcast( remap%src%grid_center_lat_units, scriplen, &
            & mpi_character, 0, mycomm, ierr )
         CALL mpi_bcast( remap%dst%grid_center_lat_units, scriplen, &
            & mpi_character, 0, mycomm, ierr )
         CALL mpi_bcast( remap%src%grid_center_lon_units, scriplen, &
            & mpi_character, 0, mycomm, ierr )
         CALL mpi_bcast( remap%dst%grid_center_lon_units, scriplen, &
            & mpi_character, 0, mycomm, ierr )
         CALL mpi_bcast( remap%src%grid_corner_lat_units, scriplen, &
            & mpi_character, 0, mycomm, ierr )
         CALL mpi_bcast( remap%src%grid_corner_lon_units, scriplen, &
            & mpi_character, 0, mycomm, ierr )
         CALL mpi_bcast( remap%dst%grid_corner_lat_units, scriplen, &
            & mpi_character, 0, mycomm, ierr )
         CALL mpi_bcast( remap%dst%grid_corner_lon_units, scriplen, &
            & mpi_character, 0, mycomm, ierr )
         CALL mpi_bcast( remap%src%grid_imask_units, scriplen, &
            & mpi_character, 0, mycomm, ierr )
         CALL mpi_bcast( remap%dst%grid_imask_units, scriplen, &
            & mpi_character, 0, mycomm, ierr )
         CALL mpi_bcast( remap%src_grid_area_units, scriplen, &
            & mpi_character, 0, mycomm, ierr )
         CALL mpi_bcast( remap%dst_grid_area_units, scriplen, &
            & mpi_character, 0, mycomm, ierr )
         CALL mpi_bcast( remap%src_grid_frac_units, scriplen, &
            & mpi_character, 0, mycomm, ierr )
         CALL mpi_bcast( remap%dst_grid_frac_units, scriplen, &
            & mpi_character, 0, mycomm, ierr )
         CALL mpi_bcast( remap%title, scriplen, &
            & mpi_character, 0, mycomm, ierr )
         CALL mpi_bcast( remap%normalization, scriplen, &
            & mpi_character, 0, mycomm, ierr )
         CALL mpi_bcast( remap%map_method, scriplen, &
            & mpi_character, 0, mycomm, ierr )
         CALL mpi_bcast( remap%history, scriplen, &
            & mpi_character, 0, mycomm, ierr )
         CALL mpi_bcast( remap%conventions, scriplen, &
            & mpi_character, 0, mycomm, ierr )
      ENDIF
      
      CALL mpi_bcast( remap%src_address, remap%num_links, &
         & mpi_integer, 0, mycomm, ierr )
      CALL mpi_bcast( remap%dst_address, remap%num_links, &
         & mpi_integer, 0, mycomm, ierr )
      CALL mpi_bcast( remap%remap_matrix, remap%num_wgts*remap%num_links, &
         & mpi_double_precision, 0, mycomm, ierr )
      CALL mpi_bcast(  remap%src%grid_imask, remap%src%grid_size, &
         & mpi_integer, 0, mycomm, ierr )
      CALL mpi_bcast(  remap%dst%grid_imask, remap%dst%grid_size, &
         & mpi_integer, 0, mycomm, ierr )

#endif
   END SUBROUTINE scripremap_read_sgl

   SUBROUTINE scripremap_write(cdfilename,remap)

      CHARACTER(len=*) :: cdfilename 
      TYPE(scripremaptype) :: remap

      INTEGER :: ncid
      INTEGER :: dimsgs,dimdgs,dimsgc,dimdgc,dimsgr,dimdgr,dimnl,dimnw
      INTEGER :: dims1(1),dims2(2)
      INTEGER :: idsgd,iddgd,idsgea,iddgea,idsgeo,iddgeo
      INTEGER :: idsgoa,idsgoo,iddgoa,iddgoo,idsgim,iddgim,idsgar,iddgar
      INTEGER :: idsgf,iddgf,idsga,iddga,idsa,idda,idrm

      CALL nchdlerr(nf90_create(TRIM(cdfilename),nf90_clobber,ncid), &
         &          __LINE__, __MYFILE__ )
      
      CALL nchdlerr(nf90_def_dim(ncid,'src_grid_size',&
         &                       remap%src%grid_size,dimsgs),&
         &          __LINE__,__MYFILE__)

      CALL nchdlerr(nf90_def_dim(ncid,'dst_grid_size',&
         &                       remap%dst%grid_size,dimdgs),&
         &          __LINE__,__MYFILE__)
      
      CALL nchdlerr(nf90_def_dim(ncid,'src_grid_corners',&
         &                       remap%src%grid_corners,dimsgc),&
         &          __LINE__,__MYFILE__)

      CALL nchdlerr(nf90_def_dim(ncid,'dst_grid_corners',&
         &                       remap%dst%grid_corners,dimdgc),&
         &          __LINE__,__MYFILE__)

      CALL nchdlerr(nf90_def_dim(ncid,'src_grid_rank',&
         &                       remap%src%grid_rank,dimsgr),&
         &          __LINE__,__MYFILE__)

      CALL nchdlerr(nf90_def_dim(ncid,'dst_grid_rank',&
         &                       remap%dst%grid_rank,dimdgr),&
         &          __LINE__,__MYFILE__)

      CALL nchdlerr(nf90_def_dim(ncid,'num_links',&
         &                       remap%num_links,dimnl),&
         &          __LINE__,__MYFILE__)

      CALL nchdlerr(nf90_def_dim(ncid,'num_wgts',&
         &                       remap%num_wgts,dimnw),&
         &          __LINE__,__MYFILE__)

      dims1(1)=dimsgr
      CALL nchdlerr(nf90_def_var(ncid,'src_grid_dims',&
         &                       nf90_int,dims1,idsgd),&
         &          __LINE__,__MYFILE__)

      dims1(1)=dimdgr
      CALL nchdlerr(nf90_def_var(ncid,'dst_grid_dims',&
         &                       nf90_int,dims1,iddgd), &
         &         __LINE__,__MYFILE__)

      dims1(1)=dimsgs
      CALL nchdlerr(nf90_def_var(ncid,'src_grid_center_lat',&
         &                       nf90_double,dims1,idsgea), &
         &         __LINE__,__MYFILE__)

      dims1(1)=dimdgs
      CALL nchdlerr(nf90_def_var(ncid,'dst_grid_center_lat',&
         &                       nf90_double,dims1,iddgea), &
         &         __LINE__,__MYFILE__)

      dims1(1)=dimsgs
      CALL nchdlerr(nf90_def_var(ncid,'src_grid_center_lon',&
         &                       nf90_double,dims1,idsgeo), &
         &         __LINE__,__MYFILE__)

      dims1(1)=dimdgs
      CALL nchdlerr(nf90_def_var(ncid,'dst_grid_center_lon',&
         &                       nf90_double,dims1,iddgeo), &
         &         __LINE__,__MYFILE__)

      dims2(1)=dimsgc
      dims2(2)=dimsgs
      CALL nchdlerr(nf90_def_var(ncid,'src_grid_corner_lat',&
         &                       nf90_double,dims2,idsgoa), &
         &         __LINE__,__MYFILE__)

      dims2(1)=dimsgc
      dims2(2)=dimsgs
      CALL nchdlerr(nf90_def_var(ncid,'src_grid_corner_lon',&
         &                       nf90_double,dims2,idsgoo), &
         &         __LINE__,__MYFILE__)
 
      dims2(1)=dimdgc
      dims2(2)=dimdgs
      CALL nchdlerr(nf90_def_var(ncid,'dst_grid_corner_lat',&
         &                       nf90_double,dims2,iddgoa), &
         &         __LINE__,__MYFILE__)

      dims2(1)=dimdgc
      dims2(2)=dimdgs
      CALL nchdlerr(nf90_def_var(ncid,'dst_grid_corner_lon',&
         &                       nf90_double,dims2,iddgoo), &
         &         __LINE__,__MYFILE__) 

      dims1(1)=dimsgs
      CALL nchdlerr(nf90_def_var(ncid,'src_grid_imask',&
         &                       nf90_int,dims1,idsgim), &
         &         __LINE__,__MYFILE__)

      dims1(1)=dimdgs
      CALL nchdlerr(nf90_def_var(ncid,'dst_grid_imask',&
         &                       nf90_int,dims1,iddgim), &
         &         __LINE__,__MYFILE__)

      dims1(1)=dimsgs
      CALL nchdlerr(nf90_def_var(ncid,'src_grid_area',&
         &                       nf90_double,dims1,idsga), &
         &         __LINE__,__MYFILE__)

      dims1(1)=dimdgs
      CALL nchdlerr(nf90_def_var(ncid,'dst_grid_area',&
         &                       nf90_double,dims1,iddga), &
         &         __LINE__,__MYFILE__)

      dims1(1)=dimsgs
      CALL nchdlerr(nf90_def_var(ncid,'src_grid_frac',&
         &                       nf90_double,dims1,idsgf), &
         &         __LINE__,__MYFILE__)

      dims1(1)=dimdgs
      CALL nchdlerr(nf90_def_var(ncid,'dst_grid_frac',&
         &                       nf90_double,dims1,iddgf), &
         &         __LINE__,__MYFILE__)

      dims1(1)=dimnl
      CALL nchdlerr(nf90_def_var(ncid,'src_address',&
         &                       nf90_int,dims1,idsa), &
         &         __LINE__,__MYFILE__)

      dims1(1)=dimnl
      CALL nchdlerr(nf90_def_var(ncid,'dst_address',&
         &                       nf90_int,dims1,idda), &
         &         __LINE__,__MYFILE__)

      dims2(1)=dimnw
      dims2(2)=dimnl
      CALL nchdlerr(nf90_def_var(ncid,'remap_matrix',&
         &                       nf90_double,dims2,idrm), &
         &         __LINE__,__MYFILE__)

      CALL nchdlerr(nf90_put_att(ncid,idsgea,'units',&
         &                       remap%src%grid_center_lat_units),&
         &          __LINE__,__MYFILE__)
      CALL nchdlerr(nf90_put_att(ncid,iddgea,'units',&
         &                       remap%dst%grid_center_lat_units),&
         &          __LINE__,__MYFILE__)
      CALL nchdlerr(nf90_put_att(ncid,idsgeo,'units',&
         &                       remap%src%grid_center_lon_units),&
         &          __LINE__,__MYFILE__)
      CALL nchdlerr(nf90_put_att(ncid,iddgeo,'units',&
         &                       remap%dst%grid_center_lon_units),&
         &          __LINE__,__MYFILE__)
      CALL nchdlerr(nf90_put_att(ncid,idsgoa,'units',&
         &                       remap%src%grid_corner_lat_units),&
         &          __LINE__,__MYFILE__)
      CALL nchdlerr(nf90_put_att(ncid,idsgoo,'units',&
         &                       remap%src%grid_corner_lon_units),&
         &          __LINE__,__MYFILE__)
      CALL nchdlerr(nf90_put_att(ncid,iddgoa,'units',&
         &                       remap%dst%grid_corner_lat_units),&
         &          __LINE__,__MYFILE__)
      CALL nchdlerr(nf90_put_att(ncid,iddgoo,'units',&
         &                       remap%dst%grid_corner_lon_units),&
         &          __LINE__,__MYFILE__)
      CALL nchdlerr(nf90_put_att(ncid,idsgim,'units',&
         &                       remap%src%grid_imask_units),&
         &          __LINE__,__MYFILE__)
      CALL nchdlerr(nf90_put_att(ncid,iddgim,'units',&
         &                       remap%dst%grid_imask_units),&
         &          __LINE__,__MYFILE__)
      CALL nchdlerr(nf90_put_att(ncid,idsga,'units',&
         &                       remap%src_grid_area_units),&
         &          __LINE__,__MYFILE__)
      CALL nchdlerr(nf90_put_att(ncid,iddga,'units',&
         &                       remap%dst_grid_area_units),&
         &          __LINE__,__MYFILE__)
      CALL nchdlerr(nf90_put_att(ncid,idsgf,'units',&
         &                       remap%src_grid_frac_units),&
         &          __LINE__,__MYFILE__)
      CALL nchdlerr(nf90_put_att(ncid,iddgf,'units',&
         &                       remap%dst_grid_frac_units),&
         &          __LINE__,__MYFILE__)

      CALL nchdlerr(nf90_put_att(ncid,nf90_global,'title',&
         &                       remap%title),&
         &          __LINE__,__MYFILE__)
      CALL nchdlerr(nf90_put_att(ncid,nf90_global,'normalization',&
         &                       remap%normalization),&
         &          __LINE__,__MYFILE__)
      CALL nchdlerr(nf90_put_att(ncid,nf90_global,'map_method',&
         &                       remap%map_method),&
         &          __LINE__,__MYFILE__)
      CALL nchdlerr(nf90_put_att(ncid,nf90_global,'history',&
         &                       remap%history),&
         &          __LINE__,__MYFILE__)
      CALL nchdlerr(nf90_put_att(ncid,nf90_global,'conventions',&
         &                       remap%conventions),&
         &          __LINE__,__MYFILE__)
      CALL nchdlerr(nf90_put_att(ncid,nf90_global,'dest_grid',&
         &                       remap%dst%title),&
         &          __LINE__,__MYFILE__)
      CALL nchdlerr(nf90_put_att(ncid,nf90_global,'source_grid',&
         &                       remap%src%title),&
         &          __LINE__,__MYFILE__)

      CALL nchdlerr(nf90_enddef(ncid),__LINE__,__MYFILE__)
      
      CALL nchdlerr(nf90_put_var(ncid,idsgd,remap%src%grid_dims),&
         &          __LINE__,__MYFILE__)

      CALL nchdlerr(nf90_put_var(ncid,iddgd,remap%dst%grid_dims),&
         &          __LINE__,__MYFILE__)

      CALL nchdlerr(nf90_put_var(ncid,idsgea,remap%src%grid_center_lat),&
         &          __LINE__,__MYFILE__)

      CALL nchdlerr(nf90_put_var(ncid,iddgea,remap%dst%grid_center_lat),&
         &          __LINE__,__MYFILE__)

      CALL nchdlerr(nf90_put_var(ncid,idsgeo,remap%src%grid_center_lon),&
         &          __LINE__,__MYFILE__)

      CALL nchdlerr(nf90_put_var(ncid,iddgeo,remap%dst%grid_center_lon),&
         &          __LINE__,__MYFILE__)

      CALL nchdlerr(nf90_put_var(ncid,idsgoa,remap%src%grid_corner_lat),&
         &          __LINE__,__MYFILE__)
      CALL nchdlerr(nf90_put_var(ncid,idsgoo,remap%src%grid_corner_lon),&
         &          __LINE__,__MYFILE__)
      CALL nchdlerr(nf90_put_var(ncid,iddgoa,remap%dst%grid_corner_lat),&
         &          __LINE__,__MYFILE__)
      CALL nchdlerr(nf90_put_var(ncid,iddgoo,remap%dst%grid_corner_lon),&
         &          __LINE__,__MYFILE__)
      CALL nchdlerr(nf90_put_var(ncid,idsgim,remap%src%grid_imask),&
         &          __LINE__,__MYFILE__)
      CALL nchdlerr(nf90_put_var(ncid,iddgim,remap%dst%grid_imask),&
         &          __LINE__,__MYFILE__)
      CALL nchdlerr(nf90_put_var(ncid,idsga,remap%src_grid_area),&
         &          __LINE__,__MYFILE__)
      CALL nchdlerr(nf90_put_var(ncid,iddga,remap%dst_grid_area),&
         &          __LINE__,__MYFILE__)
      CALL nchdlerr(nf90_put_var(ncid,idsgf,remap%src_grid_frac),&
         &          __LINE__,__MYFILE__)
      CALL nchdlerr(nf90_put_var(ncid,iddgf,remap%dst_grid_frac),&
         &          __LINE__,__MYFILE__)
      CALL nchdlerr(nf90_put_var(ncid,idsa,remap%src_address),&
         &          __LINE__,__MYFILE__)
      CALL nchdlerr(nf90_put_var(ncid,idda,remap%dst_address),&
         &          __LINE__,__MYFILE__)
      CALL nchdlerr(nf90_put_var(ncid,idrm,remap%remap_matrix),&
         &          __LINE__,__MYFILE__)

      CALL nchdlerr(nf90_close(ncid),__LINE__, __MYFILE__ )

   END SUBROUTINE scripremap_write

   SUBROUTINE scripremap_init(remap) 

      TYPE(scripremaptype) :: remap

      CALL scripgrid_init(remap%src)
      CALL scripgrid_init(remap%dst)
      remap%num_links = 0
      remap%num_wgts = 0
      remap%title=''
      remap%normalization=''
      remap%map_method=''
      remap%history=''
      remap%conventions=''
      remap%src_grid_area_units=''
      remap%dst_grid_area_units=''
      remap%src_grid_frac_units=''
      remap%dst_grid_frac_units=''

   END SUBROUTINE scripremap_init

   SUBROUTINE scripremap_alloc(remap) 

      TYPE(scripremaptype) :: remap

      IF ( (remap%num_links == 0) .OR. &
         & (remap%num_wgts == 0) ) THEN
         WRITE(*,*)'scripremaptype not initialized'
         CALL abort
      ENDIF

      CALL scripgrid_alloc(remap%src)
      CALL scripgrid_alloc(remap%dst)
      
      ALLOCATE( &
         & remap%src_grid_area(remap%src%grid_size), &
         & remap%dst_grid_area(remap%dst%grid_size), &
         & remap%src_grid_frac(remap%src%grid_size), &
         & remap%dst_grid_frac(remap%dst%grid_size), &
         & remap%src_address(remap%num_links), &
         & remap%dst_address(remap%num_links), &
         & remap%remap_matrix(remap%num_wgts, remap%num_links) &
         & )

   END SUBROUTINE scripremap_alloc

   SUBROUTINE scripremap_dealloc(remap) 

      TYPE(scripremaptype) :: remap

      DEALLOCATE( &
         & remap%src_grid_area, &
         & remap%dst_grid_area, &
         & remap%src_grid_frac, &
         & remap%dst_grid_frac, &
         & remap%src_address, &
         & remap%dst_address, &
         & remap%remap_matrix &
         & )

      CALL scripgrid_dealloc(remap%src)
      CALL scripgrid_dealloc(remap%dst)

      CALL scripremap_init(remap)

   END SUBROUTINE scripremap_dealloc

END MODULE scripremap

MODULE parinter

#if defined key_mpp_mpi
   USE mpi
#endif
   USE scripremap
   USE scrippar
   USE nctools

   IMPLICIT NONE

   ! Type to contains interpolation information
   ! (like what is in scripremaptype) and message
   ! passing information

   TYPE parinterinfo  
      ! Number of local links
      INTEGER :: num_links
      ! Destination side
      INTEGER, POINTER, DIMENSION(:) :: dst_address
      ! Source addresses and work array
      INTEGER, POINTER, DIMENSION(:) :: src_address
      ! Local remap matrix
      REAL(scripdp), POINTER, DIMENSION(:,:) :: remap_matrix
      ! Message passing information
      ! Array of local addresses for send buffer
      ! packing
      INTEGER, POINTER, DIMENSION(:) :: send_address
      ! Sending bookkeeping
      INTEGER :: nsendtot                
      INTEGER, POINTER, DIMENSION(:) :: nsend,nsdisp
      ! Receiving bookkeeping
      INTEGER :: nrecvtot
      INTEGER, POINTER, DIMENSION(:) :: nrecv,nrdisp
   END TYPE parinterinfo
   LOGICAL, PUBLIC :: lparinterp2p = .TRUE.

CONTAINS

   SUBROUTINE parinter_init( mype, nproc, mpi_comm, &
      & nsrclocpoints, nsrcglopoints, srcmask, srcgloind,  &
      & ndstlocpoints, ndstglopoints, dstmask, dstgloind, &
      & remap, pinfo, lcommout, commoutprefix, iunit )

      ! Setup interpolation based on SCRIP format weights in
      ! remap and the source/destination grids information.

      ! Procedure:

      !   1) A global SCRIP remapping file is read on all processors.
      !   2) Find local destination points in the global grid.
      !   3) Find which processor needs source data and setup buffer
      !      information for sending data.
      !   4) Construct new src remapping for buffer received

      ! All information is stored in the TYPE(parinterinfo) output
      ! data type

      ! Input arguments.

      ! Message passing information
      INTEGER, INTENT(IN) :: mype, nproc, mpi_comm
      ! Source grid local and global number of grid points
      INTEGER, INTENT(IN) :: nsrclocpoints, nsrcglopoints
      ! Source integer mask (0/1) for SCRIP compliance
      INTEGER, INTENT(IN), DIMENSION(nsrclocpoints) :: srcmask
      ! Source global addresses of each local grid point
      INTEGER, INTENT(IN), DIMENSION(nsrclocpoints) :: srcgloind
      ! Destination grid local and global number of grid points
      INTEGER, INTENT(IN) :: ndstlocpoints, ndstglopoints
      ! Destination integer mask (0/1) for SCRIP compliance
      INTEGER, INTENT(IN), DIMENSION(ndstlocpoints) :: dstmask
      ! Destination global addresses of each local grid point
      INTEGER, INTENT(IN), DIMENSION(ndstlocpoints) :: dstgloind
      ! SCRIP remapping data
      TYPE(scripremaptype) :: remap
      ! Switch for output communication patterns
      LOGICAL :: lcommout
      CHARACTER(len=*) :: commoutprefix
      ! Unit to use for output
      INTEGER :: iunit

      ! Output arguments

      ! Interpolation and message passing information
      TYPE(parinterinfo), INTENT(OUT) :: pinfo

      ! Local variable

      ! Variable for glocal <-> local address/pe information
      INTEGER, DIMENSION(nsrcglopoints) :: ilsrcmppmap, ilsrclocind
      INTEGER, DIMENSION(nsrcglopoints) :: igsrcmppmap, igsrclocind
      INTEGER, DIMENSION(ndstglopoints) :: ildstmppmap, ildstlocind
      INTEGER, DIMENSION(ndstglopoints) :: igdstmppmap, igdstlocind
      INTEGER, DIMENSION(nsrcglopoints) :: isrcpe,isrcpetmp
      INTEGER, DIMENSION(nsrcglopoints) :: isrcaddtmp
      INTEGER, DIMENSION(0:nproc-1)     :: isrcoffset
      INTEGER, DIMENSION(nproc)         :: isrcno, isrcoff, isrccur
      INTEGER, DIMENSION(nproc)         :: ircvoff, ircvcur
      INTEGER, DIMENSION(:), ALLOCATABLE :: isrctot, ircvtot

      ! Misc variable
      INTEGER :: i,n,pe
      INTEGER :: istatus
      CHARACTER(len=256) :: cdfile

      ! Check that masks are consistent.

      ! Remark: More consistency tests between remapping information
      ! and input argument could be code, but for now we settle
      ! for checking the masks.

      ! Source grid

      DO i=1,nsrclocpoints
         IF (srcmask(i)/=remap%src%grid_imask(srcgloind(i))) THEN
            WRITE(iunit,*)'Source imask is inconsistent at '
            WRITE(iunit,*)'global index = ',srcgloind(i)
            WRITE(iunit,*)'Source mask  = ',srcmask(i)
            WRITE(iunit,*)'Remap  mask  = ',remap%src%grid_imask(srcgloind(i))
            WRITE(iunit,*)'Latitude     = ',remap%src%grid_center_lat(srcgloind(i))
            WRITE(iunit,*)'Longitude    = ',remap%src%grid_center_lon(srcgloind(i))
            CALL flush(iunit)
            CALL abort
         ENDIF
      ENDDO

      ! Destination grid
      
      DO i=1,ndstlocpoints
         IF (dstmask(i)/=remap%dst%grid_imask(dstgloind(i))) THEN
            WRITE(iunit,*)'Destination imask is inconsistent at '
            WRITE(iunit,*)'global index = ',dstgloind(i)
            WRITE(iunit,*)'Destin mask  = ',dstmask(i)
            WRITE(iunit,*)'Remap  mask  = ',remap%dst%grid_imask(dstgloind(i))
            WRITE(iunit,*)'Latitude     = ',remap%dst%grid_center_lat(dstgloind(i))
            WRITE(iunit,*)'Longitude    = ',remap%dst%grid_center_lon(dstgloind(i))
            CALL flush(iunit)
            CALL abort
         ENDIF
      ENDDO

      ! Setup global to local and vice versa mappings.

      ilsrcmppmap(:)=-1
      ilsrclocind(:)=0
      ildstmppmap(:)=-1
      ildstlocind(:)=0      
      
      DO i=1,nsrclocpoints
         ilsrcmppmap(srcgloind(i))=mype
         ilsrclocind(srcgloind(i))=i
      ENDDO

      DO i=1,ndstlocpoints
         ildstmppmap(dstgloind(i))=mype
         ildstlocind(dstgloind(i))=i
      ENDDO
      
#if defined key_mpp_mpi
      CALL mpi_allreduce(ilsrcmppmap,igsrcmppmap,nsrcglopoints, &
         &               mpi_integer,mpi_max,mpi_comm,istatus)
      CALL mpi_allreduce(ilsrclocind,igsrclocind,nsrcglopoints, &
         &               mpi_integer,mpi_max,mpi_comm,istatus)
      CALL mpi_allreduce(ildstmppmap,igdstmppmap,ndstglopoints, &
         &               mpi_integer,mpi_max,mpi_comm,istatus)
      CALL mpi_allreduce(ildstlocind,igdstlocind,ndstglopoints, &
         &               mpi_integer,mpi_max,mpi_comm,istatus)
#else
      igsrcmppmap(:)=ilsrcmppmap(:)
      igsrclocind(:)=ilsrclocind(:)
      igdstmppmap(:)=ildstmppmap(:)
      igdstlocind(:)=ildstlocind(:)
#endif

      ! Optionally construct an ascii file listing what src and
      ! dest points belongs to which task

      ! Since igsrcmppmap and igdstmppmap are global data only do
      ! this for mype==0.

      IF (lcommout.AND.(mype==0)) THEN
         WRITE(cdfile,'(A,I4.4,A)')commoutprefix//'_srcmppmap_',mype+1,'.dat'
         OPEN(9,file=cdfile)
         DO i=1,nsrcglopoints
            WRITE(9,*)remap%src%grid_center_lat(i),&
               & remap%src%grid_center_lon(i), &
               & igsrcmppmap(i)+1,remap%src%grid_imask(i)
         ENDDO
         CLOSE(9)
         WRITE(cdfile,'(A,I4.4,A)')commoutprefix//'_dstmppmap_',mype+1,'.dat'
         OPEN(9,file=cdfile)
         DO i=1,ndstglopoints
            WRITE(9,*)remap%dst%grid_center_lat(i),&
               & remap%dst%grid_center_lon(i), &
               & igdstmppmap(i)+1,remap%dst%grid_imask(i)
         ENDDO
         CLOSE(9)
      ENDIF

      !
      ! Standard interpolation in serial case is
      !
      ! DO n=1,remap%num_links
      !   zdst(remap%dst_address(n)) = zdst(remap%dst_address(n)) + &
      !      & remap%remap_matrix(1,n)*zsrc(remap%src_address(n))
      ! END DO
      !

      ! In parallel we need to first find local number of links
      
      pinfo%num_links=0
      DO i=1,remap%num_links
         IF (igdstmppmap(remap%dst_address(i))==mype) &
            & pinfo%num_links=pinfo%num_links+1
      ENDDO
      ALLOCATE(pinfo%dst_address(pinfo%num_links),&
         &     pinfo%src_address(pinfo%num_links),&
         &     pinfo%remap_matrix(1,pinfo%num_links))

      ! Get local destination addresses

      n=0
      DO i=1,remap%num_links
         IF (igdstmppmap(remap%dst_address(i))==mype) THEN
            n=n+1
            pinfo%dst_address(n)=&
               & igdstlocind(remap%dst_address(i))
            pinfo%remap_matrix(:,n)=&
               & remap%remap_matrix(:,i)
         ENDIF
      ENDDO

      ! Get sending processors maps.

      ! The same data point might need to be sent to many processors
      ! so first construct a map for processors needing the data

      isrcpe(:)=-1
      DO i=1,remap%num_links
         IF (igdstmppmap(remap%dst_address(i))==mype) THEN
            isrcpe(remap%src_address(i))=&
               & igsrcmppmap(remap%src_address(i))
         ENDIF
      ENDDO

      ! Optionally write a set if ascii file listing which tasks
      ! mype needs to send to communicate with

      IF (lcommout) THEN
         ! Destination processors
         WRITE(cdfile,'(A,I4.4,A)')commoutprefix//'_dsts_',mype+1,'.dat'
         OPEN(9,file=cdfile)
         DO pe=0,nproc-1
            IF (pe==mype) THEN
               isrcpetmp(:)=isrcpe(:)
            ENDIF
#if defined key_mpp_mpi
            CALL mpi_bcast(isrcpetmp,nsrcglopoints,mpi_integer,pe,mpi_comm,istatus)
#endif
            DO i=1,nsrcglopoints
               IF (isrcpetmp(i)==mype) THEN
                  WRITE(9,*)remap%src%grid_center_lat(i),&
                     & remap%src%grid_center_lon(i), &
                     & pe+1,mype+1
               ENDIF
            ENDDO
         ENDDO
         CLOSE(9)
      ENDIF

      ! Get number of points to send to each processor

      ALLOCATE(pinfo%nsend(0:nproc-1))
      isrcno(:)=0
      DO i=1,nsrcglopoints
         IF (isrcpe(i)>=0) THEN
            isrcno(isrcpe(i)+1)=isrcno(isrcpe(i)+1)+1
         ENDIF
      ENDDO
#if defined key_mpp_mpi
      CALL mpi_alltoall(isrcno,1,mpi_integer, &
         &              pinfo%nsend(0:nproc-1),1,mpi_integer, &
         &              mpi_comm,istatus)
#else
      pinfo%nsend(0:nproc-1) = isrcno(1:nproc)
#endif
      pinfo%nsendtot=SUM(pinfo%nsend(0:nproc-1))

      ! Construct sending buffer mapping. Data is mapping in
      ! processor order.

      ALLOCATE(pinfo%send_address(pinfo%nsendtot))
      
      ! Temporary arrays for mpi all to all.

      ALLOCATE(isrctot(SUM(isrcno(1:nproc))))
      ALLOCATE(ircvtot(SUM(pinfo%nsend(0:nproc-1))))

      ! Offset for message parsing

      isrcoff(1)=0
      ircvoff(1)=0
      DO i=1,nproc-1
         isrcoff(i+1) = isrcoff(i) + isrcno(i)
         ircvoff(i+1) = pinfo%nsend(i-1) + ircvoff(i)
      ENDDO

      ! Pack indices i into a buffer

      isrccur(:)=0
      DO i=1,nsrcglopoints
         IF (isrcpe(i)>=0) THEN
            isrccur(isrcpe(i)+1)=isrccur(isrcpe(i)+1)+1
            isrctot(isrccur(isrcpe(i)+1)+isrcoff(isrcpe(i)+1)) = i
         ENDIF
      ENDDO
      
      ! Send the data

#if defined key_mpp_mpi
      CALL mpi_alltoallv(&
         & isrctot,isrccur,isrcoff,mpi_integer, &
         & ircvtot,pinfo%nsend(0:nproc-1),ircvoff,mpi_integer, &
         & mpi_comm,istatus)
#else
      ircvtot(:)=isrctot(:)
#endif

      ! Get the send address. ircvtot will at this point contain the 
      ! addresses in the global index needed for message passing

      DO i=1,pinfo%nsendtot
         pinfo%send_address(i)=igsrclocind(ircvtot(i))
      ENDDO

      ! Deallocate the mpi all to all arrays

      DEALLOCATE(ircvtot,isrctot)

      ! Get number of points to receive to each processor

      ALLOCATE(pinfo%nrecv(0:nproc-1))
      pinfo%nrecv(0:nproc-1)=0
      DO i=1,nsrcglopoints
        IF (isrcpe(i)>=0 .AND. isrcpe(i)<nproc) THEN
          pinfo%nrecv(isrcpe(i))=pinfo%nrecv(isrcpe(i))+1
        ENDIF
      ENDDO
      pinfo%nrecvtot=SUM(pinfo%nrecv(0:nproc-1))

      ! Find new src address mapping

      ! Setup local positions in the global array in a temporary array
      ! taking into accound the processor ordering

      isrcaddtmp(:)=-1
      isrcoffset(0)=0
      DO pe=1,nproc-1
        isrcoffset(pe)=isrcoffset(pe-1)+pinfo%nrecv(pe-1)
      ENDDO

      DO i=1,nsrcglopoints
        IF (isrcpe(i)>=0 .AND. isrcpe(i)<nproc) THEN
          isrcoffset(isrcpe(i))=isrcoffset(isrcpe(i))+1
          isrcaddtmp(i)=isrcoffset(isrcpe(i))
        ENDIF
      ENDDO 

      ! Find the local positions in the temporary array.

      n=0
      DO i=1,remap%num_links
         IF (igdstmppmap(remap%dst_address(i))==mype) THEN
            n=n+1
            pinfo%src_address(n)=isrcaddtmp(remap%src_address(i))
         ENDIF
      ENDDO

      ! MPI displacements for mpi_alltoallv

      ALLOCATE(pinfo%nsdisp(0:nproc-1),pinfo%nrdisp(0:nproc-1))
      pinfo%nsdisp(0)=0
      pinfo%nrdisp(0)=0
      DO pe=1,nproc-1
         pinfo%nsdisp(pe)=pinfo%nsdisp(pe-1)+pinfo%nsend(pe-1)
         pinfo%nrdisp(pe)=pinfo%nrdisp(pe-1)+pinfo%nrecv(pe-1)
      ENDDO

      ! Optionally construct an ascii file listing number of
      ! data to send and receive from which processor.

      IF (lcommout) THEN
         WRITE(cdfile,'(A,I4.4,A)')commoutprefix//'_nsend_',mype+1,'.dat'
         OPEN(9,file=cdfile)
         DO pe=0,nproc-1
            WRITE(9,*)pe+1,pinfo%nsend(pe)
         ENDDO
         CLOSE(9)
         WRITE(cdfile,'(A,I4.4,A)')commoutprefix//'_nrecv_',mype+1,'.dat'
         OPEN(9,file=cdfile)
         DO pe=0,nproc-1
            WRITE(9,*)pe+1,pinfo%nrecv(pe)
         ENDDO
         CLOSE(9)
      ENDIF

   END SUBROUTINE parinter_init

   SUBROUTINE parinter_fld( mype, nproc, mpi_comm, &
      & pinfo, nsrclocpoints, zsrc, ndstlocpoints, zdst )

      ! Perform a single interpolation from the zsrc field
      ! to zdst field based on the information in pinfk
      
      ! Input arguments

      ! Message passing information
      INTEGER, INTENT(IN) :: mype, nproc, mpi_comm
      ! Interpolation setup
      TYPE(parinterinfo), INTENT(IN) :: pinfo
      ! Source data/
      INTEGER, INTENT(IN) :: nsrclocpoints
      REAL, INTENT(IN), DIMENSION(nsrclocpoints) :: zsrc

      ! Output arguments

      ! Destination data
      INTEGER , INTENT(IN):: ndstlocpoints
      REAL, DIMENSION(ndstlocpoints) :: zdst

      ! Local variables

      ! MPI send/recv buffers
#if defined key_parinter_alloc
      REAL(scripdp) , ALLOCATABLE :: zsend(:),zrecv(:)
#else
      REAL(scripdp) :: zsend(pinfo%nsendtot),zrecv(pinfo%nrecvtot)
#endif
      ! Misc variables
      INTEGER :: i,iproc,istatus,ierr,off,itag,irq,nreqs
      INTEGER :: reqs(0:2*(nproc-1))

#if defined key_parinter_alloc
      ALLOCATE(zsend(pinfo%nsendtot),zrecv(pinfo%nrecvtot))
#endif

      ! Pack the sending buffer

      DO i=1,pinfo%nsendtot
         zsend(i)=zsrc(pinfo%send_address(i))
      ENDDO
      
      ! Do the message passing

#if defined key_mpp_mpi

      IF (lparinterp2p) THEN

         ! total num of reqs ( recv + send )
         nreqs = 2*(nproc-1)
         ! post irecvs first
         irq = 0
         DO iproc=0,nproc-1
            IF( pinfo%nrecv(iproc) > 0 .AND. iproc /= mype ) THEN
               off = pinfo%nrdisp(iproc)
               itag=100
               CALL mpi_irecv(zrecv(off+1),pinfo%nrecv(iproc),mpi_double_precision, &
                  & iproc,itag,mpi_comm,reqs(irq),ierr) 
               irq = irq + 1
               IF( irq > nreqs ) THEN
                  WRITE(0,*)'parinter_fld_mult: exceeded number of reqs when posting recvs',mype
                  CALL abort
               ENDIF
            ENDIF
         ENDDO
         ! post isends
         DO iproc=0,nproc-1
            IF( pinfo%nsend(iproc) > 0 )  THEN
               IF( iproc == mype ) THEN
                  zrecv( pinfo%nrdisp(iproc)+1:pinfo%nrdisp(iproc)+pinfo%nsend(iproc) ) = &
                     & zsend( pinfo%nsdisp(iproc)+1:pinfo%nsdisp(iproc)+pinfo%nsend(iproc) )
               ELSE
                  off = pinfo%nsdisp(iproc)
                  itag=100
                  CALL mpi_isend(zsend(off+1),pinfo%nsend(iproc),mpi_double_precision, &
                     & iproc,itag,mpi_comm,reqs(irq),ierr)
                  irq = irq + 1 
                  IF( irq > nreqs ) THEN
                     WRITE(0,*)'parinter_fld_mult: exceeded number of reqs when posting sends',mype
                     CALL abort
                  ENDIF
               ENDIF
            ENDIF
         ENDDO
         ! wait on requests
         CALL mpi_waitall(irq,reqs,MPI_STATUSES_IGNORE,ierr)

      ELSE

         CALL mpi_alltoallv(&
            & zsend,pinfo%nsend(0:nproc-1),&
            & pinfo%nsdisp(0:nproc-1),mpi_double_precision, &
            & zrecv,pinfo%nrecv(0:nproc-1), &
            & pinfo%nrdisp(0:nproc-1),mpi_double_precision, &
            & mpi_comm,istatus)

      ENDIF
#else

      zrecv(:)=zsend(:)

#endif

      ! Do the interpolation

      zdst(:)=0.0
      DO i=1,pinfo%num_links
         zdst(pinfo%dst_address(i)) = zdst(pinfo%dst_address(i)) + &
            & pinfo%remap_matrix(1,i)*zrecv(pinfo%src_address(i))
      END DO

#if defined key_parinter_alloc
      DEALLOCATE(zsend,zrecv)
#endif

   END SUBROUTINE parinter_fld

   SUBROUTINE parinter_fld_mult( nfield, &
      & mype, nproc, mpi_comm, &
      & pinfo, nsrclocpoints, zsrc, ndstlocpoints, zdst )

      ! Perform nfield interpolations from the zsrc fields
      ! to zdst fields based on the information in pinfo

      ! Input arguments

      ! Message passing information
      INTEGER, INTENT(IN) :: mype, nproc, mpi_comm, nfield
      ! Interpolation setup
      TYPE(parinterinfo), INTENT(IN) :: pinfo
      ! Source data/
      INTEGER, INTENT(IN) :: nsrclocpoints
      REAL, INTENT(IN), DIMENSION(nsrclocpoints,nfield) :: zsrc

      ! Output arguments

      ! Destination data
      INTEGER, INTENT(IN):: ndstlocpoints
      REAL, DIMENSION(ndstlocpoints,nfield) :: zdst

      INTEGER :: nsend(0:nproc-1), nrecv(0:nproc-1),nrdisp(0:nproc-1), nsdisp(0:nproc-1)

      ! Local variables

      ! MPI send/recv buffers
#if defined key_parinter_alloc
      REAL(scripdp), ALLOCATABLE, DIMENSION(:,:) :: zrecvnf
      REAL(scripdp), ALLOCATABLE, DIMENSION(:) :: zsend, zrecv
#else
      REAL(scripdp) :: zrecvnf(pinfo%nrecvtot,nfield)
      REAL(scripdp) :: zsend(pinfo%nsendtot*nfield), &
         & zrecv(pinfo%nrecvtot*nfield)
#endif
      ! Misc variables
      INTEGER :: i,istatus,ierr
      INTEGER :: nf, ibases, ibaser, np, iproc, off, itag, irq, nreqs
      INTEGER :: reqs(0:2*(nproc-1))
      INTEGER, DIMENSION(nfield,0:nproc-1) :: ibaseps, ibasepr

      ! Allocate temporary arrays on heap
      
#if defined key_parinter_alloc
      ALLOCATE(zrecvnf(pinfo%nrecvtot,nfield),&
         & zsend(pinfo%nsendtot*nfield),zrecv(pinfo%nrecvtot*nfield))
#endif

      ! Compute starts for packing

      ibases=0
      ibaser=0
      DO np=0,nproc-1
         DO nf=1,nfield
            ibaseps(nf,np) = ibases
            ibasepr(nf,np) = ibaser
            ibases = ibases + pinfo%nsend(np)
            ibaser = ibaser + pinfo%nrecv(np)
         ENDDO
      ENDDO

      ! Pack the sending buffer

      !$omp parallel default(shared) private(nf,np,i)
      !$omp do schedule(dynamic)
      DO np=0,nproc-1
         DO nf=1,nfield
            DO i=1,pinfo%nsend(np)
               zsend(i+ibaseps(nf,np))=&
                  & zsrc(pinfo%send_address(i+pinfo%nsdisp(np)),nf)
            ENDDO
         ENDDO
         nsend(np)=pinfo%nsend(np)*nfield
         nrecv(np)=pinfo%nrecv(np)*nfield
         nrdisp(np)=pinfo%nrdisp(np)*nfield
         nsdisp(np)=pinfo%nsdisp(np)*nfield
      ENDDO
      !$omp end do
      !$omp end parallel

      ! Do the message passing

#if defined key_mpp_mpi

      IF (lparinterp2p) THEN

         ! total num of reqs ( recv + send )
         nreqs = 2*(nproc-1)
         ! post irecvs first
         irq = 0
         DO iproc=0,nproc-1
            IF( nrecv(iproc) > 0 .AND. iproc /= mype ) THEN
               off = nrdisp(iproc)
               itag=100
               CALL mpi_irecv(zrecv(off+1),nrecv(iproc),mpi_double_precision, &
                  & iproc,itag,mpi_comm,reqs(irq),ierr) 
               irq = irq + 1
               IF( irq > nreqs ) THEN
                  WRITE(0,*)'parinter_fld_mult: exceeded number of reqs when posting recvs',mype
                  CALL abort
               ENDIF
            ENDIF
         ENDDO
         ! post isends
         DO iproc=0,nproc-1
            IF( nsend(iproc) > 0 )  THEN
               IF( iproc == mype ) THEN
                  zrecv( nrdisp(iproc)+1:nrdisp(iproc)+nsend(iproc) ) = &
                     & zsend( nsdisp(iproc)+1:nsdisp(iproc)+nsend(iproc) )
               ELSE
                  off = nsdisp(iproc)
                  itag=100
                  CALL mpi_isend(zsend(off+1),nsend(iproc),mpi_double_precision, &
                     & iproc,itag,mpi_comm,reqs(irq),ierr)
                  irq = irq + 1 
                  IF( irq > nreqs ) THEN
                     WRITE(0,*)'parinter_fld_mult: exceeded number of reqs when posting sends',mype
                     CALL abort
                  ENDIF
               ENDIF
            ENDIF
         ENDDO
         ! wait on requests
         CALL mpi_waitall(irq,reqs,MPI_STATUSES_IGNORE,ierr)

      ELSE

         IF(mype==0)WRITE(0,*)'lparinterp2p off'
         CALL mpi_alltoallv(&
            & zsend,nsend(0:nproc-1),&
            & nsdisp(0:nproc-1),mpi_double_precision, &
            & zrecv,nrecv(0:nproc-1), &
            & nrdisp(0:nproc-1),mpi_double_precision, &
            & mpi_comm,istatus)

      ENDIF
#else

      zrecv(:)=zsend(:)

#endif
 
      ! Unpack individual fields

      !$omp parallel default(shared) private(nf,np,i)

      !$omp do schedule (dynamic)
      DO np=0,nproc-1
         DO nf=1,nfield
            DO i=1,pinfo%nrecv(np)
               zrecvnf(i+pinfo%nrdisp(np),nf)=zrecv(i+ibasepr(nf,np))
            ENDDO
         ENDDO
      ENDDO
      !omp end do

      ! Do the interpolation

      !$omp do
      DO nf=1,nfield
         zdst(:,nf)=0.0
         DO i=1,pinfo%num_links
            zdst(pinfo%dst_address(i),nf) = zdst(pinfo%dst_address(i),nf) + &
               & pinfo%remap_matrix(1,i)*zrecvnf(pinfo%src_address(i),nf)
         END DO
      END DO
      !$omp end do

      !$omp end parallel

#if defined key_parinter_alloc
      DEALLOCATE( zrecvnf, zsend, zrecv )
#endif

   END SUBROUTINE parinter_fld_mult

   SUBROUTINE parinter_write( mype, nproc, &
      & nsrcglopoints, ndstglopoints, &
      & pinfo, cdpath, cdprefix )

      ! Write pinfo information in a netCDF file in order to
      ! be able to read it rather than calling parinter_init

      ! Input arguments.

      ! Message passing information
      INTEGER, INTENT(IN) :: mype, nproc
      ! Source grid local global number of grid points
      INTEGER, INTENT(IN) :: nsrcglopoints
      ! Destination grid global number of grid points
      INTEGER, INTENT(IN) :: ndstglopoints
      ! Interpolation and message passing information
      TYPE(parinterinfo), INTENT(IN) :: pinfo
      ! Path and file prefix
      CHARACTER(len=*) :: cdpath, cdprefix

      ! Local variable

      ! Misc variable
      CHARACTER(len=1024) :: cdfile
      INTEGER :: ncid, dimnl, dimnw, dimnst, dimnrt, dimnpr
      INTEGER :: dims1(1), dims2(2)
      INTEGER :: idda, idsa, idrm, idns, idsaa, idnr, idnrp, idnsp

      WRITE(cdfile, '(A,I0,A,2(I8.8,A),2(I4.4,A),A)') &
         & TRIM(cdpath) // '/dist_', nproc, '/'// TRIM(cdprefix)//'_', &
         & nsrcglopoints, '_', ndstglopoints, '_', mype, '_', nproc, '.nc'

      call execute_command_line('mkdir -p ' // cdpath)
      CALL nchdlerr(nf90_create(TRIM(cdfile),nf90_clobber,ncid), &
         &          __LINE__, __MYFILE__ )

      ! To avoid problems with multiple unlimited netCDF dimensions
      ! we don't write num_links, nsendtot and nrecvtot if
      ! they are 0. If they are not there we assume they are 0
      ! in parinter_read.

      IF (pinfo%num_links>0) THEN
         CALL nchdlerr(nf90_def_dim(ncid,'num_links',&
            &                       pinfo%num_links,dimnl),&
            &          __LINE__,__MYFILE__)
      ENDIF

      CALL nchdlerr(nf90_def_dim(ncid,'num_wgts',&
         &                       1,dimnw),&
         &          __LINE__,__MYFILE__)

      IF (pinfo%nsendtot>0) THEN
         CALL nchdlerr(nf90_def_dim(ncid,'nsendtot',&
            &                       pinfo%nsendtot,dimnst),&
            &          __LINE__,__MYFILE__)
      ENDIF
      
      IF (pinfo%nrecvtot>0) THEN
         CALL nchdlerr(nf90_def_dim(ncid,'nrecvtot',&
            &                       pinfo%nrecvtot,dimnrt),&
            &          __LINE__,__MYFILE__)
      ENDIF
      
      CALL nchdlerr(nf90_def_dim(ncid,'nproc',&
         &                       nproc,dimnpr),&
         &          __LINE__,__MYFILE__)

      IF (pinfo%num_links>0) THEN

         dims1(1)=dimnl
         CALL nchdlerr(nf90_def_var(ncid,'dst_address',&
            &                       nf90_int,dims1,idda),&
            &          __LINE__,__MYFILE__)
         
         dims1(1)=dimnl
         CALL nchdlerr(nf90_def_var(ncid,'src_address',&
            &                       nf90_int,dims1,idsa),&
            &          __LINE__,__MYFILE__)
         
         dims2(1)=dimnw
         dims2(2)=dimnl
         CALL nchdlerr(nf90_def_var(ncid,'remap_matrix',&
            &                       nf90_double,dims2,idrm),&
            &          __LINE__,__MYFILE__)

      ENDIF

      dims1(1)=dimnpr
      CALL nchdlerr(nf90_def_var(ncid,'nsend',&
         &                       nf90_int,dims1,idns),&
         &          __LINE__,__MYFILE__)

      IF (pinfo%nsendtot>0) THEN

         dims1(1)=dimnst
         CALL nchdlerr(nf90_def_var(ncid,'send_address',&
            &                       nf90_int,dims1,idsaa),&
            &          __LINE__,__MYFILE__)

      ENDIF

      dims1(1)=dimnpr
      CALL nchdlerr(nf90_def_var(ncid,'nrecv',&
         &                       nf90_int,dims1,idnr),&
         &          __LINE__,__MYFILE__)

      dims1(1)=dimnpr
      CALL nchdlerr(nf90_def_var(ncid,'nsdisp',&
         &                       nf90_int,dims1,idnsp),&
         &          __LINE__,__MYFILE__)

      dims1(1)=dimnpr
      CALL nchdlerr(nf90_def_var(ncid,'nrdisp',&
         &                       nf90_int,dims1,idnrp),&
         &          __LINE__,__MYFILE__)

      CALL nchdlerr(nf90_enddef(ncid),__LINE__,__MYFILE__)


      IF (pinfo%num_links>0) THEN

         CALL nchdlerr(nf90_put_var(ncid,idda,pinfo%dst_address),&
            &          __LINE__,__MYFILE__)

         CALL nchdlerr(nf90_put_var(ncid,idsa,pinfo%src_address),&
            &          __LINE__,__MYFILE__)

         CALL nchdlerr(nf90_put_var(ncid,idrm,pinfo%remap_matrix),&
            &          __LINE__,__MYFILE__)

      ENDIF

      CALL nchdlerr(nf90_put_var(ncid,idns,pinfo%nsend(0:nproc-1)),&
         &          __LINE__,__MYFILE__)

      IF (pinfo%nsendtot>0) THEN

         CALL nchdlerr(nf90_put_var(ncid,idsaa,pinfo%send_address),&
            &          __LINE__,__MYFILE__)

      ENDIF

      CALL nchdlerr(nf90_put_var(ncid,idnr,pinfo%nrecv(0:nproc-1)),&
         &          __LINE__,__MYFILE__)

      CALL nchdlerr(nf90_put_var(ncid,idnsp,pinfo%nsdisp(0:nproc-1)),&
         &          __LINE__,__MYFILE__)

      CALL nchdlerr(nf90_put_var(ncid,idnrp,pinfo%nrdisp(0:nproc-1)),&
         &          __LINE__,__MYFILE__)

      CALL nchdlerr(nf90_close(ncid),__LINE__, __MYFILE__ )

   END SUBROUTINE parinter_write

   SUBROUTINE parinter_read( mype, nproc, &
      & nsrcglopoints, ndstglopoints, &
      & pinfo, cdpath, cdprefix, lexists )

      ! Write pinfo information in a netCDF file in order to
      ! be able to read it rather than calling parinter_init

      ! Input arguments.

      ! Message passing information
      INTEGER, INTENT(IN) :: mype, nproc
      ! Source grid local global number of grid points
      INTEGER, INTENT(IN) :: nsrcglopoints
      ! Destination grid global number of grid points
      INTEGER, INTENT(IN) :: ndstglopoints
      ! Interpolation and message passing information
      TYPE(parinterinfo), INTENT(OUT) :: pinfo
      ! Does the information exists
      LOGICAL :: lexists
      ! Path and file prefix
      CHARACTER(len=*) :: cdpath, cdprefix

      ! Local variable

      ! Misc variable
      CHARACTER(len=1024) :: cdfile
      INTEGER :: ncid, dimid, varid, num_wgts

      WRITE(cdfile, '(A,I0,A,2(I8.8,A),2(I4.4,A),A)') &
         & TRIM(cdpath) // '/dist_', nproc, '/'// TRIM(cdprefix)//'_', &
         & nsrcglopoints, '_', ndstglopoints, '_', mype, '_', nproc, '.nc'


      lexists=nf90_open(TRIM(cdfile),nf90_nowrite,ncid)==nf90_noerr

      IF (lexists) THEN
        IF (mype==0) THEN
          WRITE(*,*) "FESOM ifs-interface: precomputed weight file exists, reading it in."
        ENDIF

         ! If num_links is not present we assume it to be zero.
         
         IF (nf90_inq_dimid(ncid,'num_links',dimid)==nf90_noerr) THEN
            CALL nchdlerr(nf90_inquire_dimension(ncid,dimid,&
               &                                 len=pinfo%num_links),&
               &          __LINE__,__MYFILE__)
         ELSE
            pinfo%num_links=0
         ENDIF

         CALL nchdlerr(nf90_inq_dimid(ncid,'num_wgts',dimid),&
            &          __LINE__,__MYFILE__)
         CALL nchdlerr(nf90_inquire_dimension(ncid,dimid,&
            &                                 len=num_wgts),&
            &          __LINE__,__MYFILE__)
         IF (num_wgts/=1) THEN
            WRITE(0,*)'parinter_read: num_wgts has to be 1 for now'
            CALL abort
         ENDIF

         ! If nsendtot is not present we assume it to be zero.

         IF (nf90_inq_dimid(ncid,'nsendtot',dimid)==nf90_noerr) THEN
            CALL nchdlerr(nf90_inquire_dimension(ncid,dimid,&
               &                                 len=pinfo%nsendtot),&
               &          __LINE__,__MYFILE__)
         ELSE
            pinfo%nsendtot=0
         ENDIF

         IF(nf90_inq_dimid(ncid,'nrecvtot',dimid)==nf90_noerr) THEN
            CALL nchdlerr(nf90_inquire_dimension(ncid,dimid,&
               &                                 len=pinfo%nrecvtot),&
               &          __LINE__,__MYFILE__)
         ELSE
            pinfo%nrecvtot=0
         ENDIF

         ALLOCATE(pinfo%dst_address(pinfo%num_links),&
            &     pinfo%src_address(pinfo%num_links),&
            &     pinfo%remap_matrix(num_wgts,pinfo%num_links),&
            &     pinfo%nsend(0:nproc-1),&
            &     pinfo%send_address(pinfo%nsendtot),&
            &     pinfo%nrecv(0:nproc-1),&
            &     pinfo%nsdisp(0:nproc-1),&
            &     pinfo%nrdisp(0:nproc-1))
 
         IF (pinfo%num_links>0) THEN
            CALL nchdlerr(nf90_inq_varid(ncid,'dst_address',varid),&
               &          __LINE__,__MYFILE__)
            CALL nchdlerr(nf90_get_var(ncid,varid,pinfo%dst_address),&
               &          __LINE__,__MYFILE__)
            
            CALL nchdlerr(nf90_inq_varid(ncid,'src_address',varid),&
               &          __LINE__,__MYFILE__)
            CALL nchdlerr(nf90_get_var(ncid,varid,pinfo%src_address),&
               &          __LINE__,__MYFILE__)
            
            CALL nchdlerr(nf90_inq_varid(ncid,'remap_matrix',varid),&
               &          __LINE__,__MYFILE__)
            CALL nchdlerr(nf90_get_var(ncid,varid,pinfo%remap_matrix),&
               &          __LINE__,__MYFILE__)
         ENDIF

         CALL nchdlerr(nf90_inq_varid(ncid,'nsend',varid),&
            &          __LINE__,__MYFILE__)
         CALL nchdlerr(nf90_get_var(ncid,varid,pinfo%nsend(0:nproc-1)),&
            &          __LINE__,__MYFILE__)

         IF (pinfo%nsendtot>0) THEN

            CALL nchdlerr(nf90_inq_varid(ncid,'send_address',varid),&
               &          __LINE__,__MYFILE__)
            CALL nchdlerr(nf90_get_var(ncid,varid,pinfo%send_address),&
               &          __LINE__,__MYFILE__)

         ENDIF

         CALL nchdlerr(nf90_inq_varid(ncid,'nrecv',varid),&
            &          __LINE__,__MYFILE__)
         CALL nchdlerr(nf90_get_var(ncid,varid,pinfo%nrecv(0:nproc-1)),&
            &          __LINE__,__MYFILE__)

         CALL nchdlerr(nf90_inq_varid(ncid,'nsdisp',varid),&
            &          __LINE__,__MYFILE__)
         CALL nchdlerr(nf90_get_var(ncid,varid,pinfo%nsdisp(0:nproc-1)),&
            &          __LINE__,__MYFILE__)

         CALL nchdlerr(nf90_inq_varid(ncid,'nrdisp',varid),&
            &          __LINE__,__MYFILE__)
         CALL nchdlerr(nf90_get_var(ncid,varid,pinfo%nrdisp(0:nproc-1)),&
            &          __LINE__,__MYFILE__)

         CALL nchdlerr(nf90_close(ncid),__LINE__, __MYFILE__ )

      ENDIF

   END SUBROUTINE parinter_read
   
END MODULE parinter

MODULE interinfo

   ! Parallel regridding information

   USE parinter

   IMPLICIT NONE

   SAVE

   ! IFS to NEMO

   TYPE(parinterinfo) :: gausstoT,gausstoUV

   ! NEMO to IFS

   TYPE(parinterinfo) :: Ttogauss, UVtogauss

   ! Read parinterinfo on task 0 only and broadcast.

   LOGICAL :: lparbcast = .FALSE.

   ! Use multiple fields option

   LOGICAL :: lparintmultatm = .FALSE.
   
END MODULE interinfo
