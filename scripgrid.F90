#define __MYFILE__ 'scripgrid.F90'
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
