#define __MYFILE__ 'scripremap.F90'
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
