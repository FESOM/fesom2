SUBROUTINE nemogcmcoup_coupinit( mype, npes, icomm, &
   &                             npoints, nlocmsk, ngloind )

   ! Initialize single executable coupling 
   USE parinter
   USE scripremap
   USE interinfo
   IMPLICIT NONE

   ! Input arguments

   ! Message passing information
   INTEGER, INTENT(IN) :: mype,npes,icomm
   ! Gaussian grid information   
   ! Number of points
   INTEGER, INTENT(IN) :: npoints
   ! Integer mask and global indices
   INTEGER, DIMENSION(npoints), INTENT(IN) :: nlocmsk, ngloind
   INTEGER :: iunit = 0

   ! Local variables

   ! Namelist containing the file names of the weights
   CHARACTER(len=256) :: cdfile_gauss_to_T, cdfile_gauss_to_UV, &
      &                  cdfile_T_to_gauss, cdfile_UV_to_gauss
   CHARACTER(len=256) :: cdpathdist
   LOGICAL :: lwritedist, lreaddist
   LOGICAL :: lcommout
   CHARACTER(len=128) :: commoutprefix
   NAMELIST/namnemocoup/cdfile_gauss_to_T,&
      &                 cdfile_gauss_to_UV,&
      &                 cdfile_T_to_gauss,&
      &                 cdfile_UV_to_gauss,&
      &                 cdpathdist, &
      &                 lreaddist, &
      &                 lwritedist, &
      &                 lcommout, &
      &                 commoutprefix,&
      &                 lparbcast

   ! Global number of gaussian gridpoints
   INTEGER :: nglopoints
   ! Ocean grids accessed with NEMO modules
   INTEGER :: noglopoints,nopoints
   INTEGER, ALLOCATABLE, DIMENSION(:) :: omask,ogloind
   ! SCRIP remapping data structures.
   TYPE(scripremaptype) :: remap_gauss_to_T, remap_T_to_gauss, &
      & remap_gauss_to_UV, remap_UV_to_gauss
   ! Misc variables
   INTEGER :: i,j,k,ierr
   LOGICAL :: lexists

   ! Read namelists
   
   cdfile_gauss_to_T = 'gausstoT.nc'
   cdfile_gauss_to_UV = 'gausstoUV.nc'
   cdfile_T_to_gauss = 'Ttogauss.nc'
   cdfile_UV_to_gauss = 'UVtogauss.nc'
   lcommout          = .FALSE.
   commoutprefix     = 'parinter_comm'
   cdpathdist        = './'
   lreaddist         = .FALSE.
   lwritedist        = .FALSE.

   OPEN(9,file='namnemocoup.in')
   READ(9,namnemocoup)
   CLOSE(9)

   ! Global number of Gaussian gridpoints

#if defined key_mpp_mpi
   CALL mpi_allreduce( npoints, nglopoints, 1, &
      &                mpi_integer, mpi_sum, icomm, ierr)
#else
   nglopoints=npoints
#endif

   WRITE(0,*)'Update FESOM global scalar points'
   noglopoints=126858
   IF (mype==0) THEN
      nopoints=126858
   ELSE
      nopoints=0
   ENDIF

   ! Ocean mask and global indicies
   
   ALLOCATE(omask(MAX(nopoints,1)),ogloind(MAX(nopoints,1)))

   omask(:) = 1
   IF (mype==0) THEN
      DO i=1,nopoints
         ogloind(i)=i
      ENDDO
   ENDIF

   ! Read the interpolation weights and setup the parallel interpolation
   ! from atmosphere Gaussian grid to ocean T-grid
   
   IF (lreaddist) THEN
      CALL parinter_read( mype, npes, nglopoints, noglopoints, gausstoT,  &
         & cdpathdist,'ifs_to_fesom_gridT',lexists)
   ENDIF
   IF ((.NOT.lreaddist).OR.(.NOT.lexists)) THEN
      IF (lparbcast) THEN
         CALL scripremap_read_sgl(cdfile_gauss_to_T,remap_gauss_to_T,&
            & mype,npes,icomm,.TRUE.)
      ELSE
         CALL scripremap_read(cdfile_gauss_to_T,remap_gauss_to_T)
      ENDIF
      CALL parinter_init( mype, npes, icomm, &
         & npoints, nglopoints, nlocmsk, ngloind, &
         & nopoints, noglopoints, omask, ogloind, & 
         & remap_gauss_to_T, gausstoT, lcommout, TRIM(commoutprefix)//'_gtoT', &
         & iunit )
      CALL scripremap_dealloc(remap_gauss_to_T)
      IF (lwritedist) THEN
         CALL parinter_write( mype, npes, nglopoints, noglopoints, gausstoT,  &
            & cdpathdist,'ifs_to_fesom_gridT')
      ENDIF
   ENDIF

   ! From ocean T-grid to atmosphere Gaussian grid

   IF (lreaddist) THEN
      CALL parinter_read( mype, npes, noglopoints, nglopoints, Ttogauss,  &
         & cdpathdist,'fesom_gridT_to_ifs',lexists)
   ENDIF
   IF ((.NOT.lreaddist).OR.(.NOT.lexists)) THEN
      IF (lparbcast) THEN
         CALL scripremap_read_sgl(cdfile_T_to_gauss,remap_T_to_gauss,&
            & mype,npes,icomm,.TRUE.)
      ELSE
         CALL scripremap_read(cdfile_T_to_gauss,remap_T_to_gauss)
      ENDIF

      CALL parinter_init( mype, npes, icomm, &
         & nopoints, noglopoints, omask, ogloind, & 
         & npoints, nglopoints, nlocmsk, ngloind, &
         & remap_T_to_gauss, Ttogauss, lcommout, TRIM(commoutprefix)//'_Ttog', &
         & iunit )
      CALL scripremap_dealloc(remap_T_to_gauss)
      IF (lwritedist) THEN
         CALL parinter_write( mype, npes, noglopoints, nglopoints, Ttogauss,  &
            & cdpathdist,'fesom_gridT_to_ifs')
      ENDIF
   ENDIF
   
   DEALLOCATE(omask,ogloind)

   WRITE(0,*)'Update FESOM global vector points'
   noglopoints=244659
   IF (mype==0) THEN
      nopoints=244659
   ELSE
      nopoints=0
   ENDIF

   ! Ocean mask and global indicies
   
   ALLOCATE(omask(MAX(nopoints,1)),ogloind(MAX(nopoints,1)))

   omask(:) = 1
   IF (mype==0) THEN
      DO i=1,nopoints
         ogloind(i)=i
      ENDDO
   ENDIF

   ! Read the interpolation weights and setup the parallel interpolation
   ! from atmosphere Gaussian grid to ocean UV-grid
   
   IF (lreaddist) THEN
      CALL parinter_read( mype, npes, nglopoints, noglopoints, gausstoUV,  &
         & cdpathdist,'ifs_to_fesom_gridUV',lexists)
   ENDIF
   IF ((.NOT.lreaddist).OR.(.NOT.lexists)) THEN
      IF (lparbcast) THEN
         CALL scripremap_read_sgl(cdfile_gauss_to_UV,remap_gauss_to_UV,&
            & mype,npes,icomm,.TRUE.)
      ELSE
         CALL scripremap_read(cdfile_gauss_to_UV,remap_gauss_to_UV)
      ENDIF
      CALL parinter_init( mype, npes, icomm, &
         & npoints, nglopoints, nlocmsk, ngloind, &
         & nopoints, noglopoints, omask, ogloind, & 
         & remap_gauss_to_UV, gausstoUV, lcommout, TRIM(commoutprefix)//'_gtoUV', &
         & iunit )
      CALL scripremap_dealloc(remap_gauss_to_UV)
      IF (lwritedist) THEN
         CALL parinter_write( mype, npes, nglopoints, noglopoints, gausstoUV,  &
            & cdpathdist,'ifs_to_fesom_gridUV')
      ENDIF
   ENDIF

   ! From ocean UV-grid to atmosphere Gaussian grid

   IF (lreaddist) THEN
      CALL parinter_read( mype, npes, noglopoints, nglopoints, UVtogauss,  &
         & cdpathdist,'fesom_gridUV_to_ifs',lexists)
   ENDIF
   IF ((.NOT.lreaddist).OR.(.NOT.lexists)) THEN
      IF (lparbcast) THEN
         CALL scripremap_read_sgl(cdfile_UV_to_gauss,remap_UV_to_gauss,&
            & mype,npes,icomm,.TRUE.)
      ELSE
         CALL scripremap_read(cdfile_UV_to_gauss,remap_UV_to_gauss)
      ENDIF

      CALL parinter_init( mype, npes, icomm, &
         & nopoints, noglopoints, omask, ogloind, & 
         & npoints, nglopoints, nlocmsk, ngloind, &
         & remap_UV_to_gauss, UVtogauss, lcommout, TRIM(commoutprefix)//'_UVtog', &
         & iunit )
      CALL scripremap_dealloc(remap_UV_to_gauss)
      IF (lwritedist) THEN
         CALL parinter_write( mype, npes, noglopoints, nglopoints, UVtogauss,  &
            & cdpathdist,'fesom_gridUV_to_ifs')
      ENDIF
   ENDIF
   
   DEALLOCATE(omask,ogloind)
         
END SUBROUTINE nemogcmcoup_coupinit
