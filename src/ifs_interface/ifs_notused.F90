! Routines usually provided by the library that are currently
! not implemented for FESOM2.
!
! -Original code by Kristian Mogensen, ECMWF.

SUBROUTINE nemogcmcoup_mlflds_get( mype, npes, icomm, &
   &                               nlev, nopoints, pgt3d, pgs3d, pgu3d, pgv3d )

   ! Interpolate sst, ice: surf T; albedo; concentration; thickness,
   ! snow thickness and currents from the ORCA grid to the Gaussian grid. 
   
   ! This routine can be called at any point in time since it does
   ! the necessary message passing in parinter_fld. 

   USE par_kind
   IMPLICIT NONE
   
   ! Arguments
   REAL(wpIFS), DIMENSION(nopoints,nlev) :: pgt3d, pgs3d, pgu3d, pgv3d
   ! Message passing information
   INTEGER, INTENT(IN) :: mype, npes, icomm
   ! Number Gaussian grid points
   INTEGER, INTENT(IN) :: nopoints,nlev

   ! Local variables

   WRITE(0,*)'nemogcmcoup_mlflds_get should not be called when coupling to fesom.'
   CALL abort

END SUBROUTINE nemogcmcoup_mlflds_get


SUBROUTINE nemogcmcoup_get( mype, npes, icomm, &
   &                        nopoints, pgsst, pgice, pgucur, pgvcur )

   ! Interpolate sst, ice and currents from the ORCA grid
   ! to the Gaussian grid. 
   
   ! This routine can be called at any point in time since it does
   ! the necessary message passing in parinter_fld. 

   USE par_kind

   IMPLICIT NONE

   
   ! Arguments

   ! Message passing information
   INTEGER, INTENT(IN) :: mype, npes, icomm
   ! Number Gaussian grid points
   INTEGER, INTENT(IN) :: nopoints
   ! Local arrays of sst, ice and currents
   REAL(wpIFS), DIMENSION(nopoints) :: pgsst, pgice, pgucur, pgvcur

   ! Local variables

   WRITE(0,*)'nemogcmcoup_get should not be called with FESOM'
   CALL abort

END SUBROUTINE nemogcmcoup_get


SUBROUTINE nemogcmcoup_get_1way( mype, npes, icomm )

   ! Interpolate sst, ice and currents from the ORCA grid
   ! to the Gaussian grid. 
   
   ! This routine can be called at any point in time since it does
   ! the necessary message passing in parinter_fld. 

   IMPLICIT NONE

   
   ! Arguments

   ! Message passing information
   INTEGER, INTENT(IN) :: mype, npes, icomm

   ! Local variables

   WRITE(0,*)'nemogcmcoup_get_1way should not be called when coupling to fesom.'
   CALL abort

END SUBROUTINE nemogcmcoup_get_1way


SUBROUTINE nemogcmcoup_mlinit( mype, npes, icomm, &
   &                           nlev, nopoints, pdep, pmask )

   ! Get information about the vertical discretization of the ocean model
   
   ! nlevs are maximum levels on input and actual number levels on output

   USE par_kind

   IMPLICIT NONE

   ! Input arguments

   ! Message passing information
   INTEGER, INTENT(IN) :: mype,npes,icomm
   ! Grid information
   INTEGER, INTENT(INOUT) :: nlev, nopoints
   REAL(wpIFS), INTENT(OUT), DIMENSION(nlev) :: pdep
   REAL(wpIFS), INTENT(OUT), DIMENSION(nopoints,nlev) :: pmask

   ! Local variables

   ! dummy argument with explicit INTENT(OUT) declaration needs an explicit value
   pdep=0.
   pmask=0.

   WRITE(0,*)'nemogcmcoup_mlinit should not be called when coupling to fesom.'
   CALL abort
   
END SUBROUTINE nemogcmcoup_mlinit


SUBROUTINE nemogcmcoup_update( mype, npes, icomm, &
   &                           npoints, pgutau, pgvtau, &
   &                           pgqsr, pgqns, pgemp, kt, ldebug )

   ! Update fluxes in nemogcmcoup_data by parallel
   ! interpolation of the input gaussian grid data
   
   USE par_kind

   IMPLICIT NONE

   ! Arguments

   ! MPI communications
   INTEGER, INTENT(IN) :: mype,npes,icomm
   ! Fluxes on the Gaussian grid.
   INTEGER, INTENT(IN) :: npoints
   REAL(wpIFS), DIMENSION(npoints), intent(IN) :: &
      & pgutau, pgvtau, pgqsr, pgqns, pgemp
   ! Current time step
   INTEGER, INTENT(in) :: kt
   ! Write debugging fields in netCDF
   LOGICAL, INTENT(IN) :: ldebug

   ! Local variables

   WRITE(0,*)'nemogcmcoup_update should be called with with.'
   CALL abort
   
END SUBROUTINE nemogcmcoup_update

SUBROUTINE nemogcmcoup_update_add( mype, npes, icomm, &
   &                               npoints, pgsst, pgtsk, kt, ldebug )

   ! Update addetiona in nemogcmcoup_data by parallel
   ! interpolation of the input gaussian grid data
   
   USE par_kind
   USE fesom_main_storage_module, only: fesom => f ! only: MPI_COMM_FESOM, mype (previously in g_parsup)

   IMPLICIT NONE

   ! Arguments

   ! MPI communications
   INTEGER, INTENT(IN) :: mype,npes,icomm
   ! Input on the Gaussian grid.
   INTEGER, INTENT(IN) :: npoints
   REAL(wpIFS), DIMENSION(npoints), intent(IN) :: &
      & pgsst, pgtsk
   ! Current time step
   INTEGER, INTENT(in) :: kt
   ! Write debugging fields in netCDF
   LOGICAL, INTENT(IN) :: ldebug

   ! Local variables

   if(fesom%mype==0) then
   WRITE(0,*)'In nemogcmcoup_update_add FESOM dummy routine. Proceeding...'
   !CALL abort
   endif   

END SUBROUTINE nemogcmcoup_update_add


SUBROUTINE nemogcmcoup_wam_coupinit( mype, npes, icomm, &
   &                                 nlocpoints, nglopoints, &
   &                                 nlocmsk, ngloind, iunit )

   ! Initialize single executable coupling between WAM and NEMO
   ! This is called from WAM.

   IMPLICIT NONE

   ! Input arguments

   ! Message passing information
   INTEGER, INTENT(IN) :: mype,npes,icomm
   ! WAM grid information   
   ! Number of local and global points
   INTEGER, INTENT(IN) :: nlocpoints, nglopoints
   ! Integer mask and global indices
   INTEGER, DIMENSION(nlocpoints), INTENT(IN) :: nlocmsk, ngloind
   ! Unit for output in parinter_init
   INTEGER :: iunit
   
   WRITE(0,*)'Wam coupling not implemented for FESOM'
   CALL abort

END SUBROUTINE nemogcmcoup_wam_coupinit


SUBROUTINE nemogcmcoup_wam_get( mype, npes, icomm, &
   &                            nopoints, pwsst, pwicecov, pwicethk, &
   &                            pwucur, pwvcur, licethk )

   ! Interpolate from the ORCA grid
   ! to the WAM grid. 
   
   ! This routine can be called at any point in time since it does
   ! the necessary message passing in parinter_fld. 

   USE par_kind
   IMPLICIT NONE
   
   ! Arguments

   ! Message passing information
   INTEGER, INTENT(IN) :: mype, npes, icomm
   ! Number WAM grid points
   INTEGER, INTENT(IN) :: nopoints
   ! Local arrays of sst, ice cover, ice thickness and currents
   REAL(wpIFS), DIMENSION(nopoints) :: pwsst, pwicecov, pwicethk, pwucur, pwvcur
   LOGICAL :: licethk

   ! Local variables

   WRITE(0,*)'nemogcmcoup_wam_get should not be called when coupling to fesom.'
   CALL abort

END SUBROUTINE nemogcmcoup_wam_get


SUBROUTINE nemogcmcoup_wam_update( mype, npes, icomm, &
   &                               npoints, pwswh, pwmwp, &
   &                               pwphioc, pwtauoc, pwstrn, &
   &                               pwustokes, pwvstokes, &
   &                               cdtpro, ldebug )

   ! Update fluxes in nemogcmcoup_data by parallel
   ! interpolation of the input WAM grid data
   
   USE par_kind

   IMPLICIT NONE

   ! Arguments

   ! MPI communications
   INTEGER, INTENT(IN) :: mype,npes,icomm
   ! Data on the WAM grid.
   INTEGER, INTENT(IN) :: npoints
   REAL(wpIFS), DIMENSION(npoints), INTENT(IN) :: &
      & pwswh, pwmwp, pwphioc, pwtauoc, pwstrn, pwustokes, pwvstokes
   ! Current time
   CHARACTER(len=14), INTENT(IN) :: cdtpro
   ! Write debugging fields in netCDF
   LOGICAL, INTENT(IN) :: ldebug

   ! Local variables
   
   WRITE(0,*)'nemogcmcoup_wam_update should not be called when coupling to fesom.'
   CALL abort

END SUBROUTINE nemogcmcoup_wam_update


SUBROUTINE nemogcmcoup_wam_update_stress( mype, npes, icomm, npoints, &
   &                                      pwutau, pwvtau, pwuv10n, pwphif,&
   &                                      cdtpro, ldebug )

   ! Update stresses in nemogcmcoup_data by parallel
   ! interpolation of the input WAM grid data
   
   USE par_kind

   IMPLICIT NONE

   ! Arguments

   ! MPI communications
   INTEGER, INTENT(IN) :: mype,npes,icomm
   ! Data on the WAM grid.
   INTEGER, INTENT(IN) :: npoints
   REAL(wpIFS), DIMENSION(npoints), INTENT(IN) :: &
      & pwutau, pwvtau, pwuv10n, pwphif
   ! Current time step
   CHARACTER(len=14), INTENT(IN) :: cdtpro
   ! Write debugging fields in netCDF
   LOGICAL, INTENT(IN) :: ldebug

   ! Local variables
   
   WRITE(0,*)'nemogcmcoup_wam_update_stress should not be called when coupling to fesom.'
   CALL abort

END SUBROUTINE nemogcmcoup_wam_update_stress
