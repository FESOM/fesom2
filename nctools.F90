#define __MYFILE__ 'nctools.F90'
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
