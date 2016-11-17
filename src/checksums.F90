!BOP
!
! !ROUTINE: checksums --- compute checksums to test forecast integrity
!
! !INTERFACE
SUBROUTINE checksums

! !DESCRIPTION:
! This routine compute checksums of the different fields of FESOM
! in order to check for deviations in between diferent model runs.
!
! !USES
  USE o_ARRAYS, only: eta_n, tr_arr
  USE g_PARSUP, only: myDIM_nod2d, MPI_DOUBLE_PRECISION, MPI_SUM, &
       MPI_COMM_WORLD, MPIerr, MPI_COMM_WORLD, mype
  USE o_MESH, only: nl, WP

  IMPLICIT NONE
!EOP

! Local variables
  INTEGER :: i,j                        ! Counters
  REAL(kind=WP) :: sum_ssh_p, sum_ssh   ! Sum over SSH field
  REAL(kind=WP) :: sum_temp_p, sum_temp ! Sum over temperature
  REAL(kind=WP) :: sum_salt_p, sum_salt ! Sum over salinity 

! *** SSH ***
  sum_ssh_p = 0.0
  DO i = 1, myDIM_nod2d
     sum_ssh_p = sum_ssh_p + eta_n(i)
  END DO

  CALL MPI_Reduce(sum_ssh_p, sum_ssh, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
          0, MPI_COMM_WORLD, MPIerr)


! *** Temperature ***
  sum_temp_p = 0.0
  DO i = 1, myDIM_nod2d
     DO j = 1, nl-1
        sum_temp_p = sum_temp_p + tr_arr(j, i, 1)
     END DO
  END DO

  CALL MPI_Reduce(sum_temp_p, sum_temp, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
       0, MPI_COMM_WORLD, MPIerr)

! *** Salinity ***
   sum_salt_p = 0.0
   DO i = 1, myDIM_nod2d
      DO j = 1, nl-1
         sum_salt_p = sum_salt_p + tr_arr(j, i, 2)
      END DO
   END DO
   CALL MPI_Reduce(sum_salt_p, sum_salt, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
           0, MPI_COMM_WORLD, MPIerr)


! *** Output ***
  IF (mype == 0) THEN
     write (*, *) ''
     write (*, *) '======== CHECKSUMS ========='
     write (*, *) 'ssh:         ',sum_ssh
     write (*, *) 'Temperature: ',sum_temp
     write (*, *) 'Salinity:    ',sum_salt
  END IF



END SUBROUTINE checksums
