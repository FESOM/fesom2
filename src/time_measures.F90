MODULE time_measures
!------------------------------------------------------------------------------
!
! Description:
!   This stand-alone-module provides service utilities for timing different 
!   parts a program.
!
!   Routines (module procedures) currently contained:
!
!     - timer_init:
!       initialize the variables for time-measuring and mpi if not done
!       so far
!
!     - timer:
!       Gets the time used for a special part of the program
!
!     - timer_finalize:
!       At the end of the program all timings are collected to one node and 
!       printed.
!
! Method:
!   timer has to be called for special parts of the program. The special
!   parts are characterized by numbers and optionally by strings, which are 
!   passed to timer. Before timer is called the variables necessary 
!   for time measuring have to be initialized by calling timer_init. At the
!   end of the program timer_finalize collects all measurements to
!   node 0 and prints them to stdout.
!
! Usage:
!
!   PROGRAM test
!
!     USE timer, ONLY: timer_init, timer, timer_finalize
!
!     <Initilization stuff>
!
!     CALL timer_init
!     CALL timer( 1, 'start', NAME = 'Total' )
!
!     <Doing some computation>
!
!     DO i = 0, n
!
!        <Doing some computation>
!
!        CALL timer( 2, 'start', NAME = 'Loop part1' )  
!
!        <Doing some computation>
!  
!        CALL timer( 2, 'stop', NAME = 'Loop part1' )  
!
!        <Doing some computation>
!
!     ENDDO
!
!     <Doing some more computation>
!
!     CALL timer( 1, 'stop', NAME = 'Total' )
!     CALL timer_finalize
!
!   END PROGRAM test
!
!  Author: Michael Schroeter, AWI
!  phone:  +49  471 4831 2084
!  email:  mschroeter@awi-bremerhaven.de
!
! History:
! Version    Date       Name
! ---------- ---------- ----
! 1.0        2005/06/08 Michael Schroeter
!  Initial release
! !VERSION!  !DATE!     <Your name>
!  <Modification comments>
!
!==============================================================================

!
!-- None variable must be implicit
    IMPLICIT NONE

    INCLUDE "mpif.h"
!
!-- number of desired significant digits for
!-- real variables corresponds to 8 byte real variables
    INTEGER, PARAMETER  :: ireals    = SELECTED_REAL_KIND (12,200), &
                           ntab      = 100

    INTEGER             :: cpu_start_time, cpu_my_id, cpu_n_pes
    INTEGER             :: ierr, imp_reals, irate, imax

    LOGICAL             :: cpu_initialized = .FALSE.

!
!-- 
    TYPE cpu_measurements

       CHARACTER(LEN=128) :: name
       LOGICAL            :: initialized
       REAL(KIND=ireals)  :: start_time, time, sum_time

    END TYPE cpu_measurements

    TYPE( cpu_measurements ), DIMENSION(ntab) :: timings

    INTERFACE timer
       MODULE PROCEDURE  timer_sub
    END INTERFACE

    INTERFACE timer_finalize
       MODULE PROCEDURE  timer_finalize_sub
    END INTERFACE

    SAVE

    CONTAINS

!******************************************************************************
!------------------------------------------------------------------------------

      SUBROUTINE timer_init

        INTEGER :: icount

        LOGICAL :: mpi_flag

!------------------------------------------------------------------------------
!- End of header -
!------------------------------------------------------------------------------

!
!--     get system clock settings
        CALL SYSTEM_CLOCK( COUNT      = icount, &
                           COUNT_RATE = irate,  &
                           COUNT_MAX  = imax    )
        cpu_start_time  = icount
        cpu_initialized = .TRUE.

!
!--     Initialize the timing
        timings(:)%initialized = .FALSE.
        timings(:)%start_time  = 0.0
        timings(:)%time        = 0.0
        timings(:)%sum_time    = 0.0

!
!--     If not already done, initialize mpi
        CALL MPI_INITIALIZED( mpi_flag, ierr )
        IF ( .NOT. mpi_flag ) CALL MPI_INIT( ierr )

!
!--     This module holds his own PE-informations
        CALL MPI_COMM_SIZE( MPI_COMM_WORLD, cpu_n_pes, ierr )
        CALL MPI_COMM_RANK( MPI_COMM_WORLD, cpu_my_id, ierr )

!
!--     Determine the type of REALs for MPI and other variables
!--     If the KIND-type parameters in data_parameters are changed, the
!--     variables here have to be changed accordingly.
!--     Model Real variables
!         IF ( KIND( 1.0 )       == ireals ) THEN
!            write (*,*) 'MPI_REAL, ireals', ireals
!            imp_reals = MPI_REAL
!         ELSEIF ( KIND( 1.0D0 ) == ireals ) THEN
!            write (*,*) 'MPI_DOUBLE_PRECISION'
           imp_reals = MPI_DOUBLE_PRECISION
!         ELSEIF ( KIND( 1.0 ) == 8 .AND. ireals == 4 ) THEN
!            ! it seems that this is a T3E where 4 Byte REALs are used
!            write (*,*) 'MPI_REAL4'
!            imp_reals = MPI_REAL4
!         ELSE
!            PRINT*, '### TIMER: cannot find MPI REAL Type for ireals'
!            STOP
!         ENDIF

      END SUBROUTINE timer_init

!******************************************************************************
!------------------------------------------------------------------------------

      SUBROUTINE timer_sub( itab, switch, name, barrier )

        CHARACTER(LEN=*), OPTIONAL :: barrier
        CHARACTER(LEN=*), OPTIONAL :: name
        CHARACTER(LEN=*)           :: switch
        CHARACTER(LEN=3)           :: suffix

        INTEGER            :: icount, itab

        REAL(KIND=ireals)  :: time

        LOGICAL :: wait = .FALSE.

!------------------------------------------------------------------------------
!- End of header -
!------------------------------------------------------------------------------

!
!--     Do some checks
        IF ( .NOT. cpu_initialized ) THEN
           PRINT*, '+++ Time measurement not initialized: '
           STOP
        ENDIF

!
!--     To many measurement points?
        IF ( itab > ntab ) THEN
           PRINT*, '+++ Time measurement: Increase ntab (> ', ntab, ')'
           STOP
        ENDIF

!
!--     Wrong name?
        IF ( timings(itab)%initialized ) THEN
           IF ( PRESENT( name ) .AND. &
                TRIM( name ) /= TRIM( timings(itab)%name ) ) THEN
              PRINT*, '+++ Time measurement: Incorrect name at point', itab
           ENDIF
        ENDIF

        IF ( .NOT. timings(itab)%initialized ) THEN
           IF ( PRESENT( name ) ) THEN
              timings(itab)%name = TRIM( name )
           ELSE
              WRITE( suffix, '(I3.3)' ) itab
              timings(itab)%name = 'measure point ' // suffix
           ENDIF
        ENDIF

!
!--     Check, whether system clock is present and initialize icountsold
        CALL SYSTEM_CLOCK( COUNT = icount  )
        IF ( icount >= cpu_start_time ) THEN
           time = ( REAL( icount - cpu_start_time ) ) &
                / REAL( irate, ireals )
        ELSE
           time = REAL( imax - ( cpu_start_time - icount ), ireals ) &
                / REAL( irate, ireals )
        ENDIF

        IF ( ( TRIM( switch ) == 'start' ) .OR. &
             ( TRIM( switch ) == 'START' ) ) THEN 

           timings(itab)%start_time  = time
           timings(itab)%initialized = .TRUE.
 
        ELSE IF ( ( TRIM( switch ) == 'stop' ) .OR. &
             ( TRIM( switch ) == 'STOP' ) ) THEN

           IF ( .NOT. timings(itab)%initialized ) THEN 
              PRINT*, '+++ Time measurement ', itab, ' not initialized'
              STOP
           ENDIF

           timings(itab)%time     = time - timings(itab)%start_time
           timings(itab)%sum_time = timings(itab)%sum_time + timings(itab)%time

        ELSE

           PRINT*, '+++ unknown switch of time measurement: ', switch
           STOP

        ENDIF

!
!--     Set a barrier?
        IF ( PRESENT( barrier ) ) THEN
           IF ( TRIM( barrier ) == 'no_wait' ) wait = .FALSE.
           IF ( TRIM( barrier ) == 'wait' ) wait = .TRUE.
        ENDIF

!
!--     Set the barrier
        IF ( wait ) CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )



      END SUBROUTINE timer_sub

!******************************************************************************
!------------------------------------------------------------------------------

      SUBROUTINE timer_finalize_sub( headerline )

        CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: headerline
	CHARACTER(LEN=94)      :: header

        INTEGER, DIMENSION(1)  :: imax, imin
        INTEGER                :: h_len, itab, nrecv, nroot, nsend

        REAL(KIND=ireals)  :: abs_time, max_time, min_time, rel_time, var_time
        REAL(KIND=ireals), ALLOCATABLE, DIMENSION(:,:) :: times_global

!------------------------------------------------------------------------------
!- End of header -
!------------------------------------------------------------------------------

        ALLOCATE( times_global(1:ntab,0:cpu_n_pes-1) )

        nsend = ntab
        nrecv = nsend
        nroot = 0

!
!--     Gather local times from all nodes
        IF ( cpu_n_pes > 1 ) THEN
           CALL MPI_GATHER( timings%sum_time, nsend, imp_reals, &
                            times_global,     nrecv, imp_reals, &
                            nroot, MPI_COMM_WORLD, ierr )
        ELSE
           times_global(:,0) = timings(:)%sum_time
        ENDIF

!
!--     Do some output only for PE 0
        IF ( cpu_my_id == 0 ) THEN
!
!--        Set the headerline
	   IF ( PRESENT( headerline ) ) THEN
	        header(1:94) = ' '
                h_len = LEN_TRIM( headerline )
		IF ( h_len > 94 ) h_len = 94
	        header(1:h_len) = headerline(1:h_len)
	   ELSE
	   	header = 'Time measurements'
	   ENDIF

           WRITE(*,1030) header
           WRITE(*,1020)
           WRITE(*,1025) cpu_n_pes
           WRITE(*,1020)
           WRITE(*,1010)
           itab = 1
           DO WHILE ( itab <= ntab )
              
              IF ( timings(itab)%initialized ) THEN
                 abs_time = SUM( times_global(itab,:) ) /     &
                           REAL( cpu_n_pes, ireals )
                 rel_time = abs_time /                  &
                          ( SUM( times_global(1,:) )  / &
                            REAL( cpu_n_pes, ireals ) )
                 var_time = SQRT( 1.0 / ( cpu_n_pes-1 ) *           &
                      SUM( ( times_global(itab,:) - abs_time )**2 ) &
                                )
                 max_time = MAXVAL( times_global(itab,:) )
                 min_time = MINVAL( times_global(itab,:) )
                 imax     = MAXLOC( times_global(itab,:) )
                 imin     = MINLOC( times_global(itab,:) )
                 WRITE(*,1000) itab, &
                      TRIM( timings(itab)%name ),           &
                      abs_time, rel_time * 100.0, var_time, &
                      max_time, imax(1), min_time, imin(1)
              END IF
              itab = itab + 1

           END DO
           WRITE(*,990)

        ENDIF

        DEALLOCATE( times_global )

!                  1234567890123456789012345678901234567890123456789012345
1030    FORMAT( / '#',  95('-'), '#' &
                / '| ', A94,     '|' )
!                  1020
1025    FORMAT(   '| on ', I5.5, ' PEs', 82X,       '|' )
1020    FORMAT(   '#', 95('-'),                     '#' )
1010    FORMAT(   '| Nr.  Location', 25X, 'time', 7X, '%', 8X, 'var      max  at i1      min  at i2 |' &
                / '|              ', 25X, ' (s)', 7X, '         (s)       (s)             (s)       |' )
1000    FORMAT(   '| ', I3, 1X, A25, ': ', 3X, F9.2, 3X, F5.1, 1X, F9.2, 1X, F9.2, 1X,  &
                        I5, 1X, F9.2, 1X, I5,      ' |' )
990     FORMAT(   '#', 95('-'),                    '#' /  )
        
      END SUBROUTINE timer_finalize_sub

END MODULE time_measures
