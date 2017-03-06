
!==========================================================================
SUBROUTINE write_step_info(istep)
use g_config, only: logfile_outfreq
USE o_MESH
USE o_PARAM
USE g_PARSUP
USE o_ARRAYS
USE g_comm_auto
IMPLICIT NONE

integer                                   :: n, elnodes(3), istep
real(kind=8)                              :: int_eta, int_hbar
real(kind=8)                              :: min_eta, min_hbar
real(kind=8)                              :: max_eta, max_hbar
real(kind=8)                              :: loc_eta, loc_hbar

int_eta=0.
int_hbar=0.
loc_eta=0.
loc_hbar=0.

  do n=1, myDim_nod2D
     loc_eta =loc_eta +area(1, n)*eta_n(n)
     loc_hbar=loc_hbar+area(1, n)*hbar(n)
  end do

  call MPI_AllREDUCE(loc_eta,  int_eta,  1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, MPIerr)
  call MPI_AllREDUCE(loc_hbar, int_hbar, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, MPIerr)

  loc_eta= minval(eta_n(1:myDim_nod2D))
  loc_hbar=minval(hbar(1:myDim_nod2D))

  call MPI_AllREDUCE(loc_eta,  min_eta,  1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, MPIerr)
  call MPI_AllREDUCE(loc_hbar, min_hbar, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, MPIerr)

  loc_eta= maxval(eta_n(1:myDim_nod2D))
  loc_hbar=maxval(hbar(1:myDim_nod2D))

  call MPI_AllREDUCE(loc_eta,  max_eta,  1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, MPIerr)
  call MPI_AllREDUCE(loc_hbar, max_hbar, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, MPIerr)

if (mod(istep,logfile_outfreq)==0 .and. mype==0) then
    write(*,*) 'at step             =', mstep
    write(*,*) 'integrated eta/hbar =', int_eta, int_hbar
    write(*,*) 'the error is        =', int_hbar-int_eta
    write(*,*) 'min eta/hbar        =', min_eta, min_hbar
    write(*,*) 'max eta/hbar        =', max_eta, max_hbar
endif
END SUBROUTINE write_step_info

