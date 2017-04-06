!==========================================================================
SUBROUTINE write_step_info(istep)
use g_config, only: logfile_outfreq, dt
USE o_MESH
USE o_PARAM
USE g_PARSUP
USE o_ARRAYS
USE g_comm_auto
IMPLICIT NONE

integer                                   :: n, elnodes(3), istep
real(kind=8)                              :: int_eta, int_hbar, int_deta, int_dhbar, int_wflux
real(kind=8)                              :: min_eta, min_hbar, min_wflux
real(kind=8)                              :: max_eta, max_hbar, max_wflux
real(kind=8)                              :: loc_eta, loc_hbar, loc_deta, loc_dhbar, loc_wflux

	if (mod(istep,logfile_outfreq)==0) then
		!_______________________________________________________________________
		int_eta   =0.
		int_hbar  =0.
		int_deta  =0.
		int_dhbar =0.
		int_wflux =0.
		loc_eta   =0.
		loc_hbar  =0.
		
		!_______________________________________________________________________
		do n=1, myDim_nod2D
			loc_eta   = loc_eta   + area(1, n)*eta_n(n)
			loc_hbar  = loc_hbar  + area(1, n)*hbar(n)
			loc_wflux = loc_wflux + area(1, n)*water_flux(n)*dt*(-1.0)
			loc_deta  = loc_deta  + area(1, n)*d_eta(n)
			loc_dhbar = loc_dhbar + area(1, n)*(hbar(n)-hbar_old(n))
		end do
		
		!_______________________________________________________________________
		call MPI_AllREDUCE(loc_eta  , int_eta  , 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, MPIerr)
		call MPI_AllREDUCE(loc_hbar , int_hbar , 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, MPIerr)
		call MPI_AllREDUCE(loc_wflux, int_wflux, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, MPIerr)
		call MPI_AllREDUCE(loc_deta , int_deta , 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, MPIerr)
		call MPI_AllREDUCE(loc_dhbar, int_dhbar, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, MPIerr)
		
		!_______________________________________________________________________
		loc_eta= minval(eta_n(1:myDim_nod2D))
		loc_hbar=minval(hbar(1:myDim_nod2D))
		loc_wflux=minval(water_flux(1:myDim_nod2D))
		call MPI_AllREDUCE(loc_eta  , min_eta  , 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, MPIerr)
		call MPI_AllREDUCE(loc_hbar , min_hbar , 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, MPIerr)
		call MPI_AllREDUCE(loc_wflux, min_wflux, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, MPIerr)
		
		!_______________________________________________________________________
		loc_eta= maxval(eta_n(1:myDim_nod2D))
		loc_hbar=maxval(hbar(1:myDim_nod2D))
		loc_wflux=maxval(water_flux(1:myDim_nod2D))
		call MPI_AllREDUCE(loc_eta  , max_eta  , 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, MPIerr)
		call MPI_AllREDUCE(loc_hbar , max_hbar , 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, MPIerr)
		call MPI_AllREDUCE(loc_wflux, max_wflux, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, MPIerr)
		
		!_______________________________________________________________________
		if (mype==0) then
			write(*,*) '	___global estimat of eta & hbar --> mstep=',mstep,'____'
			write(*,*) '	 int(eta), int(hbar)      =', int_eta, int_hbar
			write(*,*) '	  --> error(eta-hbar)     =', int_eta-int_hbar
			write(*,*) '	 min(eta) , max(eta)      =', min_eta, max_eta
			write(*,*) '	 max(hbar), max(hbar)     =', min_hbar, max_hbar
			write(*,*)
			write(*,*) '	 int(deta), int(dhbar)    =', int_deta, int_dhbar
			write(*,*) '	  --> error(deta-dhbar)   =', int_deta-int_dhbar
			write(*,*) '	  --> error(deta-wflux)   =', int_deta-int_wflux
			write(*,*) '	  --> error(dhbar-wflux)  =', int_dhbar-int_wflux
			write(*,*)
			write(*,*) '	 -int(wflux)*dt           =', int_wflux
			write(*,*) '	 int(deta )-int(wflux)*dt =', int_deta-int_wflux
			write(*,*) '	 int(dhbar)-int(wflux)*dt =', int_dhbar-int_wflux
			write(*,*) '	 min(wflux), max(wflux)   =', min_wflux, max_wflux
			write(*,*)
		endif
	endif
END SUBROUTINE write_step_info

