module write_step_info_interface
  interface
    subroutine write_step_info(istep,outfreq, mesh)
      use MOD_MESH
      integer								:: istep,outfreq
      type(t_mesh), intent(in)                               , target :: mesh
    end subroutine
  end interface
end module

!
!
!===============================================================================
subroutine write_step_info(istep,outfreq, mesh)
	use g_config, only: dt, use_ice
	use MOD_MESH
	use o_PARAM
	use g_PARSUP
	use o_ARRAYS
	use i_ARRAYS
	use g_comm_auto
	implicit none
	
	integer								:: n, istep,outfreq
	real(kind=WP)						:: int_eta, int_hbar, int_wflux, int_hflux, int_temp, int_salt
	real(kind=WP)						:: min_eta, min_hbar, min_wflux, min_hflux, min_temp, min_salt, &
                                           min_wvel,min_hnode,min_deta,min_wvel2,min_hnode2, &
                                           min_vvel, min_vvel2, min_uvel, min_uvel2
	real(kind=WP)						:: max_eta, max_hbar, max_wflux, max_hflux, max_temp, max_salt, &
                                           max_wvel, max_hnode, max_deta, max_wvel2, max_hnode2, max_m_ice, &
                                           max_vvel, max_vvel2, max_uvel, max_uvel2, &
                                           max_cfl_z, max_pgfx, max_pgfy, max_kv, max_av 
	real(kind=WP)						:: int_deta , int_dhbar
	real(kind=WP)						:: loc, loc_eta, loc_hbar, loc_deta, loc_dhbar, loc_wflux,loc_hflux, loc_temp, loc_salt
    type(t_mesh), intent(in)                               , target :: mesh
#include "associate_mesh.h"
	if (mod(istep,outfreq)==0) then
		
		!_______________________________________________________________________
		int_eta   =0.
		int_hbar  =0.
		int_deta  =0.
		int_dhbar =0.
		int_wflux =0.
		int_hflux =0.
		int_temp  =0.
		int_salt  =0.
		loc_eta   =0.
		loc_hbar  =0.
		loc_deta  =0.
		loc_dhbar =0.
		loc_wflux =0.
!!PS 		loc_hflux =0.
!!PS 		loc_temp  =0.
!!PS 		loc_salt  =0.
		loc       =0.
		!_______________________________________________________________________
		do n=1, myDim_nod2D
            if (ulevels_nod2D(n)>1) cycle
			loc_eta   = loc_eta   + area(ulevels_nod2D(n), n)*eta_n(n)
			loc_hbar  = loc_hbar  + area(ulevels_nod2D(n), n)*hbar(n)
			loc_deta  = loc_deta  + area(ulevels_nod2D(n), n)*d_eta(n)
			loc_dhbar = loc_dhbar + area(ulevels_nod2D(n), n)*(hbar(n)-hbar_old(n))
			loc_wflux = loc_wflux + area(ulevels_nod2D(n), n)*water_flux(n)
!!PS 			loc_hflux = loc_hflux + area(1, n)*heat_flux(n)
!!PS 			loc_temp  = loc_temp  + area(1, n)*sum(tr_arr(:,n,1))/(nlevels_nod2D(n)-1)
!!PS 			loc_salt  = loc_salt  + area(1, n)*sum(tr_arr(:,n,2))/(nlevels_nod2D(n)-1)
		end do
		
		!_______________________________________________________________________
		call MPI_AllREDUCE(loc_eta  , int_eta  , 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_FESOM, MPIerr)
		call MPI_AllREDUCE(loc_hbar , int_hbar , 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_FESOM, MPIerr)
		call MPI_AllREDUCE(loc_deta , int_deta , 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_FESOM, MPIerr)
		call MPI_AllREDUCE(loc_dhbar, int_dhbar, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_FESOM, MPIerr)
		call MPI_AllREDUCE(loc_wflux, int_wflux, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_FESOM, MPIerr)
!!PS 		call MPI_AllREDUCE(loc_hflux, int_hflux, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_FESOM, MPIerr)
!!PS 		call MPI_AllREDUCE(loc_temp , int_temp , 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_FESOM, MPIerr)
!!PS 		call MPI_AllREDUCE(loc_salt , int_salt , 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_FESOM, MPIerr)

		int_eta  = int_eta  /ocean_area
		int_hbar = int_hbar /ocean_area
		int_deta = int_deta /ocean_area
		int_dhbar= int_dhbar/ocean_area
		int_wflux= int_wflux/ocean_area
		
!!PS 		int_eta  = int_eta  /ocean_areawithcav
!!PS 		int_hbar = int_hbar /ocean_areawithcav
!!PS 		int_deta = int_deta /ocean_areawithcav
!!PS 		int_dhbar= int_dhbar/ocean_areawithcav
!!PS 		int_wflux= int_wflux/ocean_areawithcav
		
!!PS 		int_hflux= int_hflux/ocean_area
!!PS 		int_temp = int_temp /ocean_area
!!PS 		int_salt = int_salt /ocean_area
		
		!_______________________________________________________________________
		loc = minval(eta_n(1:myDim_nod2D))
		call MPI_AllREDUCE(loc , min_eta  , 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_FESOM, MPIerr)
		loc = minval(hbar(1:myDim_nod2D))
		call MPI_AllREDUCE(loc , min_hbar , 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_FESOM, MPIerr)
		loc = minval(water_flux(1:myDim_nod2D))
		call MPI_AllREDUCE(loc , min_wflux, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_FESOM, MPIerr)
		loc = minval(heat_flux(1:myDim_nod2D))
		call MPI_AllREDUCE(loc , min_hflux, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_FESOM, MPIerr)
		loc = minval(tr_arr(:,1:myDim_nod2D,1),MASK=(tr_arr(:,1:myDim_nod2D,2)/=0.0))
		call MPI_AllREDUCE(loc , min_temp , 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_FESOM, MPIerr)
		loc = minval(tr_arr(:,1:myDim_nod2D,2),MASK=(tr_arr(:,1:myDim_nod2D,2)/=0.0))
		call MPI_AllREDUCE(loc , min_salt , 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_FESOM, MPIerr)
		loc = minval(Wvel(1,1:myDim_nod2D))
		call MPI_AllREDUCE(loc , min_wvel , 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_FESOM, MPIerr)
		loc = minval(Wvel(2,1:myDim_nod2D))
		call MPI_AllREDUCE(loc , min_wvel2 , 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_FESOM, MPIerr)
		loc = minval(Unode(1,1,1:myDim_nod2D))
		call MPI_AllREDUCE(loc , min_uvel , 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_FESOM, MPIerr)
		loc = minval(Unode(1,2,1:myDim_nod2D))
		call MPI_AllREDUCE(loc , min_uvel2 , 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_FESOM, MPIerr)
		loc = minval(Unode(2,1,1:myDim_nod2D))
		call MPI_AllREDUCE(loc , min_vvel , 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_FESOM, MPIerr)
		loc = minval(Unode(2,2,1:myDim_nod2D))
		call MPI_AllREDUCE(loc , min_vvel2 , 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_FESOM, MPIerr)
		loc = minval(d_eta(1:myDim_nod2D))
		call MPI_AllREDUCE(loc , min_deta  , 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_FESOM, MPIerr)
		loc = minval(hnode(1,1:myDim_nod2D),MASK=(hnode(1,1:myDim_nod2D)/=0.0))
		call MPI_AllREDUCE(loc , min_hnode , 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_FESOM, MPIerr)
		loc = minval(hnode(2,1:myDim_nod2D),MASK=(hnode(2,1:myDim_nod2D)/=0.0))
		call MPI_AllREDUCE(loc , min_hnode2 , 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_FESOM, MPIerr)
		
		!_______________________________________________________________________
		loc = maxval(eta_n(1:myDim_nod2D))
		call MPI_AllREDUCE(loc , max_eta  , 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_FESOM, MPIerr)
		loc = maxval(hbar(1:myDim_nod2D))
		call MPI_AllREDUCE(loc , max_hbar , 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_FESOM, MPIerr)
		loc = maxval(water_flux(1:myDim_nod2D))
		call MPI_AllREDUCE(loc , max_wflux, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_FESOM, MPIerr)
		loc = maxval(heat_flux(1:myDim_nod2D))
		call MPI_AllREDUCE(loc , max_hflux, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_FESOM, MPIerr)
		loc = maxval(tr_arr(:,1:myDim_nod2D,1),MASK=(tr_arr(:,1:myDim_nod2D,2)/=0.0))
		call MPI_AllREDUCE(loc , max_temp , 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_FESOM, MPIerr)
		loc = maxval(tr_arr(:,1:myDim_nod2D,2),MASK=(tr_arr(:,1:myDim_nod2D,2)/=0.0))
		call MPI_AllREDUCE(loc , max_salt , 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_FESOM, MPIerr)
		loc = maxval(Wvel(1,1:myDim_nod2D))
		call MPI_AllREDUCE(loc , max_wvel , 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_FESOM, MPIerr)
		loc = maxval(Wvel(2,1:myDim_nod2D))
		call MPI_AllREDUCE(loc , max_wvel2 , 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_FESOM, MPIerr)
		loc = maxval(Unode(1,1,1:myDim_nod2D))
		call MPI_AllREDUCE(loc , max_uvel , 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_FESOM, MPIerr)
		loc = maxval(Unode(1,2,1:myDim_nod2D))
		call MPI_AllREDUCE(loc , max_uvel2 , 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_FESOM, MPIerr)
		loc = maxval(Unode(2,1,1:myDim_nod2D))
		call MPI_AllREDUCE(loc , max_vvel , 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_FESOM, MPIerr)
		loc = maxval(Unode(2,2,1:myDim_nod2D))
		call MPI_AllREDUCE(loc , max_vvel2 , 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_FESOM, MPIerr)
		loc = maxval(d_eta(1:myDim_nod2D))
		call MPI_AllREDUCE(loc , max_deta  , 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_FESOM, MPIerr)
		loc = maxval(hnode(1,1:myDim_nod2D),MASK=(hnode(1,1:myDim_nod2D)/=0.0))
		call MPI_AllREDUCE(loc , max_hnode , 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_FESOM, MPIerr)
		loc = maxval(hnode(2,1:myDim_nod2D),MASK=(hnode(2,1:myDim_nod2D)/=0.0))
		call MPI_AllREDUCE(loc , max_hnode2 , 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_FESOM, MPIerr)
		loc = maxval(CFL_z(:,1:myDim_nod2D))
		call MPI_AllREDUCE(loc , max_cfl_z , 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_FESOM, MPIerr)
		loc = maxval(abs(pgf_x(:,1:myDim_nod2D)))
		call MPI_AllREDUCE(loc , max_pgfx , 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_FESOM, MPIerr)
		loc = maxval(abs(pgf_y(:,1:myDim_nod2D)))
		call MPI_AllREDUCE(loc , max_pgfy , 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_FESOM, MPIerr)
                if (use_ice) then
                   loc = maxval(m_ice(1:myDim_nod2D))
                   call MPI_AllREDUCE(loc , max_m_ice , 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_FESOM, MPIerr)
                end if
		loc = maxval(abs(Av(:,1:myDim_nod2D)))
		call MPI_AllREDUCE(loc , max_av , 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_FESOM, MPIerr)
		loc = maxval(abs(Kv(:,1:myDim_nod2D)))
		call MPI_AllREDUCE(loc , max_kv , 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_FESOM, MPIerr)
		!_______________________________________________________________________
		if (mype==0) then
			write(*,*) '___CHECK GLOBAL OCEAN VARIABLES --> mstep=',mstep
			write(*,*) '	___global estimat of eta & hbar____________________'
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
			write(*,*) '	 -int(wflux)*dt           =', int_wflux*dt*(-1.0)
			write(*,*) '	 int(deta )-int(wflux)*dt =', int_deta-int_wflux*dt*(-1.0)
			write(*,*) '	 int(dhbar)-int(wflux)*dt =', int_dhbar-int_wflux*dt*(-1.0)
			write(*,*)
			write(*,*) '	___global min/max/mean  --> mstep=',mstep,'____________'
			write(*,"(A, ES10.3, A, ES10.3, A, A     )") ' 	       eta= ', min_eta  ,' | ',max_eta  ,' | ','N.A.'
			write(*,"(A, ES10.3, A, ES10.3, A, A     )") ' 	      deta= ', min_deta ,' | ',max_deta ,' | ','N.A.'
			write(*,"(A, ES10.3, A, ES10.3, A, A     )") ' 	      hbar= ', min_hbar ,' | ',max_hbar ,' | ','N.A.'
			write(*,"(A, ES10.3, A, ES10.3, A, ES10.3)") ' 	     wflux= ', min_wflux,' | ',max_wflux,' | ',int_wflux
			write(*,"(A, ES10.3, A, ES10.3, A, ES10.3)") ' 	     hflux= ', min_hflux,' | ',max_hflux,' | ',int_hflux
			write(*,"(A, ES10.3, A, ES10.3, A, ES10.3)") ' 	      temp= ', min_temp ,' | ',max_temp ,' | ',int_temp
			write(*,"(A, ES10.3, A, ES10.3, A, ES10.3)") ' 	      salt= ', min_salt ,' | ',max_salt ,' | ',int_salt
			write(*,"(A, ES10.3, A, ES10.3, A, A     )") ' 	 wvel(1,:)= ', min_wvel ,' | ',max_wvel ,' | ','N.A.'
			write(*,"(A, ES10.3, A, ES10.3, A, A     )") ' 	 wvel(2,:)= ', min_wvel2,' | ',max_wvel2,' | ','N.A.'
			write(*,"(A, ES10.3, A, ES10.3, A, A     )") ' 	 uvel(1,:)= ', min_uvel ,' | ',max_uvel ,' | ','N.A.'
			write(*,"(A, ES10.3, A, ES10.3, A, A     )") ' 	 uvel(2,:)= ', min_uvel2,' | ',max_uvel2,' | ','N.A.'
			write(*,"(A, ES10.3, A, ES10.3, A, A     )") ' 	 vvel(1,:)= ', min_vvel ,' | ',max_vvel ,' | ','N.A.'
			write(*,"(A, ES10.3, A, ES10.3, A, A     )") ' 	 vvel(2,:)= ', min_vvel2,' | ',max_vvel2,' | ','N.A.'
			write(*,"(A, ES10.3, A, ES10.3, A, A     )") ' 	hnode(1,:)= ', min_hnode,' | ',max_hnode,' | ','N.A.'
			write(*,"(A, ES10.3, A, ES10.3, A, A     )") ' 	hnode(2,:)= ', min_hnode2,' | ',max_hnode2,' | ','N.A.'
			write(*,"(A, A     , A, ES10.3, A, A     )") ' 	     cfl_z= ',' N.A.     ',' | ',max_cfl_z  ,' | ','N.A.'
			write(*,"(A, A     , A, ES10.3, A, A     )") ' 	     pgf_x= ',' N.A.     ',' | ',max_pgfx  ,' | ','N.A.'
			write(*,"(A, A     , A, ES10.3, A, A     )") ' 	     pgf_y= ',' N.A.     ',' | ',max_pgfy  ,' | ','N.A.'
			write(*,"(A, A     , A, ES10.3, A, A     )") ' 	        Av= ',' N.A.     ',' | ',max_av    ,' | ','N.A.'
			write(*,"(A, A     , A, ES10.3, A, A     )") ' 	        Kv= ',' N.A.     ',' | ',max_kv    ,' | ','N.A.'
	if (use_ice)    write(*,"(A, A     , A, ES10.3, A, A     )") ' 	     m_ice= ',' N.A.     ',' | ',max_m_ice  ,' | ','N.A.'
			write(*,*)
		endif
	endif ! --> if (mod(istep,logfile_outfreq)==0) then
end subroutine write_step_info
!
!
!===============================================================================
subroutine check_blowup(istep, mesh)
	use g_config, only: logfile_outfreq, which_ALE
	use MOD_MESH
	use o_PARAM
	use g_PARSUP
	use o_ARRAYS
	use i_ARRAYS
	use g_comm_auto
	use io_BLOWUP
	use g_forcing_arrays
	use diagnostics
	use write_step_info_interface
	implicit none
	
	integer           :: n, nz, istep, found_blowup_loc=0, found_blowup=0
	integer 		  :: el, elidx
        type(t_mesh), intent(in), target :: mesh
#include "associate_mesh.h"
	!___________________________________________________________________________
! ! 	if (mod(istep,logfile_outfreq)==0) then
! ! 		if (mype==0) then 
! ! 			write(*,*) '___CHECK FOR BLOW UP___________ --> mstep=',istep
! ! 			write(*,*)
! ! 		endif 
		do n=1, myDim_nod2d
			
			!___________________________________________________________________
			! check ssh
			if ( ((eta_n(n) /= eta_n(n)) .or. &
				eta_n(n)<-50.0 .or. eta_n(n)>50.0 .or. &
				(d_eta(n) /= d_eta(n)) ) ) then
!!PS 				eta_n(n)<-10.0 .or. eta_n(n)>10.0)) then
				found_blowup_loc=1
				write(*,*) '___CHECK FOR BLOW UP___________ --> mstep=',istep
				write(*,*) ' --STOP--> found eta_n become NaN or <-10.0, >10.0'
				write(*,*) 'mype        = ',mype
				write(*,*) 'mstep       = ',istep
				write(*,*) 'node        = ',n
				write(*,*) 'uln, nln    = ',ulevels_nod2D(n), nlevels_nod2D(n)
				write(*,*) 'glon,glat   = ',geo_coord_nod2D(:,n)/rad
				write(*,*)
				write(*,*) 'eta_n(n)    = ',eta_n(n)
				write(*,*) 'd_eta(n)    = ',d_eta(n)
				write(*,*)
				write(*,*) 'zbar_3d_n   = ',zbar_3d_n(:,n)
				write(*,*) 'Z_3d_n      = ',Z_3d_n(:,n)
				write(*,*)
				write(*,*) 'ssh_rhs = ',ssh_rhs(n),', ssh_rhs_old = ',ssh_rhs_old(n)
				write(*,*)
				write(*,*) 'hbar = ',hbar(n),', hbar_old = ',hbar_old(n)
				write(*,*)
				write(*,*) 'wflux = ',water_flux(n)
				write(*,*)
				write(*,*) 'u_wind = ',u_wind(n),', v_wind = ',v_wind(n)
				write(*,*)
				do nz=1,nod_in_elem2D_num(n)
                    write(*,*) 'stress_surf(1:2,',nz,') = ',stress_surf(:,nod_in_elem2D(nz,n))
				end do
				write(*,*)
				write(*,*) 'm_ice = ',m_ice(n),', m_ice_old = ',m_ice_old(n)
				write(*,*) 'a_ice = ',a_ice(n),', a_ice_old = ',a_ice_old(n)
!!PS 				write(*,*) 'thdgr = ',thdgr(n),', thdgr_old = ',thdgr_old(n)
!!PS 				write(*,*) 'thdgrsn = ',thdgrsn(n)
				write(*,*)
!!PS 				if (lcurt_stress_surf) then
!!PS                     write(*,*) 'curl_stress_surf = ',curl_stress_surf(n)
!!PS                     write(*,*)
!!PS 				endif 
!!PS  				do el=1,nod_in_elem2d_num(n)
!!PS  					elidx = nod_in_elem2D(el,n)
!!PS  					write(*,*) ' elem#=',el,', elemidx=',elidx
!!PS  					write(*,*) ' 	 pgf_x =',pgf_x(:,elidx)
!!PS  					write(*,*) ' 	 pgf_y =',pgf_y(:,elidx)
!!PS ! 					write(*,*) ' 	     U =',UV(1,:,elidx)
!!PS ! 					write(*,*) ' 	     V =',UV(2,:,elidx)
!!PS                     write(*,*)
!!PS  				enddo
!!PS 				write(*,*) 'Wvel(1, n)  = ',Wvel(,n)
				write(*,*) 'Wvel(:, n)  = ',Wvel(ulevels_nod2D(n):nlevels_nod2D(n),n)
				write(*,*)
				write(*,*) 'CFL_z(:,n)  = ',CFL_z(ulevels_nod2D(n):nlevels_nod2D(n),n)
				write(*,*)
!!PS 				write(*,*) 'hnode(1, n)  = ',hnode(1,n)
				write(*,*) 'hnode(:, n)  = ',hnode(ulevels_nod2D(n):nlevels_nod2D(n),n)
				write(*,*)
				
            endif
			
			!___________________________________________________________________
			! check surface vertical velocity --> in case of zlevel and zstar 
			! vertical coordinate its indicator if Volume is conserved  for 
			! Wvel(1,n)~maschine preccision
!!PS 			if ( .not. trim(which_ALE)=='linfs' .and. ( Wvel(1, n) /= Wvel(1, n)  .or. abs(Wvel(1,n))>1e-12 )) then
			if ( .not. trim(which_ALE)=='linfs' .and. ( Wvel(1, n) /= Wvel(1, n)  )) then
				found_blowup_loc=1
				write(*,*) '___CHECK FOR BLOW UP___________ --> mstep=',istep
				write(*,*) ' --STOP--> found surface layer vertical velocity becomes NaN or >1e-12'
				write(*,*) 'mype        = ',mype
				write(*,*) 'mstep       = ',istep
				write(*,*) 'node        = ',n
				write(*,*) 'uln, nln    = ',ulevels_nod2D(n), nlevels_nod2D(n)
				write(*,*) 'glon,glat   = ',geo_coord_nod2D(:,n)/rad
				write(*,*)
				write(*,*) 'Wvel(1, n)  = ',Wvel(1,n)
				write(*,*) 'Wvel(:, n)  = ',Wvel(:,n)
				write(*,*)
				write(*,*) 'hnode(1, n) = ',hnode(1,n)
				write(*,*) 'hnode(:, n) = ',hnode(:,n)
				write(*,*) 'hflux       = ',heat_flux(n)
                write(*,*) 'wflux       = ',water_flux(n)
                write(*,*)
                write(*,*) 'eta_n       = ',eta_n(n)
                write(*,*) 'd_eta(n)    = ',d_eta(n)
                write(*,*) 'hbar        = ',hbar(n)
                write(*,*) 'hbar_old    = ',hbar_old(n)
                write(*,*) 'ssh_rhs     = ',ssh_rhs(n)
                write(*,*) 'ssh_rhs_old = ',ssh_rhs_old(n)
                write(*,*)
                write(*,*) 'CFL_z(:,n)  = ',CFL_z(:,n)
                write(*,*)
                
			end if ! --> if ( .not. trim(which_ALE)=='linfs' .and. ...
				
			!___________________________________________________________________
			! check surface layer thinknesss
			if ( .not. trim(which_ALE)=='linfs' .and. ( hnode(1, n) /= hnode(1, n)  .or. hnode(1,n)< 0 )) then
				found_blowup_loc=1
				write(*,*) '___CHECK FOR BLOW UP___________ --> mstep=',istep
				write(*,*) ' --STOP--> found surface layer thickness becomes NaN or <0'
				write(*,*) 'mype        = ',mype
				write(*,*) 'mstep       = ',istep
				write(*,*) 'node        = ',n
				write(*,*)
				write(*,*) 'hnode(1, n)  = ',hnode(1,n)
				write(*,*) 'hnode(:, n)  = ',hnode(:,n)
				write(*,*)
				write(*,*) 'glon,glat   = ',geo_coord_nod2D(:,n)/rad
				write(*,*)
			end if ! --> if ( .not. trim(which_ALE)=='linfs' .and. ...
				
			
			do nz=1,nlevels_nod2D(n)-1
				!_______________________________________________________________
				! check temp
				if ( (tr_arr(nz, n,1) /= tr_arr(nz, n,1)) .or. &
					tr_arr(nz, n,1) < -5.0 .or. tr_arr(nz, n,1)>60) then
					found_blowup_loc=1
					write(*,*) '___CHECK FOR BLOW UP___________ --> mstep=',istep
					write(*,*) ' --STOP--> found temperture becomes NaN or <-5.0, >60'
					write(*,*) 'mype        = ',mype
					write(*,*) 'mstep       = ',istep
					write(*,*) 'node        = ',n
					write(*,*) 'lon,lat     = ',geo_coord_nod2D(:,n)/rad
					write(*,*) 'nz          = ',nz
					write(*,*) 'nzmin, nzmax= ',ulevels_nod2D(n),nlevels_nod2D(n)
					write(*,*) 'x=', geo_coord_nod2D(1,n)/rad, ' ; ', 'y=', geo_coord_nod2D(2,n)/rad
					write(*,*) 'z=', Z_n(nz)
					write(*,*) 'temp(nz, n) = ',tr_arr(nz, n,1)
					write(*,*) 'temp(: , n) = ',tr_arr(:, n,1)
					write(*,*) 'temp_old(nz,n)= ',tr_arr_old(nz, n,1)
					write(*,*) 'temp_old(: ,n)= ',tr_arr_old(:, n,1)
					write(*,*)
					write(*,*) 'hflux       = ',heat_flux(n)
					write(*,*) 'wflux       = ',water_flux(n)
					write(*,*)
					write(*,*) 'eta_n       = ',eta_n(n)
					write(*,*) 'd_eta(n)    = ',d_eta(n)
					write(*,*) 'hbar        = ',hbar(n)
					write(*,*) 'hbar_old    = ',hbar_old(n)
					write(*,*) 'ssh_rhs     = ',ssh_rhs(n)
					write(*,*) 'ssh_rhs_old = ',ssh_rhs_old(n)
					write(*,*)
					write(*,*) 'm_ice       = ',m_ice(n)
					write(*,*) 'm_ice_old   = ',m_ice_old(n)
					write(*,*) 'm_snow      = ',m_snow(n)
					write(*,*) 'm_snow_old  = ',m_snow_old(n)
					write(*,*)
					write(*,*) 'hnode       = ',hnode(:,n)
					write(*,*) 'hnode_new   = ',hnode_new(:,n)
					write(*,*)
					write(*,*) 'Kv          = ',Kv(:,n)
					write(*,*)
					write(*,*) 'W           = ',Wvel(:,n)
					write(*,*)
					write(*,*) 'CFL_z(:,n)  = ',CFL_z(:,n)
					write(*,*)
! 					do el=1,nod_in_elem2d_num(n)
! 						elidx = nod_in_elem2D(el,n)
! 						write(*,*) ' elem#=',el,', elemidx=',elidx
! 						write(*,*) ' 	 helem =',helem(:,elidx)
! 						write(*,*) ' 	     U =',UV(1,:,elidx)
! 						write(*,*) ' 	     V =',UV(2,:,elidx)
! 					enddo
					write(*,*)
					
				endif ! --> if ( (tr_arr(nz, n,1) /= tr_arr(nz, n,1)) .or. & ...
				
				!_______________________________________________________________
				! check salt
				if ( (tr_arr(nz, n,2) /= tr_arr(nz, n,2)) .or.  &
					tr_arr(nz, n,2) < 0 .or. tr_arr(nz, n,2)>50 ) then
					found_blowup_loc=1
					write(*,*) '___CHECK FOR BLOW UP___________ --> mstep=',istep
					write(*,*) ' --STOP--> found salinity becomes NaN or <0, >50'
					write(*,*) 'mype        = ',mype
					write(*,*) 'mstep       = ',istep
					write(*,*) 'node        = ',n
					write(*,*) 'nz          = ',nz
					write(*,*) 'nzmin, nzmax= ',ulevels_nod2D(n),nlevels_nod2D(n)
					write(*,*) 'x=', geo_coord_nod2D(1,n)/rad, ' ; ', 'y=', geo_coord_nod2D(2,n)/rad
					write(*,*) 'z=', Z_n(nz)
					write(*,*) 'salt(nz, n) = ',tr_arr(nz, n,2)
					write(*,*) 'salt(: , n) = ',tr_arr(:, n,2)
					write(*,*)
					write(*,*) 'temp(nz, n) = ',tr_arr(nz, n,1)
					write(*,*) 'temp(: , n) = ',tr_arr(:, n,1)
					write(*,*)
					write(*,*) 'hflux       = ',heat_flux(n)
					write(*,*)
                                        write(*,*) 'wflux       = ',water_flux(n)
					write(*,*) 'eta_n       = ',eta_n(n)
					write(*,*) 'd_eta(n)    = ',d_eta(n)
					write(*,*) 'hbar        = ',hbar(n)
					write(*,*) 'hbar_old    = ',hbar_old(n)
					write(*,*) 'ssh_rhs     = ',ssh_rhs(n)
					write(*,*) 'ssh_rhs_old = ',ssh_rhs_old(n)
					write(*,*)
					write(*,*) 'hnode       = ',hnode(:,n)
					write(*,*) 'hnode_new   = ',hnode_new(:,n)
					write(*,*)
					write(*,*) 'zbar_3d_n   = ',zbar_3d_n(:,n)
					write(*,*) 'Z_3d_n      = ',Z_3d_n(:,n)
					write(*,*)
					write(*,*) 'Kv          = ',Kv(:,n)
					write(*,*)
 					do el=1,nod_in_elem2d_num(n)
 						elidx = nod_in_elem2D(el,n)
 						write(*,*) ' elem#=',el,', elemidx=',elidx
 						write(*,*) ' 	 Av =',Av(:,elidx)
! 						write(*,*) ' 	 helem =',helem(:,elidx)
! 						write(*,*) ' 	     U =',UV(1,:,elidx)
! 						write(*,*) ' 	     V =',UV(2,:,elidx)
 					enddo
					write(*,*) 'Wvel        = ',Wvel(:,n)
					write(*,*)
					write(*,*) 'CFL_z(:,n)  = ',CFL_z(:,n)
					write(*,*)
					write(*,*) 'glon,glat   = ',geo_coord_nod2D(:,n)/rad
					write(*,*)
				endif ! --> if ( (tr_arr(nz, n,2) /= tr_arr(nz, n,2)) .or.  & ...
			end do ! --> do nz=1,nlevels_nod2D(n)-1
		end do ! --> do n=1, myDim_nod2d
! ! 	end if 
		
		!_______________________________________________________________________
		! check globally if one of the cpus hat a blowup situation. if its the
		! case CPU mype==0 needs to write out the stuff. Write out occurs in 
		! moment only over CPU mype==0
		call MPI_AllREDUCE(found_blowup_loc  , found_blowup  , 1, MPI_INTEGER, MPI_MAX, MPI_COMM_FESOM, MPIerr)
		if (found_blowup==1) then
			call write_step_info(istep,1,mesh)
			if (mype==0) then
				call sleep(1)
				write(*,*)
				write(*,*) '                       MODEL BLOW UP !!!'
				write(*,*) '                             ____'
				write(*,*) '                      __,-~~/~    `---.'
				write(*,*) '                    _/_,---(      ,    )'
				write(*,*) '                __ /        <    /   )  \___'
				write(*,*) '- -- ----===;;;`====------------------===;;;===---- -- -'
				write(*,*) '                   \/  ~"~"~"~"~"~\~"~)~"/'
				write(*,*) '                   (_ (   \  (     >    \)'
				write(*,*) '                    \_( _ <         >_>`'
				write(*,*) '                       ~ `-i` ::>|--"'
				write(*,*) '                           I;|.|.|'
				write(*,*) '                          <|i::|i|`'
				write(*,*) '                         (` ^`"`- ")'
				write(*,*) '                  _____.,-#%&$@%#&#~,._____'
				write(*,*)
			end if
			call blowup(istep, mesh)
			if (mype==0) write(*,*) ' --> finished writing blow up file'
			call par_ex
		endif 
end subroutine
