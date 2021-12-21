!==============================================================================
! calculates basically the new iceberg velocity; if melting is enabled, the
! iceberg dimensions are adjusted as well.
!
!   Thomas Rackow, 29.06.2010
!==============================================================================
subroutine iceberg_dyn(mesh, ib, new_u_ib, new_v_ib, u_ib, v_ib, lon,lat, depth_ib, &
                       height_ib, length_ib, width_ib, iceberg_elem, &
		       mass_ib, Ci, Ca, Co, Cda_skin, Cdo_skin, &
		       rho_ice, rho_air, rho_h2o, P_sill, conc_sill, frozen_in, &
		       file1, file2, P_ib, conci_ib, dt_ib, lastsubstep, &
		       f_u_ib_old, f_v_ib_old, l_semiimplicit, &
		       semiimplicit_coeff, AB_coeff, file3, rho_icb)

! kh 19.02.21
 use g_forcing_arrays 	!for u_wind, v_wind or u_wind_ib, v_wind_ib respectively

! kh 19.02.21
 use i_arrays		!for u_ice , v_ice or u_ice_ib, v_ice_ib respectively

! kh 17.03.21 specification of structures used
!use o_arrays, only: coriolis_ib, Tsurf_ib, Ssurf_ib
 use o_arrays, only: coriolis, Tsurf_ib, Ssurf_ib

 use o_param		!for dt
 !use o_mesh
 use iceberg_params,only: l_melt, coriolis_scale !are icebergs allowed to melt?

 USE MOD_MESH
 use g_parsup

 implicit none

 integer, intent(IN) 	:: ib !current iceberg's index
 real,    intent(OUT)	:: new_u_ib, new_v_ib
 real,    intent(IN)    :: u_ib, v_ib
 real,    intent(IN)	:: lon,lat !radiant
 real,    intent(INOUT) :: depth_ib	!inout for case of melting iceberg
 real,    intent(INOUT) :: height_ib	!inout for case of melting iceberg
 real,    intent(INOUT) :: length_ib	!inout for case of melting iceberg
 real,    intent(INOUT) :: width_ib	!inout for case of melting iceberg
 integer, intent(IN) 	:: iceberg_elem !local
 real,	  intent(IN)	:: mass_ib
 real,    intent(IN)	:: Ci, Ca, Co, Cda_skin, Cdo_skin
 real,    intent(IN)	:: rho_ice, rho_air, rho_h2o
 real,    intent(IN) 	:: P_sill, conc_sill
 real,    intent(INOUT) :: frozen_in
 character,  intent(IN)	:: file1*80, file2*80
 real, 	  intent(OUT)	:: P_ib, conci_ib
 real,    intent(IN) 	:: dt_ib
 logical, intent(IN)	:: lastsubstep
 real, 	  intent(INOUT)	:: f_u_ib_old, f_v_ib_old
 logical, intent(IN)	:: l_semiimplicit
 real, 	  intent(IN)	:: semiimplicit_coeff, AB_coeff
 
 real, dimension(3) 	:: uo_dz, vo_dz, hi_ib3, conci_ib3, uo_keel, vo_keel, T_dz,S_dz, T_keel,S_keel
 real 			:: uo_ib, vo_ib, ua_ib, va_ib, ui_ib, vi_ib, hi_ib, uo_skin_ib, vo_skin_ib
 real 			:: Ao, Aa, Ai, Ad, fcoriolis
 real 			:: au_ib, av_ib
 real, dimension(2,2) 	:: SI_matrix
 real, dimension(2)	:: SI_velo
 real 			:: u_ib_tmp, v_ib_tmp, normold, normnew, abs_omib, abs_omib_skin, ocean_drag
 integer 		:: iter_ib
 real          		:: M_b, M_v, M_e, M_bv, sst_ib, sss_ib ! meltrates (basal, lateral, erosion, lateral 'basal'), temp. & salinity
 real			:: T_ave_ib, S_ave_ib, T_keel_ib, S_keel_ib
 character,  intent(IN)	:: file3*80
 real, intent(IN)	:: rho_icb
 
 integer, dimension(3)  :: tmp_arr

type(t_mesh), intent(in) , target :: mesh
#include "associate_mesh.h"
 
 !OCEAN VELOCITIES: 
 ! - (uo_ib, vo_ib)		: integrated mean velocity at location of iceberg
 ! - (uo_skin_ib, vo_skin_ib)	: velocity below the draft of the iceberg
 ! call iceberg_avvelo_ufkeel(uo_dz,vo_dz, uo_keel,vo_keel, depth_ib,iceberg_elem)
 call iceberg_average_andkeel(mesh, uo_dz,vo_dz, uo_keel,vo_keel, T_dz,S_dz, T_keel,S_keel, depth_ib,iceberg_elem, ib)
 call FEM_3eval(mesh, uo_ib,vo_ib,lon,lat,uo_dz,vo_dz,iceberg_elem)
 call FEM_3eval(mesh, uo_skin_ib,vo_skin_ib,lon,lat,uo_keel,vo_keel,iceberg_elem)
 

 !TEMPERATURE AND SALINITY:
 ! - T_ave_ib, S_ave_ib		: Mean T & S (integrated) at location of iceberg
 ! - T_keel_ib, S_keel_ib	: T & S below the draft of the iceberg (depth_ib)
 call FEM_3eval(mesh, T_ave_ib,S_ave_ib,lon,lat,T_dz,S_dz,iceberg_elem)
 !write(*,*) "LA DEBUG: T_ave_ib = ", T_ave_ib, ", S_ave_ib = ", S_ave_ib, ", T_dz = ", T_dz, ", S_dz = ", S_dz, ", iceberg_elem = ", iceberg_elem
 call FEM_3eval(mesh, T_keel_ib,S_keel_ib,lon,lat,T_keel,S_keel,iceberg_elem)
 !write(*,*) "LA DEBUG: T_keel_ib = ", T_keel_ib, ", S_keel_ib = ", S_keel_ib, ", T_keel = ", T_keel, ", S_keel = ", S_keel, ", iceberg_elem = ", iceberg_elem


 !ATMOSPHERIC VELOCITY ua_ib, va_ib
 ! kh 11.02.21
 call FEM_eval(mesh, ua_ib,va_ib,lon,lat,u_wind_ib,v_wind_ib,iceberg_elem)      

 !ICE VELOCITY ui_ib, vi_ib
 ! kh 11.02.21
 call FEM_eval(mesh, ui_ib,vi_ib,lon,lat,u_ice_ib,v_ice_ib,iceberg_elem)    
  
 !ICE THICKNESS (CONCENTRATION) hi_ib, conci_ib
 
 ! LA debug
 tmp_arr=elem2D_nodes(:,iceberg_elem)

 hi_ib3    = m_ice_ib(tmp_arr)
 conci_ib3 = a_ice_ib(tmp_arr) 
 call FEM_3eval(mesh, hi_ib,conci_ib,lon,lat,hi_ib3,conci_ib3,iceberg_elem)
 P_ib = 20000. * hi_ib * exp(-20.*(1-conci_ib))
 
 call compute_areas(Ao, Aa, Ai, Ad, depth_ib, &
                    height_ib, length_ib, width_ib, hi_ib)


 !fcoriolis = coriolis_param_elem2D(iceberg_elem) * coriolis_scale(ib)

! kh 16.03.21 use coriolis_ib buffered values here, test only
!fcoriolis = coriolis_ib(iceberg_elem) * coriolis_scale(ib)
 fcoriolis = coriolis(iceberg_elem) * coriolis_scale(ib)

 call iceberg_acceleration( mesh, ib, au_ib, av_ib, Ao, Aa, Ai, Ad, 		&
                            uo_ib,vo_ib, uo_skin_ib, vo_skin_ib,	&
			    ua_ib,va_ib, ui_ib,vi_ib, 			&
		            u_ib, v_ib, mass_ib, fcoriolis, 		&
			    Ci, Ca, Co, Cda_skin, Cdo_skin, 		&
			    rho_ice, rho_air, rho_h2o, length_ib, 	&
			    iceberg_elem, conci_ib, file1, file2, 	&
			    lon, lat, lastsubstep, f_u_ib_old, 		&
			    f_v_ib_old, l_semiimplicit, 		&
			    semiimplicit_coeff, AB_coeff )


 !========================THERMODYNAMICS============================
 if(l_melt) then

! kh 15.03.21 use Tsurf_ib and Ssurf_ib buffered values here
  call FEM_eval(mesh, sst_ib,sss_ib,lon,lat,Tsurf_ib,Ssurf_ib,iceberg_elem)

  !if(ib==3497) then
  ! write(*,*) '#####################################################'
  ! write(*,*) 'SST:',sst_ib,'T_ave:',T_ave_ib,'T_keel:',T_keel_ib
  ! write(*,*) 'SSS:',sss_ib,'S_ave:',S_ave_ib,'S_keel:',S_keel_ib,'depth:',depth_ib
  ! write(*,*) 'uo_keel:',uo_skin_ib,'vo_keel:',vo_skin_ib
  ! write(*,*) 'uo_ave:',uo_ib,'vo_ave:',vo_ib
  !end if
 
  !call iceberg_meltrates(M_b, M_v, M_e, M_bv, u_ib,v_ib, uo_ib,vo_ib, ua_ib,va_ib, &
  !			  sst_ib, length_ib, conci_ib)
  call iceberg_meltrates(	M_b, M_v, M_e, M_bv, &
				u_ib,v_ib, uo_ib,vo_ib, ua_ib,va_ib, &
				sst_ib, length_ib, conci_ib, &
				uo_skin_ib, vo_skin_ib, T_keel_ib, S_keel_ib, depth_ib, &
				T_ave_ib, S_ave_ib, ib)

  call iceberg_newdimensions(ib, depth_ib,height_ib,length_ib,width_ib,M_b,M_v,M_e,M_bv, &
			     rho_h2o, rho_icb, file3)

 end if 
 !====================END OF THERMODYNAMICS=========================

  new_u_ib = u_ib + au_ib * dt_ib
  new_v_ib = v_ib + av_ib * dt_ib

 if (l_semiimplicit) then !a matrix multiplication is to be performed
  			  !for semiimpl. coriolis term and implicit
			  !water drag
			   
   abs_omib = sqrt( (uo_ib - u_ib)**2 + (vo_ib - v_ib)**2 )
   abs_omib_skin = sqrt( (uo_skin_ib - u_ib)**2 + (vo_skin_ib - v_ib)**2 )
   
   ocean_drag = (0.5 * Co * rho_h2o * Ao * abs_omib + rho_h2o * Cdo_skin * Ad  &
  	           * abs_omib_skin)/mass_ib
   
  
   SI_matrix(1,1) = 1. + dt_ib*ocean_drag
   SI_matrix(1,2) = dt_ib*fcoriolis*semiimplicit_coeff
   SI_matrix(2,1) =-SI_matrix(1,2)
   SI_matrix(2,2) = SI_matrix(1,1)
   SI_matrix = (1./( SI_matrix(2,2)**2 + SI_matrix(1,2)**2 )) * SI_matrix
   !new velocity
   SI_velo = MATMUL(SI_matrix, (/ new_u_ib, new_v_ib /))
      
   !for P_sill calculate "effective acceleration" in case of semi-impl. scheme by
   !(new velo - old velo)/dt for au_ib (av_ib) isn't the real acceleration
   !au_ib = (SI_velo(1) - u_ib)/dt_ib
   !av_ib = (SI_velo(2) - v_ib)/dt_ib
      
   !now the velocity can be updated
   new_u_ib = SI_velo(1)
   new_v_ib = SI_velo(2)
   
 else !compute only water drag implicitly, coriolis: AB
  
   abs_omib = sqrt( (uo_ib - u_ib)**2 + (vo_ib - v_ib)**2 )
   abs_omib_skin = sqrt( (uo_skin_ib - u_ib)**2 + (vo_skin_ib - v_ib)**2 )
   ocean_drag = (0.5 * Co * rho_h2o * Ao * abs_omib + rho_h2o * Cdo_skin * Ad  &
  	           * abs_omib_skin)/mass_ib  
   		   
   SI_matrix(1,1) = 1. + dt_ib*ocean_drag
   SI_matrix(1,2) = 0.0
   SI_matrix(2,1) = 0.0
   SI_matrix(2,2) = SI_matrix(1,1)
   SI_matrix = (1./(SI_matrix(2,2)**2)) * SI_matrix
   !new velocity
   SI_velo = MATMUL(SI_matrix, (/ new_u_ib, new_v_ib /))
      
   !now the velocity can be updated
   new_u_ib = SI_velo(1)
   new_v_ib = SI_velo(2)  
    
 end if !matrix-multiplication


 !icebergs may be frozen in: .true. or .false.
 call iceberg_frozen(frozen_in, P_sill, P_ib, conc_sill, conci_ib, ib)
      
      
 new_u_ib = (1-frozen_in) * new_u_ib + frozen_in * ui_ib
 new_v_ib = (1-frozen_in) * new_v_ib + frozen_in * vi_ib
 
end subroutine iceberg_dyn


 !***************************************************************************************************************************
 !***************************************************************************************************************************


subroutine iceberg_frozen(festgefroren, P_sill, P_ib, conc_sill, conci_ib, ib)
 
 use iceberg_params, only : l_freeze !use 'capturing mechanism' of sea ice?
 implicit none
 
 real,	  intent(OUT)	:: festgefroren
 real,	  intent(IN)	:: P_sill, P_ib
 real,    intent(IN)	:: conc_sill, conci_ib
 integer, intent(IN)	:: ib
 
 festgefroren=0.0
 
 if( l_freeze(ib) ) then
  festgefroren= factor(P_ib,P_sill,P_sill-3000.) &
              * factor(conci_ib,conc_sill,conc_sill-0.04)
 end if  
 
 contains
 
  !=================================================================
  ! ... a function for linear transition zone between "frozen in" 
  ! and "not frozen in", where sill & zone are defined as below
  !
  !	      sill
  !              |
  !	        /--------- 1
  ! 0  ________/
  !           |
  !         zone			and value is P_ib or conci_ib.
  !
  !=================================================================
 real function factor(value, sill,zone)
  implicit none
  real, intent(IN) :: value
  real, intent(IN) :: sill
  real, intent(IN) :: zone

  if(value < zone) then
   factor=0.
  else if(value > sill) then
      	factor=1.
       else
      	factor=(value-zone)/(sill-zone)
  end if	
 end function factor
 
end subroutine iceberg_frozen


 !***************************************************************************************************************************
 !***************************************************************************************************************************


subroutine iceberg_acceleration(mesh, ib, au_ib, av_ib, Ao, Aa, Ai, Ad, &
                                uo_ib,vo_ib, uo_skin_ib, vo_skin_ib,  &
				ua_ib,va_ib, ui_ib,vi_ib, &
				u_ib, v_ib, mass_ib, fcoriolis, &
				Ci, Ca, Co, Cda_skin, Cdo_skin, &
				rho_ice, rho_air, rho_h2o, length_ib, &
				iceberg_elem, conci_ib, file1, file2, &
				lon_rad, lat_rad, output, f_u_ib_old, &
				f_v_ib_old, l_semiimplicit, &
				semiimplicit_coeff, AB_coeff )

! kh 17.03.21 not really used here
! use o_arrays 	!for ssh


 use o_param 	!for g
 use g_config	!for istep
 use g_parsup	!for mype
 use g_rotate_grid,  only: vector_r2g, vector_g2r
 use iceberg_params, only: l_wave, l_tides, l_geo_out, surfslop_scale, ascii_out

use MOD_MESH

 implicit none
 
 integer,   intent(IN)	:: ib
 real, 	    intent(OUT) :: au_ib, av_ib
 real, 	    intent(IN)  :: Ao, Aa, Ai, Ad 
 real, 	    intent(IN)  :: uo_ib,vo_ib, uo_skin_ib, vo_skin_ib, ua_ib,va_ib, ui_ib,vi_ib, u_ib,v_ib
 real, 	    intent(IN)  :: mass_ib, fcoriolis
 real, 	    intent(IN)  :: Ci, Ca, Co, Cda_skin, Cdo_skin
 real, 	    intent(IN)  :: rho_ice, rho_air, rho_h2o, length_ib
 integer,   intent(IN)	:: iceberg_elem
 real, 	    intent(IN)  :: conci_ib
 character, intent(IN)	:: file1*80, file2*80
 real, 	    intent(IN)  :: lon_rad, lat_rad 
 logical,   intent(IN) 	:: output
 real, 	   intent(INOUT):: f_u_ib_old, f_v_ib_old
 logical,   intent(IN) 	:: l_semiimplicit
 real, 	    intent(IN) 	:: semiimplicit_coeff, AB_coeff
 
 real 			:: vel_atm, wave_amplitude, direction_u, direction_v
 real 			:: abs_omib, abs_amib, abs_imib, abs_omib_skin
 real 			:: ocean_drag_u, ocean_skin_u, air_drag_u, air_skin_u 
 real 			:: ice_drag_u, wave_radiation_u, surface_slope_u, slope_tides_u
 real 			:: ocean_drag_v, ocean_skin_v, air_drag_v, air_skin_v 
 real 			:: ice_drag_v, wave_radiation_v, surface_slope_v, slope_tides_v
 real, dimension(2) 	:: nablaeta
 real, dimension(8)	:: accs_u_out, accs_v_out
 real, dimension(4)	:: vels_u_out, vels_v_out
 real 			:: oneminus_AB !, test1, test2
 integer 		:: icbID, i, istep
 
type(t_mesh), intent(in) , target :: mesh
#include "associate_mesh.h"

icbID = mype+10
 
 !estimate wave height at the icebergs location (Bigg et al., 1997),
 !so wave_amplitude = 0.5 * wave_height = 0.5 * const. * abs(atm velo)**2
 vel_atm = sqrt(ua_ib**2 + va_ib**2)
 wave_amplitude = 0.5 * 0.02025 * vel_atm**2
 
 !assume that waves have same direction as the winds
 direction_u = ua_ib / vel_atm
 direction_v = va_ib / vel_atm
  
 !absolute values of relative velocities 
 abs_omib = sqrt( (uo_ib - u_ib)**2 + (vo_ib - v_ib)**2 )
 abs_amib = sqrt( (ua_ib - u_ib)**2 + (va_ib - v_ib)**2 )
 abs_imib = sqrt( (ui_ib - u_ib)**2 + (vi_ib - v_ib)**2 )  
 abs_omib_skin = sqrt( (uo_skin_ib - u_ib)**2 + (vo_skin_ib - v_ib)**2 )
 
 ! u-components
 ocean_drag_u     = (0.5 * Co * rho_h2o * Ao  &
                   * abs_omib * uo_ib)/mass_ib	!calculate part of it implicitly

 
 !ocean_skin_u     = (rho_h2o * Cdo_skin * Ad  &
 ! 	           * abs_omib * uo_ib)/mass_ib !calculate part of it implicitly
 ocean_skin_u     = (rho_h2o * Cdo_skin * Ad  &
  	           * abs_omib_skin * uo_skin_ib)/mass_ib !calculate part of it implicitly
		       
 air_drag_u       = (0.5 * Ca * rho_air * Aa  &
	           * abs_amib * (ua_ib - u_ib))/mass_ib
 air_skin_u       = (rho_air * Cda_skin * Ad   &
  	           * abs_amib * (ua_ib - u_ib))/mass_ib 
	     
 ice_drag_u       = (0.5 * Ci * rho_ice * Ai  &
	           * abs_imib * (ui_ib - u_ib))/mass_ib
		   
 if(l_wave(ib)) then
   wave_radiation_u = 1./4. * rho_h2o * g * length_ib &
                     * wave_amplitude**2 * direction_u /mass_ib
 else
   wave_radiation_u = 0.0
 end if
 
 !use gradient smoothing for surface slope term
 call mean_gradient(mesh, iceberg_elem, lon_rad, lat_rad, nablaeta)
 
 !additional surface slope due to tides
 if(l_tides) then
   slope_tides_u = 0.0  !-g * sum( opbnd_z_tide(elem2D_nodes(:,iceberg_elem)) * bafux_2D(:,iceberg_elem) )
 else
   slope_tides_u = 0.0
 end if
   
 !surface_slope_u = -g * sum( ssh(elem2D_nodes(:,iceberg_elem)) & 
 !		   * bafux_2D(:,iceberg_elem) ) !-g* nabla ssh
 surface_slope_u = (-g * nablaeta(1) + slope_tides_u) * surfslop_scale(ib) !default scaling is 1.0
 
 
 ! v-components	       
 !ocean_drag_v     = (0.5 * Co * rho_h2o * Ao  &
 !                  * abs_omib * (vo_ib - v_ib))/mass_ib
 ocean_drag_v     = (0.5 * Co * rho_h2o * Ao  &
                   * abs_omib * vo_ib)/mass_ib	!calculate part of it implicitly
 
 !ocean_skin_v     = (rho_h2o * Cdo_skin * Ad  &
 ! 	           * abs_omib * vo_ib)/mass_ib !calculate part of it implicitly
 ocean_skin_v     = (rho_h2o * Cdo_skin * Ad  &
  	           * abs_omib_skin * vo_skin_ib)/mass_ib !calculate part of it implicitly
 
	     
 air_drag_v       = (0.5 * Ca * rho_air * Aa  &
	           * abs_amib * (va_ib - v_ib))/mass_ib
 air_skin_v       = (rho_air * Cda_skin * Ad   &
  	           * abs_amib * (va_ib - v_ib))/mass_ib 
	     
 ice_drag_v       = (0.5 * Ci * rho_ice * Ai  &
	           * abs_imib * (vi_ib - v_ib))/mass_ib
 
 if(l_wave(ib)) then	       
   wave_radiation_v = 1./4. * rho_h2o * g * length_ib &
                     * wave_amplitude**2 * direction_v /mass_ib	
 else
   wave_radiation_v = 0.0
 end if
 
 !additional surface slope due to tides
 if(l_tides) then
   slope_tides_v = 0.0  !-g * sum( opbnd_z_tide(elem2D_nodes(:,iceberg_elem)) * bafuy_2D(:,iceberg_elem) )	
 else
   slope_tides_v = 0.0
 end if
    	   
 !surface_slope_v = -g * sum( ssh(elem2D_nodes(:,iceberg_elem)) &
 !		   * bafuy_2D(:,iceberg_elem) )
 surface_slope_v = (-g * nablaeta(2) + slope_tides_v) * surfslop_scale(ib) !default scaling is 1.0
 
 
 if (l_semiimplicit) then !USE (SEMI-)IMPLICIT SCHEME for coriolis term
 
  if (conci_ib .GT. 0.15) then

      au_ib = ocean_drag_u     &
	    + ocean_skin_u     &
	    + air_drag_u       &
	    + air_skin_u       &
	    + ice_drag_u       &
	    + wave_radiation_u &
	    + surface_slope_u  &
	    + (1.-semiimplicit_coeff)*fcoriolis*v_ib

      av_ib = ocean_drag_v     &
	    + ocean_skin_v     &
	    + air_drag_v       &
	    + air_skin_v       &
	    + ice_drag_v       &
	    + wave_radiation_v &
	    + surface_slope_v  &	   
	    - (1.-semiimplicit_coeff)*fcoriolis*u_ib
  else

      au_ib = ocean_drag_u     &
	    + ocean_skin_u     &
	    + air_drag_u       &
	    + air_skin_u       &
	    + wave_radiation_u &
	    + surface_slope_u  &
	    + (1.-semiimplicit_coeff)*fcoriolis*v_ib

      av_ib = ocean_drag_v     &
	    + ocean_skin_v     &
	    + air_drag_v       &
	    + air_skin_v       &  
	    + wave_radiation_v &
	    + surface_slope_v  &
	    - (1.-semiimplicit_coeff)*fcoriolis*u_ib
  end if
  
 else !USE ADAMS-BASHFORTH SCHEME for coriolis
  oneminus_AB= 1. - AB_coeff
 
  if (conci_ib .GT. 0.15) then
    
      au_ib = ocean_drag_u     &
	    + ocean_skin_u     &
	    + air_drag_u       &
	    + air_skin_u       &
	    + ice_drag_u       &
	    + wave_radiation_u &
	    + surface_slope_u  &
	    + AB_coeff*fcoriolis*v_ib+oneminus_AB*f_v_ib_old


      av_ib = ocean_drag_v     &
	    + ocean_skin_v     &
	    + air_drag_v       &
	    + air_skin_v       &
	    + ice_drag_v       &	
	    + wave_radiation_v & 
	    + surface_slope_v  &    
	    - (AB_coeff*fcoriolis*u_ib) - (oneminus_AB*f_u_ib_old)
  else
    
      au_ib = ocean_drag_u     &
	    + ocean_skin_u     &
	    + air_drag_u       &
	    + air_skin_u       &
	    + wave_radiation_u &
	    + surface_slope_u  &
	    + AB_coeff*fcoriolis*v_ib+oneminus_AB*f_v_ib_old

      av_ib = ocean_drag_v     &
	    + ocean_skin_v     &
	    + air_drag_v       &
	    + air_skin_v       &
	    + wave_radiation_v &
	    + surface_slope_v  &
	    - (AB_coeff*fcoriolis*u_ib) - (oneminus_AB*f_u_ib_old)
  end if
    
  !save f*velocity for A.B. scheme before it is updated
  f_u_ib_old=fcoriolis*u_ib
  f_v_ib_old=fcoriolis*v_ib
  
 end if !use semiimplicit scheme or AB method?
 
 !if(.false.) then
 if(output .AND. ascii_out) then
 
   accs_u_out(1) = ocean_drag_u - (0.5 * Co * rho_h2o * Ao * abs_omib * u_ib)/mass_ib
   accs_u_out(2) = ocean_skin_u - (rho_h2o * Cdo_skin * Ad * abs_omib_skin * u_ib)/mass_ib
   accs_u_out(3) = air_drag_u
   accs_u_out(4) = air_skin_u
   accs_u_out(5) = ice_drag_u
   accs_u_out(6) = wave_radiation_u
   accs_u_out(7) = surface_slope_u
   accs_u_out(8) = fcoriolis*v_ib
   vels_u_out(1) = uo_ib
   vels_u_out(2) = uo_skin_ib
   vels_u_out(3) = ua_ib
   vels_u_out(4) = ui_ib
  
   accs_v_out(1) = ocean_drag_v - (0.5 * Co * rho_h2o * Ao * abs_omib * v_ib)/mass_ib
   accs_v_out(2) = ocean_skin_v - (rho_h2o * Cdo_skin * Ad * abs_omib_skin * v_ib)/mass_ib
   accs_v_out(3) = air_drag_v
   accs_v_out(4) = air_skin_v
   accs_v_out(5) = ice_drag_v
   accs_v_out(6) = wave_radiation_v
   accs_v_out(7) = surface_slope_v
   accs_v_out(8) = -fcoriolis*u_ib
   vels_v_out(1) = vo_ib
   vels_v_out(2) = vo_skin_ib
   vels_v_out(3) = va_ib
   vels_v_out(4) = vi_ib
   
   if(l_geo_out) then
    do i=1,8
     call vector_r2g(accs_u_out(i), accs_v_out(i), lon_rad, lat_rad, 0)
    end do
    do i=1,4
     call vector_r2g(vels_u_out(i), vels_v_out(i), lon_rad, lat_rad, 0)
    end do   
   end if
   
   !test1 = -0.1
   !test2 = 0.1
   !write(*,*) 'Test rotation:', ib, test1, test2 
   !call vector_r2g(test1, test2, lon_rad, lat_rad, 0)
   !write(*,*) 'After vector_r2g:', ib, test1, test2 
   !call vector_g2r(test1, test2, lon_rad, lat_rad, 0)
   !write(*,*) 'After vector_g2r (should be the same):', ib, test1, test2 
   
   open(unit=icbID,file=file1,position='append')
   write(icbID,'(2I,12e15.7)') mype, istep, 	&
  		 accs_u_out(1),	accs_u_out(2),	&
		 accs_u_out(3),	accs_u_out(4),	&
		 accs_u_out(5),	accs_u_out(6),	&	
		 accs_u_out(7),	accs_u_out(8), 	&
		 vels_u_out(1),	vels_u_out(2),	&
		 vels_u_out(3),	vels_u_out(4)
		 
		 
		 
   close(icbID)  
 
   open(unit=icbID,file=file2,position='append')      
   write(icbID,'(2I,12e15.7)') mype,istep, 	&
  		 accs_v_out(1),	accs_v_out(2),	&
		 accs_v_out(3),	accs_v_out(4),	&
		 accs_v_out(5),	accs_v_out(6),	&	
		 accs_v_out(7),	accs_v_out(8),	&		 
		 vels_v_out(1),	vels_v_out(2),	&
		 vels_v_out(3),	vels_v_out(4)
   close(icbID)
 end if !output
  
end subroutine iceberg_acceleration


!***************************************************************************************************************************
!***************************************************************************************************************************


subroutine compute_areas(Ao, Aa, Ai, Ad, depth_ib, &
                         height_ib, length_ib, width_ib, hi_ib)
 implicit none
 
 real, intent(OUT) :: Ao, Aa, Ai, Ad
 real, intent(IN)  :: depth_ib, height_ib, length_ib, width_ib, hi_ib
  
  !area of iceberg exposed to ocean, atm, seaice, horizontal area
  Ao = abs(depth_ib) * length_ib
  Aa = (height_ib - abs(depth_ib)) * length_ib
  Ai = hi_ib * length_ib
  Ad = length_ib * width_ib
  
end subroutine compute_areas

 !***************************************************************************************************************************
 !***************************************************************************************************************************


subroutine iceberg_average_andkeel(mesh, uo_dz,vo_dz, uo_keel,vo_keel, T_dz,S_dz, T_keel,S_keel, depth_ib,iceberg_elem, ib)
  !use o_mesh
  USE MOD_MESH
  use o_param
  use i_therm_param
  use i_param
  use i_arrays
  use g_parsup
  
! kh 17.03.21 specification of structures used
  use o_arrays, only: UV_ib, tr_arr_ib, Z_3d_n_ib

  use g_clock
  use g_forcing_arrays
  use g_rotate_grid
  
  implicit none

  REAL, DIMENSION(3), INTENT(OUT) :: uo_dz
  REAL, DIMENSION(3), INTENT(OUT) :: vo_dz
  REAL, DIMENSION(3), INTENT(OUT) :: uo_keel
  REAL, DIMENSION(3), INTENT(OUT) :: vo_keel
  REAL, DIMENSION(3), INTENT(OUT) :: T_dz
  REAL, DIMENSION(3), INTENT(OUT) :: S_dz
  REAL, DIMENSION(3), INTENT(OUT) :: T_keel
  REAL, DIMENSION(3), INTENT(OUT) :: S_keel
  REAl,               INTENT(IN)  :: depth_ib
  INTEGER,            INTENT(IN)  :: iceberg_elem, ib

  real           :: lev_up, lev_low
  integer        :: m, k, n2, n_up, n_low, cavity_count
  ! depth over which is integrated (layer and sum)
  real           :: dz, ufkeel1, ufkeel2, Temkeel, Salkeel

type(t_mesh), intent(in) , target :: mesh
#include "associate_mesh.h"

  cavity_count=0

  !LOOP: over all nodes of the iceberg element
  nodeloop: do m=1, 3
   !for each 2D node of the iceberg element..
   n2=elem2D_nodes(m,iceberg_elem)

   uo_dz(m)=0.0
   vo_dz(m)=0.0
   uo_keel(m)=0.0
   vo_keel(m)=0.0
   T_dz(m)=0.0
   S_dz(m)=0.0
   T_keel(m)=0.0
   S_keel(m)=0.0

   ! LOOP: consider all neighboring pairs (n_up,n_low) of 3D nodes
   ! below n2..
   innerloop: do k=1, nl+1


    !##################################################
    ! *** LA changed to zbar_3d_n ***

! kh 18.03.21 use zbar_3d_n_ib buffered values here
    if( k==1 ) then
        lev_up = 0.0
    else
        lev_up  = Z_3d_n_ib(k-1, n2)
    end if
    lev_low = Z_3d_n_ib(k, n2)
    dz = abs( lev_low - lev_up )
    
    !if( abs(lev_up)>=abs(depth_ib) ) then
    !    ! ...icb bottom above lev_up --> no further integration
    !end if
    
    !if( (abs(coord_nod3D(3, n_low))>abs(depth_ib)) .AND. (abs(coord_nod3D(3, n_up))>abs(depth_ib)) ) then
    ! write(*,*) 'INFO, k:',k,'z_up:',coord_nod3D(3, n_up),'z_lo:',coord_nod3D(3, n_low),'depth:',depth_ib,'cavity:',(cavity_flag_nod2d(elem2D_nodes(m,iceberg_elem))==1)
    !end if

#ifdef use_cavity
    ! if cavity node ..
    if( cavity_flag_nod2d(elem2D_nodes(m,iceberg_elem))==1 .AND. abs(depth_ib)<abs(lev_up) ) then
    ! LA: Never go here for k=1, because abs(depth_ib)>=0.0 for all icebergs

      !cavity_count=cavity_count + 1
      !if(cavity_count==3) then
      !  write(*,*) '3 cavity nodes in the same element! Iceberg should not be here!'
      !end if

! kh 08.03.21 use UV_ib buffered values here
      uo_dz(m)=UV_ib(1,k-1,n2)*abs(depth_ib)
      vo_dz(m)=UV_ib(2,k-1,n2)*abs(depth_ib)
      uo_keel(m)=UV_ib(1,k-1,n2)
      vo_keel(m)=UV_ib(2,k-1,n2)

      ! *** LA changed to tr_arr(nl-1,node_size,num_tr) ***
! kh 15.03.21 use tr_arr_ib buffered values here
      T_dz(m)=tr_arr_ib(k-1,n2,1)*abs(depth_ib)
      S_dz(m)=tr_arr_ib(k-1,n2,2)*abs(depth_ib)
      T_keel(m)=tr_arr_ib(k-1,n2,1)
      S_keel(m)=tr_arr_ib(k-1,n2,2) ! check those choices with RT: OK

      exit innerloop

    ! if the lowest z coord is below the iceberg draft, exit
    !else if( abs(coord_nod3D(3, n_low))>=abs(depth_ib) .AND. abs(coord_nod3D(3, n_up))<=abs(depth_ib) ) then
    
    !****************************************************************
    ! LA 23.11.21 case if depth_ib<lev_up
    else if( abs(lev_low)>=abs(depth_ib) ) then !.AND. (abs(lev_up)<=abs(depth_ib)) ) then
#else
    !if( abs(coord_nod3D(3, n_low))>=abs(depth_ib) .AND. abs(coord_nod3D(3, n_up))<=abs(depth_ib) ) then
    if( abs(lev_low)>=abs(depth_ib) ) then !.AND. (abs(lev_up)<=abs(depth_ib)) ) then
#endif
      if( abs(lev_up)<abs(depth_ib) ) then
        dz = abs ( lev_up - depth_ib )
      else
        ! LA: Never go here, when starting with k=1
        dz = abs(depth_ib)
      end if
    !****************************************************************

      ! *** LA changed to tr_arr(nl-1,node_size,num_tr) ***
! kh 08.03.21 use UV_ib buffered values here
      if( k==1 ) then
        ufkeel1 = UV_ib(1,k,n2)
        ufkeel2 = UV_ib(2,k,n2)
        Temkeel = tr_arr_ib(k,n2,1)
        Salkeel = tr_arr_ib(k,n2,2)
        uo_dz(m)=ufkeel1*dz 
        vo_dz(m)=ufkeel2*dz
        T_dz(m)=Temkeel*dz
        S_dz(m)=Salkeel*dz
      else
        ufkeel1 = interpol1D(abs(lev_up),UV_ib(1,k-1,n2),abs(lev_low),UV_ib(1,k,n2),abs(depth_ib))
        ufkeel2 = interpol1D(abs(lev_up),UV_ib(2,k-1,n2),abs(lev_low),UV_ib(2,k,n2),abs(depth_ib))
! kh 1  5.03.21 use tr_arr_ib buffered values here
        Temkeel = interpol1D(abs(lev_up),tr_arr_ib(k-1,n2,1),abs(lev_low),tr_arr_ib(k,n2,1),abs(depth_ib))
        Salkeel = interpol1D(abs(lev_up),tr_arr_ib(k-1,n2,2),abs(lev_low),tr_arr_ib(k,n2,2),abs(depth_ib))
! kh 08.03.21 use UV_ib buffered values here
        uo_dz(m)=uo_dz(m)+ 0.5*(UV_ib(1,k-1,n2) + ufkeel1)*dz 
        vo_dz(m)=vo_dz(m)+ 0.5*(UV_ib(2,k-1,n2) + ufkeel2)*dz
! kh 15.03.21 use tr_arr_ib buffered values here
        T_dz(m)=T_dz(m)+ 0.5*(tr_arr_ib(k-1,n2,1)+ Temkeel)*dz
        S_dz(m)=S_dz(m)+ 0.5*(tr_arr_ib(k-1,n2,2)+ Salkeel)*dz
      end if

      uo_keel(m)=ufkeel1
      vo_keel(m)=ufkeel2

      T_keel(m)=Temkeel
      S_keel(m)=Salkeel
      
      if(S_dz(m)/abs(depth_ib)>70.) then
       write(*,*) 'innerloop, dz:',dz,', depth:',depth_ib,',S_dz(m):',S_dz(m),"m:",m,", k:",k,", tr_arr_ib(k-1,n2,1):",tr_arr_ib(k-1,n2,1),", tr_arr_ib(k,n2,1):", tr_arr_ib(k,n2,1),", Salkeel:",Salkeel,", lev_low:",lev_low,", lev_up:",lev_up
      end if
      
      if(T_dz(m)/abs(depth_ib)>70.) then
       write(*,*) 'innerloop, dz:',dz,', depth:',depth_ib,',T_dz(m):',T_dz(m),"m:",m,", k:",k,", tr_arr_ib(k-1,n2,2):",tr_arr_ib(k-1,n2,2),", tr_arr_ib(k,n2,2):", tr_arr_ib(k,n2,2),",Temkeel:",Temkeel,", lev_low:",lev_low,", lev_up:",lev_up
      end if

      exit innerloop
    
    !****************************************************************
    ! LA 23.11.21 case if lev_low==0
    else if(lev_low==lev_up) then
      exit innerloop
    !****************************************************************
    
    else	
      if( k==1 ) then
          !write(*,*) "LA DEBUG: lev_low = ", lev_low, ", lev_up = ", lev_up, ", depth_ib = ", depth_ib
          cycle
      end if


      ! .. and sum up the layer-integrated velocities ..
! kh 08.03.21 use UV_ib buffered values here      
      uo_dz(m)=uo_dz(m)+ 0.5*(UV_ib(1,k-1,n2)+UV_ib(1,k,n2))*dz
      vo_dz(m)=vo_dz(m)+ 0.5*(UV_ib(2,k-1,n2)+UV_ib(2,k,n2))*dz
      ! .. and T & S:

! kh 15.03.21 use tr_arr_ib buffered values here
      T_dz(m)=T_dz(m)+ 0.5*(tr_arr_ib(k-1,n2,1)+ tr_arr_ib(k,n2,1))*dz
      S_dz(m)=S_dz(m)+ 0.5*(tr_arr_ib(k-1,n2,2)+ tr_arr_ib(k,n2,2))*dz

   
      if(S_dz(m)/abs(depth_ib)>70.) then
       write(*,*) 'innerloop, dz:',dz,', depth:',depth_ib,',S_dz(m):',S_dz(m),"m:",m,", k:",k,", tr_arr_ib(k-1,n2,1):",tr_arr_ib(k-1,n2,1),", tr_arr_ib(k,n2,1):", tr_arr_ib(k,n2,1),", lev_low:",lev_low,", lev_up:",lev_up
      end if
      
      if(T_dz(m)/abs(depth_ib)>70.) then
       write(*,*) 'innerloop, dz:',dz,', depth:',depth_ib,',T_dz(m):',T_dz(m),"m:",m,", k:",k,", tr_arr_ib(k-1,n2,2):",tr_arr_ib(k-1,n2,2),", tr_arr_ib(k,n2,2):", tr_arr_ib(k,n2,2),", lev_low:",lev_low,", lev_up:",lev_up
      end if

! save the current deepest values; will be overwritten most of the time in lines 707 to 713.
! kh 08.03.21 use UV_ib buffered values here
      uo_keel(m)=UV_ib(1,k,n2)
      vo_keel(m)=UV_ib(2,k,n2)

! kh 15.03.21 use tr_arr_ib buffered values here
      T_keel(m)=tr_arr_ib(k,n2,1)
      S_keel(m)=tr_arr_ib(k,n2,2)
    end if

   end do innerloop
 
   ! divide by depth over which was integrated
   uo_dz(m)=uo_dz(m)/abs(depth_ib)
   vo_dz(m)=vo_dz(m)/abs(depth_ib)
   T_dz(m)=T_dz(m)/abs(depth_ib)
   S_dz(m)=S_dz(m)/abs(depth_ib)
   
   !if(T_dz(m)>70.) then
   ! write(*,*) 'nodeloop, dz:',dz,'z_up:',coord_nod3D(3, n_up),'z_lo:',coord_nod3D(3, n_low),'depth:',depth_ib,'T_up:',tracer(n_up,2),'T_lo:',tracer(n_low,2)
   !end if

 end do nodeloop !loop over all nodes of iceberg element
       
 contains
 
  real function interpol1D(x0,f0,x1,f1,x)
   implicit none
   real, intent(IN) :: x0,f0,x1,f1,x
   real :: frac
   
   frac = (f1 - f0)/(x1 - x0)
   interpol1D = f0 + frac * (x - x0)
  	
  end function interpol1D
  
end subroutine iceberg_average_andkeel


 !***************************************************************************************************************************
 !***************************************************************************************************************************


!subroutine iceberg_avvelo_layer3(uo_dz,vo_dz,depth_ib,iceberg_elem)
!  use g_parsup 		!for myDim_nod3D, eDim_nod3D
!  !use o_mesh		!for coord_nod3D
!  use o_arrays		!for uf
!  implicit none
!  
!  !inout
!  real, dimension(3),intent(out):: uo_dz
!  real, dimension(3),intent(out):: vo_dz
!  real,	intent(in)		:: depth_ib
!  integer, intent(in)		:: iceberg_elem
!  !local
!  real :: displ3, dz, depth, ufkeel1, ufkeel2
!  integer :: m, k, n2, n_up, n_low, l3_node
!  
!  displ3=myDim_nod3d+eDim_nod3D
!  uo_dz=0.
!  vo_dz=0.
!  
! do m=1, 3
!   n2=elem2D_nodes(m,iceberg_elem)
!   depth=0.0
!   
!   if(num_layers_below_nod2d(n2)+1 <= 3) then
!   !the following do loop will not start then
!   write(*,*) 'no ocean velocity computed'
!   call par_ex
!   stop 'NO OCEAN VELOCITY COMPUTED'
!   end if
!   
!   !start with second node below n2; neglect upper layers
!   do k=4, num_layers_below_nod2d(n2)+1 
!    ! ..consider all neighboring pairs (n_up,n_low) of 3D nodes
!    ! below n2..
!    n_up = nod3d_below_nod2d(k-1, n2) ! for k=4 we have second node below n2
!    n_low= nod3d_below_nod2d(k,   n2)
!    dz   = abs( coord_nod3D(3, n_up)-coord_nod3D(3, n_low) )
!    
!    if(k .EQ. 4) then
!      l3_node = n_up
!      !WHAT IF ICEBERG IS SO SMALL THAT THE DRAFT IS 
!      !ABOVE THE THIRD LAYER?  
!      if(abs(coord_nod3D(3, l3_node))>= abs(depth_ib)) then
!        !just take layer 3 value
!	write(*,*) 'took layer 3 value'	
!        uo_dz(m)= uf(l3_node)
!        vo_dz(m)= uf(l3_node+displ3)
!        depth=1.
!        exit
!      end if
!    end if
!    
!    ! if the lowest z coord is below the iceberg draft, exit
!    if ( abs(coord_nod3D(3, n_low))>= abs(depth_ib)) then
!      dz = abs( coord_nod3D(3, n_up)- depth_ib )
!      ufkeel1 = interpol1D(	abs(coord_nod3D(3, n_up)),UV(1,k-1,n2), 	&
!      				abs(coord_nod3D(3, n_low)),UV(1,k,n2), 	&
!				abs(depth_ib) )
!      ufkeel2 = interpol1D(	abs(coord_nod3D(3, n_up)),UV(2,k-1,n2), 	&
!      				abs(coord_nod3D(3, n_low)),UV(2,k,n2),	&
!				abs(depth_ib) ) 
!      uo_dz(m)=uo_dz(m)+0.5*(UV(1,k-1,n2)+ ufkeel1)*dz 
!      vo_dz(m)=vo_dz(m)+0.5*(UV(2,k-1,n2)+ ufkeel2)*dz
!      depth=depth+dz
!      exit
!
!    else	
!      ! .. and sum up the layer-integrated velocities:
!      uo_dz(m)=uo_dz(m)+0.5*(UV(1,k-1,n2)+UV(1,k,n2))*dz
!      vo_dz(m)=vo_dz(m)+0.5*(UV(2,k-1,n2)+UV(2,k,n2))*dz
!      depth=depth+dz
!	 
!    end if
!   end do
!   
!   ! divide by depth over which was integrated
!   uo_dz(m)=uo_dz(m)/abs(depth)
!   vo_dz(m)=vo_dz(m)/abs(depth)
!         
! end do !loop over all nodes of iceberg element
!       
! contains
! 
!  real function interpol1D(x0,f0,x1,f1,x)
!   implicit none
!   real, intent(IN) :: x0,f0,x1,f1,x
!   real :: frac
!   
!   frac = (f1 - f0)/(x1 - x0)
!   interpol1D = f0 + frac * (x - x0)
!  	
!  end function interpol1D 
!end subroutine iceberg_avvelo_layer3

 !***************************************************************************************************************************
 !***************************************************************************************************************************


subroutine iceberg_avvelo(mesh, uo_dz,vo_dz,depth_ib,iceberg_elem)
  !use o_mesh
  USE MOD_MESH
  use o_param
  use i_therm_param
  use i_param
  use i_arrays
  use g_parsup
  
! kh 17.03.21 specification of structures used
  use o_arrays, only: UV_ib, tr_arr_ib, Z_3d_n_ib

  use g_clock
  use g_forcing_arrays
  use g_rotate_grid
  
  implicit none

  REAL, DIMENSION(3), INTENT(OUT) :: uo_dz
  REAL, DIMENSION(3), INTENT(OUT) :: vo_dz
  REAl,               INTENT(IN)  :: depth_ib
  INTEGER,            INTENT(IN)  :: iceberg_elem

  real           :: lev_up, lev_low  
  integer        :: m, k, n2, n_up, n_low
  ! depth over which is integrated (layer and sum)
  real           :: dz, ufkeel1, ufkeel2
  ! variables for velocity correction
  real           :: delta_depth, u_bottom_x, u_bottom_y

type(t_mesh), intent(in) , target :: mesh
#include "associate_mesh.h"
  
  
  ! loop over all nodes of the iceberg element
  do m=1, 3
   !for each 2D node of the iceberg element..
   n2=elem2D_nodes(m,iceberg_elem)
   uo_dz(m)=0.0
   vo_dz(m)=0.0

   ! ..consider all neighboring pairs (n_up,n_low) of 3D nodes
   ! below n2..
   do k=1, nl+1

! kh 18.03.21 use zbar_3d_n_ib buffered values here
    if( k==1 ) then
        lev_up = 0.0
    else
        lev_up = Z_3d_n_ib(k-1, n2)
    end if
    lev_low = Z_3d_n_ib(k, n2)
    
    if (lev_up==lev_low) then
        exit
    end if
    dz = abs( lev_low - lev_up )
	
    if(dz < 1) then
      !write(*,*) 'z coord of up node', n_up, ':', coord_nod3D(3, n_up), 'z coord of low node', n_low, ':', coord_nod3D(3, n_low)
      call par_ex
      stop
    end if
	
    ! if the lowest z coord is below the iceberg draft, exit
    if ( abs(lev_low)>= abs(depth_ib)) then
  
      dz = abs( lev_up - depth_ib )

      if( k==1 ) then
          ufkeel1 = UV_ib(1,k,n2)
          ufkeel2 = UV_ib(2,k,n2)
          uo_dz(m)= ufkeel1*dz 
          vo_dz(m)= ufkeel2*dz
      else
! kh 08.03.21 use UV_ib buffered values here
          ufkeel1 = interpol1D(abs(lev_up),UV_ib(1,k-1,n2),abs(lev_low),UV_ib(1,k,n2),abs(depth_ib))
          ufkeel2 = interpol1D(abs(lev_up),UV_ib(2,k-1,n2),abs(lev_low),UV_ib(2,k,n2),abs(depth_ib))
! kh 08.03.21 use UV_ib buffered values here
          uo_dz(m)=uo_dz(m)+0.5*(UV_ib(1,k-1,n2)+ ufkeel1)*dz 
          vo_dz(m)=vo_dz(m)+0.5*(UV_ib(2,k-1,n2)+ ufkeel2)*dz
      end if

      exit
	 
    else	
      if( k==1 ) then
        cycle
      end if
      ! .. and sum up the layer-integrated velocities:
! kh 08.03.21 use UV_ib buffered values here
      uo_dz(m)=uo_dz(m)+0.5*(UV_ib(1,k-1,n2)+UV_ib(1,k,n2))*dz
      vo_dz(m)=vo_dz(m)+0.5*(UV_ib(2,k-1,n2)+UV_ib(2,k,n2))*dz
	 
    end if
   end do
 
   ! divide by depth over which was integrated
   uo_dz(m)=uo_dz(m)/abs(depth_ib)
   vo_dz(m)=vo_dz(m)/abs(depth_ib)
         
 end do !loop over all nodes of iceberg element
       
 contains
 
  real function interpol1D(x0,f0,x1,f1,x)
   implicit none
   real, intent(IN) :: x0,f0,x1,f1,x
   real :: frac
   
   frac = (f1 - f0)/(x1 - x0)
   interpol1D = f0 + frac * (x - x0)
  	
  end function interpol1D
  
end subroutine iceberg_avvelo

 !***************************************************************************************************************************
 !***************************************************************************************************************************


!subroutine iceberg_avvelo_ufkeel(uo_dz,vo_dz, uo_keel,vo_keel, depth_ib,iceberg_elem)
!  !use o_mesh
!  use o_param
!  use i_therm_param
!  use i_param
!  use i_arrays
!  use g_parsup
!  
!  use o_arrays         
!  use g_clock
!  use g_forcing_index
!  use g_forcing_arrays
!  use g_rotate_grid
!  
!  implicit none
!  
!  REAL, DIMENSION(3), INTENT(OUT) :: uo_dz
!  REAL, DIMENSION(3), INTENT(OUT) :: vo_dz
!  REAL, DIMENSION(3), INTENT(OUT) :: uo_keel
!  REAL, DIMENSION(3), INTENT(OUT) :: vo_keel
!  REAl,               INTENT(IN)  :: depth_ib
!  INTEGER,            INTENT(IN)  :: iceberg_elem
!  
!  integer        :: m, k, n2, n_up, n_low, displ3
!  ! depth over which is integrated (layer and sum)
!  real           :: dz, depth, ufkeel1, ufkeel2
!  ! variables for velocity correction
!  real           :: delta_depth, u_bottom_x, u_bottom_y
!  
!  displ3=myDim_nod3d+eDim_nod3D
!  
!  ! loop over all nodes of the iceberg element
!  do m=1, 3
!   !for each 2D node of the iceberg element..
!   n2=elem2D_nodes(m,iceberg_elem)
!   depth=0.0
!   uo_dz(m)=0.0
!   vo_dz(m)=0.0
!   uo_keel(m)=0.0
!   vo_keel(m)=0.0
!
!   ! ..consider all neighboring pairs (n_up,n_low) of 3D nodes
!   ! below n2..
!   do k=2, num_layers_below_nod2d(n2)+1
!    n_up = nod3d_below_nod2d(k-1, n2) ! for k=2 we have n2
!    n_low= nod3d_below_nod2d(k,   n2)
!
!    ! ..compute their distance dz ..
!    dz   = abs( coord_nod3D(3, n_up)-coord_nod3D(3, n_low) )
!	
!    if(dz < 1) then
!      write(*,*) 'dz 0!', dz
!      write(*,*) 'z coord of up node', n_up, ':', coord_nod3D(3, n_up), 'z coord of low node', n_low, ':', coord_nod3D(3, n_low)
!      call par_ex
!      stop
!    end if
!	
!    ! if the lowest z coord is below the iceberg draft, exit
!    if ( abs(coord_nod3D(3, n_low))>= abs(depth_ib)) then
!  
!      dz = abs( coord_nod3D(3, n_up)- depth_ib )
!
!      ufkeel1 = interpol1D(abs(coord_nod3D(3, n_up)),UV(1,k-1,n2),abs(coord_nod3D(3, n_low)),UV(1,k,n2),abs(depth_ib))
!      ufkeel2 = interpol1D(abs(coord_nod3D(3, n_up)),UV(2,k-1,n2),abs(coord_nod3D(3, n_low)),UV(2,k,n2),abs(depth_ib))
!
!      uo_dz(m)=uo_dz(m)+0.5*(UV(1,k-1,n2)+ ufkeel1)*dz 
!      vo_dz(m)=vo_dz(m)+0.5*(UV(2,k-1,n2)+ ufkeel2)*dz
!      uo_keel(m)=ufkeel1
!      vo_keel(m)=ufkeel2
!      exit
!	 
!    else	
!      ! .. and sum up the layer-integrated velocities:
!      uo_dz(m)=uo_dz(m)+0.5*(UV(1,k-1,n2)+UV(1,k,n2))*dz
!      vo_dz(m)=vo_dz(m)+0.5*(UV(2,k-1,n2)+UV(2,k,n2))*dz
!	 
!    end if
!   end do
! 
!   ! divide by depth over which was integrated
!   uo_dz(m)=uo_dz(m)/abs(depth_ib)
!   vo_dz(m)=vo_dz(m)/abs(depth_ib)
!         
! end do !loop over all nodes of iceberg element
!       
! contains
! 
!  real function interpol1D(x0,f0,x1,f1,x)
!   implicit none
!   real, intent(IN) :: x0,f0,x1,f1,x
!   real :: frac
!   
!   frac = (f1 - f0)/(x1 - x0)
!   interpol1D = f0 + frac * (x - x0)
!  	
!  end function interpol1D
!  
!end subroutine iceberg_avvelo_ufkeel


 !***************************************************************************************************************************
 !***************************************************************************************************************************


!subroutine init_global_tides
!  ! initialize tides: read constituents (amplitude and phase)
!  ! The global tidal data should be on a grid with longitude range 0-360 degree.
!  ! Assumption: the number of columes/rows of the grids for global tidal data
!  !             are the same for all constituents and types (U,V, z), though 
!  !             the exact lon/lat values can be (are) different for different
!  ! 	        types.
!  !
!  ! This is an adapted version of the subroutine init_tidal_opbnd for use of
!  ! global tides in iceberg case; external model data is interpolated on the
!  ! grid and can be added to the surface slope ssh computed by FESOM for process studies.
!  !
!  ! Thomas Rackow, 17.12.2010, some cleaning up still to be done
!  !
!  use o_param
!  use o_arrays
!  !use o_mesh
!  use g_rotate_grid
!  use g_parsup
!  implicit none
!  !
!  integer					:: i, j, n, num, fid, nodmin,m(1)
!  integer					:: num_lon_reg, num_lat_reg
!  integer                                       :: p_flag
!  character(2)    				:: cons_name
!  real(kind=8)					:: lon, lat, dep
!  real(kind=8)					:: x, y, x2, y2, d, dmin
!  real(kind=8), allocatable, dimension(:)	:: lon_reg_4u, lat_reg_4u
!  real(kind=8), allocatable, dimension(:)	:: lon_reg_4v, lat_reg_4v
!  real(kind=8), allocatable, dimension(:)	:: lon_reg_4z, lat_reg_4z
!  real(kind=8), allocatable, dimension(:,:)	:: amp_reg, pha_reg
!  real(kind=8), allocatable, dimension(:)	:: lon_opbnd_n2d, lat_opbnd_n2d
!  integer					:: mySize2D
!
!  mySize2D=myDim_nod2D+eDim_nod2D
!  
!  ! check
!  num=len(trim(tidal_constituent))
!  if(mod(num,2) /= 0 .or. num/2 /= nmbr_tidal_cons) then
!     write(*,*) 'wrong specification of tidal constituents in O_module'
!     call par_ex
!     stop
!  end if
!
!  ! allocate arrays
!  if(trim(tide_opbnd_type)=='Flather') then
!     allocate(opbnd_u_tide(nmbr_opbnd_t2D))
!     allocate(opbnd_v_tide(nmbr_opbnd_t2D))
!     allocate(opbnd_z_tide(nmbr_opbnd_t2D))
!     opbnd_u_tide=0.
!     opbnd_v_tide=0.
!     opbnd_z_tide=0.
!  elseif(trim(tide_opbnd_type)=='ssh') then
!     allocate(opbnd_z_tide(mySize2D))	!CHANGED: opbnd_z_tide(nmbr_opbnd_t2D)
!     allocate(opbnd_z0_tide(mySize2D))	!CHANGED: opbnd_z0_tide(nmbr_opbnd_t2D)
!     opbnd_z_tide=0.  
!     opbnd_z0_tide=0.
!  else !'vel'
!     allocate(opbnd_u_tide(nmbr_opbnd_t2D))
!     allocate(opbnd_v_tide(nmbr_opbnd_t2D))
!     opbnd_u_tide=0.
!     opbnd_v_tide=0.
!  end if
!
!  ! allocate tidal amp and pha arrays
!  allocate(tide_period_coeff(nmbr_tidal_cons))
!  tide_period_coeff=0.
!  if(trim(tide_opbnd_type)/='ssh') then
!     allocate(tide_u_amp(nmbr_opbnd_t2d, nmbr_tidal_cons))
!     allocate(tide_v_amp(nmbr_opbnd_t2d, nmbr_tidal_cons))
!     allocate(tide_u_pha(nmbr_opbnd_t2d, nmbr_tidal_cons))
!     allocate(tide_v_pha(nmbr_opbnd_t2d, nmbr_tidal_cons))
!     tide_u_amp=0.
!     tide_v_amp=0.
!     tide_u_pha=0.
!     tide_v_pha=0.
!  end if
!  if(trim(tide_opbnd_type)/='vel') then !'ssh' enters here
!     allocate(tide_z_amp(mySize2D, nmbr_tidal_cons))	!tide_z_amp(nmbr_opbnd_t2d, nmbr_tidal_cons)
!     allocate(tide_z_pha(mySize2D, nmbr_tidal_cons))	!tide_z_pha(nmbr_opbnd_t2d, nmbr_tidal_cons)
!     tide_z_amp=0.
!     tide_z_pha=0.
!  end if
!
!  ! read the regular grid (for velocity)
!  if(trim(tide_opbnd_type)/='ssh') then
!     open(101,file=trim(TideForcingPath)//'lonlat_4U_'//trim(tidemodelname)//'.dat', status='old')
!     read(101,*) num_lon_reg, num_lat_reg 
!     allocate(lon_reg_4u(num_lon_reg), lat_reg_4u(num_lat_reg))
!     do i=1,num_lon_reg
!        read(101,*) lon_reg_4u(i)
!     end do
!     do i=1,num_lat_reg
!        read(101,*) lat_reg_4u(i)
!     end do
!     close(101)
!     open(102,file=trim(TideForcingPath)//'lonlat_4V_'//trim(tidemodelname)//'.dat', status='old')
!     read(102,*) num_lon_reg, num_lat_reg 
!     allocate(lon_reg_4v(num_lon_reg), lat_reg_4v(num_lat_reg))
!     do i=1,num_lon_reg
!        read(102,*) lon_reg_4v(i)
!     end do
!     do i=1,num_lat_reg
!        read(102,*) lat_reg_4v(i)
!     end do
!     close(102)
!  end if
!
!  ! read the regular grid (for ssh)
!  if(trim(tide_opbnd_type)/='vel') then !'ssh' enters here
!     open(103,file=trim(TideForcingPath)//'lonlat_4z_'//trim(tidemodelname)//'.dat', status='old')
!     read(103,*) num_lon_reg, num_lat_reg 
!     allocate(lon_reg_4z(num_lon_reg), lat_reg_4z(num_lat_reg))
!     do i=1,num_lon_reg
!        read(103,*) lon_reg_4z(i)
!     end do
!     do i=1,num_lat_reg
!        read(103,*) lat_reg_4z(i)
!     end do
!     close(103)
!  end if
!
!  ! allocate arrays for global data on regular grids
!  allocate(amp_reg(num_lon_reg, num_lat_reg), pha_reg(num_lon_reg, num_lat_reg)) 
!
!  ! open-boundary nodes coordinates
!  !allocate(lon_opbnd_n2d(nmbr_opbnd_n2d), lat_opbnd_n2d(nmbr_opbnd_n2d))
!  !do i=1, nmbr_opbnd_n2d
!  !   n=opbnd_n2d(i)
!  !   if(rotated_grid) then
!  !      call r2g(lon, lat, coord_nod2d(1,n), coord_nod2d(2,n))
!  !      lon_opbnd_n2d(i)=lon/rad   ! change unit to degree
!  !      lat_opbnd_n2d(i)=lat/rad
!  !   else
!  !      lon_opbnd_n2d(i)=coord_nod2d(1,n)/rad
!  !      lat_opbnd_n2d(i)=coord_nod2d(2,n)/rad
!  !   end if
!  !   ! change lon range to [0 360]
!  !   if(lon_opbnd_n2d(i)<0.) lon_opbnd_n2d(i)=lon_opbnd_n2d(i) + 360.0  
!  !end do
!  allocate(lon_opbnd_n2d(mySize2D), lat_opbnd_n2d(mySize2D))
!  do i=1, mySize2D
!  !   n=opbnd_n2d(i)
!     if(rotated_grid) then
!        call r2g(lon, lat, coord_nod2d(1,i), coord_nod2d(2,i))
!        lon_opbnd_n2d(i)=lon/rad   ! change unit to degree
!        lat_opbnd_n2d(i)=lat/rad
!     else
!        lon_opbnd_n2d(i)=coord_nod2d(1,i)/rad
!        lat_opbnd_n2d(i)=coord_nod2d(2,i)/rad
!     end if
!     ! change lon range to [0 360]
!     if(lon_opbnd_n2d(i)<0.) lon_opbnd_n2d(i)=lon_opbnd_n2d(i) + 360.0  
!  end do
!  
!
!  ! read and interpolate each constituent (for U,V,z and their phase) 
!  do n=1,nmbr_tidal_cons
!     cons_name=tidal_constituent(2*n-1:2*n)
!     call tide_period(cons_name, tide_period_coeff(n))
!
!     if(trim(tide_opbnd_type)/='ssh') then
!        fid=103+n
!        open(fid,file=trim(TideForcingPath)//'U_'//cons_name//'_'//trim(tidemodelname)//'.dat', status='old')
!        do i=1, num_lon_reg
!           do j=1, num_lat_reg
!              read(fid, *) amp_reg(i,j)         
!           end do
!        end do
!        do i=1, num_lon_reg
!           do j=1, num_lat_reg
!              read(fid, *) pha_reg(i,j)         
!           end do
!        end do
!        close(fid)
!
!        p_flag=0
!        call interp_2d_field(num_lon_reg, num_lat_reg, lon_reg_4u, lat_reg_4u, amp_reg, &
!             nmbr_opbnd_n2D, lon_opbnd_n2d, lat_opbnd_n2d, tide_u_amp(1:nmbr_opbnd_n2d,n),p_flag)
!        tide_u_amp(:,n)=tide_u_amp(:,n)/opbnd_dep  ! change transport (m^2/s) to depth mean velocity (m/s)
!
!        p_flag=1
!        call interp_2d_field(num_lon_reg, num_lat_reg, lon_reg_4u, lat_reg_4u, pha_reg, &
!             nmbr_opbnd_n2D, lon_opbnd_n2d, lat_opbnd_n2d, tide_u_pha(1:nmbr_opbnd_n2d,n),p_flag)
!        tide_u_pha(:,n)=tide_u_pha(:,n)*rad   ! change units of phase from degree to radian 
!
!
!        open(fid,file=trim(TideForcingPath)//'V_'//cons_name//'_'//trim(tidemodelname)//'.dat', status='old')
!        do i=1, num_lon_reg
!           do j=1, num_lat_reg
!              read(fid, *) amp_reg(i,j)         
!           end do
!        end do
!        do i=1, num_lon_reg
!           do j=1, num_lat_reg
!              read(fid, *) pha_reg(i,j)         
!           end do
!        end do
!        close(fid)
!
!        p_flag=0
!        call interp_2d_field(num_lon_reg, num_lat_reg, lon_reg_4v, lat_reg_4v, amp_reg, &
!             nmbr_opbnd_n2D, lon_opbnd_n2d, lat_opbnd_n2d, tide_v_amp(1:nmbr_opbnd_n2d,n),p_flag)
!        tide_v_amp(:,n)=tide_v_amp(:,n)/opbnd_dep  ! change transport (m^2/s) to depth mean velocity (m/s)
!
!        p_flag=1
!        call interp_2d_field(num_lon_reg, num_lat_reg, lon_reg_4v, lat_reg_4v, pha_reg, &
!             nmbr_opbnd_n2D, lon_opbnd_n2d, lat_opbnd_n2d, tide_v_pha(1:nmbr_opbnd_n2d,n),p_flag)
!        tide_v_pha(:,n)=tide_v_pha(:,n)*rad   ! change units of phase from degree to radian 
!
!     end if
!
!
!     if(trim(tide_opbnd_type)/='vel') then !'ssh' enters here
!        fid=103+n
!        open(fid,file=trim(TideForcingPath)//'z_'//cons_name//'_'//trim(tidemodelname)//'.dat', status='old')
!        do i=1, num_lon_reg
!           do j=1, num_lat_reg
!              read(fid, *) amp_reg(i,j)         
!           end do
!        end do
!        do i=1, num_lon_reg
!           do j=1, num_lat_reg
!              read(fid, *) pha_reg(i,j)         
!           end do
!        end do
!        close(fid)
!
!        p_flag=0
!	!CHANGED:
!        !call interp_2d_field(num_lon_reg, num_lat_reg, lon_reg_4z, lat_reg_4z, amp_reg, &
!        !     nmbr_opbnd_n2D, lon_opbnd_n2d, lat_opbnd_n2d, tide_z_amp(1:nmbr_opbnd_n2d,n),p_flag)
!	call interp_2d_field(num_lon_reg, num_lat_reg, lon_reg_4z, lat_reg_4z, amp_reg, &
!             mySize2D, lon_opbnd_n2d, lat_opbnd_n2d, tide_z_amp(1:mySize2D,n),p_flag)
!	     
!        p_flag=1
!	!CHANGED:
!        !call interp_2d_field(num_lon_reg, num_lat_reg, lon_reg_4z, lat_reg_4z, pha_reg, &
!        !     nmbr_opbnd_n2D, lon_opbnd_n2d, lat_opbnd_n2d, tide_z_pha(1:nmbr_opbnd_n2d,n),p_flag)
!	call interp_2d_field(num_lon_reg, num_lat_reg, lon_reg_4z, lat_reg_4z, pha_reg, &
!             mySize2D, lon_opbnd_n2d, lat_opbnd_n2d, tide_z_pha(1:mySize2D,n),p_flag)
!        tide_z_pha(:,n)=tide_z_pha(:,n)*rad   ! change units of phase from degree to radian 
!     end if
!  end do
!
!  ! amplifying the magnitude for process studies
!  if(trim(tide_opbnd_type)/='ssh') then
!     tide_u_amp=tide_u_amp*tide_amplify_coeff
!     tide_v_amp=tide_v_amp*tide_amplify_coeff
!  end if
!  if(trim(tide_opbnd_type)/='vel') then !'ssh' enters here
!     tide_z_amp=tide_z_amp*tide_amplify_coeff
!  end if
!
!
!  ! deallocate temporary arrays
!  deallocate(lat_opbnd_n2d, lon_opbnd_n2d)
!  deallocate(pha_reg, amp_reg)
!  if(trim(tide_opbnd_type)/='ssh') deallocate(lat_reg_4v,lon_reg_4v,lat_reg_4u,lon_reg_4u)
!  if(trim(tide_opbnd_type)/='vel') deallocate(lat_reg_4z,lon_reg_4z) !for 'ssh'
!
!  if(mype==0) then 
!   write(*,*) 'Global Tidal Constituents ', trim(tidal_constituent), ' have been loaded for Iceberg Case'
!   write(*,*) '*************************************************************'
!  end if 
!end subroutine init_global_tides

 !***************************************************************************************************************************
 !***************************************************************************************************************************

!subroutine update_global_tides  
!  !
!  ! This is an adapted version of the subroutine update_tidal_opbnd;
!  ! Thomas Rackow, 17.12.2010, some clean up still to be done
!  !
!  use o_param
!  use o_arrays
!  use g_config !for istep
!  !use o_mesh
!  implicit none
!  integer		:: n
!  real(kind=8)		:: aux
!  !
!  aux=2.0*pi*istep*dt
!  if(trim(tide_opbnd_type)=='Flather') then
!     opbnd_u_tide=0.
!     opbnd_v_tide=0.
!     opbnd_z_tide=0.
!     do n=1, nmbr_tidal_cons
!        opbnd_u_tide=opbnd_u_tide + tide_u_amp(:,n)*cos(aux/tide_period_coeff(n)-tide_u_pha(:,n))
!        opbnd_v_tide=opbnd_v_tide + tide_v_amp(:,n)*cos(aux/tide_period_coeff(n)-tide_v_pha(:,n))
!        opbnd_z_tide=opbnd_z_tide + tide_z_amp(:,n)*cos(aux/tide_period_coeff(n)-tide_z_pha(:,n))
!     end do
!  elseif(trim(tide_opbnd_type)=='ssh') then
!     opbnd_z0_tide=opbnd_z_tide
!     opbnd_z_tide=0.
!     do n=1, nmbr_tidal_cons
!        opbnd_z_tide=opbnd_z_tide + tide_z_amp(:,n)*cos(aux/tide_period_coeff(n)-tide_z_pha(:,n))
!     end do
!  else
!     opbnd_u_tide=0.
!     opbnd_v_tide=0.
!     do n=1, nmbr_tidal_cons
!        opbnd_u_tide=opbnd_u_tide + tide_u_amp(:,n)*cos(aux/tide_period_coeff(n)-tide_u_pha(:,n))
!        opbnd_v_tide=opbnd_v_tide + tide_v_amp(:,n)*cos(aux/tide_period_coeff(n)-tide_v_pha(:,n))
!     end do
!  end if
!  !
!end subroutine update_global_tides
