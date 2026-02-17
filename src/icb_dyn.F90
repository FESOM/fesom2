module iceberg_dynamics
 USE MOD_MESH
 use MOD_PARTIT
 USE MOD_PARSUP
 use MOD_ICE
 USE MOD_DYN
 use iceberg_params
 use iceberg_element
 !use iceberg_step

implicit none

 public ::  iceberg_dyn
 public ::  iceberg_frozen
 public ::  iceberg_acceleration
 public ::  compute_areas
 public ::  iceberg_average_andkeel
 public ::  iceberg_avvelo

 contains

!==============================================================================
! calculates basically the new iceberg velocity; if melting is enabled, the
! iceberg dimensions are adjusted as well.
!
!   Thomas Rackow, 29.06.2010
!==============================================================================
subroutine iceberg_dyn(mesh, partit, ice, dynamics, ib, new_u_ib, new_v_ib, u_ib, v_ib, lon,lat, depth_ib, &
               height_ib, length_ib, width_ib, iceberg_elem, &
		       mass_ib, Ci, Ca, Co, Cda_skin, Cdo_skin, &
		       rho_ice, rho_air, rho_h2o, P_sill, conc_sill, frozen_in, &
		       file1, file2, P_ib, conci_ib, dt_ib, lastsubstep, &
		       f_u_ib_old, f_v_ib_old, l_semiimplicit, &
		       semiimplicit_coeff, AB_coeff, file3, rho_icb)

 use g_forcing_arrays 	!for u_wind, v_wind or u_wind_ib, v_wind_ib respectively
 use o_arrays, only: Tsurf_ib, Ssurf_ib
 use o_param		!for dt

 integer                :: ib_n_lvls, m
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

!LA 2023-03-07
 real, dimension(:), pointer    :: hi_ib3, conci_ib3, coriolis
 real, dimension(3) 	:: uo_keel, vo_keel, T_keel,S_keel, uo_dz, vo_dz, T_dz,S_dz!hi_ib3, conci_ib3, 
 real, dimension(:,:), allocatable 	:: arr_uo_dz, arr_vo_dz, arr_T_dz,arr_S_dz !hi_ib3, conci_ib3, 
 real 			:: uo_ib, vo_ib, ua_ib, va_ib, ui_ib, vi_ib, hi_ib, uo_skin_ib, vo_skin_ib
 real, dimension(:), allocatable :: arr_uo_ib, arr_vo_ib, arr_T_ave_ib, arr_S_ave_ib
 real 			:: Ao, Aa, Ai, Ad, fcoriolis
 real 			:: au_ib, av_ib
 real, dimension(2,2) 	:: SI_matrix
 real, dimension(2)	:: SI_velo
 real 			:: u_ib_tmp, v_ib_tmp, normold, normnew, abs_omib, abs_omib_skin, ocean_drag
 integer 		:: iter_ib, n, n2
 real          		:: M_b, M_v, M_e, M_bv, sst_ib, sss_ib ! meltrates (basal, lateral, erosion, lateral 'basal'), temp. & salinity
 real			:: T_ave_ib, S_ave_ib, T_keel_ib, S_keel_ib
 character,  intent(IN)	:: file3*80
 real, intent(IN)	:: rho_icb
 
type(t_ice)   , intent(inout), target :: ice
type(t_mesh), intent(in) , target :: mesh
type(t_partit), intent(inout), target :: partit
type(t_dyn)   , intent(inout), target :: dynamics
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"
  
n2=elem2D_nodes(1,iceberg_elem)
allocate(arr_uo_dz(3,nlevels_nod2D(n2)))
allocate(arr_vo_dz(3,nlevels_nod2D(n2)))
allocate(arr_T_dz(3,nlevels_nod2D(n2)))
allocate(arr_S_dz(3,nlevels_nod2D(n2)))
arr_uo_dz = 0.0
arr_vo_dz = 0.0
arr_T_dz = 0.0
arr_S_dz = 0.0

allocate(arr_uo_ib(nlevels_nod2D(n2)))
allocate(arr_vo_ib(nlevels_nod2D(n2)))
allocate(arr_T_ave_ib(nlevels_nod2D(n2)))
allocate(arr_S_ave_ib(nlevels_nod2D(n2)))
arr_uo_ib = 0.0
arr_vo_ib = 0.0
arr_T_ave_ib = 0.0
arr_S_ave_ib = 0.0

 !OCEAN VELOCITIES: 
 ! - (uo_ib, vo_ib)		: integrated mean velocity at location of iceberg
 ! - (uo_skin_ib, vo_skin_ib)	: velocity below the draft of the iceberg
 ! call iceberg_avvelo_ufkeel(uo_dz,vo_dz, uo_keel,vo_keel, depth_ib,iceberg_elem) 
 call iceberg_average_andkeel(mesh, partit, dynamics, uo_dz,vo_dz, uo_keel,vo_keel, T_dz,S_dz, T_keel,S_keel, depth_ib,iceberg_elem, ib)
 
 !OCEANIC VELOCITY uo_ib, vo_ib
 call FEM_3eval(mesh, partit, uo_ib,vo_ib,lon,lat,uo_dz,vo_dz,iceberg_elem)

 call iceberg_levelwise_andkeel(mesh, partit, dynamics, arr_uo_dz,arr_vo_dz, uo_keel,vo_keel, arr_T_dz,arr_S_dz, T_keel,S_keel, depth_ib,iceberg_elem, ib, ib_n_lvls)
 do n=1, ib_n_lvls
    call FEM_3eval(mesh, partit, arr_uo_ib(n),arr_vo_ib(n),lon,lat,arr_uo_dz(:,n),arr_vo_dz(:,n),iceberg_elem)
    call FEM_3eval(mesh, partit, arr_T_ave_ib(n),arr_S_ave_ib(n),lon,lat,arr_T_dz(:,n),arr_S_dz(:,n),iceberg_elem)
 end do 
 call FEM_3eval(mesh, partit, uo_skin_ib,vo_skin_ib,lon,lat,uo_keel,vo_keel,iceberg_elem)

 !TEMPERATURE AND SALINITY:
 ! - T_ave_ib, S_ave_ib		: Mean T & S (integrated) at location of iceberg
 ! - T_keel_ib, S_keel_ib	: T & S below the draft of the iceberg (depth_ib)
 call FEM_3eval(mesh, partit, T_keel_ib,S_keel_ib,lon,lat,T_keel,S_keel,iceberg_elem)


 !ATMOSPHERIC VELOCITY ua_ib, va_ib
 call FEM_eval(mesh, partit, ua_ib,va_ib,lon,lat,u_wind_ib,v_wind_ib,iceberg_elem)      

 !ICE VELOCITY ui_ib, vi_ib
 call FEM_eval(mesh, partit, ui_ib,vi_ib,lon,lat,ice%uice_ib,ice%vice_ib,iceberg_elem)    
  
 !ICE THICKNESS (CONCENTRATION) hi_ib, conci_ib
 hi_ib3    => ice%data(size(ice%data))%values(:) !ice%m_ice_ib(tmp_arr)
 conci_ib3 => ice%data(size(ice%data)-1)%values(:) !ice%a_ice_ib(tmp_arr) 
 call FEM_3eval(mesh, partit, hi_ib,conci_ib,lon,lat,hi_ib3,conci_ib3,iceberg_elem)
 P_ib = 20000. * hi_ib * exp(-20.*(1-conci_ib))
 
 call compute_areas(Ao, Aa, Ai, Ad, depth_ib, &
                    height_ib, length_ib, width_ib, hi_ib)

 coriolis => mesh%coriolis(:)
 fcoriolis = coriolis(iceberg_elem) * coriolis_scale(ib)

 call iceberg_acceleration( mesh, partit, dynamics, ib, au_ib, av_ib, Ao, Aa, Ai, Ad, 		&
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
  call FEM_eval(mesh, partit, sst_ib,sss_ib,lon,lat,Tsurf_ib,Ssurf_ib,iceberg_elem)

  call iceberg_meltrates(partit, mesh, M_b, M_v, M_e, M_bv, &
				u_ib,v_ib, arr_uo_ib,arr_vo_ib, ua_ib,va_ib, &
				sst_ib, length_ib, conci_ib, &
				uo_skin_ib, vo_skin_ib, T_keel_ib, S_keel_ib, depth_ib, height_ib, &
				arr_T_ave_ib, arr_S_ave_ib, ib, rho_icb, uo_ib, vo_ib, ib_n_lvls, iceberg_elem, nlevels_nod2D(n2))

  call iceberg_newdimensions(partit, ib, depth_ib,height_ib,length_ib,width_ib,M_b,M_v,M_e,M_bv, &
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
 
 !use iceberg_params, only : l_freeze !use 'capturing mechanism' of sea ice?
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


subroutine iceberg_acceleration(mesh, partit, dynamics, ib, au_ib, av_ib, Ao, Aa, Ai, Ad, &
                                uo_ib,vo_ib, uo_skin_ib, vo_skin_ib,  &
				ua_ib,va_ib, ui_ib,vi_ib, &
				u_ib, v_ib, mass_ib, fcoriolis, &
				Ci, Ca, Co, Cda_skin, Cdo_skin, &
				rho_ice, rho_air, rho_h2o, length_ib, &
				iceberg_elem, conci_ib, file1, file2, &
				lon_rad, lat_rad, output, f_u_ib_old, &
				f_v_ib_old, l_semiimplicit, &
				semiimplicit_coeff, AB_coeff )
 
 use o_param 	!for g
 use g_config	!for istep
 use MOD_PARTIT	!for mype
 use g_rotate_grid,  only: vector_r2g, vector_g2r
 !use iceberg_params, only: l_wave, l_tides, l_geo_out, surfslop_scale, ascii_out
 use MOD_MESH
 use MOD_DYN

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
 integer 		:: i, istep
 
type(t_mesh), intent(in) , target :: mesh
type(t_partit), intent(inout), target :: partit
type(t_dyn)   , intent(inout), target :: dynamics
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"

 !estimate wave height at the icebergs location (Bigg et al., 1997),
 !so wave_amplitude = 0.5 * wave_height = 0.5 * const. * abs(atm velo)**2
 vel_atm = sqrt(ua_ib**2 + va_ib**2)
 wave_amplitude = 0.5 * 0.02025 * vel_atm**2
 
 !assume that waves have same direction as the winds;
 !guard against zero wind speed to avoid division by zero
 if(vel_atm > 0.) then
   direction_u = ua_ib / vel_atm
   direction_v = va_ib / vel_atm
 else
   direction_u = 0.
   direction_v = 0.
 end if
  
 !absolute values of relative velocities 
 abs_omib = sqrt( (uo_ib - u_ib)**2 + (vo_ib - v_ib)**2 )
 abs_amib = sqrt( (ua_ib - u_ib)**2 + (va_ib - v_ib)**2 )
 abs_imib = sqrt( (ui_ib - u_ib)**2 + (vi_ib - v_ib)**2 )  
 abs_omib_skin = sqrt( (uo_skin_ib - u_ib)**2 + (vo_skin_ib - v_ib)**2 )
 
 ! u-components
 ocean_drag_u     = (0.5 * Co * rho_h2o * Ao  &
                   * abs_omib * uo_ib)/mass_ib	!calculate part of it implicitly

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
 call mean_gradient(mesh, partit, dynamics, iceberg_elem, lon_rad, lat_rad, nablaeta)
 
 !additional surface slope due to tides
 if(l_tides) then
   slope_tides_u = 0.0  !-g * sum( opbnd_z_tide(elem2D_nodes(:,iceberg_elem)) * bafux_2D(:,iceberg_elem) )
 else
   slope_tides_u = 0.0
 end if
   
 surface_slope_u = (-g * nablaeta(1) + slope_tides_u) * surfslop_scale(ib) !default scaling is 1.0
 
 
 ! v-components	       
 ocean_drag_v     = (0.5 * Co * rho_h2o * Ao  &
                   * abs_omib * vo_ib)/mass_ib	!calculate part of it implicitly
 
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
 
! !if(.false.) then
! if(output .AND. ascii_out) then
! 
!   accs_u_out(1) = ocean_drag_u - (0.5 * Co * rho_h2o * Ao * abs_omib * u_ib)/mass_ib
!   accs_u_out(2) = ocean_skin_u - (rho_h2o * Cdo_skin * Ad * abs_omib_skin * u_ib)/mass_ib
!   accs_u_out(3) = air_drag_u
!   accs_u_out(4) = air_skin_u
!   accs_u_out(5) = ice_drag_u
!   accs_u_out(6) = wave_radiation_u
!   accs_u_out(7) = surface_slope_u
!   accs_u_out(8) = fcoriolis*v_ib
!   vels_u_out(1) = uo_ib
!   vels_u_out(2) = uo_skin_ib
!   vels_u_out(3) = ua_ib
!   vels_u_out(4) = ui_ib
!  
!   accs_v_out(1) = ocean_drag_v - (0.5 * Co * rho_h2o * Ao * abs_omib * v_ib)/mass_ib
!   accs_v_out(2) = ocean_skin_v - (rho_h2o * Cdo_skin * Ad * abs_omib_skin * v_ib)/mass_ib
!   accs_v_out(3) = air_drag_v
!   accs_v_out(4) = air_skin_v
!   accs_v_out(5) = ice_drag_v
!   accs_v_out(6) = wave_radiation_v
!   accs_v_out(7) = surface_slope_v
!   accs_v_out(8) = -fcoriolis*u_ib
!   vels_v_out(1) = vo_ib
!   vels_v_out(2) = vo_skin_ib
!   vels_v_out(3) = va_ib
!   vels_v_out(4) = vi_ib
!   
!   if(l_geo_out) then
!    do i=1,8
!     call vector_r2g(accs_u_out(i), accs_v_out(i), lon_rad, lat_rad, 0)
!    end do
!    do i=1,4
!     call vector_r2g(vels_u_out(i), vels_v_out(i), lon_rad, lat_rad, 0)
!    end do   
!   end if
!   
!   !open(unit=icbID,file=file1,position='append')
!   !write(icbID,'(2I,12e15.7)') mype, istep, 	&
!   ! 	 accs_u_out(1),	accs_u_out(2),	&
!   ! 	 accs_u_out(3),	accs_u_out(4),	&
!   ! 	 accs_u_out(5),	accs_u_out(6),	&	
!   ! 	 accs_u_out(7),	accs_u_out(8), 	&
!   ! 	 vels_u_out(1),	vels_u_out(2),	&
!   ! 	 vels_u_out(3),	vels_u_out(4)
!   ! 	 
!   ! 	 
!   ! 	 
!   !close(icbID)  
! 
!   !open(unit=icbID,file=file2,position='append')      
!   !write(icbID,'(2I,12e15.7)') mype,istep, 	&
!   ! 	 accs_v_out(1),	accs_v_out(2),	&
!   ! 	 accs_v_out(3),	accs_v_out(4),	&
!   ! 	 accs_v_out(5),	accs_v_out(6),	&	
!   ! 	 accs_v_out(7),	accs_v_out(8),	&		 
!   ! 	 vels_v_out(1),	vels_v_out(2),	&
!   ! 	 vels_v_out(3),	vels_v_out(4)
!   !close(icbID)
! end if !output
  
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


subroutine iceberg_levelwise_andkeel(mesh, partit, dynamics, uo_dz,vo_dz, uo_keel,vo_keel, T_dz,S_dz, T_keel,S_keel, depth_ib,iceberg_elem, ib, ib_n_lvls)
  !--------------------------------------------------------------------------
  ! Per-level velocities and T/S at iceberg location, plus keel values.
  ! Uses Z_3d_n_ib (mid-level depths) for consistent index mapping to UV_ib.
  !--------------------------------------------------------------------------
  USE MOD_MESH
  use o_param
  use MOD_PARTIT
  use MOD_DYN

  use o_arrays, only: Tclim_ib, Sclim_ib !, UV_ib, Z_3d_n_ib

  use g_clock
  use g_forcing_arrays
  use g_rotate_grid
  
  implicit none

  REAL, DIMENSION(:,:), allocatable, INTENT(OUT) :: uo_dz
  REAL, DIMENSION(:,:), allocatable, INTENT(OUT) :: vo_dz
  REAL, DIMENSION(3), INTENT(OUT) :: uo_keel
  REAL, DIMENSION(3), INTENT(OUT) :: vo_keel
  REAL, DIMENSION(:,:), allocatable, INTENT(OUT) :: T_dz
  REAL, DIMENSION(:,:), allocatable, INTENT(OUT) :: S_dz
  REAL, DIMENSION(3), INTENT(OUT) :: T_keel
  REAL, DIMENSION(3), INTENT(OUT) :: S_keel
  REAl,               INTENT(IN)  :: depth_ib
  INTEGER                         :: ib_n_lvls_old
  INTEGER,            INTENT(IN)  :: iceberg_elem, ib
  INTEGER,            INTENT(OUT) :: ib_n_lvls
  INTEGER, dimension(3)           :: arr_ib_n_lvls
  REAL, dimension(:,:,:), pointer :: UV_ib

  real           :: lev_up, lev_low
  integer        :: m, k, n2, max_node_level_count
  ! depth over which is integrated (layer and sum)
  real           :: dz, ufkeel1, ufkeel2, Temkeel, Salkeel

type(t_mesh), intent(in) , target :: mesh
type(t_dyn), intent(in) , target :: dynamics
type(t_partit), intent(inout), target :: partit
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"

  UV_IB     => dynamics%uv_ib(:,:,:)

  do m=1,3
    if(m==1) then
        max_node_level_count = nlevels_nod2D(elem2D_nodes(m,iceberg_elem))
    else
        max_node_level_count = max(max_node_level_count, nlevels_nod2D(elem2D_nodes(m,iceberg_elem)))
    end if
  end do

  allocate(uo_dz(3,max_node_level_count))
  allocate(vo_dz(3,max_node_level_count))
  allocate(T_dz(3,max_node_level_count))
  allocate(S_dz(3,max_node_level_count))

  ib_n_lvls_old = 0
  ib_n_lvls = 0
  arr_ib_n_lvls = 0

  uo_dz     =   0.0
  vo_dz     =   0.0
  uo_keel   =   0.0
  vo_keel   =   0.0
  T_dz      =   0.0
  S_dz      =   0.0
  T_keel    =   0.0
  S_keel    =   0.0
  
  !LOOP: over all nodes of the iceberg element
  nodeloop: do m=1, 3
   !for each 2D node of the iceberg element..
   n2=elem2D_nodes(m,iceberg_elem)

   ! LOOP over mid-levels: Z_3d_n_ib(k) gives depth where UV_ib(k) lives
   innerloop: do k=1, nlevels_nod2D(n2)
    if( k==1 ) then
        lev_up = 0.0
    else
        lev_up = mesh%Z_3d_n_ib(k-1, n2)
    end if
    if( k <= nlevels_nod2D(n2)-1 ) then
        lev_low = mesh%Z_3d_n_ib(k, n2)
    else
        lev_low = mesh%zbar_n_bot(n2)          ! bottom boundary
    end if

    if (lev_up==lev_low) then
      arr_ib_n_lvls(m) = k
      exit innerloop
    end if
    dz = abs( lev_low - lev_up )

    ! Cavity check: use zbar_3d_n for geometric comparison
    if (use_cavity .AND. mesh%cavity_depth(elem2D_nodes(m,iceberg_elem)) /= 0.0 &
        .AND. abs(depth_ib) < abs(mesh%zbar_3d_n(k, n2))) then
      ! Iceberg draft is above the ocean surface at this cavity node
      uo_keel(m)=UV_ib(1,k-1,n2)
      vo_keel(m)=UV_ib(2,k-1,n2)

      T_keel(m)=Tclim_ib(k-1,n2)
      S_keel(m)=Sclim_ib(k-1,n2) ! check those choices with RT: OK

      uo_dz(m,k)=UV_ib(1,k-1,n2)
      vo_dz(m,k)=UV_ib(2,k-1,n2)
      T_dz(m,k)=Tclim_ib(k-1,n2)
      S_dz(m,k)=Sclim_ib(k-1,n2)
      exit innerloop

    ! Keel layer: mid-level k is at or below the iceberg draft
    else if( abs(lev_low)>=abs(depth_ib) ) then
      if( k==1 ) then
        ! Draft within first half-layer: use first mid-level value
        ufkeel1 = UV_ib(1,1,n2)
        ufkeel2 = UV_ib(2,1,n2)
        Temkeel = Tclim_ib(1,n2)
        Salkeel = Sclim_ib(1,n2)
      else if( k > nlevels_nod2D(n2)-1 ) then
        ! Last half-layer to bottom: no UV_ib(k), use piecewise constant
        ufkeel1 = UV_ib(1,k-1,n2)
        ufkeel2 = UV_ib(2,k-1,n2)
        Temkeel = Tclim_ib(k-1,n2)
        Salkeel = Sclim_ib(k-1,n2)
      else
        ! Interpolate keel values between mid-levels k-1 and k
        ufkeel1 = interpol1D(abs(lev_up),UV_ib(1,k-1,n2),abs(lev_low),UV_ib(1,k,n2),abs(depth_ib))
        ufkeel2 = interpol1D(abs(lev_up),UV_ib(2,k-1,n2),abs(lev_low),UV_ib(2,k,n2),abs(depth_ib))
        Temkeel = interpol1D(abs(lev_up),Tclim_ib(k-1,n2),abs(lev_low),Tclim_ib(k,n2),abs(depth_ib))
        Salkeel = interpol1D(abs(lev_up),Sclim_ib(k-1,n2),abs(lev_low),Sclim_ib(k,n2),abs(depth_ib))
      end if
      
      uo_dz(m,k)=ufkeel1 
      vo_dz(m,k)=ufkeel2
      T_dz(m,k)=Temkeel
      S_dz(m,k)=Salkeel

      uo_keel(m)=ufkeel1
      vo_keel(m)=ufkeel2
      T_keel(m)=Temkeel
      S_keel(m)=Salkeel

      arr_ib_n_lvls(m) = k
      exit innerloop

    ! Regular layer: iceberg extends deeper
    else
      if( k > nlevels_nod2D(n2)-1 ) then
        ! Last half-layer to bottom: no UV_ib(k), use piecewise constant
        uo_dz(m,k)=UV_ib(1,k-1,n2)
        vo_dz(m,k)=UV_ib(2,k-1,n2)
        T_dz(m,k)=Tclim_ib(k-1,n2)
        S_dz(m,k)=Sclim_ib(k-1,n2)
        uo_keel(m)=UV_ib(1,k-1,n2)
        vo_keel(m)=UV_ib(2,k-1,n2)
        T_keel(m)=Tclim_ib(k-1,n2)
        S_keel(m)=Sclim_ib(k-1,n2)
      else
        ! Store per-level value (UV_ib(k) lives at Z_3d_n_ib(k))
        uo_dz(m,k)=UV_ib(1,k,n2)
        vo_dz(m,k)=UV_ib(2,k,n2)
        T_dz(m,k)=Tclim_ib(k,n2)
        S_dz(m,k)=Sclim_ib(k,n2)
        ! Update keel to deepest level so far
        uo_keel(m)=UV_ib(1,k,n2)
        vo_keel(m)=UV_ib(2,k,n2)
        T_keel(m)=Tclim_ib(k,n2)
        S_keel(m)=Sclim_ib(k,n2)
      end if
    end if

   end do innerloop
 end do nodeloop !loop over all nodes of iceberg element
 
 do m=1,3
    if (arr_ib_n_lvls(m)==0) then
        cycle
    else
        ib_n_lvls_old = ib_n_lvls
        ib_n_lvls = arr_ib_n_lvls(m)
        if ((ib_n_lvls_old.ne.0) .and. ib_n_lvls_old<ib_n_lvls) ib_n_lvls=ib_n_lvls_old
    end if
 end do

 contains
 
  real function interpol1D(x0,f0,x1,f1,x)
   implicit none
   real, intent(IN) :: x0,f0,x1,f1,x
   real :: frac
   
   frac = (f1 - f0)/(x1 - x0)
   interpol1D = f0 + frac * (x - x0)
  	
  end function interpol1D
end subroutine iceberg_levelwise_andkeel

!***************************************************************************************************************************
!***************************************************************************************************************************


subroutine iceberg_average_andkeel(mesh, partit, dynamics, uo_dz,vo_dz, uo_keel,vo_keel, T_dz,S_dz, T_keel,S_keel, depth_ib,iceberg_elem, ib)
  !--------------------------------------------------------------------------
  ! Depth-averaged velocities and T/S from surface to iceberg draft, plus
  ! keel values.  Uses Z_3d_n_ib with trapezoidal integration (as iceberg_avvelo).
  !--------------------------------------------------------------------------
  USE MOD_MESH
  use o_param
  use MOD_PARTIT
  use MOD_DYN

  use o_arrays, only: Tclim_ib, Sclim_ib !, UV_ib, Z_3d_n_ib

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
  REAL, dimension(:,:,:), pointer :: UV_ib

  real           :: lev_up, lev_low
  integer        :: m, k, n2
  ! depth over which is integrated (layer and sum)
  real           :: dz, ufkeel1, ufkeel2, Temkeel, Salkeel

type(t_mesh), intent(in) , target :: mesh
type(t_dyn), intent(in) , target :: dynamics
type(t_partit), intent(inout), target :: partit
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"

  UV_IB     => dynamics%uv_ib(:,:,:)

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

   ! LOOP over mid-levels: Z_3d_n_ib(k) gives depth where UV_ib(k) lives
   innerloop: do k=1, nlevels_nod2D(n2)
    if( k==1 ) then
        lev_up = 0.0
    else
        lev_up = mesh%Z_3d_n_ib(k-1, n2)
    end if

    if( k <= nlevels_nod2D(n2)-1 ) then
        lev_low = mesh%Z_3d_n_ib(k, n2)
    else
        lev_low = mesh%zbar_n_bot(n2)          ! bottom boundary
    end if

    if (lev_up==lev_low) then
      exit innerloop
    end if
    dz = abs( lev_low - lev_up )

    ! Cavity check: use zbar_3d_n for geometric comparison
    if (use_cavity .AND. mesh%cavity_depth(elem2D_nodes(m,iceberg_elem)) /= 0.0 &
        .AND. abs(depth_ib) < abs(mesh%zbar_3d_n(k, n2))) then
      ! Iceberg draft is above the ocean surface at this cavity node
      uo_dz(m)=UV_ib(1,k-1,n2)*abs(depth_ib)
      vo_dz(m)=UV_ib(2,k-1,n2)*abs(depth_ib)
      uo_keel(m)=UV_ib(1,k-1,n2)
      vo_keel(m)=UV_ib(2,k-1,n2)

      T_dz(m)=Tclim_ib(k-1,n2)*abs(depth_ib)
      S_dz(m)=Sclim_ib(k-1,n2)*abs(depth_ib)
      T_keel(m)=Tclim_ib(k-1,n2)
      S_keel(m)=Sclim_ib(k-1,n2) ! check those choices with RT: OK

      exit innerloop

    ! Keel layer: mid-level k is at or below the iceberg draft
    else if( abs(lev_low)>=abs(depth_ib) ) then
      dz = abs( lev_up - depth_ib )

      if( k==1 ) then
        ! Draft within first half-layer: piecewise constant
        ufkeel1 = UV_ib(1,1,n2)
        ufkeel2 = UV_ib(2,1,n2)
        Temkeel = Tclim_ib(1,n2)
        Salkeel = Sclim_ib(1,n2)
        uo_dz(m)=ufkeel1*dz 
        vo_dz(m)=ufkeel2*dz
        T_dz(m)=Temkeel*dz
        S_dz(m)=Salkeel*dz
      else if( k > nlevels_nod2D(n2)-1 ) then
        ! Last half-layer to bottom: no UV_ib(k), use piecewise constant
        ufkeel1 = UV_ib(1,k-1,n2)
        ufkeel2 = UV_ib(2,k-1,n2)
        Temkeel = Tclim_ib(k-1,n2)
        Salkeel = Sclim_ib(k-1,n2)
        uo_dz(m)=uo_dz(m)+ ufkeel1*dz
        vo_dz(m)=vo_dz(m)+ ufkeel2*dz
        T_dz(m)=T_dz(m)+ Temkeel*dz
        S_dz(m)=S_dz(m)+ Salkeel*dz
      else
        ! Interpolate keel values between mid-levels k-1 and k
        ufkeel1 = interpol1D(abs(lev_up),UV_ib(1,k-1,n2),abs(lev_low),UV_ib(1,k,n2),abs(depth_ib))
        ufkeel2 = interpol1D(abs(lev_up),UV_ib(2,k-1,n2),abs(lev_low),UV_ib(2,k,n2),abs(depth_ib))
        Temkeel = interpol1D(abs(lev_up),Tclim_ib(k-1,n2),abs(lev_low),Tclim_ib(k,n2),abs(depth_ib))
        Salkeel = interpol1D(abs(lev_up),Sclim_ib(k-1,n2),abs(lev_low),Sclim_ib(k,n2),abs(depth_ib))
        ! Trapezoidal integration for partial layer to keel
        uo_dz(m)=uo_dz(m)+ 0.5*(UV_ib(1,k-1,n2) + ufkeel1)*dz
        vo_dz(m)=vo_dz(m)+ 0.5*(UV_ib(2,k-1,n2) + ufkeel2)*dz
        T_dz(m)=T_dz(m)+ 0.5*(Tclim_ib(k-1,n2) + Temkeel)*dz
        S_dz(m)=S_dz(m)+ 0.5*(Sclim_ib(k-1,n2) + Salkeel)*dz
      end if

      uo_keel(m)=ufkeel1
      vo_keel(m)=ufkeel2

      T_keel(m)=Temkeel
      S_keel(m)=Salkeel
      exit innerloop

    ! Regular layer: iceberg extends deeper
    else
      if( k==1 ) then
        ! First half-layer (surface to first mid-level): piecewise constant
        uo_dz(m)=uo_dz(m)+ UV_ib(1,1,n2)*dz
        vo_dz(m)=vo_dz(m)+ UV_ib(2,1,n2)*dz
        T_dz(m)=T_dz(m)+ Tclim_ib(1,n2)*dz
        S_dz(m)=S_dz(m)+ Sclim_ib(1,n2)*dz
        cycle
      end if

      if( k > nlevels_nod2D(n2)-1 ) then
        ! Last half-layer to bottom: no UV_ib(k), use piecewise constant
        uo_dz(m)=uo_dz(m)+ UV_ib(1,k-1,n2)*dz
        vo_dz(m)=vo_dz(m)+ UV_ib(2,k-1,n2)*dz
        T_dz(m)=T_dz(m)+ Tclim_ib(k-1,n2)*dz
        S_dz(m)=S_dz(m)+ Sclim_ib(k-1,n2)*dz
        uo_keel(m)=UV_ib(1,k-1,n2)
        vo_keel(m)=UV_ib(2,k-1,n2)
        T_keel(m)=Tclim_ib(k-1,n2)
        S_keel(m)=Sclim_ib(k-1,n2)
      else
        ! Trapezoidal integration between consecutive mid-levels
        uo_dz(m)=uo_dz(m)+ 0.5*(UV_ib(1,k-1,n2)+UV_ib(1,k,n2))*dz
        vo_dz(m)=vo_dz(m)+ 0.5*(UV_ib(2,k-1,n2)+UV_ib(2,k,n2))*dz
        T_dz(m)=T_dz(m)+ 0.5*(Tclim_ib(k-1,n2)+Tclim_ib(k,n2))*dz
        S_dz(m)=S_dz(m)+ 0.5*(Sclim_ib(k-1,n2)+Sclim_ib(k,n2))*dz
        ! Update keel to deepest level so far
        uo_keel(m)=UV_ib(1,k,n2)
        vo_keel(m)=UV_ib(2,k,n2)
        T_keel(m)=Tclim_ib(k,n2)
        S_keel(m)=Sclim_ib(k,n2)
      end if
    end if
 
   end do innerloop

   ! divide by depth over which was integrated
   uo_dz(m)=uo_dz(m)/abs(depth_ib)
   vo_dz(m)=vo_dz(m)/abs(depth_ib)
   T_dz(m)=T_dz(m)/abs(depth_ib)
   S_dz(m)=S_dz(m)/abs(depth_ib)

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

subroutine iceberg_avvelo(mesh, partit, dynamics, uo_dz,vo_dz,depth_ib,iceberg_elem)
  USE MOD_MESH
  use o_param
  use MOD_PARTIT
  use MOD_DYN

  use o_arrays, only: Tclim_ib, Sclim_ib !, UV_ib, Z_3d_n_ib

  use g_clock
  use g_forcing_arrays
  use g_rotate_grid
  
  implicit none

  REAL, DIMENSION(3), INTENT(OUT) :: uo_dz
  REAL, DIMENSION(3), INTENT(OUT) :: vo_dz
  REAl,               INTENT(IN)  :: depth_ib
  REAL, dimension(:,:,:), pointer :: UV_ib
  INTEGER,            INTENT(IN)  :: iceberg_elem

  real           :: lev_up, lev_low  
  integer        :: m, k, n2, n_up, n_low
  ! depth over which is integrated (layer and sum)
  real           :: dz, ufkeel1, ufkeel2
  ! variables for velocity correction
  real           :: delta_depth, u_bottom_x, u_bottom_y

type(t_mesh), intent(in) , target :: mesh
type(t_dyn), intent(in) , target :: dynamics
type(t_partit), intent(inout), target :: partit
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"

  UV_IB     => dynamics%uv_ib(:,:,:)
  ! loop over all nodes of the iceberg element
  do m=1, 3
   !for each 2D node of the iceberg element..
   n2=mesh%elem2D_nodes(m,iceberg_elem)
   uo_dz(m)=0.0
   vo_dz(m)=0.0

   ! ..consider all neighboring pairs (n_up,n_low) of 3D nodes
   ! below n2..
   do k=1, nl+1

! kh 18.03.21 use zbar_3d_n_ib buffered values here
    if( k==1 ) then
        lev_up = 0.0
    else
        lev_up = mesh%Z_3d_n_ib(k-1, n2)
    end if
    lev_low = mesh%Z_3d_n_ib(k, n2)
    
    if (lev_up==lev_low) then
        exit
    end if
    dz = abs( lev_low - lev_up )
	
    if(dz < 1) then
      !write(*,*) 'z coord of up node', n_up, ':', coord_nod3D(3, n_up), 'z coord of low node', n_low, ':', coord_nod3D(3, n_low)
      call par_ex (partit%MPI_COMM_FESOM, partit%mype)
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
          ufkeel1 = interpol1D(abs(lev_up),UV_ib(1,k-1,n2),abs(lev_low),UV_ib(1,k,n2),abs(depth_ib))
          ufkeel2 = interpol1D(abs(lev_up),UV_ib(2,k-1,n2),abs(lev_low),UV_ib(2,k,n2),abs(depth_ib))
          uo_dz(m)=uo_dz(m)+0.5*(UV_ib(1,k-1,n2)+ ufkeel1)*dz 
          vo_dz(m)=vo_dz(m)+0.5*(UV_ib(2,k-1,n2)+ ufkeel2)*dz
      end if

      exit
	 
    else	
      if( k==1 ) then
        uo_dz(m)=uo_dz(m)+ UV_ib(1,k,n2)*dz
        vo_dz(m)=vo_dz(m)+ UV_ib(2,k,n2)*dz
        cycle
      end if
      ! .. and sum up the layer-integrated velocities:
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
end module iceberg_dynamics
