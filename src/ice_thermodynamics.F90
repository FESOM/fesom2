!===================================================================
subroutine cut_off()
use o_param
use i_arrays
implicit none

where(a_ice>1.0_WP)
 a_ice=1.0_WP
end where

where(a_ice<0.1e-8)
 a_ice=0.0_WP
end where

where(m_ice<0.1e-8)
 m_ice=0.0_WP
end where

end subroutine cut_off
!===================================================================
! Sea-ice thermodynamics routines
!
! Coded by N. Yakovlev and S. Danilov.
! Adjusted for upgraded model physics (NCEP forcing data; parameterization
! of ice-ocean heat flux considering friction velocity, etc) 
! by Ralph Timmermann.
! Adjusted for general forcing data and NlFs option, cleaned up, bug fixing,
! by Qiang Wang, 13.01.2009
!----------------------------------------------------------------------------

subroutine thermodynamics
  !
  ! For every surface node, this subroutine extracts the information
  ! needed for computation of thermodydnamics, calls the relevant
  ! subroutine, and returns the information to the vectors of prognostic
  ! variables.
  !------------------------------------------------------------------------
  
  use o_param
  use o_mesh
  use i_therm_param
  use i_param
  use i_arrays
  use g_config
  use g_forcing_param
  use g_forcing_arrays
  use g_parsup
  use g_comm_auto
  implicit none
  real(kind=WP)  :: h,hsn,A,fsh,flo,Ta,qa,rain,snow,runo,rsss,rsf,evap_in
  real(kind=WP)  :: ug,ustar,T_oc,S_oc,h_ml,t,ch,ce,ch_i,ce_i,fw,ehf,evap
  real(kind=WP)  :: ithdgr, ithdgrsn, iflice, hflatow, hfsenow, hflwrdout
  real(kind=WP)  :: lat
  integer        :: i, j, elem
  real(kind=WP), allocatable  :: ustar_aux(:)
  real*8  lid_clo

  rsss=ref_sss

  ! u_ice and v_ice are at nodes
  ! u_w, v_w is always on elements
  ! u_wind and v_wind are always at nodes
  ! ================
  ! Friction velocity 
  ! ================
  allocate(ustar_aux(myDim_nod2D+eDim_nod2D))
    DO i=1, myDim_nod2D
       ustar=0.0_WP
           ustar=((u_ice(i)-u_w(i))**2+ &
	                     (v_ice(i)-v_w(i))**2)
           ustar_aux(i)=sqrt(ustar*Cd_oce_ice)
    END DO	
  call exchange_nod(ustar_aux) !TODO Why do we need it?
  
  ! ================
  ! end: friction velocity 
  ! ================
  
  do i=1, myDim_nod2d+eDim_nod2D
     h       = m_ice(i)
     hsn     = m_snow(i)
     A       = a_ice(i)
     fsh     = shortwave(i)
     flo     = longwave(i)
     Ta      = Tair(i)
     qa      = shum(i)  
     if(precip_data_source=='NCEP') then
        if(Ta>=0.0_WP) then
           rain=prec_rain(i)
           snow=0.0_WP
        else
           rain=0.0_WP
           snow=prec_rain(i)
        endif
        evap_in=-evaporation(i) !evap_in: positive down
     else
        rain = prec_rain(i)
        snow = prec_snow(i)
        evap_in=0.0_WP
     end if
     runo    = runoff(i)
     ug      = sqrt(u_wind(i)**2+v_wind(i)**2)
     ustar   = ustar_aux(i)
     T_oc    = T_oc_array(i)      
     S_oc    = S_oc_array(i)
     if(ref_sss_local) rsss = S_oc
     t       = t_skin(i)   
     ch	     = Ch_atm_oce_arr(i)
     ce	     = Ce_atm_oce_arr(i)
     ch_i    = Ch_atm_ice
     ce_i    = Ce_atm_ice
     h_ml    = 10.0_WP       	         ! 10.0 or 30. used previously
     fw      = 0.0_WP
     ehf     = 0.0_WP
     lid_Clo=h0
     if (coord_nod2D(2,i)>0) then !TODO 2 separate pars for each hemisphere
       lid_clo=0.5
     else
       lid_clo=1.
     endif

     call therm_ice(h,hsn,A,fsh,flo,Ta,qa,rain,snow,runo,rsss, &
          ug,ustar,T_oc,S_oc,h_ml,t,ice_dt,ch,ce,ch_i,ce_i,evap_in,fw,ehf,evap, &
          rsf, ithdgr, ithdgrsn, iflice, hflatow, hfsenow, hflwrdout,lid_clo,i)
	
	ehf=0.0_WP
	fw=0.0_WP
	if (abs(ehf)>1.5e6 .and. mstep>60) then
! 	if (abs(ehf)>1.5e6 .and. mstep>60) then
		write(*,*) ' FOUND HEATFLUX> 0.7e6'
		write(*,*) ' mype            = ',mype
		write(*,*) ' mstep           = ',mstep
		write(*,*) ' i               = ',i
		write(*,*) ' lon/lat         =',coord_nod2D(:,i)/rad
		write(*,*)
		write(*,*) '___INPUT_FORCING______________________________'
		write(*,*) ' Ta              =',Ta
		write(*,*) ' qa              =',qa
		write(*,*) ' ug              =',ug
		write(*,*) ' ustar           =',ustar
		write(*,*) ' rain            =',rain
		write(*,*) ' snow            =',snow
		write(*,*) ' runo            =',runo
		write(*,*) ' evap            =',evap_in
		write(*,*) ' fsh             =',fsh
		write(*,*) ' flo             =',flo
		write(*,*)
		write(*,*) '___OUTPUT_FORCING__________________________'
		write(*,*) ' fw              =',fw
		write(*,*) ' fresh_wa_flux   =',fresh_wa_flux(i)
		write(*,*) ' ehf             =',ehf
		write(*,*) ' net_heat_flux   =',net_heat_flux(i)
		write(*,*) ' evap            =',-evap
		write(*,*)
		write(*,*) '___INPUT_OCEAN______________________________'
		write(*,*) ' T_oc            =',T_oc
		write(*,*) ' S_oc            =',S_oc
		write(*,*) ' ch              =',ch
		write(*,*) ' ce              =',ce
		write(*,*) ' ch_i            =',ch_i
		write(*,*) ' ce_i            =',ce_i
		write(*,*) ' t_skin          =',t_skin(i)
		write(*,*)
		write(*,*) '___OUTPUT_OCEAN______________________________'
		write(*,*) ' t_skin          =',t
		write(*,*)
		write(*,*) '___INPUT_ICE________________________________'
		write(*,*) ' h_old           =', m_ice(i)
		write(*,*) ' hsn_old         =', m_snow(i)
		write(*,*) ' A_old           =', a_ice(i)
		write(*,*) ' thdgr_old       =', thdgr(i)
		write(*,*) ' thdgrsn_old     =', thdgrsn(i)
		write(*,*) ' flice_old       =', flice(i)
		write(*,*) ' olat_heat_old   =', olat_heat(i)
		write(*,*) ' osen_heat_old   =', osen_heat(i)
		write(*,*) ' olwout_old      =', olwout(i)
		write(*,*)
		write(*,*) '___OUTPUT_ICE_______________________________'
		write(*,*) ' h_new           =', h
		write(*,*) ' hsn_new         =', hsn
		write(*,*) ' A_new           =', A
		write(*,*) ' thdgr_new       =', ithdgr
		write(*,*) ' thdgrsn_new     =', ithdgrsn
		write(*,*) ' flice_new       =', iflice
		write(*,*) ' olat_heat_new   =', hflatow
		write(*,*) ' osen_heat_new   =', hfsenow
		write(*,*) ' olwout_new      =', hflwrdout
		write(*,*)
		call par_ex(1)	
	end if

     m_ice(i)         = h
     m_snow(i)        = hsn
     a_ice(i)         = A
     t_skin(i)        = t
     fresh_wa_flux(i) = fw      !positive down
     net_heat_flux(i) = ehf     !positive down
     evaporation(i)   = -evap   !positive up

     thdgr(i)         = ithdgr
     thdgrsn(i)       = ithdgrsn
     flice(i)         = iflice
     olat_heat(i)     = hflatow
     osen_heat(i)     = hfsenow
     olwout(i)        = hflwrdout

  end do
     deallocate(ustar_aux)
end subroutine thermodynamics
!
!===================================================================
!
subroutine therm_ice(h,hsn,A,fsh,flo,Ta,qa,rain,snow,runo,rsss, &
     ug,ustar,T_oc,S_oc,H_ML,t,ice_dt,ch,ce,ch_i,ce_i,evap_in,fw,ehf,evap, &
     rsf, dhgrowth, dhsngrowth, iflice, hflatow, hfsenow, hflwrdout,lid_clo,n)
  ! Ice Thermodynamic growth model     
  !
  ! Input parameters:
  !------------------
  ! h - ice mass [m]
  ! hsn - snow mass [m]
  ! A - ice compactness
  ! fsh - shortwave radiation
  ! flo - longwave radiation
  ! Ta - air temperature
  ! qa - specific humidity
  ! rain - precipitation rain
  ! snow - precipitation snow
  ! runo - runoff
  ! ug - wind speed
  ! ustar - friction velocity
  ! T_oc, S_oc - ocean temperature and salinity beneath the ice (mixed layer)
  ! H_ML - mixed layer depth - should be specified.
  ! t - temperature of snow/ice top surface
  ! ice_dt - time step [s]
  ! ch - transfer coefficient for sensible heat (for open ocean)
  ! ce - transfer coefficient for evaporation   (for open ocean)
  ! ch_i - transfer coefficient for sensible heat (for ice)
  ! ce_i - transfer coefficient for evaporation   (for ice)  
  ! lid_clo - lid closing parameter
  ! Output parameters:
  !-------------------
  ! h - ice mass
  ! hsn - snow mass
  ! A - ice compactness
  ! t - temperature of snow/ice top surface
  ! fw - freshwater flux due to ice melting [m water/ice_dt]
  ! ehf - net heat flux at the ocean surface [W/m2]        !RTnew

  use i_therm_param
  use g_forcing_param,  only: precip_data_source
  
  use o_param
  use g_parsup
  
  implicit none

  integer k,n
  real*8  h,hsn,A,fsh,flo,Ta,qa,rain,snow,runo,rsss,evap_in
  real*8  ug,ustar,T_oc,S_oc,H_ML,t,ice_dt,ch,ce,ch_i,ce_i,fw,ehf
  real*8  dhgrowth,dhsngrowth,ahf,prec,subli,subli_i,rsf
  real*8  rhow,show,rhice,shice,sh,thick,thact,lat
  real*8  rh,rA,qhst,sn,hsntmp,o2ihf,evap
  real*8  iflice,hflatow,hfsenow,hflwrdout
  real*8, external  :: TFrez  ! Sea water freeze temperature.
  real*8  lid_clo
  ! Store ice thickness at start of growth routine
  dhgrowth=h  	  

  ! determine h(i,j)/a(i,j) = actual ice thickness.
  ! if snow layer is present, add hsn weighted with quotient
  ! of conductivities of ice and snow, according to 0-layer approach
  ! of Semtner (1976).   	    
  ! thickness at the ice covered part
  thick=hsn*(con/consn)/max(A,Armin)    ! Effective snow thickness
  thick=thick+h/max(A,Armin)            ! Effective total snow-ice thickness

  ! Growth rate for ice in open ocean
  rhow=0.0
  evap=0.0
  call obudget(qa,fsh,flo,T_oc,ug,ta,ch,ce,rhow,evap,hflatow,hfsenow,hflwrdout) 
  hflatow=hflatow*(1.0-A)
  hfsenow=hfsenow*(1.0-A)
  hflwrdout=hflwrdout*(1.0-A)

  ! add heat loss at open ocean due to melting snow fall
  !rhow=rhow+snow*1000.0/rhoice !qiang
  ! ice_dt and (1-A) will be multiplied afterwards

  ! growth rate of ice in ice covered part
  ! following Hibler 1984
  ! assuming ice thickness has an euqal, 7-level distribution from zero to two times h 
  rhice=0.0                      
  subli=0.0
  if (thick.gt.hmin) then
     do k=1,iclasses			  
        thact = (2*k-1)*thick/float(iclasses)  	! Thicknesses of actual ice class
        call budget(thact,hsn,t,ta,qa,fsh,flo,ug,S_oc,ch_i,ce_i,shice,subli_i) 
        !Thick ice K-class growth rate
        rhice=rhice+shice/float(iclasses)      	! Add to average heat flux
        subli=subli+subli_i/float(iclasses)
     end do
  end if

  ! Convert growth rates [m ice/sec] into growth per time step DT.
  rhow=rhow*ice_dt
  rhice=rhice*ice_dt

  ! Multiply ice growth of open water and ice
  ! with the corresponding areal fractions of grid cell
  show =rhow*(1.0-A)
  shice=rhice*A
  sh   =show+shice

  ! Store atmospheric heat flux, average over grid cell [W/m**2]
  ahf=-cl*sh/ice_dt   

  ! precipitation (into the ocean)
  prec=rain+runo+snow*(1.0-A)  	        ! m water/s

  ! snow fall above ice
  hsn=hsn+snow*ice_dt*A*rhowat*inv_rhosno	! Add snow fall to temporary snow thickness    !!!
  dhsngrowth=hsn   		        ! Store snow thickness after snow fall 

  ! evaporation/sublimation
  if(precip_data_source=='NCEP') then
     evap=evap_in*(1.0-A)
     subli=evap_in*A
  else
     evap=evap*(1.0-A)    		! m water/s
     subli=subli*A
  end if

  ! If there is atmospheric melting, first melt any snow that is present.
  ! Atmospheric heat flux available for melting
  ! sh = MINUS atm. heat flux / specific latent heat of sea ice
  ! Note: (sh<0) for melting, (sh>0) for freezing
  hsntmp= -min(sh,0.0_WP)*rhoice*inv_rhosno

  ! hsntmp is the decrease in snow thickness due to atmospheric melting
  ! [m/DT]. Do not melt more snow than available
  hsntmp=min(hsntmp,hsn)
  hsn=hsn-hsntmp  ! Update snow thickness after atmospheric snow melt

  ! Negative atmospheric heat flux left after melting of snow
  ! Note: (sh<0) and (hsntmp>0) for melting conditions
  ! hsntmp=0 for non-snow-melting conditions
  rh=sh+hsntmp*rhosno/rhoice
  h=max(h,0.0_WP)

  ! Compute heat content qhst of mixed layer - sea ice system
  !
  ! Total heat content is the sum of
  !	h	ice thickness after calculation of dynamic effects
  !	rh	change in ice thickness due to atmospheric forcing
  ! and heat available in mixed layer, with
  !	T_oc	temperature of ocean surface layer
  !	Tfrez	freezing point of sea water
  !	H_ML	thickness of uppermost layer
  !
  !RT:
  ! There are three possibilities to do this.
  ! 1.: Assume an instantaneous adjustment of mixed layer heat content.
  !     Any heat available is then instantaneously used to melt ice.
  !     (so-called ice-bath approach)
  !     This is what used to be used in the Lemke sea ice-mixed layer model.
  ! rh=rh-(T_oc-TFrez(S_oc))*H_ML*cc/cl
  ! qhst=h+rh 
  !
  ! 2.: Parameterize the ocean-to-ice heat flux (o2ihf)
  !     as a function of temperature difference. For a first step 
  !     we can assume a constant exchange coefficient gamma_t:
  ! o2ihf= (T_oc-TFrez(S_oc))*gamma_t*cc*A     &
  !        +(T_oc-Tfrez(S_oc))*H_ML/ice_dt*cc*(1.0-A) ! [W/m2]
  ! rh=rh-o2ihf*ice_dt/cl
  ! qhst=h+rh		                      	! [m]
  !
  ! 3.  Parameterize the ocean-to-ice heat flux (o2ihf)
  !     as a function of temperature difference and the
  !     friction velocity:
  o2ihf= (T_oc-TFrez(S_oc))*0.006*ustar*cc*A  &
       +(T_oc-Tfrez(S_oc))*H_ML/ice_dt*cc*(1.0-A)  	! [W/m2]
  rh=rh-o2ihf*ice_dt/cl
  qhst=h+rh		              		! [m]

  ! Melt snow if there is any ML heat content left (qhst<0).
  ! This may be the case if advection moves ice (with snow) to regions
  ! with a warm mixed layer.
  sn=hsn+min(qhst,0.0_WP)*rhoice*inv_rhosno

  ! New temporary snow thickness must not be negative:
  sn=max(sn,0.0_WP)

  ! Update snow and ice depth
  hsn=sn
  h=max(qhst,0.0_WP)     
  if (h.lt.1E-6) h=0.        ! Avoid very small ice thicknesses

  ! heat and fresh water fluxes
  if (mype==32 .and. mstep>=50 .and. n==1529) then
		write(*,*) ' --CHECK (B)--> heatflux in subroutine therm_ice()'
		write(*,*) 'h_new      = ',h
		write(*,*) 'h_old      = ',dhgrowth
		write(*,*) 'dhgrowth   =', h-dhgrowth
		write(*,*) 'hsn_new    = ',hsn
		write(*,*) 'hsn_old    = ',dhsngrowth
		write(*,*) 'dhsngrowth =', hsn-dhsngrowth
		
  endif
  dhgrowth=h-dhgrowth        ! Change in ice thickness due to thermodynamic effects
  dhsngrowth=hsn-dhsngrowth  ! Change in snow thickness due to thermodynamic melting

  ! (without snow fall). This is a negative value (MINUS snow melt)

  dhgrowth=dhgrowth/ice_dt       ! Conversion: 'per time step' -> 'per second'
  dhsngrowth=dhsngrowth/ice_dt   ! Conversion: 'per time step' -> 'per second'
  ! (radiation+turbulent) + freezing(-melting) sea-ice&snow 
  
!   if (mype==32 .and. mstep==247 .and. n==1529) then
  if (mype==32 .and. mstep>=50 .and. n==1529) then
		write(*,*) ' --CHECK (A)--> heatflux in subroutine therm_ice()'
		write(*,*) 'ahf        = ',ahf
		write(*,*) 'cl         = ',cl
		write(*,*) 'dhgrowth   = ',dhgrowth
		write(*,*) 'dhsngrowth = ',dhsngrowth
		write(*,*) 'rhosno     = ',rhosno
		write(*,*) 'rhoice     = ',rhoice
		write(*,*) 'ice_dt     = ',ice_dt
		write(*,*) 'cl*(dhgrowth+(rhosno/rhoice)*dhsngrowth)=',cl*(dhgrowth+(rhosno/rhoice)*dhsngrowth) 
		
  endif
  ehf = ahf + cl*(dhgrowth+(rhosno/rhoice)*dhsngrowth) 
  
  
  if (mype==32 .and. mstep==247 .and. n==1529) then
		write(*,*) ' --CHECK--> heatflux in subroutine therm_ice()'
		write(*,*) 'ehf = ',ehf
! 		call par_ex(1)
  endif
  
  
  
  
  ! (prec+runoff)+evap - freezing(+melting) ice&snow
#ifdef use_fullfreesurf
  fw= prec+evap - dhgrowth*rhoice*inv_rhowat - dhsngrowth*rhosno*inv_rhowat 
  rsf= -dhgrowth*rhoice*inv_rhowat*Sice
#else
  fw= prec+evap - dhgrowth*rhoice*inv_rhowat*(rsss-Sice)/rsss - dhsngrowth*rhosno*inv_rhowat 
#endif

  ! Changes in compactnesses (equation 16 of Hibler 1979)
  rh=-min(h,-rh)   ! Make sure we do not try to melt more ice than is available
  rA= rhow - o2ihf*ice_dt/cl !Qiang: it was -(T_oc-TFrez(S_oc))*H_ML*cc/cl, changed in June 2010
  !rA= rhow - (T_oc-TFrez(S_oc))*H_ML*cc/cl*(1.0-A)
  A=A + 0.5*min(rh,0.0_WP)*A/max(h,hmin) + max(rA,0.0_WP)*(1.-A)/lid_clo  !/h0   
  !meaning:           melting                         freezing
 
  A=min(A,h*1.e6)     ! A -> 0 for h -> 0
  A=min(max(A,0.0_WP),1.) ! A >= 0, A <= 1

  ! Flooding (snow to ice conversion)
  iflice=h
  call flooding(h,hsn)     
  iflice=(h-iflice)/ice_dt

  ! to maintain salt conservation for the current model version
  !(a way to avoid producing net salt from snow-type-ice) 
#ifdef use_fullfreesurf
  rsf=rsf-iflice*rhoice*inv_rhowat*Sice
#else
  fw=fw+iflice*rhoice*inv_rhowat*Sice/rsss
#endif   

evap=evap+subli
end subroutine therm_ice
!
!=====================================================================================
!
subroutine budget (hice,hsn,t,ta,qa,fsh,flo,ug,S_oc,ch_i,ce_i,fh,subli)
  ! Thick ice growth rate [m ice/sec]
  !
  ! INPUT:
  ! hice - actual ice thickness [m]
  ! hsn - snow thickness, used for albedo parameterization [m]
  ! t - temperature of snow/ice surface [C]
  ! ta - air temperature [C]
  ! qa - specific humidity [Kg/Kg]
  ! fsh - shortwave radiation [W/m**2]
  ! flo - longwave radiation  [W/m**2]
  ! ug - wind speed [m/sec]
  ! S_oc - ocean salinity for the temperature of the ice base calculation [ppt]
  ! ch_i - transfer coefficient for sensible heat (for ice)
  ! ce_i - transfer coefficient for evaporation   (for ice) 
  !
  ! OUTPUT: fh - growth rate
  !
  ! qiang: The formular for saturated humidity was modified according to Large/Yeager2004
  ! to allow direct comparison with the CORE results (Griffies et al. 2009). The new
  ! formular does not require sea level pressure.
  ! A similar change was also made for the obudget routine.
  ! It was found through experiments that the results are quite similar to that from the
  ! original code, and the simulated ice volume is only slightly larger after modification. 
  
  use i_therm_param
  implicit none

  integer iter, imax      ! Number of iterations
  real*8  hice,hsn,t,ta,qa,fsh,flo,ug,S_oc,ch_i,ce_i,fh
  real*8  hfsen,hfrad,hflat,hftot,subli         
  real*8  alb             ! Albedo of sea ice
  real*8  q1, q2	  ! coefficients for saturated specific humidity
  real*8  A1,A2,A3,B,C, d1, d2, d3   
  real*8, external :: TFrez

  data q1 /11637800.0/, q2 /-5897.8/ 
  data imax /5/

  ! set albedo
  ! ice and snow, freezing and melting conditions are distinguished.
  if (t<0.0) then	        ! freezing condition    
     if (hsn.gt.0.0) then	!   snow cover present  
        alb=albsn         	
     else              		!   no snow cover       
        alb=albi       	
     endif
  else			        ! melting condition     
     if (hsn.gt.0.0) then	!   snow cover present  
        alb=albsnm	    	
     else			!   no snow cover       
        alb=albim		
     endif
  endif

  d1=rhoair*cpair*Ch_i
  d2=rhoair*Ce_i
  d3=d2*clhi

  ! total incoming atmospheric heat flux
  A1=(1.0-alb)*fsh + flo + d1*ug*ta + d3*ug*qa   ! in LY2004 emiss is multiplied wiht flo
  ! NEWTON-RHAPSON TO GET TEMPERATURE AT THE TOP OF THE ICE LAYER

  do iter=1,imax

     B=q1*inv_rhoair*exp(q2/(t+tmelt))		! (saturated) specific humidity over ice
     A2=-d1*ug*t-d3*ug*B &
          -emiss_ice*boltzmann*((t+tmelt)**4)	! sensible and latent heat,and outward radiation
     A3=-d3*ug*B*q2/((t+tmelt)**2)		! gradient coefficient for the latent heat part
     C=con/hice                     		! gradient coefficient for downward heat conductivity
     A3=A3+C+d1*ug & 			! gradient coefficient for sensible heat and radiation 
          +4.0*emiss_ice*boltzmann*((t+tmelt)**3)    
     C=C*(TFrez(S_oc)-t)       		! downward conductivity term

     t=t+(A1+A2+C)/A3 		        ! NEW ICE TEMPERATURE AS THE SUM OF ALL COMPONENTS
  end do

  t=min(0.0,t)
  ! heat fluxes [W/m**2]:

  hfrad= (1.0-alb)*fsh &	        ! absorbed short wave radiation
       +flo &           	        ! long wave radiation coming in  ! in LY2004 emiss is multiplied
       -emiss_ice*boltzmann*((t+tmelt)**4) 	! long wave radiation going out

  hfsen=d1*ug*(ta-t) 			! sensible heat 
  subli=d2*ug*(qa-B) 			! sublimation
  hflat=clhi*subli                     	! latent heat

  hftot=hfrad+hfsen+hflat               ! total heat

  fh= -hftot/cl                         ! growth rate [m ice/sec]
  !                                      	+: ML gains energy, ice melts
  !                                      	-: ML loses energy, ice grows
  subli=subli*inv_rhowat                    ! negative upward

  return
end subroutine budget
!
!======================================================================================
!
subroutine obudget (qa,fsh,flo,t,ug,ta,ch,ce,fh,evap,hflatow,hfsenow,hflwrdout)  
  ! Ice growth rate for open ocean [m ice/sec]
  !
  ! INPUT:
  ! t - temperature of open water [C]
  ! fsh - shortwave radiation
  ! flo - longwave radiation
  ! ta - air temperature [C]
  ! qa  - specific humidity             
  ! ug - wind speed [m/sec]
  ! ch - transfer coefficient for sensible heat
  ! ce - transfer coefficient for evaporation
  !
  ! OUTPUT: fh - growth rate
  !         evap - evaporation

  use i_therm_param
  implicit none

  real*8 qa,t,Ta,fsh,flo,ug,ch,ce,fh,evap
  real*8 hfsenow,hfradow,hflatow,hftotow,hflwrdout,b
  real*8 q1, q2 		! coefficients for saturated specific humidity
  real*8 c1, c4, c5
  logical :: standard_saturation_shum_formula = .true.


  !data c1, c4, c5 /3.8e-3, 17.67, 243.5/
  data c1, c4, c5 /3.8e-3, 17.27, 237.3/
  data q1 /640380./, q2 /-5107.4/

  ! (saturated) surface specific humidity
  if(standard_saturation_shum_formula) then
     b=c1*exp(c4*t/(t+c5))                      ! a standard one
  else
     b=0.98*q1*inv_rhoair*exp(q2/(t+tmelt)) 	! LY2004 NCAR version 
  end if

  ! heat fluxes [W/m**2]:

  hfradow= (1.0-albw)*fsh &	                ! absorbed short wave radiation
       +flo             	                ! long wave radiation coming in !put emiss/check
  hflwrdout=-emiss_wat*boltzmann*((t+tmelt)**4) ! long wave radiation going out !in LY2004 emiss=1
  hfradow=hfradow+hflwrdout

  hfsenow=rhoair*cpair*ch*ug*(ta-t)             ! sensible heat 
  evap=rhoair*ce*ug*(qa-b)                      ! evaporation kg/m2/s
  hflatow=clhw*evap                             ! latent heat W/m2

  hftotow=hfradow+hfsenow+hflatow               ! total heat W/m2

  fh= -hftotow/cl                             	! growth rate [m ice/sec]
  !                                           	+: ML gains energy, ice melts
  !                                           	-: ML loses energy, ice grows
  evap=evap*inv_rhowat 	 			! evaporation rate [m water/s],negative up !!!

  return
end subroutine obudget
!
!======================================================================================
!
subroutine flooding (h,hsn)
  use i_therm_param

  real*8 h,hsn,hdraft,hflood

  hdraft=(rhosno*hsn+h*rhoice)*inv_rhowat ! Archimedes: displaced water
  hflood=hdraft-min(hdraft,h)         ! Increase in mean ice thickness due to flooding
  h=h+hflood                          ! Add converted snow to ice volume
  hsn=hsn-hflood*rhoice*inv_rhosno        ! Subtract snow from snow layer

  !RT   This is what all AWI sea ice models do, but
  !RT   I wonder whether it really is correct for the heat budget.
  !RT   I suggest we initially keep it to allow for a comparison with BRIOS results
  !RT   and rethink it at a later stage.

  return
end subroutine flooding
!
!======================================================================================
!
function TFrez(S)
  ! Nonlinear correlation for the water freezing temperature.
  ! Millero (1978) - UNESCO. Reference - See A. Gill, 1982.
  implicit none
  real*8 S, TFrez

  TFrez= -0.0575*S+1.7105e-3 *sqrt(S**3)-2.155e-4 *S*S

end function TFrez
