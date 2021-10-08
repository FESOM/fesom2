subroutine REcoM_sms(n,Nn,state,thick,recipthick,SurfSR,sms,Temp,SinkVel,zF,PAR, mesh)

    use REcoM_declarations
    use REcoM_LocVar
    use REcoM_GloVar
    use recom_config
    use REcoM_ciso
    use g_clock
    use o_PARAM
    use g_PARSUP
    use g_rotate_grid
    use g_config
    use mod_MESH
    use i_arrays
    use o_param
    use i_param
    use o_arrays
    use g_forcing_arrays
    use g_comm_auto
    use i_therm_param
    use g_comm
    use g_support

    implicit none
    type(t_mesh), intent(in) , target                       :: mesh
    integer, intent(in)                                     :: Nn                   !< Total number of nodes in the vertical
    real(kind=8),dimension(mesh%nl-1,bgc_num),intent(inout) :: state                !< ChlA conc in phytoplankton [mg/m3]
									  	    !! should be in instead of inout

    real(kind=8),dimension(mesh%nl-1)                       :: thick                !< [m] Vertical distance between two nodes = Thickness 
    real(kind=8),dimension(mesh%nl-1)                       :: recipthick           !< [1/m] reciprocal of thick
    real(kind=8), intent(in)                                :: SurfSR               !< [W/m2] ShortWave radiation at surface

    real(kind=8),dimension(mesh%nl-1,bgc_num),intent(inout) :: sms                  !< Source-Minus-Sinks term
    real(kind=8),dimension(mesh%nl-1)        ,intent(in)    :: Temp                 !< [degrees C] Ocean temperature
    real(kind=8),dimension(mesh%nl,4)        ,intent(in)    :: SinkVel

    real(kind=8),dimension(mesh%nl),intent(in)              :: zF                   !< [m] Depth of fluxes
    real(kind=8),dimension(mesh%nl-1),intent(inout)         :: PAR

    real(kind=8)                                            :: net                  

    real(kind=8)                                            :: dt_d                 !< Size of time steps [day]
    real(kind=8)                                            :: dt_b                 !< Size of time steps [day]
    real(kind=8),dimension(mesh%nl-1)                       :: Sink
    real(kind=8)                                            :: dt_sink              !< Size of local time step
    real(kind=8)                                            :: Fc                   !< Flux of labile C into sediment, used for denitrification calculation [umolC/cm2/s]
    real(kind=8)                                            :: recip_hetN_plus      !< MB's addition to heterotrophic respiration
    real(kind=8)                                            :: recip_res_het        !< [day] Reciprocal of respiration by heterotrophs and mortality (loss to detritus)

    integer                                                 :: k,step,ii, idiags,n, aux
    real(kind=8)                                            :: & 
        DIN,     & !< Dissolved Inorganic Nitrogen 				[mmol/m3] 
        DIC,     & !< Dissolved Inorganic Carbon				[mmol/m3]
        Alk,     & !< Total Alkalinity					        [mmol/m3]
        PhyN,    & !< Intracellular conc of Nitrogen in small phytoplankton	[mmol/m3]
        PhyC,    & !< Intracellular conc of Carbon in small phytoplankton 	[mmol/m3]
        PhyChl,  & !< Current intracellular ChlA conc. 			        [mg/m3]
        DetN,    & !< Conc of N in Detritus 				        [mmol/m3]
        DetC,    & !< Conc of C in Detritus					[mmol/m3]
        HetN,    & !< Conc of N in heterotrophs				        [mmol/m3]
        HetC,    & !< Conc of C in heterotrophs				        [mmol/m3]
        DON,     & !< Dissolved organic N in the water			        [mmol/m3]
        EOC,     & !< Extracellular Organic C conc				[mmol/m3]
        DiaN,    &
        DiaC,    &
        DiaChl,  &
        DiaSi,   &
        DetSi,   &
        Si,      &
        Fe,      &
        PhyCalc, &
        DetCalc, &
!    if (REcoM_Second_Zoo) then
        Zoo2N,    &
        Zoo2C,    &
        DetZ2N,   &
        DetZ2C,   &
        DetZ2Si,  &
        DetZ2Calc,&
!    endif
        FreeFe,  &
        O2
     
#include "../associate_mesh.h"

    if (ciso) then
        Cphot_z = zero                 ! initialize vertical profiles of
        Cphot_dia_z = zero             ! daily-mean photosynthesis rates
        if (.not. ciso_calcdiss) then  ! turn off isotopic fractionation
            alpha_calc_13 = 1.d0   ! during calcification / dissolution
            alpha_calc_14 = 1.d0
            alpha_dcal_13 = 1.d0
            alpha_dcal_14 = 1.d0
        endif
    endif

    sms = zero ! double precision

    tiny_N   = tiny_chl/chl2N_max      ! 0.00001/ 3.15d0   Chl2N_max [mg CHL/mmol N] Maximum CHL a : N ratio = 0.3 gCHL gN^-1
    tiny_N_d = tiny_chl/chl2N_max_d    ! 0.00001/ 4.2d0
    tiny_C   = tiny_N  /NCmax          ! NCmax   = 0.2d0   [mmol N/mmol C] Maximum cell quota of nitrogen (N:C)
    tiny_C_d = tiny_N_d/NCmax_d        ! NCmax_d = 0.2d0 
    tiny_Si  = tiny_C_d/SiCmax         ! SiCmax = 0.8d0

    recip_res_het = 1.d0/res_het       ! res_het = 0.01d0  [1/day] Respiration by heterotrophs and mortality (loss to detritus)

!-------------------------------------------------------------------------------
!> REcoM time steps [day]
!-------------------------------------------------------------------------------

    rTref =  real(one)/recom_Tref
  
    dt_d  =  dt/SecondsPerDay     !< Size of FESOM time step [day]
    dt_b  =  dt_d/real(biostep)   !< Size of REcoM time step [day]

!-------------------------------------------------------------------------------
!Main time loop starts
    do step  = one,biostep

        kdzUpper	= 0.d0	        !< Upper light attenuation of top cell is set to zero

        if (any(abs(sms(:,:)) <= tiny)) sms(:,:) = zero      ! tiny = 2.23D-16
!if (mype==0) write(*,*) sms(:,:)

!-------------------------------------------------------------------------------
! Main vertical loop starts
        do k = one,Nn   ! nzmin, nzmax 
!    do n=1, myDim_nod2D!+eDim_nod2D 
!       Nn=nlevels_nod2D(n)-1  !nzmax
!       nzmin = ulevels_nod2D(row)
!       nzmax = nlevels_nod2D(row)
!       do k=1, Nn
            DIN    = max(tiny,state(k,idin)  	 	+ sms(k,idin  )) !< Avoids division by zero
            DIC    = max(tiny,state(k,idic)   		+ sms(k,idic  )) !! and updates Conc between
            ALK    = max(tiny,state(k,ialk)   		+ sms(k,ialk  )) !! local steps in REcoM when
            PhyN   = max(tiny_N,state(k,iphyn)  	+ sms(k,iphyn )) !! biostep > 1
            PhyC   = max(tiny_C,state(k,iphyc) 		+ sms(k,iphyc ))
            PhyChl = max(tiny_chl,state(k,ipchl)  	+ sms(k,ipchl ))
            DetN   = max(tiny,state(k,idetn)  		+ sms(k,idetn ))
            DetC   = max(tiny,state(k,idetc)  		+ sms(k,idetc ))
            HetN   = max(tiny,state(k,ihetn)  		+ sms(k,ihetn ))
            HetC   = max(tiny,state(k,ihetc)  		+ sms(k,ihetc ))
            if (REcoM_Second_Zoo) then 
                Zoo2N  = max(tiny,state(k,izoo2n)       + sms(k,izoo2n))
                Zoo2C  = max(tiny,state(k,izoo2c)       + sms(k,izoo2c))
                DetZ2N = max(tiny,state(k,idetz2n)      + sms(k,idetz2n))
                DetZ2C = max(tiny,state(k,idetz2c)      + sms(k,idetz2c))
                DetZ2Si = max(tiny,state(k,idetz2si)    + sms(k,idetz2si))
                DetZ2Calc = max(tiny,state(k,idetz2calc)+ sms(k,idetz2calc))

            endif
            DON    = max(tiny,state(k,idon)   		+ sms(k,idon  ))
            EOC    = max(tiny,state(k,idoc)   		+ sms(k,idoc  ))
            DiaN   = max(tiny_N,state(k,idian)  	+ sms(k,idian ))
            DiaC   = max(tiny_C,state(k,idiac)  	+ sms(k,idiac ))
            DiaChl = max(tiny_chl,state(k,idchl)  	+ sms(k,idchl ))
            DiaSi  = max(tiny_si,state(k,idiasi)      	 + sms(k,idiasi))
            !DiaSi  = max(tiny,state(k,idiasi) 	+ sms(k,idiasi))
            DetSi  = max(tiny,state(k,idetsi) 		 + sms(k,idetsi))
            Si     = max(tiny,state(k,isi)    		 + sms(k,isi   ))
            Fe     = max(tiny,state(k,ife)    		+ sms(k,ife   ))
            O2     = max(tiny,state(k,ioxy)             + sms(k,ioxy  ))
            FreeFe = zero
!! REcoM calcification
            PhyCalc= max(tiny,state(k,iphycal)		+ sms(k,iphycal))
            DetCalc= max(tiny,state(k,idetcal)		+ sms(k,idetcal))

            calc_diss      = calc_diss_rate * SinkVel(k,ivdet) /20.d0 ! Dissolution rate of CaCO3 scaled by the sinking velocity at the current depth 0.005714   !20.d0/3500.d0
            calc_diss2     = calc_diss_rate2  ! Dissolution rate of CaCO3 for seczoo

            quota          =  PhyN / PhyC ! include variability of the N: C ratio, cellular chemical composition 
            recipquota     =  real(one) / quota
            Chl2C          =  PhyChl  / PhyC ! Chl a:phytoplankton carbon ratio, cellular chemical composition [gCHL gC^-1]
            Chl2N          =  PhyChl  / PhyN ! Chl a:phytoplankton nitrogen ratio, cellular chemical composition [gCHL gN^-1]
            CHL2C_plast    =  Chl2C * (quota/(quota - NCmin))
    
            quota_dia      	=  DiaN / DiaC
            recipQuota_dia 	=  real(one)/quota_dia
            Chl2C_dia      	=  DiaChl / DiaC
            Chl2N_dia      	=  DiaChl / DiaN
            CHL2C_plast_dia 	=  Chl2C_dia * (quota_dia/(quota_dia - NCmin_d))
            qSiC           	=  DiaSi / DiaC
            qSiN           	=  DiaSi / DiaN

            recipQZoo      	=  HetC / HetN
            recip_hetN_plus 	= 1. / (HetN + tiny_het) ! MB's addition for more stable zoo respiration

            if (REcoM_Second_Zoo) then  
                recipQZoo2     =  Zoo2C / Zoo2N
            endif
            if (Grazing_detritus) then
                 recipDet  = DetC / DetN
                 recipDet2 = DetZ2C / DetZ2N
            end if

            if (ciso) then
!<       additional variables are declared in module REcoM_ciso
                DIC_13     = max(tiny,state(k,idic_13)    + sms(k,idic_13  ))
                DIC_14     = max(tiny,state(k,idic_14)    + sms(k,idic_14  ))
                PhyC_13    = max(tiny_C,state(k,iphyc_13) + sms(k,iphyc_13 ))
                PhyC_14    = max(tiny_C,state(k,iphyc_14) + sms(k,iphyc_14 ))
                DetC_13    = max(tiny,state(k,idetc_13)   + sms(k,idetc_13 ))
                DetC_14    = max(tiny,state(k,idetc_14)   + sms(k,idetc_14 ))
                HetC_13    = max(tiny,state(k,ihetc_13)   + sms(k,ihetc_13 ))
                HetC_14    = max(tiny,state(k,ihetc_14)   + sms(k,ihetc_14 ))
                EOC_13     = max(tiny,state(k,idoc_13)    + sms(k,idoc_13  ))
                EOC_14     = max(tiny,state(k,idoc_14)    + sms(k,idoc_14  ))
                DiaC_13    = max(tiny_C,state(k,idiac_13) + sms(k,idiac_13 ))
                DiaC_14    = max(tiny_C,state(k,idiac_14) + sms(k,idiac_14 ))
                PhyCalc_13 = max(tiny,state(k,iphycal_13) + sms(k,iphycal_13))
                PhyCalc_14 = max(tiny,state(k,iphycal_14) + sms(k,iphycal_14))
                DetCalc_13 = max(tiny,state(k,idetcal_13) + sms(k,idetcal_13))
                DetCalc_14 = max(tiny,state(k,idetcal_14) + sms(k,idetcal_14))

                calc_diss_13   = alpha_dcal_13 * calc_diss
                calc_diss_14   = alpha_dcal_14 * calc_diss

                quota_13             = PhyN / PhyC_13
                quota_14             = PhyN / PhyC_14
                recipQuota_13        = real(one) / quota_13
                recipQuota_14        = real(one) / quota_14

                quota_dia_13         = DiaN / DiaC_13
                quota_dia_14         = DiaN / DiaC_14
                recipQuota_dia_13    = real(one) / quota_dia_13
                recipQuota_dia_14    = real(one) / quota_dia_14

                recipQZoo_13         = HetC_13 / HetN 
                recipQZoo_14         = HetC_14 / HetN 
            end if ! ciso

!-------------------------------------------------------------------------------
!> Temperature dependence of rates
!------------------------------------------------------------------------------- 
!< Schourup 2013 Eq. A54
!< Temperature dependence of metabolic rate, fT, dimensionless
!< Ae: Slope of the linear region of the Arrhenius plot
!< rTloc: Inverse of local temperature in [1/Kelvin]
!< rTref=288.15 (15 degC): Reference temperature for Arrhenius equation [1/Kelvin]
!< See Figure A1
!< Other functions can be used for temperature dependency (Eppley 1972; Li 1980; Ahlgren 1987)

    rTloc   = real(one)/(Temp(k) + C2K)
    arrFunc = exp(-Ae * ( rTloc - rTref))

    if (REcoM_Second_Zoo) then 
        arrFuncZoo2 = exp(t1_zoo2/t2_zoo2 - t1_zoo2*rTloc)/(1 + exp(t3_zoo2/t4_zoo2 - t3_zoo2*rTloc))
    endif

!< Silicate temperature dependence 
!    reminSiT = min(1.32e16 * exp(-11200.d0 * rTloc),reminSi) !! arrFunc control, reminSi=0.02d0 ! Kamatani (1982)
    reminSiT = max(0.023d0 * 2.6d0**((Temp(k)-10.)/10.),reminSi)
!    reminSiT = reminSi

!-------------------------------------------------------------------------------
!> Photosynthesis section, light parameters and rates
!-------------------------------------------------------------------------------

!_______________________________________________________________________
!< Small phytoplankton
!< Schourup 2013 Appendix A6.2
!< Intracellular regulation of C uptake
!< qlimitFac, qlimitFacTmp: Factor that regulates photosynthesis
!< NMinSlope: 50.d0
!< NCmin: 0.04d0
!< quota: PhyN/PhyC
!< qlimitFac [0.0, 1.0] 
!< if quota < NCmin qlimitFac=0
!< if quota > ≈ 9 * NCmin qlimitFac=1
!< P_cm: 3.0d0 [1/day], Rate of C-specific photosynthesis 


!< pMax = The carbon-specific, light-saturated rate of photosynthesis [day^-1]
!< Nutrient limited environment
!< Small pyhtoplankton is limited by iron and nitrogen
!< Diatoms are additionally limited by silicon

    qlimitFac     	= recom_limiter(NMinSlope,NCmin,quota) ! Eqn A55
    feLimitFac  	= Fe/(k_Fe + Fe) ! Use Michaelis–Menten kinetics
    qlimitFac   	= min(qlimitFac,feLimitFac) ! Liebig law of the minimum
    pMax          	= P_cm * qlimitFac * arrFunc ! Maximum value of C-specific rate of photosynthesis
    
!_______________________________________________________________________
!< Diatoms
    qlimitFac     	= recom_limiter(NMinSlope,NCmin_d,quota_dia)
    qlimitFacTmp  	= recom_limiter(SiMinSlope,SiCmin,qSiC)
    qlimitFac     	= min(qLimitFac,qlimitFacTmp)
    feLimitFac  	= Fe/(k_Fe_d + Fe)
    qlimitFac   	= min(qlimitFac,feLimitFac)
    pMax_dia      	= P_cm_d * qlimitFac * arrFunc

!_______________________________________________________________________
!< Light
! Exponential diminition of downward irradiance
    if (k==1) then  ! ulevels_nod2D(n)==1
        PARave = max(tiny,SurfSR)
        PAR(k) = PARave
        chl_upper = (PhyChl + DiaChl)
    else    
        chl_lower = PhyChl + DiaChl
        Chlave    = (chl_upper+chl_lower)*0.5d0

        kappar         =  k_w + a_chl * (Chlave)
        kappastar      =  kappar / cosAI(n)
        kdzLower       =  kdzUpper + kappastar * thick(k-1)
        Lowerlight     =  SurfSR * exp(-kdzLower)
        Lowerlight     =  max(tiny,Lowerlight)
        PARave         =  Lowerlight
        PAR(k)         =  PARave
        chl_upper      =  chl_lower
        kdzUpper       =  kdzLower
    end if

!-------------------------------------------------------------------------------
!< Small phytoplankton photosynthesis rate
    if ( pMax .lt. tiny .OR. PARave /= PARave                  &
         .OR. CHL2C /= CHL2C) then
        Cphot       = zero
    else
        Cphot       = pMax*( real(one) &
                     - exp(-alfa * Chl2C * PARave / pMax))
    end if
    if ( Cphot .lt. tiny) Cphot = zero
    
!------------------------------------------------------------------------------
!< Diatom photosynthesis rate
    if ( pMax_dia .lt. tiny .OR. PARave /= PARave               &
         .OR. CHL2C_dia /= CHL2C_dia) then
        Cphot_dia   = zero
    else
        Cphot_dia   = pMax_dia * (real(one) &
                     - exp(-alfa_d * Chl2C_dia * PARave / pMax_dia))
    end if
    if (Cphot_dia .lt. tiny) Cphot_dia = zero

    if (ciso) then
!<       Photosynthesis rates involving 13|14C should be smaller than
!!       C_phot,_dia but they are poorly constrained for photoperiods
!!       shorter than a few hours. For this reason we disregard rate
!!       differences but treat the isotopic fractionation associated with
!!       photosynthesis as a pseudo equilibrium process at the end of time
!!       stepping (cf subroutine recom_forcing). To do so, we need to store
!!       interim values of photosynthesis rates at depth:
        Cphot_z(k)     = Cphot_z(k)     + Cphot
        Cphot_dia_z(k) = Cphot_dia_z(k) + Cphot_dia
    end if

!-------------------------------------------------------------------------------- 
!< chlorophyll degradation
    KOchl = deg_Chl
    KOchl_dia = deg_Chl_d
        
    if (use_photodamage) then
!< Phyto Chla loss
!< adding a minimum value for photodamage 
        if (pMax .lt. tiny .OR. PARave /= PARave                 &
             .OR. CHL2C_plast /= CHL2C_plast) then
            KOchl = deg_Chl*0.1d0
        else
            KOchl = deg_Chl*(1-exp(-alfa * CHL2C_plast * PARave / pMax))
            KOchl = max((deg_Chl*0.1d0), KOchl)
        end if

!< Phyto Chla loss                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
        if (pMax_dia .lt. tiny .OR. PARave /= PARave             &
            .OR. CHL2C_plast_dia /= CHL2C_plast_dia) then
            KOchl_dia = deg_Chl_d*0.1d0
        else
            KOchl_dia = deg_Chl_d * ( 1 -                         &
               exp( -alfa_d * CHL2C_plast_dia * PARave / pMax_dia ))
            KOchl_dia = max((deg_Chl_d*0.1d0), KOchl_dia)
        end if

        if (KOchl /= KOchl) then
            print*,' KOchl is ', KOchl
            print*,' deg_Chl is ', deg_Chl
            print*,' alfa is ', alfa
            print*,' CHL2C is ', CHL2C_plast
            print*,' PARave is ', PARave
            print*,' pMax is ', pMax
            stop
        end if
        if (KOchl_dia /= KOchl_dia) then
            print*,' KOchl_dia is ', KOchl_dia
            print*,' deg_Chl_d is ', deg_Chl_d
            print*,' alfa_d is ', alfa_d
            print*,' CHL2C_d is ', CHL2C_plast_dia
            print*,' PARave is ', PARave
            print*,' pMax_d is ', pMax_dia
            stop
        end if

    end if ! photodamage 
 
!-------------------------------------------------------------------------------
!> Assimilation section
!-------------------------------------------------------------------------------
!< Nitrogen and silicon part
!< Compute assimilation from Geider et al 1998
!< V_cm: Scaling factor for C-specific N uptake, dimensionless
!< NCmax: Maximum cell quota of nitrogen (N:C) [mmol N/mmol C]
!< NMaxSlope: Max slope for limiting function
!< NCuptakeRatio: Maximum uptake ratio N:C [mmol N mmol C−1]
!< SiCUptakeRatio: Maximum uptake ratio Si : C [mmol Si mmol C−1 ]
!< The N:C ratio is taken into account, as a
!! too high ratio indicates that the intracellular 
!! concentration of energy rich carbon molecules becomes too low to 
!! use energy on silicon uptake.

    V_cm           = V_cm_fact
    limitFacN      = recom_limiter(NMaxSlope,quota,NCmax)
    N_assim        = V_cm * pMax * NCuptakeRatio &                ! [mmol N / (mmol C * day)]
                      * limitFacN * (DIN/(DIN + k_din)) ! Michaelis–Menten kinetics

    V_cm           = V_cm_fact_d
    limitFacN_dia  = recom_limiter(NMaxSlope,quota_dia,NCmax_d)
    N_assim_dia    = V_cm * pMax_dia * NCUptakeRatio_d &
                      * limitFacN_dia * DIN/(DIN + k_din_d)

    limitFacSi     = recom_limiter(SiMaxSlope,qSiC,SiCmax)  &
                      * limitFacN_dia
    Si_assim       = V_cm * P_cm_d * arrFunc * SiCUptakeRatio &
                      * limitFacSi * Si/(Si + k_si)

!-------------------------------------------------------------------------------
!< Iron chemistry
 
    freeFe = iron_chemistry(Fe,totalligand,ligandStabConst)

!-------------------------------------------------------------------------------
!< Chlorophyll synthesis
!< Coupled to N uptake
!< Converted to chlorophyll units with a maximum Chl:N ratio, Chl2N_max
!< Chl2N_max: Maximum Chl:N ratio for phytoplankton [mg Chl mmol N−1 ]

    chlSynth = zero
    if (PARave .ge. tiny .AND. PARave .eq. PARave) then
        chlSynth = N_assim * Chl2N_max                    &
            * min( real(one),Cphot/(alfa * Chl2C * PARave))
    end if
    ChlSynth_dia = zero
    if (PARave .ge. tiny .AND. PARave .eq. PARave) then
        ChlSynth_dia = N_assim_dia                        &
            * Chl2N_max_d * min(real(one),                 &
            Cphot_dia /(alfa_d * Chl2C_dia * PARave))
    end if
!-------------------------------------------------------------------------------
!< Phytoplankton respiraion rate
!< res_phy: Maintenance respiration rate constant [day−1 ]
!< biosynth: The cost of biosynthesis of N [mmol C mmol N−1 ]
    phyRespRate     = res_phy * limitFacN + biosynth * N_assim
    phyRespRate_dia = res_phy_d * limitFacN_dia +        &
        biosynth * N_assim_dia + biosynthSi * Si_assim

    if (ciso) then
!       we assume that
!       phyRespRate_13,14     = phyRespRate
!       phyRespRate_dia_13,14 = phyRespRate_dia
    end if

!-------------------------------------------------------------------------------
!< Zooplankton grazing on small phytoplankton and diatoms
!< At the moment there is no preference for one or the other food. Change this!

    if (REcoM_Grazing_Variable_Preference) then
        if (Grazing_detritus) then
            if (Graz_pref_new) then

!< pzPhy: Maximum nanophytoplankton preference
!< pzDia: Maximum diatom preference
!< pzDet: Maximum small detritus prefence by first zooplankton
!< pzDetZ2: Maximum large detritus preference by first zooplankton

                varpzPhy    = pzPhy   * PhyN /(pzPhy*PhyN + pzDia*DiaN + PzDet*DetN + pzDetZ2*DetZ2N)
                varpzDia    = pzDia   * DiaN /(pzPhy*PhyN + pzDia*DiaN + PzDet*DetN + pzDetZ2*DetZ2N)
                varpzDet    = pzDet   * DetN /(pzPhy*PhyN + pzDia*DiaN + PzDet*DetN + pzDetZ2*DetZ2N)
                varpzDetZ2  = pzDetZ2 * DetN /(pzPhy*PhyN + pzDia*DiaN + PzDet*DetN + pzDetZ2*DetZ2N)
!                aux         = pzPhy*PhyN + pzDia*DiaN + PzDet*DetN + pzDetZ2*DetZ2N
!                varpzPhy    = pzPhy   * PhyN /aux
!                varpzDia    = pzDia   * DiaN /aux
!                varpzDet    = pzDet   * DetN /aux
!                varpzDetZ2  = pzDetZ2 * DetN /aux
            else
!                DiaNsq      = DiaN**2
!                PhyNsq      = PhyN**2
!                DetNsq      = DetN**2
!                DetZ2Nsq    = DetZ2N**2
                PhyNsq      = PhyN * PhyN
                varpzPhy    = pzPhy * PhyNsq /(sPhyNsq + PhyNsq)
                DiaNsq      = DiaN * DiaN
                varpzDia    = pzDia * DiaNsq /(sDiaNsq + DiaNsq)
                DetNsq      = DetN * DetN
                varpzDet    = pzDet * DetNsq /(sDetNsq + DetNsq)
                DetZ2Nsq    = DetZ2N * DetZ2N
                varpzDetZ2  = pzDetZ2 * DetZ2Nsq /(sDetZ2Nsq + DetZ2Nsq)

            endif

            fDiaN    = varpzDia   * DiaN
            fPhyN    = varpzPhy   * PhyN
            fDetN    = varpzDet   * DetN
            fDetZ2N  = varpzDetZ2 * DetZ2N

        else

            if (Graz_pref_new) then
                varpzPhy      = pzPhy * PhyN /(pzPhy*PhyN + pzDia*DiaN)
                varpzDia      = pzDia * DiaN /(pzPhy*PhyN + pzDia*DiaN)
            else
                DiaNsq        = DiaN  * DiaN
                varpzDia      = pzDia * DiaNsq /(sDiaNsq + DiaNsq)
                PhyNsq        = PhyN  * PhyN
                varpzPhy      = pzPhy * PhyNsq /(sPhyNsq + PhyNsq)
            end if

            fDiaN         = varpzDia * DiaN
            fPhyN         = varpzPhy * PhyN

        end if
    else
        fDiaN         = pzDia * DiaN
        fPhyN         = pzPhy * PhyN
        if (Grazing_detritus) then
            fDetN        = pzDet   * DetN
            fDetZ2N      = pzDetZ2 * DetZ2N
        end if
    end if

    if (Grazing_detritus) then
       food              = fPhyN + fDiaN + fDetN + fDetZ2N
       foodsq            = food * food
       grazingFlux       = (Graz_max * foodsq)/(epsilonr + foodsq) * HetN * arrFunc
       grazingFlux_phy   = grazingFlux * fphyN / food
       grazingFlux_Dia   = grazingFlux * fDiaN / food
       grazingFlux_Det   = grazingFlux * fDetN / food
       grazingFlux_DetZ2 = grazingFlux * fDetZ2N / food
    else
       food              = PhyN + fDiaN
       foodsq            = food * food
       grazingFlux       = (Graz_max * foodsq)/(epsilonr + foodsq) * HetN * arrFunc
       grazingFlux_phy   = grazingFlux * phyN / food
       grazingFlux_Dia   = grazingFlux * fDiaN / food
    endif

    if (REcoM_Second_Zoo) then
!-------------------------------------------------------------------------------
!< Second Zooplankton grazing on small phytoplankton, diatoms and heterotrophs
!< At the moment there is no preference for one or the other food. Change this!
        
    if (REcoM_Grazing_Variable_Preference) then
       if (Grazing_detritus) then
          if (Graz_pref_new) then
             varpzDia2      = pzDia2   * DiaN   /(pzPhy2 * PhyN + PzDia2 * DiaN + pzHet * HetN + pzDet2 * DetN + pzDetZ22 * DetZ2N)
             varpzPhy2      = pzPhy2   * PhyN   /(pzPhy2 * PhyN + PzDia2 * DiaN + pzHet * HetN + pzDet2 * DetN + pzDetZ22 * DetZ2N)
             varpzHet       = pzHet    * HetN   /(pzPhy2 * PhyN + PzDia2 * DiaN + pzHet * HetN + pzDet2 * DetN + pzDetZ22 * DetZ2N)  
             varpzDet2      = pzDet2   * DetN   /(pzPhy2 * PhyN + PzDia2 * DiaN + pzHet * HetN + pzDet2 * DetN + pzDetZ22 * DetZ2N)
             varpzDetZ22    = pzDetZ22 * DetZ2N /(pzPhy2 * PhyN + PzDia2 * DiaN + pzHet * HetN + pzDet2 * DetN + pzDetZ22 * DetZ2N)
!             aux            = pzPhy2 * PhyN + PzDia2 * DiaN + pzHet * HetN + pzDet2 * DetN + pzDetZ22 * DetZ2N
!             varpzDia2      = pzDia2   * DiaN   /aux
!             varpzPhy2      = pzPhy2   * PhyN   /aux
!             varpzHet       = pzHet    * HetN   /aux  
!             varpzDet2      = pzDet2   * DetN   /aux
!             varpzDetZ22    = pzDetZ22 * DetZ2N /aux
          else
             DiaNsq2        = DiaN * DiaN
             varpzDia2      = pzDia2 * DiaNsq2 /(sDiaNsq2 + DiaNsq2)
             fDiaN2         = varpzDia2 * DiaN
             PhyNsq2        = PhyN * PhyN
             varpzPhy2      = pzPhy2 * PhyNsq2 /(sPhyNsq2 + PhyNsq2)
             fPhyN2         = varpzPhy2 * PhyN
             HetNsq         = HetN * HetN
             varpzHet       = pzHet * HetNsq /(sHetNsq + HetNsq)
             fHetN          = varpzHet * HetN   
             DetNsq         = DetN * DetN
             varpzDet2      = pzDet2 * DetNsq /(sDetNsq2 + DetNsq)
             DetZ2Nsq       = DetZ2N * DetZ2N
             varpzDetZ22    = pzDetZ22 * DetZ2Nsq /(sDetZ2Nsq2 + DetZ2Nsq)
          end if
          fDiaN2         = varpzDia2 * DiaN
          fPhyN2         = varpzPhy2 * PhyN
          fHetN          = varpzHet * HetN
          fDetN2         = varpzDet2 * DetN
          fDetZ2N2       = varpzDetZ22 * DetZ2N
       else
          if (Graz_pref_new) then
             varpzDia2      = pzDia2 * DiaN /(pzPhy2 * PhyN + PzDia2 * DiaN + pzHet * HetN)
             varpzPhy2      = pzPhy2 * PhyN /(pzPhy2 * PhyN + PzDia2 * DiaN + pzHet * HetN)
             varpzHet       = pzHet * HetN /(pzPhy2 * PhyN + PzDia2 * DiaN + pzHet * HetN)
!             aux            = pzPhy2 * PhyN + PzDia2 * DiaN + pzHet * HetN
!             varpzDia2      = pzDia2 * DiaN /aux
!             varpzPhy2      = pzPhy2 * PhyN /aux
!             varpzHet       = pzHet  * HetN /aux

          else
             DiaNsq2        = DiaN * DiaN
             varpzDia2      = pzDia2 * DiaNsq2 /(sDiaNsq2 + DiaNsq2)
             fDiaN2         = varpzDia2 * DiaN
             PhyNsq2        = PhyN * PhyN
             varpzPhy2      = pzPhy2 * PhyNsq2 /(sPhyNsq2 + PhyNsq2)
             fPhyN2         = varpzPhy2 * PhyN
             HetNsq         = HetN * HetN
             varpzHet       = pzHet * HetNsq /(sHetNsq + HetNsq)
          end if
          fDiaN2         = varpzDia2 * DiaN
          fPhyN2         = varpzPhy2 * PhyN
          fHetN          = varpzHet * HetN
       end if
    else
       fDiaN2         = pzDia2 * DiaN
       fPhyN2         = pzPhy2 * PhyN
       fHetN          = pzHet * HetN
       if (Grazing_detritus) then
          fDetN2         = pzDet2 * DetN
          fDetZ2N2       = pzDetZ22 * DetZ2N
       end if
    end if

    if (Grazing_detritus) then
       food2             = fPhyN2 + fDiaN2 + fHetN + fDetN2 + fDetZ2N2
       foodsq2           = food2 * food2
       grazingFlux2     = (Graz_max2 * foodsq2)/(epsilon2 + foodsq2) * Zoo2N * arrFuncZoo2
       grazingFlux_phy2  = grazingFlux2 * fphyN2 / food2
       grazingFlux_Dia2  = grazingFlux2 * fDiaN2 / food2
       grazingFlux_het2  = grazingFlux2 * fHetN / food2
       grazingFlux_Det2  = grazingFlux2 * fDetN2 / food2
       grazingFlux_DetZ22  = grazingFlux2 * fDetZ2N2 / food2

       grazingFluxcarbonzoo2 = (grazingFlux_phy2 * recipQuota * grazEff2) &
                         + (grazingFlux_Dia2 * recipQuota_Dia * grazEff2) &
                         + (grazingFlux_het2 * recipQZoo * grazEff2)      &
                         + (grazingFlux_Det2 * recipDet * grazEff2)       &
                         + (grazingFlux_DetZ22 *recipDet2 * grazEff2)    
    else
       food2             = fPhyN2 + fDiaN2 + fHetN
       foodsq2           = food2 * food2
       grazingFlux2      = (Graz_max2 * foodsq2)/(epsilon2 + foodsq2) * Zoo2N * arrFuncZoo2
       grazingFlux_phy2  = grazingFlux2 * fphyN2 / food2
       grazingFlux_Dia2  = grazingFlux2 * fDiaN2 / food2
       grazingFlux_het2  = grazingFlux2 * fHetN / food2

       grazingFluxcarbonzoo2 = (grazingFlux_phy2 * recipQuota * grazEff2) &
                          + (grazingFlux_Dia2 * recipQuota_Dia * grazEff2) &
                          + (grazingFlux_het2 * recipQZoo * grazEff2)
    end if
    end if

!-------------------------------------------------------------------------------
! Heterotrophic respiration is assumed to drive zooplankton back to Redfield C:N
! if their C:N becomes higher than Redfield
! res_het: Timescale for zooplankton respiration [day−1 ]

    if (het_resp_noredfield) then
        HetRespFlux    = res_het *  arrFunc * HetC ! tau * f_T [HetC]
    else
        if (HetRespFlux_plus) then
            HetRespFlux    = recip_res_het * arrFunc * (hetC    * recip_hetN_plus - redfield) * HetC
        else
! default computation scheme 
            HetRespFlux    = recip_res_het * arrFunc * (recipQZoo    - redfield) * HetC  ! 1/tau * f_T * (q_zoo - q_standard)
        endif
        HetRespFlux    = max(zero,HetRespFlux)
    endif

! Next part changes results, but is needed: Otherwise heterotrophs take up 
! inorganic C when their C:N becomes lower than Redfield

!    HetRespFlux    = max(zero,HetRespFlux)

    if (ciso) then
!MB    set HetRespFlux_plus = .true. !
        HetRespFlux_13 = max(zero, recip_res_het * arrFunc * (hetC_13 * recip_hetN_plus - redfield) * HetC_13)
        HetRespFlux_14 = max(zero, recip_res_het * arrFunc * (hetC_14 * recip_hetN_plus - redfield) * HetC_14)
    end if

!-------------------------------------------------------------------------------
! Quadratic zooplanton mortality
    hetLossFlux    = loss_het * HetN * HetN

    if (REcoM_Second_Zoo) then
!-------------------------------------------------------------------------------
! if (REcoM_Second_Zoo .eq. .true.) then
! Second zooplankton respiration is assumed to drive zooplankton back to Redfield C:N
! if their C:N becomes higher than Redfield
       call krill_resp(n,mesh)
   !print*,'daynew: ', daynew
   !print*,'Latr ', Latr 
       if((grazingFluxcarbonzoo2/Zoo2C) <= 0.1)then
          res_zoo2_f = 0.1*(grazingFluxcarbonzoo2/Zoo2C*100)
       else
          res_zoo2_f = 1.
       end if  
       recip_res_zoo22 = res_zoo2*(1.+res_zoo2_f + res_zoo2_a)
       Zoo2RespFlux    = recip_res_zoo22 * Zoo2C                   
!-------------------------------------------------------------------------------
! Quadratic second zooplanton mortality
       Zoo2LossFlux    = loss_zoo2 * zoo2N * zoo2N

       if(zoo2_fecal_loss) then
          Zoo2fecalloss_n   = fecal_rate_n * grazingFlux2
          Zoo2fecalloss_c   = fecal_rate_c * grazingFluxcarbonzoo2
       else
          Zoo2fecalloss_n   = 0.0
          Zoo2fecalloss_c   = 0.0
       end if
    end if
!-------------------------------------------------------------------------------
! Phytoplankton and detritus aggregation
    if (diatom_mucus) then
       qlimitFac     = recom_limiter(NMinSlope,NCmin_d,quota_dia)
       qlimitFacTmp  = recom_limiter(SiMinSlope,SiCmin,qSiC)
       qlimitFac     = min(qLimitFac,qlimitFacTmp)
       feLimitFac  = Fe/(k_Fe_d + Fe)
       qlimitFac   = min(qlimitFac,feLimitFac)
       if (REcoM_Second_Zoo) then 
          aggregationrate = agg_PD * DetN + agg_PD * DetZ2N &
                             + agg_PP * PhyN + agg_PP * (1 - qlimitFac) * DiaN
       else
          aggregationrate = agg_PD * DetN + agg_PP * PhyN &
                             + agg_PP * (1 - qlimitFac) * DiaN
       endif
    else
       if (REcoM_Second_Zoo) then
          aggregationrate = agg_PD * DetN + agg_PD * DetZ2N &
                             + agg_PP * PhyN + agg_PP * DiaN
       else
          aggregationrate = agg_PD * DetN + agg_PP * PhyN &
                             + agg_PP * DiaN
       endif
    endif
!-------------------------------------------------------------------------------
! Phytoplankton and detritus aggregation
!    aggregationrate = agg_PD * DetN + agg_PP * PhyN &
!                    + agg_PP * DiaN
!-------------------------------------------------------------------------------
! Terms required for the formation and dissolution of CaCO3
!#ifdef REcoM_calcification
!< calc_prod_ratio: Calcite production ratio, dimensionless
    calcification = calc_prod_ratio * Cphot * PhyC   ! Z in equations
    calc_loss_agg = aggregationRate * PhyCalc

    if(REcoM_Second_Zoo)  then
!     calc_loss_gra  = (grazingFlux_phy + grazingFlux_phy2)  &
!                     * recipQuota/(PhyC + tiny)    * PhyCalc
       calc_loss_gra  =  grazingFlux_phy   &
                       * recipQuota/(PhyC + tiny)    * PhyCalc
       calc_loss_gra2 = grazingFlux_phy2           &
                       * recipQuota/(PhyC + tiny)    * PhyCalc
    else
       calc_loss_gra = grazingFlux_phy           &
                       * recipQuota/(PhyC + tiny)    * PhyCalc
    endif
 
 
    if (ciso) then
        calcification_13 = calc_prod_ratio * Cphot * PhyC_13 * alpha_calc_13
        calcification_14 = calc_prod_ratio * Cphot * PhyC_14 * alpha_calc_14
        calc_loss_agg_13 = aggregationRate * PhyCalc_13
        calc_loss_agg_14 = aggregationRate * PhyCalc_14
        calc_loss_gra_13 = grazingFlux_phy *              &
            recipQuota_13/(PhyC_13 + tiny)   * PhyCalc_13
        calc_loss_gra_14 = grazingFlux_phy *              &
            recipQuota_14/(PhyC_14 + tiny)   * PhyCalc_14
    end if

!#endif
		
!-------------------------------------------------------------------------------
! Sources minus sinks are calculated
!-------------------------------------------------------------------------------

!____________________________________________________________
!< DIN 

!< DON:            Extracellular dissolved organic nitrogen [mmolN m-3]
!< rho_N*arrFunc : Remineralization rate, temperature dependency is calculated with arrFunc [day^-1]
!< DiaC:           Intracellular carbon concentration in diatoms [mmolC m-3]
!< PhyC:           Intracellular carbon concentration in nanophytoplankton [mmolC m-3]
!< N_assim:        N assimilation rate for nanophytoplankton [mmolN mmolC-1 day-1]
!< N_assim_Dia:    N assimilation rate for diatoms [mmolN mmolC-1 day-1]
!< dt_b:           REcoM time step [day]

!! Schourup 2013 Eq. A2

    sms(k,idin)      = (                       &
        - N_assim                      * PhyC    &  ! --> N assimilation Nanophytoplankton, [mmol N/(mmol C * day)] C specific N utilization rate
        - N_assim_Dia                  * DiaC    &  ! --> N assimilation Diatoms
        + rho_N * arrFunc              * DON     &  ! --> DON remineralization, temperature dependent [day^-1 * mmol/m3]       
!        + LocRiverDIN                            & ! --> added in FESOM2 
                                             ) * dt_b + sms(k,idin)  

!____________________________________________________________
!< DIC

!< rho_C1: Temperature dependent C degradation of extracellular organic C (EOC) [day^-1]

    if(REcoM_Second_Zoo)  then
        sms(k,idic)      = (                             &
            - Cphot                         * PhyC       & ! --> Small pyhtoplankton photosynthesis 
            + phyRespRate                   * PhyC       & ! --> Small pyhtoplankton respiration 
            - Cphot_Dia                     * DiaC       & ! --> Diatom photosynthesis 
            + phyRespRate_Dia               * DiaC       & ! --> Diatom respiration
            + rho_C1 * arrFunc              * EOC        & ! --> Remineralization of DOC
            + HetRespFlux                                & ! --> Zooplankton respiration                     
            + Zoo2RespFlux                               &                    
            + calc_diss                     * DetCalc    & ! --> Calcite dissolution from detritus 
            + calc_loss_gra  * calc_diss_guts            &
            + calc_loss_gra2 * calc_diss_guts            &
            + calc_diss2                     * DetZ2Calc & 
            - calcification                              & ! --> Calcification
!     + LocRiverDIC                             &
                                                 ) * dt_b + sms(k,idic)
    else
        sms(k,idic)      = (                       &
            - Cphot                         * PhyC    &
            + phyRespRate                   * PhyC    &
            - Cphot_Dia                     * DiaC    &
            + phyRespRate_Dia               * DiaC    &
            + rho_C1 * arrFunc              * EOC     &
            + HetRespFlux                             & 
!#ifdef REcoM_calcification                     
            + calc_diss                     * DetCalc &
            + calc_loss_gra * calc_diss_guts          &
            - calcification                           &
!     + LocRiverDIC                             &
!#endif
                                                 ) * dt_b + sms(k,idic)
    endif
!  if((Latd(1)<-45.0) .and. ((state(k,idic)+sms(k,idic))>2500)) then
!     !co2flux(1)=0.0  
!      print*,'ERROR: strange dic !'
!      print*,'state(k,idic): ', state(k,idic)
!      print*,'sms Cphot: ', -Cphot*PhyC
!      print*,'sms resp: ', phyRespRate*PhyC
!      print*,'sms Cphot dia: ', -Cphot_Dia*DiaC
!      print*,'sms resp dia: ', phyRespRate_Dia * DiaC
!      print*,'sms eoc: ', rho_C1* arrFunc *EOC
!      print*,'sms het resp: ', HetRespFlux
!      print*, 'sms co2: ',  dflux(1) * recipdzF(k) * max( 2-k, 0 )
!      print*, 'sms calcdiss: ', calc_diss * DetCalc
!      print*, 'sms calc_loss: ', calc_loss_gra * calc_diss_guts
!      print*, 'sms calcification: ', -calcification
!      stop
!    endif

!____________________________________________________________
!< Alkalinity (Assumes that N:P follows a constant Redfield ratio

!< N_assimC: 1.0625 = 1/16 + 1
    if (REcoM_Second_Zoo) then
        sms(k,ialk)      = (                       &
            + 1.0625 * N_assim             * PhyC    &
            + 1.0625 * N_assim_Dia         * DiaC    &
            - 1.0625 * rho_N * arrFunc     * DON     &
!#ifdef REcoM_calcification     
            + 2.d0 * calc_diss             * DetCalc &
            + 2.d0 * calc_loss_gra * calc_diss_guts  &
            + 2.d0 * calc_loss_gra2 * calc_diss_guts &
            + 2.d0 * calc_diss2            * DetZ2Calc &
            - 2.d0 * calcification                   &
!#endif                       
                                             ) * dt_b + sms(k,ialk)
    else
        sms(k,ialk)      = (                       &
            + 1.0625 * N_assim             * PhyC    &
            + 1.0625 * N_assim_Dia         * DiaC    &
            - 1.0625 * rho_N * arrFunc     * DON     &
!#ifdef REcoM_calcification
            + 2.d0 * calc_diss             * DetCalc &
            + 2.d0 * calc_loss_gra * calc_diss_guts  &
            - 2.d0 * calcification                   &
!      + LocRiverAlk                             &
!#endif 
                                             ) * dt_b + sms(k,ialk) 
    endif
!____________________________________________________________
!< Small phytoplankton
 
!< lossN:  Phytoplankton loss of organic N compounds [day^-1]
    if (REcoM_Second_Zoo) then
        sms(k,iphyn)      = (                        &
            + N_assim                      * PhyC    & ! -- N assimilation
            - lossN * limitFacN            * PhyN    & ! -- DON excretion
            - aggregationRate              * phyN    & ! -- Aggregation loss
            - grazingFlux_phy                        & ! -- Grazing loss
            - grazingFlux_phy2                       &     
                                                   ) * dt_b + sms(k,iphyn)
    else
        sms(k,iphyn)      = (                        &
            + N_assim                      * PhyC    &
            - lossN * limitFacN            * PhyN    &
            - aggregationRate              * phyN    &
            - grazingFlux_phy                        &
                                                   ) * dt_b + sms(k,iphyn)
   endif
!____________________________________________________________
!< Small phytoplankton C

!< lossC: Phytoplankton loss of carbon [day^-1]
!< When N : C ratio becomes too high, excretion of DOC is downregulated
!< by the limiter function limitFacN
!< aggregationRate transfers C to the detritus pool

    if (REcoM_Second_Zoo) then
        sms(k,iphyc)      = (                        &
            + Cphot                        * PhyC    & ! -- Photosynthesis ---->/
            - lossC * limitFacN            * PhyC    & ! -- Excretion of DOC   / Net photosynthesis
            - phyRespRate                  * PhyC    & ! -- Respiration ----->/
            - aggregationRate              * PhyC    & ! -- Aggregation loss
            - grazingFlux_phy * recipQuota           & ! -- Grazing loss
            - grazingFlux_phy2 * recipQuota          &
                                                   ) * dt_b + sms(k,iphyc)
    else
        sms(k,iphyc)      = (                        &
            + Cphot                        * PhyC    &
            - lossC * limitFacN            * PhyC    &
            - phyRespRate                  * PhyC    &
            - aggregationRate              * PhyC    &
            - grazingFlux_phy * recipQuota           &
                                                   ) * dt_b + sms(k,iphyc)
   endif
!____________________________________________________________
! Phytoplankton ChlA

!< Chl2N: Conversion factor from mmolN to mgChla 
!< Chl2N = PhyChl/PhyN

    if (REcoM_Second_Zoo) then
        sms(k,ipchl)       = (                       &
     	    + chlSynth                     * phyC    & ! -- Chl-a synthesis
     	    - KOchl                        * PhyChl  & ! -- Degradation loss
            - aggregationRate              * PhyChl  & ! -- Aggregation loss
     	    - grazingFlux_phy * Chl2N                & ! -- Grazing loss
            - grazingFlux_phy2 * Chl2N               & 
                                                   ) * dt_b + sms(k,ipchl)
    else
        sms(k,ipchl)       = (                       &
     	    + chlSynth                     * phyC    &
     	    - KOchl                        * PhyChl  &
     	    - aggregationRate              * PhyChl  &
     	    - grazingFlux_phy * Chl2N                &
                                                   ) * dt_b + sms(k,ipchl)
    endif
!-------------------------------------------------------------------------------
! Detritus N
!   if (REcoM_Second_Zoo) then
!   if (Grazing_detritus) then
!    sms(k,idetn)       = (                       &
!        + grazingFlux_phy                        &
!        + grazingFlux_dia                        &
!        - grazingFlux * grazEff                  &
!        + grazingFlux_phy2                       &
!        + grazingFlux_dia2                       &
!        + grazingFlux_het2                       &
!        - grazingFlux2 * grazEff2                &
!        + aggregationRate              * PhyN    &
!        + aggregationRate              * DiaN    &
!        + hetLossFlux                            &
!        + Zoo2LossFlux                           &
!        - reminN * arrFunc             * DetN    &
!                                               ) * dt_b + sms(k,idetn)
!   else

!    sms(k,idetn)       = (                       &
!     	+ grazingFlux_phy                        &
!     	+ grazingFlux_dia                        &
!     	- grazingFlux * grazEff                  &
!     	+ aggregationRate              * PhyN    & 
!     	+ aggregationRate              * DiaN    &
!     	+ hetLossFlux                            &
!     	- reminN * arrFunc             * DetN    &
!                                               ) * dt_b + sms(k,idetn)
!   endif   
!-------------------------------------------------------------------------------
! Detritus C
!   if (REcoM_Second_Zoo) then
!    sms(k,idetc)       = (                             &
!        + grazingFlux_phy * recipQuota                 &
!        - grazingFlux_phy * recipQuota * grazEff       &
!        + grazingFlux_Dia * recipQuota_Dia             &
!        - grazingFlux_Dia * recipQuota_Dia * grazEff   &
!        + grazingFlux_phy2 * recipQuota                &
!        - grazingFlux_phy2 * recipQuota * grazEff2     &
!        + grazingFlux_Dia2 * recipQuota_Dia            &
!        - grazingFlux_Dia2 * recipQuota_Dia * grazEff2 &
!        + grazingFlux_het2 * recipQZoo                 &
!        - grazingFlux_het2 * recipQZoo * grazEff2      &
!        + aggregationRate              * phyC          &
!        + aggregationRate              * DiaC          &
!        + hetLossFlux * recipQZoo                      &
!         + Zoo2LossFlux * recipQZoo2                   &
!        - reminC * arrFunc             * DetC          &
!                                             )   * dt_b + sms(k,idetc)
!   else

!    sms(k,idetc)       = (                           &
!     	+ grazingFlux_phy * recipQuota               &
!     	- grazingFlux_phy * recipQuota * grazEff     &
!     	+ grazingFlux_Dia * recipQuota_Dia           &
!     	- grazingFlux_Dia * recipQuota_Dia * grazEff &
!     	+ aggregationRate              * phyC        &
!     	+ aggregationRate              * DiaC        &
!     	+ hetLossFlux * recipQZoo                    &
!     	- reminC * arrFunc             * DetC        &
!                                             )   * dt_b + sms(k,idetc)
!   endif   
!-------------------------------------------------------------------------------
! Detritus N
   if (Grazing_detritus) then
    sms(k,idetn)       = (                       &
	+ grazingFlux_phy                        &
        - grazingFlux_phy * grazEff              &
        + grazingFlux_dia                        &
        - grazingFlux_dia * grazEff              &
        - grazingFlux_Det * grazEff              & !!!Sloppy feeding is  thought because of grazing flux multiplied with grazeff 
        - grazingFlux_Det2 * grazEff             &
        + aggregationRate              * PhyN    &
        + aggregationRate              * DiaN    &
        + hetLossFlux                            &
        - reminN * arrFunc             * DetN    &
                                               ) * dt_b + sms(k,idetn)
   else
    sms(k,idetn)       = (                       &
        + grazingFlux_phy                        &
        + grazingFlux_dia                        &
        - grazingFlux * grazEff                  &
        + aggregationRate              * PhyN    &
        + aggregationRate              * DiaN    &
        + hetLossFlux                            &
        - reminN * arrFunc             * DetN    &
                                               ) * dt_b + sms(k,idetn)
   end if   
!-------------------------------------------------------------------------------
! Detritus C
   if (Grazing_detritus) then
    sms(k,idetc)       = (                            &
        + grazingFlux_phy * recipQuota                 &
        - grazingFlux_phy * recipQuota * grazEff       &
        + grazingFlux_Dia * recipQuota_Dia             &
        - grazingFlux_Dia * recipQuota_Dia * grazEff   &
        - grazingFlux_Det * recipDet * grazEff         &
        - grazingFlux_Det2 * recipDet2 * grazEff       &     
        + aggregationRate              * phyC          &
        + aggregationRate              * DiaC          &
        + hetLossFlux * recipQZoo                      &
        - reminC * arrFunc             * DetC          &
                                             )   * dt_b + sms(k,idetc)
   else
    sms(k,idetc)       = (                             &
        + grazingFlux_phy * recipQuota                 &
        - grazingFlux_phy * recipQuota * grazEff       &
        + grazingFlux_Dia * recipQuota_Dia             &
        - grazingFlux_Dia * recipQuota_Dia * grazEff   &
        + aggregationRate              * phyC          &
        + aggregationRate              * DiaC          &
        + hetLossFlux * recipQZoo                      &
        - reminC * arrFunc             * DetC          &
                                             )   * dt_b + sms(k,idetc)
   end if
!____________________________________________________________
!< Heterotrophic N
    if (REcoM_Second_Zoo) then
        sms(k,ihetn)       = (                       &
    	    + grazingFlux * grazEff                  & ! -- Grazig on phytoplankton
            - grazingFlux_het2                       &
     	    - hetLossFlux                            & ! -- Mortality
     	    - lossN_z                      * HetN    & ! -- Excretion of DON
                                                   ) * dt_b + sms(k,ihetn)
    else
        sms(k,ihetn)       = (                       &
     	    + grazingFlux * grazEff                  &
     	    - hetLossFlux                            &
     	    - lossN_z                      * HetN    &
                                                   ) * dt_b + sms(k,ihetn)
    endif   
!____________________________________________________________
!< Heterotrophic C

    if (REcoM_Second_Zoo) then
        if (Grazing_detritus) then
            sms(k,ihetc)      = (                            &
     	        + grazingFlux_phy * recipQuota * grazEff     & ! -- Grazing on small phytoplankton
     	        + grazingFlux_Dia * recipQuota_Dia * grazEff & ! -- Grazing on diatom
                + grazingFlux_Det * recipDet * grazEff       & ! -- Grazing on detritus
                + grazingFlux_DetZ2 * recipDet2 * grazEff    &
     	        - hetLossFlux * recipQZoo                    & ! -- Mortality loss
                - grazingFlux_het2 * recipQZoo               &
     	        - lossC_z                      * HetC        & ! -- Excretion loss
      	        - hetRespFlux                                & ! -- REspiration loss
                                                        ) * dt_b + sms(k,ihetc)
        else
            sms(k,ihetc)      = (                            &
     	        + grazingFlux_phy * recipQuota * grazEff     &
     	        + grazingFlux_Dia * recipQuota_Dia * grazEff &
     	        - hetLossFlux * recipQZoo                    &
                - grazingFlux_het2 * recipQZoo               &
     	        - lossC_z                      * HetC        &
     	        - hetRespFlux                                & 
                                                        ) * dt_b + sms(k,ihetc)
        endif
    else
        sms(k,ihetc)      = (                              &
            + grazingFlux_phy * recipQuota * grazEff     &
            + grazingFlux_Dia * recipQuota_Dia * grazEff &
            - hetLossFlux * recipQZoo                    &
            - lossC_z                      * HetC        &
            - hetRespFlux                                &
                                                        ) * dt_b + sms(k,ihetc)
    endif
!____________________________________________________________
   if (REcoM_Second_Zoo) then
 ! Second Zooplankton N                                                                                              
     sms(k,izoo2n)       = (                        &
         + grazingFlux2 * grazEff2                  &
         - Zoo2LossFlux                             &
         - lossN_z2                      * Zoo2N    &
         - Zoo2fecalloss_n                            & 
                                               ) * dt_b + sms(k,izoo2n)
 !-------------------------------------------------------------------------------                       
  ! Second Zooplankton C                                                                                 
                  
     if (Grazing_detritus) then
             
     sms(k,izoo2c)      = (                             &
         + grazingFlux_phy2 * recipQuota * grazEff2     &
         + grazingFlux_Dia2 * recipQuota_Dia * grazEff2 &
         + grazingFlux_het2 * recipQZoo * grazEff2      &
         + grazingFlux_Det2 * recipDet * grazEff2       &
         + grazingFlux_DetZ22 * recipDet2 * grazEff2    &
         - zoo2LossFlux * recipQZoo2                    &
         - lossC_z2                      * Zoo2C        &
         - Zoo2RespFlux                                 &
         - Zoo2fecalloss_c                              &
                                                ) * dt_b + sms(k,izoo2c)  
     else
      sms(k,izoo2c)      = (                             &
         + grazingFlux_phy2 * recipQuota * grazEff2     &
         + grazingFlux_Dia2 * recipQuota_Dia * grazEff2 &
         + grazingFlux_het2 * recipQZoo * grazEff2      &
         - zoo2LossFlux * recipQZoo2                    &
         - lossC_z2                      * Zoo2C        &
         - Zoo2RespFlux                                 &
         - Zoo2fecalloss_c                              &
                                                ) * dt_b + sms(k,izoo2c)
     end if   
!---------------------------------------------------------------------------------
  ! Second Zooplankton Detritus N
     if (Grazing_detritus) then
      sms(k,idetz2n)       = (                       &
         + grazingFlux_phy2                       &
         - grazingFlux_phy2 * grazEff2            &         
         + grazingFlux_dia2                       &
         - grazingFlux_dia2 * grazEff2            & 
         + grazingFlux_het2                       &
         - grazingFlux_het2 * grazEff2            &
         - grazingFlux_DetZ2 * grazEff2           &
         - grazingFlux_DetZ22 * grazEff2          &   
         + Zoo2LossFlux                           &
         + Zoo2fecalloss_n                          &
         - reminN * arrFunc             * DetZ2N  &
                                               ) * dt_b + sms(k,idetz2n)
     else
      sms(k,idetz2n)       = (                       &
         + grazingFlux_phy2                       &
         + grazingFlux_dia2                       &
         + grazingFlux_het2                       &
         - grazingFlux2 * grazEff2                &
         + Zoo2LossFlux                           &
         + Zoo2fecalloss_n                          &
         - reminN * arrFunc             * DetZ2N  &
                                               ) * dt_b + sms(k,idetz2n)
     end if
 !---------------------------------------------------------------------------------                                                                                            
                                                           
  ! Second Zooplankton Detritus C
     if (Grazing_detritus) then
      sms(k,idetz2c)       = (                             &
        + grazingFlux_phy2 * recipQuota                &
        - grazingFlux_phy2 * recipQuota * grazEff2     &
        + grazingFlux_Dia2 * recipQuota_Dia            &
        - grazingFlux_Dia2 * recipQuota_Dia * grazEff2 &
        + grazingFlux_het2 * recipQZoo                 &
        - grazingFlux_het2 * recipQZoo * grazEff2      &
        - grazingFlux_DetZ2 * recipDet * grazEff2      &
        - grazingFlux_DetZ22 * recipDet2 * grazEff2    &
        + Zoo2LossFlux * recipQZoo2                    &
        + Zoo2fecalloss_c                   & 
        - reminC * arrFunc             * DetZ2C          &
                                             )   * dt_b + sms(k,idetz2c)
     else
      sms(k,idetz2c)       = (                             &
        + grazingFlux_phy2 * recipQuota                &
        - grazingFlux_phy2 * recipQuota * grazEff2     &
        + grazingFlux_Dia2 * recipQuota_Dia            &
        - grazingFlux_Dia2 * recipQuota_Dia * grazEff2 &
        + grazingFlux_het2 * recipQZoo                 &
        - grazingFlux_het2 * recipQZoo * grazEff2      &
        + Zoo2LossFlux * recipQZoo2                    &
        + Zoo2fecalloss_c                    &
        - reminC * arrFunc             * DetZ2C          &
                                             )   * dt_b + sms(k,idetz2c)
     end if

 !-------------------------------------------------------------------------------------
  !Second Zooplankton  Detritus Si                                                                                                                                                                                                      
     sms(k,idetz2si)     = (                          &
         + grazingFlux_dia2 * qSiN                    &  ! --> qSin convert N to Si
         - reminSiT                        * DetZ2Si  &
                                             ) * dt_b + sms(k,idetz2si)
 !-------------------------------------------------------------------------------                                                                                              
 ! Detritus calcite   
     sms(k,idetz2calc)   = (               &
       + calc_loss_gra2                  &
       - calc_loss_gra2 * calc_diss_guts &
       - calc_diss2     * DetZ2Calc       &
                                           ) * dt_b + sms(k,idetz2calc)
   endif                                                                               
!-------------------------------------------------------------------------------
! DON (Extracellular organic N)
   if (REcoM_Second_Zoo) then
    sms(k,idon)      = (                        &
      + lossN * limitFacN              * phyN   &
      + lossN_d * limitFacN_Dia        * DiaN   &
      + reminN * arrFunc               * DetN   &
      + reminN * arrFunc               * DetZ2N &
      + lossN_z                        * HetN   &
      + lossN_z2                       * Zoo2N  &
      - rho_N * arrFunc                * DON    &
!      + LocRiverDON                             &
                                             ) * dt_b + sms(k,idon)
   else

    sms(k,idon)      = (                        &
      + lossN * limitFacN              * phyN   &
      + lossN_d * limitFacN_Dia        * DiaN   &
      + reminN * arrFunc               * DetN   &
      + lossN_z                        * HetN   &
      - rho_N * arrFunc                * DON    &
!      + LocRiverDON                             &
                                              ) * dt_b + sms(k,idon)
   endif
!-------------------------------------------------------------------------------
! EOC
   if (REcoM_Second_Zoo) then
    sms(k,idoc)       = (                       &
      + lossC * limitFacN              * phyC   &
      + lossC_d * limitFacN_dia        * DiaC   &
      + reminC * arrFunc               * DetC   &
      + reminC * arrFunc               * DetZ2C &
      + lossC_z                        * HetC   &
      + lossC_z2                       * Zoo2C  &
      - rho_c1 * arrFunc               * EOC    &
!      + LocRiverDOC                             &
                                              ) * dt_b + sms(k,idoc)	
   else 
    sms(k,idoc)       = (                       &
      + lossC * limitFacN              * phyC   &
      + lossC_d * limitFacN_dia        * DiaC   &
      + reminC * arrFunc               * DetC   &
      + lossC_z                        * HetC   &
      - rho_c1 * arrFunc               * EOC    &
!      + LocRiverDOC                             &
                                              ) * dt_b + sms(k,idoc)
   endif 		
!____________________________________________________________
!< Diatom N
 
!< lossN:  Diatom loss of organic N compounds [day^-1]
!< When N : C ratio becomes too high, excretion of DON is downregulated
!< by the limiter function limitFacN_dia
!< aggregationRate transfers N to the detritus pool
    if (REcoM_Second_Zoo) then
        sms(k,idian)      = (                      &
            + N_assim_dia                    * DiaC  & ! -- N assimilation
            - lossN_d * limitFacN_dia        * DiaN  & ! -- DON excretion
            - aggregationRate                * DiaN  & ! -- Aggregation loss
            - grazingFlux_Dia                        & ! -- Grazing loss
            - grazingFlux_Dia2                       &
                                                 ) * dt_b + sms(k,idian)
    else
        sms(k,idian)      = (                      &
            + N_assim_dia                    * DiaC  &
            - lossN_d * limitFacN_dia        * DiaN  &
            - aggregationRate                * DiaN  &
            - grazingFlux_Dia                        &
                                                 ) * dt_b + sms(k,idian)
    endif 		
!____________________________________________________________
!< Diatom C

!< lossC_d: Diatom loss of carbon [day^-1]
!< When N : C ratio becomes too high, excretion of DOC is downregulated
!< by the limiter function limitFacN_dia
!< aggregationRate transfers C to the detritus pool
    if (REcoM_Second_Zoo) then
        sms(k,idiac)      = (                      &
            + Cphot_dia                      * DiaC  & ! -- Photosynthesis ---->/
            - lossC_d * limitFacN_dia        * DiaC  & ! -- Excretion of DOC --/ Net Photosynthesis
            - phyRespRate_dia                * DiaC  & ! -- Respiration ----->/
            - aggregationRate                * DiaC  &
            - grazingFlux_dia * recipQuota_dia       &
            - grazingFlux_dia2 * recipQuota_dia      &
     	                                           ) * dt_b + sms(k,idiac)
    else
        sms(k,idiac)      = (                      &
            + Cphot_dia                      * DiaC  &
            - lossC_d * limitFacN_dia        * DiaC  &
            - phyRespRate_dia                * DiaC  &
            - aggregationRate                * DiaC  &
            - grazingFlux_dia * recipQuota_dia       &
     	                                           ) * dt_b + sms(k,idiac)
    endif
!____________________________________________________________
!< Diatom Chl

    if (REcoM_Second_Zoo) then
        sms(k,idchl)      = (                         &
            + chlSynth_dia                   * DiaC   & ! -- Chl a synthesis
            - KOchl_dia                      * DiaChl & ! -- Degradation loss
            - aggregationRate                * DiaChl & ! -- Aggregation loss
            - grazingFlux_dia * Chl2N_dia             & ! -- Grazing loss
            - grazingFlux_dia2 * Chl2N_dia            &                 
                                                   ) * dt_b + sms(k,idchl)
    else
        sms(k,idchl)      = (                         &
            + chlSynth_dia                   * DiaC   &
            - KOchl_dia                      * DiaChl &
            - aggregationRate                * DiaChl &
            - grazingFlux_dia * Chl2N_dia             &
                                                   ) * dt_b + sms(k,idchl)
    endif 
!____________________________________________________________
!< Diatom Si
!< lossN_d: Diatom loss of organic nitrogen compunds [day^-1]
!< When N : C ratio becomes too high, excretion is downregulated
!< by the limiter function limitFacN_dia
!< aggregationRate transfers Si to the detritus pool

    if (REcoM_Second_Zoo) then
        sms(k,idiasi)        = (                      &
            + Si_assim                        * DiaC  & ! -- Diatom silicon assimilation
            - lossN_d * limitFacN_dia         * DiaSi & ! -- Excretion to detritus
            - aggregationRate                 * DiaSi & ! -- Aggregation loss
            - grazingFlux_dia * qSiN                  & ! -- Grazing loss
            - grazingFlux_dia2 * qSiN                 &
                                                   ) * dt_b + sms(k,idiasi)
    else
        sms(k,idiasi)        = (                    &
            + Si_assim                        * DiaC  &
            - lossN_d * limitFacN_dia         * DiaSi &
            - aggregationRate                 * DiaSi &
            - grazingFlux_dia * qSiN                  &
                                                   ) * dt_b + sms(k,idiasi)
    endif 
!-------------------------------------------------------------------------------
! Detritus Si
   if (REcoM_Second_Zoo) then
    sms(k,idetsi)     = (                       &
      + aggregationRate                 * DiaSi &
      + lossN_d * limitFacN_dia         * DiaSi &
      + grazingFlux_dia * qSiN                  &
!      + grazingFlux_dia2 * qSiN                 & should not be here
      - reminSiT                        * DetSi &
                                             ) * dt_b + sms(k,idetsi)
   else
    sms(k,idetsi)     = (                       &
      + aggregationRate                 * DiaSi &
      + lossN_d * limitFacN_dia         * DiaSi &
      + grazingFlux_dia * qSiN                  &
      - reminSiT                        * DetSi &
                                             ) * dt_b + sms(k,idetsi)
   endif
!____________________________________________________________
!< DSi, Silicate 

!    if (mype==0) print*,' sms(k,idin) ', sms(k,idin)

!< DiaC:        Intracellular carbon concentration in diatoms [mmolC m-3]
!< DetSi:       Detritus silicon concentration [mmolSi m-3]
!< Si_assim:    Si assimilation rate for diatoms [mmolSi mmolC-1 day-1]
!< reminSiT:    Remineralization rate of silicon, temperature dependency [day-1]
!< dt_b:        REcoM time step [day]

!! Schourup 2013 Eq. A3
    if (REcoM_Second_Zoo) then
        sms(k,isi)        = (                           &
            - Si_assim                        * DiaC    &  ! --> Si assimilation of diatoms
            + reminSiT                        * DetSi   &  ! --> Remineralization of detritus, temperature dependent
            + reminSiT                        * DetZ2Si &
!             + LocRiverDSi                             &  ! --> added in FESOM2 
                                             ) * dt_b + sms(k,isi)
    else
        sms(k,isi)        = (                           &
            - Si_assim                        * DiaC    &
            + reminSiT                        * DetSi   &
!        + LocRiverDSi                                  &

                                             ) * dt_b + sms(k,isi)
   endif

!____________________________________________________________
!< Fe

!< Fe2C: Intracellular Fe : C ratio [μmol Fe mmol C^-1]
!< Fe2N: Intracellular Fe : N ratio [μmol Fe mmol N^-1] Fe2N = Fe2C * 6.625
!< PhyC: Intracellular carbon concentration in nanophytoplankton [mmolCm^-3]
!< Cphot: C-specific actual rate of photosynthesis for nanopyhtoplankton [day^-1]
!< DiaC: Intracellular carbon concentration in diatoms [mmol C m^-3 ]
!< Cphot_dia: C-specific actual rate of photosynthesis for diatom [day^-1]
!< phyRespRate: Nanopyhtoplankton respiration rate [day^-1]
!< phyRespRate_dia: Diatom respiration rate [day^-1]
!< lossC: Nanopyhtoplankton excretion of organic C [day^-1]
!< limitFacN: limiting factor
!< lossC_d: Diatom excretion of organic C [day^-1]
!< limitFacN_dia: limiting factor
!< detC: Detritus carbon concentration [mmol C m^-3]
!< reminC: Temperature dependent remineralisation rate of detritus [day^-1]
!< arrFunc: Arrhenius function
!< hetC: Zooplankton carbon concentration [mmol C m^-3 ]
!< lossC_z: Zooplankton excretion of organic C [day^-1 ]
!< hetRespFlux: Zooplankton respiration rate [day^-1]
!< kScavFe: Scavenging rate of iron [m3 mmol C^-1 day^-1]

    if (REcoM_Second_Zoo) then
        if (use_Fe2N) then
            sms(k,ife) = ( Fe2N * (                  &
                - N_assim                 * PhyC     & ! --> N assimilation Nanophytoplankton, [mmol N/(mmol C * day)] C specific N utilization rate  
                - N_assim_dia             * DiaC     & ! --> N assimilation Diatom
                + lossN*limitFacN         * PhyN     & ! --> Excretion from small pythoplankton
                + lossN_d*limitFacN_dia   * DiaN     & ! --> Excretion from diatom
                + reminN * arrFunc        * DetN     & ! --> Remineralization of detritus
                + reminN * arrFunc        * DetZ2N   &
                + lossN_z                 * HetN     & ! --> Excretion from zooplanton
                + lossN_z2                * Zoo2N    &         
                                               )     &
                - kScavFe                 * DetC   * FreeFe   & 
                - kScavFe                 * DetZ2C * FreeFe   &
                                               ) * dt_b + sms(k,ife)
        else
            sms(k,ife)      = ( Fe2C *(          &
                -  Cphot                  * PhyC     & ! Small pyhtoplankton photosynthesis ---/
                -  Cphot_dia              * DiaC     & ! Diatom photosynthesis                / -> net growth
                +  phyRespRate            * PhyC     & ! Small pyhtoplankton respiration --- /
                +  phyRespRate_dia        * DiaC     & ! Diatom respiration
                +  lossC*limitFacN        * phyC     & ! Exrcetion from small pythoplankton
                +  lossC_d*limitFacN_dia  * diaC     & ! Excretion from diatom
                +  reminC * arrFunc       * detC     & ! Remineralization of detritus
                +  reminC * arrFunc       * DetZ2C   &
                +  lossC_z                * hetC     & ! Excretion from zooplanton
                +  hetRespFlux                       & ! Zooplankton respiration
                +  lossC_z2               * Zoo2C    &
                +  zoo2RespFlux                      &
                                               )     &
                -  kScavFe                * DetC   * FreeFe   & ! Scavenging of free iron (correlated with detC)
                -  kScavFe                * DetZ2C * FreeFe   &   
                                                ) * dt_b + sms(k,ife)
        end if
   else

        if (use_Fe2N) then
            sms(k,ife) = ( Fe2N * (                   &
                - N_assim                 * PhyC      &
                - N_assim_dia             * DiaC      &
                + lossN*limitFacN         * PhyN      &
                + lossN_d*limitFacN_dia   * DiaN      &
                + reminN * arrFunc        * DetN      &
                + lossN_z                 * HetN )    &
                - kScavFe                 * DetC * FreeFe &
                                              ) * dt_b + sms(k,ife)
        else
            sms(k,ife)      = ( Fe2C *(          &
                -  Cphot                  * PhyC     &
                -  Cphot_dia              * DiaC     &
                +  phyRespRate            * PhyC     &
                +  phyRespRate_dia        * DiaC     &
                +  lossC*limitFacN        * phyC     &
                +  lossC_d*limitFacN_dia  * diaC     &
                +  reminC * arrFunc       * detC     &
                +  lossC_z                * hetC     &
                +  hetRespFlux  )                    &
                -  kScavFe                * DetC * FreeFe   & 
                                              ) * dt_b + sms(k,ife)
    end if
   endif
!____________________________________________________________
!< Small phytoplankton calcite

!#ifdef REcoM_calcification
    if (REcoM_Second_Zoo) then
        sms(k,iphycal)    = (             &
            + calcification               & ! -- Calcification
            - lossC * limitFacN * phyCalc & ! -- Excretion loss
            - phyRespRate       * phyCalc & ! -- Respiration
            - calc_loss_agg               & ! -- Aggregation loss
            - calc_loss_gra               & ! -- Grazing loss
            - calc_loss_gra2              &
                                                  ) * dt_b + sms(k,iphycal)
    else
        sms(k,iphycal)    = (             &
            + calcification               &
            - lossC * limitFacN * phyCalc &
            - phyRespRate       * phyCalc &
            - calc_loss_agg               &
            - calc_loss_gra               &
                                                  ) * dt_b + sms(k,iphycal)
    endif

!-------------------------------------------------------------------------------
! Detritus calcite
    sms(k,idetcal)   = (               &
      + lossC * limitFacN * phyCalc    &
      + phyRespRate       * phyCalc    &
      + calc_loss_agg                  &
      + calc_loss_gra                  &
      - calc_loss_gra * calc_diss_guts &
      - calc_diss     * DetCalc        &
                                           ) * dt_b + sms(k,idetcal)
!#endif
!-------------------------------------------------------------------------------
! Oxygen
   if (REcoM_Second_Zoo) then
    sms(k,ioxy)   = (               &
      + Cphot              * phyC  &
      - phyRespRate         * phyC  &
      + Cphot_dia          * diaC  &
      - phyRespRate_dia     * diaC  &
      - rho_C1  * arrFunc   * EOC   &
      - hetRespFlux                 &
      - Zoo2RespFlux                 &
                                        )*redO2C * dt_b + sms(k,ioxy)  
   else
    sms(k,ioxy)   = (               &
      + Cphot              * phyC  &
      - phyRespRate         * phyC  &
      + Cphot_dia          * diaC  &
      - phyRespRate_dia     * diaC  &
      - rho_C1  * arrFunc   * EOC   &
      - hetRespFlux                 &
                                      )*redO2C * dt_b + sms(k,ioxy)
   endif

!
!   
    if (ciso) then
!-------------------------------------------------------------------------------
! DIC_13
        sms(k,idic_13) =        (                           &
            - Cphot                         * PhyC_13       &
            + phyRespRate                   * PhyC_13       &
            - Cphot_Dia                     * DiaC_13       &
            + phyRespRate_Dia               * DiaC_13       &
            + rho_C1 * arrFunc              * EOC_13        &
            + HetRespFlux_13                                &
            + calc_diss_13                  * DetCalc_13    &
            + calc_loss_gra_13 * calc_diss_guts             &
            - calcification_13                              &
                                ) * dt_b + sms(k,idic_13)
!-------------------------------------------------------------------------------
! DIC_14
        sms(k,idic_14) =        (                           &
            - Cphot                         * PhyC_14       &
            + phyRespRate                   * PhyC_14       &
            - Cphot_Dia                     * DiaC_14       &
            + phyRespRate_Dia               * DiaC_14       &
            + rho_C1 * arrFunc              * EOC_14        &
            + HetRespFlux_14                                &
            + calc_diss_14                  * DetCalc_14    &
            + calc_loss_gra_14 * calc_diss_guts             &
            - calcification_14                              &
                                ) * dt_b + sms(k,idic_14)
!-------------------------------------------------------------------------------
! Phytoplankton C_13
            sms(k,iphyc_13)      =  (                           &
                + Cphot                        * PhyC_13        &
                - lossC * limitFacN            * PhyC_13        &
                - phyRespRate                  * PhyC_13        &
                - aggregationRate              * PhyC_13        &
                - grazingFlux_phy * recipQuota_13               &
                                    ) * dt_b + sms(k,iphyc_13)
!-------------------------------------------------------------------------------
! Phytoplankton C_14
            sms(k,iphyc_14)      =  (                           &
                + Cphot                        * PhyC_14        &
                - lossC * limitFacN            * PhyC_14        &
                - phyRespRate                  * PhyC_14        &
                - aggregationRate              * PhyC_14        &
                - grazingFlux_phy * recipQuota_14               &
                                    ) * dt_b + sms(k,iphyc_14)
!-------------------------------------------------------------------------------
! Detritus C_13
            sms(k,idetc_13)       = (                           &
                + grazingFlux_phy * recipQuota_13               &
                - grazingFlux_phy * recipQuota_13 * grazEff     &
                + grazingFlux_Dia * recipQuota_dia_13           &
                - grazingFlux_Dia * recipQuota_dia_13 * grazEff &
                + aggregationRate              * phyC_13        &
                + aggregationRate              * DiaC_13        &
                + hetLossFlux * recipQZoo_13                    &
                - reminC * arrFunc             * DetC_13        &
                                                   )   * dt_b + sms(k,idetc_13)
!-------------------------------------------------------------------------------
! Detritus C_14
            sms(k,idetc_14)       = (                           &
                + grazingFlux_phy * recipQuota_14               &
                - grazingFlux_phy * recipQuota_14 * grazEff     &
                + grazingFlux_Dia * recipQuota_dia_14           &
                - grazingFlux_Dia * recipQuota_dia_14 * grazEff &
                + aggregationRate              * phyC_14        &
                + aggregationRate              * DiaC_14        &
                + hetLossFlux * recipQZoo_14                    &
                - reminC * arrFunc             * DetC_14        &
                                                   )   * dt_b + sms(k,idetc_14)
!-------------------------------------------------------------------------------
! Heterotrophic C_13
        sms(k,ihetc_13)      = (                            &
            + grazingFlux_phy * recipQuota_13 * grazEff     &
            + grazingFlux_Dia * recipQuota_dia_13 * grazEff &
            - hetLossFlux * recipQZoo_13                    &
            - lossC_z                      * HetC_13        &
            - hetRespFlux_13                                &
                                                 ) * dt_b + sms(k,ihetc_13)
!-------------------------------------------------------------------------------
! Heterotrophic C_14
        sms(k,ihetc_14)      = (                            &
            + grazingFlux_phy * recipQuota_14 * grazEff     &
            + grazingFlux_Dia * recipQuota_dia_14 * grazEff &
            - hetLossFlux * recipQZoo_14                    &
            - lossC_z                      * HetC_14        &
            - hetRespFlux_14                                &
                                                 ) * dt_b + sms(k,ihetc_14)
!-------------------------------------------------------------------------------
! EOC_13
            sms(k,idoc_13)       = (                            &
                + lossC * limitFacN              * phyC_13      &
                + lossC_d * limitFacN_dia        * DiaC_13      &
                + reminC * arrFunc               * DetC_13      &
                + lossC_z                        * HetC_13      &
                - rho_c1 * arrFunc               * EOC_13       &
                + LocRiverDOC * alpha_iorg_13                   &
                                                    ) * dt_b + sms(k,idoc_13)
!-------------------------------------------------------------------------------
! EOC_14
            sms(k,idoc_14)       = (                            &
                + lossC * limitFacN              * phyC_14      &
                + lossC_d * limitFacN_dia        * DiaC_14      &
                + reminC * arrFunc               * DetC_14      &
                + lossC_z                        * HetC_14      &
                - rho_c1 * arrFunc               * EOC_14       &
                + LocRiverDOC * alpha_iorg_14                   &
                                                    ) * dt_b + sms(k,idoc_14)
!-------------------------------------------------------------------------------
! Diatom C_13
            sms(k,idiac_13)      = (                            &
                + Cphot_dia                      * DiaC_13      &
                - lossC_d * limitFacN_dia        * DiaC_13      &
                - phyRespRate_dia                * DiaC_13      &
                - aggregationRate                * DiaC_13      &
                - grazingFlux_dia * recipQuota_dia_13           &
                                                   ) * dt_b + sms(k,idiac_13)
!-------------------------------------------------------------------------------
! Diatom C_14
            sms(k,idiac_14)      = (                            &
                + Cphot_dia                      * DiaC_14      &
                - lossC_d * limitFacN_dia        * DiaC_14      &
                - phyRespRate_dia                * DiaC_14      &
                - aggregationRate                * DiaC_14      &
                - grazingFlux_dia * recipQuota_dia_14           &
                                                   ) * dt_b + sms(k,idiac_14)
!-------------------------------------------------------------------------------
! Small phytoplankton calcite_13
        sms(k,iphycal_13)    = (                            &
            + calcification_13                              &
            - lossC * limitFacN * phyCalc_13                &
            - phyRespRate       * phyCalc_13                &
            - calc_loss_agg_13                              &
            - calc_loss_gra_13                              &
                                              ) * dt_b + sms(k,iphycal_13)
!-------------------------------------------------------------------------------
! Small phytoplankton calcite_14
        sms(k,iphycal_14)    = (                            &
            + calcification_14                              &
            - lossC * limitFacN * phyCalc_14                &
            - phyRespRate       * phyCalc_14                &
            - calc_loss_agg_14                              &
            - calc_loss_gra_14                              &
                                              ) * dt_b + sms(k,iphycal_14)
!-------------------------------------------------------------------------------
! Detritus calcite_13
        sms(k,idetcal_13)   = (                             &
            + lossC * limitFacN * phyCalc_13                &
            + phyRespRate       * phyCalc_13                &
            + calc_loss_agg_13                              &
            + calc_loss_gra_13                              &
            - calc_loss_gra_13 * calc_diss_guts             &
            - calc_diss_13     * DetCalc_13                 &
                                             ) * dt_b + sms(k,idetcal_13)
!-------------------------------------------------------------------------------
! Detritus calcite_14
        sms(k,idetcal_14)   = (                             &
            + lossC * limitFacN * phyCalc_14                &
            + phyRespRate       * phyCalc_14                &
            + calc_loss_agg_14                              &
            + calc_loss_gra_14                              &
            - calc_loss_gra_14 * calc_diss_guts             &
            - calc_diss_14     * DetCalc_14                 &
                                             ) * dt_b + sms(k,idetcal_14)
!-------------------------------------------------------------------------------

    end if ! ciso

!-------------------------------------------------------------------------------
! Diagnostics: Averaged rates
	
	recipbiostep    = 1.d0/real(biostep)

!*** Net primary production [mmol C /(m3 * day)]
	Diags3Dloc(k,1) = Diags3Dloc(k,1) + (   &
     	+ Cphot                   * PhyC  &
     	- PhyRespRate             * PhyC  &
     	) * recipbiostep

	Diags3Dloc(k,2) = Diags3Dloc(k,2) + (   &
     	+ Cphot_dia               * DiaC  &
     	- PhyRespRate_dia         * DiaC  &
     	) * recipbiostep

!*** Gross primary production [mmol C /(m3 * day)]
	Diags3Dloc(k,3) = Diags3Dloc(k,3) + (   &
     	+ Cphot                   * PhyC  &
     	) * recipbiostep

	Diags3Dloc(k,4) = Diags3Dloc(k,4) + (   &
     	+ Cphot_dia               * DiaC  &
     	) * recipbiostep

!*** Net N-assimilation [mmol N/(m3 * day)]
	Diags3Dloc(k,5) = Diags3Dloc(k,5) + (   &
     	+ N_assim                 * PhyC  &
     	- lossN * limitFacN       * PhyN  &
     	) * recipbiostep

	Diags3Dloc(k,6) = Diags3Dloc(k,6) + (   &
     	+ N_assim_dia             * DiaC  &
     	- lossN * limitFacN_dia   * DiaN  &
     	) * recipbiostep

!*** Changed to chlorophyll degradation (commented out gross N-assimilation below)
        Diags3Dloc(k,7) = Diags3Dloc(k,7) + (   &
        + KOchl  &
        ) * recipbiostep

        Diags3Dloc(k,8) = Diags3Dloc(k,8) + (   &
        + KOchl_dia  &
        ) * recipbiostep
!	Diags3Dloc(k,7) = Diags3Dloc(k,7) + (   &
!     	+ N_assim                 * PhyC  &
!     	) * recipbiostep

!	Diags3Dloc(k,8) = Diags3Dloc(k,8) + (   &
!     	+ N_assim_dia             * DiaC  &
!     	) * recipbiostep

  end do ! Main vertikal loop ends

!-------------------------------------------------------------------------------
! Remineralization from the sediments into the bottom layer
!		kLoc = Nn      
      
!   if (mype==0) then
!      write(*,*) '____________________________________________________________'
!      write(*,*) ' --> Remineralization from the sediments'
!      write(*,*) '     into the bottom layer'
!   endif

!*** DIN ***
!< decayRateBenN: Remineralization rate for benthic N [day^-1]
!< LocBenthos(1): Vertically integrated N concentration in benthos (1 layer) [mmolN/m^2]
    decayBenthos(1) = decayRateBenN * LocBenthos(1)                                       ! flux
    LocBenthos(1)   = LocBenthos(1)   - decaybenthos(1) * dt_b ! / depth of benthos    ! remove from benthos (flux)
!    sms(Nn,idin)    = sms(Nn,idin) + decayBenthos(1) * dt_b  &
!                      / bottom_node_thickness(n) !recipthick(Nn) 			  ! add to water column (flux)

!*** DIC ***
!< decayRateBenC: Remineralization rate for benthic C [day^-1]
!< LocBenthos(2): Vertically integrated C concentration in benthos (1 layer) [mmolC/m^2]		
    decayBenthos(2) = decayRateBenC * LocBenthos(2)
    LocBenthos(2)   = LocBenthos(2)   - decaybenthos(2) * dt_b ! / depth of benthos
!    sms(Nn,idic)    = sms(Nn,idic) + decayBenthos(2) * dt_b  &
!                      / bottom_node_thickness(n) !recipthick(Nn)

!*** Si ***
!< decayRateBenSi: Remineralization rate for benthic Si [day^-1]
!< LocBenthos(3) : Vertically integrated N concentration in benthos (1 layer) [mmolSi/m^2]
		
    decayBenthos(3) = decayRateBenSi * LocBenthos(3)                                      ! [Si/day] * [mmolSi/m2] -> [mmolSi/m2/day]
    LocBenthos(3)   = LocBenthos(3)   - decaybenthos(3) * dt_b ! / depth of benthos    ! [mmolSi/m2]
!    sms(Nn,isi)     = sms(Nn,isi)  + decayBenthos(3) * dt_b  &
!                      / bottom_node_thickness(n) !recipthick(Nn)

!*** Calc: DIC, Alk ***
!#ifdef REcoM_calcification
    decayBenthos(4) = calc_diss * LocBenthos(4)
    LocBenthos(4)      = LocBenthos(4)   - decayBenthos(4) * dt_b ! / depth of benthos
!    sms(Nn,idic)    = sms(Nn,idic) + decayBenthos(4) * dt_b  &
!                      / bottom_node_thickness(n) !* recipthick(Nn) 
!    sms(Nn,ialk)    = sms(Nn,ialk) + decayBenthos(4) * dt_b  &
!                      / bottom_node_thickness(n)* 2.d0  !recipthick(Nn)* 2.d0 
!#endif
    if (ciso) then
!*** DIC_13 ***  We ignore isotopic fractionation during remineralization.
        decayBenthos(5) = alpha_dcal_13   * decayRateBenC   * LocBenthos(5)
        LocBenthos(5)   = LocBenthos(5)   - decayBenthos(5) * dt_b
!        sms(Nn,idic_13) = sms(Nn,idic_13) + decayBenthos(5) * dt_b  &
!                          * recipthick(Nn) !recipdzF(Nn)
!*** DIC_14 ***  We ignore isotopic fractionation during remineralization.
        decayBenthos(6) = alpha_dcal_14   * decayRateBenC   * LocBenthos(6)
        LocBenthos(6)   = LocBenthos(6)   - decayBenthos(6) * dt_b
!        sms(Nn,idic_14) = sms(Nn,idic_14) + decayBenthos(6) * dt_b &
!                          * recipthick(Nn) !recipdzF(Nn)
!*** Calc: DIC_13,14 ***
        decayBenthos(7) = calc_diss_13    * LocBenthos(7)
        LocBenthos(7)   = LocBenthos(7)   - decayBenthos(7) * dt_b ! / depth of benthos
!        sms(Nn,idic_13) = sms(Nn,idic_13) + decayBenthos(7) * dt_b &
!                          * recipthick(Nn) !recipdzF(Nn)
        decayBenthos(8) = calc_diss_14    * LocBenthos(8)
        LocBenthos(8)   = LocBenthos(8)   - decayBenthos(8) * dt_b ! / depth of benthos
!        sms(Nn,idic_14) = sms(Nn,idic_14) + decayBenthos(8) * dt_b &
!                          * recipthick(Nn) ! recipdzF(Nn)
    end if ! ciso
!*** DFe ***
!  if(use_Fe2N) then 
!    Ironflux          = decayRateBenN * LocBenthos(1) * Fe2N_benthos
!  else
!    Ironflux          = decayRateBenC * LocBenthos(2) * Fe2C_benthos
!  end if
!  sms(Nn,ife)       = sms(Nn,ife) + Ironflux * recipthick(Nn) * dt_b

!*** O2 ***
!   sms(Nn,ioxy)    = sms(Nn,ioxy) - decayBenthos(2) * redO2C *dt_b  &
!                      * recipthick(Nn)


  end do ! Main time loop ends

!-------------------------------------------------------------------------------
! Sinking

  dt_sink = dt_b * real(biostep)	
!  dt = dt_s / SecondsPerDay ! Physics time step converted to [day] as sinking velocity is in [m/day]

  if (VPhy .gt. 0.1) then
    call REcoM_sinking(dt_sink,Nn,SinkVel(:,ivphy),thick,recipthick,state(:,iphyn),sink,zF, mesh)
    sms(1:Nn,iphyn) = sms(1:Nn,iphyn) + sink(1:Nn)

    call REcoM_sinking(dt_sink,Nn,SinkVel(:,ivphy),thick,recipthick,state(:,iphyc),sink,zF, mesh)
    sms(1:Nn,iphyc) = sms(1:Nn,iphyc) + sink(1:Nn)	

    call REcoM_sinking(dt_sink,Nn,SinkVel(:,ivphy),thick,recipthick,state(:,ipchl),sink,zF, mesh)
    sms(1:Nn,ipchl) = sms(1:Nn,ipchl) + sink(1:Nn)

!#ifdef REcoM_calcification
    call REcoM_sinking(dt_sink,Nn,SinkVel(:,ivphy),thick,recipthick,state(:,iphycal),sink,zF, mesh)
    sms(1:Nn,iphycal) = sms(1:Nn,iphycal) + sink(1:Nn)

!#endif

  end if		

  if (VDia .gt. 0.1) then
    call REcoM_sinking(dt_sink,Nn,SinkVel(:,ivdia),thick,recipthick,state(:,idian),sink,zF, mesh)
    sms(1:Nn,idian) = sms(1:Nn,idian) + sink(1:Nn)

    call REcoM_sinking(dt_sink,Nn,SinkVel(:,ivdia),thick,recipthick,state(:,idiac),sink,zF, mesh)
    sms(1:Nn,idiac) = sms(1:Nn,idiac) + sink(1:Nn)

    call REcoM_sinking(dt_sink,Nn,SinkVel(:,ivdia),thick,recipthick,state(:,idchl),sink,zF, mesh)
    sms(1:Nn,idchl) = sms(1:Nn,idchl) + sink(1:Nn)

    call REcoM_sinking(dt_sink,Nn,SinkVel(:,ivdia),thick,recipthick,state(:,idiasi),sink,zF, mesh)
    sms(1:Nn,idiasi) = sms(1:Nn,idiasi) + sink(1:Nn)
!
  end if		

  if (VDet .gt. 0.1) then
    call REcoM_sinking(dt_sink,Nn,SinkVel(:,ivdet),thick,recipthick,state(:,idetn),sink,zF, mesh)
    sms(1:Nn,idetn) = sms(1:Nn,idetn) + sink(1:Nn)

    call REcoM_sinking(dt_sink,Nn,SinkVel(:,ivdet),thick,recipthick,state(:,idetc),sink,zF, mesh)
    sms(1:Nn,idetc) = sms(1:Nn,idetc) + sink(1:Nn)	

    call REcoM_sinking(dt_sink,Nn,SinkVel(:,ivdet),thick,recipthick,state(:,idetsi),sink,zF, mesh)
    sms(1:Nn,idetsi) = sms(1:Nn,idetsi) + sink(1:Nn)

!#ifdef REcoM_calcification
    call REcoM_sinking(dt_sink,Nn,SinkVel(:,ivdet),thick,recipthick,state(:,idetcal),sink,zF, mesh)
    sms(1:Nn,idetcal) = sms(1:Nn,idetcal) + sink(1:Nn)
!#endif


   if (REcoM_Second_Zoo) then
    call REcoM_sinking2(dt_sink,Nn,SinkVel(:,ivdetsc),thick,recipthick,state(:,idetz2n),sink,zF, mesh)
    sms(1:Nn,idetz2n) = sms(1:Nn,idetz2n) + sink(1:Nn)

    call REcoM_sinking2(dt_sink,Nn,SinkVel(:,ivdetsc),thick,recipthick,state(:,idetz2c),sink,zF, mesh)
    sms(1:Nn,idetz2c) = sms(1:Nn,idetz2c) + sink(1:Nn)

    call REcoM_sinking2(dt_sink,Nn,SinkVel(:,ivdetsc),thick,recipthick,state(:,idetz2si),sink,zF, mesh)
    sms(1:Nn,idetz2si) = sms(1:Nn,idetz2si) + sink(1:Nn)

    call REcoM_sinking2(dt_sink,Nn,SinkVel(:,ivdetsc),thick,recipthick,state(:,idetz2calc),sink,zF, mesh)
    sms(1:Nn,idetz2calc) = sms(1:Nn,idetz2calc) + sink(1:Nn)
   end if
  
  end if		
  if (ciso) then
!   sinking of 13|14C compounds
    if (VPhy .gt. 0.1) then
        call REcoM_sinking(dt_sink,Nn,SinkVel(:,ivphy),thick,recipthick,state(:,iphyc_13),sink,zF,mesh)
        sms(:,iphyc_13) = sms(:,iphyc_13) + sink

        call REcoM_sinking(dt_sink,Nn,SinkVel(:,ivphy),thick,recipthick,state(:,iphyc_14),sink,zF,mesh)
        sms(:,iphyc_14) = sms(:,iphyc_14) + sink

        call REcoM_sinking(dt_sink,Nn,SinkVel(:,ivphy),thick,recipthick,state(:,iphycal_13),sink,zF,mesh)
        sms(:,iphycal_13) = sms(:,iphycal_13) + sink

        call REcoM_sinking(dt_sink,Nn,SinkVel(:,ivphy),thick,recipthick,state(:,iphycal_14),sink,zF,mesh)
        sms(:,iphycal_14) = sms(:,iphycal_14) + sink
    end if

    if (VDia .gt. 0.1) then
        call REcoM_sinking(dt_sink,Nn,SinkVel(:,ivdia),thick,recipthick,state(:,idiac_13),sink,zF,mesh)
        sms(:,idiac_13) = sms(:,idiac_13) + sink

        call REcoM_sinking(dt_sink,Nn,SinkVel(:,ivdia),thick,recipthick,state(:,idiac_14),sink,zF,mesh)
        sms(:,idiac_14) = sms(:,idiac_14) + sink
    end if

    if (VDet .gt. 0.1) then
        call REcoM_sinking(dt_sink,Nn,SinkVel(:,ivdet),thick,recipthick,state(:,idetc_13),sink,zF,mesh)
        sms(:,idetc_13) = sms(:,idetc_13) + sink

        call REcoM_sinking(dt_sink,Nn,SinkVel(:,ivdet),thick,recipthick,state(:,idetc_14),sink,zF,mesh)
        sms(:,idetc_14) = sms(:,idetc_14) + sink

        call REcoM_sinking(dt_sink,Nn,SinkVel(:,ivdet),thick,recipthick,state(:,idetcal_13),sink,zF,mesh)
        sms(:,idetcal_13) = sms(:,idetcal_13) + sink

        call REcoM_sinking(dt_sink,Nn,SinkVel(:,ivdet),thick,recipthick,state(:,idetcal_14),sink,zF,mesh)
        sms(:,idetcal_14) = sms(:,idetcal_14) + sink
    end if
  end if ! ciso

!-------------------------------------------------------------------------------
! Benthic layer: Detritus and phytoplankton sink into the benthic layer and are 
! lost from the water column
! (but remineralized and re-released in dissolved form later

!   if (mype==0) then
!      write(*,*) '____________________________________________________________'
!      write(*,*) ' --> Benthic layer: Detritus and phytoplankton sink'
!      write(*,*) '     into the bottom layer'
!   endif
 if (REcoM_Second_Zoo) then
  if (allow_var_sinking) then
    Vben_det = Vdet_a * abs(zF(Nn)) + VDet
    Vben_det_seczoo = VDet_zoo2
    Vben_phy = Vdet_a * abs(zF(Nn)) + VPhy
    Vben_dia = Vdet_a * abs(zF(Nn)) + VDia
  else
    Vben_det = VDet
    Vben_det_seczoo = VDet_zoo2
    Vben_phy = VPhy
    Vben_dia = VDia
  endif
 else 
  if (allow_var_sinking) then
    Vben_det = Vdet_a * abs(zF(Nn)) + VDet
    Vben_phy = Vdet_a * abs(zF(Nn)) + VPhy
    Vben_dia = Vdet_a * abs(zF(Nn)) + VDia
  else
    Vben_det = VDet
    Vben_phy = VPhy
    Vben_dia = VDia
  endif
 endif


!*** Det *** (index: 1=N,2=C,3=Si,4=calc)
    wFluxDet(1)    = Vben_det * state(Nn,idetn)
    sms(Nn,idetn)  = sms(Nn,idetn)    - wFluxDet(1) * dt_b  &
                     / bottom_node_thickness(n) ! * recipthick(Nn)
    LocBenthos(1)  = LocBenthos(1)    + wFluxDet(1) * dt_b ! / thickness of benthos layer

    wFluxDet(2)    = Vben_det * state(Nn,idetc)
    Fc             = max(tiny,0.1d0 * wFluxDet(2))                                     ! Conversion of detC from [mmolC/m2/day] => [umol/cm2/day]
    if (NitrogenSS) then
      LocDenit       = -0.9543 + 0.7662 * log(Fc) - 0.2350 * log(Fc) * log(Fc) ! Denitrification factor following Middelburg et al. (1996)
      LocDenit       = exp(LocDenit * log(10.d0))                                 ! = 10^Denit, following Middelburg et al. (1996)
      LocDenit       = q_NC_Denit * LocDenit * wFluxDet(2) ! q_NC_Denit = 0.86 N:C quota of the denitrification process
      sms(Nn,idin)   = sms(Nn,idin)    -   LocDenit * dt_b * recipthick(Nn)   ! Part of benthos din is removed and lost due to denitrification
    else
      LocDenit = zero
    end if                     

    sms(Nn,idetc)  = sms(Nn,idetc)    - wFluxDet(2) * dt_b  &
                     / bottom_node_thickness(n) ! * recipthick(Nn)
    LocBenthos(2)  = LocBenthos(2)    + wFluxDet(2) * dt_b ! / thickness of benthos layer

    wFluxDet(3)    = Vben_det * state(Nn,idetsi)                        ! m/day*[mmol/m3 ] 
    sms(Nn,idetsi) = sms(Nn,idetsi)   - wFluxDet(3) * dt_b  &
                      / bottom_node_thickness(n) ! * recipthick(Nn)
    LocBenthos(3)  = LocBenthos(3)    + wFluxDet(3) * dt_b ! / thickness of benthos layer   [mmol/day/m2]*[day] 
	  
!#ifdef REcoM_calcification 
    wFluxDet(4)     = Vben_det * state(Nn,idetcal)
    sms(Nn,idetcal) = sms(Nn,idetcal) - wFluxDet(4) * dt_b &
                      / bottom_node_thickness(n) ! * recipthick(Nn) 
    LocBenthos(4)   = LocBenthos(4)   + wFluxDet(4) * dt_b ! / thickness of benthos layer
!#endif

!*** Phy *** (index: 1=N,2=C,3=calc)
    wFluxPhy(1)     = Vben_phy * state(Nn,iphyn)
    sms(Nn,iphyn)   = sms(Nn,iphyn)   - wFluxPhy(1) * dt_b  &
                      / bottom_node_thickness(n) ! * recipthick(Nn)
    LocBenthos(1)   = LocBenthos(1)   + wFluxPhy(1) * dt_b ! / thickness of benthos layer

    wFluxPhy(2)     = Vben_phy * state(Nn,iphyc)
    sms(Nn,iphyc)   = sms(Nn,iphyc)   - wFluxPhy(2) * dt_b  &
                      / bottom_node_thickness(n) ! * recipthick(Nn)
    LocBenthos(2)   = LocBenthos(2)   + wFluxPhy(2) * dt_b ! / thickness of benthos layer
!#ifdef REcoM_calcification
    wFluxPhy(3)     = Vben_phy * state(Nn,iphycal)
    sms(Nn,iphycal) = sms(Nn,iphycal) - wFluxPhy(3) * dt_b  &
                      / bottom_node_thickness(n) ! * recipthick(Nn)
    LocBenthos(4)   = LocBenthos(4)   + wFluxPhy(3) * dt_b ! / thickness of benthos layer
!#endif

    wFluxPhy(4)     = Vben_phy * state(Nn,ipchl)
    sms(Nn,ipchl)   = sms(Nn,ipchl)   - wFluxPhy(4) * dt_b  &
     	              / bottom_node_thickness(n) !* recipthick(Nn)

!------------------------------------
if (REcoM_Second_Zoo) then
    wFluxDet(5)    = Vben_det_seczoo * state(Nn,idetz2n)
    sms(Nn,idetz2n)  = sms(Nn,idetz2n)    - wFluxDet(5) * dt_b  &
                      / bottom_node_thickness(n) ! * recipthick(Nn)
    LocBenthos(1)  = LocBenthos(1)    + wFluxDet(5) * dt_b ! / thickness of benthos layer  

    wFluxDet(6)    = Vben_det_seczoo * state(Nn,idetz2c)
    sms(Nn,idetz2c)  = sms(Nn,idetz2c)    - wFluxDet(6) * dt_b  &
                      / bottom_node_thickness(n) ! * recipthick(Nn)
    LocBenthos(2)  = LocBenthos(2)    + wFluxDet(6) * dt_b ! / thickness of benthos layer  

    wFluxDet(7)    = Vben_det_seczoo * state(Nn,idetz2si)
    sms(Nn,idetz2si)  = sms(Nn,idetz2si)    - wFluxDet(7) * dt_b  &
                      / bottom_node_thickness(n) ! * recipthick(Nn)
    LocBenthos(3)  = LocBenthos(3)    + wFluxDet(7) * dt_b ! / thickness of benthos layer  

    wFluxDet(8)    = Vben_det_seczoo * state(Nn,idetz2calc)
    sms(Nn,idetz2calc)  = sms(Nn,idetz2calc)    - wFluxDet(8) * dt_b  &
                      / bottom_node_thickness(n) ! * recipthick(Nn)
    LocBenthos(4)  = LocBenthos(4)    + wFluxDet(8) * dt_b ! / thickness of benthos layer  
endif

!*** Dia *** (index: 1=N,2=C,3=Si)
    wFluxDia(1)     = Vben_dia * state(Nn,idian)
    sms(Nn,idian)   = sms(Nn,idian)   - wFluxDia(1) * dt_b  &
                      / bottom_node_thickness(n) !* recipthick(Nn)
    LocBenthos(1)   = LocBenthos(1)   + wFluxDia(1) * dt_b ! / thickness of benthos layer
 
    wFluxDia(2)     = Vben_dia * state(Nn,idiac)
    sms(Nn,idiac)   = sms(Nn,idiac)   - wFluxDia(2) * dt_b  &
                      / bottom_node_thickness(n) !* recipthick(Nn) 
    LocBenthos(2)   = LocBenthos(2)   + wFluxDia(2) * dt_b ! / thickness of benthos layer

    wFluxDia(3)     = Vben_dia * state(Nn,idiasi)
    sms(Nn,idiasi)  = sms(Nn,idiasi)  - wFluxDia(3) * dt_b  &
                      / bottom_node_thickness(n) !* recipthick(Nn)
    LocBenthos(3)   = LocBenthos(3)   + wFluxDia(3) * dt_b ! / thickness of benthos layer
    
    wFluxDia(4)     = Vben_dia * state(Nn,idchl)
    sms(Nn,idchl)   = sms(Nn,idchl)   - wFluxDia(4) * dt_b  &
     	              / bottom_node_thickness(n) !* recipthick(Nn)


 
    if (ciso) then ! reorganize the index here after secondzoopnakton and detritus
!*** Det *** (index: 5=C_13,6=C_14,7=calc_13,8=calc_14)
        wFluxDet(5)        = Vben_det * state(Nn,idetc_13)
        sms(Nn,idetc_13)   = sms(Nn,idetc_13) - wFluxDet(5) * dt_b * recipthick(Nn)
        LocBenthos(5)      = LocBenthos(5)    + wFluxDet(5) * dt_b ! / thickness of benthos layer

        wFluxDet(6)        = Vben_det * state(Nn,idetc_14)
        sms(Nn,idetc_14)   = sms(Nn,idetc_14) - wFluxDet(6) * dt_b * recipthick(Nn)
        LocBenthos(6)      = LocBenthos(6)    + wFluxDet(6) * dt_b ! / thickness of benthos layer

        wFluxDet(7)        = Vben_det * state(Nn,idetcal_13)
        sms(Nn,idetcal_13) = sms(Nn,idetcal_13) - wFluxDet(7) * dt_b * recipthick(Nn)
        LocBenthos(7)      = LocBenthos(7)   + wFluxDet(7) * dt_b ! / thickness of benthos layer

        wFluxDet(8)        = Vben_det * state(Nn,idetcal_14)
        sms(Nn,idetcal_14) = sms(Nn,idetcal_14) - wFluxDet(8) * dt_b * recipthick(Nn)
        LocBenthos(8)      = LocBenthos(8)   + wFluxDet(8) * dt_b ! / thickness of benthos layer

!*** Phy *** (index: 5=C_13,6=C_14,7=calc_13,8=calc_14)
        wFluxPhy(5)        = Vben_phy * state(Nn,iphyc_13)
        sms(Nn,iphyc_13)   = sms(Nn,iphyc_13) - wFluxPhy(5) * dt_b * recipthick(Nn)
        LocBenthos(5)      = LocBenthos(5)   + wFluxPhy(5) * dt_b ! / thickness of benthos layer

        wFluxPhy(6)        = Vben_phy * state(Nn,iphyc_14)
        sms(Nn,iphyc_14)   = sms(Nn,iphyc_14) - wFluxPhy(6) * dt_b * recipthick(Nn)
        LocBenthos(6)      = LocBenthos(6)   + wFluxPhy(6) * dt_b ! / thickness of benthos layer

        wFluxPhy(7)        = Vben_phy * state(Nn,iphycal_13)
        sms(Nn,iphycal_13) = sms(Nn,iphycal_13) - wFluxPhy(7) * dt_b * recipthick(Nn)
        LocBenthos(7)      = LocBenthos(7)   + wFluxPhy(7) * dt_b ! / thickness of benthos layer

        wFluxPhy(8)        = Vben_phy * state(Nn,iphycal_14)
        sms(Nn,iphycal_14) = sms(Nn,iphycal_14) - wFluxPhy(8) * dt_b * recipthick(Nn)
        LocBenthos(8)      = LocBenthos(8)   + wFluxPhy(8) * dt_b ! / thickness of benthos layer

!*** Dia *** (index: 5=C_13,6=C_14)
        wFluxDia(5)        = Vben_dia * state(Nn,idiac_13)
        sms(Nn,idiac_13)   = sms(Nn,idiac_13) - wFluxDia(5) * dt_b * recipthick(Nn)
        LocBenthos(5)      = LocBenthos(5)   + wFluxDia(5) * dt_b ! / thickness of benthos layer

        wFluxDia(6)        = Vben_dia * state(Nn,idiac_14)
        sms(Nn,idiac_14)   = sms(Nn,idiac_14) - wFluxDia(6) * dt_b * recipthick(Nn)
        LocBenthos(6)      = LocBenthos(6)   + wFluxDia(6) * dt_b ! / thickness of benthos layer
    end if ! ciso


end subroutine REcoM_sms

!-------------------------------------------------------------------------------
! Function for calculating limiter
!-------------------------------------------------------------------------------

function recom_limiter(slope,qa,qb)
  use recom_config
  Implicit None
  Real(kind=8) :: recom_limiter
  Real(kind=8) :: slope, qa, qb
  Real(kind=8) :: dq
	
  dq = qa - qb
  if (REcoM_Geider_limiter) then
    recom_limiter = max(min( -slope*dq, 1.d0),0.d0)
  else
    recom_limiter = 1.d0 - exp( -slope*( abs(dq)-dq )**2)
  endif
  return
  end

!-------------------------------------------------------------------------------
! Function for iron chemistry
!-------------------------------------------------------------------------------

function iron_chemistry(Fe, totalLigand, ligandStabConst)
  implicit none

  Real(kind=8) :: iron_chemistry
  Real(kind=8) :: Fe, totalLigand, ligandStabConst ! Input
  Real(kind=8) :: FreeFe                          ! Output
  Real(kind=8) :: ligand,FeL,a,b,c,discrim

! Abbrevations
  a = ligandstabConst
  b = ligandstabConst * (Fe - totalLigand) + 1.d0
  c = -totalLigand
  discrim = b*b - 4.d0 * a * c
	
  if (a .ne. 0.d0 .and. discrim .ge. 0.d0) then
    ligand = ( -b + sqrt(discrim) ) / (2.d0 * a)
    FeL    = totalLigand - ligand
    freeFe = Fe - FeL
  else ! No free iron
    freeFe = 0.d0
  end if

  iron_chemistry = freeFe

  return
  end

!===============================================================================
! Subroutine calculating sinking of Detritus, diatoms and possibly phy
!===============================================================================
subroutine REcoM_sinking(dt, Nn, wF, dzF, recipDzF, C, sink, zF, mesh)

  use recom_config
  use g_PARSUP
  use mod_MESH
 
  implicit none

  type(t_mesh), intent(in), target  :: mesh

  Real(kind=8)                      :: dt
  Integer                           :: Nn

  Integer	                    :: k,km1,km2,kp1

  Real(kind=8),dimension(mesh%nl-1) :: C
  Real(kind=8),dimension(mesh%nl)   :: wF            ! Vertical velocity for fluxes

  Real(kind=8)                      :: wLoc,wM,wPs
  Real(kind=8)                      :: Rjp,Rj,Rjm

  Real(kind=8),dimension(mesh%nl-1) :: dzF,recipdzF

  Real(kind=8)                      :: cfl, d0, d1, thetaP, thetaM, psiP, psiM
  Real(kind=8)                      :: onesixth	= 	1.d0/6.d0
  Real(kind=8)                      :: wflux
  Real(kind=8)                      :: wfluxkp1

  Real(kind=8),dimension(mesh%nl-1) :: sink	

  Real(kind=8),dimension(mesh%nl)   :: zF                   ! [m] Depth of fluxes
  Real(kind=8),dimension(mesh%nl)   :: new_wF

#include  "../associate_mesh.h"

!! zbar (depth of levels) and Z (mid depths of layers)
!! zbar is negative 
!! zbar(mesh%nl) allocate the array for storing the standard depths
!! Z(mesh%nl-1)  mid-depths of cells

    !    ----------- zbar_1, V_1 (Volume eq. to Area)
    ! Z_1 o T_1
    !    ----------- zbar_2, V_2
    ! Z_2 o T_2
    !    ----------- zbar_3, V_3
    ! Z_3 o T_3
    !    ----------- zbar_4


!-------------------------------------------------------------------------------
! 
!-------------------------------------------------------------------------------

  if (allow_var_sinking) then

!! loop over the layers
!write (*,*) "nlevels_nod2D= ", nlevels_nod2D(n), " Nn= ", Nn

    do k=1,Nn+1


      ! Vdet_a = 0.0288       ! [1/day]
      new_wF(k)= Vdet_a * abs(zF(k)) + wF(k) ! Applies to phy, dia and det.
!!                |-----------------|
!!                       |                | 
!!                       |                '--> ! constant through the water column and positive downwards ((calculated in recom_extra)
!!                       '-->	increasing sinking velocity with depth [m/day] (calculated in recom_extra)     
    enddo

  else
!! constant through the water column positive downwards ((calculated in recom_extra)
    new_wF = wF 

  endif

!-------------------------------------------------------------------------------
! 
!-------------------------------------------------------------------------------

  
  wfluxkp1   = 0.d0
  wflux      = 0.d0	
	
! mid layers
 
  do k = Nn,2,-1   

     !___________________________________________________________________
            ! if topographic depth dd is deeper than depth of deepest full cell 
            ! depth level zbar(nle)
            !       : 
            !       : 
            ! ______________ zbar(Nn-4)--------->+---->+
            !                                     |     |
            !                                     |     |
            ! -------------- Z(Nn-4)              |     |
            !                                     |     |
            !                                     |     |
            ! ______________ zbar(Nn-3)--------->+---->+   --> km2
            !                                     |     |            |
            !                                     |     |            |
            ! -------------- Z(Nn-3)              |     |            Rjm 
            !                                     |     |            |
            !                                     |     |            |
            ! ______________ zbar(Nn-2)--------->+---->+   --> km1
            !                                     |     |            |
            !                                     |     |            |
            ! -------------- Z(Nn-2)              |     |            Rj
            !                                     |     |            |
            !                                     |     |            |
            ! ______________ zbar(Nn-1)--------->+---->+   --> k       
            !                                     |     |            |
            !                                     |     |            |
            ! -------------- Z(Nn-1)              |     |            Rjp
            !                                     |     |            |
            !                                     |     |            |
            ! ______________ zbar(Nn)----------->+---->+   --> kp1
            ! / / / / / / / / / / / / / / / / / / 
           


    ! -- lower layer  
    km1      = max(1,  k-1) 
    ! -- upper layer
    km2      = max(1,  k-2) 

    kp1      = min(Nn, k+1) ! -- kp1 --> maximum depth Nn 

    wLoc     = -new_wF(k)  ! new_wF(k) is positive downwards therefore wLoc is negative

    wPs      = wLoc + abs(wLoc) ! --> Positive vertical velocity
    wM       = wLoc - abs(wLoc) ! --> Negative vertical velocity

    Rjp	     = C(k)	- C(kp1) 

    Rj       = C(km1)   - C(k)  

    Rjm	     = C(km2)   - C(km1)	

    cfl      = abs(wLoc * dt * recipDzF(k))         ! [m/day] * [day] * [1/m]

    d0       = (2.d0 - cfl)*(1.d0 - cfl)*onesixth

    d1       = (1.d0 - cfl*cfl)*onesixth
	
    thetaP   = Rjm/(1.d-20+Rj)

    psiP     = d0 + d1*thetaP
    psiP     = max(0.d0, min(min(1.d0,psiP), &
                 (1.d0-cfl)/(1.d-20+cfl)*thetaP))

    thetaM   = Rjp/(1.d-20 + Rj)	
    psiM     = d0 + d1*thetaM

    psiM     = max(0.d0, min(min(1.d0,psiM), &
		(1.d0-cfl)/(1.d-20-cfl)*thetaM))

    wflux    = ( 0.5*wPs*(C(k)  + psiM * Rj)+ &
		0.5*wM*(C(km1)+ psiP * Rj))

    sink(k)  = -(wflux-wfluxkp1)*recipDzF(k)*dt

    wfluxkp1 =	wflux	

  enddo

! surface layer 

  k          = 1
  wflux      = 0
  sink(k)    = -(wflux-wfluxkp1)*recipDzF(k)*dt
	
end subroutine REcoM_sinking
!===============================================================================
! Subroutine calculating sinking2 of secondzoo
!===============================================================================
subroutine REcoM_sinking2(dt, Nn, wF, dzF, recipDzF, C, sink, zF, mesh)

  use recom_config
  use g_PARSUP
  use mod_MESH
 
  implicit none

  type(t_mesh), intent(in), target  :: mesh

  Real(kind=8)                      :: dt
  Integer                           :: Nn

  Integer	                    :: k,km1,km2,kp1

  Real(kind=8),dimension(mesh%nl-1) :: C
  Real(kind=8),dimension(mesh%nl)   :: wF            ! Vertical velocity for fluxes

  Real(kind=8)                      :: wLoc,wM,wPs
  Real(kind=8)                      :: Rjp,Rj,Rjm

  Real(kind=8),dimension(mesh%nl-1) :: dzF,recipdzF

  Real(kind=8)                      :: cfl, d0, d1, thetaP, thetaM, psiP, psiM
  Real(kind=8)                      :: onesixth	= 	1.d0/6.d0
  Real(kind=8)                      :: wflux
  Real(kind=8)                      :: wfluxkp1

  Real(kind=8),dimension(mesh%nl-1) :: sink	

  Real(kind=8),dimension(mesh%nl)   :: zF                   ! [m] Depth of fluxes
  Real(kind=8),dimension(mesh%nl)   :: new_wF

#include  "../associate_mesh.h"


  new_wF = wF 

!-------------------------------------------------------------------------------
! 
!-------------------------------------------------------------------------------

  
  wfluxkp1   = 0.d0
  wflux      = 0.d0	
	
! mid layers
 
  do k = Nn,2,-1   

     !___________________________________________________________________
            ! if topographic depth dd is deeper than depth of deepest full cell 
            ! depth level zbar(nle)
            !       : 
            !       : 
            ! ______________ zbar(Nn-4)--------->+---->+
            !                                     |     |
            !                                     |     |
            ! -------------- Z(Nn-4)              |     |
            !                                     |     |
            !                                     |     |
            ! ______________ zbar(Nn-3)--------->+---->+   --> km2
            !                                     |     |            |
            !                                     |     |            |
            ! -------------- Z(Nn-3)              |     |            Rjm 
            !                                     |     |            |
            !                                     |     |            |
            ! ______________ zbar(Nn-2)--------->+---->+   --> km1
            !                                     |     |            |
            !                                     |     |            |
            ! -------------- Z(Nn-2)              |     |            Rj
            !                                     |     |            |
            !                                     |     |            |
            ! ______________ zbar(Nn-1)--------->+---->+   --> k       
            !                                     |     |            |
            !                                     |     |            |
            ! -------------- Z(Nn-1)              |     |            Rjp
            !                                     |     |            |
            !                                     |     |            |
            ! ______________ zbar(Nn)----------->+---->+   --> kp1
            ! / / / / / / / / / / / / / / / / / / 
           


    ! -- lower layer  
    km1      = max(1,  k-1) 
    ! -- upper layer
    km2      = max(1,  k-2) 

    kp1      = min(Nn, k+1) ! -- kp1 --> maximum depth Nn 

    wLoc     = -new_wF(k)  ! new_wF(k) is positive downwards therefore wLoc is negative

    wPs      = wLoc + abs(wLoc) ! --> Positive vertical velocity
    wM       = wLoc - abs(wLoc) ! --> Negative vertical velocity

    Rjp	     = C(k)	- C(kp1) 

    Rj       = C(km1)   - C(k)  

    Rjm	     = C(km2)   - C(km1)	

    cfl      = abs(wLoc * dt * recipDzF(k))         ! [m/day] * [day] * [1/m]

    d0       = (2.d0 - cfl)*(1.d0 - cfl)*onesixth

    d1       = (1.d0 - cfl*cfl)*onesixth
	
    thetaP   = Rjm/(1.d-20+Rj)

    psiP     = d0 + d1*thetaP
    psiP     = max(0.d0, min(min(1.d0,psiP), &
                 (1.d0-cfl)/(1.d-20+cfl)*thetaP))

    thetaM   = Rjp/(1.d-20 + Rj)	
    psiM     = d0 + d1*thetaM

    psiM     = max(0.d0, min(min(1.d0,psiM), &
		(1.d0-cfl)/(1.d-20-cfl)*thetaM))

    wflux    = ( 0.5*wPs*(C(k)  + psiM * Rj)+ &
		0.5*wM*(C(km1)+ psiP * Rj))

    sink(k)  = -(wflux-wfluxkp1)*recipDzF(k)*dt

    wfluxkp1 =	wflux	

  enddo

! surface layer 

  k          = 1
  wflux      = 0
  sink(k)    = -(wflux-wfluxkp1)*recipDzF(k)*dt
	
end subroutine REcoM_sinking2
