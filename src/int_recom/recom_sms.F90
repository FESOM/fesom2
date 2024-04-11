subroutine REcoM_sms(n,Nn,state,thick,recipthick,SurfSR,sms,Temp, Sali_depth &
        , Loc_slp, zF, PAR, Lond, Latd, ice, dynamics, tracers, partit, mesh)

    use recom_declarations
    use recom_locvar
    use recom_glovar
    use recom_config
    use recoM_ciso
    use g_clock

    use g_config
    use MOD_MESH
    use MOD_TRACER
    use MOD_DYN
    USE MOD_ICE
    use o_ARRAYS
    use o_PARAM
    USE MOD_PARTIT
    USE MOD_PARSUP

    use g_forcing_arrays
    use g_comm_auto
    use mvars
    use mdepth2press                                   
    use gsw_mod_toolbox, only: gsw_sa_from_sp,gsw_ct_from_pt,gsw_rho

    implicit none
    type(t_dyn)   , intent(inout), target :: dynamics
    type(t_tracer), intent(inout), target :: tracers
    type(t_partit), intent(inout), target :: partit
    type(t_mesh)  , intent(inout), target :: mesh
    type(t_ice)   , intent(inout), target :: ice

    integer, intent(in)                                     :: Nn                   !< Total number of nodes in the vertical
    real(kind=8),dimension(mesh%nl-1,bgc_num),intent(inout) :: state                !< ChlA conc in phytoplankton [mg/m3]
									  	    !! should be in instead of inout

    real(kind=8),dimension(mesh%nl-1)                       :: thick                !< [m] Vertical distance between two nodes = Thickness 
    real(kind=8),dimension(mesh%nl-1)                       :: recipthick           !< [1/m] reciprocal of thick
    real(kind=8),intent(in)                                 :: SurfSR               !< [W/m2] ShortWave radiation at surface

    real(kind=8),dimension(mesh%nl-1,bgc_num),intent(inout) :: sms                  !< Source-Minus-Sinks term
    real(kind=8),dimension(mesh%nl-1)        ,intent(in)    :: Temp                 !< [degrees C] Ocean temperature
    real(kind=8),dimension(mesh%nl-1)        ,intent(in)    :: Sali_depth           !< NEW MOCSY Salinity for the whole water column

    real(kind=8),dimension(mesh%nl)          ,intent(in)    :: zF                   !< [m] Depth of fluxes
    real(kind=8),dimension(mesh%nl-1),intent(inout)         :: PAR

    real(kind=8)                                            :: dt_d                 !< Size of time steps [day]
    real(kind=8)                                            :: dt_b                 !< Size of time steps [day]
    real(kind=8),dimension(mesh%nl-1)                       :: Sink
    real(kind=8)                                            :: dt_sink              !< Size of local time step

    real(kind=8)                                            :: recip_hetN_plus      !< MB's addition to heterotrophic respiration
    real(kind=8)                                            :: recip_res_het        !< [day] Reciprocal of respiration by heterotrophs and mortality (loss to detritus)
    real(kind=8)                                            :: Sink_Vel
    real(kind=8)                                            :: aux
    integer                                                 :: k,step,ii, idiags,n

    real(kind=8),                      intent(in)           :: Loc_slp              ! MOCSY [Pa] sea-level pressure
    real(kind=8)                                            :: Patm_depth(1)        ! MOCSY
    real(kind=8)                                            :: REcoM_T_depth(1)     ! MOCSY temperature for the whole water column for mocsy minimum defined as -2
    real(kind=8)                                            :: REcoM_S_depth(1)     ! MOCSY
    real(kind=8)                                            :: REcoM_DIC_depth(1)   ! MOCSY
    real(kind=8)                                            :: REcoM_Alk_depth(1)   ! MOCSY
    real(kind=8)                                            :: REcoM_Si_depth(1)    ! MOCSY
    real(kind=8)                                            :: REcoM_Phos_depth(1)  ! MOCSY
    real(kind=8),                      intent(in)           :: Latd(1)              ! latitude in degree
    real(kind=8),                      intent(in)           :: Lond(1)              ! NEW MOCSY longitude in degree 
    real(kind=8)                                            :: mocsy_step_per_day 
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
        FreeFe,  &
        O2

#include "../associate_part_def.h"
#include "../associate_mesh_def.h"
#include "../associate_part_ass.h"
#include "../associate_mesh_ass.h"

    sms = zero ! double precision

    tiny_N   = tiny_chl/chl2N_max      !< 0.00001/ 3.15d0   Chl2N_max [mg CHL/mmol N] Maximum CHL a : N ratio = 0.3 gCHL gN^-1
    tiny_N_d = tiny_chl/chl2N_max_d    !< 0.00001/ 4.2d0

    tiny_C   = tiny_N  /NCmax          !< NCmax   = 0.2d0   [mmol N/mmol C] Maximum cell quota of nitrogen (N:C)
    tiny_C_d = tiny_N_d/NCmax_d        !< NCmax_d = 0.2d0 

    tiny_Si  = tiny_C_d/SiCmax         !< SiCmax = 0.8d0


    recip_res_het = 1.d0/res_het       !< res_het = 0.01d0  [1/day] Respiration by heterotrophs and mortality (loss to detritus)

    Patm_depth    = Loc_slp/Pa2atm     ! MOCSY convert from Pa to atm.

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

!-------------------------------------------------------------------------------
! Main vertical loop starts
        do k = one,Nn   ! nzmin, nzmax 
!    do n=1, myDim_nod2D!+eDim_nod2D 
!       Nn=nlevels_nod2D(n)-1  !nzmax
!       nzmin = ulevels_nod2D(row)
!       nzmax = nlevels_nod2D(row)
!       do k=1, Nn
            DIN    = max(tiny,state(k,idin)  	     + sms(k,idin  )) !< Avoids division by zero
            DIC    = max(tiny,state(k,idic)   	     + sms(k,idic  )) !! and updates Conc between
            ALK    = max(tiny,state(k,ialk)   	     + sms(k,ialk  )) !! local steps in REcoM when
            PhyN   = max(tiny_N,state(k,iphyn)       + sms(k,iphyn )) !! biostep > 1
            PhyC   = max(tiny_C,state(k,iphyc) 	     + sms(k,iphyc ))
            PhyChl = max(tiny_chl,state(k,ipchl)     + sms(k,ipchl ))
            DetN   = max(tiny,state(k,idetn)  	     + sms(k,idetn ))
            DetC   = max(tiny,state(k,idetc)  	     + sms(k,idetc ))
            HetN   = max(tiny,state(k,ihetn)  	     + sms(k,ihetn ))
            HetC   = max(tiny,state(k,ihetc)  	     + sms(k,ihetc ))
            DON    = max(tiny,state(k,idon)          + sms(k,idon  ))
            EOC    = max(tiny,state(k,idoc)   	     + sms(k,idoc  ))
            DiaN   = max(tiny_N_d,state(k,idian)     + sms(k,idian ))
            DiaC   = max(tiny_C_d,state(k,idiac)     + sms(k,idiac ))
            DiaChl = max(tiny_chl,state(k,idchl)     + sms(k,idchl ))
            DiaSi  = max(tiny_si,state(k,idiasi)     + sms(k,idiasi)) 
            DetSi  = max(tiny,state(k,idetsi) 	     + sms(k,idetsi)) 
            Si     = max(tiny,state(k,isi)    	     + sms(k,isi   ))
            Fe     = max(tiny,state(k,ife)    	     + sms(k,ife ))
            O2     = max(tiny,state(k,ioxy)          + sms(k,ioxy))
            FreeFe = zero

            PhyCalc = max(tiny,state(k,iphycal)      + sms(k,iphycal))
            DetCalc = max(tiny,state(k,idetcal)      + sms(k,idetcal))

!!------------------------------------------------------------------------------
!< Quotas
            ! *** Small phytoplankton
            quota       =  PhyN / PhyC                     ! include variability of the N: C ratio, cellular chemical composition 
            recipquota  =  real(one) / quota
            Chl2C       =  PhyChl  / PhyC                  ! Chl a:phytoplankton carbon ratio, cellular chemical composition [gCHL gC^-1]
            Chl2N       =  PhyChl  / PhyN                  ! Chl a:phytoplankton nitrogen ratio, cellular chemical composition [gCHL gN^-1]
            CHL2C_plast =  Chl2C * (quota/(quota - NCmin))

            ! *** Diatoms
            quota_dia       =  DiaN / DiaC
            recipQuota_dia  =  real(one)/quota_dia
            Chl2C_dia       =  DiaChl / DiaC
            Chl2N_dia       =  DiaChl / DiaN
            CHL2C_plast_dia =  Chl2C_dia * (quota_dia/(quota_dia - NCmin_d))
            qSiC            =  DiaSi / DiaC
            qSiN            =  DiaSi / DiaN

            recipQZoo       = HetC / HetN
            recip_hetN_plus = 1.d0 / (HetN + tiny_het) ! MB's addition for more stable zoo respiration
            if (Grazing_detritus) recipDet  = DetC / DetN

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

!< Silicate temperature dependence 
!            reminSiT = min(1.32e16 * exp(-11200.d0 * rTloc),reminSi) !! arrFunc control, reminSi=0.02d0 ! Kamatani (1982)
!            reminSiT = reminSi
            reminSiT = max(0.023d0 * 2.6d0**((Temp(k)-10.)/10.),reminSi)

!-------------------------------------------------------------------------------
!> O2 dependence of rates
!------------------------------------------------------------------------------- 
!! O2 dependency of organic matter remineralization
!! O2Func [0.0, 1.0]
!! k_o2_remin = 15.d0 mmol m-3; Table 1 in Cram 2018 cites 
!! DeVries & Weber 2017 for a range of 0-30 mmol m-3                

            O2Func = 1.d0 ! in this case, remin. rates only depend on temperature
            if (O2dep_remin) O2Func = O2/(k_o2_remin + O2) ! O2remin

!< *** Light ***
!< *************
!! Has to be calculated here already to use the 1%PAR depth.
            if (k==1) then 
                PARave     = max(tiny,SurfSR)
                PAR(k)     = PARave

                chl_upper  = (PhyChl + DiaChl)
            else
                chl_lower  = PhyChl + DiaChl
                Chlave     = (chl_upper+chl_lower)*0.5

                kappa      =  k_w + a_chl * (Chlave)
                kappastar  =  kappa / cosAI(n)
                kdzLower   =  kdzUpper + kappastar * thick(k-1)
                Lowerlight =  SurfSR * exp(-kdzLower)
                Lowerlight =  max(tiny,Lowerlight) 
                PARave     =  Lowerlight
                PAR(k)     =  PARave
                chl_upper  =  chl_lower
                kdzUpper   =  kdzLower
            end if

!------------------------------------------------------------------------------
! Calcite dissolution dependent on OmegaC ! DISS
!------------------------------------------------------------------------------
            Sink_Vel    = Vdet_a* abs(zF(k)) + Vdet

                calc_diss = calc_diss_rate * Sink_Vel/20.d0 ! Dissolution rate of CaCO3 scaled by the sinking velocity at the current depth
            calc_diss_ben = calc_diss_rate * Sink_Vel/20.d0 ! DISS added the variable calc_diss_ben to keep the calcite dissolution in the benthos with the old formulation

!-------------------------------------------------------------------------------
!> Photosynthesis section, light parameters and rates
!-------------------------------------------------------------------------------
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

!< *** Small phytoplankton ***
!< ***************************
            qlimitFac    = recom_limiter(NMinSlope, NCmin, quota) ! Eqn A55
            feLimitFac   = Fe/(k_Fe + Fe)                         ! Use Michaelis–Menten kinetics
            qlimitFac    = min(qlimitFac, feLimitFac)             ! Liebig law of the minimum
            pMax         = P_cm * qlimitFac * arrFunc             ! Maximum value of C-specific rate of photosynthesis
    
!< *** Diatoms ***
!< ***************
            qlimitFac    = recom_limiter(NMinSlope, NCmin_d, quota_dia)
            qlimitFacTmp = recom_limiter(SiMinSlope, SiCmin, qSiC)
            qlimitFac    = min(qLimitFac, qlimitFacTmp)
            feLimitFac   = Fe/(k_Fe_d + Fe)
            qlimitFac    = min(qlimitFac, feLimitFac)
            pMax_dia     = P_cm_d * qlimitFac * arrFunc
!-------------------------------------------------------------------------------
!< *** Small phytoplankton photosynthesis rate ***
!< ***********************************************
            if (pMax .lt. tiny .OR. PARave /= PARave .OR. CHL2C /= CHL2C) then  ! OG in case of only respiration, i.e. darkness??
                Cphot = zero
            else
                Cphot = pMax*(real(one) - exp(-alfa * Chl2C * PARave / pMax))
                if (CO2lim) Cphot = Cphot * PhyCO2 ! Added the CO2 dependence
            end if
            if (Cphot .lt. tiny) Cphot = zero

!< *** Diatom photosynthesis rate ***
!< **********************************
            if ( pMax_dia .lt. tiny .OR. PARave /= PARave .OR. CHL2C_dia /= CHL2C_dia) then
                Cphot_dia = zero
            else
                Cphot_dia = pMax_dia * (real(one) - exp(-alfa_d * Chl2C_dia * PARave / pMax_dia))
                if (CO2lim) Cphot_dia = Cphot_dia  * DiaCO2 ! Added the CO2 dependence
            end if
            if (Cphot_dia .lt. tiny) Cphot_dia = zero

!------------------------------------------------------------------------------- 
!< chlorophyll degradation
!-------------------------------------------------------------------------------
            KOchl = deg_Chl
            KOchl_dia = deg_Chl_d
        
            if (use_photodamage) then
!< add a minimum value for photodamage
!< *** Phytoplankton Chla loss ***
!< *******************************
                if (pMax .lt. tiny .OR. PARave /= PARave .OR. CHL2C_plast /= CHL2C_plast) then
                    KOchl = deg_Chl*0.1d0
                else
                    KOchl = deg_Chl*(real(one) - exp(-alfa * CHL2C_plast * PARave / pMax))
                    KOchl = max((deg_Chl*0.1d0), KOchl)
                end if
!< *** Diatoms Chla loss ***                                                                                                                                
!< *************************
                if (pMax_dia .lt. tiny .OR. PARave /= PARave .OR. CHL2C_plast_dia /= CHL2C_plast_dia) then
                    KOchl_dia = deg_Chl_d*0.1d0
                else
                    KOchl_dia = deg_Chl_d * (real(one) - exp(-alfa_d * CHL2C_plast_dia * PARave / pMax_dia ))
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

            V_cm      = V_cm_fact
            limitFacN = recom_limiter(NMaxSlope, quota, NCmax)
            N_assim   = V_cm * pMax * NCuptakeRatio &                ! [mmol N / (mmol C * day)]
                         * limitFacN * (DIN/(DIN + k_din))           ! Michaelis–Menten kinetics

            V_cm          = V_cm_fact_d
            limitFacN_dia = recom_limiter(NMaxSlope, quota_dia, NCmax_d)
            N_assim_dia   = V_cm * pMax_dia * NCUptakeRatio_d &
                             * limitFacN_dia * DIN/(DIN + k_din_d)


            limitFacSi     = recom_limiter(SiMaxSlope, qSiC, SiCmax)  &
                              * limitFacN_dia
            Si_assim       = V_cm_fact_d * P_cm_d * arrFunc * SiCUptakeRatio &
                              * limitFacSi * Si/(Si + k_si)

!-------------------------------------------------------------------------------
!< *** Iron chemistry ***
!< ********************** 
! select the method to calculate freeFe
               freeFe = iron_chemistry(Fe,totalligand,ligandStabConst)
!-------------------------------------------------------------------------------
!< *** Chlorophyll synthesis ***
!< *****************************

!< Coupled to N uptake
!< Converted to chlorophyll units with a maximum Chl:N ratio, Chl2N_max
!< Chl2N_max: Maximum Chl:N ratio for phytoplankton [mg Chl mmol N−1 ]

            chlSynth = zero
            if (PARave .ge. tiny .AND. PARave .eq. PARave) then
                chlSynth = N_assim * Chl2N_max                    &
                         * min(real(one),Cphot/(alfa * Chl2C * PARave))
            end if
            ChlSynth_dia = zero
            if (PARave .ge. tiny .AND. PARave .eq. PARave) then
                ChlSynth_dia = N_assim_dia * Chl2N_max_d                         &
                             * min(real(one),Cphot_dia /(alfa_d * Chl2C_dia * PARave))
            end if
            ChlSynth_cocco = zero
!-------------------------------------------------------------------------------
!< *** Phytoplankton respiraion rate ***
!< *************************************

!< res_phy: Maintenance respiration rate constant [day−1 ]
!< biosynth: The cost of biosynthesis of N [mmol C mmol N−1 ]

            phyRespRate     = res_phy   * limitFacN     + biosynth * N_assim
            phyRespRate_dia = res_phy_d * limitFacN_dia + biosynth * N_assim_dia + biosynthSi * Si_assim

!-------------------------------------------------------------------------------
! Mesozooplankton
!-------------------------------------------------------------------------------
!< Grazing on small phytoplankton, diatoms, coccolithophore (optional), 
!< microzooplankton (optional), slow- and fast-sinking detritus

!< *** Food availability ***
!< *************************
!< pzPhy: Maximum nanophytoplankton preference
!< pzDia: Maximum diatom preference
!< pzDet: Maximum slow-sinking detritus prefence by first zooplankton
!< pzDetZ2: Maximum fast-sinking detritus preference by first zooplankton

            if (REcoM_Grazing_Variable_Preference) then   ! CHECK ONUR
                aux = pzPhy*PhyN + pzDia*DiaN              
                if (Grazing_detritus) aux = aux + PzDet*DetN
! ******************************************************************************
                varpzPhy = (pzPhy*PhyN)/aux
                varpzDia = (pzDia*DiaN)/aux
                if (Grazing_detritus) varpzDet = (pzDet*DetN)/aux
! ******************************************************************************
                fDiaN = varpzDia * DiaN
                fPhyN = varpzPhy * PhyN
                if (Grazing_detritus) fDetN = varpzDet * DetN
            else ! REcoM_Grazing_Variable_Preference = .false.
                fPhyN = pzPhy * PhyN
                fDiaN = pzDia * DiaN
                if (Grazing_detritus) fDetN = pzDet * DetN
            end if ! REcoM_Grazing_Variable_Preference

!< *** Grazing fluxes ***
!< **********************
            food = fPhyN + fDiaN
            if (Grazing_detritus) food = food + fDetN
! ******************************************************************************
            foodsq            = food**2
            grazingFlux       = (Graz_max * foodsq)/(epsilonr + foodsq) * HetN * arrFunc
            grazingFlux_phy   = grazingFlux * fphyN / food
            grazingFlux_Dia   = grazingFlux * fDiaN / food
            if (Grazing_detritus) grazingFlux_Det   = grazingFlux * fDetN / food

!< *** Grazing efficiency ***
!< **************************
            grazEff = gfin + 1/(0.2*food + 2)

            grazingFluxcarbon_mes = (grazingFlux_phy   * recipQuota       * grazEff)   &
                                  + (grazingFlux_Dia   * recipQuota_Dia   * grazEff)   

            if (Grazing_detritus) grazingFluxcarbon_mes = grazingFluxcarbon_mes        &
                                  + (grazingFlux_Det   * recipDet         * grazEff)   


!-------------------------------------------------------------------------------
!< Heterotrophic respiration is assumed to drive zooplankton back to 
!< Redfield C:N if their C:N becomes higher than Redfield
!< res_het: Timescale for zooplankton respiration [day−1 ]

            if (het_resp_noredfield) then
                HetRespFlux = res_het * arrFunc * HetC ! tau * f_T [HetC]
            else
                HetRespFlux = recip_res_het * arrFunc * (hetC * recip_hetN_plus - redfield) * HetC
                HetRespFlux = max(zero, HetRespFlux)  !!!!!!!! CHECK Judith Valid for het_resp_noredfield case as well ???????? Then move it below
            endif


!-------------------------------------------------------------------------------
!< Zooplanton mortality (Quadratic)

            hetLossFlux = loss_het * HetN * HetN


!-------------------------------------------------------------------------------
! Phytoplankton and detritus aggregation
!-------------------------------------------------------------------------------
            if (diatom_mucus) then
                qlimitFac = recom_limiter(NMinSlope, NCmin_d, quota_dia)
                qlimitFacTmp = recom_limiter(SiMinSlope, SiCmin, qSiC)
                qlimitFac = min(qLimitFac, qlimitFacTmp)
                feLimitFac= Fe/(k_Fe_d + Fe)
                qlimitFac = min(qlimitFac, feLimitFac)
                aggregationrate = agg_PP * (1 - qlimitFac) * DiaN
            else
                aggregationrate = agg_PP * DiaN
            endif
                   
            aggregationrate = aggregationrate + agg_PD * DetN + agg_PP * PhyN
    
!-------------------------------------------------------------------------------
! Calcification
!-------------------------------------------------------------------------------
! Terms required for the formation and dissolution of CaCO3
! Without this, calcification is performed by a fraction of small phytoplankton

!< calc_prod_ratio: Calcite production ratio, dimensionless
            calcification = calc_prod_ratio * Cphot * PhyC ! Z in equations

            calc_loss_agg = aggregationrate * PhyCalc

!< *** Small phytoplankton ***
!< ***************************
            aux = recipQuota/(PhyC + tiny) * PhyCalc
            calc_loss_gra  = grazingFlux_phy  * aux
		
!-------------------------------------------------------------------------------
! Sources minus sinks (SMS)
!-------------------------------------------------------------------------------

!< *** DIN *** 
!< ***********

!< N_assim:        N assimilation rate for nanophytoplankton [mmolN mmolC-1 day-1]
!< PhyC:           Intracellular carbon concentration in nanophytoplankton [mmolC m-3]
!< N_assim_Dia:    N assimilation rate for diatoms [mmolN mmolC-1 day-1]
!< DiaC:           Intracellular carbon concentration in diatoms [mmolC m-3]
!< N_assim_Cocco:  N assimilation rate for coccolithophore [mmolN mmolC-1 day-1]
!< CoccoC:         Intracellular carbon concentration in coccolithophore [mmolC m-3]
!< rho_N*arrFunc:  Remineralization rate and temperature dependency which is calculated with arrFunc [day^-1]
!< O2Func:         O2 dependency of organic matter remineralization
!< DON:            Extracellular dissolved organic nitrogen [mmolN m-3]
!< dt_b:           REcoM time step [day]

!! Schourup 2013 Eq. A2

    sms(k,idin)      = (                     &
        - N_assim                  * PhyC    &  ! --> N assimilation Nanophytoplankton, [mmol N/(mmol C * day)] C specific N utilization rate
        - N_assim_Dia              * DiaC    &  ! --> N assimilation Diatoms
        + rho_N * arrFunc * O2Func * DON     &  ! --> DON remineralization, temperature dependent [day^-1 * mmol/m3] ! O2remin
                                            ) * dt_b + sms(k,idin)  

!< *** DIC ***
!< ***********

!< rho_C1: Temperature dependent C degradation of extracellular organic C (EOC) [day^-1]

    sms(k,idic)      = (                              &
        - Cphot                           * PhyC      & ! --> Small pyhtoplankton photosynthesis 
        + phyRespRate                     * PhyC      & ! --> Small pyhtoplankton respiration 
        - Cphot_Dia                       * DiaC      & ! --> Diatom photosynthesis 
        + phyRespRate_Dia                 * DiaC      & ! --> Diatom respiration
        + rho_C1 * arrFunc * O2Func       * EOC       & ! --> Remineralization of DOC ! NEW O2remin
        + HetRespFlux                                 & ! --> Mesozooplankton respiration
        + calc_diss                       * DetCalc   & ! --> Calcite dissolution from slow-sinking detritus 
        + calc_loss_gra  * calc_diss_guts             & ! --> Additional dissolution in mesozooplankton guts
        - calcification                               & ! --> Calcification
                                                     ) * dt_b + sms(k,idic)

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

!< *** Alk ***
!< ***********

!< Alkalinity (Assumes that N:P follows a constant Redfield ratio
!< N_assimC: 1.0625 = 1/16 + 1

    sms(k,ialk)      = (                                   &
        + 1.0625 * N_assim                       * PhyC    &
        + 1.0625 * N_assim_Dia                   * DiaC    &
        - 1.0625 * rho_N * arrFunc * O2Func      * DON     & ! O2remin      
        + 2.d0 * calc_diss                       * DetCalc &
        + 2.d0 * calc_loss_gra  * calc_diss_guts           &
        - 2.d0 * calcification                             &                      
                                                          ) * dt_b + sms(k,ialk)
!< *** Small Phytoplankton ***
!< ***************************

!____________________________________________________________
!< Small phytoplankton N

!< lossN:  Phytoplankton loss of organic N compounds [day^-1]

    sms(k,iphyn)      = (                                  &
        + N_assim                                * PhyC    & ! --> N assimilation
        - lossN * limitFacN                      * PhyN    & ! --> DON excretion
        - aggregationRate                        * PhyN    & ! --> Aggregation loss
        - grazingFlux_phy                                  & ! --> Grazing loss
                                                          ) * dt_b + sms(k,iphyn)
!____________________________________________________________
!< Small phytoplankton C

!< lossC: Phytoplankton loss of carbon [day^-1]
!< When N : C ratio becomes too high, excretion of DOC is downregulated
!< by the limiter function limitFacN
!< aggregationRate transfers C to the detritus pool

    sms(k,iphyc)      = (                                  &
        + Cphot                                  * PhyC    & ! --> Photosynthesis ---->/
        - lossC * limitFacN                      * PhyC    & ! --> Excretion of DOC   / Net photosynthesis
        - phyRespRate                            * PhyC    & ! --> Respiration ----->/
        - aggregationRate                        * PhyC    & ! --> Aggregation loss          
        - grazingFlux_phy  * recipQuota                    & ! --> Grazing loss
                                                          ) * dt_b + sms(k,iphyc)
!____________________________________________________________
! Phytoplankton ChlA

!< Chl2N: Conversion factor from mmolN to mgChla 
!< Chl2N = PhyChl/PhyN

    sms(k,ipchl)       = (                                 &
        + chlSynth                               * PhyC    & ! --> Chl-a synthesis
        - KOchl                                  * PhyChl  & ! --> Degradation loss
        - aggregationRate                        * PhyChl  & ! --> Aggregation loss
     	- grazingFlux_phy  * Chl2N                         & ! --> Grazing loss
                                                          ) * dt_b + sms(k,ipchl)

!< *** Slow-sinking Detritus ***
!< *****************************

!____________________________________________________________
! Detritus N
    if (Grazing_detritus) then
        sms(k,idetn)       = (                             &
	    + grazingFlux_phy                              &
            - grazingFlux_phy   * grazEff                  &
            + grazingFlux_dia                              &
            - grazingFlux_dia   * grazEff                  &
            - grazingFlux_Det   * grazEff                  & ! Sloppy feeding is thought because of grazing flux multiplied with grazeff 
            - grazingFlux_Det2  * grazEff2                 &  !!!!!!!!!!CHECK
            + aggregationRate               * PhyN         &
            + aggregationRate               * DiaN         &
            + hetLossFlux                                  &
            - reminN * arrFunc * O2Func     * DetN         & ! O2remin
                                                          ) * dt_b + sms(k,idetn)
   else
        sms(k,idetn)       = (                             &
            + grazingFlux_phy                              &
            + grazingFlux_dia                              &
            - grazingFlux        * grazEff                 &
            + aggregationRate               * PhyN         &
            + aggregationRate               * DiaN         &
            + hetLossFlux                                  &
            - reminN * arrFunc * O2Func     * DetN         & ! O2remin
                                                          ) * dt_b + sms(k,idetn)
   end if

!____________________________________________________________
! Detritus C
    if (Grazing_detritus) then
        sms(k,idetc)       = (                                 &
            + grazingFlux_phy * recipQuota                     &
            - grazingFlux_phy * recipQuota         * grazEff   &
            + grazingFlux_Dia * recipQuota_Dia                 &
            - grazingFlux_Dia * recipQuota_Dia     * grazEff   &
            - grazingFlux_Det  * recipDet  * grazEff           &
            - grazingFlux_Det2 * recipDet2 * grazEff           & !!!!!! CHECK
            + aggregationRate              * phyC              &
            + aggregationRate              * DiaC              &
            + hetLossFlux      * recipQZoo                     &
            - reminC * arrFunc * O2Func    * DetC              & ! O2remin
                                                              ) * dt_b + sms(k,idetc)
   else
        sms(k,idetc)       = (                                 &
            + grazingFlux_phy   * recipQuota                   &
            - grazingFlux_phy   * recipQuota        * grazEff  &
            + grazingFlux_Dia   * recipQuota_Dia               &
            - grazingFlux_Dia   * recipQuota_Dia    * grazEff  &
            + aggregationRate                       * phyC     &
            + aggregationRate                       * DiaC     &
            + hetLossFlux       * recipQZoo                    &
            - reminC * arrFunc * O2Func    * DetC              & ! O2remin
                                                              ) * dt_b + sms(k,idetc)
   end if

!< *** Mesozooplankton ***
!< ***********************

!____________________________________________________________
!< Heterotrophic N
        sms(k,ihetn)       = (                                 &
    	    + grazingFlux      * grazEff                       & ! --> Grazing on phytoplankton -> okay, because of recipQuota
     	    - hetLossFlux                                      & ! --> Mortality
     	    - lossN_z                      * HetN              & ! --> Excretion of DON
                                                              ) * dt_b + sms(k,ihetn)  
!____________________________________________________________
!< Heterotrophic C
        if (Grazing_detritus) then
            sms(k,ihetc)      = (                                 &
     	        + grazingFlux_phy    * recipQuota       * grazEff & ! --> Grazing on small phytoplankton
     	        + grazingFlux_Dia    * recipQuota_Dia   * grazEff & ! --> Grazing on diatom
                + grazingFlux_Det    * recipDet         * grazEff & ! --> Grazing on detritus
     	        - hetLossFlux        * recipQZoo                  & ! --> Mortality loss
     	        - lossC_z                               * HetC    & ! --> Excretion loss
      	        - hetRespFlux                                     & ! --> REspiration loss
                                                                 ) * dt_b + sms(k,ihetc)
        else
            sms(k,ihetc)      = (                                 &
     	        + grazingFlux_phy    * recipQuota       * grazEff &
     	        + grazingFlux_Dia    * recipQuota_Dia   * grazEff &
     	        - hetLossFlux        * recipQZoo                  &
     	        - lossC_z                               * HetC    &
     	        - hetRespFlux                                     & 
                                                                 ) * dt_b + sms(k,ihetc)
        endif

!< *** DOM ***
!< ***********

!____________________________________________________________
!< DON (Extracellular organic N)

    sms(k,idon)      = (                      &
        + lossN * limitFacN         * phyN    &
        + lossN_d * limitFacN_Dia   * DiaN    &
        + reminN * arrFunc * O2Func * DetN    & 
        + lossN_z                   * HetN    &
        - rho_N * arrFunc * O2Func  * DON     & ! O2remin
                                             ) * dt_b + sms(k,idon)

!____________________________________________________________
!< EOC

    sms(k,idoc)       = (                     &
        + lossC * limitFacN         * phyC    &
        + lossC_d * limitFacN_dia   * DiaC    &
        + reminC * arrFunc * O2Func * DetC    &
        + lossC_z                   * HetC    &
        - rho_c1 * arrFunc * O2Func * EOC     & ! O2remin
                                             ) * dt_b + sms(k,idoc)	

!< *** Diatoms ***
!< ***************

!____________________________________________________________
!< Diatom N
 
!< lossN:  Diatom loss of organic N compounds [day^-1]
!< When N : C ratio becomes too high, excretion of DON is downregulated
!< by the limiter function limitFacN_dia
!< aggregationRate transfers N to the detritus pool

    sms(k,idian)      = (                        &
        + N_assim_dia                    * DiaC  & ! --> N assimilation
        - lossN_d * limitFacN_dia        * DiaN  & ! --> DON excretion
        - aggregationRate                * DiaN  & ! --> Aggregation loss
        - grazingFlux_Dia                        & ! --> Grazing loss
                                                ) * dt_b + sms(k,idian)

!____________________________________________________________
!< Diatom C

!< lossC_d: Diatom loss of carbon [day^-1]
!< When N : C ratio becomes too high, excretion of DOC is downregulated
!< by the limiter function limitFacN_dia
!< aggregationRate transfers C to the detritus pool

    sms(k,idiac)      = (                           &
        + Cphot_dia                         * DiaC  & ! -- Photosynthesis ---->/
        - lossC_d * limitFacN_dia           * DiaC  & ! -- Excretion of DOC --/ Net Photosynthesis
        - phyRespRate_dia                   * DiaC  & ! -- Respiration ----->/
        - aggregationRate                   * DiaC  &
        - grazingFlux_dia  * recipQuota_dia         &
     	                                           ) * dt_b + sms(k,idiac)

!____________________________________________________________
!< Diatom Chl

    sms(k,idchl)      = (                           &
        + chlSynth_dia                     * DiaC   & ! --> Chl a synthesis
        - KOchl_dia                        * DiaChl & ! --> Degradation loss
        - aggregationRate                  * DiaChl & ! --> Aggregation loss
        - grazingFlux_dia  * Chl2N_dia              & ! --> Grazing loss
                                                   ) * dt_b + sms(k,idchl)

!____________________________________________________________
!< Diatom Si

!< lossN_d: Diatom loss of organic nitrogen compunds [day^-1]
!< When N : C ratio becomes too high, excretion is downregulated
!< by the limiter function limitFacN_dia
!< aggregationRate transfers Si to the detritus pool

    sms(k,idiasi)        = (                        &
        + Si_assim                          * DiaC  & ! -- Diatom silicon assimilation
        - lossN_d * limitFacN_dia           * DiaSi & ! -- Excretion to detritus
        - aggregationRate                   * DiaSi & ! -- Aggregation loss
        - grazingFlux_dia  * qSiN                   & ! -- Grazing loss
                                                   ) * dt_b + sms(k,idiasi)

!< *** Silicate ***
!< ****************

!____________________________________________________________
!< Detritus Si
    sms(k,idetsi)     = (                          &
        + aggregationRate                  * DiaSi &
        + lossN_d         * limitFacN_dia  * DiaSi &
        + grazingFlux_dia * qSiN                   &
        - reminSiT                         * DetSi &
                                                  ) * dt_b + sms(k,idetsi)
!____________________________________________________________
!< DSi, Silicate 

!< DiaC:        Intracellular carbon concentration in diatoms [mmolC m-3]
!< DetSi:       Detritus silicon concentration [mmolSi m-3]
!< Si_assim:    Si assimilation rate for diatoms [mmolSi mmolC-1 day-1]
!< reminSiT:    Remineralization rate of silicon, temperature dependency [day-1]
!< dt_b:        REcoM time step [day]

!! Schourup 2013 Eq. A3

    sms(k,isi)        = (                            &
        - Si_assim                         * DiaC    &  ! --> Si assimilation of diatoms
        + reminSiT                         * DetSi   &  ! --> Remineralization of detritus, temperature dependent
                                                    ) * dt_b + sms(k,isi)
!< *** Iron ***
!< ************

!____________________________________________________________
!< Fe

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

    sms(k,ife) = ( Fe2N * (                         &
        - N_assim                         * PhyC    & ! --> N assimilation Nanophytoplankton, [mmol N/(mmol C * day)] C specific N utilization rate  
        - N_assim_dia                     * DiaC    & ! --> N assimilation Diatom
        + lossN         * limitFacN       * PhyN    & ! --> Excretion from small pythoplankton
        + lossN_d       * limitFacN_dia   * DiaN    & ! --> Excretion from diatom
        + reminN * arrFunc * O2Func * DetN          & ! --> Remineralization of detritus   ! NEW O2remin
        + lossN_z                   * HetN          & ! --> Excretion from zooplankton
                                              )     &
        - kScavFe          * DetC   * FreeFe        & 
                                                   ) * dt_b + sms(k,ife)

!< *** Calcification ***
!< *********************

!____________________________________________________________
!< Small phytoplankton calcite
    sms(k,iphycal)    = (                      &
        + calcification                        & ! --> Calcification  
        - lossC          * limitFacN * PhyCalc & ! --> Excretion loss  
        - phyRespRate                * PhyCalc & ! --> Respiration 
        - calc_loss_agg                        & ! --> Aggregation loss 
        - calc_loss_gra                        & ! --> Grazing loss
                                              ) * dt_b + sms(k,iphycal)

!____________________________________________________________
! Detritus calcite
    sms(k,idetcal)   = (                           &
        + lossC         * limitFacN      * PhyCalc &
        + phyRespRate                    * PhyCalc &
        + calc_loss_agg                            &
        + calc_loss_gra                            &
        - calc_loss_gra * calc_diss_guts           &
        - calc_diss                      * DetCalc &
                                                  ) * dt_b + sms(k,idetcal)

!____________________________________________________________
! Oxygen

    sms(k,ioxy)   = (                       &
      + Cphot                      * phyC   &
      - phyRespRate                * phyC   &
      + Cphot_dia                  * diaC   &
      - phyRespRate_dia            * diaC   &
      - rho_C1  * arrFunc * O2Func * EOC    & ! O2remin
      - hetRespFlux                         &
                                           ) * redO2C * dt_b + sms(k,ioxy)  
!   
!-------------------------------------------------------------------------------
! Diagnostics: Averaged rates
	
	recipbiostep    = 1.d0/real(biostep)

  end do ! Main vertikal loop ends

!-------------------------------------------------------------------------------
! Remineralization from the sediments into the bottom layer
!*** DIN ***
!< decayRateBenN: Remineralization rate for benthic N [day^-1]
!< LocBenthos(1): Vertically integrated N concentration in benthos (1 layer) [mmolN/m^2]
  decayBenthos(1) = decayRateBenN * LocBenthos(1)             
  LocBenthos(1)   = LocBenthos(1)   - decaybenthos(1) * dt_b ! remove from benthos (flux)

!*** DIC ***
!< decayRateBenC: Remineralization rate for benthic C [day^-1]
!< LocBenthos(2): Vertically integrated C concentration in benthos (1 layer) [mmolC/m^2]		
  decayBenthos(2) = decayRateBenC * LocBenthos(2)
  LocBenthos(2)   = LocBenthos(2)   - decaybenthos(2) * dt_b

!*** Si ***
!< decayRateBenSi: Remineralization rate for benthic Si [day^-1]
!< LocBenthos(3) : Vertically integrated N concentration in benthos (1 layer) [mmolSi/m^2]		
  decayBenthos(3) = decayRateBenSi * LocBenthos(3) ! [1/day] * [mmolSi/m2] -> [mmolSi/m2/day]
  LocBenthos(3)   = LocBenthos(3)   - decaybenthos(3) * dt_b 

!*** Calc: DIC, Alk ***  ! OG calc_diss_ben is taken from the deepest level 
  decayBenthos(4) = calc_diss_ben * LocBenthos(4)    ! NEW DISS changed calc_diss to calc_diss_ben to not make the dissolution omega dependent when using the switch OmegaC_diss
  LocBenthos(4)      = LocBenthos(4)   - decayBenthos(4) * dt_b

  end do ! Main time loop ends


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
function iron_chemistry_2ligands(fet,l1t,l2t,k1,k2)
      implicit none

      Real(kind=8) :: iron_chemistry_2ligands
      Real(kind=8) :: l1t,l2t,fet,k1,k2
      Real(kind=8) :: a3,a2,a1,a0,a,b,c,p,q,discr,rho,phi,amp,pi
      Real(kind=8) :: one3rd, one27th
      Real(kind=8) :: fe1,fe2,fe3

! coefficients of the 4th-order polynomial
      a3 = k1*k2
      a2 = ( k1*k2*(l1t + l2t - fet) + k1 + k2 )
      a1 = ( 1 - (k1 + k2)*fet + k1*l1t + k2*l2t )
      a0 = -fet

! coefficients of the normalized polynomial
      a = a2/a3
      b = a1/a3
      c = a0/a3

! some numbers that are used several times
      one3rd = 1.0/3.0
      one27th = 1.0/27.0

! now solve the polynomial stepwise
      p = b - a*a*one3rd
      q = c - a*b*one3rd + 2.0*a*a*a*one27th
      discr = q*q/4.0 + p*p*p*one27th

      rho = sqrt(-(p*p*p*one27th))
      phi = acos(-q/(2.0*rho))
      amp = 2.0*rho**one3rd
      pi = 3.1415926535897931

! the equation has three real roots
      fe1 = amp*cos(phi*one3rd) - a*one3rd
      fe2 = amp*cos((phi + 2.0*pi)*one3rd) - a*one3rd
      fe3 = amp*cos((phi + 4.0*pi)*one3rd) - a*one3rd

      iron_chemistry_2ligands = max(fe1,fe2,fe3)

end function iron_chemistry_2ligands
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

