subroutine REcoM_sms(n,Nn,state,thick,recipthick,SurfSW,sms,Temp,SinkVel,zF,PAR, mesh)

  Use REcoM_declarations
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


  Implicit none
  type(t_mesh), intent(in) , target                       :: mesh
  Integer,                           intent(in)           :: Nn                   ! Total number of nodes in the vertical
  Real(kind=8),dimension(mesh%nl-1,bgc_num),intent(inout) :: state                ! ChlA conc in phytoplankton [mg/m3]
											! should be in instead of inout

  Real(kind=8),dimension(mesh%nl-1)                       :: thick                ! [m] Vertical distance between two nodes = Thickness 
  Real(kind=8),dimension(mesh%nl-1)                       :: recipthick           ! [1/m] reciprocal of thick
  Real(kind=8),                      intent(in)           :: SurfSW               ! [W/m2] ShortWave radiation at surface

  Real(kind=8),dimension(mesh%nl-1,bgc_num),intent(inout) :: sms                  ! Source-Minus-Sinks term
  Real(kind=8),dimension(mesh%nl-1)        ,intent(in)    :: Temp                 ! [degrees C] Ocean temperature
  Real(kind=8),dimension(mesh%nl,3)        ,intent(in)    :: SinkVel

  Real(kind=8),dimension(mesh%nl)      ,intent(in)        :: zF                   ! [m] Depth of fluxes
  Real(kind=8),dimension(mesh%nl-1),intent(inout)         :: PAR

  Real(kind=8)                                            :: dt_d                   ! Size of time steps [day]
  Real(kind=8)                                            :: dt_b                   ! Size of time steps [day]
  Real(kind=8),dimension(mesh%nl-1)                       :: Sink
  Real(kind=8)                                            :: dt_sink              ! Size of local time step
  Real(kind=8)                                            :: Fc                   ! Flux of labile C into sediment, used for denitrification calculation [umolC/cm2/s]
  Real(kind=8)                                            :: recip_hetN_plus      ! MB's addition to heterotrophic respiration
  Real(kind=8)                                            :: recip_res_het        ! [day] Reciprocal of respiration by heterotrophs and mortality (loss to detritus)

  Integer                                                 :: k,step,ii, idiags,n
  Real(kind=8)                                            :: & 
    DIN,     & ! [mmol/m3] Dissolved Inorganic Nitrogen 	
    DIC,     & ! [mmol/m3] Dissolved Inorganic Carbon
    Alk,     & !? [mmol/m3] Total Alkalinity
    PhyN,    & ! [mmol/m3] Intracellular conc of Nitrogen in small phytoplankton
    PhyC,    & ! [mmol/m3] Intracellular conc of Carbon in small phytoplankton
    ChlA,    & ! [mg/m3]   Current intracellular ChlA conc
    DetN,    & ! [mmol/m3] Conc of N in Detritus
    DetC,    & ! [mmol/m3] Conc of C in Detritus
    HetN,    & ! [mmol/m3] Conc of N in heterotrophs
    HetC,    & ! [mmol/m3] Conc of C in heterotrophs
    DON,     & ! [mmol/m3] Dissolved organic N in the water
    EOC,     & !? [mmol/m3] Extracellular Organic C conc
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
!    endif
    FreeFe,  &
    O2
     ! mmol/m3
#include "../associate_mesh.h"

  sms = zero
  tiny_N   = tiny_chl/chl2N_max
  tiny_N_d = tiny_chl/chl2N_max_d
  tiny_C   = tiny_N  /NCmax    
  tiny_C_d = tiny_N_d/NCmax_d
  tiny_Si  = tiny_C_d/SiCmax

  recip_res_het = 1./res_het

!-------------------------------------------------------------------------------
! Size of REcoM time steps are calculated [day]
!-------------------------------------------------------------------------------

  rTref =  real(one)/recom_Tref
  
  dt_d	 =	dt/SecondsPerDay     ! Size of FESOM time step [day]
  dt_b	 =	dt_d/real(biostep)   ! Size of REcoM time step [day]

!-------------------------------------------------------------------------------
!Main time loop starts
  do step  = one,biostep

    kdzUpper	= 0.d0	            ! Upper light attenuation of top cell is set to zero

    if (any(abs(sms(:,:)) <= tiny)) sms(:,:) = zero

!    do k  = one,Nn
!      do ii = one,bgc_num
!        if (abs(sms(k,ii)) .le. tiny) sms(k,ii) = zero
!      end do
!    end do 

!-------------------------------------------------------------------------------
! Main vertical loop starts
  do k = one,Nn
!    do n=1, myDim_nod2D!+eDim_nod2D 
!       Nn=nlevels_nod2D(n)-1  !nzmax
!       do k=1, Nn

    DIN    = max(tiny,state(k,idin)  	 	+ sms(k,idin  )) ! Avoids division by zero
    DIC    = max(tiny,state(k,idic)   		+ sms(k,idic  )) ! and updates Conc between
    ALK    = max(tiny,state(k,ialk)   		+ sms(k,ialk  )) ! local steps in REcoM when
    PhyN   = max(tiny_N,state(k,iphyn)  	+ sms(k,iphyn )) ! biostep > 1
    PhyC   = max(tiny_C,state(k,iphyc) 		+ sms(k,iphyc ))
    ChlA   = max(tiny_chl,state(k,ipchl)  	+ sms(k,ipchl ))
    DetN   = max(tiny,state(k,idetn)  		+ sms(k,idetn ))
    DetC   = max(tiny,state(k,idetc)  		+ sms(k,idetc ))
    HetN   = max(tiny,state(k,ihetn)  		+ sms(k,ihetn ))
    HetC   = max(tiny,state(k,ihetc)  		+ sms(k,ihetc ))
if (REcoM_Second_Zoo) then 
    Zoo2N  = max(tiny,state(k,izoo2n)           + sms(k,izoo2n))
    Zoo2C  = max(tiny,state(k,izoo2c)           + sms(k,izoo2c))
endif
    DON    = max(tiny,state(k,idon)   		+ sms(k,idon  ))
    EOC    = max(tiny,state(k,idoc)   		+ sms(k,idoc  ))
    DiaN   = max(tiny_N,state(k,idian)  	+ sms(k,idian ))
    DiaC   = max(tiny_C,state(k,idiac)  	+ sms(k,idiac ))
    DiaChl = max(tiny_chl,state(k,idchl)  	+ sms(k,idchl ))
    DiaSi  = max(tiny_si,state(k,idiasi) 	+ sms(k,idiasi))
    DetSi  = max(tiny,state(k,idetsi) 		+ sms(k,idetsi))
    Si     = max(tiny,state(k,isi)    		+ sms(k,isi   ))
    Fe     = max(tiny,state(k,ife)    		+ sms(k,ife   ))
    O2     = max(tiny,state(k,ioxy)             + sms(k,ioxy   ))
    FreeFe = zero
!#ifdef REcoM_calcification		
    PhyCalc= max(tiny,state(k,iphycal)		+ sms(k,iphycal))
    DetCalc= max(tiny,state(k,idetcal)		+ sms(k,idetcal))

!YY: SinkVel(nl,3) calculated in recom_forcing (routine:depth_calculations in
!recom_extra), sinking velocity at the 1. level is 0 and we start from the
!2. level
    calc_diss      = calc_diss_rate * SinkVel(k+1,ivdet) /20.d0 ! Dissolution rate of CaCO3 scaled by the sinking velocity at the current depth
!#endif

    quota          	=  PhyN / PhyC
    recipquota     	=  real(one) / quota
    Chl2C          	=  ChlA  / PhyC
    Chl2N          	=  ChlA  / PhyN
    CHL2C_plast    	=  Chl2C * (quota/(quota - NCmin))
    
    quota_dia      	=  DiaN / DiaC
    recipQuota_dia 	=  real(one)/quota_dia
    Chl2C_dia      	=  DiaChl / DiaC
    Chl2N_dia      	=  DiaChl / DiaN
    CHL2C_plast_dia 	=  Chl2C_dia * (quota_dia/(quota_dia - NCmin_d))
    qSiC           	=  DiaSi / DiaC
    qSiN           	=  DiaSi / DiaN


    recipQZoo      	=  HetC / HetN
    recip_hetN_plus 	= 1. / (hetN + tiny_het) ! MB's addition for more stable zoo respiration
    if (REcoM_Second_Zoo) then  
      recipQZoo2     =  Zoo2C / Zoo2N
    endif
    if (ciso) then
!     additional variables are declared in module REcoM_ciso
      DIC_13              = max(tiny,state(k,idic_13)    + sms(k,idic_13  ))
      PhyC_13             = max(tiny_C,state(k,iphyc_13) + sms(k,iphyc_13 ))
      DetC_13             = max(tiny,state(k,idetc_13)   + sms(k,idetc_13 ))
      HetC_13             = max(tiny,state(k,ihetc_13)   + sms(k,ihetc_13 ))
      EOC_13              = max(tiny,state(k,idoc_13)    + sms(k,idoc_13  ))
      DiaC_13             = max(tiny_C,state(k,idiac_13) + sms(k,idiac_13 ))
      PhyCalc_13          = max(tiny,state(k,iphycal_13) + sms(k,iphycal_13))
      DetCalc_13          = max(tiny,state(k,idetcal_13) + sms(k,idetcal_13))

      calc_diss_13        = alpha_dcal_13 * calc_diss

      quota_13            = PhyN / PhyC_13
      recipQuota_13       = real(one) / quota_13

      quota_dia_13        = DiaN / DiaC_13
      recipQuota_dia_13   = real(one) / quota_dia_13

      recipQZoo_13        = HetC_13 / HetN 

      if (ciso_14) then
        DIC_14            = max(tiny,state(k,idic_14)    + sms(k,idic_14  ))
        if (ciso_organic_14) then
          PhyC_14           = max(tiny_C,state(k,iphyc_14) + sms(k,iphyc_14 ))
          DetC_14           = max(tiny,state(k,idetc_14)   + sms(k,idetc_14 ))
          HetC_14           = max(tiny,state(k,ihetc_14)   + sms(k,ihetc_14 ))
          EOC_14            = max(tiny,state(k,idoc_14)    + sms(k,idoc_14  ))
          DiaC_14           = max(tiny_C,state(k,idiac_14) + sms(k,idiac_14 ))
          PhyCalc_14        = max(tiny,state(k,iphycal_14) + sms(k,iphycal_14))
          DetCalc_14        = max(tiny,state(k,idetcal_14) + sms(k,idetcal_14))
          
          calc_diss_14      = alpha_dcal_14 * calc_diss
          
          quota_14          = PhyN / PhyC_14
          recipQuota_14     = real(one) / quota_14
          
          quota_dia_14      = DiaN / DiaC_14
          recipQuota_dia_14 = real(one) / quota_dia_14
          recipQZoo_14      = HetC_14 / HetN
        end if ! ciso_organic_14 
      end if   ! ciso_14
    end if     ! ciso

!-------------------------------------------------------------------------------
! Temperature dependence of rates
!------------------------------------------------------------------------------- 
!! Schourup 2013 Eq. A54
!! Temperature dependence of metabolic rate, fT, dimensionless
!! rTloc: inverse of local temperature in [1/Kelvin]
!! rTref=288.15 (15 degC): Reference temperature for Arrhenius function [1/Kelvin]
!!  See Figure A1 
    rTloc          	= real(one)/(Temp(k) + C2K)
    arrFunc        	= exp( -Ae * ( rTloc - rTref))
    if (REcoM_Second_Zoo) then 
      arrFuncZoo2   = exp(t1_zoo2/t2_zoo2 - t1_zoo2*rTloc)/(1 + exp(t3_zoo2/t4_zoo2 - t3_zoo2*rTloc))
    endif

!! Silicate temperature dependence
    reminSiT 		= min(1.32e16 * exp(-11200.d0 * rTloc),reminSi) !! arrFunc control

!-------------------------------------------------------------------------------
! Photosynthesis section, light parameters and rates
!-------------------------------------------------------------------------------

!   if (mype==0) then
!      write(*,*) '____________________________________________________________'
!      write(*,*) ' --> Photosynthesis section'
!   endif
   
   !_______________________________________________________________________
   !! Small phytoplankton
   !! qlimitFac, qlimitFacTmp: Factor that regulates photosynthesis
   !! NMinSlope: 50.d0
   !! NCmin: 0.04d0
   !! quota: PhyN/PhyC
   !! qlimitFac [0.0, 1.0] 
   !! if quota < NCmin qlimitFac=0
   !! if quota > ≈ 9 * NCmin qlimitFac=1
   !! P_cm: 3.0d0 [1/day], Rate of C-specific photosynthesis 
   qlimitFac     	= recom_limiter(NMinSlope,NCmin,quota)
   feLimitFac  		= Fe/(k_Fe + Fe)
   qlimitFac   		= min(qlimitFac,feLimitFac)

   pMax          	= P_cm * qlimitFac * arrFunc
    
   !_______________________________________________________________________
   ! Diatoms
    qlimitFac     	= recom_limiter(NMinSlope,NCmin_d,quota_dia)
    qlimitFacTmp  	= recom_limiter(SiMinSlope,SiCmin,qSiC)
    qlimitFac     	= min(qLimitFac,qlimitFacTmp)
      feLimitFac  	= Fe/(k_Fe_d + Fe)
      qlimitFac   	= min(qlimitFac,feLimitFac)

    pMax_dia      	= P_cm_d * qlimitFac * arrFunc

   !_______________________________________________________________________
   ! Light
   if (k==1) then 
      PARave = max(tiny,SurfSW)
      PAR(k) = PARave
      chl_upper = (ChlA + DiaChl)
   else    
      chl_lower = ChlA + DiaChl
      Chlave    = (chl_upper+chl_lower)*0.5

      kappar          =  k_w + a_chl * (Chlave)
      kappastar      =  kappar / cosAI(n)
      kdzLower       =  kdzUpper + kappastar * thick(k)
      Lowerlight     =  SurfSW * exp(-kdzLower)
      Lowerlight         =  max(tiny,Lowerlight)
      PARave         =  Lowerlight
      PAR(k)         =  PARave
      chl_upper      =  chl_lower
      kdzUpper       =  kdzLower
   end if

!-------------------------------------------------------------------------------
! Small phytoplankton photosynthesis rate

    if ( pMax .lt. tiny .OR. PARave /= PARave                  &
         .OR. CHL2C /= CHL2C) then
      Cphot       = zero
    else
      Cphot       = pMax*( real(one) &
                   - exp(-alfa * Chl2C * PARave / pMax))
    end if
    if ( Cphot .lt. tiny) Cphot = zero
    
 !------------------------------------------------------------------------------
 ! Diatom photosynthesis rate

    if ( pMax_dia .lt. tiny .OR. PARave /= PARave               &
         .OR. CHL2C_dia /= CHL2C_dia) then
      Cphot_dia   = zero
    else
      Cphot_dia   = pMax_dia * (real(one) &
           	- exp( -alfa_d * Chl2C_dia * PARave / pMax_dia))
    end if
    if (Cphot_dia .lt. tiny) Cphot_dia = zero

!-------------------------------------------------------------------------------- 
! chlorophyll degradation

    KOchl = deg_Chl
    KOchl_dia = deg_Chl_d
        
    if (use_photodamage) then
!    Phyto chla loss                                                                                                                                                      
! adding a minimum value for photodamage 
         if (pMax .lt. tiny .OR. PARave /= PARave                 &
             .OR. CHL2C_plast /= CHL2C_plast) then
           KOchl = deg_Chl*0.1d0
        else
           KOchl = deg_Chl*(1-exp(-alfa * CHL2C_plast * PARave / pMax))
           KOchl = max((deg_Chl*0.1d0), KOchl)
        end if

!    Diatoms chla loss                                                                                                                                                     
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
! Assimilation section
!-------------------------------------------------------------------------------
! Compute assimilation from Geider et al 1998

!   if (mype==0) then
!      write(*,*) '____________________________________________________________'
!      write(*,*) ' --> Assimilation section'
!   endif

    V_cm           = V_cm_fact
    limitFacN      = recom_limiter(NMaxSlope,quota,NCmax)
    N_assim        = V_cm * pMax * NCuptakeRatio &                ! [mmol N / (mmol C * day)]
                      * limitFacN * (DIN/(DIN + k_din))

    V_cm           = V_cm_fact_d
    limitFacN_dia  = recom_limiter(NMaxSlope,quota_dia,NCmax_d)
    N_assim_dia    = V_cm * pMax_dia * NCUptakeRatio_d &
                      * limitFacN_dia * DIN/(DIN + k_din_d)

    limitFacSi     = recom_limiter(SiMaxSlope,qSiC,SiCmax)  &
                      * limitFacN_dia
    Si_assim       = V_cm * P_cm_d * arrFunc * SiCUptakeRatio &
                      * limitFacSi * Si/(Si + k_si)

!-------------------------------------------------------------------------------
! Iron chemistry
 
      freeFe      = iron_chemistry(Fe,totalligand,ligandStabConst)

!-------------------------------------------------------------------------------
! Chlorophyll synthesis

    chlSynth       = zero
    if (PARave .ge. tiny .AND. PARave .eq. PARave) then
       chlSynth = N_assim * Chl2N_max                    &
          * min( real(one),Cphot/(alfa*Chl2C*PARave))
    end if
    ChlSynth_dia   = zero
    if (PARave .ge. tiny .AND. PARave .eq. PARave) then
       ChlSynth_dia = N_assim_dia                        &
          * Chl2N_max_d * min(real(one),                 &
            Cphot_dia /(alfa_d * Chl2C_dia * PARave))
    end if
!-------------------------------------------------------------------------------
! Phytoplankton respiration rate
    phyRespRate     = res_phy * limitFacN + biosynth * N_assim
    phyRespRate_dia = res_phy_d * limitFacN_dia +        &
            biosynth * N_assim_dia + biosynthSi * Si_assim
    if (ciso) then
!       we assume that
!       phyRespRate_13,14     = phyRespRate
!       phyRespRate_dia_13,14 = phyRespRate_dia
    end if
!-------------------------------------------------------------------------------
! Zooplankton grazing on small phytoplankton and diatoms
! At the moment there is no preference for one or the other food. Change this!

      if (REcoM_Grazing_Variable_Preference) then
        DiaNsq        = DiaN * DiaN
        varpzDia      = pzDia * DiaNsq /(sDiaNsq + DiaNsq)
        fDiaN         = varpzDia * DiaN
      else
        fDiaN         = pzDia * DiaN
      end if
      food            = PhyN + fDiaN
      foodsq          = food * food
      grazingFlux     = (Graz_max * foodsq)/(epsilonr + foodsq) * HetN * arrFunc
      grazingFlux_phy = grazingFlux * phyN / food
      grazingFlux_Dia = grazingFlux * fDiaN / food

     if (REcoM_Second_Zoo) then
!-------------------------------------------------------------------------------
! Second Zooplankton grazing on small phytoplankton, diatoms and heterotrophs
! At the moment there is no preference for one or the other food. Change this!
        
      if (REcoM_Grazing_Variable_Preference) then
        DiaNsq2        = DiaN * DiaN
        varpzDia2      = pzDia2 * DiaNsq2 /(sDiaNsq2 + DiaNsq2)
        fDiaN2         = varpzDia2 * DiaN
        PhyNsq2        = PhyN * PhyN
        varpzPhy2      = pzPhy2 * PhyNsq2 /(sPhyNsq2 + PhyNsq2)
        fPhyN2         = varpzPhy2 * PhyN
        HetNsq         = HetN * HetN
        varpzHet       = pzHet * HetNsq /(sHetNsq + HetNsq)
        fHetN          = varpzHet * HetN   
      else
        fDiaN2         = pzDia2 * DiaN
        fPhyN2         = pzPhy2 * PhyN
        fHetN          = pzHet * HetN
      end if
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

!-------------------------------------------------------------------------------
! Heterotrophic respiration is assumed to drive zooplankton back to Redfield C:N
! if their C:N becomes higher than Redfield
    if (HetRespFlux_plus) then
       HetRespFlux    = recip_res_het * arrFunc * (hetC    * recip_hetN_plus - redfield) * HetC
    else
! default computation scheme 
       HetRespFlux    = recip_res_het * arrFunc * (recipQZoo    - redfield) * HetC  
    endif
! Next part changes results, but is needed: Otherwise heterotrophs take up 
! inorganic C when their C:N becomes lower than Redfield

    HetRespFlux    = max(zero,HetRespFlux)
    if (ciso) then
!MB    set HetRespFlux_plus = .true. in namelist.recom
!      HetRespFlux_13   = max(zero, recip_res_het * arrFunc * (hetC_13 * recip_hetN_plus - redfield) * HetC_13)
!      Numerically safer parametrization avoiding instable results which may result from different cutoff values -- CHECK
       HetRespFlux_13     = HetRespFlux * HetC_13 / HetC 
!!     HetRespFlux_13     = HetRespFlux * (HetC_13 / HetC) **2
       if (ciso_14 .and. ciso_organic_14) then
!        HetRespFlux_14 = max(zero, recip_res_het * arrFunc * (hetC_14 * recip_hetN_plus - redfield) * HetC_14)
         HetRespFlux_14   = HetRespFlux * HetC_14 / HetC
!!       HetRespFlux_14   = HetRespFlux * (HetC_14 / HetC) **2
       end if
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
 end if
!-------------------------------------------------------------------------------
! Phytoplankton and detritus aggregation
    aggregationrate = agg_PD * DetN + agg_PP * PhyN &
                    + agg_PP * DiaN
!-------------------------------------------------------------------------------
! Terms required for the formation and dissolution of CaCO3
!#ifdef REcoM_calcification
    calcification = calc_prod_ratio * Cphot * PhyC   ! Z in equations
    calc_loss_agg = aggregationRate * PhyCalc

  if(REcoM_Second_Zoo)  then
     calc_loss_gra  = (grazingFlux_phy + grazingFlux_phy2)  &
                     * recipQuota/(PhyC + tiny)    * PhyCalc
  else
    calc_loss_gra = grazingFlux_phy           &
                    * recipQuota/(PhyC + tiny)    * PhyCalc
  endif  
    if (ciso) then
      calcification_13   = calc_prod_ratio * Cphot * PhyC_13 * alpha_calc_13
      calcification_13   = calcification   * alpha_calc_13  !! CHECK
      calc_loss_agg_13   = aggregationRate * PhyCalc_13
      calc_loss_gra_13   = grazingFlux_phy * recipQuota_13/(PhyC_13 + tiny) * PhyCalc_13
      if (ciso_14 .and. ciso_organic_14) then
        calcification_14 = calc_prod_ratio * Cphot * PhyC_14 * alpha_calc_14
        calc_loss_agg_14 = aggregationRate * PhyCalc_14
        calc_loss_gra_14 = grazingFlux_phy * recipQuota_14/(PhyC_14 + tiny) * PhyCalc_14
      end if
    end if

!#endif
		
!-------------------------------------------------------------------------------
! Sources minus sinks are calculated
!-------------------------------------------------------------------------------

!   if (mype==0) then
!      write(*,*) '____________________________________________________________'
!      write(*,*) ' --> Sources minus sinks section'
!   endif

!____________________________________________________________
! ---- DIN 

!! DON:            Extracellular dissolved organic nitrogen [mmolN m-3]
!! rho_N*arrFunc : Remineralization rate, temperature dependency is calculated with arrFunc [day-1]
!! DiaC:           Intracellular carbon concentration in diatoms [mmolC m-3]
!! PhyC:           Intracellular carbon concentration in nanophytoplankton [mmolC m-3]
!! N_assim:        N assimilation rate for nanophytoplankton [mmolN mmolC-1 day-1]
!! N_assim_Dia:    N assimilation rate for diatoms [mmolN mmolC-1 day-1]
!! dt_b:           Size of REcoM time step [day]

!! Schourup 2013 Eq. A2

    sms(k,idin)      = (                       &
      - N_assim                      * PhyC    &  ! --> N assimilation Nanophytoplankton, [mmol N/(mmol C * day)] C specific N utilization rate
      - N_assim_Dia                  * DiaC    &  ! --> N assimilation Diatoms
      + rho_N * arrFunc              * DON     &  ! --> DON remineralization, temperature dependent             
      + LocRiverDIN                            &
                                             ) * dt_b + sms(k,idin)  

!-------------------------------------------------------------------------------
! DIC
  if(REcoM_Second_Zoo)  then
    sms(k,idic)      = (                       &
     - Cphot                         * PhyC    &
     + phyRespRate                   * PhyC    &
     - Cphot_Dia                     * DiaC    &
     + phyRespRate_Dia               * DiaC    &
     + rho_C1 * arrFunc              * EOC     &
     + HetRespFlux                             &                      
     + Zoo2RespFlux                           &                    
     + calc_diss                     * DetCalc &
     + calc_loss_gra * calc_diss_guts          &
     - calcification                           &
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

!-------------------------------------------------------------------------------
! Alkalinity (Assumes that N:P follows a constant Redfield ratio
! N_assimC: 1.0625 = 1/16 + 1
    sms(k,ialk)      = (                       &
      + 1.0625 * N_assim             * PhyC    &
      + 1.0625 * N_assim_Dia         * DiaC    &
      - 1.0625 * rho_N * arrFunc     * DON     &
!#ifdef REcoM_calcification
      + 2.d0 * calc_diss             * DetCalc &
      + 2.d0 * calc_loss_gra * calc_diss_guts  &
      - 2.d0 * calcification                   &
!#endif 
                                             ) * dt_b + sms(k,ialk) 
!-------------------------------------------------------------------------------
! Phytoplankton N
   if (REcoM_Second_Zoo) then
     sms(k,iphyn)      = (                      &
      + N_assim                      * PhyC    &
      - lossN * limitFacN            * PhyN    &
      - aggregationRate              * phyN    &
      - grazingFlux_phy                        &
      - grazingFlux_phy2                       &     
  
                                             ) * dt_b + sms(k,iphyn)
   else

    sms(k,iphyn)      = (                      &
      + N_assim                      * PhyC    &
      - lossN * limitFacN            * PhyN    &
      - aggregationRate              * phyN    &
      - grazingFlux_phy                        &
                                             ) * dt_b + sms(k,iphyn)
   endif
!-------------------------------------------------------------------------------
! Phytoplankton C
   if (REcoM_Second_Zoo) then

    sms(k,iphyc)      = (                      &
      + Cphot                        * PhyC    &
      - lossC * limitFacN            * PhyC    &
      - phyRespRate                  * PhyC    &
      - aggregationRate              * PhyC    &
      - grazingFlux_phy * recipQuota           &
      - grazingFlux_phy2 * recipQuota          &
                                             ) * dt_b + sms(k,iphyc)
   else
    sms(k,iphyc)      = (                      &
      + Cphot                        * PhyC    &
      - lossC * limitFacN            * PhyC    &
      - phyRespRate                  * PhyC    &
      - aggregationRate              * PhyC    &
      - grazingFlux_phy * recipQuota           &
                                             ) * dt_b + sms(k,iphyc)
   endif
!-------------------------------------------------------------------------------
! Phytoplankton ChlA
   if (REcoM_Second_Zoo) then
   
    sms(k,ipchl)       = (                       &
     	+ chlSynth                     * phyC    &
     	- KOchl                        * ChlA     &
     	- aggregationRate              * ChlA     &
     	- grazingFlux_phy * Chl2N                &
        - grazingFlux_phy2 * Chl2N               & 
                                               ) * dt_b + sms(k,ipchl)
   else
    sms(k,ipchl)       = (                       &
     	+ chlSynth                     * phyC    &
     	- KOchl                        * ChlA    &
     	- aggregationRate              * ChlA    &
     	- grazingFlux_phy * Chl2N                &
                                               ) * dt_b + sms(k,ipchl)
   endif
!-------------------------------------------------------------------------------
! Detritus N
   if (REcoM_Second_Zoo) then
    sms(k,idetn)       = (                       &
        + grazingFlux_phy                        &
        + grazingFlux_dia                        &
        - grazingFlux * grazEff                  &
        + grazingFlux_phy2                       &
        + grazingFlux_dia2                       &
        + grazingFlux_het2                       &
        - grazingFlux2 * grazEff2                &
        + aggregationRate              * PhyN    &
        + aggregationRate              * DiaN    &
        + hetLossFlux                            &
        + Zoo2LossFlux                           &
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
   endif   
!-------------------------------------------------------------------------------
! Detritus C
   if (REcoM_Second_Zoo) then
    sms(k,idetc)       = (                             &
        + grazingFlux_phy * recipQuota                 &
        - grazingFlux_phy * recipQuota * grazEff       &
        + grazingFlux_Dia * recipQuota_Dia             &
        - grazingFlux_Dia * recipQuota_Dia * grazEff   &
        + grazingFlux_phy2 * recipQuota                &
        - grazingFlux_phy2 * recipQuota * grazEff2     &
        + grazingFlux_Dia2 * recipQuota_Dia            &
        - grazingFlux_Dia2 * recipQuota_Dia * grazEff2 &
        + grazingFlux_het2 * recipQZoo                 &
        - grazingFlux_het2 * recipQZoo * grazEff2      &
        + aggregationRate              * phyC          &
        + aggregationRate              * DiaC          &
        + hetLossFlux * recipQZoo                      &
         + Zoo2LossFlux * recipQZoo2                   &
        - reminC * arrFunc             * DetC          &
                                             )   * dt_b + sms(k,idetc)
   else

    sms(k,idetc)       = (                           &
     	+ grazingFlux_phy * recipQuota               &
     	- grazingFlux_phy * recipQuota * grazEff     &
     	+ grazingFlux_Dia * recipQuota_Dia           &
     	- grazingFlux_Dia * recipQuota_Dia * grazEff &
     	+ aggregationRate              * phyC        &
     	+ aggregationRate              * DiaC        &
     	+ hetLossFlux * recipQZoo                    &
     	- reminC * arrFunc             * DetC        &
                                             )   * dt_b + sms(k,idetc)
   endif   
!-------------------------------------------------------------------------------
! Heterotrophic N
   if (REcoM_Second_Zoo) then
    sms(k,ihetn)       = (                       &
    	+ grazingFlux * grazEff                  &
        - grazingFlux_het2                       &
     	- hetLossFlux                            &
     	- lossN_z                      * HetN    &
                                               ) * dt_b + sms(k,ihetn)
   else

    sms(k,ihetn)       = (                       &
     	+ grazingFlux * grazEff                  &
     	- hetLossFlux                            &
     	- lossN_z                      * HetN    &
                                               ) * dt_b + sms(k,ihetn)
   endif   

!-------------------------------------------------------------------------------
! Heterotrophic C
   if (REcoM_Second_Zoo) then
    sms(k,ihetc)      = (                            &
     	+ grazingFlux_phy * recipQuota * grazEff     &
     	+ grazingFlux_Dia * recipQuota_Dia * grazEff &
     	- hetLossFlux * recipQZoo                    &
        - grazingFlux_het2 * recipQZoo               &
     	- lossC_z                      * HetC        &
      	- hetRespFlux                                & 
                                                ) * dt_b + sms(k,ihetc)
   else
    sms(k,ihetc)      = (                            &
     	+ grazingFlux_phy * recipQuota * grazEff     &
     	+ grazingFlux_Dia * recipQuota_Dia * grazEff &
     	- hetLossFlux * recipQZoo                    &
     	- lossC_z                      * HetC        &
     	- hetRespFlux                                & 
                                               ) * dt_b + sms(k,ihetc)
   endif
!-------------------------------------------------------------------------------

   if (REcoM_Second_Zoo) then
 ! Second Zooplankton N                                                                                              
     sms(k,izoo2n)       = (                        &
         + grazingFlux2 * grazEff2                  &
         - Zoo2LossFlux                             &
         - lossN_z2                      * Zoo2N    &
                                               ) * dt_b + sms(k,izoo2n)
 !-------------------------------------------------------------------------------                       
  ! Second Zooplankton C                                                                                 
                                       
     sms(k,izoo2c)      = (                             &
         + grazingFlux_phy2 * recipQuota * grazEff2     &
         + grazingFlux_Dia2 * recipQuota_Dia * grazEff2 &
         + grazingFlux_het2 * recipQZoo * grazEff2      &
         - zoo2LossFlux * recipQZoo2                    &
         - lossC_z2                      * Zoo2C        &
         - Zoo2RespFlux                                 &
                                                ) * dt_b + sms(k,izoo2c)  
   endif                                                                               
!-------------------------------------------------------------------------------
! DON (Extracellular organic N)
   if (REcoM_Second_Zoo) then
    sms(k,idon)      = (                        &
      + lossN * limitFacN              * phyN   &
      + lossN_d * limitFacN_Dia        * DiaN   &
      + reminN * arrFunc               * DetN   &
      + lossN_z                        * HetN   &
      + lossN_z2                       * Zoo2N  &
      - rho_N * arrFunc                * DON    &
      + LocRiverDON                             &
                                             ) * dt_b + sms(k,idon)
   else

    sms(k,idon)      = (                        &
      + lossN * limitFacN              * phyN   &
      + lossN_d * limitFacN_Dia        * DiaN   &
      + reminN * arrFunc               * DetN   &
      + lossN_z                        * HetN   &
      - rho_N * arrFunc                * DON    &
      + LocRiverDON                             &
                                              ) * dt_b + sms(k,idon)
   endif
!-------------------------------------------------------------------------------
! EOC
   if (REcoM_Second_Zoo) then
    sms(k,idoc)       = (                       &
      + lossC * limitFacN              * phyC   &
      + lossC_d * limitFacN_dia        * DiaC   &
      + reminC * arrFunc               * DetC   &
      + lossC_z                        * HetC   &
      + lossC_z2                       * Zoo2C  &
      - rho_c1 * arrFunc               * EOC    &
      + LocRiverDOC                             &
                                              ) * dt_b + sms(k,idoc)	
   else 
    sms(k,idoc)       = (                       &
      + lossC * limitFacN              * phyC   &
      + lossC_d * limitFacN_dia        * DiaC   &
      + reminC * arrFunc               * DetC   &
      + lossC_z                        * HetC   &
      - rho_c1 * arrFunc               * EOC    &
      + LocRiverDOC                             &
                                              ) * dt_b + sms(k,idoc)
   endif 		
!-------------------------------------------------------------------------------
! Diatom N
   if (REcoM_Second_Zoo) then
    sms(k,idian)      = (                      &
      + N_assim_dia                    * DiaC  &
      - lossN_d * limitFacN_dia        * DiaN  &
      - aggregationRate                * DiaN  &
      - grazingFlux_Dia                        &
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
!-------------------------------------------------------------------------------
! Diatom C
   if (REcoM_Second_Zoo) then
    sms(k,idiac)      = (                      &
      + Cphot_dia                      * DiaC  &
      - lossC_d * limitFacN_dia        * DiaC  &
      - phyRespRate_dia                * DiaC  &
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
!-------------------------------------------------------------------------------
! Diatom Chl
   if (REcoM_Second_Zoo) then
    sms(k,idchl)      = (                       &
      + chlSynth_dia                   * DiaC   &
      - KOchl_dia                      * DiaChl &
      - aggregationRate                * DiaChl &
      - grazingFlux_dia * Chl2N_dia             &
      - grazingFlux_dia2 * Chl2N_dia            &                 
                                             ) * dt_b + sms(k,idchl)
   else
    sms(k,idchl)      = (                       &
      + chlSynth_dia                   * DiaC   &
      - KOchl_dia                      * DiaChl &
      - aggregationRate                * DiaChl &
      - grazingFlux_dia * Chl2N_dia             &
                                             ) * dt_b + sms(k,idchl)
   endif 
!-------------------------------------------------------------------------------
! Diatom Si
   if (REcoM_Second_Zoo) then
    sms(k,idiasi)        = (                    &
      + Si_assim                        * DiaC  &
      - lossN_d * limitFacN_dia         * DiaSi &
      - aggregationRate                 * DiaSi &
      - grazingFlux_dia * qSiN                  &
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
      + grazingFlux_dia2 * qSiN                 &
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
! ---- Silicate 

!    if (mype==0) print*,' sms(k,idin) ', sms(k,idin)

!! DiaC:        Intracellular carbon concentration in diatoms [mmolC m-3]
!! DetSi:       Detritus silicon concentration [mmolSi m-3]
!! Si_assim:    Si assimilation rate for diatoms [mmolSi mmolC-1 day-1]
!! reminSiT:    Remineralization rate of silicon, temperature dependency [day-1]
!! dt_b:        Size of REcoM time step [day]

!! Schourup 2013 Eq. A3

    sms(k,isi)        = (                       &
      - Si_assim                        * DiaC  &  ! --> Si assimilation of diatoms
      + reminSiT                        * DetSi &  ! --> Remineralization of detritus, temperature dependent
      + LocRiverDSi                             &
                                             ) * dt_b + sms(k,isi)
!-------------------------------------------------------------------------------
! Fe
   if (REcoM_Second_Zoo) then
    if (use_Fe2N) then
         sms(k,ife) = ( Fe2N * (                  &
            - N_assim                 * PhyC      &
            - N_assim_dia             * DiaC      &
            + lossN*limitFacN         * PhyN      &
            + lossN_d*limitFacN_dia   * DiaN      &
            + reminN * arrFunc        * DetN      &
            + lossN_z                 * HetN      &
            + lossN_z2                * Zoo2N     &         
                                              )   &
            - kScavFe                 * DetC * FreeFe &
                          ) * dt_b           + sms(k,ife)
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
          +  hetRespFlux                       &
          +  lossC_z2               * Zoo2C    &
          +  zoo2RespFlux                      &
                                          )    &
          -  kScavFe                * DetC * FreeFe   &  
                                            ) * dt_b + sms(k,ife)
    end if
   else

    if (use_Fe2N) then
         sms(k,ife) = ( Fe2N * (                  &
            - N_assim                 * PhyC      &
            - N_assim_dia             * DiaC      &
            + lossN*limitFacN         * PhyN      &
            + lossN_d*limitFacN_dia   * DiaN      &
            + reminN * arrFunc        * DetN      &
            + lossN_z                 * HetN )    &
            - kScavFe                 * DetC * FreeFe &
                          ) * dt_b           + sms(k,ife)
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
!-------------------------------------------------------------------------------
! Small phytoplankton calcite
!#ifdef REcoM_calcification

    sms(k,iphycal)    = (           &
      + calcification               &
      - lossC * limitFacN * phyCalc &
      - phyRespRate       * phyCalc &
      - calc_loss_agg               &
      - calc_loss_gra               &
                                            ) * dt_b + sms(k,iphycal)
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
    sms(k,ioxy)   = (               &
      + Cphot              * phyC  &
      - phyRespRate         * phyC  &
      + Cphot_dia          * diaC  &
      - phyRespRate_dia     * diaC  &
      - rho_C1  * arrFunc   * EOC   &
      - hetRespFlux                 &
                                        )*redO2C * dt_b + sms(k,ioxy)  
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
! Phytoplankton C_13
      sms(k,iphyc_13)      =  (                           &
                + Cphot                        * PhyC_13        &
                - lossC * limitFacN            * PhyC_13        &
                - phyRespRate                  * PhyC_13        &
                - aggregationRate              * PhyC_13        &
                - grazingFlux_phy * recipQuota_13               &
                                    ) * dt_b + sms(k,iphyc_13)
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
! Heterotrophic C_13
      sms(k,ihetc_13)      = (                            &
            + grazingFlux_phy * recipQuota_13 * grazEff     &
            + grazingFlux_Dia * recipQuota_dia_13 * grazEff &
            - hetLossFlux * recipQZoo_13                    &
            - lossC_z                      * HetC_13        &
            - hetRespFlux_13                                &
                                                 ) * dt_b + sms(k,ihetc_13)
!-------------------------------------------------------------------------------
! EOC_13
      sms(k,idoc_13)       = (                            &
                + lossC * limitFacN              * phyC_13      &
                + lossC_d * limitFacN_dia        * DiaC_13      &
                + reminC * arrFunc               * DetC_13      &
                + lossC_z                        * HetC_13      &
                - rho_c1 * arrFunc               * EOC_13       &
                + LocRiverDOC * r_iorg_13                       &
                                                    ) * dt_b + sms(k,idoc_13)
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
! Small phytoplankton calcite_13
      sms(k,iphycal_13)    = (                            &
            + calcification_13                              &
            - lossC * limitFacN * phyCalc_13                &
            - phyRespRate       * phyCalc_13                &
            - calc_loss_agg_13                              &
            - calc_loss_gra_13                              &
                                              ) * dt_b + sms(k,iphycal_13)
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
      if (ciso_14) then
!-------------------------------------------------------------------------------
        if (ciso_organic_14) then
! DIC_14
          sms(k,idic_14) =        (                         &
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
! Phytoplankton C_14
          sms(k,iphyc_14)      =  (                           &
            + Cphot                        * PhyC_14        &
            - lossC * limitFacN            * PhyC_14        &
            - phyRespRate                  * PhyC_14        &
            - aggregationRate              * PhyC_14        &
            - grazingFlux_phy * recipQuota_14               &
                                    ) * dt_b + sms(k,iphyc_14)
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
! Heterotrophic C_14
          sms(k,ihetc_14)      = (                            &
            + grazingFlux_phy * recipQuota_14 * grazEff     &
            + grazingFlux_Dia * recipQuota_dia_14 * grazEff &
            - hetLossFlux * recipQZoo_14                    &
            - lossC_z                      * HetC_14        &
            - hetRespFlux_14                                &
                                                 ) * dt_b + sms(k,ihetc_14)
!-------------------------------------------------------------------------------
! EOC_14
          sms(k,idoc_14)       = (                            &
            + lossC * limitFacN              * phyC_14      &
            + lossC_d * limitFacN_dia        * DiaC_14      &
            + reminC * arrFunc               * DetC_14      &
            + lossC_z                        * HetC_14      &
            - rho_c1 * arrFunc               * EOC_14       &
            + LocRiverDOC * r_iorg_14                       &
                                                    ) * dt_b + sms(k,idoc_14)
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
! Small phytoplankton calcite_14
          sms(k,iphycal_14)    = (                            &
            + calcification_14                              &
            - lossC * limitFacN * phyCalc_14                &
            - phyRespRate       * phyCalc_14                &
            - calc_loss_agg_14                              &
            - calc_loss_gra_14                              &
                                              ) * dt_b + sms(k,iphycal_14)
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
        else
!         "Abiotic" DIC_14, identical to DIC except for radioactive decay (-> recom_forcing)
          sms(k,idic_14) = sms(k,idic)
        end if ! ciso_organic_14
      end if   ! ciso_14
    end if     ! ciso

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

  if (use_MEDUSA .and. (sedflx_num .ne. 0)) then
!!---------------------------------------------------------------------
!! CV: calculation of the effect of sediment fluxes has been shifted to
!! S/R oce_ale_tracer.F90     
!!---------------------------------------------------------------------
!!$    sms(Nn,idic)    = sms(Nn,idic) + LocSed(1) * dt_b  &
!!$                      * recipthick(Nn) !recipdzF(Nn)
!!$    sms(Nn,ialk)    = sms(Nn,ialk) + LocSed(2) * dt_b  &
!!$                      * recipthick(Nn) !recipdzF(Nn)
!!$    sms(Nn,idin)    = sms(Nn,idin) + LocSed(3) * dt_b  &
!!$                      * recipthick(Nn) !recipdzF(Nn)
!!$    sms(Nn,isi)     = sms(Nn,isi) + LocSed(4) * dt_b  &
!!$                      * recipthick(Nn) !recipdzF(Nn)
!!$    sms(Nn,ioxy)    = sms(Nn,ioxy) + LocSed(5) * dt_b  &
!!$                      * recipthick(Nn) !recipdzF(Nn)
!!$   if(ciso) then
!!$    sms(Nn,idic_13) = sms(Nn,idic_13) + LocSed(6) * dt_b  &
!!$                      * recipthick(Nn) !recipdzF(Nn)
!!$    if(ciso_14) then 
!!$    sms(Nn,idic_14) = sms(Nn,idic_14) + LocSed(7) * dt_b  &
!!$                      * recipthick(Nn) !recipdzF(Nn)
!!$    endif
!!$   endif ! ciso
!!$
!!$!   Ying: Fe is not yet calculated in medusa and thus related to DIC or DIN
!!$   if(use_Fe2N) then
!!$    sms(Nn,ife)     = sms(Nn,ife) + LocSed(3) * dt_b  &
!!$                      * recipthick(Nn) * Fe2N_benthos
!!$   else
!!$    sms(Nn,ife)     = sms(Nn,ife) + LocSed(1) * dt_b  &
!!$                      * recipthick(Nn) * Fe2C_benthos
!!$   endif ! use_Fe2N
!!---------------------------------------------------------------------
  
  else ! not use_MEDUSA or sedflx_num = 0

!*** DIN ***		
    decayBenthos(1) = decayRateBenN * LocBenthos(1)                                       ! flux
    LocBenthos(1)      = LocBenthos(1)   - decaybenthos(1) * dt_b ! / depth of benthos    ! remove from benthos (flux)

!*** DIC ***
    decayBenthos(2) = decayRateBenC * LocBenthos(2)
    LocBenthos(2)      = LocBenthos(2)   - decaybenthos(2) * dt_b ! / depth of benthos

!*** Si ***
!if (mype==0) write (*,*) "sms1, n= " , sms(Nn,isi)

    decayBenthos(3) = decayRateBenSi * LocBenthos(3)                                      ! [Si/day] * [mmolSi/m2] -> [mmolSi/m2/day]
    LocBenthos(3)      = LocBenthos(3)   - decaybenthos(3) * dt_b ! / depth of benthos    ! [mmolSi/m2]

!*** Calc: DIC, Alk ***
    decayBenthos(4) = calc_diss * LocBenthos(4)
    LocBenthos(4)      = LocBenthos(4)   - decayBenthos(4) * dt_b ! / depth of benthos

    if (ciso) then
!*** DIC_13 ***  We ignore isotopic fractionation during remineralization.
        decayBenthos(5) = alpha_dcal_13   * decayRateBenC   * LocBenthos(5)
        LocBenthos(5)   = LocBenthos(5)   - decayBenthos(5) * dt_b
!*** Calc: DIC_13 ***
        decayBenthos(6) = calc_diss_13    * LocBenthos(6)
        LocBenthos(6)   = LocBenthos(6)   - decayBenthos(6) * dt_b ! / depth of benthos
      if (ciso_14) then
        if (ciso_organic_14) then
!*** DIC_14 ***  We ignore isotopic fractionation during remineralization.
          decayBenthos(7) = alpha_dcal_14   * decayRateBenC   * LocBenthos(7)
          LocBenthos(7)   = LocBenthos(7)   - decayBenthos(7) * dt_b
!*** Calc: DIC_14 ***
          decayBenthos(8) = calc_diss_14    * LocBenthos(8)
          LocBenthos(8)   = LocBenthos(8)   - decayBenthos(8) * dt_b ! / depth of benthos
        else
!         Do nothing here because sms(idic_14) is defined as sms(idic) further above
        end if ! ciso_organic_14
      end if   ! ciso_14
    end if ! ciso
  endif ! use_MEDUSA

  end do ! Main time loop ends

! CV this still needs to be shifted to oce_ale_tracer.F90
  
!-------------------------------------------------------------------------------
! save sinking flux of PON[1],POC[2],opal[3],Calc[4],POC13[5],POC14[6],Calc13[7],Calc14[8] 
!    if (use_MEDUSA) then
!         LocSinkFlx(1) = wFluxDet(1) + wFluxPhy(1) + wFluxDia(1)
!         LocSinkFlx(2) = wFluxDet(2) + wFluxPhy(2) + wFluxDia(2)
!         LocSinkFlx(3) = wFluxDet(3) + wFluxDia(3)
!         LocSinkFlx(4) = wFluxDet(4) + wFluxPhy(3)
!     if (ciso) then
!         LocSinkFlx(5) = wFluxDet(5) + wFluxPhy(5) + wFluxDia(5)
!         LocSinkFlx(7) = wFluxDet(6) + wFluxPhy(6)
!       if (ciso_14) then
!         LocSinkFlx(6) = wFluxDet(7) + wFluxPhy(7) + wFluxDia(7)
!         LocSinkFlx(8) = wFluxDet(8) + wFluxPhy(8)         
!       endif
!     endif
!    endif

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

