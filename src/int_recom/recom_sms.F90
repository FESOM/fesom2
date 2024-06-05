subroutine REcoM_sms(n,Nn,state,thick,recipthick,SurfSR,sms,Temp, Sali_depth &
        , CO2_watercolumn                                                    & 
        , pH_watercolumn                                                     &
        , pCO2_watercolumn                                                   &
        , HCO3_watercolumn                                                   &
        , CO3_watercolumn                                                    &
        , OmegaC_watercolumn                                                 &
        , kspc_watercolumn                                                   &
        , rhoSW_watercolumn                                                  &
        , Loc_slp, zF, PAR, Lond, Latd, mesh)

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
    use mvars                                                                       
    use mdepth2press                                                                ! Ballasting
    use gsw_mod_toolbox, only: gsw_sa_from_sp,gsw_ct_from_pt,gsw_rho                ! Ballasting

    implicit none
    type(t_mesh), intent(in) , target                       :: mesh
    integer, intent(in)                                     :: Nn                   !< Total number of nodes in the vertical
    real(kind=8),dimension(mesh%nl-1,bgc_num),intent(inout) :: state                !< ChlA conc in phytoplankton [mg/m3]
									  	    !! should be in instead of inout

    real(kind=8),dimension(mesh%nl-1)                       :: thick                !< [m] Vertical distance between two nodes = Thickness 
    real(kind=8),dimension(mesh%nl-1)                       :: recipthick           !< [1/m] reciprocal of thick
    real(kind=8),intent(in)                                 :: SurfSR               !< [W/m2] ShortWave radiation at surface

    real(kind=8),dimension(mesh%nl-1,bgc_num),intent(inout) :: sms                  !< Source-Minus-Sinks term
    real(kind=8),dimension(mesh%nl-1)        ,intent(in)    :: Temp                 !< [degrees C] Ocean temperature
    real(kind=8),dimension(mesh%nl-1)        ,intent(in)    :: Sali_depth           !< NEW MOCSY Salinity for the whole water column

    Real(kind=8),dimension(mesh%nl-1),intent(inout)         :: CO2_watercolumn      !< [mol/m3]
    Real(kind=8),dimension(mesh%nl-1),intent(inout)         :: pH_watercolumn       !< on total scale
    Real(kind=8),dimension(mesh%nl-1),intent(inout)         :: pCO2_watercolumn     !< [uatm]
    Real(kind=8),dimension(mesh%nl-1),intent(inout)         :: HCO3_watercolumn     !< [mol/m3]
    Real(kind=8),dimension(mesh%nl-1),intent(inout)         :: CO3_watercolumn      !< [mol/m3]
    Real(kind=8),dimension(mesh%nl-1),intent(inout)         :: OmegaC_watercolumn   !< calcite saturation state
    Real(kind=8),dimension(mesh%nl-1),intent(inout)         :: kspc_watercolumn     !< stoichiometric solubility product [mol^2/kg^2]
    Real(kind=8),dimension(mesh%nl-1),intent(inout)         :: rhoSW_watercolumn    !< in-situ density of seawater [kg/m3]

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

    Real(kind=8),                      intent(in)           :: Loc_slp              ! [Pa] sea-level pressure
    Real(kind=8)                                            :: Patm_depth(1)        
    Real(kind=8)                                            :: REcoM_T_depth(1)     ! temperature for the whole water column for mocsy minimum defined as -2
    Real(kind=8)                                            :: REcoM_S_depth(1)     
    Real(kind=8)                                            :: REcoM_DIC_depth(1)   
    Real(kind=8)                                            :: REcoM_Alk_depth(1)   
    Real(kind=8)                                            :: REcoM_Si_depth(1)    
    Real(kind=8)                                            :: REcoM_Phos_depth(1)  
    Real(kind=8),                      intent(in)           :: Latd(1)              ! latitude in degree
    Real(kind=8),                      intent(in)           :: Lond(1)              ! longitude in degree 
    Real(kind=8)                                            :: mocsy_step_per_day 
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
#if defined (__coccos) 
        CoccoN,  &
        CoccoC,  &
        CoccoChl,&
#endif
        Si,      &
        Fe,      &
        PhyCalc, &
        DetCalc, &
#if defined (__3Zoo2Det)
        Zoo2N,    &
        Zoo2C,    &
        DetZ2N,   &
        DetZ2C,   &
        DetZ2Si,  &
        DetZ2Calc,&
        MicZooN,  & ! 3Zoo
        MicZooC,  & ! 3Zoo
#endif
        FreeFe,  &
        O2
     
#include "../associate_mesh.h"

    sms = zero ! double precision

    tiny_N   = tiny_chl/chl2N_max      !< 0.00001/ 3.15d0   Chl2N_max [mg CHL/mmol N] Maximum CHL a : N ratio = 0.3 gCHL gN^-1
    tiny_N_d = tiny_chl/chl2N_max_d    !< 0.00001/ 4.2d0

    tiny_C   = tiny_N  /NCmax          !< NCmax   = 0.2d0   [mmol N/mmol C] Maximum cell quota of nitrogen (N:C)
    tiny_C_d = tiny_N_d/NCmax_d        !< NCmax_d = 0.2d0 

    tiny_Si  = tiny_C_d/SiCmax         !< SiCmax = 0.8d0

#if defined (__coccos) 
    tiny_N_c = tiny_chl/chl2N_max_c 
    tiny_C_c = tiny_N_c/NCmax_c
#endif

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
#if defined (__3Zoo2Det)
            Zoo2N     = max(tiny,state(k,izoo2n)     + sms(k,izoo2n))
            Zoo2C     = max(tiny,state(k,izoo2c)     + sms(k,izoo2c))
            DetZ2N    = max(tiny,state(k,idetz2n)    + sms(k,idetz2n))
            DetZ2C    = max(tiny,state(k,idetz2c)    + sms(k,idetz2c))
            DetZ2Si   = max(tiny,state(k,idetz2si)   + sms(k,idetz2si)) 
            DetZ2Calc = max(tiny,state(k,idetz2calc) + sms(k,idetz2calc))
            MicZooN   = max(tiny,state(k,imiczoon)   + sms(k,imiczoon))
            MicZooC   = max(tiny,state(k,imiczooc)   + sms(k,imiczooc))
#endif
            DON    = max(tiny,state(k,idon)          + sms(k,idon  ))
            EOC    = max(tiny,state(k,idoc)   	     + sms(k,idoc  ))
            DiaN   = max(tiny_N_d,state(k,idian)     + sms(k,idian ))
            DiaC   = max(tiny_C_d,state(k,idiac)     + sms(k,idiac ))
            DiaChl = max(tiny_chl,state(k,idchl)     + sms(k,idchl ))
            DiaSi  = max(tiny_si,state(k,idiasi)     + sms(k,idiasi)) 
            DetSi  = max(tiny,state(k,idetsi) 	     + sms(k,idetsi)) 
            Si     = max(tiny,state(k,isi)    	     + sms(k,isi   ))
#if defined (__coccos) 
            CoccoN   = max(tiny_N_c,state(k,icocn)   + sms(k,icocn ))
            CoccoC   = max(tiny_C_c,state(k,icocc)   + sms(k,icocc ))
            CoccoChl = max(tiny_chl,state(k,icchl)   + sms(k,icchl ))
#endif 
            Fe     = max(tiny,state(k,ife)    	     + sms(k,ife ))
            O2     = max(tiny,state(k,ioxy)          + sms(k,ioxy))
            FreeFe = zero

! For Mocsy
            REcoM_T_depth    = max(2.d0, Temp(k))        ! minimum set to 2 degC: K1/K2 Lueker valid between 2degC-35degC and 19-43psu
            REcoM_T_depth    = min(REcoM_T_depth, 40.d0) ! maximum set to 40 degC: K1/K2 Lueker valid between 2degC-35degC and 19-43psu
            REcoM_S_depth    = max(21.d0, Sali_depth(k)) ! minimum set to 21: K1/K2 Lueker valid between 2degC-35degC and 19-43psu, else causes trouble in regions with S between 19 and 21 and ice conc above 97%
            REcoM_S_depth    = min(REcoM_S_depth, 43.d0) ! maximum set to 43: K1/K2 Lueker valid between 2degC-35degC and 19-43psu, else causes trouble
            REcoM_DIC_depth  = max(tiny*1e-3,state(k,idic)*1e-3   + sms(k,idic  )*1e-3)     
            REcoM_Alk_depth  = max(tiny*1e-3,state(k,ialk)*1e-3   + sms(k,ialk  )*1e-3)      
            REcoM_Si_depth   = max(tiny*1e-3,state(k,isi)*1e-3    + sms(k,isi   )*1e-3)      
            REcoM_Phos_depth = max(tiny*1e-3,state(k,idin)*1e-3   + sms(k,idin  )*1e-3) /16 ! convert N to P with Redfield [mol/m3] 

            PhyCalc = max(tiny,state(k,iphycal)      + sms(k,iphycal))
            DetCalc = max(tiny,state(k,idetcal)      + sms(k,idetcal))

!-------------------------------------------------------------------------------
!< Quotas

            quota       =  PhyN / PhyC ! include variability of the N: C ratio, cellular chemical composition 
            recipquota  =  real(one) / quota
            Chl2C       =  PhyChl  / PhyC ! Chl a:phytoplankton carbon ratio, cellular chemical composition [gCHL gC^-1]
            Chl2N       =  PhyChl  / PhyN ! Chl a:phytoplankton nitrogen ratio, cellular chemical composition [gCHL gN^-1]
            CHL2C_plast =  Chl2C * (quota/(quota - NCmin))
    
            quota_dia       =  DiaN / DiaC
            recipQuota_dia  =  real(one)/quota_dia
            Chl2C_dia       =  DiaChl / DiaC
            Chl2N_dia       =  DiaChl / DiaN
            CHL2C_plast_dia =  Chl2C_dia * (quota_dia/(quota_dia - NCmin_d))
            qSiC            =  DiaSi / DiaC
            qSiN            =  DiaSi / DiaN
#if defined (__coccos) 
            quota_cocco       = CoccoN / CoccoC
            recipQuota_cocco  = real(one)/quota_cocco
            Chl2C_cocco       = CoccoChl / CoccoC
            Chl2N_cocco       = CoccoChl / CoccoN
            CHL2C_plast_cocco = Chl2C_cocco * (quota_cocco/(quota_cocco - NCmin_c))
#endif
            recipQZoo       = HetC / HetN
            recip_hetN_plus = 1.d0 / (HetN + tiny_het) ! MB's addition for more stable zoo respiration
            if (Grazing_detritus) recipDet  = DetC / DetN
#if defined (__3Zoo2Det)
            recipQZoo2 = Zoo2C / Zoo2N
            recipQZoo3 = MicZooC / MicZooN
            if (Grazing_detritus) recipDet2 = DetZ2C / DetZ2N
#endif

            if (ciso) then
!<       additional variables are declared in module REcoM_ciso
                DIC_13      = max(tiny,state(k,idic_13)    + sms(k,idic_13  ))
                PhyC_13     = max(tiny_C,state(k,iphyc_13) + sms(k,iphyc_13 ))
                DetC_13     = max(tiny,state(k,idetc_13)   + sms(k,idetc_13 ))
                HetC_13     = max(tiny,state(k,ihetc_13)   + sms(k,ihetc_13 ))
                EOC_13      = max(tiny,state(k,idoc_13)    + sms(k,idoc_13  ))
                DiaC_13     = max(tiny_C,state(k,idiac_13) + sms(k,idiac_13 ))
                PhyCalc_13  = max(tiny,state(k,iphycal_13) + sms(k,iphycal_13))
                DetCalc_13  = max(tiny,state(k,idetcal_13) + sms(k,idetcal_13))

                calc_diss_13      = alpha_dcal_13 * calc_diss

                quota_13          = PhyN / PhyC_13
                recipQuota_13     = real(one) / quota_13

                quota_dia_13      = DiaN / DiaC_13
                recipQuota_dia_13 = real(one) / quota_dia_13

                recipQZoo_13      = HetC_13 / HetN

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
#if defined (__coccos) 
            CoccoTFunc = max(0.1419d0 * Temp(k)**0.8151d0,tiny) ! Function from Fielding 2013; is based on observational GR, but range fits best to ours
#endif

#if defined (__3Zoo2Det)
            arrFuncZoo2 = exp(t1_zoo2/t2_zoo2 - t1_zoo2*rTloc)/(1 + exp(t3_zoo2/t4_zoo2 - t3_zoo2*rTloc)) ! 2Zoo
            q10_mes      = 1.0242**(Temp(k)) ! 3Zoo
            q10_mic      = 1.04**(Temp(k))   ! 3Zoo
            q10_mes_res  = 1.0887**(Temp(k)) ! 3Zoo
            q10_mic_res  = 1.0897**(Temp(k)) ! 3Zoo
#endif

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
#if defined (__coccos) 
                chl_upper = chl_upper + CoccoChl 
#endif
            else
                chl_lower  = PhyChl + DiaChl
#if defined (__coccos) 
                chl_lower  = chl_lower + CoccoChl
#endif            
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

!-------------------------------------------------------------------------------
! Depth component of Mocsy (see http://ocmip5.ipsl.jussieu.fr/mocsy/pyth.html)
!-------------------------------------------------------------------------------

! Calculate the carbonate system for the very first time step of the first year of the run
    !if (mocsy_restart==.false. .and. recom_istep==1) then    ! r_restart is defined in gen_modules_clock in fesom_cpl.
            dpos(1) = -zF(k)
            if (mstep==1) then
                call vars_sprac(ph_depth, pco2_depth, fco2_depth, co2_depth, hco3_depth, co3_depth, OmegaA_depth, OmegaC_depth, kspc_depth, BetaD_depth,  &
                       rhoSW_depth, p_depth, tempis_depth,                                                                                                &
                       REcoM_T_depth, REcoM_S_depth, REcoM_Alk_depth, REcoM_DIC_depth, REcoM_Si_depth, REcoM_Phos_depth, Patm_depth, dpos, Latd, Nmocsy,  &
                       optCON='mol/m3', optT='Tpot   ', optP='m ', optB='u74', optK1K2='l  ', optKf='dg', optGAS='Pinsitu', optS='Sprc')
                CO2_watercolumn(k)    = co2_depth(1)
                pH_watercolumn(k)     = ph_depth(1)
                pCO2_watercolumn(k)   = pco2_depth(1)
                HCO3_watercolumn(k)   = hco3_depth(1)
                CO3_watercolumn(k)    = co3_depth(1)
                OmegaC_watercolumn(k) = OmegaC_depth(1)
                kspc_watercolumn(k)   = kspc_depth(1)
                rhoSW_watercolumn(k)  = rhoSW_depth(1)
            endif 

!! Calculate carbonate system every 7 days for depths < 1%PAR, and every 30 days for the depths below.
            mocsy_step_per_day = 1/dt_b   ! NEW ms: time steps per day in recom -> is that correct? Not necessary to define in namelist? 
            logfile_outfreq_7  = mocsy_step_per_day*7
            logfile_outfreq_30 = mocsy_step_per_day*30

            if (PARave > 0.01*SurfSR .and. mod(mstep,logfile_outfreq_7)==0) then
                call vars_sprac(ph_depth, pco2_depth, fco2_depth, co2_depth, hco3_depth, co3_depth, OmegaA_depth, OmegaC_depth, kspc_depth, BetaD_depth,  & 
                       rhoSW_depth, p_depth, tempis_depth,                                                                                                &
                       REcoM_T_depth, REcoM_S_depth, REcoM_Alk_depth, REcoM_DIC_depth, REcoM_Si_depth, REcoM_Phos_depth, Patm_depth, dpos, Latd, Nmocsy,  &
                       optCON='mol/m3', optT='Tpot   ', optP='m ', optB='u74', optK1K2='l  ', optKf='dg', optGAS='Pinsitu', optS='Sprc')
                CO2_watercolumn(k)    = co2_depth(1)
                pH_watercolumn(k)     = ph_depth(1)
                pCO2_watercolumn(k)   = pco2_depth(1)
                HCO3_watercolumn(k)   = hco3_depth(1)
                CO3_watercolumn(k)    = co3_depth(1)
                OmegaC_watercolumn(k) = OmegaC_depth(1)
                kspc_watercolumn(k)   = kspc_depth(1)
                rhoSW_watercolumn(k)  = rhoSW_depth(1)

            elseif (PARave < 0.01*SurfSR .and. mod(mstep,logfile_outfreq_30)==0) then
                call vars_sprac(ph_depth, pco2_depth, fco2_depth, co2_depth, hco3_depth, co3_depth, OmegaA_depth, OmegaC_depth, kspc_depth, BetaD_depth,  &
                       rhoSW_depth, p_depth, tempis_depth,                                                                                                &
                       REcoM_T_depth, REcoM_S_depth, REcoM_Alk_depth, REcoM_DIC_depth, REcoM_Si_depth, REcoM_Phos_depth, Patm_depth, dpos, Latd, Nmocsy,  &
                       optCON='mol/m3', optT='Tpot   ', optP='m ', optB='u74', optK1K2='l  ', optKf='dg', optGAS='Pinsitu', optS='Sprc')
                CO2_watercolumn(k)    = co2_depth(1)
                pH_watercolumn(k)     = ph_depth(1)
                pCO2_watercolumn(k)   = pco2_depth(1)
                HCO3_watercolumn(k)   = hco3_depth(1)
                CO3_watercolumn(k)    = co3_depth(1)
                OmegaC_watercolumn(k) = OmegaC_depth(1)
                kspc_watercolumn(k)   = kspc_depth(1)
                rhoSW_watercolumn(k)  = rhoSW_depth(1)
            endif

!-------------------------------------------------------------------------------
! CO2 dependence of rates ! NEW CO2
!-------------------------------------------------------------------------------
! Convert pH to proton concentration
            h_depth(1) = 10.**(-ph_depth(1))
! Conversion factor Cunits not needed for [H], because in model and function derived from pH and therefore in [mol/L]

! Small phytoplankton
            PhyCO2 = a_co2_phy * HCO3_watercolumn(k) * Cunits / (b_co2_phy + HCO3_watercolumn(k) * Cunits) &
                - exp(-c_co2_phy * CO2_watercolumn(k) * Cunits) - d_co2_phy * 10.**(-pH_watercolumn(k))
            PhyCO2 = min(PhyCO2,3.d0)      ! April 2022: limitation to 3
            PhyCO2 = max(0.d0,PhyCO2)      ! July 2022: limitation to zero

! Diatoms
            DiaCO2 = a_co2_dia * HCO3_watercolumn(k) * Cunits / (b_co2_dia + HCO3_watercolumn(k) * Cunits) &
                - exp(-c_co2_dia * CO2_watercolumn(k) * Cunits) - d_co2_dia * 10.**(-pH_watercolumn(k))
            DiaCO2 = min(DiaCO2,3.d0)      ! April 2022: limitation to 3
            DiaCO2 = max(0.d0,DiaCO2)      ! July 2022: limitation to zero

#if defined (__coccos) 
! Coccolithophores
            CoccoCO2 = a_co2_cocco * HCO3_watercolumn(k) * Cunits / (b_co2_cocco + HCO3_watercolumn(k) * Cunits) & 
                - exp(-c_co2_cocco * CO2_watercolumn(k) * Cunits) - d_co2_cocco * 10.**(-pH_watercolumn(k))
            CoccoCO2 = min(CoccoCO2,3.d0)  ! April 2022: limitation to 3
            CoccoCO2 = max(0.d0,CoccoCO2)  ! July 2022: limitation to zero
#endif

!------------------------------------------------------------------------------
! Calcite dissolution
!------------------------------------------------------------------------------
            Sink_Vel    = Vdet_a* abs(zF(k)) + Vdet

            if (OmegaC_diss) then    ! Calcdiss dependent on carbonate saturation
                Ca        = (0.02128d0/40.078d0) * Sali_depth(k)/1.80655d0 ! Calcium ion concentration [mol/kg], function from varsolver.f90
                CO3_sat   = (kspc_watercolumn(k) / Ca) * rhoSW_watercolumn(k) ! Saturated carbonate ion concentration, converted to [mol/m3]
                calc_diss = calc_diss_omegac * max(zero,(1-(CO3_watercolumn(k)/CO3_sat)))**(calc_diss_exp) ! Dissolution rate scaled by carbonate ratio, after Aumont et al. 2015
#if defined (__3Zoo2Det)
                calc_diss2 = calc_diss
#endif
                calc_diss_ben = calc_diss
            else    ! Calcdiss dependent on depth
                calc_diss = calc_diss_rate * Sink_Vel/20.d0 ! Dissolution rate of CaCO3 scaled by the sinking velocity at the current depth
#if defined (__3Zoo2Det)
                calc_diss2 = calc_diss_rate2* Sink_Vel/20.d0 
#endif
                calc_diss_ben = calc_diss_rate * Sink_Vel/20.0
            endif

!#if defined (__3Zoo2Det)
!            calc_diss2    = calc_diss_rate2 ! Dissolution rate of CaCO3 scaled by the sinking velocity at the current depth seczoo
!#endif
!            calc_diss_ben = calc_diss_rate * Sink_Vel/20.0

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

!< *** Coccolithophores ***
!< ************************
#if defined (__coccos) 
            qlimitFac  = recom_limiter(NMinSlope, NCmin_c, quota_cocco) 
            feLimitFac = Fe/(k_Fe_c + Fe)                               
            qlimitFac  = min(qlimitFac, feLimitFac)                     
            pMax_cocco = P_cm_c * qlimitFac * CoccoTFunc                ! Here the T dependency is changed
#endif
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

!< *** Coccolithophore photosynthesis rate ***
!< *******************************************
#if defined (__coccos) 
            if ( pMax_cocco .lt. tiny .OR. Parave /= Parave .OR. CHL2C_cocco /= CHL2C_cocco) then 
                Cphot_cocco = zero
            else 
                Cphot_cocco = pMax_cocco * (real(one) - exp( -alfa_c * Chl2C_cocco * PARave / pMax_cocco))  
                if (CO2lim) Cphot_cocco = Cphot_cocco * CoccoCO2 ! Added the CO2 dependence 
            end if 
            if (Cphot_cocco .lt. tiny) Cphot_cocco = zero
#endif
!------------------------------------------------------------------------------- 
!< chlorophyll degradation
!-------------------------------------------------------------------------------
            KOchl = deg_Chl
            KOchl_dia = deg_Chl_d
#if defined (__coccos) 
            KOchl_cocco = deg_Chl_c
#endif
        
            if (use_photodamage) then
!< add a minimum value for photodamage
!< *** Phytoplankton Chla loss ***
!< *******************************
                if (pMax .lt. tiny .OR. PARave /= PARave .OR. CHL2C_plast /= CHL2C_plast) then
                    KOchl = deg_Chl*0.1d0
                else
                   !KOchl = deg_Chl*(real(one) - exp(-alfa * CHL2C_plast * PARave / pMax))
                    KOchl = deg_Chl * CHL2C_plast * PARave
                    KOchl = max((deg_Chl*0.1d0), KOchl)
                    KOchl = min(KOChl, 0.3d0)
                end if
!< *** Diatoms Chla loss ***                                                                                                                                
!< *************************
                if (pMax_dia .lt. tiny .OR. PARave /= PARave .OR. CHL2C_plast_dia /= CHL2C_plast_dia) then
                    KOchl_dia = deg_Chl_d*0.1d0
                else
                   !KOchl_dia = deg_Chl_d * (real(one) - exp(-alfa_d * CHL2C_plast_dia * PARave / pMax_dia ))
                    KOchl_dia = deg_Chl_d * CHL2C_plast_dia * PARave
                    KOchl_dia = max((deg_Chl_d*0.1d0), KOchl_dia)
                    KOchl_dia = min(KOChl_dia, 0.3d0)
                end if
!< *** Coccolithophores chla loss ***
!< **********************************
#if defined (__coccos) 
                if (pMax_cocco .lt. tiny .OR. PARave /= Parave .OR. CHL2C_plast_cocco /= CHL2C_plast_cocco) then 
                    KOchl_cocco = deg_Chl_c*0.1d0 
                else 
                   !KOchl_cocco = deg_Chl_c * (real(one) - exp( -alfa_c * CHL2C_plast_cocco * PARave / pMax_cocco ))
                    KOchl_cocco = deg_Chl_c * CHL2C_plast_cocco * PARave
                    KOchl_cocco = max((deg_Chl_c*0.1d0), KOchl_cocco)
                    KOchl_cocco = min(KOChl_cocco, 0.3d0)
                end if
#endif
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
#if defined (__coccos) 
                if (KOchl_cocco /= KOchl_cocco) then 
                    print*,' KOchl_cocco is ', KOchl_cocco
                    print*,' deg_Chl_c is ', deg_Chl_c
                    print*,' alfa_c is ', alfa_c
                    print*,' CHL2C_c is ', CHL2C_plast_cocco
                    print*,' PARave is ', PARave
                    print*,' pMax_c is ', pMax_cocco 
                    stop
                end if  
#endif
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

#if defined (__coccos) 
            V_cm            = V_cm_fact_c                                   
            limitFacN_cocco = recom_limiter(NMaxSlope, quota_cocco, NCmax_c)
            N_assim_cocco   = V_cm * pMax_cocco * NCUptakeRatio_c &         
                               * limitFacN_cocco * DIN/(DIN + k_din_c)      
#endif

            limitFacSi     = recom_limiter(SiMaxSlope, qSiC, SiCmax)  &
                              * limitFacN_dia
            Si_assim       = V_cm_fact_d * P_cm_d * arrFunc * SiCUptakeRatio &
                              * limitFacSi * Si/(Si + k_si)

!-------------------------------------------------------------------------------
!< *** Iron chemistry ***
!< ********************** 
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
#if defined (__coccos) 
            if (PARave .ge. tiny .AND. PARave .eq. PARave) then                             
                ChlSynth_cocco = N_assim_cocco * Chl2N_max_c                &               
                               * min(real(one),Cphot_cocco /(alfa_c * Chl2C_cocco * PARave))
            end if                                                                          
#endif
!-------------------------------------------------------------------------------
!< *** Phytoplankton respiraion rate ***
!< *************************************

!< res_phy: Maintenance respiration rate constant [day−1 ]
!< biosynth: The cost of biosynthesis of N [mmol C mmol N−1 ]

            phyRespRate     = res_phy   * limitFacN     + biosynth * N_assim
            phyRespRate_dia = res_phy_d * limitFacN_dia + biosynth * N_assim_dia + biosynthSi * Si_assim
#if defined (__coccos) 
            phyRespRate_cocco = res_phy_c * limitFacN_cocco + biosynth * N_assim_cocco
#endif

!-------------------------------------------------------------------------------
! Mesozooplankton
!-------------------------------------------------------------------------------
!< Grazing on small phytoplankton, diatoms, coccolithophore (optional), 
!< microzooplankton (optional), slow- and fast-sinking detritus

!< *** Food availability ***
!< *************************
!< pzPhy: Maximum nanophytoplankton preference
!< pzDia: Maximum diatom preference
!< pzCocco: Maximum coccolithophore preference
!< pzDet: Maximum slow-sinking detritus prefence
!< pzDetZ2: Maximum fast-sinking detritus preference
!< pzMicZoo: Maximum microzooplankton preference

            if (REcoM_Grazing_Variable_Preference) then   ! CHECK ONUR
                aux = pzPhy*PhyN + pzDia*DiaN              
                if (Grazing_detritus) aux = aux + PzDet*DetN
#if defined (__3Zoo2Det)
                if (Grazing_detritus) aux = aux + pzDetZ2*DetZ2N ! 2Det
                aux = aux + pzMicZoo*MicZooN                     ! 3Zoo
#endif
#if defined (__coccos) 
                aux = aux + pzCocco*CoccoN
#endif
! ******************************************************************************
                varpzPhy = (pzPhy*PhyN)/aux
                varpzDia = (pzDia*DiaN)/aux
                if (Grazing_detritus) varpzDet = (pzDet*DetN)/aux
#if defined (__3Zoo2Det)
                if (Grazing_detritus) varpzDetZ2 = (pzDetZ2*DetZ2N)/aux ! 2Det
                varpzMicZoo = (pzMicZoo*MicZooN)/aux                    ! 3Zoo
#endif
#if defined (__coccos) 
                varpzCocco = (pzCocco*CoccoN)/aux
#endif
! ******************************************************************************
                fDiaN = varpzDia * DiaN
                fPhyN = varpzPhy * PhyN
                if (Grazing_detritus) fDetN = varpzDet * DetN
#if defined (__3Zoo2Det)
                if (Grazing_detritus) fDetZ2N = varpzDetZ2 * DetZ2N ! 2Det
                fMicZooN = varpzMicZoo * MicZooN                    ! 3Zoo                                         
#endif
#if defined (__coccos) 
                fCoccoN = varpzCocco * CoccoN
#endif
            else ! REcoM_Grazing_Variable_Preference = .false.
                fPhyN = pzPhy * PhyN
                fDiaN = pzDia * DiaN
                if (Grazing_detritus) fDetN = pzDet * DetN
#if defined (__3Zoo2Det)
                if (Grazing_detritus) fDetZ2N = pzDetZ2 * DetZ2N ! 2Det
                fMicZooN = pzMicZoo * MicZooN                    ! 3Zoo
#endif
#if defined (__coccos) 
                fCoccoN = pzCocco * CoccoN 
#endif
            end if ! REcoM_Grazing_Variable_Preference

!< *** Grazing fluxes ***
!< **********************
            food = fPhyN + fDiaN
            if (Grazing_detritus) food = food + fDetN
#if defined (__3Zoo2Det)
            if (Grazing_detritus) food = food + fDetZ2N
            food = food + fMicZooN ! 3Zoo
#endif
#if defined (__coccos) 
            food = food + fCoccoN
#endif
! ******************************************************************************
            foodsq            = food**2
            grazingFlux       = (Graz_max * foodsq)/(epsilonr + foodsq) * HetN * arrFunc
#if defined (__3Zoo2Det)
            grazingFlux       = (Graz_max * foodsq)/(epsilonr + foodsq) * HetN * q10_mes
#endif
            grazingFlux_phy   = grazingFlux * fphyN / food
            grazingFlux_Dia   = grazingFlux * fDiaN / food
            if (Grazing_detritus) grazingFlux_Det   = grazingFlux * fDetN / food
#if defined (__3Zoo2Det)
            if (Grazing_detritus) grazingFlux_DetZ2 = grazingFlux * fDetZ2N / food
            grazingFlux_miczoo = grazingFlux * fMicZooN / food ! 3Zoo
#endif
#if defined (__coccos) 
            grazingFlux_Cocco = grazingFlux * fCoccoN / food
#endif

!< *** Grazing efficiency ***
!< **************************
            grazEff = gfin + 1/(0.2*food + 2)

            grazingFluxcarbon_mes = (grazingFlux_phy   * recipQuota       * grazEff)   &
                                  + (grazingFlux_Dia   * recipQuota_Dia   * grazEff)   

            if (Grazing_detritus) grazingFluxcarbon_mes = grazingFluxcarbon_mes        &
                                  + (grazingFlux_Det   * recipDet         * grazEff)   
#if defined (__3Zoo2Det)
            if (Grazing_detritus) grazingFluxcarbon_mes = grazingFluxcarbon_mes        &
                                  + (grazingFlux_DetZ2 * recipDet2        * grazEff)
            grazingFluxcarbon_mes = grazingFluxcarbon_mes                              &
                                  + (grazingFlux_miczoo * recipQZoo3      * grazEff) ! 3Zoo  
#endif
#if defined (__coccos) 
            grazingFluxcarbon_mes = grazingFluxcarbon_mes                              &
                                  + (grazingFlux_Cocco * recipQuota_Cocco * grazEff)         
#endif

!-------------------------------------------------------------------------------
! Second Zooplankton
!-------------------------------------------------------------------------------
!< Grazing on small phytoplankton, diatoms, coccolithophore (optional), 
!< heterotrophs, microzooplankton (optional), slow- and fast-sinking detritus

!< *** Food availability ***
!< *************************  
!< pzPhy2: Maximum nanophytoplankton preference 
!< pzDia2: Maximum diatom preference
!< pzCocco2: Maximum coccolithophore preference  
!< pzDet2: Maximum slow-sinking detritus prefence
!< pzDetZ22: Maximum fast-sinking detritus preference
!< pzHet: Maximum mesozooplankton preference            
!< pzMicZoo2: Maximum microzooplankton preference

#if defined (__3Zoo2Det)     
            if (REcoM_Grazing_Variable_Preference) then
                aux = pzPhy2 * PhyN + PzDia2 * DiaN + pzHet * HetN 
                if (Grazing_detritus) aux = aux + pzDet2 * DetN + pzDetZ22 * DetZ2N
                aux = aux + pzMicZoo2 * MicZooN ! 3Zoo
#if defined (__coccos) 
                aux = aux + pzCocco2 * CoccoN
#endif
! ******************************************************************************
                varpzPhy2 = (pzPhy2 * PhyN)/aux
                varpzDia2 = (pzDia2 * DiaN)/aux
                varpzMicZoo2 = (pzMicZoo2 * MicZooN)/aux ! 3Zoo

#if defined (__coccos) 
                varpzCocco2 = (pzCocco2 * CoccoN)/aux 
#endif 
                varpzHet = (pzHet * HetN)/aux
                if (Grazing_detritus) then
                    varpzDet2   = (pzDet2 * DetN)/aux
                    varpzDetZ22 = (pzDetZ22 * DetZ2N)/aux
                end if
! ******************************************************************************
                fDiaN2 = varpzDia2 * DiaN
                fPhyN2 = varpzPhy2 * PhyN
                fMicZooN2 = varpzMicZoo2 * MicZooN ! 3Zoo
#if defined (__coccos) 
                fCoccoN2 = varpzCocco2 * CoccoN
#endif
                fHetN = varpzHet * HetN
                if (Grazing_detritus) then
                    fDetN2   = varpzDet2 * DetN
                    fDetZ2N2 = varpzDetZ22 * DetZ2N
                end if
            else ! REcoM_Grazing_Variable_Preference = .false.

                fDiaN2 = pzDia2 * DiaN
                fPhyN2 = pzPhy2 * PhyN
                fMicZooN2 = pzMicZoo2 * MicZooN ! 3Zoo
#if defined (__coccos) 
                fCoccoN2 = pzCocco2 * CoccoN 
#endif
                fHetN = pzHet * HetN
                if (Grazing_detritus) then
                    fDetN2 = pzDet2 * DetN
                    fDetZ2N2 = pzDetZ22 * DetZ2N
                end if
            end if ! REcoM_Grazing_Variable_Preference

!< *** Grazing fluxes ***
!< **********************
            food2 = fPhyN2 + fDiaN2 + fHetN
            if (Grazing_detritus) food2 = food2 + fDetN2 + fDetZ2N2
            food2 = food2 + fMicZooN2 ! 3Zoo
#if defined (__coccos) 
            food2 = food2 + fCoccoN2
#endif
! ******************************************************************************
            foodsq2 = food2**2
            grazingFlux2 = (Graz_max2 * foodsq2)/(epsilon2 + foodsq2) * Zoo2N * arrFuncZoo2

            grazingFlux_phy2 = (grazingFlux2 * fphyN2)/food2
            grazingFlux_Dia2 = (grazingFlux2 * fDiaN2)/food2
            grazingFlux_miczoo2 = (grazingFlux2 * fMicZooN2)/food2 ! 3Zoo

#if defined (__coccos) 
            grazingFlux_Cocco2 = (grazingFlux2 * fCoccoN2)/food2
#endif
            grazingFlux_het2 = (grazingFlux2 * fHetN)/food2
            if (Grazing_detritus) then
                grazingFlux_Det2 = (grazingFlux2 * fDetN2)/food2
                grazingFlux_DetZ22 = (grazingFlux2 * fDetZ2N2)/food2
            end if

            grazingFluxcarbonzoo2 = (grazingFlux_phy2    * recipQuota       * grazEff2)    &
                                  + (grazingFlux_Dia2    * recipQuota_Dia   * grazEff2)    &
                                  + (grazingFlux_het2    * recipQZoo        * grazEff2)    
            if (Grazing_detritus) then
                grazingFluxcarbonzoo2 = grazingFluxcarbonzoo2 +                            &
                                  + (grazingFlux_Det2    * recipDet         * grazEff2)    &
                                  + (grazingFlux_DetZ22  * recipDet2        * grazEff2)
            end if
            grazingFluxcarbonzoo2 = grazingFluxcarbonzoo2 +                                &
                                  + (grazingFlux_miczoo2 * recipQZoo3       * grazEff2) ! 3Zoo
#if defined (__coccos) 
            grazingFluxcarbonzoo2 = grazingFluxcarbonzoo2 +                                &
                                  + (grazingFlux_Cocco2  * recipQuota_Cocco * grazEff2)
#endif

!-------------------------------------------------------------------------------
! Third Zooplankton (Microzooplankton)
!-------------------------------------------------------------------------------
!< Grazing on small phytoplankton, diatoms and coccolithophore (optional)

!< *** Food availability ***
!< *************************
!< pzPhy3: Maximum nanophytoplankton preference     
!< pzDia3: Maximum diatom preference 
!< pzCocco3: Maximum coccolithophore preference

            if (REcoM_Grazing_Variable_Preference) then
                aux = pzPhy3 * PhyN + pzDia3 * DiaN
#if defined (__coccos)
                aux = aux + pzCocco3 * CoccoN
#endif
! ******************************************************************************
                varpzPhy3 = (pzPhy3 * PhyN)/aux
                varpzDia3 = (pzDia3 * DiaN)/aux
#if defined (__coccos)
                varpzCocco3 = (pzCocco3 * CoccoN)/aux
#endif
! ******************************************************************************
                fPhyN3 = varpzPhy3 * PhyN
                fDiaN3 = varpzDia3 * DiaN
#if defined (__coccos)
                fCoccoN3 = varpzCocco3 * CoccoN
#endif
            else ! REcoM_Grazing_Variable_Preference = .false.

                fPhyN3 = pzPhy3 * PhyN
                fDiaN3 = pzDia3 * DiaN
#if defined (__coccos)
                fCoccoN3 = pzCocco3 * CoccoN
#endif
            endif !REcoM_Grazing_Variable_Preference

!< *** Grazing fluxes ***
!< **********************
            food3 = fPhyN3 + fDiaN3
#if defined (__coccos)
            food3 = food3 + fCoccoN3
#endif
! ******************************************************************************
            foodsq3 = food3**2
            grazingFlux3 = (Graz_max3 * foodsq3)/(epsilon3 + foodsq3) * MicZooN * q10_mic
            grazingFlux_phy3 = (grazingFlux3 * fphyN3)/food3
            grazingFlux_Dia3 = (grazingFlux3 * fDiaN3)/food3
#if defined (__coccos)
            grazingFlux_Cocco3 = (grazingFlux3 * fCoccoN3)/food3
#endif
#endif

!-------------------------------------------------------------------------------
!< Heterotrophic respiration is assumed to drive zooplankton back to 
!< Redfield C:N if their C:N becomes higher than Redfield
!< res_het: Timescale for zooplankton respiration [day−1 ]

            if (het_resp_noredfield) then
#if defined (__3Zoo2Det)
                HetRespFlux = res_het * q10_mes_res * HetC ! 3Zoo
#else
                HetRespFlux = res_het * arrFunc * HetC ! tau * f_T [HetC]
#endif
            else
                HetRespFlux = recip_res_het * arrFunc * (hetC * recip_hetN_plus - redfield) * HetC
                HetRespFlux = max(zero, HetRespFlux)  !!!!!!!! CHECK Judith Valid for het_resp_noredfield case as well ???????? Then move it below
            endif

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
!< Zooplanton mortality (Quadratic)

            hetLossFlux = loss_het * HetN * HetN

#if defined (__3Zoo2Det)
!-------------------------------------------------------------------------------
!< Second zooplankton respiration 

            call krill_resp(n,mesh)

            if((grazingFluxcarbonzoo2/Zoo2C) <= 0.1)then
                res_zoo2_f = 0.1*(grazingFluxcarbonzoo2/Zoo2C*100)
            else
                res_zoo2_f = 1.
            end if
            recip_res_zoo22 = res_zoo2*(1.+ res_zoo2_f + res_zoo2_a)
            Zoo2RespFlux = recip_res_zoo22 * Zoo2C                   
!-------------------------------------------------------------------------------
!< Second zooplankton mortality (Quadratic)

            Zoo2LossFlux = loss_zoo2 * zoo2N * zoo2N

!-------------------------------------------------------------------------------
!< Second zooplankton fecal pellets

            Zoo2fecalloss_n = fecal_rate_n * grazingFlux2
            Zoo2fecalloss_c = fecal_rate_c * grazingFluxcarbonzoo2

!-------------------------------------------------------------------------------
!< Mesozooplankton fecal pellets

            mesfecalloss_n = fecal_rate_n_mes * grazingFlux
            mesfecalloss_c = fecal_rate_c_mes * grazingFluxcarbon_mes
!------------------------------------------------------------------------------- 
! Third zooplankton, microzooplankton, respiration ! 3Zoo

            MicZooRespFlux =  res_miczoo * q10_mic_res * MicZooC
!-------------------------------------------------------------------------------
! Third zooplankton, microzooplankton, mortality  (Quadratic) ! 3Zoo

            MicZooLossFlux = loss_miczoo * MicZooN * MicZooN
#endif

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

#if defined (__3Zoo2Det)
            aggregationrate = aggregationrate + agg_PD * DetZ2N ! 2Det
#endif
#if defined (__coccos)
            aggregationrate = aggregationrate + agg_PP * CoccoN
#endif
    
!-------------------------------------------------------------------------------
! Calcification
!-------------------------------------------------------------------------------
! Without coccolithophores, calcification is performed by a fraction of small phytoplankton

#if defined (__coccos) 
            if (Temp(k) < 10.6) then                 ! (PICPOC definition after Krumhardt et al. 2017, 2019; Temp(k) because we need degC here)
                PICPOCtemp = 0.104d0 * Temp(k) - 0.108d0
            else
                PICPOCtemp = 1.0d0
            end if
            PICPOCtemp = max(tiny,PICPOCtemp)

            PICPOCCO2 = a_co2_calc * HCO3_watercolumn(k) * Cunits / (b_co2_calc + HCO3_watercolumn(k) * Cunits) - exp(-c_co2_calc * CO2_watercolumn(k) * Cunits) - d_co2_calc * 10.**(-pH_watercolumn(k))
            PICPOCCO2 = min(PICPOCCO2,3.d0) ! April 2022: limitation to 3
            PICPOCCO2 = max(0.d0,PICPOCCO2) ! July 2022: limitation to zero

            PICPOCN = -0.31 * (DIN/(DIN + k_din_c)) + 1.31
            PICPOCN = max(tiny,PICPOCN)
    
            calcification = 1.d0 * Cphot_cocco * CoccoC * PICPOCtemp * PICPOCN
            if (CO2lim) calcification = calcification * PICPOCCO2

#else
!< calc_prod_ratio: Calcite production ratio, dimensionless
            calcification = calc_prod_ratio * Cphot * PhyC ! Z in equations
#endif

            calc_loss_agg = aggregationrate * PhyCalc

#if defined (__coccos) 
!< *** Coccolithophores ***
!< ************************
            aux = recipQuota_Cocco/(CoccoC + tiny) * PhyCalc
            calc_loss_gra  = grazingFlux_Cocco  * aux 
#if defined (__3Zoo2Det)
            calc_loss_gra2 = grazingFlux_Cocco2 * aux
            calc_loss_gra3 = grazingFlux_Cocco3 * aux ! 3Zoo
#endif

#else
!< *** Small phytoplankton ***
!< ***************************
            aux = recipQuota/(PhyC + tiny) * PhyCalc
            calc_loss_gra  = grazingFlux_phy  * aux
#if defined (__3Zoo2Det)
            calc_loss_gra2 = grazingFlux_phy2 * aux
            calc_loss_gra3 = grazingFlux_phy3 * aux ! 3Zoo
#endif
#endif

            if (ciso) then
                calcification_13 = calc_prod_ratio * Cphot * PhyC_13 * alpha_calc_13
                calcification_13 = calcification   * alpha_calc_13
                calc_loss_agg_13 = aggregationRate * PhyCalc_13
                calc_loss_gra_13 = grazingFlux_phy * recipQuota_13/(PhyC_13 + tiny) * PhyCalc_13
                if (ciso_14 .and. ciso_organic_14) then
                    calcification_14 = calc_prod_ratio * Cphot * PhyC_14 * alpha_calc_14
                    calc_loss_agg_14 = aggregationRate * PhyCalc_14
                    calc_loss_gra_14 = grazingFlux_phy * recipQuota_14/(PhyC_14 + tiny) * PhyCalc_14
                end if
            end if

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
#if defined (__coccos)
        - N_assim_Cocco            * CoccoC  &  ! --> N assimilation Coccolithophore
#endif
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
#if defined (__coccos)
        - Cphot_Cocco                     * CoccoC    & ! --> Coccolithophore photosynthesis
        + phyRespRate_Cocco               * CoccoC    & ! --> Coccolithophore respiration
#endif
        + rho_C1 * arrFunc * O2Func       * EOC       & ! --> Remineralization of DOC
        + HetRespFlux                                 & ! --> Mesozooplankton respiration
#if defined (__3Zoo2Det)                     
        + Zoo2RespFlux                                & ! --> Macrozooplankton respiration            
        + MicZooRespFlux                              & ! --> Microzooplankton respiration
#endif          
        + calc_diss                       * DetCalc   & ! --> Calcite dissolution from slow-sinking detritus 
        + calc_loss_gra  * calc_diss_guts             & ! --> Additional dissolution in mesozooplankton guts
#if defined (__3Zoo2Det)
        + calc_loss_gra2 * calc_diss_guts             & ! --> Additional dissolution in macrozooplankton guts
        + calc_loss_gra3 * calc_diss_guts             & ! --> Additional dissolution in microzooplankton guts
        + calc_diss2                      * DetZ2Calc & ! --> Calcite dissolution from fast-sinking detritus
#endif
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
#if defined (__coccos)
        + 1.0625 * N_assim_Cocco                 * CoccoC  & 
#endif
        - 1.0625 * rho_N * arrFunc * O2Func      * DON     &     
        + 2.d0 * calc_diss                       * DetCalc &
        + 2.d0 * calc_loss_gra  * calc_diss_guts           &
#if defined (__3Zoo2Det)
        + 2.d0 * calc_loss_gra2 * calc_diss_guts           &
        + 2.d0 * calc_loss_gra3 * calc_diss_guts           & ! 3Zoo
        + 2.d0 * calc_diss2            * DetZ2Calc         &
#endif
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
#if defined (__3Zoo2Det)
        - grazingFlux_phy2                                 & 
        - grazingFlux_phy3                                 & ! 3Zoo    
#endif
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
#if defined (__3Zoo2Det)
        - grazingFlux_phy2 * recipQuota                    &
        - grazingFlux_phy3 * recipQuota                    & ! 3Zoo
#endif
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
#if defined (__3Zoo2Det)
        - grazingFlux_phy2 * Chl2N                         & 
        - grazingFlux_phy3 * Chl2N                         & ! 3Zoo
#endif
                                                          ) * dt_b + sms(k,ipchl)

!< *** Slow-sinking Detritus ***
!< *****************************

!____________________________________________________________
! Detritus N
    if (Grazing_detritus) then
#if defined (__3Zoo2Det)
        sms(k,idetn)       = (                             &
            + grazingFlux_phy3                             &
            - grazingFlux_phy3   * grazEff3                &
            + grazingFlux_dia3                             &
            - grazingFlux_dia3   * grazEff3                &
#if defined (__coccos)
            + grazingFlux_Cocco3                           &
            - grazingFlux_Cocco3 * grazEff3                &
            + aggregationRate               * CoccoN       & 
#endif
            - grazingFlux_Det    * grazEff                 &   ! --> grazing of first zoo (meso) on first detritus
            - grazingFlux_Det2   * grazEff2                &   ! --> grazing of second zoo on first detritus
            + aggregationRate               * PhyN         &
            + aggregationRate               * DiaN         &
            + miczooLossFlux                               &
            - reminN * arrFunc * O2Func     * DetN         & ! O2remin
                                                          ) * dt_b + sms(k,idetn)
#else
        sms(k,idetn)       = (                             &
	    + grazingFlux_phy                              &
            - grazingFlux_phy   * grazEff                  &
            + grazingFlux_dia                              &
            - grazingFlux_dia   * grazEff                  &
#if defined (__coccos)
            + grazingFlux_Cocco                            &
            - grazingFlux_Cocco * grazEff                  &
            + aggregationRate               * CoccoN       &
#endif
            - grazingFlux_Det   * grazEff                  & ! Sloppy feeding is thought because of grazing flux multiplied with grazeff 
            - grazingFlux_Det2  * grazEff2                 &
            + aggregationRate               * PhyN         &
            + aggregationRate               * DiaN         &
            + hetLossFlux                                  &
            - reminN * arrFunc * O2Func     * DetN         & ! O2remin
                                                          ) * dt_b + sms(k,idetn)
#endif
   else
#if defined (__3Zoo2Det)
        sms(k,idetn)       = (                             &
            + grazingFlux_phy3                             &
            + grazingFlux_dia3                             &
#if defined (__coccos)
            + grazingFlux_Cocco3                           &
            + aggregationRate               * CoccoN       &
#endif
            - grazingFlux        * grazEff3                &
            + aggregationRate               * PhyN         &
            + aggregationRate               * DiaN         &
            + miczooLossFlux                               &
            - reminN * arrFunc * O2Func     * DetN         & ! O2remin
                                                          ) * dt_b + sms(k,idetn)
#else
        sms(k,idetn)       = (                             &
            + grazingFlux_phy                              &
            + grazingFlux_dia                              &
#if defined (__coccos)
            + grazingFlux_Cocco                            &
            + aggregationRate               * CoccoN       & 
#endif
            - grazingFlux        * grazEff                 &
            + aggregationRate               * PhyN         &
            + aggregationRate               * DiaN         &
            + hetLossFlux                                  &
            - reminN * arrFunc * O2Func     * DetN         & ! O2remin
                                                          ) * dt_b + sms(k,idetn)
#endif
   end if

!____________________________________________________________
! Detritus C
    if (Grazing_detritus) then
#if defined (__3Zoo2Det)
        sms(k,idetc)       = (                                 &
            + grazingFlux_phy3 * recipQuota                    &
            - grazingFlux_phy3 * recipQuota         * grazEff3 &
            + grazingFlux_Dia3 * recipQuota_Dia                & 
            - grazingFlux_Dia3 * recipQuota_Dia     * grazEff3 &
#if defined (__coccos)
            + grazingFlux_Cocco3 * recipQuota_Cocco            &
            - grazingFlux_Cocco3 * recipQuota_Cocco * grazEff3 &
            + aggregationRate                       * CoccoC   &
#endif
            - grazingFlux_Det  * recipDet  * grazEff           &
            - grazingFlux_Det2 * recipDet2 * grazEff2          &
            + aggregationRate                         * PhyC   &
            + aggregationRate                         * DiaC   &
            + miczooLossFlux   * recipQZoo3                    &
            - reminC * arrFunc * O2Func   * DetC               & ! O2remin
                                                              ) * dt_b + sms(k,idetc)
#else
        sms(k,idetc)       = (                                 &
            + grazingFlux_phy * recipQuota                     &
            - grazingFlux_phy * recipQuota         * grazEff   &
            + grazingFlux_Dia * recipQuota_Dia                 &
            - grazingFlux_Dia * recipQuota_Dia     * grazEff   &
#if defined (__coccos)
            + grazingFlux_Cocco * recipQuota_Cocco             &
            - grazingFlux_Cocco * recipQuota_Cocco * grazEff   &
            + aggregationRate                      * CoccoC    &
#endif
            - grazingFlux_Det  * recipDet  * grazEff           &
            - grazingFlux_Det2 * recipDet2 * grazEff           &
            + aggregationRate              * phyC              &
            + aggregationRate              * DiaC              &
            + hetLossFlux      * recipQZoo                     &
            - reminC * arrFunc * O2Func    * DetC              & ! O2remin
                                                              ) * dt_b + sms(k,idetc)
#endif
   else
#if defined (__3Zoo2Det)
        sms(k,idetc)       = (                                 &
            + grazingFlux_phy3   * recipQuota                  &
            - grazingFlux_phy3   * recipQuota       * grazEff3 &
            + grazingFlux_Dia3   * recipQuota_Dia              &
            - grazingFlux_Dia3   * recipQuota_Dia   * grazEff3 &
#if defined (__coccos)
            + grazingFlux_Cocco3 * recipQuota_Cocco            &
            - grazingFlux_Cocco3 * recipQuota_Cocco * grazEff3 &
            + aggregationRate                       * CoccoC   &
#endif
            + aggregationRate                       * PhyC     &
            + aggregationRate                       * DiaC     &
            + miczooLossFlux     * recipQZoo3                  &
            - reminC * arrFunc * O2Func    * DetC              & ! O2remin
                                                              ) * dt_b + sms(k,idetc)
#else
        sms(k,idetc)       = (                                 &
            + grazingFlux_phy   * recipQuota                   &
            - grazingFlux_phy   * recipQuota        * grazEff  &
            + grazingFlux_Dia   * recipQuota_Dia               &
            - grazingFlux_Dia   * recipQuota_Dia    * grazEff  &
#if defined (__coccos)
            + grazingFlux_Cocco * recipQuota_Cocco             & 
            - grazingFlux_Cocco * recipQuota_Cocco  * grazEff  & 
            + aggregationRate                       * CoccoC   & 
#endif
            + aggregationRate                       * phyC     &
            + aggregationRate                       * DiaC     &
            + hetLossFlux       * recipQZoo                    &
            - reminC * arrFunc * O2Func    * DetC              & ! O2remin
                                                              ) * dt_b + sms(k,idetc)
#endif
   end if

!< *** Mesozooplankton ***
!< ***********************

!____________________________________________________________
!< Heterotrophic N
        sms(k,ihetn)       = (                                 &
    	    + grazingFlux      * grazEff                       & ! --> Grazing on phytoplankton -> okay, because of recipQuota
#if defined (__3Zoo2Det)
            - grazingFlux_het2                                 &
            - Mesfecalloss_n                                   & ! 3Zoo
#endif
     	    - hetLossFlux                                      & ! --> Mortality
     	    - lossN_z                      * HetN              & ! --> Excretion of DON
                                                              ) * dt_b + sms(k,ihetn)  
!____________________________________________________________
!< Heterotrophic C
        if (Grazing_detritus) then
            sms(k,ihetc)      = (                                 &
     	        + grazingFlux_phy    * recipQuota       * grazEff & ! --> Grazing on small phytoplankton
     	        + grazingFlux_Dia    * recipQuota_Dia   * grazEff & ! --> Grazing on diatom
#if defined (__coccos)
                + grazingFlux_Cocco  * recipQuota_Cocco * grazEff & 
#endif
#if defined (__3Zoo2Det)
                + grazingFlux_miczoo * recipQZoo3       * grazEff & ! 3Zoo
                + grazingFlux_DetZ2  * recipDet2        * grazEff &
                - grazingFlux_het2   * recipQZoo                  &
                - Mesfecalloss_c                                  & ! 3Zoo
#endif
                + grazingFlux_Det    * recipDet         * grazEff & ! --> Grazing on detritus
     	        - hetLossFlux        * recipQZoo                  & ! --> Mortality loss
     	        - lossC_z                               * HetC    & ! --> Excretion loss
      	        - hetRespFlux                                     & ! --> REspiration loss
                                                                 ) * dt_b + sms(k,ihetc)
        else
            sms(k,ihetc)      = (                                 &
     	        + grazingFlux_phy    * recipQuota       * grazEff &
     	        + grazingFlux_Dia    * recipQuota_Dia   * grazEff &
#if defined (__coccos)
                + grazingFlux_Cocco  * recipQuota_Cocco * grazEff &
#endif
#if defined (__3Zoo2Det)
                + grazingFlux_miczoo * recipQZoo3       * grazEff & ! 3Zoo
                - grazingFlux_het2   * recipQZoo                  &
                - Mesfecalloss_c                                  & ! 3Zoo
#endif
     	        - hetLossFlux        * recipQZoo                  &
     	        - lossC_z                               * HetC    &
     	        - hetRespFlux                                     & 
                                                                 ) * dt_b + sms(k,ihetc)
        endif

!< *** Macrozooplankton ***
!< ************************

#if defined (__3Zoo2Det)
!____________________________________________________________
!< Second Zooplankton N                                                                                              
    sms(k,izoo2n)       = (                        &
        + grazingFlux2     * grazEff2              &
        - Zoo2LossFlux                             &
        - lossN_z2                    * Zoo2N      &
        - Zoo2fecalloss_n                          & 
                                                  ) * dt_b + sms(k,izoo2n)

!____________________________________________________________
!< Second Zooplankton C                                                                                  
    if (Grazing_detritus) then
             
        sms(k,izoo2c)      = (                                 &
            + grazingFlux_phy2   * recipQuota       * grazEff2 &
            + grazingFlux_Dia2   * recipQuota_Dia   * grazEff2 &
#if defined (__coccos)
            + grazingFlux_Cocco2 * recipQuota_Cocco * grazEff2 &
#endif
            + grazingFlux_het2   * recipQZoo        * grazEff2 &
            + grazingFlux_miczoo2* recipQZoo3       * grazEff2 & ! 3Zoo
            + grazingFlux_Det2   * recipDet         * grazEff2 &
            + grazingFlux_DetZ22 * recipDet2        * grazEff2 &
            - zoo2LossFlux       * recipQZoo2                  &
            - lossC_z2                              * Zoo2C    &
            - Zoo2RespFlux                                     &
            - Zoo2fecalloss_c                                  &
                                                              ) * dt_b + sms(k,izoo2c)  
    else
        sms(k,izoo2c)      = (                                 &
            + grazingFlux_phy2   * recipQuota       * grazEff2 &
            + grazingFlux_Dia2   * recipQuota_Dia   * grazEff2 &
#if defined (__coccos)
            + grazingFlux_Cocco2 * recipQuota_Cocco * grazEff2 & 
#endif
            + grazingFlux_het2   * recipQZoo        * grazEff2 &
            + grazingFlux_miczoo2* recipQZoo3       * grazEff2 & ! 3Zoo
            - zoo2LossFlux       * recipQZoo2                  &
            - lossC_z2                              * Zoo2C    &
            - Zoo2RespFlux                                     &
            - Zoo2fecalloss_c                                  &
                                                              ) * dt_b + sms(k,izoo2c)
    end if

!< *** Microzooplankton ***
!< ************************

!____________________________________________________________
!< Third Zooplankton N             
    sms(k,imiczoon)       = (                      &
        + grazingFlux3        * grazEff3           &
        - grazingFlux_miczoo                       &
        - grazingFlux_miczoo2                      &
        - MicZooLossFlux                           &
        - lossN_z3                       * MicZooN &
                                                  ) * dt_b + sms(k,imiczoon)

!____________________________________________________________
!< Third Zooplankton C            
    sms(k,imiczooc)      = (                                &
        + grazingFlux_phy3    * recipQuota       * grazEff3 &
        + grazingFlux_Dia3    * recipQuota_Dia   * grazEff3 &
#if defined (__coccos)
        + grazingFlux_Cocco3  * recipQuota_Cocco * grazEff3 & 
#endif
        - MicZooLossFlux      * recipQZoo3                  &
        - grazingFlux_miczoo  * recipQZoo3                  & 
        - grazingFlux_miczoo2 * recipQZoo3                  &
        - lossC_z3                               * MicZooC  &
        - MicZooRespFlux                                    &
                                                           ) * dt_b + sms(k,imiczooc)

!< *** Fast-sinking Detritus ***
!< *****************************

!____________________________________________________________
!< Second Zooplankton Detritus N
    if (Grazing_detritus) then
        sms(k,idetz2n)       = (                     &
            + grazingFlux_phy2                       &
            - grazingFlux_phy2    * grazEff2         &
            + grazingFlux_dia2                       &
            - grazingFlux_dia2    * grazEff2         &
#if defined (__coccos)
            + grazingFlux_Cocco                      &
            - grazingFlux_Cocco   * grazEff          &
            + grazingFlux_Cocco2                     &
            - grazingFlux_Cocco2  * grazEff2         &
#endif
            + grazingFlux_het2                       &
            - grazingFlux_het2    * grazEff2         &
            + grazingFlux_miczoo2                    &
            - grazingFlux_miczoo2 * grazEff2         &
            + grazingFlux_phy                        &
            - grazingFlux_phy     * grazEff          &
            + grazingFlux_dia                        &
            - grazingFlux_dia     * grazEff          &
            + grazingFlux_miczoo                     &
            - grazingFlux_miczoo  * grazEff          &
            - grazingFlux_DetZ2   * grazEff          &
            - grazingFlux_DetZ22  * grazEff2         &
            + Zoo2LossFlux                           &
            + hetLossFlux                            &
            + Zoo2fecalloss_n                        &
            + Mesfecalloss_n                         &
            - reminN * arrFunc * O2Func   * DetZ2N   & ! O2remin
                                                    ) * dt_b + sms(k,idetz2n)
    else
        sms(k,idetz2n)       = (                     &
            + grazingFlux_phy2                       &
            + grazingFlux_dia2                       &
#if defined (__coccos)
            + grazingFlux_Cocco                      &
            + grazingFlux_Cocco2                     &
#endif
            + grazingFlux_het2                       &
            + grazingFlux_miczoo2                    &
            - grazingFlux2        * grazEff2         &
            + grazingFlux_phy                        &
            + grazingFlux_dia                        &
            + grazingFlux_miczoo                     &
            - grazingFlux         * grazEff          &
            + Zoo2LossFlux                           &
            + hetLossFlux                            &
            + Zoo2fecalloss_n                        &
            + Mesfecalloss_n                         &
            - reminN * arrFunc * O2Func   * DetZ2N   & ! O2remin
                                                    ) * dt_b + sms(k,idetz2n)
     end if

!____________________________________________________________
!< Second Zooplankton Detritus C
    if (Grazing_detritus) then
        sms(k,idetz2c)       = (                                 &
            + grazingFlux_phy2    * recipQuota                  &
            - grazingFlux_phy2    * recipQuota       * grazEff2 &
            + grazingFlux_Dia2    * recipQuota_Dia              &
            - grazingFlux_Dia2    * recipQuota_Dia   * grazEff2 &
#if defined (__coccos)
            + grazingFlux_Cocco   * recipQuota_Cocco            &
            - grazingFlux_Cocco   * recipQuota_Cocco * grazEff  &
            + grazingFlux_Cocco2  * recipQuota_Cocco            & 
            - grazingFlux_Cocco2  * recipQuota_Cocco * grazEff2 & 
#endif
            + grazingFlux_het2    * recipQZoo                   &
            - grazingFlux_het2    * recipQZoo        * grazEff2 &
            + grazingFlux_miczoo2 * recipQZoo3                  &
            - grazingFlux_miczoo2 * recipQZoo3       * grazEff2 &
            + grazingFlux_phy     * recipQuota                  &
            - grazingFlux_phy     * recipQuota       * grazEff  &
            + grazingFlux_Dia     * recipQuota_Dia              &
            - grazingFlux_Dia     * recipQuota_Dia   * grazEff  &
            + grazingFlux_miczoo  * recipQZoo3                  &
            - grazingFlux_miczoo  * recipQZoo3       * grazEff  &
            - grazingFlux_DetZ2   * recipDet2        * grazEff  &
            - grazingFlux_DetZ22  * recipDet2        * grazEff2 &
            + Zoo2LossFlux        * recipQZoo2                  &
            + hetLossFlux         * recipQZoo                   &
            + Zoo2fecalloss_c                                   &
            + Mesfecalloss_c                                    &
            - reminC * arrFunc * O2Func   * DetZ2C              & ! O2remin
                                                               ) * dt_b + sms(k,idetz2c)
    else
        sms(k,idetz2c)       = (                                &
            + grazingFlux_phy2    * recipQuota                  &
            - grazingFlux_phy2    * recipQuota       * grazEff2 &
            + grazingFlux_Dia2    * recipQuota_Dia              &
            - grazingFlux_Dia2    * recipQuota_Dia   * grazEff2 &
#if defined (__coccos)
            + grazingFlux_Cocco   * recipQuota_Cocco            &
            - grazingFlux_Cocco   * recipQuota_Cocco * grazEff  &
            + grazingFlux_Cocco2  * recipQuota_Cocco            & 
            - grazingFlux_Cocco2  * recipQuota_Cocco * grazEff2 & 
#endif  
            + grazingFlux_het2    * recipQZoo                   &
            - grazingFlux_het2    * recipQZoo        * grazEff2 &
            + grazingFlux_miczoo2 * recipQZoo3                  &
            - grazingFlux_miczoo2 * recipQZoo3       * grazEff2 &
            + grazingFlux_phy     * recipQuota                  &
            - grazingFlux_phy     * recipQuota       * grazEff  &
            + grazingFlux_Dia     * recipQuota_Dia              &
            - grazingFlux_Dia     * recipQuota_Dia   * grazEff  &
            + grazingFlux_miczoo  * recipQZoo3                  &
            - grazingFlux_miczoo  * recipQZoo3       * grazEff  &
            + Zoo2LossFlux        * recipQZoo2                  &
            + hetLossFlux         * recipQZoo                   &
            + Zoo2fecalloss_c                                   &
            + Mesfecalloss_c                                    &
            - reminC * arrFunc * O2Func    * DetZ2C             & ! O2remin
                                                               ) * dt_b + sms(k,idetz2c)
     end if

!____________________________________________________________
!< Second Zooplankton Detritus Si 
    sms(k,idetz2si)     = (                 &
        + grazingFlux_dia2 * qSiN           &  ! --> qSin convert N to Si
        + grazingFlux_dia  * qSiN           &
        - reminSiT                * DetZ2Si &
                                          ) * dt_b + sms(k,idetz2si)

!____________________________________________________________
!< Second Zooplankton Detritus calcite   
    sms(k,idetz2calc)   = (               &
        + calc_loss_gra2                  &
        - calc_loss_gra2 * calc_diss_guts &
        + calc_loss_gra                   &
        - calc_loss_gra  * calc_diss_guts &
        - calc_diss2     * DetZ2Calc      &
                                         ) * dt_b + sms(k,idetz2calc)
#endif 

!< *** DOM ***
!< ***********

!____________________________________________________________
!< DON (Extracellular organic N)

    sms(k,idon)      = (                      &
        + lossN * limitFacN         * phyN    &
        + lossN_d * limitFacN_Dia   * DiaN    &
#if defined (__coccos)
        + lossN_c * limitFacN_Cocco * CoccoN  & 
#endif
        + reminN * arrFunc * O2Func * DetN    & 
        + lossN_z                   * HetN    &
#if defined (__3Zoo2Det)
        + reminN * arrFunc * O2Func * DetZ2N  &
        + lossN_z2                  * Zoo2N   &
        + lossN_z3                  * MicZooN & ! 3Zoo
#endif 
        - rho_N * arrFunc * O2Func  * DON     & ! O2remin
                                             ) * dt_b + sms(k,idon)

!____________________________________________________________
!< EOC

    sms(k,idoc)       = (                     &
        + lossC * limitFacN         * phyC    &
        + lossC_d * limitFacN_dia   * DiaC    &
#if defined (__coccos)
        + lossC_c * limitFacN_cocco * CoccoC  &
#endif
        + reminC * arrFunc * O2Func * DetC    &
        + lossC_z                   * HetC    &
#if defined (__3Zoo2Det)
        + reminC * arrFunc * O2Func * DetZ2C  & 
        + lossC_z2                  * Zoo2C   &
        + lossC_z3                  * MicZooC & ! 3Zoo
#endif
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
#if defined (__3Zoo2Det)
        - grazingFlux_Dia2                       &
        - grazingFlux_Dia3                       & ! 3Zoo
#endif
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
#if defined (__3Zoo2Det)
        - grazingFlux_dia2 * recipQuota_dia         &
        - grazingFlux_dia3 * recipQuota_dia         & ! 3Zoo
#endif
     	                                           ) * dt_b + sms(k,idiac)

!____________________________________________________________
!< Diatom Chl

    sms(k,idchl)      = (                           &
        + chlSynth_dia                     * DiaC   & ! --> Chl a synthesis
        - KOchl_dia                        * DiaChl & ! --> Degradation loss
        - aggregationRate                  * DiaChl & ! --> Aggregation loss
        - grazingFlux_dia  * Chl2N_dia              & ! --> Grazing loss
#if defined (__3Zoo2Det)
        - grazingFlux_dia2 * Chl2N_dia              &          
        - grazingFlux_dia3 * Chl2N_dia              & ! 3Zoo
#endif
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
#if defined (__3Zoo2Det)
        - grazingFlux_dia2 * qSiN                   &
        - grazingFlux_dia3 * qSiN                   & ! 3Zoo
#endif
                                                   ) * dt_b + sms(k,idiasi)

!< *** Coccolithophore ***
!< ***********************

#if defined (__coccos)
!____________________________________________________________
!< Coccolithophore N
    sms(k,icocn)      = (                    &                        
        + N_assim_cocco             * CoccoC &                        
        - lossN_c * limitFacN_cocco * CoccoN &                        
        - aggregationRate           * CoccoN &                        
        - grazingFlux_Cocco                  &                       
#if defined (__3Zoo2Det)
        - grazingFlux_Cocco2                 &                        
        - grazingFlux_Cocco3                 &  ! 3Zoo 
#endif
                                            ) * dt_b + sms(k,icocn)    

!____________________________________________________________
!< Coccolithophore C

    sms(k,icocc)      = (                            & 
      + Cphot_cocco               * CoccoC           &
      - lossC_c * limitFacN_cocco * CoccoC           &
      - phyRespRate_cocco         * CoccoC           &
      - aggregationRate           * CoccoC           &
      - grazingFlux_cocco         * recipQuota_cocco &
#if defined (__3Zoo2Det)
      - grazingFlux_Cocco2        * recipQuota_cocco &
      - grazingFlux_Cocco3        * recipQuota_cocco & ! 3Zoo
#endif
                                                    ) * dt_b + sms(k,icocc)

   if(sms(k,icocc)>100) then   
       print*,'ERROR: strange CoccoC !'
       print*,'k= ', k
       print*,'dt= ', dt
       print*,'dt_b= ', dt_b
       print*,'state(k,icocc): ', state(k,icocc)
       print*,'sms CoccoC: ', CoccoC
       print*,'sms CoccoN: ', CoccoN
       print*,'sms Cphot cocco: ', Cphot_cocco*CoccoC
       print*,'sms lossC_c: ', lossC_c
       print*,'sms limitFacN_cocco: ', limitFacN_cocco
       print*,'sms phyRespRate_cocco: ', phyRespRate_cocco
       print*,'sms grazingFlux_cocco: ', grazingFlux_cocco
       print*,'sms grazingFlux_cocco2: ', grazingFlux_Cocco2
       print*,'sms grazingFlux_cocco3: ', grazingFlux_Cocco3
       print*,'sms recipQuota_cocco: ', recipQuota_cocco
       call par_ex
       stop
    endif

!____________________________________________________________
!< Coccolithophore Chl 

    sms(k,icchl)      = (                           & 
      + ChlSynth_cocco                   * CoccoC   & 
      - KOchl_cocco                      * CoccoChl &
      - aggregationRate                  * CoccoChl & 
      - grazingFlux_cocco  * Chl2N_cocco            &
#if defined (__3Zoo2Det)
      - grazingFlux_Cocco2 * Chl2N_cocco            & 
      - grazingFlux_Cocco3 * Chl2N_cocco            & ! 3Zoo
#endif
                                                   ) * dt_b + sms(k,icchl)
#endif

!< *** Silicate ***
!< ****************

!____________________________________________________________
!< Detritus Si
#if defined (__3Zoo2Det)
    sms(k,idetsi)     = (                          &
        + aggregationRate                  * DiaSi &
        + lossN_d          * limitFacN_dia * DiaSi &
        + grazingFlux_dia3 * qSiN                  &
        - reminSiT                         * DetSi &
                                                  ) * dt_b + sms(k,idetsi)
#else
    sms(k,idetsi)     = (                          &
        + aggregationRate                  * DiaSi &
        + lossN_d         * limitFacN_dia  * DiaSi &
        + grazingFlux_dia * qSiN                   &
        - reminSiT                         * DetSi &
                                                  ) * dt_b + sms(k,idetsi)
#endif
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
#if defined (__3Zoo2Det)
        + reminSiT                         * DetZ2Si &
#endif
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
#if defined (__coccos)
        - N_assim_cocco                   * CoccoC  & 
        + lossN_c       * limitFacN_cocco * CoccoN  & 
#endif
        + lossN         * limitFacN       * PhyN    & ! --> Excretion from small pythoplankton
        + lossN_d       * limitFacN_dia   * DiaN    & ! --> Excretion from diatom
        + reminN * arrFunc * O2Func * DetN          & ! --> Remineralization of detritus
        + lossN_z                   * HetN          & ! --> Excretion from zooplankton
#if defined (__3Zoo2Det)
        + reminN * arrFunc * O2Func * DetZ2N        & ! O2remin
        + lossN_z2                        * Zoo2N   &         
        + lossN_z3                        * MicZooN & ! 3Zoo
#endif
                                              )     &
        - kScavFe          * DetC   * FreeFe        & 
#if defined (__3Zoo2Det)
        - kScavFe          * DetZ2C * FreeFe        &
#endif
                                                   ) * dt_b + sms(k,ife)

!< *** Calcification ***
!< *********************

!____________________________________________________________
!< Small phytoplankton calcite

#if defined (__coccos)
    sms(k,iphycal)    = (                               &
        + calcification                                 & ! --> Calcification
        - lossC_c           * limitFacN_cocco * PhyCalc & ! --> Excretion loss 
        - phyRespRate_cocco                   * PhyCalc & ! --> Respiration      
        - calc_loss_agg                                 & ! --> Aggregation loss
        - calc_loss_gra                                 & ! --> Grazing loss
#if defined (__3Zoo2Det)
        - calc_loss_gra2                                &
        - calc_loss_gra3                                & ! 3Zoo
#endif
                                                       ) * dt_b + sms(k,iphycal)
#else
    sms(k,iphycal)    = (                      &
        + calcification                        & ! --> Calcification  
        - lossC          * limitFacN * PhyCalc & ! --> Excretion loss  
        - phyRespRate                * PhyCalc & ! --> Respiration 
        - calc_loss_agg                        & ! --> Aggregation loss 
        - calc_loss_gra                        & ! --> Grazing loss
#if defined (__3Zoo2Det)
        - calc_loss_gra2                       &
        - calc_loss_gra3                       & ! 3Zoo
#endif
                                              ) * dt_b + sms(k,iphycal)
#endif

!____________________________________________________________
! Detritus calcite
#if defined (__coccos)

#if defined (__3Zoo2Det)
    sms(k,idetcal)   = (                                &
        + lossC_c           * limitFacN_cocco * PhyCalc &
        + phyRespRate_cocco                   * PhyCalc & 
        + calc_loss_agg                                 &
        + calc_loss_gra3                                &
        - calc_loss_gra3    * calc_diss_guts            &
        - calc_diss                           * DetCalc &
                                                       ) * dt_b + sms(k,idetcal)

#else
    sms(k,idetcal)   = (                                &
        + lossC_c           * limitFacN_cocco * PhyCalc & 
        + phyRespRate_cocco                   * PhyCalc & 
        + calc_loss_agg                                 &
        + calc_loss_gra                                 &
        - calc_loss_gra     * calc_diss_guts            &
        - calc_diss                           * DetCalc &
                                                       ) * dt_b + sms(k,idetcal)

#endif

#else

#if defined (__3Zoo2Det)
    sms(k,idetcal)   = (                            &
        + lossC          * limitFacN      * PhyCalc &
        + phyRespRate                     * PhyCalc &
        + calc_loss_agg                             &
        + calc_loss_gra3                            &
        - calc_loss_gra3 * calc_diss_guts           &
        - calc_diss                       * DetCalc &
                                                   ) * dt_b + sms(k,idetcal)
#else
    sms(k,idetcal)   = (                           &
        + lossC         * limitFacN      * PhyCalc &
        + phyRespRate                    * PhyCalc &
        + calc_loss_agg                            &
        + calc_loss_gra                            &
        - calc_loss_gra * calc_diss_guts           &
        - calc_diss                      * DetCalc &
                                                  ) * dt_b + sms(k,idetcal)
#endif
#endif

!____________________________________________________________
! Oxygen

    sms(k,ioxy)   = (                       &
      + Cphot                      * phyC   &
      - phyRespRate                * phyC   &
      + Cphot_dia                  * diaC   &
      - phyRespRate_dia            * diaC   &
#if defined (__coccos)
      + Cphot_cocco                * CoccoC &
      - phyRespRate_cocco          * CoccoC &
#endif
      - rho_C1  * arrFunc * O2Func * EOC    & ! O2remin
      - hetRespFlux                         &
#if defined (__3Zoo2Det)
      - Zoo2RespFlux                        &
      - MicZooRespFlux                      & ! 3Zoo
#endif
                                           ) * redO2C * dt_b + sms(k,ioxy)  
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
!         "Abiotic" DIC_14, identical to DIC except for radioactive decay (->
!         recom_forcing)
          sms(k,idic_14) = sms(k,idic)
        end if ! ciso_organic_14
      end if   ! ciso_14
    end if     ! ciso
!-------------------------------------------------------------------------------
! Diagnostics: Averaged rates
	
	recipbiostep    = 1.d0/real(biostep)
if (Diags) then
!*** Net primary production [mmol C /(m3 * day)]
        vertNPPn(k) = vertNPPn(k)  + (   &
     	+ Cphot                   * PhyC  &
     	- PhyRespRate             * PhyC  &
     	) * recipbiostep

        vertNPPd(k) = vertNPPd(k)  + (   &
     	+ Cphot_dia               * DiaC  &
     	- PhyRespRate_dia         * DiaC  &
     	) * recipbiostep

#if defined (__coccos)
        vertNPPc(k) = vertNPPc(k)  + (   &
        + Cphot_cocco             * CoccoC  &
        - PhyRespRate_cocco       * CoccoC  & 
        ) * recipbiostep
#endif

!*** Gross primary production [mmol C /(m3 * day)]
        vertGPPn(k) = vertGPPn(k)  + (   &
     	+ Cphot                   * PhyC  &
     	) * recipbiostep

        vertGPPd(k) = vertGPPd(k)  + (   &
     	+ Cphot_dia               * DiaC  &
     	) * recipbiostep

#if defined (__coccos)
        vertGPPc(k) = vertGPPc(k)  + (   &
        + Cphot_cocco             * CoccoC  &
        ) * recipbiostep
#endif

!*** Net N-assimilation [mmol N/(m3 * day)]
        vertNNAn(k) = vertNNAn(k)  + (   &
     	+ N_assim                 * PhyC  &
     	- lossN * limitFacN       * PhyN  &
     	) * recipbiostep

        vertNNAd(k) = vertNNAd(k)  + (   &
     	+ N_assim_dia             * DiaC  &
     	- lossN * limitFacN_dia   * DiaN  &
     	) * recipbiostep

#if defined (__coccos)
        vertNNAc(k) = vertNNAc(k)  + (   &
        + N_assim_cocco           * CoccoC  &
        - lossN * limitFacN_cocco * CoccoN  &
        ) * recipbiostep
#endif

!*** Changed to chlorophyll degradation (commented out gross N-assimilation below)
        vertChldegn(k) = vertChldegn(k)  + (   &
        + KOchl  &
        ) * recipbiostep

        vertChldegd(k) = vertChldegd(k)  + (   &
        + KOchl_dia  &
        ) * recipbiostep

#if defined (__coccos)
        vertChldegc(k) = vertChldegc(k)  + (   &
        + KOchl_cocco & 
        ) * recipbiostep
#endif

!--------------------------------------------------------------------------------------------------------------------------------------

! GRAZING FLUXES
! Only for the case with detritus grazing, not without detritus grazing, because this output is probably anyway not needed as a default.
! diagnostics, combined from Onur and Cara, modified by Miriam

    if (Grazing_detritus) then

! **** First zooplankton / Mesozooplankton ****
! -----------------------------------------
!*** Total grazing of first zooplankton (with graz_eff, i.e. what reaches ZOO)
        vertgrazmeso_tot(k) = vertgrazmeso_tot(k) + (    &                       
        + grazingFlux_phy * recipQuota * grazEff         &
        + grazingFlux_Dia * recipQuota_Dia * grazEff     &
        + grazingFlux_Det * recipDet * grazEff           &
#if defined (__coccos)
        + grazingFlux_Cocco * recipQuota_Cocco * grazEff &
#endif
#if defined (__3Zoo2Det)
        + GrazingFlux_DetZ2 * recipDet2 * grazEff        &
        + grazingFlux_miczoo * recipQZoo3 * grazEff      &
#endif
        ) * recipbiostep                                 

!*** Grazing on small phytoplankton by first zooplankton (without grazeff, i.e. loss term for PHY)
        vertgrazmeso_n(k) = vertgrazmeso_n(k) + ( &  
        + grazingFlux_phy * recipQuota            &
        ) * recipbiostep                                 

!*** Grazing on diatoms by first zooplankton (without grazeff, i.e. loss term for DIA)
        vertgrazmeso_d(k) = vertgrazmeso_d(k) + ( &                                               
        + grazingFlux_dia * recipQuota_dia        & 
        ) * recipbiostep  

#if defined (__coccos)
!*** Grazing on cocclithophores by first zooplankton (without grazeff, i.e. loss term for COCCO)
        vertgrazmeso_c(k) = vertgrazmeso_c(k) + ( &
        + grazingFlux_Cocco * recipQuota_cocco    &
        ) * recipbiostep
#endif

!*** Grazing on first detritus by first zooplankton (without grazeff, i.e. loss term for DET)
        vertgrazmeso_det(k) = vertgrazmeso_det(k) + ( &
        + grazingFlux_Det * recipDet                  &
        ) * recipbiostep

#if defined (__3Zoo2Det)
!*** Grazing on microzooplankton by first zooplankton (without grazeff, i.e. loss term for MICZOO)
        vertgrazmeso_mic(k) = vertgrazmeso_mic(k) + ( &
        + grazingFlux_miczoo * recipQZoo3             &
        ) * recipbiostep

!*** Grazing on second detritus by first zooplankton (without grazeff, i.e. loss term for DET2)
        vertgrazmeso_det2(k) = vertgrazmeso_det2(k) + ( &
        + GrazingFlux_DetZ2 * recipDet2                 &
        ) * recipbiostep
#endif

        

! **** Second zooplankton / Macrozooplankton / Krill ****  
! -----------------------------------------
#if defined (__3Zoo2Det)
        
!*** Total grazing of second zooplankton (with graz_eff, i.e. what reaches ZOO)
        vertgrazmacro_tot(k) = vertgrazmacro_tot(k) + (    &
        + grazingFlux_phy2 * recipQuota * grazEff2         &
        + grazingFlux_Dia2   * recipQuota_Dia   * grazEff2 &
#if defined (__coccos)
        + grazingFlux_Cocco2 * recipQuota_Cocco * grazEff2 &
#endif
        + grazingFlux_het2   * recipQZoo        * grazEff2 &
        + grazingFlux_miczoo2* recipQZoo3       * grazEff2 &
        + grazingFlux_Det2   * recipDet         * grazEff2 &
        + grazingFlux_DetZ22 * recipDet2        * grazEff2 &
        ) * recipbiostep

!*** Grazing on small phytoplankton by second zooplankton (without grazeff, i.e. loss term for PHY)
        vertgrazmacro_n(k) = vertgrazmacro_n(k) + ( &
        + grazingFlux_phy2 * recipQuota             &
        ) * recipbiostep

!*** Grazing on diatoms by second zooplankton (without grazeff, i.e. loss term for DIA) 
        vertgrazmacro_d(k) = vertgrazmacro_d(k) + ( &
        + grazingFlux_Dia2 * recipQuota_Dia         &
        ) * recipbiostep
        
#if defined (__coccos)
!*** Grazing on cocclithophores by second zooplankton (without grazeff, i.e. loss term for COCCO)
        vertgrazmacro_c(k) = vertgrazmacro_c(k) + ( &
        + grazingFlux_Cocco2 * recipQuota_cocco     &
        ) * recipbiostep
#endif

!*** Grazing on mesozooplankton by second zooplankton (without grazeff, i.e. loss term for HET)
        vertgrazmacro_mes(k) = vertgrazmacro_mes(k) + ( &
        + grazingFlux_het2   * recipQZoo                &
        ) * recipbiostep
        
!*** Grazing on first detritus by second zooplankton (without grazeff, i.e. loss term for DET)
        vertgrazmacro_det(k) = vertgrazmacro_det(k) + ( &
        + grazingFlux_Det2 * recipDet                   &
        ) * recipbiostep

!*** Grazing on microzooplankton by second zooplankton (without grazeff, i.e. loss term for MICZOO)
        vertgrazmacro_mic(k) = vertgrazmacro_mic(k) + ( &
        + grazingFlux_miczoo2 * recipQZoo3              &
        ) * recipbiostep

!*** Grazing on second detritus by second zooplankton (without grazeff, i.e. loss term for DET2)
        vertgrazmacro_det2(k) = vertgrazmacro_det2(k) + ( &
        + GrazingFlux_DetZ22 * recipDet2                  &
        ) * recipbiostep

#endif


! **** Third zooplankton / Microzooplankton ****
! -----------------------------------------                                                                                                                  
#if defined (__3Zoo2Det)
        
!*** Total grazing of third zooplankton (with graz_eff, i.e. what reaches ZOO) 
        vertgrazmicro_tot(k) = vertgrazmicro_tot(k) + (     &
        + grazingFlux_phy3    * recipQuota       * grazEff3 &
        + grazingFlux_Dia3    * recipQuota_Dia   * grazEff3 &
#if defined (__coccos)
        + grazingFlux_Cocco3  * recipQuota_Cocco * grazEff3 &
#endif
        ) * recipbiostep

!*** Grazing on small phytoplankton by third zooplankton (without grazeff, i.e. loss term for PHY)
        vertgrazmicro_n(k) = vertgrazmicro_n(k) + ( &
        + grazingFlux_phy3 * recipQuota             &
        ) * recipbiostep

!*** Grazing on diatoms by third zooplankton (without grazeff, i.e. loss term for DIA)
        vertgrazmicro_d(k) = vertgrazmicro_d(k) + ( &
        + grazingFlux_Dia3 * recipQuota_Dia         &
        ) * recipbiostep

#if defined (__coccos)
!*** Grazing on cocclithophores by third zooplankton (without grazeff, i.e. loss term for COCCO)
        vertgrazmicro_c(k) = vertgrazmicro_c(k) + ( &
        + grazingFlux_Cocco3 * recipQuota_cocco     &
        ) * recipbiostep
#endif

#endif

     end if ! grazing_detritus

!-----------------------------------------------------------------------------------------------------------------
     
!*** zooplankton1 respiration
        vertrespmeso(k) = vertrespmeso(k) + (     &
        + HetRespFlux                             &
        ) * recipbiostep
#if defined (__3Zoo2Det)
!*** zooplankton2 respiration
        vertrespmacro(k) = vertrespmacro(k) + (   &
        + Zoo2RespFlux                            &
        ) * recipbiostep

!*** zooplankton3 respiration
        vertrespmicro(k) = vertrespmicro(k) + (  &
        + MicZooRespFlux                         &
        ) * recipbiostep
#endif
!*** calc_diss
        vertcalcdiss(k) = vertcalcdiss(k) + (     &
        + calc_diss * DetCalc                     &
        ) * recipbiostep

!***    aggregation by  small phytoplankton                                                                      
        vertaggn(k) = vertaggn(k) + (             & 
        + aggregationrate * PhyC                  &
        ) * recipbiostep

!***    aggregation by  diatoms                                                                                  
        vertaggd(k) = vertaggd(k) + (             &
        + aggregationrate * DiaC                  &
        ) * recipbiostep  

#if defined (__coccos)
!***    aggregation by coccolithophores
        vertaggc(k) = vertaggc(k) + (             &
        + aggregationrate * CoccoC                &    
        ) * recipbiostep
#endif                                

!*** excrection of DOC by phytoplankton
        vertdocexn(k) = vertdocexn(k) + (         &
        + lossC * limitFacN              * phyC   &
        ) * recipbiostep  
  
!*** excrection of DOC by diatoms
        vertdocexd(k) = vertdocexd(k) + (         &
        + lossC_d * limitFacN_dia        * DiaC   &
        ) * recipbiostep  

#if defined (__coccos)
!*** excretion of DOC by coccolithophores
        vertdocexc(k) = vertdocexc(k) + (         &
        + lossC_c * limitFacN_cocco      * CoccoC &  
        ) * recipbiostep
#endif

!*** calcification
        vertcalcif(k) = vertcalcif(k) + (         &
        + calcification                           &
        ) * recipbiostep  

! phy respiration
        vertrespn(k) = vertrespn(k) + (           &
     	+ PhyRespRate             * PhyC          &
     	) * recipbiostep

! dia respiration
        vertrespd(k) = vertrespd(k) + (           &
     	+ PhyRespRate_dia         * DiaC          &
     	) * recipbiostep

#if defined (__coccos)
! cocco resipration
        vertrespc(k) = vertrespc(k) + (           &
        + PhyRespRate_cocco       * CoccoC        &
        ) * recipbiostep

#endif
endif
  end do ! Main vertikal loop ends

!   if (mype==0 .and. my_fesom_group == 0) then !OG
!      write(*,*) ' --> Main vertikal loop ends'
!   endif

!-------------------------------------------------------------------------------
! Remineralization from the sediments into the bottom layer

  if (use_MEDUSA .and. (sedflx_num .ne. 0)) then
   if (mype==0 .and. my_fesom_group == 0) then !OG
      write(*,*) ' --> Sedimentary input of nutrients through MEDUSA'
   endif

  else ! not use_MEDUSA or sedflx_num = 0
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
  decayBenthos(4) = calc_diss_ben * LocBenthos(4)
  LocBenthos(4)      = LocBenthos(4)   - decayBenthos(4) * dt_b

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
!         Do nothing here because sms(idic_14) is defined as sms(idic) further
!         above
        end if ! ciso_organic_14
      end if   ! ciso_14
    end if ! ciso
  endif ! use_MEDUSA

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

