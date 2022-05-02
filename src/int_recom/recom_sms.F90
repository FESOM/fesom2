subroutine REcoM_sms(n,Nn,state,thick,recipthick,SurfSR,sms,Temp, Sali_depth &
        , CO2_watercolumn                                              & ! NEW MOCSY [mol/m3]
        , pH_watercolumn                                               & ! NEW MOCSY on total scale
        , pCO2_watercolumn                                             & ! NEW MOCSY [uatm]
        , HCO3_watercolumn                                             & ! NEW MOCSY [mol/m3]
        , CO3_watercolumn                                              & ! NEW DISS [mol/m3]
        , OmegaC_watercolumn                                           & ! NEW DISS calcite saturation state
        , kspc_watercolumn                                             & ! NEW DISS stoichiometric solubility product [mol^2/kg^2]
        , rhoSW_watercolumn                                            & ! NEW DISS in-situ density of seawater [kg/m3]
        , Nutlim_phy                                                   & ! NEWOUT nutrient limitation of small phytoplankton
        , Nutlim_dia                                                   & ! NEWOUT nutrient limitation of diatoms
        , Nutlim_cocco                                                 & ! NEWOUT nutrient limitation of coccolithophores
        , Tlim_arr                                                     & ! NEWOUT temperature limitation according to Arrhenius function
        , Tlim_cocco                                                   & ! NEWOUT temperature limitation of coccolithophores
        , Llim_phy                                                     & ! NEWOUT light limitation of small phytoplankton
        , Llim_dia                                                     & ! NEWOUT light limitation of diatoms
        , Llim_cocco                                                   & ! NEWOUT light limitation of coccolithophores
        , CO2lim_phy                                                   & ! NEWOUT CO2 limitation of small phytoplankton
        , CO2lim_dia                                                   & ! NEWOUT CO2 limitation of diatoms
        , CO2lim_cocco                                                 & ! NEWOUT CO2 limitation of coccolithophores
        , PR_phy                                                       & ! NEWOUT Photosynthesis rate of small phytoplankton
        , PR_dia                                                       & ! NEWOUT Photosynthesis rate of diatoms
        , PR_cocco                                                     & ! NEWOUT Photosynthesis rate of coccolithophores
        , Cal_Tlim                                                     & ! NEWOUT Temperature dependence of calcification
        , Cal_CO2lim                                                   & ! NEWOUT CO2 dependence of calcification
        , Cal_Nlim                                                     & ! NEWOUT Nitrate dependence of calcification
        , Cal_pure                                                     & ! NEWOUT PIC only dependent on PICPOCmax, CoccoC, T, N, CO2
        , Loc_slp, SinkVel,zF,PAR, Lond, Latd, mesh)

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
    use mvars                                                                ! NEW MOCSY


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
    real(kind=8),dimension(mesh%nl,4)        ,intent(in)    :: SinkVel
    Real(kind=8),dimension(mesh%nl-1),intent(inout)         :: CO2_watercolumn      ! NEW MOCSY
    Real(kind=8),dimension(mesh%nl-1),intent(inout)         :: pH_watercolumn       ! NEW MOCSY
    Real(kind=8),dimension(mesh%nl-1),intent(inout)         :: pCO2_watercolumn     ! NEW MOCSY
    Real(kind=8),dimension(mesh%nl-1),intent(inout)         :: HCO3_watercolumn     ! NEW MOCSY
    Real(kind=8),dimension(mesh%nl-1),intent(inout)         :: CO3_watercolumn      ! NEW DISS
    Real(kind=8),dimension(mesh%nl-1),intent(inout)         :: OmegaC_watercolumn   ! NEW DISS
    Real(kind=8),dimension(mesh%nl-1),intent(inout)         :: kspc_watercolumn     ! NEW DISS
    Real(kind=8),dimension(mesh%nl-1),intent(inout)         :: rhoSW_watercolumn    ! NEW DISS
    Real(kind=8),dimension(mesh%nl-1),intent(inout)         :: Nutlim_phy           ! NEWOUT
    Real(kind=8),dimension(mesh%nl-1),intent(inout)         :: Nutlim_dia           ! NEWOUT
    Real(kind=8),dimension(mesh%nl-1),intent(inout)         :: Nutlim_cocco         ! NEWOUT
    Real(kind=8),dimension(mesh%nl-1),intent(inout)         :: Tlim_arr             ! NEWOUT
    Real(kind=8),dimension(mesh%nl-1),intent(inout)         :: Tlim_cocco           ! NEWOUT
    Real(kind=8),dimension(mesh%nl-1),intent(inout)         :: Llim_phy             ! NEWOUT
    Real(kind=8),dimension(mesh%nl-1),intent(inout)         :: Llim_dia             ! NEWOUT
    Real(kind=8),dimension(mesh%nl-1),intent(inout)         :: Llim_cocco           ! NEWOUT
    Real(kind=8),dimension(mesh%nl-1),intent(inout)         :: CO2lim_phy           ! NEWOUT
    Real(kind=8),dimension(mesh%nl-1),intent(inout)         :: CO2lim_dia           ! NEWOUT
    Real(kind=8),dimension(mesh%nl-1),intent(inout)         :: CO2lim_cocco         ! NEWOUT
    Real(kind=8),dimension(mesh%nl-1),intent(inout)         :: PR_phy               ! NEWOUT
    Real(kind=8),dimension(mesh%nl-1),intent(inout)         :: PR_dia               ! NEWOUT
    Real(kind=8),dimension(mesh%nl-1),intent(inout)         :: PR_cocco             ! NEWOUT
    Real(kind=8),dimension(mesh%nl-1),intent(inout)         :: Cal_Tlim             ! NEWOUT
    Real(kind=8),dimension(mesh%nl-1),intent(inout)         :: Cal_CO2lim           ! NEWOUT
    Real(kind=8),dimension(mesh%nl-1),intent(inout)         :: Cal_Nlim             ! NEWOUT
    Real(kind=8),dimension(mesh%nl-1),intent(inout)         :: Cal_pure             ! NEWOUT

    real(kind=8),dimension(mesh%nl)          ,intent(in)    :: zF                   !< [m] Depth of fluxes
    real(kind=8),dimension(mesh%nl-1),intent(inout)         :: PAR

    real(kind=8)                                            :: net                  

    real(kind=8)                                            :: dt_d                 !< Size of time steps [day]
    real(kind=8)                                            :: dt_b                 !< Size of time steps [day]
    real(kind=8),dimension(mesh%nl-1)                       :: Sink
    real(kind=8)                                            :: dt_sink              !< Size of local time step
    real(kind=8)                                            :: Fc                   !< Flux of labile C into sediment, used for denitrification calculation [umolC/cm2/s]
    real(kind=8)                                            :: recip_hetN_plus      !< MB's addition to heterotrophic respiration
    real(kind=8)                                            :: recip_res_het        !< [day] Reciprocal of respiration by heterotrophs and mortality (loss to detritus)
    real(kind=8)                                            :: aux
    integer                                                 :: k,step,ii, idiags,n

!if (use_coccos) then
    Real(kind=8),                      intent(in)    :: Loc_slp              ! NEW MOCSY [Pa] sea-level pressure
    Real(kind=8)                                     :: Patm_depth(1)        ! NEW MOCSY
    Real(kind=8)                                     :: REcoM_T_depth(1)        ! NEW MOCSY temperature for the whole water column for mocsy minimum defined as -2
    Real(kind=8)                                     :: REcoM_S_depth(1)        ! NEW MOCSY
    Real(kind=8)                                     :: REcoM_DIC_depth(1)      ! NEW MOCSY
    Real(kind=8)                                     :: REcoM_Alk_depth(1)      ! NEW MOCSY
    Real(kind=8)                                     :: REcoM_Si_depth(1)       ! NEW MOCSY
    Real(kind=8)                                     :: REcoM_Phos_depth(1)     ! NEW MOCSY
    Real(kind=8),                      intent(in)    :: Latd(1)              ! latitude in degree
    Real(kind=8),                      intent(in)    :: Lond(1)              ! NEW MOCSY longitude in degree 
    Real(kind=8)                                     :: mocsy_step_per_day
!end if 

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
!    if (use_coccos) then
    CoccoN,  & ! NEW
    CoccoC,  & ! NEW
    CoccoChl,& ! NEW
!    endif
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
            alpha_calc_13 = 1.d0       ! during calcification / dissolution
            alpha_calc_14 = 1.d0
            alpha_dcal_13 = 1.d0
            alpha_dcal_14 = 1.d0
        endif
    endif

    sms = zero ! double precision

    tiny_N   = tiny_chl/chl2N_max      !< 0.00001/ 3.15d0   Chl2N_max [mg CHL/mmol N] Maximum CHL a : N ratio = 0.3 gCHL gN^-1
    tiny_N_d = tiny_chl/chl2N_max_d    !< 0.00001/ 4.2d0
    tiny_N_c = tiny_chl/chl2N_max_c    ! NEW
    tiny_C   = tiny_N  /NCmax          !< NCmax   = 0.2d0   [mmol N/mmol C] Maximum cell quota of nitrogen (N:C)
    tiny_C_d = tiny_N_d/NCmax_d        !< NCmax_d = 0.2d0
    tiny_C_c = tiny_N_c/NCmax_c       ! NEW
    tiny_Si  = tiny_C_d/SiCmax         !< SiCmax = 0.8d0

    recip_res_het = 1.d0/res_het       !< res_het = 0.01d0  [1/day] Respiration by heterotrophs and mortality (loss to detritus)

    Patm_depth       = Loc_slp/Pa2atm             ! NEW MOCSY convert from Pa to atm.

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
            DiaSi  = max(tiny_si,state(k,idiasi)      	+ sms(k,idiasi)) 
            DetSi  = max(tiny,state(k,idetsi) 		+ sms(k,idetsi)) 
            Si     = max(tiny,state(k,isi)    		+ sms(k,isi   ))
!if (use_coccos) then
            CoccoN = max(tiny,state(k,icocn)            + sms(k,icocn  )) ! NEW
            CoccoC = max(tiny,state(k,icocc)            + sms(k,icocc  )) ! NEW
            CoccoChl = max(tiny,state(k,icchl)          + sms(k,icchl)) ! NEW
!end if 
            Fe     = max(tiny,state(k,ife)    		+ sms(k,ife   ))
            O2     = max(tiny,state(k,ioxy)             + sms(k,ioxy  ))
            FreeFe = zero

!if (use_coccos) then
! For Mocsy
    REcoM_T_depth = max(2.d0, Temp(k))   ! NEW MOCSY minimum set to 2 degC: K1/K2 Lueker valid between 2degC-35degC and 19-43psu
    REcoM_T_depth = min(REcoM_T_depth, 40.d0)! NEW MOCSY maximum set to 40 degC: K1/K2 Lueker valid between 2degC-35degC and 19-43psu
    REcoM_S_depth = max(21.d0, Sali_depth(k))! NEW MOCSY minimum set to 21: K1/K2 Lueker valid between 2degC-35degC and 19-43psu, else causes trouble in regions with S between 19 and 21 and ice conc above 97%                      
    REcoM_S_depth    = min(REcoM_S_depth, 43.d0)  ! NEW MOCSY maximum set to 43: K1/K2 Lueker valid between 2degC-35degC and 19-43psu, else causes trouble   
    REcoM_DIC_depth  = max(tiny*1e-3,state(k,idic)*1e-3   + sms(k,idic  )*1e-3)      ! NEW MOCSY
    REcoM_Alk_depth  = max(tiny*1e-3,state(k,ialk)*1e-3   + sms(k,ialk  )*1e-3)      ! NEW MOCSY
    REcoM_Si_depth   = max(tiny*1e-3,state(k,isi)*1e-3    + sms(k,isi   )*1e-3)      ! NEW MOCSY
    REcoM_Phos_depth = max(tiny*1e-3,state(k,idin)*1e-3   + sms(k,idin  )*1e-3) /16  ! NEW MOCSY convert N to P with Redfield
!endif
            PhyCalc= max(tiny,state(k,iphycal)		+ sms(k,iphycal))
            DetCalc= max(tiny,state(k,idetcal)		+ sms(k,idetcal))

            !calc_diss      = calc_diss_rate * SinkVel(k,ivdet) /20.d0 ! Dissolution rate of CaCO3 scaled by the sinking velocity at the current depth 0.005714   !20.d0/3500.d0 ! NEW DISS: Moved behind Mocsy call
            !calc_diss2     = calc_diss_rate2  ! Dissolution rate of CaCO3 for seczoo ! NEW DISS: Moved behind Mocsy call

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
!if (use_coccos) then
    quota_cocco    = CoccoN / CoccoC           ! NEW
    recipQuota_cocco = real(one)/quota_cocco   ! NEW
    Chl2C_cocco    = CoccoChl / CoccoC         ! NEW
    Chl2N_cocco    = CoccoChl / CoccoN         ! NEW
    CHL2C_plast_cocco = Chl2C_cocco * (quota_cocco/(quota_cocco - NCmin_c)) ! NEW
!end if 
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
    CoccoTFunc     = max(0.1419d0 * Temp(k)**0.8151d0,tiny)        ! NEW (function from Fielding 2013; is based on observational GR, but range fits best to our arrFunc; they use T in degree celsius, so C2K is not needed)
    if (REcoM_Second_Zoo) then 
        arrFuncZoo2 = exp(t1_zoo2/t2_zoo2 - t1_zoo2*rTloc)/(1 + exp(t3_zoo2/t4_zoo2 - t3_zoo2*rTloc))
    endif

!< Silicate temperature dependence 
!    reminSiT = min(1.32e16 * exp(-11200.d0 * rTloc),reminSi) !! arrFunc control, reminSi=0.02d0 ! Kamatani (1982)
    reminSiT = max(0.023d0 * 2.6d0**((Temp(k)-10.)/10.),reminSi)
!    reminSiT = reminSi

    Tlim_arr(k)    = arrFunc                                       ! NEWOUT
    Tlim_cocco(k)  = CoccoTFunc                                    ! NEWOUT

!-------------------------------------------------------------------------------
! Light                                                                                                   ! NEW MOCSY
!-------------------------------------------------------------------------------
! Has to be calculated here already to use the 1%PAR depth.

    if (k==1) then
       PARave    = max(tiny,SurfSR)
       PAR(k)    = PARave
       chl_upper = (PhyChl + DiaChl + CoccoChl)                       ! NEW (edited term)
                                                                                                                                  
    else
       chl_lower = PhyChl + DiaChl + CoccoChl                         ! NEW (edited term)                            
                                                                                                                                    
       Chlave    = (chl_upper+chl_lower)*0.5

       kappa          =  k_w + a_chl * (Chlave)
       kappastar      =  kappa / cosAI(n)
       kdzLower       =  kdzUpper + kappastar * thick(k-1)
       Lowerlight     =  SurfSR * exp(-kdzLower)
       Lowerlight     =  max(tiny,Lowerlight) 
       PARave         =  Lowerlight
       PAR(k)         =  PARave
       chl_upper      =  chl_lower
       kdzUpper       =  kdzLower
    end if

!-------------------------------------------------------------------------------
! Depth component of Mocsy (see http://ocmip5.ipsl.jussieu.fr/mocsy/pyth.html)                                                                                                      ! NEW MOCSY
!-------------------------------------------------------------------------------

! Calculate the carbonate system for the very first time step of the first year of the run
    !if (mocsy_restart==.false. .and. recom_istep==1) then    ! r_restart is defined in gen_modules_clock in fesom_cpl.
    if (mstep==1) then
!       call vars_sprac(ph_depth, pco2_depth, fco2_depth, co2_depth, hco3_depth, co3_depth, OmegaA_depth, OmegaC_depth, kspc_depth, BetaD_depth, rhoSW_depth, p_depth, tempis_depth, & ! NEW MOCSY
!        REcoM_T_depth, REcoM_S_depth, REcoM_Alk_depth, REcoM_DIC_depth, REcoM_Si_depth, REcoM_Phos_depth, Patm_depth, zF(k), Latd, Nmocsy, & ! NEW MOCSY
!            optCON='mol/m3', optT='Tpot   ', optP='m ', optB='u74', optK1K2='l  ', optKf='dg', optGAS='Pinsitu', optS='Sprc') ! NEW MOCSY
       CO2_watercolumn(k)    = co2_depth(1)
       pH_watercolumn(k)     = ph_depth(1)
       pCO2_watercolumn(k)   = pco2_depth(1)
       HCO3_watercolumn(k)   = hco3_depth(1)
       CO3_watercolumn(k)    = co3_depth(1)    ! NEW DISS
       OmegaC_watercolumn(k) = OmegaC_depth(1) ! NEW DISS
       kspc_watercolumn(k)   = kspc_depth(1)   ! NEW DISS
       rhoSW_watercolumn(k)  = rhoSW_depth(1)  ! NEW DISS
    endif 

! Calculate carbonate system every 7 days for depths < 1%PAR, and every 30 days for the depths below.
    mocsy_step_per_day = dt ! define this in namelist.recom OG
    logfile_outfreq_7  = mocsy_step_per_day*7                    ! mocsy_step_per_day seems to be defined in the fesom namelist.
    logfile_outfreq_30 = mocsy_step_per_day*30

    if (PARave > 0.01*SurfSR .and. mod(mstep,logfile_outfreq_7)==0) then
!       call vars_sprac(ph_depth, pco2_depth, fco2_depth, co2_depth, hco3_depth, co3_depth, OmegaA_depth, OmegaC_depth, kspc_depth, BetaD_depth, rhoSW_depth, p_depth, tempis_depth, & ! NEW MOCSY
!                            REcoM_T_depth, REcoM_S_depth, REcoM_Alk_depth, REcoM_DIC_depth, REcoM_Si_depth, REcoM_Phos_depth, Patm_depth, zF(k), Latd, Nmocsy,                      & ! NEW MOCSY
!                            optCON='mol/m3', optT='Tpot   ', optP='m ', optB='u74', optK1K2='l  ', optKf='dg', optGAS='Pinsitu', optS='Sprc')                                         ! NEW MOCSY
       CO2_watercolumn(k)    = co2_depth(1)
       pH_watercolumn(k)     = ph_depth(1)
       pCO2_watercolumn(k)   = pco2_depth(1)
       HCO3_watercolumn(k)   = hco3_depth(1)
       CO3_watercolumn(k)    = co3_depth(1)    ! NEW DISS
       OmegaC_watercolumn(k) = OmegaC_depth(1) ! NEW DISS
       kspc_watercolumn(k)   = kspc_depth(1)   ! NEW DISS
       rhoSW_watercolumn(k)  = rhoSW_depth(1)  ! NEW DISS
    elseif (PARave < 0.01*SurfSR .and. mod(mstep,logfile_outfreq_30)==0) then
!       call vars_sprac(ph_depth, pco2_depth, fco2_depth, co2_depth, hco3_depth, co3_depth, OmegaA_depth, OmegaC_depth, kspc_depth, BetaD_depth, rhoSW_depth, p_depth, tempis_depth, & ! NEW MOCSY
!                           REcoM_T_depth, REcoM_S_depth, REcoM_Alk_depth, REcoM_DIC_depth, REcoM_Si_depth, REcoM_Phos_depth, Patm_depth, zF(k), Latd, Nmocsy,                       & ! NEW MOCSY
!                           optCON='mol/m3', optT='Tpot   ', optP='m ', optB='u74', optK1K2='l  ', optKf='dg', optGAS='Pinsitu', optS='Sprc')                                          ! NEW MOCSY
       CO2_watercolumn(k)    = co2_depth(1)
       pH_watercolumn(k)     = ph_depth(1)
       pCO2_watercolumn(k)   = pco2_depth(1)
       HCO3_watercolumn(k)   = hco3_depth(1)
       CO3_watercolumn(k)    = co3_depth(1)    ! NEW DISS
       OmegaC_watercolumn(k) = OmegaC_depth(1) ! NEW DISS
       kspc_watercolumn(k)   = kspc_depth(1)   ! NEW DISS
       rhoSW_watercolumn(k)  = rhoSW_depth(1)  ! NEW DISS
    endif
    
!-------------------------------------------------------------------------------
! CO2 dependence of rates ! NEW CO2
!-------------------------------------------------------------------------------
! Convert pH to proton concentration
   h_depth(1) = 10.**(-ph_depth(1))
   !if(recom_mype==1) write(*,*), 'pH =',pH_depth(1)
   !if(recom_mype==1) write(*,*), 'H converted from pH =',h_depth(1)
   !if(recom_mype==1) write(*,*), 'pCO2 =',pco2_depth(1)

! Conversion factor between [mol/m3] (model) and [umol/kg] (function): (1000 * 1000) / 1024
   Cunits = 976.5625
   ! Conversion not needed for [H], because in model and function derived from pH and therefore in [mol/L]

! Multiple driver interaction: calculation of d in the CO2 function  ! NEW inter
   d_CT_CL_phy = 1.023e+07 + m_CT_phy * (Temp(k) - T0_CT) + m_CL_phy * (PARave - L0_CL)
   d_CT_CL_dia = 2.640e+06 + m_CT_dia * (Temp(k) - T0_CT) + m_CL_dia * (PARave - L0_CL)
   d_CT_CL_coc = 9.450e+06 + m_CT_coc * (Temp(k) - T0_CT) + m_CL_coc * (PARave - L0_CL)

! Small phytoplankton
   if (inter_CT_CL) then ! NEW inter
    PhyCO2 = 1.162e+00 * HCO3_watercolumn(k) * Cunits / (4.888e+01 + HCO3_watercolumn(k) * Cunits) - exp(-2.255e-01 * CO2_watercolumn(k) * Cunits) - d_CT_CL_phy * 10.**(-pH_watercolumn(k))
   else
    PhyCO2 = 1.162e+00 * HCO3_watercolumn(k) * Cunits / (4.888e+01 + HCO3_watercolumn(k) * Cunits) - exp(-2.255e-01 * CO2_watercolumn(k) * Cunits) - 1.023e+07 * 10.**(-pH_watercolumn(k))
    !PhyCO2 = 1.278e+00 * HCO3_watercolumn(k) * Cunits / (5.377e+01 + HCO3_watercolumn(k) * Cunits) - exp(-2.481e-01 * CO2_watercolumn(k) * Cunits) - 1.125e+07 * 10.**(-pH_watercolumn(k)) ! lowSens
   endif
    !PhyCO2 = min(PhyCO2,real(one))
    PhyCO2 = min(PhyCO2,3.0)      ! April 2022: limitation to 3
    !if(recom_mype==1) write(*,*), 'PhyCO2 =',PhyCO2 
    CO2lim_phy(k)   = PhyCO2      ! NEWOUT

! Diatoms
   if (inter_CT_CL) then ! NEW inter
    DiaCO2 = 1.040e+00 * HCO3_watercolumn(k) * Cunits / (2.890e+01 + HCO3_watercolumn(k) * Cunits) - exp(-8.778e-01 * CO2_watercolumn(k) * Cunits) - d_CT_CL_dia * 10.**(-pH_watercolumn(k))
   else
    DiaCO2 = 1.040e+00 * HCO3_watercolumn(k) * Cunits / (2.890e+01 + HCO3_watercolumn(k) * Cunits) - exp(-8.778e-01 * CO2_watercolumn(k) * Cunits) - 2.640e+06 * 10.**(-pH_watercolumn(k))
    !DiaCO2 = 1.144e+00 * HCO3_watercolumn(k) * Cunits / (3.179e+01 + HCO3_watercolumn(k) * Cunits) - exp(-9.656e-01 * CO2_watercolumn(k) * Cunits) - 2.904e+06 * 10.**(-pH_watercolumn(k)) ! lowSens
   endif
    !DiaCO2 = min(DiaCO2,real(one))
    DiaCO2 = min(DiaCO2,3.0)      ! April 2022: limitation to 3
    !if(recom_mype==1) write(*,*), 'DiaCO2 =',DiaCO2
    CO2lim_dia(k)   = DiaCO2      ! NEWOUT

! Coccolithophores
   if (inter_CT_CL) then ! NEW inter
    CoccoCO2 = 1.109e+00 * HCO3_watercolumn(k) * Cunits / (3.767e+01 + HCO3_watercolumn(k) * Cunits) - exp(-3.912e-01 * CO2_watercolumn(k) * Cunits) - d_CT_CL_coc * 10.**(-pH_watercolumn(k))
   else
    CoccoCO2 = 1.109e+00 * HCO3_watercolumn(k) * Cunits / (3.767e+01 + HCO3_watercolumn(k) * Cunits) - exp(-3.912e-01 * CO2_watercolumn(k) * Cunits) - 9.450e+06 * 10.**(-pH_watercolumn(k))
    !CoccoCO2 = 1.220e+00 * HCO3_watercolumn(k) * Cunits / (4.144e+01 + HCO3_watercolumn(k) * Cunits) - exp(-4.303e-01 * CO2_watercolumn(k) * Cunits) - 1.040e+07 * 10.**(-pH_watercolumn(k)) ! lowSens
   endif
    !CoccoCO2 = min(CoccoCO2,real(one))
    CoccoCO2 = min(CoccoCO2,3.0)  ! April 2022: limitation to 3
    !if(recom_mype==1) write(*,*), 'CoccoCO2 =',CoccoCO2
    CO2lim_cocco(k) = CoccoCO2    ! NEWOUT

    !if(recom_mype==1) write(*,*), 'CoccoCO2 =',CoccoCO2

!------------------------------------------------------------------------------
! Calcite dissolution dependent on OmegaC !NEW DISS
!------------------------------------------------------------------------------

   if (OmegaC_diss) then
    Ca                 = (0.02128d0/40.078d0) * Sali_depth(k)/1.80655d0 ! Calcium ion concentration [mol/kg], function from varsolver.f90
    CO3_sat            = (kspc_watercolumn(k) / Ca) * rhoSW_watercolumn(k) ! Saturated carbonate ion concentration, converted to [mol/m3]
    calc_diss          = calc_diss_omegac * max(zero,(1-(CO3_watercolumn(k)/CO3_sat)))**(calc_diss_exp) ! Dissolution rate scaled by carbonate ratio, after Aumont et al. 2015
    !if(recom_mype==1) write(*,*), 'Max without exponent = ',max(zero,(1-(CO3_watercolumn(k)/CO3_sat)))
    !if(recom_mype==1) write(*,*), 'Max with exponent as integer = ',max(zero,(1-(CO3_watercolumn(k)/CO3_sat)))**(1)
    !if(recom_mype==1) write(*,*), 'Max with exponent variable = ',max(zero,(1-(CO3_watercolumn(k)/CO3_sat)))**(calc_diss_exp)
    !if(recom_mype==1) write(*,*), 'Max with less parentheses = ',max(zero,1-(CO3_watercolumn(k)/CO3_sat))**(calc_diss_exp)
    !CO3_watercolumn(k) = CO3_watercolumn(k) + calc_diss * DetCalc * dt  ! Add newly dissolved carbonate iones to pool inbetween mocsy timesteps to account for them in diss function
    !calc_diss      = calc_diss_omegac * max(zero,1-OmegaC_watercolumn(k))**(1) ! Dissolution rate scaled by the calcite saturation state OmegaC
   else
    calc_diss      = calc_diss_rate * SinkVel(k,ivdet) /20.d0 ! Dissolution rate of CaCO3 scaled by the sinking velocity at the current depth
   endif
    calc_diss2     = calc_diss_rate2  ! Dissolution rate of CaCO3 scaled by the sinking velocity at the current depth  seczoo

   calc_diss_ben   = calc_diss_rate * SinkVel(k,ivdet) /20.0 ! NEW DISS added the variable calc_diss_ben to keep the calcite dissolution in the benthos with the old formulation




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
    feLimitFac  	= Fe/(k_Fe + Fe)                       ! Use Michaelis–Menten kinetics
    qlimitFac   	= min(qlimitFac,feLimitFac)            ! Liebig law of the minimum

    Nutlim_phy(k) = qlimitFac                                     ! NEWOUT

    pMax          	= P_cm * qlimitFac * arrFunc           ! Maximum value of C-specific rate of photosynthesis
    
!_______________________________________________________________________
!< Diatoms
    qlimitFac     	= recom_limiter(NMinSlope,NCmin_d,quota_dia)
    qlimitFacTmp  	= recom_limiter(SiMinSlope,SiCmin,qSiC)
    qlimitFac     	= min(qLimitFac,qlimitFacTmp)
    feLimitFac  	= Fe/(k_Fe_d + Fe)
    qlimitFac   	= min(qlimitFac,feLimitFac)

    Nutlim_dia(k) = qlimitFac                                     ! NEWOUT

    pMax_dia      	= P_cm_d * qlimitFac * arrFunc

!_______________________________________________________________________
!< Coccolithophores (NEW!!!)

    qlimitFac     = recom_limiter(NMinSlope,NCmin_c,quota_cocco)  ! NEW
      feLimitFac  = Fe/(k_Fe_c + Fe)                              ! NEW
      qlimitFac   = min(qlimitFac,feLimitFac)                     ! NEW
    Nutlim_cocco(k)=qlimitFac                                     ! NEWOUT

    pMax_cocco    = P_cm_c * qlimitFac * CoccoTFunc               ! NEW (and here also the T dependency is changed)

!_______________________________________________________________________
!< Light
! Exponential diminition of downward irradiance
! Light ! NEW MOCSY (moved this part further up)
!    if (k==1) then  ! ulevels_nod2D(n)==1
!        PARave = max(tiny,SurfSR)
!        PAR(k) = PARave
!        chl_upper = (PhyChl + DiaChl)
!    else    
!        chl_lower = PhyChl + DiaChl
!        Chlave    = (chl_upper+chl_lower)*0.5d0

!        kappar         =  k_w + a_chl * (Chlave)
!        kappastar      =  kappar / cosAI(n)
!        kdzLower       =  kdzUpper + kappastar * thick(k-1)
!        Lowerlight     =  SurfSR * exp(-kdzLower)
!        Lowerlight     =  max(tiny,Lowerlight)
!        PARave         =  Lowerlight
!        PAR(k)         =  PARave
!        chl_upper      =  chl_lower
!        kdzUpper       =  kdzLower
!    end if

!-------------------------------------------------------------------------------
!< Small phytoplankton photosynthesis rate
    if ( pMax .lt. tiny .OR. PARave /= PARave                  &
         .OR. CHL2C /= CHL2C) then
        Cphot       = zero
      Llim_phy(k) = zero                                               ! NEWOUT
    else
     if (CO2lim) then
        Cphot     = pMax*( real(one) &
                   - exp(-alfa * Chl2C * PARave / pMax)) * PhyCO2      ! NEW CO2 added the CO2 dependence
      else
        Cphot       = pMax*( real(one) &
                     - exp(-alfa * Chl2C * PARave / pMax))
      endif
      Llim_phy(k) = real(one) - exp(-alfa * Chl2C * PARave / pMax)     ! NEWOUT
    end if
    if ( Cphot .lt. tiny) Cphot = zero

    PR_phy(k)     = Cphot                                              ! NEWOUT
    
!------------------------------------------------------------------------------
!< Diatom photosynthesis rate
    if ( pMax_dia .lt. tiny .OR. PARave /= PARave               &
         .OR. CHL2C_dia /= CHL2C_dia) then
        Cphot_dia   = zero
      Llim_dia(k) = zero                                               ! NEWOUT
    else
      if (CO2lim) then
        Cphot_dia = pMax_dia * (real(one) &
           	- exp( -alfa_d * Chl2C_dia * PARave / pMax_dia)) * DiaCO2 ! NEW CO2 added the CO2 dependence
      else
        Cphot_dia   = pMax_dia * (real(one) &
                     - exp(-alfa_d * Chl2C_dia * PARave / pMax_dia))
      endif
        Llim_dia(k) = real(one)-exp(-alfa_d*Chl2C_dia*PARave/pMax_dia)   ! NEWOUT
    end if
    if (Cphot_dia .lt. tiny) Cphot_dia = zero

    PR_dia(k)     = Cphot_dia                                          ! NEWOUT

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
! Coccolithophore photosynthesis rate (NEW!!!)

    if ( pMax_cocco .lt. tiny .OR. Parave /= Parave             &      ! NEW
         .OR. CHL2C_cocco /= CHL2C_cocco) then                         ! NEW
       Cphot_cocco = zero                                              ! NEW
       Llim_cocco(k)=zero                                              ! NEWOUT
    else                                                               ! NEW
       if (CO2lim) then
         Cphot_cocco = pMax_cocco * (real(one) &                       ! NEW
                 - exp( -alfa_c * Chl2C_cocco * PARave / pMax_cocco)) * CoccoCO2  ! NEW NEW CO2 added the CO2 dependence
       else
         Cphot_cocco = pMax_cocco * (real(one) &                       ! NEW
                 - exp( -alfa_c * Chl2C_cocco * PARave / pMax_cocco))  ! NEW
       endif
       Llim_cocco(k)=real(one)-exp(-alfa_c*Chl2C_cocco*PARave/pMax_cocco)! NEWOUT
    end if                                                             ! NEW
    if (Cphot_cocco .lt. tiny) Cphot_cocco = zero                      ! NEW

    PR_cocco(k)    = Cphot_cocco                                       ! NEW

!-------------------------------------------------------------------------------- 
!< chlorophyll degradation
    KOchl = deg_Chl
    KOchl_dia = deg_Chl_d
    KOchl_cocco = deg_Chl_c                                            ! NEW
        
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

!< Diatoms Chla loss                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     
        if (pMax_dia .lt. tiny .OR. PARave /= PARave             &
            .OR. CHL2C_plast_dia /= CHL2C_plast_dia) then
            KOchl_dia = deg_Chl_d*0.1d0
        else
            KOchl_dia = deg_Chl_d * ( 1 -                         &
               exp( -alfa_d * CHL2C_plast_dia * PARave / pMax_dia ))
            KOchl_dia = max((deg_Chl_d*0.1d0), KOchl_dia)
        end if

!    Coccolithophores chla loss (NEW!!!)
        if (pMax_cocco .lt. tiny .OR. PARave /= Parave           &     ! NEW
                 .OR. CHL2C_plast_cocco /= CHL2C_plast_cocco) then     ! NEW
           KOchl_cocco = deg_Chl_c*0.1d0                               ! NEW
        else                                                           ! NEW
           KOchl_cocco = deg_Chl_c * ( 1 -                       &     ! NEW
                exp( -alfa_c * CHL2C_plast_cocco * PARave / pMax_cocco )) ! NEW
           KOchl_cocco = max((deg_Chl_c*0.1d0), KOchl_cocco)           ! NEW
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
        if (KOchl_cocco /= KOchl_cocco) then                               ! NEW
           print*,' KOchl_cocco is ', KOchl_cocco                          ! NEW
           print*,' deg_Chl_c is ', deg_Chl_c                              ! NEW
           print*,' alfa_c is ', alfa_c                                    ! NEW
           print*,' CHL2C_c is ', CHL2C_plast_cocco                        ! NEW
           print*,' PARave is ', PARave                                    ! NEW
           print*,' pMax_c is ', pMax_cocco                                ! NEW
           stop                                                            ! NEW
        end if                                                             ! NEW

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
                      * limitFacN * (DIN/(DIN + k_din))           ! Michaelis–Menten kinetics

    V_cm           = V_cm_fact_d
    limitFacN_dia  = recom_limiter(NMaxSlope,quota_dia,NCmax_d)
    N_assim_dia    = V_cm * pMax_dia * NCUptakeRatio_d &
                      * limitFacN_dia * DIN/(DIN + k_din_d)

    V_cm           = V_cm_fact_c                                       ! NEW
    limitFacN_cocco= recom_limiter(NMaxSlope,quota_cocco,NCmax_c)      ! NEW
    N_assim_cocco  = V_cm * pMax_cocco * NCUptakeRatio_c &             ! NEW
                      * limitFacN_cocco * DIN/(DIN + k_din_c)          ! NEW

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
            * Chl2N_max_d * min(real(one),                &
            Cphot_dia /(alfa_d * Chl2C_dia * PARave))
    end if
    ChlSynth_cocco    = zero                                           ! NEW
    if (PARave .ge. tiny .AND. PARave .eq. PARave) then                ! NEW
       ChlSynth_cocco = N_assim_cocco                    &             ! NEW
          * Chl2N_max_c * min(real(one),                 &             ! NEW
            Cphot_cocco /(alfa_c * Chl2C_cocco * PARave))              ! NEW
    end if                                                             ! NEW

!-------------------------------------------------------------------------------
!< Phytoplankton respiraion rate
!< res_phy: Maintenance respiration rate constant [day−1 ]
!< biosynth: The cost of biosynthesis of N [mmol C mmol N−1 ]
    phyRespRate     = res_phy * limitFacN + biosynth * N_assim
    phyRespRate_dia = res_phy_d * limitFacN_dia +        &
        biosynth * N_assim_dia + biosynthSi * Si_assim
    phyRespRate_cocco = res_phy_c * limitFacN_cocco +        &         ! NEW
            biosynth * N_assim_cocco                                   ! NEW
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

!                varpzPhy    = pzPhy   * PhyN /(pzPhy*PhyN + pzDia*DiaN + PzDet*DetN + pzDetZ2*DetZ2N)
!                varpzDia    = pzDia   * DiaN /(pzPhy*PhyN + pzDia*DiaN + PzDet*DetN + pzDetZ2*DetZ2N)
!                varpzDet    = pzDet   * DetN /(pzPhy*PhyN + pzDia*DiaN + PzDet*DetN + pzDetZ2*DetZ2N)
!                varpzDetZ2  = pzDetZ2 * DetN /(pzPhy*PhyN + pzDia*DiaN + PzDet*DetN + pzDetZ2*DetZ2N)
                aux         = pzPhy*PhyN + pzDia*DiaN + pzCocco*CoccoN + PzDet*DetN + pzDetZ2*DetZ2N ! NEW added Cocco
                varpzPhy    = pzPhy   * PhyN /aux
                varpzDia    = pzDia   * DiaN /aux
                varpzCocco  = pzCocco * CoccoN /aux ! NEW
                varpzDet    = pzDet   * DetN /aux
                varpzDetZ2  = pzDetZ2 * DetN /aux
            else
                DiaNsq      = DiaN**2
                PhyNsq      = PhyN**2
                CoccoNsq      = CoccoN**2                                     ! NEW
                DetNsq      = DetN**2
                DetZ2Nsq    = DetZ2N**2
!                PhyNsq      = PhyN * PhyN
                varpzPhy    = pzPhy * PhyNsq /(sPhyNsq + PhyNsq)
!                DiaNsq      = DiaN * DiaN
                varpzDia    = pzDia * DiaNsq /(sDiaNsq + DiaNsq)
!                CoccoNsq      = CoccoN * CoccoN                              ! NEW
                varpzCocco    = pzCocco * CoccoNsq /(sCoccoNsq + CoccoNsq)    ! NEW
!                DetNsq      = DetN * DetN
                varpzDet    = pzDet * DetNsq /(sDetNsq + DetNsq)
!                DetZ2Nsq    = DetZ2N * DetZ2N
                varpzDetZ2  = pzDetZ2 * DetZ2Nsq /(sDetZ2Nsq + DetZ2Nsq)

            endif

            fDiaN    = varpzDia   * DiaN
            fPhyN    = varpzPhy   * PhyN
            fCoccoN  = varpzCocco * CoccoN                            ! NEW
            fDetN    = varpzDet   * DetN
            fDetZ2N  = varpzDetZ2 * DetZ2N

        else

            if (Graz_pref_new) then
                varpzPhy      = pzPhy * PhyN /(pzPhy*PhyN + pzDia*DiaN + pzCocco*CoccoN)        ! NEW added Coccos
                varpzDia      = pzDia * DiaN /(pzPhy*PhyN + pzDia*DiaN + pzCocco*CoccoN)        ! NEW added Coccos
                varpzCocco    = pzCocco * CoccoN /(pzPhy*PhyN + pzDia*DiaN + pzCocco*CoccoN)    ! NEW
            else
                DiaNsq        = DiaN  * DiaN
                varpzDia      = pzDia * DiaNsq /(sDiaNsq + DiaNsq)
                PhyNsq        = PhyN  * PhyN
                varpzPhy      = pzPhy * PhyNsq /(sPhyNsq + PhyNsq)
                CoccoNsq      = CoccoN * CoccoN                                ! NEW
                varpzCocco    = pzCocco * CoccoNsq /(sCoccoNsq + CoccoNsq)     ! NEW
            end if
            fDiaN         = varpzDia * DiaN
            fPhyN         = varpzPhy * PhyN
            fCoccoN       = varpzCocco * CoccoN                            ! NEW
        end if
    else
        fDiaN         = pzDia * DiaN
        fPhyN         = pzPhy * PhyN
        fCoccoN       = pzCocco * CoccoN                                 ! NEW
        if (Grazing_detritus) then
            fDetN        = pzDet   * DetN
            fDetZ2N      = pzDetZ2 * DetZ2N
        end if
    end if

    if (Grazing_detritus) then
       food              = fPhyN + fDiaN + fCoccoN + fDetN + fDetZ2N      ! NEW added Coccos
       foodsq            = food * food
       grazingFlux       = (Graz_max * foodsq)/(epsilonr + foodsq) * HetN * arrFunc
       grazingFlux_phy   = grazingFlux * fphyN / food
       grazingFlux_Dia   = grazingFlux * fDiaN / food
       grazingFlux_Cocco = grazingFlux * fCoccoN / food                 ! NEW
       grazingFlux_Det   = grazingFlux * fDetN / food
       grazingFlux_DetZ2 = grazingFlux * fDetZ2N / food
    else
       food              = fPhyN + fDiaN + fCoccoN                        ! NEW added Coccos
       foodsq            = food * food
       grazingFlux       = (Graz_max * foodsq)/(epsilonr + foodsq) * HetN * arrFunc
       grazingFlux_phy   = grazingFlux * fphyN / food
       grazingFlux_Dia   = grazingFlux * fDiaN / food
       grazingFlux_Cocco = grazingFlux * fCoccoN / food                 ! NEW
    endif

    if (REcoM_Second_Zoo) then
!-------------------------------------------------------------------------------
!< Second Zooplankton grazing on small phytoplankton, diatoms and heterotrophs
!< At the moment there is no preference for one or the other food. Change this!
        
    if (REcoM_Grazing_Variable_Preference) then
       if (Grazing_detritus) then
          if (Graz_pref_new) then
!             varpzDia2      = pzDia2   * DiaN   /(pzPhy2 * PhyN + PzDia2 * DiaN + pzHet * HetN + pzDet2 * DetN + pzDetZ22 * DetZ2N)
!             varpzPhy2      = pzPhy2   * PhyN   /(pzPhy2 * PhyN + PzDia2 * DiaN + pzHet * HetN + pzDet2 * DetN + pzDetZ22 * DetZ2N)
!             varpzHet       = pzHet    * HetN   /(pzPhy2 * PhyN + PzDia2 * DiaN + pzHet * HetN + pzDet2 * DetN + pzDetZ22 * DetZ2N)  
!             varpzDet2      = pzDet2   * DetN   /(pzPhy2 * PhyN + PzDia2 * DiaN + pzHet * HetN + pzDet2 * DetN + pzDetZ22 * DetZ2N)
!             varpzDetZ22    = pzDetZ22 * DetZ2N /(pzPhy2 * PhyN + PzDia2 * DiaN + pzHet * HetN + pzDet2 * DetN + pzDetZ22 * DetZ2N)
             aux            = pzPhy2 * PhyN + PzDia2 * DiaN + pzCocco2 * CoccoN + pzHet * HetN + pzDet2 * DetN + pzDetZ22 * DetZ2N       ! NEW added Coccos
             varpzDia2      = pzDia2   * DiaN   /aux
             varpzPhy2      = pzPhy2   * PhyN   /aux
             varpzCocco2    = pzCocco2 * CoccoN /aux      ! NEW added Coccos
             varpzHet       = pzHet    * HetN   /aux  
             varpzDet2      = pzDet2   * DetN   /aux
             varpzDetZ22    = pzDetZ22 * DetZ2N /aux
          else
             DiaNsq2        = DiaN * DiaN
             varpzDia2      = pzDia2 * DiaNsq2 /(sDiaNsq2 + DiaNsq2)
             fDiaN2         = varpzDia2 * DiaN
             PhyNsq2        = PhyN * PhyN
             varpzPhy2      = pzPhy2 * PhyNsq2 /(sPhyNsq2 + PhyNsq2)
             fPhyN2         = varpzPhy2 * PhyN
             varpzCocco2    = pzCocco2 * CoccoNsq2 /(sCoccoNsq2 + CoccoNsq2)    ! NEW
             fCoccoN2       = varpzCocco2 * CoccoN                              ! NEW
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
          fCoccoN2       = varpzCocco2 * CoccoN                              ! NEW
          fHetN          = varpzHet * HetN
          fDetN2         = varpzDet2 * DetN
          fDetZ2N2       = varpzDetZ22 * DetZ2N
       else
          if (Graz_pref_new) then
!             varpzDia2      = pzDia2 * DiaN /(pzPhy2 * PhyN + PzDia2 * DiaN + pzHet * HetN)
!             varpzPhy2      = pzPhy2 * PhyN /(pzPhy2 * PhyN + PzDia2 * DiaN + pzHet * HetN)
!             varpzHet       = pzHet * HetN /(pzPhy2 * PhyN + PzDia2 * DiaN + pzHet * HetN)
             aux            = pzPhy2 * PhyN + PzDia2 * DiaN + pzCocco2 * CoccoN + pzHet * HetN         ! NEW added Coccos
             varpzDia2      = pzDia2 * DiaN /aux
             varpzPhy2      = pzPhy2 * PhyN /aux
             varpzCocco2    = pzCocco2 * CoccoN /aux     ! NEW 
             varpzHet       = pzHet  * HetN /aux

          else
             DiaNsq2        = DiaN * DiaN
             varpzDia2      = pzDia2 * DiaNsq2 /(sDiaNsq2 + DiaNsq2)
             fDiaN2         = varpzDia2 * DiaN
             PhyNsq2        = PhyN * PhyN
             varpzPhy2      = pzPhy2 * PhyNsq2 /(sPhyNsq2 + PhyNsq2)
             fPhyN2         = varpzPhy2 * PhyN
             CoccoNsq2      = CoccoN * CoccoN                                   ! NEW
             varpzCocco2    = pzCocco2 * CoccoNsq2 /(sCoccoNsq2 + CoccoNsq2)    ! NEW
             fCoccoN2       = varpzCocco2 * CoccoN                              ! NEW
             HetNsq         = HetN * HetN
             varpzHet       = pzHet * HetNsq /(sHetNsq + HetNsq)
          end if
          fDiaN2         = varpzDia2 * DiaN
          fPhyN2         = varpzPhy2 * PhyN
          fCoccoN2       = varpzCocco2 * CoccoN                                 ! NEW
          fHetN          = varpzHet * HetN
       end if
    else
       fDiaN2         = pzDia2 * DiaN
       fPhyN2         = pzPhy2 * PhyN
       fCoccoN2       = pzCocco2 * CoccoN                                       ! NEW
       fHetN          = pzHet * HetN
       if (Grazing_detritus) then
          fDetN2         = pzDet2 * DetN
          fDetZ2N2       = pzDetZ22 * DetZ2N
       end if
    end if

    if (Grazing_detritus) then
       food2             = fPhyN2 + fDiaN2 + fCoccoN2 + fHetN + fDetN2 + fDetZ2N2    ! NEW added Coccos
       foodsq2           = food2 * food2
       grazingFlux2     = (Graz_max2 * foodsq2)/(epsilon2 + foodsq2) * Zoo2N * arrFuncZoo2
       grazingFlux_phy2  = grazingFlux2 * fphyN2 / food2
       grazingFlux_Dia2  = grazingFlux2 * fDiaN2 / food2
       grazingFlux_Cocco2 = grazingFlux2 * fCoccoN2 / food2                          ! NEW
       grazingFlux_het2  = grazingFlux2 * fHetN / food2
       grazingFlux_Det2  = grazingFlux2 * fDetN2 / food2
       grazingFlux_DetZ22  = grazingFlux2 * fDetZ2N2 / food2

       grazingFluxcarbonzoo2 = (grazingFlux_phy2 * recipQuota * grazEff2) &
                         + (grazingFlux_Dia2 * recipQuota_Dia * grazEff2) &
                         + (grazingFlux_het2 * recipQZoo * grazEff2)      &
                         + (grazingFlux_Cocco2 * recipQuota_Cocco * grazEff2) &      ! NEW
                         + (grazingFlux_Det2 * recipDet * grazEff2)       &
                         + (grazingFlux_DetZ22 *recipDet2 * grazEff2)    
    else
       food2             = fPhyN2 + fDiaN2 + fCoccoN2 + fHetN                        ! NEW
       foodsq2           = food2 * food2
       grazingFlux2      = (Graz_max2 * foodsq2)/(epsilon2 + foodsq2) * Zoo2N * arrFuncZoo2
       grazingFlux_phy2  = grazingFlux2 * fphyN2 / food2
       grazingFlux_Dia2  = grazingFlux2 * fDiaN2 / food2
       grazingFlux_Cocco2 = grazingFlux2 * fCoccoN2 / food2                          ! NEW
       grazingFlux_het2  = grazingFlux2 * fHetN / food2

       grazingFluxcarbonzoo2 = (grazingFlux_phy2 * recipQuota * grazEff2) &
                          + (grazingFlux_Dia2 * recipQuota_Dia * grazEff2) &
                          + (grazingFlux_Cocco2 * recipQuota_Cocco * grazEff2) &      ! NEW
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


!       if (REcoM_Second_Zoo) then 
!          aggregationrate = agg_PD * DetN + agg_PD * DetZ2N + agg_PP * PhyN + agg_PP * (1 - qlimitFac) * DiaN
!       else
!          aggregationrate = agg_PD * DetN                   + agg_PP * PhyN + agg_PP * (1 - qlimitFac) * DiaN
!       endif
       aggregationrate = agg_PD * DetN + agg_PP * PhyN + agg_PP * CoccoN & + agg_PP * (1 - qlimitFac) * DiaN               ! NEW added Coccos
       if (REcoM_Second_Zoo) aggregationrate = aggregationrate + agg_PD * DetZ2N


    else

!       if (REcoM_Second_Zoo) then
!          aggregationrate = agg_PD * DetN + agg_PD * DetZ2N + agg_PP * PhyN + agg_PP * DiaN
!       else
!          aggregationrate = agg_PD * DetN                   + agg_PP * PhyN + agg_PP * DiaN
!       endif
       aggregationrate = agg_PD * DetN + agg_PP * PhyN + agg_PP * CoccoN + agg_PP * DiaN                  ! NEW added Coccos
       if (REcoM_Second_Zoo) aggregationrate = aggregationrate + agg_PD * DetZ2N

    endif

!-------------------------------------------------------------------------------
! Phytoplankton and detritus aggregation
!    aggregationrate = agg_PD * DetN + agg_PP * PhyN + agg_PP * DiaN
!                    
!-------------------------------------------------------------------------------

! Terms required for the formation and dissolution of CaCO3
    if (Temp(k) < 10.6) then                                  ! VERY NEW (PICPOC definition after Krumhardt et al. 2017, 2019; Temp(k) because we need degC here)
      PICPOCtemp = 0.104d0 * Temp(k) - 0.108d0
       else
      PICPOCtemp = 1
    end if
    PICPOCtemp   = max(tiny,PICPOCtemp)
    Cal_Tlim(k)  = PICPOCtemp                                 ! NEWOUT

    PICPOCCO2     = 1.102e+00 * HCO3_watercolumn(k) * Cunits / (4.238e+01 + HCO3_watercolumn(k) * Cunits) - exp(-7.079e-01 * CO2_watercolumn(k) * Cunits) - 1.343e+07 * 10.**(-pH_watercolumn(k))
    !PICPOCCO2     = 1.212e+00 * HCO3_watercolumn(k) * Cunits / (4.662e+01 + HCO3_watercolumn(k) * Cunits) - exp(-7.787e-01 * CO2_watercolumn(k) * Cunits) - 1.477e+07 * 10.**(-pH_watercolumn(k)) ! lowSens
    !PICPOCCO2     = min(PICPOCCO2,real(one))
    PICPOCCO2     = min(PICPOCCO2,3.0)                                       ! April 2022: limitation to 3
    !if(recom_mype==1) write(*,*), 'PICPOCCO2 =',PICPOCCO2
    Cal_CO2lim(k) = PICPOCCO2                                                ! NEWOUT

    PICPOCN     = -0.31 * (DIN/(DIN + k_din_c)) + 1.31
    PICPOCN     = max(tiny,PICPOCN)
    Cal_Nlim(k) = PICPOCN                                                    ! NEWOUT
    
    if (CO2lim) then
      calcification = 1.d0 * Cphot_cocco * CoccoC * PICPOCtemp * PICPOCN * PICPOCCO2  ! VERY NEW
    else
      calcification = 1.d0 * Cphot_cocco * CoccoC * PICPOCtemp * PICPOCN
    endif
    Cal_pure(k)   = calcification                                            ! NEWOUT


!< calc_prod_ratio: Calcite production ratio, dimensionless
!    calcification = calc_prod_ratio * Cphot * PhyC   ! Z in equations    ! OG
    calc_loss_agg = aggregationrate * PhyCalc

!    if(REcoM_Second_Zoo)  then
!       calc_loss_gra  =  grazingFlux_phy   &
!                       * recipQuota/(PhyC + tiny)    * PhyCalc
!       calc_loss_gra2 = grazingFlux_phy2   &
!                       * recipQuota/(PhyC + tiny)    * PhyCalc
!    else
!
!       calc_loss_gra = grazingFlux_phy     &
!                       * recipQuota/(PhyC + tiny)    * PhyCalc
!    endif

!    calc_loss_gra = grazingFlux_phy          &
!                  * recipQuota/(PhyC + tiny) &  
!                  * PhyCalc

!    if (REcoM_Second_Zoo) calc_loss_gra2 = grazingFlux_phy2         &
!                                         * recipQuota/(PhyC + tiny) &
!                                         * PhyCalc

  if(REcoM_Second_Zoo)  then
     calc_loss_gra  =  grazingFlux_Cocco   &                                 ! NEW changed to Coccos
                     * recipQuota_Cocco/(CoccoC + tiny)    * PhyCalc         ! NEW changed to Coccos
     calc_loss_gra2 = grazingFlux_Cocco2           &                         ! NEW changed to Coccos
                     * recipQuota_Cocco/(CoccoC + tiny)    * PhyCalc         ! NEW changed to Coccos
   else
    calc_loss_gra = grazingFlux_Cocco           &                            ! NEW changed to Coccos
                    * recipQuota_Cocco/(CoccoC + tiny)    * PhyCalc          ! NEW changed to Coccos
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
        - N_assim_Cocco                * CoccoC  &          ! NEW
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
            - Cphot_Cocco                   * CoccoC  &          ! NEW
            + phyRespRate_Cocco             * CoccoC  &          ! NEW
            + rho_C1 * arrFunc              * EOC        & ! --> Remineralization of DOC
            + HetRespFlux                                & ! --> Zooplankton respiration                     
            + Zoo2RespFlux                               &                    
            + calc_diss                     * DetCalc    & ! --> Calcite dissolution from detritus 
            + calc_loss_gra  * calc_diss_guts            &
            + calc_loss_gra2 * calc_diss_guts            &
            + calc_diss2                     * DetZ2Calc & 
            - calcification                              & ! --> Calcification
                                                        ) * dt_b + sms(k,idic)
    else
        sms(k,idic)      = (                       &
            - Cphot                         * PhyC    &
            + phyRespRate                   * PhyC    &
            - Cphot_Dia                     * DiaC    &
            + phyRespRate_Dia               * DiaC    &
            - Cphot_Cocco                   * CoccoC  &          ! NEW
            + phyRespRate_Cocco             * CoccoC  &          ! NEW
            + rho_C1 * arrFunc              * EOC     &
            + HetRespFlux                             & 
                    
            + calc_diss                     * DetCalc &
            + calc_loss_gra * calc_diss_guts          &
            - calcification                           &
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
            + 1.0625 * N_assim_Cocco       * CoccoC  &           ! NEW
            - 1.0625 * rho_N * arrFunc     * DON     &     
            + 2.d0 * calc_diss             * DetCalc &
            + 2.d0 * calc_loss_gra * calc_diss_guts  &
            + 2.d0 * calc_loss_gra2 * calc_diss_guts &
            + 2.d0 * calc_diss2            * DetZ2Calc &
            - 2.d0 * calcification                   &                      
                                             ) * dt_b + sms(k,ialk)
    else
        sms(k,ialk)      = (                       &
            + 1.0625 * N_assim             * PhyC    &
            + 1.0625 * N_assim_Dia         * DiaC    &
            + 1.0625 * N_assim_Cocco       * CoccoC  &           ! NEW
            - 1.0625 * rho_N * arrFunc     * DON     &
            + 2.d0 * calc_diss             * DetCalc &
            + 2.d0 * calc_loss_gra * calc_diss_guts  &
            - 2.d0 * calcification                   &
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
        + grazingFlux_Cocco                      &              ! NEW
        - grazingFlux_Cocco * grazEff            &              ! NEW
        - grazingFlux_Det * grazEff              & !!!Sloppy feeding is  thought because of grazing flux multiplied with grazeff 
        - grazingFlux_Det2 * grazEff             &
        + aggregationRate              * PhyN    &
        + aggregationRate              * DiaN    &
        + aggregationRate              * CoccoN  &              ! NEW
        + hetLossFlux                            &
        - reminN * arrFunc             * DetN    &
                                               ) * dt_b + sms(k,idetn)
   else
    sms(k,idetn)       = (                       &
        + grazingFlux_phy                        &
        + grazingFlux_dia                        &
        + grazingFlux_Cocco                      &              ! NEW
        - grazingFlux * grazEff                  &
        + aggregationRate              * PhyN    &
        + aggregationRate              * DiaN    &
        + aggregationRate              * CoccoN  &              ! NEW
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
        + grazingFlux_Cocco * recipQuota_Cocco         &         ! NEW
        - grazingFlux_Cocco * recipQuota_Cocco * grazEff &       ! NEW
        - grazingFlux_Det * recipDet * grazEff         &
        - grazingFlux_Det2 * recipDet2 * grazEff       &     
        + aggregationRate              * phyC          &
        + aggregationRate              * DiaC          &
        + aggregationRate              * CoccoC        &         ! NEW
        + hetLossFlux * recipQZoo                      &
        - reminC * arrFunc             * DetC          &
                                             )   * dt_b + sms(k,idetc)
   else
    sms(k,idetc)       = (                             &
        + grazingFlux_phy * recipQuota                 &
        - grazingFlux_phy * recipQuota * grazEff       &
        + grazingFlux_Dia * recipQuota_Dia             &
        - grazingFlux_Dia * recipQuota_Dia * grazEff   &
        + grazingFlux_Cocco * recipQuota_Cocco         &         ! NEW
        - grazingFlux_Cocco * recipQuota_Cocco * grazEff &       ! NEW
        + aggregationRate              * phyC          &
        + aggregationRate              * DiaC          &
        + aggregationRate              * CoccoC        &         ! NEW
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
                + grazingFlux_Cocco*recipQuota_Cocco*grazEff &            ! NEW
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
                + grazingFlux_Cocco*recipQuota_Cocco*grazEff &             ! NEW
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
            + grazingFlux_Cocco*recipQuota_Cocco*grazEff &             ! NEW
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
         + grazingFlux_Cocco2*recipQuota_Cocco*grazEff2 &         ! NEW
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
         + grazingFlux_Cocco2*recipQuota_Cocco*grazEff2 &          ! NEW
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
         + grazingFlux_Cocco2                     &                ! NEW
         - grazingFlux_Cocco2 * grazEff2          &                ! NEW 
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
         + grazingFlux_Cocco2                     &                ! NEW
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
        + grazingFlux_Cocco2 * recipQuota_Cocco        &           ! NEW
        - grazingFlux_Cocco2*recipQuota_Cocco*grazEff2 &           ! NEW
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
        + grazingFlux_Cocco2 * recipQuota_Cocco        &           ! NEW
        - grazingFlux_Cocco2*recipQuota_Cocco*grazEff2 &           ! NEW
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
      + lossN_c * limitFacN_Cocco      * CoccoN &                  ! NEW
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
      + lossN_c * limitFacN_Cocco      * CoccoN &                  ! NEW
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
      + lossC_c * limitFacN_cocco      * CoccoC &                  ! NEW
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
      + lossC_c * limitFacN_cocco      * CoccoC &                  ! NEW
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
! Coccolithophore N (NEW!!!)
   if (REcoM_Second_Zoo) then                                           ! NEW
    sms(k,icocn)      = (                      &                        ! NEW
      + N_assim_cocco                 * CoccoC &                        ! NEW
      - lossN_c * limitFacN_cocco     * CoccoN &                        ! NEW
      - aggregationRate               * CoccoN &                        ! NEW
      - grazingFlux_Cocco                      &                        ! NEW
      - grazingFlux_Cocco2                     &                        ! NEW
                                             ) * dt + sms(k,icocn)      ! NEW
   else                                                                 ! NEW
    sms(k,icocn)      = (                      &                        ! NEW
      + N_assim_cocco                 * CoccoC &                        ! NEW
      - lossN_c * limitFacN_cocco     * CoccoN &                        ! NEW
      - aggregationRate               * CoccoN &                        ! NEW
      - grazingFlux_Cocco                      &                        ! NEW
                                                 ) * dt + sms(k,icocn)  ! NEW
   endif                                                                ! NEW
!-------------------------------------------------------------------------------
! Coccolithophore C (NEW!!!)
    if (REcoM_Second_Zoo) then                                          ! NEW
    sms(k,icocc)      = (                      &                        ! NEW
      + Cphot_cocco                   * CoccoC &                        ! NEW
      - lossC_c * limitFacN_cocco     * CoccoC &                        ! NEW
      - phyRespRate_cocco             * CoccoC &                        ! NEW
      - aggregationRate               * CoccoC &                        ! NEW
      - grazingFlux_cocco * recipQuota_cocco   &                        ! NEW
      - grazingFlux_Cocco2* recipQuota_cocco   &                        ! NEW
                                             ) * dt + sms(k,icocc)      ! NEW
   else                                                                 ! NEW
    sms(k,icocc)      = (                      &                        ! NEW
      + Cphot_cocco                   * CoccoC &                        ! NEW
      - lossC_c * limitFacN_cocco     * CoccoC &                        ! NEW
      - phyRespRate_cocco             * CoccoC &                        ! NEW
      - aggregationRate               * CoccoC &                        ! NEW
      - grazingFlux_cocco * recipQuota_cocco   &                        ! NEW
                                             ) * dt + sms(k,icocc)      ! NEW
   endif

!-------------------------------------------------------------------------------
! Coccolithophore Chl (NEW!!!)
    if (REcoM_Second_Zoo) then                                          ! NEW
    sms(k,icchl)      = (                         &                     ! NEW
      + ChlSynth_cocco                 * CoccoC   &                     ! NEW
      - KOchl_cocco                    * CoccoChl &                     ! NEW
      - aggregationRate                * CoccoChl &                     ! NEW
      - grazingFlux_cocco * Chl2N_cocco           &                     ! NEW
      - grazingFlux_Cocco2* Chl2N_cocco           &                     ! NEW
                                             ) * dt + sms(k,icchl)      ! NEW
   else                                                                 ! NEW
    sms(k,icchl)      = (                         &                     ! NEW
      + chlSynth_cocco                 * CoccoC   &                     ! NEW
      - KOchl_cocco                    * CoccoChl &                     ! NEW
      - aggregationRate                * CoccoChl &                     ! NEW
      - grazingFlux_cocco * Chl2N_cocco           &                     ! NEW
                                             ) * dt + sms(k,icchl)      ! NEW
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
                - N_assim_cocco           * CoccoC    &                      ! NEW
                + lossN*limitFacN         * PhyN     & ! --> Excretion from small pythoplankton
                + lossN_d*limitFacN_dia   * DiaN     & ! --> Excretion from diatom
                + lossN_c*limitFacN_cocco * CoccoN    &                      ! NEW
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
                -  Cphot_cocco            * CoccoC   &                         ! NEW
                +  phyRespRate            * PhyC     & ! Small pyhtoplankton respiration --- /
                +  phyRespRate_dia        * DiaC     & ! Diatom respiration
                +  phyRespRate_cocco      * CoccoC   &                         ! NEW
                +  lossC*limitFacN        * phyC     & ! Exrcetion from small pythoplankton
                +  lossC_d*limitFacN_dia  * diaC     & ! Excretion from diatom
                +  lossC_c*limitFacN_cocco* CoccoC   &                         ! NEW
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
                - N_assim_cocco           * CoccoC    &                      ! NEW
                + lossN*limitFacN         * PhyN      &
                + lossN_d*limitFacN_dia   * DiaN      &
                + lossN_c*limitFacN_cocco * CoccoN    &                      ! NEW
                + reminN * arrFunc        * DetN      &
                + lossN_z                 * HetN )    &
                - kScavFe                 * DetC * FreeFe &
                                              ) * dt_b + sms(k,ife)
        else
            sms(k,ife)      = ( Fe2C *(          &
                -  Cphot                  * PhyC     &
                -  Cphot_dia              * DiaC     &
                -  Cphot_cocco            * CoccoC   &                         ! NEW
                +  phyRespRate            * PhyC     &
                +  phyRespRate_dia        * DiaC     &
                +  phyRespRate_cocco      * CoccoC   &                         ! NEW
                +  lossC*limitFacN        * phyC     &
                +  lossC_d*limitFacN_dia  * diaC     &
                +  lossC_c*limitFacN_cocco* CoccoC   &                         ! NEW
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
            - lossC_c * limitFacN_cocco * phyCalc &                             ! NEW
            - phyRespRate_cocco         * phyCalc &                             ! NEW
            !- lossC * limitFacN * phyCalc & ! -- Excretion loss
            !- phyRespRate       * phyCalc & ! -- Respiration
            - calc_loss_agg               & ! -- Aggregation loss
            - calc_loss_gra               & ! -- Grazing loss
            - calc_loss_gra2              &
                                                  ) * dt_b + sms(k,iphycal)
    else
        sms(k,iphycal)    = (             &
            + calcification               &
            - lossC_c * limitFacN_cocco * phyCalc &                             ! NEW
            - phyRespRate_cocco         * phyCalc &                             ! NEW
            !- lossC * limitFacN * phyCalc &
            !- phyRespRate       * phyCalc &
            - calc_loss_agg               &
            - calc_loss_gra               &
                                                  ) * dt_b + sms(k,iphycal)
    endif

!-------------------------------------------------------------------------------
! Detritus calcite
    sms(k,idetcal)   = (               &
      + lossC_c * limitFacN_cocco * phyCalc &                            ! VERY NEW (changed from phy to cocco)
      + phyRespRate_cocco         * phyCalc &                            ! VERY NEW (changed from phy to cocco)
      !+ lossC * limitFacN * phyCalc    &
      !+ phyRespRate       * phyCalc    &
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
      + Cphot_cocco         * CoccoC&               ! NEW
      - phyRespRate_cocco   * CoccoC&               ! NEW
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
      + Cphot_cocco         * CoccoC&               ! NEW
      - phyRespRate_cocco   * CoccoC&               ! NEW
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

        Diags3Dloc(k,21) = Diags3Dloc(k,21) + (   &                      ! NEW
        + Cphot_cocco             * CoccoC  &                            ! NEW
        - PhyRespRate_cocco       * CoccoC  &                            ! NEW
        ) * recipbiostep

!*** Gross primary production [mmol C /(m3 * day)]
	Diags3Dloc(k,3) = Diags3Dloc(k,3) + (   &
     	+ Cphot                   * PhyC  &
     	) * recipbiostep

	Diags3Dloc(k,4) = Diags3Dloc(k,4) + (   &
     	+ Cphot_dia               * DiaC  &
     	) * recipbiostep

       Diags3Dloc(k,22) = Diags3Dloc(k,22) + (   &                       ! NEW                                                                                                   
       + Cphot_cocco             * CoccoC  &
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

        Diags3Dloc(k,23) = Diags3Dloc(k,23) + (   &                       ! NEW                                                      
        + N_assim_cocco           * CoccoC  &
        - lossN * limitFacN_cocco * CoccoN  &
        ) * recipbiostep

!*** Changed to chlorophyll degradation (commented out gross N-assimilation below)
        Diags3Dloc(k,7) = Diags3Dloc(k,7) + (   &
        + KOchl  &
        ) * recipbiostep

        Diags3Dloc(k,8) = Diags3Dloc(k,8) + (   &
        + KOchl_dia  &
        ) * recipbiostep

        Diags3Dloc(k,24) = Diags3Dloc(k,24) + ( &                         ! NEW
        + N_assim_cocco            * CoccoC &                             ! NEW
        ) * recipbiostep                                                  ! NEW

! diagnostics, combined from Onur and Cara
!*** Total grazing of first zooplankton (with graz_eff, i.e. what reaches ZOO)
        Diags3Dloc(k,9) = Diags3Dloc(k,9) + (   &                             
        + grazingFlux_phy * recipQuota * grazEff     &
        + grazingFlux_Dia * recipQuota_Dia * grazEff &
        + grazingFlux_Cocco * recipQuota_Cocco * grazEff &                ! NEW
        ) * recipbiostep                                 

!*** Grazing on small Phytoplankton by First Zooplankton (without grazeff, i.e. loss term for PHY)
        Diags3Dloc(k,10) = Diags3Dloc(k,10)+ (   &   
        + grazingFlux_phy * recipQuota           &
        ) * recipbiostep                                 

!*** Grazing on diatoms by First Zooplankton (without grazeff, i.e. loss term for DIA)
        Diags3Dloc(k,11) = Diags3Dloc(k,11) + (  &                                                 
        + grazingFlux_dia * recipQuota_dia       & 
        ) * recipbiostep  

!*** Grazing on cocclithophores by First Zooplankton (without grazeff, i.e. loss term for COCCO)
        Diags3Dloc(k,25) = Diags3Dloc(k,25) +(   &                         ! NEW
        + grazingFlux_Cocco * recipQuota_cocco   &                         ! NEW
        ) * recipbiostep                                                   ! NEW

!*** zooplankton1 respiration
        Diags3Dloc(k,12) = Diags3Dloc(k,12) + (   &
        + HetRespFlux                             &
        ) * recipbiostep

!*** calc_diss
        Diags3Dloc(k,13) = Diags3Dloc(k,13) + (   &
        + calc_diss * DetCalc                     &
        ) * recipbiostep

!***    aggregation by  small phytoplankton                                                                      
        Diags3Dloc(k,14) = Diags3Dloc(k,14) + (   &   
        + aggregationrate * PhyC                  &
        ) * recipbiostep

!***    aggregation by  diatoms                                                                                  
        Diags3Dloc(k,15) = Diags3Dloc(k,15) + (   &
        + aggregationrate * DiaC                  &
        ) * recipbiostep  

!***    aggregation by coccolithophores
        Diags3Dloc(k,26) = Diags3Dloc(k,26) + (   &                         ! NEW
        + aggregationrate * CoccoC                &                         ! NEW
        ) * recipbiostep                                                    ! NEW

!*** excrection of DOC by phytoplankton
        Diags3Dloc(k,16) = Diags3Dloc(k,16) + (   &
        + lossC * limitFacN              * phyC   &
        ) * recipbiostep  
  
!*** excrection of DOC by diatoms
        Diags3Dloc(k,17) = Diags3Dloc(k,17) + (   &
        + lossC_d * limitFacN_dia        * DiaC   &
        ) * recipbiostep  

!*** excretion of DOC by coccolithophores
        Diags3Dloc(k,27) = Diags3Dloc(k,27) + (   &                          ! NEW
        + lossC_c * limitFacN_cocco      * CoccoC &                          ! NEW
        ) * recipbiostep                                                     ! NEW

!*** calcification
        Diags3Dloc(k,18) = Diags3Dloc(k,18) + (   &
        + calcification                           &
        ) * recipbiostep  

! phy respiration
	Diags3Dloc(k,19) = Diags3Dloc(k,19) + (   &
     	+ PhyRespRate             * PhyC          &
     	) * recipbiostep

! dia respiration
	Diags3Dloc(k,20) = Diags3Dloc(k,20) + (   &
     	+ PhyRespRate_dia         * DiaC          &
     	) * recipbiostep

! cocco resipration
        Diags3Dloc(k,28) = Diags3Dloc(k,28) + (   &                          ! NEW
        + PhyRespRate_cocco       * CoccoC        &                          ! NEW
        ) * recipbiostep                                                     ! NEW

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

!*** Calc: DIC, Alk ***
  decayBenthos(4) = calc_diss_ben * LocBenthos(4)              ! NEW DISS changed calc_diss to calc_diss_ben to not make the dissolution omega dependent when using the switch OmegaC_diss
  LocBenthos(4)      = LocBenthos(4)   - decayBenthos(4) * dt_b

  if (ciso) then
!*** DIC_13 ***  We ignore isotopic fractionation during remineralization.
     decayBenthos(5) = alpha_dcal_13   * decayRateBenC   * LocBenthos(5)
     LocBenthos(5)   = LocBenthos(5)   - decayBenthos(5) * dt_b
!*** DIC_14 ***  We ignore isotopic fractionation during remineralization.
     decayBenthos(6) = alpha_dcal_14   * decayRateBenC   * LocBenthos(6)
     LocBenthos(6)   = LocBenthos(6)   - decayBenthos(6) * dt_b
!*** Calc: DIC_13,14 ***
     decayBenthos(7) = calc_diss_13    * LocBenthos(7)
     LocBenthos(7)   = LocBenthos(7)   - decayBenthos(7) * dt_b ! / depth of benthos
     decayBenthos(8) = calc_diss_14    * LocBenthos(8)
     LocBenthos(8)   = LocBenthos(8)   - decayBenthos(8) * dt_b ! / depth of benthos
  end if ! ciso

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

