!===============================================================================
! REcoM_Forcing
!===============================================================================
subroutine REcoM_Forcing(zNodes, n, Nn, state, SurfSW, Loc_slp , Temp, Sali, Sali_depth &
            , CO2_watercolumn                                          &
            , pH_watercolumn                                           &
            , pCO2_watercolumn                                         &
            , HCO3_watercolumn                                         &
            , CO3_watercolumn                                          &
            , OmegaC_watercolumn                                       &
            , kspc_watercolumn                                         &
            , rhoSW_watercolumn                                        &
            , PAR, ice, dynamics, tracers, partit, mesh)

    use recom_declarations
    use recom_locvar
    use recom_config
    use recom_glovar
    use gasx
    use recom_ciso
    use g_clock
    use o_PARAM
    use g_rotate_grid
    use g_config
    use mod_mesh
    USE MOD_PARTIT
    USE MOD_PARSUP
    use mod_tracer
    use MOD_DYN
    use MOD_ICE

    use o_param
    use o_arrays
    use g_forcing_arrays
    use g_comm_auto
    use g_support
    implicit none

    type(t_dyn)   , intent(inout), target :: dynamics
    type(t_ice)   , intent(inout), target :: ice
    type(t_tracer), intent(inout), target :: tracers
    type(t_partit), intent(inout), target :: partit
    type(t_mesh)  , intent(inout), target :: mesh

    real(kind=8)                              :: Latr          
    integer                                   :: n, Nn         ! Nn -Total number of nodes
    real(kind=8),dimension(mesh%nl-1)	      :: zNodes	       ! Depth of nodes   zr(1:nzmax) = Z_3d_n(1:nzmax,n)
    real(kind=8),dimension(mesh%nl-1,bgc_num) :: state         	
    real(kind=8)                              :: SurfSW        ! [W/m2] ShortWave radiation at surface
    real(kind=8)                              :: Loc_slp       ! [Pa] sea-level pressure
    real(kind=8),dimension(mesh%nl-1)         :: Temp          ! [degrees C] Ocean temperature
    real(kind=8),dimension(mesh%nl-1)         :: Sali_depth    ! Salinity for the whole water column

    !!---- Watercolumn carbonate chemistry
    real(kind=8),dimension(mesh%nl-1)         :: CO2_watercolumn
    real(kind=8),dimension(mesh%nl-1)         :: pH_watercolumn
    real(kind=8),dimension(mesh%nl-1)         :: pCO2_watercolumn
    real(kind=8),dimension(mesh%nl-1)         :: HCO3_watercolumn
    real(kind=8),dimension(mesh%nl-1)         :: CO3_watercolumn
    real(kind=8),dimension(mesh%nl-1)         :: OmegaC_watercolumn
    real(kind=8),dimension(mesh%nl-1)         :: kspc_watercolumn
    real(kind=8),dimension(mesh%nl-1)         :: rhoSW_watercolumn

    real(kind=8),dimension(mesh%nl-1)         :: PAR

    !!---- Subroutine Depth

    real(kind=8),dimension(mesh%nl)           :: zF                   ! [m] Depth of fluxes
    real(kind=8),dimension(mesh%nl,5)         :: SinkVel              ! [m/day]
    real(kind=8),dimension(mesh%nl-1)         :: thick                ! [m] Vertical distance between two nodes = Thickness 
    real(kind=8),dimension(mesh%nl-1)         :: recipthick           ! [1/m] reciprocal of thick

    !!---- Subroutine CO2Flux /mocsy
    real(kind=8)                              :: REcoM_DIC(1)         ! [mol/m3] Conc of DIC in the surface water, used to calculate CO2 flux
    real(kind=8)                              :: REcoM_Alk(1)         ! [mol/m3] Conc of Alk in the surface water, used to calculate CO2 flux
    real(kind=8)                              :: REcoM_Si(1)          ! [mol/m3] Conc of Si in the surface water, used to calculate CO2 flux
    real(kind=8)                              :: REcoM_Phos(1)        ! [mol/m3] Conc of Phos in the surface water, used to calculate the CO2 flux
    real(kind=8)                              :: Sali(1)              ! Salinity of current surface layer
    real(kind=8)                              :: Latd(1)              ! latitude in degree
    real(kind=8)                              :: Lond(1)              ! longitude in degree
    real(kind=8)                              :: REcoM_T(1)           ! temperature again, for mocsy minimum defined as -2
    real(kind=8)                              :: REcoM_S(1)           ! temperature again, for mocsy minimum defined as 21
! atm pressure, now read in as forcing!!
    !!---- atm pressure
    real(kind=8)                              :: Patm(1)              ! atmospheric pressure [atm]

    !!---- Subroutine o2flux /mocsy 
    real(kind=8)                              :: ppo(1)               ! atmospheric pressure, divided by 1 atm 
    real(kind=8)                              :: REcoM_O2(1)          ! [mmol/m3] Conc of O2 in the surface water, used to calculate O2 flux

    !!---- Subroutine REcoM_sms
    real(kind=8),dimension(mesh%nl-1,bgc_num) :: sms                  ! matrix that entail changes in tracer concentrations

    !!---- Diagnostics
    integer                                   :: idiags,k

    integer                    :: tr_num

#include "../associate_part_def.h"
#include "../associate_mesh_def.h"
#include "../associate_part_ass.h"
#include "../associate_mesh_ass.h"


    tiny_N   = tiny_chl/chl2N_max   ! 0.00001/ 3.15d0   Chl2N_max [mg CHL/mmol N] Maximum CHL a : N ratio = 0.3 gCHL gN^-1
    tiny_N_d = tiny_chl/chl2N_max_d ! 0.00001/ 4.2d0

    tiny_C   = tiny_N  /NCmax       ! NCmax   = 0.2d0   [mmol N/mmol C] Maximum cell quota of nitrogen (N:C)
    tiny_C_d = tiny_N_d/NCmax_d     ! NCmax_d = 0.2d0 

    tiny_Si  = tiny_C_d/SiCmax      ! SiCmax = 0.8d0

#if defined (__coccos)
    tiny_N_c = tiny_chl/chl2N_max_c ! 0.00001/ 3.5d0 
    tiny_C_c = tiny_N_c/NCmax_c     ! NCmax_c = 0.15d0
#endif

    call Cobeta(partit, mesh)      
    call Depth_calculations(n, Nn,SinkVel,zF,thick,recipthick, partit, mesh)

    !! *** Mocsy ***

    !!---- convert from mmol/m3 to mol/m3
    REcoM_DIC  = max(tiny*1e-3, state(one,idic)*1e-3)
    REcoM_Alk  = max(tiny*1e-3, state(one,ialk)*1e-3)
    REcoM_Si   = max(tiny*1e-3, state(one,isi) *1e-3)

    !!---- convert N to P with Redfield ratio
    REcoM_Phos = max(tiny*1e-3, state(one,idin)*1e-3) /16.

    !!---- minimum set to 2 degC: K1/K2 Lueker valid between 2degC-35degC and 19-43psu
    REcoM_T    = max(2.d0, Temp(1))
    !!---- maximum set to 40 degC: K1/K2 Lueker valid between 2degC-35degC and 19-43psu
    REcoM_T    = min(REcoM_T, 40.d0) 

    !!---- minimum set to 21: K1/K2 Lueker valid between 2degC-35degC and 19-43psu, else causes trouble in regions with S between 19 and 21 and ice conc above 97%
    REcoM_S    = max(21.d0, Sali(1)) 
    !!---- maximum set to 43: K1/K2 Lueker valid between 2degC-35degC and 19-43psu, else causes trouble   REcoM_S    = min(REcoM_S, 43.d0)  !!!!!!!!

    !!---- convert from Pa to atm.
    Patm = Loc_slp/Pa2atm  

    !!---- lon
    Lond=geo_coord_nod2D(1,n)/rad !! convert from rad to degree
    !!---- lat
    Latr=geo_coord_nod2D(2,n)
    Latd=geo_coord_nod2D(2,n)/rad !! convert from rad to degree

    !!---- calculate piston velocity kw660, which is an input to the flxco2 calculation
    !!---- pistonvel already scaled for ice-free area
    !!---- compute piston velolicty kw660 (at 25 C) from wind speed
    !!---- BUT without Schmidt number temperature correction (Sc differs each gas)
    !! ULoc: wind speed at 10-m height
    !! Loc_ice_conc: modeled sea-ice cover: fraction of grid cell, varying between 0.0 (no ice) and 1.0 (full cover)
    !! kw660: piston velocity at 25°C [m/s], uncorrected by the Schmidt number for different temperatures

    call pistonvel(ULoc, Loc_ice_conc, Nmocsy, kw660)

    !! *** check ***

    if((REcoM_DIC(1) > 10000.d0)) then               ! NEW: added this entire print statement (if to endif)
        print*, 'NEW ERROR: DIC !'  
        print*, 'pco2surf: ',pco2surf
        print*, 'co2: ',co2
        print*, 'rhoSW: ', rhoSW
        print*, 'temp: ',REcoM_T
        print*, 'tempis: ',tempis
        print*, 'REcoM_S: ', REcoM_S
        print*, 'REcoM_Alk: ', REcom_Alk
        print*, 'REcoM_DIC: ', REcoM_DIC
        print*, 'REcoM_Si: ', REcoM_Si
        print*, 'REcoM_Phos: ', REcoM_Phos
        print*, 'kw660: ',kw660
        print*, 'LocAtmCO2: ', LocAtmCO2
        print*, 'Patm: ', Patm
        print*, 'thick(One): ',thick(One)
        print*, 'Nmocsy: ', Nmocsy
        print*, 'Lond: ', Lond
        print*, 'Latd: ', Latd
        print*, 'ULoc: ', ULoc
        print*, 'Loc_ice_conc: ', Loc_ice_conc
        stop
    endif

    call flxco2(co2flux, co2ex, dpco2surf,                                                   &
                ph, pco2surf, fco2, co2, hco3, co3, OmegaA, OmegaC, BetaD, rhoSW, p, tempis, &
                REcoM_T, REcoM_S, REcoM_Alk, REcoM_DIC, REcoM_Si, REcoM_Phos, kw660, LocAtmCO2, Patm, thick(One), Nmocsy, Lond,Latd, &
                optCON='mol/m3',optT='Tpot   ',optP='m ',optB='u74',optK1K2='l  ',optKf='dg',optGAS='Pinsitu',optS='Sprc')

! changed optK1K2='l  ' to 'm10'
  if((co2flux(1)>1.e10) .or. (co2flux(1)<-1.e10)) then
!     co2flux(1)=0.0  
      print*, 'ERROR: co2 flux !'
      print*, 'pco2surf: ',pco2surf
      print*, 'co2: ',co2
      print*, 'rhoSW: ', rhoSW
      print*, 'temp: ',REcoM_T
      print*, 'tempis: ',tempis
      print*, 'REcoM_S: ', REcoM_S
      print*, 'REcoM_Alk: ', REcom_Alk
      print*, 'REcoM_DIC: ', REcoM_DIC
      print*, 'REcoM_Si: ', REcoM_Si
      print*, 'REcoM_Phos: ', REcoM_Phos
      print*, 'kw660: ',kw660
      print*, 'LocAtmCO2: ', LocAtmCO2
      print*, 'Patm: ', Patm
      print*, 'thick(One): ',thick(One) 
      print*, 'Nmocsy: ', Nmocsy
      print*, 'Lond: ', Lond
      print*, 'Latd: ', Latd   
      print*, 'ULoc: ', ULoc
      print*, 'Loc_ice_conc: ', Loc_ice_conc
      stop
    endif

! use ice-free area and also convert from mol/m2/s to mmol/m2/d
!   if(mype==0) write(*,*), 'co2flux (mol/m2/s) =',co2flux

! ice-fraction is already considered in piston-velocity, so don't apply it here
   dflux     = co2flux * 1.e3 *SecondsPerDay  !* (1.d0 - Loc_ice_conc)
!   if(mype==0) write(*,*), 'dflux (mmol/m2/d) =',dflux

   co2flux_seaicemask = co2flux * 1.e3 !  [mmol/m2/s]  * (1.d0 - Loc_ice_conc)
!   if(mype==0) write(*,*), 'co2flux_seaicemask (mmol/m2/s) =',co2flux_seaicemask

! then oxygen
   ppo = Loc_slp/Pa2atm !1 !slp divided by 1 atm
   REcoM_O2 = max(tiny*1e-3,state(one,ioxy)*1e-3) ! convert from mmol/m3 to mol/m3 for mocsy

   call  o2flux(REcoM_T, REcoM_S, kw660, ppo, REcoM_O2, Nmocsy, o2ex)
   oflux     = o2ex * 1.e3 *SecondsPerDay  !* (1.d0 - Loc_ice_conc) [mmol/m2/d]
   o2flux_seaicemask = o2ex * 1.e3 ! back to mmol here [mmol/m2/s] 

! Source-Minus-Sinks

if (recom_debug .and. mype==0) print *, achar(27)//'[36m'//'     --> REcoM_sms'//achar(27)//'[0m'

!  call REcoM_sms(n, Nn, state, thick, recipthick, SurfSW, sms, Temp ,zF, PAR, mesh)

  call REcoM_sms(n, Nn, state, thick, recipthick, SurfSW, sms, Temp, Sali_depth &
        , CO2_watercolumn                                              & ! MOCSY [mol/m3]
        , pH_watercolumn                                               & ! MOCSY on total scale
        , pCO2_watercolumn                                             & ! MOCSY [uatm]
        , HCO3_watercolumn                                             & ! MOCSY [mol/m3]
        , CO3_watercolumn                                              & ! DISS [mol/m3]
        , OmegaC_watercolumn                                           & ! DISS calcite saturation state
        , kspc_watercolumn                                             & ! DISS stoichiometric solubility product [mol^2/kg^2]
        , rhoSW_watercolumn                                            & ! DISS in-situ density of seawater [kg/m3]
        , Loc_slp                                                      &
        , zF, PAR, Lond, Latd, ice, dynamics, tracers, partit, mesh)

  state(1:nn,:)      = max(tiny,state(1:nn,:) + sms(1:nn,:))

  state(1:nn,ipchl)  = max(tiny_chl,state(1:nn,ipchl))
  state(1:nn,iphyn)  = max(tiny_N,  state(1:nn,iphyn))
  state(1:nn,iphyc)  = max(tiny_C,  state(1:nn,iphyc))
  state(1:nn,idchl)  = max(tiny_chl,state(1:nn,idchl))
  state(1:nn,idian)  = max(tiny_N_d,state(1:nn,idian))
  state(1:nn,idiac)  = max(tiny_C_d,state(1:nn,idiac))
  state(1:nn,idiasi) = max(tiny_Si, state(1:nn,idiasi))

#if defined (__coccos)
  state(1:nn,icchl)  = max(tiny_chl,state(1:nn,icchl))
  state(1:nn,icocn)  = max(tiny_N_c,state(1:nn,icocn))
  state(1:nn,icocc)  = max(tiny_C_c,state(1:nn,icocc))
#endif

#if defined (__3Zoo2Det)
  state(1:nn,imiczoon)  = max(tiny,state(1:nn,imiczoon))
  state(1:nn,imiczooc)  = max(tiny,state(1:nn,imiczooc))
#endif

if (recom_debug .and. mype==0) print *, achar(27)//'[36m'//'     --> ciso after REcoM_Forcing'//achar(27)//'[0m'

!-------------------------------------------------------------------------------
! Diagnostics
  if (Diags) then

!    logical, optional                 :: lNPPn

!    if (present(lNPPn))then
!        locNPPn = sum(diags3Dloc(1:nn,idiags) * thick(1:nn))
!    endif
     locNPPn = sum(vertNPPn(1:nn) * thick(1:nn))
     locGPPn = sum(vertGPPn(1:nn) * thick(1:nn))
     locNNAn = sum(vertNNAn(1:nn) * thick(1:nn))
     locChldegn = sum(vertChldegn(1:nn) * thick(1:nn))

     locNPPd = sum(vertNPPd(1:nn) * thick(1:nn))
     locGPPd = sum(vertGPPd(1:nn) * thick(1:nn))
     locNNAd = sum(vertNNAd(1:nn) * thick(1:nn))
     locChldegd = sum(vertChldegd(1:nn) * thick(1:nn))

     locNPPc = sum(vertNPPc(1:nn) * thick(1:nn))
     locGPPc = sum(vertGPPc(1:nn) * thick(1:nn))
     locNNAc = sum(vertNNAc(1:nn) * thick(1:nn))
     locChldegc = sum(vertChldegc(1:nn) * thick(1:nn))

  end if
end subroutine REcoM_Forcing
