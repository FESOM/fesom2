!===============================================================================
! REcoM_Forcing
!===============================================================================
subroutine REcoM_Forcing(zNodes, n, Nn, state, SurfSW, Loc_slp, Temp, Sali, PAR, mesh)

  use REcoM_declarations
  use REcoM_LocVar
  use recom_config
  use REcoM_GloVar
  use GASX
  use REcoM_ciso
  use g_clock
  use o_PARAM
  use g_PARSUP
  use g_rotate_grid
  use g_config
  use MOD_MESH
  use i_arrays
  use o_param
  use i_param
  use o_arrays
  use g_forcing_arrays
  use g_comm_auto
  use i_therm_param
  use g_comm
  use g_support

  Implicit none
  type(t_mesh), intent(in) , target         :: mesh

  Real(kind=8)                              :: Latr          
  Integer                                   :: Nn            ! Total number of nodes
  Real(kind=8),dimension(mesh%nl-1)	    :: zNodes	     ! Depth of nodes   zr(1:nzmax) = Z_3d_n(1:nzmax,n)
  Real(kind=8),dimension(mesh%nl-1,bgc_num) :: state         	
  Real(kind=8)                              :: SurfSW        ! [W/m2] ShortWave radiation at surface
  Real(kind=8)                              :: Loc_slp       ! [Pa] se-level pressure

  Real(kind=8),dimension(mesh%nl-1)         :: Temp          ! [degrees C] Ocean temperature
  real(kind=8),dimension(mesh%nl-1)         :: PAR

! Subroutine Depth

  Real(kind=8),dimension(mesh%nl)           :: zF                   ! [m] Depth of fluxes
  Real(kind=8),dimension(mesh%nl,4)         :: SinkVel              ! [m/day]
  Real(kind=8),dimension(mesh%nl-1)         :: thick                ! [m] Vertical distance between two nodes = Thickness 
  Real(kind=8),dimension(mesh%nl-1)         :: recipthick           ! [1/m] reciprocal of thick

! Subroutine CO2Flux /mocsy
  Real(kind=8)                         :: REcoM_DIC(1)            ! [mmol/m3] Conc of DIC in the surface water, used to calculate CO2 flux
  Real(kind=8)                         :: REcoM_Alk(1)            ! [mmol/m3] Conc of Alk in the surface water, used to calculate CO2 flux
  Real(kind=8)                         :: REcoM_Si(1)             ! [mol/m3] Conc of Si in the surface water, used to calculate CO2 flux
  Real(kind=8)                         :: REcoM_Phos(1)           ! [mol/m3] Conc of Phos in the surface water, used to calculate the CO2 flux
  Real(kind=8)                         :: Sali(1)                 ! Salinity of current surface layer
  Real(kind=8)                         :: Latd(1)                 ! latitude in degree
  Real(kind=8)                         :: Lond(1)                 ! longitude in degree
  Real(kind=8)                         :: REcoM_T(1)              ! temperature again, for mocsy minimum defined as -2
  Real(kind=8)                         :: REcoM_S(1)              ! temperature again, for mocsy minimum defined as 21
! atm pressure, now read in as forcing!!
  Real(kind=8)                         :: Patm(1)                 ! atmospheric pressure [atm]

! Subroutine o2flux /mocsy 
  Real(kind=8)                          :: ppo(1)                 ! atmospheric pressure, divided by 1 atm 
  Real(kind=8)                          :: REcoM_O2(1)            ! [mmol/m3] Conc of O2 in the surface water, used to calculate O2 flux

! Subroutine REcoM_sms
  Real(kind=8),dimension(mesh%nl-1,bgc_num) :: sms, aux                ! matrix that entail changes in tracer concentrations

!Diagnostics
  integer                              :: idiags,n,k
#include "../associate_mesh.h"
!  do n=1, myDim_nod2D
!     !!---- Number of vertical layers
!     nzmax=nlevels_nod2D(n)-1

    tiny_N   = tiny_chl/chl2N_max      ! 0.00001/ 3.15d0   Chl2N_max [mg CHL/mmol N] Maximum CHL a : N ratio = 0.3 gCHL gN^-1
    tiny_N_d = tiny_chl/chl2N_max_d    ! 0.00001/ 4.2d0
    tiny_C   = tiny_N  /NCmax          ! NCmax   = 0.2d0   [mmol N/mmol C] Maximum cell quota of nitrogen (N:C)
    tiny_C_d = tiny_N_d/NCmax_d        ! NCmax_d = 0.2d0 
    tiny_Si  = tiny_C_d/SiCmax         ! SiCmax = 0.8d0


  call Cobeta(mesh)        
  call Depth_calculations(n, Nn,SinkVel,zF,thick,recipthick, mesh)

!! ----- mocsy -------! 

!! convert from mmol/m3 to mol/m3
  REcoM_DIC  = max(tiny*1e-3, state(one,idic)*1e-3)
  REcoM_Alk  = max(tiny*1e-3, state(one,ialk)*1e-3)
  REcoM_Si   = max(tiny*1e-3, state(one,isi) *1e-3)

!! convert N to P with Redfield ratio
  REcoM_Phos = max(tiny*1e-3, state(one,idin)*1e-3) /16
!! minimum set to 2 degC: K1/K2 Lueker valid between 2degC-35degC and 19-43psu
  REcoM_T    = max(2.d0, Temp(1))
!! maximum set to 40 degC: K1/K2 Lueker valid between 2degC-35degC and 19-43psu
  REcoM_T    = min(REcoM_T, 40.d0) 
!! minimum set to 21: K1/K2 Lueker valid between 2degC-35degC and 19-43psu, else causes trouble in regions with S between 19 and 21 and ice conc above 97%
  REcoM_S    = max(21.d0, Sali(1)) 
!! maximum set to 43: K1/K2 Lueker valid between 2degC-35degC and 19-43psu, else causes trouble   REcoM_S    = min(REcoM_S, 43.d0)

!! convert from Pa to atm.
  Patm = Loc_slp/Pa2atm  

!! lon
  Lond=geo_coord_nod2D(1,n)/rad !! convert from rad to degree
!! lat
  Latr=geo_coord_nod2D(2,n)
  Latd=geo_coord_nod2D(2,n)/rad !! convert from rad to degree


! first calculate piston velocity kw660, which is an input to the flxco2 calculation:
! pistonvel already scaled for ice-free area:

  call pistonvel(ULoc, Loc_ice_conc, Nmocsy, kw660)
  call flxco2(co2flux, co2ex, dpco2surf,                                                    &
                  ph, pco2surf, fco2, co2, hco3, co3, OmegaA, OmegaC, BetaD, rhoSW, p, tempis,  &
                  REcoM_T, REcoM_S, REcoM_Alk, REcoM_DIC, REcoM_Si, REcoM_Phos, kw660, LocAtmCO2, Patm, thick(One), Nmocsy, Lond,Latd,     &
                  optCON='mol/m3',optT='Tpot   ',optP='m ',optB='u74',optK1K2='l  ',optKf='dg',optGAS='Pinsitu',optS='Sprc')

! changed optK1K2='l  ' to 'm10'
  if((co2flux(1)>1.e10) .or. (co2flux(1)<-1.e10)) then
!     co2flux(1)=0.0  
      print*,'ERROR: co2 flux !'
      print*,'pco2surf: ',pco2surf
      print*,'co2: ',co2
      print*,'rhoSW: ', rhoSW
      print*,'temp: ',REcoM_T
      print*,'tempis: ',tempis
      print*,'REcoM_S: ', REcoM_S
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
   dflux     = co2flux *1.e3 *SecondsPerDay  !* (1.d0 - Loc_ice_conc)
!   if(mype==0) write(*,*), 'dflux (mmol/m2/d) =',dflux
   co2flux_seaicemask = co2flux * 1.e3 !  [mmol/m2/s]  * (1.d0 - Loc_ice_conc)
!   if(mype==0) write(*,*), 'co2flux_seaicemask (mmol/m2/s) =',co2flux_seaicemask

! then oxygen
   ppo = Loc_slp/Pa2atm !1 !slp divided by 1 atm
   REcoM_O2 = max(tiny*1e-3,state(one,ioxy)*1e-3) ! convert from mmol/m3 to mol/m3 for mocsy

   call  o2flux(REcoM_T, REcoM_S, kw660, ppo, REcoM_O2, Nmocsy, o2ex)
     
   o2flux_seaicemask = o2ex * 1.e3 ! back to mmol here [mmol/m2/s] 

! Source-Minus-Sinks

if (recom_debug .and. mype==0) print *, achar(27)//'[36m'//'     --> REcoM_sms'//achar(27)//'[0m'


!  addtiny(1:nn,1) = state(1:nn,isi)
!  addtiny(1:nn,2) = state(1:nn,idetsi) 
!  addtiny(1:nn,3) = state(1:nn,idiasi)
!  addtiny(1:nn,4) = state(1:nn,idetz2si)

  call REcoM_sms(n, Nn, state, thick, recipthick, SurfSW, sms, Temp, SinkVel, zF, PAR, mesh)

!  addtiny(1:nn,1) = (state(1:nn,isi)           - aux(1:nn,isi))
!  addtiny(1:nn,2) = (state(1:nn,idetsi)        - aux(1:nn,idetsi))
!  addtiny(1:nn,3) = (state(1:nn,idiasi)        - aux(1:nn,idiasi)) 
!  addtiny(1:nn,4) = (state(1:nn,idetz2si)      - aux(1:nn,idetz2si))

!  aux=0.0d0
!  aux(1:nn,:)        = state(1:nn,:) + sms(1:nn,:)

  state(1:nn,:)      = max(tiny,state(1:nn,:) + sms(1:nn,:))
  state(1:nn,ipchl)  = max(tiny_chl,state(1:nn,ipchl))
  state(1:nn,iphyn)  = max(tiny_N,  state(1:nn,iphyn))
  state(1:nn,iphyc)  = max(tiny_C,  state(1:nn,iphyc))
  state(1:nn,idchl)  = max(tiny_chl,state(1:nn,idchl))
  state(1:nn,idian)  = max(tiny_N_d,state(1:nn,idian))
  state(1:nn,idiac)  = max(tiny_C_d,state(1:nn,idiac))
  state(1:nn,idiasi) = max(tiny_Si, state(1:nn,idiasi))

!  addtiny(1:nn,5) = (state(1:nn,isi)           - aux(1:nn,isi))
!  addtiny(1:nn,6) = (state(1:nn,idetsi)        - aux(1:nn,idetsi))
!  addtiny(1:nn,7) = (state(1:nn,idiasi)        - aux(1:nn,idiasi)) 
!  addtiny(1:nn,8) = (state(1:nn,idetz2si)      - aux(1:nn,idetz2si))

!  addtiny(1:nn,5) = state(1:nn,isi)
!  addtiny(1:nn,6) = state(1:nn,idetsi) 
!  addtiny(1:nn,7) = state(1:nn,idiasi)
!  addtiny(1:nn,8) = state(1:nn,idetz2si)

if (recom_debug .and. mype==0) print *, achar(27)//'[36m'//'     --> ciso after REcoM_Forcing'//achar(27)//'[0m'
  if (ciso) then
!    Calculcate isotopic fractionation of 13|14C, radioactive decay of 14C is calculated in oce_ale_tracer.F90
     kwco2  = kw660(1) * (660/scco2(REcoM_T(1)))**0.5 ! piston velocity (via mocsy)
     co2sat = co2flux(1) / (kwco2 + tiny) + co2(1)    ! saturation concentration of CO2 (via mocsy)

     do kj = one, nn
        recom_dic    = max(tiny*1e-3,state(kj,idic)*1e-3)    ! overriding REcoM_DIC from above!
        recom_dic_13 = max(tiny*1e-3,state(kj,idic_13)*1e-3)
        recom_dic_14 = max(tiny*1e-3,state(kj,idic_14)*1e-3)
        r_dic_13     = recom_dic_13(1) / recom_dic(1)
        r_dic_14     = recom_dic_14(1) / recom_dic(1)

!       Air-sea exchange and chemical speciation of CO2, fco3 is taken from mocsy
!       First calculate fractionation factors alpha_aq and alpha_dic. We ignore
!       vertical fco3 variations, because the effect on alpha_dic (and hence on alpha_p)
!       is much smaller than the effect of vertical temperature variations.
        call recom_ciso_airsea(temp(kj), co3(1), recom_dic(1))

        r_co2s_13 = alpha_aq_13 / alpha_dic_13 * r_dic_13  ! for biogenic fractionation further below
        r_co2s_14 = alpha_aq_14 / alpha_dic_14 * r_dic_14

        if (kj == one) then
           co2flux_13 = alpha_k_13 * alpha_aq_13 * kwco2 * (r_atm_13 * co2sat - r_dic_13 * co2(1) / alpha_dic_13)
           co2flux_14 = alpha_k_14 * alpha_aq_14 * kwco2 * (r_atm_14 * co2sat - r_dic_14 * co2(1) / alpha_dic_14)
           co2flux_seaicemask_13 = co2flux_13 * 1.e3 !  [mmol/m2/s]
           co2flux_seaicemask_14 = co2flux_14 * 1.e3 !  [mmol/m2/s]
        end if
!       Biogenic fractionation due to photosynthesis
!       Vertical variations of CO2* are estimated by a pressure correction (Weiss 1974, eq. 5)
!       involving the surface density provided by flxco2 (mind that znodes(k) < 0):
        prt = Patm(1) - 1.0e-5 * rhoSW(1) * 9.81 * znodes(kj)        ! total pressure at level k in atm
!       CO2* at depth, pressure correction following Weiss (1974, eq. 5) and mocsy
        co2s = co2(1) * exp((1 - prt) * vco2rgas / (Temp(kj) + C2K)) ! CO2* at depth
        call recom_ciso_photo(Cphot_z(kj), Cphot_dia_z(kj), co2s, r_co2s_13, r_co2s_14, rhoSW(1)) ! -> alpha_p
!       Now incorporate biogenic fractionation
!       iphyc_13|14 and idiac_13|14 are only used in REcoM_sms to calculate
!       DIC_13|14, DOC_13|14 and DetC_13|14
        state(kj,iphyc_13) = state(kj,iphyc) * r_co2s_13 / alpha_p_13
        state(kj,iphyc_14) = state(kj,iphyc) * r_co2s_14 / alpha_p_14
        state(kj,idiac_13) = state(kj,idiac) * r_co2s_13 / alpha_p_dia_13
        state(kj,idiac_14) = state(kj,idiac) * r_co2s_14 / alpha_p_dia_14
     end do
     state(:,iphyc_13)  = max(tiny_C,  state(:,iphyc_13))
     state(:,iphyc_14)  = max(tiny_C,  state(:,iphyc_14))
     state(:,idiac_13)  = max(tiny_C_d,state(:,idiac_13))
     state(:,idiac_14)  = max(tiny_C_d,state(:,idiac_14))
  end if

!-------------------------------------------------------------------------------
! Diagnostics
  if (Diags) then
	do idiags = one,8
	  LocDiags2D(idiags) = sum(diags3Dloc(1:nn,idiags) * thick(1:nn))
	end do
  end if

end subroutine REcoM_Forcing

