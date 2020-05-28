!=======================================================================
!
! This submodule contain an ocean mixed layer implementation 
!
! Author: L. Zampieri ( lorenzo.zampieri@awi.de )
!
!=======================================================================

      submodule (icedrv_main) icedrv_ocean_mix_layer

      contains

      module subroutine ocn_mixed_layer_icepack(                &
                                         alvdr_ocn, swvdr,      &
                                         alidr_ocn, swidr,      &
                                         alvdf_ocn, swvdf,      &
                                         alidf_ocn, swidf,      &
                                         sst,       flwout_ocn, &
                                         fsens_ocn, shcoef,     &
                                         flat_ocn,  lhcoef,     &
                                         evap_ocn,  flw,        &
                                         delt,      delq,       &
                                         aice,      fhocn,      &
                                         fswthru,   hmix,       &
                                         Tf,        fresh,      &
                                         frain,     fsnow,      &
                                         fhocn_tot, fresh_tot,  &
                                         frzmlt)                

          use i_therm_param, only: emiss_wat
    
          real (kind=dbl_kind), intent(in) :: &
             alvdr_ocn , & ! visible, direct   (fraction)
             alidr_ocn , & ! near-ir, direct   (fraction)
             alvdf_ocn , & ! visible, diffuse  (fraction)
             alidf_ocn , & ! near-ir, diffuse  (fraction)
             swvdr     , & ! sw down, visible, direct  (W/m^2)
             swvdf     , & ! sw down, visible, diffuse (W/m^2)
             swidr     , & ! sw down, near IR, direct  (W/m^2)
             swidf     , & ! sw down, near IR, diffuse (W/m^2)
             flw       , & ! incoming longwave radiation (W/m^2)
             Tf        , & ! freezing temperature (C)
             hmix      , & ! mixed layer depth (m)
             delt      , & ! potential temperature difference   (K)
             delq      , & ! specific humidity difference   (kg/kg)
             shcoef    , & ! transfer coefficient for sensible heat
             lhcoef    , & ! transfer coefficient for latent heat
             fswthru   , & ! shortwave penetrating to ocean (W/m^2)
             aice      , & ! ice area fraction
             sst       , & ! sea surface temperature (C)
             frain     , & ! rainfall rate (kg/m^2/s)
             fsnow         ! snowfall rate (kg/m^2/s)
    
          real (kind=dbl_kind), intent(inout) :: &
             flwout_ocn, & ! outgoing longwave radiation (W/m^2)
             fsens_ocn , & ! sensible heat flux (W/m^2)
             flat_ocn  , & ! latent heat flux   (W/m^2)
             evap_ocn  , & ! evaporative water flux (kg/m^2/s)
             fhocn     , & ! net heat flux to ocean (W/m^2)
             fresh     , & ! fresh water flux to ocean (kg/m^2/s)
             frzmlt        ! freezing/melting potential (W/m^2)

          real (kind=dbl_kind), intent(out) :: &
             fhocn_tot , & ! net total heat flux to ocean (W/m^2)
             fresh_tot     ! fresh total water flux to ocean (kg/m^2/s)
    
          real (kind=dbl_kind), parameter :: &
             frzmlt_max = c1000   ! max magnitude of frzmlt (W/m^2)
    
          real (kind=dbl_kind) :: &
             TsfK ,    &  ! surface temperature (K)
             swabs,    &  ! surface absorbed shortwave heat flux (W/m^2)
             sst_n,    &  ! temporary new sst
             fhocn_n,  &
             fresh_n
    
          character(len=*),parameter :: subname='(icepack_ocn_mixed_layer)'
    
          ! shortwave radiative flux ! Visible is absorbed by clorophil afterwards
          swabs = (c1-alidr_ocn) * swidr  + (c1-alidf_ocn) * swidf + &
                  (c1-alvdr_ocn) * swvdr  + (c1-alvdf_ocn) * swvdf
    
          ! ocean surface temperature in Kelvin
          TsfK = sst + Tffresh
    
          ! longwave radiative flux
          ! Water emissivity added to be consistent
          ! with the standard FESOM2 version
          flwout_ocn = - emiss_wat * stefan_boltzmann * TsfK**4
    
          ! downward latent and sensible heat fluxes
          fsens_ocn =  shcoef * delt
          flat_ocn  =  lhcoef * delq
          evap_ocn  = -flat_ocn / Lvap
    
          ! Compute heat change due to exchange between ocean and atmosphere
    
          fhocn_tot = fhocn + fswthru                                  &  ! these are *aice already
                    + (fsens_ocn + flat_ocn + flwout_ocn + flw + swabs &
                    + Lfresh*fsnow) * (c1-aice) + max(c0,frzmlt)*aice
    
          fresh_tot = fresh + (-evap_ocn + frain + fsnow)*(c1-aice)

      end subroutine ocn_mixed_layer_icepack

      end submodule icedrv_ocean_mix_layer 
