!===============================================================================
! For variables saved locally for each column and then used in REcoM
!-------------------------------------------------------------------------------
Module REcoM_locVar

  Real(kind=8),allocatable,dimension(:) :: LocBenthos ! Storing the values for benthos in current watercolumn: N,C,Si and Calc
  Real(kind=8) :: Hplus                     ! [mol/kg] Concentrations of H-plus ions in the surface node
  Real(kind=8) :: pCO2surf(1)                  ! [uatm] Partial pressure of CO2 in surface layer at current 2D node	
  Real(kind=8) :: dflux(1)                     ! [mmol/m2/day] Flux of CO2 into the ocean
  Real(kind=8) :: oflux(1)                     ! [mmol/m2/day] Flux of O2 into the ocean
  Real(kind=8) :: o2ex(1)                     ! [mmol/m2/s] Flux of O2 into the ocean
  Real(kind=8) :: ULoc(1)                      ! Wind strength above current 2D node, change array size if used with mocsy input vector longer than one
  Real(kind=8) :: dpCO2surf(1)              ! [uatm] difference of oceanic pCO2 minus atmospheric pCO2

! mocsy output -----------------------------------------------------------------------------------------------------------------------------
  Real(kind=8) :: co2flux(1)                   ! air-to-sea flux of CO2 [mol/(m^2 * s)]
  Real(kind=8) :: co2ex(1)                     ! time rate of change of surface CO2 due to gas exchange [mol/(m^3 * s)]
  Real(kind=8) :: dpco2(1)                     ! difference of oceanic pCO2 minus atmospheric pCO2 [uatm]
  Real(kind=8) :: ph(1)                        ! pH on total scale
  Real(kind=8) :: pco2(1)                      ! oceanic partial pressure of CO2 (uatm)
  Real(kind=8) :: fco2(1)                      ! oceanic fugacity of CO2 (uatm)
  Real(kind=8) :: co2(1)                       ! aqueous CO2 concentration [mol/m^3]
  Real(kind=8) :: hco3(1)                      ! bicarbonate (HCO3-) concentration [mol/m^3]
  Real(kind=8) :: co3(1)                       ! carbonate (CO3--) concentration [mol/m^3]
  Real(kind=8) :: OmegaA(1)                    ! Omega for aragonite, i.e., the aragonite saturation state
  Real(kind=8) :: OmegaC(1)                    ! Omega for calcite, i.e., the   calcite saturation state
  Real(kind=8) :: BetaD(1)                     ! BetaD = Revelle factor   dpCO2/pCO2 / dDIC/DIC
  Real(kind=8) :: rhoSW(1)                     ! rhoSW  = in-situ density of seawater; rhoSW = f(s, t, p)
  Real(kind=8) :: p(1)                         ! pressure [decibars]; p = f(depth, latitude) if computed from depth [m] OR p = depth if [db]
  Real(kind=8) :: tempis(1)                    ! in-situ temperature [degrees C]
  Real(kind=8) :: dpos(1)                      ! depth converted to positive values, needed in the mocsy routine
  Real(kind=8) :: kw660(1)                     ! gas transfer velocity (piston velocity) for CO2 [m/s] 
  Real(kind=8) :: K0(1)                        ! CO2 solubility
  Real(kind=8) :: co2flux_seaicemask(1)        ! air-to-sea flux of CO2 [mmol/m2/s]
  Real(kind=8) :: o2flux_seaicemask(1)         ! air-to-sea flux of CO2 [mmol/m2/s]

! mocsy output entire depth range ------------------------------------------------------------------------------------------------------------  ! NEW MOCSY
  Real(kind=8) :: ph_depth(1)                  ! NEW MOCSY pH on total scale
  Real(kind=8) :: pco2_depth(1)                ! NEW MOCSY oceanic partial pressure of CO2 (uatm)
  Real(kind=8) :: fco2_depth(1)                ! NEW MOCSY oceanic fugacity of CO2 (uatm)
  Real(kind=8) :: co2_depth(1)                 ! NEW MOCSY aqueous CO2 concentration [mol/m^3]
  Real(kind=8) :: hco3_depth(1)                ! NEW MOCSY bicarbonate (HCO3-) concentration [mol/m^3]
  Real(kind=8) :: co3_depth(1)                 ! NEW MOCSY carbonate (CO3--) concentration [mol/m^3]
  Real(kind=8) :: OmegaA_depth(1)              ! NEW MOCSY Omega for aragonite, i.e., the aragonite saturation state
  Real(kind=8) :: OmegaC_depth(1)              ! NEW MOCSY Omega for calcite, i.e., the   calcite saturation state
  Real(kind=8) :: BetaD_depth(1)               ! NEW MOCSY BetaD = Revelle factor   dpCO2/pCO2 / dDIC/DIC
  Real(kind=8) :: kspc_depth(1)                ! NEW DISS  stoichiometric solubility product of calcite (mol^2/kg^2)
  Real(kind=8) :: rhoSW_depth(1)               ! NEW MOCSY rhoSW  = in-situ density of seawater; rhoSW = f(s, t, p)
  Real(kind=8) :: p_depth(1)                   ! NEW MOCSY pressure [decibars]; p = f(depth, latitude) if computed from depth [m] OR p = depth if [db]
  Real(kind=8) :: tempis_depth(1)              ! NEW MOCSY in-situ temperature [degrees C]

  Integer      :: logfile_outfreq_7            ! NEW MOCSY helper value to calculate the timesteps for the carbonate system (every 7th day)
  Integer      :: logfile_outfreq_30           ! NEW MOCSY helper value to calculate the timesteps for the carbonate system (every 30th day)

!-------------------------------------------------------------------------------

  Real(kind=8) :: bt, dic_molal, talk_molal    ! Common block: Species
  Real(kind=8) :: k1, k2, kw, kb, ff           ! Common block: Equilibrium_constants
  Real(kind=8) :: FeDust                       ! [umol/m2/s]
  Real(kind=8) :: NDust                        ! [mmol/m2/s]
  Real(kind=8) :: Loc_ice_conc(1)              ! Used to calculate flux of DIC in REcoM 0 -> 1
  Real(kind=8) :: LocAtmCO2(1)                 ! [uatm]
  Real(kind=8) :: LocDiags2D(12)               ! (changed it from 8 to 12)
!  Real(kind=8) :: LocDenit                    ! BALL
  Real(kind=8) :: LocRiverDIN, LocRiverDON, LocRiverDOC, LocRiverDSi, LocRiverDIC, LocRiverAlk

  Real(kind=8) :: res_zoo2_a, res_zoo2_f
  Real(kind=8) :: grazingFluxcarbonzoo2     ! grazingfluxcarbon 
  Real(kind=8) :: grazingFluxcarbon_mes        ! Zoo3

  Real(kind=8) :: PICPOCtemp                   ! (added to make the calcification dependent on the temperature, after Krumhardt et al. 2017/2019)
  Real(kind=8) :: PICPOCCO2                    ! (to make calcification dependent on CO2)
  Real(kind=8) :: PICPOCN                      ! (to make calcification dependent on N-limitation)
  Real(kind=8) :: calc_prod_final              ! (added to make the calcification dependent on nutrients (N, Fe), after Krumhardt et al. 2017/2019)
  Integer      :: currentCO2year

end module REcoM_LocVar
