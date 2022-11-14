    subroutine recom_atbox(mesh)
!     Simple 0-d box model to calculate the temporal evolution of atmospheric CO2.
!     Initially the box model was part of module recom_ciso. Now it can be run also 
!     without carbon isotopes (ciso==.false.)
!     mbutzin, 2021-07-08

!     Settings are copied from subroutine bio_fluxes,
!     some of the following modules may be unnecessary here
!     use REcoM_declarations
!     use REcoM_LocVar
      use REcoM_GloVar
      use recom_config
      use recom_ciso

      use mod_MESH
      USE g_CONFIG
      use o_ARRAYS
      use i_ARRAYS
      use g_comm_auto
      use g_forcing_arrays
      use g_PARSUP
      use g_support
      use i_therm_param

      
      implicit none
      integer                           :: n, elem, elnodes(3),n1
      real(kind=WP)                     :: total_co2flux,    &     ! (mol / s) 
                                           total_co2flux_13, &     ! (mol / s) carbon-13
                                           total_co2flux_14        ! (mol / s) radiocarbon
      real(kind=WP), parameter          :: mol_allatm = 1.7726e20  ! atmospheric inventory of all compounds (mol)
      type(t_mesh), intent(in) , target :: mesh
#include "../associate_mesh.h"

!     Globally integrated air-sea CO2 flux (mol / s)
      total_co2flux    = 0.
      call integrate_nod(0.001 * GloCO2flux_seaicemask   , total_CO2flux,    mesh)

!     Atmospheric carbon budget (mol)
!     mass of the dry atmosphere = 5.1352e18 kg (Trenberth & Smith 2005, doi:10.1175/JCLI-3299.1)
!     mean density of air = 0.02897 kg / mol (https://nssdc.gsfc.nasa.gov/planetary/factsheet/earthfact.html)
!     => total molecular inventory of the dry atmosphere: moles_atm = 1.7726e20 mol == constant.
!     mol_co2atm    = mol_co2atm    - total_co2flux    * dt 
!     Atmospheric mixing ratios in ppm
!     x_co2atm(1)    = mol_co2atm    / mol_allatm * 1.e6 ! ppm
      x_co2atm(1)    = x_co2atm(1)    - total_co2flux    / mol_allatm * dt * 1.e6
      x_co2atm       = x_co2atm(1)

      if (ciso) then
!       Consider 13CO2 (and maybe also 14CO2)

!       Globally integrated air-sea 13CO2 flux (mol / s)
        total_co2flux_13 = 0.
        call integrate_nod(0.001 * GloCO2flux_seaicemask_13, total_co2flux_13, mesh)

!       Atmospheric carbon-13 budget (mol)
!       mol_co2atm_13 = mol_co2atm_13 - total_co2flux_13 * dt 
!       Budget in terms of the 13C / 12C volume mixing ratio
!       x_co2atm_13(1) = mol_co2atm_13 / mol_allatm * 1.e6
        x_co2atm_13(1) = x_co2atm_13(1) - total_co2flux_13 / mol_allatm * dt * 1.e6
        x_co2atm_13    = x_co2atm_13(1)

        if (ciso_14) then
          total_co2flux_14   = 0.  ! globally integrated air-sea 14CO2 flux (mol / s)
          call integrate_nod(0.001 * GloCO2flux_seaicemask_14, total_co2flux_14, mesh)
!         Atmospheric radiocarbon budget in mol:
!         mol_co2atm_14 = mol_co2atm_14  + dt * (cosmic_14(1) - mol_co2atm_14 * lambda_14 - total_co2flux_14)
!                       = (mol_co2atm_14  + dt * (cosmic_14(1) - total_co2flux_14)) / (1 + lambda_14 * dt)
!         Budget in terms of the 14C / 12C volume mixing ratio
          x_co2atm_14(1) = (x_co2atm_14(1) + dt * (cosmic_14(1) - total_co2flux_14) / mol_allatm * 1.e6) / &
                           (1 + lambda_14 * dt)
          x_co2atm_14    = x_co2atm_14(1)

!         Adjust cosmogenic 14C production (mol / s) in spinup runs,
          r_atm_14        = x_co2atm_14(1) / x_co2atm(1)
!         r_atm_spinup_14 is calculated once-only in subroutine recom_init
          if (atbox_spinup .and. abs(r_atm_14 - r_atm_spinup_14) > 0.001) then
            cosmic_14(1) = cosmic_14(1) * (r_atm_spinup_14 / r_atm_14)
!           cosmic_14(1) = cosmic_14(1) * (1 + 0.01 * (r_atm_14_spinup / r_atm_14))
          end if
          cosmic_14 = cosmic_14(1)
        endif
      end if

      return
    end subroutine recom_atbox

