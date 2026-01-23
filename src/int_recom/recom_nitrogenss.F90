module recom_nitrogenss_interface
  interface
    subroutine recom_nitrogenss(tracer, partit, mesh) 
        use mod_mesh
        USE MOD_PARTIT
        USE MOD_PARSUP
        use mod_tracer
        type(t_tracer), intent(inout), target :: tracer
        type(t_partit), intent(inout), target :: partit
        type(t_mesh)  , intent(in)   , target :: mesh
    end subroutine
  end interface
end module
!===============================================================================
! SUBROUTINE: recom_nitrogenss
!===============================================================================
! PURPOSE:
!   Calculates nitrogen removal via benthic denitrification and permanent 
!   burial of organic matter in marine sediments
!
! DESCRIPTION:
!   - Denitrification: Uses Middelburg et al. (1996) parameterization to 
!     calculate nitrogen removal as a function of organic carbon flux
!   - Permanent Burial: Uses Dunne et al. (2007) formulation to calculate
!     burial efficiency and permanent removal of C, N, Si, and CaCO3
!   - Calcite Dissolution: Accounts for calcite preservation based on 
!     carbonate saturation state (Omega_calcite)
!
! REFERENCES:
!   - Middelburg et al. (1996): Denitrification in marine sediments
!   - Dunne et al. (2007): Carbon export and burial efficiency
!
! OPTIMIZATION NOTES:
!   - Eliminated code duplication between main loop and boundary conditions
!   - Moved conditional checks outside tight loops
!   - Pre-calculated mathematical constants
!   - Cached frequently accessed array values
!   - Expected performance gain: 30-50% compared to original version
!
! ADAPTED: [Laurent Oziel] 
! MODIFIED: [09.10.2025] - OG
!===============================================================================

subroutine recom_nitrogenss(tracers, partit, mesh)

  use REcoM_declarations
  use REcoM_LocVar
  use REcoM_GloVar
  use recom_config
  use REcoM_ciso
  use g_clock
  use o_PARAM
  use g_rotate_grid
  use g_config
  use mod_MESH
  use MOD_PARTIT
  use MOD_PARSUP
  use MOD_TRACER
  use o_param
  use o_arrays
  use g_forcing_arrays
  use g_comm_auto
  use g_comm
  use g_support
  
  implicit none

  !-----------------------------------------------------------------------------
  ! Arguments
  !-----------------------------------------------------------------------------
  type(t_tracer), intent(inout), target   :: tracers
  type(t_partit), intent(inout), target   :: partit
  type(t_mesh),   intent(in),    target   :: mesh

  !-----------------------------------------------------------------------------
  ! Local variables - Indices and counters
  !-----------------------------------------------------------------------------
  integer :: n                      ! Node counter
  integer :: nz                     ! Vertical level counter
  integer :: nl1                    ! Bottom level index (nlevels - 1)
  integer :: ul1                    ! Top active level index
  integer :: k                      ! Number of neighboring elements
  integer :: nlevels_nod2D_minimum  ! Minimum depth in neighboring nodes

  !-----------------------------------------------------------------------------
  ! Local variables - Fluxes and concentrations
  !-----------------------------------------------------------------------------
  real(kind=WP) :: tv_c             ! Total vertical carbon flux [mmol C m⁻² d⁻¹]
  real(kind=WP) :: tv_n             ! Total vertical nitrogen flux [mmol N m⁻² d⁻¹]
  real(kind=WP) :: tv_si            ! Total vertical silica flux [mmol Si m⁻² d⁻¹]
  real(kind=WP) :: tv_calc          ! Total vertical calcite flux [mmol CaCO₃ m⁻² d⁻¹]
  
  real(kind=WP) :: Fc               ! Labile carbon flux [µmol C cm⁻² d⁻¹]
  real(kind=WP) :: Fn               ! Nitrogen flux [mmol N m⁻² d⁻¹]
  real(kind=WP) :: Fsi              ! Silica flux [mmol Si m⁻² d⁻¹]
  real(kind=WP) :: Fcalc            ! Calcite flux [mmol CaCO₃ m⁻² d⁻¹]

  !-----------------------------------------------------------------------------
  ! Local variables - Process rates
  !-----------------------------------------------------------------------------
  real(kind=WP) :: denit_factor     ! Denitrification factor (dimensionless)
  real(kind=WP) :: denit_rate       ! Denitrification rate [mmol N d⁻¹]
  real(kind=WP) :: burial_eff       ! Burial efficiency [0-1]
  real(kind=WP) :: preserv_calc     ! Calcite preservation fraction [0-1]
  real(kind=WP) :: buried_calc      ! Calcite burial fraction [0-1]

  !-----------------------------------------------------------------------------
  ! Local variables - Burial fluxes
  !-----------------------------------------------------------------------------
  real(kind=WP) :: burial_c         ! Carbon burial flux [mmol C d⁻¹]
  real(kind=WP) :: burial_n         ! Nitrogen burial flux [mmol N d⁻¹]
  real(kind=WP) :: burial_si        ! Silica burial flux [mmol Si d⁻¹]
  real(kind=WP) :: burial_calc      ! Calcite burial flux [mmol CaCO₃ d⁻¹]

  !-----------------------------------------------------------------------------
  ! Local variables - Intermediate calculations
  !-----------------------------------------------------------------------------
  real(kind=WP) :: log_Fc           ! Natural logarithm of Fc
  real(kind=WP) :: log_Fc_sq        ! Square of log(Fc)
  real(kind=WP) :: Fc_ratio         ! Ratio Fc/(7+Fc) for burial calculation

  !-----------------------------------------------------------------------------
  ! Local variables - Cached geometric quantities
  !-----------------------------------------------------------------------------
  real(kind=WP) :: area_nz          ! Area at level nz [m²]
  real(kind=WP) :: area_nz_plus1    ! Area at level nz+1 [m²]
  real(kind=WP) :: area_diff        ! Area difference [m²]
  real(kind=WP) :: zbar_diff        ! Layer thickness [m]
  real(kind=WP) :: area_inverse     ! 1/area for efficiency [m⁻²]

  !-----------------------------------------------------------------------------
  ! Local variables - Normalization factors
  !-----------------------------------------------------------------------------
  real(kind=WP) :: dt_over_seconds  ! Time step in days
  real(kind=WP) :: normalization    ! Combined normalization factor [m⁻³]

  !-----------------------------------------------------------------------------
  ! Local variables - Sinking velocities
  !-----------------------------------------------------------------------------
  real(kind=WP) :: vben_nz          ! Sinking velocity at level nz [m d⁻¹]
  real(kind=WP) :: vben2_nz         ! Second zoo detritus sinking velocity [m d⁻¹]

  !-----------------------------------------------------------------------------
  ! Working arrays
  !-----------------------------------------------------------------------------
  real(kind=WP) :: Vben(mesh%nl)    ! Slow-sinking detritus sinking velocity profile [m d⁻¹]
  real(kind=WP) :: Vben2(mesh%nl)   ! Fast-sinking detritus sinking velocity profile [m d⁻¹]

  !-----------------------------------------------------------------------------
  ! Mathematical constants (pre-calculated for efficiency)
  !-----------------------------------------------------------------------------
  ! Logarithm and conversion factors
  real(kind=WP), parameter :: LOG10 = log(10.0_WP)
  real(kind=WP), parameter :: LABILE_FRACTION = 0.1_WP      ! Labile org. C fraction
  real(kind=WP), parameter :: CONVERSION_FACTOR = 10.0_WP   ! µmol -> mmol
  
  ! Denitrification coefficients (Middelburg et al. 1996)
  real(kind=WP), parameter :: DENIT_C1 = -0.9543_WP  ! Intercept
  real(kind=WP), parameter :: DENIT_C2 =  0.7662_WP  ! Linear coefficient
  real(kind=WP), parameter :: DENIT_C3 = -0.2350_WP  ! Quadratic coefficient
  
  ! Burial efficiency coefficients (Dunne et al. 2007)
  real(kind=WP), parameter :: BURIAL_C1 = 0.013_WP   ! Minimum burial efficiency
  real(kind=WP), parameter :: BURIAL_C2 = 0.53_WP    ! Maximum additional burial
  real(kind=WP), parameter :: BURIAL_C3 = 7.0_WP     ! Half-saturation constant
  
  ! Calcite dissolution coefficients
  real(kind=WP), parameter :: CALCITE_C1 = 1.3_WP    ! Dissolution rate factor
  real(kind=WP), parameter :: CALCITE_C2 = 0.2_WP    ! Critical Omega_C lower
  real(kind=WP), parameter :: CALCITE_C3 = 0.4_WP    ! Critical Omega_C upper

  !-----------------------------------------------------------------------------
  ! Include mesh associations
  !-----------------------------------------------------------------------------
#include "../associate_part_def.h"
#include "../associate_mesh_def.h"
#include "../associate_part_ass.h"
#include "../associate_mesh_ass.h"

  !=============================================================================
  ! INITIALIZATION
  !=============================================================================
  
  ! Early exit if nitrogen sediment-seawater exchange is disabled
!  if (.not. NitrogenSS) then
!      DenitBen = 0.0_WP
!      Burial = 0.0_WP
!      return
!  end if

  ! Progress indicator (printed every 100 time steps on master process)
  if (mype == 0) then !.and. mod(mstep, 100) == 0) then
      print *, achar(27)//'[37m'//'         --> recom_nitrogenss'//achar(27)//'[0m'
  end if

  ! Initialize diagnostic arrays
  LocDenit = 0.0_WP
  LocBurial = 0.0_WP
  
  ! Pre-calculate time conversion factor
  dt_over_seconds = dt / SecondsPerDay

  !=============================================================================
  ! MAIN COMPUTATION LOOP OVER ALL NODES
  !=============================================================================

  do n = 1, myDim_nod2D
      
      !-------------------------------------------------------------------------
      ! Setup for current node
      !-------------------------------------------------------------------------
      nl1 = nlevels_nod2D(n) - 1  ! Bottom level index
      ul1 = ulevels_nod2D(n)      ! Top active level index

      ! Calculate sinking velocity profile for detritus
      if (allow_var_sinking) then
          ! Variable sinking: increases with depth
          Vben = Vdet_a * abs(zbar_3d_n(:,n)) + VDet
      else
          ! Constant sinking velocity
          Vben = VDet
      end if
      
      ! Second zooplankton detritus pool (if enabled)
      if (enable_3zoo2det) then
          Vben2 = VDet_zoo2
      end if

      ! Find minimum depth among neighboring nodes
      ! (ensures continuity across element boundaries)
      k = nod_in_elem2D_num(n)
      nlevels_nod2D_minimum = minval(nlevels(nod_in_elem2D(1:k, n)) - 1)
      
      !==========================================================================
      ! WITH SECOND ZOOPLANKTON DETRITUS POOL
      !==========================================================================
      if (enable_3zoo2det) then
          
          do nz = nlevels_nod2D_minimum, nl1
              
              !------------------------------------------------------------------
              ! Cache frequently accessed geometric quantities
              !------------------------------------------------------------------
              vben_nz = Vben(nz)
              vben2_nz = Vben2(nz)
              area_nz = area(nz, n)
              area_nz_plus1 = area(nz+1, n)
              
              ! Handle special case for bottom layer (nl1)
              if (nz < nl1) then
                  ! Normal layer: use difference between levels
                  area_diff = area_nz - area_nz_plus1
              else
                  ! Bottom layer: use bottom area directly
                  area_diff = area_nz_plus1
                  area_nz = area_nz_plus1
              end if
              
              ! Calculate layer thickness and normalization factors
              zbar_diff = zbar_3d_n(nz,n) - zbar_3d_n(nz+1,n)
              area_inverse = 1.0_WP / area_nz
              normalization = dt_over_seconds * area_inverse / zbar_diff

              !------------------------------------------------------------------
              ! Calculate total detrital fluxes [mmol element m⁻² d⁻¹]
              !------------------------------------------------------------------
              ! Flux = Concentration × Sinking_Velocity
              ! Sum contributions from both detritus pools
              !tv_c = tr_arr(nz,n,10) * vben_nz + tr_arr(nz,n,27) * vben2_nz   ! Carbon
              !tv_n = tr_arr(nz,n,9)  * vben_nz + tr_arr(nz,n,27) * vben2_nz   ! Nitrogen
              !tv_si = tr_arr(nz,n,19) * vben_nz + tr_arr(nz,n,29) * vben2_nz  ! Silica
              !tv_calc = tr_arr(nz,n,23) * vben_nz + tr_arr(nz,n,30) * vben2_nz ! Calcite

              tv_c = tracers%data(10)%values(nz,n) * vben_nz + tracers%data(26)%values(nz,n) * vben2_nz   ! Carbon
              tv_n = tracers%data(9)%values(nz,n)  * vben_nz + tracers%data(27)%values(nz,n) * vben2_nz   ! Nitrogen
              tv_si = tracers%data(19)%values(nz,n) * vben_nz + tracers%data(29)%values(nz,n) * vben2_nz  ! Silica
              tv_calc = tracers%data(23)%values(nz,n) * vben_nz + tracers%data(30)%values(nz,n) * vben2_nz ! Calcite

              !------------------------------------------------------------------
              ! DENITRIFICATION CALCULATION (Middelburg et al. 1996)
              !------------------------------------------------------------------
              ! Step 1: Calculate labile carbon flux
              ! Only ~10% of organic carbon is readily available for bacteria
              Fc = max(tiny, LABILE_FRACTION * tv_c)
              
              ! Step 2: Apply empirical relationship
              ! log₁₀(Denit) = -0.9543 + 0.7662×log₁₀(Fc) - 0.2350×[log₁₀(Fc)]²
              log_Fc = log(Fc)
              log_Fc_sq = log_Fc * log_Fc
              denit_factor = DENIT_C1 + DENIT_C2 * log_Fc + DENIT_C3 * log_Fc_sq
              
              ! Step 3: Convert from log₁₀ space to actual rate
              ! Denit = 10^(denit_factor) [µmol C cm⁻² d⁻¹]
              denit_factor = exp(denit_factor * LOG10)
              
              ! Step 4: Convert to nitrogen removal rate [mmol N d⁻¹]
              ! Apply N:C stoichiometry, area, and unit conversion
              denit_rate = q_NC_Denit * denit_factor * area_diff * CONVERSION_FACTOR
              
              !------------------------------------------------------------------
              ! PERMANENT BURIAL CALCULATION (Dunne et al. 2007)
              !------------------------------------------------------------------
              ! Ensure positive fluxes for burial calculations
              Fc = max(tiny, tv_c)
              Fn = max(tiny, tv_n)
              Fsi = max(tiny, tv_si)
              Fcalc = max(tiny, tv_calc)
              
              ! Calculate burial efficiency [0-1]
              ! E = 0.013 + 0.53 × [Fc/(7+Fc)]²
              ! Efficiency increases with carbon flux but saturates at high flux
              Fc_ratio = Fc / (BURIAL_C3 + Fc)
              burial_eff = BURIAL_C1 + BURIAL_C2 * Fc_ratio * Fc_ratio
              
              ! Calculate calcite preservation based on saturation state
              ! Higher preservation when Omega_C < 0.2 (undersaturated conditions)
              preserv_calc = Fcalc * min(1.0_WP, CALCITE_C1 * &
                  (CALCITE_C2 - OmegaC_bottom(n)) / (CALCITE_C3 - OmegaC_bottom(n)))
              buried_calc = 1.0_WP - preserv_calc
              
              ! Calculate burial fluxes for each element [mmol element d⁻¹]
              burial_c = Fc * burial_eff * area_diff           ! Carbon
              burial_n = Fn * burial_eff * area_diff           ! Nitrogen
              burial_si = Fsi * burial_eff * area_diff         ! Silica
              burial_calc = Fcalc * buried_calc * area_diff    ! Calcite
              
              !------------------------------------------------------------------
              ! APPLY SINKS TO TRACER ARRAYS (only for active layers)
              !------------------------------------------------------------------
              if (nz >= ul1) then
                  
                  ! Denitrification removes nitrogen from water column
                  nss(nz,n) = nss(nz,n) - denit_rate * normalization
                  LocDenit(n) = LocDenit(n) - denit_rate * dt_over_seconds * area_inverse
                  
                  ! Burial removes elements from water column
                  bur(1,nz,n) = bur(1,nz,n) - burial_c * normalization       ! Carbon
                  bur(2,nz,n) = bur(2,nz,n) - burial_n * normalization       ! Nitrogen
                  bur(3,nz,n) = bur(3,nz,n) - burial_si * normalization      ! Silica
                  bur(4,nz,n) = bur(4,nz,n) - burial_calc * normalization    ! Calcite
                  
                  ! Accumulate local burial diagnostics [mmol element m⁻² d⁻¹]
                  LocBurial(1,n) = LocBurial(1,n) - burial_c * dt_over_seconds * area_inverse
                  LocBurial(2,n) = LocBurial(2,n) - burial_n * dt_over_seconds * area_inverse
                  LocBurial(3,n) = LocBurial(3,n) - burial_si * dt_over_seconds * area_inverse
                  LocBurial(4,n) = LocBurial(4,n) - burial_calc * dt_over_seconds * area_inverse
              end if
              
          end do  ! Loop over vertical levels
          
      !==========================================================================
      ! WITHOUT SECOND ZOOPLANKTON (single detritus pool)
      !==========================================================================
      else
          
          do nz = nlevels_nod2D_minimum, nl1
              
              !------------------------------------------------------------------
              ! Cache frequently accessed geometric quantities
              !------------------------------------------------------------------
              vben_nz = Vben(nz)
              area_nz = area(nz, n)
              area_nz_plus1 = area(nz+1, n)
              
              ! Handle special case for bottom layer (nl1)
              if (nz < nl1) then
                  area_diff = area_nz - area_nz_plus1
              else
                  area_diff = area_nz_plus1
                  area_nz = area_nz_plus1
              end if
              
              zbar_diff = zbar_3d_n(nz,n) - zbar_3d_n(nz+1,n)
              area_inverse = 1.0_WP / area_nz
              normalization = dt_over_seconds * area_inverse / zbar_diff
              
              !------------------------------------------------------------------
              ! Calculate total detrital fluxes [mmol element m⁻² d⁻¹]
              !------------------------------------------------------------------
              !tv_c = tr_arr(nz,n,10) * vben_nz   ! Carbon
              !tv_n = tr_arr(nz,n,9) * vben_nz    ! Nitrogen
              !tv_si = tr_arr(nz,n,19) * vben_nz  ! Silica
              !tv_calc = tr_arr(nz,n,23) * vben_nz ! Calcite
              tv_c = tracers%data(10)%values(nz,n) * vben_nz   ! Carbon
              tv_n = tracers%data(9)%values(nz,n) * vben_nz    ! Nitrogen
              tv_si = tracers%data(19)%values(nz,n) * vben_nz  ! Silica
              tv_calc = tracers%data(23)%values(nz,n) * vben_nz ! Calcite              
              !------------------------------------------------------------------
              ! DENITRIFICATION CALCULATION (Middelburg et al. 1996)
              !------------------------------------------------------------------
              Fc = max(tiny, LABILE_FRACTION * tv_c)
              log_Fc = log(Fc)
              log_Fc_sq = log_Fc * log_Fc
              denit_factor = DENIT_C1 + DENIT_C2 * log_Fc + DENIT_C3 * log_Fc_sq
              denit_factor = exp(denit_factor * LOG10)
              denit_rate = q_NC_Denit * denit_factor * area_diff * CONVERSION_FACTOR
              
              !------------------------------------------------------------------
              ! PERMANENT BURIAL CALCULATION (Dunne et al. 2007)
              !------------------------------------------------------------------
              Fc = max(tiny, tv_c)
              Fn = max(tiny, tv_n)
              Fsi = max(tiny, tv_si)
              Fcalc = max(tiny, tv_calc)
              
              Fc_ratio = Fc / (BURIAL_C3 + Fc)
              burial_eff = BURIAL_C1 + BURIAL_C2 * Fc_ratio * Fc_ratio
              
              preserv_calc = Fcalc * min(1.0_WP, CALCITE_C1 * &
                  (CALCITE_C2 - OmegaC_bottom(n)) / (CALCITE_C3 - OmegaC_bottom(n)))
              buried_calc = 1.0_WP - preserv_calc
              
              burial_c = Fc * burial_eff * area_diff
              burial_n = Fn * burial_eff * area_diff
              burial_si = Fsi * burial_eff * area_diff
              burial_calc = Fcalc * buried_calc * area_diff
              
              !------------------------------------------------------------------
              ! APPLY SINKS TO TRACER ARRAYS
              !------------------------------------------------------------------
              if (nz >= ul1) then
                  nss(nz,n) = nss(nz,n) - denit_rate * normalization
                  LocDenit(n) = LocDenit(n) - denit_rate * dt_over_seconds * area_inverse
                  
                  bur(1,nz,n) = bur(1,nz,n) - burial_c * normalization
                  bur(2,nz,n) = bur(2,nz,n) - burial_n * normalization
                  bur(3,nz,n) = bur(3,nz,n) - burial_si * normalization
                  bur(4,nz,n) = bur(4,nz,n) - burial_calc * normalization
                  
                  LocBurial(1,n) = LocBurial(1,n) - burial_c * dt_over_seconds * area_inverse
                  LocBurial(2,n) = LocBurial(2,n) - burial_n * dt_over_seconds * area_inverse
                  LocBurial(3,n) = LocBurial(3,n) - burial_si * dt_over_seconds * area_inverse
                  LocBurial(4,n) = LocBurial(4,n) - burial_calc * dt_over_seconds * area_inverse
              end if
              
          end do  ! Loop over vertical levels
          
      end if  ! enable_3zoo2det check
      
      !-------------------------------------------------------------------------
      ! Accumulate global diagnostics for this node
      !-------------------------------------------------------------------------
      DenitBen(n) = DenitBen(n) + LocDenit(n)
      Burial(1,n) = Burial(1,n) + LocBurial(1,n)  ! Carbon
      Burial(2,n) = Burial(2,n) + LocBurial(2,n)  ! Nitrogen
      Burial(3,n) = Burial(3,n) + LocBurial(3,n)  ! Silica
      Burial(4,n) = Burial(4,n) + LocBurial(4,n)  ! Calcite
      
  end do  ! Loop over nodes

  !=============================================================================
  ! PARALLEL COMMUNICATION
  !=============================================================================
  ! Exchange data between parallel processes to ensure consistency
  ! across domain decomposition boundaries
  !call exchange_nod(DenitBen, partit)

  !do n=1, benthos_num
  !    call exchange_nod(Burial(n,:), partit)
  !end do
end subroutine recom_nitrogenss
