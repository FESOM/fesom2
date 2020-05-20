!=======================================================================
!
! Primary state variables in various configurations
! Note: other state variables are at the end of this...
! The primary state variable names are:
!-------------------------------------------------------------------
! for each category   aggregated over     units
!                       categories
!-------------------------------------------------------------------
! aicen(i,n)         aice(i)           ---
! vicen(i,n)         vice(i)           m
! vsnon(i,n)         vsno(i)           m
! trcrn(i,it,n)      trcr(i,it)        
!
! Area is dimensionless because aice is the fractional area
! (normalized so that the sum over all categories, including open
! water, is 1.0).  That is why vice/vsno have units of m instead of m^3.
!
! Variable names follow these rules:
!
! (1) For 3D variables (indices i,n), write 'ice' or 'sno' or
!     'sfc' and put an 'n' at the end.
! (2) For 2D variables (indices i) aggregated over all categories,
!     write 'ice' or 'sno' or 'sfc' without the 'n'.
! (3) For 2D variables (indices i) associated with an individual
!     category, write 'i' or 's' instead of 'ice' or 'sno' and put an 'n'
!     at the end: e.g. hin, hsn.  These are not declared here
!     but in individual modules (e.g., ice_therm_vertical).
!
! authors C. M. Bitz, UW
!         Elizabeth C. Hunke and William H. Lipscomb, LANL

      module icedrv_state

      use icedrv_kinds
      use icedrv_domain_size, only: nx, ncat, max_ntrcr

      implicit none
      private

      !-----------------------------------------------------------------
      ! state of the ice aggregated over all categories
      !-----------------------------------------------------------------

      real (kind=dbl_kind), dimension(nx), &
         public :: &
         aice  , & ! concentration of ice
         vice  , & ! volume per unit area of ice          (m)
         vsno      ! volume per unit area of snow         (m)

      real (kind=dbl_kind), &
         dimension(nx,max_ntrcr), public :: &
         trcr      ! ice tracers
                   ! 1: surface temperature of ice/snow   (C)

      !-----------------------------------------------------------------
      ! state of the ice for each category
      !-----------------------------------------------------------------

      real (kind=dbl_kind), dimension (nx), &
         public:: &
         aice0     ! concentration of open water

      real (kind=dbl_kind), &
         dimension (nx,ncat), public :: &
         aicen , & ! concentration of ice
         vicen , & ! volume per unit area of ice          (m)
         vsnon     ! volume per unit area of snow         (m)

      real (kind=dbl_kind), public, &
         dimension (nx,max_ntrcr,ncat) :: &
         trcrn     ! tracers
                   ! 1: surface temperature of ice/snow   (C)

      !-----------------------------------------------------------------
      ! tracer infrastructure arrays
      !-----------------------------------------------------------------

      integer (kind=int_kind), dimension (max_ntrcr), public :: &
         trcr_depend   ! = 0 for ice area tracers
                       ! = 1 for ice volume tracers
                       ! = 2 for snow volume tracers

      integer (kind=int_kind), dimension (max_ntrcr), public :: &
         n_trcr_strata ! number of underlying tracer layers

      integer (kind=int_kind), dimension (max_ntrcr,2), public :: &
         nt_strata     ! indices of underlying tracer layers

      real (kind=dbl_kind), dimension (max_ntrcr,3), public :: &
         trcr_base     ! = 0 or 1 depending on tracer dependency
                       ! argument 2:  (1) aice, (2) vice, (3) vsno

      !-----------------------------------------------------------------
      ! dynamic variables closely related to the state of the ice
      !-----------------------------------------------------------------

      real (kind=dbl_kind), dimension(nx), &
         public :: &
         uvel     , & ! x-component of velocity (m/s)
         vvel     , & ! y-component of velocity (m/s)
         divu     , & ! strain rate I component, velocity divergence (1/s)
         shear    , & ! strain rate II component (1/s)
         strength     ! ice strength (N/m)

      !-----------------------------------------------------------------
      ! ice state at start of time step, saved for later in the step 
      !-----------------------------------------------------------------

      real (kind=dbl_kind), dimension(nx), &
         public :: &
         aice_init       ! initial concentration of ice, for diagnostics

      real (kind=dbl_kind), &
         dimension(nx,ncat), public :: &
         aicen_init  , & ! initial ice concentration, for linear ITD
         vicen_init  , & ! initial ice volume (m), for linear ITD
         vsnon_init      ! initial snow volume (m), for aerosol 

!=======================================================================

      end module icedrv_state

!=======================================================================
