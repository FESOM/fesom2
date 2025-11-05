!
! synopsis: Simple driver program for toy ocean model
!           Mimics fesom_main.F90 but uses the simplified toy_fesom_module
!
! This program runs a minimal ocean-only configuration suitable for
! idealized test cases like the Soufflet channel or double gyre experiments.
!
! Excluded features:
!   - Sea ice model (use_ice=.false.)
!   - Icebergs (use_icebergs=.false.)
!   - OASIS coupling (no __oasis)
!   - OpenIFS coupling (no __oifs)
!   - RECOM biogeochemistry (no __recom)
!   - Transient tracers (use_transit=.false.)
!   - Global tides (use_global_tides=.false.)
!   - Age tracers (use_age_tracer=.false.)
!   - Land ice coupling (use_landice_water=.false.)
!   - Icepack thermodynamics (no __icepack)
!
program fesom_toy_main
  use fesom_toy_module
  implicit none

  integer :: nsteps

  ! Initialize the toy ocean model
  call toy_fesom_init(nsteps)

  ! Run the time-stepping loop
  call toy_fesom_runloop(nsteps)

  ! Finalize and output statistics
  call toy_fesom_finalize

end program fesom_toy_main
