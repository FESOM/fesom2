!=============================================================================!
!
!                 Finite Volume Sea-ice Ocean Model
!
!=============================================================================!
!                      The main driving routine
!=============================================================================!    

program main
  use fesom_module

  integer nsteps

  call fesom_init(nsteps)
  call fesom_runloop(nsteps)
  call fesom_finalize

end program
