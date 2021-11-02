! synopsis: save any derived types we initialize
!           so they can be reused after fesom_init
module fvom_types_storage_module

end module

! synopsis: main FESOM program split into 3 parts
!           this way FESOM can e.g. be used as a library with an external time loop driver
!           used with IFS-FESOM
module fvom_module
  implicit none
  public fesom_init, fesom_runloop, fesom_finalize
  private

contains
 
  subroutine fesom_init(nsteps)
    integer, intent(out) :: nsteps
    ! EO parameters
    
  end subroutine


  subroutine fesom_runloop(nsteps)
    integer, intent(in) :: nsteps 
    ! EO parameters

  end subroutine


  subroutine fesom_finalize()

  end subroutine

end module
