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
