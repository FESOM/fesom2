#if defined(__yac)
module cpl_yac_driver

  use mo_yac_finterface, only: yac_finit, yac_fdef_calendar, &
       YAC_PROLEPTIC_GREGORIAN, yac_fdef_comp, yac_fget_comp_comm, &
       yac_ffinalize
  implicit none
  save

  character(len=*), PARAMETER   :: comp_name = "fesom2"
  integer :: comp_id

contains

  subroutine cpl_yac_init( localCommunicator )
    implicit none

    integer, intent(out) :: localCommunicator

#ifdef VERBOSE
    print *, '================================================='
    print *, 'cpl_yac_init : coupler initialization for YAC'
    print *, '*************************************************'
#endif /* VERBOSE */

    CALL yac_finit()
    CALL yac_fdef_calendar(YAC_PROLEPTIC_GREGORIAN)
    CALL yac_fdef_comp(comp_name, comp_id)

    CALL yac_fget_comp_comm(comp_id, localCommunicator)

  end subroutine cpl_yac_init

  subroutine cpl_yac_finalize ()
    implicit none
#ifdef VERBOSE
    print *, '================================================='
    print *, 'cpl_yac_finalize : coupler finalization for YAC'
    print *, '*************************************************'
#endif /* VERBOSE */
    CALL yac_ffinalize()
  end subroutine cpl_yac_finalize

end module cpl_yac_driver
#endif
