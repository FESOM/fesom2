!=======================================================================
!
! Diagnostic information output during run
!
! author: Tony Craig

      module icedrv_system

      use icedrv_kinds
      use g_parsup,         only: par_ex          
      use icedrv_constants, only: ice_stderr
      use icepack_intfc,    only: icepack_warnings_flush, icepack_warnings_aborted

      implicit none

      public :: icedrv_system_abort
   
      private

!=======================================================================

      contains

!=======================================================================

      subroutine icedrv_system_abort(icell, istep, string, file, line)

      integer (kind=int_kind), intent(in), optional :: &
         icell       , & ! indices of grid cell where model aborts
         istep       , & ! time step number
         line            ! line number

      character (len=*), intent(in), optional :: string, file

      ! local variables

      character(len=*), parameter :: subname='(icedrv_system_abort)'

      write(ice_stderr,*) ' '

      call icepack_warnings_flush(ice_stderr)

      write(ice_stderr,*) ' '
      write(ice_stderr,*) subname,' ABORTED: '
      if (present(file))   write (ice_stderr,*) subname,' called from ', trim(file)
      if (present(line))   write (ice_stderr,*) subname,' line number',  line
      if (present(istep))  write (ice_stderr,*) subname,' istep =',      istep
      if (present(string)) write (ice_stderr,*) subname,' string =',     trim(string)
       
      ! Stop FESOM2

      call par_ex(1)
      stop

      end subroutine icedrv_system_abort

!=======================================================================

      end module icedrv_system

!=======================================================================

