!=======================================================================
!
! Diagnostic information output during run
!
! author: Tony Craig

      module icedrv_system

      use icedrv_kinds
      use mod_parsup,       only: par_ex
      use icedrv_constants, only: ice_stderr
      use icepack_intfc,    only: icepack_warnings_flush, icepack_warnings_aborted
      use mod_partit
      implicit none

      public :: icedrv_system_abort, icedrv_system_init
   
      private
      type(t_partit), save, pointer  :: p_partit             ! a pointer to the mesh partitioning (has been accessed via "use g_parsup" in the original code)
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

      call par_ex(p_partit%MPI_COMM_FESOM, p_partit%mype, 1)
      stop

      end subroutine icedrv_system_abort

!=======================================================================
      subroutine icedrv_system_init(partit)
      implicit none
      type(t_partit), intent(inout), target :: partit

      p_partit => partit
      end subroutine icedrv_system_init

!=======================================================================


      end module icedrv_system

!=======================================================================

