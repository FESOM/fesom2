module mod_hecuba
  use o_PARAM, only : WP
  use, intrinsic :: iso_fortran_env, only: real64, real32
  implicit none
  public:: output_hecuba
  private
contains

subroutine output_hecuba(istep, mesh)
  use g_clock
  use mod_mesh
  use g_PARSUP
  use io_gather_module
  use o_arrays
#if defined (__icepack)
  use icedrv_main,    only: init_io_icepack
#endif

  implicit none

  integer       :: istep
  logical, save :: lfirst=.true.
  type(t_mesh), intent(in) , target :: mesh
  integer :: glsize, lcsize ! =mesh%nod2D , lcsize= mesh%myDim_nod2d   
  real(real64)                  :: rtime !timestamp of the record
  character(len=100) :: root_dir, filename, dirnames, tmp
  character(len=200) :: variable_path
  character(len=5) :: mype_key
  integer  :: fileunit
  real(kind=WP),allocatable :: local_data2d(:)
  real(kind=WP),allocatable :: global_data2d(:)

  integer, save :: counter=0
  
  ! write 2d variable sst
  write(filename,'(A I0 A I0)')  "data_", istep, "_sst_pe", mype
  open(newunit=fileunit,file=trim(filename), form="unformatted", access='stream', action='write')
  write(fileunit) tr_arr(1,1:myDim_nod2d,1) 
  flush(fileunit)
  close(fileunit)

  ! write a 3d variable
  write(filename,'(A I0 A I0)')  "data_", istep, "_temp_", mype
  open(newunit=fileunit,file=trim(filename), form="unformatted", access='stream', action='write')
  write(fileunit) tr_arr(:,1:myDim_nod2d,1) 
  flush(fileunit)
  close(fileunit)

end subroutine output_hecuba 
   
end module mod_hecuba
