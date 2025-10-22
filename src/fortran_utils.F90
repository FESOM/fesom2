 ! synopsis: basic Fortran utilities, no MPI, dependencies only to INTRINSIC modules
module fortran_utils
  use mpi
  implicit none

contains


  function int_to_txt(val) result(txt)
    integer, intent(in) :: val
    character(:), allocatable :: txt
    ! EO parameters
    integer val_width

    if(val == 0) then
      val_width = 1
    else
      val_width = int(log10(real(val)))+1 ! does not work for val=0
    end if
    allocate(character(val_width) :: txt)
    write(txt,'(i0)') val
  end function int_to_txt


  function int_to_txt_pad(val, width) result(txt)
    integer, intent(in) :: val, width
    character(:), allocatable :: txt
    ! EO parameters
    integer w, val_width
    character(:), allocatable :: widthtxt

    if(val == 0) then
      val_width = 1
    else
      val_width = int(log10(real(val)))+1 ! does not work for val=0
    end if
    w = width
    if(w < val_width) w = val_width
    widthtxt = int_to_txt(w)
    allocate(character(w) :: txt)
    write(txt,'(i0.'//widthtxt//')') val
  end function int_to_txt_pad


  function mpirank_to_txt(mpicomm) result(txt)
    integer, intent(in) :: mpicomm
    character(:), allocatable :: txt
    ! EO parameters
    integer mype
    integer npes
    integer mpierr

    call MPI_Comm_Rank(mpicomm, mype, mpierr)
    call MPI_Comm_Size(mpicomm, npes, mpierr)
    txt = int_to_txt_pad(mype,int(log10(real(npes)))+1) ! pad to the width of the number of processes
  end function mpirank_to_txt

  subroutine mkdir(path)
    character(len=*), intent(in) :: path
    integer stat, cmdstat

    write(*,*) "Creating directory by calling mkdir -p : " // path
    call execute_command_line("mkdir -p " // trim(path), exitstat=stat,cmdstat=cmdstat)
    if(cmdstat /= 0)then
        write(*,'(a)')'<ERROR>:failed mkdir '
    endif
  end subroutine mkdir

end module fortran_utils
