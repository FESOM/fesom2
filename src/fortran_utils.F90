 ! synopsis: basic Fortran utilities, no MPI, dependencies only to INTRINSIC modules
module fortran_utils
  implicit none

contains


  function positiveint_to_txt(val) result(txt)
    integer, intent(in) :: val
    character(int(log10(real(val)))+1) :: txt ! does not work for val=0
    ! EO parameters
    write(txt,'(i0)') val
  end function


  function positiveint_to_txt_pad(val, width) result(txt) ! for val=0 width must be > 0
    integer, intent(in) :: val, width
    character(:), allocatable :: txt
    ! EO parameters
    integer w, val_width
    character(:), allocatable :: widthtxt

    val_width = int(log10(real(val)))+1
    w = width
    if(w < val_width) w = val_width
    widthtxt = positiveint_to_txt(w)
    allocate(character(w) :: txt)
    write(txt,'(i0.'//widthtxt//')') val
  end function

end module
