 ! synopsis: basic Fortran utilities, no MPI, dependencies only to INTRINSIC modules
module fortran_utils
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
  end function


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
  end function

end module
