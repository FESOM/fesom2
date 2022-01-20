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
  
  
  ! using EXECUTE_COMMAND_LINE to call mkdir sometimes fail (EXECUTE_COMMAND_LINE is forked as an new process, which may be the problem)
  ! try to use the C mkdir as an alternative
  subroutine mkdir(path)
    use iso_c_binding
    character(len=*), intent(in) :: path
    ! EO parameters
    integer result
    character(:,kind=C_CHAR), allocatable :: pathcopy

    interface
      function mkdir_c(path, mode) bind(c,name="mkdir")
        use iso_c_binding
        integer(c_int) mkdir_c
        character(kind=c_char,len=1) path(*)
        integer(c_int), value :: mode
      end function
    end interface
      
    pathcopy = path ! we need to pass an array of c_char to the C funcktion (this is not a correct type conversion, but Fortran characters seem to be of the same kind as c_char)
    ! result is 0 if the dir has been created from this call, otherwise -1
    ! the mode will not exactly be what we pass here, as it is subtracted by the umask bits (and possibly more)
    result = mkdir_c(pathcopy//C_NULL_CHAR, int(o'777', c_int))
  end subroutine

end module
