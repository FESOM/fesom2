module mo_roctx
  use iso_c_binding
  implicit none

  integer, private, parameter :: ROCTX_MAX_LEN = 256

  interface roctxRangePushA
     subroutine roctxRangePushA(name) bind(C, name='roctxRangePushA')
       use iso_c_binding
       character(kind=c_char) :: name(256)
     end subroutine roctxRangePushA
  end interface roctxRangePushA

  interface roctxRangePop
     subroutine roctxRangePop() bind(C, name='roctxRangePop')
     end subroutine roctxRangePop
  end interface roctxRangePop

  interface roctracer_start
     subroutine roctracer_start() bind(C, name='roctracer_start')
     end subroutine roctracer_start
  end interface roctracer_start

  interface roctracer_stop
     subroutine roctracer_stop() bind(C, name='roctracer_stop')
     end subroutine roctracer_stop
  end interface roctracer_stop

contains

  subroutine roctxStartRange(name)
    character(kind=c_char,len=*) :: name
    character :: c_name(ROCTX_MAX_LEN)
    integer:: i, str_len

    str_len = len(trim(name))    
    do i = 1, len(trim(name))
       c_name(i) = name(i:i)
    end do
    c_name(str_len+1) = C_NULL_CHAR

    call roctxRangePushA(c_name)

  end subroutine roctxStartRange

end module mo_roctx
