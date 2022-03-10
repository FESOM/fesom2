module command_line_options_module
! synopsis: read options passed to the main executable and trigger corresponding actions

  implicit none  
  public command_line_options
  private

  type :: command_line_options_type
  contains
    procedure, nopass :: parse
  end type
  type(command_line_options_type) command_line_options

contains

  subroutine parse()
    use info_module
    integer i
    character(len=:), allocatable :: arg
    integer arglength

    do i = 1, command_argument_count()    
      call get_command_argument(i, length=arglength)
      allocate(character(arglength) :: arg)
      call get_command_argument(i, value=arg)
      select case (arg)
      case('--smoketest')
        print '(g0)', 'smoketest'
      case('--info')
        print '(g0)', '# Definitions'
        call info%print_definitions()
      case default
        print *, 'unknown option: ', arg
        error stop
      end select
      deallocate(arg)
    end do

  end subroutine

end module
