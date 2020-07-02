module async_threads_module
  implicit none  
  public thread_type
  private
  

  abstract interface
    subroutine callback_interface(index_argument)
      integer, intent(in) :: index_argument
    end subroutine
  end interface

  
  type thread_type
    private
    integer :: idx = 0
    procedure(callback_interface), nopass, pointer :: run_ptr => null()
    integer :: run_arg = 0
  contains
    procedure initialize
    procedure begin
    procedure wait
  end type


  type wrapper_type
    type(thread_type), pointer :: ptr => null()
  end type
  type(wrapper_type), save, allocatable :: threads(:)


contains


  subroutine initialize(this, thread_callback_procedure, procedure_argument)
    class(thread_type), target :: this
    procedure(callback_interface) thread_callback_procedure
    integer procedure_argument
    ! EO args
    type(wrapper_type), allocatable :: tmparr(:)

    if( .not. allocated(threads)) then
      allocate(threads(1))
    else
      allocate( tmparr(size(threads)+1) )
      tmparr(1:size(threads)) = threads
      deallocate(threads)
      call move_alloc(tmparr, threads)
    end if
    threads(size(threads))%ptr => this

    this%idx = size(threads)
    this%run_ptr => thread_callback_procedure
    this%run_arg = procedure_argument
  end subroutine


  subroutine begin(this)
    class(thread_type) this
    ! EO args
    call begin_ccall(this%idx)
  end subroutine


  subroutine wait(this)
    class(thread_type) this
    ! EO args
    call end_ccall(this%idx)
  end subroutine


  subroutine execute_fcall(thread_idx) bind (C, name="execute_fcall")
    use iso_c_binding
    integer(c_int), intent(in), value :: thread_idx
    ! EO args
    type(thread_type), pointer :: t
    
    call assert(size(threads) >= thread_idx, __LINE__)
    t => threads(thread_idx)%ptr
    
    call t%run_ptr( t%run_arg )
  end subroutine


  subroutine assert(val, line)
    logical, intent(in) :: val
    integer, intent(in) :: line
    ! EO args
    if(.NOT. val) then
      print *, "error in line ",line, __FILE__
      stop 1
    end if
  end subroutine

end module
