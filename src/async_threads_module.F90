#define __FILENAME__ "async_threads_module.F90"
module async_threads_module
  implicit none  
  public thread_type
  private
  

  abstract interface
    subroutine callback_interface(index_argument)
      integer, intent(in) :: index_argument
    end subroutine callback_interface
  end interface

  
  type thread_type
    private
    integer :: idx = 0
    procedure(callback_interface), nopass, pointer :: run_ptr => null()
    integer :: run_arg = 0
    logical :: with_real_threads = .true.
  contains
    procedure initialize
    procedure run
    procedure join
    procedure disable_async
  end type thread_type


  type wrapper_type
    class(thread_type), pointer :: ptr => null()
  end type wrapper_type
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
    
#ifdef DISABLE_MULTITHREADING
    call this%disable_async()
#else
    if(this%with_real_threads) call init_ccall(this%idx)
#endif
  end subroutine initialize


  subroutine run(this)
    class(thread_type) this
    ! EO args
    if(this%with_real_threads) then
#ifndef DISABLE_MULTITHREADING
      call begin_ccall(this%idx)
#endif
    else
      call async_threads_execute_fcall(this%idx)
    end if
  end subroutine run


  subroutine join(this)
    class(thread_type) this
    ! EO args    
#ifndef DISABLE_MULTITHREADING
    if(this%with_real_threads) call end_ccall(this%idx)
#endif
  end subroutine join
  
  
  subroutine disable_async(this)
    class(thread_type) this
    ! EO args
    this%with_real_threads = .false.  
  end subroutine disable_async


  subroutine async_threads_execute_fcall(thread_idx) bind (C, name="async_threads_execute_fcall")
    use iso_c_binding
    integer(c_int), intent(in), value :: thread_idx
    ! EO args
    class(thread_type), pointer :: t
    
    call assert(size(threads) >= thread_idx, __LINE__)
    t => threads(thread_idx)%ptr
    
    call t%run_ptr( t%run_arg )
  end subroutine async_threads_execute_fcall


  subroutine assert(val, line)
    logical, intent(in) :: val
    integer, intent(in) :: line
    ! EO args
    if(.NOT. val) then
      print *, "error in line ",line, __FILENAME__
      stop 1
    end if
  end subroutine assert

end module async_threads_module
