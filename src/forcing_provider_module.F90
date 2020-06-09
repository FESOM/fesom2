module forcing_provider_module
  use forcing_lookahead_reader_module
  implicit none
  public forcing_provider
  private
  
  type forcing_provider_type
    private
    type(forcing_lookahead_reader_type), allocatable :: all_readers(:)
    contains
    procedure, public :: get_forcingdata
  end type
  type(forcing_provider_type), save :: forcing_provider ! handle to singleton instance  
  
  contains

  subroutine get_forcingdata(this, varcount, varindex, filepath, fileyear, varname, time_index, forcingdata)
    class(forcing_provider_type), intent(inout) :: this
    integer, intent(in) :: varcount ! not used here, but we want to have the same arguments as for the async forcing provider
    integer, intent(in) :: varindex ! todo: remove this arg and just use a hashmap for varname
    character(len=*), intent(in) :: filepath
    integer, intent(in) :: fileyear
    character(len=*), intent(in) :: varname
    integer, intent(in) :: time_index
    real(4), intent(out) :: forcingdata(:,:)
    ! EO args
    type(forcing_lookahead_reader_type), allocatable :: tmparr(:)
    type(forcing_lookahead_reader_type) new_reader
    
    ! init our all_readers array if not already done
    if(.not. allocated(this%all_readers)) then
      allocate(this%all_readers(varindex))
    end if
        
    if(size(this%all_readers) < varindex) then
      allocate( tmparr(varindex) )
      tmparr(1:size(this%all_readers)) = this%all_readers
      deallocate(this%all_readers)
      call move_alloc(tmparr, this%all_readers)      
    end if
    
    if(.not. this%all_readers(varindex)%is_initialized() ) then ! reader has never been initialized
      call this%all_readers(varindex)%initialize(filepath, fileyear, varname)
    else if(fileyear /= this%all_readers(varindex)%fileyear()) then
      ! close our reader and create a new one
      call this%all_readers(varindex)%finalize()
      this%all_readers(varindex) = new_reader
      call this%all_readers(varindex)%initialize(filepath, fileyear, varname)
    end if

    call this%all_readers(varindex)%yield_data(time_index, forcingdata)
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
