module io_MEANDATA

  use g_config
  use o_mesh
  use g_parsup
  use g_clock
  use g_comm_auto
  use o_ARRAYS
  implicit none
#include "netcdf.inc"
  private
  public :: def_stream3D, output
!
!--------------------------------------------------------------------------------------------
!
  type Meandata
    private
    integer                                            :: ndim
    integer                                            :: lcsize(2)
    integer                                            :: glsize(2)
    real(kind=8),  public, allocatable, dimension(:,:) :: local_values
    integer                                            :: addcounter=0
    real(kind=WP), pointer                             :: ptr2(:), ptr3(:,:)

    character(500)                                     :: filename
    character(100)                                     :: name
    character(500)                                     :: description
    character(100)                                     :: units
    character(100)                                     :: dimname(2)
    integer                                            :: ncid
    integer                                            :: rec_count=0
    integer                                            :: recID, tID
    integer                                            :: dimID(2), varID
    integer                                            :: freq=1
    character                                          :: freq_unit='m'
    integer                                            :: error_status(1000), error_count
    logical                                            :: is_in_use=.false.
  end type
!
!--------------------------------------------------------------------------------------------
!
  type(Meandata), save, allocatable, target :: io_stream(:)
  integer, save                             :: io_NSTREAMS=0
  real(kind=WP)                             :: ctime !current time in seconds from the beginning of the year
!
!--------------------------------------------------------------------------------------------
! generic interface was required to associate variables of unknown rank with the pointers of the same rank
! this allows for automatic streaming of associated variables into the netcdf file
  INTERFACE def_stream
            MODULE PROCEDURE def_stream2D, def_stream3D
  END INTERFACE
!
!--------------------------------------------------------------------------------------------
!
  contains
!
!--------------------------------------------------------------------------------------------
!
subroutine ini_mean_io
  implicit none
!					global size           local size      varname     varname long      unit        array           writeout unit (day,mon,year...) 
! call def_stream((/nl-1, elem2D/), (/nl-1, myDim_elem2D/), 'u', 'horizontal velocity', 'm/s',     uv(1,:,:),     1, 'm')  
!3D
  call def_stream((/nl-1, elem2D/), (/nl-1, myDim_elem2D/), 'u', 'horizontal velocity', 'm/s',     uv(1,:,:),     1, 'm')
  call def_stream((/nl-1, elem2D/), (/nl-1, myDim_elem2D/), 'v', 'meridional velocity', 'm/s',     uv(2,:,:),     1, 'm')
  call def_stream((/nl-1, nod2D/), (/nl-1, myDim_nod2D+eDim_nod2D/), 'temp', 'temperature', 'C',   tr_arr(:,:,1), 1, 'm')
  call def_stream((/nl-1, nod2D/), (/nl-1, myDim_nod2D+eDim_nod2D/), 'salt', 'salinity',    'psu', tr_arr(:,:,2), 1, 'm')
!2D
!   call def_stream(nod2D, myDim_nod2D+eDim_nod2D, 'ssh', 'sea surface elevation',   'm', eta_n,                                1, 'd')
!   call def_stream(nod2D, myDim_nod2D+eDim_nod2D, 'sst', 'sea surface temperature', 'C', tr_arr(1,1:myDim_nod2D+eDim_nod2D,1), 1, 'd')
!   call def_stream(nod2D, myDim_nod2D+eDim_nod2D, 'sss', 'sea surface salinity',  'psu', tr_arr(1,1:myDim_nod2D+eDim_nod2D,2), 1, 'd')
  call def_stream(nod2D, myDim_nod2D+eDim_nod2D, 'ssh', 'sea surface elevation',   'm', eta_n,                                1, 'm')
  call def_stream(nod2D, myDim_nod2D+eDim_nod2D, 'sst', 'sea surface temperature', 'C', tr_arr(1,1:myDim_nod2D+eDim_nod2D,1), 1, 'm')
  call def_stream(nod2D, myDim_nod2D+eDim_nod2D, 'sss', 'sea surface salinity',  'psu', tr_arr(1,1:myDim_nod2D+eDim_nod2D,2), 1, 'm')
end subroutine ini_mean_io
!
!--------------------------------------------------------------------------------------------
!
function get_dimname(n) result(s)
  implicit none
  integer       :: n
  character(50) :: s
  if (n==nod2D) then
     s='nod2'
  elseif (n==elem2D) then
     s='elem'
  elseif (n==nl) then
     s='nz'
  elseif (n==nl-1) then
     s='nz1'
  else
     s='unknown'
     if (mype==0) write(*,*) 'WARNING: unknown dimension in mean I/O with zise of ', n
  end if
  end function
!
!--------------------------------------------------------------------------------------------
!
subroutine create_new_file(entry)
  implicit none
  integer                       :: c, j
  character(2000)               :: att_text

  type(Meandata), intent(inout) :: entry
  ! Serial output implemented so far
  if (mype/=0) return
  c=1
  entry%error_status=0
  ! create an ocean output file
  write(*,*) 'initializing I/O file for ', trim(entry%name)

  if (restart_offset==32) then
     entry%error_status(c) = nf_create(entry%filename, nf_clobber, entry%ncid);                      c=c+1
  else
     entry%error_status(c) = nf_create(entry%filename, IOR(NF_CLOBBER,NF_64BIT_OFFSET), entry%ncid); c=c+1
  end if

  do j=1, entry%ndim
!___Create mesh related dimentions__________________________________________
     entry%error_status(c) = nf_def_dim(entry%ncid, entry%dimname(j), entry%glsize(j), entry%dimID(j)); c=c+1
  end do
!___Create time related dimentions__________________________________________
  entry%error_status(c) = nf_def_dim(entry%ncid, 'time', NF_UNLIMITED, entry%recID);                     c=c+1
!___Define the time and iteration variables_________________________________
  entry%error_status(c) = nf_def_var(entry%ncid, 'time', NF_DOUBLE, 1, entry%recID, entry%tID); c=c+1

  att_text='time'
  entry%error_status(c) = nf_put_att_text(entry%ncid, entry%tID, 'long_name', len_trim(att_text), trim(att_text)); c=c+1
  write(att_text, '(a14,I4.4,a1,I2.2,a1,I2.2,a6)'), 'seconds since ', yearold, '-', 1, '-', 1, ' 0:0:0'
  entry%error_status(c) = nf_put_att_text(entry%ncid, entry%tID, 'units', len_trim(att_text), trim(att_text)); c=c+1

  entry%error_status(c) = nf_def_var(entry%ncid, trim(entry%name), NF_DOUBLE, entry%ndim+1, &
                       (/entry%dimid(1:entry%ndim), entry%recID/), entry%varID); c=c+1
  entry%error_status(c)=nf_put_att_text(entry%ncid, entry%varID, 'description', len_trim(entry%description), entry%description); c=c+1
  entry%error_status(c)=nf_put_att_text(entry%ncid, entry%varID, 'units',       len_trim(entry%units),       entry%units);       c=c+1
  entry%error_status(c)=nf_close(entry%ncid); c=c+1
  entry%error_count=c-1
end subroutine create_new_file
!
!--------------------------------------------------------------------------------------------
!
subroutine assoc_ids(entry)
  implicit none

  type(Meandata), intent(inout) :: entry
  integer                       :: c, j, k
  real(kind=8)                  :: rtime !timestamp of the record
  ! Serial output implemented so far
  if (mype/=0) return
  c=1
  entry%error_status=0
  ! open existing netcdf file
  write(*,*) 'associating mean I/O file ', trim(entry%filename)

  entry%error_status(c) = nf_open(entry%filename, nf_nowrite, entry%ncid);

  if (entry%error_status(c) .ne. nf_noerr) then
     call create_new_file(entry) ! error status counter will be reset
     c=entry%error_count+1
     entry%error_status(c) = nf_open(entry%filename, nf_nowrite, entry%ncid); c=c+1
  end if

  do j=1, entry%ndim
!___Create mesh related dimentions__________________________________________
     entry%error_status(c) = nf_inq_dimid(entry%ncid, entry%dimname(j), entry%dimID(j)); c=c+1
  end do
!___Associate time related dimentions_______________________________________
  entry%error_status(c) = nf_inq_dimid (entry%ncid, 'time', entry%recID);          c=c+1
  entry%error_status(c) = nf_inq_dimlen(entry%ncid, entry%recID, entry%rec_count); c=c+1
!___Associate the time and iteration variables______________________________
  entry%error_status(c) = nf_inq_varid(entry%ncid, 'time', entry%tID); c=c+1
!___if the time rtime at the rec_count is larger than ctime we look for the closest record with the 
! timestamp less than ctime
  do k=entry%rec_count, 1, -1
     entry%error_status(c)=nf_get_vara_double(entry%ncid, entry%tID, k, 1, rtime, 1);
     if (ctime > rtime) then
        entry%rec_count=k+1
!       write(*,*) 'I/O '//trim(entry%name)//' : current record = ', entry%rec_count, '; ', entry%rec_count, ' records in the file;'
        exit ! a proper rec_count detected, exit the loop
     end if
     if (k==1) then
        write(*,*) 'I/O '//trim(entry%name)//' WARNING: the existing output file will be overwritten'//'; ', entry%rec_count, ' records in the file;'
        entry%rec_count=1
        exit ! no appropriate rec_count detected
     end if
  end do
  c=c+1 ! check will be made only for the last nf_get_vara_double

  entry%rec_count=max(entry%rec_count, 1)
!___Associate physical variables____________________________________________
  entry%error_status(c) = nf_inq_varid(entry%ncid, entry%name, entry%varID); c=c+1
  entry%error_status(c)=nf_close(entry%ncid); c=c+1
  entry%error_count=c-1
  write(*,*) trim(entry%name)//': current mean I/O counter = ', entry%rec_count
end subroutine assoc_ids
!
!--------------------------------------------------------------------------------------------
!
subroutine write_mean(entry)
  implicit none
  type(Meandata), intent(inout) :: entry
  real(kind=8), allocatable     :: aux1(:), aux2(:,:) 
  integer                       :: i, size1, size2
  integer                       :: c
  ! Serial output implemented so far
  if (mype==0) then
     c=1
     write(*,*) 'writing mean record for ', trim(entry%name), '; rec. count = ', entry%rec_count
     entry%error_status(c)=nf_open(entry%filename, nf_write, entry%ncid); c=c+1
     entry%error_status(c)=nf_put_vara_double(entry%ncid, entry%Tid, entry%rec_count, 1, ctime, 1); c=c+1
  end if
!_______writing 2D fields________________________________________________
     if (entry%ndim==1) then
        size1=entry%glsize(1)
        if (mype==0) allocate(aux1(size1))
        if (size1==nod2D)  call gather_nod (entry%local_values(1:entry%lcsize(1),1), aux1)
        if (size1==elem2D) call gather_elem(entry%local_values(1:entry%lcsize(1),1), aux1)
        if (mype==0) then
           entry%error_status(c)=nf_put_vara_double(entry%ncid, entry%varID, (/1, entry%rec_count/), (/size1, 1/), aux1, 1); c=c+1
        end if
        if (mype==0) deallocate(aux1)
!_______writing 3D fields________________________________________________
     elseif (entry%ndim==2) then
        size1=entry%glsize(1)
        size2=entry%glsize(2)
        if (mype==0) allocate(aux2(size1, size2))
        if (size1==nod2D  .or. size2==nod2D)  call gather_nod (entry%local_values(1:entry%lcsize(1),1:entry%lcsize(2)), aux2)
        if (size1==elem2D .or. size2==elem2D) call gather_elem(entry%local_values(1:entry%lcsize(1),1:entry%lcsize(2)), aux2)
        if (mype==0) then
           entry%error_status(c)=nf_put_vara_double(entry%ncid, entry%varID, (/1, 1, entry%rec_count/), (/size1, size2, 1/), aux2, 2); c=c+1
        end if
        if (mype==0) deallocate(aux2)
     else
        if (mype==0) write(*,*) 'not supported shape of array in mean I/O file'
           call par_ex
           stop
     end if

  if (mype==0) entry%error_count=c-1
  call was_error(entry)
  if (mype==0) entry%error_status(1)=nf_close(entry%ncid);
  entry%error_count=1
  call was_error(entry)
end subroutine write_mean
!
!--------------------------------------------------------------------------------------------
!
subroutine update_means
  implicit none
  type(Meandata), pointer :: entry
  integer                 :: n

  do n=1, io_NSTREAMS
     entry=>io_stream(n)
     if (entry%ndim==1) then 
        entry%local_values(:, 1) = entry%local_values(:, 1)+entry%ptr2(:)
     elseif (entry%ndim==2) then 
        entry%local_values = entry%local_values+entry%ptr3(:,:)
     else
        if (mype==0) write(*,*) 'not supported size in update_means'
        call par_ex
        stop
     end if
     entry%addcounter=entry%addcounter+1
  end do
end subroutine update_means
!
!--------------------------------------------------------------------------------------------
!
subroutine output(istep)
  implicit none

  integer       :: istep
  logical, save :: lfirst=.true.
  integer       :: mpierr
  integer       :: n
  logical       :: do_output
  type(Meandata), pointer :: entry

  ctime=timeold+(dayold-1.)*86400
  if (lfirst) call ini_mean_io

  call update_means

  do n=1, io_NSTREAMS
     entry=>io_stream(n)
     !check whether restart will be written
     do_output=.false.

     if (entry%freq_unit.eq.'y') then
        call annual_event(do_output)

     else if (entry%freq_unit == 'm') then 
        call monthly_event(do_output) 

     else if (entry%freq_unit == 'd') then
        call daily_event(do_output, entry%freq)

     else if (entry%freq_unit == 'h') then
        call daily_event(do_output, entry%freq)

     else if (entry%freq_unit == 's') then
        call daily_event(do_output, istep, entry%freq)

     else
        write(*,*) 'You did not specify a supported outputflag.'
        write(*,*) 'The program will stop to give you opportunity to do it.'
        call par_ex(1)
        stop
     endif

     if (.not. do_output) cycle
     entry%filename=trim(ResultPath)//trim(entry%name)//'.'//trim(runid)//'.'//cyearnew//'.nc'
     call assoc_ids(entry)
     entry%local_values = entry%local_values /float(entry%addcounter) ! compute_means
     call write_mean(entry)
     entry%local_values = 0. ! clean_meanarrays
     entry%addcounter   = 0  ! clean_meanarrays
  end do
  lfirst=.false.
end subroutine output
!
!--------------------------------------------------------------------------------------------
!
subroutine def_stream3D(glsize, lcsize, name, description, units, data, freq, freq_unit)
  implicit none
  integer                              :: glsize(2), lcsize(2)
  character(len=*),      intent(in)    :: name, description, units
  real(kind=WP), target, intent(inout) :: data(:,:)
  integer,               intent(in)    :: freq
  character,             intent(in)    :: freq_unit
  type(Meandata),        allocatable   :: tmparr(:)
  type(Meandata),        pointer       :: entry
   ! add this instance to io_stream array
  if ( .not. allocated(io_stream)) then
    allocate(io_stream(1))
  else
    allocate(tmparr(size(io_stream)+1))
    tmparr(1:size(io_stream)) = io_stream
    deallocate(io_stream)
    call move_alloc(tmparr, io_stream)
  end if
  entry=>io_stream(size(io_stream))
  entry%ptr3 => data
  allocate(entry%local_values(lcsize(1), lcsize(2)))
  entry%ndim=2
  entry%lcsize=lcsize
  entry%glsize=glsize
  entry%name = name
  entry%description = description
  entry%units = units

  entry%dimname(1)=get_dimname(glsize(1))
  entry%dimname(2)=get_dimname(glsize(2))
  entry%freq=freq
  entry%freq_unit=freq_unit
  ! clean_meanarrays
  entry%local_values = 0. 
  entry%addcounter   = 0
  entry%is_in_use=.true.
  io_NSTREAMS=io_NSTREAMS+1
end subroutine def_stream3D
!
!--------------------------------------------------------------------------------------------
!
subroutine def_stream2D(glsize, lcsize, name, description, units, data, freq, freq_unit)
  implicit none
  integer                              :: glsize, lcsize
  character(len=*),      intent(in)    :: name, description, units
  real(kind=WP), target, intent(inout) :: data(:)
  integer,               intent(in)    :: freq
  character,             intent(in)    :: freq_unit
  type(Meandata),        allocatable   :: tmparr(:)
  type(Meandata),        pointer       :: entry
   ! add this instance to io_stream array
  if ( .not. allocated(io_stream)) then
    allocate(io_stream(1))
  else
    allocate(tmparr(size(io_stream)+1))
    tmparr(1:size(io_stream)) = io_stream
    deallocate(io_stream)
    call move_alloc(tmparr, io_stream)
  end if
  entry=>io_stream(size(io_stream))
  entry%ptr2 => data
  allocate(entry%local_values(lcsize, 1))
  entry%ndim=1
  entry%lcsize=(/lcsize, 1/)
  entry%glsize=(/glsize, 1/)
  entry%name = name
  entry%description = description
  entry%units = units

  entry%dimname(1)=get_dimname(glsize)
  entry%dimname(2)='unknown'
  entry%freq=freq
  entry%freq_unit=freq_unit
  ! clean_meanarrays
  entry%local_values = 0. 
  entry%addcounter   = 0
  entry%is_in_use=.true.
  io_NSTREAMS=io_NSTREAMS+1
end subroutine def_stream2D
!
!--------------------------------------------------------------------------------------------
!
subroutine was_error(entry)
  implicit none
  type(Meandata), intent(inout) :: entry
  integer                       :: k, status, ierror

  call MPI_BCast(entry%error_count, 1,  MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)
  call MPI_BCast(entry%error_status(1), entry%error_count, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)

  do k=1, entry%error_count
     status=entry%error_status(k)
     if (status .ne. nf_noerr) then
        if (mype==0) call handle_err(status)
        call par_ex
        stop
     end if
  end do
end subroutine was_error

end module io_MEANDATA

