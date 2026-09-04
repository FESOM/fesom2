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


  !> Calendar stamp of a model instant: YYYY-MM-DD.
  !>
  !> The instant is normalised first, because FESOM's clock represents midnight
  !> two ways and both have to produce one name. The end of day d
  !> (sec_of_day == 86400) is the start of day d+1, and the end of the last day
  !> of a year is 1 January of the next. That is the same normalisation
  !> clock_finish applies before writing fesom.clock and clock_init applies
  !> after reading it back, which is what lets a writer and a later reader
  !> derive the same name from their own clock state with no bookkeeping in
  !> between. Applying it to already-normalised input changes nothing, so it is
  !> safe to call on either side.
  !>
  !> days_in_year is 365 or 366 and carries the leap flag; the month lengths are
  !> the same table g_clock keeps, repeated here so this stays free of model
  !> dependencies and can be tested on its own.
  !>
  !> sec_of_day is an input even though the stamp has no time of day, because
  !> the normalisation needs it: it is what says whether the instant is the end
  !> of a day and therefore belongs to the next one. Two instants inside the
  !> same day get the same stamp.
  pure function calendar_stamp(year, day_of_year, sec_of_day, days_in_year) result(stamp)
    integer, intent(in) :: year, day_of_year, sec_of_day, days_in_year
    character(:), allocatable :: stamp
    ! EO parameters
    integer, parameter :: month_len(12,0:1) = reshape( &
         [31,28,31,30,31,30,31,31,30,31,30,31, &
          31,29,31,30,31,30,31,31,30,31,30,31], [12,2])
    integer :: yy, mm, dd, doy, leap, i, elapsed
    character(10) :: buf

    yy   = year
    doy  = day_of_year
    leap = max(0, min(1, days_in_year-365))

    ! end of a day is the start of the next one ...
    if(sec_of_day >= 86400) doy = doy+1
    ! ... and the end of the last day of a year is 1 January
    if(doy > 365+leap) then
      doy = 1
      yy = yy+1
      leap = 0 ! 1 January is 1 January in either kind of year
    end if

    ! day of year -> month and day in month
    mm = 1
    dd = doy
    elapsed = 0
    do i=1, 12
      if(doy > elapsed .and. doy <= elapsed+month_len(i,leap)) then
        mm = i
        dd = doy-elapsed
        exit
      end if
      elapsed = elapsed+month_len(i,leap)
    end do

    write(buf, '(i4.4,"-",i2.2,"-",i2.2)') yy, mm, dd
    stamp = buf
  end function calendar_stamp


end module fortran_utils
