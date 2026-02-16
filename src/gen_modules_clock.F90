module g_clock
  !combining RT and Lars version
  !
  use g_config
  use iso_fortran_env, only: error_unit
  use mpi
  implicit none
  save
  real(kind=WP)            :: timeold, timenew     !time in a day, unit: sec
  integer                  :: dayold, daynew       !day in a year
  integer                  :: yearold, yearnew     !year before and after time step
  integer                  :: yearstart            !year when simulation started
  integer                  :: month, day_in_month  !month and day in a month
  integer                  :: fleapyear            !1 fleapyear, 0 not 
  integer                  :: ndpyr                !number of days in yearnew 
  integer                  :: num_day_in_month(0:1,12)
  character(4)             :: cyearold, cyearnew   !year as character string      
  character(2)             :: cmonth               !month as character string      
  data num_day_in_month(0,:) /31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/
  data num_day_in_month(1,:) /31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/


contains
  !
  !--------------------------------------------------------------------------------
  !
  subroutine clock

    implicit none
    integer          :: i
    real(kind=WP)    :: aux1, aux2
    !
    timeold=timenew 
    dayold=daynew
    yearold=yearnew

    ! update time
    timenew=timenew+dt          
     
    ! update day
    if (timenew>86400._WP) then  !assumed that time step is less than one day!
       daynew=daynew+1
       timenew=timenew-86400._WP
    endif

    ! update year
    if (daynew>ndpyr) then
       daynew=1
       yearnew=yearnew+1
       call check_fleapyr(yearnew, fleapyear)
       ndpyr=365+fleapyear
       write(cyearold,'(i4)') yearold
       write(cyearnew,'(i4)') yearnew
    endif

    ! find month and dayinmonth at new time step
    aux1=0
    do i=1,12
       aux2=aux1+num_day_in_month(fleapyear,i)
       if(daynew>aux1 .and. daynew<=aux2) then
          month=i
          write(cmonth, '(I2.2)') month
          day_in_month=daynew-aux1
          exit
       end if
       aux1=aux2
    end do
       
  end subroutine clock
  !
  !--------------------------------------------------------------------------------
  !
  subroutine clock_init(partit)
    USE MOD_PARTIT
    use g_config
    use mod_transit, only: ti_transit, ti_start_transit
    implicit none
    type(t_partit), intent(in), target    :: partit
    integer                               :: i, daystart
    real(kind=WP)                         :: aux1, aux2, timestart
    integer                               :: ierr
    integer                               :: file_unit
    character(512)                        :: errmsg
 
    ! init clock for this run - read clock file FIRST
    open(newunit=file_unit, file=trim(RestartInPath)//trim(runid)//'.clock', action='read', &
        status='old', iostat=ierr, iomsg=errmsg)
    if (ierr /= 0) then
      write (unit=error_unit, fmt='(3A)') &
        '### error: can not open file ', trim(RestartInPath)//trim(runid)//'.clock', &
        ', error: ' // trim(errmsg)
      call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
    end if
    read(unit=file_unit, fmt=*) timeold, dayold, yearold
    read(unit=file_unit, fmt=*) timenew, daynew, yearnew
    close(unit=file_unit)
    if(daynew==0) daynew=1
    
    ! the model initialized at - set AFTER reading clock file
    timestart=timenew
    daystart=daynew
    yearstart=yearnew
    
    ! check if this is a restart or not
    ! For initial run: clock file has same values on both lines (old=new)
    ! For restart: clock file has different old/new values
    if(yearnew==yearold .and. daynew==dayold .and. timenew==timeold) then
       r_restart=.false.
       yearold=yearnew-1 !required for checking if create new output files
    else
       r_restart=.true.
    end if
!   For simulations with transient tracer input data
    if (use_transit) ti_transit = yearnew - yearstart + ti_start_transit


    ! year as character string 
    write(cyearold,'(i4)') yearold
    write(cyearnew,'(i4)') yearnew

    ! if restart model at beginning of a day, set timenew to be zero
    if (timenew==86400._WP) then  
       timenew=0.0_WP
       daynew=daynew+1
    endif

    ! set timeold to be timenew, ready for initializing forcing fields,
    ! yearold should not be updated here, which is requird to open input files.
    ! timeold=timenew 
    ! dayold=daynew
    
    ! check fleap year
    call check_fleapyr(yearnew, fleapyear)
    ndpyr=365+fleapyear

    ! find month and dayinmonth at the new time step
    aux1=0
    do i=1,12
       aux2=aux1+num_day_in_month(fleapyear,i)
       if(daynew>aux1 .and. daynew<=aux2) then
          month=i
          day_in_month=daynew-aux1
          exit
       end if
       aux1=aux2
    end do

    if(partit%mype==0) then
        if(r_restart) then
            write(*,*)
            print *, achar(27)//'[31m'    //'____________________________________________________________'//achar(27)//'[0m'
            print *, achar(27)//'[5;7;31m'//' --> THIS IS A RESTART RUN !!!                              '//achar(27)//'[0m'
            write(*,"(A, F8.2, I4, I5)") '     > clock restarted at time:', timenew, daynew, yearnew
            write(*,"(A, I5)") '     > yearstart for annual_event:', yearstart
            write(*,*)
        else
            write(*,*)
            print *, achar(27)//'[32m'  //'____________________________________________________________'//achar(27)//'[0m'
            print *, achar(27)//'[7;32m'//' --> THIS IS A INITIALISATION RUN !!!                       '//achar(27)//'[0m'
            write(*,"(A, F8.2, I4, I5)")'     > clock initialized at time:', timenew, daynew, yearnew
            write(*,"(A, I5)") '     > yearstart for annual_event:', yearstart
            write(*,*)
        end if
    end if
  
  end subroutine clock_init
  !
  !-------------------------------------------------------------------------------
  !
  subroutine clock_finish
    implicit none
    !
    real(kind=WP)            :: dum_timenew     !time in a day, unit: sec
    integer                  :: dum_daynew       !day in a year
    integer                  :: dum_yearnew     !year before and after time step
    integer                               :: ierr
    integer                               :: file_unit
    character(512)                        :: errmsg
    
    dum_timenew = timenew
    dum_daynew  = daynew
    dum_yearnew = yearnew
    if ((dum_daynew==ndpyr) .and. (dum_timenew==86400._WP)) then
       dum_timenew=0.0_WP
       dum_daynew=1
       dum_yearnew=yearold+1
    endif

    open(newunit=file_unit, file=trim(RestartOutPath)//trim(runid)//'.clock', action='write', &
        status='unknown', iostat=ierr, iomsg=errmsg)
    if (ierr /= 0) then
      write (unit=error_unit, fmt='(3A)') &
        '### error: can not open file ', trim(RestartOutPath)//trim(runid)//'.clock', &
        ', error: ' // trim(errmsg)
      call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
    end if
    write(unit=file_unit, fmt=*) timeold, dayold, yearold
    write(unit=file_unit, fmt=*) dum_timenew, dum_daynew, dum_yearnew
    close(unit=file_unit)
  end subroutine clock_finish
  !
  !----------------------------------------------------------------------------
  !
  subroutine clock_newyear
    implicit none
    !
    if ((daynew>=ndpyr).and.(timenew==86400._WP)) then
       timenew=0.0_WP
       daynew=1
       yearnew=yearold+1
       write(cyearnew,'(i4)') yearnew
    endif
  end subroutine clock_newyear
  !
  !----------------------------------------------------------------------------
  !
  subroutine check_fleapyr(year, flag)
    implicit none
    integer, intent(in) :: year      
    integer, intent(out):: flag

    flag=0

    if(.not.include_fleapyear) return
    call is_fleapyr(year, flag)
  end subroutine check_fleapyr

  subroutine is_fleapyr(year, flag)
    implicit none
    integer, intent(in) :: year      
    integer, intent(out):: flag
    flag=0
    if ((mod(year,4)==0.and.mod(year,100)/=0) .or. mod(year,400)==0) then
       flag=1
    endif
  end subroutine is_fleapyr
  !
  !----------------------------------------------------------------------------
  !
end module g_clock