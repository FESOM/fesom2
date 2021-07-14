MODULE io_RESTART
  use restart_file_group_module
  use g_clock
  use o_arrays
  use i_arrays
  use g_cvmix_tke
  implicit none
  public :: restart, finalize_restart
  private

  integer,       save       :: globalstep=0 ! todo: remove this from module scope as it will mess things up if we use async read/write from the same process
  real(kind=WP)             :: ctime !current time in seconds from the beginning of the year

  type(restart_file_group), save :: oce_files
  type(restart_file_group), save :: ice_files
  character(:), allocatable, save :: oce_path
  character(:), allocatable, save :: ice_path

  character(:), allocatable, save :: raw_restart_dirpath


  contains
!
!--------------------------------------------------------------------------------------------
! ini_ocean_io initializes ocean_file datatype which contains information of all variables need to be written into 
! the ocean restart file. This is the only place need to be modified if a new variable is added!
subroutine ini_ocean_io(year, mesh)
  integer, intent(in)       :: year
  integer                   :: j
  character(500)            :: longname
  character(500)            :: trname, units
  character(4)              :: cyear
  type(t_mesh), intent(in) , target :: mesh
  logical, save :: has_been_called = .false.

  write(cyear,'(i4)') year
  oce_path = trim(ResultPath)//trim(runid)//'.'//cyear//'.oce.restart.nc'

  if(has_been_called) return
  has_been_called = .true.

  !===========================================================================
  !===================== Definition part =====================================
  !===========================================================================
  !___Define the netCDF variables for 2D fields_______________________________
  !___SSH_____________________________________________________________________
  call oce_files%def_node_var('ssh', 'sea surface elevation', 'm',   eta_n, mesh)
  !___ALE related fields______________________________________________________
  call oce_files%def_node_var('hbar', 'ALE surface elevation', 'm',   hbar, mesh)
!!PS   call oce_files%def_node_var('ssh_rhs', 'RHS for the elevation', '?',   ssh_rhs, mesh)
  call oce_files%def_node_var('ssh_rhs_old', 'RHS for the elevation', '?',   ssh_rhs_old, mesh)
  call oce_files%def_node_var('hnode', 'nodal layer thickness', 'm',   hnode, mesh)
  
  !___Define the netCDF variables for 3D fields_______________________________
  call oce_files%def_elem_var('u', 'zonal velocity',        'm/s', UV(1,:,:), mesh)
  call oce_files%def_elem_var('v', 'meridional velocity',   'm/s', UV(2,:,:), mesh)
  call oce_files%def_elem_var('urhs_AB', 'Adams–Bashforth for u', 'm/s', UV_rhsAB(1,:,:), mesh)
  call oce_files%def_elem_var('vrhs_AB', 'Adams–Bashforth for v', 'm/s', UV_rhsAB(2,:,:), mesh)
  
  !___Save restart variables for TKE and IDEMIX_________________________________
  if (trim(mix_scheme)=='cvmix_TKE' .or. trim(mix_scheme)=='cvmix_TKE+IDEMIX') then
        call oce_files%def_node_var('tke', 'Turbulent Kinetic Energy', 'm2/s2', tke(:,:), mesh)
  endif
  if (trim(mix_scheme)=='cvmix_IDEMIX' .or. trim(mix_scheme)=='cvmix_TKE+IDEMIX') then
        call oce_files%def_node_var('iwe', 'Internal Wave eneryy', 'm2/s2', tke(:,:), mesh)
  endif 
  if (visc_option==8) then
        call oce_files%def_elem_var('uke', 'unresolved kinetic energy', 'm2/s2', uke(:,:), mesh)
        call oce_files%def_elem_var('uke_rhs', 'unresolved kinetic energy rhs', 'm2/s2', uke_rhs(:,:), mesh)
  endif
  
  do j=1,num_tracers
     SELECT CASE (j) 
       CASE(1)
         trname='temp'
         longname='potential temperature'
         units='degC'
       CASE(2)
         trname='salt'
         longname='salinity'
         units='psu'
       CASE DEFAULT
         write(trname,'(A3,i1)') 'tra_', j
         write(longname,'(A15,i1)') 'passive tracer ', j
         units='none'
     END SELECT
     call oce_files%def_node_var(trim(trname), trim(longname), trim(units), tr_arr(:,:,j), mesh)
     longname=trim(longname)//', Adams–Bashforth'
     call oce_files%def_node_var(trim(trname)//'_AB', trim(longname), trim(units), tr_arr_old(:,:,j), mesh)
  end do
  call oce_files%def_node_var('w', 'vertical velocity', 'm/s', Wvel, mesh)
  call oce_files%def_node_var('w_expl', 'vertical velocity', 'm/s', Wvel_e, mesh)
  call oce_files%def_node_var('w_impl', 'vertical velocity', 'm/s', Wvel_i, mesh)
end subroutine ini_ocean_io
!
!--------------------------------------------------------------------------------------------
! ini_ice_io initializes ice_file datatype which contains information of all variables need to be written into 
! the ice restart file. This is the only place need to be modified if a new variable is added!
subroutine ini_ice_io(year, mesh)
  integer,      intent(in)  :: year
  character(4)              :: cyear
  type(t_mesh), intent(in) , target :: mesh
  logical, save :: has_been_called = .false.

  write(cyear,'(i4)') year
  ice_path = trim(ResultPath)//trim(runid)//'.'//cyear//'.ice.restart.nc'

  if(has_been_called) return
  has_been_called = .true.

  !===========================================================================
  !===================== Definition part =====================================
  !===========================================================================
  !___Define the netCDF variables for 2D fields_______________________________
  call ice_files%def_node_var('area', 'ice concentration [0 to 1]', '%',   a_ice, mesh)
  call ice_files%def_node_var('hice', 'effective ice thickness',    'm',   m_ice, mesh)
  call ice_files%def_node_var('hsnow', 'effective snow thickness',   'm',   m_snow, mesh)
  call ice_files%def_node_var('uice', 'zonal velocity',             'm/s', u_ice, mesh)
  call ice_files%def_node_var('vice', 'meridional velocity',        'm',   v_ice, mesh)
#if defined (__oifs)
  call ice_files%def_node_var('ice_albedo', 'ice albedo',                 '-',   ice_alb, mesh)
  call ice_files%def_node_var('ice_temp', 'ice surface temperature',  'K',   ice_temp, mesh)
#endif /* (__oifs) */

end subroutine ini_ice_io
!
!--------------------------------------------------------------------------------------------
!
subroutine restart(istep, l_read, mesh)

#if defined(__icepack)
  icepack restart not merged here ! produce a compiler error if USE_ICEPACK=ON; todo: merge icepack restart from 68d8b8b
#endif

  use fortran_utils
  implicit none
  ! this is the main restart subroutine
  ! if l_read   is TRUE the restart file will be read

  integer :: istep
  logical :: l_read
  logical :: is_restart
  type(t_mesh), intent(in) , target :: mesh
  logical dumpfiles_exist
  logical, save :: initialized = .false.
  integer cstat, estat
  character(500) cmsg ! there seems to be no documentation about the max size of this text
  
  if(.not. initialized) then
    initialized = .true.
    raw_restart_dirpath = trim(ResultPath)//"/raw_restart/np"//int_to_txt(npes)
    if(raw_restart_length_unit /= "off") then
      if(mype == 0) then
        print *,"creating raw restart directory: "//raw_restart_dirpath
        call execute_command_line("mkdir -p "//raw_restart_dirpath, exitstat=estat, cmdstat=cstat, cmdmsg=cmsg) ! sometimes does not work on aleph
        if(cstat /= 0) print *,"creating raw restart directory ERROR ", trim(cmsg)
      end if
      call MPI_Barrier(MPI_COMM_FESOM, mpierr) ! make sure the dir has been created before we continue...
    end if
  end if

  ctime=timeold+(dayold-1.)*86400
  if (.not. l_read) then
               call ini_ocean_io(yearnew, mesh)
  if (use_ice) call ini_ice_io  (yearnew, mesh)
  else
               call ini_ocean_io(yearold, mesh)
  if (use_ice) call ini_ice_io  (yearold, mesh)
  end if

  if (l_read) then
    ! determine if we can load raw restart dump files
    if(mype == 0) then
      inquire(file=trim(ResultPath)//"/raw_restart/np"//int_to_txt(npes)//".info", exist=dumpfiles_exist)
    end if
    call MPI_Bcast(dumpfiles_exist, 1, MPI_LOGICAL, 0, MPI_COMM_FESOM, MPIerr)
    call read_restart(oce_path, oce_files)
   if (use_ice) then
      call read_restart(ice_path, ice_files)
    end if
  end if

  if (istep==0) return

  !check whether restart will be written
  is_restart = is_due(restart_length_unit, restart_length, istep)

  if (.not. is_restart) return

  ! write restart
  if(mype==0) write(*,*)'Do output (netCDF, restart) ...'
  call write_restart(oce_path, oce_files, istep)
  if (use_ice) then
    call write_restart(ice_path, ice_files, istep)
  end if
  
  ! actualize clock file to latest restart point
  if (mype==0) then
		write(*,*) ' --> actualize clock file to latest restart point'
		call clock_finish  
  end if
  
end subroutine restart


subroutine write_restart(path, filegroup, istep)
  character(len=*), intent(in) :: path
  type(restart_file_group), intent(inout) :: filegroup
  integer,  intent(in)          :: istep
  ! EO parameters
  integer cstep
  integer i
  character(:), allocatable :: dirpath
  
  cstep = globalstep+istep
  
  do i=1, filegroup%nfiles
    call filegroup%files(i)%join() ! join the previous write (if required)

    if(filegroup%files(i)%is_iorank()) then
      if(filegroup%files(i)%is_attached()) call filegroup%files(i)%close_file() ! close the file from previous write
      
      dirpath = path(1:len(path)-3) ! chop of the ".nc" suffix
      if(filegroup%files(i)%path .ne. dirpath//"/"//filegroup%files(i)%varname//".nc") then
        call execute_command_line("mkdir -p "//dirpath)
        filegroup%files(i)%path = dirpath//"/"//filegroup%files(i)%varname//".nc"
        call filegroup%files(i)%open_write_create(filegroup%files(i)%path)
     else
        call filegroup%files(i)%open_write_append(filegroup%files(i)%path) ! todo: keep the file open between writes
      end if

      write(*,*) 'writing restart record ', filegroup%files(i)%rec_count()+1, ' to ', filegroup%files(i)%path
      ! todo: write iter to a separate (non-mesh-variable) file
      call filegroup%files(i)%write_var(filegroup%files(i)%iter_varindex, [filegroup%files(i)%rec_count()+1], [1], [cstep])
      ! todo: write time via the fesom_file_type
      call filegroup%files(i)%write_var(filegroup%files(i)%time_varindex(), [filegroup%files(i)%rec_count()+1], [1], [ctime])
    end if

    call filegroup%files(i)%async_gather_and_write_variables()
  end do
  
end subroutine


! join remaining threads and close all open files
subroutine finalize_restart()
  integer i

  ! join all previous writes
  ! close all restart files

  do i=1, oce_files%nfiles
    call oce_files%files(i)%join()
    if(oce_files%files(i)%is_iorank()) then
      if(oce_files%files(i)%is_attached()) call oce_files%files(i)%close_file()
    end if
  end do

  if(use_ice) then
    do i=1, ice_files%nfiles
      call ice_files%files(i)%join()
      if(ice_files%files(i)%is_iorank()) then
        if(ice_files%files(i)%is_attached()) call ice_files%files(i)%close_file()
      end if
    end do
  end if
end subroutine


subroutine read_restart(path, filegroup)
  character(len=*), intent(in) :: path
  type(restart_file_group), intent(inout) :: filegroup
  ! EO parameters
  real(kind=WP) rtime
  integer i
  character(:), allocatable :: dirpath

  do i=1, filegroup%nfiles
    if( filegroup%files(i)%is_iorank() ) then
      dirpath = path(1:len(path)-3) ! chop of the ".nc" suffix
      if(filegroup%files(i)%path .ne. dirpath//"/"//filegroup%files(i)%varname//".nc") then
        call execute_command_line("mkdir -p "//dirpath)
        filegroup%files(i)%path = dirpath//"/"//filegroup%files(i)%varname//".nc"
#ifndef DISABLE_PARALLEL_RESTART_READ
        write(*,*) 'reading restart PARALLEL for ', filegroup%files(i)%varname, ' at ', filegroup%files(i)%path
#else
        write(*,*) 'reading restart SEQIENTIAL for ', filegroup%files(i)%varname, ' at ', filegroup%files(i)%path
#endif
        call filegroup%files(i)%open_read(filegroup%files(i)%path) ! do we need to bother with read-only access?
        ! todo: print a reasonable error message if the file does not exist
      end if
    end if

    call filegroup%files(i)%async_read_and_scatter_variables()
#ifndef DISABLE_PARALLEL_RESTART_READ
  end do
  
  do i=1, filegroup%nfiles
#endif
    call filegroup%files(i)%join()

    if(filegroup%files(i)%is_iorank()) then
      write(*,*) 'restart from record ', filegroup%files(i)%rec_count(), ' of ', filegroup%files(i)%rec_count(), filegroup%files(i)%path

      ! read the last entry from the iter variable
      call filegroup%files(i)%read_var1(filegroup%files(i)%iter_varindex, [filegroup%files(i)%rec_count()], globalstep)

      ! read the last entry from the time variable
      call filegroup%files(i)%read_var1(filegroup%files(i)%time_varindex(), [filegroup%files(i)%rec_count()], rtime)
      call filegroup%files(i)%close_file()

     if (int(ctime)/=int(rtime)) then
        write(*,*) 'Reading restart: timestamps in restart and in clock files do not match'
        write(*,*) 'restart/ times are:', ctime, rtime
        write(*,*) 'the model will stop!'
        stop 1
      end if
    end if    
  end do
end subroutine


  function is_due(unit, frequency, istep) result(d)
    character(len=*), intent(in) :: unit
    integer, intent(in) :: frequency
    integer, intent(in) :: istep
    logical d
    ! EO parameters
    d = .false.
    
    if(unit.eq.'y') then
      call annual_event(d)
    else if(unit.eq.'m') then 
      call monthly_event(d) 
    else if(unit.eq.'d') then
      call daily_event(d, frequency)
    else if(unit.eq.'h') then
      call hourly_event(d, frequency)
    else if(unit.eq.'s') then
      call step_event(d, istep, frequency)
    else
      write(*,*) 'You did not specify a supported outputflag.'
      write(*,*) 'The program will stop to give you opportunity to do it.'
      stop 1
    stop
    end if
  end function

end module
