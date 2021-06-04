MODULE io_RESTART
  use restart_file_group_module
  use g_config
  use g_clock
  use g_parsup
  use g_comm_auto
  use mod_mesh
  use o_arrays
  use i_arrays
  use g_cvmix_tke
  use g_cvmix_idemix
  implicit none
#include "netcdf.inc"
!
!--------------------------------------------------------------------------------------------
! 
  type nc_dims
    integer        :: size
    character(100) :: name
    integer        :: code
  end type nc_dims
!
!--------------------------------------------------------------------------------------------
!
  type nc_vars
    character(100) :: name
    integer        :: code
    character(500) :: longname
    character(100) :: units
    integer        :: ndim
    integer        :: dims(2) !<=2; assume there are no variables with dimension more than 2xNLxT
    real(kind=WP), pointer :: pt1(:), pt2(:,:)
  end type nc_vars
!
!--------------------------------------------------------------------------------------------
!
  type nc_file
    character(500)                                :: filename
    type(nc_dims), allocatable, dimension(:) :: dim
    type(nc_vars), allocatable, dimension(:) :: var
    integer :: ndim=0, nvar=0
    integer :: rec_dimid, time_varid, iter_varid
    integer :: ncid
    integer :: rec_count=0
    integer :: error_status(250), error_count
    logical :: is_in_use=.false.
  end type nc_file


  type(nc_file), save       :: ocean_file, ice_file
  integer,       save       :: globalstep=0
  real(kind=WP)             :: ctime !current time in seconds from the beginning of the year

  PRIVATE
  PUBLIC :: restart 
  
  type(restart_file_group), save :: oce_files
  type(restart_file_group), save :: ice_files
  character(:), allocatable, save :: oce_path
  character(:), allocatable, save :: ice_path


  contains
!
!--------------------------------------------------------------------------------------------
! ini_ocean_io initializes ocean_file datatype which contains information of all variables need to be written into 
! the ocean restart file. This is the only place need to be modified if a new variable is added!
subroutine ini_ocean_io(year, mesh)
  implicit none

  integer, intent(in)       :: year
  integer                   :: ncid, j
  integer                   :: varid
  character(500)            :: longname
  character(500)            :: filename
  character(500)            :: trname, units
  character(4)              :: cyear
  type(t_mesh), intent(in) , target :: mesh

#include  "associate_mesh.h"

  write(cyear,'(i4)') year
  oce_path = trim(ResultPath)//trim(runid)//'.'//cyear//'.oce.restart.nc'
  if (ocean_file%is_in_use) return
  ocean_file%is_in_use=.true.
  call def_dim(ocean_file, 'node', nod2d)
  call def_dim(ocean_file, 'elem', elem2d)
  call def_dim(ocean_file, 'nz_1', nl-1)
  call def_dim(ocean_file, 'nz',   nl)

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
  implicit none

  integer,      intent(in)  :: year
  integer                   :: ncid, j
  integer                   :: varid
  character(500)            :: longname
  character(500)            :: filename
  character(500)            :: trname, units
  character(4)              :: cyear
  type(t_mesh), intent(in) , target :: mesh

#include  "associate_mesh.h"

  write(cyear,'(i4)') year
  ice_path = trim(ResultPath)//trim(runid)//'.'//cyear//'.ice.restart.nc'
  if (ice_file%is_in_use) return
  ice_file%is_in_use=.true.
  call def_dim(ice_file, 'node', nod2d)

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
subroutine restart(istep, l_write, l_read, mesh)
  implicit none
  ! this is the main restart subroutine
  ! if l_write  is TRUE writing restart file will be forced
  ! if l_read   is TRUE the restart file will be read

  integer :: istep
  logical :: l_write, l_read
  logical :: is_restart
  integer :: mpierr
  type(t_mesh), intent(in) , target :: mesh

  ctime=timeold+(dayold-1.)*86400
  if (.not. l_read) then
    call ini_ocean_io(yearnew, mesh)
    if (use_ice) call ini_ice_io  (yearnew, mesh)
  else
    call ini_ocean_io(yearold, mesh)
    if (use_ice) call ini_ice_io  (yearold, mesh)
  end if

  if (l_read) then
   call assoc_ids(ocean_file);          call was_error(ocean_file)
   call read_restart(ocean_file, mesh); call was_error(ocean_file)
   if (use_ice) then
      call assoc_ids(ice_file);          call was_error(ice_file)
      call read_restart(ice_file, mesh); call was_error(ice_file)
   end if
  end if

  if (istep==0) return

  !check whether restart will be written
  is_restart=.false.

  if (restart_length_unit.eq.'y') then
     call annual_event(is_restart)
  else if (restart_length_unit.eq.'m') then 
     call monthly_event(is_restart) 
  else if (restart_length_unit.eq.'d') then
     call daily_event(is_restart, restart_length)
  else if (restart_length_unit.eq.'h') then
     call hourly_event(is_restart, restart_length)
  else if (restart_length_unit.eq.'s') then
     call step_event(is_restart, istep, restart_length)
  else
     write(*,*) 'You did not specify a supported outputflag.'
     write(*,*) 'The program will stop to give you opportunity to do it.'
     call par_ex(1)
     stop
  endif

  if (l_write) is_restart=.true.  

  if (.not. is_restart) return

  ! write restart
  if(mype==0) write(*,*)'Do output (netCDF, restart) ...'
  call assoc_ids(ocean_file);                  call was_error(ocean_file)  
  call write_restart(ocean_file, istep, mesh); call was_error(ocean_file)
  if (use_ice) then
     call assoc_ids(ice_file);                  call was_error(ice_file)  
     call write_restart(ice_file, istep, mesh); call was_error(ice_file)
  end if
  
  ! actualize clock file to latest restart point
  if (mype==0) then
		write(*,*) ' --> actualize clock file to latest restart point'
		call clock_finish  
  end if
  
end subroutine restart


subroutine create_new_file(file)
  implicit none

  type(nc_file),  intent(inout) :: file
  integer                       :: c, j
  integer                       :: n, k, l, kdim, dimid(4)
  character(2000)               :: att_text
  ! Serial output implemented so far
  if (mype/=0) return
  c=1
  file%error_status=0
  ! create an ocean output file
  write(*,*) 'initializing restart file ', trim(file%filename)
  file%error_status(c) = nf_create(file%filename, IOR(NF_NOCLOBBER,IOR(NF_NETCDF4,NF_CLASSIC_MODEL)), file%ncid); c=c+1

  do j=1, file%ndim
!___Create mesh related dimentions__________________________________________
     file%error_status(c) = nf_def_dim(file%ncid, file%dim(j)%name, file%dim(j)%size, file%dim(j)%code ); c=c+1
  end do

!___Create time related dimentions__________________________________________
  file%error_status(c) = nf_def_dim(file%ncid, 'time', NF_UNLIMITED, file%rec_dimid);         c=c+1
!___Define the time and iteration variables_________________________________
  file%error_status(c) = nf_def_var(file%ncid, 'time', NF_DOUBLE, 1, file%rec_dimid, file%time_varid); c=c+1
  file%error_status(c) = nf_def_var(file%ncid, 'iter', NF_INT,    1, file%rec_dimid, file%iter_varid); c=c+1


  att_text='time'
  file%error_status(c) = nf_put_att_text(file%ncid, file%time_varid, 'long_name', len_trim(att_text), trim(att_text)); c=c+1
  write(att_text, '(a14,I4.4,a1,I2.2,a1,I2.2,a6)') 'seconds since ', yearold, '-', 1, '-', 1, ' 0:0:0'
  file%error_status(c) = nf_put_att_text(file%ncid, file%time_varid, 'units', len_trim(att_text), trim(att_text)); c=c+1

  att_text='iteration_count'
  file%error_status(c) = nf_put_att_text(file%ncid, file%iter_varid, 'long_name', len_trim(att_text), trim(att_text)); c=c+1

  do j=1, file%nvar
!___associate physical dimension with the netcdf IDs________________________
     n=file%var(j)%ndim ! shape size of the variable (exluding time)
     do k=1, n
        !k_th dimension of the variable
        kdim=file%var(j)%dims(k)
        do l=1, file%ndim ! list all defined dimensions 
           if (kdim==file%dim(l)%size) dimid(k)=file%dim(l)%code
        end do
        !write(*,*) "j",j,kdim, ' -> ', dimid(k)
     end do
     file%error_status(c) = nf_def_var(file%ncid, trim(file%var(j)%name), NF_DOUBLE, file%var(j)%ndim+1, (/dimid(1:n), file%rec_dimid/), file%var(j)%code); c=c+1
     !if (n==1) then
     !   file%error_status(c)=nf_def_var_chunking(file%ncid, file%var(j)%code, NF_CHUNKED, (/1/)); c=c+1 
     if (n==2) then
        file%error_status(c)=nf_def_var_chunking(file%ncid, file%var(j)%code, NF_CHUNKED, (/1, file%dim(1)%size/)); ! c=c+1 
     end if
     file%error_status(c)=nf_put_att_text(file%ncid, file%var(j)%code, 'description', len_trim(file%var(j)%longname), file%var(j)%longname); c=c+1
     file%error_status(c)=nf_put_att_text(file%ncid, file%var(j)%code, 'units',       len_trim(file%var(j)%units),    file%var(j)%units);    c=c+1
  end do

  file%error_status(c)=nf_close(file%ncid); c=c+1
  file%error_count=c-1
end subroutine create_new_file


subroutine def_dim(file, name, ndim)
  implicit none
  type(nc_file),    intent(inout) :: file
  character(len=*), intent(in)    :: name
  integer,          intent(in)    :: ndim
  type(nc_dims), allocatable, dimension(:) :: temp

  if (file%ndim > 0) then
     ! create temporal dimension
     allocate(temp(file%ndim)); temp=file%dim
     ! deallocate the input data array
     deallocate(file%dim)
     ! then reallocate
     file%ndim=file%ndim+1
     allocate(file%dim(file%ndim))
     ! restore the original data
     file%dim(1:file%ndim-1)=temp  
     deallocate(temp)
   else
     ! first dimension in a file
     file%ndim=1
     allocate(file%dim(file%ndim))
   end if
   file%dim(file%ndim)%name=trim(name)
   file%dim(file%ndim)%size=ndim
end subroutine def_dim


subroutine write_restart(file, istep, mesh)
  implicit none
  type(nc_file),  intent(inout) :: file
  integer,  intent(in)          :: istep
  type(t_mesh), intent(in)     , target :: mesh
  real(kind=WP), allocatable    :: aux(:), laux(:)
  integer                       :: i, lev, size1, size2, shape
  integer                       :: c
  real(kind=WP)                 :: t0, t1, t2, t3

#include  "associate_mesh.h"

  ! Serial output implemented so far
  if (mype==0) then
     c=1
     !file%rec_count=file%rec_count+1
     write(*,*) 'writing restart record ', file%rec_count
     file%error_status(c)=nf_open(file%filename, nf_write, file%ncid); c=c+1
     file%error_status(c)=nf_put_vara_double(file%ncid, file%time_varid, file%rec_count, 1, ctime, 1); c=c+1
     file%error_status(c)=nf_put_vara_int(file%ncid,    file%iter_varid, file%rec_count, 1, globalstep+istep, 1);   c=c+1
  end if

  call was_error(file); c=1

  do i=1, file%nvar
     shape=file%var(i)%ndim
!_______writing 2D fields________________________________________________
     if (shape==1) then
        size1=file%var(i)%dims(1)
        if (mype==0) allocate(aux(size1))
        t0=MPI_Wtime()
        if (size1==nod2D)  call gather_nod (file%var(i)%pt1, aux)
        if (size1==elem2D) call gather_elem(file%var(i)%pt1, aux)
        t1=MPI_Wtime()
        if (mype==0) then
           file%error_status(c)=nf_put_vara_double(file%ncid, file%var(i)%code, (/1, file%rec_count/), (/size1, 1/), aux, 1); c=c+1
        end if
        t2=MPI_Wtime()
#ifdef DEBUG
        ! Timeing information for collecting and writing restart file
        if (mype==0) write(*,*) 'nvar: ', i, 'size: ', size1, 'gather_nod: ', t1-t0
        if (mype==0) write(*,*) 'nvar: ', i, 'size: ', size1, 'nf_put_var: ', t2-t1
#endif
        if (mype==0) deallocate(aux)
!_______writing 3D fields________________________________________________
     elseif (shape==2) then
        size1=file%var(i)%dims(1)
        size2=file%var(i)%dims(2)
        if (mype==0)       allocate(aux (size2))
        if (size2==nod2D)  allocate(laux(myDim_nod2D +eDim_nod2D ))
        if (size2==elem2D) allocate(laux(myDim_elem2D+eDim_elem2D))
        do lev=1, size1
           laux=file%var(i)%pt2(lev,:)
!          if (size1==nod2D  .or. size2==nod2D)  call gather_nod (file%var(i)%pt2(lev,:), aux)
!          if (size1==elem2D .or. size2==elem2D) call gather_elem(file%var(i)%pt2(lev,:), aux)
           t0=MPI_Wtime()
           if (size1==nod2D  .or. size2==nod2D)  call gather_nod (laux, aux)
           if (size1==elem2D .or. size2==elem2D) call gather_elem(laux, aux)
           t1=MPI_Wtime()
           if (mype==0) then
              file%error_status(c)=nf_put_vara_double(file%ncid, file%var(i)%code, (/lev, 1, file%rec_count/), (/1, size2, 1/), aux, 1); c=c+1
           end if
           t2=MPI_Wtime()
#ifdef DEBUG
           ! Timeing information for collecting and writing output file
           if (mype==0) write(*,*) 'nvar: ', i, 'size: ', size2, 'lev: ', lev, 'gather_nod: ', t1-t0
           if (mype==0) write(*,*) 'nvar: ', i, 'size: ', size2, 'lev: ', lev, 'nf_put_var: ', t2-t1
#endif
        end do
        deallocate(laux)
        if (mype==0) deallocate(aux)
     else
        if (mype==0) write(*,*) 'not supported shape of array in restart file'
           call par_ex
           stop
     end if
     call was_error(file); c=1
  end do

  if (mype==0) file%error_count=c-1
  call was_error(file)
  if (mype==0) file%error_status(1)=nf_close(file%ncid);
  file%error_count=1
  call was_error(file)
end subroutine write_restart


subroutine read_restart(file, mesh, arg)
  implicit none
  type(nc_file),     intent(inout) :: file
  integer, optional, intent(in)    :: arg
  real(kind=WP), allocatable       :: aux(:), laux(:)
  integer                          :: i, lev, size1, size2, shape
  integer                          :: rec2read, c
  real(kind=WP)                    :: rtime !timestamp of the record
  logical                          :: file_exist=.False.
  type(t_mesh), intent(in)        , target :: mesh

#include  "associate_mesh.h"

  ! laux=0.
  ! Serial output implemented so far
  c=1
  if (mype==0) then
     file_exist=.False.
     inquire(file=file%filename,exist=file_exist) 
     if (file_exist) then
        write(*,*) '     reading restart file:  ', trim(file%filename)
        file%error_status(c)=nf_open(file%filename, nf_nowrite, file%ncid);                           c=c+1
        file%error_status(c)=nf_get_vara_int(file%ncid,    file%iter_varid, file%rec_count, 1, globalstep, 1); c=c+1
        file%error_status(c)=nf_get_vara_double(file%ncid, file%time_varid, file%rec_count, 1, rtime, 1);      c=c+1
     else
        write(*,*)
        print *, achar(27)//'[33m'
        write(*,*) '____________________________________________________________________'
        write(*,*) ' ERROR: could not find restart_file:',trim(file%filename),'!'    
        write(*,*) '____________________________________________________________________'
        print *, achar(27)//'[0m'
        write(*,*)
        call par_ex
     end if 
     
     if (.not. present(arg)) then
        rec2read=file%rec_count
     else
        rec2read=arg
     end if
     write(*,*) 'restart from record ', rec2read, ' of ', file%rec_count

     if (int(ctime)/=int(rtime)) then
        write(*,*) 'Reading restart: timestamps in restart and in clock files do not match'
        write(*,*) 'restart/ times are:', ctime, rtime
        write(*,*) 'the model will stop!'
        file%error_status(c)=-310; c=c+1
     end if
  end if

  call was_error(file); c=1
 
  do i=1, file%nvar
     shape=file%var(i)%ndim
     if (mype==0) write(*,*) 'reading restart for ', trim(file%var(i)%name)
!_______writing 2D fields________________________________________________
     if (shape==1) then
        size1=file%var(i)%dims(1)
        if (mype==0) then
           allocate(aux(size1))
           file%error_status(c)=nf_get_vara_double(file%ncid, file%var(i)%code, (/1, file%rec_count/), (/size1, 1/), aux, 1); c=c+1
!          write(*,*) 'min/max 2D =', minval(aux), maxval(aux)
        end if
        if (size1==nod2D)  call broadcast_nod (file%var(i)%pt1, aux)
        if (size1==elem2D) call broadcast_elem(file%var(i)%pt1, aux)
        if (mype==0) deallocate(aux)
!_______writing 3D fields________________________________________________
     elseif (shape==2) then
        size1=file%var(i)%dims(1)
        size2=file%var(i)%dims(2)
        if (mype==0)       allocate(aux (size2))
        if (size2==nod2D)  allocate(laux(myDim_nod2D +eDim_nod2D ))
        if (size2==elem2D) allocate(laux(myDim_elem2D+eDim_elem2D))        
        do lev=1, size1
           if (mype==0) then
              file%error_status(c)=nf_get_vara_double(file%ncid, file%var(i)%code, (/lev, 1, file%rec_count/), (/1, size2, 1/), aux, 1); c=c+1
!             write(*,*) 'min/max 3D ', lev,'=', minval(aux), maxval(aux)
           end if
           file%var(i)%pt2(lev,:)=0.
!          if (size1==nod2D  .or. size2==nod2D)  call broadcast_nod (file%var(i)%pt2(lev,:), aux)
!          if (size1==elem2D .or. size2==elem2D) call broadcast_elem(file%var(i)%pt2(lev,:), aux)
           if (size2==nod2D)  then
              call broadcast_nod (laux, aux)
              file%var(i)%pt2(lev,:)=laux(1:myDim_nod2D+eDim_nod2D)
           end if
           if (size2==elem2D) then
              call broadcast_elem(laux, aux)
              file%var(i)%pt2(lev,:)=laux(1:myDim_elem2D+eDim_elem2D)
           end if
        end do
        deallocate(laux)
        if (mype==0) deallocate(aux)
     else
        if (mype==0) write(*,*) 'not supported shape of array in restart file when reading restart'
           call par_ex
           stop
     end if
     call was_error(file); c=1
  end do

  if (mype==0) file%error_status(1)=nf_close(file%ncid);
  file%error_count=1
  call was_error(file)
end subroutine read_restart


subroutine assoc_ids(file)
  implicit none

  type(nc_file),  intent(inout) :: file
  character(500)                :: longname
  integer                       :: c, j, k
  real(kind=WP)                 :: rtime !timestamp of the record
  ! Serial output implemented so far
  if (mype/=0) return
  c=1
  file%error_status=0
  ! open existing netcdf file
  write(*,*) 'associating restart file ', trim(file%filename)

  file%error_status(c) = nf_open(file%filename, nf_nowrite, file%ncid)
  !if the file does not exist it will be created!
  if (file%error_status(c) .ne. nf_noerr) then
     call create_new_file(file) ! error status counter will be reset
     c=file%error_count+1
     file%error_status(c) = nf_open(file%filename, nf_nowrite, file%ncid); c=c+1
  end if

  do j=1, file%ndim
!___Associate mesh related dimentions_______________________________________
    file%error_status(c) = nf_inq_dimid(file%ncid, file%dim(j)%name, file%dim(j)%code); c=c+1
  end do
!___Associate time related dimentions_______________________________________
  file%error_status(c) = nf_inq_dimid (file%ncid, 'time', file%rec_dimid);       c=c+1
  file%error_status(c) = nf_inq_dimlen(file%ncid, file%rec_dimid, file%rec_count); c=c+1
!___Associate the time and iteration variables______________________________
  file%error_status(c) = nf_inq_varid(file%ncid, 'time', file%time_varid); c=c+1
  file%error_status(c) = nf_inq_varid(file%ncid, 'iter', file%iter_varid); c=c+1
!___if the time rtime at the rec_count does not equal ctime we look for the closest record with the 
! timestamp less than ctime
  do k=file%rec_count, 1, -1
     file%error_status(c)=nf_get_vara_double(file%ncid, file%time_varid, k, 1, rtime, 1);
     if (ctime > rtime) then
        file%rec_count=k+1
        exit ! a proper rec_count detected, ready for writing restart, exit the loop
     elseif (ctime == rtime) then
        file%rec_count=k
        exit ! a proper rec_count detected, ready for reading restart, exit the loop
     end if
     if (k==1) then
        if (mype==0) write(*,*) 'WARNING: all dates in restart file are after the current date'
        if (mype==0) write(*,*) 'reading restart will not be possible !'
        if (mype==0) write(*,*) 'the model attempted to start with the time stamp = ', int(ctime)
        file%error_status(c)=-310;
     end if
  end do
  c=c+1 ! check will be made only for the last nf_get_vara_double
  file%rec_count=max(file%rec_count, 1)
!___Associate physical variables____________________________________________
  do j=1, file%nvar
     file%error_status(c) = nf_inq_varid(file%ncid, file%var(j)%name, file%var(j)%code); c=c+1
  end do
  file%error_status(c)=nf_close(file%ncid); c=c+1
  file%error_count=c-1
  write(*,*) 'current restart counter = ',       file%rec_count
end subroutine assoc_ids
!
!--------------------------------------------------------------------------------------------
!
subroutine was_error(id)
  implicit none
  type(nc_file),  intent(inout) :: id
  integer                       :: k, status, ierror

  call MPI_BCast(id%error_count, 1,  MPI_INTEGER, 0, MPI_COMM_FESOM, ierror)
  call MPI_BCast(id%error_status(1), id%error_count, MPI_INTEGER, 0, MPI_COMM_FESOM, ierror)

  do k=1, id%error_count
     status=id%error_status(k)
     if (status .ne. nf_noerr) then
        if (mype==0) write(*,*) 'error counter=', k
        if (mype==0) call handle_err(status)
        call par_ex
        stop
     end if
  end do
end subroutine was_error
END MODULE io_RESTART
