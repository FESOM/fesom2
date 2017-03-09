MODULE io_RESTART
  use g_config
  use g_clock
  use g_parsup
  use g_comm_auto
  use o_mesh
  use o_arrays
  use i_arrays
  implicit none
#include "netcdf.inc"
 
  type nc_file_dims
    integer        :: size
    character(100) :: name
    integer        :: code
  end type nc_file_dims

  type nc_file_vars
    character(100) :: name
    integer        :: code
    character(500) :: longname
    character(100) :: units
    integer        :: ndim
    integer        :: dims(3) !<=3; assume there are no variables with dimension more than 3xT
    real(kind=WP), pointer :: pt1(:), pt2(:,:)
  end type nc_file_vars


  type nc_file
    character(500)                                :: filename
    type(nc_file_dims), allocatable, dimension(:) :: dim
    type(nc_file_vars), allocatable, dimension(:) :: var
    integer :: ndim=0, nvar=0
    integer :: rec, Tid, Iid
    integer :: ncid
    integer :: error_status(100), error_count
    integer :: rec_count=0
  end type nc_file

  type type_id
    integer :: nd, el, nz, nz1, T, rec, iter
  end type type_id


  ! id will keep the IDs of all required dimentions and variables
  type(nc_file), save       :: oid
  integer                   :: error_status(100), error_count

  PRIVATE
  PUBLIC :: check_restart
  !generic interface was required to associate variables of unknown rank with the pointers of the same rank
  !this allows for automatic streaming of associated variables into the netcdf file
  INTERFACE def_variable
            MODULE PROCEDURE def_variable_1d, def_variable_2d
  END INTERFACE

  contains

subroutine create_restart_file(l_create)
  implicit none

  logical, intent(in)       :: l_create  
  integer                   :: ncid, j
  integer                   :: varid
  character(500)            :: longname
  character(500)            :: filename
  character(500)            :: trname, units

  if (.not. l_create) return
  ! create an ocean restart file; serial output implemented so far
  oid.filename=trim(ResultPath)//trim(runid)//'.'//cyearnew//'.oce.restart.nc'
  call def_dim(oid, 'node', nod2d)
  call def_dim(oid, 'elem', elem2d)
  call def_dim(oid, 'nz_1', nl-1)
  call def_dim(oid, 'nz',   nl)

  !===========================================================================
  !===================== Definition part =====================================
  !===========================================================================
  !___Define the netCDF variables for 2D fields_______________________________
  !___SSH_____________________________________________________________________
  call def_variable(oid, 'ssh',         (/nod2D/),  'sea surface elevation', 'm', eta_n);
  !___ALE related fields______________________________________________________
  call def_variable(oid, 'hbar',        (/nod2D/), 'ALE surface elevation hbar_n+0.5', 'm', hbar);
  call def_variable(oid, 'hbar_old',    (/nod2D/), 'ALE surface elevation hbar_n-0.5', 'm', hbar_old);
  call def_variable(oid, 'ssh_rhs',     (/nod2D/), 'RHS for the elevation', '?', ssh_rhs);
  call def_variable(oid, 'ssh_rhs_old', (/nod2D/), 'RHS for the elevation', '?', ssh_rhs_old);
  !___Define the netCDF variables for 3D fields_______________________________
  call def_variable(oid, 'hnode',    (/nl-1,  nod2D/), 'ALE stuff', '?', hnode);
  call def_variable(oid, 'u',        (/nl-1, elem2D/), 'zonal velocity', 'm/s', UV(1,:,:));
  call def_variable(oid, 'v',        (/nl-1, elem2D/), 'meridional velocity', 'm/s', UV(2,:,:));
  call def_variable(oid, 'urhs_AB',  (/nl-1, elem2D/), 'Adams–Bashforth for u', 'm/s', UV_rhsAB(1,:,:));
  call def_variable(oid, 'vrhs_AB',  (/nl-1, elem2D/), 'Adams–Bashforth for v', 'm/s', UV_rhsAB(2,:,:));

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
         write(trname,'(A3,i1)') 'ptr', j
         write(longname,'(A15,i1)') 'passive tracer ', j
         units='none'
     END SELECT
     call def_variable(oid, trim(trname),       (/nl-1, nod2D/), trim(longname), trim(units), tr_arr(:,:,j));
     longname=trim(longname)//', Adams–Bashforth'
     call def_variable(oid, trim(trname)//'_AB',(/nl-1, nod2D/), trim(longname), trim(units), tr_arr_old(:,:,j));
  end do
  call def_variable(oid, 'w',      (/nl, nod2D/), 'vertical velocity', 'm/s', Wvel);
  call def_variable(oid, 'w_expl', (/nl, nod2D/), 'vertical velocity', 'm/s', Wvel_e);
  call def_variable(oid, 'w_impl', (/nl, nod2D/), 'vertical velocity', 'm/s', Wvel_i);


  call open_new_file(oid); call was_error

end subroutine create_restart_file
!
!--------------------------------------------------------------------------------------------
!
subroutine open_new_file(id)
  implicit none

  type(nc_file),  intent(inout) :: id
  character(500)                :: longname
  integer                       :: c, j
  integer                       :: n, k, l, kdim, dimid(4)
  ! Serial output implemented so far
  if (mype/=0) return
  c=1
  error_status=0
  ! create an ocean output file
  write(*,*) 'initializing restart file ', trim(id.filename)

  if (restart_offset==32) then
     error_status(c) = nf_create(id.filename, nf_clobber, id.ncid);                      c=c+1
  else
     error_status(c) = nf_create(id.filename, IOR(NF_CLOBBER,NF_64BIT_OFFSET), id.ncid); c=c+1
  end if

do j=1, id.ndim
!___Create mesh related dimentions__________________________________________
  error_status(c) = nf_def_dim(id.ncid, id.dim(j).name, id.dim(j).size, id.dim(j).code ); c=c+1
end do
!___Create time related dimentions__________________________________________
  error_status(c) = nf_def_dim(id.ncid, 'time', NF_UNLIMITED, id.rec);         c=c+1
!___Define the time and iteration variables_________________________________
  error_status(c) = nf_def_var(id.ncid, 'time', NF_DOUBLE, 1, id.rec, id.Tid); c=c+1
  error_status(c) = nf_def_var(id.ncid, 'iter', NF_INT,    1, id.rec, id.Iid); c=c+1

  longname='model time'
  error_status(c) = nf_put_att_text(id.ncid, id%Tid, 'long_name', len_trim(longname), trim(longname)); c=c+1
  error_status(c) = nf_put_att_text(id.ncid, id%Tid, 'units', 1, 's');                                 c=c+1
  longname='iteration_count'
  error_status(c) = nf_put_att_text(id.ncid, id%Iid, 'long_name', len_trim(longname), trim(longname)); c=c+1

  do j=1, id.nvar
!___associate physical dimension with the netcdf IDs________________________
     n=id.var(j).ndim ! shape size of the variable (exluding time)
     do k=1, n
        !k_th dimension of the variable
        kdim=id.var(j).dims(k)
        do l=1, id.ndim ! list all defined dimensions 
           if (kdim==id.dim(l).size) dimid(k)=id.dim(l).code
        end do
!________write(*,*) kdim, ' -> ', dimid(k)__________________________________
     end do
     error_status(c) = nf_def_var(id.ncid, trim(id.var(j).name), NF_DOUBLE, id.var(j).ndim+1, &
                       (/dimid(1:n), id.rec/), id.var(j).code); c=c+1
  end do
  error_status(c)=nf_close(id.ncid); c=c+1
  error_count=c-1
end subroutine open_new_file
!
!--------------------------------------------------------------------------------------------
!
subroutine was_error
  implicit none

  integer                   :: k, status, ierror

call MPI_BCast(error_count, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)
call MPI_BCast(error_status(1), error_count, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)

  do k=1, error_count
     status=error_status(k)
     if (status .ne. nf_noerr) then
        if (mype==0) call handle_err(status)
        call par_ex
        stop
     end if
  end do
end subroutine was_error
!
!--------------------------------------------------------------------------------------------
!
subroutine check_restart(istep, l_write, l_create)
  implicit none
  ! this is the main restart subroutine
  ! if l_write  is TRUE the restart will be forced
  ! if l_create is TRUE the new restart file will be created

  integer :: istep
  logical :: l_create, l_write
  logical :: is_restart
  logical, save :: lfirst=.true.
  integer :: mpierr

  lfirst=.false.
  call create_restart_file(l_create)
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
  call write_restart(oid, istep); call was_error
end subroutine check_restart


subroutine def_dim(id, name, ndim)
  implicit none
  type(nc_file),    intent(inout) :: id
  character(len=*), intent(in)    :: name
  integer,          intent(in)    :: ndim
  type(nc_file_dims), allocatable, dimension(:) :: temp

  if (id.ndim > 0) then
     ! create temporal dimension
     allocate(temp(id.ndim)); temp=id.dim
     ! deallocate the input data array
     deallocate(id.dim)
     ! then reallocate
     id.ndim=id.ndim+1
     allocate(id.dim(id.ndim))
     ! restore the original data
     id.dim(1:id.ndim-1)=temp  
     deallocate(temp)
   else
     ! first dimension in a file
     id.ndim=1
     allocate(id.dim(id.ndim))
   end if
   id.dim(id.ndim).name=trim(name)
   id.dim(id.ndim).size=ndim
end subroutine def_dim



subroutine def_variable_1d(id, name, dims, longname, units, data)
  implicit none
  type(nc_file),    intent(inout)        :: id
  character(len=*), intent(in)           :: name
  integer, intent(in)                    :: dims(1)
  character(len=*), intent(in), optional :: units, longname
  real(kind=8),target,     intent(inout)        :: data(:)
  integer                                :: c
  type(nc_file_vars), allocatable, dimension(:) :: temp

  if (id.nvar > 0) then
     ! create temporal dimension
     allocate(temp(id.nvar)); temp=id.var
     ! deallocate the input data array
     deallocate(id.var)
     ! then reallocate
     id.nvar=id.nvar+1
     allocate(id.var(id.nvar))
     ! restore the original data
     id.var(1:id.nvar-1)=temp  
     deallocate(temp)
   else
     ! first dimension in a file
     id.nvar=1
     allocate(id.var(id.nvar))
   end if
   id.var(id.nvar).name=trim(name)
   id.var(id.nvar).longname=trim(longname)
   id.var(id.nvar).units=trim(units)
   id.var(id.nvar).ndim=1
   id.var(id.nvar).dims(1)=dims(1)
   id.var(id.nvar).pt1=>data
end subroutine def_variable_1d

subroutine def_variable_2d(id, name, dims, longname, units, data)
  implicit none
  type(nc_file),    intent(inout)        :: id
  character(len=*), intent(in)           :: name
  integer, intent(in)                    :: dims(2)
  character(len=*), intent(in), optional :: units, longname
  real(kind=8),target,     intent(inout)        :: data(:,:)
  integer                                :: c
  type(nc_file_vars), allocatable, dimension(:) :: temp

  if (id.nvar > 0) then
     ! create temporal dimension
     allocate(temp(id.nvar)); temp=id.var
     ! deallocate the input data array
     deallocate(id.var)
     ! then reallocate
     id.nvar=id.nvar+1
     allocate(id.var(id.nvar))
     ! restore the original data
     id.var(1:id.nvar-1)=temp  
     deallocate(temp)
   else
     ! first dimension in a file
     id.nvar=1
     allocate(id.var(id.nvar))
   end if
   id.var(id.nvar).name=trim(name)
   id.var(id.nvar).longname=trim(longname)
   id.var(id.nvar).units=trim(units)
   id.var(id.nvar).ndim=2
   id.var(id.nvar).dims(1:2)=dims
   id.var(id.nvar).pt2=>data
end subroutine def_variable_2d

subroutine write_restart(id, istep)
  implicit none
  type(nc_file),  intent(inout) :: id
  integer,  intent(in)       :: istep
  real(kind=8), allocatable  :: aux1(:), aux2(:,:) 
  integer                    :: i, size1, size2, shape
  integer                    :: c
  c=1
  error_status(c)=nf_open(id.filename, nf_write, id.ncid); c=c+1
  id.rec_count=id.rec_count+1
  do i=1, id.nvar
     shape=id.var(i).ndim
!_______writing 2D fields________________________________________________
     if (shape==1) then
        size1=id.var(i).dims(1)
        allocate(aux1(size1))
        if (size1==nod2D)  call gather_nod (id.var(i).pt1, aux1)
        if (size1==elem2D) call gather_elem(id.var(i).pt1, aux1)
        if (mype==0) then
           error_status(c)=nf_put_vara_double(id.ncid, id.var(i).code, (/1, id.rec_count/), (/size1, 1/), aux1, 2); c=c+1
        end if
        deallocate(aux1)
!_______writing 3D fields________________________________________________
     elseif (shape==2) then
        size1=id.var(i).dims(1)
        size2=id.var(i).dims(2)
        allocate(aux2(size1, size2))
        if (size1==nod2D  .or. size2==nod2D)  call gather_nod (id.var(i).pt2, aux2)
        if (size1==elem2D .or. size2==elem2D) call gather_elem(id.var(i).pt2, aux2)
        if (mype==0) then
           error_status(c)=nf_put_vara_double(id.ncid, id.var(i).code, (/1, 1, id.rec_count/), (/size1, size2, 1/), aux2, 3); c=c+1
        end if
        deallocate(aux2)
     else
        if (mype==0) write(*,*) 'not supported shape of array in restart file'
           call par_ex
           stop
     end if
  end do
  error_count=c-1
  call was_error
  if (mype==0) error_status(1)=nf_close(id.ncid);
  error_count=1
  call was_error

end subroutine write_restart


END MODULE io_RESTART
