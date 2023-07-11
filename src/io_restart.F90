MODULE io_RESTART
  use g_config
  use g_clock
  use g_parsup
  use g_comm_auto
  use mod_mesh
  use o_arrays
  use i_arrays
  use g_cvmix_tke
  use g_cvmix_idemix
#if defined(__recom)
  use REcoM_GloVar
  use recom_config
  use REcoM_ciso
#endif
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
    integer :: rec, Tid, Iid
    integer :: ncid
    integer :: rec_count=0
    integer :: error_status(2000), error_count
    logical :: is_in_use=.false.
  end type nc_file
!
!--------------------------------------------------------------------------------------------
!
  type type_id
    integer :: nd, el, nz, nz1, T, rec, iter
  end type type_id
!
!--------------------------------------------------------------------------------------------
! id will keep the IDs of all required dimentions and variables
  type(nc_file), save       :: oid, iid, bid
  integer,       save       :: globalstep=0
  type(nc_file), save       :: ip_id
  real(kind=WP)             :: ctime !current time in seconds from the beginning of the year

  PRIVATE
  PUBLIC :: restart, oid, iid, bid
  PUBLIC :: ip_id, def_dim, def_variable_1d, def_variable_2d 
!
!--------------------------------------------------------------------------------------------
! generic interface was required to associate variables of unknown rank with the pointers of the same rank
! this allows for automatic streaming of associated variables into the netcdf file
  INTERFACE def_variable
            MODULE PROCEDURE def_variable_1d, def_variable_2d
  END INTERFACE
!
!--------------------------------------------------------------------------------------------
!
  contains
!
!--------------------------------------------------------------------------------------------
! ini_ocean_io initializes oid datatype which contains information of all variables need to be written into 
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
  ! create an ocean restart file; serial output implemented so far
  oid%filename=trim(ResultPath)//trim(runid)//'.'//cyear//'.oce.restart.nc'
  if (oid%is_in_use) return
  oid%is_in_use=.true.
  call def_dim(oid, 'node', nod2d)
  call def_dim(oid, 'elem', elem2d)
  call def_dim(oid, 'nz_1', nl-1)
  call def_dim(oid, 'nz',   nl)

  !===========================================================================
  !===================== Definition part =====================================
  !===========================================================================
  !___Define the netCDF variables for 2D fields_______________________________
  !___SSH_____________________________________________________________________
  call def_variable(oid, 'ssh',             (/nod2D/), 'sea surface elevation', 'm',   eta_n);
  !___ALE related fields______________________________________________________
  call def_variable(oid, 'hbar',            (/nod2D/), 'ALE surface elevation', 'm',   hbar);
!!PS   call def_variable(oid, 'ssh_rhs',         (/nod2D/), 'RHS for the elevation', '?',   ssh_rhs);
  call def_variable(oid, 'ssh_rhs_old',     (/nod2D/), 'RHS for the elevation', '?',   ssh_rhs_old);
  call def_variable(oid, 'hnode',    (/nl-1,  nod2D/), 'nodal layer thickness', 'm',   hnode);
  
  !___Define the netCDF variables for 3D fields_______________________________
  call def_variable(oid, 'u',        (/nl-1, elem2D/), 'zonal velocity',        'm/s', UV(1,:,:));
  call def_variable(oid, 'v',        (/nl-1, elem2D/), 'meridional velocity',   'm/s', UV(2,:,:));
  call def_variable(oid, 'urhs_AB',  (/nl-1, elem2D/), 'Adams–Bashforth for u', 'm/s', UV_rhsAB(1,:,:));
  call def_variable(oid, 'vrhs_AB',  (/nl-1, elem2D/), 'Adams–Bashforth for v', 'm/s', UV_rhsAB(2,:,:));
  
  !___Save restart variables for TKE and IDEMIX_________________________________
  if (trim(mix_scheme)=='cvmix_TKE' .or. trim(mix_scheme)=='cvmix_TKE+IDEMIX') then
        call def_variable(oid, 'tke',  (/nl, nod2d/), 'Turbulent Kinetic Energy', 'm2/s2', tke(:,:));
  endif
  if (trim(mix_scheme)=='cvmix_IDEMIX' .or. trim(mix_scheme)=='cvmix_TKE+IDEMIX') then
        call def_variable(oid, 'iwe',  (/nl, nod2d/), 'Internal Wave eneryy', 'm2/s2', tke(:,:));
  endif 
  if (visc_option==8) then
        call def_variable(oid, 'uke',      (/nl-1, elem2D/), 'unresolved kinetic energy', 'm2/s2', uke(:,:));
        call def_variable(oid, 'uke_rhs',  (/nl-1, elem2D/), 'unresolved kinetic energy rhs', 'm2/s2', uke_rhs(:,:));
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
!         write(trname,'(A3,i1)') 'tra_', j
!         write(longname,'(A15,i1)') 'passive tracer ', j
  	 write(trname,'(A3,i4.4)') 'tra', j		 ! OG i1 -> i4
         write(longname,'(A15,i4.4)') 'passive tracer ', j
         units='none'
     END SELECT
     call def_variable(oid, trim(trname),       (/nl-1, nod2D/), trim(longname), trim(units), tr_arr(:,:,j));
     longname=trim(longname)//', Adams–Bashforth'
     call def_variable(oid, trim(trname)//'_AB',(/nl-1, nod2D/), trim(longname), trim(units), tr_arr_old(:,:,j));
  end do
  call def_variable(oid, 'w',      (/nl, nod2D/), 'vertical velocity', 'm/s', Wvel);
  call def_variable(oid, 'w_expl', (/nl, nod2D/), 'vertical velocity', 'm/s', Wvel_e);
  call def_variable(oid, 'w_impl', (/nl, nod2D/), 'vertical velocity', 'm/s', Wvel_i);
end subroutine ini_ocean_io
!
!--------------------------------------------------------------------------------------------
! ini_ice_io initializes iid datatype which contains information of all variables need to be written into 
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
  ! create an ocean restart file; serial output implemented so far
  iid%filename=trim(ResultPath)//trim(runid)//'.'//cyear//'.ice.restart.nc'
  if (iid%is_in_use) return
  iid%is_in_use=.true.
  call def_dim(iid, 'node', nod2d)

  !===========================================================================
  !===================== Definition part =====================================
  !===========================================================================
  !___Define the netCDF variables for 2D fields_______________________________
  call def_variable(iid, 'area',       (/nod2D/), 'ice concentration [0 to 1]', '%',   a_ice);
  call def_variable(iid, 'hice',       (/nod2D/), 'effective ice thickness',    'm',   m_ice);
  call def_variable(iid, 'hsnow',      (/nod2D/), 'effective snow thickness',   'm',   m_snow);
  call def_variable(iid, 'uice',       (/nod2D/), 'zonal velocity',             'm/s', u_ice);
  call def_variable(iid, 'vice',       (/nod2D/), 'meridional velocity',        'm',   v_ice);
#if defined (__oifs)
  call def_variable(iid, 'ice_albedo', (/nod2D/), 'ice albedo',                 '-',   ice_alb);
  call def_variable(iid, 'ice_temp',(/nod2D/), 'ice surface temperature',  'K',   ice_temp);
#endif /* (__oifs) */

end subroutine ini_ice_io
#if defined(__recom)
!
!--------------------------------------------------------------------------------------------
! ini_bio_io initializes bid datatype which contains information of all variables need to be written into 
! the bio restart file. This is the only place need to be modified if a new variable is added!
subroutine ini_bio_io(year, mesh)
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
  ! create an ocean restart file; serial output implemented so far
  bid%filename=trim(ResultPath)//trim(runid)//'.'//cyear//'.bio.restart.nc'
  if (bid%is_in_use) return
  bid%is_in_use=.true.
  call def_dim(bid, 'node', nod2d)

  !===========================================================================
  !===================== Definition part =====================================
  !===========================================================================
  !___Define the netCDF variables for 2D fields_______________________________
  call def_variable(bid, 'BenN',       (/nod2D/), 'Benthos Nitrogen', 'mmol/m3',   Benthos(:,1));
  call def_variable(bid, 'BenC',       (/nod2D/), 'Benthos Carbon',   'mmol/m3',   Benthos(:,2));
  call def_variable(bid, 'BenSi',      (/nod2D/), 'Benthos Silicate', 'mmol/m3',   Benthos(:,3));
  call def_variable(bid, 'BenCalc',    (/nod2D/), 'Benthos Calcite',  'mmol/m3',   Benthos(:,4));
  call def_variable(bid, 'HPlus',      (/nod2D/), 'Conc. of H-plus ions in the surface water', 'mol/kg',   GloHplus);
  if (ciso) then
    call def_variable(bid, 'BenC_13',       (/nod2D/), 'Benthos Carbon-13',   'mmol/m3',   Benthos(:,5));
    call def_variable(bid, 'BenC_14',       (/nod2D/), 'Benthos Carbon-14',   'mmol/m3',   Benthos(:,6));
    call def_variable(bid, 'BenCalc_13',    (/nod2D/), 'Benthos Calcite-13',  'mmol/m3',   Benthos(:,7)); 
    call def_variable(bid, 'BenCalc_14',    (/nod2D/), 'Benthos Calcite-14',  'mmol/m3',   Benthos(:,8)); 
  end if
end subroutine ini_bio_io
#endif
!
!--------------------------------------------------------------------------------------------
!
subroutine restart(istep, l_write, l_read, mesh)

#if defined(__icepack)
  use icedrv_main,   only: init_restart_icepack
#endif

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
#if defined(__recom)
  if (use_REcoM) then
               call ini_bio_io  (yearnew, mesh)
  end if 
#endif
#if defined(__icepack)
  if (use_ice) call init_restart_icepack(yearnew, mesh)
#endif
  else
               call ini_ocean_io(yearold, mesh)
  if (use_ice) call ini_ice_io  (yearold, mesh)
#if defined(__recom)
  if (REcoM_restart) then
               call ini_bio_io(yearold, mesh)
  end if 
#endif
#if defined(__icepack)
  if (use_ice) call init_restart_icepack(yearold, mesh)
#endif
  end if

  if (l_read) then
   call assoc_ids(oid);          call was_error(oid)
   call read_restart(oid, mesh); call was_error(oid)
   if (use_ice) then
      call assoc_ids(iid);          call was_error(iid)
      call read_restart(iid, mesh); call was_error(iid)
#if defined(__icepack)
      call assoc_ids(ip_id);          call was_error(ip_id)
      call read_restart(ip_id, mesh); call was_error(ip_id)
#endif
    end if
#if defined(__recom)
!RECOM restart
!read here
if(mype==0)  write(*,*) 'REcoM_restart= ',REcoM_restart 
   if (REcoM_restart) then
      call assoc_ids(bid);          call was_error(bid)
      call read_restart(bid, mesh); call was_error(bid)
   end if 
#endif
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
  call assoc_ids(oid);                  call was_error(oid)  
  call write_restart(oid, istep, mesh); call was_error(oid)
  if (use_ice) then
     call assoc_ids(iid);                  call was_error(iid)  
     call write_restart(iid, istep, mesh); call was_error(iid)
#if defined(__icepack)
     call assoc_ids(ip_id);                  call was_error(ip_id)
     call write_restart(ip_id, istep, mesh); call was_error(ip_id)
#endif
  end if
#if defined(__recom)
!RECOM restart
!write here
   if (REcoM_restart .or. use_REcoM) then
     call assoc_ids(bid);                  call was_error(bid)
     call write_restart(bid, istep, mesh); call was_error(bid)
   end if
#endif     
  ! actualize clock file to latest restart point
  if (mype==0) then
		write(*,*) ' --> actualize clock file to latest restart point'
		call clock_finish  
  end if
  
end subroutine restart
!
!--------------------------------------------------------------------------------------------
!
subroutine create_new_file(id)
  implicit none

  type(nc_file),  intent(inout) :: id
  integer                       :: c, j
  integer                       :: n, k, l, kdim, dimid(4)
  character(2000)               :: att_text
  ! Serial output implemented so far
  if (mype/=0) return
  c=1
  id%error_status=0
  ! create an ocean output file
  write(*,*) 'initializing restart file ', trim(id%filename)
  id%error_status(c) = nf_create(id%filename, IOR(NF_NOCLOBBER,IOR(NF_NETCDF4,NF_CLASSIC_MODEL)), id%ncid); c=c+1

  do j=1, id%ndim
!___Create mesh related dimentions__________________________________________
     id%error_status(c) = nf_def_dim(id%ncid, id%dim(j)%name, id%dim(j)%size, id%dim(j)%code ); c=c+1
  end do

!___Create time related dimentions__________________________________________
  id%error_status(c) = nf_def_dim(id%ncid, 'time', NF_UNLIMITED, id%rec);         c=c+1
!___Define the time and iteration variables_________________________________
  id%error_status(c) = nf_def_var(id%ncid, 'time', NF_DOUBLE, 1, id%rec, id%tID); c=c+1
  id%error_status(c) = nf_def_var(id%ncid, 'iter', NF_INT,    1, id%rec, id%iID); c=c+1


  att_text='time'
  id%error_status(c) = nf_put_att_text(id%ncid, id%tID, 'long_name', len_trim(att_text), trim(att_text)); c=c+1
  write(att_text, '(a14,I4.4,a1,I2.2,a1,I2.2,a6)') 'seconds since ', yearold, '-', 1, '-', 1, ' 0:0:0'
  id%error_status(c) = nf_put_att_text(id%ncid, id%tID, 'units', len_trim(att_text), trim(att_text)); c=c+1

  att_text='iteration_count'
  id%error_status(c) = nf_put_att_text(id%ncid, id%iID, 'long_name', len_trim(att_text), trim(att_text)); c=c+1

  do j=1, id%nvar
!___associate physical dimension with the netcdf IDs________________________
     n=id%var(j)%ndim ! shape size of the variable (exluding time)
     do k=1, n
        !k_th dimension of the variable
        kdim=id%var(j)%dims(k)
        do l=1, id%ndim ! list all defined dimensions 
           if (kdim==id%dim(l)%size) dimid(k)=id%dim(l)%code
        end do
        !write(*,*) "j",j,kdim, ' -> ', dimid(k)
     end do
     id%error_status(c) = nf_def_var(id%ncid, trim(id%var(j)%name), NF_DOUBLE, id%var(j)%ndim+1, (/dimid(1:n), id%rec/), id%var(j)%code); c=c+1
     !if (n==1) then
     !   id%error_status(c)=nf_def_var_chunking(id%ncid, id%var(j)%code, NF_CHUNKED, (/1/)); c=c+1 
     if (n==2) then
        id%error_status(c)=nf_def_var_chunking(id%ncid, id%var(j)%code, NF_CHUNKED, (/1, id%dim(1)%size/)); ! c=c+1 
     end if
     id%error_status(c)=nf_put_att_text(id%ncid, id%var(j)%code, 'description', len_trim(id%var(j)%longname), id%var(j)%longname); c=c+1
     id%error_status(c)=nf_put_att_text(id%ncid, id%var(j)%code, 'units',       len_trim(id%var(j)%units),    id%var(j)%units);    c=c+1
  end do

  id%error_status(c)=nf_close(id%ncid); c=c+1
  id%error_count=c-1
end subroutine create_new_file
!
!--------------------------------------------------------------------------------------------
!
subroutine def_dim(id, name, ndim)
  implicit none
  type(nc_file),    intent(inout) :: id
  character(len=*), intent(in)    :: name
  integer,          intent(in)    :: ndim
  type(nc_dims), allocatable, dimension(:) :: temp

  if (id%ndim > 0) then
     ! create temporal dimension
     allocate(temp(id%ndim)); temp=id%dim
     ! deallocate the input data array
     deallocate(id%dim)
     ! then reallocate
     id%ndim=id%ndim+1
     allocate(id%dim(id%ndim))
     ! restore the original data
     id%dim(1:id%ndim-1)=temp  
     deallocate(temp)
   else
     ! first dimension in a file
     id%ndim=1
     allocate(id%dim(id%ndim))
   end if
   id%dim(id%ndim)%name=trim(name)
   id%dim(id%ndim)%size=ndim
end subroutine def_dim
!
!--------------------------------------------------------------------------------------------
!
subroutine def_variable_1d(id, name, dims, longname, units, data)
  implicit none
  type(nc_file),    intent(inout)        :: id
  character(len=*), intent(in)           :: name
  integer, intent(in)                    :: dims(1)
  character(len=*), intent(in), optional :: units, longname
  real(kind=WP),target,     intent(inout)        :: data(:)
  integer                                :: c
  type(nc_vars), allocatable, dimension(:) :: temp

  if (id%nvar > 0) then
     ! create temporal dimension
     allocate(temp(id%nvar)); temp=id%var
     ! deallocate the input data array
     deallocate(id%var)
     ! then reallocate
     id%nvar=id%nvar+1
     allocate(id%var(id%nvar))
     ! restore the original data
     id%var(1:id%nvar-1)=temp  
     deallocate(temp)
   else
     ! first dimension in a file
     id%nvar=1
     allocate(id%var(id%nvar))
   end if
   id%var(id%nvar)%name=trim(name)
   id%var(id%nvar)%longname=trim(longname)
   id%var(id%nvar)%units=trim(units)
   id%var(id%nvar)%ndim=1
   id%var(id%nvar)%dims(1)=dims(1)
   id%var(id%nvar)%pt1=>data
end subroutine def_variable_1d
!
!--------------------------------------------------------------------------------------------
!
subroutine def_variable_2d(id, name, dims, longname, units, data)
  implicit none
  type(nc_file),    intent(inout)        :: id
  character(len=*), intent(in)           :: name
  integer, intent(in)                    :: dims(2)
  character(len=*), intent(in), optional :: units, longname
  real(kind=WP),target,     intent(inout) :: data(:,:)
  integer                                :: c
  type(nc_vars), allocatable, dimension(:) :: temp

  if (id%nvar > 0) then
     ! create temporal dimension
     allocate(temp(id%nvar)); temp=id%var
     ! deallocate the input data array
     deallocate(id%var)
     ! then reallocate
     id%nvar=id%nvar+1
     allocate(id%var(id%nvar))
     ! restore the original data
     id%var(1:id%nvar-1)=temp  
     deallocate(temp)
   else
     ! first dimension in a file
     id%nvar=1
     allocate(id%var(id%nvar))
   end if
   id%var(id%nvar)%name=trim(name)
   id%var(id%nvar)%longname=trim(longname)
   id%var(id%nvar)%units=trim(units)
   id%var(id%nvar)%ndim=2
   id%var(id%nvar)%dims(1:2)=dims
   id%var(id%nvar)%pt2=>data
end subroutine def_variable_2d
!
!--------------------------------------------------------------------------------------------
!
subroutine write_restart(id, istep, mesh)
  implicit none
  type(nc_file),  intent(inout) :: id
  integer,  intent(in)          :: istep
  type(t_mesh), intent(in)     , target :: mesh
  real(kind=WP), allocatable    :: aux(:), laux(:)
  real(kind=WP)                 :: t0, t1, t2, t3
  integer                       :: i, lev, size1, size2, size_gen, size_lev, shape
  integer                       :: c, order

#include  "associate_mesh.h"

  ! Serial output implemented so far
  if (mype==0) then
     c=1
     !id%rec_count=id%rec_count+1
     write(*,*) 'writing restart record ', id%rec_count
     id%error_status(c)=nf_open(id%filename, nf_write, id%ncid); c=c+1
     id%error_status(c)=nf_put_vara_double(id%ncid, id%tID, id%rec_count, 1, ctime, 1); c=c+1
     id%error_status(c)=nf_put_vara_int(id%ncid,    id%iID, id%rec_count, 1, globalstep+istep, 1);   c=c+1
  end if

  call was_error(id); c=1

  do i=1, id%nvar
     shape=id%var(i)%ndim
!_______writing 2D fields________________________________________________
     if (shape==1) then
        size1=id%var(i)%dims(1)
        if (mype==0) allocate(aux(size1))
        t0=MPI_Wtime()
        if (size1==nod2D)  call gather_nod (id%var(i)%pt1, aux)
        if (size1==elem2D) call gather_elem(id%var(i)%pt1, aux)
        t1=MPI_Wtime()
        if (mype==0) then
           id%error_status(c)=nf_put_vara_double(id%ncid, id%var(i)%code, (/1, id%rec_count/), (/size1, 1/), aux, 1); c=c+1
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
        size1=id%var(i)%dims(1)
        size2=id%var(i)%dims(2)
        ! I assume that the model has always more surface nodes or elements than
        ! vertical layers => more flexibility in terms of array dimensions
        if (size1==nod2D .or. size1==elem2D) then
            size_gen=size1
            size_lev=size2
            order=1
        else if (size2==nod2D .or. size2==elem2D) then
            size_gen=size2
            size_lev=size1
            order=2
        else
            if (mype==0) write(*,*) 'the shape of the array in the restart file and the grid size are different'
            call par_ex
            stop
        end if 
        if (mype==0)          allocate(aux (size_gen))
        if (size_gen==nod2D)  allocate(laux(myDim_nod2D +eDim_nod2D ))
        if (size_gen==elem2D) allocate(laux(myDim_elem2D+eDim_elem2D))
        do lev=1, size_lev
           if (order==1) laux=id%var(i)%pt2(:,lev)
           if (order==2) laux=id%var(i)%pt2(lev,:)
           if (size_gen==nod2D)  call gather_nod (laux, aux)
           if (size_gen==elem2D) call gather_elem(laux, aux)
           if (mype==0) then
              if (order==1) id%error_status(c)=nf_put_vara_double(id%ncid, id%var(i)%code, (/1, lev, id%rec_count/), (/size_gen, 1, 1/), aux, 1); c=c+1
              if (order==2) id%error_status(c)=nf_put_vara_double(id%ncid, id%var(i)%code, (/lev, 1, id%rec_count/), (/1, size_gen, 1/), aux, 1); c=c+1
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
     call was_error(id); c=1
  end do

  if (mype==0) id%error_count=c-1
  call was_error(id)
  if (mype==0) id%error_status(1)=nf_close(id%ncid);
  id%error_count=1
  call was_error(id)
end subroutine write_restart
!
!--------------------------------------------------------------------------------------------
!
subroutine read_restart(id, mesh, arg)
  implicit none
  type(nc_file),     intent(inout) :: id
  integer, optional, intent(in)    :: arg
  real(kind=WP), allocatable       :: aux(:), laux(:)
  integer                          :: i, lev, size1, size2, size_gen, size_lev, shape
  integer                          :: rec2read, c, order
  real(kind=WP)                    :: rtime !timestamp of the record
  logical                          :: file_exist=.False.
  type(t_mesh), intent(in)        , target :: mesh

#include  "associate_mesh.h"

  ! laux=0.
  ! Serial output implemented so far
  c=1
  if (mype==0) then
     file_exist=.False.
     inquire(file=id%filename,exist=file_exist) 
     if (file_exist) then
        write(*,*) '     reading restart file:  ', trim(id%filename)
        id%error_status(c)=nf_open(id%filename, nf_nowrite, id%ncid);                           c=c+1
        id%error_status(c)=nf_get_vara_int(id%ncid,    id%iID, id%rec_count, 1, globalstep, 1); c=c+1
        id%error_status(c)=nf_get_vara_double(id%ncid, id%tID, id%rec_count, 1, rtime, 1);      c=c+1
     else
        write(*,*)
        print *, achar(27)//'[33m'
        write(*,*) '____________________________________________________________________'
        write(*,*) ' ERROR: could not find restart_file:',trim(id%filename),'!'    
        write(*,*) '____________________________________________________________________'
        print *, achar(27)//'[0m'
        write(*,*)
        call par_ex
     end if 
     
     if (.not. present(arg)) then
        rec2read=id%rec_count
     else
        rec2read=arg
     end if
     write(*,*) 'restart from record ', rec2read, ' of ', id%rec_count

     if (int(ctime)/=int(rtime)) then
        write(*,*) 'Reading restart: timestamps in restart and in clock files do not match'
        write(*,*) 'restart/ times are:', ctime, rtime
        write(*,*) 'the model will stop!'
        id%error_status(c)=-310; c=c+1
     end if
  end if

  call was_error(id); c=1
 
  do i=1, id%nvar
     shape=id%var(i)%ndim
     if (mype==0) write(*,*) 'reading restart for ', trim(id%var(i)%name)
!_______writing 2D fields________________________________________________
     if (shape==1) then
        size1=id%var(i)%dims(1)
        if (mype==0) then
           allocate(aux(size1))
           id%error_status(c)=nf_get_vara_double(id%ncid, id%var(i)%code, (/1, id%rec_count/), (/size1, 1/), aux, 1); c=c+1
!          write(*,*) 'min/max 2D =', minval(aux), maxval(aux)
        end if
        if (size1==nod2D)  call broadcast_nod (id%var(i)%pt1, aux)
        if (size1==elem2D) call broadcast_elem(id%var(i)%pt1, aux)
        if (mype==0) deallocate(aux)
!_______writing 3D fields________________________________________________
     elseif (shape==2) then
        size1=id%var(i)%dims(1)
        size2=id%var(i)%dims(2)
        ! I assume that the model has always more surface nodes or elements than
        ! vertical layers => more flexibility in terms of array dimensions
        if (size1==nod2D .or. size1==elem2D) then
            size_gen=size1
            size_lev=size2
            order=1
        else if (size2==nod2D .or. size2==elem2D) then
            size_gen=size2
            size_lev=size1
            order=2
        else
            if (mype==0) write(*,*) 'the shape of the array in the restart file and the grid size are different'
            call par_ex
            stop
        end if
        if (mype==0)          allocate(aux (size_gen))
        if (size_gen==nod2D)  allocate(laux(myDim_nod2D +eDim_nod2D ))
        if (size_gen==elem2D) allocate(laux(myDim_elem2D+eDim_elem2D))
        do lev=1, size_lev
           if (mype==0) then
              if (order==1) id%error_status(c)=nf_get_vara_double(id%ncid, id%var(i)%code, (/1, lev, id%rec_count/), (/size_gen, 1, 1/), aux, 1); c=c+1
              if (order==2) id%error_status(c)=nf_get_vara_double(id%ncid, id%var(i)%code, (/lev, 1, id%rec_count/), (/1, size_gen, 1/), aux, 1); c=c+1
           end if
           id%var(i)%pt2(lev,:)=0.
           if (size_gen==nod2D)  then
              call broadcast_nod (laux, aux)
              if (order==1) id%var(i)%pt2(:,lev)=laux(1:myDim_nod2D+eDim_nod2D)
              if (order==2) id%var(i)%pt2(lev,:)=laux(1:myDim_nod2D+eDim_nod2D)
           end if
           if (size_gen==elem2D) then
              call broadcast_elem(laux, aux)
              if (order==1) id%var(i)%pt2(:,lev)=laux(1:myDim_elem2D+eDim_elem2D)
              if (order==2) id%var(i)%pt2(lev,:)=laux(1:myDim_elem2D+eDim_elem2D)
           end if
        end do
        deallocate(laux)
        if (mype==0) deallocate(aux)
     else
        if (mype==0) write(*,*) 'not supported shape of array in restart file when reading restart'
           call par_ex
           stop
     end if
     call was_error(id); c=1
  end do

  if (mype==0) id%error_status(1)=nf_close(id%ncid);
  id%error_count=1
  call was_error(id)
end subroutine read_restart
!
!--------------------------------------------------------------------------------------------
!
subroutine assoc_ids(id)
  implicit none

  type(nc_file),  intent(inout) :: id
  character(500)                :: longname
  integer                       :: c, j, k
  real(kind=WP)                 :: rtime !timestamp of the record
  ! Serial output implemented so far
  if (mype/=0) return
  c=1
  id%error_status=0
  ! open existing netcdf file
  write(*,*) 'associating restart file ', trim(id%filename)

  id%error_status(c) = nf_open(id%filename, nf_nowrite, id%ncid)
  !if the file does not exist it will be created!
  if (id%error_status(c) .ne. nf_noerr) then
     call create_new_file(id) ! error status counter will be reset
     c=id%error_count+1
     id%error_status(c) = nf_open(id%filename, nf_nowrite, id%ncid); c=c+1
  end if

  do j=1, id%ndim
!___Associate mesh related dimentions_______________________________________
    id%error_status(c) = nf_inq_dimid(id%ncid, id%dim(j)%name, id%dim(j)%code); c=c+1
  end do
!___Associate time related dimentions_______________________________________
  id%error_status(c) = nf_inq_dimid (id%ncid, 'time', id%rec);       c=c+1
  id%error_status(c) = nf_inq_dimlen(id%ncid, id%rec, id%rec_count); c=c+1
!___Associate the time and iteration variables______________________________
  id%error_status(c) = nf_inq_varid(id%ncid, 'time', id%tID); c=c+1
  id%error_status(c) = nf_inq_varid(id%ncid, 'iter', id%iID); c=c+1
!___if the time rtime at the rec_count does not equal ctime we look for the closest record with the 
! timestamp less than ctime
  do k=id%rec_count, 1, -1
     id%error_status(c)=nf_get_vara_double(id%ncid, id%tID, k, 1, rtime, 1);
     if (ctime > rtime) then
        id%rec_count=k+1
        exit ! a proper rec_count detected, ready for writing restart, exit the loop
     elseif (ctime == rtime) then
        id%rec_count=k
        exit ! a proper rec_count detected, ready for reading restart, exit the loop
     end if
     if (k==1) then
        if (mype==0) write(*,*) 'WARNING: all dates in restart file are after the current date'
        if (mype==0) write(*,*) 'reading restart will not be possible !'
        if (mype==0) write(*,*) 'the model attempted to start with the time stamp = ', int(ctime)
        id%error_status(c)=-310;
     end if
  end do
  c=c+1 ! check will be made only for the last nf_get_vara_double
  id%rec_count=max(id%rec_count, 1)
!___Associate physical variables____________________________________________
  do j=1, id%nvar
     id%error_status(c) = nf_inq_varid(id%ncid, id%var(j)%name, id%var(j)%code); c=c+1
  end do
  id%error_status(c)=nf_close(id%ncid); c=c+1
  id%error_count=c-1
  write(*,*) 'current restart counter = ',       id%rec_count
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
