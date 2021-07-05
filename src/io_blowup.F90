MODULE io_BLOWUP
	use g_config
	use g_clock
	use g_parsup
	use g_comm_auto
        USE MOD_MESH
	use o_arrays
	use i_arrays
	implicit none
#include "netcdf.inc"
	!___________________________________________________________________________
	type nc_dims
		integer        :: size
		character(100) :: name
		integer        :: code
	end type nc_dims
	!___________________________________________________________________________
	type nc_vars
		character(100) :: name
		integer        :: code
		character(500) :: longname
		character(100) :: units
		integer        :: ndim
		integer        :: dims(2) !<=2; assume there are no variables with dimension more than 2xNLxT
		real(kind=WP), pointer :: pt1(:), pt2(:,:)
	end type nc_vars
	!___________________________________________________________________________
	type nc_file
		character(500)                                :: filename
		type(nc_dims), allocatable, dimension(:) :: dim
		type(nc_vars), allocatable, dimension(:) :: var
		integer :: ndim=0, nvar=0
		integer :: rec, Tid, Iid
		integer :: ncid
		integer :: rec_count=0
		integer :: error_status(200), error_count
		logical :: is_in_use=.false.
	end type nc_file
	!___________________________________________________________________________
	type type_id
		integer :: nd, el, nz, nz1, T, rec, iter
	end type type_id
	!___________________________________________________________________________
	! id will keep the IDs of all required dimentions and variables
	type(nc_file), save       :: bid
	integer,       save       :: globalstep=0
	real(kind=WP)             :: ctime !current time in seconds from the beginning of the year
	
	PRIVATE
	PUBLIC :: blowup, bid 
	
	!___________________________________________________________________________
	! generic interface was required to associate variables of unknown rank with the pointers of the same rank
	! this allows for automatic streaming of associated variables into the netcdf file
	INTERFACE def_variable
		MODULE PROCEDURE def_variable_1d, def_variable_2d
	END INTERFACE
	!___________________________________________________________________________
	contains
	!
	!
	!_______________________________________________________________________________
	! ini_ocean_io initializes bid datatype which contains information of all variables need to be written into 
	! the ocean restart file. This is the only place need to be modified if a new variable is added!
	subroutine ini_blowup_io(year, mesh)
		implicit none
                type(t_mesh), intent(in) , target :: mesh
		integer, intent(in)       :: year
		integer                   :: ncid, j
		integer                   :: varid
		character(500)            :: longname
		character(500)            :: filename
		character(500)            :: trname, units
		character(4)              :: cyear

#include  "associate_mesh.h"

		if(mype==0) write(*,*)' --> Init. blowpup file '
		write(cyear,'(i4)') year
		! create an ocean restart file; serial output implemented so far
		bid%filename=trim(ResultPath)//trim(runid)//'.'//cyear//'.oce.blowup.nc'
		if (bid%is_in_use) return
		bid%is_in_use=.true.
		call def_dim(bid, 'node', nod2d)
		call def_dim(bid, 'elem', elem2d)
		call def_dim(bid, 'nz_1', nl-1)
		call def_dim(bid, 'nz',   nl)
		
		!===========================================================================
		!===================== Definition part =====================================
		!===========================================================================
		!___Define the netCDF variables for 2D fields_______________________________
		!___SSH_____________________________________________________________________
		call def_variable(bid, 'eta_n'		, (/nod2D/)			, 'sea surface elevation', 'm', eta_n);
		call def_variable(bid, 'd_eta'		, (/nod2D/)			, 'change in ssh from solver', 'm', d_eta);
		!___ALE related fields______________________________________________________
		call def_variable(bid, 'hbar'		, (/nod2D/)			, 'ALE surface elevation hbar_n+0.5', 'm', hbar);
!!PS 		call def_variable(bid, 'hbar_old'	, (/nod2D/)			, 'ALE surface elevation hbar_n-0.5', 'm', hbar_old);
		call def_variable(bid, 'ssh_rhs'	, (/nod2D/)			, 'RHS for the elevation', '?', ssh_rhs);
		call def_variable(bid, 'ssh_rhs_old', (/nod2D/)			, 'RHS for the elevation', '?', ssh_rhs_old);
		!___Define the netCDF variables for 3D fields_______________________________
		call def_variable(bid, 'hnode'		, (/nl-1,  nod2D/)	, 'ALE stuff', '?', hnode);
		call def_variable(bid, 'helem'		, (/nl-1, elem2D/)	, 'Element layer thickness', 'm/s', helem(:,:));
		call def_variable(bid, 'u'			, (/nl-1, elem2D/)	, 'zonal velocity', 'm/s', UV(1,:,:));
		call def_variable(bid, 'v'			, (/nl-1, elem2D/)	, 'meridional velocity', 'm/s', UV(2,:,:));
		call def_variable(bid, 'u_rhs'			, (/nl-1, elem2D/)	, 'zonal velocity', 'm/s', UV_rhs(1,:,:));
		call def_variable(bid, 'v_rhs'			, (/nl-1, elem2D/)	, 'meridional velocity', 'm/s', UV_rhs(2,:,:));
		call def_variable(bid, 'urhs_AB'	, (/nl-1, elem2D/)	, 'Adams–Bashforth for u', 'm/s', UV_rhsAB(1,:,:));
		call def_variable(bid, 'vrhs_AB'	, (/nl-1, elem2D/)	, 'Adams–Bashforth for v', 'm/s', UV_rhsAB(2,:,:));
		call def_variable(bid, 'zbar_n_bot' , (/nod2D/)			, 'node bottom depth', 'm', zbar_n_bot);
		call def_variable(bid, 'zbar_e_bot' , (/elem2d/)		, 'elem bottom depth', 'm', zbar_e_bot);
		call def_variable(bid, 'bottom_node_thickness' , (/nod2D/)			, 'node bottom thickness', 'm', bottom_node_thickness);
		call def_variable(bid, 'bottom_elem_thickness' , (/elem2d/)		, 'elem bottom thickness', 'm', bottom_elem_thickness);
!!PS 		call def_variable(bid, 'pgf_x'	, (/nl-1, elem2D/)	, 'zonal pressure gradient force', '???', pgf_x(:,:));
!!PS 		call def_variable(bid, 'pgf_y'	, (/nl-1, elem2D/)	, 'meridional pressure gradient force', '???', pgf_y(:,:));
!!PS 		call def_variable(bid, 'density_m_rho0'	, (/nl-1, nod2D/)	, 'density minus rho0', '???', density_m_rho0(:,:));
		
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
			call def_variable(bid, trim(trname),       (/nl-1, nod2D/), trim(longname), trim(units), tr_arr(:,:,j));
!!PS 			longname=trim(longname)//', Adams–Bashforth'
!!PS 			call def_variable(bid, trim(trname)//'_AB',(/nl-1, nod2D/), trim(longname), trim(units), tr_arr_old(:,:,j));
		end do
		call def_variable(bid, 'w'			, (/nl, nod2D/)		, 'vertical velocity', 'm/s', Wvel);
		call def_variable(bid, 'w_expl'		, (/nl, nod2D/)		, 'vertical velocity', 'm/s', Wvel_e);
		call def_variable(bid, 'w_impl'		, (/nl, nod2D/)		, 'vertical velocity', 'm/s', Wvel_i);
		call def_variable(bid, 'cfl_z'		, (/nl-1, nod2D/)		, 'vertical CFL criteria', '', CFL_z);
		
		!_____________________________________________________________________________
		! write snapshot ice variables to blowup file
		call def_variable(bid, 'a_ice'		, (/nod2D/)			, 'ice concentration [0 to 1]', '%', a_ice);
		call def_variable(bid, 'm_ice'		, (/nod2D/)			, 'effective ice thickness',    'm', m_ice);
		call def_variable(bid, 'm_snow'		, (/nod2D/)			, 'effective snow thickness',   'm', m_snow);
		call def_variable(bid, 'u_ice'		, (/nod2D/)			, 'zonal velocity',    'm/s', u_ice);
		call def_variable(bid, 'v_ice'		, (/nod2D/)			, 'meridional velocity', 'm', v_ice);
!!PS  		call def_variable(bid, 'a_ice_old'	, (/nod2D/)			, 'ice concentration [0 to 1]', '%', a_ice_old); !PS
!!PS  		call def_variable(bid, 'm_ice_old'	, (/nod2D/)			, 'effective ice thickness',    'm', m_ice_old); !PS
!!PS  		call def_variable(bid, 'm_snow_old'	, (/nod2D/)			, 'effective snow thickness',   'm', m_snow_old); !PS
!!PS 		call def_variable(bid, 'u_ice_old'	, (/nod2D/)			, 'zonal velocity',    'm/s', u_ice_old);
!!PS  		call def_variable(bid, 'v_ice_old'	, (/nod2D/)			, 'meridional velocity', 'm', v_ice_old);
 		call def_variable(bid, 'heat_flux'	, (/nod2D/)			, 'heat flux ',     '?', heat_flux); !PS
		call def_variable(bid, 'water_flux'	, (/nod2D/)			, 'water flux ',    '?', water_flux); !PS
		call def_variable(bid, 'tx_sur'	    , (/elem2d/)		, 'tx_sur',         '?', stress_surf(1, :)); !PS
		call def_variable(bid, 'ty_sur'	    , (/elem2d/)		, 'ty_sur',         '?', stress_surf(2, :)); !PS
		
!!PS 		call def_variable(bid, 'fer_k'			, (/nl, nod2D/)		, 'GM diffusivity', '', fer_K);
		
		call def_variable(bid, 'Kv'			, (/nl, nod2D/)		, 'vertical diffusivity', '', Kv);
		call def_variable(bid, 'Av'			, (/nl, elem2d/)	, 'vertical viscosity', '', Av);
		call def_variable(bid, 'N2'			, (/nl, nod2D/)		, 'squared bouyancy freq', '', bvfreq);
		
	end subroutine ini_blowup_io
!
!
!_______________________________________________________________________________
	subroutine blowup(istep, mesh)
		implicit none
                type(t_mesh), intent(in) , target :: mesh		
		integer                   :: istep
		
		ctime=timeold+(dayold-1.)*86400
		call ini_blowup_io(yearnew, mesh)
		if(mype==0) write(*,*)'Do output (netCDF, blowup) ...'
		if(mype==0) write(*,*)' --> call assoc_ids(bid)'
		call assoc_ids(bid) ; call was_error(bid)  
		if(mype==0) write(*,*)' --> call write_blowup(bid, istep)'
		call write_blowup(bid, istep, mesh) ; call was_error(bid)
	
	end subroutine blowup
!
!
!_______________________________________________________________________________
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
		if(mype==0) write(*,*) 'initializing blowup file ', trim(id%filename)
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
		!________write(*,*) kdim, ' -> ', dimid(k)__________________________________
			end do
			id%error_status(c) = nf_def_var(id%ncid, trim(id%var(j)%name), NF_DOUBLE, id%var(j)%ndim+1, &
							(/dimid(1:n), id%rec/), id%var(j)%code); c=c+1
			id%error_status(c)=nf_put_att_text(id%ncid, id%var(j)%code, 'description', len_trim(id%var(j)%longname), id%var(j)%longname); c=c+1
			id%error_status(c)=nf_put_att_text(id%ncid, id%var(j)%code, 'units',       len_trim(id%var(j)%units),    id%var(j)%units);    c=c+1
		end do
		
		id%error_status(c)=nf_close(id%ncid); c=c+1
		id%error_count=c-1
	end subroutine create_new_file
!
!
!_______________________________________________________________________________
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
!
!_______________________________________________________________________________
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
!
!_______________________________________________________________________________
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
!
!_______________________________________________________________________________
	subroutine write_blowup(id, istep, mesh)
		implicit none
		type(nc_file),  intent(inout) :: id
		integer,  intent(in)          :: istep
		real(kind=WP), allocatable     :: aux1(:), aux2(:,:) 
		integer                       :: i, size1, size2, shape
		integer                       :: c
        type(t_mesh), intent(in)     , target :: mesh

#include  "associate_mesh.h"

		! Serial output implemented so far
		if (mype==0) then
			c=1
			!id%rec_count=id%rec_count+1
			write(*,*) 'writing blowup record ', id%rec_count
			id%error_status(c)=nf_open(id%filename, nf_write, id%ncid); c=c+1
			id%error_status(c)=nf_put_vara_double(id%ncid, id%tID, id%rec_count, 1, ctime, 1); c=c+1
			id%error_status(c)=nf_put_vara_int(id%ncid,    id%iID, id%rec_count, 1, globalstep+istep, 1);   c=c+1
		end if
		
		do i=1, id%nvar
			shape=id%var(i)%ndim
		!_______writing 2D fields________________________________________________
			if (shape==1) then
				size1=id%var(i)%dims(1)
				if (mype==0) allocate(aux1(size1))
				if (size1==nod2D)  call gather_nod (id%var(i)%pt1, aux1)
				if (size1==elem2D) call gather_elem(id%var(i)%pt1, aux1)
				if (mype==0) then
				id%error_status(c)=nf_put_vara_double(id%ncid, id%var(i)%code, (/1, id%rec_count/), (/size1, 1/), aux1, 1); c=c+1
				end if
				if (mype==0) deallocate(aux1)
		!_______writing 3D fields________________________________________________
			elseif (shape==2) then
				size1=id%var(i)%dims(1)
				size2=id%var(i)%dims(2)
				if (mype==0) allocate(aux2(size1, size2))
				if (size1==nod2D  .or. size2==nod2D)  call gather_nod (id%var(i)%pt2, aux2)
				if (size1==elem2D .or. size2==elem2D) call gather_elem(id%var(i)%pt2, aux2)
				if (mype==0) then
				id%error_status(c)=nf_put_vara_double(id%ncid, id%var(i)%code, (/1, 1, id%rec_count/), (/size1, size2, 1/), aux2, 2); c=c+1
				end if
				if (mype==0) deallocate(aux2)
			else
				if (mype==0) write(*,*) 'not supported shape of array in restart file'
				call par_ex
				stop
			end if
		end do
		
		if (mype==0) id%error_count=c-1
		call was_error(id)
		if (mype==0) id%error_status(1)=nf_close(id%ncid);
		id%error_count=1
		call was_error(id)
	end subroutine write_blowup
!
!
!_______________________________________________________________________________
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
		write(*,*) 'associating blowup file ', trim(id%filename)
		
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
!
!_______________________________________________________________________________
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
END MODULE io_BLOWUP
