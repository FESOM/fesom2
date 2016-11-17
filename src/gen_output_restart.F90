subroutine init_output_restart(do_init)
  use o_mesh
  use g_config
  use g_clock
  use g_parsup
  implicit none

#include "netcdf.inc"

  logical, intent(in)       :: do_init  
  integer                   :: status, ncid, dimid_rec, j
  integer                   :: dimid_2d, dimid_2de, dimid_nl1, dimid_nl
  integer                   :: dimids(2), dimid3(3)
  integer                   :: time_varid, iter_varid
  integer                   :: ssh_varid, tra_varid(num_tracer), tra_varid_ab(num_tracer)
  integer                   :: u_varid, v_varid, wpot_varid
  integer                   :: urhs_varid_AB, vrhs_varid_AB
  integer                   :: w_varid, we_varid, wi_varid
  integer                   :: area_varid, hice_varid, hsnow_varid
  integer                   :: uice_varid, vice_varid
  character(100)            :: longname
  character(100)            :: filename
  character(100)            :: trname, units

  if (.not. do_init) return
  ! Serial output implemented so far
  if (mype/=0) return

  ! create an ocean output file
  write(*,*) 'initialize new output files'
  filename=trim(ResultPath)//runid//'.'//cyearnew//'.oce.restart.nc'

  if (restart_offset==32) then
     status = nf_create(filename, nf_clobber, ncid)
  else
     status = nf_create(filename, IOR(NF_CLOBBER,NF_64BIT_OFFSET), ncid)
  end if

  if (status.ne.nf_noerr) call handle_err(status)

  ! Define the dimensions
  status = nf_def_dim(ncid, 'nodes_2d', nod2d, dimid_2d)
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_def_dim(ncid, 'elem_2d', elem2d, dimid_2de)  !!!! 
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_def_dim(ncid, 'nl_1', nl-1, dimid_nl1)
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_def_dim(ncid, 'nl', nl, dimid_nl)
  if (status .ne. nf_noerr) call handle_err(status)       !!!!

  status = nf_def_dim(ncid, 'T', NF_UNLIMITED, dimid_rec)
  if (status .ne. nf_noerr) call handle_err(status)
  ! Define the time and iteration variables
  status = nf_def_var(ncid, 'time', NF_DOUBLE, 1, dimid_rec, time_varid)
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_def_var(ncid, 'iter', NF_INT, 1, dimid_rec, iter_varid)
  if (status .ne. nf_noerr) call handle_err(status)

  ! Define the netCDF variables for 2D fields.
  ! In Fortran, the unlimited dimension must come
  ! last on the list of dimids.
  dimids(1) = dimid_2d
  dimids(2) = dimid_rec
  
  status = nf_def_var(ncid, 'ssh', NF_DOUBLE, 2, dimids, ssh_varid)
  if (status .ne. nf_noerr) call handle_err(status)

  ! Define the netCDF variables for 3D fields
  ! ocean horizontal velocity u, v
  dimid3(1) = dimid_nl1
  dimid3(2) = dimid_2de
  dimid3(3) = dimid_rec

  status = nf_def_var(ncid, 'u', NF_DOUBLE, 3, dimid3, u_varid)
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_def_var(ncid, 'v', NF_DOUBLE, 3, dimid3, v_varid)
  if (status .ne. nf_noerr) call handle_err(status)

  status = nf_def_var(ncid, 'urhs_AB', NF_DOUBLE, 3, dimid3, urhs_varid_AB)
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_def_var(ncid, 'vrhs_AB', NF_DOUBLE, 3, dimid3, vrhs_varid_AB)
  if (status .ne. nf_noerr) call handle_err(status)

     
  ! ocean vertical velocity w
  dimid3(1) = dimid_nl
  dimid3(2) = dimid_2d
  dimid3(3) = dimid_rec

  status = nf_def_var(ncid, 'w', NF_DOUBLE, 3, dimid3, w_varid)
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_def_var(ncid, 'w_expl', NF_DOUBLE, 3, dimid3, we_varid)
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_def_var(ncid, 'w_impl', NF_DOUBLE, 3, dimid3, wi_varid)
  if (status .ne. nf_noerr) call handle_err(status)

  ! scalar fields like T, S
  dimid3(1) = dimid_nl1
  dimid3(2) = dimid_2d
  dimid3(3) = dimid_rec

  do j=1,num_tracers
     SELECT CASE (j) 
       CASE(1)
         trname='temp'
       CASE(2)
         trname='salt'
       CASE DEFAULT
         write(trname,'(A3,i1)') 'ptr', j
     END SELECT
     status = nf_def_var(ncid, trim(trname), NF_DOUBLE, 3, dimid3, tra_varid(j))
     if (status .ne. nf_noerr) call handle_err(status)
     ! Adams–Bashforth part
     trname=trim(trname)//'_AB'
     status = nf_def_var(ncid, trim(trname), NF_DOUBLE, 3, dimid3, tra_varid_ab(j))
     if (status .ne. nf_noerr) call handle_err(status)
  end do

  ! Assign long_name and units attributes to variables.

  longname='model time'
  status = nf_PUT_ATT_TEXT(ncid, time_varid, 'long_name', len_trim(longname), trim(longname))
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_put_att_text(ncid, time_varid, 'units', 1, 's')
  if (status .ne. nf_noerr) call handle_err(status)
  longname='iteration_count'
  status = nf_PUT_ATT_TEXT(ncid, iter_varid, 'long_name', len_trim(longname), trim(longname))
  if (status .ne. nf_noerr) call handle_err(status)

  longname='sea surface elevation'
  status = nf_put_att_text(ncid, ssh_varid, 'description', len_trim(longname), trim(longname)) 
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_put_att_text(ncid, ssh_varid, 'units', 1, 'm')
  if (status .ne. nf_noerr) call handle_err(status)
  longname='zonal velocity'
  status = nf_put_att_text(ncid, u_varid, 'description', len_trim(longname), trim(longname)) 
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_put_att_text(ncid, u_varid, 'units', 3, 'm/s')
  if (status .ne. nf_noerr) call handle_err(status)
  longname='meridional velocity'
  status = nf_put_att_text(ncid, v_varid, 'description', len_trim(longname), trim(longname)) 
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_put_att_text(ncid, v_varid, 'units', 3, 'm/s')
  if (status .ne. nf_noerr) call handle_err(status)
  longname='Adams-Bashforth RHS for zonal velocity'
  status = nf_put_att_text(ncid, urhs_varid_AB, 'description', len_trim(longname), trim(longname)) 
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_put_att_text(ncid, urhs_varid_AB, 'units', 3, 'm3/s2')
  if (status .ne. nf_noerr) call handle_err(status)
  longname='Adams-Bashforth RHS for meridional velocity'
  status = nf_put_att_text(ncid, vrhs_varid_AB, 'description', len_trim(longname), trim(longname)) 
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_put_att_text(ncid, vrhs_varid_AB, 'units', 3, 'm3/s2')
  if (status .ne. nf_noerr) call handle_err(status)
  longname='vertical velocity'
  status = nf_put_att_text(ncid, w_varid, 'description', len_trim(longname), trim(longname)) 
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_put_att_text(ncid, w_varid, 'units', 3, 'm/s')
  if (status .ne. nf_noerr) call handle_err(status)
  longname='vertical velocity (explicit part)'
  status = nf_put_att_text(ncid, we_varid, 'description', len_trim(longname), trim(longname)) 
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_put_att_text(ncid, we_varid, 'units', 3, 'm/s')
  if (status .ne. nf_noerr) call handle_err(status)
  longname='vertical velocity (implicit part)'
  status = nf_put_att_text(ncid, wi_varid, 'description', len_trim(longname), trim(longname)) 
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_put_att_text(ncid, wi_varid, 'units', 3, 'm/s')
  if (status .ne. nf_noerr) call handle_err(status)

  do j=1,num_tracers
     SELECT CASE (j) 
       CASE(1)
         longname='potential temperature'
         units='degC'
       CASE(2)
         longname='salinity'
         units='psu'
       CASE DEFAULT
         write(longname,'(A15,i1)') 'passive tracer ', j
         units='none'
     END SELECT
     status = nf_put_att_text(ncid, tra_varid(j), 'description', len_trim(longname), trim(longname))
     if (status .ne. nf_noerr) call handle_err(status)
     status = nf_put_att_text(ncid, tra_varid(j), 'units', len_trim(units), trim(units))
     if (status .ne. nf_noerr) call handle_err(status)
     ! Adams-Bashforth part
     longname=trim(longname)//', Adams–Bashforth'
     status = nf_put_att_text(ncid, tra_varid_ab(j), 'description', len_trim(longname), trim(longname))
     if (status .ne. nf_noerr) call handle_err(status)
     status = nf_put_att_text(ncid, tra_varid_ab(j), 'units', len_trim(units), trim(units))
     if (status .ne. nf_noerr) call handle_err(status)
  end do

  status = nf_enddef(ncid)
  if (status .ne. nf_noerr) call handle_err(status)

  status=nf_close(ncid)
  if (status .ne. nf_noerr) call handle_err(status)

  ! ice part
  ! create an ice output file
  if (use_ice) then
     filename=trim(ResultPath)//runid//'.'//cyearnew//'.ice.restart.nc'
     status = nf_create(filename, nf_clobber, ncid)
     if (status.ne.nf_noerr) call handle_err(status)

     ! Define the dimensions
     status = nf_def_dim(ncid, 'nodes_2d', nod2d, dimid_2d)
     if (status .ne. nf_noerr) call handle_err(status)
     status = nf_def_dim(ncid, 'elem_2d', elem2d, dimid_2de)
     if (status .ne. nf_noerr) call handle_err(status)
     status = nf_def_dim(ncid, 'T', NF_UNLIMITED, dimid_rec)
     if (status .ne. nf_noerr) call handle_err(status)

     ! Define the time and iteration variables
     status = nf_def_var(ncid, 'time', NF_DOUBLE, 1, dimid_rec, time_varid)
     if (status .ne. nf_noerr) call handle_err(status)
     status = nf_def_var(ncid, 'iter', NF_INT, 1, dimid_rec, iter_varid)
     if (status .ne. nf_noerr) call handle_err(status)
     ! Define the netCDF variables for 2D fields.
     ! In Fortran, the unlimited dimension must come
     ! last on the list of dimids.
     dimids(1) = dimid_2d
     dimids(2) = dimid_rec

     status = nf_def_var(ncid, 'area', NF_DOUBLE, 2, dimids, area_varid)
     if (status .ne. nf_noerr) call handle_err(status)
     status = nf_def_var(ncid, 'hice', NF_DOUBLE, 2, dimids, hice_varid)
     if (status .ne. nf_noerr) call handle_err(status)
     status = nf_def_var(ncid, 'hsnow', NF_DOUBLE, 2, dimids, hsnow_varid)
     if (status .ne. nf_noerr) call handle_err(status)
  
     status = nf_def_var(ncid, 'uice', NF_DOUBLE, 2, dimids, uice_varid)
     if (status .ne. nf_noerr) call handle_err(status)
     status = nf_def_var(ncid, 'vice', NF_DOUBLE, 2, dimids, vice_varid)
     if (status .ne. nf_noerr) call handle_err(status)
  
     ! Assign long_name and units attributes to variables.
     longname='model time'
     status = nf_PUT_ATT_TEXT(ncid, time_varid, 'long_name', len_trim(longname), trim(longname)) 
     if (status .ne. nf_noerr) call handle_err(status)
     status = nf_put_att_text(ncid, time_varid, 'units', 1, 's')
     if (status .ne. nf_noerr) call handle_err(status)
     longname='iteration_count'
     status = nf_PUT_ATT_TEXT(ncid, iter_varid, 'long_name', len_trim(longname), trim(longname)) 
     if (status .ne. nf_noerr) call handle_err(status)
     longname='ice concentration [0 to 1]'
     status = nf_PUT_ATT_TEXT(ncid, area_varid, 'description', len_trim(longname), trim(longname)) 
     if (status .ne. nf_noerr) call handle_err(status)
     longname='effective ice thickness'
     status = nf_PUT_ATT_TEXT(ncid, hice_varid, 'description', len_trim(longname), trim(longname)) 
     if (status .ne. nf_noerr) call handle_err(status)
     status = nf_put_att_text(ncid, hice_varid, 'units', 1, 'm')
     if (status .ne. nf_noerr) call handle_err(status)
     longname='effective snow thickness'
     status = nf_PUT_ATT_TEXT(ncid, hsnow_varid, 'description', len_trim(longname), trim(longname)) 
     if (status .ne. nf_noerr) call handle_err(status)
     status = nf_put_att_text(ncid, hsnow_varid, 'units', 1, 'm')
     if (status .ne. nf_noerr) call handle_err(status)
     longname='zonal velocity'
     status = nf_PUT_ATT_TEXT(ncid, uice_varid, 'description', len_trim(longname), trim(longname)) 
     if (status .ne. nf_noerr) call handle_err(status)
     status = nf_put_att_text(ncid, uice_varid, 'units', 3, 'm/s')
     if (status .ne. nf_noerr) call handle_err(status)
     longname='meridional velocity'
     status = nf_PUT_ATT_TEXT(ncid, vice_varid, 'description', len_trim(longname), trim(longname)) 
     if (status .ne. nf_noerr) call handle_err(status)
     status = nf_put_att_text(ncid, vice_varid, 'units', 3, 'm/s')
     if (status .ne. nf_noerr) call handle_err(status)

     status = nf_enddef(ncid)
     if (status .ne. nf_noerr) call handle_err(status)

     status=nf_close(ncid)
     if (status .ne. nf_noerr) call handle_err(status)
  endif
end subroutine init_output_restart
!======================================================================
subroutine write_restarts(istep)

  use o_arrays
  use o_mesh
  use i_arrays
  use g_config
  use g_clock
  use g_parsup
  use g_comm_auto
  implicit none

#include "netcdf.inc" 

  integer                   :: status, ncid, j, istep
  integer                   :: time_varid, iter_varid
  integer                   :: ssh_varid, tra_varid(num_tracer), tra_varid_ab(num_tracer)
  integer                   :: u_varid, v_varid, wpot_varid
  integer                   :: w_varid, we_varid, wi_varid
  integer                   :: urhs_varid_AB, vrhs_varid_AB
  integer                   :: area_varid, hice_varid, hsnow_varid
  integer                   :: uice_varid, vice_varid, count_id
  integer                   :: start(2), count(2),start3(3), count3(3) 
  real(kind=8)              :: sec_in_year
  character(100)            :: filename
  character(100)            :: trname
  real(kind=8), allocatable :: aux2(:), aux3(:,:) 

  ! ocean part
  if (mype==0) then ! Serial output implemented so far
     ! open files
     filename=trim(ResultPath)//runid//'.'//cyearnew//'.oce.restart.nc'
     status = nf_open(filename, nf_write, ncid)
     if (status .ne. nf_noerr) call handle_err(status)

     ! inquire variable id
     status=nf_inq_varid(ncid, 'time', time_varid)
     if (status .ne. nf_noerr) call handle_err(status)
     status=nf_inq_varid(ncid, 'iter', iter_varid)
     if (status .ne. nf_noerr) call handle_err(status)
     status=nf_inq_varid(ncid, 'ssh', ssh_varid)
     if (status .ne. nf_noerr) call handle_err(status)
     status=nf_inq_varid(ncid, 'u', u_varid)
     if (status .ne. nf_noerr) call handle_err(status)
     status=nf_inq_varid(ncid, 'v', v_varid)
     if (status .ne. nf_noerr) call handle_err(status)
     status=nf_inq_varid(ncid, 'urhs_AB', urhs_varid_AB)
     if (status .ne. nf_noerr) call handle_err(status)
     status=nf_inq_varid(ncid, 'vrhs_AB', vrhs_varid_AB)
     if (status .ne. nf_noerr) call handle_err(status)
     status=nf_inq_varid(ncid, 'w', w_varid)
     if (status .ne. nf_noerr) call handle_err(status)
     status=nf_inq_varid(ncid, 'w_expl', we_varid)
     if (status .ne. nf_noerr) call handle_err(status)
     status=nf_inq_varid(ncid, 'w_impl', wi_varid)
     if (status .ne. nf_noerr) call handle_err(status)

      do j=1,num_tracers
        SELECT CASE (j) 
          CASE(1)
            trname='temp'
          CASE(2)
            trname='salt'
          CASE DEFAULT
            write(trname,'(A3,i1)') 'ptr', j
        END SELECT
        status=nf_inq_varid(ncid, trim(trname), tra_varid(j))
        if (status .ne. nf_noerr) call handle_err(status)
        ! Adams–Bashforth part
        trname=trim(trname)//'_AB'
        status=nf_inq_varid(ncid, trim(trname), tra_varid_ab(j))
        if (status .ne. nf_noerr) call handle_err(status)
     end do

     ! continue writing netcdf keeping the old records
     status = nf_inq_dimid(ncid, 'T', count_id)
     if(status .ne. nf_noerr) call handle_err(status)
     status=  nf_inq_dimlen(ncid, count_id, save_count_restart)
     if(status .ne. nf_noerr) call handle_err(status)
     save_count_restart=save_count_restart+1

     ! write variables
     sec_in_year=dt*istep
     ! time and iteration
     status=nf_put_vara_double(ncid, time_varid, save_count_restart, 1, sec_in_year)
     if (status .ne. nf_noerr) call handle_err(status)
     status=nf_put_vara_int(ncid, iter_varid, save_count_restart, 1, istep)
     if (status .ne. nf_noerr) call handle_err(status)
  end if ! mype==0
  allocate(aux2(nod2D), aux3(nl-1, elem2D)) ! auxuary arrays for broadcasting the SSH and horizontal velocities
  ! 2d fields
  call gather_nod(eta_n, aux2)
  if(mype==0) then            
     start=(/1,save_count_restart/)
     count=(/nod2d, 1/)
     status=nf_put_vara_double(ncid, ssh_varid, start, count, aux2, 4)
     if (status .ne. nf_noerr) call handle_err(status)
  end if

  ! 3d fields U, V
  call gather_elem(UV(1,:,:), aux3)
  if (mype==0) then                  
     start3=(/1, 1, save_count_restart/)
     count3=(/nl-1, elem2D, 1/)
     status=nf_put_vara_double(ncid, u_varid, start3, count3, aux3) 
     if (status .ne. nf_noerr) call handle_err(status)
  end if

  call gather_elem(UV(2,:,:),aux3)  
  if(mype==0) then                      
     status=nf_put_vara_double(ncid, v_varid, start3, count3, aux3)
     if (status .ne. nf_noerr) call handle_err(status)
  end if

  ! 3d fields UV_rhs_AB
  call gather_elem(UV_rhsAB(1,:,:), aux3)
  if (mype==0) then                  
     status=nf_put_vara_double(ncid, urhs_varid_AB, start3, count3, aux3) 
     if (status .ne. nf_noerr) call handle_err(status)
  end if

  call gather_elem(UV_rhsAB(2,:,:), aux3)  
  if(mype==0) then                      
     status=nf_put_vara_double(ncid, vrhs_varid_AB, start3, count3, aux3)
     if (status .ne. nf_noerr) call handle_err(status)
  end if


  deallocate(aux3) !reallocate for w
  allocate(aux3(nl,nod2D))

  call gather_nod(Wvel,aux3)
  if(mype==0) then                        
     count3=(/nl, nod2D, 1/)
     status=nf_put_vara_double(ncid, w_varid, start3, count3, aux3)
     if (status .ne. nf_noerr) call handle_err(status)
  end if
  
  call gather_nod(Wvel_e,aux3)
  if(mype==0) then                        
     count3=(/nl, nod2D, 1/)
     status=nf_put_vara_double(ncid, we_varid, start3, count3, aux3)
     if (status .ne. nf_noerr) call handle_err(status)
  end if

  call gather_nod(Wvel_i,aux3)
  if(mype==0) then                        
     count3=(/nl, nod2D, 1/)
     status=nf_put_vara_double(ncid, wi_varid, start3, count3, aux3)
     if (status .ne. nf_noerr) call handle_err(status)
  end if

  deallocate(aux3) !reallocate for scalar variables
  allocate(aux3(nl-1,nod2D))
  
  !T,S and passive tracers
  count3=(/nl-1, nod2D, 1/)
  do j=1,num_tracers
     call gather_nod(tr_arr(:,:,j),aux3)
     if(mype==0) then                        
        status=nf_put_vara_double(ncid, tra_varid(j), start3, count3, aux3)
        if (status .ne. nf_noerr) call handle_err(status)
     end if
     ! Adams–Bashforth part
     call gather_nod(tr_arr_old(:,:,j),aux3)
     if(mype==0) then                        
        status=nf_put_vara_double(ncid, tra_varid_ab(j), start3, count3, aux3)
        if (status .ne. nf_noerr) call handle_err(status)
     end if
  end do

  if(mype==0) then
     status=nf_close(ncid)
     if (status .ne. nf_noerr) call handle_err(status)
  end if
  deallocate(aux3)
  ! ice part
  if (use_ice) then
     if (mype==0) then ! Serial output implemented so far
        ! open files
        filename=trim(ResultPath)//runid//'.'//cyearnew//'.ice.restart.nc'
        status = nf_open(filename, nf_write, ncid)
        if (status .ne. nf_noerr) call handle_err(status)
        ! inquire variable id
        status=nf_inq_varid(ncid, 'time', time_varid)
        if (status .ne. nf_noerr) call handle_err(status)
        status=nf_inq_varid(ncid, 'iter', iter_varid)
        if (status .ne. nf_noerr) call handle_err(status)
        status=nf_inq_varid(ncid, 'area', area_varid)
        if (status .ne. nf_noerr) call handle_err(status)
        status=nf_inq_varid(ncid, 'hice', hice_varid)
        if (status .ne. nf_noerr) call handle_err(status)
        status=nf_inq_varid(ncid, 'hsnow', hsnow_varid)
        if (status .ne. nf_noerr) call handle_err(status)
        status=nf_inq_varid(ncid, 'uice', uice_varid)
        if (status .ne. nf_noerr) call handle_err(status)
        status=nf_inq_varid(ncid, 'vice', vice_varid)
        if (status .ne. nf_noerr) call handle_err(status)

        ! write variables
        ! time and iteration
        status=nf_put_vara_double(ncid, time_varid, save_count_restart, 1, sec_in_year)
        if (status .ne. nf_noerr) call handle_err(status)
        status=nf_put_vara_int(ncid, iter_varid, save_count_restart, 1, istep)
        if (status .ne. nf_noerr) call handle_err(status)

        ! 2d fields
        start=(/1,save_count_restart/)
        count=(/nod2d, 1/)
     end if ! mype==0

     call gather_nod(a_ice,aux2)                
     if (mype==0) then                           
        status=nf_put_vara_double(ncid, area_varid, start, count, aux2)
        if (status .ne. nf_noerr) call handle_err(status)
     end if
     call gather_nod(m_ice,aux2)
     if (mype==0) then                          
        status=nf_put_vara_double(ncid, hice_varid, start, count, aux2)
        if (status .ne. nf_noerr) call handle_err(status)
     end if
     call gather_nod(m_snow,aux2)
     if (mype==0) then                            
        status=nf_put_vara_double(ncid, hsnow_varid, start, count, aux2)
        if (status .ne. nf_noerr) call handle_err(status)
     end if
     call gather_nod(u_ice,aux2)
     if (mype==0) then
        status=nf_put_vara_double(ncid, uice_varid, start, count, aux2)
        if (status .ne. nf_noerr) call handle_err(status)
     end if
     call gather_nod(v_ice,aux2)
     if (mype==0) then                          
        status=nf_put_vara_double(ncid, vice_varid, start, count, aux2)
        if (status .ne. nf_noerr) call handle_err(status)
     end if
     if(mype==0) then
        status=nf_close(ncid)
        if (status .ne. nf_noerr) call handle_err(status)
     end if
  end if !use_ice
  deallocate(aux2)
end subroutine write_restarts
!
!--------------------------------------------------------------------------------------------
!
subroutine restart(directionflag, istep)

  use g_config
  use g_clock
  use g_parsup
  implicit none

  logical :: do_restart=.false.
  integer :: directionflag,istep
  !check whether we want to do output
  if (restart_length_unit.eq.'y') then
     call annual_event(do_restart)
  else if (restart_length_unit.eq.'m') then 
     call monthly_event(do_restart) 
  else if (restart_length_unit.eq.'d') then
     call daily_event(do_restart, restart_length)
  else if (restart_length_unit.eq.'h') then
     call hourly_event(do_restart, restart_length)
  else if (restart_length_unit.eq.'s') then
     call step_event(do_restart, istep, restart_length)
  else
     write(*,*) 'You did not specify a supported outputflag.'
     write(*,*) 'The program will stop to give you opportunity to do it.'
     call par_ex(1)
     stop
  endif

  if (directionflag.eq.1) do_restart=.true.  

  if (.not.do_restart) return

  ! write results
if (mype==0) write (*,*) 'Writing of restart files is disabled!'
!  if(mype==0) write(*,*)'Do output (netCDF, restart) ...'
!  call write_restarts(istep)
end subroutine restart
