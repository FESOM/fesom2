subroutine init_output_mean(do_init)
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
  integer                   :: ssh_varid, tra_varid(num_tracer)
  integer                   :: u_varid, v_varid, w_varid, wpot_varid, cflw_varid
  integer                   :: gmu_varid, gmv_varid, gmw_varid
  integer                   :: area_varid, hice_varid, hsnow_varid
  integer                   :: uice_varid, vice_varid
  character(100)            :: longname
  character(100)            :: filename
  character(100)            :: trname, units
  integer                   :: av_varid, kv_varid(num_tracer), hbl_varid

  if (.not. do_init) return
  ! Serial output implemented so far
  if (mype/=0) return

  ! create an ocean output file
  write(*,*) 'initialize new output files'
  filename=trim(ResultPath)//runid//'.'//cyearnew//'.oce.nc'

  if (output_offset==32) then
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
  
  status = nf_def_var(ncid, 'ssh', NF_FLOAT, 2, dimids, ssh_varid)
  if (status .ne. nf_noerr) call handle_err(status)

  if (hbl_diag) then
   status = nf_def_var(ncid, 'hbl', NF_FLOAT, 2, dimids, hbl_varid)
   if (status .ne. nf_noerr) call handle_err(status)
  endif

  ! Define the netCDF variables for 3D fields
  ! ocean horizontal velocity u, v
  dimid3(1) = dimid_nl1
  dimid3(2) = dimid_2de
  dimid3(3) = dimid_rec

  status = nf_def_var(ncid, 'u', NF_FLOAT, 3, dimid3, u_varid)
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_def_var(ncid, 'v', NF_FLOAT, 3, dimid3, v_varid)
  if (status .ne. nf_noerr) call handle_err(status)

  if (Fer_GM) then
     status = nf_def_var(ncid, 'GM_u', NF_FLOAT, 3, dimid3, gmu_varid)
     if (status .ne. nf_noerr) call handle_err(status)

     status = nf_def_var(ncid, 'GM_v', NF_FLOAT, 3, dimid3, gmv_varid)
     if (status .ne. nf_noerr) call handle_err(status)
  end if
  
  ! vertical viscosity coefficient
  if (AvKv) then 
     status = nf_def_var(ncid, 'Av', NF_FLOAT, 3, dimid3, av_varid)
     if (status .ne. nf_noerr) call handle_err(status)
  endif
     
  ! ocean vertical velocity w
  dimid3(1) = dimid_nl
  dimid3(2) = dimid_2d
  dimid3(3) = dimid_rec

  status = nf_def_var(ncid, 'w', NF_FLOAT, 3, dimid3, w_varid)
  if (status .ne. nf_noerr) call handle_err(status)

  ! vertical diffusivity coefficient
  if (AvKv) then 
    do j=1,num_tracers
       SELECT CASE (j) 
          CASE(1)
            trname='Kvt'
          CASE(2)
            trname='Kvs'
          CASE DEFAULT
            write(trname,'(A3,i1)') 'ptr', j
       END SELECT
       status = nf_def_var(ncid, trim(trname), NF_FLOAT, 3, dimid3, kv_varid(j))
       if (status .ne. nf_noerr) call handle_err(status)
    end do
  endif

  if (Fer_GM) then
     status = nf_def_var(ncid, 'GM_w', NF_FLOAT, 3, dimid3, gmw_varid)
     if (status .ne. nf_noerr) call handle_err(status)
  end if

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
     status = nf_def_var(ncid, trim(trname), NF_FLOAT, 3, dimid3, tra_varid(j))
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
  if(hbl_diag) then
     longname='boundary layer depth'
     status = nf_put_att_text(ncid, hbl_varid, 'description', len_trim(longname), trim(longname)) 
     if (status .ne. nf_noerr) call handle_err(status)
     status = nf_put_att_text(ncid, hbl_varid, 'units', 1, 'm')
     if (status .ne. nf_noerr) call handle_err(status)
  endif
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
  longname='vertical velocity'
  status = nf_put_att_text(ncid, w_varid, 'description', len_trim(longname), trim(longname)) 
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_put_att_text(ncid, w_varid, 'units', 3, 'm/s')
  if (status .ne. nf_noerr) call handle_err(status)
  if(AvKv) then
     longname='vertical eddy viscosity coefficient'
     status = nf_put_att_text(ncid, av_varid, 'description', len_trim(longname), trim(longname)) 
     if (status .ne. nf_noerr) call handle_err(status)
     status = nf_put_att_text(ncid, av_varid, 'units', 4, 'm2/s')
     if (status .ne. nf_noerr) call handle_err(status)
     do j=1,num_tracers
        SELECT CASE (j) 
          CASE(1)
            longname='vertical eddy diffusion coefficient (temperature)'
            units='m2/s'
          CASE(2)
            longname='vertical eddy diffusion coefficient (salinity)'
            units='m2/s'
          CASE DEFAULT
            write(longname,'(A15,i1)') 'passive tracer ', j
            units='none'
        END SELECT
        status = nf_put_att_text(ncid, kv_varid(j), 'description', len_trim(longname), trim(longname))
        if (status .ne. nf_noerr) call handle_err(status)
        status = nf_put_att_text(ncid, kv_varid(j), 'units', len_trim(units), trim(units))
        if (status .ne. nf_noerr) call handle_err(status)
     end do
  end if 
  if (Fer_GM) then
     longname='subgrid zonal velocity'
     status = nf_put_att_text(ncid, gmu_varid, 'description', len_trim(longname), trim(longname)) 
     if (status .ne. nf_noerr) call handle_err(status)
     status = nf_put_att_text(ncid, gmu_varid, 'units', 3, 'm/s')
     if (status .ne. nf_noerr) call handle_err(status)
     longname='subgrid meridional velocity'
     status = nf_put_att_text(ncid, gmv_varid, 'description', len_trim(longname), trim(longname)) 
     if (status .ne. nf_noerr) call handle_err(status)
     status = nf_put_att_text(ncid, gmv_varid, 'units', 3, 'm/s')
     if (status .ne. nf_noerr) call handle_err(status)
     longname='subgrid vertical velocity'
     status = nf_put_att_text(ncid, gmw_varid, 'description', len_trim(longname), trim(longname)) 
     if (status .ne. nf_noerr) call handle_err(status)
     status = nf_put_att_text(ncid, gmw_varid, 'units', 3, 'm/s')
     if (status .ne. nf_noerr) call handle_err(status)
  end if

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
  end do


  status = nf_enddef(ncid)
  if (status .ne. nf_noerr) call handle_err(status)

  status=nf_close(ncid)
  if (status .ne. nf_noerr) call handle_err(status)

  ! ice part
  ! create an ice output file
  if (use_ice) then
     filename=trim(ResultPath)//runid//'.'//cyearnew//'.ice.nc'
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

     status = nf_def_var(ncid, 'area', NF_FLOAT, 2, dimids, area_varid)
     if (status .ne. nf_noerr) call handle_err(status)
     status = nf_def_var(ncid, 'hice', NF_FLOAT, 2, dimids, hice_varid)
     if (status .ne. nf_noerr) call handle_err(status)
     status = nf_def_var(ncid, 'hsnow', NF_FLOAT, 2, dimids, hsnow_varid)
     if (status .ne. nf_noerr) call handle_err(status)
  
     status = nf_def_var(ncid, 'uice', NF_FLOAT, 2, dimids, uice_varid)
     if (status .ne. nf_noerr) call handle_err(status)
     status = nf_def_var(ncid, 'vice', NF_FLOAT, 2, dimids, vice_varid)
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
end subroutine init_output_mean
!======================================================================
subroutine write_means(istep)

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
  integer                   :: ssh_varid, tra_varid(num_tracer)
  integer                   :: u_varid, v_varid, w_varid, wpot_varid
  integer                   :: gmu_varid, gmv_varid, gmw_varid
  integer                   :: area_varid, hice_varid, hsnow_varid
  integer                   :: uice_varid, vice_varid, count_id
  integer                   :: start(2), count(2),start3(3), count3(3) 
  real(kind=8)              :: sec_in_year
  character(100)            :: filename
  character(100)            :: trname
  real(kind=4), allocatable :: aux2(:), aux3(:,:) 
  integer                   :: av_varid, kv_varid(num_tracer), hbl_varid

  ! ocean part
  if (mype==0) then ! Serial output implemented so far
     ! open files
     filename=trim(ResultPath)//runid//'.'//cyearnew//'.oce.nc'
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
     status=nf_inq_varid(ncid, 'w', w_varid)
     if (status .ne. nf_noerr) call handle_err(status)
     if (hbl_diag) then
        status=nf_inq_varid(ncid, 'hbl', hbl_varid)
        if (status .ne. nf_noerr) call handle_err(status)
     endif     
     if(AvKv) then
        status=nf_inq_varid(ncid, 'Av', av_varid)
        if (status .ne. nf_noerr) call handle_err(status)
        do j=1,num_tracers
           SELECT CASE (j) 
             CASE(1)
               trname='Kvt'
             CASE(2)
               trname='Kvs'
             CASE DEFAULT
               write(trname,'(A3,i1)') 'ptr', j
           END SELECT
           status = nf_inq_varid(ncid, trim(trname), kv_varid(j))
           if (status .ne. nf_noerr) call handle_err(status)
        end do
     endif

     if (Fer_GM) then
        status=nf_inq_varid(ncid, 'GM_u', gmu_varid)
        if (status .ne. nf_noerr) call handle_err(status)
        status=nf_inq_varid(ncid, 'GM_v', gmv_varid)
        if (status .ne. nf_noerr) call handle_err(status)
        status=nf_inq_varid(ncid, 'GM_w', gmw_varid)
        if (status .ne. nf_noerr) call handle_err(status)
     end if	

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
     end do
     ! continue writing netcdf keeping the old records
     status = nf_inq_dimid(ncid, 'T', count_id)
     if(status .ne. nf_noerr) call handle_err(status)
     status=  nf_inq_dimlen(ncid, count_id, save_count_mean)
     if(status .ne. nf_noerr) call handle_err(status)
     save_count_mean=save_count_mean+1
     ! write variables
     sec_in_year=dt*istep
     ! time and iteration
     status=nf_put_vara_double(ncid, time_varid, save_count_mean, 1, sec_in_year)
     if (status .ne. nf_noerr) call handle_err(status)
     status=nf_put_vara_int(ncid, iter_varid, save_count_mean, 1, istep)
     if (status .ne. nf_noerr) call handle_err(status)
  end if ! mype==0
  allocate(aux2(nod2D), aux3(nl-1, elem2D)) ! auxuary arrays for broadcasting the SSH and horizontal velocities
  ! 2d fields
  call gather_nod(eta_n_mean, aux2)
  if(mype==0) then            
     start=(/1,save_count_mean/)
     count=(/nod2d, 1/)
     status=nf_put_vara_real(ncid, ssh_varid, start, count, aux2) 
     if (status .ne. nf_noerr) call handle_err(status)
  end if
  if (hbl_diag) then
     call gather_nod(hbl_mean, aux2)
     if(mype==0) then            
        start=(/1,save_count_mean/)
        count=(/nod2d, 1/)
        status=nf_put_vara_real(ncid, hbl_varid, start, count, aux2) 
        if (status .ne. nf_noerr) call handle_err(status)
     end if
  endif

  ! 3d fields
  call gather_elem(UV_mean(1,:,:), aux3)
  if (mype==0) then                  
     start3=(/1, 1, save_count_mean/)
     count3=(/nl-1, elem2D, 1/)
     status=nf_put_vara_real(ncid, u_varid, start3, count3, aux3) 
     if (status .ne. nf_noerr) call handle_err(status)
  end if

  call gather_elem(UV_mean(2,:,:),aux3)  
  if(mype==0) then                      
     status=nf_put_vara_real(ncid, v_varid, start3, count3, aux3)
     if (status .ne. nf_noerr) call handle_err(status)
  end if
  
  if (Fer_GM) then
     call gather_elem(fer_UV_mean(1,:,:), aux3)
     if (mype==0) then                  
        start3=(/1, 1, save_count_mean/)
        count3=(/nl-1, elem2D, 1/)
        status=nf_put_vara_real(ncid, gmu_varid, start3, count3, aux3) 
        if (status .ne. nf_noerr) call handle_err(status)
     end if

     call gather_elem(fer_UV_mean(2,:,:),aux3)  
     if(mype==0) then                      
        status=nf_put_vara_real(ncid, gmv_varid, start3, count3, aux3)
        if (status .ne. nf_noerr) call handle_err(status)
     end if     
   end if

  if(AvKv) then
     call gather_elem(Av_mean,aux3)
     if(mype==0) then                        
        start3=(/1, 1, save_count_mean/)
        count3=(/nl-1, elem2D, 1/)
        status=nf_put_vara_real(ncid, av_varid, start3, count3, aux3)
        if (status .ne. nf_noerr) call handle_err(status)
     end if
  endif

  deallocate(aux3) !reallocate for w
  allocate(aux3(nl,nod2D))

  call gather_nod(Wvel_mean,aux3)
  if(mype==0) then                        
     count3=(/nl, nod2D, 1/)
     status=nf_put_vara_real(ncid, w_varid, start3, count3, aux3)
     if (status .ne. nf_noerr) call handle_err(status)
  end if

! _OG_
  if(AvKv) then
     count3=(/nl, nod2D, 1/)
     do j=1,num_tracers
        call gather_nod(Kv_mean(:,:,j),aux3)
        if(mype==0) then                        
           status=nf_put_vara_real(ncid, kv_varid(j), start3, count3, aux3)
           if (status .ne. nf_noerr) call handle_err(status)
        end if 
     end do
  end if

  if (Fer_GM) then
     call gather_nod(fer_Wvel_mean,aux3)
     if (mype==0) then                        
        count3=(/nl, nod2D, 1/)
        status=nf_put_vara_real(ncid, gmw_varid, start3, count3, aux3)
        if (status .ne. nf_noerr) call handle_err(status)
     end if     
  end if

  deallocate(aux3) !reallocate for scalar variables
  allocate(aux3(nl-1,nod2D))

  !T,S and passive tracers
  count3=(/nl-1, nod2D, 1/)
  do j=1,num_tracers
     call gather_nod(tr_arr_mean(:,:,j),aux3)
     if(mype==0) then                        
        status=nf_put_vara_real(ncid, tra_varid(j), start3, count3, aux3)
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
        filename=trim(ResultPath)//runid//'.'//cyearnew//'.ice.nc'
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
        status=nf_put_vara_double(ncid, time_varid, save_count_mean, 1, sec_in_year)
        if (status .ne. nf_noerr) call handle_err(status)
        status=nf_put_vara_int(ncid, iter_varid, save_count_mean, 1, istep)
        if (status .ne. nf_noerr) call handle_err(status)

        ! 2d fields
        start=(/1,save_count_mean/)
        count=(/nod2d, 1/)
     end if ! mype==0

     call gather_nod(a_ice_mean,aux2)                
     if (mype==0) then                           
        status=nf_put_vara_real(ncid, area_varid, start, count, aux2)
        if (status .ne. nf_noerr) call handle_err(status)
     end if
     call gather_nod(m_ice_mean,aux2)              
     if (mype==0) then                          
        status=nf_put_vara_real(ncid, hice_varid, start, count, aux2)   
        if (status .ne. nf_noerr) call handle_err(status)
     end if
     call gather_nod(m_snow_mean,aux2)              
     if (mype==0) then                            
        status=nf_put_vara_real(ncid, hsnow_varid, start, count, aux2)  
        if (status .ne. nf_noerr) call handle_err(status)
     end if
     call gather_nod(u_ice_mean,aux2)              
     if (mype==0) then
        status=nf_put_vara_real(ncid, uice_varid, start, count, aux2)
        if (status .ne. nf_noerr) call handle_err(status)
     end if
     call gather_nod(v_ice_mean,aux2)              
     if (mype==0) then                          
        status=nf_put_vara_real(ncid, vice_varid, start, count, aux2)
        if (status .ne. nf_noerr) call handle_err(status)
     end if
     if(mype==0) then
        status=nf_close(ncid)
        if (status .ne. nf_noerr) call handle_err(status)
     end if
  end if !use_ice
  deallocate(aux2)
end subroutine write_means
!
!--------------------------------------------------------------------------------------------
!
subroutine update_means(do_output)
use g_parsup
use g_config
use o_arrays
use i_arrays
use g_clock
use o_mixing_KPP_mod !_OG_
implicit none
logical, intent(IN)       :: do_output
integer                   :: nsteps
integer                   :: j
ndpyr=365+fleapyear
  if (output_length_unit.eq.'y') then                                                                 
     nsteps=ndpyr*step_per_day                                                                    
  else if (output_length_unit.eq.'m') then                                                            
     nsteps=num_day_in_month(fleapyear,month)*step_per_day                                                                   
  else if (output_length_unit.eq.'d') then                                                            
     nsteps=step_per_day                                                                     
  else if (output_length_unit.eq.'h') then                                                            
     nsteps=int(real(step_per_day)/24.0)
  else if (output_length_unit.eq.'s') then
     nsteps=output_length                                                                      
  endif
   UV_mean=UV_mean+UV ; if (do_output) UV_mean=UV_mean/dble(nsteps)                              
   Wvel_mean=Wvel_mean+Wvel; if (do_output) Wvel_mean=Wvel_mean/dble(nsteps)
   eta_n_mean=eta_n_mean+eta_n; if (do_output) eta_n_mean=eta_n_mean/dble(nsteps)
   if(hbl_diag) then
      hbl_mean=hbl_mean+hbl; if (do_output) hbl_mean=hbl_mean/dble(nsteps)
   endif
   do j=1, num_tracers
      tr_arr_mean(:,:,j)=tr_arr_mean(:,:,j)+tr_arr(:,:,j)
   end do
   if (do_output) tr_arr_mean=tr_arr_mean/dble(nsteps)
   if (use_ice) then! ice has different update rate, so this is extra work
   U_ice_mean=U_ice_mean+U_ice ; if (do_output) U_ice_mean=U_ice_mean/dble(nsteps)
   V_ice_mean=V_ice_mean+V_ice ; if (do_output) V_ice_mean=V_ice_mean/dble(nsteps)
   m_ice_mean=m_ice_mean+m_ice ; if (do_output) m_ice_mean=m_ice_mean/dble(nsteps)
   a_ice_mean=a_ice_mean+a_ice ; if (do_output) a_ice_mean=a_ice_mean/dble(nsteps)
   m_snow_mean=m_snow_mean+m_snow ; if (do_output) m_snow_mean=m_snow_mean/dble(nsteps)
   endif
   if (AvKv) then
      Av_mean=Av_mean+Av; if (do_output) Av_mean=Av_mean/dble(nsteps)
      do j=1, num_tracers
         Kv_mean(:,:,j)=Kv_mean(:,:,j)+Kv2(:,:,j)
      end do
      if (do_output) Kv_mean=Kv_mean/dble(nsteps)
   endif
   if (Fer_GM) then
      fer_UV_mean=fer_UV_mean+fer_UV ; if (do_output) fer_UV_mean=fer_UV_mean/dble(nsteps)
      fer_wvel_mean=fer_wvel_mean+fer_wvel ; if (do_output) fer_wvel_mean=fer_wvel_mean/dble(nsteps)
   end if
end subroutine update_means
!
!--------------------------------------------------------------------------------------------
!
subroutine clean_means
use g_parsup
use g_config
use o_arrays
use i_arrays
implicit none
logical do_output

   UV_mean=0d0
   Wvel_mean=0d0
   eta_n_mean=0d0
   if (hbl_diag) then
      hbl_mean=0d0
   endif
   tr_arr_mean(:,:,:)=0d0
   if (use_ice) then! ice has different update rate, so this is extra work
      U_ice_mean=0d0
      V_ice_mean=0d0
      m_ice_mean=0d0
      a_ice_mean=0d0
      m_snow_mean=0d0
   endif
   if(AvKv) then
      Av_mean=0d0
      Kv_mean=0d0
   endif
   if (Fer_GM) then
      fer_UV_mean=0.0_WP
      fer_wvel_mean=0.0_WP
   end if
end subroutine clean_means
!
!--------------------------------------------------------------------------------------------
!
subroutine output(directionflag,istep)
  ! main output routine
  !
  ! Coded by Ralph Timmermann
  ! Modified by Qiang Wang for more diagnose output
  ! Reviewed by ??
  !--------------------------------------------------------------	

  use g_config
  use g_clock
  use g_parsup
  use o_arrays
  implicit none

  logical :: do_output=.false.
  integer :: directionflag,istep
  !check whether we want to do output
  if (output_length_unit.eq.'y') then
     call annual_event(do_output)
  else if (output_length_unit.eq.'m') then 
     call monthly_event(do_output) 
  else if (output_length_unit.eq.'d') then
     call daily_event(do_output, output_length)  
  else if (output_length_unit.eq.'h') then
     call hourly_event(do_output, output_length) 
  else if (output_length_unit.eq.'s') then
     call step_event(do_output, istep, output_length) 
  else
     write(*,*) 'You did not specify a supported outputflag.'
     write(*,*) 'The program will stop to give you opportunity to do it.'
     call par_ex(1)
     stop
  endif

  if (directionflag.eq.1) do_output=.true.  

  if (use_means) call update_means(do_output)

  if (.not.do_output) return

  ! write results

  if (use_means) then
     if(mype==0) write(*,*)'Do output (netCDF, mean) ...'
     call write_means(istep)
     call clean_means
  endif
end subroutine output
!
!--------------------------------------------------------------------------------------------
!
subroutine annual_event(do_output)
  !decides whether it's time to do output
  use g_clock
  implicit none

  logical :: do_output

  if ((daynew == ndpyr) .and. (timenew==86400.)) then
     do_output=.true.
  else
     do_output=.false.
  endif

end subroutine annual_event
!
!--------------------------------------------------------------------------------------------
!
subroutine monthly_event(do_output)
  !decides whether it's time to do output
  use g_clock
  implicit none

  logical :: do_output

  if (day_in_month==num_day_in_month(fleapyear,month) .and. &
       timenew==86400.) then
     do_output=.true.
  else
     do_output=.false.
  end if

end subroutine monthly_event
!
!--------------------------------------------------------------------------------------------
!
subroutine daily_event(do_output, N)
  !decides whether it's time to do output
  use g_clock
  implicit none

  logical             :: do_output
  integer, intent(in) :: N
  if (mod(daynew, N)==0 .and. timenew==86400.) then
     do_output=.true.
  else
     do_output=.false.
  endif

end subroutine daily_event
!
!--------------------------------------------------------------------------------------------
!
subroutine hourly_event(do_output, N)
  !decides whether it's time to do output
  use g_clock
  implicit none

  logical             :: do_output
  integer, intent(in) :: N

  if (mod(timenew, 3600.*N)==0) then
     do_output=.true.
  else
     do_output=.false.
  endif

end subroutine hourly_event
!
!--------------------------------------------------------------------------------------------
!
subroutine step_event(do_output,istep, N)
  !decides whether it's time to do output
  use g_config
  implicit none

  logical             :: do_output
  integer             :: istep
  integer, intent(in) :: N

  if (mod(istep, N)==0) then
     do_output=.true.
  else
     do_output=.false.
  endif

end subroutine step_event
!
!--------------------------------------------------------------------------------------------
!
subroutine handle_err(errcode)
  use g_parsup
  implicit none
  
#include "netcdf.inc" 
  
  integer errcode
  
  write(*,*) 'Error: ', nf_strerror(errcode)
  call par_ex(1)
  stop
end subroutine handle_err
!
!--------------------------------------------------------------------------------------------
!
