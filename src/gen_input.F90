! Routines for reading inital or re-start files
! Adapted from FESOM, changes involve: (i) vertical-horizontal storage
! (ii) horizontal velocities on elements
! (iii) T, S at midlevels, w at full levels
! SD, 18.04.2012
module g_input
contains
!read_prepared_initial_ice
!read_init_ts
!ice_input_unformatted
!oce_input_unformatted
!
subroutine read_prepared_initial_ice
  ! 
  ! Coded by Qiang Wang
  ! Reviewed by ??
  !------------------------------------------------------------------
  
  use o_mesh
  use i_arrays
  use g_config
  use g_clock
  use g_PARSUP
  implicit none

#include "netcdf.inc" 

  integer                   :: status, ncid, dimid_rec, nrec
  integer                   :: area_varid, hice_varid, hsnow_varid
  integer                   :: uice_varid, vice_varid
  integer                   :: istart(2), icount(2)
  character(100)            :: filename
  real(kind=WP), allocatable :: aux2(:)

  allocate(aux2(nod2D))  

  ! open files
  filename=trim(ResultPath)//trim(runid)//'.'//'initial_ice.nc'
  status = nf_open(filename, nf_nowrite, ncid)
  if (status .ne. nf_noerr) then
     print*,'ERROR: CANNOT READ initial ice FILE CORRECTLY !'
     print*,'Error in opening netcdf file'//filename
     call par_ex(1) 
     stop
  endif
  
  ! inquire variable id
  status=nf_inq_varid(ncid, 'area', area_varid)
  if (status .ne. nf_noerr) call handle_err(status)
  status=nf_inq_varid(ncid, 'hice', hice_varid)
  if (status .ne. nf_noerr) call handle_err(status)
  status=nf_inq_varid(ncid, 'hsnow', hsnow_varid)
  if (status .ne. nf_noerr) call handle_err(status)
  
  ! read variables

  ! which record to read
  if(restartflag=='last') then
     status = nf_inq_dimid(ncid, 'T', dimid_rec)
     if(status .ne. nf_noerr) call handle_err(status)
     status = nf_inq_dimlen(ncid, dimid_rec, nrec)
     if(status .ne. nf_noerr) call handle_err(status)
  else
     read(restartflag,'(i4)') nrec
  end if

  istart=(/1,nrec/)
  icount=(/nod2d, 1/)
  status=nf_get_vara_double(ncid, area_varid, istart, icount, aux2) 
  if (status .ne. nf_noerr) call handle_err(status)
  a_ice=aux2(myList_nod2D)     
  status=nf_get_vara_double(ncid, hice_varid, istart, icount, aux2) 
  if (status .ne. nf_noerr) call handle_err(status)
  m_ice=aux2(myList_nod2D)      
  status=nf_get_vara_double(ncid, hsnow_varid, istart, icount, aux2) 
  if (status .ne. nf_noerr) call handle_err(status)
  m_snow=aux2(myList_nod2D)      
  status=nf_close(ncid)
  if (status .ne. nf_noerr) call handle_err(status)

  deallocate(aux2)   

  !cutoff
  
  where(a_ice>1.0)
     a_ice=1.0
  end where

  where(a_ice<0.0)
     a_ice=0.0
  end where
  
  where(m_ice<0.0)
     m_ice=0.0
  end where

  where(m_snow<0.0)
     m_snow=0.0
  end where
  
end subroutine read_prepared_initial_ice
!
!-----------------------------------------------------------------------------
!
subroutine read_init_ts
  ! 
  ! Coded by Qiang Wang
  ! Reviewed by ??
  !------------------------------------------------------------------

  use o_mesh
  use o_param
  use o_arrays
  use g_config
  use g_PARSUP
  implicit none
  !
  integer                     :: i, j, n, nz, fileID
  integer                     :: num_lat_reg, num_lon_reg, num_lay_reg
  real(kind=WP)               :: pp, pr, tt, ss, lon, lat, tbott, sbott
  real(kind=WP), external     :: ptheta
  real(kind=WP), allocatable  :: lon_reg(:), lat_reg(:), lay_reg(:)
  real(kind=WP), allocatable  :: raw_data(:,:,:)
  real(kind=WP), allocatable  :: temp_x(:), temp_y(:)
  character*1000              :: filename
  real(kind=WP)               :: t0, t1
  integer                     :: ierror              ! return error code

  t0=MPI_Wtime()
  fileID=10
  filename=trim(ClimateDataPath)//trim(OceClimaDataName)
  
  ! 0 proc reads the data and sends it to other procs
  if (mype==0) then
     ! open global T/S data files
     write(*,*) 'reading ', trim(ClimateDataPath)//trim(OceClimaDataName)
     open(fileID,file=trim(filename), status='old')
     ! read reg. grid
     read(fileID, *) num_lon_reg, num_lat_reg, num_lay_reg
  end if
  call MPI_BCast(num_lon_reg,    1, MPI_INTEGER, 0, MPI_COMM_FESOM, ierror)
  call MPI_BCast(num_lat_reg,    1, MPI_INTEGER, 0, MPI_COMM_FESOM, ierror)
  call MPI_BCast(num_lay_reg,    1, MPI_INTEGER, 0, MPI_COMM_FESOM, ierror)

  allocate(lon_reg(num_lon_reg))
  allocate(lat_reg(num_lat_reg))
  allocate(lay_reg(num_lay_reg))

  ! 0 proc reads the data and sends it to other procs
  if (mype==0) then
     read(fileID, *) lon_reg
     read(fileID, *) lat_reg
     read(fileID, *) lay_reg
  end if
  call MPI_BCast(lon_reg, num_lon_reg, MPI_DOUBLE_PRECISION, 0, MPI_COMM_FESOM, ierror)
  call MPI_BCast(lat_reg, num_lat_reg, MPI_DOUBLE_PRECISION, 0, MPI_COMM_FESOM, ierror)
  call MPI_BCast(lay_reg, num_lay_reg, MPI_DOUBLE_PRECISION, 0, MPI_COMM_FESOM, ierror)

  allocate(raw_data(num_lon_reg,num_lat_reg,num_lay_reg))


  ! model grid coordinates
  allocate(temp_x(myDim_nod2d+eDim_nod2D), temp_y(myDim_nod2d+eDim_nod2D)) 
  do n=1, myDim_nod2d+eDim_nod2D                    
     temp_x(n)=geo_coord_nod2D(1,n)/rad
     temp_y(n)=geo_coord_nod2D(2,n)/rad
     ! change lon range to [0 360]
     if(temp_x(n)<0.) temp_x(n)=temp_x(n) + 360.0_WP  
  end do
 !==============================
 ! temperature
 !==============================
 ! 0 proc reads the data and sends it to other procs
  if (mype==0) then
     ! read raw data and do interpolation
     do i=1, num_lay_reg
        do j=1, num_lat_reg
           read(fileID, *) raw_data(:, j, i)
        end do
     end do  
  end if
  call MPI_BCast(raw_data(1,1,1), num_lon_reg*num_lat_reg*num_lay_reg, MPI_DOUBLE_PRECISION, 0, MPI_COMM_FESOM, ierror)

  call interp_3d_field(num_lon_reg, num_lat_reg, num_lay_reg, &
       lon_reg, lat_reg, lay_reg, raw_data, nl-1, myDim_nod2D+eDim_nod2D, &
       temp_x, temp_y, Z, tr_arr(:,:,1))

 !==============================
 ! salinity
 !==============================
 ! 0 proc reads the data and sends it to other procs
  if (mype==0) then
     ! read raw data and do interpolation
     do i=1, num_lay_reg
        do j=1, num_lat_reg
           read(fileID, *) raw_data(:, j, i)
        end do
     end do  
  end if
  call MPI_BCast(raw_data(1,1,1), num_lon_reg*num_lat_reg*num_lay_reg, MPI_DOUBLE_PRECISION, 0, MPI_COMM_FESOM, ierror)

  call interp_3d_field(num_lon_reg, num_lat_reg, num_lay_reg, &
       lon_reg, lat_reg, lay_reg,raw_data, nl-1, myDim_nod2d+eDim_nod2D, &
       temp_x, temp_y, Z, tr_arr(:,:,2))

  close(fileID) 

  ! This is recommended cutoff for the cold start initialization. Uncomment if the model explodes !!!
  ! where (tr_arr(:,:,2)<=28._WP)
  !       tr_arr(:,:,2)=28._WP
  ! end where

  ! Convert in situ temperature into potential temperature
  pr=0.0_WP
  do n=1,myDim_nod2d+eDim_nod2D
     do nz=1, nlevels_nod2D(n)-1    
        tt=tr_arr(nz,n,1)
        ss=tr_arr(nz,n,2)
        pp=abs(Z(nz))
        tr_arr(nz,n,1)=ptheta(ss, tt, pp, pr)
     end do	
  end do
  deallocate(temp_y, temp_x, raw_data, lay_reg, lat_reg, lon_reg)

  t1=MPI_Wtime()
  if (mype==0) then
     write(*,*) 'read 3D climatology completed in ', t1-t0, ' seconds'
     write(*,*) '========================='
  endif
end subroutine read_init_ts
end module g_input
