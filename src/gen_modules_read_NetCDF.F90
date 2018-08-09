! A collection of modules for reading different NetCDF files
! coded by Qiang Wang for FESOM1.4
! modified by Dmitry Sidorenko for reading with 0 PE and broadcasting to others  
module g_read_CORE_NetCDF
  ! read NetCDF data on T62 NCEP/NCAR grid

contains 
subroutine read_CORE_NetCDF(file,vari,itime,ncdata)

  use, intrinsic :: ISO_FORTRAN_ENV
  use g_config
  use g_parsup
  implicit none

#include "netcdf.inc" 

  integer, parameter                  :: nci=192, ncj=94  ! T62 grid dimension
  integer, dimension(3)               :: istart, icount
  integer			      :: ncid, varid, io
  integer                             :: i, j, itime
  real(real64), dimension (nci,ncj)   :: ncdata
  character(15)                       :: vari
  character(300)                      :: file
  integer                             :: ierror           ! return error code

  if (mype==0) then
       istart = (/1,1,itime/)
       icount= (/nci,ncj,1/)
! open netcdf input file
       io=nf_open(file, nf_nowrite, ncid)
  end if

  call MPI_BCast(io, 1, MPI_INTEGER, 0, MPI_COMM_FESOM, ierror)
  if (io.ne.nf_noerr)then
     print*,'ERROR: CANNOT READ FORCING FILE CORRECTLY !!!!!'
     print*,'Error in opening netcdf file'//file
     call par_ex
     stop 'Fatal error in open_netcdf'
  endif

  if (mype==0) then
! get handle for data variable    
       io=nf_inq_varid(ncid, vari, varid)
       io=nf_get_vara_double(ncid,varid,istart,icount,ncdata)
       io=nf_close(ncid)
  end if
  call MPI_BCast(ncdata, nci*ncj, MPI_DOUBLE_PRECISION, 0, MPI_COMM_FESOM, ierror)
  return
end subroutine read_CORE_NetCDF

end module g_read_CORE_NetCDF
!
!------------------------------------------------------------------------------------
!
module g_read_NCEP_NetCDF
  ! read NetCDF data on T62 NCEP/NCAR grid
contains 

subroutine read_NCEP_NetCDF(file,vari,itime,ncdata)

  use, intrinsic :: ISO_FORTRAN_ENV

  use g_config
  use g_parsup
  implicit none

#include "netcdf.inc" 

  integer, parameter                  :: nci=192, ncj=94  ! T62 grid dimension
  integer			      :: ncid, varid, io
  integer                             :: i, j, itime
  integer, dimension(3)               :: istart, icount
  integer(int16), dimension(nci,ncj) :: iuw 
  real(real32)                        :: xscale, xoff, miss
  real(real64), dimension(nci,ncj)    :: ncdata
  character(15)                       :: vari
  character(300)                      :: file
  integer                             :: ierror           ! return error code

  if (mype==0) then

     istart = (/1,1,itime/)
     icount= (/nci,ncj,1/)

     ! open netcdf input file
     io=nf_open(file, nf_nowrite, ncid)
  end if

  call MPI_BCast(io, 1, MPI_INTEGER, 0, MPI_COMM_FESOM, ierror)
  if (io.ne.nf_noerr)then
     print*,'ERROR: CANNOT READ FORCING FILE CORRECTLY !!!!!'
     print*,'Error in opening netcdf file'//file
     stop 'Fatal error in open_netcdf'
  endif

  if (mype==0) then
     ! get handle for data variable    
     io=nf_inq_varid(ncid, vari, varid)

     ! get att
     io=nf_get_att_real(ncid,varid,'scale_factor',xscale)
     io= nf_get_att_real(ncid,varid,'add_offset',xoff)
     ! io= nf_get_att_int(ncid,varid,'missing_value',miss)

     ! get variable
     io=nf_get_vara_int2(ncid,varid,istart,icount,iuw)

     ! close file 
     io=nf_close(ncid)
  end if
  call MPI_BCast(iuw, nci*ncj, MPI_DOUBLE_PRECISION, 0, MPI_COMM_FESOM, ierror)

  ! ncdata
  do j=1,ncj
     do i=1,nci
        ncdata(i,j)=real(iuw(i,j))*xscale+xoff
     enddo
  end do
  return
end subroutine read_NCEP_NetCDF
!
!--------------------------------------------------------------
!
subroutine upside_down(array,ni,nj)

  use, intrinsic :: ISO_FORTRAN_ENV

  implicit none
  integer                         :: i, j, ji
  integer                         :: ni, nj
  real(real64), dimension(ni,nj)  :: array
  real(real64), dimension(nj)     :: ywrk

  do i=1,ni  
     do j=1,nj
        ywrk(j)=array(i,j)
     end do
     do j=1,nj
        ji =nj-j+1
        array(i,ji)=ywrk(j)            
     enddo
  enddo

  return
end subroutine upside_down

end module g_read_NCEP_NetCDF
!
!------------------------------------------------------------------------------------
!
module g_read_other_NetCDF
contains
subroutine read_other_NetCDF(file, vari, itime, model_2Darray, check_dummy)
  ! Read 2D data and interpolate to the model grid.
  ! Currently used to read runoff and SSS.
  ! First, missing values are filled in on the raw regular grid;
  ! Second, interp_2d_field does the interpolation.
  ! The check_dummy part should be modified in new applications!
  ! if check_dummy=.true.,  missing value is replaced with a meaningful value nearby
  ! if check_dummy=.false., missing value is replaced with 0.0

  use, intrinsic :: ISO_FORTRAN_ENV

  use g_config
  use o_param
  use o_mesh
  use g_parsup
  implicit none

#include "netcdf.inc" 

  integer                    :: i, j, ii, jj, k, n, num, flag, cnt
  integer                    :: itime, latlen, lonlen
  integer                    :: status, ncid, varid
  integer                    :: lonid, latid
  integer                    :: istart(3), icount(3)
  real(real64)               :: x, y, miss, aux
  real(real64), allocatable  :: lon(:), lat(:)
  real(real64), allocatable  :: ncdata(:,:), ncdata_temp(:,:)
  real(real64), allocatable  :: temp_x(:), temp_y(:)
  real(real64)               :: model_2Darray(myDim_nod2d+eDim_nod2D)   
  character(15)              :: vari
  character(300)             :: file
  logical                    :: check_dummy
  integer                    :: ierror           ! return error code

  if (mype==0) then
     ! open file
     status=nf_open(file, nf_nowrite, ncid)
  end if

  call MPI_BCast(status, 1, MPI_INTEGER, 0, MPI_COMM_FESOM, ierror)
  if (status.ne.nf_noerr)then
     print*,'ERROR: CANNOT READ runoff FILE CORRECTLY !!!!!'
     print*,'Error in opening netcdf file'//file
     call par_ex
     stop
  endif

  if (mype==0) then
     ! lat
     status=nf_inq_dimid(ncid, 'lat', latid)
     status=nf_inq_dimlen(ncid, latid, latlen)
     ! lon
     status=nf_inq_dimid(ncid, 'lon', lonid)
     status=nf_inq_dimlen(ncid, lonid, lonlen)
  end if
  call MPI_BCast(latlen, 1, MPI_INTEGER, 0, MPI_COMM_FESOM, ierror)
  call MPI_BCast(lonlen, 1, MPI_INTEGER, 0, MPI_COMM_FESOM, ierror)

  ! lat
  allocate(lat(latlen))
  if (mype==0) then
     status=nf_inq_varid(ncid, 'lat', varid)
     status=nf_get_vara_double(ncid,varid,1,latlen,lat)
  end if
  call MPI_BCast(lat, latlen, MPI_DOUBLE_PRECISION, 0, MPI_COMM_FESOM, ierror)  

  ! lon
  allocate(lon(lonlen))
  if (mype==0) then
     status=nf_inq_varid(ncid, 'lon', varid)
     status=nf_get_vara_double(ncid,varid,1,lonlen,lon)
  end if
  call MPI_BCast(lon, lonlen, MPI_DOUBLE_PRECISION, 0, MPI_COMM_FESOM, ierror)  

  ! make sure range 0. - 360.
  do n=1,lonlen
     if (lon(n)<0.0) then
        lon(n)=lon(n)+360.
     end if
  end do

  allocate(ncdata(lonlen,latlen), ncdata_temp(lonlen,latlen))

  if (mype==0) then
    ! data
     status=nf_inq_varid(ncid, vari, varid)
     istart = (/1,1,itime/)
     icount= (/lonlen,latlen,1/)
     status=nf_get_vara_double(ncid,varid,istart,icount,ncdata)

    ! missing value
     status= nf_get_att_double(ncid,varid,'missing_value',miss)
    ! close file
    status=nf_close(ncid)
  end if
  call MPI_BCast(ncdata, lonlen*latlen, MPI_DOUBLE_PRECISION, 0, MPI_COMM_FESOM, ierror)
  call MPI_BCast(miss,               1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_FESOM, ierror)
  !write(*,*)'miss', miss
  !write(*,*)'raw',minval(ncdata),maxval(ncdata)
  ncdata_temp=ncdata
  do i=1,lonlen
     do j=1,latlen
        if (ncdata(i,j)==miss .or. ncdata(i,j)==-99.0_WP) then  !!
           if (check_dummy) then
              aux=0.0
              cnt=0
              do k=1,30
                 do ii=max(1,i-k),min(lonlen,i+k)
                    do jj=max(1,j-k),min(latlen,j+k)
                       if (ncdata_temp(ii,jj)/=miss .and. ncdata_temp(ii,jj)/=-99.0_WP) then  !!
                          aux=aux+ncdata_temp(ii,jj)
                          cnt=cnt+1                         
                       end if
                    end do	!ii
                 end do	!jj
                 if (cnt>0) then
                    ncdata(i,j)=aux/cnt
                    exit
                 end if
              end do  	!k    
           else
              ncdata(i,j)=0.0
           end if
        end if
     end do
  end do
  !write(*,*) 'post',minval(ncdata), maxval(ncdata)
  ! model grid coordinates
  num=myDim_nod2d+eDim_nod2d
  allocate(temp_x(num), temp_y(num))  
  do n=1, num
     temp_x(n)=geo_coord_nod2d(1,n)/rad              
     temp_y(n)=geo_coord_nod2d(2,n)/rad             
     ! change lon range to [0 360]
     if(temp_x(n)<0.) temp_x(n)=temp_x(n) + 360.0  
  end do
  ! interpolation
  flag=0
  call interp_2d_field(lonlen, latlen, lon, lat, ncdata, num, temp_x, temp_y, & 
       model_2Darray, flag) 
  deallocate(temp_y, temp_x, ncdata_temp, ncdata, lon, lat)
end subroutine read_other_NetCDF
!
!------------------------------------------------------------------------------------
!
  subroutine read_surf_hydrography_NetCDF(file, vari, itime, model_2Darray)
    ! Read WOA (NetCDF) surface T/S and interpolate to the model grid.
    ! Currently used for surface restoring in case of ocean-alone models
    ! Calling interp_2d_field_v2 to do interpolation, which also treats the dummy value.
    !
    ! Coded by Qiang Wang
    ! Reviewed by ??

  use, intrinsic :: ISO_FORTRAN_ENV

    use g_config
    use o_param
    use o_mesh
    use g_rotate_grid
    use g_parsup
    implicit none

#include "netcdf.inc" 

  integer                       :: i, j,  n, num
  integer                       :: itime, latlen, lonlen
  integer                       :: status, ncid, varid
  integer                       :: lonid, latid
  integer                       :: istart(4), icount(4)
  real(real64)                  :: x, y, miss
  real(real64), allocatable     :: lon(:), lat(:)
  real(real64), allocatable     :: ncdata(:,:)
  real(real64), allocatable     :: temp_x(:), temp_y(:)
  real(real64)                  :: model_2Darray(myDim_nod2d+eDim_nod2D)   
  character(15)                 :: vari
  character(300)                :: file
  logical                       :: check_dummy
  integer                       :: ierror           ! return error code

  if (mype==0) then
     ! open file
     status=nf_open(file, nf_nowrite, ncid)
  end if

  call MPI_BCast(status, 1, MPI_INTEGER, 0, MPI_COMM_FESOM, ierror)
  if (status.ne.nf_noerr)then
     print*,'ERROR: CANNOT READ runoff FILE CORRECTLY !!!!!'
     print*,'Error in opening netcdf file'//file
     call par_ex
     stop
  endif

  if (mype==0) then
    ! lat
    status=nf_inq_dimid(ncid, 'lat', latid)
    status=nf_inq_dimlen(ncid, latid, latlen)
    ! lon
    status=nf_inq_dimid(ncid, 'lon', lonid)
    status=nf_inq_dimlen(ncid, lonid, lonlen)
  end if

  ! lat
  allocate(lat(latlen))
  if (mype==0) then
     status=nf_inq_varid(ncid, 'lat', varid)
     status=nf_get_vara_double(ncid,varid,1,latlen,lat)
   end if
  call MPI_BCast(lat, latlen, MPI_DOUBLE_PRECISION, 0, MPI_COMM_FESOM, ierror)  

  ! lon
  allocate(lon(lonlen))
  if (mype==0) then
     status=nf_inq_varid(ncid, 'lon', varid)
     status=nf_get_vara_double(ncid,varid,1,lonlen,lon)
  end if
  call MPI_BCast(lon, lonlen, MPI_DOUBLE_PRECISION, 0, MPI_COMM_FESOM, ierror)  

  ! make sure range 0. - 360.
  do n=1,lonlen
     if (lon(n)<0.0) then
        lon(n)=lon(n)+360.
     end if
  end do

  ! data
  allocate(ncdata(lonlen,latlen))

  if (mype==0) then
     status=nf_inq_varid(ncid, vari, varid)
     istart = (/1,1,1,itime/)
     icount= (/lonlen,latlen,1,1/)
     status=nf_get_vara_double(ncid,varid,istart,icount,ncdata)

     ! missing value
     status= nf_get_att_double(ncid,varid,'missing_value',miss)
     !write(*,*)'miss', miss
     !write(*,*)'raw',minval(ncdata),maxval(ncdata)
     !close file
     status=nf_close(ncid)
  end if
  call MPI_BCast(ncdata, lonlen*latlen, MPI_DOUBLE_PRECISION, 0, MPI_COMM_FESOM, ierror)
  call MPI_BCast(miss,               1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_FESOM, ierror)
  ! the next step is to interpolate data to model grids
  ! model grid coordinates
  num=myDim_nod2d+eDim_nod2d
  allocate(temp_x(num), temp_y(num))  
  do n=1, num                        
     if (rotated_grid) then
        call r2g(x, y, coord_nod2d(1,n), coord_nod2d(2,n))
        temp_x(n)=x/rad   ! change unit to degree  
        temp_y(n)=y/rad                             
     else
        temp_x(n)=coord_nod2d(1,n)/rad              
        temp_y(n)=coord_nod2d(2,n)/rad             
     end if
     ! change lon range to [0 360]
     if(temp_x(n)<0.) temp_x(n)=temp_x(n) + 360.0  
  end do
  ! interpolation
  call interp_2d_field_v2(lonlen, latlen, lon, lat, ncdata, miss, &
       num, temp_x, temp_y, model_2Darray) 

  deallocate(temp_y, temp_x, ncdata, lon, lat)
end subroutine read_surf_hydrography_NetCDF
!
!------------------------------------------------------------------------------------
!
subroutine read_2ddata_on_grid_NetCDF(file, vari, itime, model_2Darray)  

  use, intrinsic :: ISO_FORTRAN_ENV

  use g_config
  use o_param
  use o_mesh
  use g_rotate_grid
  use g_parsup
  implicit none

#include "netcdf.inc" 

  integer                       :: n, i
  integer                       :: itime
  integer                       :: status, ncid, varid
  integer                       :: istart(2), icount(2)
  real(real64)                  :: ncdata(nod2D)
  real(real64),   intent(out)	:: model_2Darray(myDim_nod2D+eDim_nod2D)
  character(300), intent(in) 	:: file
  character(15),  intent(in)    :: vari
  integer                       :: ierror           ! return error code

  if (mype==0) then
    ! open file
    status=nf_open(file, nf_nowrite, ncid)
  end if
  call MPI_BCast(status, 1, MPI_INTEGER, 0, MPI_COMM_FESOM, ierror)
  if (status.ne.nf_noerr)then
     print*,'ERROR: CANNOT READ runoff FILE CORRECTLY !!!!!'
     print*,'Error in opening netcdf file'//file
     call par_ex
     stop
  endif

  if (mype==0) then
  ! get variables
     status=nf_inq_varid(ncid, vari, varid)
     istart = (/1, itime/)
     icount= (/nod2D, 1/)
     status=nf_get_vara_double(ncid,varid,istart,icount,ncdata)
     status=nf_close(ncid)
  end if      
  call MPI_BCast(ncdata, nod2D, MPI_DOUBLE_PRECISION, 0, MPI_COMM_FESOM, ierror)
  model_2Darray=ncdata(myList_nod2D)     
end subroutine read_2ddata_on_grid_NetCDF
  
end module g_read_other_NetCDF


