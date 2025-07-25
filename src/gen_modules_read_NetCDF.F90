!
!------------------------------------------------------------------------------------
!
module g_read_other_NetCDF
contains
subroutine read_other_NetCDF(file, vari, itime, model_2Darray, check_dummy, do_onvert, partit, mesh)
  ! Read 2D data and interpolate to the model grid.
  ! Currently used to read runoff and SSS.
  ! First, missing values are filled in on the raw regular grid;
  ! Second, interp_2d_field does the interpolation.
  ! The check_dummy part should be modified in new applications!
  ! if check_dummy=.true.,  missing value is replaced with a meaningful value nearby
  ! if check_dummy=.false., missing value is replaced with 0.0

  use, intrinsic :: ISO_FORTRAN_ENV, only: real64
  use g_config
  use o_param
  USE MOD_MESH
  USE MOD_PARTIT
  USE MOD_PARSUP
  use netcdf
  implicit none

  type(t_mesh),   intent(in),    target :: mesh
  type(t_partit), intent(inout), target :: partit
  integer                    :: i, j, ii, jj, k, n, num, flag, cnt
  integer                    :: itime, latlen, lonlen
  integer                    :: status, ncid, varid
  integer                    :: lonid, latid
  integer                    :: istart(3), icount(3), elnodes(3)
  real(real64)               :: x, y, miss, aux, xmin, elnodes_x(3)
  real(real64), allocatable  :: lon(:), lat(:)
  real(real64), allocatable  :: ncdata(:,:), ncdata_temp(:,:)
  real(real64), allocatable  :: temp_x(:), temp_y(:)
  real(real64)               :: model_2Darray(partit%myDim_nod2d+partit%eDim_nod2D)   
  character(*)               :: vari
  character(*)               :: file
  logical                    :: check_dummy, do_onvert
  integer                    :: ierror           ! return error code

#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"


  if (mype==0) then
     ! open file
     status=nf90_open(file, nf90_nowrite, ncid)
  end if

  call MPI_BCast(status, 1, MPI_INTEGER, 0, MPI_COMM_FESOM, ierror)
  if (status.ne.nf90_noerr)then
     print*,'ERROR: CANNOT READ 2D netCDF FILE CORRECTLY !!!!!'
     print*,'Error in opening netcdf file '//file
     call par_ex(partit%MPI_COMM_FESOM, partit%mype)
     stop
  endif

  if (mype==0) then
     ! lat
     status=nf90_inq_dimid(ncid, 'lat', latid)
     status=nf90_inquire_dimension(ncid, latid, len=latlen)
     ! lon
     status=nf90_inq_dimid(ncid, 'lon', lonid)
     status=nf90_inquire_dimension(ncid, lonid, len=lonlen)
  end if
  call MPI_BCast(latlen, 1, MPI_INTEGER, 0, MPI_COMM_FESOM, ierror)
  call MPI_BCast(lonlen, 1, MPI_INTEGER, 0, MPI_COMM_FESOM, ierror)

  ! lat
  allocate(lat(latlen))
  if (mype==0) then
     status=nf90_inq_varid(ncid, 'lat', varid)
     status=nf90_get_var(ncid, varid, lat, start=(/1/), count=(/latlen/))
  end if
  call MPI_BCast(lat, latlen, MPI_DOUBLE_PRECISION, 0, MPI_COMM_FESOM, ierror)  

  ! lon
  allocate(lon(lonlen))
  if (mype==0) then
     status=nf90_inq_varid(ncid, 'lon', varid)
     status=nf90_get_var(ncid, varid, lon, start=(/1/), count=(/lonlen/))
  end if
  call MPI_BCast(lon, lonlen, MPI_DOUBLE_PRECISION, 0, MPI_COMM_FESOM, ierror)  

  ! make sure range 0. - 360.
  do n=1,lonlen
     if (lon(n)<0.0_WP) then
        lon(n)=lon(n)+360._WP
     end if
  end do

  allocate(ncdata(lonlen,latlen), ncdata_temp(lonlen,latlen))
  ncdata = 0.0_WP
  
  if (mype==0) then
    ! data
     status=nf90_inq_varid(ncid, vari, varid)
     istart = (/1,1,itime/)
     icount= (/lonlen,latlen,1/)
     status=nf90_get_var(ncid, varid, ncdata, start=istart, count=icount)

    ! missing value
     status= nf90_get_att(ncid, varid, 'missing_value', miss)
    ! close file
    status=nf90_close(ncid)
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
              aux=0.0_WP
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
              ncdata(i,j)=0.0_WP
           end if
        end if
     end do
  end do
    !write(*,*) 'post',minval(ncdata), maxval(ncdata)
    
    !___________________________________________________________________________
    ! create interpolation coordinates
    ! do data interpolation on vertices
    if (do_onvert) then
        num=myDim_nod2d+eDim_nod2d
        allocate(temp_x(num), temp_y(num))  
        do n=1, num
            temp_x(n)=geo_coord_nod2d(1,n)/rad              
            temp_y(n)=geo_coord_nod2d(2,n)/rad             
            ! change lon range to [0 360]
            if(temp_x(n)<0._WP) temp_x(n)=temp_x(n) + 360.0_WP  
        end do
        
    ! do data interpolation on element centroids
    else
        num = myDim_elem2D
        allocate(temp_x(num), temp_y(num))  
        do n=1, num
            ! compute points of element centroids in geo frame use them here for interpolation
            elnodes  = elem2D_nodes(:,n)
            elnodes_x= geo_coord_nod2D(1, elnodes)
            xmin     = minval(elnodes_x)
            do k=1,3
                if(elnodes_x(k)-xmin>=cyclic_length/2.0_WP) elnodes_x(k)=elnodes_x(k)-cyclic_length
                if(elnodes_x(k)-xmin<-cyclic_length/2.0_WP) elnodes_x(k)=elnodes_x(k)+cyclic_length
            end do
            ! compute in units [deg], in geo frame
            temp_x(n)=sum(elnodes_x)/3.0_WP/rad
            temp_y(n)=sum(geo_coord_nod2D(2,elnodes))/3.0_WP/rad
            
            ! change lon range to [0 360]
            if(temp_x(n)<0._WP) temp_x(n)=temp_x(n) + 360.0_WP  
        end do
    end if   
    
    !___________________________________________________________________________
    ! do interpolation
    flag=0
    call interp_2d_field(lonlen, latlen, lon, lat, ncdata, num, temp_x, temp_y, & 
                         model_2Darray, flag, partit) 
    deallocate(temp_y, temp_x, ncdata_temp, ncdata, lon, lat)
    
end subroutine read_other_NetCDF
!
!------------------------------------------------------------------------------------
!
subroutine read_surf_hydrography_NetCDF(file, vari, itime, model_2Darray, partit, mesh)
    ! Read WOA (NetCDF) surface T/S and interpolate to the model grid.
    ! Currently used for surface restoring in case of ocean-alone models
    ! Calling interp_2d_field_v2 to do interpolation, which also treats the dummy value.
    !
    ! Coded by Qiang Wang
    ! Reviewed by ??
    use g_config
    use o_param
    USE MOD_MESH
    USE MOD_PARTIT
    USE MOD_PARSUP
    use g_rotate_grid
    use, intrinsic :: ISO_FORTRAN_ENV, only: real64
    use netcdf
    implicit none

  type(t_mesh),   intent(in),    target :: mesh
  type(t_partit), intent(inout), target :: partit
  integer                       :: i, j,  n, num
  integer                       :: itime, latlen, lonlen
  integer                       :: status, ncid, varid
  integer                       :: lonid, latid, drain_num
  integer                       :: istart(4), icount(4)
  real(real64)                  :: x, y, miss
  real(real64), allocatable     :: lon(:), lat(:)
  real(real64), allocatable     :: ncdata(:,:)
  real(real64), allocatable     :: temp_x(:), temp_y(:)
  real(real64)                  :: model_2Darray(partit%myDim_nod2d+partit%eDim_nod2D)   
  character(15)                 :: vari
  character(300)                :: file
  logical                       :: check_dummy
  integer                       :: ierror           ! return error code

#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"

  if (mype==0) then
     ! open file
     status=nf90_open(file, nf90_nowrite, ncid)
  end if

  call MPI_BCast(status, 1, MPI_INTEGER, 0, MPI_COMM_FESOM, ierror)
  if (status.ne.nf90_noerr)then
     print*,'ERROR: CANNOT READ runoff FILE CORRECTLY !!!!!'
     print*,'Error in opening netcdf file '//file
     call par_ex(partit%MPI_COMM_FESOM, partit%mype)
     stop
  endif

  if (mype==0) then
    ! lat
    status=nf90_inq_dimid(ncid, 'lat', latid)
    status=nf90_inquire_dimension(ncid, latid, len=latlen)
    ! lon
    status=nf90_inq_dimid(ncid, 'lon', lonid)
    status=nf90_inquire_dimension(ncid, lonid, len=lonlen)
  end if

  ! lat
  allocate(lat(latlen))
  if (mype==0) then
     status=nf90_inq_varid(ncid, 'lat', varid)
     status=nf90_get_var(ncid, varid, lat, start=(/1/), count=(/latlen/))
   end if
  call MPI_BCast(lat, latlen, MPI_DOUBLE_PRECISION, 0, MPI_COMM_FESOM, ierror)  

  ! lon
  allocate(lon(lonlen))
  if (mype==0) then
     status=nf90_inq_varid(ncid, 'lon', varid)
     status=nf90_get_var(ncid, varid, lon, start=(/1/), count=(/lonlen/))
  end if
  call MPI_BCast(lon, lonlen, MPI_DOUBLE_PRECISION, 0, MPI_COMM_FESOM, ierror)  

  ! make sure range 0. - 360.
  do n=1,lonlen
     if (lon(n)<0.0_WP) then
        lon(n)=lon(n)+360._WP
     end if
  end do

  ! data
  allocate(ncdata(lonlen,latlen))
  ncdata = 0.0_WP

  if (mype==0) then
     status=nf90_inq_varid(ncid, vari, varid)
     istart = (/1,1,1,itime/)
     icount= (/lonlen,latlen,1,1/)
     status=nf90_get_var(ncid, varid, ncdata, start=istart, count=icount)

     ! missing value
     status= nf90_get_att(ncid, varid, 'missing_value', miss)
     !write(*,*)'miss', miss
     !write(*,*)'raw',minval(ncdata),maxval(ncdata)
     !close file
     status=nf90_close(ncid)
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
     if(temp_x(n)<0._WP) temp_x(n)=temp_x(n) + 360.0_WP
  end do
  ! interpolation
  call interp_2d_field_v2(lonlen, latlen, lon, lat, ncdata, miss, &
       num, temp_x, temp_y, model_2Darray, partit) 
  deallocate(temp_y, temp_x, ncdata, lon, lat)
end subroutine read_surf_hydrography_NetCDF
!
!------------------------------------------------------------------------------------
!
subroutine read_2ddata_on_grid_NetCDF(file, vari, itime, model_2Darray, partit, mesh)  

  use, intrinsic :: ISO_FORTRAN_ENV, only: real64

  use g_config
  use o_param
  USE MOD_MESH
  USE MOD_PARTIT
  USE MOD_PARSUP
  use g_rotate_grid
  use netcdf
  implicit none

  type(t_mesh),   intent(in)   , target :: mesh
  type(t_partit), intent(inout), target :: partit
  integer                       :: n, i
  integer                       :: itime
  integer                       :: status, ncid, varid
  integer                       :: istart(2), icount(2)
  real(real64)                  :: ncdata(mesh%nod2D)
  real(real64),   intent(out)	:: model_2Darray(partit%myDim_nod2D+partit%eDim_nod2D)
  character(*),  intent(in) 	:: file
  character(*),  intent(in)     :: vari
  integer                       :: ierror           ! return error code

#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"

  if (mype==0) then
    ! open file
    status=nf90_open(file, nf90_nowrite, ncid)
  end if
  call MPI_BCast(status, 1, MPI_INTEGER, 0, MPI_COMM_FESOM, ierror)
  if (status.ne.nf90_noerr)then
     print*,'ERROR: CANNOT READ runoff FILE CORRECTLY !!!!!'
     print*,'Error in opening netcdf file '//file
     call par_ex(partit%MPI_COMM_FESOM, partit%mype)
     stop
  endif

  if (mype==0) then
  ! get variables
     status=nf90_inq_varid(ncid, vari, varid)
     istart = (/1, itime/)
     icount= (/nod2D, 1/)
     status=nf90_get_var(ncid, varid, ncdata, start=istart, count=icount)
     status=nf90_close(ncid)
  end if      
  call MPI_BCast(ncdata, nod2D, MPI_DOUBLE_PRECISION, 0, MPI_COMM_FESOM, ierror)
  model_2Darray=ncdata(myList_nod2D) 
end subroutine read_2ddata_on_grid_NetCDF
end module g_read_other_NetCDF
