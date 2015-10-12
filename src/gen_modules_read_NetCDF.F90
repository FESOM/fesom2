! A collection of modules for reading different NetCDF files

module g_read_CORE_NetCDF
  ! read NetCDF data on T62 NCEP/NCAR grid
  !
  ! Coded by Qiang Wang
  ! Reviewed by ??
  
contains 

  subroutine read_CORE_NetCDF(file,vari,itime,ncdata)
    use g_config
    use g_parsup
    implicit none

#include "netcdf.inc" 
    
    integer, parameter             	:: nci=192, ncj=94  ! T62 grid dimension
    integer, dimension(3)           	:: istart, icount
    integer				:: ncid, varid, io
    integer                             :: i, j, itime
    real(kind=8), dimension (nci,ncj)   :: ncdata
    character(15)                       :: vari
    character(80)                      	:: file

    istart = (/1,1,itime/)
    icount= (/nci,ncj,1/)

    ! open netcdf input file
    io=nf_open(file, nf_nowrite, ncid)

    if (io.ne.nf_noerr)then
       print*,'ERROR: CANNOT READ FORCING FILE CORRECTLY !!!!!'
       print*,'Error in opening netcdf file'//file
       stop 'Fatal error in open_netcdf'
    endif

    ! get handle for data variable    
    io=nf_inq_varid(ncid, vari, varid)
    if(io  .ne. nf_noerr) call handle_err(io)
    io=nf_get_vara_double(ncid,varid,istart,icount,ncdata)
    if(io  .ne. nf_noerr) call handle_err(io)  
    io=nf_close(ncid)

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
    use g_config
    use g_parsup
    implicit none

#include "netcdf.inc" 

    integer, parameter             	:: nci=192, ncj=94  ! T62 grid dimension
    integer				:: ncid, varid, io
    integer                             :: i, j, itime
    integer, dimension(3)           	:: istart, icount
    integer(kind=2), dimension(nci,ncj) :: iuw 
    real(kind=4)                        :: xscale, xoff, miss
    real(kind=8), dimension(nci,ncj)    :: ncdata
    character(15)                       :: vari
    character(80)                      	:: file

    istart = (/1,1,itime/)
    icount= (/nci,ncj,1/)

    ! open netcdf input file
    io=nf_open(file, nf_nowrite, ncid)

    if (io.ne.nf_noerr)then
       print*,'ERROR: CANNOT READ FORCING FILE CORRECTLY !!!!!'
       print*,'Error in opening netcdf file'//file
       stop 'Fatal error in open_netcdf'
    endif

    ! get handle for data variable    
    io=nf_inq_varid(ncid, vari, varid)
    if(io  .ne. nf_noerr) call handle_err(io)
  
    ! get att
    io=nf_get_att_real(ncid,varid,'scale_factor',xscale)
    io= nf_get_att_real(ncid,varid,'add_offset',xoff)
    ! io= nf_get_att_int(ncid,varid,'missing_value',miss)

    ! get variable
    io=nf_get_vara_int2(ncid,varid,istart,icount,iuw)
    if(io  .ne. nf_noerr) call handle_err(io)
  
    ! close file 
    io=nf_close(ncid)

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
  subroutine read_NCEPv2_NetCDF(file,vari,itime,ncdata)
    use o_param
    use g_parsup
    implicit none

#include "netcdf.inc" 

    integer, parameter             	:: nci=192, ncj=94  ! T62 grid dimension
    integer				:: ncid, varid, io
    integer                             :: i, j, itime
    integer, dimension(4)           	:: istart, icount
    integer(kind=2), dimension(nci,ncj) :: iuw 
    real(kind=4)                        :: xscale, xoff, miss
    real(kind=8), dimension(nci,ncj)    :: ncdata
    character(15)                       :: vari
    character(80)                      	:: file

    istart = (/1,1,1,itime/)
    icount= (/nci,ncj,1,1/)

    ! open netcdf input file
    io=nf_open(file, nf_nowrite, ncid)

    if (io.ne.nf_noerr)then
       print*,'ERROR: CANNOT READ FORCING FILE CORRECTLY !!!!!'
       print*,'Error in opening netcdf file'//file
       stop 'Fatal error in open_netcdf'
    endif

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

    ! ncdata
    do j=1,ncj
       do i=1,nci
          ncdata(i,j)=real(iuw(i,j))*xscale+xoff
       enddo
    end do
    return
  end subroutine read_NCEPv2_NetCDF
  !
  !--------------------------------------------------------------
  !
  subroutine upside_down(array,ni,nj)
    implicit none

    integer                         :: i, j, ji
    integer                         :: ni, nj
    real(kind=8), dimension(ni,nj)  :: array
    real(kind=8), dimension(nj)     :: ywrk

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
  ! Read global data with NetCDF format and interpolate to the model grid.
  !
  ! Coded by Qiang Wang
  ! Reviewed by ??

contains

  subroutine read_other_NetCDF(file, vari, itime, model_2Darray, check_dummy)
    ! Read 2D data and interpolate to the model grid.
    ! Currently used to read runoff and SSS.
    ! First, missing values are filled in on the raw regular grid;
    ! Second, interp_2d_field does the interpolation.

    ! The check_dummy part should be modified in new applications!
    ! if check_dummy=.true.,  missing value is replaced with a meaningful value nearby
    ! if check_dummy=.false., missing value is replaced with 0.0

    use g_config
    use o_param
    use o_mesh
    use g_parsup
    implicit none

#include "netcdf.inc" 

    integer			:: i, j, ii, jj, k, n, num, flag, cnt
    integer			:: itime, latlen, lonlen
    integer			:: status, ncid, varid
    integer			:: lonid, latid
    integer			:: istart(3), icount(3)
    real(kind=8)		:: x, y, miss, aux
    real(kind=8), allocatable	:: lon(:), lat(:)
    real(kind=8), allocatable	:: ncdata(:,:), ncdata_temp(:,:)
    real(kind=8), allocatable	:: temp_x(:), temp_y(:)
    real(kind=8)		:: model_2Darray(myDim_nod2d+eDim_nod2D)   
    character(15)		:: vari
    character(80)              	:: file
    logical                     :: check_dummy

    ! open file
    status=nf_open(file, nf_nowrite, ncid)
    if (status.ne.nf_noerr)then
       print*,'ERROR: CANNOT READ runoff FILE CORRECTLY !!!!!'
       print*,'Error in opening netcdf file'//file
       stop
    endif

    ! lat
    status=nf_inq_dimid(ncid, 'lat', latid)
    status=nf_inq_dimlen(ncid, latid, latlen)
    allocate(lat(latlen))
    status=nf_inq_varid(ncid, 'lat', varid)
    status=nf_get_vara_double(ncid,varid,1,latlen,lat)

    ! lon
    status=nf_inq_dimid(ncid, 'lon', lonid)
    status=nf_inq_dimlen(ncid, lonid, lonlen)
    allocate(lon(lonlen))
    status=nf_inq_varid(ncid, 'lon', varid)
    status=nf_get_vara_double(ncid,varid,1,lonlen,lon)
    ! make sure range 0. - 360.
    do n=1,lonlen
       if(lon(n)<0.0) then
          lon(n)=lon(n)+360.
       end if
    end do

    ! data
    allocate(ncdata(lonlen,latlen), ncdata_temp(lonlen,latlen))
    status=nf_inq_varid(ncid, vari, varid)
    istart = (/1,1,itime/)
    icount= (/lonlen,latlen,1/)
    status=nf_get_vara_double(ncid,varid,istart,icount,ncdata)

    ! missing value
    status= nf_get_att_double(ncid,varid,'missing_value',miss)
    ncdata_temp=ncdata
    do i=1,lonlen
       do j=1,latlen
          if(ncdata(i,j)==miss .or. ncdata(i,j)==-99.0_8) then  !!
             if(check_dummy) then
                aux=0.0
                cnt=0
                do k=1,30
                   do ii=max(1,i-k),min(lonlen,i+k)
                      do jj=max(1,j-k),min(latlen,j+k)
                         if(ncdata_temp(ii,jj)/=miss .and. ncdata_temp(ii,jj)/=-99.0_8) then  !!
                            aux=aux+ncdata_temp(ii,jj)
                            cnt=cnt+1                         
                         end if
                      end do	!ii
                   end do	!jj
                   if(cnt>0) then
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

    ! close file
    status=nf_close(ncid)
    
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

    use g_config
    use o_param
    use o_mesh
    use g_rotate_grid
    use g_parsup
    implicit none

#include "netcdf.inc" 

    integer			:: i, j,  n, num
    integer			:: itime, latlen, lonlen
    integer			:: status, ncid, varid
    integer			:: lonid, latid
    integer			:: istart(4), icount(4)
    real(kind=8)		:: x, y, miss
    real(kind=8), allocatable	:: lon(:), lat(:)
    real(kind=8), allocatable	:: ncdata(:,:)
    real(kind=8), allocatable	:: temp_x(:), temp_y(:)
    real(kind=8)		:: model_2Darray(myDim_nod2d+eDim_nod2D)   
    character(15)		:: vari
    character(80)              	:: file
    logical                     :: check_dummy

    ! open file
    status=nf_open(file, nf_nowrite, ncid)
    if (status.ne.nf_noerr)then
       print*,'ERROR: CANNOT READ runoff FILE CORRECTLY !!!!!'
       print*,'Error in opening netcdf file'//file
       stop
    endif

    ! lat
    status=nf_inq_dimid(ncid, 'lat', latid)
    status=nf_inq_dimlen(ncid, latid, latlen)
    allocate(lat(latlen))
    status=nf_inq_varid(ncid, 'lat', varid)
    status=nf_get_vara_double(ncid,varid,1,latlen,lat)

    ! lon
    status=nf_inq_dimid(ncid, 'lon', lonid)
    status=nf_inq_dimlen(ncid, lonid, lonlen)
    allocate(lon(lonlen))
    status=nf_inq_varid(ncid, 'lon', varid)
    status=nf_get_vara_double(ncid,varid,1,lonlen,lon)
    ! make sure range 0. - 360.
    do n=1,lonlen
       if(lon(n)<0.0) then
          lon(n)=lon(n)+360.
       end if
    end do

    ! data
    allocate(ncdata(lonlen,latlen))
    status=nf_inq_varid(ncid, vari, varid)
    istart = (/1,1,1,itime/)
    icount= (/lonlen,latlen,1,1/)
    status=nf_get_vara_double(ncid,varid,istart,icount,ncdata)

    ! missing value
    status= nf_get_att_double(ncid,varid,'missing_value',miss)

    ! close file
    status=nf_close(ncid)

    ! the next step is to interpolate data to model grids

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
    call interp_2d_field_v2(lonlen, latlen, lon, lat, ncdata, miss, &
         num, temp_x, temp_y, model_2Darray) 

    deallocate(temp_y, temp_x, ncdata, lon, lat)

  end subroutine read_surf_hydrography_NetCDF
  !
  !------------------------------------------------------------------------------------
  !
  subroutine read_2ddata_on_grid_NetCDF(file, vari, itime, model_2Darray)  
    ! read 2D data which are already prepared on the model grid
    !
    ! Coded by Qiang Wang
    ! Reviewed by ??
    
    use g_config
    use o_param
    use o_mesh
    use g_rotate_grid
    use g_parsup
    implicit none

#include "netcdf.inc" 

    integer			:: n, i
    integer			:: itime
    integer			:: status, ncid, varid
    integer			:: istart(2), icount(2)
    real(kind=8)           	:: ncdata(nod2D)
    real(kind=8), intent(out)	:: model_2Darray(myDim_nod2d+eDim_nod2D)
    character(80), intent(in) 	:: file
    character(15), intent(in)  :: vari

    ! open file
    status=nf_open(file, nf_nowrite, ncid)
    if (status.ne.nf_noerr)then
       print*,'ERROR: CANNOT READ runoff FILE CORRECTLY !!!!!'
       print*,'Error in opening netcdf file'//file
       stop
    endif

    ! get variables
    status=nf_inq_varid(ncid, vari, varid)
    istart = (/1, itime/)
    icount= (/nod2D, 1/)
    status=nf_get_vara_double(ncid,varid,istart,icount,ncdata)
    status=nf_close(ncid)
      
    model_2Darray=ncdata(myList_nod2D)     

  end subroutine read_2ddata_on_grid_NetCDF
  
end module g_read_other_NetCDF


