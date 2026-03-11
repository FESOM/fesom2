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
    use g_interp
    use netcdf
    implicit none

    !___Input/Output____________________________________________________________
    character(*)  , intent(in)            :: file
    character(*)  , intent(in)            :: vari
    integer       , intent(in)            :: itime
    real(kind=WP) , intent(inout)         :: model_2Darray(:)
    logical       , intent(in)            :: check_dummy
    logical       , intent(in)            :: do_onvert
    type(t_mesh)  , intent(in)   , target :: mesh
    type(t_partit), intent(inout), target :: partit
    !___Local___________________________________________________________________
    integer                               :: i, j, ii, jj, k, n, num, flag, cnt
    integer                               :: latlen, lonlen
    integer                               :: status, ncid, varid
    integer                               :: lonid, latid
    integer                               :: istart(3), icount(3), elnodes(3)
    real(kind=WP)                         :: x, y, miss, aux, xmin, elnodes_x(3)
    real(kind=WP), allocatable            :: lon(:), lat(:)
    real(kind=WP), allocatable            :: ncdata(:,:), ncdata_temp(:,:)
    real(kind=WP), allocatable            :: temp_x(:), temp_y(:)
    logical                               :: found_error=.False.
    integer                               :: ierror           ! return error code

#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"

    !___________________________________________________________________________
    ! open netcdf file, check existence
    if (mype==0) then
        ! open file
        status=nf90_open(file, nf90_nowrite, ncid)
        if (status /= nf90_noerr) then
            print *, achar(27)//'[33m'
            write(*,*) '____________________________________________________________________'
            write(*,*) ' ERROR: File not found: '
            write(*,*) '        ├> file:', file
            write(*,*) '        ├> status:', status
            write(*,*) '        ├> ', nf90_strerror(status)
            write(*,*) '        └> check: namelist.*'
            write(*,*) '____________________________________________________________________'
            print *, achar(27)//'[0m'
            write(*,*)
            found_error = .True.
        end if
    end if
    
    ! broadcast found_error variable to all other processes --> create kill signal
    ! for all processes
    call MPI_BCast(found_error, 1, MPI_LOGICAL, 0, MPI_COMM_FESOM, ierror)
    if (found_error) call par_ex(partit%MPI_COMM_FESOM, partit%mype)

    !___________________________________________________________________________
    ! read dimensions from netcdf file, checkfor the naming of dimensions
    if (mype==0) then
        ! lat
        status=nf90_inq_dimid(ncid, 'lat', latid)
        if (status == NF90_EBADDIM) then
            print *, achar(27)//'[33m'
            write(*,*) '____________________________________________________________________'
            write(*,*) ' ERROR: Dimension lat not found in file: '
            write(*,*) '        ├> file:', file
            write(*,*) '        ├> status:', status
            write(*,*) '        ├> ',nf90_strerror(status)
            write(*,*) '        └> check: Dim and CoordVariable name must be lat'
            write(*,*) '____________________________________________________________________'
            print *, achar(27)//'[0m'
            write(*,*)
            found_error = .True.
        end if
        status=nf90_inquire_dimension(ncid, latid, len=latlen)
        
        ! lon
        status=nf90_inq_dimid(ncid, 'lon', lonid)
        if (status == NF90_EBADDIM) then
            print *, achar(27)//'[33m'
            write(*,*) '____________________________________________________________________'
            write(*,*) ' ERROR: Dimension lon not found in file: '
            write(*,*) '        ├> file:', file
            write(*,*) '        ├> status:', status
            write(*,*) '        ├> ',nf90_strerror(status)
            write(*,*) '        └> check: Dim and CoordVariable name must be lon'
            write(*,*) '____________________________________________________________________'
            print *, achar(27)//'[0m'
            write(*,*)
            found_error = .True.
        end if
        status=nf90_inquire_dimension(ncid, lonid, len=lonlen)
    end if
    ! broadcast found_error variable to all other processes --> create kill signal
    ! for all processes
    call MPI_BCast(found_error, 1, MPI_LOGICAL, 0, MPI_COMM_FESOM, ierror)
    if (found_error) call par_ex(partit%MPI_COMM_FESOM, partit%mype)
    
    ! broadcast len of lon/lat to all processes
    call MPI_BCast(latlen, 1, MPI_INTEGER, 0, MPI_COMM_FESOM, ierror)
    call MPI_BCast(lonlen, 1, MPI_INTEGER, 0, MPI_COMM_FESOM, ierror)

    !___________________________________________________________________________
    ! read latitude regular coordinates from netcdf file --> lat
    allocate(lat(latlen))
    if (mype==0) then
        status=nf90_inq_varid(ncid, 'lat', varid)
        if (status == NF90_ENOTVAR) then
            print *, achar(27)//'[33m'
            write(*,*) '____________________________________________________________________'
            write(*,*) ' ERROR: Variable lat not found in file: '
            write(*,*) '        ├> file:', file
            write(*,*) '        ├> status:', status
            write(*,*) '        ├> ',nf90_strerror(status)
            write(*,*) '        └> check: Dim and CoordVariable name must be lat'
            write(*,*) '____________________________________________________________________'
            print *, achar(27)//'[0m'
            write(*,*)
            found_error = .True.
        end if
        status=nf90_get_var(ncid, varid, lat, start=(/1/), count=(/latlen/))
    end if
    ! broadcast found_error variable to all other processes --> create kill signal
    ! for all processes
    call MPI_BCast(found_error, 1, MPI_LOGICAL, 0, MPI_COMM_FESOM, ierror)
    if (found_error) call par_ex(partit%MPI_COMM_FESOM, partit%mype)
    
    ! broadcast latitude array to all processes
    call MPI_BCast(lat, latlen, MPI_DOUBLE_PRECISION, 0, MPI_COMM_FESOM, ierror)  

    !___________________________________________________________________________
    ! read longitude regular coordinates from netcdf file --> lat
    allocate(lon(lonlen))
    if (mype==0) then
        status=nf90_inq_varid(ncid, 'lon', varid)
        if (status == NF90_ENOTVAR) then
            print *, achar(27)//'[33m'
            write(*,*) '____________________________________________________________________'
            write(*,*) ' ERROR: Variable lon not found in file: '
            write(*,*) '        ├> file:', file
            write(*,*) '        ├> status:', status
            write(*,*) '        ├> ',nf90_strerror(status)
            write(*,*) '        └> check: Dim and CoordVariable name must be lon'
            write(*,*) '____________________________________________________________________'
            print *, achar(27)//'[0m'
            write(*,*)
            found_error = .True.
        end if
        status=nf90_get_var(ncid, varid, lon, start=(/1/), count=(/lonlen/))
    end if
    ! broadcast found_error variable to all other processes --> create kill signal
    ! for all processes
    call MPI_BCast(found_error, 1, MPI_LOGICAL, 0, MPI_COMM_FESOM, ierror)
    if (found_error) call par_ex(partit%MPI_COMM_FESOM, partit%mype)
    
    ! broadcast longitude array to all processes
    call MPI_BCast(lon, lonlen, MPI_DOUBLE_PRECISION, 0, MPI_COMM_FESOM, ierror)  

    ! make sure range 0. - 360.
    do n=1,lonlen
        if (lon(n)<0.0_WP) then
            lon(n)=lon(n)+360._WP
        end if
    end do

    !___________________________________________________________________________
    ! read the 2d data --> check if variable name does exist
    allocate(ncdata(lonlen,latlen), ncdata_temp(lonlen,latlen))
    ncdata = 0.0_WP
    if (mype==0) then
        ! data
        status=nf90_inq_varid(ncid, vari, varid)
        if (status == NF90_ENOTVAR) then
            print *, achar(27)//'[33m'
            write(*,*) '____________________________________________________________________'
            write(*,*) ' ERROR: Variable ',vari,' not found in file: '
            write(*,*) '        ├> file:', file
            write(*,*) '        ├> status:', status
            write(*,*) '        ├> ',nf90_strerror(status)
            write(*,*) '        └> check: namelist.*   '
            write(*,*) '____________________________________________________________________'
            print *, achar(27)//'[0m'
            write(*,*)
            found_error = .True.
            call par_ex(partit%MPI_COMM_FESOM, partit%mype, 0)
        end if
        istart= (/1,1,itime/)
        icount= (/lonlen,latlen,1/)
        status=nf90_get_var(ncid, varid, ncdata, start=istart, count=icount)
        
        ! missing value
        status= nf90_get_att(ncid, varid, 'missing_value', miss)
        
        ! close file
        status=nf90_close(ncid)
    end if
    ! broadcast found_error variable to all other processes --> create kill signal
    ! for all processes
    call MPI_BCast(found_error, 1, MPI_LOGICAL, 0, MPI_COMM_FESOM, ierror)
    if (found_error) call par_ex(partit%MPI_COMM_FESOM, partit%mype)
    
    ! broadcast 2d data array to all processes
    call MPI_BCast(ncdata, lonlen*latlen, MPI_DOUBLE_PRECISION, 0, MPI_COMM_FESOM, ierror)
    call MPI_BCast(miss  ,             1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_FESOM, ierror)
    
    ncdata_temp=ncdata
    do i=1,lonlen
        do j=1,latlen
            ! check for missing value
            if (ncdata(i,j)==miss .or. ncdata(i,j)==-99.0_WP) then
                ! check_dummy=.true.,  missing value is replaced with a 
                ! meaningful value nearby make sure that near coastal 
                ! interpoaltion is valid
                if (check_dummy) then
                    aux=0.0_WP
                    cnt=0
                    ! look for up to 30 grid boxes nearby
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
                    end do ! --> do k=1,30
                else
                    ncdata(i,j)=0.0_WP
                end if ! --> if (check_dummy) then
            end if ! --> if (ncdata(i,j)==miss .or. ncdata(i,j)==-99.0_WP) then
        end do ! --> do j=1,latlen
    end do ! --> do i=1,lonlen
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
                         model_2Darray(1:num), flag, partit) 
    deallocate(temp_y, temp_x, ncdata_temp, ncdata, lon, lat)
    
end subroutine read_other_NetCDF
!
!
!_______________________________________________________________________________
subroutine read_other_NetCDF_3d(file, vname, zvname, model_3Darray, do_onvert, partit, mesh)
    ! Read data over three dimensions e2g. time or phit (mist be dimension without 
    ! topography or lsmask information) and interpolate to the model grid.
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
    use g_interp
    use netcdf
    implicit none

    !___Input/Output____________________________________________________________
    character(*)  , intent(in)            :: file
    character(*)  , intent(in)            :: vname
    character(*)  , intent(in)            :: zvname
    real(kind=WP) , intent(inout)         :: model_3Darray(:,:)
    logical       , intent(in)            :: do_onvert
    type(t_mesh)  , intent(in)   , target :: mesh
    type(t_partit), intent(inout), target :: partit
    !___Local___________________________________________________________________
    integer                               :: i, j, ii, jj, k, n, num, flag, cnt
    integer                               :: latlen, lonlen, zlen, nz
    integer                               :: status, ncid, varid
    integer                               :: lonid, latid, zid
    integer                               :: istart(3), icount(3), elnodes(3)
    real(kind=WP)                         :: x, y, missvalue, aux, xmin, elnodes_x(3)
    real(kind=WP), allocatable            :: lon(:), lat(:) !, z(:)
    real(kind=WP), allocatable            :: ncdata(:,:,:)
    real(kind=WP), allocatable            :: temp_x(:), temp_y(:)
    logical                               :: found_error=.False.
    integer                               :: ierror           ! return error code

#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"

    !___________________________________________________________________________
    ! open netcdf file, check existence
    if (mype==0) then
        ! open file
        status=nf90_open(file, nf90_nowrite, ncid)
        if (status /= nf90_noerr) then
            print *, achar(27)//'[33m'
            write(*,*) '____________________________________________________________________'
            write(*,*) ' ERROR: File not found: '
            write(*,*) '        ├> file:', file
            write(*,*) '        ├> status:', status
            write(*,*) '        ├> ', nf90_strerror(status)
            write(*,*) '        └> check: namelist.*'
            write(*,*) '____________________________________________________________________'
            print *, achar(27)//'[0m'
            write(*,*)
            found_error = .True.
        end if
    end if
    
    ! broadcast found_error variable to all other processes --> create kill signal
    ! for all processes
    call MPI_BCast(found_error, 1, MPI_LOGICAL, 0, MPI_COMM_FESOM, ierror)
    if (found_error) call par_ex(partit%MPI_COMM_FESOM, partit%mype)

    !___________________________________________________________________________
    ! read dimensions from netcdf file, checkfor the naming of dimensions
    if (mype==0) then
        ! lat
        status=nf90_inq_dimid(ncid, 'lat', latid)
        if (status == NF90_EBADDIM) then
            print *, achar(27)//'[33m'
            write(*,*) '____________________________________________________________________'
            write(*,*) ' ERROR: Dimension lat not found in file: '
            write(*,*) '        ├> file:', file
            write(*,*) '        ├> status:', status
            write(*,*) '        ├> ',nf90_strerror(status)
            write(*,*) '        └> check: Dim and CoordVariable name must be lat'
            write(*,*) '____________________________________________________________________'
            print *, achar(27)//'[0m'
            write(*,*)
            found_error = .True.
        end if
        status=nf90_inquire_dimension(ncid, latid, len=latlen)
        
        ! lon
        status=nf90_inq_dimid(ncid, 'lon', lonid)
        if (status == NF90_EBADDIM) then
            print *, achar(27)//'[33m'
            write(*,*) '____________________________________________________________________'
            write(*,*) ' ERROR: Dimension lon not found in file: '
            write(*,*) '        ├> file:', file
            write(*,*) '        ├> status:', status
            write(*,*) '        ├> ',nf90_strerror(status)
            write(*,*) '        └> check: Dim and CoordVariable name must be lon'
            write(*,*) '____________________________________________________________________'
            print *, achar(27)//'[0m'
            write(*,*)
            found_error = .True.
        end if
        status=nf90_inquire_dimension(ncid, lonid, len=lonlen)
        
        ! read 3rd z dimension in file, could be time, or phit (idemix2) ... 
        status=nf90_inq_dimid(ncid, zvname, zid)
        if (status == NF90_EBADDIM) then
            print *, achar(27)//'[33m'
            write(*,*) '____________________________________________________________________'
            write(*,*) ' ERROR: Dimension ',zvname, ' not found in file: '
            write(*,*) '        ├> file:', file
            write(*,*) '        ├> status:', status
            write(*,*) '        ├> ',nf90_strerror(status)
            write(*,*) '        └> check: Dim and CoordVariable name must be lon'
            write(*,*) '____________________________________________________________________'
            print *, achar(27)//'[0m'
            write(*,*)
            found_error = .True.
        end if
        status=nf90_inquire_dimension(ncid, zid, len=zlen)
    end if
    ! broadcast found_error variable to all other processes --> create kill signal
    ! for all processes
    call MPI_BCast(found_error, 1, MPI_LOGICAL, 0, MPI_COMM_FESOM, ierror)
    if (found_error) call par_ex(partit%MPI_COMM_FESOM, partit%mype)
    
    ! broadcast len of lon/lat to all processes
    call MPI_BCast(latlen, 1, MPI_INTEGER, 0, MPI_COMM_FESOM, ierror)
    call MPI_BCast(lonlen, 1, MPI_INTEGER, 0, MPI_COMM_FESOM, ierror)
    call MPI_BCast(zlen  , 1, MPI_INTEGER, 0, MPI_COMM_FESOM, ierror)

    !___________________________________________________________________________
    ! read latitude regular coordinates from netcdf file --> lat
    allocate(lat(latlen))
    if (mype==0) then
        status=nf90_inq_varid(ncid, 'lat', varid)
        if (status == NF90_ENOTVAR) then
            print *, achar(27)//'[33m'
            write(*,*) '____________________________________________________________________'
            write(*,*) ' ERROR: Variable lat not found in file: '
            write(*,*) '        ├> file:', file
            write(*,*) '        ├> status:', status
            write(*,*) '        ├> ',nf90_strerror(status)
            write(*,*) '        └> check: Dim and CoordVariable name must be lat'
            write(*,*) '____________________________________________________________________'
            print *, achar(27)//'[0m'
            write(*,*)
            found_error = .True.
        end if
        status=nf90_get_var(ncid, varid, lat, start=(/1/), count=(/latlen/))
    end if
    ! broadcast found_error variable to all other processes --> create kill signal
    ! for all processes
    call MPI_BCast(found_error, 1, MPI_LOGICAL, 0, MPI_COMM_FESOM, ierror)
    if (found_error) call par_ex(partit%MPI_COMM_FESOM, partit%mype)
    
    ! broadcast latitude array to all processes
    call MPI_BCast(lat, latlen, MPI_DOUBLE_PRECISION, 0, MPI_COMM_FESOM, ierror)  

    !___________________________________________________________________________
    ! read longitude regular coordinates from netcdf file --> lat
    allocate(lon(lonlen))
    if (mype==0) then
        status=nf90_inq_varid(ncid, 'lon', varid)
        if (status == NF90_ENOTVAR) then
            print *, achar(27)//'[33m'
            write(*,*) '____________________________________________________________________'
            write(*,*) ' ERROR: Variable lon not found in file: '
            write(*,*) '        ├> file:', file
            write(*,*) '        ├> status:', status
            write(*,*) '        ├> ',nf90_strerror(status)
            write(*,*) '        └> check: Dim and CoordVariable name must be lon'
            write(*,*) '____________________________________________________________________'
            print *, achar(27)//'[0m'
            write(*,*)
            found_error = .True.
        end if
        status=nf90_get_var(ncid, varid, lon, start=(/1/), count=(/lonlen/))
    end if
    ! broadcast found_error variable to all other processes --> create kill signal
    ! for all processes
    call MPI_BCast(found_error, 1, MPI_LOGICAL, 0, MPI_COMM_FESOM, ierror)
    if (found_error) call par_ex(partit%MPI_COMM_FESOM, partit%mype)
    
    ! broadcast longitude array to all processes
    call MPI_BCast(lon, lonlen, MPI_DOUBLE_PRECISION, 0, MPI_COMM_FESOM, ierror)  

    ! make sure range 0. - 360.
    do n=1,lonlen
        if (lon(n)<0.0_WP) then
            lon(n)=lon(n)+360._WP
        end if
    end do

!     !___________________________________________________________________________
!     ! read z coordiantes 
!     allocate(z(zlen))
!     if (mype==0) then
!         status=nf90_inq_varid(ncid, zvname, varid)
!         if (status == NF90_ENOTVAR) then
!             print *, achar(27)//'[33m'
!             write(*,*) '____________________________________________________________________'
!             write(*,*) ' ERROR: Variable ', zvname ,' not found in file: '
!             write(*,*) '        ├> file:', file
!             write(*,*) '        ├> status:', status
!             write(*,*) '        ├> ',nf90_strerror(status)
!             write(*,*) '        └> check: Dim and CoordVariable name must'
!             write(*,*) '____________________________________________________________________'
!             print *, achar(27)//'[0m'
!             write(*,*)
!             found_error = .True.
!         end if
!         status=nf90_get_var(ncid, varid, z, start=(/1/), count=(/zlen/))
!     end if
!     ! broadcast found_error variable to all other processes --> create kill signal
!     ! for all processes
!     call MPI_BCast(found_error, 1, MPI_LOGICAL, 0, MPI_COMM_FESOM, ierror)
!     if (found_error) call par_ex(partit%MPI_COMM_FESOM, partit%mype)
!     
!     ! broadcast longitude array to all processes
!     call MPI_BCast(z, zlen, MPI_DOUBLE_PRECISION, 0, MPI_COMM_FESOM, ierror)  
    
    !___________________________________________________________________________
    ! read the 3d data --> check if variable name does exist
    allocate(ncdata(lonlen, latlen, zlen))
    ncdata = 0.0_WP
    if (mype==0) then
        ! data
        status=nf90_inq_varid(ncid, vname, varid)
        if (status == NF90_ENOTVAR) then
            print *, achar(27)//'[33m'
            write(*,*) '____________________________________________________________________'
            write(*,*) ' ERROR: Variable ', vname, ' not found in file: '
            write(*,*) '        ├> file:', file
            write(*,*) '        ├> status:', status
            write(*,*) '        ├> ',nf90_strerror(status)
            write(*,*) '        └> check: namelist.*   '
            write(*,*) '____________________________________________________________________'
            print *, achar(27)//'[0m'
            write(*,*)
            found_error = .True.
        end if
        istart= (/     1,      1,    1/)
        icount= (/lonlen, latlen, zlen/)
        status=nf90_get_var(ncid, varid, ncdata, start=istart, count=icount)
        
        ! missing value
        status= nf90_get_att(ncid, varid, 'missing_value', missvalue)
        
        ! close file
        status=nf90_close(ncid)
    end if
    ! broadcast found_error variable to all other processes --> create kill signal
    ! for all processes
    call MPI_BCast(found_error, 1, MPI_LOGICAL, 0, MPI_COMM_FESOM, ierror)
    if (found_error) call par_ex(partit%MPI_COMM_FESOM, partit%mype)
    
    ! broadcast 2d data array to all processes
    call MPI_BCast(ncdata   , lonlen*latlen*zlen, MPI_DOUBLE_PRECISION, 0, MPI_COMM_FESOM, ierror)
    call MPI_BCast(missvalue,                  1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_FESOM, ierror)
    
    !___________________________________________________________________________
    ! size of input array should be 3d (nz, myDim_elem2D) or (nz, myDim_nod2D)
    found_error = .False.
    nz = size(model_3Darray,1)
    if (nz /= zlen) then 
        if (mype==0) then
            print *, achar(27)//'[33m'
            write(*,*) '____________________________________________________________________'
            write(*,*) ' ERROR: length of zdimension in file does not agree with len of zdim in model'
            write(*,*) '____________________________________________________________________'
            print *, achar(27)//'[0m'
            write(*,*)
            found_error = .True.
        end if
    end if 
    ! broadcast found_error variable to all other processes --> create kill signal
    ! for all processes
    call MPI_BCast(found_error, 1, MPI_LOGICAL, 0, MPI_COMM_FESOM, ierror)
    if (found_error) call par_ex(partit%MPI_COMM_FESOM, partit%mype)
    
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
    do k=1,zlen
        call interp_2d_field_v2(lonlen, latlen, &
                                lon, lat, ncdata(:,:,k), missvalue, &
                                num, temp_x, temp_y, model_3Darray(k,1:num), &
                                partit)
    end do
    deallocate(temp_y, temp_x, ncdata, lon, lat)
    
end subroutine read_other_NetCDF_3d








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
    use g_interp
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
  use g_interp
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
