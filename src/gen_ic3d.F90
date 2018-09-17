
MODULE g_ic3d
   !!===========================================================================
   !! Initial conditions for tracers:
   !!===========================================================================
   !! History: 0.1 ! 03/2016 I. Kuznetsov
   
   !! Description:
   !!   read and interpolate initial conditions for tracers on model grid,
   !!       
   !! public: 
   !!   do_ic3d   -- provides an initial 3D boundary conditions
   !!
   USE o_ARRAYS
   USE o_MESH
   USE o_PARAM
   USE g_PARSUP
   USE g_comm_auto
   USE g_support
   USE g_config, only: dummy, ClimateDataPath
   
   IMPLICIT NONE

   include 'netcdf.inc'

   public  do_ic3d, &                                   ! read and apply 3D initial conditions
           n_ic3d, idlist, filelist, varlist, nam_ic3d  ! to be read from the namelist
      
   private

! namelists
   integer, save  :: nm_ic_unit     = 103       ! unit to open namelist file
!============== namelistatmdata variables ================
   integer, parameter                           :: ic_max=10
   integer,        save                         :: n_ic3d
   integer,        save,  dimension(ic_max)     :: idlist
   character(256), save,  dimension(ic_max)     :: filelist
   character(50),  save,  dimension(ic_max)     :: varlist

   namelist / nam_ic3d / n_ic3d, idlist, filelist, varlist 

   character(256), save                         :: filename
   character(50),  save                         :: varname
   integer,        save                         :: current_tracer
   integer,        save                         :: warn       ! warning switch node/element coordinate out of forcing bounds

   ! ========== interpolation coeficients
   integer,  allocatable, save, dimension(:)     :: bilin_indx_i ! indexs i for interpolation
   integer,  allocatable, save, dimension(:)     :: bilin_indx_j ! indexs j for interpolation
!============== NETCDF ==========================================

  ! arrays of time, lon and lat in INfiles
   real(wp), allocatable, save, dimension(:)  :: nc_lon
   real(wp), allocatable, save, dimension(:)  :: nc_lat
   real(wp), allocatable, save, dimension(:)  :: nc_depth! depth (z) in netcdf file [m] , assume we have z-coordinate system

  ! lenght of arrays in INfiles 
   integer,save              :: nc_Nlon
   integer,save              :: nc_Nlat
   integer,save              :: nc_Ndepth   

!============== NETCDF ==========================================   
CONTAINS
   SUBROUTINE nc_readGrid
   ! Read time array and grid from nc file
      IMPLICIT NONE

      integer              :: iost !I/O status     
      integer              :: ncid      ! netcdf file id
      integer              :: i
      ! ID dimensions and variables:
      integer              :: id_lon
      integer              :: id_lat
      integer              :: id_lond
      integer              :: id_latd
      integer              :: id_depth
      integer              :: id_depthd      
      integer              :: nf_start(4)
      integer              :: nf_edges(4)         
      integer              :: ierror              ! return error code
    
      !open file
      if (mype==0) then
         iost = nf_open(trim(filename),NF_NOWRITE,ncid)
      end if
      call MPI_BCast(iost, 1, MPI_INTEGER, 0, MPI_COMM_FESOM, ierror)
      call check_nferr(iost,filename)

      ! get dimensions
      if (mype==0) then
         iost = nf_inq_dimid(ncid, "lat", id_latd)
      end if
      call MPI_BCast(iost, 1, MPI_INTEGER, 0, MPI_COMM_FESOM, ierror)
      call check_nferr(iost,filename)  
      if (mype==0) then 
         iost = nf_inq_dimid(ncid, "lon", id_lond)
      end if
      call MPI_BCast(iost, 1, MPI_INTEGER, 0, MPI_COMM_FESOM, ierror)
      call check_nferr(iost,filename) 
      if (mype==0) then   
         iost = nf_inq_dimid(ncid, "depth", id_depthd)
      end if
      call MPI_BCast(iost, 1, MPI_INTEGER, 0, MPI_COMM_FESOM, ierror)
      call check_nferr(iost,filename)  

      ! get variable id
      if (mype==0) then
         iost = nf_inq_varid(ncid, "lon", id_lon)
      end if
      call MPI_BCast(iost, 1, MPI_INTEGER, 0, MPI_COMM_FESOM, ierror)
      call check_nferr(iost,filename)
      if (mype==0) then
         iost = nf_inq_varid(ncid, "lat", id_lat)
      end if
      call MPI_BCast(iost, 1, MPI_INTEGER, 0, MPI_COMM_FESOM, ierror)
      call check_nferr(iost,filename)  
      if (mype==0) then   
         iost = nf_inq_varid(ncid, "depth", id_depth)
      end if
      call MPI_BCast(iost, 1, MPI_INTEGER, 0, MPI_COMM_FESOM, ierror)
      call check_nferr(iost,filename)   
      
      ! get dimensions size
      if (mype==0) then
         iost = nf_inq_dimlen(ncid, id_latd, nc_Nlat)
      end if
      call MPI_BCast(iost, 1, MPI_INTEGER, 0, MPI_COMM_FESOM, ierror)
      call check_nferr(iost,filename)   
      if (mype==0) then      
         iost = nf_inq_dimlen(ncid, id_lond, nc_Nlon)
      end if
      call MPI_BCast(iost, 1, MPI_INTEGER, 0, MPI_COMM_FESOM, ierror)
      call check_nferr(iost,filename)   
      if (mype==0) then      
         iost = nf_inq_dimlen(ncid, id_depthd, nc_Ndepth)
      end if
      call MPI_BCast(iost, 1, MPI_INTEGER, 0, MPI_COMM_FESOM, ierror)
      call check_nferr(iost,filename) 

      call MPI_BCast(nc_Nlon,   1, MPI_INTEGER, 0, MPI_COMM_FESOM, ierror)
      call MPI_BCast(nc_Nlat,   1, MPI_INTEGER, 0, MPI_COMM_FESOM, ierror)
      call MPI_BCast(nc_Ndepth, 1, MPI_INTEGER, 0, MPI_COMM_FESOM, ierror)

      ALLOCATE( nc_lon(nc_Nlon), nc_lat(nc_Nlat),&
                &       nc_depth(nc_Ndepth))

   !read variables from file
   ! coordinates
      if (mype==0) then
         nf_start(1)=1
         nf_edges(1)=nc_Nlat
         iost = nf_get_vara_double(ncid, id_lat, nf_start, nf_edges, nc_lat)
      end if
      call MPI_BCast(iost, 1, MPI_INTEGER, 0, MPI_COMM_FESOM, ierror)
      call check_nferr(iost,filename)
      if (mype==0) then
         nf_start(1)=1
         nf_edges(1)=nc_Nlon
         iost = nf_get_vara_double(ncid, id_lon, nf_start, nf_edges, nc_lon)
      end if
      call MPI_BCast(iost, 1, MPI_INTEGER, 0, MPI_COMM_FESOM, ierror)      
      call check_nferr(iost,filename)
   ! depth
      if (mype==0) then
         nf_start(1)=1
         nf_edges(1)=nc_Ndepth
         iost = nf_get_vara_double(ncid, id_depth, nf_start, nf_edges,nc_depth)
      end if
      call MPI_BCast(iost, 1, MPI_INTEGER, 0, MPI_COMM_FESOM, ierror)      
      call check_nferr(iost,filename)

      call MPI_BCast(nc_lon,   nc_Nlon,   MPI_DOUBLE_PRECISION, 0, MPI_COMM_FESOM, ierror)
      call MPI_BCast(nc_lat,   nc_Nlat,   MPI_DOUBLE_PRECISION, 0, MPI_COMM_FESOM, ierror)
      call MPI_BCast(nc_depth, nc_Ndepth, MPI_DOUBLE_PRECISION, 0, MPI_COMM_FESOM, ierror)

      if (mype==0) then
         iost = nf_close(ncid)
      end if
      call MPI_BCast(iost, 1, MPI_INTEGER, 0, MPI_COMM_FESOM, ierror)
      call check_nferr(iost,filename)
   END SUBROUTINE nc_readGrid

   
   SUBROUTINE nc_ic3d_ini
      !!---------------------------------------------------------------------
      !! ** Purpose : inizialization of ocean forcing from NETCDF file
      !!----------------------------------------------------------------------
      IMPLICIT NONE
   
      integer            :: i
      integer            :: elnodes(3)
      real(wp)           :: x, y       ! coordinates of elements
      
      warn = 0

      if (mype==0) then
         write(*,*) 'reading input tracer file for tracer ID= ', tracer_ID(current_tracer)
         write(*,*) 'input file: ', trim(filename)
         write(*,*) 'variable  : ', trim(varname)
      end if
      call nc_readGrid

      ! prepare nearest coordinates in INfile , save to bilin_indx_i/j
      do i = 1, myDim_nod2d
!         ! get coordinates of elements
         x  = geo_coord_nod2D(1,i)/rad
         y  = geo_coord_nod2D(2,i)/rad
         if (x<0.) x=x+360.
         ! find nearest
         if ( x < nc_lon(nc_Nlon) .and. x >= nc_lon(1) ) then      
            call binarysearch(nc_Nlon, nc_lon, x, bilin_indx_i(i))
         else ! NO extrapolation in space
            if ( x < nc_lon(1) ) then               
               bilin_indx_i(i)=-1         
            else
               bilin_indx_i(i)=0
            end if
         end if
         if ( y < nc_lat(nc_Nlat) .and. y >= nc_lat(1) ) then      
            call binarysearch(nc_Nlat, nc_lat, y, bilin_indx_j(i))
         else ! NO extrapolation in space
            if ( y < nc_lat(1) ) then               
               bilin_indx_j(i)=-1         
            else
               bilin_indx_j(i)=0
            end if
         end if
      end do
                         
   END SUBROUTINE nc_ic3d_ini

   SUBROUTINE getcoeffld
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE getcoeffld ***
      !!              
      !! ** Purpose : read fields from files, inetrpolate on model mesh and prepare interpolation coeffients 
      !! ** Method  : 
      !! ** Action  : 
      !!----------------------------------------------------------------------
      IMPLICIT NONE      
      
      integer              :: iost !I/O status     
      integer              :: ncid      ! netcdf file id
      ! ID dimensions and variables:
      integer              :: id_data
      integer              :: nf_start(4)
      integer              :: nf_edges(4)         
      integer              :: fld_idx, i,j,ii, ip1, jp1, k
      integer              :: d_indx, d_indx_p1  ! index of neares      
      real(wp)             :: cf_a, cf_b, delta_d, nl1
      real(wp)             :: denom, x1, x2, y1, y2, x, y, d1,d2     
      
      real(wp), allocatable, dimension(:,:,:)  :: ncdata
      real(wp), allocatable, dimension(:)      :: data1d      
      integer              :: elnodes(3)
      integer              :: ierror              ! return error code

      ALLOCATE(ncdata(nc_Nlon,nc_Nlat,nc_Ndepth), data1d(nc_Ndepth))
      tr_arr(:,:,current_tracer)=dummy
      !open NETCDF file on 0 core     
      if (mype==0) then
         iost = nf_open(filename,NF_NOWRITE,ncid)
      end if
      call MPI_BCast(iost, 1, MPI_INTEGER, 0, MPI_COMM_FESOM, ierror)
      call check_nferr(iost,filename)
      ! get variable id
      if (mype==0) then
         iost = nf_inq_varid(ncid, varname, id_data)
      end if
      call MPI_BCast(iost, 1, MPI_INTEGER, 0, MPI_COMM_FESOM, ierror)
      call check_nferr(iost,filename)   
      !read data from file
      if (mype==0) then
         nf_start(1)=1
         nf_edges(1)=nc_Nlon
         nf_start(2)=1
         nf_edges(2)=nc_Nlat
         nf_start(3)=1
         nf_edges(3)=nc_Ndepth         
         iost = nf_get_vara_double(ncid, id_data, nf_start, nf_edges, ncdata)
      end if
      call MPI_BCast(iost, 1, MPI_INTEGER, 0, MPI_COMM_FESOM, ierror)
      call check_nferr(iost,filename)
      call MPI_BCast(ncdata, nc_Nlon*nc_Nlat*nc_Ndepth, MPI_DOUBLE_PRECISION, 0, MPI_COMM_FESOM, ierror)

      ! bilinear space interpolation,  
      ! data is assumed to be sampled on a regular grid
      do ii = 1, myDim_nod2d
         nl1 = nlevels_nod2D(ii)-1
         i = bilin_indx_i(ii)
         j = bilin_indx_j(ii)
         ip1 = i + 1   
         jp1 = j + 1
         x  = geo_coord_nod2D(1,ii)/rad
         y  = geo_coord_nod2D(2,ii)/rad
         if (x<0.) x=x+360.
         if ( min(i,j)>0 ) then
         if (any(ncdata(i:ip1,j:jp1,1) > dummy*0.99)) cycle
         x1 = nc_lon(i)
         x2 = nc_lon(ip1)
         y1 = nc_lat(j)
         y2 = nc_lat(jp1)
                  
         ! if point inside forcing domain
         denom = (x2 - x1)*(y2 - y1)
         data1d(:) = ( ncdata(i,j,:)   * (x2-x)*(y2-y)   + ncdata(ip1,j,:)    * (x-x1)*(y2-y) + &
                       ncdata(i,jp1,:) * (x2-x)*(y-y1)   + ncdata(ip1, jp1, :) * (x-x1)*(y-y1)     ) / denom          
         do k= 1, nl1
            call binarysearch(nc_Ndepth,nc_depth,-Z_3d_n(k,ii),d_indx)
            if ( d_indx < nc_Ndepth ) then
               d_indx_p1 = d_indx+1
               delta_d = nc_depth(d_indx+1)-nc_depth(d_indx)
               ! values from OB data for nearest depth           
               d1 = data1d(d_indx)
               d2 = data1d(d_indx_p1)
               if ((d1<0.99*dummy) .and. (d2<0.99*dummy)) then
               ! line a*z+b coefficients calculation
               cf_a  = (d2 - d1)/ delta_d
               ! value of interpolated OB data on Z from model
               cf_b  = d1 - cf_a * nc_depth(d_indx)
               tr_arr(k,ii,current_tracer) = -cf_a * Z_3d_n(k,ii) + cf_b
               end if
            end if
         enddo
         end if
      end do !ii
      if (mype==0) then
         iost = nf_close(ncid)
      end if
      call MPI_BCast(iost, 1, MPI_INTEGER, 0, MPI_COMM_FESOM, ierror)
      DEALLOCATE( ncdata, data1d )
   END SUBROUTINE getcoeffld  
   
   SUBROUTINE do_ic3d
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE do_ic3d ***
      !!              
      !! ** Purpose : read 3D initial conditions for tracers from netcdf and interpolate on model grid
      !!----------------------------------------------------------------------
      IMPLICIT NONE
      integer                       :: n
      if (mype==0) write(*,*) "Start: Initial conditions  for tracers"

      ALLOCATE(bilin_indx_i(myDim_nod2d+eDim_nod2D), bilin_indx_j(myDim_nod2d+eDim_nod2D))
      DO n=1, n_ic3d
      filename=trim(ClimateDataPath)//trim(filelist(n))
      varname =trim(varlist(n))
      DO current_tracer=1, num_tracers
         if (tracer_ID(current_tracer)==idlist(n)) then
            ! read initial conditions for current tracer
            call nc_ic3d_ini
            ! get first coeficients for time inerpolation on model grid for all datas
            call getcoeffld
            call nc_end ! deallocate arrqays associated with netcdf file
            call extrap_nod(tr_arr(:,:,current_tracer))
            exit
         elseif (current_tracer==num_tracers) then
            if (mype==0) write(*,*) "idlist contains tracer which is not listed in tracer_id!"
            if (mype==0) write(*,*) "check your namelists!"
            call par_ex
            stop
         end if
      END DO
      END DO
      DEALLOCATE(bilin_indx_i, bilin_indx_j)
      call insitu2pot
      if (mype==0) write(*,*) "DONE:  Initial conditions for tracers"
   
   END SUBROUTINE do_ic3d
   
   SUBROUTINE err_call(iost,fname)
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE  err_call ***
      !!----------------------------------------------------------------------
   IMPLICIT NONE
      integer, intent(in)            :: iost
      character(len=256), intent(in) :: fname
      write(*,*) 'ERROR: I/O status=',iost,' file= ',fname
      call par_ex
      stop
   END SUBROUTINE err_call
   

   SUBROUTINE nc_end

    IMPLICIT NONE
 
    DEALLOCATE(nc_lon, nc_lat, nc_depth)

   END SUBROUTINE nc_end

   SUBROUTINE check_nferr(iost,fname)
   IMPLICIT NONE 
      character(len=256), intent(in) :: fname
      integer, intent(in) :: iost
      if (iost .ne. NF_NOERR) then
         write(*,*) 'ERROR: I/O status= "',trim(nf_strerror(iost)),'";',iost,' file= ', trim(fname)
         call par_ex 
         stop
      endif
   END SUBROUTINE

   SUBROUTINE binarysearch(length, array, value, ind)!, delta)
      ! Given an array and a value, returns the index of the element that
      ! is closest to, but less than, the given value.
      ! Uses a binary search algorithm.
      ! "delta" is the tolerance used to determine if two values are equal
      ! if ( abs(x1 - x2) <= delta) then
      !    assume x1 = x2
      ! endif
      !org. source from: https://github.com/cfinch/Shocksolution_Examples/blob/master/FORTRAN/BilinearInterpolation/interpolation.f90

      IMPLICIT NONE 
      integer,  intent(in)                    :: length
      real(wp), dimension(length), intent(in) :: array
      real(wp), intent(in)                    :: value
      integer, intent(out)                    :: ind
      integer                                 :: left, middle, right
      real(wp)                                :: d

      d = 1e-9
      left = 1
      right = length
      do
         if (left > right) then
            exit
         endif
         middle = nint((left+right) / 2.0)
         if ( abs(array(middle) - value) <= d) then
            ind = middle
            return
         else if (array(middle) > value) then
            right = middle - 1
         else
            left = middle + 1
         end if
      end do
      ind = right

   END SUBROUTINE binarysearch
END MODULE g_ic3d
