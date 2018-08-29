
MODULE fv_ic
   !!===========================================================================
   !! Initial conditions:
   !!===========================================================================
   !! History: 0.1 ! 03/2016 I. Kuznetsov
   
   !! Description:
   !!   read and interpolate initial conditions (T and S) on model grid,
   !!     or use constants from namelist 
   !!   
   !!           
   !!       
   !! public: 
   !!   ic_do   -- provide a sbc (surface boundary conditions) each time step
   !!
   USE o_ARRAYS
   USE o_MESH
   USE o_PARAM
   USE g_PARSUP
   USE g_comm_auto
   USE g_support
   USE g_config, only: dummy
   
   IMPLICIT NONE

   include 'netcdf.inc'

   public  ic_do   ! read and apply IC
 
   private
!

   
   ! namelists
   integer, save  :: nm_ic_unit     = 103       ! unit to open namelist file
  !============== namelistatmdata variables ================
   integer, save  :: nm_ic        = -1        ! data  1= constant, 2=from file, -1=module off
   character(len=256), save   :: nm_ict_file = 'ic' ! name of file with initial conditions for temperature   
   character(len=256), save   :: nm_ics_file = 'ic' ! name of file with initial conditions for salinity             
   
   character(len=34), save   :: nm_temp_var = 'temp' ! name of variable in file with temperature
   character(len=34), save   :: nm_salt_var = 'salt' ! name of variable in file with salinity

   
   real(wp), save :: nm_ict = 15._wp    ! constant temperature 
   real(wp), save :: nm_ics = 32._wp    ! constant salinity 
   
   integer,save            :: warn       ! warning switch node/element coordinate out of forcing bounds

   ! ========== interpolation coeficients
   integer,  allocatable, save, dimension(:)     :: bilin_indx_i ! indexs i for interpolation
   integer,  allocatable, save, dimension(:)     :: bilin_indx_j ! indexs j for interpolation

   real(wp), allocatable, save, dimension(:,:,:) :: icdata ! IC data   
!============== NETCDF ==========================================
   integer, parameter :: i_totfl = 2 ! total number of fluxes
   integer, parameter :: i_temp = 1 ! index of tempearature [degC]
   integer, parameter :: i_salt = 2 ! index of salinity
   
   
   type, public ::   flfi_type    !flux file informations
      character(len = 256) :: file_name ! file name   
      character(len = 34)  :: var_name  ! variable name in the NetCDF file
   end type flfi_type     
 
  type(flfi_type),save, dimension(i_totfl) :: ic_flfi  !array for information about flux files

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
   SUBROUTINE nc_readTimeGrid(flf)
   ! Read time array and grid from nc file
      IMPLICIT NONE

      type(flfi_type),intent(in) :: flf

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
      
!      integer              :: nf_dims(4) ! dimensions (temporal)
      integer              :: nf_start(4)
      integer              :: nf_edges(4)         
      integer              :: sbc_alloc                   !: allocation status
      

    
      !open file
      iost = nf_open(flf%file_name,NF_NOWRITE,ncid)
      call check_nferr(iost,flf%file_name)

      ! get dimensions
      iost = nf_inq_dimid(ncid, "lat", id_latd)
      call check_nferr(iost,flf%file_name)   
      iost = nf_inq_dimid(ncid, "lon", id_lond)
      call check_nferr(iost,flf%file_name)   
      iost = nf_inq_dimid(ncid, "depth", id_depthd)
      call check_nferr(iost,flf%file_name)  

      ! get variable id
!      iost = nf_inq_varid(ncid, "air", id_data)
!      call check_nferr(iost,flf%file_name)   
      iost = nf_inq_varid(ncid, "lon", id_lon)
      call check_nferr(iost,flf%file_name)   
      iost = nf_inq_varid(ncid, "lat", id_lat)
      call check_nferr(iost,flf%file_name)    
      iost = nf_inq_varid(ncid, "depth", id_depth)
      call check_nferr(iost,flf%file_name)   
      
      !  get dimensions size
      iost = nf_inq_dimlen(ncid, id_latd, nc_Nlat)
      call check_nferr(iost,flf%file_name)   
      iost = nf_inq_dimlen(ncid, id_lond, nc_Nlon)
      call check_nferr(iost,flf%file_name)   
      iost = nf_inq_dimlen(ncid, id_depthd, nc_Ndepth)
      call check_nferr(iost,flf%file_name) 

      ALLOCATE( nc_lon(nc_Nlon), nc_lat(nc_Nlat),&
                &       nc_depth(nc_Ndepth), &

                &      STAT=sbc_alloc )  
      if( sbc_alloc /= 0 )   STOP 'ic:nc_readTimeGrid: failed to allocate arrays'   
         

   !read variables from file
   ! coordinates
      nf_start(1)=1
      nf_edges(1)=nc_Nlat
      iost = nf_get_vara_double(ncid, id_lat, nf_start, nf_edges, nc_lat)
      call check_nferr(iost,flf%file_name)
      nf_start(1)=1
      nf_edges(1)=nc_Nlon
      iost = nf_get_vara_double(ncid, id_lon, nf_start, nf_edges, nc_lon)
      call check_nferr(iost,flf%file_name)
   ! depth
      nf_start(1)=1
      nf_edges(1)=nc_Ndepth
      iost = nf_get_vara_double(ncid, id_depth, nf_start, nf_edges,nc_depth)
      call check_nferr(iost,flf%file_name)


      iost = nf_close(ncid)
      call check_nferr(iost,flf%file_name)
   END SUBROUTINE nc_readTimeGrid

   
   SUBROUTINE nc_ic_ini
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE ascii_ini ***
      !!              
      !! ** Purpose : inizialization of ocean forcing from NETCDF file
      !!----------------------------------------------------------------------

      IMPLICIT NONE
   
      integer            :: i
      integer            :: sbc_alloc

      integer            :: elnodes(4) !4 nodes from one element
      integer            :: numnodes   ! number of nodes in elem (3 for triangle, 4 for ... )
      real(wp)           :: x, y       ! coordinates of elements
      
      warn = 0
      ALLOCATE( icdata(i_totfl,nl-1, myDim_nod2d+eDim_nod2D), &
                   &      STAT=sbc_alloc )  
      if( sbc_alloc /= 0 )   STOP 'nc_ic_ini: failed to allocate arrays' 
      ALLOCATE( bilin_indx_i(myDim_nod2d+eDim_nod2D),bilin_indx_j(myDim_nod2d+eDim_nod2D), &
                   &      STAT=sbc_alloc )  
                   
      if( sbc_alloc /= 0 )   STOP 'nc_ic_ini: failed to allocate arrays' 
      
      ! we assume that all NetCDF files have identical grid and time variable
      write(ic_flfi(i_temp)%file_name,*) trim(nm_ict_file),'.nc'
      ic_flfi(i_temp)%file_name=ADJUSTL(trim(ic_flfi(i_temp)%file_name))
      ic_flfi(i_temp)%var_name=ADJUSTL(trim(nm_temp_var))

      write(ic_flfi(i_salt)%file_name,*) trim(nm_ics_file),'.nc'
      ic_flfi(i_salt)%file_name=ADJUSTL(trim(ic_flfi(i_salt)%file_name))
      ic_flfi(i_salt)%var_name=ADJUSTL(trim(nm_salt_var))
      
      call nc_readTimeGrid(ic_flfi(i_temp))

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
         if (warn == 0) then
            if (bilin_indx_i(i) < 1 .or. bilin_indx_j(i) < 1) then
               WRITE(*,*) '     WARNING:  node/element coordinate out of initial conditions bounds,'
               WRITE(*,*) '        nearest value will be used'
               warn = 1
            end if   
         end if

         
         
         
      end do
                         
   END SUBROUTINE nc_ic_ini


   
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
!      integer              :: zero_year,yyyy,mm,dd
!      character(len = 256) :: att_string ! attribute              
      integer              :: fld_idx, i,j,ii, ip1, jp1, k
      integer              :: sbc_alloc
  
      integer              :: d_indx, d_indx_p1  ! index of neares      
      real(wp)             :: cf_a, cf_b, delta_d, nl1
      real(wp)             :: denom, x1, x2, y1, y2, x, y, d1,d2     
      
      real(wp), allocatable, dimension(:,:,:)  :: sbcdata1
      real(wp), allocatable, dimension(:)      :: data1
      
      integer              :: elnodes(4) !4 nodes from one element
      integer              :: numnodes   ! number of nodes in elem (3 for triangle, 4 for ... )

      ALLOCATE( sbcdata1(nc_Nlon,nc_Nlat,nc_Ndepth), &
                data1(nc_Ndepth), &                
                &      STAT=sbc_alloc )
      if( sbc_alloc /= 0 )   STOP 'getcoeffld: failed to allocate arrays'   
      icdata=dummy
      do fld_idx = 1, i_totfl
         !open file sbc_flfi      
         iost = nf_open(ic_flfi(fld_idx)%file_name,NF_NOWRITE,ncid)
         call check_nferr(iost,ic_flfi(fld_idx)%file_name)
         ! get variable id
         iost = nf_inq_varid(ncid, ic_flfi(fld_idx)%var_name, id_data)
         call check_nferr(iost,ic_flfi(fld_idx)%file_name)   
         !read data from file
         nf_start(1)=1
         nf_edges(1)=nc_Nlon
         nf_start(2)=1
         nf_edges(2)=nc_Nlat
         nf_start(3)=1
         nf_edges(3)=nc_Ndepth         
         iost = nf_get_vara_double(ncid, id_data, nf_start, nf_edges, sbcdata1)
         call check_nferr(iost,ic_flfi(fld_idx)%file_name)


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
            if (any(sbcdata1(i:ip1,j:jp1,1) > dummy*0.99)) cycle
            x1 = nc_lon(i)
            x2 = nc_lon(ip1)
            y1 = nc_lat(j)
            y2 = nc_lat(jp1)
                  
            ! if point inside forcing domain
            denom = (x2 - x1)*(y2 - y1)
            data1(:) = ( sbcdata1(i,j,:)   * (x2-x)*(y2-y)   + sbcdata1(ip1,j,:)    * (x-x1)*(y2-y) + &
                          sbcdata1(i,jp1,:) * (x2-x)*(y-y1)   + sbcdata1(ip1, jp1, :) * (x-x1)*(y-y1)     ) / denom          
            do k= 1, nl1-1
               call binarysearch(nc_Ndepth,nc_depth,-Z_3d_n(k,ii),d_indx)
               if ( d_indx < nc_Ndepth ) then
                  d_indx_p1 = d_indx+1
                  delta_d = nc_depth(d_indx+1)-nc_depth(d_indx)
                  ! values from OB data for nearest depth           
                  d1 = data1(d_indx)
                  d2 = data1(d_indx_p1)
                  if ((d1<dummy) .and. (d2<dummy)) then
                  ! line a*z+b coefficients calculation
                  cf_a  = (d2 - d1)/ delta_d
                  ! value of interpolated OB data on Z from model
                  cf_b  = d1 - cf_a * nc_depth(d_indx)
                  icdata(fld_idx,k,ii) = cf_a * Z_3d_n(k,ii) + cf_b
                  end if
               end if
            enddo
            else          
               icdata(fld_idx,1:nl1,ii) = dummy
            end if
         end do !ii
         iost = nf_close(ncid)
         call check_nferr(iost,ic_flfi(fld_idx)%file_name)
         
      end do !fld_idx      

      DEALLOCATE( sbcdata1, data1 )
                

   END SUBROUTINE getcoeffld  

   
   SUBROUTINE ic_ini
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE sbc_ini ***
      !!              
      !! ** Purpose : inizialization of initial conditions
      !! ** Method  : 
      !! ** Action  : 
      !!----------------------------------------------------------------------
      IMPLICIT NONE
     
      integer            :: iost  ! I/O status
      integer            :: sbc_alloc                   !: allocation status
      

      
      namelist/nam_ic/ nm_ic, nm_ict_file, nm_ics_file, nm_temp_var, nm_salt_var, &
                        nm_ict, nm_ics

      
      ! OPEN and read namelist for SBC
      open( unit=nm_ic_unit, file='namelist_bc.nml', form='formatted', access='sequential', status='old', iostat=iost )
      if( iost == 0 ) then
         WRITE(*,*) '     file   : ', 'namelist_bc.nml',' open ok'
      else
         WRITE(*,*) 'ERROR: --> bad opening file   : ', 'namelist_bc.nml',' ; iostat=',iost     
         STOP 'ERROR: --> ic_ini'
      endif
      READ( nm_ic_unit, nml=nam_ic, iostat=iost )
      close( nm_ic_unit)
      
   END SUBROUTINE ic_ini  
   
   SUBROUTINE ic_do
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE ic_do ***
      !!              
      !! ** Purpose : read ICfrom netcdf, interpolate on model grid
      !! ** Method  : 
      !! ** Action  : 
      !!----------------------------------------------------------------------
      IMPLICIT NONE
      write(*,*) "Start: Initial conditions."

      call ic_ini ! read namelist
      
      if( nm_ic == -1 ) then
         ! IF module not in use
         write(*,*) "Initial conditions module is OFF."      
         write(*,*) "DONE:  Initial conditions."      
         return
      endif      
      if( nm_ic == 1 ) then
         tr_arr(:,:,1) = nm_ict;
         tr_arr(:,:,2) = nm_ics;
         write(*,*) "Initial conditions are constants."      
      endif     
      if( nm_ic == 2 ) then
         call nc_ic_ini     
         ! get first coeficients for time inerpolation on model grid for all datas
         call getcoeffld
         tr_arr(:,:,1) = icdata(i_temp,:,:)
         tr_arr(:,:,2) = icdata(i_salt,:,:)
         call extrap_nod(tr_arr(:,:,1))
         call extrap_nod(tr_arr(:,:,2))
         call ic_end ! deallocate 
      endif        
      write(*,*) "DONE:  Initial conditions."      
   
   END SUBROUTINE ic_do   

   
   SUBROUTINE err_call(iost,fname)
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE  err_call ***
      !!              
      !! ** Purpose : call Error
      !! ** Method  : 
      !! ** Action  : 
      !!----------------------------------------------------------------------
   IMPLICIT NONE
      integer, intent(in)            :: iost
      character(len=256), intent(in) :: fname
      write(*,*) 'ERROR: I/O status=',iost,' file= ',fname
      STOP 'ERROR:  stop'


   END SUBROUTINE err_call
   
  
   SUBROUTINE ic_end

      IMPLICIT NONE
 
         DEALLOCATE( nc_lon, nc_lat, nc_depth,     &
                  &  bilin_indx_i, bilin_indx_j,  &
                  &  icdata )


      
   END SUBROUTINE ic_end
   
   SUBROUTINE check_nferr(iost,fname)
   IMPLICIT NONE 
      character(len=256), intent(in) :: fname
      integer, intent(in) :: iost
   
      if (iost .ne. NF_NOERR) then
         write(*,*) 'ERROR: I/O status= "',trim(nf_strerror(iost)),'";',iost,' file= ',fname
         STOP 'ERROR: stop'
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
      integer,  intent(in) :: length
      real(wp), dimension(length), intent(in) :: array
      real(wp), intent(in) :: value
   !   real, intent(in), optional :: delta

   !   integer :: binarysearch
      integer, intent(out) :: ind

      integer :: left, middle, right
      real(wp):: d

   !   if (present(delta) .eqv. .true.) then
   !      d = delta
   !   else
      d = 1e-9
   !   endif
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
   
END MODULE fv_ic
