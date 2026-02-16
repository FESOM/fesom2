
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
   USE MOD_MESH
   USE MOD_PARTIT
   USE MOD_PARSUP
   USE MOD_TRACER
   USE o_PARAM
   USE g_comm_auto
   USE g_support
   USE g_config, only: dummy, ClimateDataPath, use_cavity
   
   IMPLICIT NONE

   include 'netcdf.inc'

   public  do_ic3d, &                                       ! read and apply 3D initial conditions
           n_ic3d, idlist, filelist, varlist, tracer_init3d, & ! to be read from the namelist
           t_insitu
   private

! namelists
   integer, save  :: nm_ic_unit     = 103       ! unit to open namelist file
!============== namelistatmdata variables ================
   integer, parameter                           :: ic_max=10
   logical                                      :: ic_cyclic=.true.
   logical                                      :: t_insitu =.true.
   integer,        save                         :: n_ic3d
   integer,        save,  dimension(ic_max)     :: idlist
   character(MAX_PATH), save, dimension(ic_max) :: filelist
   character(50),  save,  dimension(ic_max)     :: varlist

   namelist / tracer_init3d / n_ic3d, idlist, filelist, varlist, t_insitu

   character(MAX_PATH), save                    :: filename
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
   SUBROUTINE nc_readGrid(partit)
   ! Read time array and grid from nc file
      IMPLICIT NONE
      type(t_partit), intent(inout) :: partit 
      integer                       :: iost !I/O status     
      integer                       :: ncid      ! netcdf file id
      integer                       :: i
      ! ID dimensions and variables:
      integer                       :: id_lon
      integer                       :: id_lat
      integer                       :: id_lond
      integer                       :: id_latd
      integer                       :: id_depth
      integer                       :: id_depthd      
      integer                       :: nf_start(4)
      integer                       :: nf_edges(4)         
      integer                       :: ierror              ! return error code
    
      !open file
      if (partit%mype==0) then
         iost = nf_open(trim(filename),NF_NOWRITE,ncid)
      end if
      call MPI_BCast(iost, 1, MPI_INTEGER, 0, partit%MPI_COMM_FESOM, ierror)
      call check_nferr(iost,filename,partit)

      ! get dimensions
      if (partit%mype==0) then
         iost = nf_inq_dimid(ncid,    "LAT",      id_latd)
         if (iost .ne. NF_NOERR) then
            iost = nf_inq_dimid(ncid, "lat",      id_latd)
         end if
         if (iost .ne. NF_NOERR) then
            iost = nf_inq_dimid(ncid, "latitude", id_latd)
         end if
      end if
      call MPI_BCast(iost, 1, MPI_INTEGER, 0, partit%MPI_COMM_FESOM, ierror)
      call check_nferr(iost,filename,partit)
      if (partit%mype==0) then 
         iost = nf_inq_dimid(ncid,    "LON",       id_lond)
         if      (iost .ne. NF_NOERR) then
            iost = nf_inq_dimid(ncid, "longitude", id_lond)
         end if
         if (iost .ne. NF_NOERR) then
            iost = nf_inq_dimid(ncid, "lon",       id_lond)
         end if
      end if
      call MPI_BCast(iost, 1, MPI_INTEGER, 0, partit%MPI_COMM_FESOM, ierror)
      call check_nferr(iost,filename,partit) 
      if (partit%mype==0) then   
         iost = nf_inq_dimid(ncid, "depth", id_depthd)
      end if
      call MPI_BCast(iost, 1, MPI_INTEGER, 0, partit%MPI_COMM_FESOM, ierror)
      call check_nferr(iost,filename,partit)  

      ! get variable id
      call MPI_BCast(iost, 1, MPI_INTEGER, 0, partit%MPI_COMM_FESOM, ierror)
      call check_nferr(iost,filename,partit)
      if (partit%mype==0) then
         iost = nf_inq_varid(ncid,    "LAT",      id_lat)
         if     (iost .ne. NF_NOERR) then
            iost = nf_inq_varid(ncid, "lat",      id_lat)
         end if
         if (iost .ne. NF_NOERR) then
            iost = nf_inq_varid(ncid, "latitude", id_lat)
         end if
      end if
      if (partit%mype==0) then
         iost = nf_inq_varid(ncid,    "LON",       id_lon)
         if      (iost .ne. NF_NOERR) then
            iost = nf_inq_varid(ncid, "longitude", id_lon)
         end if
         if (iost .ne. NF_NOERR) then
            iost = nf_inq_varid(ncid, "lon",       id_lon)
         end if
      end if
      call MPI_BCast(iost, 1, MPI_INTEGER, 0, partit%MPI_COMM_FESOM, ierror)
      call check_nferr(iost,filename,partit)  
      if (partit%mype==0) then   
         iost = nf_inq_varid(ncid, "depth", id_depth)
      end if
      call MPI_BCast(iost, 1, MPI_INTEGER, 0, partit%MPI_COMM_FESOM, ierror)
      call check_nferr(iost,filename,partit)   
      
      ! get dimensions size
      if (partit%mype==0) then
         iost = nf_inq_dimlen(ncid, id_latd, nc_Nlat)
      end if
      call MPI_BCast(iost, 1, MPI_INTEGER, 0, partit%MPI_COMM_FESOM, ierror)
      call check_nferr(iost,filename,partit)   
      if (partit%mype==0) then      
         iost = nf_inq_dimlen(ncid, id_lond, nc_Nlon)
      end if
      call MPI_BCast(iost, 1, MPI_INTEGER, 0, partit%MPI_COMM_FESOM, ierror)
      call check_nferr(iost,filename,partit)   
      if (partit%mype==0) then      
         iost = nf_inq_dimlen(ncid, id_depthd, nc_Ndepth)
      end if
      call MPI_BCast(iost, 1, MPI_INTEGER, 0, partit%MPI_COMM_FESOM, ierror)
      call check_nferr(iost,filename,partit) 
      nc_Nlon=nc_Nlon+2 !for the halo in case of periodic boundary
      call MPI_BCast(nc_Nlon,   1, MPI_INTEGER, 0, partit%MPI_COMM_FESOM, ierror)
      call MPI_BCast(nc_Nlat,   1, MPI_INTEGER, 0, partit%MPI_COMM_FESOM, ierror)
      call MPI_BCast(nc_Ndepth, 1, MPI_INTEGER, 0, partit%MPI_COMM_FESOM, ierror)
         

      ALLOCATE( nc_lon(nc_Nlon), nc_lat(nc_Nlat),&
                &       nc_depth(nc_Ndepth))
   !read variables from file
   ! coordinates
      if (partit%mype==0) then
         nf_start(1)=1
         nf_edges(1)=nc_Nlat
         iost = nf_get_vara_double(ncid, id_lat, nf_start, nf_edges, nc_lat)
      end if
      call MPI_BCast(iost, 1, MPI_INTEGER, 0, partit%MPI_COMM_FESOM, ierror)
      call check_nferr(iost,filename,partit)
      if (partit%mype==0) then
         nf_start(1)=1
         nf_edges(1)=nc_Nlon-2
         iost = nf_get_vara_double(ncid, id_lon, nf_start, nf_edges, nc_lon(2:nc_Nlon-1))
         nc_lon(1)        =nc_lon(nc_Nlon-1)
         nc_lon(nc_Nlon)  =nc_lon(2)
      end if
      call MPI_BCast(iost, 1, MPI_INTEGER, 0, partit%MPI_COMM_FESOM, ierror)      
      call check_nferr(iost,filename,partit)
   ! depth
      if (partit%mype==0) then
         nf_start(1)=1
         nf_edges(1)=nc_Ndepth
         iost = nf_get_vara_double(ncid, id_depth, nf_start, nf_edges,nc_depth)
         if (nc_depth(2) < 0.) nc_depth=-nc_depth
      end if
      call MPI_BCast(iost, 1, MPI_INTEGER, 0, partit%MPI_COMM_FESOM, ierror)      
      call check_nferr(iost,filename,partit)

      call MPI_BCast(nc_lon,   nc_Nlon,   MPI_DOUBLE_PRECISION, 0, partit%MPI_COMM_FESOM, ierror)
      call MPI_BCast(nc_lat,   nc_Nlat,   MPI_DOUBLE_PRECISION, 0, partit%MPI_COMM_FESOM, ierror)
      call MPI_BCast(nc_depth, nc_Ndepth, MPI_DOUBLE_PRECISION, 0, partit%MPI_COMM_FESOM, ierror)

      if (partit%mype==0) then
         iost = nf_close(ncid)
      end if
      call MPI_BCast(iost, 1, MPI_INTEGER, 0, partit%MPI_COMM_FESOM, ierror)
      call check_nferr(iost,filename,partit)

      if (ic_cyclic) then
         nc_lon(1)      =nc_lon(1)-360.
         nc_lon(nc_Nlon)=nc_lon(nc_Nlon)+360.
      end if 
   END SUBROUTINE nc_readGrid

   
   SUBROUTINE nc_ic3d_ini(partit, mesh)
      !!---------------------------------------------------------------------
      !! ** Purpose : inizialization of ocean forcing from NETCDF file
      !!----------------------------------------------------------------------
      IMPLICIT NONE
      type(t_mesh),   intent(in) ,    target    :: mesh
      type(t_partit), intent(inout),  target    :: partit 
      integer                                   :: i
      integer                                   :: elnodes(3)
      real(wp)                                  :: x, y       ! coordinates of elements
      real(kind=WP), allocatable,dimension(:,:) :: cav_nrst_xyz
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"
      
      warn = 0

      if (mype==0) then
         write(*,*) 'reading ',     trim(filename)
         write(*,*) 'variable  : ', trim(varname)
      end if
      
      call nc_readGrid(partit)

      ! prepare nearest coordinates in INfile , save to bilin_indx_i/j
      !_________________________________________________________________________
      ! cavity case
      if (use_cavity) then
        ! compute bilinear index
        do i = 1, myDim_nod2d
    !       ! its a cavity node use extrapolation points of closest cavity line point
            ! exchange the coordinates of the cavity node with the coordinates of the 
            ! closest cavity-line point --> use than these coordinates to estimate 
            ! bilinear interpolation index
            if (ulevels_nod2D(i)>1) then
                x = mesh%cavity_nrst_cavlpnt_xyz(1,i)/rad
                y = mesh%cavity_nrst_cavlpnt_xyz(2,i)/rad
                !!PS if (mype==0) write(*,*) 'xold, yold, xnew, ynew = ',geo_coord_nod2D(1,i),geo_coord_nod2D(2,i),mesh%cavity_nrst_cavlpnt_xyz(1,i),mesh%cavity_nrst_cavlpnt_xyz(2,i)
            ! its no cavity use normal vertice points    
            else
                x = geo_coord_nod2D(1,i)/rad
                y = geo_coord_nod2D(2,i)/rad
            end if 
            if (x<0.)   x=x+360.
            if (x>360.) x=x-360.
            ! find nearest
            if ( x <= nc_lon(nc_Nlon) .and. x >= nc_lon(1) ) then               
                call binarysearch(nc_Nlon, nc_lon, x, bilin_indx_i(i))
            else ! NO extrapolation in space
                bilin_indx_i(i)=-1         
            end if
            if ( y <= nc_lat(nc_Nlat) .and. y >= nc_lat(1) ) then      
                call binarysearch(nc_Nlat, nc_lat, y, bilin_indx_j(i))
            else ! NO extrapolation in space
                bilin_indx_j(i)=-1         
            end if
        end do
       
      !_________________________________________________________________________ 
      ! standard non-cavity case
      else ! use_cavity==.false.
        do i = 1, myDim_nod2d
    !       ! get coordinates of elements
            x  = geo_coord_nod2D(1,i)/rad
            y  = geo_coord_nod2D(2,i)/rad
            if (x<0.)   x=x+360.
            if (x>360.) x=x-360.
            ! find nearest
            if ( x <= nc_lon(nc_Nlon) .and. x >= nc_lon(1) ) then               
                call binarysearch(nc_Nlon, nc_lon, x, bilin_indx_i(i))
            else ! NO extrapolation in space
                bilin_indx_i(i)=-1         
            end if
            if ( y <= nc_lat(nc_Nlat) .and. y >= nc_lat(1) ) then      
                call binarysearch(nc_Nlat, nc_lat, y, bilin_indx_j(i))
            else ! NO extrapolation in space
                bilin_indx_j(i)=-1         
            end if
        end do
      end if   
   END SUBROUTINE nc_ic3d_ini

   SUBROUTINE getcoeffld(tracers, partit, mesh)
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE getcoeffld ***
      !!              
      !! ** Purpose : read fields from files, inetrpolate on model mesh and prepare interpolation coeffients 
      !! ** Method  : 
      !! ** Action  : 
      !!----------------------------------------------------------------------

      USE ieee_arithmetic
      IMPLICIT NONE      
      type(t_mesh),   intent(in),    target   :: mesh
      type(t_partit), intent(inout), target   :: partit 
      type(t_tracer), intent(inout), target   :: tracers      
      integer                                 :: iost !I/O status     
      integer                                 :: ncid      ! netcdf file id
      ! ID dimensions and variables:
      integer                                 :: id_data
      integer                                 :: nf_start(4)
      integer                                 :: nf_edges(4)         
      integer                                 :: fld_idx, i,j,ii, ip1, jp1, k
      integer                                 :: d_indx, d_indx_p1  ! index of neares      
      real(wp)                                :: cf_a, cf_b, delta_d
      integer                                 :: nl1, ul1
      real(wp)                                :: denom, x1, x2, y1, y2, x, y, d1,d2, aux_z           
      real(wp), allocatable, dimension(:,:,:) :: ncdata
      real(wp), allocatable, dimension(:)     :: data1d      
      integer                                 :: elnodes(3)
      integer                                 :: ierror              ! return error code
      integer                                 :: NO_FILL             ! 0=no fillval, 1=fillval
      real(wp)                                :: FILL_VALUE
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"

      ALLOCATE(ncdata(nc_Nlon,nc_Nlat,nc_Ndepth), data1d(nc_Ndepth))
      ncdata=0.0_WP
      data1d=0.0_WP
      tracers%data(current_tracer)%values(:,:)=dummy
      !open NETCDF file on 0 core     
      if (mype==0) then
         iost = nf_open(filename,NF_NOWRITE,ncid)
      end if
      call MPI_BCast(iost, 1, MPI_INTEGER, 0, MPI_COMM_FESOM, ierror)
      call check_nferr(iost,filename,partit)
      ! get variable id
      if (mype==0) then
         iost = nf_inq_varid(ncid, varname, id_data)
         iost = nf_inq_var_fill(ncid, id_data, NO_FILL, FILL_VALUE) ! FillValue defined?
         
         write (*,*) "ncid=", ncid, "id_data=", id_data, "NO_FILL=", NO_FILL, "FILL_VALUE=", FILL_VALUE

         if (NO_FILL==1) then
            print *, 'No _FillValue is set in ', trim(filename), ', trying dummy =', dummy, FILL_VALUE
         else
            print *, 'The FillValue in ', trim(filename), ' is set to ', FILL_VALUE ! should set dummy accordingly
         end if
      end if
      call MPI_BCast(iost, 1, MPI_INTEGER, 0, MPI_COMM_FESOM, ierror)
      call check_nferr(iost,filename,partit)   
      !read data from file
      if (mype==0) then
         nf_start(1)=1
         nf_edges(1)=nc_Nlon-2
         nf_start(2)=1
         nf_edges(2)=nc_Nlat
         nf_start(3)=1
         nf_edges(3)=nc_Ndepth         
         iost = nf_get_vara_double(ncid, id_data, nf_start, nf_edges, ncdata(2:nc_Nlon-1,:,:))
         ncdata(1,:,:)      =ncdata(nc_Nlon-1,:,:)
         ncdata(nc_Nlon,:,:)=ncdata(2,:,:)

         ! replace nan (or fillvalue) by dummy value
         do k=1,nc_Ndepth
            do j=1,nc_Nlat
               do i=1,nc_Nlon
                  if (ieee_is_nan(ncdata(i,j,k)) .or. (ncdata(i,j,k)==FILL_VALUE)) then
                     ncdata(i,j,k) = dummy
                  elseif (ncdata(i,j,k) < -0.99_WP*dummy .or. ncdata(i,j,k) > dummy) then 
                     ! and in case the input data has other conventions on missing values:
                     ncdata(i,j,k) = dummy
                  endif
               end do
            end do
         end do
      end if
      call MPI_BCast(iost, 1, MPI_INTEGER, 0, MPI_COMM_FESOM, ierror)
      call check_nferr(iost,filename,partit)
      call MPI_BCast(ncdata, nc_Nlon*nc_Nlat*nc_Ndepth, MPI_DOUBLE_PRECISION, 0, MPI_COMM_FESOM, ierror)
      ! bilinear space interpolation,  
      ! data is assumed to be sampled on a regular grid
      do ii = 1, myDim_nod2d
         nl1 = nlevels_nod2D(ii)-1
         ul1 = ulevels_nod2D(ii)
         i = bilin_indx_i(ii)
         j = bilin_indx_j(ii)
         ip1 = i + 1   
         jp1 = j + 1
         !______________________________________________________________________
         ! its a cavity node use extrapolation points of closest cavity line point
         ! exchange the coordinates of the cavity node with the coordinates of the 
         ! closest cavity-line point --> use than these coordinates to estimate 
         ! bilinear interpolation index
         if (ul1>1) then
            x = mesh%cavity_nrst_cavlpnt_xyz(1,ii)/rad
            y = mesh%cavity_nrst_cavlpnt_xyz(2,ii)/rad
         ! its no cavity use normal vertice points    
         else
            x = geo_coord_nod2D(1,ii)/rad
            y = geo_coord_nod2D(2,ii)/rad
         end if 
         if (x<0.)   x=x+360.
         if (x>360.) x=x-360.
         if ( min(i,j)>0 ) then
         if (any(ncdata(i:ip1,j:jp1,1) > dummy*0.99_WP)) cycle
            x1 = nc_lon(i)
            x2 = nc_lon(ip1)
            y1 = nc_lat(j)
            y2 = nc_lat(jp1)
                    
            ! if point inside forcing domain
            denom = (x2 - x1)*(y2 - y1)
            data1d(:) = ( ncdata(i,j,:)   * (x2-x)*(y2-y)   + ncdata(ip1,j,:)     * (x-x1)*(y2-y) + &
                        ncdata(i,jp1,:) * (x2-x)*(y-y1)   + ncdata(ip1, jp1, :) * (x-x1)*(y-y1)     ) / denom
            where (ncdata(i,j,:)   > 0.99_WP*dummy .OR. ncdata(ip1,j,:)   > 0.99_WP*dummy .OR. &
                    ncdata(i,jp1,:) > 0.99_WP*dummy .OR. ncdata(ip1,jp1,:) > 0.99_WP*dummy)
                data1d(:)=dummy
            end where   
            
            !___________________________________________________________________
            ! In case of cavity --> do vertical cavity extrapolation for init TS
            if (use_cavity) then
                !!PS do k= 1, nl1
                do k= ul1, nl1
                    if (ul1>1 .and. Z_3d_n(k,ii)<mesh%cavity_nrst_cavlpnt_xyz(3,ii)) then
                        aux_z = mesh%cavity_nrst_cavlpnt_xyz(3,ii)
                    else
                        aux_z = Z_3d_n(k,ii)
                    end if 
                    call binarysearch(nc_Ndepth,nc_depth,-aux_z,d_indx)
                    
                    if ( d_indx < nc_Ndepth .and. d_indx > 0) then
                        d_indx_p1 = d_indx+1
                        delta_d = nc_depth(d_indx+1)-nc_depth(d_indx)
                        ! values from OB data for nearest depth           
                        d1 = data1d(d_indx)
                        d2 = data1d(d_indx_p1)
                        if ((d1<0.99_WP*dummy) .and. (d2<0.99_WP*dummy)) then
                            ! line a*z+b coefficients calculation
                            cf_a  = (d2 - d1)/ delta_d
                            ! value of interpolated OB data on Z from model
                            cf_b  = d1 - cf_a * nc_depth(d_indx)
                            !!PS tracers%data(current_tracer)%values(k,ii) = -cf_a * Z_3d_n(k,ii) + cf_b
                            tracers%data(current_tracer)%values(k,ii) = -cf_a * aux_z + cf_b
                        end if
                    elseif (d_indx==0) then
                        tracers%data(current_tracer)%values(k,ii)=data1d(1)
                    end if
                enddo
            !___________________________________________________________________
            ! normal non-cavity case
            else
                !!PS do k= 1, nl1
                do k= ul1, nl1
                    call binarysearch(nc_Ndepth,nc_depth,-Z_3d_n(k,ii),d_indx)
                    if ( d_indx < nc_Ndepth .and. d_indx > 0) then
                        d_indx_p1 = d_indx+1
                        delta_d = nc_depth(d_indx+1)-nc_depth(d_indx)
                        ! values from OB data for nearest depth           
                        d1 = data1d(d_indx)
                        d2 = data1d(d_indx_p1)
                        if ((d1<0.99_WP*dummy) .and. (d2<0.99_WP*dummy)) then
                            ! line a*z+b coefficients calculation
                            cf_a  = (d2 - d1)/ delta_d
                            ! value of interpolated OB data on Z from model
                            cf_b  = d1 - cf_a * nc_depth(d_indx)
                            tracers%data(current_tracer)%values(k,ii) = -cf_a * Z_3d_n(k,ii) + cf_b
                        end if
                    elseif (d_indx==0) then
                        tracers%data(current_tracer)%values(k,ii)=data1d(1)
                    end if
                enddo
            end if ! --> if (use_cavity) then
         end if ! --> if ( min(i,j)>0 ) then
      end do !ii
      if (mype==0) then
         iost = nf_close(ncid)
      end if
      call MPI_BCast(iost, 1, MPI_INTEGER, 0, MPI_COMM_FESOM, ierror)
      DEALLOCATE( ncdata, data1d )
   END SUBROUTINE getcoeffld  
   
   SUBROUTINE do_ic3d(tracers, partit, mesh)
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE do_ic3d ***
      !!              
      !! ** Purpose : read 3D initial conditions for tracers from netcdf and interpolate on model grid
      !!----------------------------------------------------------------------
      USE insitu2pot_interface 

      IMPLICIT NONE
      type(t_mesh),   intent(in),    target   :: mesh
      type(t_partit), intent(inout), target   :: partit 
      type(t_tracer), intent(inout), target   :: tracers  
      integer                                 :: n, i
      real(kind=WP)                           :: locTmax, locTmin, locSmax, locSmin, glo   
      real(kind=WP)                           :: locDINmax, locDINmin, locDICmax, locDICmin, locAlkmax !OG
      real(kind=WP)                           :: locAlkmin, locDSimax, locDSimin, locDFemax, locDFemin
      real(kind=WP)                           :: locO2min,  locO2max
      real(kind=WP)                           :: locDICremax, locDICremin ! DICremin tracer (added by Sina)


      if (partit%mype==0) write(*,*) "Start: Initial conditions  for tracers"

      ALLOCATE(bilin_indx_i(partit%myDim_nod2d+partit%eDim_nod2D), bilin_indx_j(partit%myDim_nod2d+partit%eDim_nod2D))
      DO n=1, n_ic3d
      filename=trim(ClimateDataPath)//trim(filelist(n))
      varname =trim(varlist(n))
      DO current_tracer=1, tracers%num_tracers
         if (tracers%data(current_tracer)%ID==idlist(n)) then
            ! read initial conditions for current tracer
            call nc_ic3d_ini(partit, mesh)
            ! get first coeficients for time inerpolation on model grid for all datas
            call getcoeffld(tracers, partit, mesh)
            call nc_end ! deallocate arrqays associated with netcdf file
            call extrap_nod(tracers%data(current_tracer)%values(:,:), partit, mesh)
            exit
         elseif (current_tracer==tracers%num_tracers) then
            if (partit%mype==0) write(*,*) "idlist contains tracer which is not listed in tracer_id!"
            if (partit%mype==0) write(*,*) "check your namelists!"
            call par_ex(partit%MPI_COMM_FESOM, partit%mype) 
            stop
         end if
      END DO
      END DO
      DEALLOCATE(bilin_indx_i, bilin_indx_j)

      do current_tracer=1, tracers%num_tracers
         !_________________________________________________________________________
         ! set remaining dummy values from bottom topography to 0.0_WP
         where (tracers%data(current_tracer)%values > 0.9_WP*dummy)
               tracers%data(current_tracer)%values=0.0_WP
         end where
         !_________________________________________________________________________
         ! eliminate values within cavity that result from the extrapolation of 
         ! initialisation
         do n=1,partit%myDim_nod2d + partit%eDim_nod2D
            ! ensure cavity is zero
            if (use_cavity) tracers%data(current_tracer)%values(1:mesh%ulevels_nod2D(n)-1,n)=0.0_WP
            ! ensure bottom is zero
            tracers%data(current_tracer)%values(mesh%nlevels_nod2D(n):mesh%nl-1,n)=0.0_WP
         end do
      end do
      !_________________________________________________________________________
      ! convert temperature from Kelvin --> °C
      where (tracers%data(1)%values(:,:) > 100._WP)
             tracers%data(1)%values(:,:) = tracers%data(1)%values(:,:)-273.15_WP
      end where
      
      !_________________________________________________________________________
      if (t_insitu) then
         if (partit%mype==0) write(*,*) "converting insitu temperature to potential..."
         call insitu2pot(tracers, partit, mesh)
      end if
      if (partit%mype==0) write(*,*) "DONE:  Initial conditions for tracers"             
      !_________________________________________________________________________
      ! check initial fields
      locTmax = -6666
      locTmin = 6666
      locSmax = locTmax
      locSmin = locTmin

#if defined(__recom)
        locDINmax = -66666
        locDINmin = 66666
        locDICmax = locDINmax
        locDICmin = locDINmin
        locAlkmax = locDINmax
        locAlkmin = locDINmin
        locDICremax = locDINmax ! DICremin tracer (added by Sina)
        locDICremin = locDINmin
        locDSimax = locDINmax
        locDSimin = locDINmin
        locDFemax = locDINmax
        locDFemin = locDINmin
        locO2max  = locDINmax
        locO2min  = locDINmin
#endif
      do n=1, partit%myDim_nod2d
        locTmax = max(locTmax,maxval(tracers%data(1)%values(mesh%ulevels_nod2D(n):mesh%nlevels_nod2D(n)-1,n)) )
        locTmin = min(locTmin,minval(tracers%data(1)%values(mesh%ulevels_nod2D(n):mesh%nlevels_nod2D(n)-1,n)) )
        locSmax = max(locSmax,maxval(tracers%data(2)%values(mesh%ulevels_nod2D(n):mesh%nlevels_nod2D(n)-1,n)) )
        locSmin = min(locSmin,minval(tracers%data(2)%values(mesh%ulevels_nod2D(n):mesh%nlevels_nod2D(n)-1,n)) )

#if defined(__recom)
        locDINmax = max(locDINmax,maxval(tracers%data(3)%values(mesh%ulevels_nod2D(n):mesh%nlevels_nod2D(n)-1,n)) )
        locDINmin = min(locDINmin,minval(tracers%data(3)%values(mesh%ulevels_nod2D(n):mesh%nlevels_nod2D(n)-1,n)) )
        locDICmax = max(locDICmax,maxval(tracers%data(4)%values(mesh%ulevels_nod2D(n):mesh%nlevels_nod2D(n)-1,n)) )
        locDICmin = min(locDICmin,minval(tracers%data(4)%values(mesh%ulevels_nod2D(n):mesh%nlevels_nod2D(n)-1,n)) )
        locAlkmax = max(locAlkmax,maxval(tracers%data(5)%values(mesh%ulevels_nod2D(n):mesh%nlevels_nod2D(n)-1,n)) )
        locAlkmin = min(locAlkmin,minval(tracers%data(5)%values(mesh%ulevels_nod2D(n):mesh%nlevels_nod2D(n)-1,n)) )
        locDSimax = max(locDSimax,maxval(tracers%data(20)%values(mesh%ulevels_nod2D(n):mesh%nlevels_nod2D(n)-1,n)) )
        locDSimin = min(locDSimin,minval(tracers%data(20)%values(mesh%ulevels_nod2D(n):mesh%nlevels_nod2D(n)-1,n)) )
        locDFemax = max(locDFemax,maxval(tracers%data(21)%values(mesh%ulevels_nod2D(n):mesh%nlevels_nod2D(n)-1,n)) )
        locDFemin = min(locDFemin,minval(tracers%data(21)%values(mesh%ulevels_nod2D(n):mesh%nlevels_nod2D(n)-1,n)) )
        locO2max  = max(locO2max,maxval(tracers%data(24)%values(mesh%ulevels_nod2D(n):mesh%nlevels_nod2D(n)-1,n)) )
        locO2min  = min(locO2min,minval(tracers%data(24)%values(mesh%ulevels_nod2D(n):mesh%nlevels_nod2D(n)-1,n)) )
        locDICremax = max(locDICremax,maxval(tracers%data(39)%values(mesh%ulevels_nod2D(n):mesh%nlevels_nod2D(n)-1,n)) ) ! DICremin tracer (added by Sina)
        locDICremin = min(locDICremin,minval(tracers%data(39)%values(mesh%ulevels_nod2D(n):mesh%nlevels_nod2D(n)-1,n)) )
#endif
      end do
      call MPI_AllREDUCE(locTmax , glo  , 1, MPI_DOUBLE_PRECISION, MPI_MAX, partit%MPI_COMM_FESOM, partit%MPIerr)
      if (partit%mype==0) write(*,*) '  |-> gobal max init. temp. =', glo
      call MPI_AllREDUCE(locTmin , glo  , 1, MPI_DOUBLE_PRECISION, MPI_MIN, partit%MPI_COMM_FESOM, partit%MPIerr)
      if (partit%mype==0) write(*,*) '  |-> gobal min init. temp. =', glo
      call MPI_AllREDUCE(locSmax , glo  , 1, MPI_DOUBLE_PRECISION, MPI_MAX, partit%MPI_COMM_FESOM, partit%MPIerr)
      if (partit%mype==0) write(*,*) '  |-> gobal max init. salt. =', glo
      call MPI_AllREDUCE(locSmin , glo  , 1, MPI_DOUBLE_PRECISION, MPI_MIN, partit%MPI_COMM_FESOM, partit%MPIerr)
      if (partit%mype==0) write(*,*) '  `-> gobal min init. salt. =', glo      
#if defined(__recom)

      if (partit%mype==0) write(*,*) "Sanity check for REcoM variables"
      call MPI_AllREDUCE(locDINmax , glo  , 1, MPI_DOUBLE_PRECISION, MPI_MAX, partit%MPI_COMM_FESOM, partit%MPIerr)
      if (partit%mype==0) write(*,*) '  |-> gobal max init. DIN. =', glo
      call MPI_AllREDUCE(locDINmin , glo  , 1, MPI_DOUBLE_PRECISION, MPI_MIN, partit%MPI_COMM_FESOM, partit%MPIerr)
      if (partit%mype==0) write(*,*) '  |-> gobal min init. DIN. =', glo

      call MPI_AllREDUCE(locDICmax , glo  , 1, MPI_DOUBLE_PRECISION, MPI_MAX, partit%MPI_COMM_FESOM, partit%MPIerr)
      if (partit%mype==0) write(*,*) '  |-> gobal max init. DIC. =', glo
      call MPI_AllREDUCE(locDICmin , glo  , 1, MPI_DOUBLE_PRECISION, MPI_MIN, partit%MPI_COMM_FESOM, partit%MPIerr)
      if (partit%mype==0) write(*,*) '  |-> gobal min init. DIC. =', glo
      call MPI_AllREDUCE(locAlkmax , glo  , 1, MPI_DOUBLE_PRECISION, MPI_MAX, partit%MPI_COMM_FESOM, partit%MPIerr)
      if (partit%mype==0) write(*,*) '  |-> gobal max init. Alk. =', glo
      call MPI_AllREDUCE(locAlkmin , glo  , 1, MPI_DOUBLE_PRECISION, MPI_MIN, partit%MPI_COMM_FESOM, partit%MPIerr)
      if (partit%mype==0) write(*,*) '  |-> gobal min init. Alk. =', glo
      call MPI_AllREDUCE(locDICremax , glo  , 1, MPI_DOUBLE_PRECISION, MPI_MAX, partit%MPI_COMM_FESOM, partit%MPIerr)
      if (partit%mype==0) write(*,*) '  |-> gobal max init. DICremin. =', glo
      call MPI_AllREDUCE(locDICremin , glo  , 1, MPI_DOUBLE_PRECISION, MPI_MIN, partit%MPI_COMM_FESOM, partit%MPIerr)
      if (partit%mype==0) write(*,*) '  |-> gobal min init. DICremin. =', glo
      call MPI_AllREDUCE(locDSimax , glo  , 1, MPI_DOUBLE_PRECISION, MPI_MAX, partit%MPI_COMM_FESOM, partit%MPIerr)
      if (partit%mype==0) write(*,*) '  |-> gobal max init. DSi. =', glo
      call MPI_AllREDUCE(locDSimin , glo  , 1, MPI_DOUBLE_PRECISION, MPI_MIN, partit%MPI_COMM_FESOM, partit%MPIerr)
      if (partit%mype==0) write(*,*) '  |-> gobal min init. DSi. =', glo
      call MPI_AllREDUCE(locDFemax , glo  , 1, MPI_DOUBLE_PRECISION, MPI_MAX, partit%MPI_COMM_FESOM, partit%MPIerr)
      if (partit%mype==0) write(*,*) '  |-> gobal max init. DFe. =', glo
      call MPI_AllREDUCE(locDFemin , glo  , 1, MPI_DOUBLE_PRECISION, MPI_MIN, partit%MPI_COMM_FESOM, partit%MPIerr)
      if (partit%mype==0) write(*,*) '  `-> gobal min init. DFe. =', glo
      call MPI_AllREDUCE(locO2max , glo  , 1, MPI_DOUBLE_PRECISION, MPI_MAX, partit%MPI_COMM_FESOM, partit%MPIerr)
      if (partit%mype==0) write(*,*) '  |-> gobal max init. O2. =', glo
      call MPI_AllREDUCE(locO2min , glo  , 1, MPI_DOUBLE_PRECISION, MPI_MIN, partit%MPI_COMM_FESOM, partit%MPIerr)
      if (partit%mype==0) write(*,*) '  `-> gobal min init. O2. =', glo
#endif
   END SUBROUTINE do_ic3d
    
   SUBROUTINE nc_end

    IMPLICIT NONE
 
    DEALLOCATE(nc_lon, nc_lat, nc_depth)

   END SUBROUTINE nc_end

   SUBROUTINE check_nferr(iost,fname, partit)
   IMPLICIT NONE
      type(t_partit),          intent(inout) :: partit 
      character(len=MAX_PATH), intent(in)    :: fname
      integer, intent(in)                    :: iost
      if (iost .ne. NF_NOERR) then
         write(*,*) 'ERROR: I/O status= "',trim(nf_strerror(iost)),'";',iost,' file= ', trim(fname)
         call par_ex(partit%MPI_COMM_FESOM, partit%mype)
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

      d = 1e-9_WP
      left = 1
      right = length
      do
         if (left > right) then
            exit
         endif
         middle = nint((left+right) / 2.0_WP)
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
