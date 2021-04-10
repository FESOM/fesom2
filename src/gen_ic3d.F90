
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
   USE o_PARAM
   USE g_PARSUP
   USE g_comm_auto
   USE g_support
   USE g_config, only: dummy, ClimateDataPath, use_cavity
   
   IMPLICIT NONE

   include 'netcdf.inc'

   public  do_ic3d, &                                       ! read and apply 3D initial conditions
           n_ic3d, idlist, filelist, varlist, oce_init3d, & ! to be read from the namelist
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

   namelist / oce_init3d / n_ic3d, idlist, filelist, varlist, t_insitu

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
         iost = nf_inq_dimid(ncid,    "LAT",      id_latd)
         if (iost .ne. NF_NOERR) then
            iost = nf_inq_dimid(ncid, "lat",      id_latd)
         end if
         if (iost .ne. NF_NOERR) then
            iost = nf_inq_dimid(ncid, "latitude", id_latd)
         end if
      end if
      call MPI_BCast(iost, 1, MPI_INTEGER, 0, MPI_COMM_FESOM, ierror)
      call check_nferr(iost,filename)
      if (mype==0) then 
         iost = nf_inq_dimid(ncid,    "LON",       id_lond)
         if      (iost .ne. NF_NOERR) then
            iost = nf_inq_dimid(ncid, "longitude", id_lond)
         end if
         if (iost .ne. NF_NOERR) then
            iost = nf_inq_dimid(ncid, "lon",       id_lond)
         end if
      end if
      call MPI_BCast(iost, 1, MPI_INTEGER, 0, MPI_COMM_FESOM, ierror)
      call check_nferr(iost,filename) 
      if (mype==0) then   
         iost = nf_inq_dimid(ncid, "depth", id_depthd)
      end if
      call MPI_BCast(iost, 1, MPI_INTEGER, 0, MPI_COMM_FESOM, ierror)
      call check_nferr(iost,filename)  

      ! get variable id
      call MPI_BCast(iost, 1, MPI_INTEGER, 0, MPI_COMM_FESOM, ierror)
      call check_nferr(iost,filename)
      if (mype==0) then
         iost = nf_inq_varid(ncid,    "LAT",      id_lat)
         if     (iost .ne. NF_NOERR) then
            iost = nf_inq_varid(ncid, "lat",      id_lat)
         end if
         if (iost .ne. NF_NOERR) then
            iost = nf_inq_varid(ncid, "latitude", id_lat)
         end if
      end if
      if (mype==0) then
         iost = nf_inq_varid(ncid,    "LON",       id_lon)
         if      (iost .ne. NF_NOERR) then
            iost = nf_inq_varid(ncid, "longitude", id_lon)
         end if
         if (iost .ne. NF_NOERR) then
            iost = nf_inq_varid(ncid, "lon",       id_lon)
         end if
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
      nc_Nlon=nc_Nlon+2 !for the halo in case of periodic boundary
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
         nf_edges(1)=nc_Nlon-2
         iost = nf_get_vara_double(ncid, id_lon, nf_start, nf_edges, nc_lon(2:nc_Nlon-1))
         nc_lon(1)        =nc_lon(nc_Nlon-1)
         nc_lon(nc_Nlon)  =nc_lon(2)
      end if
      call MPI_BCast(iost, 1, MPI_INTEGER, 0, MPI_COMM_FESOM, ierror)      
      call check_nferr(iost,filename)
   ! depth
      if (mype==0) then
         nf_start(1)=1
         nf_edges(1)=nc_Ndepth
         iost = nf_get_vara_double(ncid, id_depth, nf_start, nf_edges,nc_depth)
         if (nc_depth(2) < 0.) nc_depth=-nc_depth
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

      if (ic_cyclic) then
         nc_lon(1)      =nc_lon(1)-360.
         nc_lon(nc_Nlon)=nc_lon(nc_Nlon)+360.
      end if 
   END SUBROUTINE nc_readGrid

   
   SUBROUTINE nc_ic3d_ini(mesh)
      !!---------------------------------------------------------------------
      !! ** Purpose : inizialization of ocean forcing from NETCDF file
      !!----------------------------------------------------------------------
      IMPLICIT NONE
   
      integer                  :: i
      integer                  :: elnodes(3)
      real(wp)                 :: x, y       ! coordinates of elements
      real(kind=WP), allocatable,dimension(:,:) :: cav_nrst_xyz
      type(t_mesh), intent(in), target :: mesh
#include "associate_mesh.h"
      
      warn = 0

      if (mype==0) then
         write(*,*) 'reading input tracer file for tracer ID= ', tracer_ID(current_tracer)
         write(*,*) 'input file: ', trim(filename)
         write(*,*) 'variable  : ', trim(varname)
      end if
      
      call nc_readGrid

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

   SUBROUTINE getcoeffld(mesh)
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
      real(wp)             :: cf_a, cf_b, delta_d
      integer              :: nl1, ul1
      real(wp)             :: denom, x1, x2, y1, y2, x, y, d1,d2, aux_z     
      
      real(wp), allocatable, dimension(:,:,:)  :: ncdata
      real(wp), allocatable, dimension(:)      :: data1d      
      integer              :: elnodes(3)
      integer              :: ierror              ! return error code

      type(t_mesh), intent(in), target :: mesh
#include "associate_mesh.h"

      ALLOCATE(ncdata(nc_Nlon,nc_Nlat,nc_Ndepth), data1d(nc_Ndepth))
      ncdata=0.0_WP
      data1d=0.0_WP
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
         nf_edges(1)=nc_Nlon-2
         nf_start(2)=1
         nf_edges(2)=nc_Nlat
         nf_start(3)=1
         nf_edges(3)=nc_Ndepth         
         iost = nf_get_vara_double(ncid, id_data, nf_start, nf_edges, ncdata(2:nc_Nlon-1,:,:))
         ncdata(1,:,:)      =ncdata(nc_Nlon-1,:,:)
         ncdata(nc_Nlon,:,:)=ncdata(2,:,:)
         where (ncdata < -0.99_WP*dummy ) ! dummy values are only positive
                ncdata = dummy
         end where
      end if
      call MPI_BCast(iost, 1, MPI_INTEGER, 0, MPI_COMM_FESOM, ierror)
      call check_nferr(iost,filename)
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
!!PS          x  = geo_coord_nod2D(1,ii)/rad
!!PS          y  = geo_coord_nod2D(2,ii)/rad
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
                            !!PS tr_arr(k,ii,current_tracer) = -cf_a * Z_3d_n(k,ii) + cf_b
                            tr_arr(k,ii,current_tracer) = -cf_a * aux_z + cf_b
                        end if
                    elseif (d_indx==0) then
                        tr_arr(k,ii,current_tracer)=data1d(1)
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
                            tr_arr(k,ii,current_tracer) = -cf_a * Z_3d_n(k,ii) + cf_b
                        end if
                    elseif (d_indx==0) then
                        tr_arr(k,ii,current_tracer)=data1d(1)
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
   
   SUBROUTINE do_ic3d(mesh)
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE do_ic3d ***
      !!              
      !! ** Purpose : read 3D initial conditions for tracers from netcdf and interpolate on model grid
      !!----------------------------------------------------------------------
      IMPLICIT NONE
      integer                       :: n, i
      type(t_mesh), intent(in)     , target :: mesh
      real(kind=WP)                 :: locTmax, locTmin, locSmax, locSmin, glo

      if (mype==0) write(*,*) "Start: Initial conditions  for tracers"

      ALLOCATE(bilin_indx_i(myDim_nod2d+eDim_nod2D), bilin_indx_j(myDim_nod2d+eDim_nod2D))
      DO n=1, n_ic3d
      filename=trim(ClimateDataPath)//trim(filelist(n))
      varname =trim(varlist(n))
      DO current_tracer=1, num_tracers
         if (tracer_ID(current_tracer)==idlist(n)) then
            ! read initial conditions for current tracer
            call nc_ic3d_ini(mesh)
            ! get first coeficients for time inerpolation on model grid for all datas
            call getcoeffld(mesh)
            call nc_end ! deallocate arrqays associated with netcdf file
            call extrap_nod(tr_arr(:,:,current_tracer), mesh)
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

      !_________________________________________________________________________
      ! set remaining dummy values from bottom topography to 0.0_WP
      where (tr_arr > 0.9_WP*dummy)
            tr_arr=0.0_WP
      end where
      
      !_________________________________________________________________________
      ! convert temperature from Kelvin --> °C
      where (tr_arr(:,:,1) > 100._WP)
         tr_arr(:,:,1)=tr_arr(:,:,1)-273.15_WP
      end where
      
      !_________________________________________________________________________
      ! eliminate values within cavity that result from the extrapolation of 
      ! initialisation
      do n=1,myDim_nod2d + eDim_nod2D
            ! ensure cavity is zero
            if (use_cavity) tr_arr(1:mesh%ulevels_nod2D(n)-1,n,:)=0.0_WP
            ! ensure bottom is zero
            tr_arr(mesh%nlevels_nod2D(n):mesh%nl-1,n,:)=0.0_WP            
      end do 
      
      !_________________________________________________________________________
      if (t_insitu) then
         if (mype==0) write(*,*) "converting insitu temperature to potential..."
         call insitu2pot(mesh)
      end if
      if (mype==0) write(*,*) "DONE:  Initial conditions for tracers"
      
      !_________________________________________________________________________
      ! Homogenous temp salt initialisation --> for testing and debuging
!!PS       do n=1,myDim_nod2d + eDim_nod2D
!!PS             tr_arr(mesh%ulevels_nod2D(n):mesh%nlevels_nod2D(n)-1,n,1) = 16.0
!!PS             tr_arr(mesh%ulevels_nod2D(n):mesh%nlevels_nod2D(n)-1,n,2) = 35.0
!!PS       end do 
        
      !_________________________________________________________________________
      ! check initial fields
      locTmax = -6666
      locTmin = 6666
      locSmax = locTmax
      locSmin = locTmin
      do n=1,myDim_nod2d
!!PS         if (any( tr_arr(mesh%ulevels_nod2D(n):mesh%nlevels_nod2D(n)-1,n,2)>0.99_WP*dummy)) then
!!PS             write(*,*) '____________________________________________________________'
!!PS             write(*,*) ' --> check init fields SALT >0.99_WP*dummy'
!!PS             write(*,*) 'mype =',mype
!!PS             write(*,*) 'n    =',n
!!PS             write(*,*) 'lon,lat               =',mesh%geo_coord_nod2D(:,n)/rad
!!PS             write(*,*) 'mesh%ulevels_nod2D(n) =',mesh%ulevels_nod2D(n)
!!PS             write(*,*) 'mesh%nlevels_nod2D(n) =',mesh%nlevels_nod2D(n)
!!PS             write(*,*) 'tr_arr(unl:lnl,n,2) =',tr_arr(mesh%ulevels_nod2D(n):mesh%nlevels_nod2D(n)-1,n,2)
!!PS             write(*,*) 'tr_arr(  1:lnl,n,2) =',tr_arr(1:mesh%nlevels_nod2D(n)-1,n,2)
!!PS         end if 
!!PS         if (any( tr_arr(mesh%ulevels_nod2D(n):mesh%nlevels_nod2D(n)-1,n,1)>0.99_WP*dummy)) then
!!PS             write(*,*) '____________________________________________________________'
!!PS             write(*,*) ' --> check init fields TEMP >0.99_WP*dummy'
!!PS             write(*,*) 'mype =',mype
!!PS             write(*,*) 'n    =',n
!!PS             write(*,*) 'lon,lat               =',mesh%geo_coord_nod2D(:,n)/rad
!!PS             write(*,*) 'mesh%ulevels_nod2D(n) =',mesh%ulevels_nod2D(n)
!!PS             write(*,*) 'mesh%nlevels_nod2D(n) =',mesh%nlevels_nod2D(n)
!!PS             write(*,*) 'tr_arr(:,n,1) =',tr_arr(mesh%ulevels_nod2D(n):mesh%nlevels_nod2D(n)-1,n,1)
!!PS         end if 
        locTmax = max(locTmax,maxval(tr_arr(mesh%ulevels_nod2D(n):mesh%nlevels_nod2D(n)-1,n,1)) )
        locTmin = min(locTmin,minval(tr_arr(mesh%ulevels_nod2D(n):mesh%nlevels_nod2D(n)-1,n,1)) )
        locSmax = max(locSmax,maxval(tr_arr(mesh%ulevels_nod2D(n):mesh%nlevels_nod2D(n)-1,n,2)) )
        locSmin = min(locSmin,minval(tr_arr(mesh%ulevels_nod2D(n):mesh%nlevels_nod2D(n)-1,n,2)) )
      end do
      call MPI_AllREDUCE(locTmax , glo  , 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_FESOM, MPIerr)
      if (mype==0) write(*,*) '  |-> gobal max init. temp. =', glo
      call MPI_AllREDUCE(locTmin , glo  , 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_FESOM, MPIerr)
      if (mype==0) write(*,*) '  |-> gobal min init. temp. =', glo
      call MPI_AllREDUCE(locSmax , glo  , 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_FESOM, MPIerr)
      if (mype==0) write(*,*) '  |-> gobal max init. salt. =', glo
      call MPI_AllREDUCE(locSmin , glo  , 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_FESOM, MPIerr)
      if (mype==0) write(*,*) '  `-> gobal min init. salt. =', glo
      
  
   END SUBROUTINE do_ic3d
   
   SUBROUTINE err_call(iost,fname)
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE  err_call ***
      !!----------------------------------------------------------------------
   IMPLICIT NONE
      integer, intent(in)            :: iost
      character(len=MAX_PATH), intent(in) :: fname
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
      character(len=MAX_PATH), intent(in) :: fname
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
