module diag_uv_curl3
  use o_param
  use o_elements
  use o_mesh
  use o_array   
  use g_config
  use g_rotate_grid
  use g_oce_2_reg   
   implicit none
   save
#include "netcdf.inc"

   INTEGER                                    	:: year, month, io, ncid
   character(100)                             	:: filein, fileout(2)
   real(kind=8)            		      	:: x, y
   real(kind=8), dimension(:,:),   allocatable	:: u, v
   real(kind=8), dimension(:,:),   allocatable	:: arr_reg
   real(kind=8), dimension(:), 	   allocatable	:: arr, vol
   character(15)                              	:: varid2
   character(100)				:: svarid_in(2), svarid_out(2)
   integer					:: NFIELDS_IN, NFIELDS_OUT

  PRIVATE
  PUBLIC :: make_diag_uv_curl3
CONTAINS 
subroutine make_diag_uv_curl3
  implicit none
  integer        :: i, lev
  integer        :: n, el(2), enodes(2), nl1, nl2, ed
  real(kind=8)   :: deltaX1, deltaY1, deltaX2, deltaY2, c1

   NFIELDS_IN =1
   NFIELDS_OUT=1
   allocate(u(nl_1, elem2D), v(nl_1, elem2D))

   allocate(arr(nod2D), vol(nod2D))
   allocate(arr_reg(reg_nx, reg_ny))
   
   svarid_in(1) ='u'
   svarid_in(2) ='v'
   svarid_out(1)='uv_norm'

   fileout(1)=trim(outpath)//'uv_curl.nc'
   fileout(2)=trim(outpath)//'uv_ncurl.nc'

   do year=year_start, year_end
   write(filein,   '(a6,I4,a7)') trim(runid)//'.', year, '.oce.nc'
   io=nf_open(trim(datapath)//trim(filein), nf_nowrite, ncid)
   if (io.ne.nf_noerr)then
      write(*,*), 'ERROR: CANNOT READ DATA FILE CORRECTLY !!!!!'
      write(*,*), 'Error in opening netcdf file'//trim(datapath)//trim(filein)
      stop 'Fatal error in open_netcdf'
   endif

   do month=1, snap_per_year
	varid2=svarid_in(1)
	call read_fesom(ncid, month, varid2, elem2d, nl_1, u)
	varid2=svarid_in(2)
	call read_fesom(ncid, month, varid2, elem2d, nl_1, v)

        do lev=1, nl_1
           arr=0.
           DO ed=1, edge2D
              enodes=edges(:,ed)
	      el=edge_tri(:,ed)
              nl1=elvls(el(1))-1
              deltaX1=edge_cross_dxdy(1,ed)
              deltaY1=edge_cross_dxdy(2,ed)
              nl2=0
              if (el(2)>0) then
                 deltaX2=edge_cross_dxdy(3,ed)
                 deltaY2=edge_cross_dxdy(4,ed)
                 nl2=elvls(el(2))-1
              end if
              if (lev<=nl1) then
                 c1=deltaX1*U(lev,el(1))+deltaY1*V(lev,el(1))
                 arr(enodes(1))=arr(enodes(1))+c1
                 arr(enodes(2))=arr(enodes(2))-c1
              END IF
              if (lev<=nl2) then
                 c1= -deltaX2*U(lev,el(2))-deltaY2*V(lev,el(2))
                 arr(enodes(1))=arr(enodes(1))+c1
                 arr(enodes(2))=arr(enodes(2))-c1
              END IF
           END DO
           DO n=1, nod2D
              arr(n)=arr(n)/area(lev,n)
           end do
           write(*,*) minval(arr), maxval(arr)
           do i=1, NFIELDS_OUT
              call do_oce_2_reg(arr, arr_reg(:, :), 1)
           end do
           call netcdf_write (fileout(1), lev, (year==year_start .AND. month==1 .AND. lev==1), lev==1)
           call netcdf_nwrite(fileout(2), lev, (year==year_start .AND. month==1 .AND. lev==1), lev==1)
        end do
   end do
   end do
   io=nf_close(ncid)   
   deallocate(arr_reg, arr, vol, v, u)
end subroutine make_diag_uv_curl3
! ============================================================
SUBROUTINE netcdf_write(filename, lev, newold, newmonth)
  USE o_param
  USE o_mesh
  USE o_elements
  use o_array  
  use g_oce_2_reg  
  IMPLICIT NONE
  save
  
#include "netcdf.inc" 

  !  IDs
  INTEGER, save :: fileid
  ! dimensions
  character(100), INTENT(IN)  :: filename
  LOGICAL,        INTENT(IN)  :: newold, newmonth
  INTEGER,        INTENT(IN)  :: lev
  INTEGER       :: i,j,k,s,n
  INTEGER       :: lon_id, lat_id, dep_id, tim_id
  INTEGER       :: var
  INTEGER, save :: month=0
  INTEGER, DIMENSION(4)       :: dimarray
  INTEGER, DIMENSION(50)      :: stat    !  status array

  if (newold) then     
  stat(1) = NF_CREATE(TRIM(filename),0,fileid); s=s+1  
  if (stat(1).ne.NF_NOERR) then
     write(*,*) 'NetCDF error while opening file - in cpl_str_nc_init'
  end if

  s=1
  stat(s) = NF_DEF_DIM(fileid,'lons',  reg_nx,       dimarray(1)); s=s+1
  stat(s) = NF_DEF_DIM(fileid,'lats',  reg_ny,       dimarray(2)); s=s+1
  stat(s) = NF_DEF_DIM(fileid,'deps',  nl_1,         dimarray(3)); s=s+1
  stat(s) = NF_DEF_DIM(fileid,'time',  NF_UNLIMITED, dimarray(4)); s=s+1    

  do i=1,s-1
     if (stat(i).ne.NF_NOERR) write(*,*) &
        'NetCDF error in dimension definitions, no.',i
  end do

  s=1
  stat(s) = NF_DEF_VAR(fileid,'lons',    NF_FLOAT,1,dimarray(1), lon_id); s=s+1
  stat(s) = NF_DEF_VAR(fileid,'lats',    NF_FLOAT,1,dimarray(2), lat_id); s=s+1
  stat(s) = NF_DEF_VAR(fileid,'deps',    NF_FLOAT,1,dimarray(3), dep_id); s=s+1
  stat(s) = NF_DEF_VAR(fileid,'time',    NF_INT,  0,dimarray(4), tim_id); s=s+1  
  
  do n=1, NFIELDS_OUT 
     stat(s)=  NF_DEF_VAR(fileid, trim(svarid_out(n)), NF_DOUBLE, 4, dimarray(1:4), var); s=s+1
  end do
  do i=1,s-1
     if (stat(i).ne.NF_NOERR) write(*,*) &
        'NetCDF error in str. variable definition, no.',i
  end do
  stat(1)=NF_ENDDEF(fileid)
  
  s=1;
  stat(s)=NF_PUT_VAR_REAL(fileid,lon_id,real(reg_lon/rad, 4)); s=s+1
  stat(s)=NF_PUT_VAR_REAL(fileid,lat_id,real(reg_lat/rad, 4)); s=s+1  
  stat(s)=NF_PUT_VAR_REAL(fileid,dep_id,real((depths(1:nl_1)+depths(2:nl_1+1))/2., 4)); s=s+1  
  stat(s)=nf_put_var_int(fileid, tim_id, 0); s=s+1
  do i=1,s-1
     if (stat(i).ne.NF_NOERR) write(*,*) &
        'NetCDF error in str. variable definition, no.',i
  end do
  stat(1)=NF_CLOSE(fileid)
  end if

  s=1
  stat(1) = NF_OPEN(TRIM(filename),  NF_WRITE, FileId); s=s+1
  stat(s) = NF_INQ_VARID   (fileid, 'time', tim_id);s=s+1
  stat(s) = nf_get_var_int (fileid,  tim_id, month); s=s+1   

  if (newmonth) then
     month=month+1
     stat(s) = nf_put_var_int (fileid, tim_id, month); s=s+1
  end if
  !------  Writing Variables -----
  do n=1, NFIELDS_OUT 
     stat(s) = NF_INQ_VARID   (fileid, trim(svarid_out(n)), var);s=s+1
     stat(s)=NF_PUT_VARA_DOUBLE(fileid,var, (/1, 1, lev, month/), (/reg_nx, reg_ny, 1, 1/), arr_reg(:,:)); s=s+1
  end do     
  do i=1,s-1
     if (stat(i).ne.NF_NOERR) write(*,*) &
        'NetCDF error in writing variable, no.', i
  end do  
  stat(1)=NF_CLOSE(fileid); s=s+1    
end SUBROUTINE netcdf_write
SUBROUTINE netcdf_nwrite(filename, lev, newold, newmonth)
  USE o_param
  USE o_mesh
  USE o_elements
  use o_array  
  use g_oce_2_reg  
  IMPLICIT NONE
  save
  
#include "netcdf.inc" 

  !  IDs
  INTEGER, save :: fileid
  ! dimensions
  character(100), INTENT(IN)  :: filename
  LOGICAL,        INTENT(IN)  :: newold, newmonth
  INTEGER,        INTENT(IN)  :: lev
  INTEGER       :: i,j,k,s,n
  INTEGER       :: lon_id, lat_id, dep_id, tim_id, n2d
  INTEGER       :: var
  INTEGER, save :: month=0
  INTEGER, DIMENSION(4)       :: dimarray
  INTEGER, DIMENSION(50)      :: stat    !  status array

  if (newold) then     
  stat(1) = NF_CREATE(TRIM(filename),0,fileid); s=s+1  
  if (stat(1).ne.NF_NOERR) then
     write(*,*) 'NetCDF error while opening file - in cpl_str_nc_init'
  end if

  s=1

  stat(s) = NF_DEF_DIM(fileid,'deps',  nl_1,         dimarray(1)); s=s+1
  stat(s) = NF_DEF_DIM(fileid,'n2d',   nod2D,        dimarray(2)); s=s+1
  stat(s) = NF_DEF_DIM(fileid,'time',  NF_UNLIMITED, dimarray(3)); s=s+1

  do i=1,s-1
     if (stat(i).ne.NF_NOERR) write(*,*) &
        'NetCDF error in dimension definitions, no.',i
  end do

  s=1
  stat(s) = NF_DEF_VAR(fileid,'lon',   NF_FLOAT,1,   dimarray(2), lon_id); s=s+1
  stat(s) = NF_DEF_VAR(fileid,'lat',   NF_FLOAT,1,   dimarray(2), lat_id); s=s+1
  stat(s) = NF_DEF_VAR(fileid,'deps',  NF_FLOAT,1,   dimarray(1), dep_id); s=s+1
  stat(s) = NF_DEF_VAR(fileid,'time',  NF_INT,  0,   dimarray(3), tim_id); s=s+1  
 
  do n=1, NFIELDS_OUT 
     stat(s)=  NF_DEF_VAR(fileid, trim(svarid_out(n)), NF_DOUBLE, 3, dimarray(1:3), var); s=s+1
  end do
  do i=1,s-1
     if (stat(i).ne.NF_NOERR) write(*,*) &
        'NetCDF error in str. variable definition, no.',i
  end do
  stat(1)=NF_ENDDEF(fileid)
  
  s=1;
  stat(s)=NF_PUT_VAR_REAL(fileid,lon_id,real(coord_nod2D(1,:),4)); s=s+1
  stat(s)=NF_PUT_VAR_REAL(fileid,lat_id,real(coord_nod2D(1,:),4)); s=s+1  
  stat(s)=NF_PUT_VAR_REAL(fileid,dep_id,real((depths(1:nl_1)+depths(2:nl_1+1))/2., 4)); s=s+1  
  stat(s)=nf_put_var_int(fileid, tim_id, 0); s=s+1
  do i=1,s-1
     if (stat(i).ne.NF_NOERR) write(*,*) &
        'NetCDF error in str. variable definition, no.',i
  end do
  stat(1)=NF_CLOSE(fileid)
  end if

  s=1
  stat(1) = NF_OPEN(TRIM(filename),  NF_WRITE, FileId); s=s+1
  stat(s) = NF_INQ_VARID   (fileid, 'time', tim_id);s=s+1
  stat(s) = nf_get_var_int (fileid,  tim_id, month); s=s+1   

  if (newmonth) then
     month=month+1
     stat(s) = nf_put_var_int (fileid, tim_id, month); s=s+1
  end if
  !------  Writing Variables -----
  do n=1, NFIELDS_OUT 
     stat(s) = NF_INQ_VARID   (fileid, trim(svarid_out(n)), var);s=s+1
     stat(s)=NF_PUT_VARA_DOUBLE(fileid,var, (/lev, 1, month/), (/1, nod2D, 1/), arr); s=s+1
  end do     
  do i=1,s-1
     if (stat(i).ne.NF_NOERR) write(*,*) &
        'NetCDF error in writing variable, no.', i
  end do  
  stat(1)=NF_CLOSE(fileid); s=s+1    
end SUBROUTINE netcdf_nwrite
end module diag_uv_curl3
