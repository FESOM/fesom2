module diag_TS3
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
   character(100)                             	:: filein, fileout
   real(kind=8)            		      	:: x, y
   real(kind=8), dimension(:),     allocatable	:: arr
   real(kind=8), dimension(:,:,:), allocatable	:: arr_reg
   character(15)                              	:: varid2
   character(100)				:: svarid_in(100), svarid_out(100)
   integer					:: NFIELDS_IN, NFIELDS_OUT

  PRIVATE
  PUBLIC :: make_diag_TS3
CONTAINS 
subroutine make_diag_TS3
  implicit none
  integer        :: i, lev
  integer        :: n

   NFIELDS_IN =2
   NFIELDS_OUT=2
   allocate(arr(nod2D))
   allocate(arr_reg(reg_nx, reg_ny, NFIELDS_IN))
   
   svarid_in(1) ='temp'
   svarid_in(2) ='salt'
   svarid_out(1)='temp'
   svarid_out(2)='salt'

   fileout=trim(outpath)//'ts3.nc'

   do year=year_start, year_end
   write(filein,   '(a6,I4,a7)') trim(runid)//'.', year, '.oce.nc'
   io=nf_open(trim(datapath)//trim(filein), nf_nowrite, ncid)
   if (io.ne.nf_noerr)then
      write(*,*), 'ERROR: CANNOT READ DATA FILE CORRECTLY !!!!!'
      write(*,*), 'Error in opening netcdf file'//trim(datapath)//trim(filein)
      stop 'Fatal error in open_netcdf'
   endif

   do month=1, snap_per_year
      do lev=1, nl_1
         do i=1, NFIELDS_IN
            varid2=svarid_in(i)
	    call read_fesom_lev(ncid, month, varid2, nod2D, lev, arr)
            write(*,*) trim(svarid_in(i)), ': ', minval(arr), maxval(arr)
            call do_oce_2_reg(arr, arr_reg(:, :, i), 1)
         end do
         call netcdf_write (fileout, lev, (year==year_start .AND. month==1 .AND. lev==1), lev==1)
      end do
   end do
   io=nf_close(ncid)   
   end do
   deallocate(arr_reg, arr)
end subroutine make_diag_TS3
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
     stat(s)=NF_PUT_VARA_DOUBLE(fileid,var, (/1, 1, lev, month/), (/reg_nx, reg_ny, 1, 1/), arr_reg(:,:,n)); s=s+1
  end do     
  do i=1,s-1
     if (stat(i).ne.NF_NOERR) write(*,*) &
        'NetCDF error in writing variable, no.', i
  end do  
  stat(1)=NF_CLOSE(fileid); s=s+1    
end SUBROUTINE netcdf_write
end module diag_TS3
