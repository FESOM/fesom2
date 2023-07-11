module diag_uv_norm3
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
   INTEGER                                    	:: i, lev, el, elnodes(3)
   character(100)                             	:: filein, fileout
   real(kind=8)            		      	:: x, y
   real(kind=8), dimension(:,:),   allocatable	:: u, v
   real(kind=8), dimension(:,:),   allocatable	:: arr_reg
   real(kind=8), dimension(:), 	   allocatable	:: arr, vol
   character(15)                              	:: varid2
   character(100)				:: svarid_in(2), svarid_out(2)
   integer					:: NFIELDS_IN, NFIELDS_OUT

  PRIVATE
  PUBLIC :: make_diag_uv_norm3
CONTAINS 
subroutine make_diag_uv_norm3
  implicit none
   NFIELDS_IN =1
   NFIELDS_OUT=1
   allocate(u(nl_1, elem2D), v(nl_1, elem2D))

   allocate(arr(nod2D), vol(nod2D))
   allocate(arr_reg(reg_nx, reg_ny))
   
   svarid_in(1) ='u'
   svarid_in(2) ='v'
   svarid_out(1)='uv_norm'

   fileout=trim(outpath)//'uv_norm.nc'

   vol	  =0.
   do el=1, elem2D
      elnodes=elem2D_nodes(:, el)
      vol(elnodes)=vol(elnodes)+voltriangle(el)
   end do

   do year=year_start, year_end
   write(filein,   '(a6,I4,a7)') trim(runid)//'.',year, '.oce.nc'
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
           do el=1, elem2D
              elnodes=elem2D_nodes(:, el)
              arr(elnodes)=arr(elnodes)+sqrt(u(lev, el)**2+v(lev, el)**2)*voltriangle(el)
           end do
           arr=arr(:)/vol
           do i=1, NFIELDS_OUT
              call do_oce_2_reg(arr, arr_reg(:, :), 1)
           end do      
           call netcdf_write(fileout, lev, (year==year_start .AND. month==1 .AND. lev==1), lev==1)
        end do
   end do
   end do
   io=nf_close(ncid)   
   deallocate(arr_reg, arr, vol, v, u)
end subroutine make_diag_uv_norm3
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
end module diag_uv_norm3
