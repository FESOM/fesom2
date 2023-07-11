module diag_moc_w
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
   INTEGER                                    	:: i, lev, pos, el, elnodes(3)
   character(100)                             	:: filein, fileout
   real(kind=8)            		      	:: x, y, vol
   real(kind=8), dimension(:,:),   allocatable	:: w, arr
   real(kind=8), dimension(:,:),   allocatable	:: moc
   character(15)                              	:: varid2
   character(100)				:: svarid_in(1), svarid_out(1)
   integer					:: NFIELDS_IN, NFIELDS_OUT


  PRIVATE
  PUBLIC :: make_diag_moc_w
CONTAINS   
subroutine make_diag_moc_w
  implicit none

   NFIELDS_OUT=1
   svarid_out(1)='MOC'

   allocate(w(max_num_layers, nod2D), arr(max_num_layers, nod2D))
   allocate(moc(reg_ny, max_num_layers))

   fileout=trim(outpath)//'diag_moc.nc'

   do year=year_start, year_end
   write(filein,   '(a6,I4,a7)') trim(runid)//'.',year, '.oce.nc'
   io=nf_open(trim(datapath)//trim(filein), nf_nowrite, ncid)
   if (io.ne.nf_noerr)then
      write(*,*), 'ERROR: CANNOT READ DATA FILE CORRECTLY !!!!!'
      write(*,*), 'Error in opening netcdf file'//trim(datapath)//trim(filein)
      write(*,*), 'io=', io
      stop 'Fatal error in open_netcdf:'
   endif
   write(*,*) 'reading file:'
   write(*,*) trim(datapath)//trim(filein)
   w=0.
   moc=0.
   do month=1, snap_per_year
	varid2='w'
	call read_fesom(ncid, month, varid2, nod2D, max_num_layers, arr)
        w=w+arr
   end do 
   w=w/real(snap_per_year)
   write(*,*) 'min/max for w are:'
   write(*,*) minval(w), maxval(w)

   do el=1, elem2D
      elnodes=elem2D_nodes(:, el)
      if (rotated_grid) then
       call r2g(x, y, sum(coord_nod2D(1, elnodes))/3., sum(coord_nod2D(2, elnodes))/3.)
      else
         x=sum(coord_nod2D(1, elnodes))/3.
         y=sum(coord_nod2D(2, elnodes))/3.
      end if
      vol=voltriangle(el)

      do i=1, reg_ny
         if (reg_lat(i) > y) then
            pos=i
            exit
         end if
      end do
      do lev=1, elvls(el)-1
         moc(pos, lev)=moc(pos, lev)+vol*sum(w(lev, elnodes))/3.*1.e-6
      end do
   end do

   do i=2, reg_ny
      moc(i, :)=moc(i, :)+moc(i-1, :)
   end do
   call netcdf_write(fileout, year==year_start)      
   end do
   io=nf_close(ncid)   
   deallocate(w, arr, moc)
end subroutine make_diag_moc_w
! ============================================================
SUBROUTINE netcdf_write(filename, newold)
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
  LOGICAL,        INTENT(IN)  :: newold
  INTEGER       :: i,j,k,s,n
  INTEGER       :: lat_id, dep_id, tim_id
  INTEGER       :: var
  INTEGER, save :: month=0
  INTEGER, DIMENSION(3)       :: dimarray
  INTEGER, DIMENSION(50)      :: stat    !  status array

  if (newold) then     
  stat(1) = NF_CREATE(TRIM(filename),0,fileid); s=s+1  
  if (stat(1).ne.NF_NOERR) then
     write(*,*) 'NetCDF error while opening file - in cpl_str_nc_init'
     write(*,*) trim(filename)
  end if

  s=1
  stat(s) = NF_DEF_DIM(fileid,'lats',  reg_ny,         dimarray(1)); s=s+1
  stat(s) = NF_DEF_DIM(fileid,'deps',  max_num_layers, dimarray(2)); s=s+1
  stat(s) = NF_DEF_DIM(fileid,'time',  NF_UNLIMITED,   dimarray(3)); s=s+1    

  do i=1,s-1
     if (stat(i).ne.NF_NOERR) write(*,*) &
        'NetCDF error in dimension definitions, no.',i
  end do

  s=1
  stat(s) = NF_DEF_VAR(fileid,'lats',    NF_FLOAT,1,dimarray(1), lat_id); s=s+1
  stat(s) = NF_DEF_VAR(fileid,'deps',    NF_FLOAT,1,dimarray(2), dep_id); s=s+1
  stat(s) = NF_DEF_VAR(fileid,'time',    NF_INT,  0,dimarray(3), tim_id); s=s+1  
  
  stat(s)=  NF_DEF_VAR(fileid, trim(svarid_out(1)), NF_DOUBLE, 3, dimarray(1:3), var); s=s+1
  do i=1,s-1
     if (stat(i).ne.NF_NOERR) write(*,*) &
        'NetCDF error in str. variable definition, no.',i
  end do
  stat(1)=NF_ENDDEF(fileid)
  
  s=1;
  stat(s)=NF_PUT_VAR_REAL(fileid,lat_id,real(reg_lat/rad, 4)); s=s+1  
  stat(s)=NF_PUT_VAR_REAL(fileid,dep_id,real(depths, 4)); s=s+1  
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

  month=month+1
  stat(s) = nf_put_var_int (fileid, tim_id, month); s=s+1
  n=1
  stat(s) = NF_INQ_VARID   (fileid, trim(svarid_out(1)), var);s=s+1
  stat(s)=NF_PUT_VARA_DOUBLE(fileid,var, (/1, 1, month/), (/reg_ny, max_num_layers, 1/), moc); s=s+1

  do i=1,s-1
     if (stat(i).ne.NF_NOERR) write(*,*) &
        'NetCDF error in writing variable, no.', i
  end do  
  stat(1)=NF_CLOSE(fileid); s=s+1    
end SUBROUTINE  netcdf_write
end module diag_moc_w
!
!----------------------------------------------------------------------------
!
