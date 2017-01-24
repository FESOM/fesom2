module grid_info
   use o_param
   use o_elements
   use o_mesh
   use g_config
   use g_rotate_grid
   use g_oce_2_reg
   implicit none
   
   integer                   :: i, j, k
   real(kind=8)              :: x, y
   real(kind=8)              :: scos, dx, dy
   real(kind=8), allocatable :: area2(:,:),   deps3(:,:,:), area3(:,:,:)
   integer, dimension(:,:),      allocatable   :: oce_2_reg_mask
   integer, dimension(:,:,:),    allocatable   :: oce_2_reg_mask3
   real(kind=8), dimension(:),   allocatable   :: oce_fld
   real(kind=8), dimension(:,:), allocatable   :: reg_fld

  PRIVATE
  PUBLIC :: make_grid_info
CONTAINS   

subroutine make_grid_info
   implicit none
   oce_2_reg_mask =0
   oce_2_reg_mask3=0

   allocate(area2(reg_nx, reg_ny), area3(reg_nx, reg_ny, max_num_layers-1), deps3(reg_nx, reg_ny, max_num_layers-1))
   allocate(oce_fld(nod2D))
   allocate(reg_fld(reg_nx, reg_ny))
   allocate(oce_2_reg_mask  (reg_nx, reg_ny))
   allocate(oce_2_reg_mask3 (reg_nx, reg_ny, max_num_layers-1))
   area2=0.
   area3=0.
   deps3=0.
! fill 2D mask
   oce_fld=1   
   call do_oce_2_reg(oce_fld, reg_fld, 1)   
   where(reg_fld>0.5)
	oce_2_reg_mask=1
   elsewhere
	oce_2_reg_mask=0
   end where

! fill 3D mask
   do k=2, max_num_layers
      where(nlvls >= k)
           oce_fld=1
      elsewhere
           oce_fld=0
      end where
      call do_oce_2_reg(oce_fld, reg_fld, 1)
      where (reg_fld > .9)      	    
	    oce_2_reg_mask3(:,:,k-1)=1
      end where
   end do

!compute 2D area
   do i=1, reg_nx
   do j=1, reg_ny
      if (oce_2_reg_mask(i,j)==0) CYCLE
      scos=cos(reg_lat(j))
      dx=RegDx*r_earth
      dy=RegDy*r_earth*scos
      area2(i,j)=dx*dy
   end do
   end do
      
!compute 3D depths
   do k=2, max_num_layers
      deps3(:,:,k-1)=0.5*abs(depths(k)+depths(k-1))
   end do   
   deps3=deps3*real(oce_2_reg_mask3, 8)

!compute 3D area
   do k=2, max_num_layers
      area3(:,:,k-1)=area2*abs(depths(k)-depths(k-1))
   end do   

   area3=area3*real(oce_2_reg_mask3, 8)   
   call netcdf_write
   deallocate(deps3, area3, area2)
   deallocate(reg_fld, oce_fld)
   deallocate(oce_2_reg_mask, oce_2_reg_mask3)
end subroutine make_grid_info
! ============================================================
SUBROUTINE netcdf_write
  USE o_param
  USE o_mesh
  USE o_elements
  use g_oce_2_reg  
  IMPLICIT NONE

#include "netcdf.inc" 

  !  IDs
  INTEGER :: fileid
  ! dimensions
 
  INTEGER :: i,j,k, s
  INTEGER :: lon_id, lat_id, dep_id
  INTEGER :: var(3)

  INTEGER, DIMENSION(3)       :: dimarray
  INTEGER, DIMENSION(reg_nx)  :: NUM_LON
  INTEGER, DIMENSION(reg_ny)  :: NUM_LAT
  INTEGER, DIMENSION(50)      :: stat    !  status array

  stat(1) = NF_CREATE(trim(outpath)//'geometry.nc',0,fileid); s=s+1
     
  if (stat(1).ne.NF_NOERR) then
     write(*,*) 'NetCDF error while opening file - in cpl_str_nc_init'
  end if

  s=1
  stat(s) = NF_DEF_DIM(fileid,'lons', reg_nx,          dimarray(1)); s=s+1
  stat(s) = NF_DEF_DIM(fileid,'lats', reg_ny,          dimarray(2)); s=s+1
  stat(s) = NF_DEF_DIM(fileid,'deps', max_num_layers-1,dimarray(3)); s=s+1  
  
  do i=1,s-1
     if (stat(i).ne.NF_NOERR) write(*,*) &
        'NetCDF error in dimension definitions, no.',i
  end do

  s=1

  stat(s) = NF_DEF_VAR(fileid,'lons',NF_FLOAT,1,dimarray(1),lon_id); s=s+1
  stat(s) = NF_DEF_VAR(fileid,'lats',NF_FLOAT,1,dimarray(2),lat_id); s=s+1
  stat(s) = NF_DEF_VAR(fileid,'deps',NF_FLOAT,1,dimarray(3),dep_id); s=s+1
     
  stat(s)= NF_DEF_VAR(fileid,'xy_area',NF_DOUBLE, 2, dimarray(1:2), var(1)); s=s+1
  stat(s)= NF_DEF_VAR(fileid,'xyz_dep',NF_DOUBLE, 3, dimarray(1:3), var(2)); s=s+1
  stat(s)= NF_DEF_VAR(fileid,'xyz_vol',NF_DOUBLE, 3, dimarray(1:3), var(3)); s=s+1   
  do i=1,s-1
     if (stat(i).ne.NF_NOERR) write(*,*) &
        'NetCDF error in str. variable definition, no.',i
  end do
  stat(1)=NF_ENDDEF(fileid)

  !------  Writing Variables -----
  stat(s)=NF_PUT_VAR_REAL(fileid,lon_id,real(reg_lon/rad, 4)); s=s+1
  stat(s)=NF_PUT_VAR_REAL(fileid,lat_id,real(reg_lat/rad, 4)); s=s+1  
  stat(s)=NF_PUT_VAR_REAL(fileid,dep_id,real(0.5*(depths(1:max_num_layers-1)+depths(2:max_num_layers)), 4)); s=s+1
  stat(s)=NF_PUT_VAR_DOUBLE(fileid,var(1), area2); s=s+1
  stat(s)=NF_PUT_VAR_DOUBLE(fileid,var(2), deps3); s=s+1  
  stat(s)=NF_PUT_VAR_DOUBLE(fileid,var(3), area3); s=s+1    
  do i=1,s-1
     if (stat(i).ne.NF_NOERR) write(*,*) &
        'NetCDF error in writing variable, no.', i
  end do  
  stat(1)=NF_CLOSE(fileid)
end SUBROUTINE  netcdf_write
end module grid_info
! ============================================================
