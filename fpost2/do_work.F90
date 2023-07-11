PROGRAM MAIN
use o_param
use o_elements
use o_mesh
use g_config
use g_rotate_grid
use g_oce_2_reg
use grid_info
use diag_moc_w
use diag_uv_norm3
use diag_uv_curl3
use diag_TS3

IMPLICIT NONE
#include "netcdf.inc"

INTEGER                                    :: year, io, ncid
real, allocatable, dimension(:)            :: temp, salt, ssh
character(len=MAX_PATH)                         :: nmlfile

nmlfile ='namelist.interp'    ! name of general configuration namelist file
open (20,file=nmlfile)
read (20,NML=config)
read (20,NML=mask)
read (20,NML=todo)
read (20,NML=fesom_mesh)
read (20,NML=regular_mesh)
close (20)

call ocean_mesh_setup
call def_exchange_mesh
write(*,*) 'elem2D=', elem2D
write(*,*) 'Total ocean area is: ', sum(voltriangle), ' m^2'
!run it only once
if (do_mesh) then
   call build_oce_2_reg
   call save_oce_2_reg
! needs a permanent update of all entried which require interpolation
elseif (do_UVnorm .or. do_UVcurl .or. do_TS3) then
   call load_oce_2_reg
end if


!call build_oce_2_reg_el
!call save_oce_2_reg_el
!call load_oce_2_reg_el

!call make_grid_info
if (do_UVnorm) call make_diag_uv_norm3
if (do_UVcurl) call make_diag_uv_curl3
if (do_MOC)    call make_diag_moc_w
if (do_TS3)    call make_diag_TS3
END PROGRAM MAIN
