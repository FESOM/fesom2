subroutine init_diag

use o_mesh
use g_config
use g_clock
use g_parsup
implicit none
#include "netcdf.inc"

write(*,*) 'Init new diag files' 
if (yearnew==yearold) return

filename=trim(ResultPath)//runid//'.'//cyearnew//'.diag.nc'

! create a file

status = nf_create(filename, nf_clobber, ncid)
if (status.ne.nf_noerr) call handle_err(status)

status = nf_def_dim(ncid, 'nodes_2d', nod2d, dimid_2d)
if (status .ne. nf_noerr) call handle_err(status)
status = nf_def_dim(ncid, 'elem_2d', elem2d, dimid_2de)  !!!! 
if (status .ne. nf_noerr) call handle_err(status)
status = nf_def_dim(ncid, 'nl_1', nl-1, dimid_nl1)
if (status .ne. nf_noerr) call handle_err(status)
status = nf_def_dim(ncid, 'nl', nl, dimid_nl)
if (status .ne. nf_noerr) call handle_err(status)       !!!!
status = nf_def_dim(ncid, 'nl', 3, dimid_gm)
if (status .ne. nf_noerr) call handle_err(status)       !!!!

status = nf_def_dim(ncid, 'T', NF_UNLIMITED, dimid_rec)
if (status .ne. nf_noerr) call handle_err(status)

! Define the time and iteration variables
status = nf_def_var(ncid, 'time', NF_DOUBLE, 1, dimid_rec, time_varid)
if (status .ne. nf_noerr) call handle_err(status)
status = nf_def_var(ncid, 'iter', NF_INT, 1, dimid_rec, iter_varid)
if (status .ne. nf_noerr) call handle_err(status)

! Define the netCDF variables for 2D fields.
! In Fortran, the unlimited dimension must come
! last on the list of dimids.
dimids(1) = dimid_2d
dimids(2) = dimid_rec


dimid4(1) = dimid_gm
dimid4(2) = dimid_nl1 !mid levels
dimid4(3) = dimid_2d
dimid4(4) = dimid_rec

status = nf_def_var(ncid, 'Redi_GM', NF_DOUBLE, 4, dimid4, gm_varid)
if (status .ne. nf_noerr) call handle_err(status)


  longname='model time'
  status = nf_PUT_ATT_TEXT(ncid, time_varid, 'long_name', len_trim(longname), trim(longname))
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_put_att_text(ncid, time_varid, 'units', 1, 's')
  if (status .ne. nf_noerr) call handle_err(status)
  longname='iteration_count'
  status = nf_PUT_ATT_TEXT(ncid, iter_varid, 'long_name', len_trim(longname), trim(longname))
  if (status .ne. nf_noerr) call handle_err(status)

  longname='Redi_GM matrix elements'
  status = nf_put_att_text(ncid, gm_varid, 'description', len_trim(longname), trim(longname)) 
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_put_att_text(ncid, gm_varid, 'units', 1, 'm^2/sec')
  if (status .ne. nf_noerr) call handle_err(status)
end subroutine init_diag


subroutine write_diag(istep)
use o_arrays
use o_mesh
use i_arrays
use g_config
use g_clock
use g_parsup
use g_comm_auto
implicit  none

allocate(aux2(nod2D), aux3(nl-1,elem2D),aux4(3,nl-1,nod2D)) 

sec_in_year=dt*istep
filename=trim(ResultPath)//runid//'.'//cyearnew//'.diag.nc'

status = nf_open(filename, nf_write, ncid)
if (status .ne. nf_noerr) call handle_err(status)

! inquire variable id
status=nf_inq_varid(ncid, 'time', time_varid)
if (status .ne. nf_noerr) call handle_err(status)
status=nf_inq_varid(ncid, 'iter', iter_varid)
if (status .ne. nf_noerr) call handle_err(status)

status=nf_inq_varid(ncid, 'Redi_GM', gm_varid)
if (status .ne. nf_noerr) call handle_err(status)

! write variables
! time and iteration
status=nf_put_vara_double(ncid, time_varid, save_count, 1, sec_in_year)
if (status .ne. nf_noerr) call handle_err(status)
status=nf_put_vara_int(ncid, iter_varid, save_count, 1, istep)
if (status .ne. nf_noerr) call handle_err(status)


call broadcast_nod(Kd(1,:,:),aux4(1,:,:))
call broadcast_nod(Kd(2,:,:),aux4(2,:,:))
call broadcast_nod(Kd(3,:,:),aux4(3,:,:))

start4=(/1,1,1,save_count/)
count4=(/3,nl, nod2D, 1/)
status=nf_put_vara_double(ncid, gm_varid, start4, count4, aux4) 
if (status .ne. nf_noerr) call handle_err(status)


end subroutine write_diag(istep)
